library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)
library(batchtools)
library(processx)
library(data.table)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

inputs <- c(
  inchis = synPluck(syn_release, "id_mapping", "canonical_inchis_raw.csv.gz")
)

# Load compound data -----------------------------------------------------------
###############################################################################T

compounds <- inputs[["inchis"]] %>%
  syn() %>%
  read_csv()

compounds_all <- compounds %>%
  # Compute chemical formula
  drop_na(raw_inchi) %>%
  mutate(
    formula = canonical_inchi %>%
      str_split_fixed(fixed("/"), 3) %>%
      {.[, 2]},
    formula_vector = formula %>%
      str_match_all("([A-Z][a-z]?)([0-9]*)") %>%
      map(~set_names(replace_na(as.integer(.x[, 3]), 1L), .x[, 2]))
  )

compounds_selected <- compounds_all %>%
  mutate(
    inchi = case_when(
      is.na(canonical_inchi) ~ raw_inchi,
      # Only use canonicalized InChI for compounds that are "drug-like" with
      # at least 3 carbons
      map_int(formula_vector, ~if ("C" %in% names(.x)) .x[["C"]] else 0L) >= 3 ~ canonical_inchi,
      TRUE ~ raw_inchi
    )
  )

write_rds(
  compounds_selected,
  file.path(dir_release, "all_compounds_canonical.rds"),
  compress = "gz"
)
# compounds_selected <- read_rds(
#   file.path(dir_release, "all_compounds_canonical.rds")
# )

# Create molecular fingerprints ------------------------------------------------
###############################################################################T

fingerprinting_fun <- function(compound_file, output_file, fingerprint_type, fingerprint_args, port = 8000) {
  library(processx)
  library(tidyverse)
  library(lspcheminf)

  lspcheminf_script <- paste0(
    "unset PYTHONPATH
    unset PYTHONHOME
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate lspcheminf_env
    gunicorn --workers=1 -b 127.0.0.1:", port, " -t 60000 lspcheminf"
  )

  lspcheminf_p <- process$new(
    "bash", c("-c", lspcheminf_script), stderr = "|", stdout = "|"
  )

  Sys.sleep(5)

  message("lspcheminf launch output", lspcheminf_p$read_error_lines(), lspcheminf_p$read_output_lines())

  fingerprint_df <- tryCatch(
    {
      compound_df <- read_csv(compound_file)
      safely(calculate_fingerprints)(
        with(
          compound_df,
          set_names(compound, id)
        ),
        fingerprint_type = fingerprint_type,
        fingerprint_args = fingerprint_args,
        url = paste0("http://127.0.0.1:", port)
      )
    },
    finally = {
      message("Cleaning up")
      lspcheminf_p$kill_tree()
    }
  )

  message("lspcheminf fingerprinting output", fingerprint_df[["error"]])

  if (!is.null(fingerprint_df[["result"]]))
    write_csv(
      fingerprint_df[["result"]],
      output_file
    )
  else
    stop("Unsuccesful")
}

wd <- file.path("/n", "scratch3", "users", "c", "ch305", "fingerprints")
dir.create(wd, showWarnings = FALSE)

reg <- makeRegistry(
  file.dir = file.path(wd, paste0("registry_fingerprints_", gsub(" ", "_", Sys.time()))),
  seed = 1
)

fingerprinting_args <- tribble(
  ~fp_name, ~fp_type, ~fp_args,
  "morgan_chiral", "morgan", list(useChirality = TRUE),
  "morgan_normal", "morgan", list(useChirality = FALSE),
  "topological_normal", "topological", NULL
)

cmpd_ids <- compounds_selected %>%
  distinct(compound = inchi) %>%
  # Some HMSL compounds don't report inchi or smiles or anything, should have removed
  # them earlier, removing them here now
  drop_na() %>%
  mutate(id = seq_len(n()))

cmpd_chunks <- cmpd_ids %>%
  # slice(1:100) %>%
  chunk_df(300, seed = 1) %>%
  enframe("index", "compounds") %>%
  mutate(
    compound_file = file.path(wd, paste0("compounds_", index, ".csv"))
  )

pwalk(
  cmpd_chunks,
  function(compounds, compound_file, ...)
    write_csv(compounds, compound_file)
)

cmpd_fingerprint_input <- cmpd_chunks %>%
  select(index, compound_file) %>%
  crossing(fingerprinting_args) %>%
  mutate(
    index_out = seq_len(n()),
    port = 9000 + index_out,
    output_file = file.path(wd, paste0("compound_fingerprints_", index_out, ".csv"))
  )

write_rds(
  cmpd_fingerprint_input,
  file.path(wd, "cmpd_fingerprint_input.rds")
)

batchMap(
  fun = fingerprinting_fun,
  compound_file = cmpd_fingerprint_input[["compound_file"]],
  output_file = cmpd_fingerprint_input[["output_file"]],
  fingerprint_type = cmpd_fingerprint_input[["fp_type"]],
  fingerprint_args = cmpd_fingerprint_input[["fp_args"]],
  port = cmpd_fingerprint_input[["port"]]
)


# with(
#   cmpd_fingerprint_input,
#   fingerprinting_fun(
#     compound_file[[1]],
#     fp_type[[1]],
#     fp_args[[1]],
#     port[[1]]
#   )
# )

job_table <- findJobs() %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table[findExpired()],
  resources = list(
    memory = "2gb",
    ncpus = 1L,
    partition = "short",
    walltime = 5*60,
    chunks.as.arrayjobs = TRUE,
    # For some reason these nodes fail to execute R because of an "illegal instruction"
    exclude = "compute-f-17-[09-25]"
  )
)

waitForJobs()


cmpd_fingerprints_chunks <- cmpd_fingerprint_input %>%
  select(
    starts_with("fp_"), output_file
  ) %>%
  group_by(
    across(starts_with("fp_"))
  ) %>%
  summarize(
    data = output_file %>%
      map(fread, colClasses = c("character", "numeric"), sep = ",", showProgress = FALSE) %>%
      rbindlist() %>%
      list(),
    .groups = "drop"
  )


setDT(cmpd_ids)

cmpd_fingerprints  <- cmpd_fingerprints_chunks %>%
  rowwise() %>%
  mutate(
    data = data[
      cmpd_ids, on = c("names" = "id"), nomatch = NULL
    ] %>%
      list()
  ) %>%
  ungroup()

write_rds(
  cmpd_fingerprints,
  file.path(dir_release, "all_compounds_fingerprints.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

fingerprint_activity <- Activity(
  name = "Calculate molecular fingerprints",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/03_fingerprints.R"
)

fp_folder <- Folder("fingerprints", syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_fingerprints.rds")
) %>%
  synStoreMany(parent = fp_folder, activity = fingerprint_activity)
