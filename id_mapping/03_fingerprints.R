library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)
library(batchtools)
library(processx)
library(data.table)
library(qs)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v29"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

inputs <- c(
  inchis = synPluck(syn_release, "canonicalization", "canonical_inchi_ids.csv.gz")
)

# Load compound data -----------------------------------------------------------
###############################################################################T

compounds <- inputs[["inchis"]] %>%
  syn() %>%
  fread()

# Create molecular fingerprints ------------------------------------------------
###############################################################################T

fingerprinting_fun <- function(compound_file, output_file, fingerprint_type, fingerprint_args, port = NULL) {
  library(processx)
  library(tidyverse)
  library(lspcheminf)

  if (is.null(port))
    port <- httpuv::randomPort()

  lspcheminf_script <- paste0(
    "source ~/.bashrc
    conda activate lspcheminf_env
    gunicorn --workers=1 -b 127.0.0.1:", port, " -t 60000 lspcheminf"
  )

  # Because of freak bug where the entire process hangs when stdout / stderr
  # is captured disablign it. I can't see anythign strange about the output
  # when I runn the process manually
  lspcheminf_p <- process$new(
    # "bash", c("-c", lspcheminf_script), stderr = "|", stdout = "|", env = character()
    "bash", c("-c", lspcheminf_script), env = character()
  )

  Sys.sleep(5)

  # message("lspcheminf launch output", lspcheminf_p$read_error_lines(), lspcheminf_p$read_output_lines())

  fingerprint_df <- tryCatch(
    {
      compound_df <- read_csv(compound_file)
      safely(calculate_fingerprints)(
        with(
          compound_df,
          set_names(canonical_inchi, inchi_id)
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
  "morgan_normal", "morgan", list(useChirality = FALSE),
  "topological_normal", "topological", NULL
)

cmpd_chunks <- compounds %>%
  # slice(1:100) %>%
  chunk_df(300, seed = 1) %>%
  enframe("index", "compounds") %>%
  mutate(
    compound_file = file.path(wd, paste0("compounds_", index, ".csv"))
  )

pwalk(
  cmpd_chunks,
  function(compounds, compound_file, ...)
    fwrite(compounds, compound_file)
)

cmpd_fingerprint_input <- cmpd_chunks %>%
  select(index, compound_file) %>%
  crossing(fingerprinting_args) %>%
  mutate(
    index_out = seq_len(n()),
    output_file = file.path(wd, paste0("compound_fingerprints_", index_out, ".csv"))
  )

qsave(
  cmpd_fingerprint_input,
  file.path(wd, "cmpd_fingerprint_input.qs"),
  preset = "balanced"
)

# cmpd_fingerprint_input <- qread(file.path(wd, "cmpd_fingerprint_input.qs"))

batchMap(
  fun = fingerprinting_fun,
  compound_file = cmpd_fingerprint_input[["compound_file"]],
  output_file = cmpd_fingerprint_input[["output_file"]],
  fingerprint_type = cmpd_fingerprint_input[["fp_type"]],
  fingerprint_args = cmpd_fingerprint_input[["fp_args"]]
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
  job_table[findErrors()],
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

cmpd_fingerprints_chunks <- copy(setDT(cmpd_fingerprint_input))[
  ,
  data := lapply(output_file, fread, colClasses = c("character", "numeric"), sep = ",", showProgress = FALSE)
]

all(map_int(cmpd_fingerprints_chunks$data, nrow) > 1)

cmpd_fingerprints_nested <- cmpd_fingerprints_chunks[
  ,
  .(data = list(rbindlist(data, use.names = TRUE, fill = TRUE))),
  by = .(fp_name, fp_type)
]

qsave(cmpd_fingerprints_nested, file.path(dir_release, "all_compounds_fingerprints_nested.qs"))

cmpd_fingerprints <- cmpd_fingerprints_nested %>%
  unnest(data)

fwrite(cmpd_fingerprints, file.path(dir_release, "all_compounds_fingerprints.csv.gz"))

# cmpd_fingerprints <- fread(file.path(dir_release, "all_compounds_fingerprints.csv.gz"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

fingerprint_activity <- Activity(
  name = "Calculate molecular fingerprints",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/03_fingerprints.R"
)

fp_folder <- synMkdir(syn_release, "fingerprints")

c(
  file.path(dir_release, "all_compounds_fingerprints_nested.qs"),
  file.path(dir_release, "all_compounds_fingerprints.csv.gz")
) %>%
  synStoreMany(parent = fp_folder, activity = fingerprint_activity, forceVersion = FALSE)
