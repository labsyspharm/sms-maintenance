# Run on O2 cluster

library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(batchtools)
library(lspcheminf)
library(sys)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

dir.create(dir_release, showWarnings = FALSE)

# Loading files ----------------------------------------------------------------
###############################################################################T

inputs <- list(
  chembl_compounds = synPluck(syn_release, "raw_data", "chembl_compounds_raw.rds"),
  hmsl_compounds = synPluck(syn_release, "raw_data", "hmsl_compounds_raw.rds"),
  emolecules_compounds = synPluck(syn_release, "id_mapping", "emolecules", "emolecules_compounds.csv.gz")
)

raw_compounds <- tribble(
  ~source, ~data,
  "chembl", inputs[["chembl_compounds"]] %>%
    syn() %>%
    read_rds() %>%
    distinct(id = chembl_id, inchi = standard_inchi),
  "hsml", inputs[["hmsl_compounds"]] %>%
    syn() %>%
    read_rds() %>%
    distinct(id = hms_id, inchi),
  "emolecules", inputs[["emolecules_compounds"]] %>%
    syn() %>%
    read_csv() %>%
    distinct(id = parent_id, inchi)
)

# Canonicalize compounds -------------------------------------------------------
###############################################################################T

wd <- file.path("/n", "scratch3", "users", "c", "ch305", "canonicalize_compounds")
dir.create(wd)

all_inchis <- raw_compounds %>%
  # mutate(
  #   data = map2(source, data, ~mutate(.y, id = paste(.x, id, sep = "_")))
  # ) %>%
  pull(data) %>%
  map(pull, inchi) %>%
  reduce(union)

all_inchis_dfs <- tibble(compound = all_inchis) %>%
  chunk_df(ceiling(nrow(.)/10000), seed = 1) %>%
  enframe(name = "chunk", value = "data") %>%
  mutate(
    input_file = file.path(wd, paste0("compounds_", chunk, ".csv")),
    data = map(data, ~mutate(.x, row = 1:n()))
  )

pwalk(
  all_inchis_dfs,
  function(data, input_file, ...)
    write_csv(data, input_file)
)

# Set up jobs
reg <- makeRegistry(
  file.dir = file.path(wd, paste0("registry_", gsub(" ", "_", Sys.time()))),
  seed = 1
)
#reg$cluster.functions <- makeClusterFunctionsSlurm(template = "slurm-simple")

run_canonicalize_job <- function(input_file, timeout = 5*60) {
  library(sys)
  script_path <- tempfile("shell_script_", fileext = ".sh")
  script <- c(
    "unset PYTHONPATH",
    "unset PYTHONHOME",
    "source ~/miniconda3/etc/profile.d/conda.sh",
    "conda activate lspcheminf_env",
    "which -a python",
    "env",
    paste(
      "tautomers canonicalize",
      "--standardizer chembl-parent",
      "--timeout", timeout,
      input_file,
      file.path(
        dirname(input_file),
        paste0(tools::file_path_sans_ext(basename(input_file)), "_canonical.csv")
      )
    )
  )
  message("SCRIPT:\n", paste(paste0("# ", script), collapse = "\n"), "\n")
  writeLines(
    script,
    con = script_path
  )
  out <- sys::exec_wait("bash", script_path)
  if (out != 0)
    stop("Canonicalization failed with error code ", out)
  out
}

batchMap(
  fun = run_canonicalize_job,
  input_file = all_inchis_dfs[["input_file"]],
  # Limiting canonicalization to 5 min per compound
  more.args = list(timeout = 5*60)
)

job_table <- findJobs() %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table,
  resources = list(
    memory = "4gb",
    ncpus = 1L,
    partition = "short",
    # Approximately 58 compounds processed per minute
    # With 10,000 compounds per batch, should take 172 min, use 250 min
    # Turns out uses much less time, 45 min enough
    walltime = 45*60,
    chunks.as.arrayjobs = TRUE,
    # For some reason these nodes fail to execute R because of an "illegal instruction"
    exclude = "compute-f-17-[09-25]"
  )
)

waitForJobs()

canonical_inchis <- all_inchis_dfs %>%
  rowwise() %>%
  transmute(
    chunk,
    data = {
      canonical <- read_csv(
        file.path(wd, paste0("compounds_", chunk, "_canonical.csv")),
        col_types = "ci"
      ) %>%
        transmute(row = row + 1, canonical_inchi = inchi)
      full_join(
        data %>%
          rename(raw_inchi = compound),
        canonical,
        by = "row"
      ) %>%
        select(-row) %>%
        list()
    }
  ) %>%
  pull(data) %>%
  bind_rows()

write_csv(
  canonical_inchis,
  file.path(dir_release, "canonical_inchis_raw.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Canonicalize compounds",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/03_canonicalize.R"
)

syn_id_mapping <- Folder("canonicalization", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "canonical_inchis_raw.csv.gz")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = activity)

