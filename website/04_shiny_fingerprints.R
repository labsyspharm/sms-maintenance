library(tidyverse)
library(here)
library(fst)
library(synapser)
library(synExtra)
library(morgancpp)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

syn_tables <- "syn20981852"

fingerprints <- syn("syn21042105") %>%
  read_rds()

# Fingerprints -----------------------------------------------------------------
###############################################################################T

fingerprint_tables <- fingerprints %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        arrange(lspci_id) %>%
        filter(fp_name == "morgan_normal") %>%
        {set_names(.[["fingerprint"]], .[["lspci_id"]])}
    )
  )

activity <- Activity(
  name = "Convert fingerprints to binary files used by morgancpp",
  used = c(
    "syn21042105"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/04_shiny_fingerprints.R"
)

pwalk(
  fingerprint_tables,
  function(data, fp_name, ...) {
    fps <- MorganFPS$new(data)
    fps$save_file(
      file.path(dir_release, paste0("shiny_fingerprints_", fp_name, ".bin"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_fingerprints_", fp_name, ".bin")),
      parent = syn_parent,
      name = "shiny_fingerprints.bin"
    ) %>%
      synStore(activity = activity)
  }
)
