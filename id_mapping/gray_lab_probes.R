library(tidyverse)
library(here)
library(rvest)
library(httr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")


dir.create(
  here("gray_lab")
)

probes_page <- xml2::read_html(
  "https://graylab.dana-farber.org/probes.html"
)

compound_names <- html_nodes(
  probes_page,
  css = "u strong"
) %>%
  html_text(trim = TRUE) %>%
  str_replace(fixed("Probe-"), "")

dl_links <- html_nodes(
  probes_page,
  css = ".wsite-multicol-col div div div a"
) %>%
  html_attr("href") %>% {
    paste0("https://graylab.dana-farber.org", .)
  }

pdf_files <- here(
  "gray_lab",
  paste(
    seq_along(compound_names),
    paste0(compound_names, ".pdf"),
    sep = "_"
  )
)

download.file(
  dl_links,
  pdf_files
)


activity <- Activity(
  name = "Download Gray lab probe data",
  used = c(
    "https://graylab.dana-farber.org/probes.html"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/gray_lab_probes.R"
)

syn_folder <- synMkdir(syn_release, "raw_data", "gray_lab")

c(
  pdf_files
) %>%
  synStoreMany(parentId = syn_folder, activity = activity)
