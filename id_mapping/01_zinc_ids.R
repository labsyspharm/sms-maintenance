library(tidyverse)
library(synapser)
library(synExtra)
library(furrr)
library(httr)
library(future)
library(rvest)
library(polite)
library(here)
library(data.table)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

chembl_zinc_mapping <- synPluck(syn_release, "id_mapping", "unichem", "chembl_zinc_mapping.csv.gz") %>%
  syn() %>%
  read_csv()

zinc_catalogs <- synPluck(syn_release, "id_mapping", "zinc", "zinc_catalogs.csv") %>%
  syn() %>%
  read_csv()

dir_vendor <- file.path(dir_release, "vendor")
dir.create(dir_vendor, showWarnings = FALSE)

# Fetch vendor info from ZINC --------------------------------------------------
###############################################################################T

# Get catalog URLs
catalog_urls <- tribble(
  ~class_id, ~class,
  "50", "in_stock",
  "40", "in_stock",
  "30", "in_stock",
  "20", "on_demand",
  "10", "boutique",
  # "7", "one_step",
  "3", "unclear"
) %>%
  mutate(
    vendor = map(
      class_id,
      ~xml2::read_html(paste0("http://files.docking.org/catalogs/", .x)) %>%
        rvest::html_node("table") %>%
        rvest::html_table(fill = TRUE) %>%
        magrittr::extract(-1) %>%
        mutate(
          Name = str_trim(Name) %>%
            str_replace("/$", "")
        ) %>%
        filter(Size == "-", Name != "", `Last modified` != "") %>%
        pull(Name)
      # rvest::html_nodes("td[valign='top'] + td a") %>%
      # magrittr::extract(-1) %>%
      # rvest::html_text(trim = TRUE) %>%
      # str_replace("/$", "")
    )
  ) %>%
  unchop(vendor)

catalog_files <- catalog_urls %>%
  crossing(
    tribble(
      ~type, ~url_template,
      "info", "https://files.docking.org/catalogs/{class_id}/{vendor}/{vendor}.info.txt",
      "source", "https://files.docking.org/catalogs/source/{vendor}.src.txt"
    )
  ) %>%
  group_by(type) %>%
  mutate(
    url = glue::glue_data(cur_data(), url_template[[1]])
  ) %>%
  ungroup()

# Download vendor libraries

download_zinc <- function(vendor, url, type, ...) {
  message(vendor)
  path <- file.path(dir_vendor, paste0(type, "_", vendor, ".txt"))
  remote_path <- url
  for (ext in c(".gz", "")) {
    p <- paste0(path, ext)
    if (file.exists(p) && file.size(p) > 5)
      return(p)
    rp <- paste0(remote_path, ext)
    res <- try(
      download.file(rp, p),
      silent = TRUE
    )
    if (class(res) != "try-error")
      return(p)
  }
  return(NA_character_)
}

catalog_files <- catalog_files %>%
  mutate(
    path = pmap_chr(
      .,
      download_zinc
    )
  )


# Manually downloading some files that don't correspond to pattern
manual_files <- tribble(
  ~vendor, ~type, ~url,
  "cdive", "source", "https://files.docking.org/catalogs/source/cdiv.src.txt"
)


catalog_files <- catalog_files %>%
  setDT() %>%
  {
    .[
      manual_files %>%
        mutate(
          path = pmap_chr(
            .,
            download_zinc
          )
        ) %>%
        setDT(),
      on = c("vendor", "type"),
      c("url", "path") := list(i.url, i.path)
    ]
  } %>%
  setDF()


# Remove "cdive" "cdivp", doesn't contain any

vendor_libraries <- catalog_files %>%
  filter(type == "info") %>%
  drop_na(path) %>%
  mutate(
    data = pmap(
      .,
      function(vendor, path, ...) {
        message(vendor)
        tryCatch(
          fread(
            path,
            header = FALSE,
            colClasses = rep_len("character", 5),
            col.names = c("vendor_id", "zinc_id", "inchi_key", "tranche", "notes")
          ),
          error = function(e) NULL
        )
      }
    )
  )

vendor_smiles <- catalog_files %>%
  filter(type == "source") %>%
  drop_na(path) %>%
  mutate(
    data = pmap(
      .,
      function(vendor, path, ...) {
        message(vendor)
        tryCatch(
          fread(
            path,
            header = FALSE,
            colClasses = rep_len("character", 2),
            col.names = c("smiles", "vendor_id")
          ),
          error = function(e) NULL
        )
      }
    )
  )

vendor_libraries_all <- vendor_libraries %>%
  select(-url_template, -url, -path) %>%
  unnest(data)


vendor_smiles_all <- vendor_smiles %>%
  select(vendor, data) %>%
  unnest(data)

vendor_smiles_zinc <- vendor_smiles_all %>%
  inner_join(
    vendor_libraries_all %>%
      select(vendor, vendor_id, zinc_id),
    by = c("vendor", "vendor_id")
  )

write_csv(
  file.path(dir_vendor, "zinc_smiles.csv.gz"),
  vendor_smiles_all
)

write_csv(
  file.path(dir_vendor, "zinc_compounds")
)

#
# # Mapping from Zinc IDs to Chembl provided by Zinc directly
# download.file("http://files.docking.org/catalogs/1/chembl23/chembl23.codemap.txt.gz", file.path(dir_vendor, "codemap_chembl23.txt.gz"))
# chembl_codemap <- read_delim(file.path(dir_vendor, "codemap_chembl23.txt.gz"), delim = " ", col_names = c("zinc_id", "chembl_id"), col_types = "cc")
#
# zinc_mapping <- chembl_zinc_mapping %>%
#   select(-uci) %>%
#   bind_rows(chembl_codemap) %>%
#   distinct() %>%
#   setDT()
#
# vendor_tables_chembl <- vendor_tables[
#   zinc_mapping, on = "zinc_id", nomatch = NULL
# ]
#
# fwrite(
#   vendor_tables_chembl,
#   file.path(dir_vendor, "zinc_commercial_compounds.csv.gz"),
#   na = "NA"
# )
#
# vendor_tables_chembl <- read_csv()

# Get compound names from vendors ----------------------------------------------
###############################################################################T

cmpd_name_funcs <- list(
  "targetmol" = function(url) {
    Sys.sleep(2)
    xml2::read_html(url) %>%
      rvest::html_node(".cpdtitle h1") %>%
      rvest::html_text() %>%
      stringr::str_trim()
  },
  "selleck" = function(url) {
    Sys.sleep(5)
    xml2::read_html(url) %>%
      rvest::html_node("td strong") %>%
      rvest::html_text() %>%
      stringr::str_trim()
  },
  "tocris" = function(url) {
    Sys.sleep(10)
    xml2::read_html(url) %>%
      rvest::html_node("div#content_column h1") %>%
      rvest::html_text() %>%
      stringr::str_trim()
  },
  "mce" = function(url) {
    Sys.sleep(2)
    xml2::read_html(url) %>%
      rvest::html_node("th.s_pro_list_name a strong") %>%
      rvest::html_text() %>%
      stringr::str_trim()
  }
)

try_again <- function(fun, n = 3, otherwise = NA_character_, ...) {
  i = 1
  while (i <= n) {
    res <- tryCatch(
      fun(...),
      error = function(e) {
        warning("Error: ", e$message, " Args: ", deparse(list(...)))
        otherwise
      }
    )
    if (!identical(res, otherwise))
      break
    i <- i + 1
  }
  res
}

vendor_url_templates <- tribble(
  ~vendor, ~url_template,
  "targetmol", "https://www.targetmol.com/search?keyword=",
  "mce", "https://www.medchemexpress.com/search.html?q=",
  "selleck", "https://www.selleckchem.com/search.html?searchDTO.searchParam=",
  "tocris", "https://www.tocris.com/products/",
  "enamine", "https://www.enaminestore.com/catalog/",
  "enaminebb", "https://www.enaminestore.com/catalog/"
)

vendor_info_unique <- vendor_tables_chembl %>%
  inner_join(vendor_url_templates, by = "vendor") %>%
  filter(vendor %in% names(cmpd_name_funcs), !is.na(chembl_id)) %>%
  mutate(vendor_url = paste0(url_template, vendor_id)) %>%
  distinct(vendor, vendor_id, vendor_url)

max_tries <- 3

# Runs really slowly
fetch_vendor_names <- function(df) {
  plan(multisession(workers = 3))
  if (!"vendor_name" %in% colnames(df)) {
    df[["vendor_name"]] <- NA_character_
  }
  vendors <- unique(df[["vendor"]])
  df <- df %>%
    mutate(
      done = !is.na(df[["vendor_name"]]),
      fu = pmap(
        list(vendor, vendor_url, seq_len(nrow(.))),
        function(x, y, i) {
          fun <- cmpd_name_funcs[[x]]
          fu <- future(
            {try_again(fun, url = y, n = n)},
            globals = list(fun = fun, y = y, try_again = try_again, n = max_tries),
            packages = "magrittr",
            lazy = TRUE
          )
          attr(fu, "idx") <- i
          fu
        }
      )
    )
  futures <- df %>%
    group_by(vendor) %>%
    summarize(futures = list(fu)) %>%
    {set_names(.[["futures"]], .[["vendor"]])}
  walk(futures, ~run(.x[[1]]))
  idx <- which(!df[["done"]])
  while(length(idx) > 0) {
    for (v in vendors) {
      fus <- futures[[v]]
      if (length(fus) < 1)
        next
      fu <- fus[[1]]
      if(resolved(fu)) {
        fu_idx <- attr(fu, "idx", exact = TRUE)
        df[[fu_idx, "vendor_name"]] <- value(fu)
        df[[fu_idx, "done"]] <- TRUE
        futures[[v]] <- fus[-1]
        if (length(fus) > 1)
          run(fus[[2]])
        message(fu_idx, "done")
      }
    }
    Sys.sleep(0.1)
    idx <- which(!df[["done"]])
    message("remaining", length(idx))
  }
  df
}

vendor_names_result <- fetch_vendor_names(
  vendor_info_unique %>%
    {.[sample(seq_len(length(.)), length(.))]}
)

vendor_names <- vendor_names_result %>%
  select(-done, -fu) %>%
  drop_na(vendor_name)

write_rds(
  vendor_names_result,
  file.path(dir_vendor, "zinc_vendor_names.rds")
)
write_csv(
  vendor_names,
  file.path(dir_vendor, "zinc_vendor_names.csv.gz")
)

vendor_names_chembl <- vendor_tables_chembl %>%
  left_join(
    vendor_names %>%
      select(vendor, vendor_id, vendor_name),
    by = c("vendor", "vendor_id")
  ) %>%
  drop_na(chembl_id)

write_csv(
  vendor_names_chembl,
  file.path(dir_vendor, "zinc_vendor_names_chembl.csv.gz"),
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

wrangle_activity <- Activity(
  name = "Wrangle commercial availability of compounds from ZINC",
  used = c(
    "syn21572663",
    "syn21994981",
    "http://files.docking.org/catalogs/1/chembl23/chembl23.codemap.txt.gz",
    "http://files.docking.org/catalogs/"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/01_zinc_ids.R"
)

syn_zinc <- Folder("zinc", parent = "syn20830877") %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_vendor, "zinc_commercial_compounds.csv.gz"),
  file.path(dir_vendor, "zinc_vendor_names.csv.gz"),
  file.path(dir_vendor, "zinc_vendor_names_chembl.csv.gz")
) %>%
  synStoreMany(parentId = syn_zinc, activity = wrangle_activity)
