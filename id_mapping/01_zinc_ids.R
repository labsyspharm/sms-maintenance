library(tidyverse)
library(synapser)
library(synExtra)
library(furrr)
library(httr)
library(future)
library(rvest)
library(polite)
library(here)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

chembl_zinc_mapping <- syn("syn21572663") %>%
  read_csv()

dir_vendor <- file.path(dir_release, "vendor")
dir.create(dir_vendor, showWarnings = FALSE)

# Fetch vendor info from ZINC --------------------------------------------------
###############################################################################T

vendor_libraries <- tribble(
  ~id, ~url, ~vendor_url,
  "targetmol", "http://files.docking.org/catalogs/40/", "https://www.targetmol.com/search?keyword=",
  "mce", "http://files.docking.org/catalogs/50/", "https://www.medchemexpress.com/search.html?q=",
  "selleck", "http://files.docking.org/catalogs/50/", "https://www.selleckchem.com/search.html?searchDTO.searchParam=",
  "tocris", "http://files.docking.org/catalogs/50/", "https://www.tocris.com/products/",
  "enamine", "http://files.docking.org/catalogs/50/", "https://www.enaminestore.com/catalog/",
  "enaminebb", "http://files.docking.org/catalogs/50/", "https://www.enaminestore.com/catalog/"
) %>%
  mutate(
    codemap_url = paste0(url, id, "/", id, ".codemap.txt.gz"),
    info_url = paste0(url, id, "/", id, ".info.txt.gz")
  )

pwalk(
  vendor_libraries,
  function(codemap_url, info_url, id, ...) {
    download.file(codemap_url, file.path(dir_vendor, paste0("codemap_", id, ".txt.gz")))
    download.file(info_url, file.path(dir_vendor, paste0("info_", id, ".txt.gz")))
  }
)

vendor_tables <- vendor_libraries %>%
  transmute(
    id,
    codemap = map(id, ~read_delim(file.path(dir_vendor, paste0("codemap_", .x, ".txt.gz")), delim = " ", col_names = c("zinc_id", "vendor_id"), col_types = "cc")),
    info = map2(
      id, vendor_url,
      ~read_tsv(
        file.path(dir_vendor, paste0("info_", .x, ".txt.gz")),
        col_names = c("vendor_id", "zinc_id", "inchi_key", "tranche", "notes"),
        col_types = "ccccc"
      ) %>%
        mutate(
          vendor_url = paste0(.y, vendor_id)
        )
    )
  )

# Mapping from Zinc IDs to Chembl provided by Zinc directly
download.file("http://files.docking.org/catalogs/1/chembl23/chembl23.codemap.txt.gz", file.path(dir_vendor, "codemap_chembl23.txt.gz"))
chembl_codemap <- read_delim(file.path(dir_vendor, "codemap_chembl23.txt.gz"), delim = " ", col_names = c("zinc_id", "chembl_id"), col_types = "cc")

zinc_mapping <- chembl_zinc_mapping %>%
  select(-uci) %>%
  bind_rows(chembl_codemap) %>%
  distinct()

vendor_tables_chembl <- vendor_tables %>%
  dplyr::select(vendor = id, info) %>%
  unnest(info) %>%
  dplyr::select(-tranche) %>%
  left_join(
    zinc_mapping, by = "zinc_id"
  )

write_csv(
  vendor_tables_chembl,
  file.path(dir_vendor, "zinc_commercial_compounds.csv.gz")
)

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

vendor_info_unique <- vendor_tables_chembl %>%
  filter(vendor %in% names(cmpd_name_funcs), !is.na(chembl_id)) %>%
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
    "http://files.docking.org/catalogs/1/chembl23/chembl23.codemap.txt.gz",
    vendor_libraries$info_url
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
