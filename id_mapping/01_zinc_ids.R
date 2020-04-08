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
  file.path(dir_release, "zinc_vendors.csv.gz")
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
    Sys.sleep(2)
    xml2::read_html(url) %>%
      rvest::html_node("td strong") %>%
      rvest::html_text() %>%
      stringr::str_trim()
  },
  "tocris" = function(url) {
    Sys.sleep(2)
    xml2::read_html(url) %>%
      rvest::html_node("div#content_column h1") %>%
      rvest::html_text() %>%
      stringr::str_trim()
  }
  # "mce" = function(url, session) {
  #   url_l <- parse_url(url)
  #   browser()
  #   session <- nod(session, url_l$path, verbose = TRUE)
  #   scrape(session, query = url_l$query, verbose = TRUE) %>%
  #     html_node("th.s_pro_list_name a strong") %>%
  #     html_text() %>%
  #     str_trim()
  # }
)

# scrape_sessions <- vendor_libraries %>%
#   {set_names(.$vendor_url, .$id)} %>%
#   map(parse_url) %>%
#   map(~build_url(.x %>% magrittr::inset2("path", NULL) %>% magrittr::inset2("query", NULL))) %>%
#   map(bow)

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
            {fun(y)},
            globals = list(fun = fun, y = y),
            packages = "magrittr",
            lazy = TRUE
          )
          attr(fu, "idx") <- i
          fu
        }
      )
    )
  futures <- list()
  idx <- which(!df[["done"]])
  while(length(idx) > 0) {
    for (i in idx) {
      v <- df[[i, "vendor"]]
      fu <- futures[[v]]
      if (is.null(fu)) {
        futures[[v]] <- df[[i, "fu"]]
        run(futures[[v]])
      }
      else if(resolved(fu)) {
        fu_idx <- attr(fu, "idx", exact = TRUE)
        df[[fu_idx, "vendor_name"]] <- value(fu)
        df[[fu_idx, "done"]] <- TRUE
        futures[[v]] <- NULL
        message(fu_idx, "done")
      }
    }
    Sys.sleep(0.1)
    idx <- which(!df[["done"]])
    message("remaining", length(idx))
  }
  df
}

vendor_names <- fetch_vendor_names(
  vendor_info_unique %>%
    {.[sample(seq_len(length(.)), length(.))]}
)

  mutate(
    vendor_name = future_map2_chr(
      vendor, vendor_url,
      function(x, y) {
        ret <- possibly(cmpd_name_funcs[[x]], NA_character_)(y)
        message(if (is.na(ret)) "fail" else "success")
        Sys.sleep(2)
        ret
      },
      .progress = TRUE
    )
  )



vendor_tables_lscpi <- compound_mapping %>%
  mutate(
    data = map(
      data,
      ~left_join(vendor_tables_chembl, select(.x, id, lspci_id), by = c("chembl_id" = "id")) %>%
        select(-inchi_key, -chembl_id, -zinc_id, -notes) %>%
        left_join(select(vendor_libraries, id, vendor_url), by = c("vendor" = "id")) %>%
        # Some Chembl IDs can't be mapped because they are no longer used
        # in the current version, removing those
        drop_na(lspci_id) %>%
        distinct()
    )
  )

write_rds(
  vendor_tables_lscpi,
  file.path(dir_release, "compound_commercial_info_zinc.rds"),
  compress = "gz"
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
  file.path(dir_release, "zinc_vendors.csv.gz")
) %>%
  synStoreMany(parentId = syn_zinc, activity = wrangle_activity)
