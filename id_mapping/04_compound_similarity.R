library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)
library(data.table)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Loading files ----------------------------------------------------------------
###############################################################################T

inputs <- c(
  fingerprints = synPluck(syn_release, "fingerprints", "all_compounds_fingerprints.rds")
)

all_compounds_fingerprints <- inputs[["fingerprints"]] %>%
  syn() %>%
  read_rds()

# Calculate similarity between all compounds -----------------------------------
###############################################################################T

wd <- file.path("/n", "scratch3", "users", "c", "ch305", "simsearch")
dir.create(wd, showWarnings = FALSE)

# Write FPS files
pwalk(
  all_compounds_fingerprints,
  function(fp_name, data, ...) {
    fps_path <- file.path(wd, paste0(fp_name, ".fps"))
    write_lines(
      c(
        "#FPS1",
        "#num_bits=2048"
      ),
      fps_path
    )
    fwrite(
      data[, .(fingerprints, names)],
      file = fps_path,
      append = TRUE,
      sep = "\t",
      col.names = FALSE
    )
  }
)

similarity_fun <- function(input_file, output_file) {
  library(processx)
  library(tidyverse)
  library(lspcheminf)

  lspcheminf_script <- paste0(
    # Tail cut off three header lines
    "tail -n +3 ", input_file, " | sort -k1,1 --parallel=8 -S15G > ", output_file
  )

  message("script run ", lspcheminf_script)

  p <- process$new(
    "bash", c("-c", lspcheminf_script), stderr = "|", stdout = "|"
  )

  on.exit({
    message("sort output ", p$read_error_lines(), p$read_output_lines())
    p$kill()
  })

  p$wait()

  if (file.size(output_file) < 1000)
    stop("Unsuccesful")
}

# similarity_search_input <- all_compounds_fingerprints %>%
  # transmute(
  #   fp_name, fp_type, fp_args,
  #   input_file = file.path(wd, paste0(fp_name, ".fps")),
  #   output_file = file.path(wd, paste0(fp_name, ".txt"))
  # )
similarity_search_input <- tibble(fp_name = c("morgan_normal", "morgan_chiral", "topological_normal")) %>%
  mutate(
    input_file = file.path(wd, paste0(fp_name, ".fps")),
    output_file = file.path(wd, paste0(fp_name, "_sorted.fps"))
  )

pwalk(
  similarity_search_input,
  function(input_file, output_file, ...) {
    similarity_fun(input_file, output_file)
  }
)

# Load similarities ------------------------------------------------------------
###############################################################################T

load_similarities <- function(input_file) {
  # n_compounds <- 20792499L
  # n_compounds_p <- run(
  #   "wc", c("-l", input_file)
  # )
  # n_componds <- str_split_fixed(n_compounds_p[["stdout"]], fixed(" "), 2)[[1]] %>%
  #   as.integer()
  # name_vec <- character(n_compounds)
  # gruop_vec <-
  # con = file(input_file, "r")
  # while (TRUE) {
  #   line = readLines(con, n = 1)
  #   if ( length(line) == 0 ) {
  #     break
  #   }
  #   print(line)
  # }
  # close(con)
  sims <- fread(
    input_file,
    sep = "\t",
    header = FALSE,
    col.names = c("fingerprint", "id")
  )
  sims[
    ,
    .(
      identity_group = (fingerprint[-1L] != fingerprint[-length(fingerprint)]) %>%
        {c(0L, cumsum(.))},
      id
    )
  ]
}

similarity_search_result <- similarity_search_input %>%
  mutate(
    data = map(
      output_file,
      load_similarities
    )
  )

write_rds(
  similarity_search_result %>%
    select(fp_name, data),
  file.path(dir_release, "all_compounds_similarity.rds"),
  compress = "gz"
)
# similarity_search_result <- read_rds(file.path(dir_release, "all_compounds_similarity.rds"))


# Calculate compound mass -----------------------------------
###############################################################################T

# To establish identity between compounds, require that the topological FP
# matches and one of the Morgan FPs and that their molecular mass is identical

library(furrr)
# library(future.apply)

plan(multicore(workers = 8))

cmpd_mass_input <- all_compounds_fingerprints$data[[1]] %>%
  select(names, compound) %>%
  chunk_df(1000, seed = 1) %>%
  enframe("chunk", "data") %>%
  chunk_df(50, seed = 1)

cmpd_mass_raw <- cmpd_mass_input %>%
  # future_map(
  #   ~molecular_mass(
  #     set_names(.x[["compound"]], .x[["names"]])
  #   ),
  #   .progress = TRUE,
  #   .options = furrr_options(globals = FALSE, chunk_size = 1L)
  # ) %>%
  # map(
  #   function(df) {
  #     df %>%
  #       pull(data) %>%
  #       future_lapply(
  #         function(x)
  #           safely(molecular_mass)(
  #             set_names(x[["compound"]], x[["names"]])
  #           ),
  #         future.globals = FALSE,
  #         future.packages = c("lspcheminf", "rlang"),
  #         future.scheduling = TRUE
  #       )
  #   }
  # )
  map(
    function(df) {
      df %>%
        pull(data) %>%
        future_map(
          function(x)
            safely(molecular_mass)(
              set_names(x[["compound"]], x[["names"]])
            ),
          .progress = TRUE
        )
    }
  )
  # map(
  #   ~molecular_mass(
  #     set_names(.x[["compound"]], .x[["names"]])
  #   )
  # ) %>%

cmpd_mass_all <- cmpd_mass_raw %>%
  unlist(recursive = FALSE)

cmpd_mass_errors <- cmpd_mass_all %>%
  map("error")

cmpd_mass_errors %>% map_lgl(is.null) %>% all()


cmpd_mass <- cmpd_mass_all %>%
  map("result") %>%
  bind_rows()

write_csv(
  cmpd_mass,
  file.path(dir_release, "compound_masses.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Calculate chemical similarity between all compounds",
  used = c(
    "syn20692501",
    "syn22080194"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_similarity.R"
)

syn_id_mapping <- Folder("id_mapping", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_similarity.rds"),
  file.path(dir_release, "compound_masses.csv.gz")
) %>%
  synStoreMany(parentId = syn_id_mapping, activity = activity)

