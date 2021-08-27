library(tidyverse)
library(data.table)
library(qs)
library(openxlsx)
library(synapser)
library(fst)

.load_functions <- list(
  `.csv` = fread,
  `.tsv` = fread,
  `.rds` = read_rds,
  `.qs` = qread,
  `.xlsx` = read.xlsx,
  `.fst` = partial(read_fst, as.data.table = TRUE)
)

# Some columns are always incorrectly read as character by fread
# https://github.com/Rdatatable/data.table/issues/4802
.column_types <- list(
  lspci_id = as.integer,
  ontarget_IC50_Q1 = as.numeric,
  offtarget_IC50_Q1 = as.numeric,
  Q1 = as.numeric
)

.load_data <- function(x) {
  message("Loading ", x)
  read_fun <- .load_functions[[
    which(str_detect(x, stringr::fixed(names(.load_functions))))
  ]]
  df <- read_fun(x)
  cols <- intersect(colnames(df), names(.column_types)) %>%
    set_names() %>%
    map(function(col) .column_types[[col]](df[[col]]))
  set(df, j = names(cols), value = unname(cols))
  df
}

load_input_data <- function(inputs, syn) {
  inputs %>%
    map(syn) %>%
    map(.load_data)
}

pluck_inputs <- function(inputs, syn_parent) {
  map_chr(
    inputs,
    ~exec(synPluck, !!!c(syn_parent, .x))
  )
}
