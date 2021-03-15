library(tidyverse)
library(data.table)
library(qs)
library(openxlsx)
library(synapser)
library(fst)


load_input_data <- function(inputs, syn) {
  # Some columns are always incorrectly read as character by fread
  # https://github.com/Rdatatable/data.table/issues/4802
  column_types <- list(
    lspci_id = as.integer,
    ontarget_IC50_Q1 = as.numeric,
    offtarget_IC50_Q1 = as.numeric,
    Q1 = as.numeric
  )
  inputs %>%
    map(syn) %>%
    map(
      function(x)
        list(
          `.csv` = fread,
          `.tsv` = fread,
          `.rds` = read_rds,
          `.qs` = qread,
          `.xlsx` = read.xlsx,
          `.fst` = read_fst
        ) %>%
        magrittr::extract2(which(str_detect(x, fixed(names(.))))) %>% {
          df <- .(x)
          cols <- intersect(colnames(df), names(column_types)) %>%
            set_names() %>%
            map(function(col) column_types[[col]](df[[col]]))
          set(df, j = names(cols), value = unname(cols))
        }
    )
}

pluck_inputs <- function(inputs, syn_parent) {
  map_chr(
    inputs,
    ~exec(synPluck, !!!c(syn_parent, .x))
  )
}
