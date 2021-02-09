library(tidyverse)
library(data.table)
library(qs)
library(synapser)

load_input_data <- function(inputs, syn) {
  inputs %>%
    map(syn) %>%
    map(
      function(x)
        list(
          `.csv` = partial(
            fread,
            colClasses = c(
              lspci_id = "integer"
            )
          ),
          `.tsv` = fread,
          `.rds` = read_rds,
          `.qs` = qread
        ) %>%
        magrittr::extract2(which(str_detect(x, fixed(names(.))))) %>%
        {.(x)}
    )
}

pluck_inputs <- function(inputs, syn_parent) {
  map(
    inputs,
    ~exec(synPluck, !!!c(syn_parent, .x))
  )
}
