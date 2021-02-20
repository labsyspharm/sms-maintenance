## this script calculates selectivity scores from biochemical data obtained with scripts '02_collecting_data_chembl.R'

library(tidyverse)
library(data.table)
library(here)
library(furrr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

`%nin%` <- Negate(`%in%`)

source(here("utils", "load_save.R"))

# set directories, import files ------------------------------------------------
###############################################################################T

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

inputs <- list(
  dose_response_measurements = c("aggregate_data", "dose_response_measurements.csv.gz"),
  target_dictionary = c("id_mapping", "target_dictionary_wide.csv.gz")
) %>%
  pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  load_input_data(syn = syn)

# set toolscore function -------------------------------------------------------
###############################################################################T

calc_toolscore<-function(data, target_lspci_target_id)
{
  example_subset<-data

  ontarget<-example_subset[lspci_target_id == target_lspci_target_id]
  offtarget<-example_subset[lspci_target_id !=  target_lspci_target_id]

  chembl_active_data<-as.numeric(as.character(ontarget$value))
  chembl_active_low_nM<-chembl_active_data[chembl_active_data<=100]
  if (length(chembl_active_low_nM) > 1) {
    chembl_strength<-7
  }else if (length(chembl_active_data[chembl_active_data<=1000]) > 4) {
    chembl_strength<-4
  }else {
    chembl_strength<-1
  }

  strength<-chembl_strength + 1 # bonus for multiple sources

  ontarget_IC50<-as.numeric(as.character(ontarget[!is.na(ontarget$value),]$value))
  ontarget_IC50_N<-length(ontarget_IC50)
  ontarget_IC50_Q1<-quantile(ontarget_IC50, probs=0.25, na.rm=TRUE)
  offtarget_IC50<-as.numeric(as.character(offtarget[!is.na(offtarget$value),]$value))
  offtarget_IC50_N<-length(offtarget_IC50)
  offtarget_IC50_Q1<-quantile(offtarget_IC50, probs=0.25, na.rm=TRUE)
  Q1_IC50_diff<-log10(offtarget_IC50_Q1)-log10(ontarget_IC50_Q1);
  if (length(ontarget_IC50)==0 || length(offtarget_IC50)==0) {
    wilcox_pval<-10;
  } else {
    w<-wilcox.test(ontarget_IC50, offtarget_IC50, alternative='less');
    if (w$p.value>=1E-15){
      wilcox_pval<-w$p.value;
    } else {
      wilcox_pval<-1E-15;
    }
  }

  IC50<-c(ontarget_IC50, offtarget_IC50)
  IC50_Q1<-quantile(IC50, probs=0.25, na.rm=TRUE)
  investigation_bias<-length(ontarget_IC50)/length(IC50)

  selectivity<-(Q1_IC50_diff/3 + (1-investigation_bias) - log10(wilcox_pval)/15)/3
  tool_score<- strength * selectivity
  #names(tool_score)<-'tool score'
  return (list(tool_score[[1]],strength,selectivity[[1]],investigation_bias,wilcox_pval,
               Q1_IC50_diff[[1]],offtarget_IC50_Q1[[1]],ontarget_IC50_Q1[[1]],ontarget_IC50_N,offtarget_IC50_N))
}

iterate_targets <- function(c.data, ...) {
  toolscore<-list()
  toolscore_index<-0
  # c.data<-as.data.frame(lspci_id_list[[index_cmpd]])
  c.lspci_id<-unique(c.data$lspci_id)
  c.targets<-unique(c.data$lspci_target_id)
  # pb <- txtProgressBar(max = length(c.targets), style = 3)
  # on.exit(close(pb))
  for(index_target in 1:length(c.targets)){
    # setTxtProgressBar(pb, index_target)
    c.df<-list()
    c.lspci_target_id<-c.targets[index_target]
    c.df$lspci_id<-c.lspci_id
    c.df$lspci_target_id<-c.lspci_target_id
    c.return<-calc_toolscore(c.data, c.lspci_target_id)
    if(!is.na(c.return[[3]])){
      toolscore_index<-toolscore_index+1
      c.df$tool_score<-c.return[[1]]
      c.df$strength<-c.return[[2]]
      c.df$selectivity<-c.return[[3]]
      c.df$investigation_bias<-c.return[[4]]
      c.df$wilcox_pval<-c.return[[5]]
      c.df$IC50_diff<-c.return[[6]]
      c.df$ontarget_IC50_Q1<-c.return[[8]]
      c.df$offtarget_IC50_Q1<-c.return[[7]]
      c.df$ontarget_IC50_N<-c.return[[9]]
      c.df$offtarget_IC50_N<-c.return[[10]]
      c.df$N_total<-c.return[[9]]+c.return[[10]]
      toolscore[[toolscore_index]]<-as.data.frame(c.df)
    }
  }
  # print(paste0(index_cmpd,"-",length(lspci_id_list)))
  # print(paste("Done", c.lspci_id))
  rbindlist(toolscore)
}


# calculate toolscores ---------------------------------------------------------
###############################################################################T

# lspci_id_list<-dlply(activities_lspci_id_geneid,.(lspci_id),c)

plan(multicore(workers = 10))
toolscore.b <- input_data[["dose_response_measurements"]][
  ,
  .(
    data = list(.SD)
  ),
  keyby = .(lspci_id)
][
  ,
  data := future_map(
    data,
    iterate_targets,
    .options = furrr_options(
      scheduling = 10L
    )
  )
]

toolscore_all <- toolscore.b %>%
  unnest(data)

fwrite(
  toolscore_all,
  file.path(dir_release, "toolscores.csv.gz")
)

# Calculating selectivity classes ----------------------------------------------
###############################################################################T

is_most_selective <- function(df) {
  with(
    df,
    tool_score >= 5 &
    IC50_diff >= 2 &
    wilcox_pval <= 0.1 &
    strength == 8 &
    investigation_bias <= 0.2
  )
}

is_semi_selective <- function(df) {
  with(
    df,
    IC50_diff >= 1 &
    wilcox_pval <= 0.1 &
    strength >= 5 &
    investigation_bias <= 0.2
  )
}

is_poly_selective <- function(df) {
  with(
    df,
    ontarget_IC50_N > 1 &
    investigation_bias <= 0.2 &
    IC50_diff >= 0 &
    ontarget_IC50_Q1 < 9000
  )
}

is_unknown_selective <- function(df) {
  with(
    df,
    IC50_diff >= 0 &
    ontarget_IC50_Q1 < 9000
  )
}

calculate_selectivity_class <- function(df) {
   case_when(
     is_most_selective(df) ~ "most_selective",
     is_semi_selective(df) ~ "semi_selective",
     is_poly_selective(df) ~ "poly_selective",
     is_unknown_selective(df) ~ "unknown_selective",
     TRUE ~ "other_selective"
   )
}

selectivity_class_order <- c(
  "most_selective",
  "semi_selective",
  "poly_selective",
  "unknown_selective",
  "other_selective"
)

selectivity_classes <- toolscore_all %>%
  mutate(
    selectivity_class = factor(
      calculate_selectivity_class(.),
      levels = selectivity_class_order
    )
  ) %>%
  left_join(
    input_data[["target_dictionary"]][
      , .(lspci_target_id, entrez_gene_id, symbol)
    ],
    by = "lspci_target_id"
  )

fwrite(
  selectivity_classes,
  file.path(dir_release, "selectivity.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

selectivity_calc_activity <- Activity(
  name = "Calculate compound-target selectivity and selectivity classes",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/04_calculating_selectivity.R"
)

syn_selectivity_folder <- synMkdir(syn_release, "selectivity")

c(
  file.path(dir_release, "toolscores.csv.gz"),
  file.path(dir_release, "selectivity.csv.gz")
) %>%
  synStoreMany(parentId = syn_selectivity_folder, activity = selectivity_calc_activity)

