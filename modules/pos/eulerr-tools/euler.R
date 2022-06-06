
## load libraries
library ("dplyr")
library("ggplot2")
library("eulerr")
library("ggvenn")
library("stringr")
library("cowplot")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only

#args[1] <-"22.alt.tsv"

#args[2] <- "22.ref.tsv"

#args[3] <- "targets.changes" # output file


## get the mirmap tsv file from args
mirna_ref_file <- args[1]

## get targetscan tsv file from args
mirna_mut_file <- args[2]

## pass to named objects
chromosome <- args[3]

## Read miRNA targets
mirna_ref.df <- read.table(file= mirna_ref_file, header = T,
                           sep = "\t", stringsAsFactors = FALSE) %>% select(target_ID, prediction_tool)

mirna_alt.df <- read.table(file= mirna_mut_file, header = T,
                           sep = "\t", stringsAsFactors = FALSE) %>% select(target_ID, prediction_tool)

## Select mirnas targets predicted by TargetScan
mirna_ref_targetscan.df <- mirna_ref.df %>% filter(prediction_tool ==  "targetscan" |
                                                   prediction_tool == "both") %>% select(target_ID)

mirna_mut_targetscan.df <- mirna_alt.df %>% filter(prediction_tool ==  "targetscan" |
                                                   prediction_tool == "both") %>% select(target_ID)


## Select mirnas targets predicted by mirmap
mirna_ref_mirmap.df <- mirna_ref.df %>% filter(prediction_tool ==  "mirmap" |
                                                 prediction_tool == "both") %>% select(target_ID)

mirna_mut_mirmap.df <- mirna_alt.df %>% filter(prediction_tool ==  "mirmap" |
                                               prediction_tool == "both") %>% select(target_ID)



## Make a vector with the mirnas targets predicted by both tools
mirna_ref_targetscan.v <- mirna_ref_targetscan.df %>%  pull(target_ID) %>%  unique()

mirna_ref_mirmap.v <- mirna_ref_mirmap.df %>%  pull(target_ID) %>%  unique()

mirna_mut_targetscan.v <- mirna_mut_targetscan.df %>%  pull(target_ID) %>%  unique()

mirna_mut_mirmap.v <- mirna_mut_mirmap.df %>%  pull(target_ID) %>%  unique()

## Sort the ids list within a list for ggvenn
Venn_list <- list(
  A = mirna_ref_targetscan.v,
  B = mirna_ref_mirmap.v,
  C = mirna_mut_targetscan.v,
  D = mirna_mut_mirmap.v)

## Name the source of the ids
names(Venn_list) <- c("REF_TargetScan","REF_miRmap","MUT_TargetScan","MUT_miRmap")

## Ṕlot a Venn diagram
miRNAs_Venn.p <- ggvenn(Venn_list, fill_color = c("#D9ED92", "#99D98C", "#168AAD", "#1E6091"),
                        stroke_size = 0.5, set_name_size = 4 , text_size = 4)


## Save plot
ggsave( filename = str_interp("${chromosome}_changes3.png"),
        plot = miRNAs_Venn.p,
        device = "png",
        height = 7, width = 15,
        units = "in")
