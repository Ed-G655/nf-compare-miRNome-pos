
## load libraries
library ("dplyr")
library("ggplot2")
library("eulerr")
library("ggvenn")
library("stringr")
library("cowplot")
library("vroom")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only

#args[1] <-"targets.changes.tsv"

#args[2] <- "changes"


All_targets <- args[1] 

mirna_changes <- args[2]

All_targets.df <- vroom(All_targets)

## Make a vector with the mirnas targets predicted by both tools
mirna_ref_intersect.v <- All_targets.df %>% filter(target ==  "lost" | 
                                                     target ==  "remained" ) %>% pull(target_ID)

mirna_alt_intersect.v <- All_targets.df %>% filter(target ==  "gained" | 
                                                     target ==  "remained" ) %>% pull(target_ID)

## Sort the ids list within a list for ggvenn
Venn_list <- list(
  A = mirna_ref_intersect.v,
  B = mirna_alt_intersect.v)

## Name the source of the ids
names(Venn_list) <- c("Lost targets","Gain targets")

## Ṕlot a Venn diagram
miRNAs_Venn.p <- ggvenn(Venn_list, fill_color = c("#FF595E", "#007F5F"),
                        stroke_size = 0.5, set_name_size = 4 , text_size = 4)

## Save plot
ggsave( filename = str_interp("${mirna_changes}.png"),
        plot = miRNAs_Venn.p,
        device = "png",
        height = 7, width = 14,
        units = "in")

## Make eulerr plot
microRNAs_euler <- euler(Venn_list)

microRNAs_euler.p <- plot( x = microRNAs_euler,
                           quantities = TRUE,               
                           main = "microRNAS iDs",
                           fill = c("#FF595E", "#007F5F") )                 

# save plot
ggsave( filename = str_interp("${mirna_changes}_2.png"),        
        plot = microRNAs_euler.p,                
        device = "png",                 
        height = 7,                     
        width = 14,
        units = "in",
        dpi = 300 )                    