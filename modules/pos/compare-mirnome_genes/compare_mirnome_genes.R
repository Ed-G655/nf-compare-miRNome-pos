
## load libraries
library ("dplyr")
library("ggplot2")
library("eulerr")
library("ggvenn")
library("stringr")
library("cowplot")
library("vroom")
library("tidyr")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only

#args[1] <-"21.ref.tsv"

#args[2] <- "21.alt.tsv"

#args[3] <- "21"

## get the mirmap tsv file from args
mirna_ref_file <- args[1]

## get targetscan tsv file from args
mirna_mut_file <- args[2]

## pass to named objects
chromosome <- args[3]

## Read miRNA targets
mirna_ref.df <- read.table(file= mirna_ref_file, header = T,
                          sep = "\t", stringsAsFactors = FALSE)

mirna_alt.df <-read.table(file= mirna_mut_file, header = T,
                           sep = "\t", stringsAsFactors = FALSE)

## Select mirnas targets predicted by both tools
mirna_ref_intersect.df <- mirna_ref.df %>% filter(prediction_tool ==  "both") %>% 
  select(a_Gene_ID, miRNA_ID) %>%  unite(col = "gene_mirna", sep = ";")

mirna_alt_intersect.df <- mirna_alt.df %>% filter(prediction_tool ==  "both") %>% 
  select(a_Gene_ID, miRNA_ID) %>%  unite(col = "gene_mirna", sep = ";")


mirna_ref.v <- mirna_ref_intersect.df %>% pull(gene_mirna) %>% unique()

mirna_alt.v <- mirna_alt_intersect.df %>% pull(gene_mirna) %>% unique()


#Save venn data on dataframe
Venn_data <- data.frame(  chromosome = c(str_interp("${chromosome}")),
                         remain_mirna_genes = c(mirna_alt.v %>%  intersect(mirna_ref.v) %>% length()),
                         lost_mirna_genes = c(mirna_ref.v %>%  setdiff(mirna_alt.v) %>%  length()),
                         gain_mirna_genes = c(mirna_alt.v %>%  setdiff(mirna_ref.v) %>%  length()))

write.table(Venn_data, 
            file = str_interp("${chromosome}_venndata_mirna_genes.tsv"), 
            sep = "\t", row.names = F, 
            col.names = T)

## Sort the ids list within a list for ggvenn
Venn_list <- list(
  A = mirna_ref.v,
  B = mirna_alt.v)

## Name the source of the ids
names(Venn_list) <- c("Lost miRNA-gene pairs","Gain miRNA-gene pairs")

## Ṕlot a Venn diagram
miRNAs_Venn.p <- ggvenn(Venn_list, fill_color = c("#FF595E", "#007F5F"),
                        stroke_size = 0.5, set_name_size = 4 , text_size = 4)

## Save plot
ggsave( filename = str_interp("${chromosome}_mirna_genes.png"),
        plot = miRNAs_Venn.p,
        device = "png",
        height = 7, width = 14,
        units = "in")

## Make eulerr plot
microRNAs_euler <- euler(Venn_list)

microRNAs_euler.p <- plot( x = microRNAs_euler,
                           quantities = TRUE,               
                           main = "miRNA-gene pairs",
                           fill = c("#FF595E", "#007F5F") )                 

# save plot
ggsave( filename = str_interp("${chromosome}_mirna_genes2.png"),        
        plot = microRNAs_euler.p,                
        device = "png",                 
        height = 7,                     
        width = 14,
        units = "in",
        dpi = 300 )                

############## Compare Genes ##############

genes_ref.v <- mirna_ref.df %>% filter(prediction_tool ==  "both") %>% pull(a_Gene_ID) %>% unique()

genes_alt.v <- mirna_alt.df %>% filter(prediction_tool ==  "both") %>% pull(a_Gene_ID) %>% unique()


## Sort the ids list within a list for ggvenn
Venn_list <- list(
  A = genes_ref.v,
  B = genes_alt.v)

## Name the source of the ids
names(Venn_list) <- c("Lost target genes","Gain target genes")

## Ṕlot a Venn diagram
Genes_Venn.p <- ggvenn(Venn_list, fill_color = c("#FF595E", "#007F5F"),
                        stroke_size = 0.5, set_name_size = 4 , text_size = 4)

## Save plot
ggsave( filename = str_interp("${chromosome}_changes_genes.png"),
        plot = Genes_Venn.p,
        device = "png",
        height = 7, width = 14,
        units = "in")

## Make eulerr plot
Genes_euler <- euler(Venn_list)

Genes_euler.p <- plot( x = Genes_euler,
                           quantities = TRUE,               
                           main = "microRNA target genes",
                           fill = c("#FF595E", "#007F5F") )                 

# save plot
ggsave( filename = str_interp("${chromosome}_changes2_genes.png"),        
        plot = Genes_euler.p,                
        device = "png",                 
        height = 7,                     
        width = 14,
        units = "in",
        dpi = 300 )  



#### Save dataframe
Venn_data_genes <- data.frame(  chromosome = c(str_interp("${chromosome}")),
                          remain_genes = c(genes_alt.v %>%  intersect(genes_ref.v) %>% length()),
                          lost_genes = c(genes_ref.v %>%  setdiff(genes_alt.v) %>%  length()),
                          gain_genes = c(genes_alt.v %>%  setdiff(genes_ref.v) %>%  length()))

write.table(Venn_data_genes, 
            file = str_interp("${chromosome}_venndata_genes.tsv"), 
            sep = "\t", row.names = F, 
            col.names = T)

