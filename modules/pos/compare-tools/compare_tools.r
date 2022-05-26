
## load libraries
library("dplyr")
library("ggplot2")
library("ggvenn")
library("vroom")
library("stringr")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only

#args[1] <- "All_targets_ref.tsout"

#args[2] <- "All_targets_ref.mirmapout" 

#args[3] <- "sample.bed" 

#args[4] <- "output" # output file base name


## pass to named objects\
targetscan  <- args[1]

     mirmap <- args[2]

miRNA_bed   <- args[3]

output_Name <- args[4]

# import targetscan files 
targetscan.df <- vroom(targetscan)

# import targetscan files 
mirmap.df <- vroom(mirmap)


#Open miRNA BED
miRNA_bed.df <- vroom(file = miRNA_bed, col_names = F)

names(miRNA_bed.df)[1] <- "chrom"
names(miRNA_bed.df)[2] <- "start"
names(miRNA_bed.df)[3] <- "end"
names(miRNA_bed.df)[4] <- "mir"
names(miRNA_bed.df)[5] <- "score"
names(miRNA_bed.df)[6] <- "strand"
names(miRNA_bed.df)[8] <- "type"


miRNA_bed.df <- miRNA_bed.df %>%  mutate(miRNA_ID = str_c(miRNA_bed.df$mir,
                                         "(", 
                                         miRNA_bed.df$strand,
                                         ")") )


miRNA_bed.df <-miRNA_bed.df %>%  select(chrom, miRNA_ID) %>% unique()

targetscan_by_chrom.df <- inner_join(targetscan.df, miRNA_bed.df, by = "miRNA_ID") %>% 
  select(a_Gene_ID, miRNA_ID, UTR_start, UTR_end, Site_type, chrom)

mirmap_by_chrom.df <- inner_join(mirmap.df, miRNA_bed.df, by = "miRNA_ID")

names(mirmap_by_chrom.df)[1] <- "a_Gene_ID"

targetscan_by_chrom.df <- targetscan_by_chrom.df %>%  mutate( target_ID = str_c(targetscan_by_chrom.df$a_Gene_ID,
                                                                                ";",
                                                                                targetscan_by_chrom.df$miRNA_ID,
                                                                                ";",
                                                                                targetscan_by_chrom.df$UTR_start,
                                                                                ";",
                                                                                targetscan_by_chrom.df$UTR_end,
                                                                                ";",
                                                                                targetscan_by_chrom.df$Site_type))


mirmap_by_chrom.df <- mirmap_by_chrom.df %>%  mutate( target_ID = str_c(mirmap_by_chrom.df$a_Gene_ID,
                                                                        ";",
                                                                        mirmap_by_chrom.df$miRNA_ID,
                                                                        ";",
                                                                        mirmap_by_chrom.df$UTR_start,
                                                                        ";",
                                                                        mirmap_by_chrom.df$UTR_end,
                                                                        ";",
                                                                        mirmap_by_chrom.df$Site_type))


## Make a list with the unique IDs
IDs_mirmap.v <- mirmap_by_chrom.df  %>%  pull(target_ID) %>%  unique()

IDs_targetscan.v <- targetscan_by_chrom.df %>%  pull(target_ID) %>%  unique()

## Sort the ids list within a list for ggvenn
Venn_list <- list(
  A = IDs_mirmap.v,
  B = IDs_targetscan.v)

## Name the source of the ids
names(Venn_list) <- c("miRmap","TargetScan")

## Ṕlot a Venn diagram
miRNAs_Venn.p <- ggvenn(Venn_list, fill_color = c("#EE6352", "#59CD90"),
                        stroke_size = 0.5, set_name_size = 4 , text_size = 4)

## Save plot
ggsave( filename =str_interp("${output_Name}.png"),
        plot = miRNAs_Venn.p,
        device = "png",
        height = 7, width = 14,
        units = "in")

## Get mirna targets present in both tools
targets_intersect.df <- intersect(targetscan_by_chrom.df, mirmap_by_chrom.df)
intersect(targetscan_by_chrom.df, mirmap_by_chrom.df) %>% pull() %>%  length()

##Get miRNA targets ids that differ
mirmap_differ.df <- mirmap_by_chrom.df %>% setdiff(targetscan_by_chrom.df)
targetscan_differ.df <- targetscan_by_chrom.df %>% setdiff(mirmap_by_chrom.df)

## Define the source of miRNA ID
targets_intersect.df <-targets_intersect.df %>%  mutate(prediction_tool = "both")

mirmap_differ.df <- mirmap_differ.df %>%  mutate(prediction_tool = "mirmap")

targetscan_differ.df <- targetscan_differ.df %>%  mutate(prediction_tool = "targetscan")

## Merge the miRNA targets ids that differ into a single dataframe
targets_differ.df <- full_join(x = mirmap_differ.df, y = targetscan_differ.df,
                               by = c("a_Gene_ID", "miRNA_ID", "UTR_start",
                                      "UTR_end",  "Site_type", "chrom", 
                                      "target_ID", "prediction_tool") )

## Merge all miRNA targets ids into a single dataframe
All_targets.df <- full_join(x = targets_intersect.df, y = targets_differ.df,
                            by = c("a_Gene_ID", "miRNA_ID", "UTR_start",
                                   "UTR_end",  "Site_type", "chrom", 
                                   "target_ID", "prediction_tool") )

## Save dataframe
write.table(All_targets.df, file = str_interp("${output_Name}.tsv"), sep = "\t", na = "NA", quote = F, row.names = F)

