
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

#args[1] <-"targets.ref.tsv"

#args[2] <- "targets.alt.tsv"

#args[3] <- "targets.changes" # output file


## get the mirmap tsv file from args
mirna_ref_file <- args[1]

## get targetscan tsv file from args
mirna_mut_file <- args[2]

## pass to named objects
mirna_changes <- args[3]

## Read miRNA targets
mirna_ref.df <- read.table(file= mirna_ref_file, header = T,
                           sep = "\t", stringsAsFactors = FALSE)

mirna_alt.df <- read.table(file= mirna_mut_file, header = T,
                           sep = "\t", stringsAsFactors = FALSE)

## Select mirnas targets predicted by both tools
mirna_ref_intersect.df <- mirna_ref.df %>% filter(prediction_tool ==  "both")

mirna_alt_intersect.df <- mirna_alt.df %>% filter(prediction_tool ==  "both")

## Select mirnas targets predicted by any tool
mirna_ref_all.df <- mirna_ref.df 

mirna_alt_all.df <- mirna_alt.df 

## Get Lost target mirna pairs
lost_targets <- mirna_ref_intersect.df %>% setdiff(mirna_alt_all.df)

## Get Gain target mirna pairs
gain_targets <- mirna_alt_intersect.df %>% setdiff(mirna_ref_all.df)

## Get remained targets r
remained_targets.df <-mirna_ref_intersect.df %>%  intersect(mirna_alt_intersect.df)


## Define if one target is lost, gained o remained
lost_targets <- lost_targets %>%  mutate(target = "lost") %>% 
  select(-prediction_tool)
gain_targets <- gain_targets  %>%  mutate(target = "gained") %>% 
  select(-prediction_tool)
remained_targets.df <- remained_targets.df  %>% 
  mutate(target = "remained") %>% 
  select(-prediction_tool)


## Merge the miRNA targets gained and lost into a single dataframe
target_changes.df <- full_join(x = lost_targets, y = gain_targets,
                               by = c("GeneID","miRNA_ID", "UTR_start",
                                      "UTR_end", "Site_type", "target", "target_ID",
                                      "chrom") )

## Merge all miRNA targets ids into a single dataframe
All_targets.df <- full_join(x = target_changes.df, y = remained_targets.df,
                            by = c("GeneID","miRNA_ID", "UTR_start",
                                   "UTR_end", "Site_type", "target", "target_ID", 
                                   "chrom") )

## Save dataframe
write.table(All_targets.df, file = str_interp("${mirna_changes}.tsv"), sep = "\t", na = "NA", quote = F, row.names = F)






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
names(Venn_list) <- c("miRNAs REF","miRNAs ALT")

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
ggsave( filename = str_interp("${mirna_changes}_3.png"),
        plot = miRNAs_Venn.p,
        device = "png",
        height = 7, width = 15,
        units = "in")

count_changes.df <- All_targets.df %>% group_by(miRNA_ID, chrom, target) %>%
  summarise(Number_of_Targets = n())

count_lost_gains.df <- subset(count_changes.df, target == "lost") %>%  mutate(Number_of_Targets = Number_of_Targets*-1) %>%
  rbind(subset(count_changes.df, target == "gained")) %>%  arrange(chrom)


count_lost_gains.df <- count_lost_gains.df %>% arrange(-Number_of_Targets)

paleta <- c("lost" =  "#F94144",
            "gained" = "springgreen3") 

min <- count_lost_gains.df$Number_of_Targets %>% min()
max <- count_lost_gains.df$Number_of_Targets %>% max()
quartile <- abs(count_lost_gains.df$Number_of_Targets) %>%  max()/4
quartile <- ceiling(quartile)

chromosomes.v <- count_lost_gains.df %>% ungroup() %>%  pull(chrom) %>% unique()

lapply(chromosomes.v, function(chroms) {
  
  piramide.p <- ggplot(count_lost_gains.df, aes(x = miRNA_ID, y = Number_of_Targets, fill = target )) + 
    geom_col(data = subset(count_lost_gains.df, target == "lost" & chrom == chroms), 
             width = 0.5, fill = "#F94144") + 
    geom_col(data = subset(count_lost_gains.df, target ==  "gained" & chrom == chroms), 
             width = 0.5, fill = "springgreen3") +
    coord_flip() + scale_y_continuous(
      breaks = c(seq(min, 0, by = quartile), 
                 seq(0, max, by = quartile)),
      labels = c(seq(min, 0, by = quartile)*-1, 
                 seq(0, max, by = quartile))
    ) + labs(y= "Numero de pares miRNA/blanco", x = "miRNA", color = "Legend") +
    scale_color_manual(values = paleta) +
    labs(title = "Sitos blanco por miRNA y sus cambios debido a mutaciones en el miRNA") +
    theme_minimal() 
  
  ggsave( filename = str_interp("${chroms}_changes.png"), 
          plot = piramide.p,
          device = "png",
          height = 7, width = 14,
          units = "in")
})


# plot gain, lost and remain targets  
lapply(chromosomes.v, function(chroms) {
gain_and_lost.p <- ggplot(subset(count_changes.df,  chrom == chroms), aes(x = miRNA_ID, 
                                                y = Number_of_Targets, 
                                                     fill = target)) + 
  geom_bar(position = "stack", stat = "identity") + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 0, vjust = 0 )) + 
  scale_y_continuous(expand = c(0,0)) +
  ylab("Numero de pares miRNA/blanco") +
  xlab("miRNA") +
  labs(title = "Sitos blancos por miRNA y sus cambios debido a mutaciones en el miRNA") +
  labs(fill="Target change") +
  scale_fill_manual(values = c("lost" = "#c87570",
                               "gained" = "#70c875",
                               "remained" = "#7570c8"))

ggsave( filename =str_interp("${chroms}_barplot.png"),
        plot = gain_and_lost.p,
        device = "png",
        height = 14, width = 28,
        units = "in")
})

