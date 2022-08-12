
## load libraries
library ("dplyr")
library("ggplot2")
library("eulerr")
library("ggvenn")
library("stringr")
library("cowplot")
library("tidyr")
library("scales")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only

#args[1] <-"21.ref.tsv"

#args[2] <- "21.alt.tsv"

#args[3] <- "21" # output file


## get the mirmap tsv file from args
mirna_ref_file <- args[1]

## get targetscan tsv file from args
mirna_mut_file <- args[2]

## pass to named objects
chromosome <- args[3]

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
                               by = c("a_Gene_ID","miRNA_ID", "UTR_start",
                                      "UTR_end", "Site_type", "target", "target_ID",
                                      "chrom") )

## Merge all miRNA targets ids into a single dataframe
All_targets.df <- full_join(x = target_changes.df, y = remained_targets.df,
                            by = c("a_Gene_ID","miRNA_ID", "UTR_start",
                                   "UTR_end", "Site_type", "target", "target_ID", 
                                   "chrom") )

count_changes.df <- All_targets.df %>% group_by(miRNA_ID, chrom, target) %>%
  summarise(Number_of_Targets = n())


count_changes_long.df <- count_changes.df %>% spread(key = target, value = Number_of_Targets)

count_changes_long.df <- count_changes_long.df %>%  mutate(across(everything(), .fns = ~replace_na(.,0))) 


count_changes_long.df <- count_changes_long.df %>% mutate(total_ref_targets = lost + remained)

count_changes_long.df <- count_changes_long.df %>%  mutate( percent_lost = (1/total_ref_targets) * -lost) %>% 
  mutate( percent_gain = (1/total_ref_targets)*gained) %>% 
  mutate( percent_remain = (1/total_ref_targets)*remained)


count_changes_wide.df <- count_changes_long.df %>% select(-lost, -remained, -total_ref_targets, -gained) %>% 
  gather(key = "target", value = percent, percent_lost , percent_gain, percent_remain )


paleta <- c("lost" =  "#F94144",
            "gained" = "springgreen3") 

  piramide.p <- ggplot(count_changes_wide.df, aes(x = miRNA_ID, y = percent, fill = target )) + 
    geom_col(data = subset(count_changes_wide.df, target == "percent_lost"), 
             width = 0.5 , fill = "#F94144", color = "black", alpha = 0.5) + 
    geom_col(data = subset(count_changes_wide.df, target ==  "percent_gain"), 
             width = 0.5 , fill = "springgreen3", color = "black", alpha = 0.5) +
    coord_flip() +  scale_y_continuous(labels = label_percent()) +
    labs(y= "Numero de pares miRNA/blanco", x = "miRNA", color = "Legend") +
    scale_color_manual(values = paleta) +
    labs(title = "Sitos blanco por miRNA y sus cambios debido a mutaciones en el miRNA") +
    theme_minimal_grid() + theme(axis.text.y  = element_text(face="bold", size=5, angle= 30))
  
  ggsave( filename = str_interp("${chromosome}_percent.png"), 
          plot = piramide.p,
          device = "png",
          height = 7, width = 14,
          units = "in")
  
  count_changes_long2.df <- count_changes.df %>% spread(key = target, value = Number_of_Targets)
  
  count_changes_long2.df <- count_changes_long.df %>%  mutate(across(everything(), .fns = ~replace_na(.,0))) 
  
  count_changes_long2.df <- count_changes_long.df %>% mutate(total_targets = lost + gained + remained)
  
  count_changes_long2.df <- count_changes_long2.df %>%  mutate( percent_lost = (1 / total_targets)* lost) %>% 
    mutate( percent_gain = (1/total_targets)*gained) %>% 
    mutate( percent_remain = (1/total_targets)*remained)
  
  count_changes_wide2.df <- count_changes_long2.df %>% select(-lost, -remained, -total_ref_targets, -gained) %>% 
    gather(key = "target", value = percent, percent_lost , percent_gain, percent_remain )
  
  write.table(x = count_changes_long2.df, 
              file = str_interp("${chromosome}_percent.tsv"), 
              sep = "\t", 
              row.names =  F, 
              col.names = T)
  
# plot gain, lost and remain targets  
gain_and_lost.p <- ggplot(subset(count_changes_wide2.df), aes(x = miRNA_ID, 
                                                y = percent, 
                                                     fill = target)) + 
  geom_bar(position = "stack", stat = "identity") + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 0, vjust = 0 )) + 
  scale_y_continuous(labels = label_percent(), expand = c(0,0)) +
  ylab("Porcentaje de pares miRNA/blanco") +
  xlab("miRNA") +
  labs(title = "Sitos blancos por miRNA y sus cambios debido a mutaciones en el miRNA") +
  labs(fill="Target change") +
  scale_fill_manual(values = c("percent_lost" = "#c87570",
                               "percent_gain" = "#70c875",
                               "percent_remain" = "#7570c8"))

ggsave( filename =str_interp("${chromosome}_barplot_percent.png"),
        plot = gain_and_lost.p,
        device = "png",
        height = 14, width = 28,
        units = "in")

