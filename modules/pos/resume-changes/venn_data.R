# Cargamos pacman
library( "pacman" )

# cargamos mas paquetes
p_load( "vroom",
        "dplyr",
        "ggplot2",
        "purrr")

# Cargamos todos los tsv en un solo dataframe
filenames <- list.files( pattern="*.tsv" )

venn_data.df <- map_df( filenames, vroom)
#Sum colums
Total_venn_data.df <- venn_data.df %>%  summarise( total_remain_targets = sum(remain_targets)) %>% 
        bind_cols(venn_data.df %>% summarise( total_lost_targets = sum(lost_targets))) %>% 
        bind_cols(venn_data.df %>%  summarise( total_gain_targets = sum(gain_targets)))

#Write dataframe of data by chromosome
write.table(venn_data.df, 
            file = "changes_by_chr.tsv", 
            sep = "\t", 
            row.names = F, 
            col.names = T)

#Write dataframe of total data
write.table(Total_venn_data.df, 
            file = "total_changes.tsv", 
            sep = "\t", 
            row.names = F, 
            col.names = T)
