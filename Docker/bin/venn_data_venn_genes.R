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
Total_venn_data.df <- venn_data.df %>%  summarise( total_remain_genes = sum(remain_genes)) %>% 
        bind_cols(venn_data.df %>% summarise( total_lost_genes = sum(lost_genes))) %>% 
        bind_cols(venn_data.df %>%  summarise( total_gain_genes = sum(gain_genes)))

#Write dataframe of data by chromosome
write.table(venn_data.df, 
            file = "changes_by_chr_genes.tsv", 
            sep = "\t", 
            row.names = F, 
            col.names = T)

#Write dataframe of total data
write.table(Total_venn_data.df, 
            file = "total_changes_genes.tsv", 
            sep = "\t", 
            row.names = F, 
            col.names = T)

