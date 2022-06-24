# Cargamos pacman
library( "pacman" )

# cargamos mas paquetes
p_load( "vroom",
        "dplyr",
        "ggplot2",
        "purrr",
        "stringr",
        "cowplot",
        "tidyr",
        "scales")

# Cargamos todos los tsv en un solo dataframe
filenames <- list.files( pattern="*.tsv" )

changes_data.df <- map_df( filenames, vroom)


#Write dataframe of data by chromosome
write.table(changes_data.df, 
            file = "all_percent_changes.tsv", 
            sep = "\t", 
            row.names = F, 
            col.names = T)

histogram <- changes_data.df %>%  
  ggplot( aes(x = percent, fill = target))  +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = c("percent_lost" = "#c87570",
                               "percent_gain" = "#70c875",
                               "percent_remain" = "#7570c8")) +
  labs( x = "percent of miRNA/targets pairs changes") +
  scale_x_continuous(labels = label_percent()) + labs(fill="") +
  theme_minimal_grid()

ggsave( filename = "miRnome_percent_histogram.png", 
        plot = histogram,
        device = "png",
        height = 7, width = 14,
        units = "in")

# Density plot
Density <- ggplot(data=changes_data.df, aes(x=percent, group=target, fill=target)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs( x = "percent of miRNA/targets pairs changes") +
  scale_fill_manual(values = c("percent_lost" = "#c87570",
                               "percent_gain" = "#70c875",
                               "percent_remain" = "#7570c8")) +
  scale_x_continuous(labels = label_percent()) +
  theme_minimal_grid()

ggsave( filename = "miRnome_density.png", 
        plot = Density,
        device = "png",
        height = 7, width = 14,
        units = "in")



