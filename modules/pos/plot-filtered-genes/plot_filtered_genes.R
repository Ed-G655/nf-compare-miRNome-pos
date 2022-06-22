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

paleta <- c("lost" =  "#F94144",
            "gained" = "springgreen3") 

piramide.p <- ggplot(changes_data.df, aes(x = miRNA_ID, y = percent, fill = target )) + 
  geom_col(data = subset(changes_data.df, target == "percent_lost"), 
           width = 0.5 , fill = "#F94144", color = "black", alpha = 0.5) + 
  geom_col(data = subset(changes_data.df, target ==  "percent_gain"), 
           width = 0.5 , fill = "springgreen3", color = "black", alpha = 0.5) +
  coord_flip() +  scale_y_continuous(labels = label_percent()) +
  labs(y= "Numero de pares miRNA/blanco", x = "miRNA", color = "Legend") +
  scale_color_manual(values = paleta) +
  labs(title = "Genes blanco por miRNA y sus cambios debido a mutaciones en el miRNA") +
  theme_minimal_grid() + theme(axis.text.y  = element_text(face="bold", size=5, angle= 30))

ggsave( filename = str_interp("miRnome_percent_filtered.png"), 
        plot = piramide.p,
        device = "png",
        height = 7, width = 14,
        units = "in")

#Write dataframe of data by chromosome
write.table(changes_data.df, 
            file = "changes_filtered.tsv", 
            sep = "\t", 
            row.names = F, 
            col.names = T)
