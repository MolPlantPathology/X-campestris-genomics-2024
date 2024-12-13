library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(stringr)
library(patchwork)
library(readxl)
library(ggbeeswarm)
library(RColorBrewer)

theme_mp <- theme_bw() + theme(axis.text = element_text(color = "black", size = 6)) +
  theme(axis.title = element_text(color = "black", size = 8)) +
  theme(panel.grid.major = element_line(linewidth = 0.3), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) + 
  theme(panel.background = element_rect(colour = "black")) +
  theme(axis.ticks = element_line(colour = "black")) + 
  theme(strip.background=element_rect(fill="grey95"), strip.text = element_text(size = 8, face = 'bold')) 


setwd("/Users/misha/surfdrive/03 Xanthomonas genomics/Genomics/2023_transposons-T3E")

df <- read_tsv("distance_total.tsv", col_names = FALSE)
metadata <- read_tsv("metadata_AtH.txt", col_names = TRUE)

# some data cleanup
df <- df %>% mutate(effector = str_split(X4,";", simplify = T)[,2])
df <- df %>%  mutate(effector_clean = str_replace(effector,"Name=","")) %>% filter(!grepl("inference=",effector_clean))

# connect pathovar metadata
df_metadata <- merge(df, metadata, by.x = "X10", by.y="barcode")
colours = c("Xcc" = scales::alpha("#E69F00",0.3), "Xcr" = scales::alpha("#56B4E9", 0.3))

df_metadata$distance_adjusted <- df_metadata$X9 + 1  ## to enable log transformation
df_metadata <- df_metadata %>% filter(effector_clean != "CM002755.1.898")

## remove duplicates
duplicates <- c("barcode38", "barcode51")
df_metadata <- df_metadata %>% filter(!(X10 %in% duplicates))


plot <- ggplot(data = df_metadata %>% filter(pathovar != 'Not_pathogenic'), aes(x = reorder(effector_clean, distance_adjusted, FUN = median), y = distance_adjusted)) + 
                annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = 0, ymax = 2500, colour = "grey85", alpha = 0.1, linewidth = 0) +
                geom_boxplot(aes(fill = pathovar), alpha = 0.4, outlier.shape = NA) + 
                geom_jitter(aes(colour = pathovar), width = 0.3, height = 0, alpha = 0.3, size = 1) +
                theme_mp + 
                theme(axis.text.x = element_text(angle = 45 , vjust = 1, hjust=1)) + 
                theme(axis.text.y = element_text(angle = 0 , vjust = 0.5, hjust=1)) +
                scale_y_continuous(labels = scales::label_comma(), trans = "log2") +
                facet_wrap(.~pathovar)  + 
                coord_flip() + 
                scale_color_manual(values = colours) +
                scale_fill_manual(values = colours) +
                theme(legend.position = "none") + 
                labs(y = "Distance to nearest insertion sequence (bp)", x = "") 


ggsave("distances_boxplt_v2.pdf", plot, width = 10, height = 18, units = "cm", dpi = 300)



plot <- ggplot(data = df, aes(x = effector, y = X9)) + geom_boxplot() + theme(axis.text.x = element_text(angle =90, vjust = 0.5, hjust=1))
plot <- ggplot(data = df_metadata, aes(x = X9)) + geom_density(aes(color = pathovar)) + facet_grid(pathovar ~ effector_clean, scales = "free")
 
ggsave("small_multiples.png", plot, width = 1000, height = 260, units = "mm", dpi = 300)

