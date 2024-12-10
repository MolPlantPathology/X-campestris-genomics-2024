setwd("/Users/misha/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_8_CRISPR_target_plasmids") # replace with you working directory

# getting and connecting all the data
matrix <- read.csv('xop_selection_presence_absence.csv')
species_list <- read.csv('name_list.csv')
long_df <- matrix %>% gather(Plasmid, Presence_absence, -c(Class, Gene, Annotation))
long_df <- merge(long_df, species_list, by = "Plasmid")

# reordering facets
new_order <- c("X. campestris", "X. hortorum", "X. citri", "X. phaseoli", "other")
long_df$species_simplified <- factor(long_df$species_simplified , levels = new_order)

new_order2 <- c("Plasmid", "T3E (Plasmid encoded)", "T3E (Core)")
long_df$Class <- factor(long_df$Class, levels = new_order2)

# read the tree (or in this case, rather a hierarchical clustering) to reorder the plasmids according to the clustering
tree <- read.tree("accessory_binary_genes.fa.newick") # with outgroup now
p <- ggtree(tree, size = 0.3, colour = "grey30") 
tree_order <- get_taxa_name(p)
long_df$Plasmid <- factor(long_df$Plasmid, levels = tree_order)

# setting the colors
colours2 <- c("1" = "#D55E00", "0" = "Grey90")

plot <- ggplot(data = long_df, aes(x = Plasmid, y = Gene, fill = as.factor(Presence_absence))) + geom_tile(colour = "white", linewidth = 0.1) + facet_grid(Class ~ species_simplified, scales = "free", space = "free") +
  theme_bw() +
  theme(panel.spacing = unit(0.04, "cm")) + 
  scale_fill_manual(values = colours2) +
  theme(legend.position = "none") + 
  theme(text = element_text(size = 8, family="sans",colour = "black")) +  
  theme(axis.text.y = element_text(size = 8, family="sans",colour = "black")) +
  theme(axis.text.x = element_blank()) + 
  labs(x = "Plasmids targetted by Xcr spacers", y = "") +
  theme(axis.ticks = element_line(linewidth = 0.25)) +
  theme(strip.text.y = element_text(size = 8, family="sans",colour = "black"))+
  theme(strip.text.x = element_text(size = 8, family="sans",colour = "black", face = "italic"))
  
# Fig S5
ggsave("panelA_heatmap.pdf", plot, width = 19, height = 13, units = "cm", dpi = 300)



# insertion seqeunces per KB calculation
library(readr)
IS <- read.delim("plasmids_IS_data.tsv", header = FALSE, na.strings = c("NA"," "))
IS_count <-IS %>% group_by(V1) %>% tally()

plasmid_length <- read.table('contig_lengths.tsv', header = TRUE) %>% select(c(file, min_len))
plasmid_length_species <- merge(plasmid_length, species_list, by.x = "file", by.y = "Plasmid") 

all_data <- merge(plasmid_length_species, IS_count, by.x="file", by.y="V1", all = TRUE) %>%  mutate(n = replace_na(n, 0))
all_data$IS_per_10kb <- all_data$n / all_data$min_len * 10000

## compare this measure to Xcc8004 and Xcr756
xcc <- 85 / 5150550 * 10000
xcr <- 49 / 4941800 * 10000

new_order <- c("X. campestris", "X. hortorum", "X. citri", "X. phaseoli", "other")
all_data$species_simplified <- factor(all_data$species_simplified , levels = new_order)

plot2 <- ggplot(data = all_data, aes(x = species_simplified, y = IS_per_kb, colour = species_simplified)) + 
          geom_boxplot(aes(fill = species_simplified), alpha = 0.3, colour = "black", outlier.shape = NA, width = 0.5) + geom_jitter(height=  0, width = 0.1, alpha = 0.7) +
          theme_bw() + 
          labs(x = "Plasmids targetted by Xcr spacers (by species)", y = "Number of ISs per 10.000 bp") +
          scale_fill_brewer(palette = "Set2") + 
          scale_colour_brewer(palette = "Set2") + 
          theme(legend.position = "none") +
          theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
          theme(axis.text.y = element_text(size = 8, family="sans",colour = "black"))+
          theme(axis.text.x = element_text(size = 8, family="sans",colour = "black", face = "italic")) +
            theme(axis.title.y = element_text(size = 8, family="sans",colour = "black"))+
            theme(axis.title.x = element_text(size = 8, family="sans",colour = "black"))+
          geom_hline(yintercept = xcc, linewidth = 0.5, linetype = "dashed", color = "#E69F00") +
          geom_hline(yintercept = xcr, linewidth = 0.5, linetype = "dashed", color = "#56B4E9")

ggsave("panelB_boxplot.pdf", plot2, width = 12, height = 7, units = "cm", dpi = 300)