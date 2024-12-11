library(ggplot2)
library(tidyr)
library(dplyr)
library(here)

theme_mp <- theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  theme(strip.text = element_text(size = 6,  margin = margin(t = 1.5, b = 1.5))) +
  theme(strip.background = element_rect(fill = "grey90")) +
  theme(axis.text=element_text(colour="black")) +
  theme(legend.key.size = unit(0.1, "cm")) +
  #theme(legend.position = "right") +
  theme(text = element_text(size = 10, family="sans", colour = "black")) 

here()
setwd("/Users/misha/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_1_v2_XC_collection")
bioassays <- read.csv(here("Fig1/data/M28M29M30_simplified.csv"), sep=";")

colours = c("Xcc" = "#E69F00", "Xcr" = "#56B4E9", "not_pathogenic" = "grey50", "not_available" = "grey90")

bioassays_long <- bioassays %>% gather(key = "Bioassay", value = "Result",starts_with("M")) %>% mutate(number = reorder(number, -evidence_level)) 
custom_order <- c("M28", "M29", "M30")
bioassays_long$Bioassay <- factor(bioassays_long$Bioassay, levels = custom_order)

custom_order2 <- c("not_available", "not_pathogenic", "Xcr", "Xcc")
bioassays_long$Result <- factor(bioassays_long$Result, levels = custom_order2)

custom_order3 <- c("Xcc", "Xcr", "not_pathogenic")
bioassays_long$pathovar_designation <- factor(bioassays_long$pathovar_designation, levels = custom_order3)

## all strains, vertical heatmap edition
heatmap <- ggplot(data = bioassays_long) + geom_tile(aes(x = Bioassay, y = number, fill = Result), linewidth = 0.2, colour = "white") +
                  facet_grid(pathovar_designation ~ Bioassay, scales = "free", space = "free") + 
                  scale_fill_manual(values = colours, breaks = rev(levels(bioassays_long$Result))) +
                  theme_bw() +
                  theme(panel.grid.major = element_blank()) +
                  scale_y_discrete(labels = NULL) +
                  labs(x = "Experiment", y = "Strains", fill = "Phenotype:") +
                  theme(axis.ticks.x = element_blank()) +
                  theme(axis.ticks.y = element_line(linewidth = 0.4)) +
                  theme(strip.text = element_blank(),      
                  strip.background = element_blank()) +
                  theme(axis.text=element_text(colour="black")) +
                  theme(text = element_text(size = 10, family="sans", colour = "black")) +
                   theme(legend.position = "none")

ggsave("bioassays_heatmap.pdf", heatmap, width = 5, height = 22, units = "cm", dpi = 300)

########## tree
library(ggtree)
library(ape)
library(ggtreeExtra)
library(phytools)
library(ggnewscale)

tree <- read.tree("core_gene_tree_dna_outgroup.newick")
metadata <- read.csv('data/sequenced_isolates_metadata.csv', na.strings=c(""," ","NA"))

## remove duplicate samples
duplicates <- c("Xcc_barcode38", "Xcr_barcode51")
tree <- drop.tip(tree, duplicates)

outgroup_tip <- which(tree$tip.label == "Xfragariae")
tree_root <- phytools::reroot(tree, outgroup_tip, 0.002)
tree_root$edge.length[length(tree_root$edge.length)] = 0.02 

######### Half fan
p <- ggtree(tree_root, layout = "fan", open.angle = 180, size = 0) +
                    geom_highlight(node = 140, fill = "#E69F00", alpha = 0.5,  extendto = 0.04) +
                    geom_highlight(node = 175, fill = "#56B4E9", alpha = 0.5,  extendto = 0.04) +
                    geom_highlight(node = 148, fill = "#56B4E9", alpha = 0.5,  extendto = 0.04) +
                    geom_highlight(node = 145, fill = "#56B4E9", alpha = 0.5, extendto = 0.04) + 
                    geom_tree(size = 0.4) +
                    geom_tiplab2(align = TRUE, size = 0, offset = 0.004, linesize = 0.3, linetype = "solid", colour = "grey80") +
                    geom_fruit(data = metadata, 
                              geom = geom_bar,
                              aes(y = pathovar_id, x = 2, fill = continent), 
                              colour = "black",
                              size = 0.1,
                              pwidth = 0.05,
                              offset = 0.047,
                              stat = 'identity', 
                              orientation = 'y') +
                  scale_fill_brewer(palette = "Set2", na.value = "grey80",
                                    guide = guide_legend(frame.colour = "black", barwidth = 2, barheight = 0.5, ticks.colour = "black")) +
                  new_scale_fill() +
                  geom_fruit(data = metadata, 
                              geom=geom_bar, 
                              aes(y = pathovar_id, x = 2, fill = year), 
                              colour = "black", size = 0.1,
                              pwidth = 0.05,
                              offset = 0.02,
                              stat = 'identity', 
                              orientation = 'y') +
                    theme(legend.position = "bottom", 
                          legend.text=element_text(size = 8), 
                          legend.title = element_text(size = 8)) +
                    scale_fill_distiller(palette = "Greens", direction = 1, na.value = "grey80", 
                                         guide = guide_colorbar(frame.colour = "black", barwidth = 4, barheight = 0.5, ticks.colour = "black", label.theme=element_text(angle = -90, size = 6, vjust = 0.5)))

                    # geom_text(aes(label=node), hjust=0, size = 0.5) # if you want to know internal node numbers

        
ggsave("haf_fan_with_metadata_v3.pdf", p, width = 15, height = 7, units = "cm", dpi = 300)


########## Number of T3E plot (panel E of Figure 1)

dat <- read.csv('data/effectors.csv')

# remove duplicate barcodes
duplicates <- c("barcode38", "barcode51")
dat <- dat %>% filter(!(barcode  %in% duplicates))

plot <- ggplot(data = dat, aes(x = pathovar, y = n_T3E)) +
                  geom_point(aes(colour = factor(pathovar)), alpha = 0.5, shape = 19, size = 1.2, stroke = 0, position = position_jitterdodge()) + 
                  geom_boxplot(aes(fill = factor(pathovar)), colour = "black", alpha = 0.7, show.legend = FALSE, outlier.shape = NA, lwd = 0.3, width = 0.7) +
                  theme_mp + 
                  scale_fill_manual(values = colours) +
                  scale_color_manual(values = colours) + 
                  theme(legend.position = "none") + 
                  ylim(c(0,40)) +
                  labs(x = "Pathovar", y = "Type III effectors") +
                  theme(text = element_text(size = 8, family="sans", colour = "black")) +
                  theme(axis.text.x = element_text(size = 8, family="sans", colour = "black")) +
                  theme(axis.text.y = element_text(size = 8, family="sans", colour = "black")) 

ggsave("T3E_plot.pdf", plot, width = 5, height = 5, units = "cm", dpi = 300)
                
dat %>% filter(pathovar == "Xcc") %>% summarise(m = mean(n_T3E))
dat %>% filter(pathovar == "Xcr") %>% summarise(m = mean(n_T3E))

dat <- dat %>% filter(pathovar != "np") 

wilcox.test(n_T3E ~ pathovar, data = dat)
t.test(n_T3E ~ pathovar, data = dat)

colours <- c("Xcc" = "#E69F00", "Xcr" = "#56B4E9", "np" = "black")
 
 

