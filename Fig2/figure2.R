library(ggtree)
library(phytools)
library(ggplot2)
library(tidyr)
library(ggtreeExtra)
library(ggnewscale)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(patchwork)
library(stringr)

## Documentation
## https://github.com/YuLab-SMU/ggtreeExtra/issues/2
## https://yulab-smu.top/treedata-book/chapter10.html

setwd("/Users/misha/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_2_T3E")

# load the metadata 
dat <- data.frame(read.csv("metadata.csv", header = TRUE))

# read the tree 
tree <- read.tree("core_gene_tree_dna_outgroup.newick") # with outgroup now

## remove duplicate samples
duplicates <- c("Xcc_barcode38", "Xcr_barcode51")
tree <- drop.tip(tree, duplicates)

# reroot at outgroup
outgroup_tip <- which(tree$tip.label == "Xfragariae")
tree_root <- phytools::reroot(tree, outgroup_tip, 0.002)
tree_root$edge.length[length(tree_root$edge.length)] = 0.02 # Hacky way to reduce the edge length of the outgroup. TODO: in adobe illustrator, make a little hatch in this edge

## plotting initial tree
p <- ggtree(tree_root, size = 0.3, colour = "grey30") 
p2 <- p + geom_highlight(node = 140, fill = "#E69F00", alpha = 0.5,  extendto = 0.04, to.bottom = TRUE) +
    geom_highlight(node = 175, fill = "#56B4E9", alpha = 0.5,  extendto = 0.04, to.bottom = TRUE) +
    geom_highlight(node = 148, fill = "#56B4E9", alpha = 0.5,  extendto = 0.04, to.bottom = TRUE) +
    geom_highlight(node = 145, fill = "#56B4E9", alpha = 0.5, extendto = 0.04, to.bottom = TRUE) +
    geom_tiplab(align = TRUE, linesize = 0.2, size = 0, offset = 0.005)

p2

colours <- c("Xcc" = "#E69F00", "Xcr" = "#56B4E9")
colours2 <- c("1" = "#D55E00", "0" = "Grey90")

## read T3E data
t3e <- data.frame(read.csv("T3E_binary.csv", header = TRUE))

# subset the presence absence data
t3ex <- t3e[,1:71] %>% gather(T3E, presence_absence, -id)
t3ex <- merge(t3ex, dat, by.x = 'id', by.y = 'barcode')

## This code makes a new variable called t3e_class where different alleles/copies of the same effector are summarised
t3ex <- t3ex %>% mutate(t3e_class = case_when(
                                    str_detect(T3E, "TALE") ~ "TALE",
                                    str_detect(T3E, "XopAL") ~ "XopAL",
                                    str_detect(T3E, "XopAZ") ~ "XopAZ",
                                    str_detect(T3E, "XopD") ~ "XopD",
                                    str_detect(T3E, "XopE") ~ "XopE",
                                    str_detect(T3E, "XopH") ~ "XopH",
                                    str_detect(T3E, "XopK") ~ "XopK",
                                    str_detect(T3E, "XopL") ~ "XopL",
                                    str_detect(T3E, "XopM") ~ "XopM",
                                    str_detect(T3E, "XopR") ~ "XopR",
                                    str_detect(T3E, "XopX") ~ "XopX",
                                    str_detect(T3E, "XopZ") ~ "XopZ",
                                    str_detect(T3E, "AvrBs2") ~ "AvrBs2",
                                    str_detect(T3E, "AvrXcc") ~ "AvrXcc",
                                    str_detect(T3E, "Rip") ~ "Rip",
                                    str_detect(T3E, "Hop") ~ "Hop",
                                    TRUE ~ T3E  # If none of the conditions are met
                                  ))

# This is to define different categories, used to sort the heatmap in our desired order
core <- c("XopA_Harpin", "XopAZ", "XopP", "XopF1", "AvrXcc", "XopAL", "XopR", "XopM")
core_xcr <- c("XopAD")
var_both <- c("XopAC","XopZ", "XopE")
core_xcc_var_xcr <- c("XopL", "XopN")
core_xcc <- c("XopQ", "XopAG", "XopAM", "XopAY", "AvrBs2","XopK","XopX")
var_xcc <- c("XopD", "XopH", "XopJ5", "AvrBs1", "XopAH", "XopG1", "XopAE", "TALE")
misc <- c("Rip", "Hop")

new_order <- c(core, core_xcr, var_both,  core_xcc_var_xcr, core_xcc, var_xcc, misc)
t3ex$t3e_class <- factor(t3ex$t3e_class, levels = new_order)

## remove duplicate samples
duplicates <- c("barcode38", "barcode51")
t3ex <- t3ex %>% filter(!(id %in% duplicates))

## we make the heatmap seperate from the phylogenetic tree
tree_order <- get_taxa_name(p)
tree_order_clean <- sub(".*_", "", tree_order)
t3ex$id <- factor(t3ex$id, levels = rev(tree_order_clean))

new_order2 <- c("Xcc","Xcr","np")

t3ex$pathovar <- factor(t3ex$pathovar, levels = new_order2)


heatmap <- ggplot(data = t3ex, aes(y= id, x = T3E, fill = as.factor(presence_absence))) +
                          geom_tile(linewidth = 0.2, colour = "white") + 
                          facet_grid(pathovar~t3e_class, scales = "free", space = "free") + ## This is key: facetting by pathovar/class with scales and space "free" gives us nicely seperated blocks
                          theme_bw() + 
                          theme(panel.spacing = unit(0.06, "cm")) +
                          scale_fill_manual(values = colours2) +
                          theme(legend.position = 'none')  +
                          theme(
                            axis.text.y = element_blank(),
                            strip.text = element_blank(),        # Hide the text labels
                            strip.background = element_blank()) + # Hide the background
                          theme(
                            axis.text.x = element_text(
                              size = 7,         # Font size
                              color = "black",    # Font color
                              angle = 90,        # Text angle (rotation)
                              hjust = 1,          # Horizontal justification (0.5 = center)
                              vjust = 0.5)) +  # Vertical justification (0.5 = center)
                            labs(x = "Putative type III effectors", y = NULL) +
                            theme(panel.border = element_rect(color = "black", linewidth = 0.3)) +
                            theme(axis.ticks = element_line(linewidth = 0.2))


library(patchwork)
collage <- (p2| heatmap) + plot_layout(widths = c(1,8))                       
                                                       
ggsave("heatmap_t3E_reordered_v2.pdf",collage, width = 26, height = 16, units = "cm", dpi = 300)

# Because we introduced a small separation between Xcc and Xcr in the heatmap, we need to manually adjust the phylogenetic tree slightly in adobe illustrator!










