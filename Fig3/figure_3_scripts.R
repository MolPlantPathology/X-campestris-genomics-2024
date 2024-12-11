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

## Documentation
## https://github.com/YuLab-SMU/ggtreeExtra/issues/2
## https://yulab-smu.top/treedata-book/chapter10.html

setwd("/Users/misha/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_3_MGE")


tree <- read.tree("data/core_gene_tree_dna_outgroup.newick") # with outgroup now

# remove duplicates
duplicates <- c("Xcc_barcode38", "Xcr_barcode51")
tree <- drop.tip(tree, duplicates)

# reroot
outgroup_tip <- which(tree$tip.label == "Xfragariae")
tree_root <- phytools::reroot(tree, outgroup_tip, 0.002)
tree_root$edge.length[length(tree_root$edge.length)] = 0.02 # Hacky way to reduce the edge length of the outgroup. TODO: in adobe illustrator, make a little hatch in this edge


## plotting initial tree
p <- ggtree(tree_root, size = 0.45, colour = "grey30") 
p2 <- p + geom_highlight(node = 140, fill = "#E69F00", alpha = 0.5,  extendto = 0.04, to.bottom = TRUE) +
  geom_highlight(node = 175, fill = "#56B4E9", alpha = 0.5,  extendto = 0.04, to.bottom = TRUE) +
  geom_highlight(node = 148, fill = "#56B4E9", alpha = 0.5,  extendto = 0.04, to.bottom = TRUE) +
  geom_highlight(node = 145, fill = "#56B4E9", alpha = 0.5, extendto = 0.04, to.bottom = TRUE) +
  geom_tiplab(align = TRUE, linesize = 0.2, size = 0, offset = 0.011) 

 # load the metadata 
dat <- data.frame(read.csv("data/effectors.csv", header = TRUE))

cols <- c("Xcc" = "#E69F00", "Xcr" = "#56B4E9", "np" = "grey80")
cols2 <- c("Yes" = "black", "No" = "grey80")
cols3 <- c("#E69F00", "#56B4E9", "#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000","#808080","#101010")
cols_tol <- c('#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00', '#DDDDDD', 'black')


### Plasmid
p3 <- p2 + new_scale_fill() + geom_fruit(data = dat,
                                         geom = geom_bar,
                                         mapping = aes(y=id, x=1, fill = Plasmid),
                                         width = 0.75,
                                         pwidth = 0.25, 
                                         offset = 0.045,
                                         stat = 'identity',
                                          ) +
                                          scale_fill_manual(values = cols2)

###### transposons 

####### Transposons per family
trans <- data.frame(read.table("data/IS_per_family.tsv", FALSE))
trans$barcode <- substr(trans$V1, 1, 9)
trans$IS_family <- trans$V2

dat2 <- merge(dat, trans, by.x = "barcode", by.y = "barcode")


### transposons
p4<- p3 + new_scale_fill() + geom_fruit(data = dat2[which(dat2$IS_family != "total"),],
                                  geom = geom_bar,
                                    mapping = aes(y=id, x=V3, fill = IS_family),
                                    width = 0.75,
                                    pwidth = 2.5, 
                                    offset = 0.04,
                                    stat = 'identity', 
                                    axis.param = list(axis="xy", nbreak=5, vjust = 1, text.size = 2.8, line.color = "black"),
                                    grid.params = list(color = "black", size = 0.1),
                                    alpha = 1) +
                                    scale_fill_manual(values = cols_tol) +
                                    ylim(0,98) +
                                    theme(legend.position = "None")

ggsave("IS_per_family_tree_v3.pdf", p4, width = 12, height = 18, units = "cm", dpi = 300)


### Regions of genome plasticity
plastic <-data.frame(read.table("data/plastic_regions_v2.tsv", header = TRUE))

# little function to extract the last n characters from a vector, used to extract barcode from the organism labels
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# count the number of genes in plastic regions, total length in bp, and the number of regions of plasticity.
plastic_summary <- plastic %>% group_by(organism) %>% summarise(n_genes = sum(genes), length=sum(stop - start), n_plastic_regions = n())
plastic_summary$barcode <- substrRight(plastic_summary$organism, 9)

dat <- merge(dat, plastic_summary, by.x = "barcode", by.y = "barcode")

############## comparative boxplots
boxplotter <- function(df, y_axis, xlabel, ylabel)
{
  plot <- ggplot(data = df, aes(y = y_axis, x = factor(pathovar)), environment = environment()) +
    geom_point(aes(colour = factor(pathovar)), alpha = 0.5, shape = 19, size = 1.2, stroke = 0, position = position_jitterdodge()) + 
    geom_boxplot(aes(fill = factor(pathovar)), colour = "black", alpha = 0.7, show.legend = FALSE, outlier.shape = NA, lwd = 0.3, width = 0.5) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) + 
    theme_bw() +
    theme(text = element_text(size = 8, family="sans")) + 
    theme(axis.text = element_text(colour = "black", size = 8)) + 
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.5)) + 
    theme(axis.ticks = element_line(colour = "black")) + 
    theme(panel.grid.minor = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    # scale_y_continuous(breaks = seq(64.5, 65.6, by = 0.2), limits = c(64.5, 65.5), expand = c(0,0)) +
    theme(legend.position = "None") +
    xlab(xlabel) +
    ylab(ylabel) +
    scale_y_continuous(limits = c(0, max(y_axis) + (max(y_axis) * 0.1)))
}


pl1 <- boxplotter(dat[which(dat$pathovar != "np"),], dat[which(dat$pathovar != "np"),]$n_plastic_regions, "Pathovar", "Regions of genome plasticity") 
pl2 <- boxplotter(dat[which(dat$pathovar != "np"),], dat[which(dat$pathovar != "np"),]$n_genes, "Pathovar", "ORFs in regions of genome plasticity") 

collage <- pl1 | pl2

ggsave("ppanggolin_data.pdf", collage, width = 8, height = 5, units = "cm")

wilcox.test(n_plastic_regions ~ pathovar, data = dat[which(dat$pathovar != "np"), ]) # p-value = 1.967e-05
wilcox.test(n_genes ~ pathovar, data = dat[which(dat$pathovar != "np"), ]) # p-value = 2.738e-13











