library(ggplot2)
library(cowplot)
library(readr)
library(patchwork)
library(stringr)
setwd("~/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_S1_DNA_sequencing_assembly")
data_exp <- read_delim("~/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_S1_DNA_sequencing_assembly/data/mastertable_samplenames_barcodes_read_statistics.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

data_exp$coverage = as.numeric(str_replace(data_exp$total.gigabases, ",",".")) / 0.005

colours = c("Xcc" = "#E69F00", "Xcr" = "#56B4E9", "Xc" = "black")

# remove duplicate samples
to_remove <- c("barcode38", "barcode51")
data_exp <- data_exp %>% filter(!(barcode %in% to_remove)) 



p1 <- ggplot(data_exp, aes(x = N50.length, y = coverage)) +
  geom_point(aes(colour = pathovar), shape = 19, size = 1.2, alpha = 1, stroke = 0) +
  xlab("N50 read length (bp)") +
  ylab("Coverage (x genome)") + 
  theme_bw() +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) + 
  theme(axis.ticks = element_line(colour = "black")) + 
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  scale_y_continuous(breaks = seq(0, 700, by = 100), limits = c(0, 710), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0, 45000, by = 5000), limits = c(0,47000), expand = c(0,0)) + coord_flip() +
  scale_colour_manual(values = colours) +
  theme(plot.margin = unit(c(-1, -1.2, -1.2, -1.5), "cm")) + 
  theme(legend.position = "None") +
  annotate("text", x = 42500, y = 260, label = "Xc", size = 2) +
  theme(axis.text.y = element_text(size = 8, family="sans",colour = "black"))+
  theme(axis.text.x = element_text(size = 8, family="sans",colour = "black")) +
  theme(axis.title.y = element_text(size = 8, family="sans",colour = "black"))+
  theme(axis.title.x = element_text(size = 8, family="sans",colour = "black"))
       
p1

# density plots with the same data
x_plot <- ggplot(data_exp, aes(y = N50.length)) + geom_density(aes(fill = pathovar, colour = pathovar), alpha = 0.5, lwd = 0.3) +
          theme_void() + scale_y_continuous(breaks = seq(0, 45000, by = 5000), limits = c(0,47000), expand = c(0,0)) +
          scale_fill_manual(values = colours) +
          scale_colour_manual(values = colours) +
          theme(legend.position = "None")

x_plot

y_plot <- ggplot(data_exp, aes(y = coverage)) + geom_density(aes(fill = pathovar, colour = pathovar), alpha = 0.5, lwd = 0.3) +
          theme_void() +
          scale_y_continuous(breaks = seq(0, 700, by = 50), limits = c(0, 710), expand = c(0,0)) +
          coord_flip() +
          scale_fill_manual(values = colours) + 
          scale_colour_manual(values = colours) + 
          theme(legend.position = "None") +
          annotate("text", x = 0.0041, y = 230, label = "Xcr", size = 2, colour = "#56B4E9") +
          annotate("text", x = 0.003, y = 380, label = "Xcc", size = 2, colour = "#E69F00")   
y_plot

# assemble the plots using patchwork
combined_plots <- y_plot + plot_spacer() + p1 + x_plot + plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

combined_plots

# save
ggsave("readdescriptives_v2.pdf", combined_plots, width = 8, height = 8, unit = "cm", dpi = 300)


#### BUSCO scores

data <- read.table("~/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_S1_DNA_sequencing_assembly/data/combined_busco_homopolish.tsv", quote="\"", comment.char="")
### TODO: get new version here with homopolish data

# some parsing
data$intermediate <- sapply(strsplit(as.character(data$V1),"/"), "[", 1)
data$method <- sapply(strsplit(as.character(data$intermediate), "-"), "[", 1)
data$intermediate2 <- sapply(strsplit(as.character(data$V1), "-"), "[", 2)
data$barcode <- sapply(strsplit(as.character(data$intermediate2), "/"), "[", 1)

data$V1 = NULL
data$intermediate = NULL

# make percentual
data$Proportion_BUSCOs = data$V2 / 1152
# make percentual
data$Proportion_BUSCOs = data$V2 / 1152 * 100

# subset for complete busco's
subset_complete <- data[which(data$V3 == "Complete"),]

# add the barcode info to the busco dataframe (together with a lot of other info that's not needed and not used)
combined <- merge(subset_complete, data_exp, by.x = "barcode", by.y = "barcode")

## TODO: remove not pathogenic isolate
combined <- combined[which(combined$barcode != "barcode90"),]

combined$method <- factor(combined$method, levels = c("flye", "medaka", "homopolish"))

# remove duplicate samples
to_remove <- c("barcode38", "barcode51")
combined <- combined %>% filter(!(barcode %in% to_remove)) 

plot <- ggplot(data = combined, aes(y = Proportion_BUSCOs, x = method, fill = factor(pathovar))) +
        geom_hline(yintercept = 99.6, colour = "grey80", linetype = "dashed") +
        #geom_hline(yintercept = 100, colour = "black")+ ## TODO: Add BUSCO score of reference assembly here
        geom_point(aes(colour = factor(pathovar)), alpha = 0.5, shape = 19, size = 1.2, stroke = 0, position = position_jitterdodge()) + 
        geom_boxplot(aes(colour = factor(pathovar)), colour = "black", alpha = 0.7, show.legend = FALSE, outlier.shape = NA, lwd = 0.3) +
        scale_fill_manual(values = colours) +
        scale_colour_manual(values = colours) + 
        theme_bw() +
        theme(text = element_text(size = 10, family="sans")) + 
        theme(axis.text = element_text(colour = "black")) + 
        theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.5)) + 
        theme(axis.ticks = element_line(colour = "black")) + 
        theme(panel.grid.minor = element_blank()) +
        theme(panel.grid.major.x = element_blank()) +
        scale_y_continuous(breaks = seq(80, 100, by = 5), limits = c(79, 102), expand = c(0,0)) +
        theme(legend.position = "None") +
        xlab("Assembly/polishing stage") +
        ylab("Genome completeness\n(% complete BUSCO genes)") +
        theme(axis.text.y = element_text(size = 8, family="sans",colour = "black"))+
        theme(axis.text.x = element_text(size = 8, family="sans",colour = "black")) +
        theme(axis.title.y = element_text(size = 8, family="sans",colour = "black"))+
        theme(axis.title.x = element_text(size = 8, family="sans",colour = "black"))
        
plot

ggsave("BUSCO_scores_v2.pdf", plot, width = 8, height = 8, unit = "cm", dpi = 300)

wilcox.test(Proportion_BUSCOs ~ pathovar, data = combined, subset = method == 'flye') # p-value = 0.1179
wilcox.test(Proportion_BUSCOs ~ pathovar, data = combined, subset = method == 'medaka') # p-value = 0.08098
wilcox.test(Proportion_BUSCOs ~ pathovar, data = combined, subset = method == 'homopolish')  #  p-value = 0.05001




## GC content etc

GC <- read.table("~/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_S1_DNA_sequencing_assembly/data/genome_size_GC_content.tsv", quote="\"", comment.char="", header = TRUE)
CDS <- read.table("~/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_S1_DNA_sequencing_assembly/data/prokka_stats.tsv", quote="\"", comment.char="", header = TRUE)

        
combined <- merge(GC, CDS, by.x = "barcode", by.y = "barcode")
combined <- merge(combined, data_exp, by.x = "barcode", by.y = "barcode")

combined$GC_content <- combined$GC_content * 100

combined <- combined[which(combined$barcode != "barcode90"),]


boxplotter <- function(df, y_axis, xlabel, ylabel)
{
  plot <- ggplot(data = df, aes(y = y_axis, x = factor(pathovar)), environment = environment()) +
    geom_point(aes(colour = factor(pathovar)), alpha = 0.5, shape = 19, size = 1.2, stroke = 0, position = position_jitterdodge()) + 
    geom_boxplot(aes(fill = factor(pathovar)), colour = "black", alpha = 0.7, show.legend = FALSE, outlier.shape = NA, lwd = 0.3, width = 0.5) +
    scale_fill_manual(values = colours) +
    scale_colour_manual(values = colours) + 
    theme_bw() +
    theme(text = element_text(size = 8, family="sans")) + 
    theme(axis.text.y = element_text(size = 8, family="sans",colour = "black"))+
    theme(axis.text.x = element_text(size = 8, family="sans",colour = "black")) +
    theme(axis.title.y = element_text(size = 8, family="sans",colour = "black"))+
    theme(axis.title.x = element_text(size = 8, family="sans",colour = "black")) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.5)) + 
    theme(axis.ticks = element_line(colour = "black")) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    # scale_y_continuous(breaks = seq(64.5, 65.6, by = 0.2), limits = c(64.5, 65.5), expand = c(0,0)) +
    theme(legend.position = "None") +
    xlab(xlabel) +
    ylab(ylabel) +
    scale_y_continuous(expand = c(0,0))
}


pl1 <- boxplotter(combined, combined$GC_content, "Pathovar", "GC content (%)") 
pl2 <- boxplotter(combined, combined$genome_size/1000000, "Pathovar", "Genome size (Mbp)")
pl3 <- boxplotter(combined, combined$tRNA, "Pathovar", "tRNAs")
pl4 <- boxplotter(combined, combined$CDS, "Pathovar", "Coding sequences")


wilcox.test(GC_content ~ pathovar, data = combined) # p-value = 9.161e-16
wilcox.test(genome_size ~ pathovar, data = combined) # p-value < 2.2e-16
wilcox.test(CDS ~ pathovar, data = combined)  # p-value = 7.128e-13


#### save collage plot.
library(patchwork)
#arranged <- plot_grid(pl1, pl2, pl3, pl4,labels = c('E', 'F', 'G', 'H'), align  = 'v')  # no longer used
#arranged

arranged <- (plot | pl1 | pl2 | pl4) + plot_layout(widths = c(1.5, 1,1,1))
ggsave("genome_characteristics_v2.pdf", arranged, width = 19, height = 6, unit = "cm", dpi = 300)


################
# Phylogenetic tree with pathovar incanae genomes
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(readr)
library(phytools)


setwd("~/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_S1_DNA_sequencing_assembly")
tree <- read.tree("data/core_genome_tree_plus_incanae.newick")

## remove duplicate samples
duplicates <- c("Xcc_barcode38", "Xcr_barcode51")
tree <- drop.tip(tree, duplicates)

# reroot
outgroup_tip <- which(tree$tip.label == "Xfragariae")
tree_root <- phytools::reroot(tree, outgroup_tip, 0.002)
tree_root$edge.length[length(tree_root$edge.length)] = 0.02

# read name table
names <- read.csv("data/name_table.csv")

colours <- c("Xcc" = "#E69F00", "Xcr" = "#56B4E9", "np" = "grey70", "Xf"= "#D55E00", "Xci" = "#CC79A7")

t1 <- ggtree(tree_root, linewidth = 0.5) + geom_tiplab(size = 0, linesize = 0.3, align = TRUE)

t2 <- t1 + geom_fruit(data = names, geom = geom_bar, mapping = aes(y=tip_label, x = 0.5, fill = pathovvar), stat = 'identity', width = 0.8, alpha = 0.8, pwidth = 0.1) + scale_fill_manual(values = colours)
t3 <- t2 + geom_fruit(data = names, geom = geom_text, mapping = aes(y = tip_label, label = pathovar_strainname), hjust = "left", size = 1.9, pwidth = 10) + theme(legend.position = "bottom") + xlim(0,0.08)

ggsave("tree_with_incanae_v2.pdf", t3, width = 10, height = 19, unit = "cm", dpi = 300)





