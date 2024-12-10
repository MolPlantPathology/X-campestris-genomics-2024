library(ggtree)
library(ggtreeExtra)
library(readr)
library(phytools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

# common theme for plots
theme_mp <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  theme(strip.text = element_text(size = 6,  margin = margin(t = 1.5, b = 1.5))) +
  theme(strip.background = element_rect(fill = "grey90")) +
  theme(axis.text=element_text(colour="black")) +
  theme(legend.key.size = unit(0.1, "cm")) +
  theme(legend.position = "right") +
  theme(text = element_text(size = 8, family="sans", colour = "black")) 

theme_mp_stacked_bar <- theme_bw() +
  theme(strip.text = element_text(size = 6,  margin = margin(t = 1.5, b = 1.5))) +
  theme(strip.background = element_rect(fill = "grey90")) +
  theme(axis.text=element_text(colour="black")) +
  theme(legend.key.size = unit(0.1, "cm")) +
  theme(legend.position = "right") +
  theme(text = element_text(size = 8, family="sans", colour = "black")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# read tree
setwd("/Users/misha/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_5_crispr")
tree <- read.tree("core_gene_tree_dna_outgroup.newick")

## remove duplicate samples
duplicates <- c("Xcc_barcode38", "Xcr_barcode51")
tree <- drop.tip(tree, duplicates)

outgroup_tip <- which(tree$tip.label == "Xfragariae")
tree_root <- phytools::reroot(tree, outgroup_tip, 0.002)
tree_root$edge.length[length(tree_root$edge.length)] = 0.02

ggtree(tree_root)

# read barcode metadata from two files, merge this data
datx <- data.frame(read.csv("effectors.csv", header = TRUE))
dat2 <- data.frame(read.table('sequenced_isolates_metadata_taqman.txt', header = TRUE))
dat3 <- merge(dat2, datx, by.x = "id", by.y="barcode")

# read CRISPRCasFinder results --> filtered on evidence level 4 arrays
CRISPRS <- data.frame(read_tsv("results_with_outgroup.tsv", col_names = TRUE))

# read CRISPR presence absence data
CRISPRS_PA <- read.csv("crispr_presence_absence_outgroup.csv")

CRISPRx <- merge(CRISPRS, CRISPRS_PA, by = "barcode")

# connect metadata with CRISPR data table
dat4 <- merge(dat3, CRISPRx, by.x="id", by.y = "barcode", all.y = TRUE)
dat4[dat4$id == "Xfragariae", "id.y"] <- "Xfragariae"

dat5 <- dat4 %>%  pivot_longer(cols = c("CRISPR.Cas", "XopN", "IS_in_spacer_array"), names_to = "Gene", values_to="Presence_absence")
dat5$Gene <- factor(dat5$Gene, levels = c("CRISPR.Cas", "XopN", "IS_in_spacer_array"))
dat5$Presence_absence <- factor(dat5$Presence_absence, levels = c("Yes", "No", "nvt"))

######## Tree + CRISPR presence absence, + number of spacers = FIGURE 5B
p<- ggtree(tree_root, color = "grey30", linewidth = 0.3) + geom_tiplab(align = TRUE, linesize = 0.2, size = 0, offset = 0.0015) +
          geom_highlight(node = 140, fill = "#E69F00", alpha = 0.25,  extendto = 0.04, to.bottom = TRUE) +
          geom_highlight(node = 175, fill = "#56B4E9", alpha = 0.25,  extendto = 0.04, to.bottom = TRUE) +
          geom_highlight(node = 148, fill = "#56B4E9", alpha = 0.25,  extendto = 0.04, to.bottom = TRUE) +
          geom_highlight(node = 145, fill = "#56B4E9", alpha = 0.25, extendto = 0.04, to.bottom = TRUE) 
       
p1<- p + geom_fruit(data = dat5, geom=geom_tile, mapping=aes(y = id.y, x = Gene, fill = Gene, alpha = Presence_absence), 
                    width = 0.007, colour= "white", offset = 0.07, size = 0.1, pwidth = 0.4, 
                    axis.params = list()) + ylim(0,98) + theme(legend.position = "None")  +
                    scale_alpha_manual(values = c(1,0.2,0)) + 
                    scale_fill_manual(values = c("#79bfc7", "#fdc13f", "#b5fa55"))
                                                                                                
p2 <- p1 + geom_fruit(data = dat4,
                  geom = geom_bar,
                  mapping = aes(y=id.y, x=CRISPRspacers),
                  fill = "grey30",
                  width = 0.5,
                  pwidth = 1, 
                  offset = 0.20,
                  stat = 'identity',
                  axis.params = list(axis = "xy",
                                     line.color = "black",
                                     text.size = 2, vjust = 1.5)) +
                  ylim(0,98)   # a bit hacky, to make space for the axis labels (bar graph)

ggsave('test_plus_outgroupje_V2.pdf', p2, dpi = 300, width = 6, height = 11, units = "cm")

to_remove <- c("barcode38", "barcode51", "Xfragariae")
CRISPRx %>% filter(!(barcode %in% to_remove)) %>% filter(CRISPRspacers > 0) %>% summarise(mean = mean(CRISPRspacers), sum = sum(CRISPRspacers))

############ Data from CRISPRTarget run:
### total number of evidence level 4 spacers: 3740
### number of nonredundant spacers: 877

### Read CRISPRTarget results
targets <- data.frame(read_tsv("pan_xcr_v2_targets.txt", col_names = TRUE))

# remove hits on duplicate genome that needs to be removed.
targets <- targets[targets$Protospacer_description != "barcode38_contig_1", ]

## add class data
targets <- targets %>% mutate(class = ifelse(grepl("phage|virus", targets$Protospacer_description, ignore.case = TRUE), "Phage",ifelse(grepl("plasmid", targets$Protospacer_description, ignore.case = TRUE), "Plasmid", ifelse(grepl("barcode", targets$Protospacer_description, ignore.case = TRUE), "Xcc genome", "Unclassified"))))

targets <- targets %>% mutate(class = case_when(str_detect(class, "Unclassified") ~ "Phage", TRUE ~ class)) 

targets_HQ <- targets %>% filter(Score > 28)

count_per_class <- targets %>% filter(Score > 28) %>% count(Spacer_index, class) %>% count(class)
sum(count_per_class$n) # total number of spacers with known targets
count_per_class <- rbind(count_per_class, c("Unknown", 877 - sum(count_per_class$n)))  # subtract that sum from the total number of nonredundant spacers to find out how many are still unknown
count_per_class$n <- as.integer(count_per_class$n)
count_per_class$threshold = "> 28"

cc <- count_per_class
count_per_class <- rbind(count_per_class, cc)

################### target distribution stacked bar graph FIGURE 5D
colours = c("Phage" = "#CC79A7", "Xcc genome" = "#E69F00", "Plasmid" = "#D0E442", "Unknown" = "grey80")

p0 <- ggplot(data = count_per_class, aes(fill = reorder(class, n), y = n, x = threshold)) + 
      geom_bar(position = "stack", stat = "identity", colour = "white", lwd = 0.1) + theme_mp + scale_y_continuous(breaks = breaks_pretty(n = 6), expand = c(0.005,0.005)) + 
      # scale_x_discrete(labels = NULL, breaks = NULL) +
      scale_fill_manual(values = colours) +
      labs(y = "Number of spacers", x = "Score threshold", fill = "Spacer target") +
      theme(legend.position = "right") +
      theme(legend.key.width = unit(0.35, "cm"))

ggsave('target_distribution_score28.pdf', p0, dpi = 300, width = 6, height = 5, units = "cm")


## calculate number of targets per spacer, and number of spacers per target
targets_per_spacer <- targets %>% filter(Score > 28) %>% group_by(Spacer_index, class) %>% count(Spacer_index, name = "n_targets") %>% group_by(Spacer_index) %>% top_n(1, n_targets)
spacers_per_target <- targets %>% filter(Score > 28) %>% group_by(Protospacer_seq_id, class) %>% count(Protospacer_seq_id, name = "n_spacers") %>%  group_by(Protospacer_seq_id) %>% top_n(1, n_spacers)

################# plot targets per spacer FIGURE 5E
plot1 <- ggplot(data = targets_per_spacer) + geom_boxplot(aes(x = class, y = n_targets), outlier.shape=NA) + 
    geom_jitter(aes(x = class, y = n_targets, colour = class), shape = 16, width = 0.2, alpha = 0.85, size = 0.8) + theme_mp + 
    theme(axis.title.x=element_blank()) +
    #axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
    scale_y_continuous(breaks = breaks_pretty(n = 6)) +
    scale_colour_manual(values = colours) +
    theme(legend.position = "none") +
    ylab("Number of distinct targets\nper spacer") 

ggsave('targets_per_spacer2.pdf', plot1, dpi = 300, width = 6, height = 5, units = "cm")

################# plot spacers per target FIGURE 5F
plot2 <- ggplot(data = spacers_per_target) + geom_boxplot(aes(x = class, y = n_spacers), outlier.shape=NA) + 
  geom_jitter(aes(x = class, y = n_spacers, colour = class), shape = 16, width = 0.2, alpha = 0.85, size = 0.8) + theme_mp + 
  theme(axis.title.x=element_blank()) +
 #, axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
  scale_x_discrete(guide_axis(n.dodge = 2)) +
  labs(y = "Number of distinct spacers\nmapping to each target") +
  scale_colour_manual(values = colours) +
  theme(legend.position = "none") 
  
ggsave('spacers_per_target2.pdf', plot2, dpi = 300, width = 6, height = 5, units = "cm")

# long dplyr clause to find the 'main class' targetted by each spacer. I.e., it removes double entries where a spacer would match both on Plasmids, and on Xcc genome database (which includes those plasmids.)
targets_per_spacer <- targets %>% filter(Score > 28) %>% group_by(Spacer_index, class) %>% count(Spacer_index, name = "n_targets") %>% group_by(Spacer_index) %>% top_n(1, n_targets)

############### # plot prevalence or conservation of spacers
conservation <- data.frame(read_tsv("presence_list.tsv", col_names = TRUE))

# remove duplicate samples
duplicates <- c("barcode38", "barcode51")
conservation <- conservation %>% filter(!(Barcode_present %in% duplicates))

conservation_per_spacer <- conservation %>% count(spacer_id, name = "prevalence")
merged_data <- merge(targets_per_spacer, conservation_per_spacer, by.x="Spacer_index", by.y="spacer_id", all = TRUE) %>% mutate_all(~replace(., is.na(.), 0))

### FIGURE 5C
plot3 <- ggplot(data = merged_data) + geom_violin(aes(x = as.factor(1), y = prevalence)) + geom_jitter(aes(x = as.factor(1), y = prevalence), size = 0.8, shape = 16, width = 0.2, alpha = 0.2) + 
  theme_mp + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0,40)) +
  ylab("Prevalence of each spacer \n(out of 40 CRISPR+ genomes)")
  
ggsave('spacers_prevalence2.pdf', plot3, dpi = 300, width = 4, height = 5, units = "cm")

## Stitching the subpanels together it together
library(patchwork)
collage <- plot3 | p0 | plot1 | plot2
ggsave("CRISPR_count_figures_v2.pdf", width = 20, height = 6, units = "cm")


