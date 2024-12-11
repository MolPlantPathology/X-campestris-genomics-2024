library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(patchwork)
library(stringr)
library(readr)

setwd("/Users/misha/surfdrive/06 Writing/Papers/Ch4_Xcc Genomics/figure_1D_fastANI")

# load the metadata 
dat <- data.frame(read.csv("metadata.csv", header = TRUE))
duplicates <- c("barcode38","barcode51")

# load the ANI Data
ANI <- data.frame(read_tsv("all_v_all.txt"))

ANI_cleaned <- ANI %>%  mutate(q = str_replace(q, "\\.fasta", "")) %>% mutate(r = str_replace(r, "\\.fasta", ""))

ANI_metadata <- merge(ANI_cleaned, dat, by.x = "r", by.y = "barcode") %>% mutate(r_pathovar = pathovar) %>% mutate(pathovar = NULL) %>%
                merge(dat, by.x = "q", by.y = "barcode") %>% mutate(q_pathovar = pathovar) %>% mutate(pathovar = NULL) %>% mutate(combination = paste(q_pathovar, r_pathovar, sep = "_")) %>% 
                filter(!grepl("np", combination))  %>% #remove combinations with non pathogenic isolate
                filter(combination != "Xcr_Xcc")  %>% #remove Xcr_Xcc combinations (as they are the same as Xcc_Xcr combinations)
                filter(ANI != 100) %>% #remove ANI 100 (== self combinations)
                filter(!(q %in% duplicates)) %>%
                filter(!(r %in% duplicates))
              

cols <- c("Xcc_Xcc" = "#E69F00", "Xcr_Xcr" = "#56B4E9", "Xcc_Xcr" = "grey60")
my_order <- c("Xcc_Xcc", "Xcr_Xcr", "Xcc_Xcr")

ANI_metadata$combination <- factor(ANI_metadata$combination, levels = my_order)


density_plots <- ggplot(dat = ANI_metadata, aes(x = ANI)) + geom_density(aes(colour = combination, fill = combination), alpha = 0.5, linewidth = 0.2) + facet_wrap(.~combination, scales = "free_y", ncol = 1) +
                        theme_bw() + 
                        scale_colour_manual(values = cols) +
                        scale_fill_manual(values = cols) +
                        scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
                        theme(panel.grid.major.y = element_blank()) +
                        theme(panel.grid.minor = element_blank()) +
                        theme(strip.text = element_text(size = 6)) +
                        theme(strip.background = element_rect(fill = "grey90")) +
                        theme(axis.text=element_text(colour="black")) +
                        theme(legend.key.size = unit(0.3, "cm")) +
                        theme(legend.position = "None") +
                        theme(text = element_text(size = 8, family="sans", colour = "black")) +
                        theme(
                          strip.text = element_blank(),       # Hide the text labels
                          strip.background = element_blank()  # Hide the background
                        ) +
                      scale_y_continuous(labels = NULL) +
                      labs(x = "Average Nucleotide Identity", y ="Density") +
                      theme(axis.ticks.y = element_blank())


ggsave("density_plot_ANI.pdf", density_plots, dpi = 300, width = 8, height = 6, units = "cm")

ANI_metadata %>% filter(combination == "Xcc_Xcc") %>% summarise(m = mean(ANI))
ANI_metadata %>% filter(combination == "Xcc_Xcr") %>% summarise(m = mean(ANI))
ANI_metadata %>% filter(combination == "Xcr_Xcr") %>% summarise(m = mean(ANI))












