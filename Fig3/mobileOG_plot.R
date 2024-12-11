library(readr)
library(tidyr)
library(scales)
library(dplyr)

setwd("/Users/misha/surfdrive/03 Xanthomonas genomics/Genomics/2023_09_08_mobileOG")

data <- read.csv("mobileOG_total.csv")
metadata <- read.delim('metadata.tsv')

data <- data %>% mutate(barcode = substr(Query.Title, 1, 9)) # add barcode column

data_merge <-  merge(data, metadata, by = "barcode") # add pathovar column

# dplyr magic
data_dedup <- data_merge %>% select(barcode, pathovar, Major.mobileOG.Category, Contig.ORF.Name) %>% group_by(barcode, pathovar, Contig.ORF.Name, Major.mobileOG.Category) %>% summarize() %>% ungroup() %>% group_by(barcode, pathovar, Major.mobileOG.Category) %>% summarize(n = n())
          
colours = c("Xcc" = "#E69F00", "Xcr" = "#56B4E9")

data_dedup$Major.mobileOG.Category <- gsub("/", ", ", data_dedup$Major.mobileOG.Category)

# remove duplicate samples
duplicates <- c("barcode38", "barcode51")
data_dedup <- data_dedup %>% filter(!(barcode %in% duplicates))

       
plotje <- ggplot(data = data_dedup, aes(x = pathovar, y = n)) +
              geom_boxplot(aes(fill = factor(pathovar)), colour = "black", alpha = 0.7, show.legend = FALSE, outlier.shape = NA, lwd = 0.3, width = 0.7) +
              geom_jitter(aes(colour = factor(pathovar)), alpha = 0.5, shape = 19, size = 1.2, stroke = 0, width = 0.1, height = 0) +
              facet_wrap(~Major.mobileOG.Category, scales = "free_y", ncol = 5, labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
              scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
              expand_limits(y = 0) +
              theme_bw() +
              theme(text = element_text(size = 10, family="sans")) +
              theme(strip.text = element_text(size = 8, family = 'sans')) + 
              theme(axis.text = element_text(colour = "black")) + 
              theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) + 
              theme(axis.ticks = element_line(colour = "black")) + 
              theme(panel.grid.minor = element_blank()) +
              theme(panel.grid.major.x = element_blank()) +
              theme(legend.position = "None") +
              xlab("Pathovar") +
              ylab("ORFs in category\n(mobileOG-db)") +
              scale_fill_manual(values=colours) +
              scale_colour_manual(values=colours)
              

ggsave("major_category_v2.pdf", plotje, width = 19, height = 6, units = "cm")


### significance testing
wilcox.test(n ~ pathovar, data = data_dedup[which(data_dedup$pathovar != "np"), ], subset = Major.mobileOG.Category == 'integration, excision')  #  p-value < 2.2e-16
wilcox.test(n ~ pathovar, data = data_dedup[which(data_dedup$pathovar != "np"), ], subset = Major.mobileOG.Category == 'phage') # p-value = 2.926e-09
wilcox.test(n ~ pathovar, data = data_dedup[which(data_dedup$pathovar != "np"), ], subset = Major.mobileOG.Category == 'replication, recombination, repair') # p-value = 1.109e-10
wilcox.test(n ~ pathovar, data = data_dedup[which(data_dedup$pathovar != "np"), ], subset = Major.mobileOG.Category == 'stability, transfer, defense') # p-value = 0.04302
wilcox.test(n ~ pathovar, data = data_dedup[which(data_dedup$pathovar != "np"), ], subset = Major.mobileOG.Category == 'transfer') # p-value = 1.821e-13


