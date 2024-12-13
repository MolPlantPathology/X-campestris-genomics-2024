

theme_mp <- theme_bw() + theme(axis.text = element_text(color = "black", size = 8)) +
  theme(axis.title = element_text(color = "black", size = 8)) +
  theme(panel.grid.major = element_blank() , panel.grid.minor = element_blank()) +
  theme(panel.grid.major.y = element_line(linewidth = 0.3)) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) + 
  theme(panel.background = element_rect(colour = "black")) +
  theme(axis.ticks = element_line(colour = "black")) + 
  theme(strip.background=element_rect(fill="grey95"), strip.text = element_text(size = 8, face = 'bold')) 

setwd("/Users/misha/surfdrive/03 Xanthomonas genomics/Genomics/2023_09_08_SVs/2023_11_29_IS_inv_bp_overlap_newdataset")

data <- read.csv("IS_inversions_overlap.csv")

data_longer <- read.csv("IS_inversions_overlap_longDF.csv")

colours = c("Xcc" = "#E69F00", "Xcr" = "#56B4E9")

plot <- ggplot(data = data_longer %>% filter(pathovar != "np"), aes(x = treatment, y = percent_overlap, fill = pathovar)) + 
              geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
              geom_point(aes(group = pathovar), position=position_dodge(width=0.75), shape = 21, size = 1, alpha = 0.5, stroke = 0) + theme_mp + 
              labs(y="Insertion sequences overlapping\nwith inversion breakpoint (%)", x = "Insertion sequence dataset") + 
              scale_fill_manual(values = colours) + 
              theme(legend.position = "none")


data_longer %>% group_by(treatment, pathovar) %>% summarise(m = mean(percent_overlap))

ggsave("enrichment_IS_inversion_BP.pdf", plot, width = 6.1, height = 6, units = "cm", dpi = 300)
