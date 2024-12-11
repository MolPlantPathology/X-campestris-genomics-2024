library(readr)
library(readxl)
library(here)
library(dplyr)
library(ggplot2)


setwd("/Users/misha/surfdrive/05 Lab/09_CRISPR_experiment/CRISPR_exp3")


data <- read_xlsx(here("Fig6/conjugation3_cfu_counts_figure.xlsx"))

data_filtered <-  data %>% filter(tech_replicate == "mean_of_tech_replicates") %>% filter(plasmid != "no_plasmid")
data_mean <- data_filtered %>% group_by(strain, plasmid_name) %>% summarise(n = n(), mean_transformation_efficiency = mean(transformation_efficiency), sd = sd(transformation_efficiency), se = sd/sqrt(n))
  
colours = c("Xcc" = "#E69F00", "Xcr" = "#56B4E9")

plot <- ggplot() + 
            geom_errorbar(data = data_mean, aes(x = plasmid_name, ymin = mean_transformation_efficiency - se, ymax = mean_transformation_efficiency + se), width = 0.2, linewidth = 0.2) +  # Error bars
            geom_point(data = data_mean, aes(x = plasmid_name, y = mean_transformation_efficiency, fill = strain), size = 2, shape = 22)  + # mean
            geom_jitter(data = data_filtered, aes(x = plasmid_name, y = transformation_efficiency, fill = strain), alpha = 0.6, size = 1.5, width = 0.1, shape = 21) + # individual observations
            scale_fill_manual(values = colours) +
            facet_grid(.~strain) + theme_bw() +
            scale_y_continuous(trans = "log10") +
            annotation_logticks(base = 10, sides = "l", linewidth = 0.3, short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
            theme(panel.grid.minor = element_blank(), panel.grid.major.x=element_blank()) +
            theme(text = element_text(size = 9, family="sans", colour = "black")) +
            theme(axis.text=element_text(colour="black")) +
            theme(legend.position = "None") +
            labs(x = "pBBR1 plasmid", y = "Transformation efficiency\n(tranformants/recipients)") 
          
ggsave("conjugation_3.pdf", plot, width = 7, height = 7, units = "cm", dpi = 1500)

xcr_target <- data_filtered %>% filter(strain == "Xcr") %>% filter(plasmid_name == "target+") %>% select(transformation_efficiency)
xcr_notarget <- data_filtered %>% filter(strain == "Xcr") %>% filter(plasmid_name == "target-") %>% select(transformation_efficiency)

t.test(xcr_target$transformation_efficiency, xcr_notarget$transformation_efficiency, var.equal = TRUE) # p = 0.006 (same as t-test in excel)
t.test(xcr_target$transformation_efficiency, xcr_notarget$transformation_efficiency, var.equal = FALSE) # p = 0.02 (maybe this one is officially better. However, still significantly different)
