library(readr)
library(tidyr)
library(scales)
library(ggplot2)
library(dplyr)
library(stringr)

# set WD
setwd("/Users/misha/surfdrive/03 Xanthomonas genomics/Genomics/2023_12_06_defensefinder") # change to local working directory

# read the data
data <- read.table("defensefinder_data.tsv", header = TRUE,, sep = "\t")
metadata <- read.delim('metadata.tsv')

# add barcode column and pathovar column
data <- data %>% mutate(barcode = substr(sys_id, 1, 9)) # add barcode column
data <-  merge(data, metadata, by = "barcode") # add pathovar column

# transform some columns to factor
data$pathovar <- as.factor(data$pathovar)
data$sys_id <- as.factor(data$sys_id)
data$type <- as.factor(data$type)
data$subtype <- as.factor(data$subtype)

# remove duplicate samples
duplicates <- c("barcode38", "barcode51")
data <- data %>% filter(!(barcode %in% duplicates))

# count the number of defense systems per barcode
count_per_barcode <- data %>% group_by(barcode) %>% summarize(n_systems = n())
count_per_barcode <- merge(count_per_barcode, metadata, by = "barcode")

colours <- c("Xcc" = "#E69F00", "Xcr" = "#56B4E9", "np" = "grey70")

#### FIGURE 5H
plot1 <- ggplot(data = count_per_barcode, aes(x = pathovar, y = n_systems, fill = pathovar)) + geom_boxplot(alpha = 0.4, outlier.shape = NA) + 
          geom_jitter(aes(colour = pathovar), alpha = 0.5, height = 0, width = 0.2, size = 0.8) + ylim(c(0,18)) + 
          scale_fill_manual(values = colours) + scale_color_manual(values = colours) + 
          theme_bw() +
            theme(axis.text.y = element_text(size = 8, family="sans",colour = "black"))+
            theme(axis.text.x = element_text(size = 8, family="sans",colour = "black")) +
            theme(axis.title.y = element_text(size = 8, family="sans",colour = "black"))+
            theme(axis.title.x = element_text(size = 8, family="sans",colour = "black")) + 
          theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) + 
          theme(axis.ticks = element_line(colour = "black")) + 
          theme(panel.grid.minor = element_blank()) +
          theme(panel.grid.major.x = element_blank()) +
          ylab("Defense systems\n(number per genome)") +
          xlab("Pathovar") +
          theme(legend.position = "None" )

ggsave("defensefinder_systems_per_genome.pdf", plot1, width = 5, height = 5, units = "cm")

# significance testing for figure 5H
wilcox.test(n_systems ~ pathovar, data = count_per_barcode[which(count_per_barcode$pathovar != "np"), ], ) # p-value = 6.287e-06

## count presence absence for all defense systems per pathovar
count_per_system <- data %>% group_by(subtype, barcode, .drop = FALSE) %>% tally() %>% ungroup %>% mutate_if(is.numeric, ~1 * (. > 0)) 
count_per_system <- merge(count_per_system, metadata, by = "barcode") 

count_per_system <- count_per_system %>% group_by(subtype, pathovar) %>% tally(name = "presence") 
count_per_system$pathovar <- as.factor(count_per_system$pathovar)

count_per_system <- count_per_system %>% complete(pathovar, fill = list(presence = 0)) 

count_per_system <- count_per_system %>%
                      group_by(subtype) %>%
                      mutate(total_presence = sum(presence))

## add absence:
count_per_system <- count_per_system %>% filter(pathovar != "np") %>% mutate(absent = case_when(pathovar == "Xcc" ~ 49 - presence, pathovar == "Xcr" ~ 44 - presence))
cont_per_system <- count_per_system  %>% pivot_longer(c(presence, absent), names_to = "status", values_to = "count") %>% mutate(status_2 = paste(status, pathovar, sep = "_")) 

## reorder based to total presence
cont_per_system$subtype <- reorder(cont_per_system$subtype, -cont_per_system$total_presence)

# some color things
colours <- c("Xcc" = "#E69F00", "Xcr" = "#56B4E9", "np" = "grey70")
colours <- c("absent" = "grey90", "presence" = "#56B4E9") 
colours_alpha <- c("absent_Xcc" =  c(alpha("#E69F00", 0.25)), "absent_Xcr" =  c(alpha("#56B4E9", 0.25)), "presence_Xcc" = "#E69F00", "presence_Xcr" = "#56B4E9")

# renaming some to a shorter name for space reasons
cont_per_system <- cont_per_system %>%
  mutate(subtype = case_when(
    subtype == "CAS_Class1-Subtype-I-F" ~ "CRISPRCas IF",
    subtype == "NLR_like_bNACHT01" ~ "bNACHT01",
    subtype == "Lamassu-Cap4_nuclease" ~ "Lamassu-Cap4",
    subtype == "Old_exonuclease" ~ "Exonuclease",
    subtype == "Mokosh_Type_I_A" ~ "Mokosh IA",
    subtype == "Mokosh_TypeII" ~ "Mokosh II",
    # Add more conditions as needed
    TRUE ~ as.character(subtype)  # Keep other entries unchanged
  ))

cont_per_system$subtype <- reorder(cont_per_system$subtype, -cont_per_system$total_presence)

#### Figure S6
plot2 <- ggplot(data = cont_per_system, aes(x = pathovar, y = count)) +
              geom_bar(aes(fill = status_2), stat = "identity", position = "fill", show.legend = FALSE, width = 0.6) +
              facet_wrap(~subtype, ncol = 8) +
              scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0.005,0.005)) +
              theme_bw() +
              theme(strip.text = element_text(size = 8, family = 'sans', margin = margin(0, 0, 0, 0))) +
              theme(strip.background = element_rect(fill = "white", colour = NA)) +
              theme(panel.border = element_rect(fill = NA, colour = "grey50", linewidth = 0.5)) + 
              theme(axis.ticks = element_line(colour = "grey50")) + 
              theme(panel.grid.minor = element_blank()) +
              theme(panel.grid.major = element_blank()) +
              theme(legend.position = "None") +
              xlab("Pathovar") +
              ylab("Proportion of genomes encoding defense system") +
              scale_fill_manual(values=colours_alpha) +
              theme(axis.text.y = element_text(size = 8, family="sans",colour = "black"))+
              theme(axis.text.x = element_text(size = 8, family="sans",colour = "black")) +
              theme(axis.title.y = element_text(size = 8, family="sans",colour = "black"))+
              theme(axis.title.x = element_text(size = 8, family="sans",colour = "black"))
           
ggsave("defensefinder_all_systems.pdf", plot2, width = 17, height = 20, units = "cm")