library(ggplot2)
library(readr)
library(ggrepel)
library(patchwork)
# library(ggridges)

theme_mp <- theme_bw() + theme(axis.text = element_text(color = "black", size = 8)) +
            theme(axis.title = element_text(color = "black", size = 8)) +
            theme(panel.grid.major = element_line(linewidth = 0.5), panel.grid.minor = element_blank()) +
            theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) + 
            theme(panel.background = element_rect(colour = "black")) +
            theme(axis.ticks = element_line(colour = "black")) + 
            theme(strip.background=element_rect(fill="grey95"), strip.text = element_text(size = 8, face = 'bold')) 


GC_content_bc71 <- read_csv("C:/Users/PhD/surfdrive/03 Xanthomonas genomics/Genomics/2022_02_15_GC_content_per_gene/GC_content_bc71.csv")


setwd("~/surfdrive/03 Xanthomonas genomics/Genomics/2022_02_15_GC_content_per_gene/")
GC_content_bc71 <- read_csv("~/surfdrive/03 Xanthomonas genomics/Genomics/2022_02_15_GC_content_per_gene/GC_content_bc71.csv")
GC_content_bc01 <- read_csv("~/surfdrive/03 Xanthomonas genomics/Genomics/2022_02_15_GC_content_per_gene/GC_content_bc01.csv")

GC_content_bc71$pathovar = "Xcc"
GC_content_bc01$pathovar = "Xcr"

dat <- rbind(GC_content_bc71, GC_content_bc01)

## to make the categories!

xcr <- paste0("bc01_0", seq(1343, 1358))
xcc <- paste0("bc71_0", seq(3040,3055))


dat <- dat %>% mutate(category = case_when(
  str_detect(Description, "Xop") ~ "Type III effector",
  str_detect(Description, "Rip") ~ "Type III effector",
  str_detect(Description, "Avr") ~ "Type III effector",
  Name %in% xcr ~ "Type III secretion system",
  Name %in% xcc ~ "Type III secretion system",
  TRUE ~ "no category"  # If none of the conditions are met
))

## remove tRNAs based on gene sisze 
dat <- dat %>% filter(Length > 96)


dat_xcr <- dat %>% filter(pathovar == "Xcr")
dat_xcc <- dat %>% filter(pathovar == "Xcc")

xcc_qCAI <- quantile(dat_xcc$CAI, probs = seq(0,1,1/10), na.rm = TRUE, names = FALSE)
xcc_qGC <- quantile(dat_xcc$GC, probs = seq(0,1,1/10), na.rm = TRUE, names = FALSE)

xcr_qCAI <- quantile(dat_xcr$CAI, probs = seq(0,1,1/10), na.rm = TRUE, names = FALSE)
xcr_qGC <- quantile(dat_xcr$GC, probs = seq(0,1,1/10), na.rm = TRUE, names = FALSE)


pathovar_to_draw = "Xcc"

# plot simple
pl1 <- ggplot() +  
            geom_rect(aes(xmin = 0, xmax = xcc_qCAI[2], ymin = 0, ymax =  xcc_qGC[2]), fill = "grey85",  alpha = 0.5) +
            geom_point(data = dat %>% filter(pathovar == pathovar_to_draw & category == "no category"), aes(x = CAI, y = GC), alpha = 0.1, size = 1.5, colour = "black") + 
            geom_point(data = dat %>% filter(pathovar == pathovar_to_draw & category == "Type III secretion system"), aes(x = CAI, y = GC), shape = 22, colour = "white", stroke = 0.5, alpha = 1, size = 2, fill = "#228833") + 
            geom_point(data = dat %>% filter(pathovar == pathovar_to_draw & category == "Type III effector"), aes(x = CAI, y = GC), shape = 21, colour = "white", stroke = 0.5, alpha = 1, size = 2, fill = "#EE6677") + 
            geom_label_repel(data = dat %>% filter(pathovar == pathovar_to_draw & category == 'Type III effector' & CAI < xcc_qCAI[2] & GC < xcc_qGC[2]), aes(x = CAI, y = GC, label = Description), label.padding = 0.15, max.overlaps=Inf, min.segment.length = 0, size = 2.8, segment.colour = '#EE6677') +
            theme_mp + 
            facet_wrap(.~pathovar) + 
            theme(strip.background=element_rect(fill= scales::alpha("#E69F00", 0.2))) + 
            coord_cartesian(ylim = c(35,75), xlim = c(0.2, 0.9)) + 
            xlab("Codon Adaptation Index") +
            ylab("GC content") 
           

pathovar_to_draw = "Xcr"

pl2 <- ggplot() +  
  geom_rect(aes(xmin = 0, xmax = xcr_qCAI[2], ymin = 0, ymax =  xcr_qGC[2]), fill = "grey85",  alpha = 0.5) +
  geom_point(data = dat %>% filter(pathovar == pathovar_to_draw & category == "no category"), aes(x = CAI, y = GC), alpha = 0.1, size = 1.5, colour = "black") + 
  geom_point(data = dat %>% filter(pathovar == pathovar_to_draw & category == "Type III secretion system"), aes(x = CAI, y = GC), shape = 22, colour = "white", stroke = 0.5, alpha = 1, size = 2, fill = "#228833") + 
  geom_point(data = dat %>% filter(pathovar == pathovar_to_draw & category == "Type III effector"), aes(x = CAI, y = GC), shape = 21, colour = "white", stroke = 0.5, alpha = 1, size = 2, fill = "#EE6677") + 
  geom_label_repel(data = dat %>% filter(pathovar == pathovar_to_draw & category == 'Type III effector' & CAI < xcr_qCAI[2] & GC < xcr_qGC[2]), aes(x = CAI, y = GC, label = Description), label.padding = 0.15, max.overlaps=Inf, min.segment.length = 0, size = 2.8, segment.colour = '#EE6677') +
  theme_mp + 
  facet_wrap(.~pathovar) + 
  theme(strip.background=element_rect(fill= scales::alpha("#56B4E9", 0.2))) +
  coord_cartesian(ylim = c(35,75), xlim = c(0.2, 0.9)) + 
  xlab("Codon Adaptation Index") +
  ylab("GC content")

  

pl <- (pl1/pl2)

ggsave('GC_CAI_v3.pdf', pl, width = 9, height = 16, units = 'cm')


### no longer used but maybe useful at some point
colours <- c("no category" = "black", "Type III effector" = "#0072B2", "Type III secretion system" = "#CC79A7")
alphas <- c("no category" = 0.05, "Type III effector" = 1, "Type III secretion system" = 1)
sizes <- c("no category" = 1, "Type III effector" = 2, "Type III secretion system" = 2)
order <- c( "Type III effector", "no category", "Type III secretion system")
dat$category <- factor(dat$category, levels = order)