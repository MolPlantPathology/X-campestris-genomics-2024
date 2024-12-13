library(here)
library(ggplot2)
library(dplyr)
library(forcats)

# read data and add an informative header
dat <- read.delim("closest_result.tsv", header = FALSE)
colnames(dat) <- c("Barcode", "Dataset", "Iteration", "T3E", "Distance")


distribution <- ggplot(data = dat) + 
                      geom_histogram(aes(x = Distance), bins = 100) +
                      facet_wrap(~Barcode~Dataset, scales = "free_y") + theme_bw() +
                      xlab("Distance between T3E and insertion sequence (bp)") +
                      theme(text = element_text(size = 8, family="sans")) + 
                      theme(axis.text = element_text(colour = "black", size = 8))

ggsave("distribution.pdf", distribution, width = 10, height = 10, units = "cm")


 # Okay, so the distributions look pretty skewed, both for Xcc and Xcr, and for both randomized and observed datasets

# Lets summarize the dataset with their mean and median
summary_dat <- dat %>% group_by(Barcode, Dataset) %>%
      summarise(mean = mean(Distance), median = median(Distance), n = n()) %>%
      mutate(T3E_class = case_when(Dataset == "Randomized" ~ "Randomized",
                                   TRUE ~ "All"))

# Here we add a new column to check whether an T3E is accessory or other, in Xcc, based on fig 2 or 3 of the paper
dat_xcc_observed <- dat %>% 
  filter(Barcode == "barcode71") %>% 
  filter(Dataset == "Observed") %>% 
  mutate(T3E_class = case_when(
    T3E %in% c("XopH1 ", "XopJ5 ","AvrBs1 ", "XopAH ", "XopG1 ","XopAE ","XopAL2 ","XopE2 ", "XopAC ", "AvrXccA2 ") ~ "Xcc accessory",
    TRUE ~ "Other"
  )
)

summary_dat_xcc_observed <- dat_xcc_observed %>% group_by(Barcode, Dataset, T3E_class) %>%
  summarise(mean = mean(Distance), median = median(Distance), n = n())

## add data frames together and relevel the factor T3E class, so the order looks nice in the graph
summary_dat <- bind_rows(summary_dat, summary_dat_xcc_observed) %>%
    mutate(T3E_class = fct_relevel(T3E_class, "Randomized", "All", "Xcc accessory", "Other"))
# gives us a warning but seems to work.

# plot
cols <- c("Xcc" = "#E69F00", "Xcr" = "#56B4E9", "np" = "grey80")

xcc_plot <- ggplot(data = summary_dat %>% filter(Barcode == "barcode71" & T3E_class != "Other")) + 
                geom_col(aes(x = T3E_class, y = median/1000, fill = T3E_class)) + 
                geom_text(aes(x = T3E_class, y = (median/1000) + 3.5, label = paste(median,"bp")), size = 3) +
                scale_fill_manual(values = c("Randomized" = "grey80", "All" = "#E69F00", "Xcc accessory" = "#E69F00", "Other" = "#E69F00")) +
                theme_bw() +
                theme(text = element_text(size = 8, family="sans")) + 
                theme(axis.text = element_text(colour = "black", size = 8)) + 
                theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) + 
                theme(axis.ticks = element_line(colour = "black")) + 
                theme(panel.grid.minor = element_blank()) +
                theme(panel.grid.major.x = element_blank()) +
                scale_y_continuous(limits = c(0, 42)) +
                theme(legend.position = "None") +
                xlab("T3E dataset") +
                ylab("Median distance to nearest IS (kbp)") 

xcc_plot
ggsave("xcc_plot_h4.pdf", xcc_plot, width = 6, height = 4, units = "cm")


xcr_plot <- ggplot(data = summary_dat %>% filter(Barcode == "barcode01")) + 
  geom_col(aes(x = T3E_class, y = median, fill = T3E_class)) + 
  scale_fill_manual(values = c("Randomized" = "grey80", "All" = "#56B4E9")) +
  theme_bw() +
  theme(text = element_text(size = 8, family="sans")) + 
  theme(axis.text = element_text(colour = "black", size = 8)) + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) + 
  theme(axis.ticks = element_line(colour = "black")) + 
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  scale_y_continuous(limits = c(0, 50000)) +
  theme(legend.position = "None") +
  xlab("T3E dataset") +
  ylab("Median distance to nearest IS") 

ggsave("xcr_plot.pdf", xcr_plot, width = 5, height = 6, units = "cm")

### KS tests

# Test 1: all observed Xcc T3Es vs Xcc randomized
x <- dat %>% filter(Barcode == "barcode71" & Dataset == "Observed") %>% 
  select(Distance)

y <- dat %>% filter(Barcode == "barcode71") %>% filter(Dataset == "Randomized") %>% select(Distance)

ks.test(x$Distance, y$Distance) # D = 0.32664, p-value = 0.001521

# Test 2: Xcc accessories vs Xcc randomized
x <- dat_xcc_observed %>% filter(T3E_class == "Xcc accessory") %>% 
  select(Distance)

ks.test(x$Distance, y$Distance) # D = 0.9081, p-value = 0.00199

# Test 3: Xcc other vs Xcc randomized
x <- dat_xcc_observed %>% filter(T3E_class == "Other") %>% 
  select(Distance)

ks.test(x$Distance, y$Distance) # D = 0.22626, p-value = 0.1165


# Test 4: Xcr T3Es vs Xcr randomized
# Test 1: all observed Xcc T3Es vs Xcc randomized
x <- dat %>% filter(Barcode == "barcode01" & Dataset == "Observed") %>% 
  select(Distance)

y <- dat %>% filter(Barcode == "barcode01") %>% filter(Dataset == "Randomized") %>% select(Distance)

ks.test(x$Distance, y$Distance) # D = 0.43656, p-value = 0.00477
