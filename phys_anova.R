# anova of phaeo physiology

# setwd("C:\\Users\\Scott\\Google Drive\\Projects\\phaeo-miao\\")

# to do

# reorder factor levels for light from low to high
# to tukeys and anova for light levels

library(tidyverse)
library(edgeR)
library(xlsx)
library(readxl)
library(ggfortify)
library(testthat)
library(broom)
library(multcompView)
library(gridExtra)
library(cowplot)

# reading data in and transforming column name and factor levels
phys_data <- read.csv('data/original physiological data-tidy.csv')
levels(phys_data$treatment) <- c("High Fe High Mn", "High Fe Low Mn", "Low Fe High Mn", "Low Fe Low Mn")
phys_data <- phys_data %>% 
  mutate(light = fct_relevel(phys_data$light, "low", "opt", "high"))
levels(phys_data$light) <- c("Low", "Optimum", "High")

phys_data_low <- phys_data %>% dplyr::filter(light == 'Low')
phys_data_other <- phys_data %>% dplyr::filter(treatment == 'High Fe High Mn')

phys_data$treat_sum <- paste(phys_data$treatment %>% as.character(), ', ', phys_data$light %>% as.character(), ' Light', sep = "")
phys_data$treat_sum <- as.factor(phys_data$treat_sum)
phys_data$treat_sum <- fct_relevel(phys_data$treat_sum, "High Fe High Mn, High Light", 
                                   "High Fe High Mn, Optimum Light", "High Fe High Mn, Low Light")

# anova and tukey for low light and different metal treatments
chla_aov <- aov(chla ~ treatment, data = phys_data_low)
chla_tuk <- TukeyHSD(chla_aov) %>% tidy()
  
fvfm_aov <- aov(fvfm ~ treatment, data = phys_data_low)
fvfm_tuk <- TukeyHSD(fvfm_aov) %>% tidy()

gr_aov <- aov(gr ~ treatment, data = phys_data_low)
gr_tuk <- TukeyHSD(gr_aov) %>% tidy()

size_aov <- aov(size ~ treatment, data = phys_data_low)
size_tuk <- TukeyHSD(size_aov) %>% tidy()

###  anova and tukey for light
chla_aov_l <- aov(chla ~ light, data = phys_data_other)
chla_tuk_l <- TukeyHSD(chla_aov_l) %>% tidy()
  
fvfm_aov_l <- aov(fvfm ~ light, data = phys_data_other)
fvfm_tuk_l <- TukeyHSD(fvfm_aov_l) %>% tidy()

gr_aov_l <- aov(gr ~ light, data = phys_data_other)
gr_tuk_l <- TukeyHSD(gr_aov_l) %>% tidy()

size_aov_l <- aov(size ~ light, data = phys_data_other)
size_tuk_l <- TukeyHSD(size_aov_l) %>% tidy()

# writing anova results to excel file:

# Write the TukeyHSD results from low light into an excel file:
openxlsx::write.xlsx(chla_tuk, file="data/intermediate-data/phys_anova_low_light.xlsx",
           sheetName="chla_tuk", append=FALSE)
# Add a second data set in a new worksheet
openxlsx::write.xlsx(fvfm_tuk, file="data/intermediate-data/phys_anova_low_light.xlsx", 
           sheetName="fvfm_tuk", append=TRUE)
# Add a third data set
openxlsx::write.xlsx(gr_tuk, file="data/intermediate-data/phys_anova_low_light.xlsx", 
           sheetName="gr_tuk", append=TRUE)
openxlsx::write.xlsx(size_tuk, file="data/intermediate-data/phys_anova_low_light.xlsx", 
           sheetName="size_tuk", append=TRUE)


# Write the TukeyHSD results from high metal different light into an excel file:
openxlsx::write.xlsx(chla_tuk_l, file="data/intermediate-data/phys_anova_light_treat.xlsx",
           sheetName="chla_tuk_l", append=FALSE)
# Add a second data set in a new worksheet
openxlsx::write.xlsx(fvfm_tuk_l, file="data/intermediate-data/phys_anova_light_treat.xlsx", 
           sheetName="fvfm_tuk_l", append=TRUE)
# Add a third data set
openxlsx::write.xlsx(gr_tuk_l, file="data/intermediate-data/phys_anova_light_treat.xlsx", 
           sheetName="gr_tuk_l", append=TRUE)
openxlsx::write.xlsx(size_tuk_l, file="data/intermediate-data/phys_anova_light_treat.xlsx", 
           sheetName="size_tuk_l", append=TRUE)


## Composite plots:

p10 <- ggplot(phys_data, aes(x = treat_sum, y = chla)) + 
  geom_point(alpha = 0.6, size = 3) + 
  theme_bw()+
  # annotate("text", x = 0.65, y = 700, label = 'A)', size = 5) +
  ylim(0, 800) + 
  
  annotate("text", x = 1, y = 760, label = 'italic(a)', parse = TRUE, colour = 'grey50') + 
  annotate("text", x = 2, y = 760, label = 'italic(b)', parse = TRUE, colour = 'grey50') + 
  annotate("text", x = 3, y = 760, label = 'italic(c)', parse = TRUE, colour = 'grey50') + 
  
  annotate("text", x = 3, y = 695, label = 'a') + 
  annotate("text", x = 4, y = 695, label = 'ab') + 
  annotate("text", x = 5, y = 695, label = 'bc') + 
  annotate("text", x = 6, y = 695, label = 'c') + 
  ylab("Chlorophyll a \n(fg/cell)") + 
  xlab("") + 
  theme(axis.text.x = element_blank());p10
  # theme(axis.text.x = element_text(size = 7, angle = 15, hjust = 1))

p11 <- ggplot(phys_data, aes(x = treat_sum, y = fvfm)) + 
  
  geom_point(alpha = 0.6, size = 3) + 
  theme_bw()+
  
  ylim(0.1, 0.85) +
  
  
  annotate("text", x = 1, y = 0.82, label = 'italic(a)', parse = TRUE, colour = 'grey50') + 
  annotate("text", x = 2, y = 0.82, label = 'italic(a)', parse = TRUE, colour = 'grey50') + 
  annotate("text", x = 3, y = 0.82, label = 'italic(b)', parse = TRUE, colour = 'grey50') + 
  
  # annotate("text", x = 0.65, y = 0.8, label = 'B)', size = 5) +
  annotate("text", x = 3, y = 0.75, label = 'a') +
  annotate("text", x = 4, y = 0.75, label = 'b') +
  annotate("text", x = 5, y = 0.75, label = 'c') +
  annotate("text", x = 6, y = 0.75, label = 'c') +
  
  ylab("Fv/Fm") + 
  xlab("") + 
  theme(axis.text.x = element_blank());p11
  # theme(axis.text.x = element_text(size = 7, angle = 15, hjust = 1))
p12 <- ggplot(phys_data, aes(x = treat_sum, y = gr)) + 
  
  geom_point(alpha = 0.6, size = 3) + 
  theme_bw()+
  ylim(0.15, 0.5) +
  
  # annotate("text", x = 1, y = 0.45, label = 'bold(a)', parse = TRUE) + 
  # annotate("text", x = 2, y = 0.45, label = 'bold(a)', parse = TRUE) + 
  # annotate("text", x = 3, y = 0.45, label = 'bold(b)', parse = TRUE) + 
  
  # annotate("text", x = 0.65, y = 0.48, label = 'C)', size = 5) +
  annotate("text", x = 3, y = 0.45, label = 'ab') +
  annotate("text", x = 4, y = 0.45, label = 'a') +
  annotate("text", x = 5, y = 0.45, label = 'b') +
  annotate("text", x = 6, y = 0.45, label = 'ab') +
  
  ylab(expression(paste("Growth Rate (d"^-1, ")"))) +
  xlab("") + 
  theme(axis.text.x = element_blank());p12
  # theme(axis.text.x = element_text(size = 7, angle = 15, hjust = 1))
p13 <- ggplot(phys_data, aes(x = treat_sum, y = size)) +
  
  geom_point(alpha = 0.6, size = 3) + 
  theme_bw()+
  ylim(520000, 885000) +
  
  # annotate("text", x = 0.65, y = 860000, label = 'D)', size = 5) +
  
  annotate("text", x = 1, y = 870000, label = 'italic(a)', parse = TRUE, colour = 'grey50') + 
  annotate("text", x = 2, y = 870000, label = 'italic(b)', parse = TRUE, colour = 'grey50') + 
  annotate("text", x = 3, y = 870000, label = 'italic(c)', parse = TRUE, colour = 'grey50') + 
  
  annotate("text", x = 3, y = 835000, label = 'a') +
  annotate("text", x = 4, y = 835000, label = 'a') +
  annotate("text", x = 5, y = 835000, label = 'ab') +
  annotate("text", x = 6, y = 835000, label = 'b') +
  ylab("Forward Scatter \nValues") + 
  xlab("") + 
  # theme(axis.text.x = element_text(size = 7, angle = 15, hjust = 1))
  theme(axis.text.x = element_blank());p13

dev.off()

tiff("figures/main-all-phys-data.tiff", width=23*0.75, height=27.94*0.65, units = 'cm', res=400)

plot_grid(p10, p11, p12, p13, align = 'v', nrow = 4)

dev.off()




