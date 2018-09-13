# plotting spectra results

library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)

# change your working directory computer
# setwd("C:\\Users\\Scott\\Google Drive\\Projects\\phaeo-miao\\")

pep1 <- read.csv("data/ms_data/AKPNFYVK_SSU2_sample5A.csv")
pep1$sequence <- rep('AKPNFYVK', nrow(pep1))

pep2 <- read.csv("data/ms_data/AWIAQIK_Fld1_sample5A.csv")
pep2$sequence <- rep('AWIAQIK', nrow(pep2))

pep3 <- read.csv("data/ms_data/GDSITWINNK_PC1_sample5A.csv")
pep3$sequence <- rep('GDSITWINNK', nrow(pep3))

pep4 <- read.csv("data/ms_data/GGPHNVVFVEDAIPK_PC2_sample5A.csv")
pep4$sequence <- rep('GGPHNVVFVEDAIPK', nrow(pep4))

pep5 <- read.csv("data/ms_data/QIQYALNK_SSU2_sample5A.csv")
pep5$sequence <- rep('QIQYALNK', nrow(pep5))

pep1b <- melt(pep1, id.vars = c('time_min', 'sequence'))
pep2b <- melt(pep2, id.vars = c('time_min', 'sequence'))
pep3b <- melt(pep3, id.vars = c('time_min', 'sequence'))
pep4b <- melt(pep4, id.vars = c('time_min', 'sequence'))
pep5b <- melt(pep5, id.vars = c('time_min', 'sequence'))

finale <- rbind(pep1b, pep2b, pep3b, pep4b, pep5b)

names(finale)[3] <- c('Fragment Ion')

finale %>%
  filter(time_min > 8 & time_min < 20) %>%
  group_by(sequence) %>% 
  mutate(int_scaled = rescale(value, to = c(0, 100))) %>% 
  ggplot(aes(x = time_min, y = int_scaled, colour = `Fragment Ion`)) + 
  geom_line(lwd = 0.2) +
  # geom_point() +
  facet_wrap(~sequence, nrow = 3) + 
  theme_bw() + 
  ylab('Scaled Ion Intensity') + 
  xlab('Retention Time')
