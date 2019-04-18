# IRS Normalization, plotting, and Differential Expression Analysis

# tmt normalization phaeocystis mn fe experiment

# Note: TMT IRS normalization done as in Plubell, D.L., Wilmarth, P.A., Zhao, Y., Fenton, A.M., Minnier, J., Reddy, A.P., Klimek, J., Yang, X., David, L.L. and Pamir, N., 2017. Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high fat diet specific proteome changes in mouse epididymal adipose tissue. Molecular & Cellular Proteomics, 16(5), pp.873-890.

# Code in this script is modelled after https://github.com/pwilmart/IRS_normalization


############## NOTE

# this script imputes values if one TMT channel is missing. The imputed value is equal to 0.5* lowest value observed

##############

library(tidyverse)
library(edgeR)
library(readxl)
library(gridExtra)
library(ggfortify)
library(testthat)

'%!in%' <- function(x,y)!('%in%'(x,y))

# checks for directory structure

if(!dir.exists('data/')){
  warning('no directory for data')
}

if(!dir.exists("figures/")){
  warning('no directory for figures')
}

if(!dir.exists("data/intermediate-data/")){
  warning('no directory for intermediate data')
}


# Reading in PD data ------------------------------------------------------

## reading in new database searches with repaired gene model for PSI

just_annot2_a <- read_excel(path = 'data/annotation_Pa1374_PsaABchange.asc.xlsx', sheet = 1)
just_annot_appended <- read_excel(path = 'data/annotation_Pa1374_contigs_appended.xlsx')

# there were 61 contigs that were not included in the transcriptome annotation file
# this was because they were low abundant transcripts. We identified 3 of them in the 
# proteome (across all 20 or 19 channels), and manually annotated them with BLASTP.
just_annot2 <- rbind(just_annot2_a, just_annot_appended)

names(just_annot2)[1] <- "protein_id"

# reading in data from proteome discoverer
pd_2_w_crap <- read.csv('data/170529_0692_097_repeat2019_Proteins.txt', sep = '\t')

# converting columns to appropriate ones
pd_2_w_crap$Checked <- as.logical(pd_2_w_crap$Checked)

pd_2_w_crap$Protein.FDR.Confidence.Combined <- as.character(pd_2_w_crap$Protein.FDR.Confidence.Combined)
pd_2_w_crap$Master <- as.character(pd_2_w_crap$Master)
pd_2_w_crap$Accession <- as.character(pd_2_w_crap$Accession)
pd_2_w_crap$Description <- as.character(pd_2_w_crap$Description)

# remove the contaminant matches ('cRAP_')
pd_2_no_impute <- pd_2_w_crap[which(!grepl(pattern = 'cRAP_', x = pd_2_w_crap$Accession)), ]

# imputing function for pd_output (note that it's specific for this file)
impute_values <- function(pd_output_df){
  
  # pd_output_df <- tester
  
  for(i in 1:nrow(pd_output_df)){
    
    # i <- 87
    number_of_nas <- sum(is.na(pd_output_df[i,c(18:37)]))
    
    
    if(number_of_nas == 1){ # if there is only one na in the row, then replace it with half minimum
      na_index <- which(is.na(pd_output_df[i, c(18:37)])) + 17 # determine where that NA is
      half_lowest_value <- 0.5*min(pd_output_df[i, c(18:37)], na.rm = TRUE) # calculate the imputed value
      pd_output_df[i, na_index] <- half_lowest_value # replace it! :) 
    } else if(number_of_nas != 1){
      next # if there is more or less than 1 NA, then move onto the next one
    }
  }
  
  return(pd_output_df)
  
}

# impute channels with one value missing
pd_2 <- impute_values(pd_2_no_impute)

# subset pd1 as S01-S10
pd1_sub_w_cv <- pd_2[c(which(pd_2$`Protein.FDR.Confidence.Combined` == 'High')),
                     c(which(names(pd_2) == 'Accession'),
                       which(grepl(pattern = 'S1.10', x = names(pd_2), fixed = TRUE)))]
pd1_sub <- pd1_sub_w_cv[, which(!grepl(pattern = 'CV.in.Percent', names(pd1_sub_w_cv)))]


# subset pd_2 as S11-S20
pd2_sub_w_cv <- pd_2[c(which(pd_2$`Protein.FDR.Confidence.Combined` == 'High'))
                     ,c(which(names(pd_2) == 'Accession'),
                        which(grepl(pattern = 'S11', x = names(pd_2))))]
pd2_sub <- pd2_sub_w_cv[, which(!grepl(pattern = 'CV.in.Percent', names(pd2_sub_w_cv)))]


# changing the column names in the files
pd1_names <- c('protein_id', 'ol_highmn_highfe_1_1', 'ol_highmn_highfe_2_1', 'ol_highmn_highfe_3_1', 'll_highmn_highfe_1_1', 'll_highmn_highfe_2_1', 'll_highmn_highfe_3_1', 'll_lowmn_highfe_1_1', 'll_lowmn_highfe_2_1', 'hl_highmn_highfe_1_1', 'hl_highmn_highfe_2_1')
pd2_names <- c('protein_id', 'll_lowmn_highfe_3_2', 'll_highmn_lowfe_1_2', 'll_highmn_lowfe_2_2', 'll_highmn_lowfe_3_2', 'll_lowmn_lowfe_1_2', 'll_lowmn_lowfe_2_2', 'll_lowmn_lowfe_3_2', 'hl_highmn_highfe_3_2', 'hl_highmn_highfe_1_2', 'hl_highmn_highfe_2_2')

names(pd1_sub) <- pd1_names
names(pd2_sub) <- pd2_names

# combining data files

pd <- dplyr::full_join(pd2_sub, pd1_sub, by = 'protein_id')

# subsetting out the 'protein_id' column
pd_pa <- pd[,2:ncol(pd)]

# a zero should be interpreted as not observed, hence NA
pd[pd == 0] <- NA

# subsetting only proteins found in all treatments
pd_all <- pd[complete.cases(pd),]

# removing the protein_id column
pd_raw <- pd_all[, 2:ncol(pd_all)]


# Plotting raw total intensities ------------------------------------------

png("figures/s2-all-normalization-boxplots-imputed-values.png", width=23*0.75, height=27.94*0.75, units="cm", res=1200)

par(mfrow = c(2, 2), mar = c(7, 6, 4, 2))

boxplot(log2(pd_raw), 
        col = rep(c('firebrick', 'darkblue'), each = 10), 
        xaxt = 'n', 
        xlab = '', main = 'No Normalization')
axis(1, 
     labels = FALSE)
labels <- names(pd_raw)
text(x = seq_along(labels), 
     y = par("usr")[3] - 1, 
     srt = 45, 
     adj = 1,
     cex = 0.7,
     labels = labels, 
     xpd = TRUE)

# plotDensities(log2(pd_raw), col = rep(c("red", "green", "blue"), 6), main = "Raw data")

# format(round(colSums(pd_raw), digits = 0), big.mark = ",")


# Sample loading normalization --------------------------------------------

exp1_raw <- pd_raw[, 11:20]
exp2_raw <- pd_raw[, 1:10]

# this adjusts for differences in sampling loading and reaction efficiency of the labels
norm_facs1 <- mean(colSums(exp1_raw)) / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs1, FUN = "*")

norm_facs2 <- mean(colSums(exp2_raw)) / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs2, FUN = "*")

# sample loading normalization dataframe

data_sl <- cbind(exp1_sl, exp2_sl)

boxplot(log2(data_sl), 
        col = rep(c('firebrick', 'darkblue'), each = 10), 
        xaxt = 'n', 
        xlab = '', main = 'SL Normalization')
axis(1, 
     labels = FALSE)
labels <- names(pd_raw)
text(x = seq_along(labels), 
     y = par("usr")[3] - 0.5, 
     srt = 45, 
     adj = 1,
     cex = 0.7,
     labels = labels, 
     xpd = TRUE)

# sample loading correction helped the total ion intensity reported among different injections
# format(round(colSums(data_sl), digits = 0), big.mark = ",")

# plotDensities(log2(data_sl), col = rep(c("red", "green", "blue"), 6), main = "SL normalization")


# IRS normalization -------------------------------------------------------

# normalize with internal standard channels, across separate TMT experiments. 

irs1 <- data_sl %>% 
  dplyr::select(hl_highmn_highfe_1_1, hl_highmn_highfe_2_1)
irs2 <- data_sl %>% 
  dplyr::select(hl_highmn_highfe_1_2, hl_highmn_highfe_2_2)

irs_sums <- tibble(rowSums(irs1), rowSums(irs2))
names(irs_sums) <- c('sum1', 'sum2')

irs_sums$average <- (irs_sums$sum1 + irs_sums$sum2 )/ 2

# compute scaling factors for each protein

irs_sums$scaling_fac1 <- irs_sums$average / irs_sums$sum1
irs_sums$scaling_fac2 <- irs_sums$average / irs_sums$sum2

#make new df with normalized data, by applying TMT experiment/protein specific scaling factors

data_irs1 <- exp1_sl * irs_sums$scaling_fac1
data_irs2 <- exp2_sl * irs_sums$scaling_fac2

data_irs <- cbind(data_irs1, data_irs2)

# what does the IRS and sample loading scaling look like?
boxplot(log2(data_irs), 
        col = rep(c('firebrick', 'darkblue'), each = 10), 
        xaxt = 'n', 
        xlab = '', main = 'SL, IRS, Normalization')
axis(1, 
     labels = FALSE)
labels <- names(pd_raw)
text(x = seq_along(labels), 
     y = par("usr")[3] - 0.5, 
     srt = 45, 
     adj = 1,
     cex = 0.7,
     labels = labels, 
     xpd = TRUE)

# plotDensities(log2(data_irs), col = rep(c("red", "green", "blue"), 6), main = "IRS data")

# format(round(colSums(data_irs), digits = 0), big.mark = ",")

# now apply TMM normalization

irs_tmm <- calcNormFactors(data_irs)
data_irs_tmm <- sweep(data_irs, 2, irs_tmm, FUN = "/") # this is data after SL, IRS, and TMM on original scale

boxplot(log2(data_irs_tmm), 
        col = rep(c('firebrick', 'darkblue'), each = 10), 
        xaxt = 'n', 
        xlab = '', main = 'SL, IRS, TMM Normalization')
axis(1, 
     labels = FALSE)
labels <- names(pd_raw)
text(x = seq_along(labels), 
     y = par("usr")[3] - 0.5, 
     srt = 45, 
     adj = 1,
     cex = 0.7,
     labels = labels, 
     xpd = TRUE)

dev.off()


## exploratory pca comparing TMM additionally normalized and not TMM additionally normalized

t_data_irs_tmm <- t(data_irs_tmm)
test <- cbind(rownames(t_data_irs_tmm), t_data_irs_tmm)
pca_tmm <- prcomp(t_data_irs_tmm, center = TRUE, scale. = TRUE)
autoplot(pca_tmm, data = test, label = TRUE, title = 'IRS and TMM',  xlim = c(-0.5, 0.5))

t_data <- t(data_irs)
expression_matrix_vars <- cbind(rownames(t_data), t_data)
pca_1 <- prcomp(t_data, center = TRUE, scale. = TRUE)
autoplot(pca_1, data = expression_matrix_vars, label = T, title = 'IRS Only', xlim = c(-0.5, 0.5))


# making a nice PCA of sample-loaded, IRS, TMM, normalized expression

pca_data <- prcomp(t(data_irs_tmm), center = TRUE, scale. = TRUE)
pca_data_no_tmm <- prcomp(t(data_irs), center = TRUE, scale. = TRUE)

# making descriptive characteristics dataframe

sample_data <- t(as.data.frame(strsplit(names(data_irs_tmm), "_")))
rownames(sample_data) <- NULL
sample_df <- as.data.frame(sample_data)
names(sample_df) <- c('light', 'mn', 'fe', 'replicate', 'injection')
sample_df$metal <- paste(sample_df$mn, sample_df$fe, sep = "_")


# plotting pca

injection_pca <- autoplot(pca_data, data = sample_df, colour = 'injection', size = 3)
injection_pca2 <-  injection_pca + theme_bw() + scale_colour_discrete(name = 'TMT Experiment')

ggsave(injection_pca2, filename = "figures/sx-tmt-pca-imputed-values.png", width = 23*0.65, height = 27.94*0.5, units = "cm")

metal_pca <- autoplot(pca_data, 
                      data = sample_df, 
                      colour = 'metal', 
                      size = 3, 
                      shape = 'light',
                      alpha = 0.7)
metal_pca_notmm <- autoplot(pca_data_no_tmm, 
                            data = sample_df, 
                            colour = 'metal', 
                            size = 3, 
                            shape = 'light',
                            alpha = 0.7)
metal_pca2 <- metal_pca + scale_colour_manual(values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) + theme_bw() + ggtitle('SL, IRS, TMM Normalized')

# What did TMM normalization do?
metal_pca2_notmm <- metal_pca_notmm + 
  scale_colour_manual(values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) + 
  theme_bw() + 
  ggtitle('SL, IRS Normalized')

# grid.arrange(metal_pca2, metal_pca2_notmm)

all_pca <- metal_pca + scale_colour_manual(labels = c("High Mn, High Fe", "High Mn, Low Fe", "Low Mn, High Fe", "Low Mn, Low Fe"), 
                                           values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 
                                           name = "Metal") + 
  theme_bw() + scale_shape_discrete(labels = c("High Light", "Low Light", "Ordinary Light"), name = "Light")

ggsave(all_pca, filename = "figures/sx-all-pca-imputed_values.png", width = 23*0.65, height = 27.94*0.5, units = "cm")

## It looks like there are separate responses between treatments, but there are only 280 proteins because we cannot normalize what we cannot see. 
# So what happens when we isolate just the metal response, and just the light response

# isolating only the low light treatments with different metals
metal_columns <- grep(x = names(pd), pattern = "ll")
normalization_columns <- c(10, 11, 20, 21)

pd_metal <- pd[, c(metal_columns, normalization_columns)]
pd_metal_all <- pd_metal[complete.cases(pd_metal),]

pd_metal_prots <- pd[, c(1, metal_columns, normalization_columns)]
pd_metal_all_prots <- pd_metal_prots[complete.cases(pd_metal_prots),]

# dim(pd_metal_all) # 283 common proteins


# isolating only the light treatments all with high mn and high fe
light_columns <- grep(x = names(pd), pattern = 'hl')
light_columns2 <- grep(x = names(pd), pattern = 'ol')
light_columns3 <- grep(x = names(pd), pattern = 'll_highmn_highfe')

pd_light <- pd[, c(light_columns, light_columns2, light_columns3)]
pd_light_all <- pd_light[complete.cases(pd_light),]

pd_light_prots <- pd[, c(1, light_columns, light_columns2, light_columns3)]
pd_light_all_prots <- pd_light_prots[complete.cases(pd_light_prots),]

# dim(pd_light_all) # 625 commons proteins


# pd_light SL normalization -----------------------------------------------

sl_normalization <- function(protein_df, 
                             tmt1, 
                             tmt2, 
                             box_labels){
  
  expect_is(tmt1, 'integer')
  expect_is(tmt2, 'integer')
  expect_is(protein_df, 'data.frame')
  
  exp1_vals <- protein_df[,tmt1]
  exp2_vals <- protein_df[,tmt2]
  
  normalization_factor1 <- mean(colSums(exp1_vals)) / colSums(exp1_vals)
  exp1_vals_sl <- sweep(exp1_vals, 2, normalization_factor1, FUN = "*")
  
  normalization_factor2 <- mean(colSums(exp2_vals)) / colSums(exp2_vals)
  exp2_vals_sl <- sweep(exp2_vals, 2, normalization_factor2, FUN = "*")
  
  exp_data_sl <- cbind(exp1_vals_sl, exp2_vals_sl)
  
  col_blank <- rep('firebrick', length(box_labels))
  col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  
  boxplot(log2(exp_data_sl), 
          col = col_blank,
          xaxt = 'n', 
          xlab = '',
          main = 'B) SL Normalization',
          ylab = 'Intensity')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(box_labels), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = box_labels, 
       xpd = TRUE)
  
  list_data_sl <- list(exp1_vals_sl, exp2_vals_sl)
  
  return(list_data_sl)
}

irs_normalization <- function(sl_normalized_list, 
                              tmt1_common_channel, 
                              tmt2_common_channel,
                              tmt1,
                              tmt2,
                              box_labels){
  
  # sl_normalized_list <- tester
  # tmt1_common_channel <- c("hl_highmn_highfe_1_1", "hl_highmn_highfe_2_1")
  # tmt2_common_channel <- c("hl_highmn_highfe_1_2", "hl_highmn_highfe_2_2")
  
  sl_normalized_df_1 <- sl_normalized_list[[1]]
  irs1 <- sl_normalized_df_1[tmt1_common_channel]
  
  sl_normalized_df_2 <- sl_normalized_list[[2]]
  irs2 <- sl_normalized_df_2[tmt2_common_channel]
  # normalize with internal standard channels, across separate TMT experiments. 
  irs_sums <- tibble(rowSums(irs1), rowSums(irs2))
  names(irs_sums) <- c('sum1', 'sum2')
  
  irs_sums$average <- (irs_sums$sum1 + irs_sums$sum2) /2
  
  # compute scaling factors for each protein
  
  irs_sums$scaling_fac1 <- irs_sums$average / irs_sums$sum1
  irs_sums$scaling_fac2 <- irs_sums$average / irs_sums$sum2
  
  #make new df with normalized data, by applying TMT experiment/protein specific scaling factors
  
  data_irs1 <- sl_normalized_df_1 * irs_sums$scaling_fac1
  data_irs2 <- sl_normalized_df_2 * irs_sums$scaling_fac2
  
  data_irs <- cbind(data_irs1, data_irs2)  
  
  col_blank <- rep('firebrick', length(box_labels))
  col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  
  boxplot(log2(data_irs), 
          col = col_blank,
          xaxt = 'n', 
          xlab = '',
          main = 'C) SL, IRS Normalization',
          ylab = 'Intensity')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(box_labels), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = box_labels, 
       xpd = TRUE)
  
  return(data_irs)
  
}

tmm_normalization <- function(irs_normalized_df, 
                              tmt1, 
                              tmt2, 
                              box_labels){
  
  irs_tmm <- edgeR::calcNormFactors(irs_normalized_df)
  data_irs_tmm <- sweep(irs_normalized_df, 2, irs_tmm, FUN = "/")
  
  col_blank <- rep('firebrick', length(box_labels))
  col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  
  boxplot(log2(data_irs_tmm), 
          # col = c(rep('firebrick', length(tmt1)),
          # rep('darkblue', length(tmt2))),
          col = col_blank,
          xaxt = 'n', 
          xlab = '',
          ylab = 'Intensity',
          main = 'D) SL, IRS, TMM Normalization')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(box_labels), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = box_labels, 
       xpd = TRUE)
  
  return(data_irs_tmm)
}

# transformations for light treatments only (high mnfe)

dev.off()
# 
# # jpeg("figures/Metal-Normalization.jpeg", width=170, height=210, units="mm", res=850)
# 
# # below is an old figure not included in the MS, which shows normalization of just the light treatments
# 
# png("figures/s1-light-normalization.png", width=23*0.75, height=27.94*0.75, units="cm", res=1200)
# 
# par(mfrow = c(2, 2), mar = c(7, 6, 4, 2))
# 
# boxplot(log2(pd_light_all),
#         col = rep(c('firebrick', 'darkblue'), each = 10),
#         xaxt = 'n',
#         xlab = '', main = 'A) No Normalization')
# axis(1,
#      labels = FALSE)
# labels <- names(pd_light_all)
# text(x = seq_along(labels),
#      y = par("usr")[3] - 1,
#      srt = 45,
#      adj = 1,
#      labels = labels,
#      xpd = TRUE)
# 
sl_light <- sl_normalization(protein_df = pd_light_all,
                             tmt1 = c(4:11),
                             tmt2 = c(1:3),
                             box_labels = names(pd_light))

irs_light <- irs_normalization(sl_normalized_list = sl_light,
                               tmt1_common_channel = c("hl_highmn_highfe_1_1", "hl_highmn_highfe_2_1"),
                               tmt2_common_channel = c("hl_highmn_highfe_1_2", "hl_highmn_highfe_2_2"),
                               tmt1 = c(4:11),
                               tmt2 = c(1:3),
                               box_labels = names(pd_light))

tmm_light <- tmm_normalization(irs_normalized_df = irs_light,
                               tmt1 = c(4:11),
                               tmt2 = c(1:3),
                               box_labels = names(pd_light))
# dev.off()

tmm_light_prots <- data.frame(protein_id = pd_light_all_prots$protein_id, tmm_light)

names(sl_light[[1]]) <- paste0(names(sl_light[[1]]), "_sl")
names(sl_light[[2]]) <- paste0(names(sl_light[[2]]), "_sl")
names(irs_light) <- paste0(names(irs_light), "_irs")
names(tmm_light) <- paste0(names(tmm_light), "_tmm")

tmm_light_prots2 <- data.frame(protein_id = pd_light_all_prots$protein_id, tmm_light, irs_light, sl_light[[1]], sl_light[[2]])

# normalizations for low light, and different metal stressors only


dev.off()

# jpeg("figures/Metal-Normalization.jpeg", width=170, height=210, units="mm", res=850)
# png("figures/s1-metal-normalization.png", width=23*0.75, height=27.94*0.75, units="cm", res=1200)
# 
# par(mfrow = c(2, 2), mar = c(7, 6, 4, 2))
# 
# 
# col_blank <- rep('firebrick', labels %>% length())
# col_blank[which(grepl(pattern = "_2$", x = labels))] <- 'darkblue'
# 
# boxplot(log2(pd_metal_all), 
#         col = col_blank, 
#         xaxt = 'n', 
#         xlab = '', main = 'A) No Normalization', ylab = 'Intensity')
# axis(1, 
#      labels = FALSE)
# labels <- names(pd_metal_all)
# text(x = seq_along(labels), 
#      y = par("usr")[3] - 1, 
#      srt = 45, 
#      adj = 1,
#      labels = labels, 
#      xpd = TRUE)
# 
sl_metal <- sl_normalization(protein_df = pd_metal_all,
                             tmt1 = c(8:12, 15:16),
                             tmt2 = c(1:7, 13:14),
                             box_labels = names(pd_metal))
irs_metal <- irs_normalization(sl_normalized_list = sl_metal,
                               tmt1_common_channel = c("hl_highmn_highfe_1_1", "hl_highmn_highfe_2_1"),
                               tmt2_common_channel = c("hl_highmn_highfe_1_2", "hl_highmn_highfe_2_2"),
                               tmt1 = c(8:12, 15:16),
                               tmt2 = c(1:7, 13:14),
                               box_labels = names(pd_metal))
tmm_metal <- tmm_normalization(irs_normalized_df = irs_metal,
                               tmt1 = c(8:12, 15:16),
                               tmt2 = c(1:7, 13:14),
                               box_labels = names(pd_metal))
# 
# dev.off()

tmm_metal_prots <- data.frame(protein_id = pd_metal_all_prots$protein_id, tmm_metal)

names(sl_metal[[1]]) <- paste0(names(sl_metal[[1]]), "_sl")
names(sl_metal[[2]]) <- paste0(names(sl_metal[[2]]), "_sl")
names(irs_metal) <- paste0(names(irs_metal), "_irs")
names(tmm_metal) <- paste0(names(tmm_metal), "_tmm")

tmm_metal_prots2 <- data.frame(protein_id = pd_metal_all_prots$protein_id, tmm_metal, irs_metal, sl_metal[[1]], sl_metal[[2]])
#edgeR differential expression analysis for pd_metal and pd_light


# light-PCA ---------------------------------------------------------------

# PCA for light only treatments:

pca_light_data <- prcomp(t(tmm_light), center = TRUE, scale. = TRUE)

# making descriptive characteristics dataframe  
light_sample_data <- t(as.data.frame(strsplit(names(tmm_light), "_")))
rownames(light_sample_data) <- NULL
light_sample_df <- as.data.frame(light_sample_data)
names(light_sample_df) <- c('light', 'mn', 'fe', 'replicate', 'injection', 'normalization')
light_sample_df$metal <- paste(light_sample_df$mn, light_sample_df$fe, sep = "_")

# plotting pca
light_pca_plot <- ggplot2::autoplot(pca_light_data, 
                                    data = light_sample_df, colour = 'light', size = 3)
light_pca_plot2 <- light_pca_plot + theme_bw() + 
  scale_colour_manual(name = 'Light Treatment',
                      labels = c('High Light', 'Low Light', 'Ordinary Light'),
                      values = c("#F0E442", "#0072B2", "#D55E00"))

ggsave(light_pca_plot2, filename = "figures/sx-light-pca-imputed-values.png", width = 23*0.65, height = 27.94*0.5, units = "cm")


# edgeR DE analysis for light ---------------------------------------------

# setting up the GLMs
light_ <- factor(c("ol", "ol", "ol", "ll", "ll", "ll", "hl", "hl", "hl"))

light_dge <- DGEList(counts = tmm_light[,-c(1:2)], 
                     group = light_, 
                     genes = tmm_light_prots$protein_id)

design <- model.matrix(~0 +light_)

light_disp <- estimateDisp(light_dge, design, robust = TRUE)

fit <- glmQLFit(light_disp, design, robust = TRUE)

# Defining the contrasts to compare factor levels to 
light_contrasts1 <- makeContrasts(light_ol - light_ll, levels = design)
light_contrasts2 <- makeContrasts(light_hl - light_ll, levels = design)

qlf_1 <- glmQLFTest(fit, contrast = light_contrasts1)
topTags(qlf_1, n = 500)

qlf_2 <- glmQLFTest(fit, contrast = light_contrasts2)
topTags(qlf_2, n = 500)

# dev.off()
# png("figures/main-light-de-ma-plot.png", width=23*0.75, height=27.94*0.75, units="cm", res=800)
# 
# par(mfrow = c(2, 1))
# plotMD(qlf_1, 
#        main = 'Ordinary Light vs. Low Light', 
#        hl.col = c('firebrick', 'darkblue'), 
#        xlab = "")
# plotMD(qlf_2, 
#        main = 'High Light vs. Low Light',
#        hl.col = c('firebrick', 'darkblue'), 
#        legend = FALSE, xlab = 'Mean log(Intersity)')
dev.off()


de_light_1 <- decideTestsDGE(qlf_1)
de_light_2 <- decideTestsDGE(qlf_2)

# setting up df objects for following DE ggplots

de_light_1_tab <- qlf_1$table
de_light_1_tab$dge <- de_light_1[,1]

de_light_2_tab <- qlf_2$table
de_light_2_tab$dge <- de_light_2[,1]

# plotting of differential expression of light treatments

de_light_1_tab_p <- de_light_1_tab %>% 
  ggplot(aes(y = logFC, x = logCPM, colour = logFC)) + 
  geom_point(alpha = 0.7, shape = 19, aes(colour = dge %>% as.factor()), size = 2) +
  # geom_point(data = de_metal_1_tab %>% filter(dge == 1), alpha = 0.9) + 
  # scale_colour_gradientn(colours = c('yellow', 'yellow4', 'black', 'cyan4', 'cyan'), values = c(0, 0.15, 0.3, 0.55, 1)) +
  scale_colour_manual(values = c('yellow3', 'grey17', 'cyan2')) +
  labs(y = "Log(Fold Change)", x = '') +
  ggtitle('A) Optimal Light vs. Low Light') +
  geom_abline(slope = 0, intercept = 0) +
  theme_classic() + 
  ylim(-5, 5) +
  theme(plot.title = element_text(size = 10), legend.position = "none");de_light_1_tab_p

de_light_2_tab_p <- de_light_2_tab %>% 
  ggplot(aes(y = logFC, x = logCPM, colour = logFC)) + 
  geom_point(alpha = 0.7, shape = 19, aes(colour = dge %>% as.factor()), size = 2) +
  # geom_point(data = de_metal_1_tab %>% filter(dge == 1), alpha = 0.9) + 
  # scale_colour_gradientn(colours = c('yellow', 'yellow4', 'black', 'cyan4', 'cyan'), values = c(0, 0.15, 0.3, 0.55, 1)) +
  scale_colour_manual(values = c('yellow3', 'grey17', 'cyan2')) +
  labs(y = "", x = 'Mean log(Intensity)') +
  ggtitle('B) High Light vs. Low Light') +
  labs(y = "Log(Fold Change)", x = 'Mean log(Intensity)') +
  geom_abline(slope = 0, intercept = 0) +
  theme_classic() + 
  ylim(-5, 5) +
  theme(plot.title = element_text(size = 10), legend.position = "none");de_light_2_tab_p

tiff("figures/s5-de-light-plot-imputed-values.tiff", width=23*0.75, height=27.94*0.75, units="cm", res=800)

grid.arrange(de_light_1_tab_p, de_light_2_tab_p, nrow = 2)

dev.off()


# summary(de_light_1)
# summary(de_light_2)

# metal-PCA ---------------------------------------------------------------

# PCA for metal only treatments:

high_light_metal <- grepl(pattern = "ll", x = names(tmm_metal))
tmm_metal2 <- tmm_metal[,high_light_metal]

pca_metal_data <- prcomp(t(tmm_metal2), center = TRUE, scale. = TRUE)

# making descriptive characteristics dataframe  
metal_sample_data <- t(as.data.frame(strsplit(names(tmm_metal2), "_")))
rownames(metal_sample_data) <- NULL
metal_sample_df <- as.data.frame(metal_sample_data)
names(metal_sample_df) <- c('light', 'mn', 'fe', 'replicate', 'injection', 'normalization')
metal_sample_df$metal <- paste(metal_sample_df$mn, metal_sample_df$fe, sep = "_")

# plotting pca
metal_pca_plot <- ggplot2::autoplot(pca_metal_data, 
                                    data = metal_sample_df, colour = 'metal', size = 3, alpha = 0.8)
metal_pca_plot2 <- metal_pca_plot + theme_bw() + 
  scale_colour_manual(name = 'Metal',
                      labels = c("High Mn, High Fe", "High Mn, Low Fe", "Low Mn, High Fe", "Low Mn, Low Fe"),
                      values = c("#F0E442", "#0072B2", "#D55E00", "grey50"))

ggsave(metal_pca_plot2, filename = "figures/sx-metal-pca-imputed-values.png", width = 23*0.65, height = 27.94*0.5, units = "cm")


# edgeR for metals --------------------------------------------------------

tmm_metal_sub <- tmm_metal[-c(6, 7, 15, 16)]

metal_ <- factor(c("highmn_highfe", "highmn_highfe", "highmn_highfe", "lowmn_highfe", "lowmn_highfe", "lowmn_highfe", 
                   "highmn_lowfe", "highmn_lowfe", "highmn_lowfe", "lowmn_lowfe", "lowmn_lowfe", "lowmn_lowfe"))

metal_dge <- DGEList(counts = tmm_metal_sub, 
                     group = metal_, 
                     genes = tmm_metal_prots$protein_id)

design_metal <- model.matrix(~0 + metal_)

metal_disp <- estimateDisp(metal_dge, design_metal, robust = TRUE)

fit_metal <- glmQLFit(metal_disp, design_metal, robust = TRUE)

# Defining the contrasts to compare factor levels to 
metal_contrasts1_low <- makeContrasts(metal_highmn_highfe - metal_lowmn_lowfe, levels = design_metal)
metal_contrasts2_low <- makeContrasts(metal_highmn_lowfe - metal_lowmn_lowfe, levels = design_metal)
metal_contrasts3_low <- makeContrasts(metal_lowmn_highfe - metal_lowmn_lowfe, levels = design_metal)

qlf_metal_1_low <- glmQLFTest(fit_metal, contrast = metal_contrasts1_low)
qlf_metal_2_low <- glmQLFTest(fit_metal, contrast = metal_contrasts2_low)
qlf_metal_3_low <- glmQLFTest(fit_metal, contrast = metal_contrasts3_low)

metal_contrasts1 <- makeContrasts(metal_lowmn_lowfe - metal_highmn_highfe, levels = design_metal)
metal_contrasts2 <- makeContrasts(metal_highmn_lowfe - metal_highmn_highfe, levels = design_metal)
metal_contrasts3 <- makeContrasts(metal_lowmn_highfe - metal_highmn_highfe, levels = design_metal)

qlf_metal_1 <- glmQLFTest(fit_metal, contrast = metal_contrasts1)
# topTags(qlf_metal_1, n = 500)

qlf_metal_2 <- glmQLFTest(fit_metal, contrast = metal_contrasts2)
# topTags(qlf_metal_2, n = 500)

qlf_metal_3 <- glmQLFTest(fit_metal, contrast = metal_contrasts3)
# tester <- topTags(qlf_metal_3, n = 500)

dev.off()

up_colour <- '#00CCCC'
down_colour <- '#FFFF00'

# png("figures/s4-metal-de-ma-plot.png", width=23*0.60, height=27.94*0.60, units="cm", res=800)
# 
# 
# par(mfrow = c(3, 1))
# plotMD(qlf_metal_1,
#        main = 'Low Mn + Low Fe vs. High Mn + High Fe',
#        hl.col = c(up_colour, down_colour),
#        xlab = "")
# plotMD(qlf_metal_2,
#        main = 'High Mn + Low Fe vs. High Mn + High Fe',
#        hl.col = c(up_colour, down_colour),
#        legend = FALSE, xlab = "")
# plotMD(qlf_metal_3,
#        main = 'Low Mn + High Fe vs. High Mn + High Fe',
#        hl.col = c(up_colour, down_colour),
#        legend = FALSE, xlab = 'Mean log(Intensity)')
# dev.off()

de_metal_1 <- decideTestsDGE(qlf_metal_1)
de_metal_2 <- decideTestsDGE(qlf_metal_2)
de_metal_3 <- decideTestsDGE(qlf_metal_3)


de_metal_1_low <- decideTestsDGE(qlf_metal_1_low)
de_metal_2_low <- decideTestsDGE(qlf_metal_2_low)
de_metal_3_low <- decideTestsDGE(qlf_metal_3_low)

summary(de_metal_1)
summary(de_metal_2)
summary(de_metal_3)

# composite DE figure of metals, four panel

de_metal_1_tab <- qlf_metal_1$table
de_metal_1_tab$dge <- de_metal_1[,1]

de_metal_2_tab <- qlf_metal_2$table
de_metal_2_tab$dge <- de_metal_2[,1]

de_metal_3_tab <- qlf_metal_3$table
de_metal_3_tab$dge <- de_metal_3[,1] 

de_metal_2_low_tab <- qlf_metal_2_low$table
de_metal_2_low_tab$dge <- de_metal_2_low[,1]

de_metal_1_tab_p <- de_metal_1_tab %>% 
  ggplot(aes(y = logFC, x = logCPM, colour = logFC)) + 
  geom_point(alpha = 0.7, shape = 19, aes(colour = dge %>% as.factor()), size = 2) +
  # geom_point(data = de_metal_1_tab %>% filter(dge == 1), alpha = 0.9) + 
  # scale_colour_gradientn(colours = c('yellow', 'yellow4', 'black', 'cyan4', 'cyan'), values = c(0, 0.15, 0.3, 0.55, 1)) +
  scale_colour_manual(values = c('yellow3', 'grey17', 'cyan2')) +
  labs(y = "Log(Fold Change)", x = '') +
  ggtitle('A) Low Mn + Low Fe vs. High Mn + High Fe') +
  geom_abline(slope = 0, intercept = 0) +
  theme_classic() + 
  ylim(-5, 5) +
  theme(plot.title = element_text(size = 10), legend.position = "none");de_metal_1_tab_p

de_metal_2_tab_p <- de_metal_2_tab %>% 
  ggplot(aes(y = logFC, x = logCPM, colour = logFC)) + 
  geom_point(alpha = 0.7, shape = 19, aes(colour = dge %>% as.factor()), size = 2) +
  # geom_point(data = de_metal_1_tab %>% filter(dge == 1), alpha = 0.9) + 
  # scale_colour_gradientn(colours = c('yellow', 'yellow4', 'black', 'cyan4', 'cyan'), values = c(0, 0.15, 0.3, 0.55, 1)) +
  scale_colour_manual(values = c('yellow3', 'grey17', 'cyan2')) +
  labs(y = "", x = '') +
  ggtitle('B) High Mn + Low Fe vs. High Mn + High Fe') +
  geom_abline(slope = 0, intercept = 0) +
  theme_classic() + 
  ylim(-5, 5) +
  theme(plot.title = element_text(size = 10), legend.position = "none");de_metal_2_tab_p

de_metal_3_tab_p <- de_metal_3_tab %>% 
  ggplot(aes(y = logFC, x = logCPM, colour = logFC)) + 
  geom_point(alpha = 0.7, shape = 19, aes(colour = dge %>% as.factor()), size = 2) +
  # geom_point(data = de_metal_1_tab %>% filter(dge == 1), alpha = 0.9) + 
  # scale_colour_gradientn(colours = c('yellow', 'yellow4', 'black', 'cyan4', 'cyan'), values = c(0, 0.15, 0.3, 0.55, 1)) +
  scale_colour_manual(values = c("grey17", 'cyan2')) +
  labs(y = "Log(Fold Change)", x = 'Mean log(Intensity)') +
  ggtitle('C) Low Mn + High Fe vs. High Mn + High Fe') +
  geom_abline(slope = 0, intercept = 0) +
  theme_classic() +
  ylim(-5, 5) +
  theme(plot.title = element_text(size = 10), legend.position = "none");de_metal_3_tab_p

de_metal_2_low_p <- de_metal_2_low_tab %>% 
  ggplot(aes(y = logFC, x = logCPM, colour = logFC)) + 
  geom_point(alpha = 0.7, shape = 19, aes(colour = dge %>% as.factor()), size = 2) +
  # geom_point(data = de_metal_1_tab %>% filter(dge == 1), alpha = 0.9) + 
  # scale_colour_gradientn(colours = c('yellow', 'yellow4', 'black', 'cyan4', 'cyan'), values = c(0, 0.15, 0.3, 0.55, 1)) +
  scale_colour_manual(values = c('yellow3', 'grey17', 'cyan2')) +
  # scale_colour_gradientn(colours = orange_palette(20)) + 
  labs(y = "", x = 'Mean log(Intensity)') +
  ggtitle("D) High Mn + Low Fe vs. Low Mn + Low Fe") +
  geom_abline(slope = 0, intercept = 0) +
  theme_classic() + 
  ylim(-5, 5) +
  theme(plot.title = element_text(size = 10),
        legend.position = "none",
        text = element_text(family = ''));de_metal_2_low_p


tiff("figures/s4-de-plot-imputed-values.tiff", width=23*0.75, height=27.94*0.75, units="cm", res=800)

grid.arrange(de_metal_1_tab_p, de_metal_2_tab_p, de_metal_3_tab_p, de_metal_2_low_p)

dev.off()

png("figures/s3-de-venn-diagrams-imputed-values.png", width=23*0.75, height=27.94*0.75, units="cm", res=800)

par(mfrow = c(2, 1), 
    mai = c(0, 0, 0, 0))

# making light venn diagram

v_presence_light <- cbind(de_light_1[,"-1*light_ll 1*light_ol"], 
                          de_light_2[,"1*light_hl -1*light_ll"])

vcounts_light <- vennCounts(v_presence_light)
vennDiagram(vcounts_light, 
            names = c('High Light', 
                      'Ordinary Light'),
            cex = 1)

# making metal venn diagram

v_presence_metal <- cbind(de_metal_2[,"-1*metal_highmn_highfe 1*metal_highmn_lowfe"], 
                          de_metal_3[,"-1*metal_highmn_highfe 1*metal_lowmn_highfe"],
                          de_metal_1[,"-1*metal_highmn_highfe 1*metal_lowmn_lowfe"])
vcounts <- vennCounts(v_presence_metal)
vennDiagram(vcounts, 
            names = c('High Mn + Low Fe', 'Low Mn + High Fe', 'Low Mn + Low Fe'),
            cex = 1)

dev.off()

v_presence_metal_low <- cbind(de_metal_1_low[, "1*metal_highmn_highfe -1*metal_lowmn_lowfe"],
                              de_metal_2_low[, "1*metal_highmn_lowfe -1*metal_lowmn_lowfe"],
                              de_metal_3_low[, "1*metal_lowmn_highfe -1*metal_lowmn_lowfe"])


## exporting differential expression results

metal_names <- c('protein_id', 'highmn_lowfe_repl', 'lowmn_highfe_repl', 'lowmn_lowfe_repl', 'highmn_highfe_depl', 'highmn_lowfe_depl', 'lowmn_highfe_depl')
de_metal <- data.frame(qlf_metal_1$genes, v_presence_metal, v_presence_metal_low)
names(de_metal) <- metal_names

light_names <- c('protein_id', 'ol', 'hl')
de_light <- data.frame(qlf_1$genes, v_presence_light)
names(de_light) <- light_names

de_metal2 <- inner_join(de_metal, just_annot2, by = 'protein_id')
de_metal3 <- inner_join(de_metal2, tmm_metal_prots2, by = 'protein_id')

de_light2 <- inner_join(de_light, just_annot2, by = 'protein_id')
de_light3 <- inner_join(de_light2, tmm_light_prots2, by = 'protein_id')

write.csv(de_light3, 'data/intermediate-data/de-light-imputed-values.csv')
write.csv(de_metal3, 'data/intermediate-data/de-metal-imputed-values.csv')


# all normalized values plus the DE values from the metal treatments:

all_norm_pd <- data.frame(protein_id = pd_all$protein_id, data_irs_tmm)
# adding in the annotation file
all_norm_annot <- inner_join(all_norm_pd, just_annot2, by = 'protein_id')

# adding in the DE values from just the metal
# subset of DE vals from metal de test

de_metal3_sub <- de_metal3[, c(1:7)]
all_norm_annot_de <- inner_join(all_norm_annot, de_metal3_sub, by = "protein_id")

de_light3_sub <- de_light3[, c(1:3)]
all_norm_annot_de_lightmetal <- inner_join(all_norm_annot_de, de_light3_sub, by = "protein_id")

write.csv(all_norm_annot_de_lightmetal, file = 'data/intermediate-data/all_prots_normalized-imputed-values.csv')





