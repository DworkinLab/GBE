# Analysis of the mouse subset of the GBE data 

# Libraries 
library(lme4)
library(car)
library(MASS)
library(effects)
library(multcomp)
library(dplyr)
library(ggplot2)
library(knitr)
library(data.table)
library(MCMCglmm)
library(tidyverse)
library(brms)
library(lme4)
library(tidybayes)
library(bayesplot)
library(GGally)
library(randomForest)

# Working directory
setwd("~/Dropbox/Background Effects Review DworkinMerritt labs/GenBackMetaAnalysis/data")   
# Working directory for big Mac
setwd("/Volumes/Nirupama")

# read in data
gb_dat <- read.csv("gb_dat_updated.csv", na.strings= "", header = T)

# Data checks
str(gb_dat)
summary(gb_dat)
dim(gb_dat)

# Create a study variable.
gb_dat$study <- with(gb_dat,interaction(author, year, drop=TRUE))

# replace empty fields in age.specific with "none"
gb_dat$age.specific <- as.factor(gb_dat$age.specific)
levels <- levels(gb_dat$age.specific)
levels[length(levels) + 1] <- "None"
gb_dat$age.specific <- factor(gb_dat$age.specific, levels = levels)
gb_dat$age.specific[is.na(gb_dat$age.specific)] <- "None"

###############  Cleaning up the unspecific sample size column and sd column ############################
# There are some weird date values in the sample.size.unspecific, replacing those with NA
# Although coalesce and converting to numeric will not keep these anyway 
gb_dat$sample.size.unspecific <- str_replace_all(gb_dat$sample.size.unspecific,
                                                 c("Jun","Oct","Aug"), NA_character_)

# Seem to be three cases of unspecific sample size

#  First: >, i.e lower bound of sample size, remove the > and use as the sample size
gb_dat$sample.size.unspecific <- str_replace_all(gb_dat$sample.size.unspecific,
                                                 ">","")

# Second: average sample size used is given, remove avg. and use as the sample size
gb_dat$sample.size.unspecific <- str_replace_all(gb_dat$sample.size.unspecific,
                                                 "avg.","")

# Range of sample sizes: split columns based on "-" and use minimum
gb_dat <- separate(data = gb_dat,
                   col = sample.size.unspecific, 
                   into = c("Min.Size","Max.Size"), 
                   sep = "-")

# Can now remove max size column
gb_dat <- gb_dat %>% select(-Max.Size)

# Convert min size to numeric column (to combine the two sample size colummns)
gb_dat$Min.Size <- as.numeric(gb_dat$Min.Size)

# Coalesce sample size unspecific and specific columns together to make one sample size column
gb_dat <- gb_dat %>% 
  # To replace NA values in sample.size by Min.Size
  mutate( Sample.Size = coalesce(sample.size, Min.Size)) %>%
  select(-sample.size, -Min.Size)

# To remove the studies that did not report any form of sample size 
# If we want to remove these
gb_dat <- gb_dat %>%
 drop_na(Sample.Size)

summary(gb_dat$Sample.Size)

dim(gb_dat)

# Converting the standard errors to standard deviations 
gb_dat$se <- with(gb_dat, se * sqrt(Sample.Size))

# Combining the se and sd columns for one SD column 
gb_dat <- gb_dat %>% 
  mutate(SD = coalesce(se,sd)) %>%
  select(- se, -sd)

# There seem to be some negative values in the standard deviation column as well, assuming the neg sign was accidentally entered
gb_dat$SD <- with(gb_dat, abs(SD))


summary(gb_dat$SD)
#######################################################################################
# Adding a column with info about the number of strains used in each study
strain_counts <- gb_dat %>%
  group_by(study) %>%
  summarise (n_distinct(strain))

colnames(strain_counts) <- c("study","strain.count")

gb_dat <- left_join(gb_dat, strain_counts, by = "study")

# Converting to factor 
gb_dat[,c(5,8,9,16,17,18,19,20,21,27)] <- lapply(gb_dat[,c(5,8,9,16,17,18,19,20,21,27)], factor)
############################ Some  summary information ###############################
# Total number of studies (147 in total)
tot_studies <- gb_dat %>%
  group_by(study) %>%
  summarize(count = n())

dim(tot_studies)

# Frequency of strain numbers 
strain_counts_freq <- strain_counts %>%
  group_by(number) %>%
  count(number) %>%
  ungroup()

# Combined table 
study_summary <- left_join(tot_studies, strain_counts, by = "study")


# 1656 effect sizes were not reported with a sample size 
sum(is.na(gb_dat$Sample.Size))

# Studies that were not reported with a sample size 
no_sample_size <- gb_dat %>%
  filter(is.na(Sample.Size)) %>%
   group_by(study) %>%
  summarize(count = n()) 

# Number of studies with time-dependent values 
time_recorded_studies <- time_recorded_values %>%
  group_by(study) %>%
  tally() 

# Number of studies with expression data 
expression_dat <- gb_dat %>%
  filter(phenotype.type == "expression")  %>%
  group_by(study) %>%
  tally() 

# Number of studies with hybrid strains
hybrid_strains <- gb_dat %>%
  filter(hybrid == "y") %>%
  group_by(study) %>%
  tally() 
  

# Looking at how many studies reported neither standard deviation or standard error estimates
no_se_sd <- gb_dat %>% 
  filter(is.na(se) & is.na(sd)) %>%
  group_by(study) %>%
  tally() 

dim(expression_dat)

rm(tot_studies, strain_counts_freq, study_summary, no_sample_size, time_recorded_studies,
   expresion_dat, hybrid_strains, no_se_sd)

########################## Re-classification of alleletypes ##################################
# Changing line.type name to zygosity (as it is slightly confusing)
colnames(gb_dat)[colnames(gb_dat) == "line.type"] <- "zygosity"

# Changing the amorphs to null
levels(gb_dat$alleletype) <- c("null" ,"antimorph" ,"epimutation","gain of function" ,"hypermorph" ,
                              "hypomorph" ,"neomorph" ,"null" ,"unspecified")

# Re-classifying some of the unspecified and gain of function studies (usable ones which do have wildtype values)

# Dworkin-b: re-classify to hypermorph (only gain of function study)
gb_dat$alleletype[gb_dat$author == "dworkin-b"] = "hypermorph"

# alleletype studies to re-classify
# Only studies which have a wildtype value out of these (or have the wildtype value incorporated)
reclass_studies <- gb_dat %>%
  filter(alleletype %in% c("unspecified", "gain of function","epimutation","antimorph")) %>%
  filter(zygosity == "wildtype"| wt.incorporated == "yes") %>%
  group_by(study) %>%
  tally()

# Re-classifying the unspecified studies which do have wt values
gb_dat$alleletype[gb_dat$author == "al-saktawi et al."] = "hypomorph"

gb_dat$alleletype[gb_dat$author == "cendelin et al."] = "hypermorph"

gb_dat$alleletype[gb_dat$author == "chung et al."] = "hypomorph"

gb_dat$alleletype[gb_dat$author == "lachance et al."] = "hypomorph"

gb_dat$alleletype[gb_dat$study == "riedl et al. .2005"] = "hypomorph"

gb_dat$alleletype[gb_dat$author == "sheehan-rooney et al."] = "hypomorph"

gb_dat$alleletype[gb_dat$author == "spencer & promislow"] = "hypomorph"

gb_dat$alleletype[gb_dat$author == "de belle & heisenberg"] = "hypomorph"

gb_dat$alleletype[gb_dat$author == "de belle & heisenberg"] = "hypomorph"

gb_dat$alleletype[gb_dat$author == "torabi & kruglyak"] = "null"

gb_dat$alleletype[gb_dat$author == "wanat et al."] = "hypomorph"

gb_dat$alleletype[gb_dat$study == "wang et al..2013"] = "null"

gb_dat$alleletype[gb_dat$author == "xing et al."] = "hypomorph"

# Remove unusable studies (these studies were looked at again for re-classification and deemed unusable)

gb_dat <-gb_dat[!(gb_dat$author =="gerke et al."),]

gb_dat <-gb_dat[!(gb_dat$author =="vogwill et al."),]

# Remove all other studies with unspecified (as these do not have wt values and will not be used)
gb_dat <-gb_dat[!(gb_dat$alleletype =="unspecified"),]

# Droplevels 
gb_dat$alleletype <- droplevels(gb_dat$alleletype)

# Checks
dim(gb_dat)

table(gb_dat$alleletype)

rm(reclass_studies)

################### Creating grouping variables so that traits can be separated #####################
# With all the traits 
gb_dat$study_phenotype <- with(gb_dat, 
                               interaction(author, year, mutation, phenotype, age.specific, environment, sex, strain,
                                           zygosity, drop=TRUE))
# So that line.types are not grouped together 
gb_dat$study_phenotype_2 <- with(gb_dat, 
                                 interaction(author, year, mutation, phenotype, age.specific, environment, sex, strain,
                                             drop=TRUE))
# So that strains and line.types are not grouped together 
gb_dat$study_phenotype_3 <- with(gb_dat, 
                                 interaction(author, year, mutation, phenotype, age.specific, environment, sex,
                                             drop=TRUE))
# So that only strains are not grouped together 
gb_dat$study_phenotype_4 <- with(gb_dat, 
                                 interaction(author, year, mutation, phenotype, age.specific, environment, sex, zygosity,
                                             drop=TRUE))

################### Looking into wildtype-incorporated studies #######################
# Filtering
wildtype_incorporated <- gb_dat %>%
  filter(wt.incorporated == "yes")

dim(wildtype_incorporated)

# No. of studies (15)
no_wildtype_incorporated <- wildtype_incorporated %>%
  group_by(study) %>%
  tally()

table(wildtype_incorporated$phenotype.units)

# Removing them from the data-set for now and doing separate calculations for the ratios
gb_dat <- gb_dat %>%
  filter(wt.incorporated == "no")

# 
wildtype_incorporated_1 <- wildtype_incorporated  %>%
  filter(phenotype.units == "na")

# Doing mouse and non-mouse studies separately as mouse has the additional strain-level analysis 
wildtype_incorporated_no_mus <- wildtype_incorporated %>%
  filter(species != "mouse") 

wildtype_incorporated_mus <- wildtype_incorporated %>%
  filter(species == "mouse")

# Non-mouse studies 

# Relative studies (so can only be used for logRR ratio)
wildtype_relative_no_mus <- wildtype_incorporated_no_mus %>%
  filter(phenotype.units == "relative") 

# Calculating weighted phenotype measures 
wildtype_relative_no_mus$weighted_phenotype <- with(wildtype_relative_no_mus, 
                                                 (phenotype.measure)*Sample.Size)

#  Calculate weighted mean 
wildtype_relative_no_mus <- wildtype_relative_no_mus %>%
  group_by(study_phenotype_3) %>%
  mutate(wt_sum = sum(weighted_phenotype)) %>%
  mutate(total_wt_SS = sum(Sample.Size)) %>%
  mutate(logRR = log2(wt_sum/total_wt_SS)) %>%
  select(-wt_sum, -strain) %>%
  filter(logRR != "Inf" & logRR != "-Inf" & logRR != "NaN") %>%
  distinct(logRR, .keep_all = TRUE) %>%
  ungroup()

# Ratio studies (4: can't be used)
wildtype_ratio_no_mus <- wildtype_incorporated_no_mus %>%
  filter(phenotype.units == c("ratio","log-transformed ratio")) %>%
  group_by(study) %>%
  tally()

# Percent difference (1: can't be used)
wildtype_percdiff_no_mus <- wildtype_incorporated_no_mus %>%
  filter(phenotype.units == "percent difference") %>%
  group_by(study) %>%
  tally()

# Mouse studies 

# Averaged across strains

# Relative studies (so can only be used for logRR ratio)
wildtype_relative_mus <- wildtype_incorporated_mus %>%
  filter(phenotype.units == "relative") 

# Calculating weighted phenotype measures 
wildtype_relative_mus$weighted_phenotype <- with(wildtype_relative_mus, 
                                                    (phenotype.measure)*Sample.Size)

#  Calculate weighted mean 
wildtype_relative_mus <- wildtype_relative_mus %>%
  group_by(study_phenotype_3) %>%
  mutate(wt_sum = sum(weighted_phenotype)) %>%
  mutate(total_wt_SS = sum(Sample.Size)) %>%
  mutate(logRR = log2(wt_sum/total_wt_SS)) %>%
  select(-wt_sum, -strain) %>%
  filter(logRR != "Inf" & logRR != "-Inf" & logRR != "NaN") %>%
  distinct(logRR, .keep_all = TRUE) %>%
  ungroup()

  
# Percent studies (with removing wildtype values as they were extracted from figures but represent 100%)
wildtype_perc_mus <- wildtype_incorporated_mus %>%
  filter(phenotype.units == c("percent","percent of wildtype")) %>%
  filter(line.type != "wildtype")
  

# Back from percentage
wildtype_perc_mus$phenotype.measure <- with(wildtype_perc_mus, (phenotype.measure)/100)

# Calculating weighted phenotype measures 
wildtype_perc_mus$weighted_phenotype <- with(wildtype_perc_mus, 
                                                 (phenotype.measure)*Sample.Size)

#  Calculate weighted mean 
wildtype_perc_mus<- wildtype_perc_mus %>%
  group_by(study_phenotype_3) %>%
  mutate(wt_sum = sum(weighted_phenotype)) %>%
  mutate(total_wt_SS = sum(Sample.Size)) %>%
  mutate(logRR = log2(wt_sum/total_wt_SS)) %>%
  select(-wt_sum, -strain) %>%
  filter(logRR != "Inf" & logRR != "-Inf" & logRR != "NaN") %>%
  distinct(logRR, .keep_all = TRUE) %>%
  ungroup()

####################### Calculating correlations ##############################
# Only using the studies which used 5 or more strains 
correlation_data <- gb_dat %>%
  filter(strain.count > 5)

# Number of studies (18)
num_correlation_studies <- correlation_data %>%
  group_by(study) %>%
  tally()

dim(num_correlation_studies)


# Using only the studies with wildtype values (6 studies)
correlation_data_wt <- correlation_data %>%
  filter(zygosity == "wildtype") %>%
  group_by(study) %>%
  tally()

# Subsetting these 
correlation_data <- subset(correlation_data, study %in% correlation_data_wt$study)

# Removing the studies which do not have the same number of wildtypes and mutants 
correlation_data <- correlation_data[!(correlation_data$author =="milloz et al."),]

correlation_data <- correlation_data[!(correlation_data$author =="paaby & schmidt"),]

# Calculating the correlation for each of the six studies 
correlation_data <- correlation_data %>%
  group_by(study) %>%
  summarize(study_cor=cor(phenotype.measure[zygosity == "wildtype"],
                          phenotype.measure[zygosity != "wildtype"]))


rm(num_correlation_studies, correlation_data_wt)
####################### Mouse Subset Analysis ###################################
# Mouse subset
mouse_subset <- gb_dat %>%
  filter(species == "mouse") %>%
  droplevels()

dim(mouse_subset)

# Not using hybrids for the mouse subset 
mouse_subset <- mouse_subset %>%
  filter(hybrid == "n")

# Total number of studies in mouse_subset 
mouse_subset_no <- mouse_subset %>%
  group_by(study) %>%
  tally() %>%
  ungroup()

dim(mouse_subset_no)

###################### Strain consolidation for the mouse subset ############################
# For the mouse subset, since the same strains are used commonly in a lot of studies, the effect sizes will be calculated per strain 
# Convert to factor
mouse_subset$strain <- factor(mouse_subset$strain)

levels(mouse_subset$strain)
table(mouse_subset$strain)
nlevels(mouse_subset$strain)

# Consolidating certain sub-strains under "umbrella" strains 

# 129 Strains 
levels(mouse_subset$strain) <- gsub("129.*$", "129", 
                                    levels(mouse_subset$strain))

# c57bl/6 strains: classify into c57bl/6 
# c57bl/ksj cannot be combined as was compared to c57bl/6 in study
levels(mouse_subset$strain) <- gsub("c57bl/6.*$", "c57bl/6", 
                                    levels(mouse_subset$strain))
# Different notations
levels(mouse_subset$strain) <- gsub("c57bl6.*$", "c57bl/6", 
                                    levels(mouse_subset$strain))

levels(mouse_subset$strain) <- gsub("c57bl/10.*$", "c57bl/6", 
                                    levels(mouse_subset$strain))
# For c57bl/sj
levels(mouse_subset$strain) <- gsub("c57bl/sj.*$", "c57bl/6", 
                                    levels(mouse_subset$strain))

table(mouse_subset$strain)

# Classify all strains except 129 and c57bl/6 as other
fun <- function(z) {
  z[z == "129"] <- "129"
  z[z == "c57bl/6"] <- "c57bl/6"
  z[!(z %in% c("129", "c57bl/6"))] <- "Other"
  z
}
mouse_subset$strain <- fct_relabel(mouse_subset$strain, fun)

table(mouse_subset$strain)

##################### Effect Sizes Calculations for the Mouse Subset ########################
# Creating separate columns for the means, sds and sample sizes 
mouse_subset <- mouse_subset %>%
  group_by(study_phenotype_2) %>%
  mutate(wt_mean = ifelse(length(phenotype.measure[zygosity == "wildtype"]) > 0,
                          ((phenotype.measure[zygosity == "wildtype"])), NA)) %>%
  mutate(wt_sd = ifelse(length(phenotype.measure[zygosity == "wildtype"]) > 0,
                        ((SD[zygosity == "wildtype"])), NA)) %>%
  mutate(wt_SS = ifelse(length(phenotype.measure[zygosity == "wildtype"]) > 0,
                        (mean(Sample.Size[zygosity == "wildtype"])), NA)) %>%
  ungroup() %>%
  group_by(study_phenotype) %>%
  mutate(mutant_mean = ifelse(length(phenotype.measure[zygosity != "wildtype"]) > 0,
                              ((phenotype.measure[zygosity != "wildtype"])), NA)) %>%
  mutate(mutant_sd = ifelse(length(phenotype.measure[zygosity != "wildtype"]) > 0,
                            ((SD[zygosity != "wildtype"])), NA)) %>%
  mutate(mutant_SS = ifelse(length(phenotype.measure[zygosity != "wildtype"]) > 0,
                            ((Sample.Size[zygosity != "wildtype"])), NA)) %>%
  ungroup()

dim(mouse_subset)

# logRR ratio
mouse_subset_RR <- mouse_subset %>%
  mutate(logRR = (mutant_mean/wt_mean)) %>%
  filter(logRR != "Inf" & logRR != "-Inf" & logRR != "NaN")  

# Taking the abs value because we are relatively more interested in the magnitude   
mouse_subset_RR$logRR <- with(mouse_subset_RR, abs(logRR))

# Replacing zeroes to fit a Gamma distribution 
mouse_subset_RR$logRR <- with(mouse_subset_RR, ifelse(logRR == 0, 0.00000000000000001, logRR)) 

# log VR ratio
mouse_subset_VR <- mouse_subset %>%
  mutate(logVR = log2(mutant_sd/wt_sd) + 1/(2*(mutant_SS - 1)) - 1/(2*(wt_SS - 1))) %>%
  filter(logVR != "Inf" & logVR != "-Inf" & logVR != "NaN") 

# Taking the abs value because we are relatively more interested in the magnitude   
mouse_subset_VR$logVR <- with(mouse_subset_VR, abs(logVR))

# Replacing 0 values with very small values (to fit Gamma model)
mouse_subset_VR$logVR <- with(mouse_subset_VR, ifelse(logVR == 0, 0.00000000000000001, logVR)) 

# log CVR ratio
mouse_subset_CVR <- mouse_subset %>%
  mutate(logCVR = log2((mutant_sd/mutant_mean)/(wt_sd/wt_mean)) + 1/(2*(mutant_SS - 1)) 
         - 1/(2*(wt_SS - 1))) %>%
  filter(logCVR != "Inf" & logCVR != "-Inf" & logCVR != "NaN") 

# Taking the abs value because we are relatively more interested in the magnitude   
mouse_subset_CVR$logCVR <- with(mouse_subset_CVR, abs(logCVR))

# Replacing 0 values with very small values (to fit Gamma model)
mouse_subset_CVR$logCVR <- with(mouse_subset_CVR, ifelse(logCVR == 0, 0.00000000000000001, 
                                                         logCVR)) 

# Number of studies for mouse_subset logRR 
no_studies_RR <- mouse_subset_RR %>%
  group_by(study) %>%
  tally() %>%
  ungroup()

dim(no_studies_RR)

# Number of studies for mouse_subset logVR 
no_studies_VR <- mouse_subset_VR %>%
  group_by(study) %>%
  tally() %>%
  ungroup()

dim(no_studies_VR)

# Number of studies for mouse_subset logCVR 
no_studies_CVR <- mouse_subset_CVR %>%
  group_by(study) %>%
  tally() %>%
  ungroup()

dim(no_studies_CVR)

# Weights 

# logRR ratio
mouse_subset_RR <- mouse_subset_RR %>%
  mutate(logRR_var = (wt_sd)^2/((wt_SS)*(wt_mean)^2) + 
           (mutant_sd)^2/((mutant_SS)*(mutant_mean)^2))

# logVR ratio
mouse_subset_VR <- mouse_subset_VR %>%
  mutate(logVR_var = 1/(2*((mutant_SS) - 1)) + 1/(2*((wt_SS) - 1)))

# logCVR ratio 

# Correlation between wildtype mean and wildtype SD
wt_cor <- cor(mouse_subset$phenotype.measure[mouse_subset$zygosity == "wildtype"], 
    mouse_subset$SD[mouse_subset$zygosity == "wildtype"],
    use = "complete.obs")

# Correlation between mutant mean and mutant SD
mutant_cor <- cor(mouse_subset$phenotype.measure[mouse_subset$zygosity != "wildtype"], 
              mouse_subset$SD[mouse_subset$zygosity != "wildtype"],
              use = "complete.obs")

mouse_subset_CVR <- mouse_subset_CVR %>%
  mutate(logCVR_var = (wt_sd)^2/((wt_SS)*(wt_mean)^2) + 1/(2*((mutant_SS) - 1)) + 
           (mutant_sd)^2/((mutant_SS)*(mutant_mean)^2) + 1/(2*((wt_SS) - 1)) - 
           2*wt_cor*sqrt(((wt_sd)^2/((wt_SS)*(wt_mean)^2))*(1/(2*((wt_SS) - 1)))) -
           2*mutant_cor*sqrt(((mutant_sd)^2/((mutant_SS)*(mutant_mean)^2))*(1/(2*((mutant_SS) - 1)))))

# Some distribution plots
logRR_density <- density(mouse_subset_RR$logRR) 
plot(logRR_density) 

logVR_density <- density(mouse_subset_VR$logVR) 
plot(logVR_density) 

logCVR_density <- density(mouse_subset_CVR$logCVR) 
plot(logCVR_density) 

########################  For entire data-set ####################################
# Here we are averaging across strains as some studies have used multiple 

# Calculating weighted phenotype measures 
gb_dat$weighted_phenotype <- with(gb_dat, (phenotype.measure)*Sample.Size)

#  Calculate weighted wt mean and wt sd
gb_dat <-gb_dat %>%
  group_by(study_phenotype_3) %>%
  mutate(wt_sum = ifelse(length(weighted_phenotype[zygosity == "wildtype"]) > 0, 
                        (sum(weighted_phenotype[zygosity == "wildtype"])), NA)) %>%
  mutate(tot_wt_SS = ifelse(length(phenotype.measure[zygosity== "wildtype"]) > 0, 
                        sum(Sample.Size[zygosity == "wildtype"]), NA)) %>%
   mutate(wt_SS = ifelse(length(phenotype.measure[zygosity== "wildtype"]) > 0, 
                         mean(Sample.Size[zygosity == "wildtype"]), NA)) %>%
 mutate(wt_mean = wt_sum/tot_wt_SS) %>%
   mutate(wt_sd = ifelse(length(phenotype.measure[zygosity == "wildtype"]) > 0, 
                             sd(phenotype.measure[zygosity == "wildtype"]), NA)) %>%
  select(-wt_sum, -tot_wt_SS) %>%
  ungroup()

# Calculate weighted mutant mean and sd for each study (averaged for all linetypes here)
gb_dat <- gb_dat %>%
  group_by(study_phenotype_4) %>%
  mutate(mutant_sum = ifelse(length(weighted_phenotype[zygosity != "wildtype"]) > 0, 
                         (sum(weighted_phenotype[zygosity != "wildtype"])), NA)) %>%
  mutate(tot_mutant_SS = ifelse(length(weighted_phenotype[zygosity != "wildtype"]) > 0, 
                            (sum(Sample.Size[zygosity != "wildtype"])), NA)) %>%
  mutate(mutant_SS = ifelse(length(phenotype.measure[zygosity != "wildtype"]) > 0, 
                            mean(Sample.Size[zygosity != "wildtype"]), NA)) %>%
  mutate(mutant_mean = mutant_sum/tot_mutant_SS) %>%
  mutate(mutant_sd = ifelse(length(phenotype.measure[zygosity != "wildtype"]) > 0, 
                        sd(phenotype.measure[zygosity != "wildtype"]), NA)) %>%
  select(-mutant_sum, - tot_mutant_SS) %>%
  ungroup() 

dim(mouse_subset)


# Calculate logRR ratio, keeping distinct values and removing strain column
 gb_dat_RR <- gb_dat %>%
 mutate(logRR = log2(mutant_mean/wt_mean)) %>%
  ungroup() %>%
  filter(logRR != "Inf" & logRR != "-Inf" & logRR != "NaN") %>%
   distinct(logRR, .keep_all = TRUE) %>%
   select(-strain)
 
 # Taking the abs value because we are relatively more interested in the magnitude   
gb_dat_RR$logRR <- with(gb_dat_RR, abs(logRR))
 
 # Replacing zeroes to fit a Gamma distribution 
gb_dat_RR$logRR <- with(gb_dat_RR, ifelse(logRR == 0, 0.00000000000000001, logRR)) 

# Number of studies in mouse_subset_RR 
studies_RR <- gb_dat_RR %>%
  group_by(study) %>%
  tally()

dim(studies_RR)

dim(gb_dat_RR)

# logVR ratio
gb_dat_VR <- gb_dat %>%
# With the average sample sizes
  mutate(logVR = log2(mutant_sd/wt_sd) + 1/(2*(mutant_SS*sqrt(strain.count) - 1))
         - 1/(2*(wt_SS*sqrt(strain.count) - 1))) %>%
  filter(logVR != "Inf" & logVR != "-Inf" & logVR != "NaN") %>%
  distinct(logVR, .keep_all = TRUE) %>%
  select(-strain)
  
  
dim(gb_dat_VR)

# Taking the abs value because we are relatively more interested in the magnitude   
gb_dat_VR$logVR <- with(gb_dat_VR, abs(logVR))

# Replacing zeroes to fit a Gamma distribution 
gb_dat_VR$logVR <- with(gb_dat_VR, ifelse(logVR == 0, 0.00000000000000001, logVR)) 

# Number of studies in mouse_subset_VR 
studies_VR <- gb_dat_VR %>%
  group_by(study) %>%
  tally()

dim(studies_VR)

# logCVR ratio
gb_dat_CVR <- gb_dat %>%
  # With the average sample sizes
  mutate(logCVR = log2((mutant_sd/mutant_mean)/(wt_sd/wt_mean)) + 1/(2*(mutant_SS*sqrt(strain.count) - 1))
         - 1/(2*(wt_SS*sqrt(strain.count) - 1))) %>%
  filter(logCVR != "Inf" & logCVR != "-Inf" & logCVR != "NaN") %>%
  distinct(logCVR, .keep_all = TRUE) %>%
  select(-strain)

# Taking the abs value because we are relatively more interested in the magnitude   
gb_dat_CVR$logCVR <- with(gb_dat_CVR, abs(logCVR))

# Replacing zeroes to fit a Gamma distribution 
gb_dat_CVR$logCVR <- with(gb_dat_CVR, ifelse(logCVR == 0, 0.00000000000000001, logCVR)) 


# Number of studies in mouse_subset_CVR 
studies_CVR <- gb_dat_CVR %>%
  group_by(study) %>%
  tally()

dim(studies_CVR)

# Weights 

# Sampling variance of the logRR ratio
gb_dat_RR <- gb_dat_RR %>%
  mutate(logRR_var = (wt_sd)^2/((wt_SS*sqrt(strain.count))*(wt_mean)^2) + 
           (mutant_sd)^2/((mutant_SS*sqrt(strain.count))*(mutant_mean)^2))

# For the logVR ratio
gb_dat_VR <- gb_dat_VR %>%
  mutate(logVR_var = (1/(2*(mutant_SS*sqrt(strain.count) - 1)))
         + (1/(2*(wt_SS*sqrt(strain.count)  - 1))))

# For the logCVR ratio 


# Correlation between wildtype mean and wildtype SD
gb_wt_cor <- cor(gb_dat$phenotype.measure[mouse_subset$zygosity == "wildtype"], 
              gb_dat$SD[mouse_subset$zygosity == "wildtype"],
              use = "complete.obs")

# Correlation between mutant mean and mutant SD
gb_mutant_cor <- cor(gb_dat$phenotype.measure[mouse_subset$zygosity != "wildtype"], 
                  gb_dat$SD[mouse_subset$zygosity != "wildtype"],
                  use = "complete.obs")

mouse_subset_CVR <- mouse_subset_CVR %>%
  mutate(logCVR_var = (wt_sd)^2/((wt_SS*sqrt(strain.count))*(wt_mean)^2) + 1/(2*((mutant_SS*sqrt(strain.count)) - 1)) + 
           (mutant_sd)^2/((mutant_SS*sqrt(strain.count))*(mutant_mean)^2) + 1/(2*((wt_SS*sqrt(strain.count)) - 1)) - 
           2*wt_cor*sqrt(((wt_sd)^2/((wt_SS*sqrt(strain.count))*(wt_mean)^2))*(1/(2*((wt_SS*sqrt(strain.count)) - 1)))) -
           2*mutant_cor*sqrt(((mutant_sd)^2/((mutant_SS*sqrt(strain.count))*(mutant_mean)^2))*(1/(2*((mutant_SS*sqrt(strain.count)) - 1)))))

# Some distribution plots
logRR_density_gbe <- density(gb_dat_RR$logRR) 
plot(logRR_density_gbe) 

logVR_density_gbe <- density(gb_dat_VR$logVR) 
plot(logVR_density_gbe) 

logCVR_density_gbe <- density(gb_dat_CVR$logCVR) 
plot(logCVR_density_gbe) 


################################# Exploratory Plots ########################################
# Mouse subset 
# logRR
ggplot(mouse_subset_RR, aes(x = phenotype.type, y = logRR)) + geom_point(aes(colour = alleletype, shape = zygosity)) + 
  


 ############################# Imputing the weights ############################################
 
########################### Mouse subset ##################################
# logRR ratio 

# Specifying priors
priors <- c(prior(normal(0,5), class = Intercept),
            prior(cauchy(1,0.7), class = sd))

# Replacing 0 values to fit Gamma 
mouse_subset_RR$logRR_var <- with(mouse_subset_RR, 
                                  ifelse(logRR_var == 0, 0.00000000000000001, logRR_var)) 

# logRR weights 
weights_RR_fit1 <- brm(data = mouse_subset_RR, family = Gamma(link = "identity"),
            logRR_var ~ 1 + SD + phenotype.measure + (1|study) + (1|study:phenotype)
            , iter = 10000, prior = priors, warmup = 1000, cores = 2, chains = 1,
            seed = 111, inits  = 0) 

summary(weights_RR_fit1)

# pp_check
pp_check(weights_RR_fit1)

# Predicted values
predict(weights_RR_fit1, newdata =)

head(fitted_values)

## plot expected predictions against actual response
exp_vs_actual <- as.data.frame(cbind(Y = standata(weights_RR_fit1)$Y, fitted_values))
ggplot(exp_vs_actual) + geom_point(aes(x = Estimate, y = Y))

########################### Model fitting for mouse subset ########################################
# check which parameters can have priors
get_prior(logRR|weights(sqrt(logRR_var))  ~ 1 + alleletype + zygosity + (1|study) + (1|study:phenotype)
          + (1|study:sex) + (1|study:environment),
          data = mouse_subset_RR, Gamma(link = "log"))

# Specifying priors
priors <- c(prior(normal(0,5), class = Intercept),
            prior(cauchy(1,0.7), class = sd), prior(normal(0,5), class = b))

# Prior predictive checks
prior_logRR_fit <- brm(data = mouse_subset_RR, family = Gamma(link = "log"),
    logRR|weights(sqrt(logRR_var))  ~ 1 + strain + zygosity
    + alleletype + phenotype.type + alleletype:zygosity + (1|study) + (1|study:phenotype) + (1|study:phenotype.type)
    ,prior = priors, sample_prior='only') 

plot(prior_logRR_fit) 

pp_check(prior_logRR_fit) + xlim(0,10)

ppc_dens_overlay(y = mouse_subset_RR$logRR, 
                 yrep = posterior_predict(prior_logRR_fit, draws=25))

# Using weights specified by the weights argument 
logRR_fit1 <- brm(data = mouse_subset_RR, family = Gamma(link = "log"),
            logRR| weights(sqrt(logRR_var))  ~ 1 + strain + zygosity
            + alleletype + phenotype.type + alleletype:zygosity + (1|study) 
             ,iter = 25000, prior = priors, warmup = 1000, cores = 2, chains = 1,
            seed = 111, inits  = 'random') 

summary(logRR_fit1)

# Prior samples 
prior_samples(logRR_fit1)

# Posterior predictive check
pp_check(logRR_fit1)

# y and yrep for more plots
y <- mouse_subset_RR$logRR
yrep <- posterior_predict(logRR_fit1, draws = 1000)

#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  
ppc_stat(y = y, yrep = yrep,stat = "skew")

# observed and simulated mean and median
ppc_stat(y = y, yrep = yrep,stat = "mean") 
ppc_stat(y = y, yrep = yrep,stat = "median") 

# Comparing PIT values to uniform distribution 
model_loo <- loo(logRR_fit1, save_psis = TRUE, cores=5)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
ppc_loo_pit_overlay(y, yrep, lw = w)

# Pareto k-values disgnostic plot
plot(model_loo)
# Between 0.7-1 influential points

# Prettier plot 
k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Specifying a model for the between study-level variance as well 
logRR_fit2 <- brm(data = mouse_subset_RR, family = Gamma(link = "inverse"),
                  bf(logRR  ~ 1 + strain + zygosity
                  + alleletype + phenotype.type + (1|study) + (1|study:phenotype) + (1|study:phenotype.type),
                  shape ~ 1 + SD + (1|study:phenotype)), iter = 30000, prior = priors, warmup = 1000,control = list(adapt_delta = 0.95), cores = 5, chains = 2,
                  seed = 111, inits  = 'random') 

# Prior samples 
prior_samples(logRR_fit2)

# Posterior predictive check
pp_check(logRR_fit2)

# y and yrep for more plots
y <- mouse_subset_RR$logRR
yrep <- posterior_predict(logRR_fit2, draws = 1000)

#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  
ppc_stat(y = y, yrep = yrep,stat = "skew")

# observed and simulated mean and median
ppc_stat(y = y, yrep = yrep,stat = "mean") 
ppc_stat(y = y, yrep = yrep,stat = "median") 

# Comparing PIT values to uniform distribution 
model_loo <- loo(logRR_fit2, save_psis = TRUE, cores=5)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
ppc_loo_pit_overlay(y, yrep, lw = w)

# Pareto k-values disgnostic plot
plot(model_loo)
# Between 0.7-1 influential points

# Prettier plot 
k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Comparing the two models 
loo(fit1, fit2)

#  For better comparison 
kfold(fit1, fit2, compare = TRUE, K = 10)

# Posterior plots
plot(fit1, pars = c("alleletype","line.type","phenotype.type")) 

# Effects plots
plot(conditional_effects(logRR_fit1, effects = "phenotype.type:alleletype"))

plot(conditional_effects(logRR_fit2, effects = "phenotype.type:alleletype"))

######################################  Log VR Ratio #######################################
# Specifying priors
priors <- c(prior(normal(0,5), class = Intercept),
            prior(cauchy(0,2), class = sd), prior(student_t(3,0,2.5), class = b))

prior_logRR_fit <- brm(data = mouse_subset_VR, family = Gamma(link = "log"),
                       logVR|weights(sqrt(logVR_var))  ~ 1 + strain + zygosity
                       + alleletype + phenotype.type + alleletype:zygosity + (1|study) + (1|study:phenotype) + (1|study:phenotype.type)
                       ,prior = priors, sample_prior='only') 

pp_check(prior_logRR_fit) + xlim(0,5)

shape ~ 1 + offset(log(xdisp))


logVR_fit1 <- brm(data = mouse_subset_VR, family = Gamma(link = "inverse"),
                  logVR|weights(sqrt(logVR_var))  ~ 1 + strain + zygosity
                  + alleletype + phenotype.type + alleletype:zygosity + (1|study) + (1|study:phenotype) + (1|study:phenotype.type)
                  ,iter = 25000, prior = priors, warmup = 1000,control = list(adapt_delta = 0.97), cores = 2, chains = 1,
                  seed = 111, inits  = 'random') 

summary(logVR_fit1)

# Prior samples 
prior_samples(logRR_fit2)

# Posterior predictive check
pp_check(logRR_fit2)

# y and yrep for more plots
y <- mouse_subset_RR$logRR
yrep <- posterior_predict(logVR_fit2, draws = 1000)

#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  
ppc_stat(y = y, yrep = yrep,stat = "skew")

# observed and simulated mean and median
ppc_stat(y = y, yrep = yrep,stat = "mean") 
ppc_stat(y = y, yrep = yrep,stat = "median") 

# Comparing PIT values to uniform distribution 
model_loo <- loo(logVR_fit2, save_psis = TRUE, cores=5)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
ppc_loo_pit_overlay(y, yrep, lw = w)

# Pareto k-values disgnostic plot
plot(model_loo)
# Between 0.7-1 influential points

# Prettier plot 
k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

logVR_fit2 <- brm(data = mouse_subset_VR, family = Gamma(link = "inverse"),
                  bf(logVR  ~ 1 + strain + zygosity
                     + alleletype + phenotype.type + (1|study) + (1|study:phenotype) + (1|study:phenotype.type),
                     sigma ~ 1 + SD + (1|study:phenotype)), iter = 30000, 
                  prior = priors, warmup = 1000,control = list(adapt_delta = 0.95), cores = 5, chains = 2,
                  seed = 111, inits  = 'random') 

# Prior samples 
prior_samples(logRR_fit2)

# Posterior predictive check
pp_check(logRR_fit2)

# y and yrep for more plots
y <- mouse_subset_RR$logRR
yrep <- posterior_predict(logRR_fit2, draws = 1000)

#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  
ppc_stat(y = y, yrep = yrep,stat = "skew")

# observed and simulated mean and median
ppc_stat(y = y, yrep = yrep,stat = "mean") 
ppc_stat(y = y, yrep = yrep,stat = "median") 

# Comparing PIT values to uniform distribution 
model_loo <- loo(logRR_fit2, save_psis = TRUE, cores=5)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
ppc_loo_pit_overlay(y, yrep, lw = w)

# Pareto k-values disgnostic plot
plot(model_loo)
# Between 0.7-1 influential points

# Prettier plot 
k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Comparing the two models 
loo(fit1, fit2)

#  For better comparison 
kfold(fit1, fit2, compare = TRUE, K = 10)

# Posterior plots
plot(fit1, pars = c("alleletype","line.type","phenotype.type")) 

# Effects plots
plot(conditional_effects(logRR_fit1, effects = "phenotype.type:alleletype"))

plot(conditional_effects(logRR_fit2, effects = "phenotype.type:alleletype"))

################################## logCVR ########################################


#############################  For the entire data-set ############################


#############################  logRR ##################################3
# check which parameters can have priors
get_prior(logRR|weights(sqrt(logRR_var))  ~ 1 + alleletype + zygosity + (1|study) + (1|study:phenotype)
          + (1|study:sex) + (1|study:environment),
          data = mouse_subset_RR, Gamma(link = "inverse"))

# Specifying priors
priors <- c(prior(normal(0,5), class = Intercept),
            prior(cauchy(1,0.7), class = sd), prior(student_t(5,0,10), class = b))

# Using weights specified by the weights argument (where sigma is also estimated)
logRR_fit3 <- brm(data = gb_dat_RR, family = Gamma(link = "inverse"),
                  logRR|weights(sqrt(logRR_var))  ~ 1 + strain + zygosity
                  + alleletype + phenotype.type + alleletype:zygosity + (1|study) + (1|study:phenotype) + (1|study:phenotype.type) + (1|study:species)
                  ,iter = 25000, prior = priors, warmup = 1000,control = list(adapt_delta = 0.95), cores = 5, chains = 2,
                  seed = 111, inits  = 'random') 

summary(logRR_fit3)

# Prior samples 
prior_samples(logRR_fit3)

# Posterior predictive check
pp_check(logRR_fit3)

# y and yrep for more plots
y <- mouse_subset_RR$logRR
yrep <- posterior_predict(logRR_fit3, draws = 1000)

#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  
ppc_stat(y = y, yrep = yrep,stat = "skew")

# observed and simulated mean and median
ppc_stat(y = y, yrep = yrep,stat = "mean") 
ppc_stat(y = y, yrep = yrep,stat = "median") 

# Comparing PIT values to uniform distribution 
model_loo <- loo(logRR_fit3, save_psis = TRUE, cores=5)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
ppc_loo_pit_overlay(y, yrep, lw = w)

# Pareto k-values disgnostic plot
plot(model_loo)
# Between 0.7-1 influential points

# Prettier plot 
k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Specifying a model for the between study-level variance as well 
logRR_fit4 <- brm(data = gb_dat_RR, family = Gamma(link = "inverse"),
                  bf(logRR  ~ 1 + strain + zygosity
                     + alleletype + phenotype.type + (1|study) + (1|study:phenotype) + (1|study:phenotype.type) + (1|study:species),
                     sigma ~ 1 + SD + (1|study:phenotype)), iter = 30000, prior = priors, warmup = 1000,control = list(adapt_delta = 0.95), cores = 5, chains = 2,
                  seed = 111, inits  = 'random') 

# Prior samples 
prior_samples(logRR_fit4)

# Posterior predictive check
pp_check(logRR_fit4)

# y and yrep for more plots
y <- mouse_subset_RR$logRR
yrep <- posterior_predict(logRR_fit4, draws = 1000)

#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  
ppc_stat(y = y, yrep = yrep,stat = "skew")

# observed and simulated mean and median
ppc_stat(y = y, yrep = yrep,stat = "mean") 
ppc_stat(y = y, yrep = yrep,stat = "median") 

# Comparing PIT values to uniform distribution 
model_loo <- loo(logRR_fit4, save_psis = TRUE, cores=5)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
ppc_loo_pit_overlay(y, yrep, lw = w)

# Pareto k-values disgnostic plot
plot(model_loo)
# Between 0.7-1 influential points

# Prettier plot 
k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Comparing the two models 
loo(fit1, fit2)

#  For better comparison 
kfold(fit1, fit2, compare = TRUE, K = 10)

# Posterior plots
plot(fit1, pars = c("alleletype","line.type","phenotype.type")) 

# Effects plots
plot(conditional_effects(logRR_fit1, effects = "phenotype.type:alleletype"))

plot(conditional_effects(logRR_fit2, effects = "phenotype.type:alleletype"))

######################################  Log VR Ratio #######################################
# Specifying priors
priors <- c(prior(normal(0,5), class = Intercept),
            prior(cauchy(1,0.7), class = sd), prior(student_t(5,0,10), class = b))

# logVR ratio
logVR_fit3 <- brm(data = mouse_subset_VR, family = Gamma(link = "inverse"),
                  logRR|weights(sqrt(logRR_var))  ~ 1 + strain + zygosity
                  + alleletype + phenotype.type + alleletype:zygosity + (1|study) + (1|study:phenotype) + (1|study:phenotype.type) + (1|study:species)
                  ,iter = 25000, prior = priors, warmup = 1000,control = list(adapt_delta = 0.97), cores = 5, chains = 2,
                  seed = 111, inits  = 'random') 

summary(logVR_fit3)

# Prior samples 
prior_samples(logRR_fit2)

# Posterior predictive check
pp_check(logRR_fit2)

# y and yrep for more plots
y <- mouse_subset_RR$logRR
yrep <- posterior_predict(logVR_fit2, draws = 1000)

#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  
ppc_stat(y = y, yrep = yrep,stat = "skew")

# observed and simulated mean and median
ppc_stat(y = y, yrep = yrep,stat = "mean") 
ppc_stat(y = y, yrep = yrep,stat = "median") 

# Comparing PIT values to uniform distribution 
model_loo <- loo(logVR_fit2, save_psis = TRUE, cores=5)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
ppc_loo_pit_overlay(y, yrep, lw = w)

# Pareto k-values disgnostic plot
plot(model_loo)
# Between 0.7-1 influential points

# Prettier plot 
k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

logVR_fit2 <- brm(data = mouse_subset_VR, family = Gamma(link = "inverse"),
                  bf(logVR  ~ 1 + strain + zygosity
                     + alleletype + phenotype.type + (1|study) + (1|study:phenotype) + (1|study:phenotype.type),
                     sigma ~ 1 + SD + (1|study:phenotype)), iter = 30000, 
                  prior = priors, warmup = 1000,control = list(adapt_delta = 0.95), cores = 5, chains = 2,
                  seed = 111, inits  = 'random') 

logVR_fit3 <- brm(data = mouse_subset_VR, family = Gamma(link = "inverse"),
                  bf(logVR  ~ 1 + strain + zygosity
                     + alleletype + phenotype.type + (1|study) + (1|study:phenotype) + (1|study:phenotype.type),
                     sigma ~ 1 + sqrt(logRR_var)), iter = 30000, 
                  prior = priors, warmup = 1000,control = list(adapt_delta = 0.95), cores = 5, chains = 2,
                  seed = 111, inits  = 'random') 

# Prior samples 
prior_samples(logRR_fit2)

# Posterior predictive check
pp_check(logRR_fit2)

# y and yrep for more plots
y <- mouse_subset_RR$logRR
yrep <- posterior_predict(logRR_fit2, draws = 1000)

#  Fisher-Pearson Skew function based on NIST definition
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

color_scheme_set("gray")

# observed and simulated SKEW metric  
ppc_stat(y = y, yrep = yrep,stat = "skew")

# observed and simulated mean and median
ppc_stat(y = y, yrep = yrep,stat = "mean") 
ppc_stat(y = y, yrep = yrep,stat = "median") 

# Comparing PIT values to uniform distribution 
model_loo <- loo(logRR_fit2, save_psis = TRUE, cores=5)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
ppc_loo_pit_overlay(y, yrep, lw = w)

# Pareto k-values disgnostic plot
plot(model_loo)
# Between 0.7-1 influential points

# Prettier plot 
k_rintercept <- model_loo$psis_object$diagnostics$pareto_k

df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)

ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

# Comparing the two models 
loo(fit1, fit2)

#  For better comparison 
kfold(fit1, fit2, compare = TRUE, K = 10)

# Posterior plots
plot(fit1, pars = c("alleletype","line.type","phenotype.type")) 

# Effects plots
plot(conditional_effects(logRR_fit1, effects = "phenotype.type:alleletype"))

plot(conditional_effects(logRR_fit2, effects = "phenotype.type:alleletype"))

################################## logCVR ########################################





#################################################################################
# log CVR assumes a linear relationship, checking that 

# Wildtype strains

png(file="mouse_subset_logwt_meanvssd.png",width=20,height=15,units="cm",res=300)

ggplot(data = mouse_subset, aes(x = log(wt_mean), y = log(wt_sd))) +
  geom_point() + labs( x = "Log Wildtype Mean", y = "Log Wildtype Standard Deviation") +
  theme_classic() + ggtitle("Logarithmic mean-variance relationship of wildtype mouse values")

dev.off()

cor(gb_dat$wt_mean, gb_dat$wt_sd, use="complete.obs")

# Mutant strains

png(file="mouse_subset_mutant_meanvssd.png",width=20,height=15,units="cm",res=300)

ggplot(data = mouse_subset, aes(x = log(mutant_mean), y = log(mutant_sd))) +
  geom_point() + labs( x = "Log Mutant Mean", y = "Log Mutant Standard Deviation") +
  theme_classic() + ggtitle("Logarithmic mean-variance relationship of mutant mouse values")

dev.off()

cor(gb_dat$mutant_mean, gb_dat$mutant_sd, use="complete.obs")


# LogRR vs logSD ratios

png(file="mouse_subset_logRRvslogSD.png",width=20,height=15,units="cm",res=300)

ggplot(data = mouse_subset, aes(x = logRR, y = logSD)) +
  geom_point() + labs( x = "Log RR", y = "Log SD") +
  theme_classic()

dev.off()

# Some informative boxplots 

# With the logRR ratio 

# With alleletype 
png(file="mouse_subset_logRR_alleletype.png",width=20,height=15,units="cm",res=300)

ggplot(mouse_subset, aes(alleletype, logRR)) + geom_boxplot() + theme_classic()

dev.off()

# With line.type
png(file="mouse_subset_logRR_linetype.png",width=20,height=15,units="cm",res=300)

ggplot(mouse_subset, aes(line.type, logRR)) + geom_boxplot() + theme_classic()

dev.off()

# With phenotype.type 
png(file="mouse_subset_logRR_phenotype.png",width=20,height=15,units="cm",res=300)

ggplot(mouse_subset, aes(phenotype.type, logRR)) + geom_boxplot() + theme_classic()

dev.off()

# With the logSD ratio 

# With alleletype 

png(file="logSD_alleletype.png",width=20,height=15,units="cm",res=300)

ggplot(gb_dat, aes(alleletype, logSD)) + geom_boxplot() + theme_classic()

dev.off()

# With line.type
png(file="mouse_subset_logSD_linetype.png",width=20,height=15,units="cm",res=300)

ggplot(mouse_subset, aes(line.type, logSD)) + geom_boxplot() + theme_classic()

dev.off()

# With phenotype.type 
png(file="mouse_subset_logSD_alleletype.png",width=20,height=15,units="cm",res=300)

ggplot(mouse_subset, aes(phenotype.type, logSD)) + geom_boxplot() + theme_classic()

dev.off()

# With the logcCVR ratio 

# With alleletype 
png(file="mouse_subset_logCVR_alleletype.png",width=20,height=15,units="cm",res=300)

ggplot(mouse_subset, aes(alleletype, logCVR)) + geom_boxplot() + theme_classic()

dev.off()

# With line.type
png(file="mouse_subset_logCVR_linetype.png",width=20,height=15,units="cm",res=300)

ggplot(gb_dat, aes(line.type, logCVR)) + geom_boxplot() + theme_classic()

dev.off()

# With phenotype.type 
png(file="mouse_subset_logCVR_phenotype.png",width=20,height=15,units="cm",res=300)

ggplot(gb_dat, aes(phenotype.type, logCVR)) + geom_boxplot() + theme_classic()

dev.off()

# Trying code to make forest plots 
study.draws <- spread_draws(fit1, r_study[study,], b_Intercept) %>% 
  mutate(b_Intercept = r_study[,Intercept] + b_Intercept)

pooled.effect.draws <- spread_draws(fit1, b_Intercept)
