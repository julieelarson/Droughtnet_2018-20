#' ---
#' title: "Analyses for Rainfall X Grazing Experiment, Boulder CO"
#' author: "Julie Larson"
#' date: "6 August 2025"
#' output:  html_document
#' ---
#' 

#' 
#' **Welcome!**
#'
#' This R script contains code to conduct analyses and create figures included 
#' in the following manuscript in preparation:
#' 
#' Larson, J., B. Anacker, T. Merchant, K. Suding. (In Prep.) Experimental 
#' grazing and rainfall treatments reveal key contingencies underlying 
#' rangeland resistance
#' 
#' 
#' In this analysis, data are analyzed from the first three years (2018-2020) 
#' of an ongoing field experiment manipulating rainfall (dry, wet, control) 
#' and grazing (growing season, dormant season, ungrazed) in a grassland (Boulder, CO, USA).
#' 
#' This ongoing project is conducted in collaboration with the *University of Colorado* 
#' and *City of Boulder Open Space and Mountain Parks (Boulder, CO, USA)*.
#' 



#' 
#' **Load packages, functions, and options**
#' 

# Packages
#+ results=FALSE, message=FALSE, warning=FALSE
library(tidyverse)
library(gridExtra)
library(vegan)
library(lme4)
library(car)
library(lmerTest)
library(performance)
library(scales)
library(codyn)
library(GGally)
library(patchwork)
library(ggcorrplot)
library(magick)
library(ggeffects)
library(emmeans)
library(ggExtra)
library(Hmisc)

# Load source code for Oridicenter() function
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:ordicenter?do=export_code&codeblock=0')
# Ordihull() source code
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:orglhull?do=export_code&codeblock=1')

# Load function to calculate standard errors
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
} 

# Set option for contrasts (Type III sums of squares (lm's))
options(contrasts = c("contr.sum","contr.poly"))

# Set option for tibble printing
options(tibble.print_max = 50, tibble.print_min = 50)




#' 
#' **Load data**
#' 

#' Set user's path
gdrive <- "C:/Users/larsonju/OneDrive - UW/Documents/R/Droughtnet_git/Published_Files_Zenodo_2025/"

#' Soil volumetric water content (plot-level data) though experiment
vwc_path <-   paste0(gdrive, "Soil_VWC_2018-20.csv")
vwc <- read.csv(vwc_path)


#' Plant data (biomass, annuals, and bare ground data)
plant_path <-  paste0(gdrive, "Plant_Plot_Metrics_2018-23.csv")
plant <- read.csv(plant_path)

#' Forage Quality data 
forage_path <-  paste0(gdrive, "ForageQuality_2020.csv")
forage <- read.csv(forage_path)

#' Plant Cover data
cov_path <- paste0(gdrive, "Plant_Plot_Cover_AreaClipped_2018-20.csv")
cov <- read.csv(cov_path)

#' Phenology data
pheno_path <- paste0(gdrive, "PhenologyGreennessExtraction_2018-20.csv")
pheno <- read.csv(pheno_path)

pheno_mets_path <- paste0(gdrive, "GreennessCurveMetrics_compiled.csv")
pheno_mets <- read.csv(pheno_mets_path)








#####
#'
#'
#'  **Soil Moisture**
#' 
#' 



#' *Overview:*
#' We sampled soil volumetric water content periodically in half of the 
#' plots (n=36). Sampling occurred at a permanent location in the center of
#' each plot, to a depth of ~13cm.

#' Read in data 
vwc_path <-   paste0(gdrive, "Soil_VWC_2018-20.csv")
vwc_all <- read.csv(vwc_path)

#' View data structure
str(vwc_all)



#' *Prepare data*
#'
#' Filter out observations that were taken when rainfall shelters were not installed
vwc <- vwc_all %>% filter (rain_shelters_on == "Y") 

#' Create a factor capturing unique plots (block-plot combination) 
vwc$b_p <- factor(paste(vwc$block, vwc$plot,sep='-'))
 
#' Set rainfall, grazing, and year variables as factors
vwc$rain_treat <- factor(vwc$rain_treat, levels=c("D","C","W"))
vwc$graze_treat <- factor(vwc$graze_treat, levels=c("UG","GS","DS"))
vwc$year <- factor(vwc$year)



#' *Check data distributions for outliers*
#'
#' In 2018, water treatments should be equal across grazing groups, but 
#' there appears to be a positive skew for dormant-season grazing under dry conditions
#' and growing-season grazing under wet conditions.
#' 
#' Plot *3-5* (wet-growing season) and *4-11* (dry dormant-season) appear to be 
#' driving these skews, sometimes to illogical levels (>40%).

#' Filter outlier data
vwc_outlier_list <- c('3-5', '4-11')
vwc_outliers <- vwc %>% filter(b_p %in% vwc_outlier_list)

#' View all data, with outlying plots shown in red
ggplot() + 
  geom_point(aes(x=rain_treat,y=vwc_mid_FS, color=rain_treat), alpha=0.4, data=vwc)+
  geom_boxplot(aes(x=rain_treat,y=vwc_mid_FS, fill=rain_treat), alpha=0.4,data=vwc)+
  labs(y="Soil VWC\n", x="Grazing Treatment", color="Rainfall Treatment", fill="Rainfall Treatment")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  geom_point(aes(x=rain_treat,y=vwc_mid_FS, color=rain_treat), alpha=0.8, color='red', data=vwc_outliers)+
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_treat,scale='free')

#' With just two plots driving some misleading trends about grazing effects in 2018,
#' the choice moving forward is to remove these two plots from VWC analysis.
#' 
#' Remove outlying plots prior to analysis:
vwc_out <- vwc %>% filter(b_p != '3-5') %>% filter(b_p != '4-11')



#' *Prepare for analysis*
#' 
#' Create several model dataframes to view different sets of contrasts 
#' 
#' Set 1
vwc_mod <- vwc_out
vwc_mod$rain_treat <- factor(vwc_mod$rain_treat, levels=c("D","W","C"))
vwc_mod$graze_treat <- factor(vwc_mod$graze_treat, levels=c("DS","GS","UG"))
vwc_mod$year <- factor(vwc_mod$year, levels=c("2020","2019","2018"))

#' Set 2 
vwc_mod_2 <- vwc_mod
vwc_mod_2$year <- factor(vwc_mod_2$year, levels=c("2018","2019","2020"))
vwc_mod_2$rain_treat <- factor(vwc_mod_2$rain_treat, levels=c("C","D","W"))
vwc_mod_2$graze_treat <- factor(vwc_mod_2$graze_treat, levels=c("UG","DS","GS"))

#' Set 3
vwc_mod_3 <- vwc_mod
vwc_mod_3$year <- factor(vwc_mod_3$year, levels=c("2018","2019","2020"))
vwc_mod_3$rain_treat <- factor(vwc_mod_3$rain_treat, levels=c("W","C","D"))
vwc_mod_3$graze_treat <- factor(vwc_mod_3$graze_treat, levels=c("UG","GS","DS"))

#' Set 4
vwc_mod_4 <- vwc_mod
vwc_mod_4$year <- factor(vwc_mod_4$year, levels=c("2020","2019","2018"))
vwc_mod_4$rain_treat <- factor(vwc_mod_4$rain_treat, levels=c("C","D","W"))
vwc_mod_4$graze_treat <- factor(vwc_mod_4$graze_treat, levels=c("UG","GS","DS"))



#' *Analysis*
#'

#' Specify model 
vwc_mod_all <- lmer( vwc_mid_FS ~ rain_treat * graze_treat * year + (1|b_p) + (1|day_of_year), data=vwc_mod)

#' **TABLE S3 -- Model Results**
anova(vwc_mod_all)
summary(vwc_mod_all)

#' Check model diagnostics
#' check for *singularity*
tt <- getME(vwc_mod_all ,"theta")
ll <- getME(vwc_mod_all ,"lower")
min(tt[ll==0])
#' check for *R2* 
lmm_vwc_mod_all_r2 <- data.frame(unlist(r2(vwc_mod_all)))
lmm_vwc_mod_all_r2

#' Get other model contrasts
vwc_mod_all2 <- lmer( vwc_mid_FS ~ rain_treat * graze_treat * year + (1|b_p) + (1|day_of_year), data=vwc_mod_2)
vwc_mod_all3 <- lmer( vwc_mid_FS ~ rain_treat * graze_treat * year + (1|b_p) + (1|day_of_year), data=vwc_mod_3)
vwc_mod_all4 <- lmer( vwc_mid_FS ~ rain_treat * graze_treat * year + (1|b_p) + (1|day_of_year), data=vwc_mod_4)

#' Save univariate fixed effects from all models 
vwc_fix <- data.frame(round(fixef(vwc_mod_all),4)) %>% rename(fixef=round.fixef.vwc_mod_all...4.)
vwc_ci  <- data.frame(round(confint(vwc_mod_all),4))
vwc_fixef<- merge(vwc_fix, vwc_ci, by=0, all=TRUE) 
vwc_fixef$response <- "vwc"
vwc_fixef$Row.names <- gsub("rain_treat1", "rain_D", vwc_fixef$Row.names)
vwc_fixef$Row.names <- gsub("rain_treat2", "rain_W", vwc_fixef$Row.names)
vwc_fixef$Row.names <- gsub("graze_treat1", "graze_DS", vwc_fixef$Row.names)
vwc_fixef$Row.names <- gsub("graze_treat2", "graze_GS", vwc_fixef$Row.names)
vwc_fixef$Row.names <- gsub("year1", "year_2020", vwc_fixef$Row.names)
vwc_fixef$Row.names <- gsub("year2", "year_2019", vwc_fixef$Row.names)
vwc_fix2 <- data.frame(round(fixef(vwc_mod_all2),4))%>% rename(fixef=round.fixef.vwc_mod_all2...4.)
vwc_ci2  <- data.frame(round(confint(vwc_mod_all2),4))
vwc_fixef2<- merge(vwc_fix2, vwc_ci2, by=0, all=TRUE) 
vwc_fixef2$response <- "vwc"
vwc_fixef2$Row.names <- gsub("rain_treat1", "rain_C", vwc_fixef2$Row.names)
vwc_fixef2$Row.names <- gsub("rain_treat2", "rain_D", vwc_fixef2$Row.names)
vwc_fixef2$Row.names <- gsub("graze_treat1", "graze_UG", vwc_fixef2$Row.names)
vwc_fixef2$Row.names <- gsub("graze_treat2", "graze_DS", vwc_fixef2$Row.names)
vwc_fixef2$Row.names <- gsub("year1", "year_2018", vwc_fixef2$Row.names)
vwc_fixef2$Row.names <- gsub("year2", "year_2019", vwc_fixef2$Row.names)
vwc_fix3 <- data.frame(round(fixef(vwc_mod_all3),4))%>% rename(fixef=round.fixef.vwc_mod_all3...4.)
vwc_ci3  <- data.frame(round(confint(vwc_mod_all3),4))
vwc_fixef3<- merge(vwc_fix3, vwc_ci3, by=0, all=TRUE) 
vwc_fixef3$response <- "vwc"
vwc_fixef3$Row.names <- gsub("rain_treat1", "rain_W", vwc_fixef3$Row.names)
vwc_fixef3$Row.names <- gsub("rain_treat2", "rain_C", vwc_fixef3$Row.names)
vwc_fixef3$Row.names <- gsub("graze_treat1", "graze_UG", vwc_fixef3$Row.names)
vwc_fixef3$Row.names <- gsub("graze_treat2", "graze_GS", vwc_fixef3$Row.names)
vwc_fixef3$Row.names <- gsub("year1", "year_2018", vwc_fixef3$Row.names)
vwc_fixef3$Row.names <- gsub("year2", "year_2019", vwc_fixef3$Row.names)
vwc_fix4 <- data.frame(round(fixef(vwc_mod_all4),4))%>% rename(fixef=round.fixef.vwc_mod_all4...4.)
vwc_ci4  <- data.frame(round(confint(vwc_mod_all4),4))
vwc_fixef4<- merge(vwc_fix4, vwc_ci4, by=0, all=TRUE) 
vwc_fixef4$response <- "vwc"
vwc_fixef4$Row.names <- gsub("rain_treat1", "rain_C", vwc_fixef4$Row.names)
vwc_fixef4$Row.names <- gsub("rain_treat2", "rain_D", vwc_fixef4$Row.names)
vwc_fixef4$Row.names <- gsub("graze_treat1", "graze_UG", vwc_fixef4$Row.names)
vwc_fixef4$Row.names <- gsub("graze_treat2", "graze_GS", vwc_fixef4$Row.names)
vwc_fixef4$Row.names <- gsub("year1", "year_2020", vwc_fixef4$Row.names)
vwc_fixef4$Row.names <- gsub("year2", "year_2019", vwc_fixef4$Row.names)

#' view 95% confidence intervals for all main effects and interactions*
vwc_fixef_all <- unique(rbind(vwc_fixef, vwc_fixef2, vwc_fixef3, vwc_fixef4)) %>% arrange(Row.names)
vwc_fixef_all



#'
#' **FIGURE S3:  VWC over time in each year **
#' 

vwc_fig <- vwc_out

# Recode grazing treatment labels
vwc_fig$graze_treat <- gsub("UG","Ungrazed", vwc_fig$graze_treat )
vwc_fig$graze_treat <- gsub("GS","Growing-Season Grazing", vwc_fig$graze_treat )
vwc_fig$graze_treat <- gsub("DS","Dormant-Season Grazing", vwc_fig$graze_treat )

# Set grazing levels
vwc_fig$graze_treat <- factor(vwc_fig$graze_treat, levels=c('Growing Season','Dormant Season','Ungrazed'))
levels(vwc_fig$graze_treat)

#' Create separate dataframes for different years
vwc_18 <- vwc_fig %>% filter( year == 2018)
vwc_19 <- vwc_fig %>% filter( year == 2019)
vwc_20 <- vwc_fig %>% filter( year == 2020) 


#'
#'2018 figure
#'

#'Calculate Means and SEs for VWC in 2018 (rain treatment only)
vwc_18_means <- vwc_18 %>% 
  select(day_of_year, rain_treat,graze_treat, vwc_mid_FS) %>%
  group_by(day_of_year, rain_treat, graze_treat) %>%
  summarise(vwc_mean = mean(vwc_mid_FS, na.rm=T),
            vwc_se = se(vwc_mid_FS, na.rm=T))
vwc_18_means$rain_treat <- factor(vwc_18_means$rain_treat , levels = c("D","C","W"))
vwc_18_means$graze_treat <- factor(vwc_18_means$graze_treat , levels = c("Growing-Season Grazing","Dormant-Season Grazing","Ungrazed"))

#' Create dataframe of labels to indicate the number of days since watering event
vwc_18_water <- vwc_18 %>% 
  select(day_of_year, day_since_wateradd) %>%
  group_by(day_of_year) %>%
  summarise(days_since_water = max(day_since_wateradd))
vwc_18_water$y <- 0

#' Make figure for 2018
vwc_2018<- ggplot() + 
  geom_point(data = vwc_18_means, aes( x=day_of_year, y=vwc_mean,lty=graze_treat, shape=graze_treat,color=rain_treat), cex=2.5) +
  geom_line(data = vwc_18_means, aes( x=day_of_year, y=vwc_mean, lty=graze_treat, color=rain_treat), lwd=0.8)+
  geom_errorbar(data = vwc_18_means, aes(color=rain_treat, x=day_of_year, ymin=vwc_mean-vwc_se, ymax=vwc_mean+vwc_se), width=0.25)+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  xlim(c(120,270)) + 
  ylim(c(-5,50)) + 
  geom_text(aes(x=125, y=42, label= "2018"), cex=4) +
  geom_text(data=vwc_18_water, aes(x=day_of_year, y=y, label=days_since_water), color="gray40", cex = 2.75)+
  geom_text(aes(x=140, y=30, label="Note that grazing treatments\nwere not yet implemented in 2018"),color='gray40',cex=3)+
  labs(x="Day of Year", y=" Soil volumetric \n water content (%)", color="Rainfall \nTreatment") +
  theme_bw()+
  theme(legend.position="none")
vwc_2018



#'
#'2019
#'

#'Calculate Means and SEs for VWC in 2019 (rainfall and grazing)
vwc_19_means <- vwc_19 %>% 
  select(day_of_year, rain_treat, graze_treat, vwc_mid_FS) %>%
  group_by(day_of_year, rain_treat, graze_treat) %>%
  summarise(vwc_mean = mean(vwc_mid_FS, na.rm=T),
            vwc_se = se(vwc_mid_FS, na.rm=T))
vwc_19_means$rain_treat <- factor(vwc_19_means$rain_treat , levels = c("D","C","W"))
vwc_19_means$graze_treat <- factor(vwc_19_means$graze_treat , levels = c("Growing-Season Grazing","Dormant-Season Grazing","Ungrazed"))

# Create dataframe of labels to indicate number of days since watering event
vwc_19_water <- vwc_19 %>% 
  select(day_of_year, day_since_wateradd) %>%
  group_by(day_of_year) %>%
  summarise(days_since_water = max(day_since_wateradd))
vwc_19_water$y <- 0

vwc_2019 <- ggplot() + 
  geom_point(data = vwc_19_means, aes( x=day_of_year, y=vwc_mean, color=rain_treat, shape=graze_treat), cex=2.5) +
  geom_line(data = vwc_19_means, aes( x=day_of_year, y=vwc_mean, color=rain_treat, lty=graze_treat), lwd=.8)+
  geom_errorbar(data = vwc_19_means, aes(color=rain_treat, x=day_of_year, ymin=vwc_mean-vwc_se, ymax=vwc_mean+vwc_se), width=0.25)+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  xlim(c(120,270)) + 
  ylim(c(-5,50)) + 
  geom_text(aes(x=125, y=42, label= "2019"), cex=4) +
  geom_text(data=vwc_19_water, aes(x=day_of_year, y=y, label=days_since_water), color="gray40", cex = 2.75)+
  labs(x="Day of Year", y="Soil voilumetric \n water conent (%)", color = "Rainfall \nTreatment", shape="Graze \nTreatment", lty="Graze \nTreatment") +
  theme_bw() 
vwc_2019


#'
#'2020
#'

#'Calculate Means and SEs for VWC in 2019 (rainfall and grazing)
vwc_20_means <- vwc_20 %>% 
  select(day_of_year, rain_treat, graze_treat, vwc_mid_FS) %>%
  group_by(day_of_year, rain_treat, graze_treat) %>%
  summarise(vwc_mean = mean(vwc_mid_FS, na.rm=T),
            vwc_se = se(vwc_mid_FS, na.rm=T))
vwc_20_means$rain_treat <- factor(vwc_20_means$rain_treat , levels = c("D","C","W"))
vwc_20_means$graze_treat <- factor(vwc_20_means$graze_treat , levels = c("Growing-Season Grazing","Dormant-Season Grazing","Ungrazed"))

# Create dataframe of labels to indicate number of days since watering event
vwc_20_water <- vwc_20 %>% 
  select(day_of_year, day_since_wateradd) %>%
  group_by(day_of_year) %>%
  summarise(days_since_water = max(day_since_wateradd))
vwc_20_water$y <- c(-2)

vwc_2020 <- ggplot() + 
  geom_point(data = vwc_20_means, aes( x=day_of_year, y=vwc_mean, color=rain_treat, shape=graze_treat), cex=2.5) +
  geom_line(data = vwc_20_means, aes( x=day_of_year, y=vwc_mean, color=rain_treat, lty=graze_treat), lwd=.8)+
  geom_errorbar(data = vwc_20_means, aes(color=rain_treat, x=day_of_year, ymin=vwc_mean-vwc_se, ymax=vwc_mean+vwc_se), width=0.25)+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  xlim(c(120,270)) +
  ylim(c(-5,50)) +
  geom_text(data=vwc_20_water, aes(x=day_of_year, y=y, label=days_since_water), color="gray40", cex = 2.75)+
  geom_text(aes(x=125, y=42, label= "2020"), cex=4) +
  labs(x="Day of Year", y="Soil voilumetric \nwater conent (%)", color = "Rainfall \nTreatment", shape="Graze \nTreatment", lty="Graze \nTreatment") +
  theme_bw() +
  theme(legend.position="none")
vwc_2020


#'  
#'  Compile figures together with labeled fixed effects
#'  
vwc_fig_all_yr_model <- (vwc_2018/ vwc_2019 / vwc_2020) + plot_annotation(
  title = 'Fixed Effects:',
  subtitle = 'Rainfall p<0.001   Year p<0.001   Rain*Year p=0.001 \n')
vwc_fig_all_yr_model






#####
#'
#'
#'   **Plant Community Composition**
#'              
#'  


#'  **Overview**
#'  
#'  This script explores how rainfall and grazing treatments affect community
#'  composition over time, including species composition (NMDS and PerMANOVA),
#'  as well as the relative abundances of common plant functional groups
#'  (linear mixed models) 



#' 
#' **Read in data**
#' 

#' Set user's path
drive <- "C:/Users/larsonju/OneDrive - UW/Documents/R/Droughtnet_git/Published_Files_Zenodo_2025/"

comp_path <-  paste0(gdrive, "Plant_Plot_Cover_AreaClipped_2018-20.csv")
comp_raw <- read.csv(comp_path)

#' View dataset structure:
str(comp_raw)



#' 
#' **Prepare dataframes**
#'  

#' Set several grouping variables as factors
comp_raw$year <- as.factor(comp_raw$year)
comp_raw$block <- as.factor(comp_raw$block)
comp_raw$plot <- as.factor(comp_raw$plot)
comp_raw$species_code <- as.factor(comp_raw$species_code)

#' Remove non-plant objects, unknown species that could not be uniquely identified,
#' and early-phenology species that could not be reliably detected 
comp <- comp_raw %>%  filter(species_code !="cow" & species_code!= "rock") %>%  # Remove cowpies and rocks
                      filter(species_code !="uk_g" & species_code!= "uk_f") %>%  # Remove unknown grasses and forbs
                      filter(species_code !="dr", species_code!="castilleja_sp", species_code!="mertensia_sp", species_code!="lo" ) # Remove early-phenology species (not reliably detected)

#' Recode to pool species that were challenging to differentiate but functionally-similar
levels(comp$species_code)[levels(comp$species_code)=="ch"] <- "cg" # Chondrosum (syn. Bouteloua) gracilis and C. hirsuta
levels(comp$species_code)[levels(comp$species_code)=="pl"] <- "td" # Podospermum lacinatum and Tragopogon dubius

#' Check the revised species' levels to make sure this has worked
sort(unique(comp$species_code))

#'  Create long matrix of summed aerial cover for each species within 
#'  each 0.5m2 plot and year.
comp_long_pres <- comp %>% 
  select(year, block, plot, species_code, area_clipped) %>%
  group_by(year, block, plot, species_code) %>%
  summarise(cover = sum(area_clipped)) 

#' Incorporate 'present' species: 
#' 
#'   In each plot, plants were mapped within a 0.5m2 area, but any other plants 
#'   within the larger 1m2 area were also noted. These are designated by richness_only == "Y"
#'   in the 'comp' dataframe, but have currently have a cover value of NA in the 
#'   long matrix (comp_long_pres) because they were not mapped within the 0.5m2 area.
#'   
#'   We will include these species in compositional analyses by giving them a standard,
#'   minimum cover value in the 'area_clipped' column. The lowest cover value for 
#'   mapped species is 0.095, so we will use a rounded value of 0.1 for these species:
#'    
min(comp_long_pres$cover,na.rm=T)   #minimum cover of mapped species
comp_long_pres <- comp_long_pres %>%  replace(is.na(.), 0.1)  #assign minimum, standard cover value

#' Create Wide matrix of plot- and species-level aerial cover
comp_wide <- comp_long_pres %>% spread( species_code, cover)  %>%  replace(is.na(.), 0)

#' Revised long matrix of total species cover
comp_long <- comp_wide %>% gather( af:vt, key="species", value="cover") 

# Create separate wide dataframes for each year
comp_2018 <- comp_wide %>% filter( year =="2018")
comp_2019 <- comp_wide %>% filter( year =="2019")
comp_2020 <- comp_wide %>% filter( year =="2020")

#' Create separate dataframes to store plot and treatment information 
#' 
#' Create info dataframe sorted by year, block, plot
comp_wide_info <- unique(comp[c("year","block","plot","rain_trt","graze_trt")])
comp_wide_info <- comp_wide_info[with(comp_wide_info, order(year, block, plot)),]
#  Convert several variables to factors
comp_wide_info$block <- as.factor(comp_wide_info$block)
comp_wide_info$year <- as.factor(comp_wide_info$year)
#  Set reference levels for fixed effects (rain=C, graze=UG, year=2018)
comp_wide_info$rain_trt <- factor(comp_wide_info$rain_trt, levels=c("D","W","C"))
comp_wide_info$graze_trt <- factor(comp_wide_info$graze_trt, levels=c("DS","GS","UG"))
comp_wide_info$year <- factor(comp_wide_info$year, levels=c("2019","2020","2018"))
# Create separate info dataframes for each year
comp_2018_info <- comp_wide_info %>% filter(year==2018)
comp_2019_info <- comp_wide_info %>% filter(year==2019)
comp_2020_info <- comp_wide_info %>% filter(year==2020)

#' Create separate compositional dataframes to store species abundances only
#' 
#' Absolute, untransformed cover values for all years (comp_wide) and separate years
comp_wide <- comp_wide %>% ungroup() %>% select(-year, - block, -plot)
comp_2018 <- comp_2018 %>% ungroup() %>% select(-year, -block, -plot)
comp_2019 <- comp_2019 %>% ungroup() %>% select(-year, -block, -plot)
comp_2020 <- comp_2020 %>% ungroup() %>% select(-year, -block, -plot)

#' Absolute, squareroot-transformed cover values for separate years
comp_2018_sqrt <- sqrt(comp_2018)
comp_2019_sqrt <- sqrt(comp_2019)
comp_2020_sqrt <- sqrt(comp_2020)

#' Convert absolute abundances to relative abundances
spp_df_rel <- decostand(comp_wide, "total")  # All years
spp_2018_rel <- decostand(comp_2018, "total")
spp_2019_rel <- decostand(comp_2019, "total")
spp_2020_rel <- decostand(comp_2020, "total")

#' SQRT-transform relative covers for analysis
spp_df_rel_sqrt <- decostand(sqrt(comp_wide),"total")
spp_2018_rel_sqrt <- decostand(sqrt(comp_2018),"total")
spp_2019_rel_sqrt <- decostand(sqrt(comp_2019),"total")
spp_2020_rel_sqrt <- decostand(sqrt(comp_2020),"total")

#' Summarize species' average relative abundances 
#'
#' Across all years -SQRT transformed
spp_mean_rel_sqrt <- spp_df_rel_sqrt %>% summarise_all(mean) %>% t()
spp_mean_rel_sqrt <- data.frame(spp_mean_rel_sqrt)
colnames(spp_mean_rel_sqrt)[1] <- 'sqrt_rel_abun_allyrs'
spp_mean_rel_sqrt$sqrt_rel_abun_allyrs  <- 100*spp_mean_rel_sqrt$sqrt_rel_abun_allyrs
spp_mean_rel_sqrt$species <- row.names(spp_mean_rel_sqrt)
# 2018 only - SQRT transformed relative cover
spp_mean_rel_2018_sqrt <- spp_2018_rel_sqrt %>% summarise_all(mean) %>% t()
spp_mean_rel_2018_sqrt <- data.frame(spp_mean_rel_2018_sqrt)
spp_mean_rel_2018_sqrt$spp_mean_rel_2018_sqrt  <- 100*spp_mean_rel_2018_sqrt$spp_mean_rel_2018_sqrt
# 2018 only - Relative cover
spp_mean_rel_2018 <- spp_2018_rel %>% summarise_all(mean) %>% t()
spp_mean_rel_2018 <- data.frame(spp_mean_rel_2018)
spp_mean_rel_2018$spp_mean_rel_2018  <- 100*spp_mean_rel_2018

#' View list of species in order from most to least common in 2018
#'    based on raw relative abundances
spp_mean_rel_2018 %>% arrange(-spp_mean_rel_2018)


# Create long dataframe of sqrt-transformed relative abundances with plot identifier
spp_df_rel_sqrt_long_full <- cbind(comp_wide_info, spp_df_rel_sqrt)
spp_df_rel_sqrt_long_full$b_p <- paste(spp_df_rel_sqrt_long_full$block, spp_df_rel_sqrt_long_full$plot, sep="_")
spp_df_rel_sqrt_long  <- spp_df_rel_sqrt_long_full  %>% 
  select(-block, - plot, -rain_trt, -graze_trt) %>% 
  gather(af:vt, key="species", value="rel_cov")

# Create long dataframes of sqrt-transformed absolute abundances in each year
# 2018
comp_2018_sqrt_long_full <- cbind(comp_2018_info, comp_2018_sqrt)
comp_2018_sqrt_long_full$b_p <- paste(comp_2018_sqrt_long_full$block, comp_2018_sqrt_long_full$plot, sep="_")
comp_2018_sqrt_long  <- comp_2018_sqrt_long_full  %>% 
  select(-block, - plot, -rain_trt, -graze_trt) %>% 
  gather(af:vt, key="species", value="abs_cov")
# 2019
comp_2019_sqrt_long_full <- cbind(comp_2019_info, comp_2019_sqrt)
comp_2019_sqrt_long_full$b_p <- paste(comp_2019_sqrt_long_full$block, comp_2019_sqrt_long_full$plot, sep="_")
comp_2019_sqrt_long  <- comp_2019_sqrt_long_full  %>% 
  select(-block, - plot, -rain_trt, -graze_trt) %>% 
  gather(af:vt, key="species", value="abs_cov")
# 2020
comp_2020_sqrt_long_full <- cbind(comp_2020_info, comp_2020_sqrt)
comp_2020_sqrt_long_full$b_p <- paste(comp_2020_sqrt_long_full$block, comp_2020_sqrt_long_full$plot, sep="_")
comp_2020_sqrt_long  <- comp_2020_sqrt_long_full  %>% 
  select(-block, - plot, -rain_trt, -graze_trt) %>% 
  gather(af:vt, key="species", value="abs_cov")
str(comp_2020_sqrt_long)

# Create an info dataframe (plot and treatment identifiers) to bind with any comp data from any year
spp_df_rel_long_info <- spp_df_rel_sqrt_long_full %>% filter(year == "2018") %>% select(block, b_p, graze_trt, rain_trt)
spp_df_rel_long_info$b_p <- as.factor(spp_df_rel_long_info$b_p)



#' *Estimate composition metrics*
#' 


#' *1. Richness, Evenness, & Shannon Diversity* 
#' 

#' Select and prepare sqrt-transformed relative abundances
comm_df_sqrt <- spp_df_rel_sqrt_long 
comm_df_sqrt$year <- factor(comm_df_sqrt$year, levels=c("2018","2019","2020"))
comm_df_sqrt$year_num <- as.integer(as.character(comm_df_sqrt$year))

#'Estimate Richness & Evenness
structure_sqrt <- community_structure(comm_df_sqrt, time.var="year",
                                      replicate.var = "b_p",abundance.var = "rel_cov") # for Evar evenness measure
#' Shannon's diversity
div_sqrt <- community_diversity(comm_df_sqrt,  time.var="year",
                                replicate.var = "b_p",abundance.var = "rel_cov") # for Shannon's diversity measure

#' Create a single dataframe for diversity metrics
structure_sqrt$shannon_div <- div_sqrt$Shannon
structure_raw_sqrt <- left_join(structure_sqrt, spp_df_rel_long_info, by = 'b_p') %>% rename(rich = richness, even = Evar, shan = shannon_div)

#' Scale diversity metrics for analysis
structure_all_sqrt <-  structure_raw_sqrt  %>% ungroup() %>% mutate_if(is.numeric, scale)
structure_all_sqrt$b_p <- as.factor(structure_all_sqrt$b_p)
str(structure_all_sqrt, give.attr=F)


#' Format dataframe for figures
structure_sum_sqrt <- structure_raw_sqrt %>% select(-b_p, -block) %>%
  gather(c(rich, even, shan), key="metric", value="value") %>%
  group_by(rain_trt, graze_trt, year, metric) %>%
  summarise( mean = mean(value),
             se= se(value))
structure_sum_sqrt$graze_trt <- gsub("UG", "Ungrazed", structure_sum_sqrt$graze_trt)
structure_sum_sqrt$graze_trt <- gsub("DS", "Dormant-Season Grazing", structure_sum_sqrt$graze_trt)
structure_sum_sqrt$graze_trt <- gsub("GS", "Growing-Season Grazing", structure_sum_sqrt$graze_trt)
structure_sum_sqrt$year_num <- as.numeric(as.character(structure_sum_sqrt$year))
structure_sum_sqrt$metric <-  factor(structure_sum_sqrt$metric,
                                     levels = c("rich","even","shan"),
                                     labels = c("Richness","Evenness","Diversity"))
structure_sum_sqrt$rain_trt <- factor(structure_sum_sqrt$rain_trt, levels=c("D","C","W"))
structure_sum_sqrt$graze_trt <- factor(structure_sum_sqrt$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))


#' Format dataframes for analysis
#' 
#' Create several sets of unique treatment contrasts for models
#' 
#  Set 1 - View, then reset contrasts for year
contrasts(structure_all_sqrt$rain_trt)
contrasts(structure_all_sqrt$graze_trt)
contrasts(structure_all_sqrt$year)
structure_all_sqrt$year <- factor(structure_all_sqrt$year, levels=c("2020","2019","2018"))
#' Set 2 
structure_all_sqrt_2 <- structure_all_sqrt
structure_all_sqrt_2$year <- factor(structure_all_sqrt_2$year, levels=c("2018","2019","2020"))
structure_all_sqrt_2$rain_trt <- factor(structure_all_sqrt_2$rain_trt, levels=c("C","D","W"))
structure_all_sqrt_2$graze_trt <- factor(structure_all_sqrt_2$graze_trt, levels=c("UG","GS","DS"))
#' Set 3
structure_all_sqrt_3 <- structure_all_sqrt
structure_all_sqrt_3$year <- factor(structure_all_sqrt_3$year, levels=c("2018","2019","2020"))
structure_all_sqrt_3$rain_trt <- factor(structure_all_sqrt_3$rain_trt, levels=c("W","C","D"))
structure_all_sqrt_3$graze_trt <- factor(structure_all_sqrt_3$graze_trt, levels=c("UG","GS","DS"))
#' Set 4
structure_all_sqrt_4 <- structure_all_sqrt
structure_all_sqrt_4$year <- factor(structure_all_sqrt_4$year, levels=c("2020","2019","2018"))
structure_all_sqrt_4$rain_trt <- factor(structure_all_sqrt_4$rain_trt, levels=c("C","D","W"))
structure_all_sqrt_4$graze_trt <- factor(structure_all_sqrt_4$graze_trt, levels=c("UG","GS","DS"))



#' *2. Functional Groups* 
#'  

#' Starting dataframe
fun_grp_df <- comp %>% select(species_code, fun_grp) %>% 
  rename("species"="species_code")

#' View unique functional groups
unique(fun_grp_df$fun_grp)

#' Get dataframe of species' functional groups and relative abundances 
spp_info <- arrange(unique(fun_grp_df[c("species","fun_grp")]),species)
spp_info <- merge(spp_info, spp_mean_rel_sqrt, by="species")
spp_info

#' Join functional group dataframe with each relative abundances dataframe
spp_df_rel_sqrt_long_fgs <- left_join(spp_df_rel_sqrt_long, spp_info, by="species")

#' Join functional group dataframe with 2018, 2019 and 2020 dataframe
comp_2020_sqrt_long_fgs <- left_join(comp_2020_sqrt_long, spp_info, by="species")
comp_2019_sqrt_long_fgs <- left_join(comp_2019_sqrt_long, spp_info, by="species")
comp_2018_sqrt_long_fgs <- left_join(comp_2018_sqrt_long, spp_info, by="species")

#' Summarize plot-level relative abundance in each year by functional group
# All Years
spp_df_rel_sqrt_long_fg <- spp_df_rel_sqrt_long_fgs %>% 
  group_by(year, b_p, fun_grp) %>%  summarise(fg_cov = sum(rel_cov)) 
spp_df_rel_sqrt_long_fg <- left_join(spp_df_rel_sqrt_long_fg, spp_df_rel_long_info, by="b_p")
# 2018
comp_2018_sqrt_long_fg <- comp_2018_sqrt_long_fgs %>% 
  group_by(year, b_p, fun_grp) %>%  summarise(fg_cov = sum(abs_cov)) 
comp_2018_sqrt_long_fg <- left_join(comp_2018_sqrt_long_fg, spp_df_rel_long_info, by="b_p")
#2019
comp_2019_sqrt_long_fg <- comp_2019_sqrt_long_fgs %>% 
  group_by(year, b_p, fun_grp) %>%  summarise(fg_cov = sum(abs_cov)) 
comp_2019_sqrt_long_fg <- left_join(comp_2019_sqrt_long_fg, spp_df_rel_long_info, by="b_p")
#2020
comp_2020_sqrt_long_fg <- comp_2020_sqrt_long_fgs %>% 
  group_by(year, b_p, fun_grp) %>%  summarise(fg_cov = sum(abs_cov)) 
comp_2020_sqrt_long_fg <- left_join(comp_2020_sqrt_long_fg, spp_df_rel_long_info, by="b_p")

#' Save 2018 relative covers to add as extra info in dataframes
#' 
#  Squareroot-transformed relative abundances - SCALED
c4_18_sqrt <- spp_df_rel_sqrt_long_fg %>% filter(year =="2018") %>% filter(fun_grp == "G_C4")%>% 
  ungroup() %>% mutate_if(is.numeric, scale)%>% select(b_p,fg_cov) %>% rename(c4_18 = fg_cov)
c3_18_sqrt <- spp_df_rel_sqrt_long_fg %>% filter(year =="2018") %>% filter(fun_grp == "G_C3") %>% 
  ungroup() %>%  mutate_if(is.numeric, scale)%>%select(b_p,fg_cov) %>% rename(c3_18 = fg_cov)
ann_18_sqrt <- spp_df_rel_sqrt_long_fg %>% filter(year =="2018") %>% filter(fun_grp == "Ann_Bi") %>% 
  ungroup() %>% mutate_if(is.numeric, scale)%>%  select(b_p,fg_cov) %>% rename(ann_18 = fg_cov)
f_18_sqrt <- spp_df_rel_sqrt_long_fg %>% filter(year =="2018") %>% filter(fun_grp == "F")%>% 
  ungroup() %>%  mutate_if(is.numeric, scale) %>%select(b_p,fg_cov) %>% rename(f_18 = fg_cov)

#' Save absolute covers in each year
#' 
# 2018 - Squareroot-transformed absolute abundances 
c4_18_abs <- comp_2018_sqrt_long_fg %>% filter(fun_grp == "G_C4")%>% 
  ungroup() %>% select(b_p,fg_cov) %>% rename(c4_abs = fg_cov)
c3_18_abs <- comp_2018_sqrt_long_fg %>% filter(fun_grp == "G_C3") %>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(c3_abs = fg_cov)
ann_18_abs <- comp_2018_sqrt_long_fg %>% filter(fun_grp == "Ann_Bi") %>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(ann_abs = fg_cov)
f_18_abs <- comp_2018_sqrt_long_fg %>% filter(fun_grp == "F")%>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(f_abs = fg_cov)
# Save all functional groups in one joined dataframe
fg_abs_18 <- left_join(c4_18_abs, c3_18_abs, by='b_p')
fg_abs_18 <- left_join(fg_abs_18, ann_18_abs, by='b_p')
fg_abs_18 <- left_join(fg_abs_18, f_18_abs, by='b_p')
fg_abs_full_18 <- left_join(fg_abs_18, spp_df_rel_long_info, by="b_p")
fg_abs_full_18$year <- 2018
str(fg_abs_full_18)

# 2019 - Squareroot-transformed absolute abundances 
c4_19_abs <- comp_2019_sqrt_long_fg %>% filter(fun_grp == "G_C4")%>% 
  ungroup() %>% select(b_p,fg_cov) %>% rename(c4_abs = fg_cov)
c3_19_abs <- comp_2019_sqrt_long_fg %>% filter(fun_grp == "G_C3") %>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(c3_abs = fg_cov)
ann_19_abs <- comp_2019_sqrt_long_fg %>% filter(fun_grp == "Ann_Bi") %>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(ann_abs = fg_cov)
f_19_abs <- comp_2019_sqrt_long_fg %>% filter(fun_grp == "F")%>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(f_abs = fg_cov)
fg_abs_19 <- left_join(c4_19_abs, c3_19_abs, by='b_p')
fg_abs_19 <- left_join(fg_abs_19, ann_19_abs, by='b_p')
fg_abs_19 <- left_join(fg_abs_19, f_19_abs, by='b_p')
str(fg_abs_19)
fg_abs_full_19 <- left_join(fg_abs_19, spp_df_rel_long_info, by="b_p")
fg_abs_full_19$year <- 2019

# 2020 - Squareroot-transformed absolute abundances 
c4_20_abs <- comp_2020_sqrt_long_fg %>% filter(fun_grp == "G_C4")%>% 
  ungroup() %>% select(b_p,fg_cov) %>% rename(c4_abs = fg_cov)
c3_20_abs <- comp_2020_sqrt_long_fg %>% filter(fun_grp == "G_C3") %>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(c3_abs = fg_cov)
ann_20_abs <- comp_2020_sqrt_long_fg %>% filter(fun_grp == "Ann_Bi") %>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(ann_abs = fg_cov)
f_20_abs <- comp_2020_sqrt_long_fg %>% filter(fun_grp == "F")%>% 
  ungroup() %>%  select(b_p,fg_cov) %>% rename(f_abs = fg_cov)
fg_abs <- left_join(c4_20_abs, c3_20_abs, by='b_p')
fg_abs <- left_join(fg_abs, ann_20_abs, by='b_p')
fg_abs <- left_join(fg_abs, f_20_abs, by='b_p')
str(fg_abs)
fg_abs_full <- left_join(fg_abs, spp_df_rel_long_info, by="b_p")
fg_abs_full$year <- 2020

#' Combine all years
fg_abs_full_allyears <- data.frame(rbind(fg_abs_full_18, fg_abs_full_19, fg_abs_full))
fg_abs_full_allyears$year <- as.factor(fg_abs_full_allyears$year)


#' Format dataframe for figures - Filter to relevant functional groups and summarize means / SEs
spp_df_rel_sqrt_long_fg_sum <- spp_df_rel_sqrt_long_fg %>%
  ungroup() %>% select(-b_p, -block) %>%
  filter(fun_grp != "G_NA") %>%
  filter(fun_grp != "F_NA") %>%
  filter(fun_grp != "S") %>%
  group_by(year, rain_trt, graze_trt, fun_grp) %>%
  summarise( mean = mean(fg_cov),
             se= se(fg_cov))
spp_df_rel_sqrt_long_fg_sum$year <- as.numeric(as.character(spp_df_rel_sqrt_long_fg_sum$year))
spp_df_rel_sqrt_long_fg_sum$graze_trt <- gsub("UG", "Ungrazed", spp_df_rel_sqrt_long_fg_sum$graze_trt)
spp_df_rel_sqrt_long_fg_sum$rain_trt <- factor(spp_df_rel_sqrt_long_fg_sum$rain_trt, levels=c("C","D","W"))
spp_df_rel_sqrt_long_fg_sum$graze_trt <- factor(spp_df_rel_sqrt_long_fg_sum$graze_trt, levels=c("Ungrazed","GS","DS"))
spp_df_rel_sqrt_long_fg_sum <- spp_df_rel_sqrt_long_fg_sum %>% rename (metric = fun_grp) %>%
  select (year, rain_trt, graze_trt, metric, mean, se)


#' Format dataframes for analysis 
#'
#' First set generic dataframe
fg_df_sqrt <- spp_df_rel_sqrt_long_fg 
#' Create several sets of unique treatment contrasts for models
#  Set 1 - View, then reset contrasts for year
contrasts(fg_df_sqrt$rain_trt)
contrasts(fg_df_sqrt$graze_trt)
contrasts(fg_df_sqrt$year)
fg_df_sqrt$year <- factor(fg_df_sqrt$year, levels=c("2020","2019","2018"))
#' Set 2 (to view other comparisons)
fg_df_sqrt_2 <- fg_df_sqrt
fg_df_sqrt_2$year <- factor(fg_df_sqrt_2$year, levels=c("2018","2019","2020"))
fg_df_sqrt_2$rain_trt <- factor(fg_df_sqrt_2$rain_trt, levels=c("C","D","W"))
fg_df_sqrt_2$graze_trt <- factor(fg_df_sqrt_2$graze_trt, levels=c("UG","GS","DS"))
#' Set 3
fg_df_sqrt_3 <- fg_df_sqrt
fg_df_sqrt_3$year <- factor(fg_df_sqrt_3$year, levels=c("2018","2019","2020"))
fg_df_sqrt_3$rain_trt <- factor(fg_df_sqrt_3$rain_trt, levels=c("W","C","D"))
fg_df_sqrt_3$graze_trt <- factor(fg_df_sqrt_3$graze_trt, levels=c("UG","GS","DS"))
#' Set 4
fg_df_sqrt_4 <- fg_df_sqrt
fg_df_sqrt_4$year <- factor(fg_df_sqrt_4$year, levels=c("2020","2019","2018"))
fg_df_sqrt_4$rain_trt <- factor(fg_df_sqrt_4$rain_trt, levels=c("C","D","W"))
fg_df_sqrt_4$graze_trt <- factor(fg_df_sqrt_4$graze_trt, levels=c("UG","GS","DS"))


#' Create separate functional group dataframes that are scaled
#' Square root transformed (including separate contrasts)
g_c4_sqrt <- fg_df_sqrt %>% filter(fun_grp == "G_C4") %>% ungroup() %>% mutate_if(is.numeric, scale)
g_c4_sqrt_2 <- fg_df_sqrt_2 %>% filter(fun_grp == "G_C4") %>% ungroup() %>% mutate_if(is.numeric, scale)
g_c4_sqrt_3 <- fg_df_sqrt_3 %>% filter(fun_grp == "G_C4") %>% ungroup() %>% mutate_if(is.numeric, scale)
g_c4_sqrt_4 <- fg_df_sqrt_4 %>% filter(fun_grp == "G_C4") %>% ungroup() %>% mutate_if(is.numeric, scale)

g_c3_sqrt <- fg_df_sqrt %>% filter(fun_grp == "G_C3")%>% ungroup()  %>% mutate_if(is.numeric, scale)
g_c3_sqrt_2 <- fg_df_sqrt_2 %>% filter(fun_grp == "G_C3") %>% ungroup() %>% mutate_if(is.numeric, scale)
g_c3_sqrt_3 <- fg_df_sqrt_3 %>% filter(fun_grp == "G_C3") %>% ungroup() %>% mutate_if(is.numeric, scale)
g_c3_sqrt_4 <- fg_df_sqrt_4 %>% filter(fun_grp == "G_C3") %>% ungroup() %>% mutate_if(is.numeric, scale)

ann_sqrt <- fg_df_sqrt %>% filter(fun_grp == "Ann_Bi")%>% ungroup()  %>% mutate_if(is.numeric, scale)
ann_sqrt_2 <- fg_df_sqrt_2 %>% filter(fun_grp == "Ann_Bi")%>% ungroup()  %>% mutate_if(is.numeric, scale)
ann_sqrt_3 <- fg_df_sqrt_3 %>% filter(fun_grp == "Ann_Bi")%>% ungroup()  %>% mutate_if(is.numeric, scale)
ann_sqrt_4 <- fg_df_sqrt_4 %>% filter(fun_grp == "Ann_Bi")%>% ungroup()  %>% mutate_if(is.numeric, scale)

f_sqrt <- fg_df_sqrt %>% filter(fun_grp == "F") %>% ungroup() %>% mutate_if(is.numeric, scale)
f_sqrt_2 <- fg_df_sqrt_2 %>% filter(fun_grp == "F") %>% ungroup() %>% mutate_if(is.numeric, scale)
f_sqrt_3 <- fg_df_sqrt_3 %>% filter(fun_grp == "F") %>% ungroup() %>% mutate_if(is.numeric, scale)
f_sqrt_4 <- fg_df_sqrt_4 %>% filter(fun_grp == "F") %>% ungroup() %>% mutate_if(is.numeric, scale)



#' *3. Combined community metric dataframes*
#'

# Prep diversity summary metrics
structure_sum_sqrt <- structure_sum_sqrt %>% select (year, rain_trt, graze_trt, metric, mean, se) 
spp_df_rel_sqrt_long_fg_sum <- spp_df_rel_sqrt_long_fg_sum %>% select (year, rain_trt, graze_trt, metric, mean, se)
spp_df_rel_sqrt_long_fg_sum$graze_trt <- gsub("DS", "Dormant-Season Grazing", spp_df_rel_sqrt_long_fg_sum$graze_trt)
spp_df_rel_sqrt_long_fg_sum$graze_trt <- gsub("GS", "Growing-Season Grazing", spp_df_rel_sqrt_long_fg_sum$graze_trt)
spp_df_rel_sqrt_long_fg_sum$graze_trt <- factor(spp_df_rel_sqrt_long_fg_sum$graze_trt, levels=c( "Ungrazed", "Growing-Season Grazing","Dormant-Season Grazing"))
spp_df_rel_sqrt_long_fg_sum$metric <- as.factor(spp_df_rel_sqrt_long_fg_sum$metric)
structure_sum_sqrt$year <- as.numeric(as.character(structure_sum_sqrt$year))

# Combine diversity summary metrics with summary functional group metrics
community_sum_combined <- rbind(structure_sum_sqrt, spp_df_rel_sqrt_long_fg_sum)
community_sum_combined$metric <- gsub("F", "Forb", community_sum_combined$metric)
community_sum_combined$metric <- factor(community_sum_combined$metric, 
                                        levels = c("Diversity","Richness", "Evenness","G_C4", "G_C3", "Ann_Bi", "Forb"))
community_sum_combined$rain_trt <- factor(community_sum_combined$rain_trt, levels=c("D", "C", "W"))

#' Set appropriate vector settings
community_sum_combined$year <- as.numeric(as.character(community_sum_combined$year))
community_sum_combined$rain_trt <- factor(community_sum_combined$rain_trt, levels=c("D", "C", "W"))
community_sum_combined$graze_trt <- factor(community_sum_combined$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))
community_sum_combined$treat <- as.factor(paste(community_sum_combined$graze_trt, community_sum_combined$rain_trt,sep="-"))
levels(community_sum_combined$treat)

#' Create separate mean and standard error dataframes
community_sum_combined_mean <- community_sum_combined %>% select(-se) %>% spread(key=metric, value=mean) %>% arrange(graze_trt, year)
community_sum_combined_se <- community_sum_combined %>% select(-mean) %>% spread(key=metric, value=se) %>% arrange(graze_trt, year)


# For figures below, create a special dataframe for x-axis labels
year_seq <- data.frame(
  DS = rep(c('2018','2019','2020'),each=3),
  GS = rep(c('2018.15','2019.15','2020.15'),each=3),
  Ungrazed = rep(c('2017.85','2018.85','2019.85'),each=3))
year_seq2 <- year_seq %>% gather(DS:Ungrazed, key=graze_trt, value=year_seq)
year_seq2$graze_trt <- factor(year_seq2$graze_trt, levels=c("Ungrazed","GS","DS"))
year_seq2 <- year_seq2 %>% arrange(graze_trt, year_seq)
community_sum_combined_mean$year_seq <- as.numeric(year_seq2$year_seq)

# Create text dataframe for grazing labels
graze_labs <- data.frame(
  x=c(2017.75,2018,2018.25,2018.75,2019,2019.25,2019.75,2020,2020.25), 
  y=0, 
  y1=0.2,
  y2=.15,
  y3=15,
  label=c("UnG.","DS","GS"))

#' Create separate functional gorup dataframes
c3_sum <- community_sum_combined_mean %>% select(year:graze_trt, G_C3, year_seq)
c3_sum$se <-  community_sum_combined_se$G_C3
c4_sum <- community_sum_combined_mean %>% select(year:graze_trt, G_C4, year_seq)
c4_sum$se <-  community_sum_combined_se$G_C4
f_sum <- community_sum_combined_mean %>% select(year:graze_trt, Forb, year_seq)
f_sum$se <-  community_sum_combined_se$Forb
ann_sum <- community_sum_combined_mean %>% select(year:graze_trt, Ann_Bi, year_seq)
ann_sum$se <-  community_sum_combined_se$Ann_Bi
rich_sum <- community_sum_combined_mean %>% select(year:graze_trt, Richness, year_seq)
rich_sum$se <-  community_sum_combined_se$Richness
even_sum <- community_sum_combined_mean %>% select(year:graze_trt, Evenness, year_seq)
even_sum$se <-  community_sum_combined_se$Evenness
div_sum <- community_sum_combined_mean %>% select(year:graze_trt, Diversity, year_seq)
even_sum$se <-  community_sum_combined_se$Diversity



#'
#' **ANALYSIS: Species composition**
#'  

#' *PerMANOVA - Table S6* 
#' 
spp_perm_all_sqrt <- adonis2(spp_df_rel_sqrt ~ rain_trt * graze_trt * year, strata=comp_wide_info$block, data=comp_wide_info, permutations = 999, method="bray",by="terms")
spp_perm_all_sqrt


#' *NMDS*  
#'    
#'  Run across all years -- SQRT-transformed relative abundances
nmds_allyrs_sqrt <- metaMDS(spp_df_rel_sqrt, k=3, trymax=100, distance = "bray")
nmds_allyrs_sqrt
stressplot(nmds_allyrs_sqrt)


#' Save data for NMDS visualization
#' 
# Plot scores
nmds_scores <- data.frame(nmds_allyrs_sqrt$points)
plot_scores_all_sqrt_nmds <-  cbind(comp_wide_info, nmds_scores)
plot_scores_all_sqrt_nmds$b_p <- paste(plot_scores_all_sqrt_nmds$block, plot_scores_all_sqrt_nmds$plot, sep="_")
# Species scores
nmds_species <- data.frame(nmds_allyrs_sqrt$species)
nmds_species$species <- row.names(nmds_species)
nmds_species <- merge(nmds_species, spp_info, by='species')
# Separate by year to show scores through time)
plot_scores_all_2018 <- plot_scores_all_sqrt_nmds %>% filter(year == 2018) %>% 
  select( MDS1, MDS2) %>% rename(MDS1_2018 = MDS1, MDS2_2018 = MDS2)
plot_scores_all_2019 <- plot_scores_all_sqrt_nmds %>% filter(year == 2019) %>% 
  select( MDS1, MDS2) %>% rename(MDS1_2019 = MDS1, MDS2_2019 = MDS2)
plot_scores_all_2020 <- plot_scores_all_sqrt_nmds %>% filter(year == 2020) %>% 
  select( MDS1, MDS2) %>% rename(MDS1_2020 = MDS1, MDS2_2020 = MDS2)
plot_scores_all_info <- plot_scores_all_sqrt_nmds %>% filter(year==2020) %>% select(block, plot, graze_trt, rain_trt)
plot_scores_2018_20 <- cbind(plot_scores_all_info, plot_scores_all_2018, plot_scores_all_2019, plot_scores_all_2020)

#' Adjust final formatting for figure
plot_scores_all_sqrt_no2019 <- plot_scores_all_sqrt_nmds %>% filter ( year != '2019')
plot_scores_all_sqrt_no2019$graze_trt <- gsub("UG","Ungrazed", plot_scores_all_sqrt_no2019$graze_trt)
plot_scores_all_sqrt_no2019$graze_trt <- gsub("GS","Growing-Season Grazing", plot_scores_all_sqrt_no2019$graze_trt)
plot_scores_all_sqrt_no2019$graze_trt <- gsub("DS","Dormant-Season Grazing", plot_scores_all_sqrt_no2019$graze_trt)
plot_scores_all_sqrt_no2019$graze_trt <- factor(plot_scores_all_sqrt_no2019$graze_trt, levels=c("Growing-Season Grazing","Dormant-Season Grazing","Ungrazed"))
plot_scores_all_sqrt_no2019$year <- factor(plot_scores_all_sqrt_no2019$year, levels=c("2018","2020"))
plot_scores_all_sqrt_no2019$rain_trt <- factor(plot_scores_all_sqrt_no2019$rain_trt, levels=c("D","C","W"))
plot_scores_2018_20$rain_trt <- factor(plot_scores_2018_20$rain_trt, levels=c("D","C","W"))
plot_scores_2018_20$graze_trt <- gsub("UG","Ungrazed", plot_scores_2018_20$graze_trt)
plot_scores_2018_20$graze_trt <- gsub("GS","Growing-Season Grazing", plot_scores_2018_20$graze_trt)
plot_scores_2018_20$graze_trt <- gsub("DS","Dormant-Season Grazing", plot_scores_2018_20$graze_trt)
plot_scores_2018_20$graze_trt <- factor(plot_scores_2018_20$graze_trt, levels=c("Growing-Season Grazing","Dormant-Season Grazing","Ungrazed"))


#' *NMDS FIGURE - Figure S6*
#' 
#' Composition over 3 years, from a single NMDS (3 axes) with all years
#'    and vectors showing change over time 

fig_nmds_sqrt_all_vec <- ggplot() +
  #geom_text(aes(x=MDS1, y=MDS2, label=rownames(plot_species_all_sqrt)), cex=3, data=plot_species_all_sqrt) +
  geom_point(aes(x=MDS1, y=MDS2, col=rain_trt, shape=graze_trt, alpha=year), cex=2.5, data=plot_scores_all_sqrt_no2019) +
  geom_segment(aes(x=MDS1_2018, xend=MDS1_2019, y=MDS2_2018, yend=MDS2_2019, color=rain_trt, lty=graze_trt), size=0.6, data=plot_scores_2018_20) +
  geom_segment(aes(x=MDS1_2019, xend=MDS1_2020, y=MDS2_2019, yend=MDS2_2020, color=rain_trt, lty=graze_trt), size=0.6, data=plot_scores_2018_20) +
  xlim(-0.9,0.9) + ylim(-0.75,0.9) +
  geom_text(aes(x=0.35, y=0.75, label="PerMANOVA: \n  Rainfall p<=0.001\n  Rainfall:Graze p<=0.001 \n  Year<=0.001"), hjust = 0, col="black", cex=3)+
  geom_text(aes(x=-0.85, y=0.85, label="2018-2020"), fontface=2,cex=4,hjust=0) +
  labs(x="NMDS 1 - Species Composition", y="NMDS 2 - Species Composition", col="Rain trt.", shape="Graze trt.", lty="Graze trt.", alpha="Year") +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_alpha_manual(values=c(0.2,0.9)) +
  #scale_color_manual(values=c("gray85","gray15", "gray50" , "tan1", "tan3", "tan4","lightgreen", "limegreen", "green4")) +
  theme_bw()
fig_nmds_sqrt_all_vec

#' Same NDMS, but showing species labels
fig_nmds_sqrt_species <- ggplot(data = nmds_species) +
  geom_text(aes(x=-1.5, y=1.1, label="Species weighted by\nmean rel. abundance\n(sqrt-transformed)\n2018-2020"), fontface=2,cex=3,hjust=0) +
  geom_text(aes( x=MDS1, y=MDS2, size= sqrt_rel_abun_allyrs, label=species)) +
  geom_point(aes(x=MDS1, y=MDS2, col=rain_trt, shape=graze_trt), alpha=0.2, cex=2.5, data=plot_scores_all_sqrt_no2019) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_size_continuous(range=c(2,5))+
  labs(x="NMDS 1 - Species Composition", y="NMDS 2 - Species Composition", col="Rain trt.", shape="Graze trt.", lty="Graze trt.", alpha="Year") +
  theme_bw() +
  # ylim(-2,1.5) +
  theme(legend.position = 'none')
fig_nmds_sqrt_species





#' **ANALYSIS: Linear Mixed Models**
#'  
#'  Compositional Responses:
#'  
#'  --Shannon Diversity
#'     (Supporting - Evenness, Richness)
#'  --C4 grasses
#'  --C3 grasses
#'  --Annual-Biennials
#'  --Forbs
#'  




#'  *Shannon Diversity*
#'  

#' Check data distributions
structure_all_sqrt_check <- structure_all_sqrt
structure_all_sqrt_check$graze_trt <- gsub("DS", "Dormant-Season Grazing", structure_all_sqrt_check$graze_trt)
structure_all_sqrt_check$graze_trt <- gsub("GS", "Growing-Season Grazing", structure_all_sqrt_check$graze_trt)
structure_all_sqrt_check$graze_trt <- gsub("UG", "Ungrazed", structure_all_sqrt_check$graze_trt)
structure_all_sqrt_check$graze_trt <- factor(structure_all_sqrt_check$graze_trt, levels=c("Ungrazed","Dormant-Season Grazing", "Growing-Season Grazing"))
structure_all_sqrt_check$rain_trt <- factor(structure_all_sqrt_check$rain_trt, levels=c("D","C","W"))
structure_all_sqrt_check$year <- factor(structure_all_sqrt_check$year, levels=c("2018","2019","2020"))
ggplot(data=structure_all_sqrt_check ) + 
  geom_point(aes(x=rain_trt,y=shan, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=shan, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  labs(y="Shannon Diversity\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

#' Specify models (same model, different contrasts)
div_mod_all  <- lmer(shan ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=structure_all_sqrt)
div_mod_all2 <- lmer(shan ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=structure_all_sqrt_2)
div_mod_all3 <- lmer(shan ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=structure_all_sqrt_3)
div_mod_all4 <- lmer(shan ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=structure_all_sqrt_4)

# Check residuals
plot(div_mod_all)
qqnorm(residuals(div_mod_all))

#' Model summaries
anova(div_mod_all)
r2(div_mod_all)
summary(div_mod_all)

# Save univariate fixed effects from both all contrasts
div_fix <- data.frame(round(fixef(div_mod_all),4)) %>% rename(fixef=round.fixef.div_mod_all...4.)
div_ci  <- data.frame(round(confint(div_mod_all),4))
div_fixef<- merge(div_fix, div_ci, by=0, all=TRUE) 
div_fixef$response <- "div"
div_fixef$Row.names <- gsub("rain_trt1", "rain_D", div_fixef$Row.names)
div_fixef$Row.names <- gsub("rain_trt2", "rain_W", div_fixef$Row.names)
div_fixef$Row.names <- gsub("graze_trt1", "graze_DS", div_fixef$Row.names)
div_fixef$Row.names <- gsub("graze_trt2", "graze_GS", div_fixef$Row.names)
div_fixef$Row.names <- gsub("year1", "year_2020", div_fixef$Row.names)
div_fixef$Row.names <- gsub("year2", "year_2019", div_fixef$Row.names)
div_fix2 <- data.frame(round(fixef(div_mod_all2),4))%>% rename(fixef=round.fixef.div_mod_all2...4.)
div_ci2  <- data.frame(round(confint(div_mod_all2),4))
div_fixef2<- merge(div_fix2, div_ci2, by=0, all=TRUE) 
div_fixef2$response <- "div"
div_fixef2$Row.names <- gsub("rain_trt1", "rain_C", div_fixef2$Row.names)
div_fixef2$Row.names <- gsub("rain_trt2", "rain_D", div_fixef2$Row.names)
div_fixef2$Row.names <- gsub("graze_trt1", "graze_UG", div_fixef2$Row.names)
div_fixef2$Row.names <- gsub("graze_trt2", "graze_DS", div_fixef2$Row.names)
div_fixef2$Row.names <- gsub("year1", "year_2018", div_fixef2$Row.names)
div_fixef2$Row.names <- gsub("year2", "year_2019", div_fixef2$Row.names)
div_fix3 <- data.frame(round(fixef(div_mod_all3),4))%>% rename(fixef=round.fixef.div_mod_all3...4.)
div_ci3  <- data.frame(round(confint(div_mod_all3),4))
div_fixef3<- merge(div_fix3, div_ci3, by=0, all=TRUE) 
div_fixef3$response <- "div"
div_fixef3$Row.names <- gsub("rain_trt1", "rain_W", div_fixef3$Row.names)
div_fixef3$Row.names <- gsub("rain_trt2", "rain_C", div_fixef3$Row.names)
div_fixef3$Row.names <- gsub("graze_trt1", "graze_UG", div_fixef3$Row.names)
div_fixef3$Row.names <- gsub("graze_trt2", "graze_GS", div_fixef3$Row.names)
div_fixef3$Row.names <- gsub("year1", "year_2018", div_fixef3$Row.names)
div_fixef3$Row.names <- gsub("year2", "year_2019", div_fixef3$Row.names)
div_fix4 <- data.frame(round(fixef(div_mod_all4),4))%>% rename(fixef=round.fixef.div_mod_all4...4.)
div_ci4  <- data.frame(round(confint(div_mod_all4),4))
div_fixef4<- merge(div_fix4, div_ci4, by=0, all=TRUE) 
div_fixef4$response <- "div"
div_fixef4$Row.names <- gsub("rain_trt1", "rain_C", div_fixef4$Row.names)
div_fixef4$Row.names <- gsub("rain_trt2", "rain_D", div_fixef4$Row.names)
div_fixef4$Row.names <- gsub("graze_trt1", "graze_UG", div_fixef4$Row.names)
div_fixef4$Row.names <- gsub("graze_trt2", "graze_GS", div_fixef4$Row.names)
div_fixef4$Row.names <- gsub("year1", "year_2020", div_fixef4$Row.names)
div_fixef4$Row.names <- gsub("year2", "year_2019", div_fixef4$Row.names)

div_fixef_all <- unique(rbind(div_fixef, div_fixef2, div_fixef3, div_fixef4)) %>% arrange(Row.names)
div_fixef_all


#' Shannon Diversity Fig. 3
#' 
structure_sum_div <- structure_sum_sqrt %>% filter(metric=="Diversity")
structure_sum_div$graze_trt <- factor(structure_sum_div$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))
structure_sum_div$rain_trt <- factor(structure_sum_div$rain_trt, levels=c("D","C","W"))
structure_sum_div$year <- as.numeric(as.character(structure_sum_div$year))

fig_div <- 
  ggplot(data=structure_sum_div, aes(x=year, y=mean, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.3)) +
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  geom_errorbar(aes (ymin=mean-se, ymax=mean+se), size=1,width=0,alpha=0.45,position=position_dodge(0.3))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  labs(x="Year", y="Shannon Diversity\n", col="Rain trt.", shape="Graze trt.") +
  ylim(2.35,3)+  
  facet_wrap(~graze_trt) +
  ggtitle('Rain*    Year***   Graze:Year**')+
  theme_bw()+
  theme(legend.position = "none", panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.title.y = element_text(size = 10))
fig_div



#'  
#'  * Richness *
#'  

#' Check data distributions
ggplot(data=structure_all_sqrt_check) + 
  geom_point(aes(x=rain_trt,y=rich, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=rich, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  labs(y="Richness\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

#' Specify model
rich_mod_all <- lmer(rich ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=structure_all_sqrt)

# Check residuals
plot(rich_mod_all)
qqnorm(residuals(rich_mod_all))

#; Model summaries
anova(rich_mod_all)
r2(rich_mod_all)
summary(rich_mod_all)
summary(rich_mod_all2)


#' Richness Fig.S5
#' 
structure_sum_rich <- structure_sum_sqrt %>% filter(metric=="Richness")
structure_sum_rich$graze_trt <- factor(structure_sum_rich$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))
structure_sum_rich$rain_trt <- factor(structure_sum_rich$rain_trt, levels=c("D","C","W"))
structure_sum_rich$year <- as.numeric(as.character(structure_sum_rich$year))

fig_rich <- 
  ggplot(data=structure_sum_rich, aes(x=year, y=mean, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.35)) +
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  geom_errorbar(aes (ymin=mean-se, ymax=mean+se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  labs(x="Year", y="Species Richness\n", col="Rain trt.", shape="Graze trt.") +
  ylim(17.5,28)+  
  xlim(2017.5,2020.5) +
  facet_wrap(~graze_trt) +
  ggtitle('Rain*    Year***')+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.title.y = element_text(size = 10))
fig_rich



#'  
#'  * Evenness *
#'  

#' Check data distributions
ggplot(data=structure_all_sqrt_check) + 
  geom_point(aes(x=rain_trt,y=even, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=even, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  labs(y="Evenness\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

#' Specify model
even_mod_all <- lmer(even~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=structure_all_sqrt)

# Check residuals
plot(even_mod_all)
qqnorm(residuals(even_mod_all))

# Model summaries
anova(even_mod_all)
r2(even_mod_all)
summary(even_mod_all)
summary(even_mod_all2)


#' Evenness Fig.S5
#' 
structure_sum_even <- structure_sum_sqrt %>% filter(metric=="Evenness")
structure_sum_even$graze_trt <- factor(structure_sum_even$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing" ))
structure_sum_even$rain_trt <- factor(structure_sum_even$rain_trt, levels=c("D","C","W"))
structure_sum_even$year <- as.numeric(as.character(structure_sum_even$year))

fig_even <- 
  ggplot(data=structure_sum_even, aes(x=year, y=mean, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.3)) +
  geom_errorbar(aes (ymin=mean-se, ymax=mean+se), size=1,width=0,alpha=0.45,position=position_dodge(0.3))+
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  labs(x="Year", y="Species Evenness\n", col="Rain trt.", shape="Graze trt.") +
  ylim(0.2,0.5) + facet_wrap(~graze_trt) +
  ggtitle('Rain**    Year**')+
  theme_bw()+
  theme(legend.position = "none", panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.title.y = element_text(size = 10))
fig_even



#'  
#'  * C4 grasses *
#'  

#' Check data distributions
g_c4_sqrt_check <-g_c4_sqrt
g_c4_sqrt_check$graze_trt <- gsub("DS", "Dormant-Season Grazing", g_c4_sqrt_check$graze_trt)
g_c4_sqrt_check$graze_trt <- gsub("GS", "Growing-Season Grazing", g_c4_sqrt_check$graze_trt)
g_c4_sqrt_check$graze_trt <- gsub("UG", "Ungrazed", g_c4_sqrt_check$graze_trt)
g_c4_sqrt_check$graze_trt <- factor(g_c4_sqrt_check$graze_trt, levels=c("Ungrazed","Dormant-Season Grazing", "Growing-Season Grazing"))
g_c4_sqrt_check$rain_trt <- factor(g_c4_sqrt_check$rain_trt, levels=c("D","C","W"))
g_c4_sqrt_check$year <- factor(g_c4_sqrt_check$year, levels=c("2018","2019","2020"))
ggplot(data=g_c4_sqrt_check) + 
  geom_point(aes(x=rain_trt,y=fg_cov, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=fg_cov, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  labs(y="C4 Perennial Grass\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

#' Specify models with different contrasts
c4_mod_all <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=g_c4_sqrt)
c4_mod_all2 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=g_c4_sqrt_2)
c4_mod_all3 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=g_c4_sqrt_3)
c4_mod_all4 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=g_c4_sqrt_4)

# Check residuals
plot(c4_mod_all)
qqnorm(residuals(c4_mod_all))

# Model summaries
anova(c4_mod_all)
r2(c4_mod_all)
summary(c4_mod_all)

# Save univariate fixed effects from both models (all contrasts)
c4_fix <- data.frame(round(fixef(c4_mod_all),4)) %>% rename(fixef=round.fixef.c4_mod_all...4.)
c4_ci  <- data.frame(round(confint(c4_mod_all),4))
c4_fixef<- merge(c4_fix, c4_ci, by=0, all=TRUE) 
c4_fixef$response <- "c4"
c4_fixef$Row.names <- gsub("rain_trt1", "rain_D", c4_fixef$Row.names)
c4_fixef$Row.names <- gsub("rain_trt2", "rain_W", c4_fixef$Row.names)
c4_fixef$Row.names <- gsub("graze_trt1", "graze_DS", c4_fixef$Row.names)
c4_fixef$Row.names <- gsub("graze_trt2", "graze_GS", c4_fixef$Row.names)
c4_fixef$Row.names <- gsub("year1", "year_2020", c4_fixef$Row.names)
c4_fixef$Row.names <- gsub("year2", "year_2019", c4_fixef$Row.names)
c4_fix2 <- data.frame(round(fixef(c4_mod_all2),4))%>% rename(fixef=round.fixef.c4_mod_all2...4.)
c4_ci2  <- data.frame(round(confint(c4_mod_all2),4))
c4_fixef2<- merge(c4_fix2, c4_ci2, by=0, all=TRUE) 
c4_fixef2$response <- "c4"
c4_fixef2$Row.names <- gsub("rain_trt1", "rain_C", c4_fixef2$Row.names)
c4_fixef2$Row.names <- gsub("rain_trt2", "rain_D", c4_fixef2$Row.names)
c4_fixef2$Row.names <- gsub("graze_trt1", "graze_UG", c4_fixef2$Row.names)
c4_fixef2$Row.names <- gsub("graze_trt2", "graze_DS", c4_fixef2$Row.names)
c4_fixef2$Row.names <- gsub("year1", "year_2018", c4_fixef2$Row.names)
c4_fixef2$Row.names <- gsub("year2", "year_2019", c4_fixef2$Row.names)
c4_fix3 <- data.frame(round(fixef(c4_mod_all3),4))%>% rename(fixef=round.fixef.c4_mod_all3...4.)
c4_ci3  <- data.frame(round(confint(c4_mod_all3),4))
c4_fixef3<- merge(c4_fix3, c4_ci3, by=0, all=TRUE) 
c4_fixef3$response <- "c4"
c4_fixef3$Row.names <- gsub("rain_trt1", "rain_W", c4_fixef3$Row.names)
c4_fixef3$Row.names <- gsub("rain_trt2", "rain_C", c4_fixef3$Row.names)
c4_fixef3$Row.names <- gsub("graze_trt1", "graze_UG", c4_fixef3$Row.names)
c4_fixef3$Row.names <- gsub("graze_trt2", "graze_GS", c4_fixef3$Row.names)
c4_fixef3$Row.names <- gsub("year1", "year_2018", c4_fixef3$Row.names)
c4_fixef3$Row.names <- gsub("year2", "year_2019", c4_fixef3$Row.names)
c4_fix4 <- data.frame(round(fixef(c4_mod_all4),4))%>% rename(fixef=round.fixef.c4_mod_all4...4.)
c4_ci4  <- data.frame(round(confint(c4_mod_all4),4))
c4_fixef4<- merge(c4_fix4, c4_ci4, by=0, all=TRUE) 
c4_fixef4$response <- "c4"
c4_fixef4$Row.names <- gsub("rain_trt1", "rain_C", c4_fixef4$Row.names)
c4_fixef4$Row.names <- gsub("rain_trt2", "rain_D", c4_fixef4$Row.names)
c4_fixef4$Row.names <- gsub("graze_trt1", "graze_UG", c4_fixef4$Row.names)
c4_fixef4$Row.names <- gsub("graze_trt2", "graze_GS", c4_fixef4$Row.names)
c4_fixef4$Row.names <- gsub("year1", "year_2020", c4_fixef4$Row.names)
c4_fixef4$Row.names <- gsub("year2", "year_2019", c4_fixef4$Row.names)

c4_fixef_all <- unique(rbind(c4_fixef, c4_fixef2, c4_fixef3, c4_fixef4)) %>% arrange(Row.names)
c4_fixef_all


#' C4 Grasses - Fig. 4
#' 
fig_c4 <- 
  ggplot(data=c4_sum, aes(x=year, y=G_C4, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.35)) +
  geom_errorbar(aes (ymin=G_C4-se, ymax=G_C4+se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  labs(x="Year", y="C4 Perennial Grasses\n(rel. cover proportion)", col="Rain trt.", shape="Graze trt.") +
  ylim(0.2,0.425) +  
  xlim(2017.5,2020.5) +
  facet_wrap(~graze_trt) +
  ggtitle('Year***')+
  theme_bw()+
  theme(legend.position = "none", panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8),  axis.title.y = element_text(size = 10))
fig_c4



#' 
#' *C3 perennial grasses*
#' 

#' Check data distributions
g_c3_sqrt_check <- g_c3_sqrt
g_c3_sqrt_check$graze_trt <- gsub("DS", "Dormant-Season Grazing", g_c3_sqrt_check$graze_trt)
g_c3_sqrt_check$graze_trt <- gsub("GS", "Growing-Season Grazing", g_c3_sqrt_check$graze_trt)
g_c3_sqrt_check$graze_trt <- gsub("UG", "Ungrazed", g_c3_sqrt_check$graze_trt)
g_c3_sqrt_check$graze_trt <- factor(g_c3_sqrt_check$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing" ))
g_c3_sqrt_check$rain_trt <- factor(g_c3_sqrt_check$rain_trt, levels=c("D","C","W"))
g_c3_sqrt_check$year <- factor(g_c3_sqrt_check$year, levels=c("2018","2019","2020"))
ggplot(data=g_c3_sqrt_check) + 
  geom_point(aes(x=rain_trt,y=fg_cov, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=fg_cov, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  labs(y="C3 Perennial Grass\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

#' Specify models showing different contrasts in summary
c3_mod_all <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=g_c3_sqrt)
c3_mod_all2 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=g_c3_sqrt_2)
c3_mod_all3 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=g_c3_sqrt_3)
c3_mod_all4 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=g_c3_sqrt_4)

# Check residuals
plot(c3_mod_all)
qqnorm(residuals(c3_mod_all))

# Model summaries
anova(c3_mod_all)
r2(c3_mod_all)
summary(c3_mod_all)

# Save univariate fixed effects from both models (all contrasts)
c3_fix <- data.frame(round(fixef(c3_mod_all),4)) %>% rename(fixef=round.fixef.c3_mod_all...4.)
c3_ci  <- data.frame(round(confint(c3_mod_all),4))
c3_fixef<- merge(c3_fix, c3_ci, by=0, all=TRUE) 
c3_fixef$response <- "c3"
c3_fixef$Row.names <- gsub("rain_trt1", "rain_D", c3_fixef$Row.names)
c3_fixef$Row.names <- gsub("rain_trt2", "rain_W", c3_fixef$Row.names)
c3_fixef$Row.names <- gsub("graze_trt1", "graze_DS", c3_fixef$Row.names)
c3_fixef$Row.names <- gsub("graze_trt2", "graze_GS", c3_fixef$Row.names)
c3_fixef$Row.names <- gsub("year1", "year_2020", c3_fixef$Row.names)
c3_fixef$Row.names <- gsub("year2", "year_2019", c3_fixef$Row.names)
c3_fix2 <- data.frame(round(fixef(c3_mod_all2),4))%>% rename(fixef=round.fixef.c3_mod_all2...4.)
c3_ci2  <- data.frame(round(confint(c3_mod_all2),4))
c3_fixef2<- merge(c3_fix2, c3_ci2, by=0, all=TRUE) 
c3_fixef2$response <- "c3"
c3_fixef2$Row.names <- gsub("rain_trt1", "rain_C", c3_fixef2$Row.names)
c3_fixef2$Row.names <- gsub("rain_trt2", "rain_D", c3_fixef2$Row.names)
c3_fixef2$Row.names <- gsub("graze_trt1", "graze_UG", c3_fixef2$Row.names)
c3_fixef2$Row.names <- gsub("graze_trt2", "graze_DS", c3_fixef2$Row.names)
c3_fixef2$Row.names <- gsub("year1", "year_2018", c3_fixef2$Row.names)
c3_fixef2$Row.names <- gsub("year2", "year_2019", c3_fixef2$Row.names)
c3_fix3 <- data.frame(round(fixef(c3_mod_all3),4))%>% rename(fixef=round.fixef.c3_mod_all3...4.)
c3_ci3  <- data.frame(round(confint(c3_mod_all3),4))
c3_fixef3<- merge(c3_fix3, c3_ci3, by=0, all=TRUE) 
c3_fixef3$response <- "c3"
c3_fixef3$Row.names <- gsub("rain_trt1", "rain_W", c3_fixef3$Row.names)
c3_fixef3$Row.names <- gsub("rain_trt2", "rain_C", c3_fixef3$Row.names)
c3_fixef3$Row.names <- gsub("graze_trt1", "graze_UG", c3_fixef3$Row.names)
c3_fixef3$Row.names <- gsub("graze_trt2", "graze_GS", c3_fixef3$Row.names)
c3_fixef3$Row.names <- gsub("year1", "year_2018", c3_fixef3$Row.names)
c3_fixef3$Row.names <- gsub("year2", "year_2019", c3_fixef3$Row.names)
c3_fix4 <- data.frame(round(fixef(c3_mod_all4),4))%>% rename(fixef=round.fixef.c3_mod_all4...4.)
c3_ci4  <- data.frame(round(confint(c3_mod_all4),4))
c3_fixef4<- merge(c3_fix4, c3_ci4, by=0, all=TRUE) 
c3_fixef4$response <- "c3"
c3_fixef4$Row.names <- gsub("rain_trt1", "rain_C", c3_fixef4$Row.names)
c3_fixef4$Row.names <- gsub("rain_trt2", "rain_D", c3_fixef4$Row.names)
c3_fixef4$Row.names <- gsub("graze_trt1", "graze_UG", c3_fixef4$Row.names)
c3_fixef4$Row.names <- gsub("graze_trt2", "graze_GS", c3_fixef4$Row.names)
c3_fixef4$Row.names <- gsub("year1", "year_2020", c3_fixef4$Row.names)
c3_fixef4$Row.names <- gsub("year2", "year_2019", c3_fixef4$Row.names)

c3_fixef_all <- unique(rbind(c3_fixef, c3_fixef2, c3_fixef3, c3_fixef4)) %>% arrange(Row.names)
c3_fixef_all 


#' C3 Grasses - Fig. 4
#' 

fig_c3 <- 
  ggplot(data=c3_sum, aes(x=year, y=G_C3, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.35)) +
  geom_errorbar(aes (ymin=G_C3-se, ymax=G_C3+se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  labs(x="Year", y="C3 Perennial Grasses\n(rel. cover proportion)", col="Rain trt.", shape="Graze trt.") +
  ylim(0.2,0.45) + 
  xlim(2017.5,2020.5) +
  facet_wrap(~graze_trt) +
  ggtitle('Year**     Graze:Year**')+
  theme_bw()+
  theme(legend.position = "none", 
        panel.grid.minor.x = element_blank(), 
        plot.title = element_text(hjust = 0.5, size=8), 
        axis.title.y = element_text(size = 10))
fig_c3



#' 
#' *Annual-Biennials*
#' 

#' Check data distributions
#' 
#'    There is at least one outlier here to check.
#'    
ann_sqrt_check <- ann_sqrt
ann_sqrt_check$graze_trt <- gsub("DS", "Dormant-Season Grazing", ann_sqrt_check$graze_trt)
ann_sqrt_check$graze_trt <- gsub("GS", "Growing-Season Grazing", ann_sqrt_check$graze_trt)
ann_sqrt_check$graze_trt <- gsub("UG", "Ungrazed", ann_sqrt_check$graze_trt)
ann_sqrt_check$graze_trt <- factor(ann_sqrt_check$graze_trt, levels=c( "Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))
ann_sqrt_check$rain_trt <- factor(ann_sqrt_check$rain_trt, levels=c("D","C","W"))
ann_sqrt_check$year <- factor(ann_sqrt_check$year, levels=c("2018","2019","2020"))
ggplot(data=ann_sqrt_check) + 
  geom_point(aes(x=rain_trt,y=fg_cov, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=fg_cov, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  labs(y="Annual-Biennials\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')


#'  Based on qqnorm plot from the full model, there's one outlier that creates
#'  skew in the residuals, 3-4. 
#'  
#'  Create outlier dataframe and check against full model 
ann_sqrt_out <- ann_sqrt %>% filter (b_p != "3_4") 
#'
#'  Outlier model
ann_mod_all_out  <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|b_p), data=ann_sqrt_out)
#'  
#'  Full model
ann_mod_all <- lmer(fg_cov ~ rain_trt * graze_trt * year +  (1|b_p), data=ann_sqrt)
#' 
#' Compare outlier and full model results
anova(ann_mod_all)
anova(ann_mod_all_out)
#'
#' Compare fitted vs residuals
plot(ann_mod_all)
plot(ann_mod_all_out)
#'
#' Compare QQ plot of residuals
qqnorm(residuals(ann_mod_all))
qqnorm(residuals(ann_mod_all_out))


#' Conclusions: Remove outlier 3-4
#' 
#' Removing this outlier improves model diagnostics, though only slightly 
#' modifies model output.


#' Specify models
#' 
#' Create new dataframes that remove the outlier
ann_sqrt_out <- ann_sqrt %>% filter (b_p != "3_4") 
ann_sqrt_2_out <- ann_sqrt_2 %>% filter (b_p != "3_4") 
ann_sqrt_3_out <- ann_sqrt_3 %>% filter (b_p != "3_4") 
ann_sqrt_4_out <- ann_sqrt_4 %>% filter (b_p != "3_4") 

#' Specify full models (Identical models but showing different contrasts in summary)
#'     Note: Removed the random effect of block to singularity
ann_mod_all <- lmer(fg_cov ~ rain_trt * graze_trt * year +  (1|b_p), data=ann_sqrt_out)
ann_mod_all2 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|b_p), data=ann_sqrt_2_out)
ann_mod_all3 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|b_p), data=ann_sqrt_3_out)
ann_mod_all4 <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|b_p), data=ann_sqrt_4_out)

# Model summaries
Anova(ann_mod_all)
anova(ann_mod_all)
r2(ann_mod_all)
summary(ann_mod_all)

# Save univariate fixed effects from both models (all contrasts)
ann_fix <- data.frame(round(fixef(ann_mod_all),4)) %>% rename(fixef=round.fixef.ann_mod_all...4.)
ann_ci  <- data.frame(round(confint(ann_mod_all),4))
ann_fixef<- merge(ann_fix, ann_ci, by=0, all=TRUE) 
ann_fixef$response <- "ann"
ann_fixef$Row.names <- gsub("rain_trt1", "rain_D", ann_fixef$Row.names)
ann_fixef$Row.names <- gsub("rain_trt2", "rain_W", ann_fixef$Row.names)
ann_fixef$Row.names <- gsub("graze_trt1", "graze_DS", ann_fixef$Row.names)
ann_fixef$Row.names <- gsub("graze_trt2", "graze_GS", ann_fixef$Row.names)
ann_fixef$Row.names <- gsub("year1", "year_2020", ann_fixef$Row.names)
ann_fixef$Row.names <- gsub("year2", "year_2019", ann_fixef$Row.names)
ann_fix2 <- data.frame(round(fixef(ann_mod_all2),4))%>% rename(fixef=round.fixef.ann_mod_all2...4.)
ann_ci2  <- data.frame(round(confint(ann_mod_all2),4))
ann_fixef2<- merge(ann_fix2, ann_ci2, by=0, all=TRUE) 
ann_fixef2$response <- "ann"
ann_fixef2$Row.names <- gsub("rain_trt1", "rain_C", ann_fixef2$Row.names)
ann_fixef2$Row.names <- gsub("rain_trt2", "rain_D", ann_fixef2$Row.names)
ann_fixef2$Row.names <- gsub("graze_trt1", "graze_UG", ann_fixef2$Row.names)
ann_fixef2$Row.names <- gsub("graze_trt2", "graze_DS", ann_fixef2$Row.names)
ann_fixef2$Row.names <- gsub("year1", "year_2018", ann_fixef2$Row.names)
ann_fixef2$Row.names <- gsub("year2", "year_2019", ann_fixef2$Row.names)
ann_fix3 <- data.frame(round(fixef(ann_mod_all3),4))%>% rename(fixef=round.fixef.ann_mod_all3...4.)
ann_ci3  <- data.frame(round(confint(ann_mod_all3),4))
ann_fixef3<- merge(ann_fix3, ann_ci3, by=0, all=TRUE) 
ann_fixef3$response <- "ann"
ann_fixef3$Row.names <- gsub("rain_trt1", "rain_W", ann_fixef3$Row.names)
ann_fixef3$Row.names <- gsub("rain_trt2", "rain_C", ann_fixef3$Row.names)
ann_fixef3$Row.names <- gsub("graze_trt1", "graze_UG", ann_fixef3$Row.names)
ann_fixef3$Row.names <- gsub("graze_trt2", "graze_GS", ann_fixef3$Row.names)
ann_fixef3$Row.names <- gsub("year1", "year_2018", ann_fixef3$Row.names)
ann_fixef3$Row.names <- gsub("year2", "year_2019", ann_fixef3$Row.names)
ann_fix4 <- data.frame(round(fixef(ann_mod_all4),4))%>% rename(fixef=round.fixef.ann_mod_all4...4.)
ann_ci4  <- data.frame(round(confint(ann_mod_all4),4))
ann_fixef4<- merge(ann_fix4, ann_ci4, by=0, all=TRUE) 
ann_fixef4$response <- "ann"
ann_fixef4$Row.names <- gsub("rain_trt1", "rain_C", ann_fixef4$Row.names)
ann_fixef4$Row.names <- gsub("rain_trt2", "rain_D", ann_fixef4$Row.names)
ann_fixef4$Row.names <- gsub("graze_trt1", "graze_UG", ann_fixef4$Row.names)
ann_fixef4$Row.names <- gsub("graze_trt2", "graze_GS", ann_fixef4$Row.names)
ann_fixef4$Row.names <- gsub("year1", "year_2020", ann_fixef4$Row.names)
ann_fixef4$Row.names <- gsub("year2", "year_2019", ann_fixef4$Row.names)

ann_fixef_all <- unique(rbind(ann_fixef, ann_fixef2, ann_fixef3, ann_fixef4)) %>% arrange(Row.names)
ann_fixef_all 


#' Annual-Biennials - Fig. 4
#' 

#' First, remove outlier 3_4 from the figure dataframe 
ann_sum_out <- spp_df_rel_sqrt_long_fg %>%
  ungroup() %>%
  filter(b_p != "3_4") %>%
  select(-b_p, -block) %>%
  filter(fun_grp == "Ann_Bi") %>%
  group_by(year, rain_trt, graze_trt) %>%
  summarise( Ann_Bi = mean(fg_cov),
             se= se(fg_cov))
ann_sum_out$year <- as.numeric(as.character(ann_sum_out$year))
ann_sum_out$graze_trt <- gsub("UG", "Ungrazed", ann_sum_out$graze_trt)
ann_sum_out$graze_trt <- gsub("DS", "Dormant-Season Grazing", ann_sum_out$graze_trt)
ann_sum_out$graze_trt <- gsub("GS", "Growing-Season Grazing", ann_sum_out$graze_trt)
ann_sum_out$graze_trt <- factor(ann_sum_out$graze_trt, levels=c("Ungrazed",  "Growing-Season Grazing","Dormant-Season Grazing"))
ann_sum_out$rain_trt <- factor(ann_sum_out$rain_trt, levels=c("D","C","W"))
ann_sum_out <- ann_sum_out %>% arrange(graze_trt, year)

#' Create the figure
fig_ann <- 
  ggplot(data=ann_sum_out, aes(x=year, y=Ann_Bi, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.35)) +
  geom_errorbar(aes (ymin=Ann_Bi-se, ymax=Ann_Bi+se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  labs(x="Year", y="Annual-Biennials\n(rel. cover proportion)", col="Rain trt.", shape="Graze trt.") +
  ylim(0,0.125) +  
  xlim(2017.5,2020.5) +
  facet_wrap(~graze_trt) +
  ggtitle('Rain**     Graze**    Year***     Rain:Year***')+
  theme_bw()+
  theme(legend.position = "none", 
        panel.grid.minor.x = element_blank(), 
        plot.title = element_text(hjust = 0.5, size=8), 
        axis.title.y = element_text(size = 10))
fig_ann




#' 
#' *Forbs*
#' 

#' Check data distributions
f_sqrt_check <- f_sqrt
f_sqrt_check$graze_trt <- gsub("FG", "Dormant-Season Grazing", f_sqrt_check$graze_trt)
f_sqrt_check$graze_trt <- gsub("SG", "Growing-Season Grazing", f_sqrt_check$graze_trt)
f_sqrt_check$graze_trt <- gsub("C", "Ungrazed", f_sqrt_check$graze_trt)
f_sqrt_check$graze_trt <- factor(f_sqrt_check$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing" ))
f_sqrt_check$rain_trt <- factor(f_sqrt_check$rain_trt, levels=c("D","C","W"))
f_sqrt_check$year <- factor(f_sqrt_check$year, levels=c("2018","2019","2020"))
ggplot(data=f_sqrt_check) + 
  geom_point(aes(x=rain_trt,y=fg_cov, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=fg_cov, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  labs(y="Perennial Forbs\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

#' Specify full models (Identical models but showing different contrasts in summary)
f_mod_all <- lmer(fg_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=f_sqrt)
f_mod_all2 <- lmer(fg_cov ~ rain_trt * graze_trt * year+ (1|block) + (1|b_p), data=f_sqrt_2)
f_mod_all3 <- lmer(fg_cov ~ rain_trt * graze_trt * year+ (1|block) + (1|b_p), data=f_sqrt_3)
f_mod_all4 <- lmer(fg_cov ~ rain_trt * graze_trt * year+ (1|block) + (1|b_p), data=f_sqrt_4)

# Check residuals
plot(f_mod_all)
qqnorm(residuals(f_mod_all))

# Model summaries
Anova(f_mod_all)
anova(f_mod_all)
r2(f_mod_all)
summary(f_mod_all)
summary(f_mod_all2)

# Save univariate fixed effects from both models (all contrasts)
f_fix <- data.frame(round(fixef(f_mod_all),4)) %>% rename(fixef=round.fixef.f_mod_all...4.)
f_ci  <- data.frame(round(confint(f_mod_all),4))
f_fixef<- merge(f_fix, f_ci, by=0, all=TRUE) 
f_fixef$response <- "f"
f_fixef$Row.names <- gsub("rain_trt1", "rain_D", f_fixef$Row.names)
f_fixef$Row.names <- gsub("rain_trt2", "rain_W", f_fixef$Row.names)
f_fixef$Row.names <- gsub("graze_trt1", "graze_DS", f_fixef$Row.names)
f_fixef$Row.names <- gsub("graze_trt2", "graze_GS", f_fixef$Row.names)
f_fixef$Row.names <- gsub("year1", "year_2020", f_fixef$Row.names)
f_fixef$Row.names <- gsub("year2", "year_2019", f_fixef$Row.names)
f_fix2 <- data.frame(round(fixef(f_mod_all2),4))%>% rename(fixef=round.fixef.f_mod_all2...4.)
f_ci2  <- data.frame(round(confint(f_mod_all2),4))
f_fixef2<- merge(f_fix2, f_ci2, by=0, all=TRUE) 
f_fixef2$response <- "f"
f_fixef2$Row.names <- gsub("rain_trt1", "rain_C", f_fixef2$Row.names)
f_fixef2$Row.names <- gsub("rain_trt2", "rain_D", f_fixef2$Row.names)
f_fixef2$Row.names <- gsub("graze_trt1", "graze_UG", f_fixef2$Row.names)
f_fixef2$Row.names <- gsub("graze_trt2", "graze_DS", f_fixef2$Row.names)
f_fixef2$Row.names <- gsub("year1", "year_2018", f_fixef2$Row.names)
f_fixef2$Row.names <- gsub("year2", "year_2019", f_fixef2$Row.names)
f_fix3 <- data.frame(round(fixef(f_mod_all3),4))%>% rename(fixef=round.fixef.f_mod_all3...4.)
f_ci3  <- data.frame(round(confint(f_mod_all3),4))
f_fixef3<- merge(f_fix3, f_ci3, by=0, all=TRUE) 
f_fixef3$response <- "f"
f_fixef3$Row.names <- gsub("rain_trt1", "rain_W", f_fixef3$Row.names)
f_fixef3$Row.names <- gsub("rain_trt2", "rain_C", f_fixef3$Row.names)
f_fixef3$Row.names <- gsub("graze_trt1", "graze_UG", f_fixef3$Row.names)
f_fixef3$Row.names <- gsub("graze_trt2", "graze_GS", f_fixef3$Row.names)
f_fixef3$Row.names <- gsub("year1", "year_2018", f_fixef3$Row.names)
f_fixef3$Row.names <- gsub("year2", "year_2019", f_fixef3$Row.names)
f_fix4 <- data.frame(round(fixef(f_mod_all4),4))%>% rename(fixef=round.fixef.f_mod_all4...4.)
f_ci4  <- data.frame(round(confint(f_mod_all4),4))
f_fixef4<- merge(f_fix4, f_ci4, by=0, all=TRUE) 
f_fixef4$response <- "f"
f_fixef4$Row.names <- gsub("rain_trt1", "rain_C", f_fixef4$Row.names)
f_fixef4$Row.names <- gsub("rain_trt2", "rain_D", f_fixef4$Row.names)
f_fixef4$Row.names <- gsub("graze_trt1", "graze_UG", f_fixef4$Row.names)
f_fixef4$Row.names <- gsub("graze_trt2", "graze_GS", f_fixef4$Row.names)
f_fixef4$Row.names <- gsub("year1", "year_2020", f_fixef4$Row.names)
f_fixef4$Row.names <- gsub("year2", "year_2019", f_fixef4$Row.names)

f_fixef_all <- unique(rbind(f_fixef, f_fixef2, f_fixef3, f_fixef4)) %>% arrange(Row.names)
f_fixef_all 


#' Forbs - Fig. 4
#' 
fig_forb <- 
  ggplot(data=f_sum, aes(x=year, y=Forb, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.35)) +
  geom_errorbar(aes (ymin=Forb-se, ymax=Forb+se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  labs(x="Year", y="Perennial Forbs\n(rel. cover proportion)", col="Rain trt.", shape="Graze trt.") +
  ylim(0.2,0.41) +  
  xlim(2017.5,2020.5) +
  facet_wrap(~graze_trt) +
  ggtitle('Year***     Graze:Year*')+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank(), 
        plot.title = element_text(hjust = 0.5, size=8), 
        axis.title.y = element_text(size = 10))
fig_forb







#' 
#' ** Supporting data figure: Absolute relative covers for each functional groups**
#' 

#' Review dataframe
str(fg_abs_full_allyears)

#' Remove annual-biennial values for one outlier plot (3-4, in row 31, 103 and )
fg_abs_full_allyears_out <- fg_abs_full_allyears
outliers <- fg_abs_full_allyears_out$b_p == "3_4"
fg_abs_full_allyears_out$ann_abs[outliers] <- NA

#' Estimate sums
fg_abs_sum_dat <- fg_abs_full_allyears_out %>% gather(c4_abs:f_abs, key=fg, value=abs_cov) %>%
  group_by(fg, year, graze_trt, rain_trt) %>%
  summarise (abs_cov_mean= mean(abs_cov, na.rm=T),
             abs_cov_se= se(abs_cov, na.rm=T))
fg_abs_sum_dat$graze_trt <- factor(fg_abs_sum_dat$graze_trt, levels=c("GS","DS","UG"))
fg_abs_sum_dat$rain_trt <- factor(fg_abs_sum_dat$rain_trt, levels=c("D","C","W"))
fg_abs_sum_dat$year <- as.numeric(as.character(fg_abs_sum_dat$year))
fg_abs_sum_dat$graze_trt <- gsub("DS", "Dormant-Season Grazing", fg_abs_sum_dat$graze_trt)
fg_abs_sum_dat$graze_trt <- gsub("GS", "Growing-Season Grazing", fg_abs_sum_dat$graze_trt)
fg_abs_sum_dat$graze_trt <- gsub("UG", "Ungrazed", fg_abs_sum_dat$graze_trt)
fg_abs_sum_dat$fg <- gsub("ann_abs", "Annual-Biennials", fg_abs_sum_dat$fg)
fg_abs_sum_dat$fg <- gsub("c3_abs", "C3 Perennial Grasses", fg_abs_sum_dat$fg)
fg_abs_sum_dat$fg <- gsub("c4_abs", "C4 Perennial Grasses", fg_abs_sum_dat$fg)
fg_abs_sum_dat$fg <- gsub("f_abs", "Perennial Forbs", fg_abs_sum_dat$fg)
fg_abs_sum_dat$fg <- factor(fg_abs_sum_dat$fg, levels=c("Annual-Biennials","Perennial Forbs","C3 Perennial Grasses", "C4 Perennial Grasses"))
fg_abs_sum_dat$graze_trt <- factor(fg_abs_sum_dat$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))

fig_fg_abs <- 
  ggplot(data=fg_abs_sum_dat, aes(x=year, y=abs_cov_mean, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.35)) +
  #geom_line(position=position_dodge(0.2),alpha=0.4) +
  geom_errorbar(aes (ymin=abs_cov_mean-abs_cov_se, ymax=abs_cov_mean+abs_cov_se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  #geom_ribbon(alpha=0.2, aes(ymin=abs_cov_mean-abs_cov_se, ymax=abs_cov_mean+abs_cov_se, fill=rain_trt)) +
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
 # geom_vline(aes(xintercept=2018),lty=2, size=0.65,col="gray70")+
  #geom_text(aes(x=2018.1, y=0.06, label="baseline"),angle=90, col="gray60",cex=3)+ 
  labs(x="Year", y="Absolute Cover (%)", col="Rain trt.", shape="Graze trt.") +
  facet_grid(fg ~ graze_trt, scale='free') +
  theme_bw()+
  theme( panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8))
fig_fg_abs







#####
#'    
#'    
#'    
#'    [Part 2] : *FORAGE SERVICES RESPONSES*
#'                   to rainfall and grazing interactions
#'  
#'  
#'  [A] : Compile raw dataframes
#'  
#'  [B] : Prep model and figure dataframes
#'  
#'  [C] : Explore rainfall and grazing interactions on forage response
#'  
#'  [D] : Identify key indices of plant-driven functions and services
#'  
#'  
#'  
#'###         





##
#' 
#' 
#' [A]: **RAW DATA PREPARATION**
#' 
#'                Compile dataframes to assess rainfall and grazing effects on
#'                     plant communities and services
#' 
#'
#'##


#'
#'  **PLANT BIOMASS and BARE GROUND - Prep data**
#'     

#' 
#' *Data preparation*
#' 

#' Load dataframe with plant biomass (production), bare ground, and counts for some annual species
plant_path <-  paste0(gdrive, "Plant_Plot_Metrics_2018-23.csv")
plant <- read.csv(plant_path)

#' View structure of data and remove sample_ID column
str(plant)

#' Subset dataframe to replicate plots only (no exterior) and relevant variables
biomass <- plant %>%
  filter(exterior == "N") %>%
  select(year, block, plot, rain_trt, graze_trt, bare_visual_cov, biomass_g)

#' Convert biomass from [g per 0.05m2] to [pounds per acre] as follows:
#'            g per 0.05m2 * 20 [0.05m2 to 1m2] * .002205 [lb per g] * 4046.86 [1m2 to acre]
#' 
#biomass$biomass_converted <- biomass$biomass_g*20*.002205*4046.86  

#' Convert biomass from [g per 0.05m2] to [g per 0.05m2] as follows:
biomass$biomass_converted <- biomass$biomass_g*20

#' Create exterior dataframe
exterior <- plant %>%
  filter(exterior == "Y") %>%
  select(year, block, plot, plot_size_m2, biomass_g)
# Convert biomass to pounds per acre (see code block above)
exterior$biomass_converted <- exterior$biomass_g*20

#' Set rainfall and grazing treatment levels
biomass$rain_trt <- factor(biomass$rain_trt , levels = c("D","W","C"))  # See D & W as rain_trt 1 and 2 relative to grand mean
biomass$graze_trt <- factor(biomass$graze_trt , levels = c("FG","SG","C")) # See FG & SG as graze trt 1 and 2 relative to grand mean
contrasts(biomass$rain_trt)
contrasts(biomass$graze_trt)

#' Create numeric year variable
biomass$year_num <- biomass$year

#' Create a block-plot grouping variable to join datasets later
biomass$b_p <- paste(biomass$block, biomass$plot, sep="_")
str(biomass)



#'
#'  **FORAGE QUALITY METRICS - Prep data*
#'    (see full section below for full details)
#'    

#' Forage Quality data 
forage_path <-  paste0(gdrive, "ForageQuality_2020.csv")
forage <- read.csv(forage_path)

#' View structure of data; make block a factor and remove sample_ID column
str(forage)
forage$block <- as.factor(forage$block)
forage <- forage %>%
  select(-sample_id)

# Separate sample plots from known functional group samples
forage_plots <- forage %>%
  filter(type=="all_veg")

#' *Explore data and remove highly correlated variables* to subset to smaller number of responses
forage_qual <- forage_plots %>%
  select(crude_protein:total_n) %>%
  ggpairs(., labeller = labeller(facet_category = label_wrap_gen(width = 8)),
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.3)), 
          upper = list(continuous = wrap("cor", size = 3))) +
  theme_grey(base_size = 8)
forage_qual 

#' Save scatterplot 
tiff(filename="forage_quality_correlations.tiff", res=600, width=9, height = 9, units = "in")
forage_qual
dev.off()

#' Save table of correlations
library(Hmisc)
forage_qual_corr <- forage_plots %>%
  select(crude_protein:total_n) 
forage_qual_table <- rcorr(as.matrix(forage_qual_corr))
forage_qual_table$r
write.csv(forage_qual_table$r, 'forage_quality_table.csv')


#'  Wow, super highly correlated! We will reduce to only include:
#'    crude protein (~total_n), 
#'    relative feed value (~neut detergent fiber), and 
#'    digestible_energy (all others) 

forage_plots_sub <- forage_plots %>%
  select(block, plot, rain_trt, graze_trt, crude_protein, digestible_energy, relative_feed_value)
forage_plots_sub %>%
  select(crude_protein, digestible_energy, relative_feed_value)%>%
  pairs()
str(forage_plots_sub)

#' Add block/plot and year joining variables
forage_plots_sub$b_p <- paste(forage_plots_sub$block, forage_plots_sub$plot, sep="_")
forage_plots_sub$year <- as.factor(2020)



#'
#'  *PHENOLOGY METRICS - Prep data*
#'  
pheno_mets_path <- paste0(gdrive, "GreennessCurveMetrics_compiled.csv")
pheno_mets <- read.csv(pheno_mets_path)


#' Store a general treatment and grouping variables dataframe (72 rows x 4 columns)
trt_grp <- sb %>%
  filter(year == "2019", season == "S") %>%
  select(block, plot, rain_trt, graze_trt)


str(pheno_mets)
str(trt_grp)
trt_grp$block <- as.factor(trt_grp$block)
trt_grp$plot <- as.factor(trt_grp$plot)

#' Filter down to only metrics estimated using a global greenness threshold of gi.av=0.3405 
# See phenology section below for more info
pheno_mets_global <- pheno_mets %>% filter( thresh == "gi3405") %>% select(-thresh, -thresh_type) %>%
  merge(trt_grp, by=c("block", "plot"))
str(pheno_mets_global)

#' Add common joining variable
pheno_mets_global$b_p <- paste(pheno_mets_global$block, pheno_mets_global$plot, sep="_")


#' Select only necessary variables (unique variables; see full phenology section below for details)
#' Originally was going to keep total greenness (gs_integral) and growing season length (gs_length)
#'   but note that they are highly correlated *r=0.795*
#' Since growing season length is a more unique phenological variable than total greenness (which
#' could be simply highly correlated to biomass production), let's keep
#' growing season length and date of peak greenness
pheno_mets_global_add <- pheno_mets_global %>% select(b_p, max_gi_doy, gs_length)
pheno_mets_global_add$year  <- as.factor(2020)

#' Looking at the data, several plots had greenness values so low that they never
#' exceeded the global threshold, so greenup, senecesence and growing season length
#' could not be estimated: 
#'   2-10 and 3-4
#'   
#'   With greenness values so slow, the estimated day of peak greenness is also likely
#'   to be uncertain, so we will replace max_gi_doy with 'NA' for those plots as well

pheno_outliers <- pheno_mets_global_add$b_p %in% c("2_10","3_4")
pheno_mets_global_add$max_gi_doy[pheno_outliers] <- NA

#' Add to outliers list
outliers_list[4,1] <- '3_4'
outliers_list[4,2] <- 'max_gi_doy'
outliers_list[4,3] <- 'greenness too low to estimate peak green'
outliers_list[5,1] <- '2_10'
outliers_list[5,2] <- 'max_gi_doy'
outliers_list[5,3] <- 'greenness too low to estimate peak green'

#'
#'  *JOIN DATAFRAMES:  biomass/bare, forage quality, phenology, and total cover *
#'   

#'  Review dataframes

#' 3 years of data (2018-2020, 216 rows)
#' 
#' Biomass & Bare ground
str(biomass) 
biomass$year <- as.factor(biomass$year)
biomass$block <- as.factor(biomass$block)

 
#' 1 year of data (2020 only, 72 rows max)
#' 
#' Forage Quality variables (crude protein, digestible energy, relative feed value)
str(forage_plots_sub) 
#' Vegetation phenology variables (day of year for peak greenness, growing season length)
str(pheno_mets_global_add) 


#' Select only necessary variables
biomass_join <- biomass %>% select(year, block, b_p, rain_trt, graze_trt, biomass_converted, bare_visual_cov)
forage_join <- forage_plots_sub %>% select(year, b_p, crude_protein, digestible_energy, relative_feed_value )


#' Join dataframes
eco_vars_pre <- left_join(biomass_join, forage_join, by=c("b_p","year"))
eco_vars <- left_join(eco_vars_pre, pheno_mets_global_add, by=c("b_p","year")) %>%
  filter( year %in% set_years)

#' Take a first look at correlations in 2020
eco_vars %>% filter(year=="2020") %>% ungroup() %>% select(biomass_converted:gs_length) %>%  
  ggpairs(., lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.3))) +
  theme_grey(base_size = 8)


#' Save final correlation figure
tiff(filename="forage_service_correlations_2020.tiff", res=600, width=7, height = 7, units = "in")
eco_vars %>% filter(year=="2020") %>% ungroup() %>% select(biomass_converted, bare_visual_cov, crude_protein:gs_length) %>%  
  ggpairs(., lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.3))) +
  theme_grey(base_size = 8)
dev.off()






#'###
#'  
#'  
#'  
#' [B] : ** Prep model and figure dataframes **
#'  
#'  for variables in 2020 only (1 yr of data) and 2018-2020 (3 yrs of data)
#'  
#'  
#'### 

#'  Revisit structure of eco-variables dataframe
str(eco_vars)


#' Add plant diversity as forage metric
str(structure_raw_sqrt, give.attr=F)
eco_vars_diversity <- structure_raw_sqrt %>% select(b_p,year, shan)
eco_vars <- left_join(eco_vars, eco_vars_diversity, by=c("b_p","year"))

#'  For 2020 only variables, a 2020 only dataframe:
eco_vars_2020 <- eco_vars %>% filter(year=="2020")

#' Create replicated dataframes to modify for any outliers
eco_vars_out <- eco_vars
eco_vars_out_2020 <- eco_vars_2020

#'  For multi-year variables, create a 2019-2020 dataframe with 2018 covariate column:
# Get 2018 covariate columns
covariates_2018 <- eco_vars %>% filter(year=="2018") %>% select(b_p, biomass_converted, bare_visual_cov) %>%
  rename(biomass_2018 = biomass_converted, bare_2018=bare_visual_cov)
# Get 2019-2020 only dataframe
eco_vars_2019_20_prep <- eco_vars %>% filter(year!="2018")
# Join the 2019-2020 dataframe with the 2018 covariate data
eco_vars_2019_20  <- left_join(eco_vars_2019_20_prep, covariates_2018, by="b_p")
#' View dataframes
str(eco_vars_2020)    # 72 rows for 1 year
str(eco_vars_2019_20) # 144 rows for 2 years


#' Scale each dataframe for models
eco_vars_s <- eco_vars  %>%  mutate_if(is.numeric, scale)
eco_vars_2019_20_s <- eco_vars_2019_20  %>%  mutate_if(is.numeric, scale)
eco_vars_2020_s <- eco_vars_2020 %>%  mutate_if(is.numeric, scale)



#' In each model dataframe, set contrasts:
#' 
#'  Note that we use deviant coding that sums to 0, meaning that any level effects are 
#'  tested and reported in comparison to the grand mean (not in comparison to the 
#'  'control' treatment). This decision was made because -- with respect to grazing in particular --
#'  there is not necessarily a clear 'control'. One could consider spring 'active season' grazing to 
#'  be a control as a normal grazing regime in this system, historically. Alternatively, no grazing 
#'  could be considered as a more traditional control, since it excludes the variable of interest (grazing) -- 
#'  even though this marks a CHANGE from the historic grazing regime.
#'  
#'  We maintain a more neutral stance by coding contrasts to compare each treatment to the grand mean.
#'  
#'  However, the way in which contrast levels are ordered will alter which treatmens are shown in the output. 
#'  Here, we set contrasts to show D and W levels for rainfall; SG and FG levels for grazing; 2020 level for year;  

#' Set contrasts 2020 only
eco_vars_2020_s$rain_trt <- factor(eco_vars_2020_s$rain_trt, levels=c("D","W","C"))
eco_vars_2020_s$graze_trt <- factor(eco_vars_2020_s$graze_trt, levels=c("FG","SG","C"))
eco_vars_2020_s2 <- eco_vars_2020_s
eco_vars_2020_s2$rain_trt <- factor(eco_vars_2020_s2$rain_trt, levels=c("C","D","W"))
eco_vars_2020_s2$graze_trt <- factor(eco_vars_2020_s2$graze_trt, levels=c("C","SG","FG"))
eco_vars_2020_s3 <- eco_vars_2020_s
eco_vars_2020_s3$rain_trt <- factor(eco_vars_2020_s3$rain_trt, levels=c("W","C","D"))
eco_vars_2020_s3$graze_trt <- factor(eco_vars_2020_s3$graze_trt, levels=c("C","SG","FG"))

#' Create some duplicates for outlier removal if needed
eco_vars_2020_s_out <- eco_vars_2020_s
eco_vars_2020_s2_out <- eco_vars_2020_s2
eco_vars_2020_s3_out <- eco_vars_2020_s3


#' *Set contrasts for model dataframes*
#' 
#' First set of contrasts
eco_vars_s$rain_trt <- factor(eco_vars_s$rain_trt, levels=c("D","W","C"))
eco_vars_s$graze_trt <- factor(eco_vars_s$graze_trt, levels=c("FG","SG","C"))
eco_vars_s$year <- factor(eco_vars_s$year, levels=c("2020","2019","2018"))
eco_vars_s$b_p <- as.factor(eco_vars_s$b_p)

#' Create a second set of different contrasts 
eco_vars_s2<- eco_vars_s 
eco_vars_s2$graze_trt <- factor(eco_vars_s2$graze_trt, levels=c("C","SG","FG"))
eco_vars_s2$rain_trt <- factor(eco_vars_s2$rain_trt, levels=c("C","D","W"))
eco_vars_s2$year <- factor(eco_vars_s2$year, levels=c("2018","2019","2020"))
eco_vars_s2$b_p <- as.factor(eco_vars_s2$b_p)

#' Create third set of set of contrasts
eco_vars_s3 <- eco_vars_s
eco_vars_s3$graze_trt <- factor(eco_vars_s3$graze_trt, levels=c("C","SG","FG")) 
eco_vars_s3$rain_trt <- factor(eco_vars_s3$rain_trt, levels=c("W","C","D")) 
eco_vars_s3$year <- factor(eco_vars_s3$year, levels=c("2018","2019","2020")) 
eco_vars_s3$b_p <- as.factor(eco_vars_s3$b_p)

#' Create fourth set of set of contrasts
eco_vars_s4 <- eco_vars_s
eco_vars_s4$graze_trt <- factor(eco_vars_s4$graze_trt, levels=c("C","SG","FG")) 
eco_vars_s4$rain_trt <- factor(eco_vars_s4$rain_trt, levels=c("C","D","W")) 
eco_vars_s4$year <- factor(eco_vars_s4$year, levels=c("2020","2019","2018")) 
eco_vars_s4$b_p <- as.factor(eco_vars_s4$b_p)

#' Create some duplicated outlier dataframes if needed
eco_vars_s_out <- eco_vars_s
eco_vars_s2_out <- eco_vars_s2
eco_vars_s3_out <- eco_vars_s3
eco_vars_s4_out <- eco_vars_s4


#'  
#'  
#'  *Prep FIGURE dataframes* for variables in 2020 only (1 yr of data) and 2018-2020 (3 yrs of data)
#'  
#'  
 

#' We will make figures using the raw (unscaled) data
str(eco_vars)


#' For 2018-2020 variables, we will experiment with a *'percent change'* variable 
#'             (capturing change from 2018 to 2019 and 2020)
# Create 2018 dataframe
ecovars_2018 <- eco_vars %>% filter (year==2018)
# Estimate plot-level percent change for plant biomass
eco_vars$biomass_2018 <- rep(ecovars_2018$biomass_converted, times=3)
eco_vars$biomass_change <- (eco_vars$biomass_converted - eco_vars$biomass_2018)/eco_vars$biomass_2018
# Estimate plot-level percent change for bare ground
eco_vars$bare_2018 <- rep(ecovars_2018$bare_visual_cov, times=3)
eco_vars$bare_change <- (eco_vars$bare_visual_cov - eco_vars$bare_2018)/eco_vars$bare_2018
# Estimate plot-level percent change for Shannon diversity
eco_vars$shan_2018 <- rep(ecovars_2018$shan, times=3)
eco_vars$shan_change <- (eco_vars$shan - eco_vars$shan_2018)/eco_vars$shan_2018


#' For 2018-2020 variables, we will also use raw treatment means in each year
#' 
#' Summarize means and standard errors - for biomass and total cover, 2018 to 2020
eco_vars_means <- eco_vars %>% select(year, rain_trt, graze_trt, biomass_converted, shan, bare_visual_cov, biomass_change, bare_change, shan_change) %>%
  group_by(year, rain_trt, graze_trt) %>%
  summarise_all(list(mean=mean,se=se), na.rm=T) %>%
  arrange(graze_trt, year)
#' Alter variable types and labels for figure-making
eco_vars_means$year <- as.numeric(as.character(eco_vars_means$year))
eco_vars_means$rain_trt <- factor(eco_vars_means$rain_trt, levels=c("D", "C", "W"))
eco_vars_means$graze_trt <- gsub("C", "Ungrazed",eco_vars_means$graze_trt)
eco_vars_means$graze_trt <- gsub("FG", "Dormant-Season Grazing",eco_vars_means$graze_trt)
eco_vars_means$graze_trt <- gsub("SG", "Growing-Season Grazing",eco_vars_means$graze_trt)
eco_vars_means$graze_trt <- factor(eco_vars_means$graze_trt, levels=c( "Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))


#' For 2020 only variables, we will use raw treatment means in that year
#' 
#' Summarize means and standard errors - for phenology and forage quality variables
eco_vars_means_20 <- eco_vars %>% filter(year=="2020") %>%
  select(rain_trt, graze_trt, gs_length, max_gi_doy, crude_protein, digestible_energy, relative_feed_value) %>%
  group_by( rain_trt, graze_trt) %>%
  summarise_all(list(mean=mean,se=se),na.rm=T)
#' Alter variable types and labels for figure-making
eco_vars_means_20$graze_trt <- gsub("C", "Ungrazed",eco_vars_means_20$graze_trt)
eco_vars_means_20$graze_trt <- gsub("FG", "Dormant-Season Grazing",eco_vars_means_20$graze_trt)
eco_vars_means_20$graze_trt <- gsub("SG", "Growing-Season Grazing",eco_vars_means_20$graze_trt)
eco_vars_means_20$graze_trt <- factor(eco_vars_means_20$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))
eco_vars_means_20$rain_trt <- factor(eco_vars_means_20$rain_trt, levels=c("D", "C", "W"))


#' View rainfall treatment means for growing season length only
eco_vars %>% filter(year=="2020") %>%
  select(rain_trt,gs_length) %>%
  group_by( rain_trt) %>%
  summarise_all(list(mean=mean,se=se),na.rm=T)



#' 
#'    [C] : **Rainfall X Grazing Interactions **
#'    
#'    Models and figures for individual response variables
#'    
#'    

#'         
#'  For all linear mixed models (a majority of the models fit here), we expect and 
#'  test for interactions between rainfall and grazing. We use Type III sums of squares
#'  and deviant coding that sums to 0 (see contrasts, just above)
#'  



#'  
#'  *Plant Biomass Production*
#'  

#' *Check data distributions*
#' 
#'    Mostly looks good. 
#'    One outlier across years:      4-7
#'    One major outlier only in 2018:   5-6
#'
eco_vars_s_check <- eco_vars_s
eco_vars_s_check$graze_trt <- gsub("FG", "Dormant-Season Grazing",eco_vars_s_check$graze_trt)
eco_vars_s_check$graze_trt <- gsub("SG", "Growing-Season Grazing", eco_vars_s_check$graze_trt)
eco_vars_s_check$graze_trt <- gsub("C", "Ungrazed", eco_vars_s_check$graze_trt)
eco_vars_s_check$graze_trt <- factor(eco_vars_s_check$graze_trt, levels=c("Ungrazed","Dormant-Season Grazing", "Growing-Season Grazing"))
eco_vars_s_check$rain_trt <- factor(eco_vars_s_check$rain_trt, levels=c("D","C","W"))
eco_vars_s_check$year <- factor(eco_vars_s_check$year, levels=c("2018","2019","2020"))
ggplot(data=eco_vars_s_check) + 
  geom_point(aes(x=rain_trt,y=biomass_converted, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=biomass_converted, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  #geom_violin(aes(x=rain_trt,y=biomass_converted, fill=rain_trt), alpha=0.4)+
  labs(y="Standing biomass\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')


#' *Specify Models*
#'    Full Models (Identical models but showing different contrasts in summary) 
#'    
#'   **Couldn't include random BLOCK effect due to singularity**
#'   
biomass_mod_all  <- lmer(biomass_converted ~ rain_trt * graze_trt * year + (1|b_p), data=eco_vars_s)
biomass_mod_all2 <- lmer(biomass_converted ~ rain_trt * graze_trt * year + (1|b_p), data=eco_vars_s2)
biomass_mod_all3 <- lmer(biomass_converted ~ rain_trt * graze_trt * year + (1|b_p), data=eco_vars_s3)
biomass_mod_all4 <- lmer(biomass_converted ~ rain_trt * graze_trt * year + (1|b_p), data=eco_vars_s4)


# Create outlier dataframe and model
eco_vars_s_out_biomass <- eco_vars_s %>% filter (b_p != "4_7") %>% filter (b_p != "5_6") 
biomass_mod_all_out  <- lmer(biomass_converted ~ rain_trt * graze_trt * year + (1|b_p), data=eco_vars_s_out_biomass)

# Compare outlier and full models
#    Conclusions: No major impact from the two identified outliers on inferences.
#                 There is a stronger effect of two plots with extremely high 
#                 biomass values on diagnostic plots, 6-4 and 6-11, which cause the 
#                 funnel shape in the residual plot. However, these values make complete
#                 ecological sense (high with no grazing, wet treatment, and wet year 2019),
#                 and they are within a normally-distributed set of values for that combination
#                 of treatments. Removing these two outliers doesn't change inferences either
#                 (not shown here), so we will retain all data in analyses.
anova(biomass_mod_all_out)
anova(biomass_mod_all)
plot(biomass_mod_all_out)
plot(biomass_mod_all)
qqnorm(residuals(biomass_mod_all_out))
qqnorm(residuals(biomass_mod_all))

#' Model summaries
Anova(biomass_mod_all)
anova(biomass_mod_all)
r2(biomass_mod_all)
summary(biomass_mod_all)
summary(biomass_mod_all2)
summary(biomass_mod_all3)


#' Save univariate fixed effects from both models (all contrasts) for future figure
biomass_fix <- data.frame(round(fixef(biomass_mod_all),4)) %>% rename(fixef=round.fixef.biomass_mod_all...4.)
biomass_ci  <- data.frame(round(confint(biomass_mod_all),4))
biomass_fixef<- merge(biomass_fix, biomass_ci, by=0, all=TRUE) 
biomass_fixef$response <- "biomass"
biomass_fixef$Row.names <- gsub("rain_trt1", "rain_D", biomass_fixef$Row.names)
biomass_fixef$Row.names <- gsub("rain_trt2", "rain_W", biomass_fixef$Row.names)
biomass_fixef$Row.names <- gsub("graze_trt1", "graze_FG", biomass_fixef$Row.names)
biomass_fixef$Row.names <- gsub("graze_trt2", "graze_SG", biomass_fixef$Row.names)
biomass_fixef$Row.names <- gsub("year1", "year_2020", biomass_fixef$Row.names)
biomass_fixef$Row.names <- gsub("year2", "year_2019", biomass_fixef$Row.names)
biomass_fix2 <- data.frame(round(fixef(biomass_mod_all2),4))%>% rename(fixef=round.fixef.biomass_mod_all2...4.)
biomass_ci2  <- data.frame(round(confint(biomass_mod_all2),4))
biomass_fixef2<- merge(biomass_fix2, biomass_ci2, by=0, all=TRUE) 
biomass_fixef2$response <- "biomass"
biomass_fixef2$Row.names <- gsub("rain_trt1", "rain_C", biomass_fixef2$Row.names)
biomass_fixef2$Row.names <- gsub("rain_trt2", "rain_D", biomass_fixef2$Row.names)
biomass_fixef2$Row.names <- gsub("graze_trt1", "graze_C", biomass_fixef2$Row.names)
biomass_fixef2$Row.names <- gsub("graze_trt2", "graze_FG", biomass_fixef2$Row.names)
biomass_fixef2$Row.names <- gsub("year1", "year_2018", biomass_fixef2$Row.names)
biomass_fixef2$Row.names <- gsub("year2", "year_2019", biomass_fixef2$Row.names)
biomass_fix3 <- data.frame(round(fixef(biomass_mod_all3),4))%>% rename(fixef=round.fixef.biomass_mod_all3...4.)
biomass_ci3  <- data.frame(round(confint(biomass_mod_all3),4))
biomass_fixef3<- merge(biomass_fix3, biomass_ci3, by=0, all=TRUE) 
biomass_fixef3$response <- "biomass"
biomass_fixef3$Row.names <- gsub("rain_trt1", "rain_W", biomass_fixef3$Row.names)
biomass_fixef3$Row.names <- gsub("rain_trt2", "rain_C", biomass_fixef3$Row.names)
biomass_fixef3$Row.names <- gsub("graze_trt1", "graze_C", biomass_fixef3$Row.names)
biomass_fixef3$Row.names <- gsub("graze_trt2", "graze_SG", biomass_fixef3$Row.names)
biomass_fixef3$Row.names <- gsub("year1", "year_2018", biomass_fixef3$Row.names)
biomass_fixef3$Row.names <- gsub("year2", "year_2019", biomass_fixef3$Row.names)
biomass_fix4 <- data.frame(round(fixef(biomass_mod_all4),4))%>% rename(fixef=round.fixef.biomass_mod_all4...4.)
biomass_ci4  <- data.frame(round(confint(biomass_mod_all4),4))
biomass_fixef4<- merge(biomass_fix4, biomass_ci4, by=0, all=TRUE) 
biomass_fixef4$response <- "biomass"
biomass_fixef4$Row.names <- gsub("rain_trt1", "rain_C", biomass_fixef4$Row.names)
biomass_fixef4$Row.names <- gsub("rain_trt2", "rain_D", biomass_fixef4$Row.names)
biomass_fixef4$Row.names <- gsub("graze_trt1", "graze_C", biomass_fixef4$Row.names)
biomass_fixef4$Row.names <- gsub("graze_trt2", "graze_SG", biomass_fixef4$Row.names)
biomass_fixef4$Row.names <- gsub("year1", "year_2020", biomass_fixef4$Row.names)
biomass_fixef4$Row.names <- gsub("year2", "year_2019", biomass_fixef4$Row.names)

biomass_fixef_all <- unique(rbind(biomass_fixef, biomass_fixef2, biomass_fixef3, biomass_fixef4)) %>% arrange(Row.names)
biomass_fixef_all


#' 
#' FIGURES
#' 

#' Raw biomass - paneled by grazing
fig_biomass <- 
  ggplot(data=eco_vars_means, aes(x=year, y=biomass_converted_mean, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.35)) +
  #geom_line(position=position_dodge(0.2),alpha=0.4) +
  geom_errorbar(aes (ymin=biomass_converted_mean-biomass_converted_se, ymax=biomass_converted_mean+biomass_converted_se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  #geom_ribbon(alpha=0.2, aes(ymin=biomass_converted_mean-biomass_converted_se, ymax=biomass_converted_mean+biomass_converted_se, fill=rain_trt)) +
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  #geom_vline(aes(xintercept=2018),lty=2, col="gray50")+
  #geom_text(aes(x=2018.1, y=3000, label="baseline"),angle=90, col="gray50",cex=3)+ 
  labs(x="Year", y=expression(atop(End~of~Season~Biomass,(g~per~m^{2})),""), col="Rain trt.", shape="Graze trt.") +
  ylim(30,350)+
    xlim(2017.5,2020.5) +
  facet_wrap(~graze_trt) +
  ggtitle('Rain***   Graze*    Year***   Rain:Year**    Graze:Year***')+
  theme_bw()+
  theme(legend.position = "none",  panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.title.y = element_text(size = 10))
fig_biomass

fig_biomass_vert <- fig_biomass + theme(axis.title.x = element_blank(),axis.text.x = element_blank())


eco_var_means2 <- eco_vars_means
eco_var_means2$rain_trt <- gsub("D", "Dry", eco_var_means2$rain_trt)
eco_var_means2$rain_trt <- gsub("W", "Wet", eco_var_means2$rain_trt)
eco_var_means2$rain_trt <- gsub("C", "Control", eco_var_means2$rain_trt)
eco_var_means2$rain_trt <- factor(eco_var_means2$rain_trt, levels=c("Dry","Control","Wet"))
fig_biomass2 <- 
  ggplot(data=eco_var_means2, aes(x=year, y=biomass_converted_mean, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=3.2, alpha=0.8,position=position_dodge(0.2)) +
 # geom_line(aes(lty=graze_trt),position=position_dodge(0.2),alpha=0.4) +
  geom_errorbar(aes (ymin=biomass_converted_mean-biomass_converted_se, ymax=biomass_converted_mean+biomass_converted_se), size=1.5,width=0,alpha=0.35,position=position_dodge(0.2))+
  #geom_ribbon(alpha=0.2, aes(ymin=biomass_converted_mean-biomass_converted_se, ymax=biomass_converted_mean+biomass_converted_se, fill=rain_trt)) +
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  # scale_color_manual(values=c("forestgreen","goldenrod2","gray20")) + 
  #geom_vline(aes(xintercept=2018),lty=2, col="gray50")+
  #geom_text(aes(x=2018.1, y=3000, label="baseline"),angle=90, col="gray50",cex=3)+ 
    labs(x="Year", y=expression(atop(End~of~Season~Biomass,(g~per~m^{2})),""), col="Rain trt.", shape="Graze trt.") +
    ylim(0,350)+
  facet_wrap(~rain_trt) +
  ggtitle('Rain***   Graze*    Year***   Rain:Year**    Graze:Year***')+
  theme_bw()+
  theme(legend.position = "none", panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8),
        axis.text.x = element_text(angle = 90), axis.title.y = element_text(size = 10))
fig_biomass2


biomass_6yrs <- ggplot(data=eco_vars_means, aes(x=year, y=biomass_converted_mean, col=rain_trt, shape=graze_trt)) +
   geom_point()+
   geom_line() + 
    xlim(2017,2024) +
    labs(x="Year", y=expression(atop(End~of~Season~Biomass,(g~per~m^{2})),""), col="Rain trt.", shape="Graze trt.") +   
    scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
   theme_bw() +
    facet_wrap(~graze_trt)
biomass_6yrs 
tiff(filename="biomass_2018-2023.jpg", res=600, width=9, height = 6, units = "in")
biomass_6yrs 
dev.off()

#'  
#'  *Bare ground*
#'  

#' Check data distribution
#' 
#' Seems pretty good, there are just a couple of one-time outliers:
#'   2-8, 1-3, 3-4
#'   
#' Plus one plot that, for whatever reason, is totally throwing off the qqplot:
#'    1-5 (Ungrazed, wet)
#' 
ggplot(data=eco_vars_s_check) + 
  geom_point(aes(x=rain_trt,y=bare_visual_cov, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=bare_visual_cov, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  #geom_violin(aes(x=rain_trt,y=bare_visual_cov, fill=rain_trt), alpha=0.4)+
  labs(y="Bare ground\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')


#' *Specify models*
#'

#' Full Model
#'   
bare_mod_all  <- lmer(bare_visual_cov ~ rain_trt * graze_trt * year + (1|b_p) + (1|block), data=eco_vars_s)

#' Outlier dataframe and model -- Remove only the one outlier that is throwing of diagnostic plots 
#'                                           (removing others does not make a large difference)
eco_vars_s_out_bare <- eco_vars_s %>% filter (b_p != "1_5") 
bare_mod_all_out  <- lmer(bare_visual_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data=eco_vars_s_out_bare)

# Compare outlier and full models
#    Conclusions: Remove the one outlier (1-5) that is throwing of diagnostic plots 
#                    (removing others does not make a large difference)
anova(bare_mod_all_out)
anova(bare_mod_all)
plot(bare_mod_all_out)
plot(bare_mod_all)
qqnorm(residuals(bare_mod_all_out))
qqnorm(residuals(bare_mod_all))

#'
#' * Prep for final model with 1-5 removed*
#'
bare_outliers <- eco_vars_s_out$b_p %in% c("1_5")
eco_vars_s_out$bare_visual_cov[bare_outliers] <- NA
eco_vars_s2_out$bare_visual_cov[bare_outliers] <- NA
eco_vars_s3_out$bare_visual_cov[bare_outliers] <- NA
eco_vars_s4_out$bare_visual_cov[bare_outliers] <- NA

bare_outliers_2020 <- eco_vars_out_2020$b_p %in% c("1_5")
eco_vars_out_2020$bare_visual_cov[bare_outliers_2020] <- NA

#' Note outlier in list
outliers_list[6,1] <- '1_5'
outliers_list[6,2] <- 'bare soil'
outliers_list[6,3] <- 'high value in 2018 especially'

#' *Specify models*
#'
#' Full Models (Identical models but showing different contrasts in summary) 
#'   
bare_mod_all  <- lmer(bare_visual_cov ~ rain_trt * graze_trt * year + (1|b_p) + (1|block), data=eco_vars_s_out)
bare_mod_all2 <- lmer(bare_visual_cov ~ rain_trt * graze_trt * year + (1|b_p) + (1|block), data=eco_vars_s2_out)
bare_mod_all3 <- lmer(bare_visual_cov ~ rain_trt * graze_trt * year + (1|b_p) + (1|block), data=eco_vars_s3_out)
bare_mod_all4 <- lmer(bare_visual_cov ~ rain_trt * graze_trt * year + (1|b_p) + (1|block), data=eco_vars_s4_out)

#'
#' Model summaries
#' 
Anova(bare_mod_all)
anova(bare_mod_all)
r2(bare_mod_all)
summary(bare_mod_all)
summary(bare_mod_all2)
summary(bare_mod_all3)


#' Save univariate fixed effects from both models (all contrasts) for future figure
bare_fix <- data.frame(round(fixef(bare_mod_all),4)) %>% rename(fixef=round.fixef.bare_mod_all...4.)
bare_ci  <- data.frame(round(confint(bare_mod_all),4))
bare_fixef<- merge(bare_fix, bare_ci, by=0, all=TRUE) 
bare_fixef$response <- "bare"
bare_fixef$Row.names <- gsub("rain_trt1", "rain_D", bare_fixef$Row.names)
bare_fixef$Row.names <- gsub("rain_trt2", "rain_W", bare_fixef$Row.names)
bare_fixef$Row.names <- gsub("graze_trt1", "graze_FG", bare_fixef$Row.names)
bare_fixef$Row.names <- gsub("graze_trt2", "graze_SG", bare_fixef$Row.names)
bare_fixef$Row.names <- gsub("year1", "year_2020", bare_fixef$Row.names)
bare_fixef$Row.names <- gsub("year2", "year_2019", bare_fixef$Row.names)
bare_fix2 <- data.frame(round(fixef(bare_mod_all2),4))%>% rename(fixef=round.fixef.bare_mod_all2...4.)
bare_ci2  <- data.frame(round(confint(bare_mod_all2),4))
bare_fixef2<- merge(bare_fix2, bare_ci2, by=0, all=TRUE) 
bare_fixef2$response <- "bare"
bare_fixef2$Row.names <- gsub("rain_trt1", "rain_C", bare_fixef2$Row.names)
bare_fixef2$Row.names <- gsub("rain_trt2", "rain_D", bare_fixef2$Row.names)
bare_fixef2$Row.names <- gsub("graze_trt1", "graze_C", bare_fixef2$Row.names)
bare_fixef2$Row.names <- gsub("graze_trt2", "graze_FG", bare_fixef2$Row.names)
bare_fixef2$Row.names <- gsub("year1", "year_2018", bare_fixef2$Row.names)
bare_fixef2$Row.names <- gsub("year2", "year_2019", bare_fixef2$Row.names)
bare_fix3 <- data.frame(round(fixef(bare_mod_all3),4))%>% rename(fixef=round.fixef.bare_mod_all3...4.)
bare_ci3  <- data.frame(round(confint(bare_mod_all3),4))
bare_fixef3<- merge(bare_fix3, bare_ci3, by=0, all=TRUE) 
bare_fixef3$response <- "bare"
bare_fixef3$Row.names <- gsub("rain_trt1", "rain_W", bare_fixef3$Row.names)
bare_fixef3$Row.names <- gsub("rain_trt2", "rain_C", bare_fixef3$Row.names)
bare_fixef3$Row.names <- gsub("graze_trt1", "graze_C", bare_fixef3$Row.names)
bare_fixef3$Row.names <- gsub("graze_trt2", "graze_SG", bare_fixef3$Row.names)
bare_fixef3$Row.names <- gsub("year1", "year_2018", bare_fixef3$Row.names)
bare_fixef3$Row.names <- gsub("year2", "year_2019", bare_fixef3$Row.names)
bare_fix4 <- data.frame(round(fixef(bare_mod_all4),4))%>% rename(fixef=round.fixef.bare_mod_all4...4.)
bare_ci4  <- data.frame(round(confint(bare_mod_all4),4))
bare_fixef4<- merge(bare_fix4, bare_ci4, by=0, all=TRUE) 
bare_fixef4$response <- "bare"
bare_fixef4$Row.names <- gsub("rain_trt1", "rain_C", bare_fixef4$Row.names)
bare_fixef4$Row.names <- gsub("rain_trt2", "rain_D", bare_fixef4$Row.names)
bare_fixef4$Row.names <- gsub("graze_trt1", "graze_C", bare_fixef4$Row.names)
bare_fixef4$Row.names <- gsub("graze_trt2", "graze_SG", bare_fixef4$Row.names)
bare_fixef4$Row.names <- gsub("year1", "year_2020", bare_fixef4$Row.names)
bare_fixef4$Row.names <- gsub("year2", "year_2019", bare_fixef4$Row.names)

bare_fixef_all <- unique(rbind(bare_fixef, bare_fixef2, bare_fixef3, bare_fixef4)) %>% arrange(Row.names)
bare_fixef_all



#' 
#' FIGURES
#' 

#' Re-estimate means and standard errors without outliers
eco_vars_out$bare_visual_cov[bare_outliers] <- NA
eco_vars_means_out <- eco_vars_out %>% select(year, rain_trt, graze_trt, bare_visual_cov) %>%
  group_by(year, rain_trt, graze_trt) %>%
  summarise_all(list(bare_visual_cov_mean=mean,bare_visual_cov_se=se), na.rm=T) %>%
  arrange(graze_trt, year)
#' Alter variable types and labels for figure-making
eco_vars_means_out$year <- as.numeric(as.character(eco_vars_means_out$year))
eco_vars_means_out$rain_trt <- factor(eco_vars_means_out$rain_trt, levels=c("D", "C", "W"))
eco_vars_means_out$graze_trt <- gsub("C", "Ungrazed",eco_vars_means_out$graze_trt)
eco_vars_means_out$graze_trt <- gsub("FG", "Dormant-Season Grazing",eco_vars_means_out$graze_trt)
eco_vars_means_out$graze_trt <- gsub("SG", "Growing-Season Grazing",eco_vars_means_out$graze_trt)
eco_vars_means_out$graze_trt <- factor(eco_vars_means_out$graze_trt, levels=c("Ungrazed", "Growing-Season Grazing", "Dormant-Season Grazing"))

#' Raw cover of bare ground - paneled by grazing
fig_bare <- 
  ggplot(data=eco_vars_means_out, aes(x=year, y=bare_visual_cov_mean, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=2.5, alpha=0.8,position=position_dodge(0.35)) +
  #geom_line(position=position_dodge(0.2),alpha=0.4) +
  geom_errorbar(aes (ymin=bare_visual_cov_mean-bare_visual_cov_se, ymax=bare_visual_cov_mean+bare_visual_cov_se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  #geom_ribbon(alpha=0.2, aes(ymin=bare_visual_cov_mean-bare_visual_cov_se, ymax=bare_visual_cov_mean+bare_visual_cov_se, fill=rain_trt)) +
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  #geom_vline(aes(xintercept=2018),lty=2, col="gray50")+
  #geom_text(aes(x=2018.1, y=5, label="baseline"),angle=90, col="gray50",cex=3)+ 
  labs(x="Year", y="Bare Soil Cover (% aerial)\n", col="Rain trt.", shape="Graze trt.") +
  ylim(8,50)+
    xlim(2017.5,2020.5) +
  facet_wrap(~graze_trt) +
  ggtitle('Rain***   Rain:Year*  Graze:Year**')+
  theme_bw()+
  theme( legend.text = element_text(size=8.5), legend.title = element_text(size=8.5), panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8),  axis.title.y = element_text(size = 10))
fig_bare  

fig_bare_vert <- fig_bare + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
fig_bare_vert2 <- fig_bare_vert + theme(legend.position = 'none')


fig_bare2 <- 
  ggplot(data=eco_vars_means_out, aes(x=year, y=bare_visual_cov_mean, col=rain_trt, shape=graze_trt)) +
  geom_point( cex=3.2, alpha=0.8,position=position_dodge(0.2)) +
  #geom_line(aes(lty=graze_trt),size=1,position=position_dodge(0.2),alpha=0.4) +
  geom_errorbar(aes (ymin=bare_visual_cov_mean-bare_visual_cov_se, ymax=bare_visual_cov_mean+bare_visual_cov_se), size=2,width=0,alpha=0.35,position=position_dodge(0.2))+
  #geom_ribbon(alpha=0.2, aes(ymin=bare_visual_cov_mean-bare_visual_cov_se, ymax=bare_visual_cov_mean+bare_visual_cov_se, fill=rain_trt)) +
  scale_x_continuous(breaks=c(2018, 2019, 2020), limits=c(2017.75,2020.25)) + 
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  #scale_color_manual(values=c("forestgreen","goldenrod2","gray20")) + 
 #  geom_vline(aes(xintercept=2018),lty=2, col="gray50")+
  #geom_text(aes(x=2018.1, y=5, label="baseline"),angle=90, col="gray50",cex=3)+ 
  labs(x="Year", y="Bare Soil Cover (% aerial)\n", col="Rain trt.", lty="Graze trt.", shape="Graze trt.") +
  ylim(0,50)+
  facet_wrap(~rain_trt) +
  ggtitle('Rain***   Rain:Year*  Graze:Year**')+
  theme_bw()+
  theme( legend.text = element_text(size=8.5), legend.title = element_text(size=8.5), 
         panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(angle = 90), axis.title.y = element_text(size = 10))
fig_bare2



bare_6yrs <- ggplot(data=eco_vars_means, aes(x=year, y=bare_visual_cov_mean, col=rain_trt, shape=graze_trt)) +
    geom_point()+
    geom_line() + 
    xlim(2017,2023) +
    labs(x="Year", y="Bare Soil Cover (% aerial)\n", col="Rain trt.", lty="Graze trt.", shape="Graze trt.") +
    scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
    theme_bw() +
    facet_wrap(~graze_trt)
bare_6yrs 
tiff(filename="bare_2018-2022.jpg", res=600, width=9, height = 6, units = "in")
bare_6yrs 
dev.off()


#'
#' *Growing season length*
#'
#' Check data distribution
#' 
#' Seems pretty good.
#' 
eco_vars_s_2020_check <- eco_vars_2020_s
eco_vars_s_2020_check$graze_trt <- gsub("FG", "Dormant-Season Grazing",eco_vars_s_2020_check$graze_trt)
eco_vars_s_2020_check$graze_trt <- gsub("SG", "Growing-Season Grazing", eco_vars_s_2020_check$graze_trt)
eco_vars_s_2020_check$graze_trt <- gsub("C", "Ungrazed", eco_vars_s_2020_check$graze_trt)
eco_vars_s_2020_check$raze_trt <- factor(eco_vars_s_2020_check$graze_trt, levels=c("Growing-Season Grazing","Dormant-Season Grazing", "Ungrazed"))
eco_vars_s_2020_check$rain_trt <- factor(eco_vars_s_2020_check$rain_trt, levels=c("D","C","W"))
eco_vars_s_2020_check$year <- factor(eco_vars_s_2020_check$year, levels=c("2018","2019","2020"))
ggplot(data=eco_vars_s_2020_check) + 
  geom_point(aes(x=rain_trt,y=gs_length, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=gs_length, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  #geom_violin(aes(x=rain_trt,y=gs_length, fill=rain_trt), alpha=0.4)+
  labs(y="Growing Season Length\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

# Models (Identical models but showing different contrasts in summary)
gsl_mod <- lmer(gs_length ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s)
gsl_mod2 <- lmer(gs_length ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s2)
gsl_mod3 <- lmer(gs_length ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s3)


# Summaries
plot(gsl_mod)
qqnorm(residuals(gsl_mod))
Anova(gsl_mod)
anova(gsl_mod)
r2(gsl_mod)
summary(gsl_mod)
summary(gsl_mod2)

# Check block effect on GSL
ggplot(data=eco_vars_2020_s) + geom_point(aes(x=block, y=gs_length,col=rain_trt))


# Save univariate fixed effects from both models (all contrasts)
gsl_fix <- data.frame(round(fixef(gsl_mod),4)) %>% rename(fixef=round.fixef.gsl_mod...4.)
gsl_ci  <- data.frame(round(confint(gsl_mod),4))
gsl_fixef<- merge(gsl_fix, gsl_ci, by=0, all=TRUE) 
gsl_fixef$response <- "gsl"
gsl_fixef$Row.names <- gsub("rain_trt1", "rain_D", gsl_fixef$Row.names)
gsl_fixef$Row.names <- gsub("rain_trt2", "rain_W", gsl_fixef$Row.names)
gsl_fixef$Row.names <- gsub("graze_trt1", "graze_FG", gsl_fixef$Row.names)
gsl_fixef$Row.names <- gsub("graze_trt2", "graze_SG", gsl_fixef$Row.names)
gsl_fix2 <- data.frame(round(fixef(gsl_mod2),4))%>% rename(fixef=round.fixef.gsl_mod2...4.)
gsl_ci2  <- data.frame(round(confint(gsl_mod2),4))
gsl_fixef2<- merge(gsl_fix2, gsl_ci2, by=0, all=TRUE) 
gsl_fixef2$response <- "gsl"
gsl_fixef2$Row.names <- gsub("rain_trt1", "rain_C", gsl_fixef2$Row.names)
gsl_fixef2$Row.names <- gsub("rain_trt2", "rain_D", gsl_fixef2$Row.names)
gsl_fixef2$Row.names <- gsub("graze_trt1", "graze_C", gsl_fixef2$Row.names)
gsl_fixef2$Row.names <- gsub("graze_trt2", "graze_FG", gsl_fixef2$Row.names)
gsl_fixef_all <- rbind(gsl_fixef, gsl_fixef3)
gsl_fix3 <- data.frame(round(fixef(gsl_mod3),4))%>% rename(fixef=round.fixef.gsl_mod3...4.)
gsl_ci3  <- data.frame(round(confint(gsl_mod3),4))
gsl_fixef3<- merge(gsl_fix3, gsl_ci3, by=0, all=TRUE) 
gsl_fixef3$response <- "gsl"
gsl_fixef3$Row.names <- gsub("rain_trt1", "rain_W", gsl_fixef3$Row.names)
gsl_fixef3$Row.names <- gsub("rain_trt2", "rain_C", gsl_fixef3$Row.names)
gsl_fixef3$Row.names <- gsub("graze_trt1", "graze_C", gsl_fixef3$Row.names)
gsl_fixef3$Row.names <- gsub("graze_trt2", "graze_SG", gsl_fixef3$Row.names)

gsl_fixef_all <- unique(rbind(gsl_fixef, gsl_fixef2, gsl_fixef3)) %>% arrange(Row.names)
gsl_fixef_all


#' FIGURE
#' 
fig_gsl <-ggplot( data = eco_vars_means_20, aes(x=graze_trt, y=gs_length_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8,position=position_dodge(0.35)) + 
  geom_errorbar(aes (ymin=gs_length_mean-gs_length_se, ymax=gs_length_mean+gs_length_se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=220, label="Rain***"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=220, label="a"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Growing Season Length \n(days above greenness threshold)", col="Rain trt.", shape="Graze trt.") +
  ylim(0,230) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 10), axis.text.x = element_text(angle = 45))
fig_gsl

fig_gsl2 <-ggplot( data = eco_vars_means_20, aes(x=rain_trt, y=gs_length_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8,position=position_dodge(0.2)) + 
  geom_errorbar(aes (ymin=gs_length_mean-gs_length_se, ymax=gs_length_mean+gs_length_se), size=2,width=0,alpha=0.35,position=position_dodge(0.2))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=220, label="Rain***"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=220, label="a"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Growing Season Length \n(days above greenness threshold)", col="Rain trt.", shape="Graze trt.") +
  ylim(0,230) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 10))
fig_gsl2




#'
#' *Peak greenness*
#'

#' **CHECK THE PHENO CURVES ON THESE OUTLIERS AS JUSTIFICATION **


#' Check data distributions
#'
#' Five possible outliers skewing the model:
#'  1-11 - early date is not misleading; greenness curve suggests slightly early date
#'  2-4 - later date does not appear to be misleading; matches greenness curve
#'  3-6 - later date is misleading; has a very shallow bi-modal greenness curve 
#'  4-11 - later date is misleading; greenness curve is very shallow and flat, max doy should not be interpreted
#'  6-3  -   greenness curve is very shallow and flat, max doy should not be interpreted
#' 
#' Justified removal of 3-6, 4-11, 6-3
#' 
ggplot(data=eco_vars_s_2020_check) + 
  geom_point(aes(x=rain_trt,y=max_gi_doy, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=max_gi_doy, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  #geom_violin(aes(x=rain_trt,y=max_gi_doy, fill=rain_trt), alpha=0.4)+
  labs(y="Day of Peak Greenness\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod3","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')


#' Remove three outliers and update dataframes
doy_outliers <- eco_vars_2020$b_p %in% c("3_6","4_11","6_3")
eco_vars_out_2020$max_gi_doy[doy_outliers] <- NA
max_gi_doy_s <- scale(eco_vars_out_2020$max_gi_doy)
eco_vars_2020_s_out$max_gi_doy <- max_gi_doy_s 
eco_vars_2020_s2_out <- eco_vars_2020_s_out
eco_vars_2020_s3_out <- eco_vars_2020_s_out
#' Reset contrasts just in case
eco_vars_2020_s_out$rain_trt  <- factor(eco_vars_2020_s$rain_trt, levels=c("D","W","C"))
eco_vars_2020_s_out$graze_trt <- factor(eco_vars_2020_s$graze_trt, levels=c("FG","SG","C"))
eco_vars_2020_s2_out$rain_trt <- factor(eco_vars_2020_s2$rain_trt, levels=c("C","D","W"))
eco_vars_2020_s2_out$graze_trt <- factor(eco_vars_2020_s2$graze_trt, levels=c("C","SG","FG"))
eco_vars_2020_s3_out$rain_trt  <- factor(eco_vars_2020_s3$rain_trt, levels=c("W","C","D"))
eco_vars_2020_s3_out$graze_trt  <- factor(eco_vars_2020_s3$graze_trt, levels=c("C","SG","FG"))

# Full Model 
#
peakg_mod <- lmer(max_gi_doy ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s)
# Outlier dataframe and model -- Remove outliers with scaled values greater than
peakg_mod_out  <- lmer(max_gi_doy ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s_out)

# Compare outlier and full models
#    Conclusions: Removing these outliers changes model  inferences but also 
#               improves diagnostics.. although there are a few other outliers
#                evident in the qqplot, these seem ecologically valid after reviewing
#                their greenness curves, and removing them does not further change inferences
anova(peakg_mod_out)
anova(peakg_mod)
plot(peakg_mod_out)
plot(peakg_mod)
qqnorm(residuals(peakg_mod_out))
qqnorm(residuals(peakg_mod))

#' Add to outliers list
outliers_list[7,1] <- '3_6'
outliers_list[7,2] <- 'peakg'
outliers_list[7,3] <- 'greenness curve not reliable to estimate'
outliers_list[8,1] <- '4_11'
outliers_list[8,2] <- 'peakg'
outliers_list[8,3] <- 'greenness curve not reliable to estimate'
outliers_list[9,1] <-  '6_3'
outliers_list[9,2] <- 'peakg'
outliers_list[9,3] <- 'greenness curve not reliable to estimate'


#' *Specify models*
#'
#' Full Models (Identical models but showing different contrasts in summary) 
#'   
peakg_mod  <- lmer(max_gi_doy ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s_out)
peakg_mod2 <- lmer(max_gi_doy ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s2_out)
peakg_mod3 <- lmer(max_gi_doy ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s3_out)


# Summaries
Anova(peakg_mod)
anova(peakg_mod)
r2(peakg_mod)
summary(peakg_mod)
summary(peakg_mod3)

# Save univariate fixed effects from both models (all contrasts)
peakg_fix <- data.frame(round(fixef(peakg_mod),4)) %>% rename(fixef=round.fixef.peakg_mod...4.)
peakg_ci  <- data.frame(round(confint(peakg_mod),4))
peakg_fixef<- merge(peakg_fix, peakg_ci, by=0, all=TRUE) 
peakg_fixef$response <- "peakg"
peakg_fixef$Row.names <- gsub("rain_trt1", "rain_D", peakg_fixef$Row.names)
peakg_fixef$Row.names <- gsub("rain_trt2", "rain_W", peakg_fixef$Row.names)
peakg_fixef$Row.names <- gsub("graze_trt1", "graze_FG", peakg_fixef$Row.names)
peakg_fixef$Row.names <- gsub("graze_trt2", "graze_SG", peakg_fixef$Row.names)
peakg_fix2 <- data.frame(round(fixef(peakg_mod2),4))%>% rename(fixef=round.fixef.peakg_mod2...4.)
peakg_ci2  <- data.frame(round(confint(peakg_mod2),4))
peakg_fixef2<- merge(peakg_fix2, peakg_ci2, by=0, all=TRUE) 
peakg_fixef2$response <- "peakg"
peakg_fixef2$Row.names <- gsub("rain_trt1", "rain_C", peakg_fixef2$Row.names)
peakg_fixef2$Row.names <- gsub("rain_trt2", "rain_D", peakg_fixef2$Row.names)
peakg_fixef2$Row.names <- gsub("graze_trt1", "graze_C", peakg_fixef2$Row.names)
peakg_fixef2$Row.names <- gsub("graze_trt2", "graze_FG", peakg_fixef2$Row.names)
peakg_fixef_all <- rbind(peakg_fixef, peakg_fixef3)
peakg_fix3 <- data.frame(round(fixef(peakg_mod3),4))%>% rename(fixef=round.fixef.peakg_mod3...4.)
peakg_ci3  <- data.frame(round(confint(peakg_mod3),4))
peakg_fixef3<- merge(peakg_fix3, peakg_ci3, by=0, all=TRUE) 
peakg_fixef3$response <- "peakg"
peakg_fixef3$Row.names <- gsub("rain_trt1", "rain_W", peakg_fixef3$Row.names)
peakg_fixef3$Row.names <- gsub("rain_trt2", "rain_C", peakg_fixef3$Row.names)
peakg_fixef3$Row.names <- gsub("graze_trt1", "graze_C", peakg_fixef3$Row.names)
peakg_fixef3$Row.names <- gsub("graze_trt2", "graze_SG", peakg_fixef3$Row.names)

peakg_fixef_all <- unique(rbind(peakg_fixef, peakg_fixef2, peakg_fixef3)) %>% arrange(Row.names)
peakg_fixef_all


#' FIGURE

#' Re-estimate means and standard errors without outliers
eco_vars_means_out_2020 <- eco_vars_out_2020 %>% select(year, rain_trt, graze_trt, max_gi_doy) %>%
  group_by(year, rain_trt, graze_trt) %>%
  summarise_all(list(max_gi_doy_mean=mean,max_gi_doy_se=se), na.rm=T) %>%
  arrange(graze_trt, year)
#' Alter variable types and labels for figure-making
eco_vars_means_out_2020$year <- as.numeric(as.character(eco_vars_means_out_2020$year))
eco_vars_means_out_2020$rain_trt <- factor(eco_vars_means_out_2020$rain_trt, levels=c("D", "C", "W"))
eco_vars_means_out_2020$graze_trt <- gsub("C", "Ungrazed",eco_vars_means_out_2020$graze_trt)
eco_vars_means_out_2020$graze_trt <- gsub("FG", "Dormant-Season Grazing",eco_vars_means_out_2020$graze_trt)
eco_vars_means_out_2020$graze_trt <- gsub("SG", "Growing-Season Grazing",eco_vars_means_out_2020$graze_trt)
eco_vars_means_out_2020$graze_trt <- factor(eco_vars_means_out_2020$graze_trt, levels=c("Ungrazed","Growing-Season Grazing", "Dormant-Season Grazing" ))

fig_peakg <-ggplot( data = eco_vars_means_out_2020, aes(x=graze_trt, y=max_gi_doy_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8,position=position_dodge(0.35)) + 
  geom_errorbar(aes (ymin=max_gi_doy_mean-max_gi_doy_se, ymax=max_gi_doy_mean+max_gi_doy_se), size=1,width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=200, label="Graze*"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=200, label="b"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Peak Greenness \n(day of year)", col="Rain trt.", shape="Graze trt.") +
  ylim(140,205) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 10), axis.text.x = element_text(angle = 45))
fig_peakg

fig_peakg2 <-ggplot( data = eco_vars_means_out_2020, aes(x=rain_trt, y=max_gi_doy_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8,position=position_dodge(0.2)) + 
  geom_errorbar(aes (ymin=max_gi_doy_mean-max_gi_doy_se, ymax=max_gi_doy_mean+max_gi_doy_se), size=2,width=0,alpha=0.35,position=position_dodge(0.2))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=200, label="Graze*"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=200, label="b"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Peak Greenness \n(day of year)", col="Rain trt.", shape="Graze trt.") +
  ylim(140,205) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 10))
fig_peakg2


#'
#' *Crude protein*
#'

#' Check data distributions
#'
#'   Looks good
#' 
ggplot(data=eco_vars_s_2020_check) + 
  geom_point(aes(x=rain_trt,y=crude_protein, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=crude_protein, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  #geom_violin(aes(x=rain_trt,y=crude_protein, fill=rain_trt), alpha=0.4)+
  labs(y="Crude Protein\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

# Models (Identical models but showing different contrasts in summary)
protein_mod <- lmer(crude_protein ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s)
protein_mod2 <- lmer(crude_protein ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s2)
protein_mod3 <- lmer(crude_protein ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s3)

# Summaries
plot(protein_mod)
qqnorm(residuals(protein_mod))
Anova(protein_mod)
anova(protein_mod)
r2(protein_mod)
summary(protein_mod)
summary(protein_mod2)
summary(protein_mod3)

# Check block effect on crude protein
ggplot(data=eco_vars_2020_s) + geom_point(aes(x=block, y=crude_protein,col=rain_trt))

# Save univariate fixed effects from both models (all contrasts)
protein_fix <- data.frame(round(fixef(protein_mod),4)) %>% rename(fixef=round.fixef.protein_mod...4.)
protein_ci  <- data.frame(round(confint(protein_mod),4))
protein_fixef<- merge(protein_fix, protein_ci, by=0, all=TRUE) 
protein_fixef$response <- "protein"
protein_fixef$Row.names <- gsub("rain_trt1", "rain_D", protein_fixef$Row.names)
protein_fixef$Row.names <- gsub("rain_trt2", "rain_W", protein_fixef$Row.names)
protein_fixef$Row.names <- gsub("graze_trt1", "graze_FG", protein_fixef$Row.names)
protein_fixef$Row.names <- gsub("graze_trt2", "graze_SG", protein_fixef$Row.names)
protein_fix2 <- data.frame(round(fixef(protein_mod2),4))%>% rename(fixef=round.fixef.protein_mod2...4.)
protein_ci2  <- data.frame(round(confint(protein_mod2),4))
protein_fixef2<- merge(protein_fix2, protein_ci2, by=0, all=TRUE) 
protein_fixef2$response <- "protein"
protein_fixef2$Row.names <- gsub("rain_trt1", "rain_C", protein_fixef2$Row.names)
protein_fixef2$Row.names <- gsub("rain_trt2", "rain_D", protein_fixef2$Row.names)
protein_fixef2$Row.names <- gsub("graze_trt1", "graze_C", protein_fixef2$Row.names)
protein_fixef2$Row.names <- gsub("graze_trt2", "graze_FG", protein_fixef2$Row.names)
protein_fixef_all <- rbind(protein_fixef, protein_fixef3)
protein_fix3 <- data.frame(round(fixef(protein_mod3),4))%>% rename(fixef=round.fixef.protein_mod3...4.)
protein_ci3  <- data.frame(round(confint(protein_mod3),4))
protein_fixef3<- merge(protein_fix3, protein_ci3, by=0, all=TRUE) 
protein_fixef3$response <- "protein"
protein_fixef3$Row.names <- gsub("rain_trt1", "rain_W", protein_fixef3$Row.names)
protein_fixef3$Row.names <- gsub("rain_trt2", "rain_C", protein_fixef3$Row.names)
protein_fixef3$Row.names <- gsub("graze_trt1", "graze_C", protein_fixef3$Row.names)
protein_fixef3$Row.names <- gsub("graze_trt2", "graze_SG", protein_fixef3$Row.names)

protein_fixef_all <- unique(rbind(protein_fixef, protein_fixef2, protein_fixef3)) %>% arrange(Row.names)
protein_fixef_all



#' FIGURE
#' 
#' Revised
fig_protein <-ggplot( data = eco_vars_means_20, aes(x=graze_trt, y=crude_protein_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8, position=position_dodge(0.35)) + 
  geom_errorbar(aes (ymin=crude_protein_mean-crude_protein_se, ymax=crude_protein_mean+crude_protein_se),size=1, width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=10.2, label="NS"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=10.2, label="c"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Crude protein \n (%, forage quality analysis)", col="Rain trt.", shape="Graze trt.") +
  ylim(5, 10.5) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 10), axis.text.x = element_text(angle = 45))
fig_protein

fig_protein2 <-ggplot( data = eco_vars_means_20, aes(x=rain_trt, y=crude_protein_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8, position=position_dodge(0.2)) + 
  geom_errorbar(aes (ymin=crude_protein_mean-crude_protein_se, ymax=crude_protein_mean+crude_protein_se),size=2, width=0,alpha=0.35,position=position_dodge(0.2))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=10.2, label="NS"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=10.2, label="c"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Crude protein \n (%, forage quality analysis)", col="Rain trt.", shape="Graze trt.") +
  ylim(5, 10.5) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 10))
fig_protein2



#'
#' *Digestible energy*
#'

#' Check data distributions
#'
#'   Looks good
#' 
ggplot(data=eco_vars_s_2020_check) + 
  geom_point(aes(x=rain_trt,y=digestible_energy, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=digestible_energy, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  #geom_violin(aes(x=rain_trt,y=digestible_energy, fill=rain_trt), alpha=0.4)+
  labs(y="Digestible Energy\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

# Models (Identical models but showing different contrasts in summary)
energy_mod  <- lmer(digestible_energy ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s)
energy_mod2 <- lmer(digestible_energy ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s2)
energy_mod3 <- lmer(digestible_energy ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s3)

# Summaries
plot(energy_mod)
qqnorm(residuals(energy_mod))
Anova(energy_mod)
anova(energy_mod)
r2(energy_mod)
summary(energy_mod)
summary(energy_mod2)

# Save univariate fixed effects from both models (all contrasts)
energy_fix <- data.frame(round(fixef(energy_mod),4)) %>% rename(fixef=round.fixef.energy_mod...4.)
energy_ci  <- data.frame(round(confint(energy_mod),4))
energy_fixef<- merge(energy_fix, energy_ci, by=0, all=TRUE) 
energy_fixef$response <- "energy"
energy_fixef$Row.names <- gsub("rain_trt1", "rain_D", energy_fixef$Row.names)
energy_fixef$Row.names <- gsub("rain_trt2", "rain_W", energy_fixef$Row.names)
energy_fixef$Row.names <- gsub("graze_trt1", "graze_FG", energy_fixef$Row.names)
energy_fixef$Row.names <- gsub("graze_trt2", "graze_SG", energy_fixef$Row.names)
energy_fix2 <- data.frame(round(fixef(energy_mod2),4))%>% rename(fixef=round.fixef.energy_mod2...4.)
energy_ci2  <- data.frame(round(confint(energy_mod2),4))
energy_fixef2<- merge(energy_fix2, energy_ci2, by=0, all=TRUE) 
energy_fixef2$response <- "energy"
energy_fixef2$Row.names <- gsub("rain_trt1", "rain_C", energy_fixef2$Row.names)
energy_fixef2$Row.names <- gsub("rain_trt2", "rain_D", energy_fixef2$Row.names)
energy_fixef2$Row.names <- gsub("graze_trt1", "graze_C", energy_fixef2$Row.names)
energy_fixef2$Row.names <- gsub("graze_trt2", "graze_FG", energy_fixef2$Row.names)
energy_fixef_all <- rbind(energy_fixef, energy_fixef3)
energy_fix3 <- data.frame(round(fixef(energy_mod3),4))%>% rename(fixef=round.fixef.energy_mod3...4.)
energy_ci3  <- data.frame(round(confint(energy_mod3),4))
energy_fixef3<- merge(energy_fix3, energy_ci3, by=0, all=TRUE) 
energy_fixef3$response <- "energy"
energy_fixef3$Row.names <- gsub("rain_trt1", "rain_W", energy_fixef3$Row.names)
energy_fixef3$Row.names <- gsub("rain_trt2", "rain_C", energy_fixef3$Row.names)
energy_fixef3$Row.names <- gsub("graze_trt1", "graze_C", energy_fixef3$Row.names)
energy_fixef3$Row.names <- gsub("graze_trt2", "graze_SG", energy_fixef3$Row.names)

energy_fixef_all <- unique(rbind(energy_fixef, energy_fixef2, energy_fixef3)) %>% arrange(Row.names)
energy_fixef_all


#' FIGURE
#' 
#' 
#' Revised
fig_energy <-ggplot( data = eco_vars_means_20, aes(x=graze_trt, y=digestible_energy_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8, position=position_dodge(0.35)) + 
  geom_errorbar(aes (ymin=digestible_energy_mean-digestible_energy_se, ymax=digestible_energy_mean+digestible_energy_se),size=1, width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=1.32, label="NS"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=1.32, label="d"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Digestible energy \n (Mcal/lb, forage quality analysis)", col="Rain trt.", shape="Graze trt.") +
  ylim(1.15, 1.33) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 10), axis.text.x = element_text(angle = 45))
fig_energy

fig_energy2 <-ggplot( data = eco_vars_means_20, aes(x=rain_trt, y=digestible_energy_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8, position=position_dodge(0.2)) + 
  geom_errorbar(aes (ymin=digestible_energy_mean-digestible_energy_se, ymax=digestible_energy_mean+digestible_energy_se),size=2, width=0,alpha=0.35,position=position_dodge(0.2))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=1.32, label="NS"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=1.32, label="d"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Digestible energy \n (Mcal/lb, forage quality analysis)", col="Rain trt.", shape="Graze trt.") +
  ylim(1.15, 1.33) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 10))
fig_energy2



#'
#' *Relative Feed Value*
#'

#' Check data distributions
#'
#'   There's one bigger outlier skewing the residual qqplot (further down): 1-5
#'   
ggplot(data=eco_vars_s_2020_check) + 
  geom_point(aes(x=rain_trt,y=relative_feed_value, color=rain_trt), alpha=0.4, position=position_jitter(width=0.12))+
  geom_boxplot(aes(x=rain_trt,y=relative_feed_value, fill=rain_trt), alpha=0.4, outlier.shape=8, outlier.colour = 'darkred')+
  #geom_violin(aes(x=rain_trt,y=relative_feed_value, fill=rain_trt), alpha=0.4)+
  labs(y="Digestible Energy\n")+
  scale_fill_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  theme(axis.text = element_blank())+
  theme_bw() +
  facet_grid(~year ~graze_trt,scale='free')

# Models (Identical models but showing different contrasts in summary)
rfv_mod <- lmer(relative_feed_value ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s)
rfv_mod2 <- lmer(relative_feed_value ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s2)
rfv_mod3 <- lmer(relative_feed_value ~ rain_trt * graze_trt + (1|block), data=eco_vars_2020_s3)

# Create outlier dataframe and model 
eco_vars_s_out_rfv <- eco_vars_2020_s %>% filter (b_p != "1_5") 
rfv_mod_out  <- lmer(relative_feed_value ~ rain_trt * graze_trt + (1|block), data=eco_vars_s_out_rfv)

# Compare outlier and full models
#    Conclusions: Removing this outlier doesn't seem to make a major difference for inference,
#                 so we will run models using all data. There is a bit more variability around
#                 higher fitted values of digestible energy because there's a bit of a positive skew, 
#                 especially in the ungrazed treatment.. but if these were removed, it would not likely
#                 reveal any missed effects
anova(rfv_mod_out)
anova(rfv_mod)
plot(rfv_mod_out)
plot(rfv_mod)
qqnorm(residuals(rfv_mod_out))
qqnorm(residuals(rfv_mod))

# Summaries
Anova(rfv_mod)
anova(rfv_mod)
r2(rfv_mod)
summary(rfv_mod)
summary(rfv_mod2)
contrasts(eco_vars_2020_s2$rain_trt)
contrasts(eco_vars_2020_s2$graze_trt)

# Save univariate fixed effects from both models (all contrasts)
rfv_fix <- data.frame(round(fixef(rfv_mod),4)) %>% rename(fixef=round.fixef.rfv_mod...4.)
rfv_ci  <- data.frame(round(confint(rfv_mod),4))
rfv_fixef<- merge(rfv_fix, rfv_ci, by=0, all=TRUE) 
rfv_fixef$response <- "rfv"
rfv_fixef$Row.names <- gsub("rain_trt1", "rain_D", rfv_fixef$Row.names)
rfv_fixef$Row.names <- gsub("rain_trt2", "rain_W", rfv_fixef$Row.names)
rfv_fixef$Row.names <- gsub("graze_trt1", "graze_FG", rfv_fixef$Row.names)
rfv_fixef$Row.names <- gsub("graze_trt2", "graze_SG", rfv_fixef$Row.names)
rfv_fix2 <- data.frame(round(fixef(rfv_mod2),4))%>% rename(fixef=round.fixef.rfv_mod2...4.)
rfv_ci2  <- data.frame(round(confint(rfv_mod2),4))
rfv_fixef2<- merge(rfv_fix2, rfv_ci2, by=0, all=TRUE) 
rfv_fixef2$response <- "rfv"
rfv_fixef2$Row.names <- gsub("rain_trt1", "rain_C", rfv_fixef2$Row.names)
rfv_fixef2$Row.names <- gsub("rain_trt2", "rain_D", rfv_fixef2$Row.names)
rfv_fixef2$Row.names <- gsub("graze_trt1", "graze_C", rfv_fixef2$Row.names)
rfv_fixef2$Row.names <- gsub("graze_trt2", "graze_FG", rfv_fixef2$Row.names)
rfv_fixef_all <- rbind(rfv_fixef, rfv_fixef3)
rfv_fix3 <- data.frame(round(fixef(rfv_mod3),4))%>% rename(fixef=round.fixef.rfv_mod3...4.)
rfv_ci3  <- data.frame(round(confint(rfv_mod3),4))
rfv_fixef3<- merge(rfv_fix3, rfv_ci3, by=0, all=TRUE) 
rfv_fixef3$response <- "rfv"
rfv_fixef3$Row.names <- gsub("rain_trt1", "rain_W", rfv_fixef3$Row.names)
rfv_fixef3$Row.names <- gsub("rain_trt2", "rain_C", rfv_fixef3$Row.names)
rfv_fixef3$Row.names <- gsub("graze_trt1", "graze_C", rfv_fixef3$Row.names)
rfv_fixef3$Row.names <- gsub("graze_trt2", "graze_SG", rfv_fixef3$Row.names)

rfv_fixef_all <- unique(rbind(rfv_fixef, rfv_fixef2, rfv_fixef3)) %>% arrange(Row.names)
rfv_fixef_all


#' FIGURE
#'
#' Revised
fig_rfv <-ggplot( data = eco_vars_means_20, aes(x=graze_trt, y=relative_feed_value_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8, position=position_dodge(0.35)) + 
  geom_errorbar(aes (ymin=relative_feed_value_mean-relative_feed_value_se, ymax=relative_feed_value_mean+relative_feed_value_se), size=1, width=0,alpha=0.45,position=position_dodge(0.35))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=160, label="NS"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=160, label="e"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Relative Feed Value", col="Rain trt.", shape="Graze trt.") +
  ylim(50, 170) +
  theme_bw() +  theme(legend.text = element_text(size=8.5), axis.text.x = element_text(angle = 45),  legend.title = element_text(size=8.5), axis.title.y = element_text(size = 10))
fig_rfv


fig_rfv2 <-ggplot( data = eco_vars_means_20, aes(x=rain_trt, y=relative_feed_value_mean, col=rain_trt, shape=graze_trt)) +
  geom_point(cex=3.2, alpha=0.8, position=position_dodge(0.2)) + 
  geom_errorbar(aes (ymin=relative_feed_value_mean-relative_feed_value_se, ymax=relative_feed_value_mean+relative_feed_value_se), size=2, width=0,alpha=0.35,position=position_dodge(0.2))+
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) + 
  geom_text(aes(x=2, y=160, label="NS"), cex=3, col="black")+
  geom_text(aes(x=0.75, y=160, label="e"),fontface="bold", cex=4, col="black")+
  labs(x=NULL, y="Relative Feed Value", col="Rain trt.", shape="Graze trt.") +
  ylim(50, 170) +
  theme_bw() +  theme(legend.text = element_text(size=8.5),  legend.title = element_text(size=8.5), axis.title.y = element_text(size = 10))
fig_rfv2




#'
#' **ADDITIONAL FIGURES**
#' 
#'    Relationship between plant biomass and diversity, since these seem to have opposite responses
#'    

str(eco_vars_2020, give.attr=F)
str(eco_vars, give.attr=F)
eco_vars_change <- eco_vars %>% filter( year==2020 ) %>% select(rain_trt, graze_trt, biomass_change, shan_change)

ggplot(data=eco_vars_change) + 
  geom_point( aes( x=biomass_change, y=shan_change, col=rain_trt, shape=graze_trt) ) +
  geom_smooth(method='lm', aes(x=biomass_change, y=shan_change), se=F)

ggplot(data=eco_vars) + 
  geom_point( aes( x=biomass_converted, y=shan, col=rain_trt, shape=graze_trt) ) +
  geom_smooth(method='lm', aes(x=biomass_converted, y=shan), se=F)

ggplot(data=eco_vars_2020) + 
  geom_point( aes( x=biomass_converted, y=shan, col=rain_trt, shape=graze_trt) ) +
  geom_smooth(method='lm', aes(x=biomass_converted, y=shan), se=F)





#' *Figure*:  Raw responses of 2020 forage metrics 
#' 
#' Faceted by grazing
v1 <- ( fig_gsl | fig_peakg ) /  plot_spacer() / (fig_protein | fig_energy) /  plot_spacer() / ( fig_rfv | plot_spacer()) + plot_layout(heights = c(4.5, 0.25 , 4.5, 0.25, 4.5))
v1
#' Faceted by rain
v2 <- ( fig_gsl2 | fig_peakg2 ) /  plot_spacer() / (fig_protein2 | fig_energy2) /  plot_spacer() / ( fig_rfv2 | plot_spacer()) + plot_layout(heights = c(4.5, 0.25 , 4.5, 0.25, 4.5))
v2

tiff(filename="forage_service_2020_graze-facet_Jul24.tiff", res=600, width=8, height = 12, units = "in")
v1
dev.off()

tiff(filename="forage_service_2020_rain-facet.tiff", res=600, width=6, height = 8, units = "in")
v2
dev.off()



#' *Figure*:  Responses of all-year metrics
#' 
#' Faceted by grazing
v3 <-  fig_biomass + fig_bare + fig_div + guide_area() + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')') & 
    theme(plot.tag = element_text(size = 10))
v3 

tiff(filename="biomass_bare_div_graze-facet_Jun24.tiff", res=600, width=9, height = 7, units = "in")
v3 
dev.off()



v4 <-  (fig_biomass_vert / fig_bare_vert / fig_div)  + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')') & 
    theme(plot.tag = element_text(size = 10))
v4    

tiff(filename="biomass_bare_div_rain-facet_vert_Jul24.tiff", res=600, width=8, height = 7.5, units = "in")
v4 
dev.off()   





#' Faceted by rain
v5 <- fig_biomass2 + fig_div2 + fig_bare2 + guide_area() + plot_layout(guides = 'collect') + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
v5


#' All reported variables
fig_all <- (fig_biomass_vert / fig_bare_vert2 / fig_div / plot_spacer()) |
    (fig_ann_vert/ fig_forb_vert / fig_c3_vert / fig_c4) +
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')') & 
    theme(plot.tag = element_text(size = 10))
fig_all

tiff(filename="all_responses_rain_grazing.tiff", res=600, width=8, height = 10, units = "in")
fig_all 
dev.off()  



#' Save outliers list to summarize, and for future efforts:
outliers_list
write.csv(outliers_list, "Outliers_excluded_from_linear_models_2018-20.csv")




####
#'    
#'    [D] : **Identify key indices of plant functions and services** 
#'    
#'                  PCA of response variables in 2020
#'    
#'    
#'###


#'    
#' * Ordination of forage service variables **
#'    

#' View structure of dataframe with 8 variables used for ordination
#'
#' This already has outliers removed from analyses above 
#' (bare ground, 1 outlier; max_gi_doy/peakg, 5 outliers)
str(eco_vars_out_2020)
eco_vars_out_2020 %>% select(biomass_converted:shan) %>% ggpairs()
outliers_list

#' *Substitutions*
#' 
#' Note that PCA does not accept missing values, so for PCA only, we will try
#' to substitute some missing forage service variables with reasonable approximations.
#' 
#' View all rows with some missing data
eco_vars_out_2020[!complete.cases(eco_vars_out_2020), c(3:13)]

#' For crude protein and digestible energy, several plots are missing because we had 
#' to pool low-biomass plots for analysis, and assigned the data to just one plot. Here 
#' we will re-integrate the duplicated values for both plots:
#'    1-3 now also used for 1-7
#'    3-2 now also used for 3-3
#'    3-7 now also used for 3-4 (but this has other missing data and might be removed)
#'    5-5 now also used for 5-1
#'    5-3 now also used for 5-7
#'    
fq_1_7 <- eco_vars_out_2020$b_p %in% c("1_7")
eco_vars_out_2020$crude_protein[fq_1_7] <- eco_vars_out_2020 %>% filter( b_p == "1_3") %>% select(crude_protein) %>% c(1)[1]
eco_vars_out_2020$digestible_energy[fq_1_7] <- eco_vars_out_2020 %>% filter( b_p == "1_3") %>% select(digestible_energy) %>% c(1)[1]
eco_vars_out_2020$relative_feed_value[fq_1_7] <- eco_vars_out_2020 %>% filter( b_p == "1_3") %>% select(relative_feed_value) %>% c(1)[1]

fq_3_3 <- eco_vars_out_2020$b_p %in% c("3_3")
eco_vars_out_2020$crude_protein[fq_3_3] <- eco_vars_out_2020 %>% filter( b_p == "3_2") %>% select(crude_protein) %>% c(1)[1]
eco_vars_out_2020$digestible_energy[fq_3_3] <- eco_vars_out_2020 %>% filter( b_p == "3_2") %>% select(digestible_energy) %>% c(1)[1]
eco_vars_out_2020$relative_feed_value[fq_3_3] <- eco_vars_out_2020 %>% filter( b_p == "3_2") %>% select(relative_feed_value) %>% c(1)[1]

fq_3_4 <- eco_vars_out_2020$b_p %in% c("3_4")
eco_vars_out_2020$crude_protein[fq_3_4] <- eco_vars_out_2020 %>% filter( b_p == "3_7") %>% select(crude_protein) %>% c(1)[1]
eco_vars_out_2020$digestible_energy[fq_3_4] <- eco_vars_out_2020 %>% filter( b_p == "3_7") %>% select(digestible_energy) %>% c(1)[1]
eco_vars_out_2020$relative_feed_value[fq_3_4] <- eco_vars_out_2020 %>% filter( b_p == "3_7") %>% select(relative_feed_value) %>% c(1)[1]

fq_5_1 <- eco_vars_out_2020$b_p %in% c("5_1")
eco_vars_out_2020$crude_protein[fq_5_1] <- eco_vars_out_2020 %>% filter( b_p == "5_5") %>% select(crude_protein) %>% c(1)[1]
eco_vars_out_2020$digestible_energy[fq_5_1] <- eco_vars_out_2020 %>% filter( b_p == "5_5") %>% select(digestible_energy) %>% c(1)[1]
eco_vars_out_2020$relative_feed_value[fq_5_1] <- eco_vars_out_2020 %>% filter( b_p == "5_5") %>% select(relative_feed_value) %>% c(1)[1]

fq_5_7 <- eco_vars_out_2020$b_p %in% c("5_7")
eco_vars_out_2020$crude_protein[fq_5_7] <- eco_vars_out_2020 %>% filter( b_p == "5_3") %>% select(crude_protein) %>% c(1)[1]
eco_vars_out_2020$digestible_energy[fq_5_7] <- eco_vars_out_2020 %>% filter( b_p == "5_3") %>% select(digestible_energy) %>% c(1)[1]
eco_vars_out_2020$relative_feed_value[fq_5_7] <- eco_vars_out_2020 %>% filter( b_p == "5_3") %>% select(relative_feed_value) %>% c(1)[1]

#' 
#' For one missing bare ground value (4-8), we used total cover (separate metric) to predict the missing value (R2 = 0.578).
#' bare = (-0.435*totalcov) + 57.551 = 21.64205894
#' 
fq_4_8 <- eco_vars_out_2020$b_p %in% c("4_8")
eco_vars_out_2020$bare_visual_cov[fq_4_8] <- 21.64205894


#' 
#' For missing Growing Season Length values (2-10, 3-4), we used Total Greeness (total area under greenness curves)
#' to predict missing GSL values. Total greeness was estimatable for all 72 plots, and was
#' relate to GSL (R2 = 0.6319; GSL = 4.5351(totG) - 34.037
#' 
#' For one missing bare ground value (4-11), we used total cover to predict the missing value (R2 = 0.578).
#' bare = (-0.435*totalcov) + 57.551 = 21.64205894

#' 
#' The remaining plots were removed above as outliers, so we will keep them out of the ordination
#' 



#' PCA dataframe - with Shannon Diversity
eco_vars_pca_dat <- eco_vars_out_2020 %>%  select(biomass_converted:shan) 
eco_vars_pca_dat <- eco_vars_pca_dat %>% rename(Standing_Biomass=biomass_converted , Bare_Soil = bare_visual_cov,
                                    Crude_Protein = crude_protein, Digestible_Energy = digestible_energy,
                                    Peak_Greeness = max_gi_doy, Growing_Season_Len=gs_length, RFV = relative_feed_value, Diversity = shan)

#' Scale the data (but not centered around 0, to maintain positive data structure)
eco_vars_pca_dat_s <- data.frame(scale(eco_vars_pca_dat, center=FALSE))



#' 
#' *Get initial correlations* across forage service variables
#' 
eco_vars_corrs <- data.frame(round(cor(eco_vars_pca_dat_s , use="complete.obs"),3))
ggcorrplot(eco_vars_corrs, type="lower",lab="T", colors=c("navyblue","white","goldenrod4"))
# Save correlation plot
tiff(filename="forage_service_correlation_table_2020.tiff", res=300, width=7, height = 7, units = "in")
ggcorrplot(eco_vars_corrs, type="lower",lab="T", colors=c("navyblue","white","goldenrod4"))
dev.off()



#' 
#' *Conduct PCA*
#'
#' Reduce to complete observations
eco_vars_pca_dat_complete <- eco_vars_out_2020[complete.cases(eco_vars_out_2020), c(6:13)]
eco_vars_pca_dat_complete %>% ggpairs(use="complete.obs")
eco_vars_pca_dat_complete_ALL <- eco_vars_out_2020[complete.cases(eco_vars_out_2020), c(1:13)]
#' PCA
eco_vars_pca <- vegan::rda(eco_vars_pca_dat_complete, scale=T)
summary(eco_vars_pca)
rda_plot <- plot(eco_vars_pca)


#'
#' *Save PCA and supp data for figure*
#' 
rda_scores <- data.frame(rda_plot$sites)  #Site scores
rownames(rda_scores) <- eco_vars_pca_dat_complete_ALL$b_p
rda_scores$b_p <- eco_vars_pca_dat_complete_ALL$b_p
rda_scores <- left_join(rda_scores, eco_vars_pca_dat_complete_ALL, by='b_p')
rda_scores$rain_trt <- factor(rda_scores$rain_trt, levels=c("D","C","W"))

summary(aov(PC1 ~ rain_trt * graze_trt, data=rda_scores))
summary(aov(PC2 ~ rain_trt * graze_trt, data=rda_scores))

#'   Create forage service vectors dataframe
rda_vecs <- data.frame(rda_plot$species)  #Env vectors
rownames(rda_vecs) <- c("End of Season\nBiomass","Bare Soil", "Crude Protein", "Dig.Energy",
                        "RFV", "Peak Green", "GSL","Shannon Diversity")

#' Add forage service CORRELATIONS with each axis to dataframe (in addition to raw loadings)
pc_corrs_dat <- rda_scores %>% select(PC1, PC2, biomass_converted:shan)
pc_corrs <- data.frame(cor(pc_corrs_dat))
pc_corrs <- pc_corrs %>% select (PC1, PC2) %>% filter (!row.names(pc_corrs) %in% c('PC1','PC2'))
pc_corrs
rda_vecs$PC1_corr <- pc_corrs$PC1
rda_vecs$PC2_corr <- pc_corrs$PC2
rda_vecs$xlab_corr <- c(0.2,1.1,1.19,0.45,0.95,-0.15,-.74,1.0)
rda_vecs$ylab_corr <- c(-1.06,0.3,0.17,-0.73,-0.6,-0.54,-0.43,0.07)

#' Create a second vector dataframe with just diversity and plant production to highlight those variables
rda_vecs2 <- rda_vecs[c(1,8),]

#' Fix coding for graze treatment
rda_scores$graze_trt <- gsub("C", "Ungrazed", rda_scores$graze_trt)
rda_scores$graze_trt <- gsub("FG", "Dormant-Season Grazing", rda_scores$graze_trt)
rda_scores$graze_trt <- gsub("SG", "Growing-Season Grazing", rda_scores$graze_trt)
rda_scores$graze_trt <- factor(rda_scores$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))

rda_scores_gs <- rda_scores %>% filter(graze_trt == "Growing-Season Grazing")

#' FIGURE
#'
#'   Grazing and Rainfall  emphasis on PCA of forage variables
#'   
pca_graze_fig <- ggplot () +
  geom_point(aes(x=PC1, y=PC2, color=rain_trt, shape=graze_trt),  alpha=0.6, cex=2.7, data = rda_scores) +
  geom_segment(mapping=aes(x=0, y=0, xend=PC1_corr, yend=PC2_corr), size=1, data=rda_vecs) +
  geom_segment(mapping=aes(x=0, y=0, xend=PC1_corr, yend=PC2_corr), size=1, col='red4', data=rda_vecs2) +
  geom_text(mapping=aes(x=xlab_corr, y=ylab_corr, label = rownames(rda_vecs)),cex=3, data=rda_vecs) +
  geom_text(mapping=aes(x=xlab_corr, y=ylab_corr, label = rownames(rda_vecs2)),cex=3,col='red4',data=rda_vecs2) +
  scale_color_manual(values=c("goldenrod2","forestgreen","navyblue")) +
  #scale_shape_manual (values=c(16,2,0)) +
  #stat_ellipse(aes(x=PC1, y=PC2, lty=graze_trt),level=0.95, type = "t", data=rda_scores) +
  #stat_ellipse(aes(x=PC1, y=PC2, color=rain_trt),level=0.95, type = "t", data=rda_scores) +
  labs(col="Rain trt.", shape="Graze trt.", lty='Graze trt.', x="PC1 - 28.7%", y='PC2 - 24.4%') +
  ylim(-1.75,1.75) + xlim(-2,2) +
  geom_text(aes(x=-0.9,y=1.7, label="Rangeland Functions in 2020:\nTradeoffs and Synergies"), fontface=2, cex=3.5) +
  scale_alpha_manual(guide = guide_legend(override.aes = list( color = "gray") ) ) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
pca_graze_fig

tiff(filename="forage_services_PCA_Jun24_complete_cases.tiff", res=300, width=6.5, height = 4.75, units = "in")
pca_graze_fig
dev.off()



#' **Analysis: Linear mixed model of PC scores**
#' 

#' PC 1
#' 
pc1_mod <- lmer(PC1 ~ rain_trt * graze_trt + (1|block), data=rda_scores)

# Diagnostics
plot(pc1_mod)
qqnorm(residuals(pc1_mod))

# Results
anova(pc1_mod)
summary(pc1_mod)
performance(pc1_mod)



#' PC 2
#' 
pc2_mod <- lmer(PC2 ~ rain_trt * graze_trt + (1|block), data=rda_scores)

# Diagnostics
plot(pc2_mod)
qqnorm(residuals(pc2_mod))

# Results
anova(pc2_mod)
summary(pc2_mod)
performance(pc2_mod)






#####
#'   
#'   
#'   
#'    [Part 3] : SYNTHESIS: EARLY INDICATORS OF DROUGHT RESISTENCE 
#'    
#'
#'   *Community composition and grazing as mediators of forage service response to drought*
#'   
#'   
#'   [A] : Explore drought trajectories of two key response variables: diversity and biomass
#'   
#'   [B] : Merge composition and forage dataframes for analysis
#'   
#'   [C] : Models and figures of drought resistance
#'   
#'   
#'   
#'   

#'    
#'   [A] ** Movement trajectories of drought plots from 2018 to 2020 along two key forage axes: **
#'                       
#'                Diversity and Plant Production
#'    

#' PCA suggests that diversity and plant biomass are not strongly correlated (r=-0.03),
#'      with each representing two main sets of correlated forage services:
#'  
#' Along Axis 1, species diversity loaded intermediately and captured a tradeoff between
#'          high levels of diversity, crude protein and bare soil VS a longer growing season
#' Along Axis 2, plant biomass production loaded strongly and captured a suite of positive
#'          associations between high biomass, relative feed value, and digestible energy
#'          
#'
#' Since both plant production and diversity are available for 3 years with full datasets 
#' (bare ground has missing data points), we will use these as proxies to quantify plot-level 
#' change from 2018 to 2020 in core forage metrics (last part of analysis)
#' 
#' 
#' Here we will make a figure showing the trajectories of drought plots as they changed in
#' plant biomass and diversity over time.

#' Outliers
#' 
#' Note that as far as outliers go, we didn't have any plots to remove in terms of diversity
#' or plant standing biomass... however, many plots that 3-4, 3-6, 4-11, 2-10, 6-3 are outliers
#' for other variables... but only one of those was an outlier because of extreme differences
#' that could affect things (rather than simply not being able to estimate greenness values)..
#' 
#' This is plot 3-4 with lots of annuals. So, we might want to check it.
outliers_list

#' Create dataframe with separate biomass diversity columns for 2020 and 2018 for vectors
str(eco_vars)
biomass_div <- eco_vars %>% filter(year == "2020") %>% filter( rain_trt =="D") %>% 
  select (b_p, rain_trt, graze_trt, biomass_converted, shan)
biomass_div_18 <- eco_vars %>% filter (year == "2018") %>% filter (rain_trt == "D") %>% select (b_p, biomass_converted, shan)
biomass_div$biomass_converted_18 <- biomass_div_18$biomass_converted
biomass_div$shan_18 <- biomass_div_18$shan

# Fix coding for graze treatment
biomass_div$graze_trt <- gsub("C", "Ungrazed", biomass_div$graze_trt)
biomass_div$graze_trt <- gsub("FG", "Fall Graze", biomass_div$graze_trt)
biomass_div$graze_trt <- gsub("SG", "Spring Graze", biomass_div$graze_trt)
biomass_div$graze_trt <- factor(biomass_div$graze_trt, levels=c("Spring Graze","Fall Graze","Ungrazed"))


#' Create a dataframe with raw values for 2020 and 2018 only for points
biomass_div_18_20 <- eco_vars %>% filter (year != "2019") %>% filter(rain_trt =="D")
# Fix coding for graze treatment
biomass_div_18_20$graze_trt <- gsub("C", "Ungrazed", biomass_div_18_20$graze_trt)
biomass_div_18_20$graze_trt <- gsub("FG", "Fall Graze", biomass_div_18_20$graze_trt)
biomass_div_18_20$graze_trt <- gsub("SG", "Spring Graze", biomass_div_18_20$graze_trt)
biomass_div_18_20$graze_trt <- factor(biomass_div_18_20$graze_trt, levels=c("Spring Graze","Fall Graze","Ungrazed"))


#' *Figure* showing trajectories in biomass and diversity space
#' 
biomass_div_change_fig <- ggplot () + 
  geom_point(aes(x=biomass_converted, y=shan, shape=graze_trt, alpha=year), cex=3, color="goldenrod3", data=biomass_div_18_20) +
  geom_segment(aes(x=biomass_converted_18, y=shan_18, xend=biomass_converted, yend=shan,lty=graze_trt), data=biomass_div) +
  scale_alpha_manual(values=c(0.35,0.90))+
  ylim(2.2,3.2)+
  xlim(0,300)+
  labs(x="Standing Plant Biomass (g per m2)", y="Shannon Diversity", shape="Graze trt.",lty="Graze trt.", alpha="Year") + 
  theme_bw()
biomass_div_change_fig 

tiff(filename="drought_change_biplot_diversity_biomass.tiff", res=600, width=6, height = 5, units = "in")
biomass_div_change_fig 
dev.off()


#' Save figure alongside PCA
pca_graze_fig / biomass_div_change_fig + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')


tiff(filename="PCA_and_DroughtResponse_Biplot.tiff", res=600, width=7, height = 10, units = "in")
pca_graze_fig / biomass_div_change_fig + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
dev.off()




#'###
#' 
#' 
#' 
#'   [B] : Merge composition and forage dataframes for analysis        
#'
#'        Combine **Community metricS (FG, rich, even, div)** with 
#'                **rangeland functions** in a single dataframe
#'  
#'###



#' ###
#' 
#'  *Merge datasets*
#'  
#' ### 

#' Functional groups (remove 'G' which has no attributed cover)
fg_wide_sqrt <- spp_df_rel_sqrt_long_fg %>% filter(fun_grp != "G_NA") %>% 
  filter(fun_grp != "F_NA") %>% filter(fun_grp != "S") %>% spread(fun_grp, fg_cov) %>%
  select(year, b_p, F, G_C3, G_C4, Ann_Bi)
str(fg_wide_sqrt, give.attr=F)



#' Join with biomass and diversity variables
eco_vars_resist <- eco_vars %>% select (block, b_p, year, rain_trt, graze_trt, biomass_converted, shan)
comm_join2 <- left_join(eco_vars_resist, fg_wide_sqrt, by=c('b_p','year'))
comm_join2$plot <- sub(".*_", "", comm_join2$b_p)
#' Add in major species relative abundances
str(major_spp)
major_spp$b_p <- as.factor(major_spp$b_p)
major_spp2 <- major_spp %>% select(-block,-plot,-rain_trt,-graze_trt)
comm_join <- left_join(comm_join2, major_spp2, by=c('b_p','year'))
#' Check raw dataframe
str(comm_join, give.attr=FALSE)


#' 
#' *Explre initial correlations* across community variables
#' 
#' 
comm_pca_dat <- comm_join %>% ungroup() %>% rename(Forbs=F, Annual_Biennial=Ann_Bi, C3_Grass_Perennial=G_C3,
                                                   C4_Grass_Perennial=G_C4, Diversity=shan)
comm_pca_dat_all <-  comm_pca_dat %>% select(Forbs:Annual_Biennial)
comm_pca_2018_dat <- comm_pca_dat %>% filter (year == "2018") %>% select(Forbs:Annual_Biennial)
comm_pca_2020_dat <-comm_pca_dat %>% filter (year == "2020") %>% select(Forbs:Annual_Biennial)
comm_pca_2018_dry_dat <- comm_pca_dat %>% filter (year == "2018") %>% filter (rain_trt == "D") %>% select(Forbs:Annual_Biennial)
comm_pca_2020_dry_dat <-comm_pca_dat %>% filter (year == "2020") %>% filter (rain_trt == "D") %>% select(Forbs:Annual_Biennial)

# All years
comm_corrs <- data.frame(round(cor(comm_pca_dat_all, use="pairwise.complete.obs"),3))
pairs(comm_pca_dat_all)
ggcorrplot(comm_corrs, type="lower",lab="T", colors=c("navyblue","white","goldenrod4"))

# 2018
comm_corrs_2018 <- data.frame(round(cor(comm_pca_2018_dat, use="pairwise.complete.obs"),3))
pairs(comm_pca_2018_dat)
ggcorrplot(comm_corrs_2018, type="lower",lab="T", colors=c("navyblue","white","goldenrod4"))
# 2020
comm_corrs_2020 <- data.frame(round(cor(comm_pca_2020_dat, use="pairwise.complete.obs"),3))
pairs(comm_pca_2020_dat)
ggcorrplot(comm_corrs_2020, type="lower",lab="T", colors=c("navyblue","white","goldenrod4"))
# 2018 DRY only
comm_corrs_2018_dry <- data.frame(round(cor(comm_pca_2018_dry_dat, use="pairwise.complete.obs"),3))
pairs(comm_pca_2018_dry_dat)
ggcorrplot(comm_corrs_2018_dry, type="lower",lab="T", colors=c("navyblue","white","goldenrod4"))
# 2020 DRY only
comm_corrs_2020_dry <- data.frame(round(cor(comm_pca_2020_dry_dat, use="pairwise.complete.obs"),3))
pairs(comm_pca_2020_dry_dat)
ggcorrplot(comm_corrs_2020_dry, type="lower",lab="T", colors=c("navyblue","white","goldenrod4"))

#' ggpairs figure across all years
ggpairs(comm_pca_dat, columns = 8:11, ggplot2::aes(colour = year,alpha=0.3))




#' ####
#' 
#' 
#' 
#'   [C] : Models and figures of drought resistance
#' 
#' 
#' ####


#' *Check dataset structure*, make sure it is arranged identically across years
str(comm_join, give.attr=F)
comm_join <- comm_join %>% arrange (year, block, plot)


#' *Create a 'change' dataframe to explore drought resistance*
#' 
#' Goal: identical dataframes for 2018 and 2020, including only variables of interest to
#' assess change over time: 
#'  
#'    Key drought responses: biomass_converted, bare_visual_cov, total_cov
#'    Key change predictors: F, G_Ann, G_C3, G_C4, rich, even
#'    
resist_dat_2018 <- comm_join %>% filter( year=='2018') %>% ungroup() %>%
  select(F, Ann_Bi, G_C3, G_C4, mm, ag, cg, vo, td, bj, shan, biomass_converted)
resist_dat_2020 <- comm_join %>% filter( year=='2020') %>% ungroup() %>%
  select(F, Ann_Bi, G_C3, G_C4, mm, ag, cg, vo, td, bj, shan, biomass_converted)
resist_dat_bp <- comm_join %>% filter( year=='2020') %>% ungroup() %>%
    select(b_p)
    
#' Create the 'change' dataframe 
resist_dat_change_LRR <- log(resist_dat_2020 / resist_dat_2018) 

#' Change column names from 2018 and 2020 to be distinguishable when merging
colnames(resist_dat_change_LRR) <- paste(colnames(resist_dat_change_LRR),"LRR",sep="_")
colnames(resist_dat_2018) <- paste(colnames(resist_dat_2018),"2018",sep="_")
colnames(resist_dat_2020)  <- paste(colnames(resist_dat_2020),"2020",sep="_")

#' Reduce 2020 columns to just raw shan and biomass data
resist_dat_2020 <- resist_dat_2020 %>% select(shan_2020, biomass_converted_2020)
resist_dat_2020$b_p <- resist_dat_bp$b_p

#' MERGE 'CHANGE' variables with treatment info and predictors from 2018:
comm_join_info <- comm_join %>% ungroup() %>% filter (year=="2018") %>% select(block, plot, rain_trt, graze_trt)
resist_dat3 <- cbind(comm_join_info, resist_dat_2018, resist_dat_change_LRR)
resist_dat3$b_p <- paste(resist_dat3$block, resist_dat3$plot, sep="_")
resist_dat2 <- left_join (resist_dat3, fg_abs_18, by="b_p")
resist_dat <- left_join(resist_dat2, resist_dat_2020, by="b_p")

#' Filter to *drought only*
resist_dat_dry <- resist_dat %>% filter (rain_trt == "D")
str(resist_dat_dry, give.attr=F)

#' Add raw differences in biomass and shan to check against LRR
#' Positively correlated... the lower the LRR, the more that a metric decreased from 2018 to 2020
resist_dat_dry$biomass_diff <- resist_dat_dry$biomass_converted_2020- resist_dat_dry$biomass_converted_2018
resist_dat_dry$shan_diff <- resist_dat_dry$shan_2020- resist_dat_dry$shan_2018
plot(resist_dat_dry$biomass_diff, resist_dat_dry$biomass_converted_LRR)
plot(resist_dat_dry$shan_diff, resist_dat_dry$shan_LRR)

#' ** Check removing the one annual plot outlier, 3-4**
resist_dat_dry_out <- resist_dat_dry %>% filter(b_p != "3_4")



#'
#' *Check distributions of model variables*
#' 
#'  
resist_dat_dry_dist <- resist_dat_dry %>% ungroup() %>% select(F_2018:bj_2018,shan_LRR:biomass_converted_LRR)
ggpairs(resist_dat_dry_dist)

resist_dat_dry_out_dist <- resist_dat_dry_out %>% ungroup() %>% select(F_2018:cg_2018,shan_LRR:biomass_converted_LRR)
ggpairs(resist_dat_dry_out_dist)

#' Check log-transformation of one predictor variables:  annual-biennial cover
resist_dat_dry_dist$ln_Ann_Bi_2018 <- log(resist_dat_dry_dist$Ann_Bi_2018)
ggpairs(resist_dat_dry_dist)
resist_dat_dry_out_dist$ln_Ann_Bi_2018 <- log(resist_dat_dry_out_dist$Ann_Bi_2018)
ggpairs(resist_dat_dry_out_dist)

#' Log-transformation seems preferable for *annual grass cover*
#'     Add to the main dataframe
resist_dat_dry$log_Ann_Bi_2018  <- log(resist_dat_dry_dist$Ann_Bi_2018)
resist_dat_dry_out$log_Ann_Bi_2018  <- log(resist_dat_dry_out_dist$Ann_Bi_2018) 


#'
#' Scale 2018 biomass, bare, and total cover values for modeling
#'
resist_dat_dry$biomass_converted_2018_s <- scale(resist_dat_dry$biomass_converted_2018)
resist_dat_dry$shan_2018_s <- scale(resist_dat_dry$shan_2018)

resist_dat_dry_out$biomass_converted_2018_s <- scale(resist_dat_dry_out$biomass_converted_2018)
resist_dat_dry_out$shan_2018_s <- scale(resist_dat_dry_out$shan_2018)

#' 
#' *Set contrasts for models*
#'
resist_dat_dry$graze_trt <- factor(resist_dat_dry$graze_trt, levels=c("C","SG","FG"))
resist_dat_dry_out$graze_trt <- factor(resist_dat_dry_out$graze_trt, levels=c("C","SG","FG"))

#'
#' * Check final dataframe*
#' 
str(resist_dat_dry, give.attr=F)
str(resist_dat_dry_out, give.attr=F)



#'  *Question: Is an increase in Shannon diversity linked to an increase in annual-biennials?*
#'  
#'  
ggplot(data=resist_dat_dry, aes(x=Ann_Bi_LRR, y=shan_LRR, color=ag_2018)) +
    geom_point(cex=3) + 
    geom_smooth(method="lm", se=F, color='black') + 
        theme_bw()


ggplot(data=resist_dat_dry, aes(x=Ann_Bi_LRR, y=shan_LRR, color=ag_2018)) +
    geom_point(cex=3) + 
    geom_smooth(method="lm", se=F, color='black') + 
    facet_wrap(~graze_trt) + 
    theme_bw()


ggplot(data=resist_dat_dry, aes(y=Ann_Bi_LRR, x=ag_2018)) +
    geom_point(cex=3) + 
    geom_smooth(method="lm", se=F, color='black') + 
    theme_bw()

ggplot(data=resist_dat_dry, aes(y=Ann_Bi_LRR, x=ag_2018)) +
    geom_point(cex=3) + 
    geom_smooth(method="lm", se=F, color='black') + 
    facet_wrap(~graze_trt)+
    theme_bw()


#'  #########
#'  
#'  Drought resistance:  **CORRELATIONS**  
#'  
#'  #########
#'

resist_dat_dry %>% select( shan_LRR:biomass_converted_LRR) %>% ggpairs()




#'  #########
#'
#'  NOTE on *Linear Mixed Effects Models*
#' 
#'  #########
#'    
#'    Tried running linear mixed effects models first including block as a random effect
#'    
#'    Everything is singular, **AND** because we include starting biomass in the model, block isn't
#'    really a necessary grouping variable anyway -- we're already controlling for outside influences
#'    in a much more specific way.
#'    
#'    So, no random effects will be included in these models.
#'    

#'  
#'  ** Figure / Data prep**
#'  

#' 
#' *create dataframes for figures with and without the potential outlier plot* 
#' 
#' Reset figure contrasts to use Spring graze as baseline
resist_dat_dry_fig <- resist_dat_dry
resist_dat_dry_fig$graze_trt <- gsub("FG", "Dormant-Season Grazing", resist_dat_dry_fig$graze_trt)
resist_dat_dry_fig$graze_trt <- gsub("SG", "Growing-Season Grazing", resist_dat_dry_fig$graze_trt)
resist_dat_dry_fig$graze_trt <- gsub("C", "Ungrazed", resist_dat_dry_fig$graze_trt)
resist_dat_dry_fig$graze_trt <- factor(resist_dat_dry_fig$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))

resist_dat_dry_fig_out <- resist_dat_dry_out
resist_dat_dry_fig_out$graze_trt <- gsub("FG", "Dormant-Season Grazing", resist_dat_dry_fig_out$graze_trt)
resist_dat_dry_fig_out$graze_trt <- gsub("SG", "Growing-Season Grazing", resist_dat_dry_fig_out$graze_trt)
resist_dat_dry_fig_out$graze_trt <- gsub("C", "Ungrazed", resist_dat_dry_fig_out$graze_trt)
resist_dat_dry_fig_out$graze_trt <- factor(resist_dat_dry_fig_out$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))

#'  *Biomass ~ Shan Div  FIGURE*
#'  
#'   Showing the two response metrics we will look at with models below 
#'  
#' Overall biomass and shannon diversity figure
biomass_shan_resist <- ggplot( data=resist_dat_dry_fig, aes(x=biomass_converted_LRR, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_vline(aes(xintercept=0), col='gray80', size=1)+
    geom_point(fill='goldenrod3', cex=3,alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    labs(y="LRR - Shannon Diversity\nless resistant                           more resistant", 
         x="LRR - End of Season Biomass\nless resistant                           more resistant", shape="Graze trt.") +
    theme_bw()
biomass_shan_resist

tiff(filename="Resistance_biomass_shan.tiff", res=600, width=5.5, height = 3.5, units = "in")
biomass_shan_resist
dev.off()


#'  #########
#'  
#'  Drought resistance:  **BIOMASS PRODUCTION**  
#'  
#'  #########

#' *Linear Models* 
#' 
#'    When all data are retained, only C4 annual grasses have added effects on top of grazing
#'   
#'    Outlier removal:
#'    We look at what happens when removing one plot (3-4) that had outlying values for
#'    annual-biennial cover in the mixed model for annual-biennial response above.
#'    To be clear, there are several other plots with similarly high cover of annual-biennial grasses,
#'    however, I believe this is one of the only that occurred in the spring graze treatment.
#'    
#'    In the model above, plot 3-4 caused the greatest amount of skew, and removing it improved 
#'    diagnostic residual plots.
#'    
#'    Below, we are running completely different models, and while we look at removing plot 3-4,
#'    making a decision is tricky.
#'    
#'    On the one hand, whether or not this plot is included DOES seem to affect outcomes. 
#'    We have created figures to show this. 
#'    
#'    On the other hand, none of the diagnostic or residual plots indicate a problem with
#'    keeping this 'outier' in the first place... it's not necessarily an 'outlier' in 
#'    these models. What it does indicate is that with our small sample sizes, plots that have
#'    more extreme cover values (including other plots that were NOT removed) tend to have
#'    an outsized influence on relationships.
#'    
#' 
#' *FORBS  --  GRAZING ONLY *
#' 
#' Highly sensitive to removal of plot 3_4
#' 
#' With outlier removed, F p=0.007 (positive effect of F on resistance)
#' 
resist_bm_f <- lm(biomass_converted_LRR ~ graze_trt * F_2018 + biomass_converted_2018_s, data = resist_dat_dry) 
resist_bm_f_out <- lm(biomass_converted_LRR ~ graze_trt * F_2018 + biomass_converted_2018_s, data = resist_dat_dry_out) 

summary(resist_bm_f)
summary(resist_bm_f_out)
anova(resist_bm_f)
anova(resist_bm_f_out)

plot(resist_bm_f, which=c(1,1),ask=F)
plot(resist_bm_f_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_bm_f))
qqnorm(residuals(resist_bm_f_out))


#' *G_Ann  --   GRAZING ONLY*
#' 
#' Highly sensitive to removal of plot 3_4
#' 
#' With outlier removed, graze*Ann p=0.09 (lower ann-bi cover = more resistance with spring grazing)
#' 

resist_bm_ann <- lm(biomass_converted_LRR ~ graze_trt * log_Ann_Bi_2018 + biomass_converted_2018_s, data = resist_dat_dry) 
resist_bm_ann_out <- lm(biomass_converted_LRR ~ graze_trt * log_Ann_Bi_2018 + biomass_converted_2018_s, data = resist_dat_dry_out) 

summary(resist_bm_ann)
summary(resist_bm_ann_out)
anova(resist_bm_ann)
anova(resist_bm_ann_out)

plot(resist_bm_ann, which=c(1,1),ask=F)
plot(resist_bm_ann_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_bm_ann))
qqnorm(residuals(resist_bm_ann_out))


#' *G_C3  --  GRAZING ONLY*
resist_bm_c3 <- lm(biomass_converted_LRR ~ graze_trt * G_C3_2018 + biomass_converted_2018_s, data = resist_dat_dry) 
resist_bm_c3_out <- lm(biomass_converted_LRR ~ graze_trt * G_C3_2018 + biomass_converted_2018_s, data = resist_dat_dry_out)

summary(resist_bm_c3)
summary(resist_bm_c3_out)
anova(resist_bm_c3)
anova(resist_bm_c3_out)

plot(resist_bm_c3, which=c(1,1),ask=F)
plot(resist_bm_c3_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_bm_c3))
qqnorm(residuals(resist_bm_c3_out))


#' *G_C4  --    GRAZING & WEAK INTERACTION*
#' 
#' Highly sensitive to removal of plot 3_4
#' 
#' With outlier removed, interaction is lost and there is an overall negative effect of C4 grass cover
#' 
resist_bm_c4 <- lm(biomass_converted_LRR  ~ graze_trt * G_C4_2018 + as.vector(biomass_converted_2018_s), data = resist_dat_dry) 
resist_bm_c4_out <- lm(biomass_converted_LRR  ~ graze_trt * G_C4_2018 + biomass_converted_2018_s, data = resist_dat_dry_out) 

summary(resist_bm_c4)
summary(resist_bm_c4_out)
anova(resist_bm_c4)
anova(resist_bm_c4_out)

plot(resist_bm_c4, which=c(1,1),ask=F)
plot(resist_bm_c4_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_bm_c4))
qqnorm(residuals(resist_bm_c4_out))


#' *Andger -- GRAZING ONLY*
resist_bm_ag <- lm(biomass_converted_LRR ~ graze_trt * ag_2018+ biomass_converted_2018_s, data = resist_dat_dry) 
resist_bm_ag_out <- lm(biomass_converted_LRR ~ graze_trt * ag_2018+ biomass_converted_2018_s, data = resist_dat_dry_out) 
summary(resist_bm_ag)
summary(resist_bm_ag_out)
anova(resist_bm_ag)
anova(resist_bm_ag_out)
plot(resist_bm_ag, which=c(1,1),ask=F)
plot(resist_bm_ag_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_bm_ag))
qqnorm(residuals(resist_bm_ag_out))

#'*Chogra -- GRAZING ONLY*
resist_bm_cg_out <- lm(biomass_converted_LRR ~ graze_trt * cg_2018+ biomass_converted_2018_s, data = resist_dat_dry_out) 
summary(resist_bm_cg)
summary(resist_bm_cg_out)
anova(resist_bm_cg)
anova(resist_bm_cg_out)
plot(resist_bm_cg, which=c(1,1),ask=F)
plot(resist_bm_cg_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_bm_cg))
qqnorm(residuals(resist_bm_cg_out))

#'*Muhmon -- GRAZING** Biomass** Graze:MM* *
resist_bm_mm <- lm(biomass_converted_LRR ~ graze_trt * mm_2018+ biomass_converted_2018_s, data = resist_dat_dry) 
resist_bm_mm_out <- lm(biomass_converted_LRR ~ graze_trt * mm_2018+ biomass_converted_2018_s, data = resist_dat_dry_out) 
summary(resist_bm_mm)
summary(resist_bm_mm_out)
anova(resist_bm_mm)
anova(resist_bm_mm_out)
plot(resist_bm_mm, which=c(1,1),ask=F)
plot(resist_bm_mm_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_bm_mm))
qqnorm(residuals(resist_bm_mm_out))




#' 
#' **Biomass resistance FIGURES**
#'
#


#' Forb figure
biomass_resist_fig_f <- ggplot( data=resist_dat_dry_fig, aes( x=F_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F, col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-2.3,1.3)+
    labs(x="Forbs\nInitial relative cover (%)", y="LRR - Standing Biomass\nless resistant                more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.35, y=1, label="Graze**\nInitial biomass**"), size=3) + 
    theme_bw()
biomass_resist_fig_f_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=F_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F,col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-2.3,1.3)+
    labs(x="Forbs\nInitial relative cover (%)", y="LRR - Standing Biomass\nless resistant                more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.35, y=1, label="Graze*\nForb**`\nInitial biomass*"), size=3) + 
    theme_bw()

biomass_resist_fig_f
biomass_resist_fig_f_out

#' Ann-bi figure
biomass_resist_fig_ann <- ggplot( data=resist_dat_dry_fig, aes( x=log_Ann_Bi_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F, col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    xlim(-5.2,-1.5)+ylim(-2.3,1.3) +
    labs(x="Annual-Biennial\nInitial relative cover (%)", y="LRR - Standing Biomass\nless resistant                more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=-2, y=1, label="Graze*\nInitial biomass**"), size=3) + 
    theme_bw()

biomass_resist_fig_ann_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=log_Ann_Bi_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F, col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    xlim(-5.2,-1.5)+ylim(-2.3,1.3) +
    labs(x="Annual-Biennial\nInitial relative cover (%)", y="LRR - Standing Biomass\nless resistant                more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=-2, y=1, label="Graze*\nGraze : Annual-Biennial`\nInitial biomass**"), size=3) + 
    theme_bw()

biomass_resist_fig_ann
biomass_resist_fig_ann_out 

#' C3 grass figure
biomass_resist_fig_c3 <- ggplot( data=resist_dat_dry_fig, aes( x=G_C3_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F, col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-2.3,1.3)+
    labs(x="C3 Perennial Grasses\nInitial relative cover (%)", y="LRR - Standing Biomass\nless resistant                more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.35, y=1, label="Graze*\nInitial biomass**"), size=3) + 
    theme_bw()

biomass_resist_fig_out_c3 <- ggplot( data=resist_dat_dry_fig_out, aes( x=G_C3_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt),se=F, col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-2.3,1.3)+
    labs(x="C3 Perennial Grasses\nInitial relative cover (%)", y="LRR - Standing Biomass\nless resistant                more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.35, y=1, label="Graze\nInitial biomass**"), size=3) + 
    theme_bw()

biomass_resist_fig_c3
biomass_resist_fig_c3_out

###  C4 perennial grasses -- Strong Grazing effect, Weak grazing x C4 interaction (p= 0.05104)
#' C4 grass figure
biomass_resist_fig <- ggplot( data=resist_dat_dry_fig, aes( x=G_C4_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
  geom_smooth( method="lm", aes(lty=graze_trt),se=F, col='black', size=0.6, alpha=0.15) + 
  geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
    ylim(-2.3,1.3)+
  labs(x="C4 Perennial Grasses\nInitial relative cover (%)", y="LRR - Standing Biomass\nless resistant                more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
  geom_text( aes( x=0.35, y=1, label="Graze**\nGrazing : C4 Grass`\nInitial biomass**"), size=3) + 
  theme_bw()

biomass_resist_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=G_C4_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt),se=F, col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
        ylim(-2.3,1.3)+
    labs(x="C4 Perennial Grasses\nInitial relative cover (%)", y="LRR - Standing Biomass\nless resistant                more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.35, y=1, label="Graze*\nC4 Grass**\nInitial biomass*"), size=3) + 
    theme_bw()

biomass_resist_fig
biomass_resist_fig_out


#' Show all figures with significant differences when including outliers
biomass_resist_8_out <- 
    biomass_resist_fig_f  + biomass_resist_fig_f_out +
    biomass_resist_fig_ann +biomass_resist_fig_ann_out +   
    biomass_resist_fig_c3  +  biomass_resist_fig_out_c3 +
    biomass_resist_fig +  biomass_resist_fig_out + 
    guide_area() + plot_layout(ncol=2, guides = 'collect', widths=c(1,1)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
biomass_resist_8_out


#' Show figures only for datasets with NO outliers removed
biomass_resist_4 <- biomass_resist_fig_f  + biomass_resist_fig_ann+  guide_area() +   
    biomass_resist_fig_c3  +   biomass_resist_fig  + plot_layout(ncol=3, guides = 'collect', widths=c(1,1,0.5)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
biomass_resist_4


tiff(filename="Resistance_biomass_all_FGs_Outliers.tiff", res=600, width=7, height = 15, units = "in")
biomass_resist_8_out
dev.off()

tiff(filename="Resistance_biomass_all_FGs_noOutliers.tiff", res=600, width=10, height = 7, units = "in")
biomass_resist_4
dev.off()





#' ** Supporting species figures**
#' 
#' Experimental Andger figure
biomass_ag_fig <- ggplot( data=resist_dat_dry_fig, aes( x=ag_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
  geom_smooth( method="lm", aes(lty=graze_trt), se=F,col='black', size=0.6, alpha=0.15) + 
  geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
  scale_linetype(guide="none") +
  ylim(-2.5,1.5)+
  labs(x="Andropogon gerardii\nInitial relative cover (%)", y="LRR - Standing Biomass", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
  geom_text( aes( x=0.125, y=1, label="Graze*\nInitial biomass**"), size=3) + 
  theme_bw()
biomass_ag_fig

biomass_ag_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=ag_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt),se=F, col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    scale_linetype(guide="none") +
    ylim(-2.5,1.5)+
    labs(x="Andropogon gerardii\nInitial relative cover (%)", y="LRR - Standing Biomass", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.125, y=1, label="Graze*\nAndropogon`\nInitial biomass**"), size=3) + 
    theme_bw()
biomass_ag_fig_out

#' Experimental Chogra figure
biomass_cg_fig <- ggplot( data=resist_dat_dry_fig, aes( x=cg_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
  geom_smooth( method="lm", aes(lty=graze_trt),se=F, col='black', size=0.6, alpha=0.15) + 
  geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
  scale_linetype(guide="none") +
  ylim(-2.5,1.5)+
  labs(x="Bouteloua gracilis\nInitial relative cover (%)", y="LRR - Standing Biomass", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
  geom_text( aes( x=0.125, y=1, label="Graze*\nInitial biomass*"), size=3) + 
  theme_bw()
biomass_cg_fig

biomass_cg_fig_out <- ggplot( data=resist_dat_dry_fig, aes( x=cg_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F,col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    scale_linetype(guide="none") +
    ylim(-2.5,1.5)+
    labs(x="Bouteloua gracilis\nInitial relative cover (%)", y="LRR - Standing Biomass", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.125, y=1, label="Graze*\nInitial biomass**"), size=3) + 
    theme_bw()
biomass_cg_fig_out

#' Experimental muhmon figure
biomass_mm_fig <- ggplot( data=resist_dat_dry_fig, aes( x=mm_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
  geom_smooth( method="lm", aes(lty=graze_trt),se=F, col='black', size=0.6, alpha=0.15) + 
  geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
  scale_linetype(guide="none") +
  ylim(-2.5,1.5)+
  labs(x="Muhlenbergia montana\nInitial relative cover (%)", y="LRR - Standing Biomass", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
  geom_text( aes( x=0.125, y=1, label="Graze**\n Graze : Muhmon* \nInitial biomass**"), size=3) + 
  theme_bw()
biomass_mm_fig

biomass_mm_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=mm_2018, y=biomass_converted_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt),se=F, col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=biomass_converted_2018), fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    scale_linetype(guide="none") +
    ylim(-2.5,1.5)+
    labs(x="Muhlenbergia montana\nInitial relative cover (%)", y="LRR - Standing Biomass", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.125, y=1, label="Graze*\nInitial biomass**"), size=3) + 
    theme_bw()
biomass_mm_fig_out


biomass_resist_c4_sp_outliers <- biomass_ag_fig + biomass_cg_fig + biomass_mm_fig + guide_area() +
    biomass_ag_fig_out + biomass_cg_fig_out + biomass_mm_fig_out + 
    plot_layout(ncol=4, guides = 'collect', widths=c(1,1,1,0.5)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
biomass_resist_c4_sp_outliers


biomass_resist_c4_sp <- biomass_ag_fig + biomass_cg_fig + biomass_mm_fig + guide_area() +
    plot_layout(ncol=2, guides = 'collect', widths=c(1,1)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
biomass_resist_c4_sp


tiff(filename="Resistance_biomass_c4_species_Outliers.tiff", res=600, width=12, height = 7, units = "in")
biomass_resist_c4_sp_outliers
dev.off()

tiff(filename="Resistance_biomass_c4_species.tiff", res=600, width=8, height = 8, units = "in")
biomass_resist_c4_sp
dev.off()


#' 
#'  #########
#'  
#'  Drought resistance:  **Shannon Diversity**  
#'  
#'  #########
#'

#' 
#' *Linear Models* 
#' 
#' *FORBS  -- Graze p=0.06,  Graze x Forb  p=0.08*
#' 
#' Removing outlier reduces effect 
#' 
resist_shan_f <- lm(shan_LRR ~ graze_trt * F_2018 + shan_2018_s, data = resist_dat_dry) 
resist_shan_f_out <- lm(shan_LRR ~ graze_trt * F_2018 + shan_2018_s, data = resist_dat_dry_out) 

summary(resist_shan_f)
summary(resist_shan_f_out)
anova(resist_shan_f)
anova(resist_shan_f_out)

plot(resist_shan_f, which=c(1,1),ask=F)
plot(resist_shan_f_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_shan_f))
qqnorm(residuals(resist_shan_f_out))

#' *G_Ann  --   Graze p=0.04, log_ann p=0.004*
#' 
#' Removing outlier doesn't influence results
#' 
resist_shan_ann <- lm(shan_LRR ~ graze_trt * log_Ann_Bi_2018 + shan_2018_s, data = resist_dat_dry) 
resist_shan_ann_out <- lm(shan_LRR ~ graze_trt * log_Ann_Bi_2018 + shan_2018_s, data = resist_dat_dry_out) 

summary(resist_shan_ann)
summary(resist_shan_ann_out)
anova(resist_shan_ann)
anova(resist_shan_ann_out)

plot(resist_shan_ann, which=c(1,1),ask=F)
plot(resist_shan_ann_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_shan_ann))
qqnorm(residuals(resist_shan_ann_out))

#' *G_C3  --  Graze p=0.09*
#' 
#' Removing outlier doesn't influence results
resist_shan_c3 <- lm(shan_LRR ~ graze_trt * G_C3_2018 + shan_2018_s, data = resist_dat_dry) 
resist_shan_c3_out <- lm(shan_LRR ~ graze_trt * G_C3_2018 + shan_2018_s, data = resist_dat_dry_out) 

summary(resist_shan_c3)
summary(resist_shan_c3_out)
anova(resist_shan_c3)
anova(resist_shan_c3_out)

plot(resist_shan_c3, which=c(1,1),ask=F)
plot(resist_shan_c3_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_shan_c3))
qqnorm(residuals(resist_shan_c3_out))


#' *G_C4  --    Graze p=0.038, Graze:C4 p=0.02*
#' 
#' Removing outlier reduces results
resist_shan_c4 <- lm(shan_LRR  ~ graze_trt * G_C4_2018 + shan_2018_s, data = resist_dat_dry) 
resist_shan_c4_out <- lm(shan_LRR  ~ graze_trt * G_C4_2018 + shan_2018_s, data = resist_dat_dry_out) 

summary(resist_shan_c4)
summary(resist_shan_c4_out)
anova(resist_shan_c4)
anova(resist_shan_c4_out)

plot(resist_shan_c4, which=c(1,1),ask=F)
plot(resist_shan_c4_out, which=c(1,1),ask=F)
qqnorm(residuals(resist_shan_c4))
qqnorm(residuals(resist_shan_c4_out))


#' *Andger  --    Graze p=0.009  ag p=0.07  graze:ag p<0.001*
resist_shan_ag <- lm(shan_LRR  ~ graze_trt * ag_2018 + shan_2018_s, data = resist_dat_dry) 
resist_shan_ag_out<- lm(shan_LRR  ~ graze_trt * ag_2018 + shan_2018_s, data = resist_dat_dry_out) 
summary(resist_shan_ag)
summary(resist_shan_ag_out)
anova(resist_shan_ag)
anova(resist_shan_ag_out)

#' *Chogra  --   Graze p=0.059*
resist_shan_cg <- lm(shan_LRR  ~ graze_trt * cg_2018 + shan_2018_s, data = resist_dat_dry) 
resist_shan_cg_out <- lm(shan_LRR  ~ graze_trt * cg_2018 + shan_2018_s, data = resist_dat_dry_out) 
summary(resist_shan_cg)
summary(resist_shan_cg_out)
anova(resist_shan_cg)
anova(resist_shan_cg_out)

#' *Muhmon --   NS*
resist_shan_mm <- lm(shan_LRR  ~ graze_trt * mm_2018 + shan_2018_s, data = resist_dat_dry) 
resist_shan_mm_out <- lm(shan_LRR  ~ graze_trt * mm_2018 + shan_2018_s, data = resist_dat_dry_out)
summary(resist_shan_mm)
summary(resist_shan_mm_out)
anova(resist_shan_mm)
anova(resist_shan_mm_out)



#' * Little test of annual-biennials *
#' 
#' We will not 
#'   Most common annual-biennial species
spp_info
spp_info %>% filter( fun_grp =="Ann_Bi") %>% arrange( sqrt_rel_abun_allyrs)
resist_dat_dry$ln_vo_2018 <- log(resist_dat_dry$vo_2018 +0.003)
resist_dat_dry$ln_bj_2018 <- log(resist_dat_dry$bj_2018 +0.001)

resist_shan_vo <- lm(shan_LRR  ~ graze_trt * ln_vo_2018 + shan_2018_s, data = resist_dat_dry) 
resist_shan_bj <- lm(shan_LRR  ~ graze_trt * ln_bj_2018 + shan_2018_s, data = resist_dat_dry) 
anova(resist_shan_vo)
anova(resist_shan_bj)
qqnorm(residuals(resist_shan_vo))
qqnorm(residuals(resist_shan_bj))

#' Experimental vulpia figure
div_vo_fig <- ggplot( data=resist_dat_dry, aes( x=ln_vo_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=shan_2018),fill='goldenrod3',alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-0.1,0.12)+
    labs(size="Initial\ndiversity",x="vuplia\nInitial relative cover (%)", y="LRR - Shannon diversity", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=-5, y=0.1, label="Graze* \nVulpia*\nGrazing : Vulpia`\nDiversity*"), size=3) +
    theme_bw()
div_vo_fig

#' Experimental vulpia figure
div_bj_fig <- ggplot( data=resist_dat_dry, aes( x=ln_bj_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=shan_2018),fill='goldenrod3',alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-0.1,0.12)+
    labs(size="Initial\ndiversity",x="Bromus\nInitial relative cover (%)", y="LRR - Shannon diversity", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=-5, y=0.1, label="Graze' \nBromus*"), size=3) +
    theme_bw()
div_bj_fig



#' ** FIGURES**
#' 
#'  

#'   *Forbs - Diversity*
#' 
div_forb_fig <- ggplot( data=resist_dat_dry_fig, aes( x=F_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( data=resist_dat_dry_fig, method="lm", se=F,aes(x=F_2018, y=shan_LRR, lty=graze_trt), col='black', size=.6, alpha=0.15) + 
    geom_point( aes(size=shan_2018),fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    labs(size="Initial\ndiversity",x="Forbs\nInitial relative cover (%) ", y="LRR - Shannon Diversity\nless resistant                      more resistant", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.5, y=0.065, label="Graze`\nGraze : Forb`"), size=3) + 
    theme_bw()

div_forb_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=F_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( data=resist_dat_dry_fig_out,se=F, method="lm", aes(x=F_2018, y=shan_LRR, lty=graze_trt), col='black', size=.6, alpha=0.15) + 
    geom_point(aes(size=shan_2018), fill='goldenrod3',alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    labs(size="Initial\ndiversity",x="Forbs\nInitial relative cover (%) ", y="LRR - Shannon Diversity\nless resistant                      more resistant", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.5, y=0.065, label="Graze*\nInitial diversity`"), size=3) + 
    theme_bw()

div_forb_fig
div_forb_fig_out

#'   *Annual-Biennial - Diversity*
#' 
#' GET RID OF 3 SEPARATE LINES
div_ann_fig <- ggplot( data=resist_dat_dry_fig, aes( x=log_Ann_Bi_2018, y=shan_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
   geom_smooth( data=resist_dat_dry_fig, method="lm",se=F, aes(x=log_Ann_Bi_2018, y=shan_LRR, lty=graze_trt), col='black', size=.6, alpha=0.15) + 
  geom_point(aes(size=shan_2018), fill='goldenrod3',  alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
   labs(size="Initial\ndiversity",x="Annual-Biennials\nInitial relative cover (%) ", y="LRR - Shannon Diversity\nless resistant                      more resistant", shape="Graze trt.", lty="Graze trt.") +
  geom_text( aes( x=-3, y=0.085, label="Graze*\nAnn-Biennial**"), size=3) + 
  theme_bw()

div_ann_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=log_Ann_Bi_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( data=resist_dat_dry_fig_out, method="lm", se=F,aes(x=log_Ann_Bi_2018, y=shan_LRR, lty=graze_trt), col='black', size=.6, alpha=0.15) + 
    geom_point( aes(size=shan_2018),fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    labs(size="Initial\ndiversity",x="Annual-Biennials\nInitial relative cover (%) ", y="LRR - Shannon Diversity\nless resistant                      more resistant", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=-3, y=0.085, label="Graze*\nAnn-Biennial*"), size=3) + 
    theme_bw()

div_ann_fig
div_ann_fig_out

#'  Figure
#'    *C3 grass - Diversity*
#' 
div_c3_fig <- ggplot( data=resist_dat_dry_fig, aes( x=G_C3_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( data=resist_dat_dry_fig, method="lm", se=F,aes(x=G_C3_2018, y=shan_LRR, lty=graze_trt), size=0.6, col='black',alpha=0.15) + 
    geom_point(aes(size=shan_2018), fill='goldenrod3',  alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    labs(size="Initial\ndiversity",x="C3 Perennial Grasses\nInitial relative cover (%)", y="LRR - Shannon Diversity\nless resistant                   more resistant", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.35, y=0.08, label="Graze`"), size=3) + 
    theme_bw()

div_c3_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=G_C3_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( data=resist_dat_dry_fig_out, method="lm", se=F,aes(x=G_C3_2018, y=shan_LRR, lty=graze_trt), size=0.6, col='black',alpha=0.15) + 
    geom_point( aes(size=shan_2018),fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    labs(size="Initial\ndiversity",x="C3 Perennial Grasses\nInitial relative cover (%)", y="LRR - Shannon Diversity\nless resistant                   more resistant", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.35, y=0.08, label="Graze*"), size=3) + 
    theme_bw()
div_c3_fig
div_c3_fig_out

#' FIGURE:   *C4 grass - Diversity*
#' 
#' 
div_c4_fig <- ggplot( data=resist_dat_dry_fig, aes( x=G_C4_2018, y=shan_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
  geom_smooth( data=resist_dat_dry_fig, method="lm", se=F, aes(x=G_C4_2018, y=shan_LRR, lty=graze_trt), size=0.6, col='black',alpha=0.15) + 
  geom_point( aes(size=shan_2018),fill='goldenrod3', alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
  labs(size="Initial\ndiversity",x="C4 Perennial Grasses\nInitial relative cover (%)", y="LRR - Shannon Diversity\nless resistant                   more resistant", shape="Graze trt.", lty="Graze trt.") +
  geom_text( aes( x=0.35, y=0.08, label="Graze*\nGraze : C4*"), size=3) + 
  theme_bw()

div_c4_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=G_C4_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( data=resist_dat_dry_fig_out, method="lm", se=F,aes(x=G_C4_2018, y=shan_LRR, lty=graze_trt), size=0.6, col='black',alpha=0.15) + 
    geom_point( aes(size=shan_2018),fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    labs(size="Initial\ndiversity",x="C4 Perennial Grasses\nInitial relative cover (%)", y="LRR - Shannon Diversity\nless resistant                   more resistant", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.35, y=0.08, label="Graze*\nInitial diversity`"), size=3) + 
    theme_bw()

div_c4_fig
div_c4_fig_out


#' Show all figures with significant differences when including outliers
div_resist_8_out <- div_forb_fig +div_forb_fig_out+ 
    div_ann_fig + div_ann_fig_out + 
    div_c3_fig + div_c3_fig_out + 
    div_c4_fig + div_c4_fig_out  +
    guide_area() +
    plot_layout(ncol=2, guides = 'collect', widths=c(1,1)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
div_resist_8_out

tiff(filename="Resistance_diversity_all_FGs_Outliers.tiff", res=600, width=7, height = 15, units = "in")
div_resist_8_out
dev.off()

div_resist_4 <- div_forb_fig + div_ann_fig + guide_area() + div_c3_fig + div_c4_fig +  
     plot_layout(ncol=3, guides = 'collect', widths=c(1,1,0.5)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
div_resist_4

tiff(filename="Resistance_diversity_all_FGs_noOutliers.tiff", res=600, width=10, height = 7, units = "in")
div_resist_4
dev.off()



#' **Figures for c4 grass species **
#' 

#' Experimental Andger figure
div_ag_fig <- ggplot( data=resist_dat_dry_fig, aes( x=ag_2018, y=shan_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
  geom_smooth( method="lm", aes(lty=graze_trt), se=F,col='black', size=0.6, alpha=0.15) + 
  geom_point(aes(size=shan_2018),fill='goldenrod3',alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
  ylim(-0.1,0.12)+
  labs(size="Initial\ndiversity",x="Andropogon gerardii\nInitial relative cover (%)", y="LRR - Shannon diversity", shape="Graze trt.", lty="Graze trt.") +
  geom_text( aes( x=0.1, y=0.1, label="Graze** \nAndropogon`\nGrazing : Andropogon***"), size=3) +
  theme_bw()
div_ag_fig

div_ag_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=ag_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F,col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=shan_2018),fill='goldenrod3',  alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-0.1,0.12)+
    labs(size="Initial\ndiversity",x="Andropogon gerardii\nInitial relative cover (%)", y="LRR - Shannon diversity", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.1, y=0.1, label="Graze** \nInitial diversity*\nGrazing : Andropogon**"), size=3) +
    theme_bw()
div_ag_fig_out

#' Experimental Chogra figure
div_cg_fig <- ggplot( data=resist_dat_dry_fig, aes( x=cg_2018, y=shan_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
  geom_smooth( method="lm", aes(lty=graze_trt),se=F, col='black', size=0.6, alpha=0.15) + 
  geom_point(aes(size=shan_2018),fill='goldenrod3',  alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
  ylim(-0.1,0.12)+
  labs(size="Initial\ndiversity",x="Bouteloua gracilis\nInitial relative cover (%)", y="LRR - Shannon diversity", shape="Graze trt.", lty="Graze trt.") +
  geom_text( aes( x=0.125, y=0.11, label="Graze`\nInitial diversity`"), size=3) + 
  theme_bw()
div_cg_fig

div_cg_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=cg_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F,col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=shan_2018),fill='goldenrod3',  alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-0.1,0.12)+
    labs(size="Initial\ndiversity",x="Bouteloua gracilis\nInitial relative cover (%)", y="LRR - Shannon diversity", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.125, y=0.11, label="Graze*\nInitial diversity`"), size=3) + 
    theme_bw()
div_cg_fig_out

#' Experimental Chogra figure
div_mm_fig <- ggplot( data=resist_dat_dry_fig, aes( x=mm_2018, y=shan_LRR, shape=graze_trt)) + 
  geom_hline(aes(yintercept=0), col='gray80', size=1)+
  geom_smooth( method="lm", aes(lty=graze_trt), se=F,col='black', size=0.6, alpha=0.15) + 
  geom_point(aes(size=shan_2018),fill='goldenrod3',  alpha=0.7) + 
  scale_shape_manual(values=c(21,24,22))+
  ylim(-0.1,0.12)+
  labs(size="Initial\ndiversity",x="Muhlenbergia montana\nInitial relative cover (%)", y="LRR - Shannon diversity", shape="Graze trt.", lty="Graze trt.") +
  geom_text( aes( x=0.125, y=0.11, label="Graze`"), size=3) + 
  theme_bw()
div_mm_fig

div_mm_fig_out <- ggplot( data=resist_dat_dry_fig_out, aes( x=mm_2018, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_smooth( method="lm", aes(lty=graze_trt), se=F,col='black', size=0.6, alpha=0.15) + 
    geom_point(aes(size=shan_2018),fill='goldenrod3', alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    ylim(-0.1,0.12)+
    labs(size="Initial\ndiversity",x="Muhlenbergia montana\nInitial relative cover (%)", y="LRR - Shannon diversity", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.125, y=0.11, label="Graze*"), size=3) + 
    theme_bw()
div_mm_fig_out


div_resist_c4_sp_outliers <- div_ag_fig + div_cg_fig + div_mm_fig + guide_area() +
    div_ag_fig_out + div_cg_fig_out + div_mm_fig_out + 
    plot_layout(ncol=4, guides = 'collect', widths=c(1,1,1,0.5)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
div_resist_c4_sp_outliers

div_resist_c4_sp <- div_ag_fig + div_cg_fig + div_mm_fig + guide_area() +
    plot_layout(ncol=2, guides = 'collect', widths=c(1,1)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
div_resist_c4_sp


tiff(filename="Resistance_diversity_c4_species_Outliers.tiff", res=600, width=12, height = 7, units = "in")
div_resist_c4_sp_outliers
dev.off()

tiff(filename="Resistance_diversity_c4_species.tiff", res=600, width=8, height = 8, units = "in")
div_resist_c4_sp
dev.off()









#' **Final figures for manuscript**
#' 
#' With updated modeled relationships
#' 

#' View dataframe -- should have updated grazing labels
str(resist_dat_dry_fig)
resist_dat_dry_fig$graze_trt <- factor(resist_dat_dry_fig_out$graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))


#' *Biomass ~ Diversity*
#' 
#'    No changes to this figure

biomass_shan_resist <- ggplot( data=resist_dat_dry_fig, aes(x=biomass_converted_LRR, y=shan_LRR, shape=graze_trt)) + 
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    geom_vline(aes(xintercept=0), col='gray80', size=1)+
    geom_point(fill='goldenrod3', cex=3,alpha=0.7) + 
    scale_shape_manual(values=c(21,24,22))+
    labs(y="LRR - Shannon Diversity\nless resistant                           more resistant", 
         x="LRR - End of Season Biomass\nless resistant                           more resistant", shape="Graze trt.") +
    theme_bw() + 
    theme(legend.position="none")
biomass_shan_resist



#'  *Fix variables*
resist_dat_dry_fig$biomass_converted_2018_s <-   as.vector(resist_dat_dry_fig$biomass_converted_2018_s)
resist_dat_dry_fig$shan_2018_s <-   as.vector(resist_dat_dry_fig$shan_2018_s)


#' *Biomass - C4 grass figure*
#' 
#' Redo model with updated grazing treatment labels
resist_bm_c4_fig <- lm(biomass_converted_LRR ~ graze_trt * G_C4_2018 + biomass_converted_2018_s, data = resist_dat_dry_fig) 

#' View model
resist_bm_c4_fig 
anova(resist_bm_c4_fig)  # matches above


#' Make predicted relationships based on C4 grass cover and grazing
y_hat <- ggpredict(resist_bm_c4_fig, terms = c("G_C4_2018","graze_trt"))
y_hat <-  y_hat %>% 
    rename(graze_trt = group)
y_hat


#' Figure
biomass_resist_fig <- ggplot(data = y_hat, aes(x = x, y = predicted, linetype=graze_trt)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
     # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=G_C4_2018, y=biomass_converted_LRR, shape=graze_trt, size=biomass_converted_2018),
               fill='goldenrod3', alpha=0.7) +
    scale_shape_manual(values=c(21,24,22))+
    labs(x="C4 Perennial Grasses\nInitial relative cover (%)", y="LRR - End of Season Biomass\nless resistant                           more resistant", shape="Graze trt.", lty="Graze trt.", size="Initial\nbiomass") +
    geom_text( aes( x=0.5, y=0.8, label="Graze**\nGrazing : C4 Grass`\nInitial biomass**"), size=3) + 
    theme_bw()+ 
    guides(shape = guide_legend(override.aes = list(size = 3)))
biomass_resist_fig



#' Faceted figure
resist_text <- data.frame(x=c(0.15,0.15,0,0),
                          y=c(1.4,-1.90,0,0),
                          graze_trt=c("Ungrazed","Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"),
                          label=c("More resistant","Less resistant","",""))
resist_text$graze_trt <- factor(resist_text $graze_trt, levels=c("Ungrazed","Growing-Season Grazing","Dormant-Season Grazing"))


biomass_resist_fig_facet <- ggplot(data = y_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1) +
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=G_C4_2018, y=biomass_converted_LRR, shape=graze_trt), cex=3, fill='goldenrod3', alpha=0.7) +
    scale_shape_manual(values=c(21,24,22))+
    ylim(-2.1,1.5) +
    labs(x="Initial C4 Perennial Grasses\n(rel. cover)", 
         y="Change in End of Season Biomass\nDecrease                               Increase", 
         shape="Graze trt.", 
         lty="Graze trt.") +
    geom_text(data=resist_text, aes( x=x, y=y, label=label), size=3, color='gray30') + 
    ggtitle('Graze**   C4 Grass (NS)   Graze:C4 Grass`   Initial biomass**')+
    theme_bw()+ 
    theme( plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9), legend.position = 'none') +
    facet_wrap(~graze_trt) +
    guides(shape = guide_legend(override.aes = list(size = 3)))
biomass_resist_fig_facet




#' *Diversity - C4 grass figure*
#' 

#' Redo model with updated grazing treatment labels
resist_shan_c4_mod <- lm(shan_LRR ~ graze_trt * G_C4_2018 + shan_2018_s, data = resist_dat_dry_fig) 

#' View model
resist_shan_c4_mod 
anova(resist_shan_c4_mod)  # matches above


#' Make predicted relationships based on C4 grass cover and grazing
z_hat <- ggpredict(resist_shan_c4_mod, terms = c("G_C4_2018","graze_trt"))
z_hat <-  z_hat %>% 
    rename(graze_trt = group)
z_hat


#' Figure
resist_div_c4_fig <- ggplot(data = z_hat, aes(x = x, y = predicted, linetype=graze_trt)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=G_C4_2018, y=shan_LRR, shape=graze_trt), fill='goldenrod3', alpha=0.7,cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    labs(x="C4 Perennial Grasses\nInitial relative cover (%)", y="LRR - Shannon Diversity\nless resistant                           more resistant", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=0.35, y=0.12, label="\nGraze*\nGraze : C4*"), size=3) + 
    theme_bw()+ 
    guides(shape = guide_legend(override.aes = list(size = 3)))
resist_div_c4_fig

title

#' Figure - Faceted
resist_div_c4_fig_facet <- ggplot(data = z_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=G_C4_2018, y=shan_LRR, shape=graze_trt), fill='goldenrod3', alpha=0.7,cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    labs(x="Initial C4 Perennial Grasses\n(rel. cover)", 
         y="Change in Shannon Diversity\nDecrease                               Increase", 
         shape="Graze trt.", 
         lty="Graze trt.") +
    #geom_text( aes( x=0.1, y=0.085, label="Graze*\nGraze : C4*"), size=3) + 
    ggtitle("Graze*   C4 Grass (NS)   Graze:C4*   Initial diversity (NS)") +
    theme_bw()+ 
    facet_wrap( ~graze_trt)+
    theme( plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9), legend.position = 'none') +
    guides(shape = guide_legend(override.aes = list(size = 3)))
resist_div_c4_fig_facet



#' *Diversity -- Annual-Biennials Figure*
#' 


#' Redo model with updated grazing treatment labels
resist_div_ann_mod <- lm(shan_LRR ~ graze_trt * log_Ann_Bi_2018 + shan_2018_s, data = resist_dat_dry_fig) 

#' View model
resist_div_ann_mod
anova(resist_div_ann_mod)  # matches above


#' Make predicted relationships based on annual-biennial cover and grazing
x_hat <- ggpredict(resist_div_ann_mod, terms = c("log_Ann_Bi_2018","graze_trt"))
x_hat <-  x_hat %>% 
    rename(graze_trt = group)
x_hat


#' Figure
div_ann_resist_fig <- ggplot(data = x_hat, aes(x = x, y = predicted, linetype=graze_trt)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=log_Ann_Bi_2018, y=shan_LRR, shape=graze_trt),
               fill='goldenrod3', alpha=0.7, cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    labs(x="Annual-Biennials\nInitial relative cover (%)", y="LRR - Shannon Diversity\nless resistant                           more resistant", shape="Graze trt.", lty="Graze trt.") +
    geom_text( aes( x=-2.6, y=0.06, label="Graze*\nAnn-Biennial**"), size=3) + 
    theme_bw()+ 
    guides(shape = guide_legend(override.aes = list(size = 3)))
div_ann_resist_fig


#' Figure - Facet
div_ann_resist_fig_facet <- ggplot(data = x_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=log_Ann_Bi_2018, y=shan_LRR, shape=graze_trt), fill='goldenrod3', alpha=0.7, cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    labs(x="Initial Annual-Biennials\n(ln (rel. cover))", 
         y="Change in Shannon Diversity\nDecrease                               Increase", 
         shape="Graze trt.", 
         lty="Graze trt.") +
    #geom_text( aes( x=-2.7, y=0.08, label="Graze*\nAnn-Biennial**"), size=3) + 
    facet_wrap(~graze_trt) +
    ggtitle("Graze*   Ann-Bi**   Graze:Ann-Bi (NS)   Initial diversity (NS)") +
    theme_bw()+ 
    theme(  plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9), legend.position = 'none') +
    guides(shape = guide_legend(override.aes = list(size = 3)))
div_ann_resist_fig_facet



#'  **Compile final figure**
#' 
#'

resist_fig <- biomass_resist_fig +  div_ann_resist_fig+ resist_div_c4_fig + guide_area() + 
    plot_layout(ncol=2, guides = 'collect') + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')') & 
    theme(plot.tag = element_text(size = 10))
resist_fig 

tiff(filename="drought_resist_fig_Jun24.tiff", res=600, width=10.5, height = 8.5, units = "in")
resist_fig
dev.off()



resist_fig_facet <- biomass_resist_fig_facet / plot_spacer() / resist_div_c4_fig_facet / plot_spacer() / div_ann_resist_fig_facet + 
    plot_layout(ncol=1, heights=c(1.7,0.1,1.7,0.1,1.7)) + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')') & 
    theme(plot.tag = element_text(size = 9, face='bold'))
resist_fig_facet

tiff(filename="drought_resist_fig_facet_Jul24.tiff", res=600, width=6.5, height = 10, units = "in")
resist_fig_facet
dev.off()




#' *GRASS SPECIES MODELS*


#' Biomass - Andger
#' 
#' 

#' Redo model with updated grazing treatment labels
resist_bm_ag_mod <- lm(biomass_converted_LRR ~ graze_trt * ag_2018 + as.vector(biomass_converted_2018_s), data = resist_dat_dry_fig) 

#' View model
resist_bm_ag_mod
anova(resist_bm_ag_mod)  # matches above

#' Make predicted relationships based on C4 grass cover and grazing
a_hat <- ggpredict(resist_bm_ag_mod, terms = c("ag_2018","graze_trt"))
a_hat <-  a_hat %>% 
    rename(graze_trt = group)
a_hat

#' Figure
bm_ag_resist_fig <- ggplot(data = a_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=ag_2018, y=biomass_converted_LRR, shape=graze_trt),
               fill='goldenrod3', alpha=0.7,cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    facet_wrap(~graze_trt) +
    xlim(-0.02,0.17)+
    labs(x="Andropogon gerardii\nInitial relative cover", 
         y="Change in End of Season Biomass\nDecrease                               Increase", 
         shape="Graze trt.") +
    ggtitle("Graze*   Andropogon (NS)   Graze:Andropogon (NS)   Initial biomass**") + 
    theme_bw()+ 
    theme(  plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9),
            axis.text.x=element_text(angle=90), legend.position = 'none') +
    guides(shape = guide_legend(override.aes = list(size = 3)))
bm_ag_resist_fig



#' Diversity - Andger
#' 
#' 

#' Redo model with updated grazing treatment labels
resist_shan_ag_mod <- lm(shan_LRR  ~ graze_trt * ag_2018 + as.vector(shan_2018_s), data = resist_dat_dry_fig) 

#' View model
resist_shan_ag_mod
anova(resist_shan_ag_mod)  # matches above

#' Make predicted relationships based on C4 grass cover and grazing
ag_hat <- ggpredict(resist_shan_ag_mod, terms = c("ag_2018","graze_trt"))
ag_hat <-  ag_hat %>% 
    rename(graze_trt = group)
ag_hat

#' Figure
shan_ag_resist_fig <- ggplot(data = ag_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=ag_2018, y=shan_LRR, shape=graze_trt),
               fill='goldenrod3', alpha=0.7, cex=3) +
    facet_wrap(~graze_trt) +
    xlim(-0.02,0.17)+
    scale_shape_manual(values=c(21,24,22))+
    labs(x="Andropogon gerardii\nInitial relative cover", 
         y="Change in Shannon Diversity\nDecrease                               Increase", 
         shape="Graze trt.") +
    ggtitle( "Graze**   Andropogon`   Graze:Andropogon***   Initial diversity (NS)") + 
    theme_bw()+ 
    theme(  plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9), legend.position='none',
            axis.text.x=element_text(angle=90)) +
    guides(shape = guide_legend(override.aes = list(size = 3)))
shan_ag_resist_fig




#' Biomass - Chogra
#' 
#' 
resist_bm_cg_mod <- lm(biomass_converted_LRR ~ graze_trt * cg_2018+ as.vector(biomass_converted_2018_s), data = resist_dat_dry_fig) 

#' View model
resist_bm_cg_mod
anova(resist_bm_cg_mod)  # matches above

#' Make predicted relationships based on C4 grass cover and grazing
c_hat <- ggpredict(resist_bm_cg_mod, terms = c("cg_2018","graze_trt"))
c_hat <-  c_hat %>% 
    rename(graze_trt = group)
c_hat

#' Figure
bm_cg_resist_fig <- ggplot(data = c_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=cg_2018, y=biomass_converted_LRR, shape=graze_trt),
               fill='goldenrod3', alpha=0.7,cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    facet_wrap(~graze_trt) +
    xlim(-0.02,0.22)+
    labs(x="Bouteloua gracilis\nInitial relative cover",
         y="Change in End of Season Biomass\nDecrease                               Increase", 
         shape="Graze trt.") +
    ggtitle("Graze*   Bouteloua (NS)   Graze:Bouteloua (NS)   Initial biomass*") + 
    theme_bw()+ 
     theme(  plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9),legend.position='none',
             axis.text.x=element_text(angle=90)) +
    guides(shape = guide_legend(override.aes = list(size = 3)))
bm_cg_resist_fig


#' Diversity - Chogra
#' 
resist_shan_cg_mod <- lm(shan_LRR  ~ graze_trt * cg_2018 + as.vector(shan_2018_s), data = resist_dat_dry_fig) 

#' View model
resist_shan_cg_mod
anova(resist_shan_cg_mod)  # matches above

#' Make predicted relationships based on C4 grass cover and grazing
cg_hat <- ggpredict(resist_shan_cg_mod, terms = c("cg_2018","graze_trt"))
cg_hat <-  cg_hat %>% 
    rename(graze_trt = group)
cg_hat

#' Figure
shan_cg_resist_fig <- ggplot(data = cg_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
   # geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=cg_2018, y=shan_LRR, shape=graze_trt),
               fill='goldenrod3', alpha=0.7, cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    facet_wrap(~graze_trt) +
    ylim(-.125,.075) +
    xlim(-0.02,0.22)+
    labs(x="Bouteloua gracilis\nInitial relative cover",
         y="Change in Shannon Diversity\nDecrease                               Increase", 
         shape="Graze trt.") +
    ggtitle( "Graze`   Bouteloua (NS)   Graze:Bouteloua (NS)   Initial diversity`") + 
    theme_bw()+ 
    theme(  plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9), legend.position='none',
            axis.text.x=element_text(angle=90)) +
    guides(shape = guide_legend(override.aes = list(size = 3)))
shan_cg_resist_fig



#' Biomass - Muhmon
#' 
#' 
resist_bm_mm_mod <- lm(biomass_converted_LRR ~ graze_trt * mm_2018+ as.vector(biomass_converted_2018_s), data = resist_dat_dry_fig) 

#' View model
resist_bm_mm_mod
anova(resist_bm_mm_mod)  # matches above

#' Make predicted relationships based on C4 grass cover and grazing
m_hat <- ggpredict(resist_bm_mm_mod, terms = c("mm_2018","graze_trt"))
m_hat <-  m_hat %>% 
    rename(graze_trt = group)
m_hat

#' Figure
bm_mm_resist_fig <- ggplot(data = m_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=mm_2018, y=biomass_converted_LRR, shape=graze_trt),
               fill='goldenrod3', alpha=0.7,cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    facet_wrap(~graze_trt) +
    xlim(-0.02,0.17)+
    labs(x="Muhlenbergia montana\nInitial relative cover", 
         y="Change in End of Season Biomass\nDecrease                               Increase", 
         shape="Graze trt.") +
    ggtitle("Graze**   Muhlenbergia (NS)   Graze:Muhlenbergia*   Initial biomass**") + 
    theme_bw()+ 
    theme(  plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9), legend.position = 'none',
            axis.text.x=element_text(angle=90)) +
        guides(shape = guide_legend(override.aes = list(size = 3)))
bm_mm_resist_fig



#' Diversity - muhmon
#' 
#' 
resist_shan_mm_mod <- lm(shan_LRR  ~ graze_trt * mm_2018 + as.vector(shan_2018_s), data = resist_dat_dry_fig) 

#' View model
resist_shan_mm_mod
anova(resist_shan_mm_mod)  # matches above

#' Make predicted relationships based on C4 grass cover and grazing
mm_hat <- ggpredict(resist_shan_mm_mod, terms = c("mm_2018","graze_trt"))
mm_hat <-  mm_hat %>% 
    rename(graze_trt = group)
mm_hat

#' Figure
shan_mm_resist_fig <- ggplot(data = mm_hat, aes(x = x, y = predicted)) +
    # plot a line at y=0
    geom_hline(aes(yintercept=0), col='gray80', size=1)+
    # plot the fitted line
    geom_line() +
    # plot the confidence intervals
    #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill='gray70') +
    # plot the raw data
    geom_point(data = resist_dat_dry_fig, aes(x=mm_2018, y=shan_LRR, shape=graze_trt),
               fill='goldenrod3', alpha=0.7, cex=3) +
    scale_shape_manual(values=c(21,24,22))+
    facet_wrap(~graze_trt) +
    xlim(-0.02,0.17)+
    labs(x="Muhlenbergia montana\nInitial relative cover", 
         y="Change in Shannon Diversity\nDecrease                               Increase",
         shape="Graze trt.") +
    ggtitle("Graze`   Muhlenbergia (NS)   Graze:Muhlenbergia (NS)   Initial diversity (NS)") + 
    theme_bw()+ 
    theme(  plot.title = element_text(size = 8, hjust=0.5), text=element_text(size=9), legend.position='none',
            axis.text.x=element_text(angle=90)) +
        guides(shape = guide_legend(override.aes = list(size = 3)))
shan_mm_resist_fig





#'  **Compile final figures**
#' 
#'

#' Final C4 species figure for supplement
spp_fig <- bm_ag_resist_fig + shan_ag_resist_fig + bm_cg_resist_fig + shan_cg_resist_fig + bm_mm_resist_fig + shan_mm_resist_fig + 
    plot_layout(ncol=2, guides = 'collect') + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
spp_fig 

#' Diversity only figure
spp_div <- shan_ag_resist_fig + shan_cg_resist_fig + shan_mm_resist_fig + 
    plot_layout(ncol=1, guides = 'collect') + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
spp_div 

#' Biomass only figure
spp_bm <- bm_ag_resist_fig  + bm_cg_resist_fig+ bm_mm_resist_fig  +
    plot_layout(ncol=1, guides = 'collect') + 
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '.', tag_suffix = ')')
spp_bm


tiff(filename="drought_resist_fig_SPECIES_Jul24.tiff", res=600, width=11, height = 10, units = "in")
spp_fig 
dev.off()

tiff(filename="drought_resist_fig_SPECIES_June24_biomass.tiff", res=600, width=6.5, height = 8, units = "in")
spp_bm
dev.off()

tiff(filename="drought_resist_fig_SPECIES_June24_div.tiff", res=600, width=6.5, height = 8, units = "in")
spp_div
dev.off()
