
## Plant and Microbial Data Analyses
getwd()
setwd("my_dir")

# requires
require(dplyr); require(ggplot2); require(cowplot); require(lme4)
require(rstan)


# csv files 
allPlant <- read.csv("PulseChase_VegDataALL.csv") # does not have biomass values 
allPlantMass <- read.csv("PlantData_MassALL.csv") # has biomass values 

# miscellaneous functions 

# converts from delta into atom % 
atom_Cent <- function(delta_val) {
  AR <- 0.0111803
  atom_cent <- (100*AR*(delta_val/1000 + 1))/(1 + AR*(delta_val/1000 + 1))
  return(atom_cent/100)
}

# standardization function
std_func <- function(x){
  x_std <- rep(0,length(x));
  for(i in 1:length(x_std)){
    x_std[i] <- (x[i] - mean(x,na.rm=T))/sd(x,na.rm=T)
  }
  return(x_std);
}

#### Working with Plant Data 

# preparing to merge 
allPlant$Plot <- as.factor(allPlant$Plot)
allPlantMass$Plot <- as.factor(allPlantMass$Plot)
levels(allPlantMass$Pool) <- c("Ag Veg", "Ag Veg", "Rh", "Root")

# merging 
allPlant_all <- merge(allPlant, allPlantMass, by.x = c("Harvest","Treatment", "Pool", "Sub", "Plot"),
                      by.y = c("Harvest","Treatment", "Pool", "Sub", "Plot"))

# setting up a biomass C variable 
allPlant2 <- select(allPlant_all, c(1:9,11))
colnames(allPlant2)[10] <- "Biomass"
allPlant3 <- mutate(allPlant2, massC = Biomass*centC*0.01)

# scaling measured biomass to whole plot values
# belowground biomass: 
# Conversion of 4 cores to whole plot volume 
#1/(4*pi*(2.5^2)/(100*100))= 127.34
# aboveground biomass: 
# conversion of AG harvest strips to whole plot volume 
#100/7.5 = 13.33

# new dataframe 
ap4 <- mutate(allPlant3, massCcorr = ifelse(Pool == "Ag Veg", massC*13.33, massC*127.34),
              excess13C = 1000*massCcorr*(atom_Cent(d13C) - atom_Cent(-13.5)))


ap4$Harvests <- factor(ap4$Harvest, levels(ap4$Harvest)[c(1,2,4,3)])
ap4$HarvestC <- as.numeric(ap4$Harvests)


#####Multiple line graphs

shoot_root <- filter(ap4, Pool != "Rh")
shoot_root_df <- shoot_root %>%
  group_by(Treatment, Harvests,Pool) %>%
  summarise(CmassM = mean(massCcorr,na.rm = T), CmassSE = sd(massCcorr,na.rm = T)/sqrt(4),
            d13CM = mean(d13C,na.rm = T), d13CSE = sd(d13C,na.rm = T)/sqrt(4),
            excess13CM = mean(excess13C,na.rm = T), excess13CSE = sd(excess13C,na.rm = T)/sqrt(4))

# creating variables with preferred orderings 
shoot_root_df$Trt <- factor(shoot_root_df$Treatment, levels = c("Grazed", "Exc"))
levels(shoot_root_df$Trt) <- c("Grazed", "Exclosure")
shoot_root_df$day <- c(0,2,2,7,7,32,32,0,2,2,7,7,32,32)
levels(shoot_root_df$Pool)[1] <- "Shoots"

massC_shoots <- ggplot(data = shoot_root_df[which(shoot_root_df$Pool == "Shoots"),]) + 
  geom_point(aes(x = day, y = CmassM, shape = Trt), size = 3) +
  geom_line(aes(x = day, y = CmassM, group = Trt)) +
  geom_errorbar(aes(x = day, ymin = CmassM - CmassSE,
                         ymax = CmassM + CmassSE), width = 0.5) +
  xlab("day") + ylab("Biomass C (gm2)") +
  theme_bw() + 
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 

massC_roots <- ggplot(data = shoot_root_df[which(shoot_root_df$Pool == "Root"),]) + 
  geom_point(aes(x = day, y = CmassM, shape = Trt), size = 3) +
  geom_line(aes(x = day, y = CmassM, group = Trt)) +
  geom_errorbar(aes(x = day, ymin = CmassM - CmassSE,
                    ymax = CmassM + CmassSE),width = 0.5) +
  xlab("day") + ylab("") +
  theme_bw() + 
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 


d13C_shoots <- ggplot(data = shoot_root_df[which(shoot_root_df$Pool == "Shoots"),]) +
  geom_point(aes(x = day, y = d13CM, shape = Trt), size = 3) +
  geom_line(aes(x = day, y = d13CM, group = Trt)) +
  geom_errorbar(aes(x = day, ymin = d13CM - d13CSE,
                ymax = d13CM + d13CSE), width = 0.5) +
    xlab("day") + ylab("delta 13C") +
  theme_bw()+ 
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 

d13C_roots <- ggplot(data = shoot_root_df[which(shoot_root_df$Pool == "Root"),]) +
  geom_point(aes(x = day, y = d13CM, shape = Trt), size = 3) +
  geom_line(aes(x = day, y = d13CM, group = Trt)) +
  geom_errorbar(aes(x = day, ymin = d13CM - d13CSE,
                    ymax = d13CM + d13CSE), width = 0.5) +
  xlab("day") + ylab("") +
  theme_bw() + 
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 
  
excess13C_roots <- ggplot(data = shoot_root_df[which(shoot_root_df$Pool == "Root"),]) +
  geom_point(aes(x = day, y = excess13CM, shape = Trt), size = 3) +
  geom_line(aes(x = day, y = excess13CM, group = Trt)) +
  geom_errorbar(aes(x = day, ymin = excess13CM - excess13CSE,
                    ymax = excess13CM + excess13CSE), width = 0.5) +
  xlab("day") + ylab("") +
  theme_bw() + 
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 


excess13C_shoots <- ggplot(data = shoot_root_df[which(shoot_root_df$Pool == "Shoots"),]) +
  geom_point(aes(x = day, y = excess13CM, shape = Trt), size = 3) +
  geom_line(aes(x = day, y = excess13CM, group = Trt)) +
  geom_errorbar(aes(x = day, ymin = excess13CM - excess13CSE,
                    ymax = excess13CM + excess13CSE),width = 0.5) +
  xlab("day") + ylab("excess 13C (mg)") +
  theme_bw() + 
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 



########### Working with Microbial Data 

setwd("my_dir")

iso_dat <- read.csv("PulseChase_MicrobialData.csv")
str(iso_dat)
# creating a factor out of the harvest integers 
iso_dat2 <- within(iso_dat, {
  Harvest <- as.factor(Harvest)
})

# renaming harvest levels
levels(iso_dat2$Fumigated) <- c("no", "no", "yes")
levels(iso_dat2$Harvest) <- c("32day", "2day", "7day")

# reordering levels
iso_dat2$Harvest <- factor(iso_dat2$Harvest,
                    levels(iso_dat2$Harvest)[c(2,3,1)])

# renaming columns 
colnames(iso_dat2)[7:8] <- c("mgC_gSoil", "d13C")

# separating values obtained from fumigated and unfumigated samples 
iso_dat3b <- filter(iso_dat2, Fumigated == "yes")
iso_dat3c <- filter(iso_dat2, Fumigated == "no")

# calculating Microbial biomass (mg/g dry soil) as difference 
# btw fumigated and unfumigated carbon concentrations 
iso_dat3b$MB <- iso_dat3b$mgC_gSoil - iso_dat3c$mgC_gSoil

# adding the extractable C and its d13C to the iso_dat3b dataframe 
iso_dat3b$e_C <- iso_dat3c$mgC_gSoil
iso_dat3b$d13C_e <- iso_dat3c$d13C


# best estimate of background Microbial d13C is -13.5 
# leading to an atom % of 0.01090905. This value will be propagated through
# all of the microbial data from pulse-chase plots to calculate atom % excess
# and 'excess 13C' recovered in microbial biomass pool 

# first step is to deconvolute the d13C signature of microbial biomass pool 
# following mass-balance approach of Dijkstra et al. (2006). I then use the atom_Cent
# function (found below) to convert this d13C into atom %. We then calculate 
# excess label recovered as the product of microbial biomass and the atom % excess, 
# multiplied by a factor of 40 to convert to estimated value in milligrams per plot 
# (i.e. 1 m^2)  

iso_dat3d <- mutate(iso_dat3b, MB_d13C = (mgC_gSoil*d13C - e_C*d13C_e)/MB, MB_plot = MB*0.8*5*100*100, 
                    excess13C_mg = MB_plot*(atom_Cent(MB_d13C) - atom_Cent(-13.42)))



# grouping by harvest and treatment 
col_dfMIC13 <- iso_dat3d %>% 
  group_by(Trt, Harvest) %>%
  summarise(meanMIC = mean(MB,na.rm = T),
            seMIC = sd(MB,na.rm = T)/sqrt(4),
            mean_MB_d13C = mean(MB_d13C, na.rm = T),
            se_MB_d13C = sd(MB_d13C, na.rm = T)/sqrt(4), 
            mean_excess13C = mean(excess13C_mg, na.rm = T),
            se_excess13C = sd(excess13C_mg, na.rm = T)/sqrt(4))


col_dfMIC13$t <- c(2,7,32,2,7,32)
col_dfMIC13$Trt <- factor(c(rep("Exclosure", 3), rep("Grazed", 3)),
                          levels = c("Grazed", "Exclosure"))

# adding in baseline values
trt2 <- factor(c("Exclosure", "Grazed"), levels = c("Grazed", "Exclosure"))
Harvest2 <- factor(c("0day", "0day"), levels = c("0day", "2day", "7day", "32day"))
t2 <- c(0,0)
meanMIC2 <- c(0,0)
seMIC2 <- c(0,0)
mean_MB_d13C2 <- c(-13.42, -13.42)
se_MB_d13C2 <- c(0.2, 0.2)
mean_excess13C2 <- c(0,0)
se_excess13C2 <- c(0,0)

mod_df <- data.frame(Trt = trt2, Harvest = Harvest2, meanMIC = meanMIC2,
                     seMIC = seMIC2, mean_MB_d13C = mean_MB_d13C2, se_MB_d13C = se_MB_d13C2,
                     mean_excess13C = mean_excess13C2, se_excess13C = se_excess13C2, t = t2)

identical(colnames(mod_df),colnames(col_dfMIC13)) # check



col_dfMIC13C2 <- rbind.data.frame(col_dfMIC13, mod_df) 

col_dfMIC13C2$Trt <- factor(c("Exclosure", "Exclosure", "Exclosure", "Grazed", "Grazed", "Grazed", "Exclosure", "Grazed"), 
                            levels = c("Grazed", "Exclosure"))


                            
MBd13C_twoLines <- ggplot() +
  geom_point(data = col_dfMIC13C2, aes(x = t, y = mean_MB_d13C, shape = Trt), size = 3) +
  geom_line(data = col_dfMIC13C2, aes(x = t, y = mean_MB_d13C, group = Trt)) +
  geom_errorbar(data = col_dfMIC13C2, aes(x = t, ymin = mean_MB_d13C - se_MB_d13C,
                                        ymax = mean_MB_d13C + se_MB_d13C),
                width = 0.5) + 
 xlab("day") + ylab("") +
  theme_bw() +
  scale_y_continuous(limits = c(-15, 0)) +
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 


MB_excess13C <- ggplot() +
  geom_point(data = col_dfMIC13, aes(x = t, y = mean_excess13C, shape = Trt), size = 3) +
  geom_line(data = col_dfMIC13, aes(x = t, y = mean_excess13C, group = Trt)) +
  geom_errorbar(data = col_dfMIC13, aes(x = t, ymin = mean_excess13C - se_excess13C,
                                          ymax = mean_excess13C + se_excess13C),
                width = 0.5) + 
  xlab("day") + ylab("") + 
  theme_bw() + 
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 



MB_twoLines <- ggplot() +
  geom_point(data = col_dfMIC13, aes(x = t, y = meanMIC, shape = Trt), size = 3) +
  geom_line(data = col_dfMIC13, aes(x = t, y = meanMIC, group = Trt)) +
  geom_errorbar(data = col_dfMIC13, aes(x = t, ymin = meanMIC - seMIC,
                                        ymax = meanMIC + seMIC),
                width = 0.5) + 
  xlab("day") + ylab("") + 
  theme_bw() +
  scale_x_continuous(limits = c(0,33),breaks = c(0,2,7,32)) 





#################### Figure 1 
## Plot grid with all the responses including microbial 
plot_grid3 <- plot_grid(massC_shoots, massC_roots, MB_twoLines,
                        d13C_shoots, d13C_roots, MBd13C_twoLines,
                         excess13C_shoots, excess13C_roots,
                         MB_excess13C, nrow = 3)


ggsave("Fig1.pdf",plot_grid3, width = 12, height = 6, dpi = 450)





################## Figure 2, roots and microbes at 48 hours 
twoday_mic <- filter(iso_dat3d, Harvest == "2day")
ap4_root <- filter(ap4, Pool == "Root")

## Computing exudation
twoday_mic %>% group_by(Trt) %>% summarise(exudation = mean(excess13C_mg)*(1/0.5)*(1/0.6))

## Matching column names, etc. 
levels(ap4_root$Harvest) <- c("0day","2day","7day","32day")
colnames(twoday_mic)[2] <- "Treatment"
levels(ap4_root$Treatment) <- c("Exc", "G")

twoday_root <- filter(ap4_root, Harvest == "2day")
twoday_rootMic <- merge(twoday_root, twoday_mic, by = c("Treatment","Harvest", "Sub", "Plot"))

# reordering levels
twoday_rootMic$Treatment <- factor(twoday_rootMic$Treatment,
                           levels(twoday_rootMic$Treatment)[c(2,1)])

######### lme4 analyses
library(MuMIn) # for obtaining R^2 values 
## Exudation as function of root biomass 
rootMass_exudate <- lmer(excess13C_mg ~ massCcorr + Treatment + (1|Plot:Treatment),twoday_rootMic)
summary(rootMass_exudate) # 
Betas1 <- fixef(rootMass_exudate)
# Obtaining marginal and conditional R^2 values 
r.squaredGLMM(rootMass_exudate)
# R^2 = 0.78 (marginal and conditional same since random effect variance
# has been estimated as zero)

## Exudation as function of root C allocation 
rootAllocate_exudate <- lmer(excess13C_mg ~ excess13C + Treatment + (1|Plot:Treatment),twoday_rootMic)
summary(rootAllocate_exudate) # Root enrichment appears most significant 
rootAllocate_exudate2 <- lmer(excess13C_mg ~ excess13C*Treatment + (1|Plot:Treatment),twoday_rootMic)
summary(rootAllocate_exudate2) # root enrichment estimate is robust 

# Checking root enrichment by itself for purposes of plotting 
rootAllocate_exudatePLOT <- lmer(excess13C_mg ~ excess13C + (1|Plot:Treatment),twoday_rootMic)
summary(rootAllocate_exudatePLOT)
# Obtaining marginal and conditional R^2 
r.squaredGLMM(rootAllocate_exudatePLOT)
# Marginal (fixed effect only): 0.62 
# Conditional (including random effects): 0.74 

Betas2 <- fixef(rootAllocate_exudatePLOT)

##### Plotting Figure 2
# Root Biomass as a predictor of label recovery in microbial pool 
rootMassplot <- ggplot(twoday_rootMic, aes(x = massCcorr, y = excess13C_mg, shape = Treatment)) + 
  geom_point() + theme_bw() + geom_abline(aes(intercept = Betas1[1] + Betas1[3], slope = Betas1[2])) +
  ylab("13C-C label in Microbes (mg/m2)") + xlab("Root Biomass C") +
  ggtitle("13C-C label at 48 hours scales with root biomass")

# Root allocation of C as a predictor of label recovery 
rootAllocateplot <- ggplot(twoday_rootMic, aes(x = excess13C, y = excess13C_mg, shape = Treatment)) +
  geom_point() + theme_bw() + geom_abline(aes(intercept = Betas2[1], slope = Betas2[2])) +
  ylab("13C-C label in Microbes (mg/m2)") + xlab("13C-C label in Roots (mg/m2)") +
  ggtitle("Linking root C allocation to microbial C at 48 hours")

plotby2 <- plot_grid(rootMassplot, rootAllocateplot, labels = c("A", "B"),
                     rel_heights = 0.4)

# saving with better dimensions 
ggsave("Fig2.pdf",plotby2, width = 8, height = 3, dpi = 450)



#### Working with Raw Lignin Data 
setwd("my_dir")
lignin <- read.csv("chw_output.csv") # contains lignin data 
EA_IRMS <- read.csv("Soils_EAIRMS.csv")
colnames(EA_IRMS)[6:7] <- c("centN", "centC")
EA_IRMS$Pasture[24:26] <- "Skuhn"

lignin_soil <- filter(lignin, Pool != "plant")


lignin_soil$Sample.ID <- lignin_soil$sample_id
lignin_soil$Pasture = lignin_soil$Site
lignin_soil$Trt <- lignin_soil$treatment


# merging with EA_IRMS data frame below 
levels(EA_IRMS$Sample.ID) <- c("E1-P", "E2-P", "E3-P", "E4-P", "G1-P", "G2-P", "G3-P",
                               "G4-P", "Griffith_E1-B", "Griffith_E2-B", "Griffith_E3-B",
                               "Griffith_G1-B", "Griffith_G2-B", "Griffith_G3-B",
                               "S.Kuhn_E1-B", "S.Kuhn_E2-B", "S.Kuhn_E3-B",
                               "S.Kuhn_G1-B", "S.Kuhn_G2-B", "S.Kuhn_G3-B",
                               "S4_E1-B", "S4_E2-B", "S4_E3-B",
                               "S4_G1-B", "S4_G2-B", "S4_G3-B")

allSoil_df <- merge(lignin_soil, EA_IRMS, by = "Sample.ID")


### Quick check on predictive value of relationship btw lignin and SOC 
vsc.sed <- lmer(centC ~ vsc.sed*Back + (1|Pasture), allSoil_df)
summary(vsc.sed)
r.squaredGLMM(vsc.sed) # 0.9
coefs_vsc <- fixef(vsc.sed)
################### Figure 3, SOC/lignin relationship 

# Plotting scatterplot of relationship 
vsc_sed_centC <- ggplot(allSoil_df, aes(x = vsc.sed, y = 10*centC, shape = Back,
                                        group = treatment)) + 
  geom_point() + theme_bw() + 
  xlab("Lignin Phenols (mg/gSoil)") + ylab("SOC (mg/gSoil)")

ggsave("Fig3.pdf",vsc_sed_centC, height = 3, width = 4,dpi=450)


#################### Analyses of LMM in Stan and lme4, example work-flow 


library(rstan) # loading rstan 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd(" ")

# concatenate data into list, being careful with specification of baseline values
# e.g. 
ap4_abov <- filter(ap4, Pool == "Ag Veg", Harvest != 0)
y <- ap4_abov$excess13C[-11] 

# creating a plot variable that has separate labels for grazed and excluded
plot2 <- factor(interaction(ap4_abov$Treatment, ap4_abov$Plot),
                levels = c("Exc.1","Exc.2","Exc.3","Exc.4","Grazed.1","Grazed.2","Grazed.3","Grazed.4"))

# Note that for aboveground data, we must remove NA on 11th row 
stan_datList <- list(N = as.integer(38), K = as.integer(8), M = as.integer(3),
                     y = y, G = -(as.integer(ap4_abov$Treatment[-11])-2),
                     plot = as.integer(plot2)[-11], day = as.integer(ap4_abov$Harvests[-11])-1,
                     sigma_scales = sd(y)/2)


# Root data 
ap4_root <- filter(ap4, Pool == "Root")
y <- ap4_root$excess13C

# creating a plot variable that has separate labels for grazed and excluded
plot2 <- factor(interaction(ap4_root$Treatment, ap4_root$Plot),
                levels = c("Exc.1","Exc.2","Exc.3","Exc.4","Grazed.1","Grazed.2","Grazed.3","Grazed.4"))

stan_datList <- list(N = as.integer(40), K = as.integer(8), M = as.integer(3),
                     y = y, G = -(as.integer(ap4_root$Treatment)-2),
                     plot = as.integer(plot2), day = as.integer(ap4_root$Harvest)-1,
                     sigma_scales = sd(y)/2)

 
 ##### Microbes 

 y <- iso_dat3d$excess13C_mg

 plot2 <- factor(interaction(iso_dat3d$Trt, iso_dat3d$Plot),
                 levels = c("Exc.1","Exc.2","Exc.3","Exc.4","G.1","G.2","G.3","G.4"))
 
 stan_datList <- list(N = as.integer(40), K = as.integer(8), M = as.integer(3),
                      y = y, G = -(as.integer(iso_dat3d$Trt)-2),
                      plot = as.integer(plot2), day = as.integer(ap4_root$Harvest)-1,
                      sigma_scales = sd(y)/2)

# First run models with with few chains and iterations to compile and 
# debug 
                        
stan_plantMic <- stan(file = "BIRexclosures_Bayes.stan", data = stan_datList,
                        iter = 500, chains = 3);
                        
# Now re-run with pre-compiled models from above, but switching contents of 
# data list for different responses
stan_plantMicFIT <- stan(fit = stan_plantMic, data = stan_datList,
                        iter = 5000, chains = 12);

                        
# check model summaries
print(stan_plantMicFIT) # check particularly for r-hat < 1.1 (closer to 1.0 is better) 

# check model diagnostics
library(shinystan)
launch_shinystan(stan_soilFIT) # web API for visual inspection of chains, etc. 

# extracting coefficients 
AG_massCcorr <- extract(stan_plantMicFIT, par = c("B[2]","sigma_plot[1]", "sigma_plot[2]",
                                  "mean_G","mean_EXC")) # extracts the grazing coefficient 

AG_d13C <- extract(stan_plantMicFIT, par = c("B[2]","sigma_plot[1]", "sigma_plot[2]",
                                             "mean_G","mean_EXC"))


Root_massCcorr <- extract(stan_plantMicFIT, par = c("B[2]","sigma_plot[1]", "sigma_plot[2]",
                                                    "mean_G","mean_EXC"))

Root_d13C <- extract(stan_plantMicFIT, par = c("B[2]","sigma_plot[1]", "sigma_plot[2]",
                                               "mean_G","mean_EXC"))

# concatenate and save all into list 

# Comparing shoot response results to inference from "maximum likelihood" 
library(lme4)
ap4_abov$PLOT <- plot2
maxLike_shootExcess13C <- lmer(excess13C ~ Treatment + (1|Harvests) + (1|PLOT),
                               ap4_abov)
summary(maxLike_shootExcess13C)

maxLike_shootMassC <- lmer(massCcorr ~ Treatment + (1|Harvests) + (1|PLOT),
                               ap4_abov)
summary(maxLike_shootMassC)

maxLike_shootd13C <- lmer(d13C ~ Treatment + (1|Harvests) + (1|PLOT),
                               ap4_abov)
summary(maxLike_shootd13C)


plot2 <- factor(interaction(ap4_root$Treatment, ap4_root$Plot),
                levels = c("Exc.1","Exc.2","Exc.3","Exc.4","Grazed.1","Grazed.2","Grazed.3","Grazed.4"))

ap4_root$PLOT <- plot2

maxLike_rootExcess13C <- lmer(excess13C ~ Treatment + (1|Harvests) + (1|PLOT),
                              ap4_root)
summary(maxLike_rootExcess13C)
### As can be seen, these results are extremely similar. The Bayes models have a bit more
### uncertainty, and thus are a bit more conservative. 

######## Working with EA-IRMS Data

Back_df <- filter(EA_IRMS, Back == "yes")
Pulse_df <- filter(EA_IRMS, Back == "no")

str(EA_IRMS)
## Setting up data for stan model
y <- EA_IRMS$d15N # Toggle response here
G <- -(as.numeric(EA_IRMS$Trt)-2) # grazed is baseline 
Back <- as.numeric(EA_IRMS$Back)-1 # Pulse-chase plots are baseline 

EA_IRMS_df <- list(N = as.integer(26), K = as.integer(3), 
                   y = y, G = G, Back = Back, 
                   Past = as.integer(EA_IRMS$Pasture),
                   sigma_scales = sd(y)/2, sigma_location = mean(y))


EA_IRMS_m1 <- stan(file = "EA_IRMS_mod.stan", data = EA_IRMS_df, 
                   chains = 2, iter = 200)

# Run models for all responses above 
EAIRMS2 <- stan(fit = EA_IRMS_m1, data = EA_IRMS_df, chains = 12, iter = 5000);


# Analyzing treatment effects on lignin parameters 
y <- allSoil_df$adal.v # Toggle response here
G <- -(as.numeric(allSoil_df$Trt.y)-2) # grazed is baseline 
Back <- as.numeric(allSoil_df$Back)-1 # Pulse-chase plots are baseline 

EA_IRMS_df <- list(N = as.integer(26), K = as.integer(3), 
                   y = y, G = G, Back = Back, 
                   Past = as.integer(allSoil_df$Pasture),
                   sigma_scales = sd(y)/2, sigma_location = mean(y))
EAIRMS2 <- stan(fit = EA_IRMS_m1, data = EA_IRMS_df, chains = 12, iter = 5000);





 