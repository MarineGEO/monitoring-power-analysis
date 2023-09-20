library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)

# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/Coral Disease")


# Import the survey data
cover1 <- read.csv("annotations_CBC_pwr.csv")
cover2 <- read.csv("annotations_CCMI_pwr.csv")

cover2$SiteName <- recode(cover2$SiteName, "Snap Shot" = "Snapshot",
                          "Sailfin Reef" = "Sailfin")

cover <- rbind(cover1,cover2)
  
split <- strsplit(cover$Date, "/")
split2 <- matrix(unlist(split), ncol=3, byrow=TRUE)
split2 <- as.data.frame(split2)
cover <- cbind(cover, split2)

cover <- cover %>% rename("Year" = "V3") %>%
  select(Location, Name, Date, SiteName, Habitat, Depth, Year, Annotator, Label, Row, Column)

cover$Label <- as.factor(cover$Label)
levels(cover$Label)

cover$Label <- recode(cover$Label, "CNAt" = "CNAT",
                       "SMIC" = "SINT")


cover <- cover %>% mutate('Label_General' = 'x')
cover <- within(cover, Label_General[Label == "THAL"] <- "Thalassia")
cover <- within(cover, Label_General[Label == "AAGA"|Label == "APAL"|Label == "ATEN"|
                                       Label == "Coral"|Label == "MALC"|Label == "MCAV"|
                                       Label == "MCOM"|Label == "OANN"|Label == "OFAV"|
                                       Label == "PAST"|Label == "PPOR"|Label == "PSTRI"|
                                       Label == "SINT"|Label == "SSID"|Label == "PCLI"|
                                       Label == "CNAT"|Label == "ACER"|Label == "AAGA"|
                                       Label == "ATEN"|Label == "DLAB"|Label == "MCOM"|
                                       Label == "SRAD"|Label == "MLAM"|Label == "AGAR"|
                                       Label == "HCUC"|Label == "OFRA"|Label == "MDEC"|
                                       Label == "EFAS"|Label == "DSTO"] <- "Stony Coral")
cover <- within(cover, Label_General[Label == "ANTI"|Label == "GORGO"|Label == "GVEN"|
                                       Label == "MURO"|Label == "PTER"|Label == "EUNI"|
                                       Label == "BRIA"|Label == "PSEUGORG"|Label == "MURI"]<- "Octocoral")
cover <- within(cover, Label_General[Label == "Hal_spp"|Label == "Macro"|Label == "Turbin"|
                                       Label == "Lob_spp"|Label == "Dicsp"|Label == "BRAN-CALC"] <- "Macroalgae")
cover <- within(cover, Label_General[Label == "Sand"|Label == "LSUB_RUB"] <- "Unconsolidated (Sand/Rubble)")
cover <- within(cover, Label_General[Label == "Tk Tf"|Label == "Tn Tf"] <- "Turf Algae")
cover <- within(cover, Label_General[Label == "Sand"|Label == "LSUB_RUB"] <- "Unconsolidated (Sand/Rubble)")
cover <- within(cover, Label_General[Label == "SpgOth"|Label == "SPONGR"|Label == "SPTU"|
                                       Label == "ENSP"|Label == "CLIONI"|Label == "SPONGB"|
                                       Label == "SPONGV"] <- "Sponge")
cover <- within(cover, Label_General[Label == "PALY"|Label == "Other Inv"|Label == "A"] <- "Other Encrusting Invert")

cover <- within(cover, Label_General[Label == "CCA 1"|Label == "LIME"] <- "Hard Substrate")
cover <- within(cover, Label_General[Label == "TAPE"|Label == "Unk"|Label == "SHAD"] <- "Unidentified")

cover <- within(cover, Label_General[Label == "CYAN"|Label == "Cyan red"] <- "Cyanobacteria")

cover$SiteName <- recode(cover$SiteName, "CBC Lagoo" = "CBC Lagoon")

write.csv(cover, "coral_cover_pwr.csv") 
###############################
#power analysis: number of images per location & transect
###############################

cover2 <- cover %>% 
  group_by(Location, Habitat, SiteName, Year, Name, Label_General) %>%
  summarize(n = length(Label_General)) %>%
  mutate(cover = (n/40)*100) %>%
  ungroup() %>%
  complete(Label_General, nesting(Location, Habitat, SiteName, Year,  Name), 
           fill = list(n = 0, cover = 0)) %>%
  arrange(Location, Name, Habitat, SiteName, Year, Label_General) %>%
  select(Location, Name, Habitat, SiteName, Year, Label_General, n, cover)

#Look at stony coral cover first (can plug other cover types here)
covtype <- cover2 %>% subset(Label_General == "Stony Coral")


###how much replication per location (e.g. Belize, Little Cayman, to detect differences in cover by location)
  
  #reduce sample to first available year from each project to get rid of repeated sites
covtype1 <- covtype %>% subset(Location == "Belize" & Year == "2019"|Location == "Cayman" & Year == "2022"|
  Location == "Belize" & Year == "2022" & SiteName == "Tobacco Reef") %>%
  subset(!(Location == "Cayman" & SiteName == "CCMI Lagoon")) #this site is a low-cover outlier
 
 #create balanced groups by location
covtype2 <- covtype1 %>%                              
  group_by(Location, SiteName, Year) %>%
  mutate(img_number = row_number())

max <- covtype2 %>% group_by(Location, SiteName, Year) %>% summarize(max = max(img_number))
#each group has at least 20 photos 

covtype3 <- covtype2 %>% subset(img_number < 21)
  
t.test(cover ~ Location, data = covtype3)

  #calculate effect size (Cohen's D)
library(effsize)
library(effectsize)
library(pwr)

group1 <- covtype3 %>% subset(Location == "Belize")
group1$cover <- as.numeric(group1$cover)
group2 <- covtype3 %>% subset(Location == "Cayman")
group2$cover <- as.numeric(group2$cover)

cohen.d(group1$cover, group2$cover)

  #what is our power with current sample size of images per group?
pwr.t.test(n=120,d=0.23,sig.level=.05,alternative="greater") #power = 0.55

  #what sample size do we need if we want a power of .80
pwr.t.test(power = 0.8,d=0.23,sig.level=.05,alternative="greater") #n = 124000!!! (for stony coral)
                                                                  #n =  234 for octocoral

###percent stony coral cover, difference by year in Cayman (t-test) for sites with repeated surveys
cay <- covtype %>% subset(Location == "Cayman") %>% subset(Label_General == "Stony Coral") %>%
  subset(SiteName != "Meadows" & SiteName != "Pizza Hut" & SiteName != "Henrys Beach" &
         SiteName != "Tibbetts Top" & SiteName != "Splash House")

  #balance replication between years
cay1 <- cay %>%                              
  group_by(SiteName, Year) %>%
  mutate(img_number = row_number())

max <- cay1 %>% group_by(SiteName, Year) %>% summarize(max = max(img_number))

cay <- cay1 %>% subset(img_number < 22)

group1 <- cay %>% subset(Year == "2022")
group1$cover <- as.numeric(group1$cover)
group2 <- cay %>% subset(Year == "2023")
group2$cover <- as.numeric(group2$cover)

cohen.d(group1$cover, group2$cover)

  #what is our power with current sample size of images per group?
pwr.t.test(n=147,d=0.38,sig.level=.05,alternative="greater") #power = 0.95
pwr.t.test(power = 0.8 ,d=0.38,sig.level=.05,alternative="greater") #power = 0.95


  #Aiming for ~100 photos/year may be good enough for detecting change in stony coral cover over time
    #but 2023's images were scored by coralnet robot which may have underestimated cover

t.test(cover ~ Year, data = cay)


###How many photos per site/transect to detect variability among sites/transects?
  
  #create balanced groups

cay2 <- cay %>%                              
  group_by(SiteName, Year) %>%
  mutate(img_number = row_number())

max <- cay2 %>% group_by(SiteName, Year) %>% summarize(max = max(img_number))
  #each group has at least 21 photos

cay2 <- cay2 %>% subset(img_number < 22)

  #calculate effect size
mod1 <- aov(cover ~ SiteName, data = cay2)
summary(mod1)
eta_squared(mod1, partial = FALSE) #effect size is 0.09

pwr.anova.test(f=0.09,k=7,n=21,sig.level=0.05) #power is only 0.1!

pwr.anova.test(f=0.09,k=7,power=0.8,sig.level=0.05) #opimally would have 241 photos/transect to get differences among sites



