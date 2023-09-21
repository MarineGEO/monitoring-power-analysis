library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)

# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/Coral Disease")


# Import the survey data
cover <- read.csv("coral_cover_pwr.csv")

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
library(broom)

group1 <- covtype3 %>% subset(Location == "Belize")
group1$cover <- as.numeric(group1$cover)
group2 <- covtype3 %>% subset(Location == "Cayman")
group2$cover <- as.numeric(group2$cover)

cohen.d(group1$cover, group2$cover)

#what is our power with current sample size of images per group?
pwr.t.test(n=120,d=0,sig.level=.05,alternative="greater") #power = 0.55

#what sample size do we need if we want a power of .80
pwr.t.test(power = 0.8,d=0.01,sig.level=.05,alternative="greater") #n = 124000!!! (for stony coral)
#n =  234 for octocoral

df <- data.frame()

for(var in 10:10000) {
  ptt <- pwr.t.test(n=var,d=0.01,sig.level=.05,alternative="greater")
  ptt_df <- as.data.frame(tidy(ptt))
  df <- df %>%
    bind_rows(ptt_df)
}

location_pwr <- ggplot(df, aes(x = n, y = power)) +
  geom_line(color = "black") +
  scale_y_continuous("Power") +
  scale_x_continuous("Number of Images per Location") +
  ggtitle("T-test: comparison of locations") +
  theme(plot.title = element_text(size = 13,hjust = 0, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),  
        axis.title = element_text(colour = "black", size = 10),        
        legend.position = "none")
location_pwr


###percent stony coral cover, difference by year in Cayman (t-test) for sites with repeated surveys
cay <- covtype %>% subset(Location == "Cayman") %>% subset(Label_General == "Stony Coral") %>%
  subset(SiteName != "Meadows" & SiteName != "Pizza Hut" & SiteName != "Henrys Beach" &
           SiteName != "Tibbetts Top" & SiteName != "Spl ash House")

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
pwr.t.test(n=147,d=0.37,sig.level=.05,alternative="greater") #power = 0.95
pwr.t.test(power = 0.8 ,d=0.38,sig.level=.05,alternative="greater") #power = 0.95

df <- data.frame()

for(var in 10:500) {
  ptt <- pwr.t.test(n=var,d=0.37,sig.level=.05,alternative="greater")
  ptt_df <- as.data.frame(tidy(ptt))
  df <- df %>%
    bind_rows(ptt_df)
}

year_pwr <- ggplot(df, aes(x = n, y = power)) +
  geom_line(color = "black") +
  scale_y_continuous("Power") +
  scale_x_continuous("Number of Images per Year") +
  ggtitle("T-test: consecutive years w/in location") +
  theme(plot.title = element_text(size = 13,hjust = 0, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),  
        axis.title = element_text(colour = "black", size = 10),        
        legend.position = "none")
year_pwr

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

pwr.anova.test(f=0.08,k=7,n=21,sig.level=0.05) #power is only 0.1!

pwr.anova.test(f=0.08,k=7,power=0.8,sig.level=0.05) #opimally would have 241 photos/transect to get differences among sites

df <- data.frame()

for(var in 4:500) {
  pat <- pwr.anova.test(f=0.08,k=7,n=var,sig.level=0.05)
  pat_df <- as.data.frame(tidy(pat))
  df <- df %>%
    bind_rows(pat_df)
}

transect_pwr <- ggplot(df, aes(x = n, y = power)) +
  geom_line(color = "black") +
  scale_y_continuous("Power") +
  scale_x_continuous("Number of Images per Transect") +
  ggtitle("Anova: transects within a location") +
  theme(plot.title = element_text(size = 13,hjust = 0, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),  
        axis.title = element_text(colour = "black", size = 10),        
        legend.position = "none")
transect_pwr


