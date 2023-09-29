library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)
library(vegan)
library(simr)

# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/Coral Disease")


# Import the survey data
cover <- read.csv("coral_cover_pwr.csv")
# cover <- read.csv("data/raw_coral_cover_pwr.csv")

DateSplit <- cover$Date
split <- strsplit(DateSplit, "/")
split2 <- as.data.frame(matrix(unlist(split), ncol=3, byrow=TRUE))
split3 <- cbind(DateSplit, split2)
split3 <- split3 %>% distinct()
cover <- merge(cover, split3, by.x=c("Date"),by.y = c("DateSplit"),all.x=TRUE)

cover <- cover %>% 
  unite("MonthYear", c("V1", "V3"), sep = "_") %>%
  select(-c(V2))

levels(as.factor(cover$MonthYear))

cover <- cover %>% mutate(time_point = "x") %>%
  mutate(time_point = case_when(
    MonthYear == "1_20" ~ "one",
    TRUE ~ as.character(time_point))) %>% 
  mutate(time_point = case_when(
    MonthYear == "10_19" ~ "one",
    TRUE ~ as.character(time_point))) %>% 
  mutate(time_point = case_when(
    MonthYear == "12_22" ~ "three",
    TRUE ~ as.character(time_point))) %>%
  mutate(time_point = case_when(
    MonthYear == "5_22" ~ "two",
    TRUE ~ as.character(time_point))) %>%
  mutate(time_point = case_when(
    MonthYear == "6_23" ~ "two",
    TRUE ~ as.character(time_point))) %>%
  mutate(time_point = case_when(
    MonthYear == "7_22" ~ "one",
    TRUE ~ as.character(time_point)))


###############################
#power analysis: number of images per location & transect
###############################

cover2 <- cover %>% 
  group_by(Location, time_point, Date, Habitat, SiteName, Year, Name, Label_General) %>%
  summarize(n = length(Label_General)) %>%
  mutate(cover = (n/40)*100) %>%
  ungroup() %>%
  complete(Label_General, nesting(Location, Date, Habitat, SiteName, Year,  Name), 
           fill = list(n = 0, cover = 0)) %>%
  arrange(Location, time_point, Date, Name, Habitat, SiteName, Year, Label_General) %>%
  select(Location, time_point, Date, Name, Habitat, SiteName, Year, Label_General, n, cover)


#Look at stony coral cover first (can plug other cover types here)
covtype <- cover2 %>% subset(Label_General == "Stony Coral")


###how much replication per location (e.g. Belize, Little Cayman, to detect differences in cover by location)

#reduce sample to first available year from each project to get rid of repeated sites
covtype1 <- covtype %>% subset(Location == "Belize" & time_point == "one"|Location == "Cayman" & time_point == "one")
 # subset(!(Location == "Cayman" & SiteName == "CCMI Lagoon")) #this site is a low-cover outlier

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
pwr.t.test(n=133,d=0.62,sig.level=.05,alternative="greater") #power = 0.55

#what sample size do we need if we want a power of .80
pwr.t.test(power = 0.8,d=0.62,sig.level=.05,alternative="greater") #n = 30,000!!! (for stony coral)
#n =  234 for octocoral

df <- data.frame()

for(var in 2:200) {
  ptt <- pwr.t.test(n=var,d=0.62,sig.level=.05,alternative="greater")
  ptt_df <- as.data.frame(tidy(ptt))
  df <- df %>%
    bind_rows(ptt_df)
}


location_pwr <- ggplot(df, aes(x = n, y = power)) +
  geom_line(color = "black") +
  scale_y_continuous("Power") +
  scale_x_continuous("Number of Images per Location", sec.axis = sec_axis(~ . /20, name = "Number of 20-Image Transects")) +
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
pwr.t.test(n=147,d=0.37,sig.level=.05,alternative="greater") #power = 0.95
pwr.t.test(power = 0.8,d=0.38,sig.level=.05,alternative="greater") #power = 0.95

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
  scale_x_continuous("Number of Images per Year",sec.axis = sec_axis(~ . /20, name = "Number of 20-Image Transects")) +
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

###percent stony coral cover, difference by year in Belize (glm) for sites with repeated surveys
bel <- covtype %>% subset(Location == "Belize") %>% subset(Label_General == "Stony Coral") %>%
  subset(SiteName != "Tobacco Reef")

#balance replication between years
bel1 <- bel %>%                              
  group_by(SiteName, time_point) %>%
  mutate(img_number = row_number())

#try glm instead:
bel2 <- bel1 %>% mutate(propcov = cover/100) %>%
  mutate(countcov = (cover*40/100))


mod2 <- glm(countcov ~ time_point, data=bel2, family = "poisson")
summary(mod2)
hist(resid(mod2))
rsquared(mod2)
Anova(mod2)

pwr.f2.test(sig.level = 0.05, u = 1, f2 = 0.19/(1 - 0.19), power = 0.8) #n = 33 + 1 + 1 = 35
p <- pwr.f2.test(u=1, v=100, f2=0.19/(1-0.19), sig.level=.05)

model1 <- glmer(countcov ~ time_point + (1|SiteName), family="poisson",
                data=bel2)
summary(model1)
fixef(model1)["time_point"]


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

#repeat by transects: glm instead of anova
#calculate effect size

library(sensemakr)

mod2 <- glm(cover ~ SiteName, data = cay2)
summary(mod2)
Anova(mod2)
partial_r2(mod2)
eta_squared(mod1, partial = FALSE) #effect size is 0.09


############################
#Species Accumulation Curves
############################
cover3 <- cover %>%
  group_by(SiteName) %>%
  mutate(PicNum = as.numeric(as.factor(Name)))

sum <- cover3 %>% group_by(SiteName) %>%
  summarize(max = max(PicNum))

#reducing to 21 images per site ensures that we aren't including repeat surveys,
#which confound the accumulation curves (unlikely to find new species in quadrats you surveyed last year)
cover4 <- cover3 %>% subset(PicNum < 22)

sto <- within(cover4, Label[Label_General != "Stony Coral"] <- "Other")

levels(as.factor(sto$Label))
sto$Label <- recode(sto$Label, "PSTRI" = "PSTR")

sto <- sto %>% mutate("Count" = 1)

stocast <- sto %>%
  pivot_wider(names_from = Label, values_from = Count)  %>% select(-c(X, Date, Habitat, Depth, Annotator, Row, Column, Label_General, PicNum)) %>%
  replace(is.na(.), 0)


stomat = stocast[,7:ncol(stocast)]
stomat <- as.matrix(stomat)

accurvepoints<-specaccum(stomat, method="random", permutations=100)

PointAccum <- ggplot() +
  geom_line(aes(x = accurvepoints$sites, y = accurvepoints$richness), color = "black", size = 1) +
  labs(y = "Stony Coral Species Richness") +
  scale_x_continuous("Total Number of Points", sec.axis = sec_axis(~ . /25, name = "Number of 25-point Photos")) +
  theme(plot.title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 12),
        legend.key = element_rect(fill = "white"))
PointAccum  

tiff("PointAccum",width = 7, height = 5, units = "in", res = 600)
PointAccum
dev.off()

#What do the accumulation curves look like for each photo added if we do 10, 20, 30, 40 points per photo?
num_simulations <- 100

pt10_list <- lapply(1:num_simulations, function(i){
  
  pt10 <- cover4 %>%
    group_by(Name) %>%
    slice_sample(n = 10, replace = FALSE) %>%
    mutate(Count = 1)
  
  pt10 <- within(pt10, Label[Label_General != "Stony Coral"] <- "Other")
  
  cast10 <- pt10 %>% group_by(Location, Year, SiteName, Name, Label) %>%
    summarize(Count = sum(Count)) %>%
    pivot_wider(names_from = Label, values_from = Count) %>%
    replace(is.na(.), 0)
  
  mat10 = cast10[,5:ncol(cast10)]
  mat10 <- as.matrix(mat10)
  
  specaccum(mat10, method="random", permutations=500)$richness
  
})

accurvepics10_richness <- Reduce("+",pt10_list)/length(pt10_list)

#20 points

pt20_list <- lapply(1:num_simulations, function(i){
  
  pt20 <- cover4 %>%
    group_by(Name) %>%
    slice_sample(n = 20, replace = FALSE) %>%
    mutate(Count = 1)
  
  pt20 <- within(pt20, Label[Label_General != "Stony Coral"] <- "Other")
  
  cast20 <- pt20 %>% group_by(Location, Year, SiteName, Name, Label) %>%
    summarize(Count = sum(Count)) %>%
    pivot_wider(names_from = Label, values_from = Count) %>%
    replace(is.na(.), 0)
  
  mat20 = cast20[,5:ncol(cast20)]
  mat20 <- as.matrix(mat20)
  
  specaccum(mat20, method="random", permutations=500)$richness
  
})

accurvepics20_richness <- Reduce("+",pt20_list)/length(pt20_list)

#25 points
pt25_list <- lapply(1:num_simulations, function(i){
  
  pt25 <- cover4 %>%
    group_by(Name) %>%
    slice_sample(n = 25, replace = FALSE) %>%
    mutate(Count = 1)
  
  pt25 <- within(pt25, Label[Label_General != "Stony Coral"] <- "Other")
  
  cast25 <- pt25 %>% group_by(Location, Year, SiteName, Name, Label) %>%
    summarize(Count = sum(Count)) %>%
    pivot_wider(names_from = Label, values_from = Count) %>%
    replace(is.na(.), 0)
  
  mat25 = cast25[,5:ncol(cast25)]
  mat25 <- as.matrix(mat25)
  
  specaccum(mat25, method="random", permutations=500)$richness
})

accurvepics25_richness <- Reduce("+",pt25_list)/length(pt25_list)
#30 Points

pt30_list <- lapply(1:num_simulations, function(i){
  
  pt30 <- cover4 %>%
    group_by(Name) %>%
    slice_sample(n = 30, replace = FALSE) %>%
    mutate(Count = 1)
  
  pt30 <- within(pt30, Label[Label_General != "Stony Coral"] <- "Other")
  
  cast30 <- pt30 %>% group_by(Location, Year, SiteName, Name, Label) %>%
    summarize(Count = sum(Count)) %>%
    pivot_wider(names_from = Label, values_from = Count) %>%
    replace(is.na(.), 0)
  
  mat30 = cast30[,5:ncol(cast30)]
  mat30 <- as.matrix(mat30)
  
  specaccum(mat30, method="random", permutations=500)$richness
  
})

accurvepics30_richness <- Reduce("+",pt30_list)/length(pt30_list)

#35 Points
pt35_list <- lapply(1:num_simulations, function(i){
  
  pt35 <- cover4 %>%
    group_by(Name) %>%
    slice_sample(n = 35, replace = FALSE) %>%
    mutate(Count = 1)
  
  pt35 <- within(pt35, Label[Label_General != "Stony Coral"] <- "Other")
  
  cast35 <- pt35 %>% group_by(Location, Year, SiteName, Name, Label) %>%
    summarize(Count = sum(Count)) %>%
    pivot_wider(names_from = Label, values_from = Count) %>%
    replace(is.na(.), 0)
  
  mat35 = cast35[,5:ncol(cast35)]
  mat35 <- as.matrix(mat35)
  
  specaccum(mat35, method="random", permutations=500)$richness
})

accurvepics35_richness <- Reduce("+",pt35_list)/length(pt35_list)

#40 Points
pt40_list <- lapply(1:num_simulations, function(i){
  
  pt40 <- cover4 %>%
    group_by(Name) %>%
    slice_sample(n = 40, replace = FALSE) %>%
    mutate(Count = 1)
  
  pt40 <- within(pt40, Label[Label_General != "Stony Coral"] <- "Other")
  
  cast40 <- pt40 %>% group_by(Location, Year, SiteName, Name, Label) %>%
    summarize(Count = sum(Count)) %>%
    pivot_wider(names_from = Label, values_from = Count) %>%
    replace(is.na(.), 0)
  
  mat40 = cast40[,5:ncol(cast40)]
  mat40 <- as.matrix(mat40)
  
  specaccum(mat40, method="random", permutations=500)$richness
  
})

accurvepics40_richness <- Reduce("+",pt40_list)/length(pt40_list)

AccumPlot <- ggplot() +
  geom_line(aes(x = 1:503, y = accurvepics10_richness, color = "10"), size = 1) +
  geom_line(aes(x = 1:503, y = accurvepics20_richness, color = "20"), size = 1) +
  geom_line(aes(x = 1:503, y = accurvepics25_richness, color = "25"), size = 1) +
  geom_line(aes(x = 1:503, y = accurvepics30_richness, color = "30"), size = 1) + 
  geom_line(aes(x = 1:503, y = accurvepics35_richness, color = "35"), size = 1) + 
  geom_line(aes(x = 1:503, y = accurvepics40_richness, color = "40"), size = 1) +
  scale_x_continuous(breaks = seq(0,1000, by = 100)) +
  #scale_y_continuous(limits = c(0,100)) +
  labs(x = "Number of Photos",
       y = "Stony Coral Species Richness",
       color = "Points per Photo") +
  scale_color_manual(values = c("10" = "black", "20" = "purple4", "25" = "magenta3",
                                "30" = "deeppink", "35" = "salmon1","40" = "khaki")) +
  theme(plot.title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 12),
        legend.key = element_rect(fill = "white"))

tiff("AccumPicsByPointNum.tif",width = 7, height = 5, units = "in", res = 300)
AccumPlot
dev.off()


#######################
#Simulate cover measurement for different numbers of photos and points
#######################
cover5 <- cover4 %>% mutate()

num_simulations <- 100

# Silence dplyr warnings
options(dplyr.summarise.inform = FALSE)

df <- data.frame()

for(sim in 1:num_simulations){
  for(var in 5:40) {
    
    sample <- cover4 %>% mutate(Count = 1) %>%
      group_by(SiteName, Name) %>% 
      slice_sample(n = var, replace = FALSE) %>%
      group_by(SiteName, Name, Label_General) %>%
      summarize(cover = (sum(Count))/var*100) %>%
      ungroup() %>%
      complete(Label_General, nesting(SiteName, Name), 
               fill = list(cover = 0)) %>%
      group_by(Label_General) %>%
      summarize(n = var, mean = mean(cover), se = sd(cover)/sqrt(var)) %>% #se(cover)) %>%
      subset(Label_General == "Stony Coral" | Label_General == "Octocoral"|
               Label_General == "Macroalgae") %>%
      mutate(simulation = sim)
    
    df <- df %>%
      bind_rows(sample)
  }
}

df_output <- df %>%
  group_by(Label_General, n) %>%
  summarize(mean = mean(mean), se = mean(se))

MeanCovPlot <- ggplot(df_output, aes(x = n, y = mean)) +
  geom_line(aes(color = Label_General), size = 1) +
  geom_point(aes(fill = Label_General), pch = 21, size = 2.5, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se, color = Label_General), width = 0.2, alpha = 0.4,
                position = position_dodge(width = 0.9)) + 
  #scale_x_continuous(breaks = seq(0,1000, by = 100)) +
  #scale_y_continuous(limits = c(0,100)) +
  labs(x = "Number of Points",
       y = "Mean Cover",
       color = "Label Type") +
  scale_color_manual(values = c( "Octocoral" = "purple4", "Macroalgae" = "seagreen",
                                 "Stony Coral" = "salmon1")) +
  scale_fill_manual(values = c( "Octocoral" = "purple4", "Macroalgae" = "seagreen",
                                 "Stony Coral" = "salmon1")) +
  theme(plot.title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 12),
        legend.key = element_rect(fill = "white"))

tiff("CoverMeansByPointNum.tif",width = 7, height = 5, units = "in", res = 300)
MeanCovPlot
dev.off()

  
