library(tidyverse)
# Simulation to support power analysis

# 1. Calculate Mean and SD for real, observed transects ####

## Data issues
# Inconsistant site names
# Missing site names

df_raw <- read_csv(file = "data/raw_coral_cover_pwr.csv")

# Verify 40 points per image
verify_points <- df_raw %>%
  count(Name)

df <- df_raw %>%
  mutate(SiteName = case_when(
    SiteName == "CBC Lagoo" ~ "CBC Lagoon",
    T ~ SiteName
  )) %>%
  mutate(transect = paste0(SiteName, Year))

verify_transects <- df %>%
  count(transect, Name) %>%
  count(transect)

## 16 to 42 photos per transect

# Create numeric id for transects
transect_ids <- df %>%
  count(transect) %>%
  rownames_to_column("id") %>%
  mutate(id = as.numeric(id))

df_final <- left_join(df, transect_ids)

# Calculate photo-level mean Stony Coral cover
sites <- df_final %>%
  count(id, SiteName) %>%
  select(-n)

stony_coral_cover <- df_final %>%
  filter(Label_General == "Stony Coral") %>%
  group_by(id, Name) %>%
  summarize(coral_cover = (n()/40) * 100) %>%
  left_join(sites)

## Positively skewed, long tail to the right

ggplot(stony_coral_cover, aes(coral_cover)) + geom_histogram()

ggplot(stony_coral_cover, aes(coral_cover)) + 
  geom_histogram() +
  facet_wrap(~SiteName)

## Slightly more normal when plotted by Site

params <- stony_coral_cover %>%
  group_by(id) %>%
  summarize(mean = mean(coral_cover),
            sd = sd(coral_cover)) %>%
  left_join(sites)

ggplot(params, aes(mean, sd, color = SiteName)) + 
  geom_point()

## Positive correlation between mean and sd
# Suggests very patchy reef

# 2. Set mean and sd levels for simulation ####

ggplot(params, aes(mean, sd, color = SiteName)) + 
  geom_point()

minimum_mean <- min(params$mean)
maximum_mean <- max(params$mean)

minimum_sd <- min(params$sd)
maximum_sd <- max(params$sd)

# Observed mean / sd ratios
# sd = ratio * mean
params$sd_mean_ratio <- params$mean / params$sd

min_ratio <- min(params$sd_mean_ratio)
max_ratio <- max(params$sd_mean_ratio)

ggplot(params, aes(sd_mean_ratio)) + geom_histogram()

sim_mean <- c(5, 10, 20, 30)
# sim_ratio_sd_to_mean <- c(.5, 1.25, 2)
sim_ratio_sd_to_mean <- c(.2, .4, .6)

# 3. Run simulation ####

# Number of images per transect: 
num_replicates <- 21
# Number of points per grid:
num_points <- 10000
# Number of simulations:
# (Number of cover values drawn from a normal distribution)
num_simulations <- 500
# Maximum number of points to sample per frame 
max_sampled_points <- 40

results <- setNames(
  lapply(sim_mean, function(mean){
    
    setNames(
      lapply(sim_ratio_sd_to_mean, function(ratio){
        
        # Multiply ratio times mean to get standard deviation
        sd <- ratio * mean
        
        # Generate number of cover values equal to simulation rounds
        # Rounded to 2 to create integer number of points per grid
        # Having to use abs() needs to be resolved
        simulated_cover_values <- abs(
          round(rnorm(num_simulations, mean, sd), 2)
        )
        
        setNames(
          lapply(simulated_cover_values, function(simulation_number){
            
            num_points_with_cover <- simulation_number * 100
            
            # setNames(
            lapply(1:max_sampled_points, function(num_points_to_sample){
              
              # setNames(
              lapply(1:num_replicates, function(replicate){
                
                # Get random number of points from grid
                sampled_points <- sample(1:num_points, num_points_to_sample)
                
                # Determine which have coral
                values_with_cover <- which(sampled_points <= num_points_with_cover)
                
                if(length(values_with_cover) == 0){
                  0
                } else {
                  length(values_with_cover) / num_points_to_sample
                }
                
              })
              
            })
            
          }), paste0("sim", simulated_cover_values)
        )
      }),
      paste0("sd", sim_ratio_sd_to_mean * mean)
    )
    
  }),
  paste0("mean", sim_mean)
)

# 4. Evaluate outputs ####

# The sublist target will need to be updated if simulation is rerun and/or simulation parameters are changed
sub_list <- results$mean20$sd4$sim19.66

example_df <- bind_rows(
  lapply(1:40, function(x){
    
    bind_rows(
      
      lapply(1:21, function(y){
        
        tibble(number_points_sampled = x,
               replicate = y,
               cover = sub_list[[x]][[y]])
        
      })
      
    )
    
  })
) %>%
  mutate(cover = cover * 100)

example_df %>%
  group_by(number_points_sampled) %>%
  summarize(percent_cover = mean(cover)) %>%
  ggplot(aes(number_points_sampled, percent_cover)) + 
  geom_point() + 
  geom_line() +
  geom_hline(yintercept = 19.66) +
  ylim(0,50)


wilcoxen_test_df <- example_df %>%
  group_by(number_points_sampled) %>%
  summarize(percent_cover = mean(cover)) %>%
  mutate(parameter_mean = 19.66)

hist(wilcoxen_test_df$percent_cover)

wilcox.test(wilcoxen_test_df$parameter_mean, wilcoxen_test_df$percent_cover, exact = F)
