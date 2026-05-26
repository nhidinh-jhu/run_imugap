# Test run ImuGap on California
library(imuGAP)
library(tidyverse)
library(parallel)
library(bayesplot)
library(ggplot2)
library(arrow)
# Set-up some options 
optimal_cores <- parallel::detectCores() - 2
stans_spec_init <- rep(list(list(beta_bs = rep(-2.94, 5),
                                 sigma_cnty = 1,
                                 sigma_sch = 2,
                                 off_cnty = rep(0,58),# Number of counties
                                 off_sch = rep(0, 7877), # Number of schools
                                 lambda_raw = c(3,3))), 4)
## Load data specific for imuGAP

### LOCATIONS
locations <- readRDS("output/locations.rds")
# Update column names to fit imuGAP requirements
locations <- locations %>%
  rename(loc_id = id)

### OBSERVATIONS
observations <- readRDS("output/obs.rds")
# Update column names to fit imuGAP requirements
observations <- observations %>%
  rename(obs_id = row_id)

### POPULATIONS
populations <- readRDS("output/obs_populations.rds")
# Update column names to fit imuGAP requirements
populations <- populations %>%
  rename(loc_id = location)
# Normalize birth cohort to start with 1
populations <- populations %>% 
  mutate(cohort = cohort - 1988)

## Canonicalize the data
locations <- canonicalize_locations(locations)
observations <- canonicalize_observations(observations, drop_extra = TRUE)
populations <- canonicalize_populations(populations,
                                        observations,
                                        locations,
                                        #max_cohort,
                                        #max_age,
                                        max_dose = 2L)
# Weight check
weight_check <- populations %>%
  group_by(obs_id) %>%
  summarise(weight_sum = sum(weight)) %>%
  filter(abs(weight_sum - 1) > 0.001)
## GOOD

## COHORT:
cohort_min <- min(populations$cohort)
populations <- populations %>%
  mutate(cohort = cohort - cohort_min + 1L)

# Verify
range(populations$cohort)  # should be 1 to 36

# How many observations per cohort year?
populations %>%
  group_by(cohort) %>%
  summarise(n = n()) %>%
  arrange(cohort)

## Run ImuGap
ca_outputs <- imuGAP::sampling(
  observations,
  populations,
  locations,
  imugap_opts = imugap_options(df = 5L, dose_schedule = c(1, 4), object = c("default")),
  stan_opts = stan_options(chains = 4, 
                           iter = 2000,
                           init = stans_spec_init,
                           cores = optimal_cores)
)

# Save the specific object to a file
saveRDS(ca_outputs, file = "ca_outputs.rds")

# Check output
class(ca_outputs)
names(ca_outputs)
head(names(ca_outputs), 50)
print(ca_outputs, pars = names(ca_outputs)[1:10])

# Diagnostic 
# 1. Trace Plots
## Extract the raw stanfit object to make things easier
fit <- ca_outputs$stanfit

## Target Trace Plots
# A. Check global spline coefficients (beta_bs) and overall variances
mcmc_trace(fit, pars = c("beta_bs[1]", "beta_bs[2]", "sigma_cnty", "sigma_sch"))

# B. Check a specific subset of county effects (e.g., counties 30 through 35)
mcmc_trace(fit, pars = paste0("off_cnty[", 30:35, "]"))

# C. Look at the vaccine uptake rate hyperparameter (lambda_raw)
mcmc_trace(fit, pars = c("lambda_raw[1]", "lambda_raw[2]"))

# D. Look at parameters created
mcmc_trace(fit, pars = c("p_obs[1]", "p_obs[2]"))
# Print a sample of the parameter names to verify "p_obs" is the true name
print(grep("^p_", names(fit), value = TRUE)[1:10])
# 2. Parameters check
# Look at the trace plots for the first 5 observations

## Predictions
### Create target population for prediction
target_population <- imuGAP::create_target(fit = ca_outputs, 
                                  location = populations$loc_id,
                                  age = populations$age, 
                                  cohort = populations$cohort,
                                  dose = populations$dose)

ca_prediction <- predict(ca_outputs$stanfit, target_population)
# Trouble-shoot mismatch from previous version 
all_objects <- ls(getNamespace("imuGAP"), all.names = TRUE)
grep("model_", all_objects, value = TRUE)

rstantools_model_impute_school_coverage_process_v6 <- get(
  "rstantools_model_impute_school_coverage_process_v6", 
  envir = getNamespace("imuGAP")
)

.__C__Rcpp_rstantools_model_impute_school_coverage_process_v6 <- get(
  ".__C__Rcpp_rstantools_model_impute_school_coverage_process_v6", 
  envir = getNamespace("imuGAP")
)                 

ca_prediction <- predict(ca_outputs, target_population)

## Export to Parquet format
library(arrow)
# Convert the prediction to a data frame if it's not already
if (!is.data.frame(ca_prediction)) {
  ca_prediction <- as.data.frame(ca_prediction)
}
# Write the prediction to a Parquet file
write_parquet(ca_prediction, "ca_prediction.parquet")

                                  