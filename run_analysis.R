source("model_functions.R")

setwd('/storage/vaccination/influenza_202308/')

param_file    <- "parameter_files/parameters_vaccination.xlsx"
output_folder <- "./results/"

# Average duration of immunity (days) - single parameter for simple exponential waning
T_avg <- 200

# Observed vaccination uptake by age group (9 groups: 0-9, 10-19, ..., 80+)
uptake.2223 <- c(0.01, .08, .18, .18, .18, .18, .40, .63, .63)
uptake.1718 <- c(10.5, 10.5, 10.5, 9.5, 8.9, 11.2, 21.9, 38.1, 46) / 100
uptake.1819 <- c(15, 15, 15.3, 12.5, 11, 15.9, 33.4, 48.7, 50.5) / 100

# +2% vaccination in 60+ age groups (indices 7-9: 60-69, 70-79, 80+)
boost_60plus <- c(0, 0, 0, 0, 0, 0, 0.02, 0.02, 0.02)

scenarios <- list(
  list(
    name                  = "22/23 Season",
    data_file             = "nat_hosp_flu_2223.csv",
    initial_vaccination   = uptake.2223,
    T_avg                 = T_avg,
    season                = 0.5,
    severity              = 1,
    include_new_variant   = NULL,
    include_vaccination   = FALSE,
    vaccination_uptake    = uptake.2223,
    vax_scenarios = list(
      list(name = "Baseline",      vaccination_uptake = uptake.2223),
      list(name = "+2% in 60+",    vaccination_uptake = pmin(1, uptake.2223 + boost_60plus))
    )
  ),
  list(
    name                  = "17/18 Season",
    data_file             = "nat_hosp_flu_1718.csv",
    initial_vaccination   = uptake.1718,
    T_avg                 = T_avg,
    season                = 0.5,
    severity              = 1,
    include_new_variant   = NULL,
    include_vaccination   = FALSE,
    vaccination_uptake    = uptake.1718,
    vax_scenarios = list(
      list(name = "Baseline",      vaccination_uptake = uptake.1718),
      list(name = "+2% in 60+",    vaccination_uptake = pmin(1, uptake.1718 + boost_60plus))
    )
  ),
  list(
    name                  = "18/19 Season",
    data_file             = "nat_hosp_flu_1819.csv",
    initial_vaccination   = uptake.1819,
    T_avg                 = T_avg,
    season                = 0.5,
    severity              = 1,
    include_new_variant   = NULL,
    include_vaccination   = FALSE,
    vaccination_uptake    = uptake.1819,
    vax_scenarios = list(
      list(name = "Baseline",      vaccination_uptake = uptake.1819),
      list(name = "+2% in 60+",    vaccination_uptake = pmin(1, uptake.1819 + boost_60plus))
    )
  )
)

normal <- lapply(scenarios, function(x) {
  run_scenarios(x, sims_per_best = 1, n_threads = 1, L = 239) %>% mutate(season = x$name)
})
save(normal, file = paste0(output_folder, "runs.RData"))
