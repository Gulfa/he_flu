library(dplyr)
library(data.table)
library(metapop)
library(metapopnorge)
library(ggplot2)
library(zoo)

# Initial VE values (at time of vaccination)
VE_inf_max       <- 0.5   # VE against infection
VE_hosp_cond_max <- 0.2   # VE against hospitalisation given infection (conditional)


# Compute inter-compartment waning times for a simple exponential decay.
# Protection decays as VE(t) = VE_0 * exp(-t / T_avg).
# Returns the time spent in each compartment (days), ordered to match rr_hosp.
get_waning_times <- function(rr_hosp, T_avg) {
  VE <- rev(1 - rr_hosp)  # Decreasing from VE_hosp_cond_max to 0
  times <- vapply(seq_len(length(VE) - 1), function(i) {
    if (VE[i + 1] <= 0) Inf else T_avg * log(VE[i] / VE[i + 1])
  }, numeric(1))
  rev(times)
}


calc_seas <- function(date,
                      amount,
                      filename = "parameter_files/norm_pred.csv",
                      L = 365) {
  dat <- fread(filename)
  dat$out <- 1 - dat$seas * amount
  dat <- rbind(dat, dat)

  i_start <- as.numeric(date - lubridate::make_date(year(date), 1, 1))
  selection <- dat[i_start:(i_start + 365), ]
  selection$out <- selection$out / selection[1, out]

  if (L > 365) {
    selection <- rbind(selection, selection, selection)
    selection <- selection[1:L, ]
  }

  return(selection$out)
}


get_params <- function(initial_S_dist,
                       waning_time = rep(20, 11),
                       seas = 0,
                       rr_hosp = seq(0.95, 0.10, by = -0.05),
                       severity = 1,
                       include_new_variant = NULL,
                       dt = 2,
                       L = 365,
                       include_vaccination = FALSE,
                       vaccination_uptake = .7,
                       first_date = as.Date('2022-10-01')) {
  params <- read_param_file_simplified(param_file)
  n_vac <- length(rr_hosp)
  n_strain <- 1
  N <- 9
  seasonality <- calc_seas(first_date, amount = seas, L = L)
  seasonality <- rep(seasonality, each = 1 / dt)

  print(first_date, quote = F)

  # Estimating approximate starting conditions based on 2 new hospitalisations on first date
  inc <- 3 / 0.0025 / severity / 4

  I_ini  <- 0.6 * inc * 3
  P_ini  <- 0.6 * inc * 2
  A_ini  <- 0.4 * inc * 5
  Ea_ini <- 0.4 * inc * 2
  Es_ini <- 0.6 * inc * 2

  age_groups <- get_age_groups()

  S_ini <- round(age_groups * initial_S_dist)

  beta_strain <- rep(1, n_strain)

  # rr_inf derived analytically from rr_hosp assuming same exponential decay rate
  rr_inf <- 1 - (VE_inf_max / VE_hosp_cond_max) * (1 - rr_hosp)
  rr_inf[length(rr_inf)] <- 0  # just-vaccinated compartment: full protection against infection

  rr_death <- rr_hosp

  if (!is.null(include_new_variant)) {
    n_strain <- 2
    beta_strain <- c(1, include_new_variant$beta)
    rr_inf   <- c(rr_inf,   include_new_variant$rr_inf)
    rr_hosp  <- c(rr_hosp,  rep(include_new_variant$severity, 10))
    rr_death <- c(rr_death, rep(include_new_variant$severity, 10))
  }

  I_ini  <- array(round(I_ini / 90),  dim = c(N, n_vac, n_strain))
  Ea_ini <- array(round(Ea_ini / 90), dim = c(N, n_vac, n_strain))
  Es_ini <- array(round(Es_ini / 90), dim = c(N, n_vac, n_strain))
  A_ini  <- array(round(A_ini / 90),  dim = c(N, n_vac, n_strain))
  P_ini  <- array(round(P_ini / 90),  dim = c(N, n_vac, n_strain))

  if (n_strain == 2) {
    I_ini[,, 2]  <- round(I_ini[,, 2]  * include_new_variant$initial_frac)
    Ea_ini[,, 2] <- round(Ea_ini[,, 2] * include_new_variant$initial_frac)
    Es_ini[,, 2] <- round(Es_ini[,, 2] * include_new_variant$initial_frac)
    A_ini[,, 2]  <- round(A_ini[,, 2]  * include_new_variant$initial_frac)
    P_ini[,, 2]  <- round(P_ini[,, 2]  * include_new_variant$initial_frac)
  }

  vac_pars <- list(
    rr_inf       = rr_inf,
    rr_hosp      = rr_hosp,
    rr_death     = rr_death,
    rr_icu       = rr_hosp,
    rr_los_hosp  = rep(1, n_vac),
    rr_inf_asymp = rep(1, n_vac),
    rr_trans     = rep(1, n_vac)
  )

  vaccinations <- array(0, dim = c(L / dt, N, n_vac))
  if (include_vaccination) {
    # Fall vaccination start on 2023-10-01; 1/8 in each of next 8 weeks
    for (t in seq(
      as.integer(365 - (first_date - as.Date(paste0(year(first_date), "-10-01")))),
      as.integer(414 - (first_date - as.Date(paste0(year(first_date), "-10-01")))),
      by = 7
    ))
      vaccinations[t / dt, 1:9, 1] <- round(get_age_groups() * 0.125 * vaccination_uptake)
  }

  T_waning <- matrix(waning_time, nrow = N, ncol = n_vac, byrow = T)
  params <- c(params,
    list(
      N_steps          = L / dt,
      n_vac            = n_vac,
      n_strain         = n_strain,
      dt               = dt,
      T_waning         = T_waning,
      vaccinations     = vaccinations,
      beta_day         = matrix(1, ncol = N, nrow = L / dt),
      beta_strain      = beta_strain,
      cross_protection = matrix(0, ncol = n_strain, nrow = n_strain),
      n                = 9,
      S_ini            = S_ini,
      import_vec       = array(1, dim = c(L / dt, N, n_vac, n_strain)),
      I_ini            = I_ini,
      I_imp_ini        = array(0, dim = c(N, n_vac, n_strain)),
      Ea_ini           = Ea_ini,
      Es_ini           = Es_ini,
      A_ini            = A_ini,
      P_ini            = P_ini,
      H_ini            = array(0, dim = c(N, n_vac, n_strain)),
      ICU_H_ini        = array(0, dim = c(N, n_vac, n_strain)),
      ICU_R_ini        = array(0, dim = c(N, n_vac, n_strain)),
      ICU_P_ini        = array(0, dim = c(N, n_vac, n_strain)),
      B_D_ini          = array(0, dim = c(N, n_vac, n_strain)),
      B_D_H_ini        = array(0, dim = c(N, n_vac, n_strain)),
      B_D_ICU_ini      = array(0, dim = c(N, n_vac, n_strain)),
      R_ini            = array(0, dim = c(N, n_vac, n_strain)),
      D_ini            = array(0, dim = c(N, n_vac, n_strain)),
      tot_infected_ini = array(0, dim = c(N, n_vac, n_strain)),
      tot_hosp_ini     = array(0, dim = c(N, n_vac, n_strain)),
      tot_resp_ini     = array(0, dim = c(N, n_vac, n_strain)),
      tot_vac_ini      = array(0, dim = c(N, n_vac)),
      tot_vac_adm_ini  = array(0, dim = c(N, n_vac)),
      beta_norm        = age_groups,
      reg_pop_long     = age_groups,
      N_regions        = 1,
      waning_immunity_vax = array(1000, dim = c(N, n_vac, n_strain)),
      waning_inf       = 5,
      age_groups       = N,
      beta_mode        = 1,
      include_waning   = 2,
      vac_struct_length = n_vac,
      vax_type         = 2
    )
  )

  basic_params <- fix_params(params, N, n_vac, n_strain, vac_pars)
  params <- fix_beta_mode_params(basic_params)

  return(params)
}


redistribute_initial <- function(init_cond_file, rr_inf) {
  init_cond <- read.csv(init_cond_file)
  setDT(init_cond)

  n_vac       <- length(rr_inf)
  n_age_groups <- max(init_cond$age_group) - min(init_cond$age_group)

  new_cond <- data.table()
  new_cond[, age_group         := rep(min(init_cond$age_group):max(init_cond$age_group), each = n_vac)]
  new_cond[, protection_rounded := rep(1 - rr_inf, 1 + n_age_groups)]

  for (i in unique(init_cond$age_group)) {
    bins     <- seq(0, 0.999, 0.001)
    old_prop <- rep(init_cond[age_group == i, proportion / 100], each = 100)
    old_prot <- bins
    new_prop <- c()
    for (j in 1:(length(rr_inf) - 1)) {
      rr <- (1 - rev(rr_inf))[j]
      new_prop <- c(new_prop, sum(old_prop[which(old_prot >= rr)]))
      old_prop[which(old_prot >= rr)] <- 0
    }
    new_prop <- c(new_prop, 1 - sum(new_prop))
    new_cond[age_group == i, proportion := rev(new_prop)]
  }

  return(new_cond)
}


update_severity <- function(params, severity, change_icu_prob = 3.8, change_los_hosp = 2.3, base_sev = 10) {
  params$hosp_prob      <- params$hosp_prob * severity
  params$icu_prob       <- params$icu_prob * min((1 + (change_icu_prob - 1) * (severity - 1) / (base_sev - 1)), change_icu_prob)
  params$length_hosp    <- params$length_hosp * min((1 + (change_los_hosp - 1) * (severity - 1) / (base_sev - 1)), change_los_hosp)
  params$prob_death_hosp <- params$prob_death_icu <- params$prob_death_non_hosp <- params$prob_death_non_hosp * severity
  return(params)
}


run_scenarios <- function(scenario, n = 300, n_threads = 4, n_best = 10, n_parts = 10,
                          sims_per_best = 10, dt = 0.5, L = 219) {
  print(scenario$name)

  rr_hosp      <- seq(1, 0.8, -0.02)
  waning_time  <- c(1e10, get_waning_times(rr_hosp, T_avg = scenario$T_avg))
  rr_hosp      <- c(rr_hosp, 1)
  waning_time  <- c(waning_time, 4 * 31)

  initial_S <- matrix(0, nrow = 9, ncol = length(rr_hosp))
  initial_S[, 1]                    <- 1 - scenario$initial_vaccination
  initial_S[, length(rr_hosp) - 1] <- scenario$initial_vaccination

  new_params <- get_params(
    initial_S,
    waning_time          = waning_time,
    L                    = L,
    dt                   = dt,
    seas                 = scenario$season,
    severity             = scenario$severity,
    rr_hosp              = rr_hosp,
    include_new_variant  = scenario$include_new_variant,
    include_vaccination  = scenario$include_vaccination,
    vaccination_uptake   = scenario$vaccination_uptake,
    first_date           = min(fread(scenario$data_file)$date)
  )

  new_params <- update_severity(new_params, scenario$severity)

  param_sets <- fit_beta_filter(new_params, scenario$data_file, n_parts = 600, n_threads = 150, n_samples = 50)

  scenario$vax_scenarios[[length(scenario$vax_scenarios) + 1]] <-
    list(name = "Baseline", vaccination_uptake = scenario$initial_vaccination)

  runs <- list()
  for (vax_scenario in scenario$vax_scenarios) {
    i         <- 1
    initial_S <- matrix(0, nrow = 9, ncol = length(rr_hosp))
    initial_S[, 1]                    <- 1 - vax_scenario$vaccination_uptake
    initial_S[, length(rr_hosp) - 1] <- vax_scenario$vaccination_uptake
    age_groups <- get_age_groups()
    S_ini      <- round(age_groups * initial_S)
    new_ps     <- list()
    for (ps in param_sets) {
      new_ps[[length(new_ps) + 1]] <- modifyList(ps, list(S_ini = S_ini, name = paste("run", i)))
      i <- i + 1
    }
    r <- run_param_sets(
      new_ps,
      L                   = L,
      N_particles         = sims_per_best,
      N_threads_internal  = sims_per_best,
      n_threads,
      silent              = FALSE
    ) %>%
      mutate(
        sim  = paste(sim, name),
        name = vax_scenario$name,
        date = time + min(fread(scenario$data_file)$date)
      )
    runs[[length(runs) + 1]] <- r
  }

  r <- rbindlist(runs) %>% filter(time %% 1 == 0)

  if (!is.null(scenario$include_new_variant)) {
    r <- r %>% mutate(
      variant_beta         = scenario$include_new_variant$beta,
      variant_severity     = scenario$include_new_variant$severity,
      variant_initial_frac = scenario$include_new_variant$initial_frac,
      variant_rr           = min(scenario$include_new_variant$rr_inf)
    )
  } else {
    r <- r %>% mutate(
      variant_beta         = NA,
      variant_severity     = NA,
      variant_rr           = NA,
      variant_initial_frac = NA
    )
  }

  r$hosp_incidence <- apply(
    r[, .(`hosp_inc[1]`, `hosp_inc[2]`, `hosp_inc[3]`, `hosp_inc[4]`,
          `hosp_inc[5]`, `hosp_inc[6]`, `hosp_inc[7]`, `hosp_inc[8]`, `hosp_inc[9]`)],
    1, sum
  )

  return(subset(r, select = sub.vars))
}


calculate_costs <- function(r, t_ini = 65, asymp_frac = c(0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25)) {
  out_vars <- c(
    paste0("tot_infected_age_", 1:9),
    paste0("tot_hosp_age_", 1:9),
    paste0("tot_resp_age_", 1:9),
    paste0("D_age_", 1:9)
  )

  small_df <- r %>%
    filter(time %in% c(t_ini, max(r$time))) %>%
    select(c("time", all_of(out_vars)))

  diff     <- small_df %>% filter(time == max(r$time)) - small_df %>% filter(time == t_ini)
  diff$sim <- 1:nrow(diff)

  b <- reshape2::melt(diff, id.vars = "sim", measure.vars = out_vars, value.name = "value") %>%
    mutate(age_group = rep(rep(c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
                               each = nrow(diff)), 4))
  b$target <- stringr::str_sub(b$variable, 1, -3)
  data.table::setDT(b)

  b[target == "tot_infected_age", value := value * rep((1 - asymp_frac), each = nrow(diff))]

  c <- reshape2::dcast(b, sim + age_group ~ target)
  c <- c %>% mutate(cum_D = D_age, cum_I = tot_infected_age, cum_ICU = tot_resp_age, cum_H = tot_hosp_age)

  qaly_loss      <- calc_qaly_inner(c, severity = 1)
  qaly_loss$sim  <- c$sim

  ld                   <- rep(c(2, 3, 3, 3, 4, 4, 4.5, 6, 6), 2)
  loss_per_inc_illness <- (c(565 + 285, 0 + 570, 1500, 1500, 1500, 1500, 1500, 0, 0) + 285) * 1e-9 * ld

  prod_loss_by_age <- t(as.matrix(diff %>% select(starts_with("tot_infected_age")))) * (1 - asymp_frac) * loss_per_inc_illness
  prod_loss        <- colSums(prod_loss_by_age)
  prod_loss_by_age <- rowMeans(prod_loss_by_age) %>% as.vector()

  tot_qaly_per_sim <- qaly_loss %>%
    group_by(sim) %>%
    summarize(lost_qaly = sum(lost_QALYs_death_and_disease), lost_value_qaly = sum(value_lost_QALYs_death_and_disease)) %>%
    mutate(prod_loss = prod_loss, tot_loss = lost_value_qaly + prod_loss)

  overall <- tot_qaly_per_sim %>%
    reshape2::melt(id.vars = "sim") %>%
    data.table() %>%
    (function(x) x[, j = .(value = mean(value)), by = .(variable)])

  tot_qaly_by_age <- qaly_loss %>%
    group_by(age_group) %>%
    summarize(lost_qaly = mean(lost_QALYs_death_and_disease), lost_value_qaly = mean(value_lost_QALYs_death_and_disease)) %>%
    mutate(prod_loss = prod_loss_by_age, tot_loss = lost_value_qaly + prod_loss_by_age) %>%
    reshape2::melt(id.vars = "age_group") %>%
    data.table()

  tot_qaly_by_age$age_group <- paste0(tot_qaly_by_age$age_group, " years")

  tot_qaly_per_sim <- tot_qaly_per_sim %>%
    reshape2::melt(id.vars = "sim") %>%
    data.table()

  return(rbind(
    overall[, j = .(aggregation = "overall", strata = "overall", variable, value)],
    tot_qaly_by_age[, j = .(aggregation = "age", strata = age_group, variable, value)],
    tot_qaly_per_sim[, j = .(aggregation = "sim", strata = sim, variable, value)]
  ))
}


sub.vars <- c(
  "name", "time", "sim", "date",
  "tot_hosp_inc", "tot_infected", "tot_hosp", "tot_resp", "D",
  paste0("hosp_inc[", 1:9, "]"),
  paste0("tot_infected_age_", 1:9),
  paste0("tot_hosp_age_", 1:9),
  paste0("tot_resp_age_", 1:9),
  paste0("D_age_", 1:9),
  paste0("tot_vac_adm_age_", 1:9)
)

sub.vars.extended <- c(
  "name", "time", "sim",
  "hosp_incidence", "tot_infected", "tot_hosp", "tot_resp", "D",
  "tot_N", "I", "A", "incidence", "log_beta",
  paste0("hosp_inc[", 1:9, "]"),
  paste0("tot_infected_age_", 1:9),
  paste0("tot_hosp_age_", 1:9),
  paste0("tot_resp_age_", 1:9),
  paste0("D_age_", 1:9),
  paste0("tot_vac_adm_age_", 1:9),
  paste0("S_age_", 1:9), paste0("S_vac_", 1:9), paste0("N_vac_", 1:9),
  paste0("I[", 1:90, "]"), paste0("A[", 1:90, "]")
)


fixed_lag_smc <- function(filter, params, tot_steps, lag = 1) {
  state <- filter$run_begin(save_history = TRUE, pars = params)
  hists <- array(dim = c(6, filter$n_particles, tot_steps + 1 - lag))
  for (i in (lag + 1):(tot_steps + 1)) {
    state$step(i - 1)
    history_value <- state$history$value
    history_order <- state$history$order
    history_index <- state$history$index

    if (i == (lag + 1)) {
      hist  <- history_single(history_value, history_order, history_index, NULL, 1:(i - 1))
      hists <- array(dim = c(dim(hist)[1:2], tot_steps - lag))
      hists[,, 1:i - 1] <- hist
    } else {
      hist          <- history_single(history_value, history_order, history_index, NULL, i - 1 - lag)
      hists[,, i - 1 - lag] <- hist
    }
  }
  return(hists)
}


history_single <- function(history_value, history_order, history_index, index_particle, times = 1) {
  ny <- nrow(history_value)

  if (is.null(history_order)) {
    if (is.null(index_particle)) {
      ret <- history_value
    } else {
      ret <- history_value[, index_particle, , drop = FALSE]
    }
  } else {
    if (is.null(index_particle)) index_particle <- seq_len(ncol(history_value))

    np  <- length(index_particle)
    nt  <- length(times)
    idx <- matrix(NA_integer_, np, nt)

    for (i in rev(times)) {
      index_particle <- idx[, i - times[1] + 1] <- history_order[index_particle, i]
    }

    cidx <- cbind(seq_len(ny), rep(idx, each = ny), rep(seq_len(nt) + times[1] - 1, each = ny * np))
    ret  <- array(history_value[cidx], c(ny, np, nt))
  }
  rownames(ret) <- names(history_index)
  ret
}


fit_beta_filter <- function(params, input_file, n_parts = 500, n_threads = 100, n_samples = 30) {
  data     <- fread(input_file)
  hosp_inc <- data
  hosp_inc[, t := 1:nrow(hosp_inc)]
  beta_1   <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta = params$beta_day[1, ], use_eig = TRUE)
  params$beta_mode        <- 3
  params$rand_beta_sd     <- 0.17
  params$rand_beta_days   <- 7
  params$rand_beta_factors <- rep(beta_1, 9) * 0.5
  filter <- get_filter_tot_hosp(hosp_inc$hosp, params, n_particles = n_parts, n_threads = n_threads)

  filter$run(save_history = TRUE, pars = params)
  h <- filter$history()

  param_sets <- list()
  for (i in 1:n_samples) {
    beta <- exp(rep(h[9, i, ], each = 1 / params$dt)) * beta_1 * 0.5
    if (params$N_steps > length(beta)) {
      beta <- c(beta, rep(beta[length(beta)], params$N_steps - length(beta)))
    } else if (params$N_steps < length(beta)) {
      beta <- beta[1:params$N_steps]
    }
    param_sets[[i]] <- modifyList(params, list(beta_day = matrix(rep(beta, 9), ncol = 9), beta_mode = 1))
  }

  return(param_sets)
}
