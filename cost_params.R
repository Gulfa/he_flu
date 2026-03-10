cost_params <- list(

  # Time step at which to start accumulating costs (burn-in period)
  t_ini = 65,

  # Fraction of infections that are asymptomatic, by age group
  # (0-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+)
  asymp_frac = c(0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25),

  # Duration of symptomatic illness by age group (days)
  illness_days = c(2, 3, 3, 3, 4, 4, 4.5, 6, 6),

  # Daily income loss per symptomatic case by age group (NOK)
  # 0-9:  parental sick leave (565 absence + 285 part-time loss)
  # 10-19: part-time work (570)
  # 20-69: full-time work (1500)
  # 70+:  retired (0)
  daily_income_loss_by_age = c(565 + 285, 0 + 570, 1500, 1500, 1500, 1500, 1500, 0, 0),

  # Additional daily overhead cost per illness across all ages (NOK)
  # e.g. healthcare contacts, out-of-pocket costs
  daily_overhead_nok = 285,

  # Conversion factor: NOK -> BNOK
  nok_to_bnok = 1e-9,

  # Severity multiplier passed to calc_qaly_inner (1 = baseline influenza)
  qaly_severity = 1
)
