library(dplyr)
library(data.table)
library(ggplot2)
source("cost_params.R")
source("QALY.R")
source("model_functions.R")

output_dir <- "results/"
load(paste0(output_dir, "runs.RData"))
r <- rbindlist(normal)

# Observed hospitalization data for each season
obs <- rbindlist(list(
  fread("nat_hosp_flu_2223.csv")[, season := "22/23 Season"],
  fread("nat_hosp_flu_1718.csv")[, season := "17/18 Season"],
  fread("nat_hosp_flu_1819.csv")[, season := "18/19 Season"]
))
obs[, date := as.Date(date)]

level <- 0.95
age_labels <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")


# ── 1. Hospital incidence over time ──────────────────────────────────────────
# Ribbon (stochastic) + line (deterministic) vs observed, by season and scenario

sum_time <- r %>%
  group_by(season, name, run_type, date) %>%
  summarise(
    hosp_inc_l = quantile(tot_hosp_inc, (1 - level) / 2),
    hosp_inc_m = median(tot_hosp_inc),
    hosp_inc_h = quantile(tot_hosp_inc, 1 - (1 - level) / 2),
    .groups = "drop"
  )

p <- ggplot(sum_time %>% filter(run_type == "stochastic"),
            aes(x = date, fill = name, colour = name)) +
  geom_ribbon(aes(ymin = hosp_inc_l, ymax = hosp_inc_h), alpha = 0.25, colour = NA) +
  geom_line(aes(y = hosp_inc_m), linewidth = 0.8)

if ("deterministic" %in% unique(r$run_type)) {
  p <- p + geom_line(
    data = sum_time %>% filter(run_type == "deterministic"),
    aes(y = hosp_inc_m), linetype = "dashed", linewidth = 0.8
  )
}

p <- p +
  geom_point(data = obs, aes(x = date, y = hosp), inherit.aes = FALSE,
             colour = "black", size = 1, alpha = 0.6) +
  facet_wrap(~season, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = c("Baseline" = "#3C78D8", "+2% in 60+" = "#E06C00")) +
  scale_colour_manual(values = c("Baseline" = "#3C78D8", "+2% in 60+" = "#E06C00")) +
  labs(x = NULL, y = "Daily hospital admissions",
       fill = "Scenario", colour = "Scenario",
       caption = "Ribbon: 95% CI (stochastic). Dashed: deterministic. Points: observed.") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

ggsave(paste0(output_dir, "hosp_incidence_time.png"), p, width = 14, height = 12)


# ── 2. Age-specific hospital incidence at peak ────────────────────────────────
# Total cumulative hospitalisations by age group, scenario, season

hosp_age_vars <- paste0("tot_hosp_age_", 1:9)

cum_age <- r %>%
  group_by(season, name, run_type, sim) %>%
  slice_max(time) %>%
  ungroup() %>%
  select(season, name, run_type, sim, all_of(hosp_age_vars))

cum_age_long <- cum_age %>%
  tidyr::pivot_longer(all_of(hosp_age_vars), names_to = "age_var", values_to = "hosp") %>%
  mutate(age_group = age_labels[as.integer(sub("tot_hosp_age_", "", age_var))])

cum_age_sum <- cum_age_long %>%
  group_by(season, name, run_type, age_group) %>%
  summarise(
    hosp_l = quantile(hosp, (1 - level) / 2),
    hosp_m = median(hosp),
    hosp_h = quantile(hosp, 1 - (1 - level) / 2),
    .groups = "drop"
  ) %>%
  mutate(age_group = factor(age_group, levels = age_labels))

p2 <- ggplot(cum_age_sum %>% filter(run_type == "stochastic"),
             aes(x = age_group, y = hosp_m, fill = name)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = hosp_l, ymax = hosp_h),
                position = position_dodge(0.8), width = 0.3) +
  facet_wrap(~season, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Baseline" = "#3C78D8", "+2% in 60+" = "#E06C00")) +
  labs(x = "Age group", y = "Cumulative hospitalisations",
       fill = "Scenario", caption = "Bars: median. Error bars: 95% CI.") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

ggsave(paste0(output_dir, "hosp_by_age.png"), p2, width = 12, height = 12)


# ── 3. Cumulative outcomes summary ────────────────────────────────────────────
# Total infections, hospitalisations, ICU, deaths by scenario and season

cum_overall <- r %>%
  group_by(season, name, run_type, sim) %>%
  slice_max(time) %>%
  ungroup() %>%
  select(season, name, run_type, sim, tot_infected, tot_hosp, tot_resp, D)

cum_long <- cum_overall %>%
  tidyr::pivot_longer(c(tot_infected, tot_hosp, tot_resp, D),
                      names_to = "outcome", values_to = "value") %>%
  mutate(outcome = recode(outcome,
    tot_infected = "Infections",
    tot_hosp     = "Hospitalisations",
    tot_resp     = "ICU admissions",
    D            = "Deaths"
  )) %>%
  mutate(outcome = factor(outcome,
    levels = c("Infections", "Hospitalisations", "ICU admissions", "Deaths")))

cum_sum <- cum_long %>%
  group_by(season, name, run_type, outcome) %>%
  summarise(
    val_l = quantile(value, (1 - level) / 2),
    val_m = median(value),
    val_h = quantile(value, 1 - (1 - level) / 2),
    .groups = "drop"
  )

p3 <- ggplot(cum_sum %>% filter(run_type == "stochastic"),
             aes(x = season, y = val_m, fill = name)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = val_l, ymax = val_h),
                position = position_dodge(0.8), width = 0.3) +
  facet_wrap(~outcome, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Baseline" = "#3C78D8", "+2% in 60+" = "#E06C00")) +
  labs(x = NULL, y = "Cumulative count", fill = "Scenario",
       caption = "Bars: median. Error bars: 95% CI.") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(paste0(output_dir, "cumulative_outcomes.png"), p3, width = 12, height = 10)


# ── 4. Impact of +2% vaccination: absolute difference ─────────────────────────
# (Baseline - +2% in 60+) median difference for each outcome and season

diff_sum <- cum_sum %>%
  filter(run_type == "stochastic") %>%
  select(season, name, outcome, val_l, val_m, val_h) %>%
  tidyr::pivot_wider(names_from = name,
                     values_from = c(val_l, val_m, val_h)) %>%
  mutate(
    diff_m = `val_m_Baseline`              - `val_m_+2% in 60+`,
    diff_l = `val_l_Baseline`              - `val_h_+2% in 60+`,
    diff_h = `val_h_Baseline`              - `val_l_+2% in 60+`
  )

p4 <- ggplot(diff_sum, aes(x = season, y = diff_m)) +
  geom_col(fill = "#2A9D8F", width = 0.6) +
  geom_errorbar(aes(ymin = diff_l, ymax = diff_h), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~outcome, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Cases averted (Baseline \u2212 +2% in 60+)",
       caption = "Positive = cases averted by increased vaccination. Bars: median. Error bars: 95% CI.") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(paste0(output_dir, "impact_plus2pct.png"), p4, width = 12, height = 10)


# ── 5. Cost difference table ──────────────────────────────────────────────────
# Run calculate_costs for each season × scenario × run_type combination,
# then compute per-sim differences (matched on sim ID) and summarise.

var_labels <- c(
  lost_qaly       = "QALY loss (QALYs)",
  lost_value_qaly = "Value of QALY loss (BNOK)",
  prod_loss       = "Productivity loss (BNOK)",
  tot_loss        = "Total economic loss (BNOK)"
)

fmt <- function(med, lo, hi, digits = 3) {
  sprintf("%s [%s, %s]",
          formatC(med, digits = digits, format = "f"),
          formatC(lo,  digits = digits, format = "f"),
          formatC(hi,  digits = digits, format = "f"))
}

cost_rows <- list()
for (s in unique(r$season)) {
  for (rt in unique(r$run_type)) {
    costs_base <- calculate_costs(r[season == s & name == "Baseline"    & run_type == rt])
    costs_plus <- calculate_costs(r[season == s & name == "+2% in 60+" & run_type == rt])

    # Per-sim costs (matched by strata = sim ID)
    sim_base <- costs_base[aggregation == "sim", .(strata, variable, base  = value)]
    sim_plus <- costs_plus[aggregation == "sim", .(strata, variable, plus  = value)]
    sim_diff <- merge(sim_base, sim_plus, by = c("strata", "variable"))
    sim_diff[, diff := base - plus]

    summary <- sim_diff[, .(
      base_m  = median(base),  base_l  = quantile(base, (1 - level) / 2),  base_h  = quantile(base, 1 - (1 - level) / 2),
      plus_m  = median(plus),  plus_l  = quantile(plus, (1 - level) / 2),  plus_h  = quantile(plus, 1 - (1 - level) / 2),
      diff_m  = median(diff),  diff_l  = quantile(diff, (1 - level) / 2),  diff_h  = quantile(diff, 1 - (1 - level) / 2)
    ), by = variable]

    summary[, `:=`(
      season   = s,
      run_type = rt,
      outcome  = var_labels[as.character(variable)],
      Baseline          = fmt(base_m, base_l, base_h),
      `+2% in 60+`      = fmt(plus_m, plus_l, plus_h),
      `Difference (averted)` = fmt(diff_m, diff_l, diff_h)
    )]
    cost_rows[[length(cost_rows) + 1]] <- summary
  }
}

cost_table <- rbindlist(cost_rows)[
  , .(season, run_type, outcome, Baseline, `+2% in 60+`, `Difference (averted)`)
][order(season, run_type, outcome)]

fwrite(cost_table, paste0(output_dir, "cost_differences.csv"))

# Print to console
cat("\n=== Cost differences: Baseline vs +2% vaccination in 60+ ===\n")
print(cost_table, row.names = FALSE)

message("Plots and table saved to ", output_dir)
