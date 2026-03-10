

options(dplyr.summarise.inform = FALSE)


#' Calculate productivity loss
#' 
#' @export
calc_prod_loss <- function(df, params, tisk_on=1){


  ld <- rep(c(2, 3, 3, 3, 4, 4, 4.5, 6, 6),2) # days of illness
  loss_per_inc_illness <- (c(565 + 285, 0 + 570, 1500,1500, 1500, 1500,1500,0,0) + 285)*1e-9*ld
  detection_prob_sympt <- 0.7
  detection_prob_asympt <- 0.3
  N_contacts <- 3
  day_iso <- 7
  day_quarantine <- 7
  loss_per_day_iso <- rep(c(565 + 285+285, 0 + 570+285, 750+285,750+285, 750+285, 750+285,750+285,285,285)*1e-9,2)



  df$tisked_daily_inc <- df$incidence*tisk_on

  frac <- df %>% group_by(sim) %>% summarize(fra = sum(tisked_daily_inc)/sum(incidence)) %>% pull(fra)


  frac_sympt <- c(params$sympt_frac[,1,1], params$sympt_frac[,2,1]*params$susceptibility_symp[,2,1] / params$susceptibility_asymp[,2,1])
#  print(t(as.matrix(df %>%filter(t==max(t)) %>% select(starts_with("tot_infected["))))*frac_sympt*loss-per_inc_illness)
  ret_df <- data.frame(sim=1:max(df$sim),
                       prod_loss_illness=colSums(t(as.matrix(df %>%filter(t==max(t)) %>% select(starts_with("tot_infected["))))*frac_sympt*loss_per_inc_illness),
                       prod_loss_isolation=(colSums(t(as.matrix(df %>%filter(t==max(t)) %>% select(starts_with("tot_infected["))))*(1-frac_sympt)*loss_per_day_iso*day_iso*detection_prob_asympt) + colSums(t(as.matrix(df %>%filter(t==max(t)) %>% select(starts_with("tot_infected["))))*frac_sympt*loss_per_day_iso*(day_iso-ld)*detection_prob_sympt))*frac,
                       prod_loss_quarantine=N_contacts*frac*colSums(t(as.matrix(df %>%filter(t==max(t)) %>% select(starts_with("tot_infected["))))*(
                         (1-frac_sympt)*detection_prob_asympt + frac_sympt*detection_prob_sympt)*loss_per_day_iso*day_quarantine)
                       )

  return(ret_df)
}


#' Calculate quality loss
#' @export
calc_qual_post_loss <- function(inc, prev, N_beds, N_incidence_acute_patients, N_incidence_elective_patients,
                                time_acute, time_elective, daily_loss_elective, quality_reduction_function, dt,res,
                                acute_age_dist=rep(0.2, 5),
                                covid_age_dist=rep(0,2, 5),
                                loss_per_life_year=1.4e-3,
                                loss_per_QALY=1.5e-3,
                                ret_full=F){

  time_acute <- time_acute/dt
  time_elective <- time_elective/dt
  N_incidence_acute_patients <- N_incidence_acute_patients*dt
  N_incidence_elective_patients <- N_incidence_elective_patients*dt
  i <- 1
  E_post <- 0
  A <-  N_incidence_acute_patients*time_acute
  E <- N_incidence_elective_patients*time_elective
  cost <- 0
  in_hosp <- A+E
  while(E_post[i] > 0 | i <= length(inc)){
    A <- c(A, A[i] + N_incidence_acute_patients - A[i]/time_acute)
    prev_cov <- 0
    inc_cov <- 0
    if(i < length(prev)) {
      prev_cov <- prev[i]
      inc_cov <- inc[i]
    }
    free_spaces <- N_beds - A[i] - prev_cov
    E_admitted <- N_incidence_elective_patients
    if(free_spaces < E_admitted*time_elective & free_spaces > 0){
      E_admitted <- floor(free_spaces/time_elective)
    }else if(free_spaces < 0) E_admitted <- 0
    E_post_admitted <- 0
    if(E_post[i] > 0 & free_spaces  > N_incidence_elective_patients*time_elective){
      E_post_admitted <- floor((free_spaces - N_incidence_elective_patients*time_elective)/time_elective)
      if(E_post_admitted > E_post[i]) E_post_admitted <- E_post[i]
    }
    E <- c(E, E[i] + E_admitted - E[i]/time_elective + E_post_admitted)
    E_post <- c(E_post, E_post[i] + (N_incidence_elective_patients - E_admitted) - E_post_admitted)
    in_hosp <- c(in_hosp, E[i] + A[i] + prev_cov)
    
    cost <- cost + E_post[i]*daily_loss_elective*dt*loss_per_QALY + (inc_cov*quality_reduction_function((in_hosp[i])/N_beds, age_dist=covid_age_dist) + N_incidence_acute_patients*quality_reduction_function((in_hosp[i])/N_beds, age_dist=acute_age_dist))*loss_per_life_year
    i <- i + 1
    if(i > 5000) break
  }
  if(ret_full){
    return(data.frame(1:i, E_post=E_post, prev_cov=c(prev, rep(0, i - length(prev))), in_hosp=in_hosp, cost=cost))
  }
  
  return(cost)
}

quality_redction_function_template <- function(x, max_loss=0.1, max_loss_frac=2){
  if(x < 1) return(0)
  if(x > max_loss_frac) return(max_loss)
  return( (x-1)/(max_loss_frac - 1)*max_loss)
  
}


calc_ll_quality <- function(df, los_other_resp=8, dt=1){

  hosp <- df %>% group_by(scenario, R, severity, vac_scenario, sim) %>% summarize(tot_hosp_days=sum(hosp)*dt) %>%
    mutate(ll_hosp=0.05*tot_hosp_days*(0.2 + tot_hosp_days/720000),
           ll_hosp_value=1.5e-3*ll_hosp)
  
  resp_in <- fread("input_files/ll_quality.csv")


  vals_u6 <- -resp_in[status=="Covid-vaksinert", helsetap_tap_per_dag]
  vals_o6 <- -resp_in[status=="Covid-uvaksinert", helsetap_tap_per_dag]
  if(max(df[, resp]) > length(vals_u6)){
    vals_u6 <- c(vals_u6, rep(max(vals_u6), max(df[, resp]) - max(resp_in[, prevalens_covid])))
    vals_o6 <- c(vals_o6, rep(max(vals_o6), max(df[, resp]) - max(resp_in[, prevalens_covid])))
  }




  data.table::setDT(df)

  df[severity <6, los_per_resp:=vals_u6[resp+1]]
  df[severity >=6, los_per_resp:=vals_o6[resp+1]]
  df[, Ft:=pmin(900, 5*resp-50)]
  df[Ft<0, Ft:=0]
  df[, ll_resp:=(resp + 140)*los_per_resp*dt]
  df[, ll_resp_value:=1.4e-3*ll_resp]
  icu <- df %>% group_by(scenario, R, severity, vac_scenario, sim) %>% summarize(ll_resp_value=sum(ll_resp_value),
                                                                                 tot_postponed=sum(Ft)*dt) %>%
    mutate(
      postponed_qaly=0.2*tot_postponed*(0.2 + tot_postponed/60000),
      icu_postponed_value=postponed_qaly*1.5e-3)
  
  return(hosp %>% left_join(icu, on=c("scenario"="scenario", "R"="R", "severity"="severity", "vac_scenario"="vac_scenario", "sim"="sim"))) 
}


scale_with_severity <- function(severity, max_severity, max_value){
  return(1 +  min((max_value - 1)*(severity-1)/(max_severity-1), max_value-1))
}


calc_qaly_inner = function(df_model, severity=1) {
  # Input data frame of the following format:
  # df_model_test = data.frame(
  #   "age_group" = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
  #   "cum_I" = c(4844.62,11098.25,12906.51,11275.82,10369.64,6542.56,3808.49,1382.99,540.57),
  #   "cum_H" = c(5.10   , 17.02  , 95.55  , 191.20 , 227.96 , 291.54 , 356.65 , 183.51 , 131.46 ),
  #   "cum_ICU" = c(0.00  , 0.65  , 3.07  , 14.78 , 24.92 , 45.78 , 71.09 , 41.04 , 13.00),
  #   "cum_D" = c(0.00, 0.03, 0.18, 0.79, 3.40, 11.81, 32.28, 53.88, 84.29)
  # )
                                        # Put it together
  max_severity <- 10
  severity_light_disease_max <-5
  severity_hosp_max <-3
  severity_icu_max <-1.5
  fraction_asymptomatic = 0.4 # Model assumption
  VSL = 33.35       / 1000 # BNOK
  VSLY = 1.47       / 1000 # BNOK
  value_QALY = 1.53 / 1000 # BNOK
  df_par = data.frame(
    "age_group" = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
    "length_ld" = c(2, 3, 3, 3, 4, 4, 4.5, 6, 6),
    "life_quality_loss_ld" = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.2, 0.25),
    "time_to_hosp" = c(2, 3, 3, 3, 4, 4, 4.5, 5, 5),
    "life_quality_loss_prehosp" = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.2, 0.25),
    "time_in_hosp" = c(1.5, 1.5, 2, 2, 2, 2.5, 2.5, 3, 3),
    "life_quality_loss_hosp" = c(0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4),
    "time_ill_posthosp" =  c(2, 3, 3, 3, 4, 4, 4.5, 5, 5),
    "life_quality_loss_posthosp" = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.2, 0.25),
    "icu_time_to_hosp" = c(2, 3, 3, 3, 4, 4, 4.5, 5, 5),
    "icu_life_quality_loss_prehosp" = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.2, 0.25),
    "icu_time_in_hosp_preicu" = c(1, 1, 2, 2, 2, 3, 3, 3, 3),
    "icu_life_quality_loss_hosp_preicu" = c(0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5),
    "icu_time_in_icu" = c(1, 1, 2, 5, 9, 8, 9, 8, 5),
    "icu_life_quality_loss_icu" = c(0.6, 0.6, 0.6, 0.6, 0.7, 0.7, 0.7, 0.7, 0.7),
    "icu_time_in_hosp_posticu" = c(1, 1, 6, 3, 6, 6, 6, 6, 1),
    "icu_life_quality_loss_hosp_posticu" = c(0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5),
    "icu_time_ill_posthosp" = c(2, 3, 4, 5, 6, 6, 6, 6, 4),
    "icu_life_quality_loss_ill_posthosp" = c(0.2, 0.25, 0.25, 0.25, 0.3, 0.3, 0.35, 0.35, 0.35),
    "frac_longcovid_ld" = c(0.04, 0.04, 0.05, 0.05, 0.04, 0.04, 0.5, 0.06, 0.07),
    "frac_longcovid_hosp" = c(0.04, 0.04, 0.08, 0.08, 0.05, 0.05, 0.06, 0.07, 0.08),
    "frac_longcovid_icu" = c(0.05, 0.08, 0.09, 0.09, 0.06, 0.07, 0.07, 0.08, 0.09),
    "lc_length_ld" = c(40, 30, 40, 40, 45, 50, 50, 70, 80),
    "lc_life_quality_loss_ld" = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.10, 0.12),
    "lc_time_in_hosp" = c(70, 60, 60, 70, 75, 80, 90, 100, 110),
    "lc_life_quality_loss_hosp" = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.12, 0.12, 0.12),
    "lc_time_in_icu" = c(80, 80, 80, 90, 95, 100, 120, 130, 140),
    "lc_life_quality_loss_icu" = c(0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4),
    "expected_remaining_life_years" = c(78.17, 68.22, 58.45, 48.69, 39.01, 29.59, 20.77, 12.91, 6.67),
    "expected_remaining_QALYs" = c(66.5, 57.3, 48.3, 39.6, 31.3, 23.5, 15.6, 9.8, 4.9)
  )
  ##### Calculations #####
  
  ## 1. Death
  # Select relevant columns from model output
  df_d = df_model %>% select("age_group", "cum_D")
  df_d$age_group = as.character(df_d$age_group)
  df_d = df_d %>% rename(death = cum_D)
  # Calculate value of lost lives (all numbers in BNOK)
  df_d$value_lost_lives = df_d$death * VSL
  df_d$lost_life_years = df_d$death * df_par$expected_remaining_life_years
  df_d$value_lost_life_years = df_d$lost_life_years * VSLY
  df_d$lost_QALYs_death = df_d$death * df_par$expected_remaining_QALYs
  df_d$value_lost_QALYs_death = df_d$lost_QALYs * value_QALY
  # # Add column sums as bottom row
  # df_d = df_d %>%
  #   bind_rows(summarise(.,
  #                       across(where(is.numeric), sum),
  #                       across(where(is.character), ~"Totalt")))
  
  
  ## 2. Disease
  df_s = df_model %>% select("age_group")
  df_s$sympt_inf = df_model$cum_I # Assume symptomatic 
  df_s$hosp = df_model$cum_H - df_model$cum_ICU # Only hospitalised outside ICU
  df_s$ICU = df_model$cum_ICU
  # Short term QALY loss:
  df_s$QALY_loss_light_disease = (df_s$sympt_inf - df_s$hosp - df_s$ICU) * df_par$length_ld * df_par$life_quality_loss_ld / 365 *
    scale_with_severity(severity, max_severity, severity_light_disease_max)
  df_s$QALY_loss_hosp = df_s$hosp * (df_par$time_to_hosp * df_par$life_quality_loss_prehosp 
                                    + df_par$time_in_hosp * df_par$life_quality_loss_hosp
                                    + df_par$time_ill_posthosp * df_par$life_quality_loss_posthosp
                                    ) / 365*
                                    scale_with_severity(severity, max_severity, severity_hosp_max)
  df_s$QALY_loss_ICU = df_s$ICU * (df_par$icu_time_to_hosp * df_par$icu_life_quality_loss_prehosp
                                   + df_par$icu_time_in_hosp_preicu * df_par$icu_life_quality_loss_hosp_preicu
                                   + df_par$icu_time_in_icu * df_par$icu_life_quality_loss_icu
                                   + df_par$icu_time_in_hosp_posticu * df_par$icu_life_quality_loss_hosp_posticu
                                   + df_par$icu_time_ill_posthosp * df_par$icu_life_quality_loss_ill_posthosp
                                   ) / 365*
    scale_with_severity(severity, max_severity, severity_icu_max)
  df_s$total_short_term_QALY_loss = df_s$QALY_loss_light_disease + df_s$QALY_loss_hosp + df_s$QALY_loss_ICU
  # Long term QALY loss ("long covid"):
  df_s$longcov_QALY_loss_light_disease = (df_s$sympt_inf - df_s$hosp - df_s$ICU) * df_par$frac_longcovid_ld * df_par$lc_length_ld * df_par$lc_life_quality_loss_ld / 365* scale_with_severity(severity, max_severity, severity_light_disease_max)
  df_s$longcov_QALY_loss_hosp = df_s$hosp * df_par$frac_longcovid_hosp * df_par$lc_time_in_hosp * df_par$lc_life_quality_loss_hosp / 365*
    scale_with_severity(severity, max_severity, severity_hosp_max)
  df_s$longcov_QALY_loss_icu = df_s$ICU * df_par$frac_longcovid_icu * df_par$lc_time_in_icu * df_par$lc_life_quality_loss_icu / 365*
    scale_with_severity(severity, max_severity, severity_icu_max)
  df_s$total_long_term_QALY_loss = df_s$longcov_QALY_loss_light_disease + df_s$longcov_QALY_loss_hosp + df_s$longcov_QALY_loss_icu
  df_s$lost_QALYs_disease = df_s$total_short_term_QALY_loss + df_s$total_long_term_QALY_loss
  df_s$value_lost_QALYs_disease = df_s$lost_QALYs_disease * value_QALY
  
  ## 3. Put together main results like bottom of excel sheet
  df_r = df_d %>% select("age_group", "lost_QALYs_death", "value_lost_QALYs_death", "death")
  df_r = cbind(df_r, df_s %>% select("lost_QALYs_disease", "value_lost_QALYs_disease", "sympt_inf", "hosp", "ICU",
                                     "QALY_loss_light_disease", "QALY_loss_hosp", "QALY_loss_ICU", "longcov_QALY_loss_light_disease",
                                     "longcov_QALY_loss_hosp", "longcov_QALY_loss_icu"))
  df_r$lost_QALYs_death_and_disease = df_r$lost_QALYs_death + df_r$lost_QALYs_disease
  df_r$value_lost_QALYs_death_and_disease = df_r$value_lost_QALYs_death + df_r$value_lost_QALYs_disease
  df_r$lost_QALYs_death_and_disease_per_symptomatic_case = df_r$lost_QALYs_death_and_disease / df_r$sympt_inf
  df_r$value_lost_QALYs_death_and_disease_per_symptomatic_case = df_r$value_lost_QALYs_death_and_disease / df_r$sympt_inf
  # TODO add columns for "QALY-tap per person som rammes av ulike alvorlighetsgrader" (2 tables in excel doc)
  
  # # Add column sums as bottom row
  ## df_r = df_r %>%
  ##   bind_rows(summarise(.,
  ##                       across(where(is.numeric), sum),
  ##                       across(where(is.character), ~"Totalt")))
  
  
  return(df_r)
}


#' Calculate QALYs
#' @export
calc_qaly <- function(d){

  asymp_frac <- c(0.47, 0.47, 0.32, 0.32, 0.32, 0.32, 0.20, 0.20, 0.20)
  out_vars <- c("tot_infected_age_1",
                "tot_infected_age_2",
                "tot_infected_age_3",
                "tot_infected_age_4",
                "tot_infected_age_5",
                "tot_infected_age_6",
                "tot_infected_age_7",
                "tot_infected_age_8",
                "tot_infected_age_9",
                "tot_hosp_age_1",
                "tot_hosp_age_2",
                "tot_hosp_age_3",
                "tot_hosp_age_4",
                "tot_hosp_age_5",
                "tot_hosp_age_6",
                "tot_hosp_age_7",
                "tot_hosp_age_8",
                "tot_hosp_age_9",
                "tot_resp_age_1",
                "tot_resp_age_2",
                "tot_resp_age_3",
                "tot_resp_age_4",
                "tot_resp_age_5",
                "tot_resp_age_6",
                "tot_resp_age_7",
                "tot_resp_age_8",
                "tot_resp_age_9",
                "D_age_1",
                "D_age_2",
                "D_age_3",
                "D_age_4",
                "D_age_5",
                "D_age_6",
                "D_age_7",
                "D_age_8",
                "D_age_9")
  scenario_vars <- c("scenario","R", "severity", "vac_scenario", "sim")
  
  a <- d %>% filter(t==max(d$t)) %>% select(c(all_of(out_vars), all_of(scenario_vars)))
  
  d[t==200, .(scenario,
              R,
              severity,
              vac_scenario
              )]


  
  b <- reshape2::melt(a, id.vars=c(all_of(scenario_vars)), measure.vars=out_vars, value.name="value") %>% mutate(age_group = rep(rep(c("0-9",
                                                                                                                                       "10-19",
                                                                                                                                       "20-29",
                                                                                                                                       "30-39",
                                                                                                                                       "40-49",
                                                                                                                                       "50-59",
                                                                                                                                       "60-69",
                                                                                                                                       "70-79",
                                                                                                                                       "80+"), each=nrow(a)),4))
  b$target = stringr::str_sub(b$variable, 1, -3)
  data.table::setDT(b)
  b[target=="tot_infected", value:=value*rep((1-asymp_frac), each=nrow(a))]
  
  c = reshape2::dcast(b ,  R+ severity + scenario + sim + vac_scenario+ age_group ~ target )
  c <- c %>% mutate(cum_D=D_age, cum_I=tot_infected_age, cum_ICU=tot_resp_age, cum_H=tot_hosp_age)
  
  out <- list()
  for( s in unique(c$severity)){
    sub <- c %>% filter(severity==s)
    out[[length(out) + 1]] <- cbind(sub%>% select(!age_group), calc_qaly_inner(sub, severity=s), by=c("age_group"="age_group"))
  }
  
  r <- data.table::rbindlist(out)
  tot <- r %>% group_by(R, severity, scenario, vac_scenario, sim) %>% summarize(lost_qualy=sum(lost_QALYs_death_and_disease), lost_value=sum(value_lost_QALYs_death_and_disease), cum_H=sum(cum_H), cum_D=sum(cum_D), cum_I=sum(cum_I), cum_R=sum(cum_ICU))

  return(tot)
}

expected_qaly_loss <- function(params, severity=1, vac_index=1){
  use_severity=min(ceiling(severity/8*4), 4)
  d <- data.frame(age_group = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
             cum_I=params$asympt_frac[,vac_index,],
             cum_H= params$hosp_prob[,vac_index,],
             cum_ICU=params$hosp_prob[,vac_index,]*params$icu_prob[,vac_index,],
             cum_D=(1 - params$hosp_prob[,vac_index,])*params$prob_death_non_hosp[,vac_index,] +
              params$hosp_prob[,vac_index,]*(1 - params$icu_prob[,vac_index,])*params$prob_death_hosp[,vac_index,] +
              params$hosp_prob[,vac_index,]*params$icu_prob[,vac_index,]*params$prob_death_icu[,vac_index,] 
             )

              
  return(calc_qaly_inner(d, severity = use_severity)$value_lost_QALYs_death_and_disease)
    

}


### Verification by using same data as in excel file -- seems to give correct results!
# df_model_test = data.frame(
#   "age_group" = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
#   "cum_I" = c(4844.62,11098.25,12906.51,11275.82,10369.64,6542.56,3808.49,1382.99,540.57),
#   "cum_H" = c(5.10   , 17.02  , 95.55  , 191.20 , 227.96 , 291.54 , 356.65 , 183.51 , 131.46 ),
#   "cum_ICU" = c(0.00  , 0.65  , 3.07  , 14.78 , 24.92 , 45.78 , 71.09 , 41.04 , 13.00),
#   "cum_D" = c(0.00, 0.03, 0.18, 0.79, 3.40, 11.81, 32.28, 53.88, 84.29)
# )
# df_out = calc_qaly(df_model_test)


# ##### Run through all scenarios and both models
# ### Load scenarios file ###############
# df_scenarios = read.csv("../../0-pre-processing/scenarios/scenarios_csv/scenarios_N.csv", sep=";")
# #######################################
# df_out = NA
# for (i in 1:nrow(df_scenarios)) {
#   Name_curr = df_scenarios[i, Name_curr]
#   print(Name_curr)
#   
# }



# write.csv(df_out, file="QALY_loss.csv")
