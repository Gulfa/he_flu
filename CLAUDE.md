# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Health economics and epidemiological modeling project evaluating influenza vaccination strategies in Norway. Models the impact of different vaccination uptake rates on hospital incidence, infections, costs, and QALYs lost.

## Setup

Install R dependencies (uses custom drat repository for FHI packages):
```r
source("install_reqs.R")
```

Key R packages: `metapop`, `metapopnorge` (custom FHI packages), `mcstate`, `dust`, `odin.dust`, `data.table`, `tidyverse`

Python dependencies (standard): `pandas`, `numpy`, `scipy`, `matplotlib`

## Running the Analysis

**Full workflow (sequential steps):**

1. **Generate immunity distributions** (Python):
   ```bash
   python simple_immunity_model/run_simple_immunity_model_scenarios.py
   ```
   This produces `parameter_files/{low,medium,high}_proportion_in_protection_compartments_by_age_group.csv`

2. **Run main epidemiological analysis** (R):
   ```r
   source("run_analysis.R")
   ```
   Fits transmission parameter (beta) via particle MCMC, then runs vaccination scenario simulations.

3. **Health economic calculations** (R):
   ```r
   source("QALY.R")
   ```

4. **Generate plots** (R):
   ```r
   source("plot_results.R")
   ```

## Architecture

### Data Pipeline

```
Vaccination data (SYSVAK - privacy protected, template only in repo)
     ↓
simple_immunity_model/ (Python) → protection distributions by age group
     ↓
run_analysis.R (R) → fits beta → runs 6 uptake scenarios × 3 seasons
     ↓
QALY.R → health economic outputs
     ↓
plot_results.R → publication figures
```

### Epidemiological Model

- SEIR-like compartments: S, E, A (asymptomatic), I (symptomatic), P, H (hospital), ICU, R, D
- **9 age groups:** 0-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
- **3 influenza seasons:** 2022/23, 2017/18, 2018/19 (calibration data: `nat_hosp_flu_*.csv`)
- **6 vaccination uptake scenarios per season:** baseline, target, US-HHS, 95%, 90%, 80%
- Waning immunity modeled as movement through protection compartments
- Statistical inference via Sequential Monte Carlo + particle MCMC (`mcstate` package)

### Immunity Model (`simple_immunity_model/`)

- Models population immunity from vaccination and prior infection
- Exponential decay waning: ~8-9%/month (hospitalization), ~8.5%/month (infection)
- Initial VE: 50-60% against infection, 43-60% against hospitalization
- Outputs 3 scenario CSVs (low/medium/high initial conditions) used by `run_analysis.R`
- **Note:** Real SYSVAK vaccination data excluded for privacy — only a template CSV is included

### Key Parameters

All vaccination parameters are in `parameter_files/parameters_vaccination.xlsx`.
Seasonality curve: `parameter_files/norm_pred.csv` (20% seasonal amplitude applied to transmission).
Population age distribution: `simple_immunity_model/age_FHI_2022.txt`.

### Health Economics (`QALY.R`)

- VSL (value of statistical life) = 33.35 BNOK
- Value per QALY = 1.53 BNOK
- Calculates: QALYs lost, productivity loss (work days), healthcare postponement impact
- Includes long COVID fraction by severity level
