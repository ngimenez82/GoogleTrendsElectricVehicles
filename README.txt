README — Replication Package

Paper
Title: Do Google Searches Help Forecast Electric Vehicle Registrations? Evidence from European Union Countries
Authors: José Ignacio Giménez-Nadal, Antonio Gutiérrez-Lythgoe, José Alberto Molina

This folder contains the dataset and the code required to reproduce all tables and figures in the paper.

CONTENTS

Data
  dataset.dta   Country–month panel (26 European countries, January 2018 – June 2025).
                The estimation sample contains 22 EU countries (CH, IS, NO and TR lack the consumer confidence index).
                Key variables:
                  ev_bev            monthly NEW registrations of battery electric vehicles (EAFO)
                  ln_bev, d_ln_bev  log of new BEV registrations and its monthly change (dependent variable)
                  lr_elec_vs_base   log-ratio of Google searches for electric cars vs conventional cars
                  d_att             monthly change in lr_elec_vs_base (attention indicator)
                  registered_cars   all new passenger car registrations (used for the BEV share)
                  unem_rate, cci, cpi_all, d_rel_price   macroeconomic controls (entered with one lag)
                Some variables are not used in the final estimations but are provided for completeness.

Python code (reference implementation; reproduces every number in the paper)
  python/common.py        Data loading, fixed-effects estimator with country-clustered SEs (identical to
                          Stata xtreg, fe vce(cluster)), dynamic CCE mean-group estimator (Chudik & Pesaran 2015).
  python/insample.py      Main estimates and robustness (wild cluster bootstrap, year and month-by-year FE,
                          CCE lags, jackknife, BEV share, placebo search indicators)   -> out_insample.json
  python/oos.py           Pseudo out-of-sample forecasts, h = 1, 3, 6 (Clark-West, Diebold-Mariano) -> out_oos.json
  python/oos_variants.py  Parsimonious specifications and alternative evaluation windows -> out_oos_variants.json
  python/tables.py        In-sample tables (summary statistics, FE/CCE estimates, appendix tables) -> out_tables.json
  python/r2gain.py        Country x month FE column of Table 2 and R-squared gain from attention    -> out_r2gain.json
  python/oos2.py          Main forecasting evaluation (Tables 4-5): benchmarks with alternative seasonality,
                          search-augmented models, Clark-West, Diebold-Mariano, Pesaran-Timmermann tests,
                          rolling window, alternative origins, h = 3, BEV share, subsidy months  -> out_oos2.json,
                                                                                                   oos2_h1_forecasts.csv
  python/figure1.py       Figure 1                                                                   -> figures/figure1.png
  python/figure2.py       Figure 2 (cumulative squared forecast error differences)                  -> figures/figure2.png

  Requirements: Python 3.10+, pandas, numpy, scipy, statsmodels, matplotlib.
  Run from the python/ folder, in this order:
      python insample.py
      python oos.py
      python oos_variants.py
      python tables.py
      python r2gain.py
      python oos2.py
      python figure1.py
      python figure2.py

Stata code
  00_master.do, 01_analysis.do, 02_robustness.do
  Reproduce Tables 1-3, A1 and A2 (requires outreg2, xtdcce2, winsor2).
  Set the working directory to this folder and run:  do 00_master.do
  Fixed-effects results coincide exactly with the Python implementation. CCE mean-group estimates from
  xtdcce2 differ slightly from the Python implementation because of differences in how the estimation
  sample for the cross-sectional averages is defined; the paper reports the Python estimates.

DATA CONSTRUCTION NOTES

Digital attention measures were obtained from Google Trends using a Python pseudo-API. Raw query data were
processed to ensure cross-country comparability, harmonized keyword translations and consistent scaling.
Time series were seasonally adjusted using X13-ARIMA-SEATS. The code used for raw data retrieval is not
included because it generates multiple intermediate files and requires a specific data-processing pipeline.

CHANGELOG (October 2026)
  - The dependent variable is described as the growth of NEW BEV registrations (not the installed stock).
  - Joint significance of the six attention lags is now tested with testparm (the previous code tested
    their sum).
  - Winsorization in 02_robustness.do now uses winsor2 ..., suffix(_w); the previous call left the
    variable un-winsorized.
  - Relative paths in 00_master.do.
  - Added Python implementation, out-of-sample forecasting evaluation, additional robustness checks and Figure 1.
  - Forecasting version: extended evaluation (oos2.py), Table 2 column with country x month FE (r2gain.py), Figure 2.
    Note: Google Trends series in dataset.dta were seasonally adjusted on the full sample, so the pseudo
    out-of-sample exercise is not fully real-time; this look-ahead favours the search-augmented models.
