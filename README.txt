README — Replication Package

Paper
Title: Google Searches and the Adoption of Electric Vehicles
Authors: José Ignacio Giménez-Nadal, Antonio Gutiérrez-Lythgoe, José Alberto Molina

This folder contains the dataset and Stata scripts required to reproduce the empirical results reported in the paper.

CONTENTS

Data
File: dataset.dta
Final country–month panel dataset used in all estimations.

The dataset includes:

EV adoption variables

Digital attention measures (Google Trends–based indices)

Macroeconomic controls

Fixed-effects identifiers

Additional auxiliary variables

Some variables included in the dataset are not used in the final estimations but are provided for completeness and transparency.

Code

master.do
Runs the full replication sequence.

01_analysis.do
Main empirical analysis:

Baseline models

Main regression tables

02_robustness.do
Robustness checks:

Alternative specifications

Robustness tables

HOW TO REPLICATE

Open Stata

Set the working directory to the replication folder

Run:
do master.do

All tables and figures will be reproduced automatically.

DATA CONSTRUCTION NOTES

The final dataset was constructed from multiple sources.

Digital attention measures were obtained from Google Trends using a Python pseudo-API. Raw query data were processed to ensure:

Cross-country comparability

Harmonized keyword translations

Consistent scaling

Time series were seasonally adjusted using the X13-ARIMA-SEATS methodology.

The Python code used for raw data retrieval is not included in the replication package because the procedure generates multiple intermediate files and requires a specific data-processing pipeline. 

Seasonal adjustment was performed using X13-ARIMA-SEATS. 