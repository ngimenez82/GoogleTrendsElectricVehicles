clear all
set more off
* Set this to the folder that contains the replication package
global scripts "."
global results "./results"
cap mkdir "$results"

do "$scripts/01_analysis.do"
do "$scripts/02_robustness.do"
