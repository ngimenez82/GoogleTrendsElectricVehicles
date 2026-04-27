clear all
set more off

global project "C:\Users\usuario\Desktop\projects\objetivo1_dga_final"
global scripts "$project\replication_GimenezNadaletal"

do "$scripts\01_analysis.do"
do "$scripts\02_robustness.do"