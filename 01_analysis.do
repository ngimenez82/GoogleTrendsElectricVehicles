*******************************************************
clear all
set more off

*******************************************************

* Open data 
 use "$scripts/dataset.dta", clear
 
* Choose output path  
local outpath "$results"
******************************************************
*============================================================*
* Table 1: DESCRIPTIVES
*============================================================*

outreg2 using "`outpath'\Table1_Descriptives.xls", ///
    replace word ///
    sum(log) ///
    keep(d_ln_bev d_att unem_rate cci cpi_all d_rel_price) ///
    label dec(3) ///
    title("Table 1. Summary Statistics")
xtset id_iso2 tm

*============================================================*
* Table 2: Diffusion of EVs – FE vs CCE
*============================================================*

local T2 "`outpath'\Table2.xls"

*-----------------------------*
* (1) FE: dinámica + macros
*-----------------------------*
xtreg d_ln_bev L.d_ln_bev ///
    L.unem_rate L.cci L.cpi_all L.d_rel_price ///
    mdum2-mdum12 ///
    , fe vce(cluster id_iso2)

outreg2 using "`T2'", replace word label dec(3) ///
    drop(mdum*) ///
    ctitle("(1) FE: Macro") ///
    addtext("Country FE","Yes","Month FE","Yes","Attention","No","CCE","No")

*-----------------------------*
* (2) FE: + atención
*-----------------------------*
xtreg d_ln_bev L.d_ln_bev ///
    L(1/6).d_att ///
    L.unem_rate L.cci L.cpi_all L.d_rel_price ///
    mdum2-mdum12 ///
    , fe vce(cluster id_iso2)

* FIX (Oct 2026): joint significance must test the six lags jointly, not their sum
testparm L(1/6).d_att
scalar p_joint = r(p)

lincom L1.d_att + L2.d_att + L3.d_att + L4.d_att + L5.d_att + L6.d_att
scalar b_sum = r(estimate)
scalar p_sum = r(p)

outreg2 using "`T2'", append word label dec(3) ///
    drop(mdum*) ///
    ctitle("(2) FE: Macro + Attention") ///
    addtext("Country FE","Yes","Month FE","Yes","Attention","Yes","CCE","No") ///
    addstat("p(joint attn lags)", p_joint, "Sum attn lags", b_sum, "p(sum)", p_sum)

*-----------------------------*
* (3) CCE: dinámica + macros
*-----------------------------*
xtdcce2 d_ln_bev L.d_ln_bev ///
    L.unem_rate L.cci L.cpi_all L.d_rel_price ///
    mdum2-mdum12 ///
    , crosssectional(d_ln_bev) cr_lags(2)

outreg2 using "`T2'", append word label dec(3) ///
    drop(mdum*) ///
    ctitle("(3) CCE: Macro") ///
    addtext("Country FE","Yes","Month FE","Yes","Attention","No","CCE","Yes","CR lags","2")

*-----------------------------*
* (4) CCE: + atención (baseline)
*-----------------------------*
xtdcce2 d_ln_bev L.d_ln_bev ///
    L(1/6).d_att ///
    L.unem_rate L.cci L.cpi_all L.d_rel_price ///
    mdum2-mdum12 ///
    , crosssectional(d_ln_bev d_att) cr_lags(2)

* FIX (Oct 2026): joint significance must test the six lags jointly, not their sum
testparm L(1/6).d_att
scalar p_joint = r(p)

lincom L1.d_att + L2.d_att + L3.d_att + L4.d_att + L5.d_att + L6.d_att
scalar b_sum = r(estimate)
scalar p_sum = r(p)

outreg2 using "`T2'", append word label dec(3) ///
    drop(mdum*) ///
    ctitle("(4) CCE: Macro + Attention") ///
    addtext("Country FE","Yes","Month FE","Yes","Attention","Yes","CCE","Yes","CR lags","2") ///
    addstat("p(joint attn lags)", p_joint, "Sum attn lags", b_sum, "p(sum)", p_sum)


	
*============================================================*
* Table 3: Robustness - Fuel demand & N2O emissions
*============================================================*

* Requisitos:
* ssc install xtdcce2, replace

local T3 "`outpath'\Table3_Robustness_Fuel_N2O.xls"

*-----------------------------*
* (0) Construcción de variables (con cuidado con ceros)
*-----------------------------*
* Fuel demand
cap drop l_gasdies d_l_gasdies
gen double l_gasdies = ln(gasdies_ktons) if gasdies_ktons>0
gen double d_l_gasdies = D.l_gasdies

* N2O emissions
cap drop l_n2o d_l_n2o
gen double l_n2o = ln(n2o_gg) if n2o_gg>0
gen double d_l_n2o = D.l_n2o

*-----------------------------*
* (1) Baseline CCE: Macro + Attention
*-----------------------------*
xtdcce2 d_ln_bev L.d_ln_bev ///
    L(1/6).d_att ///
    L.unem_rate L.cci L.cpi_all L.d_rel_price ///
    mdum2-mdum12 ///
    , crosssectional(d_ln_bev d_att) cr_lags(2)

* FIX (Oct 2026): joint significance must test the six lags jointly, not their sum
testparm L(1/6).d_att
scalar p_joint = r(p)

lincom L1.d_att + L2.d_att + L3.d_att + L4.d_att + L5.d_att + L6.d_att
scalar b_sum = r(estimate)
scalar p_sum = r(p)

outreg2 using "`T3'", replace word label dec(3) ///
    drop(mdum*) ///
    ctitle("(1) Baseline CCE") ///
    addtext("Country intercepts","Yes","Month FE","Yes","CCE","Yes","CR lags","2") ///
    addstat("p(joint attn lags)", p_joint, "Sum attn lags", b_sum, "p(sum)", p_sum)

*-----------------------------*
* (2) + Fuel demand (Δ ln gas/diesel)
*-----------------------------*
xtdcce2 d_ln_bev L.d_ln_bev ///
    L(1/6).d_att ///
    L.unem_rate L.cci L.cpi_all L.d_rel_price ///
    L.d_l_gasdies ///
    mdum2-mdum12 ///
    , crosssectional(d_ln_bev d_att) cr_lags(2)

* FIX (Oct 2026): joint significance must test the six lags jointly, not their sum
testparm L(1/6).d_att
scalar p_joint = r(p)

lincom L1.d_att + L2.d_att + L3.d_att + L4.d_att + L5.d_att + L6.d_att
scalar b_sum = r(estimate)
scalar p_sum = r(p)

outreg2 using "`T3'", append word label dec(3) ///
    drop(mdum*) ///
    ctitle("(2) + Fuel demand") ///
    addtext("Country intercepts","Yes","Month FE","Yes","CCE","Yes","CR lags","2") ///
    addstat("p(joint attn lags)", p_joint, "Sum attn lags", b_sum, "p(sum)", p_sum)

*-----------------------------*
* (3) + N2O emissions (Δ ln N2O)
*-----------------------------*
xtdcce2 d_ln_bev L.d_ln_bev ///
    L(1/6).d_att ///
    L.unem_rate L.cci L.cpi_all L.d_rel_price ///
    L.d_l_n2o ///
    mdum2-mdum12 ///
    , crosssectional(d_ln_bev d_att) cr_lags(2)

* FIX (Oct 2026): joint significance must test the six lags jointly, not their sum
testparm L(1/6).d_att
scalar p_joint = r(p)

lincom L1.d_att + L2.d_att + L3.d_att + L4.d_att + L5.d_att + L6.d_att
scalar b_sum = r(estimate)
scalar p_sum = r(p)

outreg2 using "`T3'", append word label dec(3) ///
    drop(mdum*) ///
    ctitle("(3) + N2O emissions") ///
    addtext("Country intercepts","Yes","Month FE","Yes","CCE","Yes","CR lags","2") ///
    addstat("p(joint attn lags)", p_joint, "Sum attn lags", b_sum, "p(sum)", p_sum)

*============================================================*
* End
*============================================================*
	