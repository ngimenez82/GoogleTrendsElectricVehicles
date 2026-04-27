*============================================================*
* Winsorization setup (1-99) + Table 1 FE vs CCE
*============================================================*

* Needs:
* ssc install winsor2, replace
* ssc install outreg2, replace
* ssc install xtdcce2, replace   (if not already)

*-----------------------------*
* 0) Create winsorized variables (do NOT overwrite originals)
*-----------------------------*
preserve

* Keep only variables needed (optional but safer)
* keep id_iso2 tm d_ln_bev d_att unem_rate cci cpi_all d_rel_price mdum2-mdum12
local outpath "C:\Users\usuario\Desktop\projects\objetivo1_dga_final\results"

* Winsorize (1st and 99th percentiles)
cap which winsor2
if _rc ssc install winsor2, replace

cap drop d_ln_bev_w d_att_w cpi_all_w
gen d_ln_bev_w = d_ln_bev
winsor2 d_ln_bev_w, cuts(1 99)
* Optional (uncomment if you want)
* winsor2 d_att,     cuts(1 99) gen(d_att_w)
* winsor2 cpi_all,   cuts(1 99) gen(cpi_all_w)

* If you winsorize d_att / cpi_all, switch names below:
* local ATT d_att_w
* local CPI cpi_all_w
local ATT d_att
local CPI cpi_all

* Output file
local T1W "`outpath'\TableA1_diffusion_results_winsor.xls"

*============================================================*
* Table 1 (Winsor): Diffusion of EVs – FE vs CCE
*============================================================*

*-----------------------------*
* (1) FE: dynamics + macros (winsor DV)
*-----------------------------*
xtreg d_ln_bev_w L.d_ln_bev_w ///
    L.unem_rate L.cci L.`CPI' L.d_rel_price ///
    mdum2-mdum12 ///
    , fe vce(cluster id_iso2)

outreg2 using "`T1W'", replace word label dec(3) ///
    drop(mdum*) ///
    ctitle("(1) FE: Macro (winsor DV)") ///
    addtext("Country FE","Yes","Month FE","Yes","Attention","No","CCE","No","Winsor","DV 1-99")

*-----------------------------*
* (2) FE: + attention (winsor DV)
*-----------------------------*
xtreg d_ln_bev_w L.d_ln_bev_w ///
    L(1/6).`ATT' ///
    L.unem_rate L.cci L.`CPI' L.d_rel_price ///
    mdum2-mdum12 ///
    , fe vce(cluster id_iso2)

test L1.`ATT' + L2.`ATT' + L3.`ATT' + L4.`ATT' + L5.`ATT' + L6.`ATT' = 0
scalar p_joint = r(p)

lincom L1.`ATT' + L2.`ATT' + L3.`ATT' + L4.`ATT' + L5.`ATT' + L6.`ATT'
scalar b_sum = r(estimate)
scalar p_sum = r(p)

outreg2 using "`T1W'", append word label dec(3) ///
    drop(mdum*) ///
    ctitle("(2) FE: Macro + Attention (winsor DV)") ///
    addtext("Country FE","Yes","Month FE","Yes","Attention","Yes","CCE","No","Winsor","DV 1-99") ///
    addstat("p(joint attn lags)", p_joint, "Sum attn lags", b_sum, "p(sum)", p_sum)

*-----------------------------*
* (3) CCE: dynamics + macros (winsor DV)
*-----------------------------*
xtdcce2 d_ln_bev_w L.d_ln_bev_w ///
    L.unem_rate L.cci L.`CPI' L.d_rel_price ///
    mdum2-mdum12 ///
    , crosssectional(d_ln_bev_w) cr_lags(2)

outreg2 using "`T1W'", append word label dec(3) ///
    drop(mdum*) ///
    ctitle("(3) CCE: Macro (winsor DV)") ///
    addtext("Country intercepts","Yes","Month FE","Yes","Attention","No","CCE","Yes","CR lags","2","Winsor","DV 1-99")

*-----------------------------*
* (4) CCE: + attention (winsor DV)
*-----------------------------*
xtdcce2 d_ln_bev_w L.d_ln_bev_w ///
    L(1/6).`ATT' ///
    L.unem_rate L.cci L.`CPI' L.d_rel_price ///
    mdum2-mdum12 ///
    , crosssectional(d_ln_bev_w `ATT') cr_lags(2)

test L1.`ATT' + L2.`ATT' + L3.`ATT' + L4.`ATT' + L5.`ATT' + L6.`ATT' = 0
scalar p_joint = r(p)

lincom L1.`ATT' + L2.`ATT' + L3.`ATT' + L4.`ATT' + L5.`ATT' + L6.`ATT'
scalar b_sum = r(estimate)
scalar p_sum = r(p)

outreg2 using "`T1W'", append word label dec(3) ///
    drop(mdum*) ///
    ctitle("(4) CCE: Macro + Attention (winsor DV)") ///
    addtext("Country intercepts","Yes","Month FE","Yes","Attention","Yes","CCE","Yes","CR lags","2","Winsor","DV 1-99") ///
    addstat("p(joint attn lags)", p_joint, "Sum attn lags", b_sum, "p(sum)", p_sum)

restore
*============================================================*
* End winsor block
*============================================================*

*-----------------------------*
* Table A2 FE: + attention leads & lags (winsor DV)
*-----------------------------*
* FE baseline (no winsor): leads + lags attention
xtreg d_ln_bev L.d_ln_bev ///
    F(1/6).d_att L(1/6).d_att ///
    L.unem_rate L.cci L.cpi_all L.d_rel_price ///
    mdum2-mdum12 ///
    , fe vce(cluster id_iso2)

local outpath "C:\Users\usuario\Desktop\projects\objetivo1_dga_final\results"
local T1 "`outpath'\TableA2_forward_robustness.xls"

testparm F(1/6).d_att
scalar p_leads = r(p)

testparm L(1/6).d_att
scalar p_lags = r(p)

lincom L1.d_att + L2.d_att + L3.d_att + L4.d_att + L5.d_att + L6.d_att
scalar p_sum_lags = r(p)

outreg2 using "`T1'", replace word label dec(3) ///
    drop(mdum*) ///
    ctitle("(5) FE: Macro + Leads/Lags Attention (no winsor)") ///
    addtext("Country FE","Yes","Month FE","Yes","Attention","Leads+Lags","CCE","No") ///
    addstat("p(joint leads)", p_leads, "p(joint lags)", p_lags, ///
            "Sum lag attn",  p_sum_lags )

