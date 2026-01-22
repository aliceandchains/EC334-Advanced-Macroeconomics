*******************************************************
* Load data + preprocess once
*******************************************************

clear all
set more off

import delimited "alldata.csv", clear varnames(1)

* Time index
cap drop tq
gen tq = yq(year, quarter)
format tq %tq
sort tq
tsset tq

* Quick missingness audit (optional)
misstable summarize spf_inflation_1year michigan_1y_median swap_1year swap_5year swap_10year cpilfesl pceexcludingfoodandenergy pceall unrate nrou v_u poilbreusdm import_def

* Inflation measures (annualised q/q)
cap drop pi_core_cpi pi_core_pce pi_pce
gen pi_core_cpi = 400*(ln(cpilfesl) - ln(L.cpilfesl))
gen pi_core_pce = 400*(ln(pceexcludingfoodandenergy) - ln(L.pceexcludingfoodandenergy))
gen pi_pce      = 400*(ln(pceall) - ln(L.pceall))
label var pi_core_cpi "Core CPI inflation (annualised q/q, %)"
label var pi_core_pce "Core PCE inflation (annualised q/q, %)"
label var pi_pce      "Headline PCE inflation (annualised q/q, %)"

* Slack / tightness measures
cap drop u_gap l_vu
gen u_gap = unrate - nrou
gen l_vu  = ln(v_u)
label var u_gap "Unemployment gap (u - u*)"
label var l_vu  "Log tightness ln(v/u)"

* Supply / external price terms (annualised q/q)
cap drop d_oil pi_import
gen d_oil     = 400*(ln(poilbreusdm) - ln(L.poilbreusdm))
gen pi_import = 400*(ln(import_def) - ln(L.import_def))
label var d_oil     "Oil price inflation (annualised q/q, %)"
label var pi_import "Import price inflation (annualised q/q, %)"

* If you want a post dummy for later sections, keep it; otherwise drop it
cap drop post2020
gen post2020 = (tq >= yq(2020,1))
label var post2020 "Post-2020 dummy"

* Summary stats check (optional)
summarize pi_core_cpi pi_core_pce pi_pce u_gap l_vu spf_inflation_1year michigan_1y_median v_u d_oil pi_import


*******************************************************
* Baseline expectations-augmented Phillips Curve
* pi_t = b * E_t pi_{t+1} + g * gap_t + e_t
* Train: 1984Q1–2020Q2 ; Forecast: 2020Q3 onwards
*******************************************************

* Cutoffs
scalar tq1984q1 = yq(1984,1)
scalar tq2020q2 = yq(2020,2)
scalar tq2020q3 = yq(2020,3)

cap drop train_pc test_pc
gen train_pc = (tq >= tq1984q1 & tq <= tq2020q2)
gen test_pc  = (tq >= tq2020q3)

* Common sample for baseline PC
cap drop ok_pc
gen ok_pc = !missing(pi_pce, michigan_1y_median, u_gap)

* Estimate baseline PC on training sample
reg pi_pce michigan_1y_median u_gap if train_pc & ok_pc, robust
est store baseline_pc_8420q2

* Predicted values (fit + forecast)
cap drop pi_hat_pc pi_fit_pc pi_fc_pc
predict pi_hat_pc if ok_pc, xb
gen pi_fit_pc = pi_hat_pc if train_pc & ok_pc
gen pi_fc_pc  = pi_hat_pc if test_pc  & ok_pc

* Forecast errors in test window (>=2020Q3)
cap drop fe_pc se_pc ae_pc
gen fe_pc = pi_pce - pi_hat_pc if test_pc & ok_pc
gen se_pc = fe_pc^2 if test_pc & !missing(fe_pc)
gen ae_pc = abs(fe_pc) if test_pc & !missing(fe_pc)

quietly summarize se_pc
display "Baseline PC RMSE (forecast from 2020Q3) = " sqrt(r(mean))

quietly summarize ae_pc
display "Baseline PC MAE  (forecast from 2020Q3) = " r(mean)

* Graph: actual + fit + forecast (final style; no title)
twoway (line pi_pce tq if ok_pc & tq>=tq1984q1, lcolor(midblue) lwidth(medthick)) (line pi_fit_pc tq if !missing(pi_fit_pc) & tq>=tq1984q1, lcolor(cranberry) lwidth(medthick)) (line pi_fc_pc tq if !missing(pi_fc_pc) & tq>=tq1984q1, lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(gs10) lwidth(thin)) ytitle("Annualised quarterly inflation (pp)", size(small) color(black)) xtitle("") xscale(range(`=tq1984q1' .)) xlabel(`=tq1984q1'(40)`=tq(2024,4)', format(%tqCCYY!Qq) labsize(small) labcolor(black)) ylabel(, labsize(small) labcolor(black) glcolor(gs16) glwidth(vthin)) legend(order(1 "Actual headline PCE inflation" 2 "Fit (1984Q1–2020Q2)" 3 "Forecast (>=2020Q3)") cols(1) position(6) ring(1) size(small) color(black) region(lstyle(none))) graphregion(color(white)) plotregion(color(white)) name(pc_base, replace)

graph export "fig_baseline_pc_oos_forecast.pdf", replace








*******************************************************
* Hazell-style Phillips Curve
* Train: 1984Q1–2020Q2
* Forecast: 2020Q3 onwards
*******************************************************

* ----------------------------
* Time cutoffs
* ----------------------------
scalar tq1984q1 = yq(1984,1)
scalar tq2020q2 = yq(2020,2)
scalar tq2020q3 = yq(2020,3)

cap drop train_hz test_hz
gen train_hz = (tq >= tq1984q1 & tq <= tq2020q2)
gen test_hz  = (tq >= tq2020q3)

* ----------------------------
* Headline PCE inflation (annualised q/q)
* ----------------------------
cap drop pi_pce
gen pi_pce = 400*(ln(pceall) - ln(L.pceall))
label var pi_pce "Headline PCE inflation (annualised q/q, %)"

* ----------------------------
* Energy inflation (PCE if available, oil proxy otherwise)
* ----------------------------
cap drop pi_energy
capture confirm variable pce_energy
if _rc==0 {
    gen pi_energy = 400*(ln(pce_energy) - ln(L.pce_energy))
    label var pi_energy "PCE energy inflation (annualised q/q, %)"
}
else {
    gen pi_energy = 400*(ln(poilbreusdm) - ln(L.poilbreusdm))
    label var pi_energy "Oil price inflation proxy (annualised q/q, %)"
}

* ----------------------------
* Common sample
* ----------------------------
cap drop ok_hz
gen ok_hz = !missing(pi_pce, michigan_1y_median, u_gap, pi_energy)

* ----------------------------
* Estimate Phillips Curve (Hazell et al. style)
* ----------------------------
reg pi_pce michigan_1y_median u_gap pi_energy if train_hz & ok_hz, robust
est store hazell_8420q2

* ----------------------------
* Fitted values and forecasts
* ----------------------------
cap drop pi_hat_hz pi_fit_hz pi_fc_hz
predict pi_hat_hz if ok_hz, xb
gen pi_fit_hz = pi_hat_hz if train_hz & ok_hz
gen pi_fc_hz  = pi_hat_hz if test_hz  & ok_hz

* ----------------------------
* Forecast accuracy (out-of-sample)
* ----------------------------
cap drop fe_hz se_hz ae_hz
gen fe_hz = pi_pce - pi_hat_hz if test_hz & ok_hz
gen se_hz = fe_hz^2 if test_hz & !missing(fe_hz)
gen ae_hz = abs(fe_hz) if test_hz & !missing(fe_hz)

quietly summarize se_hz
display "RMSE (forecast from 2020Q3) = " sqrt(r(mean))

quietly summarize ae_hz
display "MAE  (forecast from 2020Q3) = " r(mean)

* ----------------------------
* Final dissertation-quality graph
* ----------------------------
twoway (line pi_pce tq if ok_hz & tq>=tq1984q1, lcolor(midblue) lwidth(medthick)) (line pi_fit_hz tq if !missing(pi_fit_hz) & tq>=tq1984q1, lcolor(cranberry) lwidth(medthick)) (line pi_fc_hz tq if !missing(pi_fc_hz) & tq>=tq1984q1, lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(gs10) lwidth(thin)) ytitle("Annualised quarterly inflation (pp)", size(small) color(black)) xtitle("") xscale(range(`=tq1984q1' .)) xlabel(`=tq1984q1'(40)`=tq(2024,4)', format(%tqCCYY!Qq) labsize(small) labcolor(black)) ylabel(, labsize(small) labcolor(black) glcolor(gs16) glwidth(vthin)) legend(order(1 "Actual headline PCE inflation" 2 "Fit (1984Q1–2020Q2)" 3 "Forecast (>=2020Q3)") cols(1) position(6) ring(1) size(small) color(black) region(lstyle(none))) graphregion(color(white)) plotregion(color(white)) name(hz20, replace)

graph export "fig_hazell_oos_forecast.pdf", replace






*******************************************************
* Robustness checks: Baseline PC with alternative inflation + expectations
* 2x2: (PCE vs CPI) x (Michigan 1y vs SPF 1y)
* Train: 1984Q1–2020Q2 ; Forecast: 2020Q3 onwards
* Output: 2x2 panels + bottom legend (NO extra plotted line)
* Style: midblue actual, cranberry fit, cranberry dashed forecast, black split line
* Panel labels: (a) (b) (c) (d) only
* Adds coefficient text: beta (expectations) and gamma (u_gap)
*******************************************************

tsset tq

* CPI inflation (using cpilfesl; replace with headline CPI if you have it)
cap drop pi_cpi
gen pi_cpi = 400*(ln(cpilfesl) - ln(L.cpilfesl))
label var pi_cpi "CPI inflation (annualised q/q, %)"

* Cutoffs
scalar tq1984q1 = yq(1984,1)
scalar tq2020q2 = yq(2020,2)
scalar tq2020q3 = yq(2020,3)

cap drop train_pc test_pc
gen train_pc = (tq >= tq1984q1 & tq <= tq2020q2)
gen test_pc  = (tq >= tq2020q3)

* X-axis formatting
local xlbl `=tq1984q1'(40)`=tq(2024,4)'
local xfmt %tqCCYY!Qq

* Where to print coefficients (x-position)
scalar tq_text = yq(1996,1)

*******************************************************
* Helper: compute y-position for coefficient text within each panel
*******************************************************
cap program drop _coefpos
program define _coefpos, rclass
    syntax varname [if]
    quietly summarize `varlist' `if', meanonly
    return scalar ytxt = r(max) - 0.08*(r(max)-r(min))
end

*******************************************************
* (a) PCE × Michigan
*******************************************************
cap drop ok1 pi_hat1 pi_fit1 pi_fc1
gen ok1 = !missing(pi_pce, michigan_1y_median, u_gap) & tq>=tq1984q1

reg pi_pce michigan_1y_median u_gap if train_pc & ok1, robust
est store m1_pce_mich
local b1 : display %5.2f _b[michigan_1y_median]
local g1 : display %5.2f _b[u_gap]

predict pi_hat1 if ok1, xb
gen pi_fit1 = pi_hat1 if train_pc & ok1
gen pi_fc1  = pi_hat1 if test_pc  & ok1

_coefpos pi_pce if ok1
scalar ytxt1 = r(ytxt)

twoway (line pi_pce tq if ok1, lcolor(midblue) lwidth(medthick)) (line pi_fit1 tq if !missing(pi_fit1), lcolor(cranberry) lwidth(medthick)) (line pi_fc1 tq if !missing(pi_fc1), lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(black) lwidth(thin)) title("(a)", size(small) color(black)) text(`=ytxt1' `=tq_text' "{it:β}=`b1'   {it:γ}=`g1'", size(small) color(black) placement(ne)) ytitle("Annualised quarterly inflation (pp)", size(small) color(black)) xtitle("") xscale(range(`=tq1984q1' .)) xlabel(`xlbl', format(`xfmt') labsize(small) labcolor(black)) ylabel(, labsize(small) labcolor(black) glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(Ga, replace)

*******************************************************
* (b) PCE × SPF
*******************************************************
cap drop ok2 pi_hat2 pi_fit2 pi_fc2
gen ok2 = !missing(pi_pce, spf_inflation_1year, u_gap) & tq>=tq1984q1

reg pi_pce spf_inflation_1year u_gap if train_pc & ok2, robust
est store m2_pce_spf
local b2 : display %5.2f _b[spf_inflation_1year]
local g2 : display %5.2f _b[u_gap]

predict pi_hat2 if ok2, xb
gen pi_fit2 = pi_hat2 if train_pc & ok2
gen pi_fc2  = pi_hat2 if test_pc  & ok2

_coefpos pi_pce if ok2
scalar ytxt2 = r(ytxt)

twoway (line pi_pce tq if ok2, lcolor(midblue) lwidth(medthick)) (line pi_fit2 tq if !missing(pi_fit2), lcolor(cranberry) lwidth(medthick)) (line pi_fc2 tq if !missing(pi_fc2), lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(black) lwidth(thin)) title("(b)", size(small) color(black)) text(`=ytxt2' `=tq_text' "{it:β}=`b2'   {it:γ}=`g2'", size(small) color(black) placement(ne)) ytitle("", size(small) color(black)) xtitle("") xscale(range(`=tq1984q1' .)) xlabel(`xlbl', format(`xfmt') labsize(small) labcolor(black)) ylabel(, labsize(small) labcolor(black) glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(Gb, replace)

*******************************************************
* (c) CPI × Michigan
*******************************************************
cap drop ok3 pi_hat3 pi_fit3 pi_fc3
gen ok3 = !missing(pi_cpi, michigan_1y_median, u_gap) & tq>=tq1984q1

reg pi_cpi michigan_1y_median u_gap if train_pc & ok3, robust
est store m3_cpi_mich
local b3 : display %5.2f _b[michigan_1y_median]
local g3 : display %5.2f _b[u_gap]

predict pi_hat3 if ok3, xb
gen pi_fit3 = pi_hat3 if train_pc & ok3
gen pi_fc3  = pi_hat3 if test_pc  & ok3

_coefpos pi_cpi if ok3
scalar ytxt3 = r(ytxt)

twoway (line pi_cpi tq if ok3, lcolor(midblue) lwidth(medthick)) (line pi_fit3 tq if !missing(pi_fit3), lcolor(cranberry) lwidth(medthick)) (line pi_fc3 tq if !missing(pi_fc3), lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(black) lwidth(thin)) title("(c)", size(small) color(black)) text(`=ytxt3' `=tq_text' "{it:β}=`b3'   {it:γ}=`g3'", size(small) color(black) placement(ne)) ytitle("Annualised quarterly inflation (pp)", size(small) color(black)) xtitle("") xscale(range(`=tq1984q1' .)) xlabel(`xlbl', format(`xfmt') labsize(small) labcolor(black)) ylabel(, labsize(small) labcolor(black) glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(Gc, replace)

*******************************************************
* (d) CPI × SPF
*******************************************************
cap drop ok4 pi_hat4 pi_fit4 pi_fc4
gen ok4 = !missing(pi_cpi, spf_inflation_1year, u_gap) & tq>=tq1984q1

reg pi_cpi spf_inflation_1year u_gap if train_pc & ok4, robust
est store m4_cpi_spf
local b4 : display %5.2f _b[spf_inflation_1year]
local g4 : display %5.2f _b[u_gap]

predict pi_hat4 if ok4, xb
gen pi_fit4 = pi_hat4 if train_pc & ok4
gen pi_fc4  = pi_hat4 if test_pc  & ok4

_coefpos pi_cpi if ok4
scalar ytxt4 = r(ytxt)

twoway (line pi_cpi tq if ok4, lcolor(midblue) lwidth(medthick)) (line pi_fit4 tq if !missing(pi_fit4), lcolor(cranberry) lwidth(medthick)) (line pi_fc4 tq if !missing(pi_fc4), lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(black) lwidth(thin)) title("(d)", size(small) color(black)) text(`=ytxt4' `=tq_text' "{it:β}=`b4'   {it:γ}=`g4'", size(small) color(black) placement(ne)) ytitle("", size(small) color(black)) xtitle("") xscale(range(`=tq1984q1' .)) xlabel(`xlbl', format(`xfmt') labsize(small) labcolor(black)) ylabel(, labsize(small) labcolor(black) glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(Gd, replace)

*******************************************************
* Build 2x2 panel figure
*******************************************************
graph combine Ga Gb Gc Gd, cols(2) imargin(2 2 2 2) graphregion(color(white)) name(Panels, replace)

*******************************************************
* Combine 4 panels and place the legend at the very bottom
*******************************************************
graph combine Ga Gb Gc Gd, cols(2) imargin(2 2 2 2) graphregion(color(white) margin(b+16)) name(fig_baseline_pc_robust_final, replace)

graph export "fig_baseline_pc_robust_final.pdf", replace
graph export "fig_baseline_pc_robust_final.png", replace width(3200)

graph combine Panels Leg, cols(1) imargin(0 0 0 0) graphregion(color(white)) name(fig_baseline_pc_robust_final, replace)

graph export "fig_baseline_pc_robust_final.pdf", replace
graph export "fig_baseline_pc_robust_final.png", replace width(3200)

*******************************************************
* Table of estimates (LaTeX)
*******************************************************
cap which esttab
if _rc ssc install estout
esttab m1_pce_mich m2_pce_spf m3_cpi_mich m4_cpi_spf using "tab_baseline_pc_robust.tex", replace se b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) label stats(N r2, labels("Observations" "R-squared")) compress









* Stage 2 - Baseline Phillips Curve: pre-2020 fit + post-2020 out-of-sample forecast *

* (Assumes Stage 1 already created: tq, train, test, pi_core_cpi, spf_inflation_1year, u_gap, and scalar tq2020q1)

* 1) Estimate baseline Phillips Curve on PRE-2020 (training) sample
reg pi_core_cpi spf_inflation_1year u_gap if train & !missing(spf_inflation_1year), robust

* 2) Generate fitted/forecast values using the same pre-2020 coefficients
predict pi_hat_all, xb

* 3) Forecast errors (available wherever RHS variables are non-missing)
gen fe_all = pi_core_cpi - pi_hat_all

* 4) Summary: how errors differ pre vs post
sum fe_all if train & !missing(pi_hat_all)
sum fe_all if test  & !missing(pi_hat_all)

* 5) Numeric forecast evaluation for post-2020 (RMSE + MAE)
gen se_oos = fe_all^2 if test & !missing(fe_all)
gen ae_oos = abs(fe_all) if test & !missing(fe_all)

sum se_oos
display "Post-2020 RMSE (baseline SPF PC, estimated pre-2020) = " sqrt(r(mean))

sum ae_oos
display "Post-2020 MAE  (baseline SPF PC, estimated pre-2020) = " r(mean)

* 6) Graph: actual inflation + pre-2020 fit (solid) + post-2020 forecast (dashed)
twoway (line pi_core_cpi tq if !missing(pi_core_cpi), lwidth(medthick)) (line pi_hat_all tq if train & !missing(pi_hat_all), lwidth(medthick)) (line pi_hat_all tq if test & !missing(pi_hat_all), lpattern(dash) lwidth(medthick)), xline(`=tq2020q1') legend(order(1 "Actual core CPI inflation" 2 "PC fit (pre-2020, in-sample)" 3 "PC forecast (post-2020, out-of-sample)")) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Phillips Curve: fit before COVID and forecast failure after") name(fig_pc_fit_forecast, replace)

graph export "fig_pc_fit_and_forecast.png", replace

* 7) Graph: forecast errors (full sample) with COVID cutoff
twoway (line fe_all tq if !missing(fe_all), lwidth(medthick)), xline(`=tq2020q1') yline(0) ytitle("Error: actual - model (pp)") xtitle("Quarter") title("Phillips Curve errors using pre-2020 estimates") name(fig_pc_errors, replace)

graph export "fig_pc_errors_prepost.png", replace










* Stage 2B - Robustness: alternative Phillips Curves (pre-2020 estimates, post-2020 forecasts)

* Define cutoffs
local cut = yq(2020,1)
local cut2 = yq(2020,2)

* Helper: drop 2020q2 if you want (Beaudry often excludes 2020q2 in residual plots)
* We will keep it for now, but you can activate this by adding: & tq!=`cut2' in the if conditions.

* -------------------------------
* SPEC 1: CPI + SPF + u_gap  (baseline)
* -------------------------------
cap drop pi_hat_s1 fe_s1 fit_s1 fc_s1
reg pi_core_cpi spf_inflation_1year u_gap if tq<`cut' & !missing(pi_core_cpi, spf_inflation_1year, u_gap), robust
predict pi_hat_s1 if !missing(pi_core_cpi, spf_inflation_1year, u_gap)
gen fe_s1 = pi_core_cpi - pi_hat_s1 if !missing(pi_core_cpi, spf_inflation_1year, u_gap)
gen fit_s1 = pi_hat_s1 if tq<`cut'
gen fc_s1  = pi_hat_s1 if tq>=`cut'

twoway (line pi_core_cpi tq if !missing(pi_core_cpi, spf_inflation_1year, u_gap), lwidth(medthick)) (line fit_s1 tq, lwidth(medthick)) (line fc_s1 tq, lpattern(dash) lwidth(medthick)), xline(`cut') legend(order(1 "Actual core CPI inflation" 2 "PC fit (pre-2020)" 3 "PC forecast (post-2020)")) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Robustness S1: CPI + SPF + u_gap")
graph export "robust_s1_cpi_spf_ugap_fit_forecast.png", replace

twoway (line fe_s1 tq if !missing(fe_s1), lwidth(medthick)), xline(`cut') yline(0) ytitle("Error: actual - model (pp)") xtitle("Quarter") title("Robustness S1 errors: CPI + SPF + u_gap")
graph export "robust_s1_cpi_spf_ugap_errors.png", replace

* -------------------------------
* SPEC 2: PCE + SPF + u_gap
* -------------------------------
cap drop pi_hat_s2 fe_s2 fit_s2 fc_s2
reg pi_core_pce spf_inflation_1year u_gap if tq<`cut' & !missing(pi_core_pce, spf_inflation_1year, u_gap), robust
predict pi_hat_s2 if !missing(pi_core_pce, spf_inflation_1year, u_gap)
gen fe_s2 = pi_core_pce - pi_hat_s2 if !missing(pi_core_pce, spf_inflation_1year, u_gap)
gen fit_s2 = pi_hat_s2 if tq<`cut'
gen fc_s2  = pi_hat_s2 if tq>=`cut'

twoway (line pi_core_pce tq if !missing(pi_core_pce, spf_inflation_1year, u_gap), lwidth(medthick)) (line fit_s2 tq, lwidth(medthick)) (line fc_s2 tq, lpattern(dash) lwidth(medthick)), xline(`cut') legend(order(1 "Actual core PCE inflation" 2 "PC fit (pre-2020)" 3 "PC forecast (post-2020)")) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Robustness S2: PCE + SPF + u_gap")
graph export "robust_s2_pce_spf_ugap_fit_forecast.png", replace

twoway (line fe_s2 tq if !missing(fe_s2), lwidth(medthick)), xline(`cut') yline(0) ytitle("Error: actual - model (pp)") xtitle("Quarter") title("Robustness S2 errors: PCE + SPF + u_gap")
graph export "robust_s2_pce_spf_ugap_errors.png", replace

* -------------------------------
* SPEC 3: CPI + Michigan + u_gap
* -------------------------------
cap drop pi_hat_s3 fe_s3 fit_s3 fc_s3
reg pi_core_cpi michigan_1y_median u_gap if tq<`cut' & !missing(pi_core_cpi, michigan_1y_median, u_gap), robust
predict pi_hat_s3 if !missing(pi_core_cpi, michigan_1y_median, u_gap)
gen fe_s3 = pi_core_cpi - pi_hat_s3 if !missing(pi_core_cpi, michigan_1y_median, u_gap)
gen fit_s3 = pi_hat_s3 if tq<`cut'
gen fc_s3  = pi_hat_s3 if tq>=`cut'

twoway (line pi_core_cpi tq if !missing(pi_core_cpi, michigan_1y_median, u_gap), lwidth(medthick)) (line fit_s3 tq, lwidth(medthick)) (line fc_s3 tq, lpattern(dash) lwidth(medthick)), xline(`cut') legend(order(1 "Actual core CPI inflation" 2 "PC fit (pre-2020)" 3 "PC forecast (post-2020)")) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Robustness S3: CPI + Michigan + u_gap")
graph export "robust_s3_cpi_mich_ugap_fit_forecast.png", replace

twoway (line fe_s3 tq if !missing(fe_s3), lwidth(medthick)), xline(`cut') yline(0) ytitle("Error: actual - model (pp)") xtitle("Quarter") title("Robustness S3 errors: CPI + Michigan + u_gap")
graph export "robust_s3_cpi_mich_ugap_errors.png", replace

* -------------------------------
* SPEC 4: CPI + SPF + v_u  (tightness proxy)
* -------------------------------
cap drop pi_hat_s4 fe_s4 fit_s4 fc_s4
reg pi_core_cpi spf_inflation_1year v_u if tq<`cut' & !missing(pi_core_cpi, spf_inflation_1year, v_u), robust
predict pi_hat_s4 if !missing(pi_core_cpi, spf_inflation_1year, v_u)
gen fe_s4 = pi_core_cpi - pi_hat_s4 if !missing(pi_core_cpi, spf_inflation_1year, v_u)
gen fit_s4 = pi_hat_s4 if tq<`cut'
gen fc_s4  = pi_hat_s4 if tq>=`cut'

twoway (line pi_core_cpi tq if !missing(pi_core_cpi, spf_inflation_1year, v_u), lwidth(medthick)) (line fit_s4 tq, lwidth(medthick)) (line fc_s4 tq, lpattern(dash) lwidth(medthick)), xline(`cut') legend(order(1 "Actual core CPI inflation" 2 "PC fit (pre-2020)" 3 "PC forecast (post-2020)")) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Robustness S4: CPI + SPF + v/u")
graph export "robust_s4_cpi_spf_vu_fit_forecast.png", replace

twoway (line fe_s4 tq if !missing(fe_s4), lwidth(medthick)), xline(`cut') yline(0) ytitle("Error: actual - model (pp)") xtitle("Quarter") title("Robustness S4 errors: CPI + SPF + v/u")
graph export "robust_s4_cpi_spf_vu_errors.png", replace











*******************************************************
* Stage 3 - Stock & Watson style changing PC slope
* (Full stage: constructs variables + estimates slopes + pretty plot)
* NO /// used anywhere
*******************************************************

* Ensure time series is set (Stage 1 should already do this)
* tsset tq

* 3.1 Construct Stock-Watson style inflation measure

* 4-quarter inflation rate (log change over 4 quarters, percent)
cap drop pi4_corecpi
gen pi4_corecpi = 100*(ln(cpilfesl) - ln(L4.cpilfesl))
label var pi4_corecpi "4-quarter core CPI inflation (%)"

* Backward-looking 4-quarter moving average of 4-quarter inflation
cap drop ma4_pi4_corecpi
gen ma4_pi4_corecpi = (pi4_corecpi + L.pi4_corecpi + L2.pi4_corecpi + L3.pi4_corecpi)/4
label var ma4_pi4_corecpi "MA(4) of 4-quarter inflation"

* Year-over-year change in smoothed inflation (4-quarter difference)
cap drop d_pi_sw
gen d_pi_sw = ma4_pi4_corecpi - L4.ma4_pi4_corecpi
label var d_pi_sw "YoY change in MA(4) inflation (pp)"

* 3.2 Unemployment gap (and MA(4))
cap drop u_gap ma4_ugap
gen u_gap = unrate - nrou
label var u_gap "Unemployment gap (u - u*)"

cap drop ma4_ugap
gen ma4_ugap = (u_gap + L.u_gap + L2.u_gap + L3.u_gap)/4
label var ma4_ugap "MA(4) unemployment gap"

* Valid observations (needs lags used in MA and in L4 of MA)
cap drop sw_ok
gen sw_ok = !missing(d_pi_sw, ma4_ugap, year)

* 3.3 Define subsamples (aligned to your data start in 1968)
cap drop period68_83 period84_99 period00_19
gen period68_83 = (year>=1968 & year<=1983)
gen period84_99 = (year>=1984 & year<=1999)
gen period00_19 = (year>=2000 & year<=2019)

* 3.4 Estimate reduced-form slopes ("Phillips correlation") and store
reg d_pi_sw ma4_ugap if sw_ok & period68_83, robust
est store s68_83

reg d_pi_sw ma4_ugap if sw_ok & period84_99, robust
est store s84_99

reg d_pi_sw ma4_ugap if sw_ok & period00_19, robust
est store s00_19

* Display slopes
est restore s68_83
display "Phillips correlation 1968–1983 = " _b[ma4_ugap]

est restore s84_99
display "Phillips correlation 1984–1999 = " _b[ma4_ugap]

est restore s00_19
display "Phillips correlation 2000–2019 = " _b[ma4_ugap]

* 3.5 Pretty plot (scatter + fitted lines)
set scheme s2color

local xmin = -2
local xmax = 6
local ymin = -4.5
local ymax = 6

local plots "(scatter d_pi_sw ma4_ugap if sw_ok & period68_83, msymbol(Oh) msize(medsmall) mcolor(navy) mlcolor(navy) mfcolor(none) mlwidth(medthick))"
local plots "`plots' (lfit d_pi_sw ma4_ugap if sw_ok & period68_83, lcolor(navy) lwidth(thick) lpattern(dash))"
local plots "`plots' (scatter d_pi_sw ma4_ugap if sw_ok & period84_99, msymbol(Sh) msize(medsmall) mcolor(forest_green) mlcolor(forest_green) mfcolor(none) mlwidth(medthick))"
local plots "`plots' (lfit d_pi_sw ma4_ugap if sw_ok & period84_99, lcolor(forest_green) lwidth(thick) lpattern(shortdash))"
local plots "`plots' (scatter d_pi_sw ma4_ugap if sw_ok & period00_19, msymbol(Dh) msize(medsmall) mcolor(black) mlcolor(black) mfcolor(none) mlwidth(medthick))"
local plots "`plots' (lfit d_pi_sw ma4_ugap if sw_ok & period00_19, lcolor(black) lwidth(thick) lpattern(solid))"

twoway `plots', legend(order(1 "1968–1983" 3 "1984–1999" 5 "2000–2019") position(3) ring(0) region(lstyle(none))) title("Changing Phillips Correlation (Reduced Form)", size(medsmall)) subtitle("Stock–Watson-style smoothing; subsample fits", size(small)) xtitle("Unemployment gap (MA(4), pp)", size(small)) ytitle("Year-over-year change in inflation (pp)", size(small)) xlabel(`xmin'(2)`xmax', grid) ylabel(`ymin'(2)`ymax', grid) xscale(range(`xmin' `xmax')) yscale(range(`ymin' `ymax')) graphregion(color(white)) plotregion(color(white)) note("Reduced-form accelerationist Phillips Curve (expectations not included). y = 4q change in MA(4) of 4q inflation; x = MA(4) unemployment gap.", size(vsmall)) name(fig_stage3_sw_pc_pretty, replace)

graph export "stage3_stockwatson_phillips_correlation_pretty.png", replace width(2400)

* 3.6 Showing that inflation expectations changed (decreased) over time acorss household, professional and market-based measures

set scheme s2color

twoway (line michigan_1y_median tq if !missing(michigan_1y_median), lcolor(navy) lwidth(thick)) (line spf_inflation_1year tq if !missing(spf_inflation_1year), lcolor(forest_green) lwidth(thick)) (line swap_1year tq if !missing(swap_1year), lcolor(black) lwidth(medthick) lpattern(dash)) (line swap_5year tq if !missing(swap_5year), lcolor(black) lwidth(thick) lpattern(solid)) (line swap_10year tq if !missing(swap_10year), lcolor(black) lwidth(medthick) lpattern(dot)), legend(order(1 "Michigan 1y median (households)" 2 "SPF 1y (professional forecasters)" 3 "1y inflation swap" 4 "5y inflation swap" 5 "10y inflation swap") position(6) ring(0) region(lstyle(none))) title("Inflation Expectations Over Time", size(medsmall)) subtitle("Surveys and market-based term structure", size(small)) ytitle("Expected inflation (%)", size(small)) xtitle("Quarter", size(small)) xline(`=yq(2020,1)', lpattern(dash) lcolor(gs8)) graphregion(color(white)) plotregion(color(white)) note("Market-based expectations reflect risk-neutral pricing. Vertical line marks 2020Q1.", size(vsmall)) name(fig_expectations_term_structure, replace)

graph export "fig_expectations_term_structure.png", replace width(2400)

*******************************************************
* 3.7 Counterfactual decomposition: holding blocks constant
* Model estimated pre-2020: pi_core_cpi = a + b*Epi + g*u_gap + so*d_oil + si*pi_import + e
* Then simulate fitted inflation holding one block constant at its pre-2020 mean
*******************************************************

* Define cutoff
local cut = yq(2020,1)

cap drop train3 test3 est3_ok
gen train3 = (tq < `cut')
gen test3  = (tq >= `cut')

* Choose expectations measure for this counterfactual exercise
local expvar spf_inflation_1year

* Ensure a clean common sample
gen est3_ok = !missing(pi_core_cpi, `expvar', u_gap, d_oil, pi_import)

* Estimate expectations-augmented PC using pre-2020 only
reg pi_core_cpi `expvar' u_gap d_oil pi_import if train3 & est3_ok, robust
est store stg3_pc_pre2020

* Store coefficients as scalars (so we can build counterfactual fits safely)
scalar b0 = _b[_cons]
scalar bE = _b[`expvar']
scalar bU = _b[u_gap]
scalar bO = _b[d_oil]
scalar bI = _b[pi_import]

* Pre-2020 means (constants)
quietly summarize `expvar' if train3 & est3_ok
scalar exp_bar = r(mean)

quietly summarize u_gap if train3 & est3_ok
scalar ugap_bar = r(mean)

quietly summarize d_oil if train3 & est3_ok
scalar oil_bar = r(mean)

quietly summarize pi_import if train3 & est3_ok
scalar imp_bar = r(mean)

* Baseline fitted inflation using actual series
cap drop pi_hat_base
gen pi_hat_base = .
replace pi_hat_base = b0 + bE*`expvar' + bU*u_gap + bO*d_oil + bI*pi_import if est3_ok

* Counterfactual A: hold expectations constant at exp_bar
cap drop pi_hat_expconst
gen pi_hat_expconst = .
replace pi_hat_expconst = b0 + bE*exp_bar + bU*u_gap + bO*d_oil + bI*pi_import if est3_ok

* Counterfactual B: hold slack constant at ugap_bar
cap drop pi_hat_slackconst
gen pi_hat_slackconst = .
replace pi_hat_slackconst = b0 + bE*`expvar' + bU*ugap_bar + bO*d_oil + bI*pi_import if est3_ok

* Counterfactual C: hold supply shocks constant at oil_bar and imp_bar
cap drop pi_hat_supplyconst
gen pi_hat_supplyconst = .
replace pi_hat_supplyconst = b0 + bE*`expvar' + bU*u_gap + bO*oil_bar + bI*imp_bar if est3_ok


*******************************************************
* 3.8 Numerical evaluation (post-2020): RMSE and MAE by counterfactual
*******************************************************

cap drop fe_base fe_exp fe_slack fe_supply
gen fe_base  = pi_core_cpi - pi_hat_base if test3 & est3_ok
gen fe_exp   = pi_core_cpi - pi_hat_expconst if test3 & est3_ok
gen fe_slack = pi_core_cpi - pi_hat_slackconst if test3 & est3_ok
gen fe_supply= pi_core_cpi - pi_hat_supplyconst if test3 & est3_ok

cap drop se_base se_exp se_slack se_supply ae_base ae_exp ae_slack ae_supply
gen se_base  = fe_base^2  if test3 & !missing(fe_base)
gen se_exp   = fe_exp^2   if test3 & !missing(fe_exp)
gen se_slack = fe_slack^2 if test3 & !missing(fe_slack)
gen se_supply= fe_supply^2 if test3 & !missing(fe_supply)

gen ae_base  = abs(fe_base)  if test3 & !missing(fe_base)
gen ae_exp   = abs(fe_exp)   if test3 & !missing(fe_exp)
gen ae_slack = abs(fe_slack) if test3 & !missing(fe_slack)
gen ae_supply= abs(fe_supply) if test3 & !missing(fe_supply)

quietly summarize se_base
scalar rmse_base = sqrt(r(mean))
quietly summarize se_exp
scalar rmse_exp = sqrt(r(mean))
quietly summarize se_slack
scalar rmse_slack = sqrt(r(mean))
quietly summarize se_supply
scalar rmse_supply = sqrt(r(mean))

quietly summarize ae_base
scalar mae_base = r(mean)
quietly summarize ae_exp
scalar mae_exp = r(mean)
quietly summarize ae_slack
scalar mae_slack = r(mean)
quietly summarize ae_supply
scalar mae_supply = r(mean)

display "=== Post-2020 forecast performance (Stage 3 counterfactuals) ==="
display "RMSE baseline (all actual)              = " rmse_base
display "RMSE hold expectations constant         = " rmse_exp
display "RMSE hold slack constant                = " rmse_slack
display "RMSE hold supply shocks constant        = " rmse_supply
display "MAE  baseline (all actual)              = " mae_base
display "MAE  hold expectations constant         = " mae_exp
display "MAE  hold slack constant                = " mae_slack
display "MAE  hold supply shocks constant        = " mae_supply

display "=== Improvements relative to baseline (lower is better) ==="
display "RMSE improvement: hold expectations = " (rmse_base - rmse_exp)
display "RMSE improvement: hold slack        = " (rmse_base - rmse_slack)
display "RMSE improvement: hold supply       = " (rmse_base - rmse_supply)


*******************************************************
* 3.9 Optional: one diagnostic graph (ONLY export if you want it as a figure)
* (Given your 5-figure limit, you may prefer NOT to export this.)
*******************************************************

cap drop pi_hat_base_in pi_hat_base_oos
gen pi_hat_base_in = pi_hat_base if train3 & est3_ok
gen pi_hat_base_oos = pi_hat_base if test3 & est3_ok

twoway (line pi_core_cpi tq if est3_ok, lwidth(medthick)) (line pi_hat_base_in tq if !missing(pi_hat_base_in), lwidth(medthick)) (line pi_hat_base_oos tq if !missing(pi_hat_base_oos), lpattern(dash) lwidth(medthick)) (line pi_hat_expconst tq if est3_ok, lpattern(shortdash) lwidth(medthick)) (line pi_hat_supplyconst tq if est3_ok, lpattern(dot) lwidth(medthick)), xline(`cut', lpattern(dash) lcolor(gs8)) legend(order(1 "Actual" 2 "Fit (pre-2020)" 3 "Forecast (post-2020)" 4 "Hold expectations const" 5 "Hold supply const") position(6) ring(0) region(lstyle(none))) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Counterfactuals: expectations vs supply (pre-2020 coefficients)") name(fig_stage3_counterfactuals, replace)

* Uncomment ONLY if you want to use this as one of your 5 figures:
* graph export "fig_stage3_counterfactuals.png", replace width(2400)











*******************************************************
* Stage 4 - Expectations: Volcker credibility and anchoring (SPF SR vs SPF LR)
* Core series:
*   Short-run expectations: spf_inflation_1year
*   Long-run expectations:  spf_inflation   (SPF "over one year")
* NO /// used anywhere
*******************************************************

* 4.0 Sanity: confirm long-run SPF exists (your file has spf_inflation)
capture confirm variable spf_inflation
if _rc!=0 {
    display "ERROR: spf_inflation not found. Check variable names with: describe spf*"
    exit 198
}

label var spf_inflation_1year "SPF expected inflation (1y, SR)"
label var spf_inflation        "SPF expected inflation (over 1y, LR)"

* 4.1 Define key regime dates (Volcker appointment is 1979Q3; disinflation credibility by mid-80s)
scalar tq1979q3 = yq(1979,3)
scalar tq1984q1 = yq(1984,1)
scalar tq2020q1 = yq(2020,1)

cap drop pre_volcker post_volcker post2020
gen pre_volcker  = (tq < tq1979q3)
gen post_volcker = (tq >= tq1984q1)
gen post2020     = (tq >= tq2020q1)

label var pre_volcker  "Pre-Volcker (before 1979Q3)"
label var post_volcker "Post-Volcker credibility era (>=1984Q1)"
label var post2020     "Post-2020 period (>=2020Q1)"

* 4.2 Construct anchoring diagnostics
* (a) Expectations gap: SR minus LR
cap drop exp_gap
gen exp_gap = spf_inflation_1year - spf_inflation
label var exp_gap "Expectations gap: SPF 1y minus SPF LR (pp)"

* (b) Rolling volatility of LR expectations (anchoring proxy): 5-year rolling SD (20 quarters)
cap drop lr_sd20
gen lr_sd20 = .
forvalues i=20/`=_N' {
    quietly summarize spf_inflation in `=`i'-19'/`i'
    replace lr_sd20 = r(sd) in `i'
}
label var lr_sd20 "Rolling SD(20q) of SPF LR expectations"

* (c) Rolling volatility of SR expectations (also useful)
cap drop sr_sd20
gen sr_sd20 = .
forvalues i=20/`=_N' {
    quietly summarize spf_inflation_1year in `=`i'-19'/`i'
    replace sr_sd20 = r(sd) in `i'
}
label var sr_sd20 "Rolling SD(20q) of SPF SR expectations"

* 4.3 Figure: SR vs LR expectations through time (this is your main expectations figure)
set scheme s2color
twoway (line spf_inflation_1year tq if !missing(spf_inflation_1year), lwidth(thick) lcolor(forest_green)) (line spf_inflation tq if !missing(spf_inflation), lwidth(thick) lcolor(navy)), legend(order(1 "SPF 1y (short-run)" 2 "SPF over 1y (long-run)") position(6) ring(0) region(lstyle(none))) title("SPF Inflation Expectations: Short-run vs Long-run", size(medsmall)) subtitle("Long-run anchoring after Volcker", size(small)) ytitle("Expected inflation (%)", size(small)) xtitle("Quarter", size(small)) xline(`=tq1979q3', lpattern(dash) lcolor(gs8)) xline(`=tq1984q1', lpattern(dash) lcolor(gs8)) xline(`=tq2020q1', lpattern(dash) lcolor(gs8)) note("Vertical lines: 1979Q3 (Volcker), 1984Q1 (post-disinflation regime), 2020Q1 (COVID).", size(vsmall)) graphregion(color(white)) plotregion(color(white)) name(fig_spf_sr_lr, replace)
graph export "fig_spf_sr_lr.png", replace width(2400)

* 4.4 Figure: expectations gap (SR - LR). This shows de-anchoring episodes vs re-anchoring.
twoway (line exp_gap tq if !missing(exp_gap), lwidth(medthick) lcolor(black)), yline(0, lcolor(gs10)) title("Expectations Gap: SR minus LR (SPF)", size(medsmall)) subtitle("Shock episodes widen the gap; anchoring implies mean reversion", size(small)) ytitle("SPF 1y - SPF LR (pp)", size(small)) xtitle("Quarter", size(small)) xline(`=tq1979q3', lpattern(dash) lcolor(gs8)) xline(`=tq1984q1', lpattern(dash) lcolor(gs8)) xline(`=tq2020q1', lpattern(dash) lcolor(gs8)) graphregion(color(white)) plotregion(color(white)) name(fig_exp_gap, replace)
graph export "fig_exp_gap.png", replace width(2400)

* 4.5 Table-style evidence: mean and volatility before/after Volcker (report in text or as a small table)
* Means
quietly summarize spf_inflation if pre_volcker & !missing(spf_inflation)
display "LR SPF mean pre-Volcker = " r(mean)
quietly summarize spf_inflation if post_volcker & !missing(spf_inflation)
display "LR SPF mean post-Volcker = " r(mean)

quietly summarize spf_inflation_1year if pre_volcker & !missing(spf_inflation_1year)
display "SR SPF mean pre-Volcker = " r(mean)
quietly summarize spf_inflation_1year if post_volcker & !missing(spf_inflation_1year)
display "SR SPF mean post-Volcker = " r(mean)

* Volatility (SD)
quietly summarize spf_inflation if pre_volcker & !missing(spf_inflation)
display "LR SPF SD pre-Volcker = " r(sd)
quietly summarize spf_inflation if post_volcker & !missing(spf_inflation)
display "LR SPF SD post-Volcker = " r(sd)

quietly summarize spf_inflation_1year if pre_volcker & !missing(spf_inflation_1year)
display "SR SPF SD pre-Volcker = " r(sd)
quietly summarize spf_inflation_1year if post_volcker & !missing(spf_inflation_1year)
display "SR SPF SD post-Volcker = " r(sd)

* 4.6 Anchoring test: does LR expectations respond less to current inflation after Volcker?
* Use lagged 4-quarter inflation from your Stage 3 construction if available; otherwise create quickly here.
capture confirm variable pi4_corecpi
if _rc!=0 {
    cap drop pi4_corecpi
    gen pi4_corecpi = 100*(ln(cpilfesl) - ln(L4.cpilfesl))
    label var pi4_corecpi "4-quarter core CPI inflation (%)"
}

* Regression: LR expectations on lagged inflation, pre vs post
reg spf_inflation L.pi4_corecpi if pre_volcker & !missing(spf_inflation, L.pi4_corecpi), robust
est store lr_pre

reg spf_inflation L.pi4_corecpi if post_volcker & !missing(spf_inflation, L.pi4_corecpi), robust
est store lr_post

est restore lr_pre
display "LR expectations sensitivity to lagged inflation (pre-Volcker) = " _b[L.pi4_corecpi]
est restore lr_post
display "LR expectations sensitivity to lagged inflation (post-Volcker) = " _b[L.pi4_corecpi]

* Optional pooled interaction test (clean evidence in one line)
cap drop DpostV
gen DpostV = post_volcker
reg spf_inflation L.pi4_corecpi c.L.pi4_corecpi#i.DpostV if (pre_volcker | post_volcker) & !missing(spf_inflation, L.pi4_corecpi), robust
test 1.DpostV#c.L.pi4_corecpi








*******************************************************
* Stage 5 - Expectations-augmented Phillips Curve (lecture-consistent)
* Goal: estimate NKPC-style regression with ONE expectations proxy at a time,
*       test stability of k, and compare post-2020 forecast performance
* NO /// used anywhere
*******************************************************

local cut = yq(2020,1)
scalar tq1984q1 = yq(1984,1)

cap drop train5 test5 post84 post20
gen train5 = (tq < `cut')
gen test5  = (tq >= `cut')
gen post84 = (tq >= tq1984q1)
gen post20 = (tq >= `cut')

label var post84 "Post-1984 dummy"
label var post20 "Post-2020 dummy"

* Supply controls (as in lectures when discussing cost-push shifts)
cap drop ok5
gen ok5 = !missing(pi_core_cpi, u_gap, d_oil, pi_import)

*******************************************************
* 5.1 SPEC A: Expectations = SPF 1-year (short-run)
*******************************************************

cap drop okA
gen okA = ok5 & !missing(spf_inflation_1year)

reg pi_core_cpi spf_inflation_1year u_gap d_oil pi_import if train5 & okA, robust
est store pc_sr_pre2020

cap drop pi_hat_sr fe_sr se_sr ae_sr
predict pi_hat_sr if okA, xb
gen fe_sr = pi_core_cpi - pi_hat_sr if okA

gen se_sr = fe_sr^2 if test5 & !missing(fe_sr)
gen ae_sr = abs(fe_sr) if test5 & !missing(fe_sr)

quietly summarize se_sr
scalar rmse_sr = sqrt(r(mean))
quietly summarize ae_sr
scalar mae_sr = r(mean)

display "SR-spec post-2020 RMSE = " rmse_sr
display "SR-spec post-2020 MAE  = " mae_sr

* k stability tests (within pre-2020 sample and at 2020 break)
reg pi_core_cpi spf_inflation_1year u_gap c.u_gap#i.post84 d_oil pi_import if train5 & okA, robust
est store pc_sr_k84
test 1.post84#c.u_gap

reg pi_core_cpi spf_inflation_1year u_gap c.u_gap#i.post20 d_oil pi_import if okA, robust
est store pc_sr_k20
test 1.post20#c.u_gap

*******************************************************
* 5.2 SPEC B: Expectations = SPF over one year (long-run)
*******************************************************

cap drop okB
gen okB = ok5 & !missing(spf_inflation)

reg pi_core_cpi spf_inflation u_gap d_oil pi_import if train5 & okB, robust
est store pc_lr_pre2020

cap drop pi_hat_lr fe_lr se_lr ae_lr
predict pi_hat_lr if okB, xb
gen fe_lr = pi_core_cpi - pi_hat_lr if okB

gen se_lr = fe_lr^2 if test5 & !missing(fe_lr)
gen ae_lr = abs(fe_lr) if test5 & !missing(fe_lr)

quietly summarize se_lr
scalar rmse_lr = sqrt(r(mean))
quietly summarize ae_lr
scalar mae_lr = r(mean)

display "LR-spec post-2020 RMSE = " rmse_lr
display "LR-spec post-2020 MAE  = " mae_lr

* k stability tests
reg pi_core_cpi spf_inflation u_gap c.u_gap#i.post84 d_oil pi_import if train5 & okB, robust
est store pc_lr_k84
test 1.post84#c.u_gap

reg pi_core_cpi spf_inflation u_gap c.u_gap#i.post20 d_oil pi_import if okB, robust
est store pc_lr_k20
test 1.post20#c.u_gap

*******************************************************
* 5.3 Small results table (RMSE/MAE) exported as CSV
*******************************************************

preserve
clear
set obs 2
gen spec = ""
replace spec = "Expectations = SPF 1y (SR)" in 1
replace spec = "Expectations = SPF over 1y (LR)" in 2
gen rmse_post2020 = .
gen mae_post2020  = .
replace rmse_post2020 = rmse_sr in 1
replace mae_post2020  = mae_sr  in 1
replace rmse_post2020 = rmse_lr in 2
replace mae_post2020  = mae_lr  in 2
format rmse_post2020 mae_post2020 %9.3f
list, noobs
export delimited using "table_stage5_forecast_compare.csv", replace
restore

*******************************************************
* 5.4 Optional figure: actual vs forecast (pick ONE spec if you need a figure)
*******************************************************

cap drop pi_hat_sr_in pi_hat_sr_oos
gen pi_hat_sr_in  = pi_hat_sr if train5 & okA
gen pi_hat_sr_oos = pi_hat_sr if test5  & okA

twoway (line pi_core_cpi tq if okA, lwidth(medthick)) (line pi_hat_sr_in tq if !missing(pi_hat_sr_in), lwidth(medthick)) (line pi_hat_sr_oos tq if !missing(pi_hat_sr_oos), lpattern(dash) lwidth(medthick)), xline(`cut', lpattern(dash) lcolor(gs8)) legend(order(1 "Actual" 2 "Fit (pre-2020)" 3 "Forecast (post-2020)") position(6) ring(0) region(lstyle(none))) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Expectations-augmented PC (SPF 1y): fit and post-2020 forecast") name(fig_stage5_sr_fit_forecast, replace)

graph export "fig_stage5_sr_fit_forecast.png", replace width(2400)


