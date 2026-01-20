* Stage 1 - understanding the data and preparing for analysis *

clear all
set more off

import delimited "alldata.csv", clear varnames(1)

misstable summarize spf_inflation_1year michigan_1y_median swap_1year swap_5year swap_10year cpilfesl pceexcludingfoodandenergy pceall unrate nrou v_u poilbreusdm import_def

gen tq = yq(year, quarter)
format tq %tq
sort tq
tsset tq

gen pi_core_cpi = 400*(ln(cpilfesl) - ln(L.cpilfesl))
gen pi_core_pce = 400*(ln(pceexcludingfoodandenergy) - ln(L.pceexcludingfoodandenergy))
label var pi_core_cpi "Core CPI inflation (annualised q/q, %)"
label var pi_core_pce "Core PCE inflation (annualised q/q, %)"

gen u_gap = unrate - nrou
label var u_gap "Unemployment gap (u - u*)"

gen d_oil = 400*(ln(poilbreusdm) - ln(L.poilbreusdm))
gen pi_import = 400*(ln(import_def) - ln(L.import_def))
label var d_oil "Oil price inflation (annualised q/q, %)"
label var pi_import "Import price inflation (annualised q/q, %)"

gen post2020 = (tq >= yq(2020,1))
label var post2020 "Post-2020 dummy"

gen l_vu = ln(v_u)
label var l_vu "Log tightness ln(v/u)"

gen train = (tq <= yq(2019,4))
gen test  = (tq >= yq(2020,1))
label var train "Training sample (<=2019Q4)"
label var test "Test sample (>=2020Q1)"

scalar tq2020q1 = yq(2020,1)

summarize pi_core_cpi pi_core_pce u_gap l_vu spf_inflation_1year michigan_1y_median v_u d_oil pi_import











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















* Stage 3 – Baseline Phillips Curve with added supply shocks (incorporating oil prices and import prices as well) *



* Augmented Phillips Curve with supply shocks
reg pi_core_cpi spf_inflation_1year u_gap d_oil pi_import if !missing(spf_inflation_1year), robust

* Fitted values from augmented Phillips Curve
predict pi_hat_pc_supply if e(sample)

* Forecast errors
gen fe_pc_supply = pi_core_cpi - pi_hat_pc_supply if e(sample)

* Graph 3 – redrawing of the graph 1 *
twoway (line pi_core_cpi tq if e(sample)) (line pi_hat_pc_supply tq if e(sample)), xline(`=q2020q1') legend(order(1 "Actual core CPI inflation" 2 "Phillips Curve fitted (with supply shocks)")) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Actual vs fitted inflation: augmented Phillips Curve")

* Graph 4 – redrawing of the graph 2 *
twoway (line fe_pc_supply tq if e(sample)), xline(`=q2020q1') yline(0) ytitle("Forecast error: actual - fitted (pp)") xtitle("Quarter") title("Phillips Curve forecasting errors (with supply shocks)")



* Stage 4 – Proving that PC has not changed structurally *


* Mehtod 1

gen ugap_post = u_gap*post2020
reg pi_core_cpi spf_inflation_1year u_gap ugap_post if !missing(spf_inflation_1year), robust
test ugap_post
lincom u_gap + ugap_post

* Method 2

reg pi_core_cpi spf_inflation_1year u_gap if tq<yq(2020,1) & !missing(spf_inflation_1year), robust
est store pre
reg pi_core_cpi spf_inflation_1year u_gap if tq>=yq(2020,1) & !missing(spf_inflation_1year), robust
est store post

* Method 3

preserve
keep if !missing(spf_inflation_1year)
rolling b_ugap=_b[u_gap] se_ugap=_se[u_gap], window(40) saving(roll_spf, replace): reg pi_core_cpi spf_inflation_1year u_gap
use roll_spf, clear
gen tq = _n
tsset tq
gen ub = b_ugap + 1.96*se_ugap
gen lb = b_ugap - 1.96*se_ugap
twoway (line b_ugap tq) (line ub tq) (line lb tq), yline(0) title("Rolling Phillips Curve slope (u_gap)") xtitle("Rolling window index") ytitle("Slope estimate")
restore



