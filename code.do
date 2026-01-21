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















