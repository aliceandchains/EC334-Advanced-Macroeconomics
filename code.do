* Stage 1 - understanding the data and preparing for analysis *

clear all
set more off

* Importing the data from the file *
import delimited "alldata.csv", clear varnames(1)

* Understanding missingness across variables *
misstable summarize spf_inflation_1year michigan_1y_median swap_1year swap_5year swap_10year cpilfesl pceexcludingfoodandenergy pceall unrate nrou v_u poilbreusdm import_def

* Constructing the time variable and tsset *
gen tq = yq(year, quarter)
format tq %tq
tsset tq

* Constructing inflation measures – Core CPI and Core PCE inflation *
gen pi_core_cpi = 400*(ln(cpilfesl) - ln(L.cpilfesl))
gen pi_core_pce = 400*(ln(pceexcludingfoodandenergy) - ln(L.pceexcludingfoodandenergy))

* Constructing Slack measures – unemployment gap *
gen u_gap = unrate - nrou

* Constructing Supply shocks variables – oil price inflation and import price inflation *
gen d_oil = 400*(ln(poilbreusdm) - ln(L.poilbreusdm))
gen pi_import = 400*(ln(import_def) - ln(L.import_def))

* Constructing a dummy variable for Pandemic period *
gen post2020 = (tq >= yq(2020,1))

* Create a numeric quarter cutoff for graphs (avoids xline(yq()) errors) *
gen tq_2020q1 = yq(2020,1)

* Summarising new variables *
summarize pi_core_cpi pi_core_pce u_gap spf_inflation_1year michigan_1y_median v_u d_oil pi_import




* Stage 2 - Baseline Phillips Curve regression (lecture-consistent OLS) *
reg pi_core_cpi spf_inflation_1year u_gap if !missing(spf_inflation_1year)



* Robust SE version (optional robustness) *
reg pi_core_cpi spf_inflation_1year u_gap if !missing(spf_inflation_1year), robust

* Generate fitted values and forecasting errors using the baseline regression sample *
reg pi_core_cpi spf_inflation_1year u_gap if !missing(spf_inflation_1year)
predict pi_hat_pc if e(sample)
gen fe_pc = pi_core_cpi - pi_hat_pc if e(sample)

* Graph 1: Actual vs fitted inflation, with vertical line at 2020Q1 *
twoway (line pi_core_cpi tq if e(sample)) (line pi_hat_pc tq if e(sample)), xline(`=yq(2020,1)') legend(order(1 "Actual core CPI inflation" 2 "Phillips Curve fitted")) ytitle("Annualised quarterly inflation (pp)") xtitle("Quarter") title("Actual vs fitted inflation: baseline Phillips Curve")

* Optional: save Graph 1 as PNG for dissertation *
graph export "fig_actual_vs_fitted_pc.png", replace

* Graph 2: Forecast errors (actual - fitted), showing systematic post-2020 errors *
twoway (line fe_pc tq if e(sample)), xline(`=yq(2020,1)') yline(0) ytitle("Forecast error: actual - fitted (pp)") xtitle("Quarter") title("Phillips Curve forecasting errors")

* Optional: save Graph 2 as PNG for dissertation *
graph export "fig_forecast_errors_pc.png", replace



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



