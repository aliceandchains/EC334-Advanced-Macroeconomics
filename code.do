*******************************************************
* Load data + preprocess once
*******************************************************

clear all
set more off

import delimited "Final Data.csv", clear varnames(1)

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







*******************************************************
* (a) PCE × Michigan + oil
*******************************************************
reg pi_pce michigan_1y_median u_gap d_oil if train_pc & ok1, robust
estimates store m1_pce_mich_oil
local b1 : display %5.2f _b[michigan_1y_median]
local g1 : display %5.2f _b[u_gap]

*******************************************************
* (b) PCE × SPF + oil
*******************************************************
reg pi_pce spf_inflation_1year u_gap d_oil if train_pc & ok2, robust
estimates store m2_pce_spf_oil
local b2 : display %5.2f _b[spf_inflation_1year]
local g2 : display %5.2f _b[u_gap]

*******************************************************
* (c) CPI × Michigan + oil
*******************************************************
reg pi_cpi michigan_1y_median u_gap d_oil if train_pc & ok3, robust
estimates store m3_cpi_mich_oil
local b3 : display %5.2f _b[michigan_1y_median]
local g3 : display %5.2f _b[u_gap]

*******************************************************
* (d) CPI × SPF + oil
*******************************************************
reg pi_cpi spf_inflation_1year u_gap d_oil if train_pc & ok4, robust
estimates store m4_cpi_spf_oil
local b4 : display %5.2f _b[spf_inflation_1year]
local g4 : display %5.2f _b[u_gap]









*******************************************************
* FEVD (VAR): inflation, expectations, slack, supply
* Output: 3 separate FEVD PNGs
* MAX compatibility version: no ///, no plotopts(), no saving()
*******************************************************

clear all
set more off

* Global look
set scheme s2color
graph set window fontface "Helvetica"

local tsize medsmall
local lsize small
local nsize vsmall
local WPNG  3200

import delimited "Final Data.csv", clear varnames(1)

cap drop tq
gen tq = yq(year, quarter)
format tq %tq
sort tq
tsset tq

cap drop pi_pce u_gap d_oil
gen pi_pce = 400*(ln(pceall) - ln(L.pceall))
gen u_gap  = unrate - nrou
gen d_oil  = 400*(ln(poilbreusdm) - ln(L.poilbreusdm))

label var pi_pce "Headline PCE inflation (annualised q/q)"
label var u_gap  "Unemployment gap (u-u*)"
label var d_oil  "Oil price inflation (annualised q/q)"

local expvar michigan_1y_median
label var `expvar' "1y-ahead inflation expectations (Michigan median)"

scalar tq1984q1 = yq(1984,1)
scalar tq2023q4 = yq(2023,4)

cap drop ok_var
gen ok_var = (tq>=tq1984q1 & tq<=tq2023q4) & !missing(pi_pce, `expvar', u_gap, d_oil)

* VAR
var pi_pce `expvar' u_gap d_oil if ok_var, lags(1/4)
varstable

* IRF horizon set here
irf set irf_pc_var, replace
irf create base, step(16) replace order(d_oil `expvar' u_gap pi_pce)

*******************************************************
* CLEAN FEVD PLOTS (manual twoway from irf table output)
* No ///, no plotopts(), no saving()
*******************************************************

* Aesthetics you can reuse
set scheme s2color
graph set window fontface "Helvetica"
local tsize medsmall
local lsize small
local nsize vsmall
local WPNG  3200

local TWBASE graphregion(color(white)) plotregion(color(white)) legend(off) ylabel(0(.02).10, angle(horizontal) labsize(`lsize') nogrid) xlabel(0(4)16, labsize(`lsize')) xtitle("Horizon (quarters)", size(`lsize')) ytitle("Share of FE variance", size(`lsize')) title(, size(`tsize')) note(, size(`nsize')) graphregion(margin(6 6 6 6)) plotregion(margin(4 4 2 2))

* Helper: build dataset from matrix returned by irf table
cap program drop fevd_make
program define fevd_make
    args impulsevar respvar outname col
    preserve
    tempname M
    quietly irf table fevd, impulse(`impulsevar') response(`respvar')
    matrix `M' = r(table)
    clear
    svmat double `M', names(col)
    gen h = _n-1
    * The FEVD series is typically in the first column of r(table) for old irf table fevd.
    * We robustify by taking the first numeric column produced by svmat.
    ds, has(type numeric)
    local numvars `r(varlist)'
    local first : word 1 of `numvars'
    gen share = `first'
    keep h share
    order h share
    twoway line share h, lcolor(`col') lwidth(medthick) `TWBASE' title("FEVD of headline PCE inflation: `outname'", size(`tsize')) note("VAR(4), sample 1984Q1–2023Q4. Cholesky order: d_oil, exp, u_gap, pi_pce.", size(`nsize'))
    graph export "fig_fevd_`outname'_pi_pce.png", replace width(`WPNG')
    restore
end

* Make the three clean plots
fevd_make d_oil  pi_pce supply_oil  blue
fevd_make `expvar' pi_pce expectations red
fevd_make u_gap  pi_pce ugap        black










*******************************************************
* Stock–Watson accelerationist PC (textbook-style, clean)
* Δπ_t = α + γ (u_t - u*_t) + v_t
* Subsamples:
*   1968Q1–1983Q4
*   1984Q1–1999Q4
*   2000Q1–2019Q4   (cutoff fixed)
* Plot window:
*   x in [-2.5, 5]
*   y in [-3, 6]
* IMPORTANT: estimation uses all available points in each subsample
*            plotting trims to the window so the figure looks clean
*******************************************************

tsset tq

* --- Cutoffs ---
scalar tq1968q1 = yq(1968,1)
scalar tq1983q4 = yq(1983,4)
scalar tq1984q1 = yq(1984,1)
scalar tq1999q4 = yq(1999,4)
scalar tq2000q1 = yq(2000,1)
scalar tq2019q4 = yq(2019,4)

* --- Accelerationist inflation ---
cap drop d_pi_pce ok_sw s1 s2 s3
gen d_pi_pce = D.pi_pce
label var d_pi_pce "Change in inflation (pp)"

gen ok_sw = !missing(d_pi_pce, u_gap) & inrange(tq, tq1968q1, tq2019q4)

gen byte s1 = inrange(tq, tq1968q1, tq1983q4)
gen byte s2 = inrange(tq, tq1984q1, tq1999q4)
gen byte s3 = inrange(tq, tq2000q1, tq2019q4)

* --- Plot window (YOU asked for this) ---
local xmin = -2.5
local xmax = 5
local ymin = -3
local ymax = 6

* Plot-only sample: trim extreme points so the figure isn’t distorted
cap drop ok_plot
gen ok_plot = ok_sw & inrange(u_gap, `xmin', `xmax') & inrange(d_pi_pce, `ymin', `ymax')

* --- Estimate α and γ on FULL subsample (not trimmed) ---
local a1 = .
local g1 = .
local a2 = .
local g2 = .
local a3 = .
local g3 = .

quietly count if ok_sw & s1
if r(N) > 10 {
    reg d_pi_pce u_gap if ok_sw & s1, robust
    local a1 = _b[_cons]
    local g1 = _b[u_gap]
}

quietly count if ok_sw & s2
if r(N) > 10 {
    reg d_pi_pce u_gap if ok_sw & s2, robust
    local a2 = _b[_cons]
    local g2 = _b[u_gap]
}

quietly count if ok_sw & s3
if r(N) > 10 {
    reg d_pi_pce u_gap if ok_sw & s3, robust
    local a3 = _b[_cons]
    local g3 = _b[u_gap]
}

* --- Gamma labels (build as single strings) ---
local lab1 "γ = " + string(`g1', "%6.3f")
local lab2 "γ = " + string(`g2', "%6.3f")
local lab3 "γ = " + string(`g3', "%6.3f")

* --- Place gamma labels near top-left but inside plot window ---
local xtext = `xmin' + 0.25*(`xmax' - `xmin')
local y1 = `a1' + `g1'*`xtext'
local y2 = `a2' + `g2'*`xtext'
local y3 = `a3' + `g3'*`xtext'

* Nudge so labels don’t overlap
local dy = 0.35

* --- Plot: textbook markers + long lines + clean window ---
twoway (scatter d_pi_pce u_gap if ok_plot & s1, mcolor(gs12) msymbol(Oh) msize(medsmall)) (scatter d_pi_pce u_gap if ok_plot & s2, mcolor(gs12) msymbol(Sh) msize(medsmall)) (scatter d_pi_pce u_gap if ok_plot & s3, mcolor(gs12) msymbol(Dh) msize(medsmall)) (function y = `a1' + `g1'*x, range(`xmin' `xmax') lcolor(cranberry) lwidth(medthick)) (function y = `a2' + `g2'*x, range(`xmin' `xmax') lcolor(midblue) lwidth(medthick)) (function y = `a3' + `g3'*x, range(`xmin' `xmax') lcolor(black) lwidth(medthick)), xscale(range(`xmin' `xmax')) yscale(range(`ymin' `ymax')) xlabel(-2(1)5, labsize(small) labcolor(black)) ylabel(-3(1)6, labsize(small) labcolor(black) glcolor(gs16) glwidth(vthin)) xtitle("Unemployment gap (u - u*)", size(small) color(black)) ytitle("Change in inflation (pp)", size(small) color(black)) legend(order(4 "Fit 1968Q1–1983Q4" 5 "Fit 1984Q1–1999Q4" 6 "Fit 2000Q1–2019Q4") cols(1) position(6) ring(1) size(small) color(black) region(lstyle(none))) text(`=`y1'+`dy'' `xtext' "`lab1'", size(small) color(cranberry) placement(e)) text(`y2' `xtext' "`lab2'", size(small) color(midblue) placement(e)) text(`=`y3'-`dy'' `xtext' "`lab3'", size(small) color(black) placement(e)) graphregion(color(white)) plotregion(color(white)) name(sw_pc_clean_3subs, replace)

graph export "fig_sw_pc_clean_3subs_2019q4.pdf", replace
graph export "fig_sw_pc_clean_3subs_2019q4.png", replace width(3200)







*******************************************************
* Counterfactual analysis (Beaudry-style accounting)
* Baseline PC: pi_t = a + b*Exp_t + g*u_gap_t + d*d_oil_t + e_t
* Train: 1984Q1–2020Q2 ; Counterfactual paths plotted through 2023Q4
* Counterfactuals: hold each regressor at its pre-2020 mean
*******************************************************

clear all
set more off

import delimited "Final Data.csv", clear varnames(1)

cap drop tq
gen tq = yq(year, quarter)
format tq %tq
sort tq
tsset tq

*******************************************************
* Build variables
*******************************************************

cap drop pi_pce u_gap d_oil
gen pi_pce = 400*(ln(pceall) - ln(L.pceall))
label var pi_pce "Headline PCE inflation (annualised q/q)"

gen u_gap = unrate - nrou
label var u_gap "Unemployment gap"

gen d_oil = 400*(ln(poilbreusdm) - ln(L.poilbreusdm))
label var d_oil "Oil price inflation"

*******************************************************
* Sample cutoffs
*******************************************************

scalar tq1984q1 = yq(1984,1)
scalar tq2020q2 = yq(2020,2)
scalar tq2020q1 = yq(2020,1)
scalar tq2023q4 = yq(2023,4)

cap drop train_pc plot_win ok_pc
gen train_pc = (tq >= tq1984q1 & tq <= tq2020q2)
gen plot_win = (tq >= tq1984q1 & tq <= tq2023q4)
gen ok_pc = !missing(pi_pce, michigan_1y_median, u_gap, d_oil) & plot_win

*******************************************************
* Pre-2020 means
*******************************************************

quietly summarize u_gap if ok_pc & tq < tq2020q1, meanonly
scalar ugap_bar = r(mean)

quietly summarize michigan_1y_median if ok_pc & tq < tq2020q1, meanonly
scalar exp_bar = r(mean)

quietly summarize d_oil if ok_pc & tq < tq2020q1, meanonly
scalar oil_bar = r(mean)

*******************************************************
* Baseline Phillips Curve
*******************************************************

reg pi_pce michigan_1y_median u_gap d_oil if ok_pc & train_pc, robust

scalar b0 = _b[_cons]
scalar bE = _b[michigan_1y_median]
scalar bU = _b[u_gap]
scalar bO = _b[d_oil]

*******************************************************
* Fitted and counterfactual paths
*******************************************************

cap drop pi_hat_full pi_hat_gapfix pi_hat_expfixed pi_hat_oilfixed
gen pi_hat_full = b0 + bE*michigan_1y_median + bU*u_gap + bO*d_oil if ok_pc
gen pi_hat_gapfix = b0 + bE*michigan_1y_median + bU*ugap_bar + bO*d_oil if ok_pc
gen pi_hat_expfixed = b0 + bE*exp_bar + bU*u_gap + bO*d_oil if ok_pc
gen pi_hat_oilfixed = b0 + bE*michigan_1y_median + bU*u_gap + bO*oil_bar if ok_pc

*******************************************************
* Common plot settings
*******************************************************

local yttl "Annualised quarterly inflation (pp)"
local xlbl `=tq1984q1'(40)`=tq2023q4'
local xfmt %tqCCYY!Qq

*******************************************************
* (1) Gap fixed at pre-2020 mean
*******************************************************

twoway (line pi_pce tq if ok_pc, lcolor(midblue) lwidth(medthick)) (line pi_hat_full tq if ok_pc, lcolor(cranberry) lwidth(medthick)) (line pi_hat_gapfix tq if ok_pc, lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q1', lpattern(dash) lcolor(gs10)) ytitle("`yttl'") xtitle("") xlabel(`xlbl', format(`xfmt')) legend(order(1 "Actual inflation" 2 "Predicted (full PC)" 3 "Gap fixed") cols(1)) graphregion(color(white)) plotregion(color(white))

graph export "fig_cf_gap_fixed.png", replace width(3200)

*******************************************************
* (2) Expectations fixed at pre-2020 mean
*******************************************************

twoway (line pi_pce tq if ok_pc, lcolor(midblue) lwidth(medthick)) (line pi_hat_full tq if ok_pc, lcolor(cranberry) lwidth(medthick)) (line pi_hat_expfixed tq if ok_pc, lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q1', lpattern(dash) lcolor(gs10)) ytitle("") xtitle("") xlabel(`xlbl', format(`xfmt')) legend(order(1 "Actual inflation" 2 "Predicted (full PC)" 3 "Expectations fixed") cols(1)) graphregion(color(white)) plotregion(color(white))

graph export "fig_cf_exp_fixed.png", replace width(3200)

*******************************************************
* (3) Oil fixed at pre-2020 mean
*******************************************************

twoway (line pi_pce tq if ok_pc, lcolor(midblue) lwidth(medthick)) (line pi_hat_full tq if ok_pc, lcolor(cranberry) lwidth(medthick)) (line pi_hat_oilfixed tq if ok_pc, lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q1', lpattern(dash) lcolor(gs10)) ytitle("`yttl'") xtitle("") xlabel(`xlbl', format(`xfmt')) legend(order(1 "Actual inflation" 2 "Predicted (full PC)" 3 "Oil fixed") cols(1)) graphregion(color(white)) plotregion(color(white))

graph export "fig_cf_oil_fixed.png", replace width(3200)










*******************************************************
* Fiscal impulse: RDPI gap (real DPI per capita vs pre-2020 trend)
* INPUT: rdpi_pc_real = BEA Table 2.1 line 39 (per capita, chained 2017$)
*******************************************************

* Replace rdpi_pc_real with your actual column name
cap drop l_rdpi t_index rdpi_trend rdpi_gap
gen l_rdpi = ln(rdpi_pc_real)
label var l_rdpi "Log real DPI per capita"

gen t_index = _n

scalar tq2019q4 = yq(2019,4)

quietly reg l_rdpi t_index if tq>=tq1984q1 & tq<=tq2019q4, robust
predict rdpi_trend, xb
gen rdpi_gap = l_rdpi - rdpi_trend
label var rdpi_gap "Fiscal impulse: log real DPIpc minus pre-2020 trend"

* Sanity plot: should spike in 2020-21 and mean-revert
twoway (line rdpi_gap tq if tq>=tq1984q1, lwidth(medthick) lcolor(black)), xline(`=yq(2020,1)', lpattern(dash) lcolor(gs10)) ytitle("Log deviation from pre-2020 trend") xtitle("") graphregion(color(white)) plotregion(color(white)) name(fig_fiscal_gap, replace)

graph export "fig_fiscal_gap.pdf", replace






