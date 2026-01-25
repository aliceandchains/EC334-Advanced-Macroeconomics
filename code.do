```stata
*******************************************************
* EC334 Dissertation — Production Master Do-file (NO ///)
* One import, no redundancies, final section = fiscal
*******************************************************

clear all
set more off

* Global style
set scheme s2color
graph set window fontface "Helvetica"
local WPNG 3200
local tsize medsmall
local lsize small
local nsize vsmall

*******************************************************
* 0) Import once + time index
*******************************************************

import delimited "Final Data.csv", clear varnames(1)

cap drop tq
gen tq = yq(year, quarter)
format tq %tq
sort tq
tsset tq

*******************************************************
* 1) Cutoffs used throughout
*******************************************************

scalar tq1968q1 = yq(1968,1)
scalar tq1984q1 = yq(1984,1)
scalar tq2019q4 = yq(2019,4)
scalar tq2020q1 = yq(2020,1)
scalar tq2020q2 = yq(2020,2)
scalar tq2020q3 = yq(2020,3)
scalar tq2023q4 = yq(2023,4)
scalar tq2024q4 = yq(2024,4)

*******************************************************
* 2) Construct variables (once)
*******************************************************

* Inflation (annualised q/q)
cap drop pi_core_cpi pi_core_pce pi_pce pi_cpi
gen pi_core_cpi = 400*(ln(cpilfesl) - ln(L.cpilfesl))
gen pi_core_pce = 400*(ln(pceexcludingfoodandenergy) - ln(L.pceexcludingfoodandenergy))
gen pi_pce      = 400*(ln(pceall) - ln(L.pceall))
gen pi_cpi      = 400*(ln(cpilfesl) - ln(L.cpilfesl))

label var pi_core_cpi "Core CPI inflation (annualised q/q, pp)"
label var pi_core_pce "Core PCE inflation (annualised q/q, pp)"
label var pi_pce      "Headline PCE inflation (annualised q/q, pp)"
label var pi_cpi      "CPI inflation (annualised q/q, pp)"

* Slack / tightness
cap drop u_gap l_vu
gen u_gap = unrate - nrou
gen l_vu  = ln(v_u)
label var u_gap "Unemployment gap (u - u*)"
label var l_vu  "Log tightness ln(v/u)"

* Supply / external prices
cap drop d_oil pi_import
gen d_oil     = 400*(ln(poilbreusdm) - ln(L.poilbreusdm))
gen pi_import = 400*(ln(import_def) - ln(L.import_def))
label var d_oil     "Oil price inflation (annualised q/q, pp)"
label var pi_import "Import price inflation (annualised q/q, pp)"

* Train/test flags
cap drop train_pc test_pc
gen train_pc = (tq >= tq1984q1 & tq <= tq2020q2)
gen test_pc  = (tq >= tq2020q3)

* Plot x-axis labels
local xlbl = tq1984q1
local xlbl2 = tq2024q4

*******************************************************
* 3) Baseline expectations-augmented Phillips Curve (OOS)
*    Dep var: headline PCE inflation (pi_pce)
*******************************************************

cap drop ok_pc
gen ok_pc = !missing(pi_pce, michigan_1y_median, u_gap) & tq>=tq1984q1

reg pi_pce michigan_1y_median u_gap if train_pc & ok_pc, robust
est store pc_base

cap drop pi_hat_pc pi_fit_pc pi_fc_pc
predict pi_hat_pc if ok_pc, xb
gen pi_fit_pc = pi_hat_pc if train_pc & ok_pc
gen pi_fc_pc  = pi_hat_pc if test_pc  & ok_pc

cap drop fe_pc se_pc ae_pc
gen fe_pc = pi_pce - pi_hat_pc if test_pc & ok_pc
gen se_pc = fe_pc^2 if test_pc & !missing(fe_pc)
gen ae_pc = abs(fe_pc) if test_pc & !missing(fe_pc)

quietly summarize se_pc
display "Baseline PC RMSE (>=2020Q3) = " sqrt(r(mean))
quietly summarize ae_pc
display "Baseline PC MAE  (>=2020Q3) = " r(mean)

twoway (line pi_pce tq if ok_pc & tq>=tq1984q1, lcolor(midblue) lwidth(medthick)) (line pi_fit_pc tq if !missing(pi_fit_pc) & tq>=tq1984q1, lcolor(cranberry) lwidth(medthick)) (line pi_fc_pc tq if !missing(pi_fc_pc) & tq>=tq1984q1, lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(gs10) lwidth(thin)) ytitle("Annualised quarterly inflation (pp)", size(`lsize')) xtitle("") xlabel(`=tq1984q1'(40)`=tq2024q4', format(%tqCCYY!Qq) labsize(`lsize')) ylabel(, labsize(`lsize') glcolor(gs16) glwidth(vthin)) legend(order(1 "Actual headline PCE inflation" 2 "Fit (1984Q1–2020Q2)" 3 "Forecast (>=2020Q3)") cols(1) position(6) ring(1) size(`lsize') region(lstyle(none))) graphregion(color(white)) plotregion(color(white)) name(fig_pc_base_oos, replace)

graph export "fig_baseline_pc_oos_forecast.pdf", replace

*******************************************************
* 4) Robustness: 2x2 PCs (inflation x expectations)
*******************************************************

cap program drop _coefpos
program define _coefpos, rclass
    syntax varname [if]
    quietly summarize `varlist' `if', meanonly
    return scalar ytxt = r(max) - 0.08*(r(max)-r(min))
end

local xlbl_list `=tq1984q1'(40)`=tq2024q4'
local xfmt %tqCCYY!Qq
scalar tq_text = yq(1996,1)

* (a) PCE × Michigan
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
twoway (line pi_pce tq if ok1, lcolor(midblue) lwidth(medthick)) (line pi_fit1 tq if !missing(pi_fit1), lcolor(cranberry) lwidth(medthick)) (line pi_fc1 tq if !missing(pi_fc1), lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(black) lwidth(thin)) title("(a)", size(`lsize')) text(`=ytxt1' `=tq_text' "{it:β}=`b1'   {it:γ}=`g1'", size(`lsize') placement(ne)) ytitle("Annualised quarterly inflation (pp)", size(`lsize')) xtitle("") xlabel(`xlbl_list', format(`xfmt') labsize(`lsize')) ylabel(, labsize(`lsize') glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(Ga, replace)

* (b) PCE × SPF
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
twoway (line pi_pce tq if ok2, lcolor(midblue) lwidth(medthick)) (line pi_fit2 tq if !missing(pi_fit2), lcolor(cranberry) lwidth(medthick)) (line pi_fc2 tq if !missing(pi_fc2), lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(black) lwidth(thin)) title("(b)", size(`lsize')) text(`=ytxt2' `=tq_text' "{it:β}=`b2'   {it:γ}=`g2'", size(`lsize') placement(ne)) ytitle("") xtitle("") xlabel(`xlbl_list', format(`xfmt') labsize(`lsize')) ylabel(, labsize(`lsize') glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(Gb, replace)

* (c) CPI × Michigan
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
twoway (line pi_cpi tq if ok3, lcolor(midblue) lwidth(medthick)) (line pi_fit3 tq if !missing(pi_fit3), lcolor(cranberry) lwidth(medthick)) (line pi_fc3 tq if !missing(pi_fc3), lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(black) lwidth(thin)) title("(c)", size(`lsize')) text(`=ytxt3' `=tq_text' "{it:β}=`b3'   {it:γ}=`g3'", size(`lsize') placement(ne)) ytitle("Annualised quarterly inflation (pp)", size(`lsize')) xtitle("") xlabel(`xlbl_list', format(`xfmt') labsize(`lsize')) ylabel(, labsize(`lsize') glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(Gc, replace)

* (d) CPI × SPF
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
twoway (line pi_cpi tq if ok4, lcolor(midblue) lwidth(medthick)) (line pi_fit4 tq if !missing(pi_fit4), lcolor(cranberry) lwidth(medthick)) (line pi_fc4 tq if !missing(pi_fc4), lcolor(cranberry) lpattern(dash) lwidth(medthick)), xline(`=tq2020q3', lpattern(dash) lcolor(black) lwidth(thin)) title("(d)", size(`lsize')) text(`=ytxt4' `=tq_text' "{it:β}=`b4'   {it:γ}=`g4'", size(`lsize') placement(ne)) ytitle("") xtitle("") xlabel(`xlbl_list', format(`xfmt') labsize(`lsize')) ylabel(, labsize(`lsize') glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(Gd, replace)

graph combine Ga Gb Gc Gd, cols(2) imargin(2 2 2 2) graphregion(color(white)) name(fig_pc_robust_2x2, replace)
graph export "fig_baseline_pc_robust_final.pdf", replace
graph export "fig_baseline_pc_robust_final.png", replace width(`WPNG')

cap which esttab
if _rc ssc install estout
esttab m1_pce_mich m2_pce_spf m3_cpi_mich m4_cpi_spf using "tab_baseline_pc_robust.tex", replace se b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) label stats(N r2, labels("Observations" "R-squared")) compress

*******************************************************
* 5) FEVD (VAR): inflation, expectations, slack, supply
*******************************************************

cap drop ok_var
gen ok_var = (tq>=tq1984q1 & tq<=tq2023q4) & !missing(pi_pce, michigan_1y_median, u_gap, d_oil)

var pi_pce michigan_1y_median u_gap d_oil if ok_var, lags(1/4)
irf set irf_pc_var, replace
irf create base, step(16) replace order(d_oil michigan_1y_median u_gap pi_pce)

local TWBASE graphregion(color(white)) plotregion(color(white)) legend(off) ylabel(0(.02).10, angle(horizontal) labsize(`lsize') nogrid) xlabel(0(4)16, labsize(`lsize')) xtitle("Horizon (quarters)", size(`lsize')) ytitle("Share of FE variance", size(`lsize')) note("VAR(4), sample 1984Q1–2023Q4. Cholesky order: oil, exp, u_gap, inflation.", size(`nsize'))

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
    ds, has(type numeric)
    local numvars `r(varlist)'
    local first : word 1 of `numvars'
    gen share = `first'
    keep h share
    twoway line share h, lcolor(`col') lwidth(medthick) `TWBASE' title("FEVD of headline PCE inflation: `outname'", size(`tsize'))
    graph export "fig_fevd_`outname'_pi_pce.png", replace width(`WPNG')
    restore
end

fevd_make d_oil pi_pce supply_oil blue
fevd_make michigan_1y_median pi_pce expectations red
fevd_make u_gap pi_pce ugap black

*******************************************************
* 6) Stock–Watson accelerationist PC (three subsamples)
*******************************************************

cap drop d_pi_pce ok_sw s1 s2 s3 ok_plot
gen d_pi_pce = D.pi_pce
gen ok_sw = !missing(d_pi_pce, u_gap) & inrange(tq, tq1968q1, tq2019q4)
gen byte s1 = inrange(tq, yq(1968,1), yq(1983,4))
gen byte s2 = inrange(tq, yq(1984,1), yq(1999,4))
gen byte s3 = inrange(tq, yq(2000,1), yq(2019,4))

local xmin = -2.5
local xmax = 5
local ymin = -3
local ymax = 6

gen ok_plot = ok_sw & inrange(u_gap, `xmin', `xmax') & inrange(d_pi_pce, `ymin', `ymax')

quietly reg d_pi_pce u_gap if ok_sw & s1, robust
local a1 = _b[_cons]
local g1 = _b[u_gap]
quietly reg d_pi_pce u_gap if ok_sw & s2, robust
local a2 = _b[_cons]
local g2 = _b[u_gap]
quietly reg d_pi_pce u_gap if ok_sw & s3, robust
local a3 = _b[_cons]
local g3 = _b[u_gap]

local xtext = `xmin' + 0.25*(`xmax' - `xmin')
local y1 = `a1' + `g1'*`xtext'
local y2 = `a2' + `g2'*`xtext'
local y3 = `a3' + `g3'*`xtext'

twoway (scatter d_pi_pce u_gap if ok_plot & s1, mcolor(gs12) msymbol(Oh) msize(medsmall)) (scatter d_pi_pce u_gap if ok_plot & s2, mcolor(gs12) msymbol(Sh) msize(medsmall)) (scatter d_pi_pce u_gap if ok_plot & s3, mcolor(gs12) msymbol(Dh) msize(medsmall)) (function y = `a1' + `g1'*x, range(`xmin' `xmax') lcolor(cranberry) lwidth(medthick)) (function y = `a2' + `g2'*x, range(`xmin' `xmax') lcolor(midblue) lwidth(medthick)) (function y = `a3' + `g3'*x, range(`xmin' `xmax') lcolor(black) lwidth(medthick)), xscale(range(`xmin' `xmax')) yscale(range(`ymin' `ymax')) xlabel(-2(1)5, labsize(`lsize')) ylabel(-3(1)6, labsize(`lsize') glcolor(gs16) glwidth(vthin)) xtitle("Unemployment gap (u - u*)", size(`lsize')) ytitle("Change in inflation (pp)", size(`lsize')) legend(order(4 "Fit 1968Q1–1983Q4" 5 "Fit 1984Q1–1999Q4" 6 "Fit 2000Q1–2019Q4") cols(1) position(6) ring(1) size(`lsize') region(lstyle(none))) text(`y1' `xtext' "γ = " + string(`g1', "%6.3f"), size(`lsize') color(cranberry) placement(e)) text(`y2' `xtext' "γ = " + string(`g2', "%6.3f"), size(`lsize') color(midblue) placement(e)) text(`y3' `xtext' "γ = " + string(`g3', "%6.3f"), size(`lsize') color(black) placement(e)) graphregion(color(white)) plotregion(color(white)) name(fig_sw_pc, replace)

graph export "fig_sw_pc_clean_3subs_2019q4.pdf", replace
graph export "fig_sw_pc_clean_3subs_2019q4.png", replace width(`WPNG')



*******************************************************
* 7) Fiscal mechanism section (2001Q1–2023Q4)
* Outputs: RDPI gap, ECI wage inflation, unemployment gap plots
* Then regressions: RDPI gap -> core PCE inflation with controls
*******************************************************

* ---------- Style (match PC figures) ----------
local lsize small
local WPNG 3200

* ---------- Ensure time series setup ----------
cap confirm variable tq
if _rc {
    di as error "Variable tq not found. Create tq first (yq(year,quarter)) and tsset tq."
    exit 198
}
tsset tq

* ---------- Time cutoffs ----------
scalar tq2001q1 = yq(2001,1)
scalar tq2020q2 = yq(2020,2)
scalar tq2023q4 = yq(2023,4)

* ---------- Ensure a linear time index exists ----------
cap confirm variable t_index
if _rc {
    gen t_index = _n
}

*******************************************************
* 7.1 RDPI gap: log(dpipercapita) minus pre-2020Q2 trend
*******************************************************

cap confirm variable dpipercapita
if _rc {
    di as error "Variable dpipercapita not found. Check your data."
    exit 198
}

cap drop l_rdpi rdpi_trend rdpi_gap
gen l_rdpi = ln(dpipercapita)
label var l_rdpi "Log real DPI per capita"

quietly reg l_rdpi t_index if inrange(tq, tq2001q1, tq2020q2) & !missing(l_rdpi), vce(robust)
predict rdpi_trend if inrange(tq, tq2001q1, tq2023q4) & !missing(l_rdpi), xb
label var rdpi_trend "Pre-2020Q2 trend"

gen rdpi_gap = l_rdpi - rdpi_trend if inrange(tq, tq2001q1, tq2023q4) & !missing(l_rdpi, rdpi_trend)
label var rdpi_gap "RDPI gap: log deviation from pre-2020Q2 trend"

twoway (line rdpi_gap tq if inrange(tq, tq2001q1, tq2023q4) & !missing(rdpi_gap), lcolor(midblue) lwidth(medthick)), xline(`=tq2020q2', lpattern(dash) lcolor(black) lwidth(thin)) ytitle("Log deviation from pre-2020Q2 trend", size(`lsize')) xtitle("") xlabel(`=tq2001q1'(40)`=tq2023q4', format(%tqCCYY!Qq) labsize(`lsize')) ylabel(, labsize(`lsize') glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(fig_fiscal_gap, replace)

graph export "fig_fiscal_gap.pdf", replace
graph export "fig_fiscal_gap.png", replace width(`WPNG')

*******************************************************
* 7.2 Wage inflation proxy: ECI (annualised q/q, pp)
*******************************************************

cap confirm variable eci
if _rc {
    di as error "Variable eci not found. Skipping ECI plot + wage regressions."
}

cap drop pi_eci ok_eci
gen pi_eci = 400*(ln(eci) - ln(L.eci))
label var pi_eci "ECI wage inflation (annualised q/q, pp)"

gen ok_eci = inrange(tq, tq2001q1, tq2023q4) & !missing(pi_eci)

twoway (line pi_eci tq if ok_eci, lcolor(midblue) lwidth(medthick)), xline(`=tq2020q2', lpattern(dash) lcolor(black) lwidth(thin)) ytitle("Annualised quarterly ECI inflation (pp)", size(`lsize')) xtitle("") xlabel(`=tq2001q1'(40)`=tq2023q4', format(%tqCCYY!Qq) labsize(`lsize')) ylabel(, labsize(`lsize') glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(fig_eci_infl, replace)

graph export "fig_eci_infl.pdf", replace
graph export "fig_eci_infl.png", replace width(`WPNG')

*******************************************************
* 7.3 Unemployment gap (u - u*)
*******************************************************

cap confirm variable u_gap
if _rc {
    di as error "Variable u_gap not found. Skipping unemployment gap plot."
}
cap drop ok_ugap
gen ok_ugap = inrange(tq, tq2001q1, tq2023q4) & !missing(u_gap)

twoway (line u_gap tq if ok_ugap, lcolor(midblue) lwidth(medthick)), xline(`=tq2020q2', lpattern(dash) lcolor(black) lwidth(thin)) ytitle("Unemployment gap (u - u*)", size(`lsize')) xtitle("") xlabel(`=tq2001q1'(40)`=tq2023q4', format(%tqCCYY!Qq) labsize(`lsize')) ylabel(, labsize(`lsize') glcolor(gs16) glwidth(vthin)) legend(off) graphregion(color(white)) plotregion(color(white)) name(fig_ugap, replace)

graph export "fig_ugap.pdf", replace
graph export "fig_ugap.png", replace width(`WPNG')

*******************************************************
* 7.4 Fiscal impulse -> inflation regressions
* Dep var: pi_core_pce (annualised q/q, pp)
* Key regressor: L1.rdpi_gap (log points)
* Controls: L1.u_gap, L1.michigan_1y_median, d_oil, pi_import, L1.pi_core_pce
*******************************************************

cap confirm variable pi_core_pce
if _rc {
    di as error "Variable pi_core_pce not found. Cannot run fiscal->inflation regressions."
    exit 198
}
cap confirm variable michigan_1y_median
if _rc {
    di as error "Variable michigan_1y_median not found. Cannot run fiscal->inflation regressions."
    exit 198
}
cap confirm variable d_oil
if _rc {
    di as error "Variable d_oil not found. Cannot run fiscal->inflation regressions."
    exit 198
}
cap confirm variable pi_import
if _rc {
    di as error "Variable pi_import not found. Cannot run fiscal->inflation regressions."
    exit 198
}

cap drop ok_fisc_reg ok_fisc_reg_w
gen ok_fisc_reg = inrange(tq, tq2001q1, tq2023q4) & !missing(pi_core_pce, rdpi_gap, u_gap, michigan_1y_median, d_oil, pi_import)

* Spec A: total effect + persistence
reg pi_core_pce L1.rdpi_gap L1.u_gap L1.michigan_1y_median d_oil pi_import L1.pi_core_pce if ok_fisc_reg, vce(robust)
estimates store SpecA_total
display "Spec A: 1% higher RDPI gap -> " 0.01*_b[L1.rdpi_gap] " pp higher core PCE inflation (annualised q/q)."

* Spec B: add wage inflation as mechanism (only if pi_eci exists)
gen ok_fisc_reg_w = ok_fisc_reg & !missing(pi_eci)

capture noisily reg pi_core_pce L1.rdpi_gap L1.u_gap L1.michigan_1y_median d_oil pi_import pi_eci L1.pi_core_pce if ok_fisc_reg_w, vce(robust)
if _rc==0 {
    estimates store SpecB_wages
    display "Spec B: 1% higher RDPI gap -> " 0.01*_b[L1.rdpi_gap] " pp higher core PCE inflation (annualised q/q)."
}

* Newey-West versions (4 lags)
capture noisily newey pi_core_pce L1.rdpi_gap L1.u_gap L1.michigan_1y_median d_oil pi_import L1.pi_core_pce if ok_fisc_reg, lag(4)
if _rc==0 estimates store SpecA_total_NW

capture noisily newey pi_core_pce L1.rdpi_gap L1.u_gap L1.michigan_1y_median d_oil pi_import pi_eci L1.pi_core_pce if ok_fisc_reg_w, lag(4)
if _rc==0 estimates store SpecB_wages_NW

* Timing robustness: 2-quarter lag
reg pi_core_pce L2.rdpi_gap L1.u_gap L1.michigan_1y_median d_oil pi_import L1.pi_core_pce if ok_fisc_reg, vce(robust)
estimates store SpecA_L2
display "Lag2: 1% higher RDPI gap -> " 0.01*_b[L2.rdpi_gap] " pp higher core PCE inflation (annualised q/q)."

*******************************************************
* END
*******************************************************
