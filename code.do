* Stage 1 - understanding the data and preparing for analysis*


* Importing the data from the file*

import delimited "alldata.csv", clear varnames(1)

  
* Undestanding the data – columns have missing values that start at different time periods. It's important to see which periods correspond to which variables *

misstable summarize spf_inflation_1year michigan_1y_median swap_1year swap_5year swap_10year cpilfesl pceexcludingfoodandenergy pceall unrate nrou v_u poilbreusdm import_def

  
* Constructing the time variable and tsset *

gen tq = yq(year, quarter)
format tq %tq
tsset tq

  
* Constructing inflation measures – Core CPI and Core PCE inflation *

gen pi_core_cpi = 400*(ln(cpilfesl) - ln(L.cpilfesl))
gen pi_core_pce = 400*(ln(pceexcludingfoodandenergy) - ln(L.pceexcludingfoodandenergy))

* Constructing Slack measures – unemployment rate *

gen u_gap = unrate - nrou


* Constructing Supply shocks variables – oil price inflation and import price inflation *

gen d_oil = 400*(ln(poilbreusdm) - ln(L.poilbreusdm))
gen pi_import = 400*(ln(import_def) - ln(L.import_def))

  
* Constructing a dummy variable for Pandemic time perios *

gen post2020 = (tq >= yq(2020,1))


* Summarising new variables in a table to get better understanding of them *

summarize pi_core_cpi pi_core_pce u_gap spf_inflation_1year michigan_1y_median v_u d_oil pi_import




