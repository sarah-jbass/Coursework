/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 710B
Task: PS10
Date created: 4/12/21
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS10\PS10_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS10"

*Question 23.8
use "PSS2017.dta"

gen lny = log(EG_total)
rename EC_c_alt x1
rename EC_d_alt x2
drop if lny==. | x1==. | x2==.

nl (lny = {beta} + ({nu}/{rho})*log({alpha}*x1^{rho} + (1-{alpha})*x2^{rho})), initial(beta 1.5 nu 1 rho 0.5 alpha 0.5)

*Question 24.14
clear
use "cps09mar.dta"

keep if female == 1 & hisp == 1
drop if education <11
gen wage 	= earnings / (hours*week)
gen lnwage	= log(wage)

foreach q in .1 .3 .5 .7 .9 {
	qui qreg lnwage education, q(`q')
    local i = `q'*10
	local b0_`i' = _b[_cons]
	local b1_`i' = _b[education]
}

twoway function y = `b0_1' + `b1_1'*x, range(education) || function y = `b0_3' + `b1_3'*x, range(education) || ///
function y = `b0_5' + `b1_5'*x, range(education) || function y = `b0_7' + `b1_7'*x, range(education)  || ///
function y = `b0_9' + `b1_9'*x, range(education) legend(label(1 "10%") label(2 "30%") label(3 "50%") label(4 "70%") label(5 "90%") cols(1)) ///
xlabel(11 12 13 14 15 16 17 18 19 20) title("Quantile Regression Results")

graph export q24_14.png, replace