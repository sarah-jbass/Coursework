/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 710B
Task: PS9
Date created: 4/1/21
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS9\PS9_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS9"

*Question 20.11
use "cps09mar.dta"

gen logwage = log(earnings/hours*week)
quietly summ education
replace education = (education - r(min))/(r(max) - r(min))

*lpoly logwage education, degree(6) ci noscatter legend(off) title("Polynomial Regression")
*graph export "20-11.png", as(png) name("Graph") replace

foreach i in 2 3 4 5 6 {
    gen educ_`i' = education^`i'
}

reg logwage education educ_2 educ_3 educ_4 educ_5 educ_6, r
predict yhat, xb
predict stdy, stdp
gen y_ub = yhat + 1.96*stdy
gen y_lb = yhat - 1.96*stdy

twoway rarea y_ub y_lb education, sort xtitle("Rescaled Education") ytitle("Log(wage)") title("Polynomial Regression") legend(off)|| line yhat education, sort legend(off)
graph export "20-11.png", as(png) name("Graph") replace

*Question 20.15
clear
use "RR2010.dta", clear

gen laggdp = l.gdp
gen lagdebt = l.debt
foreach i in 40 60 80 {
    gen lagdebt`i' = 0
    replace lagdebt`i' = lagdebt-`i' if lagdebt>=`i'
}

*Linear Regression
reg gdp laggdp lagdebt
predict yhat, xb
gen y0 = yhat
estat ic
drop yhat

*Spline, knot at 60
reg gdp laggdp lagdebt lagdebt60
predict yhat, xb
predict stdy, stdp
gen y_ub = yhat + 1.96*stdy
gen y_lb = yhat - 1.96*stdy
gen y1 = yhat
estat ic
drop yhat

*Spline, knots at 40, 80
reg gdp laggdp lagdebt lagdebt40 lagdebt80
predict yhat, xb
estat ic
gen y2 = yhat


line y0 y1 y2 lagdebt, sort xtitle("Lagged Debt") ytitle("Growth") title("Regression Estimates") legend(cols(1))
graph export "20-15a.png", as(png) name("Graph") replace

twoway rarea y_ub y_lb lagdebt, sort xtitle("Lagged Debt") ytitle("Growth") title("One Knot Spline with CI") legend(off) || line y1 lagdebt, sort legend(off)
graph export "20-15b.png", as(png) name("Graph") replace

clear