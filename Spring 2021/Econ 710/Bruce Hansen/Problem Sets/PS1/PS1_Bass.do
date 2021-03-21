/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 710B
Task: PS1
Date created: 3/21/21
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS1\PS1_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS1"

*Question 13.28
use "Card1995.dta"

gen lwage = lwage76
gen edu = ed76
gen exp = age76 - edu - 6
gen exp2per = exp^2/100
gen south = reg76r
gen urban = smsa76r
gen public = nearc4a
gen private = nearc4b
gen pubage = nearc4a*age76
gen pubage2 = nearc4a*age76^2/100

foreach x in 2sls gmm {
    ivregress `x' lwage exp exp2per south black urban (edu = public private), r
}
estat overid

foreach x in 2sls gmm {
    ivregress `x' lwage exp exp2per south black urban (edu = public private pubage pubage2), r
}
estat overid

*Question 17.15
clear
use "AB1991.dta"

xtabond k, lags(1) vce(robust)

xtdpdsys k, lags(1) vce(robust)