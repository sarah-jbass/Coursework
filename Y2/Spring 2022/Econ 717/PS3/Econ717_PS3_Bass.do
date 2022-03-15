/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 717
Task: PS3
Date created: 3/15/22
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS3\Econ717_PS3_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
capture log close
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS3"
*log using "Econ717_PS3_Bass.log", replace

use "Economics 717 Miron and Tetelbaum Data.dta"
rename _all, lower

*Question 1
    xtset state year

*Question 2
    gen mlda21 = (mlda==21) 

*Question 3
    reg rate18_20ht mlda21, robust
        outreg2 using q3_table, tex(frag) replace

*Question 4 
    *State FE
    reg rate18_20ht mlda21 i.state, robust
        *How to export results? this table is way too long
        outreg2 using q4_table, tex(frag) replace 

    *Year FE
    reg rate18_20ht mlda21 i.year, robust
        outreg2 using q4_table, tex(frag) append

*Question 5
    *State and year FE, state cluster
    reg rate18_20ht mlda21 i.state i.year, cluster(state) robust
        outreg2 using q5_table, tex(frag) replace

*Question 6
    *State and year FE
    reg rate18_20ht mlda21 i.state i.year, robust
        outreg2 using q6_table, tex(frag) replace

*Question 7
    preserve

    *Only using data from before 1990
    drop if year>1990 
    
    *State and year FE, state cluster
    reg rate18_20ht mlda21 i.state i.year, cluster(state) robust
        outreg2 using q7_table, tex(frag) replace 

    restore

*Question 8
    *Placebo treatment indicator
    gen placebo82 = (mldayr == 1987 & year >= 1982) 

    preserve    
    keep if year <= 1987 & (mlda == 21 | mldayr == 1987)

    *State and year FE, state cluster
    reg rate18_20ht mlda21 i.state i.year, cluster(state) robust
        outreg2 using q5_table, tex(frag) replace

    restore 

*Question 9


*log close