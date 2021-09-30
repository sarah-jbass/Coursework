/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 710B
Task: PS12
Date created: 4/21/21
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS12\PS12_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS12"

*Question 27.9
use "CHJ2004.dta" 

replace transfers = transfers/1000
replace transfers=0 if transfers<=0 // Assuming this is what I'm supposed to do 
replace income = income/1000
gen Dincome = (income-1) if income>1
replace Dincome = 0 if income <=1

*a
reg transfers income Dincome 

*b
summ transfers

count 
local n_obs = r(N)
count if transfers==0
display r(N)/`n_obs'

*c
drop if transfers == 0
reg transfers income Dincome 

*d 
tobit transfers income Dincome, ll

*e
clad transfers income Dincome, ll

*Question 28.12
clear
use "cps09mar.dta" 

keep if female==1
keep if hisp==1
count

gen wage = earnings/(week*hours)
gen lnwage = log(wage)
gen exp = age -education-6
foreach x in 2 3 4 5 6 {
	gen exp`x' = exp^`x'
}
gen college = (education >= 16)
gen educ_spline = (education >= 9)*(education - 9)
foreach x in 12 13 14 16 18 20 {
	gen educ`x' = (education == `x')
}
gen married =  (marital == 1 | marital == 2 | marital == 3)

mat criterion = J(9,2,.)
mat rownames criterion = "Model 1" "Model 2" "Model 3" "Model 4" "Model 5" "Model 6" "Model 7" "Model 8" "Model 9"
mat colnames criterion = "AIC" "BIC"

*Model 1
reg lnwage married i.region exp exp2 college
estat ic
mat criterion[1,1] = r(S)[1,5]
mat criterion[1,2] = r(S)[1,6]

*Model 2
reg lnwage married i.region exp exp2 educ_spline
estat ic
mat criterion[2,1] = r(S)[1,5]
mat criterion[2,2] = r(S)[1,6]

*Model 3
reg lnwage married i.region exp exp2 educ12 educ13 educ14 educ16 educ18 educ20
estat ic
mat criterion[3,1] = r(S)[1,5]
mat criterion[3,2] = r(S)[1,6]

*Model 4
reg lnwage married i.region exp exp2 exp3 exp4 college
estat ic
mat criterion[4,1] = r(S)[1,5]
mat criterion[4,2] = r(S)[1,6]

*Model 5
reg lnwage married i.region exp exp2 exp3 exp4 educ_spline
estat ic
mat criterion[5,1] = r(S)[1,5]
mat criterion[5,2] = r(S)[1,6]

*Model 6
reg lnwage married i.region exp exp2 exp3 exp4 educ12 educ13 educ14 educ16 educ18 educ20
estat ic
mat criterion[6,1] = r(S)[1,5]
mat criterion[6,2] = r(S)[1,6]

*Model 7
reg lnwage married i.region exp exp2 exp3 exp4 exp5 exp6 college
estat ic
mat criterion[7,1] = r(S)[1,5]
mat criterion[7,2] = r(S)[1,6]

*Model 8
reg lnwage married i.region exp exp2 exp3 exp4 exp5 exp6 educ_spline
estat ic
mat criterion[8,1] = r(S)[1,5]
mat criterion[8,2] = r(S)[1,6]

*Model 9
reg lnwage married i.region exp exp2 exp3 exp4 exp5 exp6 educ12 educ13 educ14 educ16 educ18 educ20
estat ic
mat criterion[9,1] = r(S)[1,5]
mat criterion[9,2] = r(S)[1,6]

mat list criterion