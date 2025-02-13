/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 710B
Task: PS11
Date created: 4/12/21
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS11\PS11_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS11"

*Question 25.15
use "cps09mar.dta"
preserve 

gen y = (union==1)
gen black = (race==2) // Not including mixed race with black in this sample
drop if female==1
probit y age education hisp black

*Question 25.17
restore
drop if education < 16 | female!=1
gen y = (marital == 1 | marital == 2 | marital == 3)

*Splines every 5
foreach x in 25 30 35 40 45 50 55 60 {
	gen age`x' = (age >= `x')*(age - `x')
}

probit y age age25 age30 age35 age40 age45 age50 age55 age60

collapse (max) age25 age30 age35 age40 age45 age50 age55 age60 (count) n_obs = female (sum) n_married = y, by(age)

*Actual percent of each age married
gen y_avg = n_married / n_obs

*Estimated percent of each age married
predict y_pred

twoway line y_pred age || scatter y_avg age, msize(tiny) mcol(black) 	///
	xlab(20(10)80) ylab(0(.2)1, angle(horizontal)) xtitle("Age")  ytitle("")	///
	title("Probability of Being Married: Women with College Degrees")		///
	leg(lab(1 "Fitted values") lab(2 "Empirical proportion"))
			
graph export q25_17.png, replace