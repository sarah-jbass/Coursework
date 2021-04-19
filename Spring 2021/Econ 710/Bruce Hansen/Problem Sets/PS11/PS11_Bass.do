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
gen age25 = (age >= 25)*(age - 25)
gen age30 = (age >= 30)*(age - 30)
gen age35 = (age >= 35)*(age - 35)
gen age40 = (age >= 40)*(age - 40)
gen age45 = (age >= 45)*(age - 45)
gen age50 = (age >= 50)*(age - 50)
gen age55 = (age >= 55)*(age - 55)
gen age60 = (age >= 60)*(age - 60)

probit y age age25 age30 age35 age40 age45 age50 age55 age60

collapse (count) n_obs = female (sum) n_married = y, by(age)

gen y_avg = n_married / n_obs // percent of eage age married

*Splines every 5
gen age25 = (age >= 25)*(age - 25)
gen age30 = (age >= 30)*(age - 30)
gen age35 = (age >= 35)*(age - 35)
gen age40 = (age >= 40)*(age - 40)
gen age45 = (age >= 45)*(age - 45)
gen age50 = (age >= 50)*(age - 50)
gen age55 = (age >= 55)*(age - 55)
gen age60 = (age >= 60)*(age - 60)

predict y_pred
twoway line y_pred age || scatter y_avg age, msize(tiny) mcol(black) 	///
	xlab(20(10)80) ylab(0(.2)1, angle(horizontal) grid nogextend) xtitle("Age") 	///
	plotregion(margin(tiny) lcolor(none)) ytitle("")						///
	title("Probability of Being Married: Women with College Degrees")		///
	leg(lab(1 "Fitted values") lab(2 "Empirical proportion") 	///
			region(lstyle(none)))
			
graph export q25_17.png, replace