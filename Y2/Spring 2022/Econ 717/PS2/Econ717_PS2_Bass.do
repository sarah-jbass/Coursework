/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 717
Task: PS2
Date created: 2/28/22
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS2\Econ717_PS2_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on
ssc install outreg2

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS2"
*log using "Econ717_PS2_Bass.log", replace

use "Economics 717 Spring 2022 NSW Data.dta"
rename _all, lower

*Question 0
drop if sample == 3

*Question 1
gen age2 = age^2
reg re78 treated age educ black hisp married nodegree re74 re75
outreg2 using q1_table, tex(frag) replace

*Question 2
drop if treated == 1 & sample==1 

*Question 3
gen treated2 = (sample ==1) // creating a variable for our definition of the treatment group

probit treated2 age age2 educ black hisp married nodegree 
outreg2 using q3_table, tex(frag) replace
predict pscorea

probit treated2 age age2 educ black hisp married nodegree re74 re75
outreg2 using q3_table, tex(frag) append
predict pscoreb

*Question 4
summ pscorea 
summ pscoreb

*Question 5
histogram pscorea, title(Propensity Score A)
graph export q4_pscorea.png

histogram pscoreb, title(Propensity Score B)
graph export q4_pscoreb.png

/*
*This is Jeff's histogram Code:

egen bins=cut(pscore), at(0(.05)1) icodes 
graph bar (count) pscore, over(treated2) over(bins, label(nolab)) asyvars
graph export q4_scoreb.png
*/

*Question 6

*Nearest neighbor w/o replacement, w/ common support
psmatch2 treated2 age age2 educ black hisp married nodegree, noreplacement common pscore(pscore6a)

*Compared to nearest neighbor w/o replacement, w/o common support
count if _support==1
by _support: summarize age age2 educ black hisp married nodegree
outreg2 using q6_table, tex(frag) replace

*How close are resulting estimates to the experimental impact estimate (Question 3)? 
egen bins=cut(pscore), at(0(.05)1) icodes 
graph bar (count) pscore, over(d) over(bins, label(nolab)) asyvars

*Do the rich scores that contain pre-treatment earnings perform better than coarse scores that do not? 
/* *Rich scores that contain pre-treatment earnings
psmatch2 treated age age2 educ black hisp married nodegree re74 re75, noreplacement common pscore(pscore6a)
 */

































