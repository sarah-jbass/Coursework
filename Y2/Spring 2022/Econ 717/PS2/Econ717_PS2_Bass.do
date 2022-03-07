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
*ssc install outreg2

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS2"
*log using "Econ717_PS2_Bass.log", replace

use "Economics 717 Spring 2022 NSW Data.dta"
rename _all, lower

*Question 0
    drop if sample == 3

*Question 1
    gen age2 = age^2

    *Regression of earnings using the original treatment definition
    reg re78 treated age age2 educ black hisp married nodegree re74 re75
    outreg2 using q1_table, tex(frag) replace

*Question 2
    drop if treated == 1 & sample==1 

*Question 3
    gen treated2 = (sample ==1) // creating a variable for our definition of the treatment group

    *Coarse propensity scores
    probit treated2 age age2 educ black hisp married nodegree 
    outreg2 using q3_table, tex(frag) replace
    predict pscorea

    *Rich propensity scores
    probit treated2 age age2 educ black hisp married nodegree re74 re75
    outreg2 using q3_table, tex(frag) append
    predict pscoreb

*Question 4
    *Look at distributions of each propensity scores
    bysort treated2: summ pscorea, detail
    bysort treated2: summ pscoreb, detail

*Question 5
    *Histogram for coarse scores
    egen binsa=cut(pscorea), at(0(.05)1) icodes 
    graph bar (count) pscorea if binsa>0, over(treated2) over(binsa, label(nolab)) asyvars
    graph export q4_scorea.png, replace

    *Histogram for rich scores
    egen binsb=cut(pscoreb), at(0(.05)1) icodes 
    graph bar (count) pscoreb if binsb>0, over(treated2) over(binsb, label(nolab)) asyvars
    graph export q4_scoreb.png, replace

*Question 6
    *Nearest neighbor w/o replacement, w/ common support
    eststo: psmatch2 treated2, noreplacement outcome(re78) pscore(pscorea) neighbor(1) common 
        count if _support==0 // How many are not on common support? 
        scalar diff_att = r(att)
        scalar se_att = r(seatt)

    *"Rich" matching estimates w/o replacement, w/ common support
    eststo: psmatch2 treated2, noreplacement outcome(re78) pscore(pscoreb) neighbor(1) common 
        count if _support==0 // How many are not on common support w/ the rich scores?

    ********************
    *ADD TABLE FOR THIS
    ********************

*Question 7
    *Nearest neighbor w/ replacement, w/ common support
    psmatch2 treated2, outcome(re78) pscore(pscorea) neighbor(1) common 
    count if _support==0 // How many are not on common support? 

    *"Rich" matching estimates w/ replacement, w/ common support
    psmatch2 treated2, outcome(re78) pscore(pscoreb) neighbor(1) common 
    count if _support==0 // How many are not on common support w/ the rich scores?

*Question 8
    *Raw data standard difference 
    stddiff re74 re75, by(treated2)
    
    *Std diff for re74 using weights = diff/std err
    psmatch2 treated2, outcome(re74) pscore(pscoreb) neighbor(1) common 

    *Std diff for re75 using weights = diff/std err
    psmatch2 treated2, outcome(re75) pscore(pscoreb) neighbor(1) common 

*Question 9
    /* psmatch2 treated2, kernel kerneltype(normal) outcome(re78) pscore(pscoreb) bwidth(0.02) common 
    psmatch2 treated2, kernel kerneltype(normal) outcome(re78) pscore(pscoreb) bwidth(0.2) common 
    psmatch2 treated2, kernel kerneltype(normal) outcome(re78) pscore(pscoreb) bwidth(2) common  */

*Question 10
    /* psmatch2 treated2, llr outcome(re78) pscore(pscoreb) bwidth(0.02) common 
    psmatch2 treated2, llr outcome(re78) pscore(pscoreb) bwidth(0.2) common 
    psmatch2 treated2, llr outcome(re78) pscore(pscoreb) bwidth(2) common  */

*Question 11
    *Reg controlling for whether in treatment group
    reg re78 treated2 age age2 educ black hisp married nodegree re74 re75
    outreg2 using q11_table, tex(frag) replace
    outreg2 using q12_table, tex(frag) replace

    *Predicted y-hat w/ control for treatment group
    predict re78hat

*Question 12
    *Reg only using comparison sample
    reg re78 age age2 educ black hisp married nodegree re74 re75 if treated2==0
    outreg2 using q12_table, tex(frag) append

    *Predicted y-hat using only comparison sample
    predict re78hat_oos

    *Look at differences in y-hat for treated2==1
    bysort treated2: summarize re78hat re78hat_oos 

*Question 13
    *Generating n0 n1
    count if treated2==0
    scalar n_0 = r(N)

    count if treated2==1
    scalar n_1 = r(N)

    *Generating p-hat
    summarize treated2
    scalar p-hat = r(mean)

    *Generating effects
    gen YiDi = re78 * treated2
    gen YiDj = 
    
