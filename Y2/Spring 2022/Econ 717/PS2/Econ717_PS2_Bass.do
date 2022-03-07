/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 717
Task: PS2
Date created: 2/28/22
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS2\Econ717_PS2_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
capture log close
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS2"
log using "Econ717_PS2_Bass.log", replace

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
    predict pscorea
        outreg2 using q3_table, tex(frag) replace
    

    *Rich propensity scores
    probit treated2 age age2 educ black hisp married nodegree re74 re75
    predict pscoreb
        outreg2 using q3_table, tex(frag) append
    

*Question 4
    *Look at distributions of each propensity scores
    bysort treated2: summ pscorea, detail
    bysort treated2: summ pscoreb, detail

*Question 5
    /* *Histogram for coarse scores
    egen binsa=cut(pscorea), at(0(.05)1) icodes 
    graph bar (count) pscorea if binsa>0, over(treated2) over(binsa, label(nolab)) asyvars title("Coarse Propensity Scores")
        graph export q4_scorea.png, replace

    *Histogram for rich scores
    egen binsb=cut(pscoreb), at(0(.05)1) icodes 
    graph bar (count) pscoreb if binsb>0, over(treated2) over(binsb, label(nolab)) asyvars title("Rich Propensity Scores")
        graph export q4_scoreb.png, replace */

*Question 6
    global table_number=6 // all of this global stuff is just to make tables

    *Nearest neighbor w/o replacement, w/ common support
    est clear
    eststo: psmatch2 treated2, noreplacement outcome(re78) pscore(pscorea) neighbor(1) common 
        global att_coarse = r(att)
        global att_coarse_se = r(seatt)
        count if _support==0 // How many are not on common support? 

    esttab, se
        global unm_coef = r(coefs)[1, 1]
        global unm_se = r(coefs)[1, 2]

    *"Rich" matching estimates w/o replacement, w/ common support
    psmatch2 treated2, noreplacement outcome(re78) pscore(pscoreb) neighbor(1) common 
        global att_fine = r(att)
        global att_fine_se = r(seatt)
        count if _support==0 // How many are not on common support w/ the rich scores?
   
   texdoc do table_maker

*Question 7
    global table_number=7

    *Nearest neighbor w/ replacement, w/ common support
    est clear
    eststo: psmatch2 treated2, outcome(re78) pscore(pscorea) neighbor(1) common 
        global att_coarse = r(att)
        global att_coarse_se = r(seatt)
        count if _support==0 // How many are not on common support? 
    
    esttab, se
        global unm_coef = r(coefs)[1, 1]
        global unm_se = r(coefs)[1, 2]

    *"Rich" matching estimates w/ replacement, w/ common support
    psmatch2 treated2, outcome(re78) pscore(pscoreb) neighbor(1) common 
        global att_fine = r(att)
        global att_fine_se = r(seatt)
        count if _support==0 // How many are not on common support w/ the rich scores?

   texdoc do table_maker

*Question 8
    *Raw data standard difference 
    stddiff re74 re75, by(treated2)
    
    *Std diff for re74 using weights = diff/std err
    psmatch2 treated2, outcome(re74) pscore(pscoreb) neighbor(1) common 

    *Std diff for re75 using weights = diff/std err
    psmatch2 treated2, outcome(re75) pscore(pscoreb) neighbor(1) common 

*Question 9
    global table_number = 9

    psmatch2 treated2, kernel kerneltype(normal) outcome(re78) pscore(pscoreb) bwidth(0.02) common 
        global att_low = r(att)
        global att_low_se = r(seatt)

    psmatch2 treated2, kernel kerneltype(normal) outcome(re78) pscore(pscoreb) bwidth(0.2) common 
        global att_med = r(att)
        global att_med_se = r(seatt)

    psmatch2 treated2, kernel kerneltype(normal) outcome(re78) pscore(pscoreb) bwidth(2) common  
        global att_hi = r(att)
        global att_hi_se = r(seatt)

    texdoc do table_maker_2

*Question 10
    global table_number = 10

    psmatch2 treated2, llr outcome(re78) pscore(pscoreb) bwidth(0.02) common 
        global att_low = r(att)
        global att_low_se = r(seatt)

    psmatch2 treated2, llr outcome(re78) pscore(pscoreb) bwidth(0.2) common 
        global att_med = r(att)
        global att_med_se = r(seatt)

    psmatch2 treated2, llr outcome(re78) pscore(pscoreb) bwidth(2) common 
        global att_hi = r(att)
        global att_hi_se = r(seatt)

    texdoc do table_maker_2
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
    scalar p_hat = r(mean)

    *Find y * d
    gen y_d = re78 * treated2
    summarize y_d
    scalar y_d_sum = r(sum)

    *Find second term wo rescaling
    gen wo_rescale = (1 - p_hat) / p_hat * pscoreb * re78 * (1 - treated2)/(1 - pscoreb)
    summarize wo_rescale
    scalar wo_rescale_sum = r(sum)

    *Find second term w rescaling
    gen rescaling = pscoreb * (1-treated2)/(1-pscoreb) 
    summarize rescaling
    scalar rescaling = 1/n_0 * r(sum)
    gen w_rescale = (1/rescaling) * pscoreb * re78 * (1 - treated2)/(1-pscoreb)
    summarize w_rescale
    scalar w_rescale_sum = r(sum)

    *Treatment effect on treated ipw w and wo rescaling
    display 1/n_1 * y_d_sum - 1/n_0 * wo_rescale_sum
    display 1/n_1 * y_d_sum - 1/n_0 * w_rescale_sum

log close