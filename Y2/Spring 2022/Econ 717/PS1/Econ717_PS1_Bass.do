/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 710B
Task: PS12
Date created: 4/21/21
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS1\Econ717_PS1_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on
ssc install outreg2

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS1"
log using "Econ717_PS1_Bass.log", replace

use "Field et al (2010) Analysis Sample.dta"
rename _all, lower

*Question 1
    drop if miss_client_age==1 | miss_client_married==1 | miss_client_education==1 | hh_income==. // Note there are households w/ zero income in the sample

*Question 2 - Linear Probability Model
    local covariates "client_age client_married client_education hh_size hh_income muslim hindu_sc_kat treated"
    reg taken_new `covariates'
    outreg2 using q2_table, tex(frag) replace

    foreach x in _cons `covariates' {
        local b_`x' = _b[`x'] // Saving the coefficients as locals for question 4
    }

*Question 3 - Linear Probability Model with Robust SE
    reg taken_new `covariates', robust
    outreg2 using q3_table, tex(frag) replace

*Question 4
    reg taken_new `covariates'

    predict taken_new_hat_lpm
    count if taken_new_hat_lpm<0 

*Question 5 - Weighted Least Squares
    *vwls taken_new `covariates' // No groups with sufficient observations
    *outreg2 using q5_table, tex(frag) replace

*Question 6 - Probit and Logit
    probit taken_new `covariates'
    outreg2 using q6_table, tex(frag) replace addtext(Model, Probit)

    logit taken_new `covariates'
    outreg2 using q6_table, tex(frag) append addtext(Model, Logit)

*Question 7 - Marginal effects

*A
    dprobit taken_new `covariates'
    outreg2 using q7_table, tex(frag) replace

*B
    probit taken_new `covariates'

    predict taken_new_hat_xb, xb 
    gen client_age_partial = normalden(taken_new_hat_xb) * e(b)[1,1] // Formula from notes

    summ client_age_partial

*C
    predict taken_new_hat_probit

    *Change client age by epsilon and recalculate
    gen client_age_eps = client_age + 0.0001
    probit taken_new client_age_eps client_married client_education hh_size hh_income muslim hindu_sc_kat treated
    predict taken_new_hat_probit_eps

    *Solve for difference
    gen client_age_partial_eps = taken_new_hat_probit_eps - taken_new_hat_probit

    summ client_age_partial_eps

*D
    probit taken_new `covariates' 
    margins, dydx(client_age) atmeans 

*Question 8 - LPM with age quartic
    forvalues i = 2/4 {
        gen client_age_`i' = client_age^`i' // Generating higher order terms
    }

    reg taken_new `covariates' client_age_2 client_age_3 client_age_4 
    outreg2 using q8_table, tex(frag) replace

    predict taken_new_hat_lpmq

    *Change client age by epsilon and recalculate
    reg taken_new client_age_eps client_married client_education hh_size hh_income muslim hindu_sc_kat treated client_age_2 client_age_3 client_age_4 
    predict taken_new_hat_lpmq_eps

    *Solve for difference
    gen client_age_partial_eps2 = taken_new_hat_lpmq_eps - taken_new_hat_lpmq

    summ client_age_partial_eps2 

*Question 9 - Calculate LRI
    probit taken_new `covariates'

    gen LRI = 1 - (e(ll)/e(ll_0))

*Question 10 - Calculate correct prediction rate
    * = 1 if p>0.5, = 0 if p<0.5
    gen p_over_50 = (taken_new_hat_probit >=0.5) if taken_new_hat_probit  !=.

    *Number of correct predictions
    count if (p_over_50 == 1 & taken_new == 1) | (p_over_50 == 0 & taken_new == 0)
    local num_right = r(N)

    *Number of incorrect predictions
    count if (p_over_50 == 0 & taken_new == 1) | (p_over_50 == 1 & taken_new == 0)
    local num_wrong = r(N)   

    *Calculate correct rate
    gen correct_rate = `num_right'/`num_wrong'

*Question 11
    probit taken_new `covariates' if imidlineid<1400
    outreg2 using q11_table, tex(frag) replace

    predict taken_new_hat_outsample

    * = 1 if p>0.5, = 0 if p<0.5
    gen p_over_50_outsample = (taken_new_hat_outsample >=0.5) if taken_new_hat_outsample  !=.

    *Number of correct predictions
    count if (p_over_50_outsample == 1 & taken_new == 1) | (p_over_50_outsample == 0 & taken_new == 0)
    local num_right_outsample = r(N)

    *Number of incorrect predictions
    count if (p_over_50_outsample == 0 & taken_new == 1) | (p_over_50_outsample == 1 & taken_new == 0)
    local num_wrong_outsample = r(N)   

    *Calculate correct rate
    gen correct_rate_outsample = `num_right_outsample'/`num_wrong_outsample'

*Question 12 - Probit with married*Muslim interaction
    gen married_muslim = client_married * muslim

    probit taken_new `covariates' married_muslim
    outreg2 using q12_table, tex(frag) replace

*Question 13 - Mean finite differences 
    margins, dydx(married_muslim)

    *Compute by hand
    predict index_hat, xb

    gen index_hat_0 = index_hat - client_married * e(b)[1,2] - muslim * e(b)[1,6] - married_muslim *e(b)[1, 9] // not muslim or married
    gen index_hat_01 = index_hat_0 + e(b)[1,2] // married but not muslim
    gen index_hat_02 = index_hat_0 + e(b)[1,6] // muslim but not married 
    gen index_hat_012 = index_hat_0 + e(b)[1,2] + e(b)[1,6] + e(b)[1,9] // married and muslim 
    
    gen finite_difference = (normal(index_hat_012) - normal(index_hat_02)) - (normal(index_hat_01) - normal(index_hat_0))
    
    summ finite_difference

*Question 14 - Standard deviation
    *No code for this question, see table

*question 15 - Heteroskedasticity test
    reg taken_new `covariates'

    *Find residuals
    predict residuals, residuals
    gen residuals2 = residuals^2 

    *Regress on residuals on X
    reg residuals2 `covariates'
    outreg2 using q15_table, tex(frag) replace

*Question 16 - Het prob
hetprob taken_new `covariates', het(client_age client_education)

log close 




