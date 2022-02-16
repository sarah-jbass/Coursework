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

*Question 3 - Linear Probability Model with Robust SE
    reg taken_new `covariates'
    outreg2 using q3_table, tex(frag) replace addtext(Robust S.E., No)

    reg taken_new `covariates', robust
    outreg2 using q3_table, tex(frag) append addtext(Robust S.E., Yes)

*Question 4 - Y-hat
    reg taken_new `covariates'

    predict taken_new_hat_lpm
    count if taken_new_hat_lpm<0 

*Question 5 - Weighted Least Squares
    *vwls taken_new `covariates' // No groups with sufficient observations
    *outreg2 using q5_table, tex(frag) replace

*Question 6 - Probit and Logit
    reg taken_new `covariates'
    outreg2 using q6_table, tex(frag) replace addtext(Model, LPM)

    probit taken_new `covariates'
    outreg2 using q6_table, tex(frag) append addtext(Model, Probit)

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

    *Estimating impact of tiny change in age
    gen taken_new_hat_probit_epsilon = normal(taken_new_hat_xb + 0.001*e(b)[1,1])
    gen client_age_partial_eps = (taken_new_hat_probit_epsilon - taken_new_hat_probit) / 0.001

    summarize client_age_partial_eps

*D
    probit taken_new `covariates' 
    margins, dydx(client_age) atmeans 

*Question 8 - LPM with age polynomial
    reg taken_new `covariates'
    outreg2 using q8_table, tex(frag) replace

    forvalues i = 2/4 {
        gen client_age_`i' = client_age^`i' // Generating higher order terms
    }

    reg taken_new `covariates' client_age_2 client_age_3 client_age_4 
    outreg2 using q8_table, tex(frag) append

    predict taken_new_hat_lpmq

    *Generate age + epsilon
    gen client_age_eps = client_age + 0.0001
    forvalues i = 2/4 {
        gen client_age_eps_`i' = client_age_eps^`i' // Generating higher order terms
    }

    reg taken_new client_age_eps client_married client_education hh_size hh_income muslim hindu_sc_kat treated client_age_eps_2 client_age_eps_3 client_age_4 
    predict taken_new_hat_lpmq_eps

    *Solve for difference
    gen client_age_partial_eps2 = taken_new_hat_lpmq_eps - taken_new_hat_lpmq

    summ client_age_partial_eps2  

*Question 9 - Calculate LRI
    probit taken_new `covariates'

    gen LRI = 1 - (e(ll)/e(ll_0))

*Question 10 - Calculate correct prediction rate
    *With cutoff = 0.5
    gen p_over_50 = (taken_new_hat_probit >=0.5) if taken_new_hat_probit  !=.
    tab taken_new p_over_50
    tab taken_new p_over_50, nofreq cell

    *With cutoff = 0.1673
    gen p_over_16 = (taken_new_hat_probit >=0.1673) if taken_new_hat_probit  !=.
    tab taken_new p_over_16
    tab taken_new p_over_16, nofreq cell   
    

*Question 11 - Calculate out of sample impact
    probit taken_new `covariates'
    outreg2 using q11_table, tex(frag) replace addtext(Full Sample, Yes)

    probit taken_new `covariates' if imidlineid<1400
    outreg2 using q11_table, tex(frag) append addtext(Full Sample, Out of Sample)

    predict taken_new_hat_outsample

    *With cutoff = 0.5
    gen p2_over_50 = (taken_new_hat_outsample >=0.5) if taken_new_hat_outsample  !=.
    tab taken_new p2_over_50 if imidlineid < 1400
    tab taken_new p2_over_50 if imidlineid < 1400, nofreq cell

    *With cutoff = 0.1673
    gen p2_over_16 = (taken_new_hat_outsample >=0.1673) if taken_new_hat_outsample  !=.
    tab taken_new p2_over_16 if imidlineid < 1400
    tab taken_new p2_over_16 if imidlineid < 1400, nofreq cell   


*Question 12 - Probit with married*Muslim interaction
    probit taken_new `covariates'
    outreg2 using q12_table, tex(frag) replace addtext(Interaction, No)
    gen married_muslim = client_married * muslim

    probit taken_new `covariates' married_muslim
    outreg2 using q12_table, tex(frag) append addtext(Interaction, Yes)

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
    *No code for this question, see margins table from 13

*Question 15 - Heteroskedasticity test
    reg taken_new `covariates'

    *Find residuals
    predict residuals, residuals
    gen residuals2 = residuals^2 

    *Regress on residuals on X
    reg residuals2 `covariates'
    outreg2 using q15_table, tex(frag) replace

*Question 16 - Het prob
    hetprob taken_new `covariates', het(client_age client_education) 
    outreg2 using q16_table, tex(frag) replace

log close 




