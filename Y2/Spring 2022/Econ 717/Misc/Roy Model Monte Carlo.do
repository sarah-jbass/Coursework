clear all

**********************************************
*Model:
**********************************************

*Normal Roy model with 1 market

*(log_S_0,log_S_1) \sim N(mu,Sigma)

*W_j = pi_j * S_j = S_j, pi_j = 1 for j = 0,1

*D = 1 if W_1 > W_0, D = 0 otherwise

*W = W_1*D + (1-D)*W_0

*5 free parameters in mu, Sigma

**********************************************
*Step 0: set model parameters
**********************************************

local sigma_log_S_0 = 0.5
local sigma_log_S_1 = 1.5
local corr = 0.5 // instead of parameterizing the covariance
local mu0 = 1
local mu1 = 1.5 // Population treatment effect = mu1 - mu0 = 0.5
local obs = 10000

**********************************************;
*Step 1: generate S0,S1;
**********************************************;

set seed 456
matrix C = (1, `corr' \ `corr', 1) 
matrix sd = (`sigma_log_S_0',`sigma_log_S_1')
matrix m = (`mu0',`mu1')
drawnorm log_S_0 log_S_1, n(`obs') corr(C) sds(sd) means(m)

*check descriptive stats--does this match our parameters?
sum log_S_0 log_S_1
correlate log_S_0 log_S_1

**********************************************
*Step 2: generate W and D
**********************************************

*potential wages
gen log_W_0 = log_S_0
gen log_W_1 = log_S_1

*occupation choice
gen D = (log_W_1 >= log_W_0) // D = 1 if W_1 > W_0, D = 0 otherwise

*observed wage
gen log_W = D*log_W_1 + (1-D)*log_W_0 // W = W_1*D + (1-D)*W_0

*what is "treatment effect" of occupation 1 vs. occupation 0?
gen Delta = log_W_1 - log_W_0 // diff between wages
sum Delta

*selectivity into occupation
sum Delta if D == 1
sum Delta if D == 0

*twoway connected Delta1 Delta0 xvar1
kdensity Delta if D == 1, gen(xvar1 Delta1) nogr
kdensity Delta if D == 0, at(xvar1) gen(xvar0 Delta0) nogr

**********************************************
*summary stats, biased regressions
**********************************************

sum log_W 
regress log_W D // Effect of working in occupation 1 on wages. Not the ATE
sum Delta // Mean is the ATE. Way lower than OLS coef estimate! 

**********************************************
*what if D was randomly assigned?
**********************************************

local obs_half = `obs'/2

gen D_rand = 1 in 1/`obs_half' 
replace D_rand = 0 if D_rand == .

gen log_W_rand = D_rand*log_W_1 + (1-D_rand)*log_W_0

regress log_W_rand D_rand // This is much closer to the true ATE (0.5) after randomizing


**********************************************
*estimate mu1 only, assuming we know mu0, sigma0, sigma1, and corr
**********************************************

*Step 1: Compute sample moment
    sum log_W if D == 1 // mean wage for those who select occupation 1. Our "target"
    local sample_moment = r(mean)

*Step 2: Compute simulated data for each possible mu1 value

    *blank dataset to collect results
    clear
    set obs 1
    gen temp = .
    save results.dta, replace

    *random number draws, this stays fixed across different parameter tries
    local R = 1000 // number of simulation draws using, this is about computation time
    set seed 847

    forvalues mu1_try = -2(0.1)3 { // grid search

        clear

        *Simulated data: 
        matrix C = (1, `corr' \ `corr', 1)
        matrix sd = (`sigma_log_S_0',`sigma_log_S_1')
        matrix m = (`mu0',`mu1_try')
        drawnorm log_S_0_try log_S_1_try, n(`R') corr(C) sds(sd) means(m) // draw skills

        *The model: 
        gen log_W_0_try = log_S_0_try
        gen log_W_1_try = log_S_1_try
        gen D_try = (log_W_1_try >= log_W_0_try)
        gen log_W_try = D_try*log_W_1_try + (1-D_try)*log_W_0_try

    *Step 3: Compute simulated moment
        sum log_W_try if D_try == 1
        gen sim_moment = r(mean) 

        *collect results
        keep in 1
        gen mu1_try = `mu1_try'
        keep mu1_try sim_moment
        append using results.dta
        save results.dta, replace

    }

    *graph results
    gen sample_moment = `sample_moment'
    twoway connected sim_moment sample_moment mu1_try, xline(`mu1') // simulated moment changes under mu1
        // intersection point: identification (e.g. moments match)
        // use intersection to identify mu1

    gen diff = sim_moment - `sample_moment'
    gen diff_sq = diff^2
    twoway connect diff_sq mu1_try, xline(`mu1') // quadratic of the difference 
        // minimum point: identification (e.g. minimized diff between moments)
        // use min to identify mu1






