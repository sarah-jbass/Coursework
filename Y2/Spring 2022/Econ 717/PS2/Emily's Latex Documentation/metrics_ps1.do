/*
	applied metrics spring 2022
	problem set 1, due february 16
	emily case 
	
	note: the tables created in this do file use eststo and esttab
*/


clear all
cd "/Users/emilycase/Desktop/metrics sp22/ps1"

use field_sample_data //i renamed the data
rename *, lower 	  //aesthetic preferences matter
keep imidlineid taken_new client_age client_married client_education hh_size hh_income muslim hindu_sc_kat treated miss_client_age miss_client_education miss_client_married

* make shorter variable labels, for LaTeX tables' sake
la var client_age "Client age"
la var client_married "Client marital status"
la var client_education "Client education"
la var hh_size "Household size"
la var hh_income "Household income"
la var muslim "Indicator for Muslim"
la var hindu_sc_kat "Indicator for Hindi"
la var treated "Treated"
la var taken_new "Client took a loan"

* global list of variables, starting with dependent (less typing = good)
global vars "taken_new client_age client_married client_education hh_size hh_income muslim hindu_sc_kat treated"


*****************
*	problem 1   *
*****************

drop if miss_client_married   == 1
drop if miss_client_education == 1
drop if miss_client_age		  == 1
drop if hh_income 			  == .


*****************
*	problem 2   *
*****************
/*
	Estimate a linear probability model with loan 
	take-up as the dependent variable and client age, client marital 
	status, client years of schooling, client household size, client 
	household income, indicators for Muslim and Hindu scheduled caste, 
	and experimental treatment status. Estimate the model without the 
	robust option. Report and discuss the coefficient estimates.
*/

eststo clear //start from a blank slate
eststo lpm: reg $vars 

esttab using "results/q2.tex", replace ///tell stata where to save the file
	   mtitles("LPM") ///rename the columns: one regression = one col
	   title(Problem 2\label{q2}) ///table title and LaTeX label
	   b(4) se(a2) /// 4 decimal places for coeffs(b) & minimum of 2 sig figs in se. also, by typing se at all, i am telling it that the secondary thing to report is standard errors instead of t statistics. the main thing it's reporting is b (the estimate coefficients)
	   nostar label booktabs nonumbers //formats it a little better

//note: i am not clearing eststo yet, because in the future problems, i want to be able to have tables that compare my results. 
	   
*****************
*	problem 3   *
*****************
/*
	Repeat the estimation in Problem 2 with the robust option.
*/

eststo lpmrobust: reg $vars , robust 
//check what is in eststo by typing eststo dir . this will show that now there is 2 estimates stored in there, so whhen i use esttab it will have 2 columns

esttab using "results/q3.tex", replace ///
	   mtitles("LPM" "LPM, robust SE") ///added column name
	   title(Problem 3\label{q3}) ///
	   b(4) se(a2) /// 4 decimal places for coeffs & min 2 sig figs in se
	   nostar label booktabs nonumbers //formats it a little better



*****************
*	problem 4   *
*****************
/*
	Generate predicted probabilities of loan take-up using the estimated 
	coefficients from Problem 2. Do any of these probabilities lie 
	outside [0, 1]? If so, do the observations corresponding to 
	these values show any particular patterns in the values of the 
	variables?
*/

* repeat the regression from problem 2
reg $vars 

* predict probabilities based on most recent regression
predict ptaken_new, xb

* check if any of them are outside [0,1]
sum ptaken_new


*****************
*	problem 5   *
*****************
/*
	Estimate the model by weighted least squares using Stata’s vwls 
	(variance weighted least squares) command. How much do the 
	coefficients differ from those obtained without weighting?
	How much do the standard errors differ from those produced by the 
	robust option?
*/

* estimate the variance of the dependent variable taken_new
egen s = sd(taken_new) 

* variance-weighted least squares, using s
eststo wls: vwls $vars , sd(s)

esttab using "results/q5.tex", replace ///
  	   mtitles("LPM" "LPM, robust SE" "vWLS") ///
	   title(Problem 5\label{q5}) ///
	   b(4) se(a2) /// 4 decimal places for coeffs & min 2 sig figs in se
	   nostar label booktabs nonumbers //formats it a little better


*****************
*	problem 6   *
*****************
/*
	Estimate probit and logit models of loan take-up using the same 
	independent variables as in Problem 2. Are the probit and logit 
	coefficient estimates similar to one another and to the LPM 
	estimates?
*/

eststo pro: probit $vars

eststo log: logit $vars

eststo drop wls 
esttab using "results/q6.tex", replace ///
	   mtitles("LPM" "LPM, robust SE" "Probit" "Logit") ///
	   title(Problem 6\label{q6}) ///
	   b(4) se(a2) /// 4 decimal places for coeffs & min 2 sig figs in se
	   nostar label booktabs nonumbers //formats it a little better
	   
	   
*****************
*	problem 7   *
*****************
/*
	Calculate the mean partial derivatives (a.k.a. marginal effects or 
	average partial effects) of the conditional probabilities of loan 
	take-up with respect to client age for the LPM, logit and probit
	models.
*/


* first doing all three with margins command
	* age is sort of continuous so don't need to worry about specifying 
		* factor variables
	* MARGINS MAY NOT WORK FOR LOGIT?
	* margins command is overkill for reg, but whatever. (because the derivative
		* is just equal to the coefficient, lol)
foreach mod in reg logit probit  {	
	* run the model:
	`mod' $vars
	* calculate the derivative using margins
	margins, dydx(client_age)
	quietly{ 
		matrix m = r(b)			//store dydx results in matrix
		local me_`mod' = m[1,1] //have to put in local alone first.. 
							//me = marginal effects
		global me_`mod' = "\texttt{`mod'}" + " & " + "`me_`mod''" + " \\" 
							//get that thang ready for a LaTex table
	}
}



*** doing extra for probit model (continuing problem 7) ***

* (a) Calculate analytic derivatives evaluated at the mean of the covariates.
eststo clear // clear the stored 
eststo dprob: dprobit $vars

sca deriv1 = round(e(dfdx)[1,1], .0000001) // store in scalar for table


* (b) Calculate mean analytic derivatives by hand. That is, use the formula to calculate the partial derivative for each observation and then summarize these to get the mean.
probit $vars // redo probit model
sca b_client_age = _b[client_age] //store the coefficent from probit est

predict probit_yhat //, xb //, index // get the predicted x*beta 
 
gen deriv2obs = normalden(probit_yhat)* b_client_age //follows formula

sum deriv2obs // get the mean over all variables 

sca deriv2 = round(_result(3), .0000001) // store in scalar for table


* (c) Calculate mean numerical derivatives by calculating predicted probabilities for each observation and then changing the client age value for each observation by a little bit (e.g. 0.001) and recalculating the predicted probabilities.

* small deviation in client age
sca dage = .0001 //even with .01 it was such a small change I got nothing...

replace client_age = client_age + dage

predict probit_yhat1 //, xb

* use the formula for derivative: 			
gen deriv3obs = (normal(probit_yhat1) - normal(probit_yhat))/dage
sum deriv3obs 
sca deriv3 = round(_result(3), .0000001)

* put client age back 
replace client_age = client_age - dage

* (d) Calculate mean derivatives using the margins command.
probit $vars

margins, dydx(client_age)
sca deriv4 = round(_result(3), .0000001)



*****************
*	problem 8   *
*****************
/*
	Re-estimate the LPM including a quartic in client age (i.e. including 
	linear, squared, cubed and fourth power terms), rather than just a 
	linear term. Calculate the average derivative numerically, as in part 
	(c), and compare it to that obtained from the probit model. Does the 
	LPM do better with greater flexibility?
*/

gen client_age2 = client_age^2
gen client_age3 = client_age^3
gen client_age4 = client_age^4

global vars2 "taken_new client_age* client_married client_education hh_size hh_income muslim hindu_sc_kat treated"

eststo clear
eststo lpm_quartic: reg $vars2 

* get the predicted values before making small change
predict lpm_xb0 

* change client age by a small amount 

replace client_age = client_age + dage
replace client_age2 = client_age^2
replace client_age3 = client_age^3
replace client_age4 = client_age^4

* predict new xb values
predict lpm_xb1

* calculate derivative at all observations 
gen lpm_deriv = (lpm_xb1 - lpm_xb0) / dage
qui sum lpm_deriv
sca lpm_deriv = _result(3)
di lpm_deriv 
di deriv3

* replace client age back to normal
replace client_age = client_age - dage



*****************
*	problem 9   *
*****************
/*
	For the probit model estimated in Problem 6, calculate the value of 
	the LRI “by hand” by obtaining the likelihood values for the model in 
	Problem 6 and a model with only an intercept and plugging them into 
	the LRI formula. Interpret the value you obtain.
*/

* repeat probit regression
probit $vars

* store scalars
sca ll = e(ll)
sca ll0 = e(ll_0)

* plug into LRI formula
sca lri = 1 - (ll / ll0)



*****************
*	problem 10   *
*****************
/*
	Calculate the correct prediction rate for the probit model estimated 
	in Problem 6 using both 0.5 and the sample fraction taking up a loan 
	as cutoff values. In each case, take the equally weighted average of 
	the correct prediction rates for those who take up a loan and those 
	who do not. Indicate and discuss how the model performs in each case. 
	Explain any differences between the two and indicate which measure of 
	goodness of fit you prefer and why. [Note: the equal weights refer to 
	the prediction rates, not the individual observations.]
*/

* get the sample fraction of taking up a loan (which is the mean)
sum taken_new
local loanfrac = _result(3)

* predicted probabilities are already done, in probit_yhat
pause on 
* create a variable for if predicted probabilities are > cutoff
forval i = 1/2{
	* define the cutoff rate
	if `i'== 1  local cutoff = 0.5
	if `i'== 2  local cutoff = `loanfrac'

	gen cutoff`i' = 0
	
	replace cutoff`i' =1 if probit_yhat> `cutoff'
	
	sum cutoff`i'
	sca cpr`i' = _result(3)
}

di cpr1
di cpr2

*****************
*	problem 11   *
*****************
/*
	Examine the difference between in-sample and out-of-sample predictive success by estimating the model using the usual covariates but only those observations with imidlineid < 1400. Compare the predictive performance of the model using the same criteria as in the preceding question for observations in and out of the estimation sample. What differences, if any, do you find? Does what you find accord with your expectations? Explain.
*/

*************
* write the tables
texdoc do metrics_ps1_tables
