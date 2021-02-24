clear all
set more off
macro drop _all
pause on
capture log close

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 709\Jack Porter\Problem Sets\PS3"

/*-------------------------------------------------------------------------------------------------------------------------------
Econ 709, Second Quarter, Problem Set 3
Sarah Bass
11/24/20

do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 709\Jack Porter\Problem Sets\PS3\PS3_Bass.do"
--------------------------------------------------------------------------------------------------------------------------------*/

* Load the data
use cps09mar.dta

* Generate transformations
gen wage = ln(earnings/(hours*week))
gen experience = age - education - 6
gen exp2 = (experience^2)/100

* Create indicator for subsamples
gen mbf = (race == 2) & (marital <= 2) & (female == 1) // married black female
gen mbf12 = (mbf == 1) & (experience == 12) // mbf with 12 yrs exp
gen sam = (race == 4) & (marital == 7) & (female == 0) //single asian male

keep if sam==1 & experience<45

*3.24a) Estimate equation 3.50 and compute R2 and sum of squared errors
    reg wage education experience exp2
    outreg2 using table1b.tex, replace ct("log(wage)") tex(fragment) lab adds(Sum-of-squared Errors, e(rss))

*3.24b) Re-estimate the slope on education using the residual regression approach
    *ln(wage) on experience
    reg wage experience exp2
    predict e_wage, residual 	// save residuals from wage-experience regression 
    outreg2 using table1b.tex, append ct("log(wage)") tex(fragment) lab adds(Sum-of-squared Errors, e(rss))

    *education on experience
    reg education experience exp2
    predict e_educ, residual 	// save residuals from education-experience regression 
    outreg2 using table1b.tex, append ct("education") tex(fragment) lab	adds(Sum-of-squared Errors, e(rss))

    *wage residuals on education residuals 
    reg e_wage e_educ 
    outreg2 using table1b.tex, append ct("$\hat{\varepsilon}_{wage}$") tex(fragment) lab adds(Sum-of-squared Errors, e(rss))

*3.25)
    reg wage education experience exp2 
    predict yhat, xb
    predict ehat, residual

    *3.25a)
    quietly summ ehat
    local a = round(r(sum), .001)

    *3.25b)
    gen x1_ehat = education * ehat
    quietly summ x1_ehat
    local b = round(r(sum), .001)

    *3.25c)
    gen x2_ehat = experience * ehat
    quietly summ x2_ehat
    local c = round(r(sum), .001)

    *3.25d)
    gen x1sq_ehat = (education^2) * ehat
    quietly summ x1sq_ehat
    local d = round(r(sum), .001)

    *3.25e)
    gen x2sq_ehat = (experience^2) * ehat
    quietly summ x2sq_ehat
    local e = round(r(sum), .001)

    *3.25f)
    gen yhat_ehat = yhat * ehat
    quietly summ yhat_ehat
    local f = round(r(sum), .001)

    *3.25g)
    gen ehat_sq = ehat^2
    quietly summ ehat_sq
    local g = round(r(sum), .001)

    file open resultsfile using "3.25_results.tex", write replace
    file write resultsfile															///
        "\begin{itemize}" 												_newline	///
            _tab "\item[(a)] $\sum_{i=1}^n\hat{e}_i = `a'$"				_newline	///
            _tab "\item[(b)] $\sum_{i=1}^nX_{1i}\hat{e}_i = `b'$"		_newline	///
            _tab "\item[(c)] $\sum_{i=1}^nX_{2i}\hat{e}_i = `c'$"		_newline	///
            _tab "\item[(d)] $\sum_{i=1}^nX_{1i}^2\hat{e}_i = `d'$"		_newline	///
            _tab "\item[(e)] $\sum_{i=1}^nX_{2i}^2\hat{e}_i = `e'$"		_newline	///
            _tab "\item[(f)] $\sum_{i=1}^n\hat{Y}_i\hat{e}_i = `f'$"	_newline	///
            _tab "\item[(g)] $\sum_{i=1}^n\hat{e}^2_i = `g'$"			_newline	///
        "\end{itemize}"
    file close resultsfile
