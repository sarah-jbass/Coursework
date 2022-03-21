/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 717
Task: PS3
Date created: 3/15/22
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS3\Econ717_PS3_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
capture log close
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS3"
log using "Econ717_PS3_Bass.log", replace

use "Economics 717 Miron and Tetelbaum Data.dta"
rename _all, lower

*Question 1
    xtset state year

*Question 2
    gen mlda21 = (mlda==21) 

*Question 3
    reg rate18_20ht mlda21, robust
        outreg2 using table1.tex, replace `opts' ctitle(Q3)	///
		addtext(Year FE, No, State FE, No, Cluster, No) keep(mlda21)

*Question 4 
    *State FE
    reg rate18_20ht mlda21 i.state, robust
        outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q4a)	///
		addtext(Year FE, Yes, State FE, No, Cluster, No)

    *Year FE
    reg rate18_20ht mlda21 i.year, robust
        outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q4b)	///
		addtext(Year FE, No, State FE, Yes, Cluster, No)

*Question 5
    *State and year FE, state cluster
    reg rate18_20ht mlda21 i.state i.year, cluster(state) robust
        outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q5)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes)

*Question 6
    *State and year FE
    reg rate18_20ht mlda21 i.state i.year, robust
        outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q6)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, No)

*Question 7
    *State and year FE, state cluster
    reg rate18_20ht mlda21 i.state i.year if year<=1990, cluster(state) robust
        outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q7)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes)


*Question 8
    *Placebo treatment indicator
    gen placebo82 = (mldayr == 1987 & year >= 1982) 
    egen min_mlda = min(mlda), by(state)

    preserve    
    keep if year <= 1987 & (mlda == 21 | mldayr == 1987)

    *State and year FE, state cluster
    reg rate18_20ht placebo82 i.state i.year, cluster(state) robust
        outreg2 using table1.tex, append `opts' ctitle(Q8)	sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) keep(placebo82)

    restore 

*Question 9

	reg rate18_20ht mlda21 i.state i.year if state == 21 | min_mlda == 21, robust cl(state)
		outreg2 using table2.tex, replace `opts' ctitle(MD Only) sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) keep(mlda21)

	reg rate18_20ht mlda21 i.state i.year if state == 23 | min_mlda == 21, robust cl(state)
		outreg2 using table2.tex, append `opts' ctitle(MI Only) sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) keep(mlda21)

*Question 10
	*Generate separate treatment indicators
	bys state mlda (year): 	gen mlda21_14 = (year - year[1] <= 3 & mlda == 21 & min_mlda == 18)
	bys state mlda (year): 	gen mlda_later = (year - year[1] > 3 & mlda == 21 & min_mlda == 18)

	
	*un analysis
	reg rate18_20ht mlda21_14 mlda_later i.state i.year if state == 21 | min_mlda == 21, robust cl(state)
		outreg2 using table2.tex, append `opts' ctitle(MD Only) sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) 	///
		keep(mlda21_14 mlda_later)

	reg rate18_20ht mlda21_14 mlda_later i.state i.year if state == 23 | min_mlda == 21, robust cl(state)
		outreg2 using table2.tex, append `opts' ctitle(MI Only) sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) 	///
		keep(mlda21_14 mlda_later)

*Bonus question
    *Goodman-Bacon decomposition
	bacondecomp rate18_20ht mlda21, ddetail nograph cl(state) robust
				
	*Save results from decomposition and load them
	foreach x in b V dd wt sumdd{
        mat `x' = e(`x')
    }
    
    mat dat = J(`=colsof(dd)', 2, .)
    mat rown dat = `: coln dd'
    mat coln dat = dd wt

    forval r = 1/`=colsof(dd)'{
        mat dat[`r', 1] = dd[1, `r']
        mat dat[`r', 2] = wt[1, `r']
    }
	
	clear
	svmat dat, n(col)
	gen category = ""
	forval i = 1/`=colsof(dd)'{
		replace category = "`: word `i' of `: coln dd''" if _n == `i'
	}
	gen b = `=b[1,1]'
	gen se = sqrt(`=V[1,1]')
	 
    *Comparison groups from category variable
	gen first = substr(category, 1, strpos(category, "_") - 1)
	gen last  = substr(category, strpos(category, "_") + 1, .)
	gen	group = "Treated vs. Always 21" if last == "Always"
	
	destring first, g(fyear) force
	destring last, g(lyear) force
	
	replace group = "Pre vs. Post" if fyear < lyear & !missing(lyear)
	replace group = "Post vs. Pre" if fyear > lyear & !missing(lyear)
	
	*Plot the estimates
	loc note "Dashed lines indicate the 95% confidence interval"
	loc note "`note' for the TWFEDD estimate"
	tw	///
		scatter dd wt if group == "Treated vs. Always 21", m(t)	||	///
		scatter dd wt if group == "Pre vs. Post", m(oh)			||	///
		scatter dd wt if group == "Post vs. Pre", m(x)				///
			graphregion(color(white)) ylab(-15(5)15, angle(horizontal))		///
			xlab(0(0.1)0.6) title("2x2 DiD Estimates vs. TWFEDD Weight")	///
			xtitle("Weight") ytitle("") note("`note'")				///
			yline(`=b[1]', lc(red)) yline(0, lw(thin) lc(black))	///
			yline(`=b[1] + 1.96*se[1]', lp(dash) lc(red) lw(vthin))	///
			yline(`=b[1] - 1.96*se[1]', lp(dash) lc(red) lw(vthin))	///
			text(`=b[1] + se[1]' 0.3 								///
				"TWFEDD Estimate = `=round(b[1], .01)'", 			///
				color(red) size(small))								///
			leg(lab(1 "Treated vs. Always 21") lab(2 "Pre vs. Post") 	///
				nobox r(1) lab(3 "Post vs. Pre") region(color(white)))
				
	graph export figure1.png, replace

log close