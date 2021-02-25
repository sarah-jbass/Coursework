/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 710A
Task: Calculate IV Estimates
Date created: 2/3/21
Programmer: Sarah Bass
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Spring 2021\Econ 710\Mikkel Solvsten\Problem Sets\PS2\Econ 714A PS2 Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Spring 2021\Econ 710\Mikkel Solvsten\Problem Sets\PS2"
use "AE80.dta"

local control agem agefstm boy1st boy2nd blackm hispm othracem
local ys workedm weeksm1 hourswm
local ym workedd weeksd1 hourswd 

*Reduced form
reg morekids samesex `control'

foreach y in `ys'{
    *Reg 1
    reg `y' morekids `control' // col 1
	matrix row=r(table)
	local b1_`y'=row[1,1]
	local se1_`y'=row[2,1]

    *Reg 2
	ivregress 2sls `y' (morekids=samesex) `control' // col 2
	matrix row=r(table)
	local b2_`y'=row[1,1]
	local se2_`y'=row[2,1]

    *Reg 3
	ivregress 2sls `y' (morekids=samesex) `control' if msample==1 // col 5
	matrix row=r(table)
	local b3_`y'=row[1,1]
	local se3_`y'=row[2,1]

    *Rounding
	local b1_`y' = round(`b1_`y'', .001)
	local b2_`y' = round(`b2_`y'', .001)
	local b3_`y' = round(`b3_`y'', .001)
	local se1_`y' = round(`se1_`y'', .001)
	local se2_`y' = round(`se2_`y'', .001)
	local se3_`y' = round(`se3_`y'', .001)
}

foreach y in `ym'{
    *Reg 1
    reg `y' morekids `control' if msample==1 // col 7 
	matrix row=r(table)
	local b4_`y'=row[1,1]
	local se4_`y'=row[2,1]

    *Reg 2
	ivregress 2sls `y' (morekids=samesex) `control' if msample ==1 // col 8
	matrix row=r(table)
	local b5_`y'=row[1,1]
	local se5_`y'=row[2,1]

    *Rounding
	local b4_`y' = round(`b4_`y'', .001)
	local b5_`y' = round(`b5_`y'', .001)
	local se4_`y' = round(`se4_`y'', .001)
	local se5_`y' = round(`se5_`y'', .001)
}


file open resultsfile using "ps2_results.tex", write replace
    file write resultsfile                                                          ///
        "\begin{tabular}{c | c c c c c}"                                               _newline    ///
		            "\hline"             _newline    ///
            "& All W., OLS & All W., 2SLS & M. W., 2SLS & H. of M.W., OLS & H. of M.W., 2SLS \\"             _newline    ///
		            "\hline"             _newline    ///
            "Worked for pay & `b1_workedm' & `b2_workedm' & `b3_workedm' & `b4_workedd' & `b5_workedd' \\"             _newline    ///
						  " & (`se1_workedm') & (`se2_workedm') & (`se3_workedm') & (`se4_workedd') & (`se5_workedd')\\"       _newline    ///
            "Weeks worked & `b1_weeksm1' & `b2_weeksm1' & `b3_weeksm1' & `b4_weeksd1' & `b5_weeksd1' \\"             _newline    ///
						  " & (`se1_weeksm1') & (`se2_weeksm1') & (`se3_weeksm1') & (`se4_weeksd1') & (`se5_weeksd1') \\"       _newline    ///
            "Hours per week & `b1_hourswm' & `b2_hourswm' & `b3_hourswm' & `b4_hourswd' & `b5_hourswd' \\"             _newline    ///
						  " & (`se1_hourswm') & (`se2_hourswm') & (`se3_hourswm') & (`se4_hourswd') & (`se5_hourswd')"       _newline    ///
        "\end{tabular}"
    file close resultsfile