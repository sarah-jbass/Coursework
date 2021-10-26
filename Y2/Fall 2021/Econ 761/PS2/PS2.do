/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 761
Task: IO Problem Set 2
Date created: 10/11/21
Programmer: Sarah Bass
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Fall 2021\Econ 761\PS2\PS2.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Fall 2021\Econ 761\PS2"

clear all
set more off
set obs 1000

set seed 543186
gen unif = uniform()

set seed 155133
gen unif2 = uniform() - 0.5

gen Num = int(unif*10+1)
gen collude = (_n <= 500 & N <= 8)

*/
gen F = 1
gen c0 = 1
gen c1 = .9
gen b0 = 1
gen b1 = 0
gen z = 0
gen h = 0
gen xi = 0 
gen a0 = 3
gen a1 = 1 
gen nu = 0 
gen eta = 0

file open table1 using "table1.tex", write replace

foreach eqn in 3 1{
	file write table1 	///
		_newline _tab "&&&&& \\ \textit{Equation (`eqn')} & & & & & \\"
	
	if (`eqn' == 3){
		gen LernerCN`eqn' = c1/Num
		gen LernerM`eqn' = c1
	}
	if (`eqn' == 1){
		gen LernerCN`eqn' = (a0+nu-b0-eta)/(a0+nu+Num*(b0+eta))
		gen LernerM`eqn' = (a0+nu-b0-eta)/(a0+nu+b0+eta)
	}

 	gen Herf`eqn' 	= 1/Num
	gen lnHerf`eqn' = ln(Herf`eqn')
	 
	gen Lerner`eqn' = collude * LernerM`eqn' + (1-collude)*LernerCN`eqn'
	gen lnLernerObs`eqn' = ln(Lerner`eqn') + .1*unif2
	
	// generate list of sample descriptions
	loc samp_descs Cournot Collusion Pooled
	
	// run analysis both pooled and separate
	forval pooled = 0/2{
		
		// save sample description
		loc samp `: word `=`pooled'+1' of `samp_descs''

		// SCP regression
		if (`pooled' == 2){
			reg lnLernerObs`eqn' lnHerf`eqn'
			loc N = e(N)
		}
		else if(`pooled' == 1) {
			reg lnLernerObs`eqn' lnHerf`eqn' if collude == 1
			loc N = e(N)
		}
        else if(`pooled' == 0){
            reg lnLernerObs`eqn' lnHerf`eqn' if collude == 0
			loc N = e(N)
        }
	
		loc beta : 	di %4.3f _b[lnHerf`eqn']
		loc se : 	di %4.3f _se[lnHerf`eqn']
		
		test _b[lnHerf`eqn'] = 1
		
		loc F : di %7.0fc r(F)
		loc p : di %4.3f r(p)
		
		file write table1 	///
			_newline _tab "`samp' & `beta' & `se' & `F' & `p' & `N' \\"

	} // pool loop
	
}
file close table1