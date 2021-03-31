/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 710B
Task: PS8
Date created: 3/29/21
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS8\PS8_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Bruce Hansen\Problem Sets\PS8"

*Question 19.9
use "Invest1993.dta"

gen I = inva
gen Q = vala

*Nadaraya-Watson
lpoly I Q if Q<=5, ci noscatter title("Nadaraya-Watson 95% CI")
graph export "9-1.png", as(png) name("Graph") replace

*Local Linear
lpoly I Q if Q<=5, degree(1) ci noscatter title("Local Linear 95% CI")
graph export "9-2.png", as(png) name("Graph") replace

*Question 19.11
clear
use "FRED-QD.dta"

gen y = 100*((gdpc1/l.gdpc1)^4-1)
gen ylag = l.y

*Nadaraya-Watson
lpoly y ylag, ci noscatter title("Nadaraya-Watson 95% CI")
graph export "11-1.png", as(png) name("Graph") replace

*Local Linear
lpoly y ylag, degree(1) ci noscatter title("Local Linear 95% CI")
graph export "11-2.png", as(png) name("Graph") replace

clear all