/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 717
Task: PS2
Date created: 2/28/22
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS2\Econ717_PS2_Bass.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all
set more off
set maxvar 32000
macro drop _all
pause on

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Spring 2022\Econ 717\PS2"

texdoc init "q6_table1.tex", replace

/*tex
\begin{table}[h!]
\caption{Problem 7}
	\label{q7a}
\begin{center}
\begin{tabular}{rc}
%\multicolumn{2}{c}{Marginal effects of the conditional probabilities of loan 
%	take-up with respect to client age} \\
\toprule
Estimation method & derivative (using margins) \\
\hline
tex*/

local row1 = 