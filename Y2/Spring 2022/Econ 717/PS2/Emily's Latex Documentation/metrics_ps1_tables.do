/*
	applied metrics spring 2022
	problem set 1 tables
	emily case 
*/

quietly{
*************************
*	problem 7: tables   *
*************************

* initialize tex doc for a table comparing ME of LPM, probit, and logit: 
texdoc init "results/q7a.tex", replace

* write the table heading:
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

* write the table contents:
	texdoc write $me_reg
	texdoc write $me_probit
	texdoc write $me_logit

* write the table footing:
/*tex
\bottomrule
\end{tabular}
\end{center}
\end{table}
tex*/

* finish this table/close the tex doc
texdoc close 


* table for parts a through d
texdoc init "results/q7b.tex", replace

* write the table heading:
/*tex
\begin{table}[h!]
\caption{Problem 7}
	\label{q7b}
\begin{center}
\begin{tabular}{lll}
\multicolumn{3}{c}{Probit model derivatives, 4 ways} \\
\toprule
& method & derivative \\
\hline
tex*/
local row1 = "(a) & \texttt{dprobit} & "
local row2 = "(b) & using formula \& \texttt{summarize} & "
local row3 = "(c) &  numerical derivatives & "
local row4 = "(d) & \texttt{margins} & "


forval i = 1/4  {	
	local deriv = round(deriv`i', .0000001)
	local row`i' = "`row`i''" + "`deriv'" + "\\" //perhaps an issue with sig figs
	texdoc write `row`i''
}

* write the table footing:
/*tex
\bottomrule
\end{tabular}
\end{center}
\end{table}
tex*/
texdoc close 


*************************
*	problem 8: tables   *
*************************

texdoc init "results/q8.tex", replace
/*tex
\begin{table}[h!]
\caption{Problem 8: comparing numerical derivatives}
	\label{q8}
\begin{center}
\begin{tabular}{cc}
\toprule
Probit & LPM (with quartic client age) \\
\toprule
tex*/
local probitderiv = round(deriv3, .0000001)
local lpmderiv = round(lpm_deriv, .0000001)

local row1 = "`probitderiv'" + "&" + "`lpmderiv'" + "\\"
macro list

texdoc write `row1' 

/*tex
\bottomrule
\end{tabular}
\end{center}
\end{table}
tex*/
texdoc close 

}
