// code for questions 2 and 3 of hw2

// clear workspace
clear

// 2a (setup)

// for each city, number of firms uniformly distributed on {1,2,...,10}
set obs 1000
gen unif = runiform()
gen num_firms = int(unif*10+1)

// in 500 cities, firms collude perfectly when num_firms <= 8
gen collude = 0
replace collude = 1 if _n <= 500 & num_firms <= 8

// 2b (construct L_i, H, e)

// demand function parameter initialization
gen c0 = 1
gen c1 = 0.9
gen xi = 0

// cost function parameter initialization
gen F = 1
gen b0 = 1
gen b1 = 0
gen eta = 0

// construct Lerner index, Herfindahl index, demand elasticity
gen lerner_cournot = c1/num_firms
gen lerner_monopoly = c1
gen lerner = collude*lerner_monopoly + (1-collude)*lerner_cournot
gen observed_lerner = ln(lerner) + 0.1*(unif - 0.5)
gen herfindahl = 1/num_firms
gen elasticity = 1/c1

// 2c (regressions and tests)

// structure-conduct-performance paradigm regressions
gen ln_herfindahl = ln(herfindahl)
regress observed_lerner ln_herfindahl if _n <= 500  // collusion is possible
test ln_herfindahl = 1
regress observed_lerner ln_herfindahl if _n > 500  // no collusion
test ln_herfindahl = 1
regress observed_lerner ln_herfindahl  // pooled sample
test ln_herfindahl = 1

// 2d (repeat 2b and 2c for linear demand)

// demand function parameter initialization
gen a0 = 3
gen a1 = 1
gen nu = 0

// construct Lerner index, Herfindahl index, demand elasticity
gen lerner_cournot2 = (a0+nu-b0-eta)/(a0+nu+num_firms*(b0+eta))
gen lerner_monopoly2 = (a0+nu-b0-eta)/(a0+nu+b0+eta)
gen lerner2 = collude*lerner_monopoly2 + (1-collude)*lerner_cournot2
gen observed_lerner2 = ln(lerner2) + 0.1*(unif - 0.5)
gen herfindahl2 = 1/num_firms
gen elasticity2 = (a0+nu+b0+eta)/(a0+nu-b0-eta)

// structure-conduct-performance paradigm regressions
gen ln_herfindahl2 = ln(herfindahl2)
regress observed_lerner2 ln_herfindahl2 if _n <= 500  // collusion is possible
test ln_herfindahl2 = 1
regress observed_lerner2 ln_herfindahl2 if _n > 500  // no collusion
test ln_herfindahl2 = 1
regress observed_lerner2 ln_herfindahl2  // pooled sample
test ln_herfindahl2 = 1

// 3a (setup, regressions, construct num_firms, L_i, H)

// clear workspace
clear

// 1000 cities
set obs 1000
gen unif = runiform()

// demand function parameter initialization
gen a0 = 5
gen a1 = 1
gen nu = 2*(unif - 0.5)

// cost function parameter initialization
gen F = 1
gen b0 = 1
gen b1 = 0
gen eta = 0

// firms enter until profits are zero
gen num_firms = (a0+nu-b0-eta-sqrt(F*a1))/(sqrt(F*a1))

// construct Lerner index, Herfindahl index, demand elasticity
gen lerner = (a0+nu-b0-eta)/(a0+nu+num_firms*(b0+eta))
gen observed_lerner = ln(lerner) + 0.1*(unif - 0.5)
gen herfindahl = 1/num_firms
gen elasticity = (a0+nu+b0+eta)/(a0+nu-b0-eta)

// structure-conduct-performance paradigm regression
gen ln_herfindahl = ln(herfindahl)
regress observed_lerner ln_herfindahl

// 3b (repeat 3a for new eta and nu)

// parameter initialization
drop nu
drop eta
gen nu = 0
gen eta = 2*(unif - 0.5)

// firms enter until profits are zero
gen num_firms2 = (a0+nu-b0-eta-sqrt(F*a1))/(sqrt(F*a1))

// construct Lerner index, Herfindahl index, demand elasticity
gen lerner2 = (a0+nu-b0-eta)/(a0+nu+num_firms2*(b0+eta))
gen observed_lerner2 = ln(lerner2) + 0.1*(unif - 0.5)
gen herfindahl2 = 1/num_firms2
gen elasticity2 = (a0+nu+b0+eta)/(a0+nu-b0-eta)

// structure-conduct-performance paradigm regression
gen ln_herfindahl2 = ln(herfindahl2)
regress observed_lerner2 ln_herfindahl2
