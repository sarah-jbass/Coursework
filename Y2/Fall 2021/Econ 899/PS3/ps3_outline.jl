#=
PS3 Outline

1. Set primitives and results structs
    a. Asset values: min (0) max (30)
    b. 3D matrices:
        Value function: age (66) x asset grid (1000) x productivity (2)
        Policy function for assets: age (66) x asset grid (1000) x productivity (2)
        Policy function for labor: age (66) x labor grid (1000) x productivity (2)
        Distribution μ: age (66) x asset grid (1000) x productivity (2)
2. Use backwards induction for the value function
    a. In last period you consume everything. Know value at period N in life.
    b. Knowing the value corresponding to each a' in the final period of life, solve for max value function to determine a'.
    c. Solve for labor value that corresponds to that a'.
3. Calculate stationary distribution
    a. Pretty much the same but now there is age too. Remember population size is decreasing in age.
4. Calculate market clearing
    a. Pause at some capital and labor demand, from the firms problem that implies prices, check labor/capital supply from household, then rerun value function iteration

Other comments:
1. Can reshape 3D matrix to get a 2D slice
2. 
