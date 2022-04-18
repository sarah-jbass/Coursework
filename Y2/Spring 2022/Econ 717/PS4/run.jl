rootpath = pwd()
outpath = joinpath(pwd(),"tex files/")

using Parameters, Distributions, Random, Optim, Statistics, KernelDensity, DataFrames
using Latexify, LaTeXStrings

# set parameters
θ_0 = parameters(1.0, 1.0, 1.15, 1.0, 0.4, 0.5, 0.25)

#simulate data
data = sim_data(θ_0)

# estimate log likelihood
b_0 = [1.15, 0.25] # = [μ_1, ρ] initial estimates
res = optimize(ll, b_0)

# output results to latex
mu_hat = latexify(res.minimizer[1]; fmt = "%2.3f")
rho_hat = latexify(res.minimizer[2]; fmt = "%2.3f")

# open(outpath*"e.p_estimates.tex", "a+") do file 
#     write(outpath*"e.p_estimates.tex", L"Estimating the parameters with MLE, I find $\hat\mu_1 =$"*mu_hat*L" and $\hat\rho =$"*rho_hat)