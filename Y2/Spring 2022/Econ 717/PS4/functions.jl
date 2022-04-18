

using Parameters, Distributions, Random, Optim, Statistics, KernelDensity, DataFrames
using Latexify


# create empty parameters 
mutable struct parameters      # mutable means we can go in and change them  
    π_1::Float64               # later
    π_2::Float64
    μ_1::Float64               # Float64 allows for decimals, have to say what
    μ_2::Float64               # type of parameters these are
    σ_1::Float64
    σ_2::Float64
    ρ::Float64
end



function simdata(θ::parameters; seed = 8, N = 1000)
    # this function (1) creates an empty dataframe to be populated, (2) simulates some data given parameters and (3) estimates individual likelihoods 

    # create an empty dataframe for simulated data
        data = DataFrame(
            ϵ_1 = Float64[],
            ϵ_2 = Float64[],
            s_1 = Float64[],
            s_2 = Float64[],
            w_1 = Float64[],
            w_2 = Float64[],
            d   = Int[],
            w   = Float64[]
        )

    # create ϵ-distribution

        # covariance matrix
        σ_11 = θ.σ_1 * θ.σ_1
        σ_22 = θ.σ_2 * θ.σ_2
        σ_12 = θ.ρ * θ.σ_1 * θ.σ_2

        ϵ_dist = MvNormal([0, 0], [σ_11 σ_12; σ_12 σ_22])
        ϵ1_d   = Normal(0, σ_11)
        ϵ2_d   = Normal(0, σ_22)

    for i in 1:N
        e1, e2 = rand(ϵ_dist)

        # create the skill variables 
        s1 = exp.(θ.μ_1 + e1)
        s2 = exp.(θ.μ_2 + e2)

        # create wages
        w1 = θ.π_1 * s1
        w2 = θ.π_2 * s2

        # dummy for if occ 1 is chosen 
        d = w1>=w2

        # observed wages
        w = d*w1 + (1-d)*w2

        # push these numbers to an observation in the dataframe 
        push!(data,[e1,e2,s1,s2,w1,w2,d,w])
    end
    return data
end

function ll(θ::parameters, df)
    ll = 0 
    N = nrow(df)
    for i in 1:N
        w = df[i, :w] 
        e2 = df[i, :ϵ_2]
        d = df[i, :d]

        A  = ((θ.σ_1 / θ.σ_2)*e2 - θ.ρ*(log(w/θ.π_1) -θ.μ_1)) / (θ.σ_1 * sqrt(1-θ.ρ^2))
        B = (1/θ.σ_1)*(log(w/θ.π_1) - θ.μ_1)
        C = ((θ.σ_2 / θ.σ_1)*(log(w/θ.π_1) -θ.μ_1) - θ.ρ*(e2)) / (θ.σ_2 * sqrt(1-θ.ρ^2))
        D = e2/θ.σ_2

        lli = d*(log(cdf(Normal(), A)) + log(pdf(Normal(),B))) + (1-d)* (log(cdf(Normal(),C)))
        ll = ll +lli
    end
return ll 
end


