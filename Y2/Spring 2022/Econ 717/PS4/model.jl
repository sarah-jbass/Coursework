# Sarah Bass
# Econ 717 - Problem Set 4 
# 4/12/2022

using Parameters, Distributions, Random, Optim, Statistics, KernelDensity, DataFrames, Plots

# parameters structure
mutable struct parameters
    π_1::Float64
    π_2::Float64
    μ_1::Float64
    μ_2::Float64
    σ_1::Float64
    σ_2::Float64
    ρ::Float64
end


function sim_data(θ::parameters; seed = 8, N = 1000)
    data = DataFrame(
        ϵ_1 = Float64[],
        ϵ_2 = Float64[],
        s_1 = Float64[],
        s_2 = Float64[],
        w_1 = Float64[],
        w_2 = Float64[],
        d = Int[],
        w = Float64[]
    )

    # covariance matrix
    σ_11 = θ.σ_1 * θ.σ_1
    σ_22 = θ.σ_2 * θ.σ_2
    σ_12 = θ.ρ * θ.σ_1 * θ.σ_2

    # create distribution
    ϵ_dist = MvNormal([0, 0], [σ_11 σ_12; σ_12 σ_22])

    # fill in our data frame with N=1000 simulated observations
    for i in 1:N
        # draw epsilon
        ϵ_1, ϵ_2 = rand(ϵ_dist)

        # compute skills
        s_1 = exp.(θ.μ_1 + ϵ_1)
        s_2 = exp.(θ.μ_2 + ϵ_2)

        # compute wages
        w_1 = θ.π_1 * s_1
        w_2 = θ.π_2 * s_2

        # occupation choice
        d = (w_1 >= w_2)

        # observed wages
        w = d * w_1 .+ (1 - d) * w_2 

        push!(data, [ϵ_1, ϵ_2, s_1, s_2, w_1, w_2, d, w])   
    end
    return data
end

function ll(b=[0.0,0.0]; θ::parameters=θ_0, df=data) # b is our guess of [μ_1,ρ]
    ll = 0 
    N = nrow(df)
    μ_1 = b[1]
    ρ = b[2]

    for i in 1:N
        # pull wage, epsilon, and d from the data (these are identified)
        w = data[i, :w] 
        # e2 = data[i, :ϵ_2]
        d = data[i, :d]

        # these are parts of the log likelihood function, broken up so that lli is cleaner below
        A = ((θ.σ_1 / θ.σ_2)*(log(w/θ.π_2) - θ.μ_2) - ρ*(log(w/θ.π_1) -μ_1)) / (θ.σ_1 * sqrt(1-ρ^2))
        B = (1/θ.σ_1)*(log(w/θ.π_1) - μ_1)
        C = ((θ.σ_2 / θ.σ_1)*(log(w/θ.π_1) - μ_1) - ρ*(log(w/θ.π_2) - θ.μ_2)) / (θ.σ_2 * sqrt(1-ρ^2))
        D = (1/θ.σ_2)*(log(w/θ.π_2) - θ.μ_2)

        # estimating our log likelihood for individual i, summing over all i
        lli = d*(log(cdf(Normal(), A)) + log(pdf(Normal(),B))) + (1-d)* (log(cdf(Normal(),C)) + log(pdf(Normal(),D))) 
        ll = ll + lli
    end
    ll = - ll
    return ll 
end


# # 2d identification plot for μ_1 and ρ
# function id_plot_3d(μ_hat::Float64, ρ_hat::Float64; θ::parameters=θ_0)

#     # grids
#     μ_1 = 1.0:0.05:1.25
#     ρ = 0.0:0.05:0.5
#     ll_hat = Array{Float64}(undef, length(μ_1), length(ρ))

#     for i in 1:length(μ_1)
#         for j in 1:length(ρ)

#         # function for plotting
#         ll_hat[i,j] =  ll([μ_1[i], ρ[j]])

#         end
#     end

#     # create plot
#     heatmap(μ_1, ρ, ll_hat)
#     scatter!([θ_0.μ_1], [θ_0.ρ], label = "(μ_1_0, ρ_0)")
#     scatter!([μ_hat], [ρ_hat], label = "(μ_1_hat, ρ_hat)")
#     xlabel!("μ_1")
#     ylabel!("ρ")
# end

# # identification plot for μ_1
# function id_plot_μ_1(μ_1_hat::Float64, ρ_hat::Float64, θ_0::parameters, d::data; method = "builtin")

#     # μ_1 grid
#     μ_1 = 1.0:0.01:1.5

#     # initialize object for log-likelihood vector
#     ll_μ_1_ρ_0 = zeros(length(μ_1))
#     ll_μ_1_ρ_hat = zeros(length(μ_1))

#     # evaluate log-likelihood across grid
#     for i in 1:length(μ_1)
#         ll_μ_1_ρ_0[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1[i], θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, θ_0.ρ), d; method=method)
#         ll_μ_1_ρ_hat[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1[i], θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_hat), d; method=method)
#     end

#     # plot log-likelihood
#     plot(μ_1, ll_μ_1_ρ_0, label = "Log Likelihood for ρ = ρ_0")
#     plot!(μ_1, ll_μ_1_ρ_hat, label = "Log Likelihood for ρ = ρ_hat")
#     vline!([μ_1_hat], label = "μ_1_hat")
#     vline!([θ_0.μ_1], label = "μ_1_0", legend = :bottomleft)
#     xlabel!("μ_1")
#     ylabel!("Log Likelihood")
# end

# # identification plot for ρ
# function id_plot_ρ(μ_1_hat::Float64, ρ_hat::Float64, θ_0::parameters, d::data; method = "builtin")

#     # ρ grid
#     ρ = 0.0:0.01:0.5

#     # initialize object for log-likelihood vector
#     ll_ρ_μ_1_0 = zeros(length(ρ))
#     ll_ρ_μ_1_hat = zeros(length(ρ))

#     # evaluate log-likelihood across grid
#     for i in 1:length(ρ)
#         ll_ρ_μ_1_0[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, θ_0.μ_1, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ[i]), d; method=method)
#         ll_ρ_μ_1_hat[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1_hat, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ[i]), d; method=method)
#     end

#     # plot log-likelihood
#     plot(ρ, ll_ρ_μ_1_0, label = "Log-likelihood for μ_1 = μ_1_0")
#     plot!(ρ, ll_ρ_μ_1_hat, label = "Log-likelihood for μ_1 = μ_1_hat")
#     vline!([ρ_hat], label = "ρ_hat")
#     vline!([θ_0.ρ], label = "ρ_0", legend = :bottomleft)
#     xlabel!("ρ")
#     ylabel!("Log Likelihood")
# end