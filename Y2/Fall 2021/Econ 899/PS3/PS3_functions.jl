using Parameters, CSV, Tables, DataFrames

#keyword-enabled structure to hold model primitives
@with_kw struct Primitives #structure with all of our constants in it
    # Life cycle params
    N::Int64 = 66
    n::Float64 = 0.011
    Jᴿ::Int64 = 46

    #Preference params
    σ::Float64 = 2.0
    γ::Float64 = 0.42
    β::Float64 = 0.97

    #Productivity params
    η::Array{Float64, 1} = map(x->parse(Float64,x), readlines("ef.txt"))
    nz::Int64 = 2
    Π::Array{Float64,2} = [0.9261 0.0739; 0.0189 0.9811]
    Π₀::Array{Float64,2} = [0.2037 0.7963]

    #Production params
    α::Float64 = 0.36
    δ::Float64 = 0.06

    #Asset grid
    a_min::Float64 = 0 #capital lower bound
    a_max::Float64 = 30 #capital upper bound
    na::Int64 = 1000 #number of capital grid points, use Int64 for integers
    a_grid::Array{Float64,1} = collect(range(a_min, length = na, stop = a_max))
end

#structure that holds model results
mutable struct Results #preallocating the results, Results is the type of structure
    θ::Float64
    z::Array{Float64, 1}
    e::Array{Float64, 2}
    value_function::Array{Float64}
    policy_function::Array{Float64}
    policy_function_i::Array{Float64}
    labor_supply::Array{Float64}
    μ::Array{Float64}
    k::Float64
    l::Float64
    w::Float64
    r::Float64
    b::Float64
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives

    θ = 0.11
    z = [3.0; 0.5] #2x1
    e = prim.η * z'#45 x 2

    # 1st dim is age, 2nd dim is asset holding, 3rd dim is productivity state
    value_function = zeros(prim.N, prim.na, prim.nz) #66x1000x2
    policy_function = zeros(prim.N, prim.na, prim.nz)
    policy_function_i = zeros(prim.N, prim.na, prim.nz)
    labor_supply = zeros(prim.Jᴿ-1, prim.na, prim.nz) #45x1000x2
    μ = ones(prim.N, prim.na, prim.nz) / sum(ones(prim.N, prim.na, prim.nz))

    # Initial values based on representative agent model
    k = 0.0
    l = 0.0
    w = 1.05
    r = 0.05
    b = 0.2

    res = Results(θ, z, e, value_function, policy_function, policy_function_i, labor_supply, μ, k, l, w, r, b)
    prim, res #return deliverables
end

function Retired_Utility(a::Float64, b::Float64, ap::Float64, r::Float64, γ::Float64, σ::Float64)
    c = (1+r)a + b - ap
    if (c>0)
        (c^((1-σ) * γ))/(1 - σ)
    else
        -Inf
    end
end

function Retired_Bellman(prim::Primitives, res::Results)
    @unpack r, b = res #unpack value functions
    @unpack γ, σ, β, N, Jᴿ, nz, a_min, a_max, na, a_grid  = prim
    res.value_function[N,:,1] = Retired_Utility.(a_grid, b, 0, r, γ, σ)

    for j=(N-1):-1:Jᴿ

        for a = 1:na
            a = a_grid[a] #value of a
            candidate_max = -Inf #bad candidate max

            for ap in 1:na#choice_lower:na #loop over possible selections of k', exploiting monotonicity of policy function
                ap = a_grid[ap]
                u = Retired_Utility(a, b, ap, r, γ, σ) #consumption given a' selection
                val = u + β*res.value_function[j+1,ap,1] #compute value
                if val>candidate_max #check for new max value
                    candidate_max = val #update max value
                    res.value_function[j,a,:] = val
                    res.policy_function[j,a,:] = a_grid[ap] #update policy function
                    res.policy_function_i[j,a,:] = ap #a_grid index that corresponds
                end
            end
        end
    end
    #Productivity doesn't matter in retirement so these are all the same
end

function Labor_Supply(γ::Float64, θ::Float64, e::Float64, w::Float64, r::Float64, a::Float64, ap::Float64)
    l = (γ*(1-θ)*e*w - (1-γ)*((1+r)*a - a_p)) / ((1-θ)*w*e)
    if l<0 || l>1
        println("INVALID LABOR SUPPLY")
    end
end

function Working_Utility(l::Float64, a::Float64, b::Float64, ap::Float64, e::Float64, r::Float64, w::Float64, θ::Float64, γ::Float64, σ::Float64)
    c = w*(1-θ)*e*l + (1+r)*a - ap
    if (c>0)
        (((c^γ) * ((1-l)^(1-γ))) ^ (1-σ)) / (1 - σ)
    else
        -Inf
    end
end

function Working_Bellman(prim::Primitives, res::Results)
    @unpack θ, e, w, r, b = res #unpack value functions
    @unpack γ, σ, β, N, Jᴿ, nz, a_min, a_max, na, a_grid  = prim

    for j=(Jᴿ-1):-1:1

        for z = 1:nz

            for a = 1:na
                a = a_grid[a] #value of a
                candidate_max = -Inf #bad candidate max

                for ap in 1:na#choice_lower:na #loop over possible selections of k', exploiting monotonicity of policy function
                    ap = a_grid[ap]
                    l = Labor_Supply(γ, θ, e[j,z], w, r, a, ap)
                    u = Working_Utility(l, a, b, ap, e[j,z], r, w, θ, γ, σ)
                    val = u + β*(Π[z,1]*res.value_function[j+1,ap,1] + Π[z,2]*res.value_function[j+1,ap,2])#compute value
                    if val>candidate_max #check for new max value
                        candidate_max = val #update max value
                        res.value_function[j,a,z] = val
                        res.policy_function[j,a,z] = a_grid[ap] #update policy function
                        res.policy_function_i[j,a,z] = ap #a_grid index that corresponds

                    end
                end
            end
        end
    end
end

function Solve_HH_problem(prim::Primitives, res::Results)
    Retired_Bellman(prim, res)
    Working_Bellman(prim, res)
end

function Solve_μ(prim::Primitives, res::Results)
    @unpack N, n, na, nz = prim
    tol = 1e-5
    init_μ = res.μ
    converged = 0

    if converged == 0
        res.μ = zeros(N, na, nz) #set PMF=0
        res.μ[1, 1, :] = Π₀
        for j = 1:(N-1)
            for a = 1:na
                for z = 1:nz
                ap = res.policy_function_i[j,a,z]
                for zp = 1:nz
                    res.μ[j+1,ap,zp] = res.μ[j+1,ap,zp] + prim.Π[z,zp]*init_μ[j,a,z]
                end
            end
            μ_diff = maximum(abs.(res.μ .- init_μ))
            if μ_diff<tol
                converged=1
            else
                init_μ = res.μ
            end
        end
    end

    age_weights = ones(N)

    for  i = 1:(N-1)
        age_weights[i+1] = age_weights[i]/(1+n)
    end

    age_weight = reshape(repeat(age_weights/sum(age_weights), na * nz), N, na, nz)
    res.μ = age_weight .* res.μ
end

function Update_prices(prim::Primitives, res::Results, k::Float64, l::Float64)
    @unpack α, δ, Jᴿ, N = prim
    res.w = (1 - α) * k ^ α * l ^ (- α)
    res.r = α * k ^ (α - 1) * l ^ (1-α) - δ
    res.b = (res.θ * res.w * l) / sum(res.μ[Jᴿ:N, :, :])
end

# #Value function iteration
# function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err_e::Float64 = 100.0, err_u::Float64 = 100.0)
#     n = 0 #counter
#
#     while err_e>tol || err_u>tol#begin iteration
#         v_next_e, v_next_u = Bellman(prim, res) #spit out new vectors, next guess
#
#         err_e = abs.(maximum(v_next_e.-res.val_func_e))/abs(v_next_e[prim.na, 1]) #reset error level, initial res.val_func is just 0
#         err_u = abs.(maximum(v_next_u.-res.val_func_u))/abs(v_next_u[prim.na, 1])
#
#         res.val_func_e = v_next_e
#         res.val_func_u = v_next_u #update value function
#         n+=1
#     end
#     println("Value function converged in ", n, " iterations.")
#
#     #T* operator - finding stationary distribution
#     pmf = prim.a_pmf
#     #pmf_prime = zeros(1000,2) #is this syntax wrong?
#     pmf_converged = 0
#
#     if pmf_converged == 0
#         res.pmf_prime = zeros(prim.na,2)
#          for i_a = 1:prim.na
#             i_ap_e = res.pol_func_i_e[i_a]
#             i_ap_u = res.pol_func_i_u[i_a]
#             res.pmf_prime[i_ap_e,1] = res.pmf_prime[i_ap_e,1] + pmf[i_a,1]*prim.pee
#             res.pmf_prime[i_ap_e,2] = res.pmf_prime[i_ap_e,2] + pmf[i_a,1]*prim.peu
#             res.pmf_prime[i_ap_u,1] = res.pmf_prime[i_ap_u,1] + pmf[i_a,2]*prim.pue
#             res.pmf_prime[i_ap_u,2] = res.pmf_prime[i_ap_u,2] + pmf[i_a,2]*prim.puu
#         end
#         if sum(res.pmf_prime) > 1+10^(-9)
#             println("Warning: pmf sums to ", sum(res.pmf_prime), '.')
#         end
#         pmf_diff = maximum(abs.(res.pmf_prime .- pmf))
#         if pmf_diff<tol
#             pmf_converged=1
#         else
#             pmf = res.pmf_prime
#         end
#     end
#     prim.a_pmf = pmf
# end
#
# #solve the model
# function Solve_model(prim::Primitives, res::Results)
#     tol::Float64 = 1e-2
#     tol2::Float64 = 10e-6
#     q_min::Float64 = prim.β
#     q_max::Float64 = 1.0
#     mc::Float64 = 1.0
#     V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
#     mc = sum(res.pol_func_e .* res.pmf_prime[:,1] + res.pol_func_u .* res.pmf_prime[:,2])
#     mc1::Float64 = 0
#     mc2::Float64 = 0
#     mc3::Float64 = 0
#     while abs(mc) > tol && (q_max - q_min)/2 > tol2
#         #Find 3 q values to test
#         q1 = (3q_min + q_max)/4
#         q2 = (q_min + q_max)/2
#         q3 = (q_min + 3q_max)/4
#         #Find 3 asset clearing values
#         for i = 1:3
#             if i == 1
#                 prim.q = q1
#                 V_iterate(prim, res)
#                 mc1 = sum(res.pol_func_e .* res.pmf_prime[:,1] + res.pol_func_u .* res.pmf_prime[:,2])
#                 println("MC1 is ", mc1)
#             elseif i == 2
#                 prim.q = q2
#                 V_iterate(prim, res)
#                 mc2 = sum(res.pol_func_e .* res.pmf_prime[:,1] + res.pol_func_u .* res.pmf_prime[:,2])
#                 println("MC2 is ", mc2)
#             else
#                 prim.q = q3
#                 V_iterate(prim, res)
#                 mc3 = sum(res.pol_func_e .* res.pmf_prime[:,1] + res.pol_func_u .* res.pmf_prime[:,2])
#                 println("MC3 is ", mc3)
#             end
#         end
#         #Replace mc value, q, q_min, q_max with new values
#         if abs(mc1) - minimum([abs(mc1) abs(mc2) abs(mc3)]) < tol2
#             mc = mc1
#             prim.q = q1
#             q_max = q2
#         elseif  abs(mc2) - minimum([abs(mc1) abs(mc2) abs(mc3)]) < tol2
#             mc = mc2
#             prim.q = q2
#             q_min = q1
#             q_max = q3
#         else abs(mc3) - minimum([abs(mc1) abs(mc2) abs(mc3)]) < tol2
#             mc = mc3
#             prim.q = q3
#             q_min = q2
#         end
#         println("Weighted sum of assets is ", mc)
#     end
# end
#
# function calculate_wealth_distribution(results::Results)
#     @unpack a_grid, ee, eu, na = Primitives()
#
#     w = zeros(2, na)
#
#     for i_s = 1:2
#         for i_a = 1:na
#
#             i_w = argmin(abs.(a_grid[i_a] .+ S[i_s] .- a_grid))
#
#             w[i_w, i_s] = results.μ[i_a, i_s]
#         end
#     end
#     w
# end
#
# function calculate_lorenz_curve(w::Array{Float64, 2})
#     @unpack a_grid, na = Primitives()
#
#     x = cumsum(w[:,1] .+ w[:,2])
#     y = cumsum((w[:,1] .+ w[:,2]) .* a_grid)
#
#     unique([x/x[a_length] y/y[a_length]]; dims = 1)
# end
#
# # https://en.wikipedia.org/wiki/Gini_coefficient
# function calculate_gini(l::Array{Float64, 2})
#     widths = diff(l[:,1])
#     heights = ((l[1:end-1,1] .+ l[2:end,1])./2 .- (l[1:end-1,2] .+ l[2:end,2])./2)
#     a = sum(widths .* heights)
#
#     l_pos = l[l[:,2].>0, :]
#     widths = diff(l_pos[:,1])
#     heights = (l_pos[1:end-1,2] .+ l_pos[2:end,2])./2
#     b = sum(widths .* heights)
#
#     a/(a+b)
# end
#
# ##############################################################################
