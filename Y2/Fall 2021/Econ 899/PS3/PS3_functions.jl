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
    η::Array{Float64, 1} = map(x->parse(Float64,x), readlines(open("ef.txt")))
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
    policy_function_i::Array{Int64}
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
    μ = ones(prim.N, prim.na, prim.nz) ./ sum(ones(prim.N, prim.na, prim.nz))

    # Initial values based on representative agent model
    k = 0.0
    l = 0.0
    w = 1.05
    r = 0.05
    b = 0.2

    res = Results(θ, z, e, value_function, policy_function, policy_function_i, labor_supply, μ, k, l, w, r, b)
    prim, res #return deliverables
end

function Retired_Utility(a::Float64, ap::Float64, prim::Primitives, res::Results)
    @unpack r, b = res #unpack value functions
    @unpack γ, σ = prim
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

    for ia=1:na
        res.value_function[N,ia,:] .= Retired_Utility(a_grid[ia], 0.0, prim, res)
    end

    for j=(N-1):-1:Jᴿ

        for ia = 1:na
            a = a_grid[ia] #value of a
            candidate_max = -Inf #bad candidate max

            for iap in 1:na#choice_lower:na #loop over possible selections of k', exploiting monotonicity of policy function
                ap = a_grid[iap]
                u = Retired_Utility(a, ap, prim, res) #consumption given a' selection
                val = u + β*res.value_function[j+1,iap,1] #compute value
                if val>candidate_max #check for new max value
                    candidate_max = val #update max value
                    res.value_function[j,ia,:] .= val
                    res.policy_function[j,ia,:] .= a_grid[iap] #update policy function
                    res.policy_function_i[j,ia,:] .= iap #a_grid index that corresponds
                end
            end
        end
    end
    #Productivity doesn't matter in retirement so these are all the same
end

function Labor_Supply(γ::Float64, θ::Float64, e::Float64, w::Float64, r::Float64, a::Float64, ap::Float64)
    interior_solution = (γ*(1-θ)*e*w - (1-γ)*((1+r)*a - ap)) / ((1-θ)*w*e)
    min(1, max(0, interior_solution))
end

function Working_Utility(l::Float64, a::Float64, ap::Float64, e::Float64, prim::Primitives, res::Results)
    @unpack θ, w, r, b = res #unpack value functions
    @unpack γ, σ  = prim

    c = w*(1-θ)*e*l + (1+r)*a - ap
    if (c>0)
        (((c^γ) * ((1-l)^(1-γ))) ^ (1-σ)) / (1 - σ)
    else
        -Inf
    end
end

function Working_Bellman(prim::Primitives, res::Results)
    @unpack θ, e, w, r, b = res #unpack value functions
    @unpack Π, γ, σ, β, N, Jᴿ, nz, a_min, a_max, na, a_grid  = prim

    for j=(Jᴿ-1):-1:1

        for z = 1:nz

            for ia = 1:na
                a = a_grid[ia] #value of a
                candidate_max = -Inf #bad candidate max

                for iap in 1:na#choice_lower:na #loop over possible selections of k', exploiting monotonicity of policy function
                    ap = a_grid[iap]
                    l = Labor_Supply(γ, θ, e[j,z], w, r, a, ap)
                    if l<0 || l>1
                        println("INVALID LABOR SUPPLY")
                    end
                    u = Working_Utility(l, a, ap, e[j,z], prim, res)
                    val = u + β*(Π[z,1]*res.value_function[j+1,iap,1] + Π[z,2]*res.value_function[j+1,iap,2])#compute value
                    if val>candidate_max #check for new max value
                        candidate_max = val #update max value
                        res.value_function[j,ia,z] = val
                        res.policy_function[j,ia,z] = a_grid[iap] #update policy function
                        res.policy_function_i[j,ia,z] = iap #a_grid index that corresponds
                        res.labor_supply[j,ia,z] = l
                    end
                end
            end
        end
    end
end

function Solve_HH_Problem(prim::Primitives, res::Results)
    Retired_Bellman(prim, res)
    Working_Bellman(prim, res)
end

# function Solve_μ(prim::Primitives, res::Results)
#     @unpack N, n, na, nz, Π, Π₀ = prim
#     tol = 1e-5
#     init_μ = res.μ
#     converged = 0
#
#     res.μ = zeros(N, na, nz) #set PMF=0
#     res.μ[1, 1, :] = Π₀
#
#     if converged == 0
#         for j = 1:(N-1)
#             for ia = 1:na
#                 for z = 1:nz
#                     iap = res.policy_function_i[j,ia,z]
#                     for zp = 1:nz
#                         res.μ[j+1,iap,zp] = res.μ[j+1,iap,zp] + Π[z,zp]*init_μ[j,ia,z]
#                     end
#                 end
#             end
#         end
#
#         age_weights = ones(N)
#
#         for  i = 1:(N-1)
#             age_weights[i+1] = age_weights[i]/(1+n)
#         end
#
#         age_weight = reshape(repeat(age_weights/sum(age_weights), na * nz), N, na, nz)
#         res.μ = age_weight .* res.μ
#
#         μ_diff = maximum(abs.(res.μ .- init_μ))
#         if μ_diff<tol
#             converged=1
#         else
#             init_μ = res.μ
#         end
#     end
# end

function Solve_μ(prim::Primitives, res::Results; progress::Bool = false)
    @unpack a_grid, N, n, nz, Π₀, Π, na = Primitives()

    # sets distribution to zero.
    res.μ = zeros(N, na, nz)

    # Fills in model-age 1 with erodgic distribution of producitivities.
    res.μ[1, 1, :] = Π₀

    for j = 1:(N-1) # Iterates through model-ages
        if progress
            println(j)
        end
        for i_a in 1:na # Iterates through asset levels
            for i_z = 1:nz
                if res.μ[j, i_a, i_z] == 0 # skips if no mass at j, i_a, i_z
                    continue
                end
                # finds index of assets tomorrow
                i_a_p = argmax(a_grid .== res.policy_function[j, i_a, i_z])
                for i_z_p = 1:nz # iterates over productivity levels tomorrow
                    res.μ[j+1, i_a_p, i_z_p] += Π[i_z, i_z_p] * res.μ[j, i_a, i_z]
                end
            end
        end
    end

    # sum(res.μ) should equal N at this point (i.e. 1 for each row)
    # Now adjusts for population growth, so that sum(res.μ) = 1

    age_weights_temp = ones(N)

    for i = 1:(N-1)
        age_weights_temp[i + 1] = age_weights_temp[i]/(1+n)
    end

    age_weight = reshape(repeat(age_weights_temp/sum(age_weights_temp), na * nz), N, na, nz)

    res.μ = age_weight .* res.μ
end

function Update_prices(prim::Primitives, res::Results, k::Float64, l::Float64)
    @unpack α, δ, Jᴿ, N = prim
    res.w = (1 - α) * k ^ α * l ^ (- α)
    res.r = α * k ^ (α - 1) * l ^ (1-α) - δ
    res.b = (res.θ * res.w * l) / sum(res.μ[Jᴿ:N, :, :])
end

function Calculate_labor_supply(prim::Primitives, res::Results)
    @unpack Jᴿ, na, nz = prim
    e_3d = reshape(repeat(res.e, na), Jᴿ -1, na, nz)
    sum(res.μ[1:(Jᴿ - 1),:,:] .* res.labor_supply .* e_3d)
end

function Calculate_capital_supply(prim::Primitives, res::Results)
    @unpack a_grid, N, nz, na, N = prim
    a_grid_3d = permutedims(reshape(repeat(a_grid, N * nz), na, N, nz), (2, 1, 3))
    sum(res.μ .* a_grid_3d)
end

function Solve_model(prim::Primitives, res::Results; θ::Float64 = 0.11, z::Array{Float64, 1} = [3.0, 0.5], γ::Float64 = 0.42, λ::Float64 = 0.5)
    # @unpack γ = prim
    # @unpack θ = res

    k_0 = 3.3
    l_0 = 0.3

    Update_prices(prim, res, k_0, l_0)

    ε = 0.001  # tolerance
    i = 0      # counter
    diff = 1

    while diff > ε
        i += 1
        println("Iteration #", i)
        println("Capital demand: ", k_0)
        println("Labor demand: ", l_0)

        Solve_HH_Problem(prim, res)
        Solve_μ(prim, res)

        k_1 = Calculate_capital_supply(prim, res)
        l_1 = Calculate_labor_supply(prim, res)

        println("Capital supply: ", k_1)
        println("Labor supply: ", l_1)

        diff = abs(k_0 - k_1) + abs(l_0 - l_1)

        println("Absolute difference: ", diff)
        println("************************************")

        k_0 = λ * k_1 + (1 - λ) * k_0
        l_0 = λ * l_1 + (1 - λ) * l_0
        Update_prices(prim, res, k_0, l_0)
    end
    res.k = k_0
    res.l = l_0
    res
end
#
# function Process_results(res::Results)
#     # calculate total welfare
#     welfare = res.value_function .* res.μ
#     welfare = sum(welfare[isfinite.(welfare)])
#
#     # calculate coefficient of variation of wealth
#     @unpack a_grid, N, nz, na, N = Primitives()
#     a_grid_3d = permutedims(reshape(repeat(a_grid, N * nz), na, N, nz), (2, 1, 3))
#     wealth_mean = sum(res.μ .* a_grid_3d)
#     wealth_second_moment = sum(res.μ .* a_grid_3d .^ 2)
#     wealth_second_central_moment = wealth_second_moment - wealth_mean^2
#     cv = wealth_mean / sqrt(wealth_second_central_moment)
#
#     # create vector of summary statistics
#     [res.θ, res.γ, res.z[1], res.k, res.l, res.w, res.r, res.b, welfare, cv]
# end
#
# function Create_table(results_vector::Array{Results})
#     table = DataFrames.DataFrame(Tables.table(reduce(hcat,process_res.(results_vector))'))
#     rename!(table, [:theta, :gamma, :z_H, :k, :l, :w, :r, :b, :welfare, :cv])
# end
