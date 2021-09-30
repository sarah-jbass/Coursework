#keyword-enabled structure to hold model primitives
@with_kw struct Primitives #structure with all of our constants in it
    β::Float64 = 0.99 #discount rate, Float64 used for decimals
    δ::Float64 = 0.025 #depreciation rate
    α::Float64 = 0.36 #capital share
    k_min::Float64 = 0.01 #capital lower bound
    k_max::Float64 = 75.0 #capital upper bound
    z_min::Float64 = 0.2
    z_max::Float64 = 1.25
    pgg::Float64 = 0.977
    pgb::Float64 = 0.023
    pbg::Float64 = 0.074
    pbb::Float64 = 0.926
    nk::Int64 = 1000 #number of capital grid points, use Int64 for integers
    k_grid::Array{Float64,1} = collect(range(k_min, length = nk, stop = k_max)) #capital grid
end


#structure that holds model results
mutable struct Results #preallocating the results, Results is the type of structure
    val_func_g::Array{Float64, 1} #value function, there is one long vector of unspecified length
    val_func_b::Array{Float64, 1}
    pol_func_g::Array{Float64, 1} #policy function
    pol_func_b::Array{Float64, 1}
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func_g = zeros(prim.nk) #initial value function guess, nk is the number of grid points
    val_func_b = zeros(prim.nk)
    pol_func_g = zeros(prim.nk) #initial policy function guess
    pol_func_b = zeros(prim.nk)
    res = Results(val_func_g, val_func_b, pol_func_g, pol_func_b) #initialize results struct
    prim, res #return deliverables
end

#Bellman Operator
function Bellman(prim::Primitives,res::Results)
    @unpack val_func_g, val_func_b = res #unpack value functions
    @unpack k_grid, β, δ, α, nk, z_min, z_max, pgg, pgb, pbg, pbb = prim #unpack model primitives
    v_next_g = zeros(nk) #next guess of value function to fill, don't need to use prim.nk since we unpacked prim
    v_next_b = zeros(nk)

    #for exploiting monotonicity of policy function
    for i_z = 1:2
        choice_lower = 1
        if i_z == 1
            z = z_min
            pg = pbg
            pb = pbb
        else
            z = z_max
            pg = pgg
            pb = pgb
        end

        for k_index = 1:nk
            k = k_grid[k_index] #value of k (capital)
            candidate_max = -Inf #bad candidate max
            budget = z*k^α + (1-δ)*k #budget

            for kp_index in 1:nk#choice_lower:nk #loop over possible selections of k', exploiting monotonicity of policy function
                c = budget - k_grid[kp_index] #consumption given k' selection
                if c>0 #check for positivity
                    val = log(c) + β*(pg*val_func_g[kp_index] + pb*val_func_b[kp_index]) #compute value
                    if val>candidate_max #check for new max value
                        candidate_max = val #update max value
                        if i_z == 1
                            res.pol_func_b[k_index] = k_grid[kp_index] #update policy function, not unpacked because then we'd have to repack
                            choice_lower = kp_index #update lowest possible choice
                        else
                            res.pol_func_g[k_index] = k_grid[kp_index] #update policy function, not unpacked because then we'd have to repack
                            choice_lower = kp_index #update lowest possible choice
                        end
                    end
                end
            end
            if i_z == 1
                v_next_b[k_index] = candidate_max #update value function
            else
                v_next_g[k_index] = candidate_max
            end
        end
        #if z==z_min
            #v_next_b #return next guess of value function
        #else
            #v_next_g
        #end
    end
    v_next_b, v_next_g
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err_g::Float64 = 100.0, err_b::Float64 = 100.0)
    n = 0 #counter

    while err_g>tol || err_b>tol#begin iteration
        v_next_b, v_next_g = Bellman(prim, res) #spit out new vectors, next guess

        err_g = abs.(maximum(v_next_g.-res.val_func_g))/abs(v_next_g[prim.nk, 1]) #reset error level, initial res.val_func is just 0
        err_b = abs.(maximum(v_next_b.-res.val_func_b))/abs(v_next_b[prim.nk, 1])

        res.val_func_g = v_next_g
        res.val_func_b = v_next_b #update value function
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end

#solve the model
function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
end
##############################################################################
