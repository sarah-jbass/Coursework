#keyword-enabled structure to hold model primitives
@with_kw mutable struct Primitives #structure with all of our constants in it
    β::Float64 = 0.9932 #discount rate, Float64 used for decimals
    q::Float64 = 0.9932
    α::Float64 = 1.5 #capital share
    ee::Int64 = 1 #employment earnings
    eu::Float64 = 0.5 #unemployment earnings
    pee::Float64 = 0.97
    peu::Float64 = 0.03
    pue::Float64 = 0.5
    puu::Float64 = 0.5
    a_min::Float64 = -2 #capital lower bound
    a_max::Float64 = 5 #capital upper bound
    na::Int64 = 1000 #number of capital grid points, use Int64 for integers
    a_grid::Array{Float64,1} = collect(range(a_min, length = na, stop = a_max))
    a_pmf::Array{Float64,2} = fill(1/(na*2), (1000,2)) #original distribution
end

#structure that holds model results
mutable struct Results #preallocating the results, Results is the type of structure
    val_func_e::Array{Float64, 1} #value function, there is one long vector of unspecified length
    val_func_u::Array{Float64, 1}
    pol_func_e::Array{Float64, 1} #policy function
    pol_func_u::Array{Float64, 1}
    pol_func_i_e::Array{Int64, 1}
    pol_func_i_u::Array{Int64, 1}
    pmf_prime::Array{Float64,2}
    μ::Array{Float64, 2}
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func_e = zeros(prim.na) #initial value function guess, na is the number of grid points
    val_func_u = zeros(prim.na)
    pol_func_e = zeros(prim.na) #initial policy function guess
    pol_func_u = zeros(prim.na)
    pol_func_i_e = zeros(prim.na)
    pol_func_i_u = zeros(prim.na)
    pmf_prime = zeros(prim.na,2)
    res = Results(val_func_e, val_func_u, pol_func_e, pol_func_u, pol_func_i_e, pol_func_i_u, pmf_prime) #initialize results struct
    prim, res #return deliverables
end

#Bellman Operator
function Bellman(prim::Primitives,res::Results)
    @unpack val_func_e, val_func_u = res #unpack value functions
    @unpack a_grid, β, α, na, ee, eu, pee, peu, pue, puu, q = prim #unpack model primitives
    v_next_e = zeros(na) #next guess of value function to fill, don't need to use prim.na since we unpacked prim
    v_next_u = zeros(na)

    for i_e = 1:2
        choice_lower = 1
        if i_e == 1
            s = ee
            pe = pee
            pu = peu
        else
            s = eu
            pe = pue
            pu = puu
        end

        #This loop solve for the optimal level of a' for each possible initial value of a
        for a_index = 1:na
            a = a_grid[a_index] #value of a
            candidate_max = -Inf #bad candidate max

            for ap_index in 1:na#choice_lower:na #loop over possible selections of k', exploiting monotonicity of policy function
                ap = a_grid[ap_index]
                c = s + a - q*ap #consumption
                if c>0 #check for positivity
                    u = (c^(1-α)-1)/(1-α) #consumption given a' selection
                    val = u + β*(pe*val_func_e[ap_index] + pu*val_func_u[ap_index]) #compute value
                    if val>candidate_max #check for new max value
                        candidate_max = val #update max value
                        if i_e == 1
                            res.pol_func_e[a_index] = a_grid[ap_index] #update policy function, not unpacked because then we'd have to repack
                            res.pol_func_i_e[a_index] = ap_index #a_grid index that corresponds
                            choice_lower = ap_index #update lowest possible choice
                        else
                            res.pol_func_u[a_index] = a_grid[ap_index] #update policy function, not unpacked because then we'd have to repack
                            res.pol_func_i_u[a_index] = ap_index
                            choice_lower = ap_index #update lowest possible choice
                        end
                    end
                end
            end
            if i_e == 1
                v_next_e[a_index] = candidate_max #update value function
            else
                v_next_u[a_index] = candidate_max
            end
        end
        #if z==z_min
            #v_next_b #return next guess of value function
        #else
            #v_next_g
        #end
    end
    v_next_e, v_next_u
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err_e::Float64 = 100.0, err_u::Float64 = 100.0)
    n = 0 #counter

    while err_e>tol || err_u>tol#begin iteration
        v_next_e, v_next_u = Bellman(prim, res) #spit out new vectors, next guess

        err_e = abs.(maximum(v_next_e.-res.val_func_e))/abs(v_next_e[prim.na, 1]) #reset error level, initial res.val_func is just 0
        err_u = abs.(maximum(v_next_u.-res.val_func_u))/abs(v_next_u[prim.na, 1])

        res.val_func_e = v_next_e
        res.val_func_u = v_next_u #update value function
        n+=1
    end
    println("Value function converged in ", n, " iterations.")

    #T* operator - finding stationary distribution
    pmf = prim.a_pmf
    #pmf_prime = zeros(1000,2) #is this syntax wrong?
    pmf_converged = 0

    if pmf_converged == 0
        res.pmf_prime = zeros(prim.na,2)
         for i_a = 1:prim.na
            i_ap_e = res.pol_func_i_e[i_a]
            i_ap_u = res.pol_func_i_u[i_a]
            res.pmf_prime[i_ap_e,1] = res.pmf_prime[i_ap_e,1] + pmf[i_a,1]*prim.pee
            res.pmf_prime[i_ap_e,2] = res.pmf_prime[i_ap_e,2] + pmf[i_a,1]*prim.peu
            res.pmf_prime[i_ap_u,1] = res.pmf_prime[i_ap_u,1] + pmf[i_a,2]*prim.pue
            res.pmf_prime[i_ap_u,2] = res.pmf_prime[i_ap_u,2] + pmf[i_a,2]*prim.puu
        end
        if sum(res.pmf_prime) > 1+10^(-9)
            println("Warning: pmf sums to ", sum(res.pmf_prime), '.')
        end
        pmf_diff = maximum(abs.(res.pmf_prime .- pmf))
        if pmf_diff<tol
            pmf_converged=1
        else
            pmf = res.pmf_prime
        end
    end
    prim.a_pmf = pmf
end

#solve the model
function Solve_model(prim::Primitives, res::Results)
    tol::Float64 = 1e-2
    tol2::Float64 = 10e-6
    q_min::Float64 = prim.β
    q_max::Float64 = 1.0
    mc::Float64 = 1.0
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
    mc = sum(res.pol_func_e .* res.pmf_prime[:,1] + res.pol_func_u .* res.pmf_prime[:,2])
    mc1::Float64 = 0
    mc2::Float64 = 0
    mc3::Float64 = 0
    while abs(mc) > tol && (q_max - q_min)/2 > tol2
        #Find 3 q values to test
        q1 = (3q_min + q_max)/4
        q2 = (q_min + q_max)/2
        q3 = (q_min + 3q_max)/4
        #Find 3 asset clearing values
        for i = 1:3
            if i == 1
                prim.q = q1
                V_iterate(prim, res)
                mc1 = sum(res.pol_func_e .* res.pmf_prime[:,1] + res.pol_func_u .* res.pmf_prime[:,2])
                println("MC1 is ", mc1)
            elseif i == 2
                prim.q = q2
                V_iterate(prim, res)
                mc2 = sum(res.pol_func_e .* res.pmf_prime[:,1] + res.pol_func_u .* res.pmf_prime[:,2])
                println("MC2 is ", mc2)
            else
                prim.q = q3
                V_iterate(prim, res)
                mc3 = sum(res.pol_func_e .* res.pmf_prime[:,1] + res.pol_func_u .* res.pmf_prime[:,2])
                println("MC3 is ", mc3)
            end
        end
        #Replace mc value, q, q_min, q_max with new values
        if abs(mc1) - minimum([abs(mc1) abs(mc2) abs(mc3)]) < tol2
            mc = mc1
            prim.q = q1
            q_max = q2
        elseif  abs(mc2) - minimum([abs(mc1) abs(mc2) abs(mc3)]) < tol2
            mc = mc2
            prim.q = q2
            q_min = q1
            q_max = q3
        else abs(mc3) - minimum([abs(mc1) abs(mc2) abs(mc3)]) < tol2
            mc = mc3
            prim.q = q3
            q_min = q2
        end
        println("Weighted sum of assets is ", mc)
    end
end

function calculate_wealth_distribution(results::Results)
    @unpack a_grid, ee, eu, na = Primitives()

    w = zeros(2, na)

    for i_s = 1:2
        for i_a = 1:na

            i_w = argmin(abs.(a_grid[i_a] .+ S[i_s] .- a_grid))

            w[i_w, i_s] = results.μ[i_a, i_s]
        end
    end
    w
end

function calculate_lorenz_curve(w::Array{Float64, 2})
    @unpack a_grid, na = Primitives()

    x = cumsum(w[:,1] .+ w[:,2])
    y = cumsum((w[:,1] .+ w[:,2]) .* a_grid)

    unique([x/x[a_length] y/y[a_length]]; dims = 1)
end

# https://en.wikipedia.org/wiki/Gini_coefficient
function calculate_gini(l::Array{Float64, 2})
    widths = diff(l[:,1])
    heights = ((l[1:end-1,1] .+ l[2:end,1])./2 .- (l[1:end-1,2] .+ l[2:end,2])./2)
    a = sum(widths .* heights)

    l_pos = l[l[:,2].>0, :]
    widths = diff(l_pos[:,1])
    heights = (l_pos[1:end-1,2] .+ l_pos[2:end,2])./2
    b = sum(widths .* heights)

    a/(a+b)
end

##############################################################################
