using Parameters, Plots, CSV, Tables, DataFrames
include("PS3_functions.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
Solve_HH_Problem(prim, res)
Solve_μ(prim, res)
# @elapsed Solve_model(prim, res) #solve the model! @elapsed will time the functions, @time name_of_file.jl also works
@unpack value_function, policy_function, labor_supply, μ = res  #separates the value and policy functions out of the results
@unpack a_grid = prim

##############Make plots
#value function
Plots.plot(a_grid, value_function[50,:,1], xlabel="Assets", ylabel="Value", title="Value Function")
Plots.savefig("Value_Function.png")

#policy functions
Plots.plot(a_grid, policy_function[20,:,:], xlabel="Assets", ylabel="New Assets", title="Policy Function")
Plots.savefig("Policy_Function.png")

#savings decision
policy_function_δ = copy(policy_function[20,:,:]).-a_grid
Plots.plot(a_grid, policy_function_δ[20,:,:], xlabel="Assets", ylabel="Savings", title="Savings Decision")
Plots.savefig("Labor_Supply.png")

#labor supply
Plots.plot(a_grid, labor_supply[20,:,:], xlabel="Assets", ylabel="Labor Supply", title="Labor Supply")
Plots.savefig("Labor_Supply.png")

################################

# Benchmark
@elapsed bm_ss = Solve_model()  # converges in ~9 iterations
@elapsed bm_no_ss = Solve_model(θ = 0.0)  # converges in ~11 iterations

# No productivity risk
@elapsed riskless_ss = Solve_model(z = [0.5, 0.5]) # converges in ~12 iterations
@elapsed riskless_no_ss = Solve_model(θ = 0.0, z = [0.5, 0.5], λ = 0.1)  # converges in ~52 iterations

# Inelastic labor supply
@elapsed inelastic_l_ss = Solve_model(γ = 1.0, λ = 0.8) # converges in ~6 iterations
@elapsed inelastic_l_no_ss = Solve_model(θ = 0.0, γ = 1.0, λ = 0.8) # converges in ~7 iterations

# table 1
table_1 = create_table([bm_ss, bm_no_ss,
                        riskless_ss, riskless_no_ss,
                        inelastic_l_ss, inelastic_l_no_ss])

CSV.write("table_1.csv", table_1)

# Figure 4
bm_ss.value_function[isinf.(bm_ss.value_function)] .= -1/eps();
bm_no_ss.value_function[isinf.(bm_no_ss.value_function)] .= -1/eps();

bm_ss_line = reshape(sum(bm_ss.value_function .* bm_ss.μ; dims = [2, 3]), 66, 1)
bm_no_ss_line = reshape(sum(bm_no_ss.value_function .* bm_no_ss.μ; dims = [2, 3]), 66, 1)

plot([bm_ss_line bm_no_ss_line],
     legend=:bottomright,
     labels = ["SS" "No SS"],
     title="Figure 6")

savefig("value_bm_ss_no_ss.png")

sum(bm_ss_line .< bm_no_ss_line)


# Figure 5
bm_ss_line = reshape(sum(bm_ss.value_function .* bm_ss.μ; dims = [1, 3]), 5000, 1)

bm_no_ss_line = reshape(sum(bm_no_ss.value_function .* bm_no_ss.μ; dims = [1, 3]), 5000, 1)

plot([bm_ss_line[2:5000] bm_no_ss_line[2:5000]],
     legend=:bottomright,
     labels = ["SS" "No SS"],
     title="Figure 4")

savefig("value_bm_ss_no_ss.png")

sum(bm_ss_line .< bm_no_ss_line)
