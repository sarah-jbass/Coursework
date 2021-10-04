using Parameters, Plots #import the packages we want
include("PS3_functions.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model! @elapsed will time the functions, @time name_of_file.jl also works
@unpack val_func_e, val_func_u, pol_func_e, pol_func_u, μ = res  #separates the value and policy functions out of the results
@unpack a_grid = prim


##############Make plots
#value function
Plots.plot(a_grid, val_func_e, xlabel="K", ylabel="Value", label="Employed")
Plots.plot!(a_grid, val_func_u, label="Unemployed", title="Value Function")
Plots.savefig("02_Value_Function.png")

#policy functions
Plots.plot(a_grid, pol_func_e, xlabel="K", ylabel="K'", label="Employed")
Plots.plot!(a_grid, pol_func_u, label="Unemployed", title="Policy Function")
Plots.savefig("02_Policy_Function.png")

#changes in policy function
pol_func_δe = copy(pol_func_e).-a_grid
pol_func_δu = copy(pol_func_u).-a_grid
Plots.plot(a_grid, pol_func_δe, xlabel="K", ylabel="K'-K", label="Employed")
Plots.plot!(a_grid, pol_func_δu, label="Unemployed", title="Policy Function Changes")
Plots.savefig("02_Policy_Functions_Changes.png")

println("All done!")
################################
