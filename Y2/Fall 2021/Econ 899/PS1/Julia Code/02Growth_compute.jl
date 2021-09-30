using Parameters, Plots #import the packages we want
include("02Growth_model.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model! @elapsed will time the functions, @time name_of_file.jl also works
@unpack val_func_g, val_func_b, pol_func_g, pol_func_b = res  #separates the value and policy functions out of the results
@unpack k_grid = prim

##############Make plots
#value function
Plots.plot(k_grid, val_func_g, xlabel="K", ylabel="Value", label="High Productivity")
Plots.plot!(k_grid, val_func_b, label="Low Productivity", title="Value Function")
Plots.savefig("02_Value_Function.png")

#policy functions
Plots.plot(k_grid, pol_func_g, xlabel="K", ylabel="K'", label="High Productivity")
Plots.plot!(k_grid, pol_func_b, label="Low Productivity", title="Policy Function")
Plots.savefig("02_Policy_Function.png")

#changes in policy function
pol_func_δg = copy(pol_func_g).-k_grid
pol_func_δb = copy(pol_func_b).-k_grid
Plots.plot(k_grid, pol_func_δg, xlabel="K", ylabel="K'-K", label="High Productivity")
Plots.plot!(k_grid, pol_func_δb, label="Low Productivity", title="Policy Function Changes")
Plots.savefig("02_Policy_Functions_Changes.png")

println("All done!")
################################
