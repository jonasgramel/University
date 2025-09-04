using Pkg
# Pkg.add("Plots")
# Pkg.add("ProgressMeter")
# Pkg.add("Statistics")
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("GLM)
# Pkg.add("DelimitedFiles")
# Pkg.add("LaTeXStrings")
# Pkg.add("StatsBase")

using Plots
using ProgressMeter
using Statistics
using DataFrames
using CSV
using GLM
using DelimitedFiles
using LaTeXStrings
using StatsBase

println("Starting plotmaker 2000...")
df_3_1_1 = CSV.read("DataFrames/H_zero/p_zero/Problem_3_1_1.csv", DataFrame)

df_3_1_1_E = plot(df_3_1_1.E_T21, label="T=2.1", xlabel="Steps", ylabel="Energy, E",legend=:topright, grid=true, fontsize=40, guidefont=font(14), legendfont=14)
plot!(df_3_1_1_E, df_3_1_1.E_T23, label="T=2.3")
plot!(df_3_1_1_E, df_3_1_1.E_T25, label="T=2.5")
savefig(df_3_1_1_E, "figures/Problem_3_1_1_E.png")


df_3_1_1_M = plot(df_3_1_1.M_T21, label="T=2.1", xlabel="Steps", ylabel="Magnetization, M",legend=:topright, grid=true, fontsize=40, guidefont=font(14), legendfont=14)
plot!(df_3_1_1_M, df_3_1_1.M_T23, label="T=2.3")
plot!(df_3_1_1_M, df_3_1_1.M_T25, label="T=2.5")
savefig(df_3_1_1_M, "figures/Problem_3_1_1_M.png")

println("Task 1 complete")

df_3_1_2 = CSV.read("DataFrames/H_zero/p_zero/Problem_3_1_2.csv", DataFrame)

df_3_1_2_E = plot(df_3_1_2.T, df_3_1_2.E, xlabel="Temperature, T", ylabel="Energy, E", grid=true, fontsize=40, guidefont=font(14), legend=false)
savefig(df_3_1_2_E, "figures/Problem_3_1_2_E.png")

df_3_1_2_M = scatter(df_3_1_2.T, abs.(df_3_1_2.M), xlabel="Temperature, T", ylabel="Magnetization, M", grid=true, fontsize=40, guidefont=font(14), legend=false)
savefig(df_3_1_2_M, "figures/Problem_3_1_2_M.png")

Problem_3_1_2_grids = heatmap(layout=(1,3),size=(1200, 400))

grid_218 = Matrix(CSV.read("grids/H_zero/p_zero/Problem_3_1_2_T_2.18.txt", DataFrame; delim=' '))
heatmap!(Problem_3_1_2_grids, grid_218, subplot=1, title="T=2.18",color=:bluesreds, grid=false, colorbar=false, titlesize=30, tickfontsize=14)

grid_278 = Matrix(CSV.read("grids/H_zero/p_zero/Problem_3_1_2_T_2.78.txt", DataFrame; delim=' '))
heatmap!(Problem_3_1_2_grids, grid_278, subplot=2, title="T=2.78", color=:bluesreds, grid=false, colorbar=false, titlesize=30, tickfontsize=14)

grid_338 = Matrix(CSV.read("grids/H_zero/p_zero/Problem_3_1_2_T_3.38.txt", DataFrame; delim=' '))   
heatmap!(Problem_3_1_2_grids, grid_338, subplot=3, title="T=3.38", color=:bluesreds, grid=false, colorbar=false, titlesize=30, tickfontsize=14)

savefig(Problem_3_1_2_grids, "figures/Problem_3_1_2_T_grids.png")
println("Task 2 complete")

df_dict = Dict{Int, DataFrame}()
L_arr = collect(10:5:60)
for (i,L) in enumerate(L_arr)
    df_dict[L] = CSV.read("DataFrames/H_zero/p_zero/Problem_3_1_3_L_$L.csv", DataFrame)
end

Problem_3_1_3_C = plot()
for L in L_arr
    plot!(Problem_3_1_3_C, df_dict[L].T, df_dict[L].C, label="L=$L", xlabel="Temperature, T", ylabel="Heat Capacity, C", lw=2, fontsize=40, guidefont=font(14), tickfontsize=14, legendfontsize=14)
end
savefig(Problem_3_1_3_C, "figures/Problem_3_1_3_C.png")

Problem_3_1_3_M = plot()
for L in L_arr
   scatter!(Problem_3_1_3_M, df_dict[L].T, df_dict[L].M, label="L=$L", xlabel="Temperature, T", ylabel="Magnetization, M", lw=2, fontsize=40, guidefont=font(14), tickfontsize=14, legendfontsize=14)
end
savefig(Problem_3_1_3_M, "figures/Problem_3_1_3_M.png")

Problem_3_1_3_E = plot()
for L in L_arr
    plot!(Problem_3_1_3_E, df_dict[L].T, df_dict[L].E, label="L=$L", xlabel="Temperature, T", ylabel="Energy, E", lw=2, fontsize=40, guidefont=font(14), tickfontsize=14, legendfontsize=14)
end
savefig(Problem_3_1_3_E, "figures/Problem_3_1_3_E.png")

p_arr = [0.02, 0.15, 0.3]
df_3_1_4_dict = Dict{Float64, DataFrame}() # Create a dictionary to store the DataFrames
Problem_3_1_4_E = plot()
Problem_3_1_4_M = scatter()
Problem_3_1_4_C = plot()

for p in p_arr
    df_3_1_4_dict[p] = CSV.read("DataFrames/H_zero/p_non_zero/Problem_3_1_4_p_$(p).csv", DataFrame)
    plot!(Problem_3_1_4_E, df_3_1_4_dict[p].T, df_3_1_4_dict[p].E, label="p=$(p)", xlabel="Temperature, T", ylabel="Energy, E", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
    scatter!(Problem_3_1_4_M, df_3_1_4_dict[p].T, df_3_1_4_dict[p].M, label="p=$(p)", xlabel="Temperature, T", ylabel="Magnetization, M", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
    plot!(Problem_3_1_4_C, df_3_1_4_dict[p].T, df_3_1_4_dict[p].C, label="p=$(p)", xlabel="Temperature, T", ylabel="Specific heat, C", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
end
savefig(Problem_3_1_4_E, "figures/Problem_3_1_4_E.png")
savefig(Problem_3_1_4_M, "figures/Problem_3_1_4_M.png")
savefig(Problem_3_1_4_C, "figures/Problem_3_1_4_C.png")

println("Task 4 complete")
"Here GitHub Copilot helped with proper formatting for subplots"
df_3_1_5_dict = Dict{Float64, DataFrame}() # Create a dictionary to store the DataFrames
Problem_3_1_5_E = plot()
Problem_3_1_5_M = scatter()
Problem_3_1_5_C = plot()

optimized_grid_plot = plot(layout=(3,3))
for p in p_arr
    df_3_1_5_dict[p] = CSV.read("DataFrames/H_zero/p_non_zero/Problem_3_1_5_p_$(p).csv", DataFrame)
    
    for i in 1:3
        plot!(Problem_3_1_5_E, df_3_1_5_dict[p].T, df_3_1_5_dict[p][!, Symbol("E_$i")], label="p=$(p), i=$(i)", xlabel="Temperature, T", ylabel="Energy, E", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
        scatter!(Problem_3_1_5_M, df_3_1_5_dict[p].T, df_3_1_5_dict[p][!, Symbol("M_$i")], label="p=$(p), i=$(i)", xlabel="Temperature, T", ylabel="Magnetization, M", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
        
        println("E_p_$(p)_i_$(i): ", df_3_1_5_dict[p][!, Symbol("E_$i")][end])
        println("M_p_$(p)_i_$(i): ", df_3_1_5_dict[p][!, Symbol("M_$i")][end])
        
        impurities = Matrix(CSV.read("grids/H_zero/p_non_zero/impurities_p_$(p)_$(i).txt", DataFrame; delim=' '))
        original_grid = Matrix(CSV.read("grids/H_zero/p_non_zero/original_grid_p_$(p)_$(i).txt", DataFrame; delim=' '))
        optimized_grid = Matrix(CSV.read("grids/H_zero/p_non_zero/optimized_grid$(p)_$(i).txt", DataFrame; delim=' '))

        original_plot = heatmap(impurities .+ original_grid, xlabel="x", ylabel="y", color=:bluesreds, grid=false, xticks=collect(10:5:40), yticks=collect(10:5:40), fontsize=40,  guidefont=font(14), tickfontsize=14)
        savefig(original_plot, "figures/SA/Problem_3_1_5_p_$(p)_i_$(i)_original_grid.png")

        subplot_index = (p == p_arr[1] ? 0 : findfirst(x -> x == p, p_arr) - 1) * 3 + i
        
        heatmap!(optimized_grid_plot, optimized_grid, color=:bluesreds, xlabel="", ylabel="", colorbar=false, title="p=$(p), i=$(i)", subplot=subplot_index)
    end
end

savefig(optimized_grid_plot, "figures/SA/Problem_3_1_5_optimized_grid.png")
savefig(Problem_3_1_5_E, "figures/SA/Problem_3_1_5_E.png")
savefig(Problem_3_1_5_M, "figures/SA/Problem_3_1_5_M.png")

println("Task 5 complete")

df_3_2_1 = CSV.read("DataFrames/H_non_zero/Problem_3_2_1.csv", DataFrame)

Plot_3_2_1_E = plot(df_3_2_1.E_h, xlabel="Steps", ylabel="Energy, E", grid=true, guidefont=font(14),  tickfontsize=14, legend=false)
savefig(Plot_3_2_1_E, "figures/Problem_3_2_1_E.png")

Plot_3_2_1_M = plot(df_3_2_1.M_h, xlabel="Steps", ylabel="Magnetization, M", grid=true, guidefont=font(14),  tickfontsize=14, legend=false)
savefig(Plot_3_2_1_M, "figures/Problem_3_2_1_M.png")

grid_3_2_1_inital = Matrix(CSV.read("grids/H_non_zero/Problem_3_2_1_initial.txt", DataFrame; delim=' '))
grid_3_2_1_final = Matrix(CSV.read("grids/H_non_zero/Problem_3_2_1_final.txt", DataFrame; delim=' '))

Plot_3_2_1_grids = plot(layout=(1,2), size=(800, 400))

heatmap!(Plot_3_2_1_grids, grid_3_2_1_inital, subplot=1, color=:bluesreds, grid=false, title="Initial",  guidefont=font(14), tickfontsize=14, colorbar=false)
heatmap!(Plot_3_2_1_grids, grid_3_2_1_final, subplot=2, color=:bluesreds, grid=false, title="Final", guidefont=font(14), tickfontsize=14, colorbar=false)
savefig(Plot_3_2_1_grids, "figures/Problem_3_2_1_grids.png")

println("Task 2.1 complete")
df_3_2_2_dict = Dict{Float64, DataFrame}() # Create a dictionary to store the DataFrames
H_arr = collect(0.01:0.01:0.05)
Plot_3_2_2_E = plot()
Plot_3_2_2_M = plot()
Plot_3_2_2_C = plot()
Plot_3_2_2_chi = plot()

for H in H_arr
    df_3_2_2_dict[H] = CSV.read("DataFrames/H_non_zero/Problem_3_2_2_H_$(H).csv", DataFrame)
    plot!(Plot_3_2_2_E, df_3_2_2_dict[H].T, df_3_2_2_dict[H].E, label="H=$(H)", xlabel="Temperature, T", ylabel="Energy, E", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
    scatter!(Plot_3_2_2_M, df_3_2_2_dict[H].T, abs.(df_3_2_2_dict[H].M), label="H=$(H)", xlabel="Temperature, T", ylabel="Magnetization, M", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
    plot!(Plot_3_2_2_C, df_3_2_2_dict[H].T, df_3_2_2_dict[H].C, label="H=$(H)", xlabel="Temperature, T", ylabel="Specific heat, C", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
    plot!(Plot_3_2_2_chi, df_3_2_2_dict[H].T, df_3_2_2_dict[H].chi, label="H=$(H)", xlabel="Temperature, T", ylabel="Susceptibility, χ", grid=true, fontsize=40, guidefont=font(14), legendfont=14)
end

savefig(Plot_3_2_2_E, "figures/Problem_3_2_2_E.png")
savefig(Plot_3_2_2_M, "figures/Problem_3_2_2_M.png")
savefig(Plot_3_2_2_C, "figures/Problem_3_2_2_C.png")
savefig(Plot_3_2_2_chi, "figures/Problem_3_2_2_chi.png")
println("Task 2.2 complete")