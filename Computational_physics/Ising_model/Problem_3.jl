using Pkg
# Pkg.add("Plots")
# Pkg.add("ProgressMeter")
# Pkg.add("Statistics")
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("GLM)
# Pkg.add("DelimitedFiles")

using Plots
using ProgressMeter
using Statistics
using DataFrames
using CSV
using GLM
using DelimitedFiles

L = 40    # Lattice size
J = 1.0
mu = 1.0  # Magnetic moment
k_b = 1.0
num_sweeps = 10000  # Number of sweeps
T_arr = collect(1:0.02:3.5)  # Temperature range

function grid_initializer(L)
    # Initialize the grid with random values
    grid = rand([1,-1], (L, L))
    return grid
end

function energy_updater(grid, i, j, L, H)
    """
    Inspired by lecture notebook of MC
    Input: 
        - grid: 2D array of spins
        - i, j: coordinates of the spin to be flipped
    Output:
        - dE: Change of energy from flipping
    """
    new_spin = -grid[i, j]

    left = i > 1 ? grid[i-1, j] : grid[end, j]
    right = i < L ? grid[i+1, j] : grid[1, j]
    up = j > 1 ? grid[i, j-1] : grid[i, end]
    down = j < L ? grid[i, j+1] : grid[i, 1]

    dE_H = 2*mu*H*grid[i,j]
    dE_spin = - 2*J*new_spin*(left + right + up + down)  # Change in energy from flipping
    
    return dE_spin + dE_H
end

function MC_sweep(grid, T, E, M, L, p, impurity_coordinates, H)
    E_arr = [E]  # Store energy value
    M_arr = [M]  # Store magnetization value
    
    for k in 1:L^2
        # Coordinates
        i = rand(1:L)
        j = rand(1:L)
        if p > 0 && (i, j) in impurity_coordinates
            push!(E_arr, E)      # Store energy value
            push!(M_arr, M)      # Store magnetization value
        else
            dE = energy_updater(grid, i, j, L, H)  # Change in energy from flipping
        
            if dE < 0 || rand() < exp(-dE / T)
                grid[i, j] = -grid[i, j]
                E += dE
                M += 2 * grid[i, j]  # Update magnetization
            end
            push!(E_arr, E)      # Store energy value
            push!(M_arr, M)      # Store magnetization value
        end
    end
    return grid, E_arr, M_arr
end

function energy_calc(grid, L, H)
    E_spin = 0.0
    E_h = 0.0
    for i in 1:L
        for j in 1:L
            left = i > 1  ? grid[i-1, j] : grid[L, j]
            right = i < L ? grid[i+1, j] : grid[1, j]
            up = j > 1 ? grid[i, j-1] : grid[i, L]
            down = j < L ? grid[i, j+1] : grid[i, 1]

            E_spin -= 0.5*J*grid[i,j]*(left + right + up + down)
            E_h -= mu*H*grid[i,j]  # Magnetic field energy
        end
    end

    return E_spin + E_h  # Total energy
end
function impurity_production(p, L)
    """ 
    Generates impurities coordinates
    """

    N_impurities = round(Int, p * L^2)  # Number of impurities to be added
    impurity_coordinates = Set((rand(1:L), rand(1:L)) for _ in 1:N_impurities)
       
    return impurity_coordinates
end

function metropolis(num_sweeps, T; L=40, p=0, impurity_coordinates=Set(), grid=nothing, H=0.0)
    """
    Input:
        - grid: 2D array of spins
        - num_sweeps: number of sweeps to perform
        - T: temperature
    Output:
        - grid: updated grid after performing the sweeps
    """
    if isnothing(grid)
        grid = grid_initializer(L)  # Initialize the grid with random values
    end
    
    M_0 = sum(grid)  # Initial magnetization
    
    E_0 = energy_calc(grid, L, H)  # Calculate initial energy
    # println("Initial energy: ", E_0)
    E_arr = [E_0]  # Store initial energy
    M_arr = [M_0]  # Store initial magnetization

    E = E_0  # Initialize energy
    M = M_0  # Initialize magnetization
   
    for _ in 1:num_sweeps
        grid, E, M = MC_sweep(grid, T, E, M, L, p, impurity_coordinates, H)
        append!(E_arr, E)  # Store energy value
        append!(M_arr, M)  # Store magnetization value
        E = E_arr[end]  # Update energy
        M = M_arr[end]  # Update magnetization
    end
    if p > 0
        return grid, E_arr, M_arr/L^2, impurity_coordinates
    else
        return grid, E_arr, M_arr/L^2
    end
end
   
function heat_capacity(E, T)
    """
    Input:
        - E: energy array
        - T: temperature
    Output:
        - C: heat capacity
    """
    
    E_avg = mean(E)         # Average energy
    E_2_avg = mean(E.^2)  # Average of energy squared
    delta_E = E_2_avg - E_avg^2  # Fluctuation in energy
    C = (1/(k_b*T^2)) * delta_E   # Heat capacity
    return C
end

function SA(p::Float64, n_sweeps::Int; T_max::Float64=10.0, T_min::Float64=0.1, cool_rate::Float64=0.96, L::Int=40)
    """
    Inspired by flowchart from Baeldung: https://www.baeldung.com/cs/simulated-annealing
    """
    T = T_max
    x_0 = grid_initializer(L)  # Initialize the grid with random values
    E = energy_calc(x_0, L,0)  # Calculate initial energy
    M = sum(x_0)/L^2  # Initial magnetization
    coord_imp = impurity_production(p, L)  # Generate impurity coordinates
    E_arr = [E]  # Store initial energy
    M_arr = [M]  # Store initial magnetization
    T_arr = [T]  # Store temperature values

    x = copy(x_0)
    
    while T > T_min
        # Perform a sweep
        x_new, E_new, M_new, _ = metropolis(n_sweeps, T, p=p, impurity_coordinates=coord_imp, grid=x)
        
        dE = E_new[end] - E

        if dE < 0 || rand() < exp(-dE / T)
            x = copy(x_new)
            E = E_new[end]
            M = M_new[end]
        end

        push!(E_arr, E)  # Store energy value
        push!(M_arr, M)  # Store magnetization value
        push!(T_arr, T)  # Store temperature value

        # Cool down the temperature
        T *= cool_rate
    end

    return x_0, x, T_arr, E_arr, M_arr, coord_imp
end

function coordinates_to_heatmap(L, coordinates)
    """
    Converts coordinates to a heatmap
    """
    map_for_impurities = zeros(Int, L, L)
    for (x,y) in coordinates
        map_for_impurities[x, y] = 1
    end

    return map_for_impurities
end

function susceptibility(M,T)
    M_avg = mean(M)  # Average magnetization
    M_2_avg = mean(M.^2)  # Average of magnetization squared
    delta_M = M_2_avg - M_avg^2  # Fluctuation in magnetization

    chi = (1/(k_b*T)) * delta_M   # Magnetic susceptibility
    return chi
end

# _, E_T21, M_T21 = metropolis(num_sweeps, 2.1)
# _, E_T23, M_T23 = metropolis(num_sweeps, 2.3)
# _, E_T25, M_T25 = metropolis(num_sweeps, 2.5)

@time begin
    # Task 3.1.1
    println("Task 3.1.1")
    start_time = time()
    df_3_1_1 = DataFrame(
        E_T21=E_T21,
        E_T23=E_T23,
        E_T25=E_T25,
        M_T21=M_T21,
        M_T23=M_T23,
        M_T25=M_T25
    )

    CSV.write("DataFrames/H_zero/p_zero/Problem_3_1_1.csv", df_3_1_1)
    println("Time spent - Task 3.1.1: ", round(time() - start_time, digits=2), " seconds")
        
    # Task 3.1.2
    println("Task 3.1.2")
    start_time = time()

    E_asf_T = zeros(length(T_arr))  # Store energy values
    M_asf_T = zeros(length(T_arr))  # Store magnetization values
    C_asf_T = zeros(length(T_arr))  # Store heat capacity values

    grid_0 = grid_initializer(L)  # Using same initial grid for alle T

    @showprogress  for (i, T) in enumerate(T_arr)
        grid, E_arr, M_arr = metropolis(num_sweeps, T; grid=copy(grid_0))
        cut_off = round(Int,0.15*length(E_arr))
        E_asf_T[i] = mean(E_arr[end-cut_off:end])
        M_asf_T[i] = mean(M_arr[end-cut_off:end])
        C_asf_T[i] = heat_capacity(E_arr[end-cut_off:end], T)
            # Last 15% of the data
        if i % 30 == 0
            writedlm("grids/H_zero/p_zero/Problem_3_1_2_T_$(T).txt", grid, ' ')
        end
    end

    df_3_1_2 = DataFrame(T=T_arr, E=E_asf_T, M=M_asf_T, C=C_asf_T)
    CSV.write("DataFrames/H_zero/p_zero/Problem_3_1_2.csv", df_3_1_2)
    println("Time spent - Task 3.1.2: ", round(time() - start_time, digits=2), " seconds")

    # Task 3.1.3
    println("Task 3.1.3")
    start_time = time()

    L_arr = collect(10:5:60)
    E_asf_L = zeros(length(L_arr), length(T_arr))
    M_asf_L = zeros(length(L_arr), length(T_arr))
    C_asf_L = zeros(length(L_arr), length(T_arr))
    @showprogress for (i, L) in enumerate(L_arr)
        for (j, T) in enumerate(T_arr)
            grid, E_arr, M_arr = metropolis(num_sweeps, T, L=L)

            cut_off = round(Int, 0.15*length(E_arr))    # Last 15% of the data
            E_asf_L[i, j] = mean(E_arr[end-cut_off:end])
            M_asf_L[i, j] = mean(M_arr[end-cut_off:end])

            C_asf_L[i,j] = heat_capacity(E_arr[end-cut_off:end], T)
        end
        df = DataFrame(T=T_arr, E=E_asf_L[i,:], M=M_asf_L[i,:], C=C_asf_L[i,:])
        CSV.write("DataFrames/H_zero/p_zero/Problem_3_1_3_L_$L.csv", df)
    end
    println("Time spent - Task 3.1.3: ", round(time() - start_time, digits=2), " seconds")

    # Task 3.1.4
    println("Task 3.1.4")
    start_time = time()

    p_arr = (0.02, 0.15, 0.30)
    T_arr = (collect(2.0:0.02:3))
    E_asf_p = zeros(length(p_arr), length(T_arr))
    M_asf_p = zeros(length(p_arr), length(T_arr))
    C_asf_p = zeros(length(p_arr), length(T_arr))

    # grid_iter = grid_initializer(L)  # Placeholder to be able to save later
    
    @showprogress for (i, p) in enumerate(p_arr)
        grid_iter = zeros(L,L)  # Empty placeholder
        for (j, T) in enumerate(T_arr)
            impurity_coordinates = impurity_production(p, L)  # 10% impurities
            grid_iter, E_arr, M_arr, _ = metropolis(num_sweeps, T, p=p, impurity_coordinates=impurity_coordinates)
            # println(sum(grid_iter))
            cut_off = round(Int, 0.15*length(E_arr))    # Last 15% of the data
            E_asf_p[i, j] = mean(E_arr[end-cut_off:end])
            M_asf_p[i, j] = mean(M_arr[end-cut_off:end])
            C_asf_p[i,j] = heat_capacity(E_arr[end-cut_off:end], T)
        end

        writedlm("grids/H_zero/p_non_zero/Problem_3_1_4_p_$(p).txt", copy(grid_iter), ' ')
        df = DataFrame(T=T_arr, E=E_asf_p[i,:], M=M_asf_p[i,:], C=C_asf_p[i,:])
        CSV.write("DataFrames/H_zero/p_non_zero/Problem_3_1_4_p_$(p).csv", df)
    end
    println("Time spent - Task 3.1.4: ", round(time() - start_time, digits=2), " seconds")

    # Task 3.1.5
    println("Task 3.1.5")
    start_time = time()
    p_arr = (0.02, 0.15, 0.30)
    @showprogress for p in p_arr
        df = DataFrame()
        for i in 1:3
            # Here ran for 500 sweeps
            original_grid, optimized_grid, T_arr_SA, E_opt_arr, M_opt_arr, impurity_coordinates = SA(p, 500)

            mapping = coordinates_to_heatmap(L, impurity_coordinates)

            writedlm("grids/H_zero/p_non_zero/original_grid_p_$(p)_$(i).txt", original_grid, ' ')
            writedlm("grids/H_zero/p_non_zero/optimized_grid$(p)_$(i).txt", optimized_grid, ' ')
            writedlm("grids/H_zero/p_non_zero/impurities_p_$(p)_$(i).txt", mapping, ' ')
            df[!, "T"] = T_arr_SA
            df[!, "E_$(i)"] = E_opt_arr
            df[!, "M_$(i)"] = M_opt_arr

        end
        
        
        CSV.write("DataFrames/H_zero/p_non_zero/Problem_3_1_5_p_$(p).csv", df)
    end

    println("Time spent - Task 3.1.5: ", round(time() - start_time, digits=2), " seconds")

#     # Task 3.2.1
    println("Task 3.2.1")
    start_time = time()

    grid_0 = grid_initializer(L)  # Initialize the grid with random values
    grid_h, E_h, M_h = metropolis(num_sweeps, 2.3, H=0.01, grid=copy(grid_0))  # Perform the simulation with a small magnetic field

    writedlm("grids/H_non_zero/Problem_3_2_1_Initial.txt", grid_0, ' ')
    writedlm("grids/H_non_zero/Problem_3_2_1_final.txt", grid_h, ' ')

    df_3_2_1 = DataFrame(
        E_h=E_h,
        M_h=M_h
    )

    CSV.write("DataFrames/H_non_zero/Problem_3_2_1.csv", df_3_2_1)
    println("Time spent - Task 3.2.1: ", round(time() - start_time, digits=2), " seconds")

    # Task 3.2.2
    println("Task 3.2.2")
    start_time = time()

#     H_arr = collect(0.01:0.01:0.05)
#     E_h_arr = zeros(length(H_arr), length(T_arr))
#     M_h_arr = zeros(length(H_arr), length(T_arr))
#     C_h_arr = zeros(length(H_arr), length(T_arr))
#     chi_h_arr = zeros(length(H_arr), length(T_arr))
#     @showprogress for (i, H) in enumerate(H_arr)
#         for (j, T) in enumerate(T_arr)
#             grid, E_arr, M_arr = metropolis(num_sweeps, T, H=H)  

#             cut_off = round(Int, 0.15*length(E_arr))    # Last 15% of the data
#             E_h_arr[i,j] = mean(E_arr[end-cut_off:end])
#             M_h_arr[i,j] = mean(M_arr[end-cut_off:end])
#             C_h_arr[i,j] = heat_capacity(E_arr[end-cut_off:end], T)
#             chi_h_arr[i,j] = susceptibility(M_arr[end-cut_off:end], T)
#         end
#         df = DataFrame(T=T_arr, E=E_h_arr[i,:], M=M_h_arr[i,:], C=C_h_arr[i,:], chi=chi_h_arr[i,:])
#         CSV.write("DataFrames/H_non_zero/Problem_3_2_2_H_$(H).csv", df)
#     end
#     println("Time spent - Task 3.2.2: ", round(time() - start_time, digits=2), " seconds")

end