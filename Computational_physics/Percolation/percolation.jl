using Pkg
#  Pkg.add("Random")
#  Pkg.add("Plots")
#  Pkg.add("LaTeXStrings")
#  Pkg.add("Statistics")
#  Pkg.add("CSV")
#  Pkg.add("JLD2")
using Random
using Plots
using LaTeXStrings
using Statistics
using CSV
using DataFrames
using JLD2

L = parse(Int, ARGS[1])
println("L:", L)

calc_binom = false

@time begin
    # Parameters
    N = L^2  # Number of nodes
    type = "square" # Type of lattice

    # RNG generator
    seed = 123
    rng = Xoshiro(seed)

    function get_root_node(i::Int, network)
        if network[i] < 0  # Check if i is the root node
            return i
        else
            network[i] = get_root_node(network[i], network)
            return network[i]
        end
    end

    # function activate_bond(k, network, lattice_matrix, cluster_dict, greatest_cluster_ref, avg_s)
    #     i = lattice_matrix[k, 1]  # Node 1
    #     j = lattice_matrix[k, 2]  # Node 2

    #     root_i = get_root_node(i, network)
    #     root_j = get_root_node(j, network)

    #     if root_i != root_j
    #         # Ensure the larger cluster remains the root
    #         if abs(network[root_i]) < abs(network[root_j])
    #             root_i, root_j = root_j, root_i  # Swap to make root_i the larger cluster
    #         end

    #         avg_s = avg_s - network[root_i]^2 - network[root_j]^2 + (network[root_i] + network[root_j])^2  # New size of the merged cluster
    #         # Merge root_j into root_i
    #         network[root_i] -= abs(network[root_j])  # Update size of the merged cluster

    #         # Update all nodes in cluster_j to point to root_i
    #         if haskey(cluster_dict, root_j)
    #             for node in cluster_dict[root_j]
    #                 network[node] = root_i  # Update parent
    #             end
    #             # Merge the node lists
    #             cluster_dict[root_i] = vcat(cluster_dict[root_i], cluster_dict[root_j])
    #             delete!(cluster_dict, root_j)  # Remove the old cluster entry
    #         end

    #         # Update the greatest cluster reference if necessary
    #         if abs(network[root_i]) > abs(network[greatest_cluster_ref[]])
    #             greatest_cluster_ref[] = root_i
    #         end
                    
    #     end
    #     return avg_s
    # end

    # New and faster
    function activate_bond(k, network, lattice_matrix, greatest_cluster_ref, avg_s)
        i = lattice_matrix[k, 1]  # Node 1
        j = lattice_matrix[k, 2]  # Node 2
    
        root_i = get_root_node(i, network)  # Find root of cluster i
        root_j = get_root_node(j, network)  # Find root of cluster j
    
        if root_i != root_j
            # Ensure the larger cluster becomes the new root
            if abs(network[root_i]) < abs(network[root_j])
                root_i, root_j = root_j, root_i  # Swap to ensure root_i is larger
            end
    
            # Update average cluster size correctly (using positive sizes)
            avg_s = avg_s - (abs(network[root_i]))^2 - (abs(network[root_j]))^2 + (abs(network[root_i]) + abs(network[root_j]))^2
    
            # Merge clusters
            network[root_i] -= abs(network[root_j])  # Update size (stored as negative value)
            network[root_j] = root_i  # Point root_j directly to root_i (only root_j is updated)
    
            # Update greatest cluster if needed
            if abs(network[root_i]) > abs(network[greatest_cluster_ref[]])
                greatest_cluster_ref[] = root_i
            end
        end
    
        return avg_s
    end
    
    

    function lattice_loader(file_name::String)
        if !isfile(file_name)
            error("File $file_name does not exist")
        end
        lattice = []
        open(file_name, "r") do file_name
            for line in eachline(file_name)
                push!(lattice, parse(Int, line))
            end
        end
        lattice_matrix = reshape(lattice, 2, :)'
    end

    function shuffle_bonds(lattice_matrix, rng)
        M = size(lattice_matrix, 1)
        for i in 2:M  
            r = rand(rng, i:M) 
            lattice_matrix[[i-1, r], :] = lattice_matrix[[r, i-1], :]
        end
        return lattice_matrix
    end

    # const log_fact_dict = Dict{Int, Float64}()  # Store log factorials for memoization
    # const binom_coeff = Dict{Int, Float64}()  # Store binomial coefficients for memoization
    # const log_fact_lock = ReentrantLock()
    # const binom_lock = ReentrantLock()

    # function log_factorial(n::Int)
    #     if n < 0
    #         throw(ArgumentError("n must be non-negative"))
    #     elseif haskey(log_fact_dict, n)
    #         return log_fact_dict[n]
    #     else
    #         result = sum(log.(1:n))
            
    #         lock(log_fact_lock) do
    #             log_fact_dict[n] = result
    #         end 
    #         return result
    #     end
    # end

    # function logbinom(n::Int, k::Int)
    #     if k < 0 || k > n
    #         return -Inf   # Invalid coeff
    #     end
    #     return log_factorial(n) - log_factorial(k) - log_factorial(n-k)
    # end

    # function precompute_log_factorials(max_n::Int)
    #     # log_values = accumulate(+, log.(1:max_n))
    #     # for i in 1:max_n
    #     #     log_factorial_dict[i] = log_values[i]
    #     # end
    #     # log_factorial_dict[0] = 0.0
    #     @threads for i in 0:max_n
    #         log_factorial(i)
    #     end
    # end

    # function compute_binom_coeffs(M::Int)
    #     @threads for n in 1:M
    #         result = logbinom(M, n)
            
    #         # Thread-safe dictionary update
    #         lock(binom_lock) do
    #             binom_coeff[n] = result
    #         end
    #     end
    # end

    function binom_coeff_row(m::Int)
        # row = Vector{BigInt}(undef, m+1)  # For non-log
        # row[1] = BigInt(1)  # C(0, 0) = 1
        # for n in 1:m
        #     row[n+1] = row[n] * (m - n + 1) ÷ n
        # end
        row = Vector{Float64}(undef, m+1)  # For log
        row[1] = 0.0  # log(1) = 0
        for n in 1:m
            row[n+1] = row[n] + log(m - n + 1) - log(n)  # For log
        end

        return row
    end
    
    function save_binomial_jld2(filename::String, ms::Vector{Int})
        JLD2.jldopen(filename, "w") do file  # Open the file in write mode
            file["ms"] = ms  # Save the list of M values
            for m in ms
                row = binom_coeff_row(m)
                file["row_$m"] = row  # Save each row under a unique key
            end
        end
    end
    
    function load_binomial_row(filename::String, m::Int)
        return JLD2.load(filename, "row_$m")
    end

    function key_exists_in_jld2(filename::String, key::String)
        JLD2.jldopen(filename, "r") do file  # Open the file in read mode
            return key in keys(file)  # Check if the key exists
        end
    end
    

    # Main function

    function percolater(lattice_matrix, N, M)
        avg_s_val = N                    # Initial value for avgerage_s (sum of the size squared for all clusters)
        network = fill(Int(-1), N)       # Initialize network with all nodes as root nodes
        shuffle_bonds(lattice_matrix, rng)    # Shuffle bonds
        #print(lattice_matrix)
        # prepare cluster dictionary (dict with root node indeces as keys and list of nodes in cluster as values)
        # cluster_dict = Dict{Int, Vector{Int}}()
        # for i in 1:N
        #     cluster_dict[i] = [i]
        # end

        # Pointer for saving values from loop of activate_bonds
        greatest_cluster_ref = Ref(1)
        k_max = Ref(1)

        # Arrays to save quantities of interest
        P_inf = zeros(M)
        avg_s = ones(M)*N

        # Iterating over all bonds
        for k in 1:M 
            
            # Activating bonds
            # avg_s_val = activate_bond(k, network, lattice_matrix, cluster_dict, greatest_cluster_ref, avg_s_val)
            avg_s_val = activate_bond(k, network, lattice_matrix, greatest_cluster_ref, avg_s_val)
            P_inf[k] = abs(network[greatest_cluster_ref[]]) / N
            avg_s[k] = avg_s_val
            
            # # Save the biggest cluster evently ten times evenly to plot
            # T = M/10
            # if k%(T) == 0
            #     greatest_cluster = greatest_cluster_ref[]
                # cluster_nodes = Int[]
                # for n in 1:N
                #     if get_root_node(n, network) == greatest_cluster
                #         push!(cluster_nodes, n)
                #     end
                # end
                # push!(pictures, cluster_nodes)
            # end

            # Cut off if all nodes are in the same cluster
            if abs(network[greatest_cluster_ref[]]) == N
                k_max[] = k    # Cut off, when all sites are in the same cluster
                # println(k_max[])
                break
            end
        end 

        if k_max[] < M
            P_inf[k_max[]+1:end] .= P_inf[k_max[]]    # Cut off rest of the array
            avg_s[k_max[]+1:end] .= avg_s[k_max[]]    # Cut off rest of the array
        end
        # Fill elements after the cut off with the last value, since no change will occur
        
        average_cluster_size = (avg_s .- (N*P_inf).^2) ./ (N*(1 .- P_inf))    # <s>

        indicies_to_change = findall(P_inf .> 0.9999)      # Avoid <s>->\infty when P_inf = 1
        average_cluster_size[indicies_to_change] .= 0

        return P_inf, average_cluster_size
    end

    # Generating lattice


    # Function for dynamical file name registration
    function format_with_underscores(n::Int)
        s = string(n)
        reversed_chunks = reverse.(Iterators.partition(reverse(s), 3))
        return join(reverse(reversed_chunks), '_')
    end

    N_str = format_with_underscores(N) * "_$type"
    file_name = "Lattices/lattice_$N_str.txt"
    println(file_name)

    # Load lattice and reshape into matrix
    lattice_matrix = lattice_loader(file_name)

    M = size(lattice_matrix)[1]  # Number of bonds


    # Calculate the binomial coefficients
    println("computing binomials:")
    start_time = time()
    if calc_binom  # Only has to be computed once
        # Precompute binomial coefficients for all M values
        N_arr = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000].^2
        M_square = 2*N_arr
        M_triangle = 3*N_arr
        M_honeycomb = round.(Int,(3/2)*N_arr)
        M_arr = vcat(M_square, M_triangle, M_honeycomb)
        M_arr = sort(M_arr, rev=false)  # Sort in ascending order

        save_binomial_jld2("binomials.jld2", M_arr)  # Save the binomial coefficients to a file
    end

    # Check that the correct binomial coefficients are calculated
    if !key_exists_in_jld2("binomials.jld2", "row_$M")
        save_binomial_jld2("binomials.jld2", [M])
    end
    println("Time spent - Computing binomials: ", round(time() - start_time, digits=2), " seconds")
    
    # Main

    num_it = 100   # Number of iterations of the percolater for averaging

    pictures = []  # Array for saving plots of the biggest cluster

    P_arr = zeros(num_it, M)                # Array for saving P_inf    
    average_cluster_size = zeros(num_it, M) # Array for saving <s>

    println("percolating:")
    start_time = time()
    # Multiple iterations of the percolater
    for l in 1:num_it
        
        P_inf, avg_s = percolater(lattice_matrix, N, M)
        P_arr[l, :] = P_inf
        average_cluster_size[l, :] = avg_s   

    end


    P_avg = vec(mean(P_arr, dims=1))
    P_2_avg = vec(mean(P_arr.^2, dims=1))

    # Add values before activatign bonds
    pushfirst!(P_avg, 0.0)  # Add intial value to the start of the array
    pushfirst!(P_2_avg, 0.0)  # Add intial value to the start of the array

    chi = N*sqrt.(abs.(P_2_avg .- P_avg.^2))
    # p = (1:M) / M   # probability of activating a node, with cut-off

    average_cluster_size_avg = vec(mean(average_cluster_size, dims=1))
    pushfirst!(average_cluster_size_avg, 1.0)  # Add intial value to the start of the array
    
    println("Time spent - Percolation: ", round(time() - start_time, digits=2), " seconds")
    
    println("Convolving:")
    start_time = time()
    # Convolution
    num_it_q = 10_000
    q = range(0.0, 1.0, length=num_it_q)

    P_conv = zeros(length(q))
    P_2_conv = zeros(length(q))
    average_cluster_size_conv = zeros(length(q))

    log_q = log.(q.+1e-12)  # Avoid log(0) by adding a small value
    log1_q = log1p.(-q.+1e-12) # Avoid log1p(-1) by adding a small value

    binom_coeff_array = load_binomial_row("binomials.jld2", M)  # Load binom_coefficients from file
    for i in 1:length(q)
        # P_conv[i] = 0.0
        # average_cluster_size_conv[i] = N
        # for n in 1:M
        #     log_P = binom_coeff[n] + n*log(q[i]) + (M-n)*log1p(- q[i]) + log(P_avg[1, n]) #log1p for better stability when 1-q[i]->1
        #     log_P_2 = binom_coeff[n] + n*log(q[i]) + (M-n)*log1p(- q[i]) + log(P_2_avg[1, n])
        #     log_s = binom_coeff[n] + n*log(q[i]) + (M-n)*log1p(- q[i]) + log(average_cluster_size_avg[1, n])

        #     P_conv[i] += exp(log_P)
        #     P_2_conv[i] += exp(log_P_2)
        #     average_cluster_size_conv[i] += exp(log_s)
        # end
        P_conv[i] = sum(exp.(binom_coeff_array[1:end] .+ (0:M) .* log_q[i] .+ (M .- (0:M)) .* log1_q[i] .+ log.(P_avg[1:end])))
        P_2_conv[i] = sum(exp.(binom_coeff_array[1:end] .+ (0:M) .* log_q[i] .+ (M .- (0:M)) .* log1_q[i] .+ log.(P_2_avg[1:end])))
        average_cluster_size_conv[i] = sum(exp.(binom_coeff_array[1:end] .+ (0:M) .* log_q[i] .+ (M .- (0:M)) .* log1_q[i] .+ log.(average_cluster_size_avg[1:end])))

    end

    argument = P_2_conv .- P_conv.^2
    argument[findall(argument .< 0)] .= 0      # Avoid all zero values

    chi_conv = N*sqrt.(argument)

    df = DataFrame(q=q, P=P_conv, chi=chi_conv, average_cluster_size=average_cluster_size_conv)
    N_str = format_with_underscores(N)
    q_str = format_with_underscores(num_it_q)
    CSV.write("DataFrames/N$(N_str)_q$(q_str)_$(type).csv", df) 
    println("Time spent - convolution: ", round(time() - start_time, digits=2), " seconds")
end
println("done")