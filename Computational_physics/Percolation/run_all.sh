#!/bin/bash

# How to run:
# 1) Want to calculate all binomial coefficients again?
#      - Set calc_binom = true in percolation.jl
# 2) Change the values of args to match the L you want
# 3) Change "type" in percolation.jl to: "square",  "triangular", "honeycom"
# 4) Run this script with the command: bash run_all.sh
# Array of arguments
args=(100 200)

# Loop through each argument and run percolation.jl
for arg in "${args[@]}"; do
    julia percolation.jl "$arg"
done