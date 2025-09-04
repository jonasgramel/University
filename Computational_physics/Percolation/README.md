# Project: Percolation

Percolation theory is an important method in many different areas of science. In this paper I have simulated clusters by connecting lattice sites through bonds for studying bond percolation. I found that the number of nearest neighbors in a lattice affects the critical probabilities, although the critical behavior itself is similar for two-dimensional lattices. The results indicates also that the triangular lattice is the dual lattice of the honeycomb lattice. 

### The files are organized as follows:
- Lattices: Lists of the two lattice sites in each bond. Every two subsequent lines are connected. The number is the number of bonds.
- DataFrames: CSV files with data for the giant component, average cluster size and the susceptibility. N is the number of bonds, q is the number of probability values
- figures: Plots
- Lattices.ipynb: Code for generating lattices
- Percolation.ipynb/percolation.jl: Code for analyzing percolation
- run_all.sh: script to automatic run percolation for multiple lattices
