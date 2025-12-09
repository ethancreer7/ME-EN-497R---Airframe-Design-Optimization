using Plots
using CSV, DataFrames, Statistics

path = joinpath(@__DIR__, "wakelength", "wakelength3.111", "wakelength3.111_convergence.csv")
df = CSV.read(path, DataFrame)
x_axis = 
CL_vals = df[:, :CL]
CD_vals = df[:, :CD]

wakelength_plot_CL = plot()