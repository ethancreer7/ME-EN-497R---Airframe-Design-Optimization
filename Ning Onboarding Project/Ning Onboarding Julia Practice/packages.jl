using CSV
using DataFrames
using LinearAlgebra
using Interpolations

# data = CSV.read("/Users/ethancreer/Library/CloudStorage/OneDrive-BrighamYoungUniversity/Junior Year/Research/Ning Onboarding Julia Practice/E203.csv", DataFrame)
# println(data)

xs = 1:0.2:5
A = log.(xs)


itp = interpolate((xs,), A, Gridded(Linear()))
println(itp(3.5))