using Plots
using VortexLattice
using DataFrames, CSV
using LaTeXStrings
using Statistics


stations = 3
mathematical_midpoint = 1.8024674446122

"""
    function get_Sref(chord, yle)

Calculates Sref by calculating the areas of the two trapezoids formed by the wing chord geometry. Note that area1 is the area of a triangle when chord[3] is 0, which is the case for 2 of the 3 variable sweeps.
"""
function get_Sref(chord, yle)
    area1 = (chord[3] + chord[2]) * (yle[3] - yle[2])/2
    area2 = (chord[2] + chord[1]) * (yle[2] - yle[1])/2
    return (area1 + area2) * 2
end


"""
    function build_system(...)

Builds a wing system and returns the CD value, which we are aiming to minimize
"""
function build_system(chord;
    xle = zeros(stations), yle = [0, 3.25, 7.5], zle = zeros(stations),
    theta = zeros(stations), phi = zeros(stations), fc = fill((xc) -> 0, stations),
    Sref = get_Sref(chord, yle), ns = 20, nc = 6, mirror = true, spacing_s = Uniform(), spacing_c = Uniform(),
    Vinf = 50, bref = 7.5, symmetric = false, alpha = 5.0 * pi/180, beta = 0.0, omega = [0.0; 0.0; 0.0], fs = Freestream(Vinf, alpha, beta, omega))

    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc; fc=fc, spacing_s = spacing_s, spacing_c = spacing_c, mirror = mirror)
    grids = [grid]
    ratios = [ratio]
    system = System(grids; ratios = ratios)

    cref = mean(chord)
    rref = [0.25 * cref, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rref, Vinf)

    steady_analysis!(system, ref, fs; symmetric = symmetric)

    CF, CM = body_forces(system; frame=Wind())
    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    dCFs, dCMs = stability_derivatives(system)

    println("CL = ", CL, ", CD = ", CDiff)
    return system, CDiff, dCFs, dCMs

end

"""
    function make_paraview(system)

Function writes ParaView Files for visualization.
"""
function make_paraview(system)
    path = joinpath(@__DIR__, "ParaView Files")
    mkpath(path)
    write_vtk(path, system)
end

"""
    function make_empty_results()

Simply generates an empty DataFrame to be used as the collection point for the data generated from our variable sweep.
"""
function make_empty_results()
    results = DataFrame(
        chordValue = Float64[], CDiff = Float64[]
    )
    return results
end


function append_result!(df::DataFrame, param_value, CDiff)
    row = Dict(
    "chordValue" => param_value,
    "CDiff" => CDiff
    )
    push!(df, row)

end


results_c3 = make_empty_results()
for i in range(start = 0.05, stop = 2.0, length = 51)
    system, CDiff, dCFs, dCMs = build_system([2, mathematical_midpoint, i])
    append_result!(results_c3, i, CDiff)
end

results_c2 = make_empty_results()
for i in range(start = 0.05, stop = 2.0, length = 51)
    system, CDiff, dCFs, dCMs = build_system([2, i, 0.05])
    append_result!(results_c2, i, CDiff)
end


results_c1 = make_empty_results()
for i in range(start = 0.05, stop = 2.0, length = 51)
    system, CDiff, dCFs, dCMs = build_system([i, mathematical_midpoint, 0.05])
    append_result!(results_c1, i, CDiff)
end

# print(results_c3)
# print(results_c2)
# print(results_c1)

