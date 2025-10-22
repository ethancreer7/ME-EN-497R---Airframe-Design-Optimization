using VortexLattice
using DataFrames, CSV
using Plots
using LaTeXStrings

"""
    volume_Ratio(stabilizer_S, x, wing_s, length_val)

Returns volume ratio (either horizontal or vertical depending on passed in arguments)
"""
function volume_Ratio(stabilizer_S, x, wing_s, length_val)
    return (stabilizer_S * x) / (wing_s * length_val)
end

"""
    function find_x(cg, mac, xle, x_shift=4.0)

Returns the value of x_h or x_v to be used in the volume_Ratio function
"""
function find_x(cg, mac, xle, x_shift=4.0)
    return xle[1] + x_shift + 0.25 * mac - cg
end

"""
    function find_mac(chord)

Returns the mean aerodynamic chord length of the given wing parameters
"""
function find_mac(chord)
    return 2 / 3 * (chord[1] + chord[2] - (chord[1] * chord[2])/(chord[1] + chord[2]))
end

"""
    function get_Stabilizer_Area(root_chord, tip_chord, half_span)

Approximates S_ref for vertical and horizontal stabilizer as a trapezoid. Note that this computes only half of the area for the horizontal stabilizer, and the whole area for the vertical stabilizer
"""
function get_Stabilizer_Area(root_chord, tip_chord, half_span)
    return (root_chord + tip_chord) * half_span / 2
end


"""
    function get_V_h(chord_h, chord, cg, xle_h, yle_h, Sref)

Combines functions to get the horizontal tail volume ratio
"""
function get_V_h(chord_h, chord, cg, xle_h, yle_h, Sref)
    mac_h = find_mac(chord_h)
    mac_w = find_mac(chord)
    x = find_x(cg, mac_h, xle_h)
    S_h = 2 * get_Stabilizer_Area(chord_h[1], chord_h[2], yle_h[2])
    V_h = volume_Ratio(S_h, x, Sref, mac_w)    
    return V_h
end

"""
    function get_V_v(chord_v, bref, cg, xle_v, yle_v, Sref)

Combines functions to get the Vertical Tail Volume Ratio
"""
function get_V_v(chord_v, bref, cg, xle_v, zle_v, Sref)
    mac_v = find_mac(chord_v)
    b_w = 2 * bref # The bref passed in is just half the wing span.
    x = find_x(cg, mac_v, xle_v)
    S_v = get_Stabilizer_Area(chord_v[1], chord_v[2], zle_v[2])
    V_v = volume_Ratio(S_v, x, Sref, b_w)
    return V_v
end

"""
    function get_sweep(xle, yle)
    
Calculates the sweep angle that can be recorded in the csv file. Rounds to 3 decimal points.
"""
function get_sweep(xle, yle)
    x = yle[2] - yle[1]
    y = xle[2] - xle[1]
    sweep_angle = atand(y,x)
    return round(sweep_angle; digits=3)
end

"""
    function get_dihedral()
        
Returns the dihedral angle to be recorded in the csv file. Rounds to 3 decimal Points.
"""
function get_dihedral(phi)
    dihedral_angle = rad2deg(phi[2])
    return round(dihedral_angle; digits=3)
end


"""
    save_results(df, filename)

Saves the results dataframe to a csv file with name filename
"""
function save_results(df, filename)

    save_path = joinpath(@__DIR__, "CSV Files", filename)
    CSV.write(save_path, df)
end

"""
    function append_result!(df::DataFrame, name, coeffs, derivs; param_name = nothing, param_value = nothing)

Appends one row of aerodynamic coefficients and stability derivatives to a specified dataframe.
If supplied, the parameter being varied (and its value) is pushed to the dataframe as well (e.g. "VolumeRatio", 0.45).
"""
function append_result!(df::DataFrame, name, coeffs, derivs; param_name = nothing, param_value = nothing)
        row = Dict(
        "Name" => name,
        "Parameter" => param_name,
        "Value" => param_value,

        "CD" => coeffs[1], "CY" => coeffs[2], "CL" => coeffs[3],
        "Cl" => coeffs[4], "Cm" => coeffs[5], "Cn" => coeffs[6],

        "CDa" => derivs[1], "CYa" => derivs[2], "CLa" => derivs[3],
        "Cla" => derivs[4], "Cma" => derivs[5], "Cna" => derivs[6],

        "CDb" => derivs[7], "CYb" => derivs[8], "CLb" => derivs[9],
        "Clb" => derivs[10], "Cmb" => derivs[11], "Cnb" => derivs[12],

        "CDp" => derivs[13], "CYp" => derivs[14], "CLp" => derivs[15],
        "Clp" => derivs[16], "Cmp" => derivs[17], "Cnp" => derivs[18],

        "CDq" => derivs[19], "CYq" => derivs[20], "CLq" => derivs[21],
        "Clq" => derivs[22], "Cmq" => derivs[23], "Cnq" => derivs[24],

        "CDr" => derivs[25], "CYr" => derivs[26], "CLr" => derivs[27],
        "Clr" => derivs[28], "Cmr" => derivs[29], "Cnr" => derivs[30]
    )
    push!(df, row)
end

"""
    function make_empty_results()

Simply generates an empty DataFrame to be used as to carry the results before writing the csv file.
"""
function make_empty_results()
    results = DataFrame(
    Name = String[], Parameter = String[], Value = Float64[],
    CD = Float64[], CY = Float64[], CL = Float64[],
    Cl = Float64[], Cm = Float64[], Cn = Float64[],
    CDa = Float64[], CYa = Float64[], CLa = Float64[],
    Cla = Float64[], Cma = Float64[], Cna = Float64[],
    CDb = Float64[], CYb = Float64[], CLb = Float64[],
    Clb = Float64[], Cmb = Float64[], Cnb = Float64[],
    CDp = Float64[], CYp = Float64[], CLp = Float64[],
    Clp = Float64[], Cmp = Float64[], Cnp = Float64[],
    CDq = Float64[], CYq = Float64[], CLq = Float64[],
    Clq = Float64[], Cmq = Float64[], Cnq = Float64[],
    CDr = Float64[], CYr = Float64[], CLr = Float64[],
    Clr = Float64[], Cmr = Float64[], Cnr = Float64[]
)
    return results
end

"""
    function perform_analysis!(system, ref, fs; surface_id)

Performs the fundamental analysis on stability derivatives and coefficients by calling VortexLattice functions
"""
function perform_analysis!(system, ref, fs; symmetric, surface_id)
    steady_analysis!(system, ref, fs; symmetric=symmetric, surface_id=surface_id)

    CF, CM = body_forces(system; frame=Wind())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    dCFs, dCMs = stability_derivatives(system)

    # traditional names for each stability derivative
    CDa, CYa, CLa = dCFs.alpha
    Cla, Cma, Cna = dCMs.alpha
    CDb, CYb, CLb = dCFs.beta
    Clb, Cmb, Cnb = dCMs.beta
    CDp, CYp, CLp = dCFs.p
    Clp, Cmp, Cnp = dCMs.p
    CDq, CYq, CLq = dCFs.q
    Clq, Cmq, Cnq = dCMs.q
    CDr, CYr, CLr = dCFs.r
    Clr, Cmr, Cnr = dCMs.r

    coefficients = [CD, CY, CL, Cl, Cm, Cn]
    stability_derivs = [CDa, CYa, CLa, Cla, Cma, Cna, CDb, CYb, CLb, Clb, Cmb, Cnb, CDp, CYp, CLp, Clp, Cmp, Cnp, CDq, CYq, CLq, Clq, Cmq, Cnq, CDr, CYr, CLr, Clr, Cmr, Cnr]

    return coefficients, stability_derivs

end

"""
    function build_system()

Returns a fully assembled wing system so that changed parameters can be passed in without needing to write so many lines of code for each system

"""
function build_system(;
    xle = [0.0, 0.2], yle = [0.0, 5.0], zle = [0.0, 1.0],
    chord = [1.0, 0.6], theta = [2.0 * pi/180, 2.0*pi/180], phi = [0.0, 0.0],
    fc = fill((xc) -> 0, 2), ns = 8, nc = 4, spacing_s = Uniform(), spacing_c = Uniform(), mirror = true,
    xle_h = [0.0, 0.14], yle_h = [0.0, 1.25], zle_h = [0.0, 0.0],
    chord_h = [0.7, 0.42], theta_h = [0.0, 0.0], phi_h = [0.0, 0.0],
    fc_h = fill((xc) -> 0, 2), ns_h = 5, nc_h = 3, spacing_s_h = Uniform(), spacing_c_h = Uniform(), mirror_h = true,
    xle_v = [0.0, 0.14], yle_v = [0.0, 0.0], zle_v = [0.0, 1.0], chord_v = [0.7, 0.42], theta_v = [0.0, 0.0], phi_v = [0.0, 0.0],
    fc_v = fill((xc) -> 0, 2), ns_v = 4, nc_v = 3, spacing_s_v = Uniform(), spacing_c_v = Uniform(), mirror_v = false,
    Sref = 9.0, cref = 0.9, bref = 10.0, rref = [0.5, 0.0, 0.0], Vinf = 1.0, ref = Reference(Sref, cref, bref, rref, Vinf),
    alpha = 5.0*pi/180, beta = 0.0, Omega = [0.0; 0.0; 0.0], fs = Freestream(Vinf, alpha, beta, Omega),
    symmetric = [false, false, false]
    ) # All of the above default parameters are set for the "base wing system" --- they're pulled directly fromt he example on the VortexLattice.jl docs. Our iterations will call this function repeatedly to vary the system for plotting.

    wgrid, wratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc, mirror=mirror, fc = fc, spacing_s = spacing_s, spacing_c = spacing_c)
    hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, fc = fc_h, spacing_s = spacing_s_h, spacing_c = spacing_c_h)
    VortexLattice.translate!(hgrid, [4.0, 0.0, 0.0])
    vgrid, vratio = wing_to_grid(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, fc = fc_v, spacing_s = spacing_s_v, spacing_c = spacing_c_v)
    VortexLattice.translate!(vgrid, [4.0, 0.0, 0.0])

    grids = [wgrid, hgrid, vgrid]
    ratios = [wratio, hratio, vratio]
    surface_id = [1, 2, 3]

    system = System(grids; ratios)
    return system, ref, fs, surface_id, chord, chord_h, chord_v, xle_h, xle_v, Sref, rref[1], yle_h, yle_v, bref, zle_v, phi, xle, yle, symmetric
end

"""Vary Horizontal Tail Volume Ratios - Here we vary the S_h parameter by varying the yle_h geometry of the horizontal stabilizer -- essentially making the stabilizer longer in the y-direction"""
results_Vh = make_empty_results()
for (count, i) in enumerate(range(start = 1.0, stop = 1.50, length = 51))
    system, ref, fs, surface_id, chord, chord_h, chord_v, xle_h, xle_v, Sref, cg, yle_h, yle_v, bref, zle_v, phi, xle, yle, symmetric = build_system(yle_h = [0.0, i])
    V_h = get_V_h(chord_h, chord, cg, xle_h, yle_h, Sref)
    coeffs, derivs = perform_analysis!(system, ref, fs; symmetric = symmetric, surface_id = surface_id)
    append_result!(results_Vh, "Iteration $count", coeffs, derivs; param_name="Horizontal Tail Volume Ratio", param_value=V_h)
end
save_results(results_Vh, "HorizontalTailRatioVary.csv")


"""Vary Vertical Tail Volume Ratios"""
results_Vv = make_empty_results()
for (count, i) in enumerate(range(start = 1.75, stop = 2.25, length = 51))
    system, ref, fs, surface_id, chord, chord_h, chord_v, xle_h, xle_v, Sref, cg, yle_h, yle_v, bref, zle_v, phi, xle, yle, symmetric = build_system(zle_v = [0.0, i])
    V_v = get_V_v(chord_v, bref, cg, xle_v, zle_v, Sref)
    coeffs, derivs = perform_analysis!(system, ref, fs; symmetric = symmetric, surface_id = surface_id)
    append_result!(results_Vv, "Iteration $count", coeffs, derivs; param_name="Vertical Tail Volume Ratio", param_value=V_v)
end
save_results(results_Vv, "VerticalTailRatioVary.csv")



"""Vary Sweep"""
results_sweep = make_empty_results()
for (count, i) in enumerate(range(start = -0.1, stop = 3, length = 51)) # A xle range iteration that models a wing from forward sweep to traditional backwards sweep
    system, ref, fs, surface_id, chord, chord_h, chord_v, xle_h, xle_v, Sref, cg, yle_h, yle_v, bref, zle_v, phi, xle, yle, symmetric = build_system(xle = [0.0, i])
    sweep_angle = get_sweep(xle, yle)
    coeffs, derivs = perform_analysis!(system, ref, fs; symmetric = symmetric, surface_id = surface_id)
    append_result!(results_sweep, "Iteration $count", coeffs, derivs; param_name = "Sweep Angle", param_value=sweep_angle)
end
save_results(results_sweep, "SweepVary.csv")

"""Vary Dihedral"""
results_dihedral = make_empty_results()
for (count, i) in enumerate(range(start = -0.08, stop=0.18, length = 51)) # A rough range of -5 degrees to 10 degrees for the dihedral/anhedral angle
    system, ref, fs, surface_id, chord, chord_h, chord_v, xle_h, xle_v, Sref, cg, yle_h, yle_v, bref, zle_v, phi, xle, yle, symmetric = build_system(phi = [0.0, i])    
    dihedral_angle = get_dihedral(phi)
    coeffs, derivs = perform_analysis!(system, ref, fs; symmetric = symmetric, surface_id = surface_id)
    append_result!(results_dihedral, "Iteration $count", coeffs, derivs; param_name = "Dihedral Angle", param_value = dihedral_angle)
end
save_results(results_dihedral, "DihedralVary.csv")


"""
    function make_df(filename)

Takes the provided csv file and converts it into a plotable dataframe 
"""
function make_df(filename)
    path = joinpath(@__DIR__, "CSV Files")
    mkpath(path)
    file_path = joinpath(@__DIR__, "CSV Files", filename)
    return CSV.read(file_path, DataFrame)
end

latex_map = Dict( # This dictionary lets us translate between notation used in the csv file to latex math commands in our plots
    # Coefficients
    "CD" => L"C_D",
    "CY" => L"C_Y",   
    "CL" => L"C_L",
    "Cl" => L"C_l",
    "Cm" => L"C_m",
    "Cn" => L"C_n",
    
    # Angle of Attack derivatives (alpha)
    "CDa" => L"C_{D,\alpha}",
    "CYa" => L"C_{Y,\alpha}",
    "CLa" => L"C_{L,\alpha}",
    "Cla" => L"C_{l,\alpha}",
    "Cma" => L"C_{m,\alpha}",
    "Cna" => L"C_{n,\alpha}",

    # Sideslip derivatives (beta)
    "CDb" => L"C_{D,\beta}",
    "CYb" => L"C_{Y,\beta}",
    "CLb" => L"C_{L,\beta}",
    "Clb" => L"C_{l,\beta}",
    "Cmb" => L"C_{m,\beta}",
    "Cnb" => L"C_{n,\beta}",

    # Roll rate derivatives (p)
    "CDp" => L"C_{D,p}",
    "CYp" => L"C_{Y,p}",
    "CLp" => L"C_{L,p}",
    "Clp" => L"C_{l,p}",
    "Cmp" => L"C_{m,p}",
    "Cnp" => L"C_{n,p}",

    # Pitch rate derivatives (q)
    "CDq" => L"C_{D,q}",
    "CYq" => L"C_{Y,q}",
    "CLq" => L"C_{L,q}",
    "Clq" => L"C_{l,q}",
    "Cmq" => L"C_{m,q}",
    "Cnq" => L"C_{n,q}",

    # Yaw rate derivatives (r)
    "CDr" => L"C_{D,r}",
    "CYr" => L"C_{Y,r}",
    "CLr" => L"C_{L,r}",
    "Clr" => L"C_{l,r}",
    "Cmr" => L"C_{m,r}",
    "Cnr" => L"C_{n,r}"
)


df_Vh = make_df("HorizontalTailRatioVary.csv")
rename!(df_Vh, Dict(col => get(latex_map, col, col) for col in names(df_Vh)))
df_Vv = make_df("VerticalTailRatioVary.csv")
rename!(df_Vv, Dict(col => get(latex_map, col, col) for col in names(df_Vv)))
df_sweep = make_df("SweepVary.csv")
rename!(df_sweep, Dict(col => get(latex_map, col, col) for col in names(df_sweep)))
df_dihedral = make_df("DihedralVary.csv")
rename!(df_dihedral, Dict(col => get(latex_map, col, col) for col in names(df_dihedral)))

"""
-------------------------------------
BELOW BEGINS THE PLOTTING OF THE DATA
-------------------------------------
"""



"""
    function plot_Vh(df, var, color, label, ylabel)

Plots a stability derivative versus the given dataframes varied parameter
"""
function plot_df(df, var, color, label, ylabel, xlabel)
    base_plot = plot(df[!, "Value"], df[!, var], linewidth = 1.5, linestyle=:solid, color=color, label = label, xlabel = xlabel, ylabel = ylabel, ylabelfontsize = 12, grid=true, dpi = 600, xticks =:auto, yticks=:auto)
    # display(base_plot)
    return base_plot
end

"""
    function save_plot(filename, plot)

Saves the provided plot into the current directory's "Figures" subfolder with the given filename
"""
function save_plot(filename, plot, subfolder)
    display(plot)
    direct_path = joinpath(@__DIR__, "Figures", subfolder)
    mkpath(direct_path)
    save_path = joinpath(direct_path, filename)
    savefig(plot, save_path)
end

dfs_data_to_plot = [
    (df_Vh, "Horizontal Tail Plots", "Horizontal Tail Volume Ratio", :orange),
    (df_Vv, "Vertical Tail Plots", "Vertical Tail Volume Ratio", :purple),
    (df_sweep, "Sweep Plots", "Sweep Angle (Degrees)", :green),
    (df_dihedral, "Dihedral Plots", "Dihedral Angle (Degrees)", :brown)
]

for (df, subfolder, xlabel, color) in dfs_data_to_plot
    y_cols = names(df)[4:end]
    for var in y_cols
        plot_title = "$(var) vs. $(xlabel)"
        to_plot = plot_df(df, var, color, plot_title, var, xlabel)
        matches = [k for (k, v) in latex_map if v == var] # Reverse engineer the mapping to get the filename
        varname = length(matches) > 0 ? matches[1] : string(var)
        save_plot("$(varname) vs $(xlabel).png", to_plot, subfolder)
    end
end

# Below creates a simple stable airframe using the specified parameters and writes the derivatives to a csv file for analysis. It also writes the system to a ParaView file for visual analysis.
optimized = make_empty_results()
system, ref, fs, surface_id, chord, chord_h, chord_v, xle_h, xle_v, Sref, cg, yle_h, yle_v, bref, zle_v, phi, xle, yle, symmetric = build_system(xle = [0.0, 1.5], phi = [0.0, 0.1], yle_h = [0.0, 1.4], zle_v = [0.0, 1.5])
coeffs, derivs = perform_analysis!(system, ref, fs; symmetric = symmetric, surface_id = surface_id)
append_result!(optimized, "Optimized", coeffs, derivs; param_name = "optimized", param_value = 0)
save_results(optimized, "Optimized.csv")
save_path = joinpath(@__DIR__, "ParaView Files", "OptimizedWing")
write_vtk(save_path, system)
