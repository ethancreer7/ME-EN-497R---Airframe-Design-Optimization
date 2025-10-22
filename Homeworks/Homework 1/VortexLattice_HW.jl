using VortexLattice
using Plots
using LaTeXStrings

# Below is a list of initialized wing design parameters. It includes data for the leading edge, Chord Length, Twist Angle, and a function for the camberline
xle = [0.0, 0.4] # leading edge x-position --- x direction is in the direction of the freestream vector. 0.0m at the root (origin), and 0.4 at the tip (provides a sweept back wing)
yle = [0.0, 7.5] # leading edge y-position --- y direction is in the spanwise direction of the aircraft. 0.0m at the root (origin), and 7.5m at the tip
zle = [0.0, 0.0] # leading edge z-position --- z direction is the vertical direction. 0.0m at the root (origin), and 0.0m at the tip. This indicates the wing has no designed dihedral
chord = [2.2, 1.8] # chord length goes from 2.2m to 1.8m over the span of the wing (again indicating a swept wing design)
theta = [2.0*pi/180, 2.0*pi/180] # twist (in Radians) - constant over the span of the wing
phi = [0.0, 0.0] # section rotation about the x-axis (dihedral is zero across the span of the wing)
fc = fill((xc) -> 0, 2) # camberline function for each section (y/c = f(x/c)) - Anonymous function dictates camberline is always flat (0). Assigns this as a vector for the 2 points specified above.


# The below section begins the discretization process across the wing. In this example, we divide the wing into 72 panels and develop
ns = 12 # number of spanwise panels
nc = 6  # number of chordwise panels
spacing_s = Sine()
spacing_c = Uniform() 


grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true) #This function puts everything together from the previous two sections and sets up the panel grid for mathematical analysis in the subsequent steps. It outputs a grid of panels, and ratios for where the control points will be placed along the chord of each panel


grids = [grid] # This line vectorizes the grid created in the above function
ratios = [ratio] # This line vectorizes the ratios created in the above function
system = System(grids; ratios) #Creates a system of the vectorized grid and ratios so that they can be operated on in future functions

#Below initializes reference dimensions of the wing geometry including area, chord, span, rotation locations, and velocity. It's put all into a 'ref' variable via the Reference() function
Sref = 30.0 # reference area
cref = 2.0  # reference chord
bref = 15.0 # reference span
rref = [0.50, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity (magnitude)
ref = Reference(Sref, cref, bref, rref, Vinf)

#Below initializes and compiles the freestream flow charactaristics. It references the aircraft's angle of attack, it's sideslip angle, and rotational velocities.
#It's put all into a 'fs' variable via the Freestream() function
alpha = 1.0*pi/180 # angle of attack
beta = 0.0 # sideslip angle, for simplicity is kept at zero, meaning the aircraft is moving directly into the freestream velcoity vector
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
fs = Freestream(Vinf, alpha, beta, Omega)

symmetric = false # Here we decide to not model our flow conditions using symmetry, but rather by mirroring the conditions across the x-z plane. If we were to use the symmetric keyword, we would get incorrect results on certain lateral stability derivatives

#This function is the root of the VLM analysis. It performs steady analysis, meaning a model of stable air flow, and mutates the system variable to an updated system variable that is referenced
#later to access all of the circulation distribution calculations. According to the documentation for this function, 'system' is in the place of the 'surfaces' argument, which makes sense as the previously defined system is
#panel grid and ratio vectors that define the surface of the wing. The naming convention of including an '!' following the function name indicates that the function is mutating one of the arguments.
#Here we are performing a near field anaylsis to find the force distribution across the panels
steady_analysis!(system, ref, fs; symmetric)

CF, CM = body_forces(system; frame=Wind()) #Here we assign the force (CF) and moment (CM) coefficients based on the analysis performed in the above function call. A near field analysis is needed to be performed before running this function (which is what we did in the steady_analysis! function). 
# Here we are using the Wind frame as opposed to the default body frame to return the CF and CM values

# extract aerodynamic forces
CD, CY, CL = CF #Here we assign variable names associated with the Force coefficients. They are respectively the Drag Coefficient, Side-force Coefficient (represents a force perpendicular to the freestream in the sideslip direction), and Lift Coefficient
Cl, Cm, Cn = CM #Here we assign variable names associated with the Moment coefficients. They are 3 ways in which a plane can rotate about axes, and are the rolling moment coefficient, ptiching moment coefficient, and yawing moment coefficient

CDiff = far_field_drag(system) # Here we compute a far field analysis to better analyze drag on the system. The near field analysis is often associated with numerical noise that can adversely affect the calculations

# calculate lifting line coefficients
cf, cm = lifting_line_coefficients(system; frame=Body()) # Now we perform analysis on the surface for lifting line coefficients. cf and cm are both vectors with a length equal to the number of surfaces (in this case 1), 
#and each element is an individual matrix with size (3, ns). It contains the x, y, and z direction force and moment coefficients (per unit span) for each spanwise segment. The lifting line is a 1D representation of the wing along the span.

dCFb, dCMb = body_derivatives(system) # Taken from the body frame

# traditional names for each body derivative (u - forward velocity,v - side velocity, w - upward veloicty, p - roll rate, q - pitch rate, r - yaw rate)
CXu, CYu, CZu = dCFb.u # How Forces change with forward velocity
CXv, CYv, CZv = dCFb.v # How forces change with side velocity
CXw, CYw, CZw = dCFb.w # How forces change with vertical velocity
CXp, CYp, CZp = dCFb.p # How forces change with roll rate
CXq, CYq, CZq = dCFb.q # How forces change with pitch rate
CXr, CYr, CZr = dCFb.r # How forces change with yaw rate
Clu, Cmu, Cnu = dCMb.u # How moments change with forward velocity
Clv, Cmv, Cnv = dCMb.v # HOw moments change with side velocity
Clw, Cmw, Cnw = dCMb.w # How moments change with upward velocity
Clp, Cmp, Cnp = dCMb.p # How moments change with roll rate
Clq, Cmq, Cnq = dCMb.q # How moments change with pitch rate
Clr, Cmr, Cnr = dCMb.r # How moments change with yaw rate

dCFs, dCMs = stability_derivatives(system) # Taken from the wind frame

# traditional names for each stability derivative (alpha - angle of attack, beta - sideslip angle, p - roll rate, q - pitch rate, r - yaw rate)
CDa, CYa, CLa = dCFs.alpha # How coefficients of drag, sideslip, and lift change with angle of attack
Cla, Cma, Cna = dCMs.alpha # How coefficients of rolling moment, pitching moment, and yawing moment change with angle of attack
CDb, CYb, CLb = dCFs.beta # How coefficients of drag, sideslip, and lift change with sideslip angle
Clb, Cmb, Cnb = dCMs.beta # How coefficients of rolling moment, pitching moment, and yawing moment change with sideslip angle
CDp, CYp, CLp = dCFs.p # How coefficients of drag, sideslip, and lift change with roll rate
Clp, Cmp, Cnp = dCMs.p # How coefficients of rolling moment, pitching moment, and yawing moment change with roll rate
CDq, CYq, CLq = dCFs.q # How coefficients of drag, sideslip, and lift change with pitch rate
Clq, Cmq, Cnq = dCMs.q # How coefficients of rolling moment, pitching moment, and yawing moment change with pitch rate
CDr, CYr, CLr = dCFs.r # How coefficients of drag, sideslip, and lift change with yaw rate
Clr, Cmr, Cnr = dCMs.r # How coefficients of rolling moment, pitching moment, and yawing moment change with yaw rate

write_vtk("simplewing", system) # This function takes everything we do above and writes it in a way that can be visualized using ParaView


# Begin plotting force and moment polars
alphas = range(-pi/12, 45 * pi/180, step = pi/360) # creates an array of alphas from 0 to pi/4 Radians that will be iterated over to find the force and moment polars
# Next 6 lines initializes arrays for each of the coefficients we're investigating
CL_array_alphas = Float64[]
CD_array_near = Float64[]
CD_array_far = Float64[]
Cl_array_alphas = Float64[]
Cm_array_alphas = Float64[]
Cn_array_alphas = Float64[]

for a in alphas # This for loop is what iterates through each alpha and creates a new system to analyze. Once analyzed, the results of coefficients are stored in the appropriate arrays before iterating to the next alpha value
    fs_alphas = Freestream(Vinf, a, beta, Omega)
    system_alphas = System(grids; ratios)
    steady_analysis!(system_alphas, ref, fs_alphas; symmetric)
    CD_far = far_field_drag(system_alphas)
    CF_alphas, CM_alphas = body_forces(system_alphas; frame=Wind())
    CD_alpha, CY_alpha, CL_alpha = CF_alphas
    Cl_alpha, Cm_alpha, Cn_alpha = CM_alphas
    push!(CL_array_alphas, CL_alpha) # Here we push the individual coefficient values to their relevant arrays. They update with each new alpha that is looped through.
    push!(CD_array_near, CD_alpha)
    push!(CD_array_far, CD_far)
    push!(Cl_array_alphas, Cl_alpha)
    push!(Cm_array_alphas, Cm_alpha)
    push!(Cn_array_alphas, Cn_alpha)
end


tick_labels = (-pi/12:π/12:π/4, [L"-\pi/12",L"0", L"\pi/12", L"\pi/6", L"\pi/4"]) # tick marks for x axis - here we start at 0, step in increments of pi/12, and max at pi/4 (or 45 degrees)

# Coefficient of Lift with varying Angle of Attack Plot
lift = plot(alphas,CL_array_alphas,color=:orange,linewidth=1.5,linestyle=:solid, label="Lift Coefficient",
    xlabel="Angle of Attack, \$\\alpha\$, in Radians",
    ylabel="Coefficient of Lift, \$C_L\$",
    xtickfont=font(12),
    ytickfont=font(10),
    xticks=tick_labels,
    dpi=600) # Various adjustments made to the plot to improve readability and asthetics.
    filename = "lift_Aoa.png" # Saved file name
    save_path = joinpath(@__DIR__, filename) # Specifies path to the current working folder
    savefig(lift, save_path) # Saves image
display(lift) # Displays plot.

# The rest of the plots below follow the same general pattern as above while occaisionally adding multiple lines to a graph using the plot!() function
# Plot displaying the relationship between far and near field induced drag and angle of attack
near_and_far_drag = plot(alphas, CD_array_far, color=:purple, linewidth=1.5, linestyle=:solid, label="Far-Field Induced Drag Coefficient",
        xlabel="Angle of Attack, \$\\alpha\$, in Radians",
        xtickfont=font(12),
        ytickfont=font(10),
        ylabel="Induced Drag Coefficient, \$C_D\$",
        xticks=tick_labels,
        dpi=600)
plot!(alphas, CD_array_near, color=:orange, linewidth=1.5, linestyle=:solid, label="Near-Field Induced Drag Coefficient")
filename = "near_far_drag_Aoa.png"
save_path = joinpath(@__DIR__, filename)
savefig(near_and_far_drag, save_path)
display(near_and_far_drag) 

# Moment Coefficients vs Angle of Attack Plot
roll = plot(alphas, Cl_array_alphas, linewidth=1.5, linestyle=:solid, color=:orange, label="Rolling Moment Coefficient", legendfontsize=10,
    xlabel="Angle of Attack, \$\\alpha\$, in Radians",
    ylabel="Moment Coefficient",
    xtickfont=font(14),
    ytickfont=font(10),
    xticks=tick_labels,
    dpi=600)
filename = "roll_Aoa.png"
save_path = joinpath(@__DIR__, filename)
savefig(roll, save_path)
display(roll)

pitch = plot(alphas, Cm_array_alphas, color=:green, linewidth=1.5, linestyle=:solid, label="Pitching Moment Coefficient",legend=:topright, legendfontsize=10,
    xlabel="Angle of Attack, \$\\alpha\$, in Radians",
    ylabel="Moment Coefficient",
    xtickfont=font(14),
    ytickfont=font(10),
    xticks=tick_labels,
    dpi=600)
filename = "pitch_Aoa.png"
save_path = joinpath(@__DIR__, filename)
savefig(pitch, save_path)
display(pitch)

yaw = plot(alphas, Cn_array_alphas, linewidth=1.5, linestlye=:solid, color=:purple, label="Yaw Moment Coefficient", legendfontsize=9,
    xlabel="Angle of Attack, \$\\alpha\$, in Radians",
    ylabel="Moment Coefficient",
    xtickfont=font(14),
    ytickfont=font(8),
    xticks=tick_labels,
    dpi=600)
filename = "yaw_Aoa.png"
save_path = joinpath(@__DIR__, filename)
savefig(yaw, save_path)
display(yaw)

#Lift vs Drag (Drag Polar)
drag_polar = plot(CL_array_alphas, CD_array_far, linewidth=1.5, linestyle=:solid, color=:orange, label="Drag Polar",
    xlabel="Coefficient of Lift, \$C_L\$",
    ylabel="Coefficient of Drag, \$C_D\$ (Far-Field)",
    ytickfont=font(12),
    xtickfont=font(8),
    dpi=600)
filename = "drag_polar.png"
save_path = joinpath(@__DIR__, filename)
savefig(drag_polar, save_path)
display(drag_polar)

#Lift/Drag vs Angle of Attack
lift_drag_alphas = plot(alphas[30:end], CL_array_alphas[30:end] ./ CD_array_far[30:end], color=:orange, linewidth=1.5, linestyle=:solid, label="Lift/Drag Ratio",
    xlabel="Angle of Attack, \$\\alpha\$, in Radians",
    ylabel="Lift/Drag Coefficient Ratio",
    xtickfont=font(12),
    ytickfont=font(10),
    xticks=(0:π/12:π/4,[L"0", L"\pi/12", L"\pi/6", L"\pi/4"]),
    dpi=600)
# plot!(alphas,CD_array_far,color=:purple,linewidth=1.5,linestyle=:solid, label="Far Field Drag Coefficient")
filename = "lift_drag_Aoa.png"
save_path = joinpath(@__DIR__, filename)
savefig(lift_drag_alphas, save_path)
display(lift_drag_alphas)


println("Lift Curve Slope = $(CLa)") # The variable CLa is defined as the change in coefficient of lift with angle of attack, so we just print it here