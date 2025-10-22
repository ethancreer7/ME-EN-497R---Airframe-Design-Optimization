using VortexLattice
using Plots
using LaTeXStrings
using Statistics

# Below is a list of initialized wing design parameters. It includes data for the leading edge, Chord Length, Twist Angle, and a function for the camberline
xle = [0.0, 0.4] # leading edge x-position --- x direction is in the direction of the freestream vector. 0.0m at the root (origin), and 0.4 at the tip (provides a swept back wing)
yle = [0.0, 7.5] # leading edge y-position --- y direction is in the spanwise direction of the aircraft. 0.0m at the root (origin), and 7.5m at the tip
zle = [0.0, 0.0] # leading edge z-position --- z direction is the vertical direction. 0.0m at the root (origin), and 0.0m at the tip. This indicates the wing has no designed dihedral
chord = [2.2, 1.8] # chord length goes from 2.2m to 1.8m over the span of the wing (again indicating a swept wing design)
theta = [0, 0] # twist (in Radians) - constant over the span of the wing
phi = [0.0, 0.0] # section rotation about the x-axis (dihedral is zero across the span of the wing)
fc = fill((xc) -> 0, 2) # camberline function for each section (y/c = f(x/c)) - Anonymous function dictates camberline is always flat (0). Assigns this as a vector for the 2 points specified above.


# The below section begins the discretization process across the wing. In this example, we divide the wing into 72 panels and develop
ns = 120 # number of spanwise panels
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
Vinf = 100 # reference velocity (magnitude)
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
#Here we are performing a near field analysis to find the force distribution across the panels
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

save_path = joinpath(@__DIR__, "Figures", "simplewing")
write_vtk(save_path, system) # This function takes everything we do above and writes it in a way that can be visualized using ParaView

# All of the above is from the getting started guid of the Julia VLM package


function getLift(rho, V, dA, cf) # returns the lifting force for each panel that is passed into it
    return 0.5 .* rho .* V .* V .* dA .* cf
end

function getArea(i, grid) # returns the trapezoidal area of each panel
    chord1 = maximum(grid[1,:,i]) - minimum(grid[1,:,i])
    chord2 = maximum(grid[1,:,i+1]) - minimum(grid[1,:,i+1])
    dy = abs(mean(grid[2,:,i + 1]) - mean(grid[2, :, i]))
    return 1/2 * (chord1 + chord2) * dy
end
# Add arrows function to a specified plot
function addArrows(xvals, yvals, plot_name, color_name, ends, lower) # Adds arrows to individual wing plots to aid in interpretability
    for i in ends:(Int(length(xvals)/10)-ends)
        x_pos = xvals[i*10]
        y_pos = yvals[i*10]
        #Add arrows
        plot!(plot_name, [x_pos, x_pos], [0, y_pos - lower], # Plots the arrow to just below the function line (so no overlapping)
            arrow=arrow(),
            color=color_name,
            alpha=0.3,
            linewidth=1.2,
            label="")
    end
end

rho = 1.225 # standard air density at sea level (units of kg/m^3) - used in getLift function

# Begin tapered chord wing system operations
cf_w, cm_w = lifting_line_coefficients(system; frame=Wind()) # Here we use the Wind() frame so that we don't have to do trigonometry that would be necessary when using the Body() frame to determine the true lift direction (perpindicular to the freestream)
c_lift_vals = cf_w[1][3,:] # Obtain the coefficient of lift values at spanwise panel section
lift_vals_array = Float64[]
iter = 1
for i in c_lift_vals
    dA = getArea(iter, grid)
    lift_val = getLift(rho, Vinf, dA, i) # calls getLift function to convert the c_l values into L' values
    push!(lift_vals_array,lift_val) # pushes each lift_val into its appropriate array
    iter += 1
end


#Obtain necessary discretized geometry for the tapered chord wing. Here we calculate the mid span point of each spanwise panel section.
span_grid = grid[2,1,:]
span_panel_centers = [(span_grid[i] + span_grid[i+1]) / 2.0 for i in 1:length(lift_vals_array)]

tick_marks = (-7.5:1.5:7.5)
cumulative = plot(span_panel_centers, lift_vals_array, linewidth=1.5, linestyle =:solid, color =:orange, label = "Tapered Chord",
    xlabel = "Span Position (m)",
    ylabel = "Panel Lifting Force (N)",
    xticks = tick_marks,
    legendfontsize = 7,
    dpi = 600) # Create a cumulative Lifting force plot
cumulative_cl = plot(span_panel_centers, c_lift_vals, linewidth = 1.5, linestyle =:solid, color =:orange, label = "Tapered Chord",
    xlabel = "Span Position (m)",
    ylabel = L"Coefficient of Lift $c_l$",
    xticks = tick_marks,
    dpi = 600,
    legendfontsize = 7,
    legend=:topright) # Create a cumulative coefficient of lift plot


tapered_L = plot(span_panel_centers,lift_vals_array, linewidth = 1.5, linestyle =:solid, color =:orange, label = "Tapered Chord",
    xlabel = "Span Position (m)",
    ylabel = "Panel Lifting Force (N)",
    xticks = tick_marks,
    dpi = 600) # Create an individual lift plot for the tapered wing
addArrows(span_panel_centers, lift_vals_array, tapered_L, :darkorange2, 3, 5) # adds the distribution arrows
filename = "tapered_chord_distribution_Lift.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(tapered_L, save_path)
display(tapered_L)

tapered_cl = plot(span_panel_centers, c_lift_vals, linewidth = 1.5, linestyle =:solid, color =:orange, label = "Tapered Chord",
    xlabel = "Span Position (m)",
    ylabel = L"Coefficient of Lift $c_l$",
    xticks = tick_marks,
    dpi = 600,
    legend=:topright) # Create an individual coefficient of lift plot for the tapered wing
filename = "tapered_chord_distribution_cL.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(tapered_cl, save_path)
display(tapered_cl)


#Create constant chord wing geometry and its grid
constant_chord = [2.0, 2.0] # This selected chord length ensures the same area is used (30 as the span length times the average chord length was used in our tapered chord system)
xle_constant = [0.0,0.0]
constant_grid, constant_ratio = wing_to_grid(xle_constant, yle, zle, constant_chord, theta, phi, ns, nc;
fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)

constant_grids = [constant_grid]
constant_ratios = [constant_ratio]
constant_system = System(constant_grids; ratios = constant_ratios) # Create the constant chord system for analysis
steady_analysis!(constant_system,ref,fs)
const_CF, const_CM = body_forces(constant_system; frame=Wind()) # for checking to see if sum and actual total are close
const_CD, const_CY, const_CL = const_CF  # for checking to see if sum and actual total are close

const_cf, const_cm = lifting_line_coefficients(constant_system; frame=Wind()) # Again using the Wind() frame for better analysis
const_lift_vals = const_cf[1][3,:]
const_lift_vals_array = Float64[]
const_iter = 1
for i in const_lift_vals
    dA = getArea(const_iter, constant_grid)
    const_lift_val = getLift(rho, Vinf, dA, i) # Again calling the getLift value to obtain the Panel Lifting force for the constant chord wing
    push!(const_lift_vals_array, const_lift_val)
    const_iter += 1
end

# Plotting of the constant chord wing, similar to the tapered chord wing plotting
const_span_grid = constant_grid[2,1,:]
const_span_panel_centers = [(const_span_grid[i] + const_span_grid[i+1]) / 2.0 for i in 1:length(const_lift_vals_array)]
constant_L = plot(const_span_panel_centers, const_lift_vals_array, linestyle =:solid, color =:cyan4, label = "Constant Chord",
    xlabel = "Span Position (m)",
    ylabel = "Panel Lifting Force (N)",
    xticks = tick_marks,
    dpi = 600)
addArrows(const_span_panel_centers, const_lift_vals_array,constant_L, :darkcyan, 3, 5)
filename = "constant_chord_dist_L.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(constant_L, save_path)
display(constant_L)

constant_cl = plot(const_span_panel_centers, const_lift_vals, linewidth = 1.5, linestyle =:solid, color =:cyan4, label = "Constant Chord",
    xlabel = "Span Position (m)",
    ylabel = L"Coefficient of Lift $c_l$",
    xticks = tick_marks,
    dpi = 600,
    legend=:topright)
filename = "const_chord_dist_cl.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(constant_cl, save_path)
display(constant_cl)


# Begin defining the parameters of an elliptic wing:
stations = 20
half_span = bref/2.0
root_chord = 4 * Sref / (pi * bref) # Calculate the root chord to meet the proper area constraint


xle_elliptic = zeros(stations)
yle_elliptic = collect(LinRange(0.0, half_span, stations))
zle_elliptic = zeros(stations)

function elliptic_chord(y,root,half) # Uses the equation of an ellipse to calculate an appropriate chord length across the span of the aircraft
    return root * sqrt(1-(y/half)^2)
end

# Defines the elliptic wing parameters and creates a grid for the elliptic wing
chord_elliptic = elliptic_chord.(yle_elliptic,root_chord,half_span)
theta_elliptic = zeros(stations)
phi_elliptic = zeros(stations)
fc_elliptic = fill((xc) -> 0, stations)

elliptic_grid, elliptic_ratio = wing_to_grid(xle_elliptic,yle_elliptic,zle_elliptic,chord_elliptic,theta_elliptic,phi_elliptic,ns,nc;
fc = fc_elliptic, spacing_s = spacing_s, spacing_c = spacing_c, mirror = true)

elliptic_grids = [elliptic_grid]
elliptic_ratios = [elliptic_ratio]
elliptic_system = System(elliptic_grids; ratios=elliptic_ratios) # Create elliptic wing system

elliptic_cref = (8/(3pi)) * root_chord # We use this as the reference chord length in our reference variable (the mean aerodynamic chord for an elliptic wing)
elliptic_rref = [elliptic_cref/4,0.0,0.0]
elliptic_ref = Reference(Sref, elliptic_cref, bref, elliptic_rref, Vinf)

steady_analysis!(elliptic_system, elliptic_ref, fs; symmetric)

elliptic_CF, elliptic_CM = body_forces(elliptic_system; frame=Wind()) 

elliptic_CD, elliptic_CY, elliptic_CL = elliptic_CF

elliptic_cf, elliptic_cm = lifting_line_coefficients(elliptic_system; frame = Wind()) # Again using the Wind() frame
elliptic_coef_lift_vals = elliptic_cf[1][3,:]
elliptic_lift_vals_array = Float64[]
elliptic_iter = 1
for i in elliptic_coef_lift_vals
    dA = getArea(elliptic_iter, elliptic_grid)
    lift_val = getLift(rho, Vinf, dA, i) #Again calling the getLift value to obtain the Panel Lifting force for the elliptic wing
    push!(elliptic_lift_vals_array,lift_val)
    elliptic_iter += 1
end

# Plotting of the elliptic chord wing, similar to the tapered and constant chord wing plotting
elliptic_wing_grid = elliptic_grid[2,1,:]
elliptic_wing_centers = [(elliptic_wing_grid[i] + elliptic_wing_grid[i+1]) / 2.0 for i in 1:length(elliptic_coef_lift_vals)]
elliptic_wing_L = plot(elliptic_wing_centers, elliptic_lift_vals_array, linestyle=:solid, color=:brown, label = "Elliptic Wing",
    xlabel = "Span Position (m)",
    ylabel = "Panel Lifting Force (N)",
    xticks = tick_marks,
    dpi = 600)
addArrows(elliptic_wing_centers, elliptic_lift_vals_array, elliptic_wing_L, :brown, 3, 4)
filename = "elliptic_wing_L.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(elliptic_wing_L, save_path)
display(elliptic_wing_L)

elliptic_wing_cl = plot(elliptic_wing_centers[begin+10:end-10], elliptic_coef_lift_vals[begin+10:end-10], linewidth = 1.5, linestyle =:solid,color=:brown, label = "Elliptic Wing",
    xlabel = "Span Position (m)",
    ylabel = L"Coefficient of Lift $c_l$",
    xticks = tick_marks,
    dpi = 600)
filename = "elliptic_wing_cl.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(elliptic_wing_cl, save_path)
display(elliptic_wing_cl)


# Begin plotting ideal elliptic lift distribution
function elliptic(gamma, y, b) # function that plots an ellipse with a provided maximum
return gamma * sqrt(1 - ((y/(b/2))^2))
end
discretized_span = -7.5:(15/120):7.5 # Gets x values for the ellipse
ideal_elliptic_lift_array = Float64[]
maxes = [maximum(const_lift_vals_array), maximum(lift_vals_array), maximum(elliptic_lift_vals_array)] # Collects the maxima of each wing's panel lifting force
gamma = mean(maxes) # Averages the maxima from the previous line to use as the maximum for the ideal elliptic lift distribution
for i in discretized_span
    push!(ideal_elliptic_lift_array, elliptic(gamma,i, bref)) # Calculates the y-axis values for the ellipse
end

# Plots the ideal elliptic lift distribution individually
elliptic_dist = plot(discretized_span, ideal_elliptic_lift_array, linewidth = 1.5, linestyle =:solid, color =:mediumpurple, label = "Elliptic Lift Distribution",
    xlabel = "Span Position (m)",
    ylabel = "Panel Lifting Force (N)",
    xticks = tick_marks,
    dpi = 600)
filename = "elliptic_dist.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(elliptic_dist, save_path)
display(elliptic_dist)

# Create ideal c_l
ideal_cl_array = fill(elliptic_CL, length(discretized_span)) # Creates a constant ideal elliptic lift coefficient array with a length of the discretized_span array.


# Create a cumulative plot for comparisons of coefficient of lift distribution
plot!(cumulative_cl, const_span_panel_centers, const_lift_vals, linewidth = 1.5, linestyle =:solid, color =:cyan4, label = "Constant Chord")
plot!(cumulative_cl, elliptic_wing_centers[begin+10:end-10], elliptic_coef_lift_vals[begin+10:end-10],linewidth = 1.5, linestyle =:solid, color =:brown, label = "Elliptic Wing") # Cut off the numerical artifacts at the first and last 10 points on the curve.
plot!(cumulative_cl, discretized_span, ideal_cl_array, linewidth = 1.5, linestyle =:solid, color =:mediumpurple, label = "Ideal Coefficient of Lift")
filename = "cumulative_cl.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(cumulative_cl, save_path)
display(cumulative_cl)

#Create a cumulative plot for comparisons of Lift distributions
plot!(cumulative, const_span_panel_centers, const_lift_vals_array, linewidth=1.5, linestyle =:solid, color =:cyan4, label = "Constant Chord")
plot!(cumulative, elliptic_wing_centers, elliptic_lift_vals_array,  linestyle=:solid, color=:brown, label = "Elliptic Wing")
plot!(cumulative, discretized_span, ideal_elliptic_lift_array, linewidth = 1.5, linestyle =:solid, color =:mediumpurple, label = "Ideal Lift")

filename = "cumulative_lift.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(cumulative, save_path)
display(cumulative)


#Begin plotting the wing's planform shape
function plot_wing_edge(xle, yle, chord, color, plot; label="")
    xte = xle .+ chord # Calculates the wings trailing edge as the leading edge position plus the chord length

    # Create full span geometry by mirroring the right half
    full_yle = vcat(-reverse(yle[2:end]), yle)
    full_xle = vcat(reverse(xle[2:end]), xle)
    full_xte = vcat(reverse(xte[2:end]), xte)

    # Creates the full continuous shape
    x_coords = vcat(full_yle, reverse(full_yle))
    y_coords = vcat(full_xle, reverse(full_xte))

    push!(x_coords, x_coords[1])
    push!(y_coords, y_coords[1])


    plot!(plot, x_coords, y_coords, linewidth=1.5, label=label, color=color) # Plots the outline
end

# Set up cumulative planform wing plot
planform_wing = plot(xlabel="Spanwise Direction (m)", ylabel="Freestream Direction (m)",
    xlims=(-10,10),
    ylims=(-5,5),
    aspect_ratio=:equal,
    yflip = true)

#Tapered Wing Outline
plot_wing_edge(xle, yle, chord, :orange, planform_wing; label="Tapered Wing")

#Constant Chord Wing Outline
plot_wing_edge(xle_constant, yle, constant_chord, :cyan4, planform_wing; label="Constant Chord")


#Elliptic Wing Outline
plot_wing_edge(xle_elliptic, yle_elliptic, chord_elliptic, :brown, planform_wing; label="Elliptic Wing")

filename = "cumulative_planform.png"
save_path = joinpath(@__DIR__, "figures", filename)
savefig(planform_wing, save_path)
display(planform_wing)

# Write vtk files for visualization in ParaView
write_vtk("simplewing", system)
write_vtk("constant_chord", constant_system)
write_vtk("elliptic_chord", elliptic_system)