#=##############################################################################
# DESCRIPTION
    45° swept-back wing at an angle of attack of 4.2°. This wing has an aspect
    ratio of 5.0, a RAE 101 airfoil section with 12% thickness, and no dihedral,
    twist, nor taper. This test case matches the experimental setup of Weber,
    J., and Brebner, G., “Low-Speed Tests on 45-deg Swept-Back Wings, Part I,”
    Tech. rep., 1951. The same case is used in a VLM calculation in Bertin's
    Aerodynamics for Engineers, Example 7.2, pp. 343.

# AUTHORSHIP
  * Author          : Eduardo J. Alvarez (edoalvarez.com)
  * Email           : Edo.AlvarezR@gmail.com
  * Created         : Feb 2023
  * Last updated    : Feb 2023
  * License         : MIT
=###############################################################################

import FLOWUnsteady as uns
import FLOWVLM as vlm
using DataFrames, CSV
using Statistics
using Plots

results = DataFrame(
n_val = Float64[], CL = Float64[], CD = Float64[])
prev_CL, prev_CD = nothing, nothing
tolerance = 0.1/100


function append_result!(df::DataFrame, tuple)
  result = Dict(
  "n_val" => tuple[1],
  "CL" => tuple[2], 
  "CD" => tuple[3]
  )
  push!(df, result)
end


for (count, i) in enumerate(range(10, 150, step = 10))
  run_name        = "wing$(i)"            # Name of this simulation

  save_path       = joinpath(@__DIR__, "wings", run_name)               # Where to save this simulation
  paraview        = true                                                # Whether to visualize with Paraview


  # ----------------- SIMULATION PARAMETERS --------------------------------------
  AOA             = 4.2                       # (deg) angle of attack
  magVinf         = 49.7                      # (m/s) freestream velocity
  rho             = 0.93                      # (kg/m^3) air density
  qinf            = 0.5*rho*magVinf^2         # (Pa) static pressure

  Vinf(X, t)      = magVinf*[cosd(AOA), 0.0, sind(AOA)]  # Freestream function


  # ----------------- GEOMETRY PARAMETERS ----------------------------------------
  b               = 2.489                     # (m) span length
  ar              = 5.0                       # Aspect ratio b/c_tip
  tr              = 1.0                       # Taper ratio c_tip/c_root
  twist_root      = 0.0                       # (deg) twist at root
  twist_tip       = 0.0                       # (deg) twist at tip
  lambda          = 45.0                      # (deg) sweep
  gamma           = 0.0                       # (deg) dihedral

  # Discretization
  n               = i                        # Number of spanwise elements per side
  r               = 10.0                      # Geometric expansion of elements
  central         = false                     # Whether expansion is central

  # NOTE: A geometric expansion of 10 that is not central means that the last
  #       element is 10 times wider than the first element. If central, the
  #       middle element is 10 times wider than the peripheral elements.

  # ----------------- SOLVER PARAMETERS ------------------------------------------
  # Time parameters
  wakelength      = 2.75*b                    # (m) length of wake to be resolved
  ttot            = wakelength/magVinf        # (s) total simulation time
  nsteps          = 20                         # Number of time steps

  # VLM and VPM parameters
  p_per_step      = 1                         # Number of particle sheds per time step

  lambda_vpm      = 2.0                       # VPM core overlap
  sigma_vpm_overwrite = lambda_vpm * magVinf * (ttot/nsteps)/p_per_step # Smoothing core size
  sigma_vlm_solver= -1                        # VLM-on-VLM smoothing radius (deactivated with <0)
  sigma_vlm_surf  = 0.05*b                    # VLM-on-VPM smoothing radius

  shed_starting   = true                      # Whether to shed starting vortex
  vlm_rlx         = 0.7                       # VLM relaxation




  # ----------------- 1) VEHICLE DEFINITION --------------------------------------
  println("Generating geometry...")

  # Generate wing
  wing = vlm.simpleWing(b, ar, tr, twist_root, lambda, gamma;
                          twist_tip=twist_tip, n=n, r=r, central=central);

  # NOTE: `FLOWVLM.simpleWing` is an auxiliary function in FLOWVLM for creating a
  #       VLM wing with constant sweep, dihedral, and taper ratio, and a linear
  #       twist between the root and the wing tips

  println("Generating vehicle...")

  # Generate vehicle
  system = vlm.WingSystem()                   # System of all FLOWVLM objects
  vlm.addwing(system, "Wing", wing)

  vlm_system = system;                        # System solved through VLM solver
  wake_system = system;                       # System that will shed a VPM wake

  vehicle = uns.VLMVehicle(   system;
                              vlm_system=vlm_system,
                              wake_system=wake_system
                          );

  # NOTE: `FLOWUnsteady.VLMVehicle` creates a vehicle made out of multiple VLM and
  #       rotor subsystems. The argument `system` represents the vehicle as a
  #       whole which will be translated and rotated with the kinematics
  #       prescribed by the maneuver. The subsystem `vlm_system` will be solved
  #       with the VLM solver. The subsystems `rotor_systems` are solved with
  #       blade elements (none in this case). The subsystem `wake_system` will
  #       shed the VPM wake.




  # ------------- 2) MANEUVER DEFINITION -----------------------------------------

  Vvehicle(t) = zeros(3)                      # Translational velocity of vehicle over time
  anglevehicle(t) = zeros(3)                  # Angle of the vehicle over time

  angle = ()                                  # Angle of each tilting system (none)
  RPM = ()                                    # RPM of each rotor system (none)

  maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)

  # NOTE: `FLOWUnsteady.KinematicManeuver` defines a maneuver with prescribed
  #       kinematics. `Vvehicle` defines the velocity of the vehicle (a vector)
  #       over time. `anglevehicle` defines the attitude of the vehicle over time
  #       (a vector with inclination angles with respect to each axis of the
  #       global coordinate system). `angle` defines the tilting angle of each
  #       tilting system over time (none in this case). `RPM` defines the RPM of
  #       each rotor system over time (none in this case).
  #       Each of these functions receives a nondimensional time `t`, which is the
  #       simulation time normalized by the total time `ttot`, from 0 to
  #       1, beginning to end of simulation. They all return a nondimensional
  #       output that is then scaled up by either a reference velocity (`Vref`) or
  #       a reference RPM (`RPMref`). Defining the kinematics and controls of the
  #       maneuver in this way allows the user to have more control over how fast
  #       to perform the maneuver, since the total time, reference velocity and
  #       RPM are then defined in the simulation parameters shown below.




  # ------------- 3) SIMULATION DEFINITION ---------------------------------------

  Vref = 0.0                                  # Reference velocity to scale maneuver by
  RPMref = 0.0                                # Reference RPM to scale maneuver by
  Vinit = Vref*Vvehicle(0)                    # Initial vehicle velocity
  Winit = pi/180*(anglevehicle(1e-6) - anglevehicle(0))/(1e-6*ttot)  # Initial angular velocity

                                              # Maximum number of particles
  max_particles = (nsteps+1)*(vlm.get_m(vehicle.vlm_system)*(p_per_step+1) + p_per_step)

  simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                      Vinit=Vinit, Winit=Winit);




  # ------------- 4) MONITORS DEFINITIONS ----------------------------------------

  # NOTE: Monitors are functions that are called at every time step to perform
  #       some secondary computation after the solution of that time step has been
  #       obtained. In this case, the wing monitor uses the circulation and
  #       induced velocities to compute aerodynamic forces and decompose them
  #       into lift and drag. The monitor then plots these forces at every time
  #       step while also saving them under a CSV file in the simulation folder.

  # Generate function that calculates aerodynamic forces
  # NOTE: We exclude skin friction since we want to compare to the experimental
  #       data reported in Weber 1951 that was measured with pressure taps
  calc_aerodynamicforce_fun = uns.generate_calc_aerodynamicforce(;
                                      add_parasiticdrag=true,
                                      add_skinfriction=false,
                                      airfoilpolar="xf-rae101-il-1000000.csv"
                                      )

  D_dir = [cosd(AOA), 0.0, sind(AOA)]         # Direction of drag
  L_dir = uns.cross(D_dir, [0,1,0])           # Direction of lift

  figs, figaxs = [], []                       # Figures generated by monitor

  # Generate wing monitor
  monitor_wing = uns.generate_monitor_wing(wing, Vinf, b, ar,
                                              rho, qinf, nsteps;
                                              calc_aerodynamicforce_fun=calc_aerodynamicforce_fun,
                                              L_dir=L_dir,
                                              D_dir=D_dir,
                                              out_figs=figs,
                                              out_figaxs=figaxs,
                                              save_path=save_path,
                                              run_name=run_name,
                                              disp_plot = false,
                                              figname="wing monitor",
                                              );




  # ------------- 5) RUN SIMULATION ----------------------------------------------
  println("Running simulation...")

  uns.run_simulation(simulation, nsteps;
                      # ----- SIMULATION OPTIONS -------------
                      Vinf=Vinf,
                      rho=rho,
                      # ----- SOLVERS OPTIONS ----------------
                      p_per_step=p_per_step,
                      max_particles=max_particles,
                      sigma_vlm_solver=sigma_vlm_solver,
                      sigma_vlm_surf=sigma_vlm_surf,
                      sigma_rotor_surf=sigma_vlm_surf,
                      sigma_vpm_overwrite=sigma_vpm_overwrite,
                      shed_starting=shed_starting,
                      vlm_rlx=vlm_rlx,
                      extra_runtime_function=monitor_wing,
                      # ----- OUTPUT OPTIONS ------------------
                      save_path=save_path,
                      run_name=run_name
                      );
  path = joinpath(@__DIR__, "wings", run_name, run_name * "_convergence.csv")                    
  result_df = CSV.read(path, DataFrame)
  CL_mean = mean(result_df[end-3:end, :CL])
  CD_mean = mean(result_df[end-3:end, :CD])
  tuple = [i, CL_mean, CD_mean]
  append_result!(results, tuple)

  if prev_CL !== nothing && prev_CD !== nothing
    percent_change_CL = abs(CL_mean - prev_CL) / abs(prev_CL)
    percent_change_CD = abs(CD_mean - prev_CD) / abs(prev_CD)
    println("Percent Change in CL = $(round(percent_change_CL*100, digits=5))%")
    println("Percent Change in CD = $(round(percent_change_CD*100, digits=5))%")

    if percent_change_CL < tolerance && percent_change_CD < tolerance
      println("Converged at n = $n (Percent Change in CL = $(percent_change_CL*100)%, Percent Change in CD = $(percent_change_CD*100)%)")
      println("Final averaged values: CL = $CL_mean, CD = $CD_mean at n = $n")
      break
    end
  end
  prev_CL, prev_CD = CL_mean, CD_mean
end

convergence_plot = plot(results[!, "n_val"], results[!, "CL"], linewidth = 1.5, linestyle=:solid, color = :orange, label = "CL Value", ylabel= "CL", xlabel = "n_value")
plot!(twinx(), results[!, "n_val"], results[!, "CD"], linewidth = 1.5, linestyle = :solid, color = :purple, label = "CD Value", ylabel = "CD Value")
display(convergence_plot)

direct_path = joinpath(@__DIR__, "Figures")
mkpath(direct_path)
save_path_fig = joinpath(direct_path, "convergence_plot.png")
savefig(convergence_plot, save_path_fig)