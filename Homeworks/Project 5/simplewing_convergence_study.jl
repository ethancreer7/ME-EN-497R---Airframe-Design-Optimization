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
prev_cl, prev_cd = nothing, nothing
tolerance = 0.25/100


"""
  function get_cv(df::DataFrame)

Returns the coefficient of variation value for the given dataframe for both the coefficient of lift and the coefficient of drag
"""
function get_cv(df::DataFrame)
  cl_tail_vals = df[end-4:end, :CL]
  cd_tail_vals = df[end-4:end, :CD]

  cl_sigma = std(cl_tail_vals)
  cl_mean = mean(cl_tail_vals)
  cl_CV = cl_sigma / cl_mean
  cd_sigma = std(cd_tail_vals)
  cd_mean = mean(cd_tail_vals)
  cd_CV = cd_sigma / cd_mean

  return cl_CV, cd_CV
end

"""
  function printstats(cl_CV, cd_CV, variable_name, variable_val)

Returns print statements summarizing the convergence
"""
function printstats(cl_CV, cd_CV, variable_name, variable_val)
  println("cl_CV = $cl_CV")
  println("cd_CV = $cd_CV")
  println("$variable_name = $variable_val")
  println("Converged within .3% variation")
end

mkpath(joinpath(@__DIR__, "wakelength"))
mkpath(joinpath(@__DIR__, "nsteps"))
mkpath(joinpath(@__DIR__, "n_sweep"))
mkpath(joinpath(@__DIR__, "p_per_sweep"))
mkpath(joinpath(@__DIR__, "verification"))


# Define Global Variables for referencing converged values
global wakelength_converged = nothing
global nsteps_converged = nothing
global n_converged = nothing
global p_per_step_min = nothing
global p_per_step = nothing

# --------------- Wakelength Convergence Study ----------------
for (count, i) in enumerate(range(0.75, 6.5, step = 0.25))
  run_name        = "wakelength$(round(i*2.489, digits=3))"             # Name of this simulation

  save_path       = joinpath(@__DIR__, "wakelength", run_name)          # Where to save this simulation
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
  n               = 20                        # Number of spanwise elements per side
  r               = 10.0                      # Geometric expansion of elements
  central         = false                     # Whether expansion is central

  # NOTE: A geometric expansion of 10 that is not central means that the last
  #       element is 10 times wider than the first element. If central, the
  #       middle element is 10 times wider than the peripheral elements.

  # ----------------- SOLVER PARAMETERS ------------------------------------------
  # Time parameters
  wakelength      = i*b                       # (m) length of wake to be resolved
  ttot            = wakelength/magVinf        # (s) total simulation time
  nsteps          = 20                        # Number of time steps

  # VLM and VPM parameters
  lambda_vpm      = 2.0                       # VPM core overlap
  # p_per_step      = 2                         # Number of particle sheds per time step

  p_per_step = Int(ceil( (2 * n * lambda_vpm * magVinf * ttot) / (nsteps * b)))

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
  
    
  path = joinpath(@__DIR__, "wakelength", run_name, run_name * "_convergence.csv")                    
  result_df = CSV.read(path, DataFrame)                    
  cl_CV, cd_CV = get_cv(result_df)

  if cl_CV < 0.002 && cd_CV < 0.002
    printstats(cl_CV, cd_CV, "Wakelength", wakelength)
    global wakelength_converged = wakelength
    break
  end  
end


# --------------- nsteps Convergence Study --------------------
for (count, i) in enumerate(range(10, 120, step = 10))
  run_name        = "nsteps$(i)"            # Name of this simulation

  save_path       = joinpath(@__DIR__, "nsteps", run_name)              # Where to save this simulation
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
  n               = 70                        # Number of spanwise elements per side
  r               = 10.0                      # Geometric expansion of elements
  central         = false                     # Whether expansion is central

  
  # ----------------- SOLVER PARAMETERS ------------------------------------------
  # Time parameters
  wakelength      = wakelength_converged      # (m) length of wake to be resolved -- Referenced from the previous convergence block
  ttot            = wakelength/magVinf        # (s) total simulation time
  nsteps          = i                         # Number of time steps

  # VLM and VPM parameters
  lambda_vpm      = 2.0                       # VPM core overlap
  p_per_step = Int(ceil( (2 * n * lambda_vpm * magVinf * ttot) / (nsteps * b)))

  
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



  # ------------- 2) MANEUVER DEFINITION -----------------------------------------

  Vvehicle(t) = zeros(3)                      # Translational velocity of vehicle over time
  anglevehicle(t) = zeros(3)                  # Angle of the vehicle over time

  angle = ()                                  # Angle of each tilting system (none)
  RPM = ()                                    # RPM of each rotor system (none)

  maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)



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
  path = joinpath(@__DIR__, "nsteps", run_name, run_name * "_convergence.csv")                    
  result_df = CSV.read(path, DataFrame)
  cl_CV, cd_CV = get_cv(result_df)

  if cl_CV < 0.002 && cd_CV < 0.002
    printstats(cl_CV, cd_CV, "nsteps", nsteps)
    global nsteps_converged = nsteps
    break
  end

end


# --------------- n Convergence Study -------------------------
for (count, i) in enumerate(range(10, 150, step = 10))
  run_name        = "n_sweep$(i)"            # Name of this simulation

  save_path       = joinpath(@__DIR__, "n_sweep", run_name)             # Where to save this simulation
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
  n               = i                         # Number of spanwise elements per side
  r               = 10.0                      # Geometric expansion of elements
  central         = false                     # Whether expansion is central

  # ----------------- SOLVER PARAMETERS ------------------------------------------
  # Time parameters
  wakelength      = wakelength_converged                    # (m) length of wake to be resolved
  ttot            = wakelength/magVinf                      # (s) total simulation time
  nsteps          = nsteps_converged                        # Number of time steps

  # VLM and VPM parameters
   lambda_vpm      = 2.0                                     # VPM core overlap
  p_per_step = Int(ceil( (2 * n * lambda_vpm * magVinf * ttot) / (nsteps_converged * b)))


  sigma_vpm_overwrite = lambda_vpm * magVinf * (ttot/nsteps)/p_per_step # Smoothing core size
  sigma_vlm_solver= -1                                      # VLM-on-VLM smoothing radius (deactivated with <0)
  sigma_vlm_surf  = 0.05*b                                  # VLM-on-VPM smoothing radius

  shed_starting   = true                                    # Whether to shed starting vortex
  vlm_rlx         = 0.7                                     # VLM relaxation




  # ----------------- 1) VEHICLE DEFINITION --------------------------------------
  println("Generating geometry...")

  # Generate wing
  wing = vlm.simpleWing(b, ar, tr, twist_root, lambda, gamma;
                          twist_tip=twist_tip, n=n, r=r, central=central);


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




  # ------------- 2) MANEUVER DEFINITION -----------------------------------------

  Vvehicle(t) = zeros(3)                      # Translational velocity of vehicle over time
  anglevehicle(t) = zeros(3)                  # Angle of the vehicle over time

  angle = ()                                  # Angle of each tilting system (none)
  RPM = ()                                    # RPM of each rotor system (none)

  maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)

  
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
  path = joinpath(@__DIR__, "n_sweep", run_name, run_name * "_convergence.csv")                    
  result_df = CSV.read(path, DataFrame)
  cl_mean = mean(result_df[end-2:end, :CL])
  cd_mean = mean(result_df[end-2:end, :CD])

  if prev_cl !== nothing && prev_cd !== nothing
    percent_change_cl = abs(cl_mean - prev_cl) / abs(prev_cl)
    percent_change_cd = abs(cd_mean - prev_cd) / abs(prev_cd)
    println("Percent Change in CL = $(round(percent_change_cl*100, digits=3))%")
    println("Percent Change in CD = $(round(percent_change_cd*100, digits=3))%")

    if percent_change_cl < tolerance && percent_change_cd < tolerance
      println("Converged at n = $n (Percent Change in CL = $(percent_change_cl*100)%, Percent Change in CD = $(percent_change_cd*100)%)")
      println("Final averaged values: CL = $cl_mean, CD = $cd_mean at n = $n")
      global cl_mean_converged = cl_mean
      global cd_mean_converged = cd_mean
      global n_converged = n

      global p_per_step_min = Int(ceil( (2 * n_converged * lambda_vpm * magVinf * ttot) / (nsteps_converged * b)))
      println("Minimum required p_per_step = $p_per_step_min to satisfy sigma < b/(2n)")

      break
    end
  end
  prev_cl, prev_cd = cl_mean, cd_mean
end


if p_per_step < p_per_step_min
  println("Need to run another test to make sure particles don't overlap")
end



println("\n========= CONVERGENCE RESULTS =========")
println("Wake length converged: $(round(wakelength_converged, digits=3)) m")
println("nsteps converged: $nsteps_converged")
println("n converged: $n_converged")
println("p_per_step_min: $p_per_step_min")
println("=======================================\n")





println("\n========== VERIFICATION =============")
println("\n====== TEST 1: DOUBLING NSTEPS ======")
nsteps_verify = nsteps_converged * 2
run_name = "verification_nsteps$(nsteps_verify)"
save_path = joinpath(@__DIR__, "verification", run_name)


println("Running verification with nsteps = $nsteps_verify and adjusting for minimum p_per_step")


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
  n               = n_converged               # Number of spanwise elements per side
  r               = 10.0                      # Geometric expansion of elements
  central         = false                     # Whether expansion is central

  # ----------------- SOLVER PARAMETERS ------------------------------------------
  # Time parameters
  wakelength      = wakelength_converged                    # (m) length of wake to be resolved
  ttot            = wakelength/magVinf                      # (s) total simulation time
  nsteps          = nsteps_verify                           # Number of time steps

  # VLM and VPM parameters
  lambda_vpm      = 2.0                                     # VPM core overlap
  p_per_step = Int(ceil( (2 * n_converged * lambda_vpm * magVinf * ttot) / (nsteps_verify * b)))

  
  sigma_vpm_overwrite = lambda_vpm * magVinf * (ttot/nsteps)/p_per_step # Smoothing core size
  sigma_vlm_solver= -1                                      # VLM-on-VLM smoothing radius (deactivated with <0)
  sigma_vlm_surf  = 0.05*b                                  # VLM-on-VPM smoothing radius

  shed_starting   = true                                    # Whether to shed starting vortex
  vlm_rlx         = 0.7                                     # VLM relaxation




  # ----------------- 1) VEHICLE DEFINITION --------------------------------------
  println("Generating geometry...")

  # Generate wing
  wing = vlm.simpleWing(b, ar, tr, twist_root, lambda, gamma;
                          twist_tip=twist_tip, n=n, r=r, central=central);


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




  # ------------- 2) MANEUVER DEFINITION -----------------------------------------

  Vvehicle(t) = zeros(3)                      # Translational velocity of vehicle over time
  anglevehicle(t) = zeros(3)                  # Angle of the vehicle over time

  angle = ()                                  # Angle of each tilting system (none)
  RPM = ()                                    # RPM of each rotor system (none)

  maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)

  
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

path = joinpath(@__DIR__, "verification", run_name, run_name * "_convergence.csv")
verify_df = CSV.read(path, DataFrame)

cl_verify = mean(verify_df[end-2:end, :CL])
cd_verify = mean(verify_df[end-2:end, :CD])

println("Verification results:")
println("CL = $cl_verify , CD = $cd_verify")
println("Original converged CL = $cl_mean_converged, CD = $cd_mean_converged")

percent_change_cl = abs(cl_verify - cl_mean_converged) / abs(cl_mean_converged) * 100
percent_change_cd = abs(cd_verify - cd_mean_converged) / abs(cd_mean_converged) * 100

println("Percent Change CL = $(round(percent_change_cl, digits=4))%, Percent Change CD = $(round(percent_change_cd, digits=4))%")

if percent_change_cl < 0.5 && percent_change_cd < 0.5
    println("Passed: Doubling nsteps does not change CL or CD significantly")
else
    println("Failed: CL/CD changed more than expected")
end
println("=======================================\n")


cl_plot = plot(verify_df[:, :T], verify_df[:, :CL], linewidth = 1.5, xlabel = "Time (s)", ylabel = "Coefficient of Lift", grid = false)
cd_plot = plot(verify_df[:, :T], verify_df[:, :CD], linewidth = 1.5, xlabel = "Time (s)", ylabel = "Coefficient of Drag", grid = false)

mkpath(joinpath(@__DIR__, "Figures"))
save_path = joinpath(@__DIR__, "Figures", "CL_vs_time.png")
savefig(cl_plot, save_path)
save_path = joinpath(@__DIR__, "Figures", "CD_vs_time.png")
savefig(cd_plot, save_path)