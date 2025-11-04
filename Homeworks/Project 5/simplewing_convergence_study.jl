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

# Define Global Variables for referencing converged values
global wakelength_converged = nothing
global nsteps_converged = nothing
global n_converged = nothing
global p_per_step_min = nothing
global p_per_step = nothing

# Make subfolders for outputs
mkpath(joinpath(@__DIR__, "wakelength"))
mkpath(joinpath(@__DIR__, "nsteps"))
mkpath(joinpath(@__DIR__, "n_sweep"))
mkpath(joinpath(@__DIR__, "verification"))
mkpath(joinpath(@__DIR__, "Figures"))

# ----------------- GLOBAL CONSTANTS AND PARAMETERS-------------------------------------------
# --- Simulation parameters ---
AOA = 4.2                     # (deg) Angle of attack
magVinf = 49.7                # (m/s) Freestream velocity
rho = 0.93                    # (kg/m^3) Air density
lambda_vpm = 2.0              # VPM core overlap

# --- Wing geometry ---
b = 2.489                     # (m) wing span
ar = 5.0                      # Aspect ratio
tr = 1.0                      # Taper ratio
twist_root = 0.0              # (deg) twist at root
twist_tip = 0.0               # (deg) twist at tip
lambda = 45.0                 # (deg) sweep angle
gamma = 0.0                   # (deg) dihedral angle


qinf = 0.5*rho*magVinf^2      # (Pa) Freestream dynamic pressure
r = 10.0
central = false
sigma_vlm_solver = -0.1
sigma_vlm_surf = 0.05 * b
shed_starting = true
vlm_rlx = 0.7
paraview = true
lambda_vpm = 2.0



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

"""
  function(lots_of_stuff)

Does all of the heavy stuff that's called repetitively.
"""
function set_parameters_run_study(study, run_name, n, wakelength, nsteps; ttot = wakelength/magVinf,   p_per_step = Int(ceil((2 * n * lambda_vpm * magVinf * ttot) / (nsteps * b))))

  ttot = wakelength/magVinf
  sigma_vpm_overwrite = lambda_vpm * magVinf * (ttot/nsteps)/p_per_step

  Vinf(X, t) = magVinf*[cosd(AOA), 0.0, sind(AOA)]
  save_path = joinpath(@__DIR__, study, run_name)

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

end

"""
  function verification_stats(verify_df::DataFrame)

Calculates and prints verification confirmation stats
"""
function verification_stats(verify_df::DataFrame, variable::String)
  cl_verify = mean(verify_df[end-2:end, :CL])
  cd_verify = mean(verify_df[end-2:end, :CD])

  println("Verification results:")
  println("CL = $cl_verify , CD = $cd_verify")
  println("Original converged CL = $cl_mean_converged, CD = $cd_mean_converged")

  percent_change_cl = abs(cl_verify - cl_mean_converged) / abs(cl_mean_converged) * 100
  percent_change_cd = abs(cd_verify - cd_mean_converged) / abs(cd_mean_converged) * 100

  println("Percent Change CL = $(round(percent_change_cl, digits=4))%, Percent Change CD = $(round(percent_change_cd, digits=4))%")

  if percent_change_cl < 0.5 && percent_change_cd < 0.5
      println("Passed: Doubling $variable does not change CL or CD significantly")
  else
      println("Failed: CL/CD changed more than expected")
  end
  println("=======================================\n")
end


# --------------- Wakelength Convergence Study ----------------
for (count, i) in enumerate(range(0.75, 5.0, step = 0.25))
  run_name        = "wakelength$(round(i*b, digits=3))"             # Name of this simulation
  set_parameters_run_study("wakelength", run_name, 30, i*2.489, 20)
    
  path = joinpath(@__DIR__, "wakelength", run_name, run_name * "_convergence.csv")                    
  result_df = CSV.read(path, DataFrame)                    
  cl_CV, cd_CV = get_cv(result_df)

  if cl_CV < 0.002 && cd_CV < 0.002
    printstats(cl_CV, cd_CV, "Wakelength", i * b)
    global wakelength_converged = i * b
    break
  end  
end


# --------------- nsteps Convergence Study --------------------
for (count, i) in enumerate(range(10, 100, step = 10))
  run_name        = "nsteps$(i)"            # Name of this simulation
  set_parameters_run_study("nsteps", run_name, 70, wakelength_converged, i)

  path = joinpath(@__DIR__, "nsteps", run_name, run_name * "_convergence.csv")                    
  result_df = CSV.read(path, DataFrame)
  cl_CV, cd_CV = get_cv(result_df)

  if cl_CV < 0.002 && cd_CV < 0.002
    printstats(cl_CV, cd_CV, "nsteps", i)
    global nsteps_converged = i
    break
  end

end


# --------------- n Convergence Study -------------------------
for (count, i) in enumerate(range(10, 100, step = 10))
  run_name = "n_sweep$(i)"

  set_parameters_run_study("n_sweep", run_name, i, wakelength_converged, nsteps_converged)

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
      println("Converged at n = $i (Percent Change in CL = $(percent_change_cl*100)%, Percent Change in CD = $(percent_change_cd*100)%)")
      println("Final averaged values: CL = $cl_mean, CD = $cd_mean at n = $i")
      global cl_mean_converged = cl_mean
      global cd_mean_converged = cd_mean
      global n_converged = i
      ttot = wakelength_converged / magVinf
      global p_per_step_min = Int(ceil( (2 * n_converged * lambda_vpm * magVinf * ttot) / (nsteps_converged * b)))
      println("Minimum required p_per_step = $p_per_step_min to satisfy sigma < b/(2n)")
      break
    end
  end
  prev_cl, prev_cd = cl_mean, cd_mean
end


println("\n========= CONVERGENCE RESULTS: SUGGESTED PARAMETERS =========")
println("Wake length converged: $(round(wakelength_converged, digits=3)) m")
println("nsteps converged: $nsteps_converged")
println("n converged: $n_converged")
println("p_per_step_min: $p_per_step_min")
println("=======================================\n")



println("\n========== VERIFICATION OF SUGGESTED PARAMETERS =============")
println("\n====== TEST 1: DOUBLING NSTEPS ======")
nsteps_verify = nsteps_converged * 2
run_name = "verification_nsteps$(nsteps_verify)"
println("Running verification with increased nsteps = $nsteps_verify and adjusting for minimum p_per_step")

set_parameters_run_study("verification", run_name, n_converged, wakelength_converged, nsteps_verify)

path = joinpath(@__DIR__, "verification", run_name, run_name * "_convergence.csv")
verify_df = CSV.read(path, DataFrame)

verification_stats(verify_df, "nsteps")


cl_plot_v1 = plot(verify_df[:, :T], verify_df[:, :CL], linewidth = 1.5, xlabel = "Time (s)", ylabel = "Coefficient of Lift", grid = false)
cd_plot_v1 = plot(verify_df[:, :T], verify_df[:, :CD], linewidth = 1.5, xlabel = "Time (s)", ylabel = "Coefficient of Drag", grid = false)
save_path = joinpath(@__DIR__, "Figures", "CL_vs_time.png")
savefig(cl_plot_v1, save_path)
save_path = joinpath(@__DIR__, "Figures", "CD_vs_time.png")
savefig(cd_plot_v1, save_path)



println("\n====== TEST 2: INCREASING P_PER_STEP_MIN ======")
p_per_step_verify = Int(ceil(p_per_step_min * 1.5))
run_name = "verification_p_per_step$(p_per_step_verify)"
println("Running verification with increased p_per_step = $p_per_step_verify")
set_parameters_run_study("verification", run_name, n_converged, wakelength_converged, nsteps_converged; p_per_step = p_per_step_verify)
path = joinpath(@__DIR__, "verification", run_name, run_name * "_convergence.csv")
verify_df = CSV.read(path, DataFrame)
verification_stats(verify_df, "p_per_step")



println("\n======= Test 3: INCREASING WAKELENGTH ========")
wakelength_verify = wakelength_converged * 1.5
run_name = "verification_wakelength$(wakelength_verify)"
println("Running verification with increased wakelength = $wakelength_verify")
set_parameters_run_study("verification", run_name, n_converged, wakelength_verify, nsteps_converged)
path = joinpath(@__DIR__, "verification", run_name, run_name * "_convergence.csv")
verify_df = CSV.read(path, DataFrame)
verification_stats(verify_df, "wakelength")