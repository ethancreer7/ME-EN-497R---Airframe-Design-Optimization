using Plots
using LaTeXStrings
using Printf

# ==========================================
# 1. DATA INPUT
# ==========================================
nsteps = [20, 40, 80, 160, 320]
CT_simulated = [0.08116769, 0.08273388, 0.08416609, 0.08549414, 0.08629362]
CT_extrapolated = 0.0871768767
inv_nsteps = 1.0 ./ nsteps

# ==========================================
# 2. PLOTTING
# ==========================================
gr() 

p = plot(
    inv_nsteps, 
    CT_simulated,
    seriestype = :scatter,
    # Removed \phantom. If you still need vertical space, we can tweak font size.
    label = "Simulated \$C_T\$\n", 
    markersize = 5,           # Slightly larger dots
    markerstrokecolor = :black,
    markerstrokewidth = 0.5,
    grid = false,
    size = (600, 450),        # Slight height increase helps spacing
    dpi = 300,
    fontfamily = "Computer Modern",
    
    # --- LEGEND SETTINGS ---
    # (x, y) coordinates are relative to the graph area (0.0 to 1.0).
    # (0.55, 0.75) puts it in the top-right area, but lowered.
    legend = (0.57, 0.85),
    legend_font_pointsize = 8, # Slightly smaller font helps prevent overlap
    background_color_legend = :white
)

# Add the Extrapolated Line
hline!(
    [CT_extrapolated],
    label = L"Richardson Extrapolated ($N_{steps} \to \infty$)",
    linestyle = :dash,
    color = :red,
    linewidth = 2
)

# ==========================================
# 3. FORMATTING LABELS
# ==========================================
xlabel!(L"Inverse Time Steps ($1 / N_{steps}$)")
ylabel!(L"Thrust Coefficient ($C_T$)")

# Margins to prevent label cutoff
plot!(left_margin=6Plots.mm, bottom_margin=6Plots.mm)

savepath = joinpath(@__DIR__, "Figures", "FINAL320_nsteps_richardson.png")
savefig(p, savepath)

display(p)