using CSV
using DataFrames
using Statistics
using Printf

# --- CONFIGURATION ---
filename = joinpath(@__DIR__, "nrevs", "nrevs8_convergence.csv")
ct_col_name = "CT_1"              # Replace with the exact name of your CT column header
                                # (You can check column names with `names(df)`)

# --- 1. LOAD DATA ---
df = CSV.read(filename, DataFrame)

# --- 2. CLEAN DATA ---
# Get the angle column (Assuming it is the first column as stated)
angle_col = names(df)[1]

# Fix rounding artifacts:
# We round to the nearest integer to ensure 359.999 becomes 360.0
# This prevents boundary points from falling into the wrong revolution bucket.
df.clean_angle = round.(df[!, angle_col])

# --- 3. ASSIGN REVOLUTION NUMBERS ---
# Old line:
# df.rev_num = floor.(Int, df.clean_angle ./ 360) .+ 1

# New robust line:
# We subtract a tiny epsilon (like 1e-6) so 2880.0 becomes 2879.999... 
# which floors correctly into Rev 8.
df.rev_num = floor.(Int, (df.clean_angle .- 1e-6) ./ 360) .+ 1

# --- 4. CALCULATE STATISTICS ---
# Group by revolution number and calculate the mean CT for each
results = combine(groupby(df, :rev_num), ct_col_name => mean => :mean_CT)

# Sort just in case the data wasn't sequential (though it likely is)
sort!(results, :rev_num)

# --- 5. CALCULATE PERCENT CHANGE ---
# Initialize the column with 0.0 or NaN
results.pct_change = zeros(Float64, nrow(results))

# Loop starting from the 2nd revolution to compare against the previous
for i in 2:nrow(results)
    current_val = results.mean_CT[i]
    prev_val    = results.mean_CT[i-1]
    
    # Formula: ((New - Old) / Old) * 100
    results.pct_change[i] = ((current_val - prev_val) / prev_val) * 100
end

# --- 6. DISPLAY OUTPUT ---
println("\nConvergence Analysis (Revolution over Revolution):")
println("----------------------------------------------------")
@printf("%-10s | %-15s | %-15s\n", "Rev #", "Mean CT", "% Change")
println("----------------------------------------------------")

for row in eachrow(results)
    # Handle the first row (no % change) separately for cleaner printing
    pct_str = row.rev_num == 1 ? "N/A" : @sprintf("%+.4f%%", row.pct_change)
    @printf("%-10d | %-15.6f | %s\n", row.rev_num, row.mean_CT, pct_str)
end

println("----------------------------------------------------")