using Printf
include("thickness_function.jl")
include("control_flow.jl")

function UpperAndLower(c, p, t, x)
    c = c/100
    t = t/100
    p = p/10
    return camber.(x, p, c) .+ thickness.(x, t)./2, camber.(x, p, c) .- thickness.(x, t)./2
end


    x = collect(range(0.0, 1.0, step=0.001))
    UpperVals, LowerVals = UpperAndLower(2,4,12,x)
    data = hcat(x, UpperVals, LowerVals)
    println("\n\e[4m  x   |      Upper     |       Lower      \e[0m")
    for i in eachindex(x)
        @printf("%5.2f | %14.8f | %14.8f\n", x[i], UpperVals[i], LowerVals[i])
    end