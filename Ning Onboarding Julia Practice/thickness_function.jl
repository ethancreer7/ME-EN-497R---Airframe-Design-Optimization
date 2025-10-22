"""
thickness(x,m)
This function calculates the thickness of a NACA 4-series airfoil given an x-position along the chord and an maximum thickness percentage of the chord (m)
"""
function thickness(x,m)
    return 10*m*(0.2969*sqrt(x) - 0.1260*x - 0.3537*x^2 + 0.2843*x^3 - 0.1015*x^4)
end
# println("\e[4m  x | t  \e[0m")
# for i in range(0.0,1.0,step=0.1)
#     t = thickness(i, 0.1)
#     println("$i | $t")
# end