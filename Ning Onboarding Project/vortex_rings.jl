using LinearAlgebra
using Plots

#Initial Conditions
d = 1.0
dt = 0.05
time = 0

tao1 = [0, 0, 1] # Used for Vortex 2 and 3
tao2 = [0, 0, -1] # Used for Vortex 3 and 4

P1 = [0.0, -0.5, 0.0] # Working position for Vortex 1
P2 = [0.0, 0.5, 0.0]  # Working position for Vortex 2
P3 = [1.0, 0.5, 0.0] # Working position for Vortex 3
P4 = [1.0, -0.5, 0.0] # Working Position for Vortex 4

vortex1 = [0.0 -0.5 0.0] # Cumulative array for Vortex 1 Position
vortex2 = [0.0 0.5 0.0] # Cumulative array for Vortex 2 Position
vortex3 = [1.0 0.5 0.0] # Cumulative array for Vortex 3 Position
vortex4 = [1.0 -0.5 0.0] # Cumulative array for Vortex 4 Position

"""
Returns the distance between two vortexes
"""
function RadiusCalc(p1,p2)

    radius = sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2)
    return radius
end

"""
Returns the induced velocity from vortex A on Vortex B
"""
function VelocityInfluence(p1,p2)
    rad = p2 .- p1
    if p1[2] < 0
        cross_vector = cross(tao2, rad)
    elseif p1[2] > 0
        cross_vector = cross(tao1, rad)
    else
        println("Sorry, your tao logic is incorrect.")
    end
    r = RadiusCalc(p1, p2)
    return cross_vector / (2 * pi * r * r)
end

"""
Returns the new position of Vortex 1, following 
the induced velocities from the other vortexes
"""
function Vortex1()
    inf2 = VelocityInfluence(P2,P1)
    inf3 = VelocityInfluence(P3, P1) 
    inf4 = VelocityInfluence(P4, P1) 
    final = inf2 .+ inf3 .+ inf4
    return final .* dt .+ P1
end

"""
Returns the new position of Vortex 2, following 
the induced velocities from the other vortexes
"""
function Vortex2()
   
    inf1 = VelocityInfluence(P1,P2)
    inf3 = VelocityInfluence(P3, P2)
    inf4 = VelocityInfluence(P4, P2)
    final = inf1 .+ inf3 .+ inf4
    return final .* dt .+ P2
end

"""
Returns the new position of Vortex 2, following 
the induced velocities from the other vortexes
"""
function Vortex3()
    inf1 = VelocityInfluence(P1,P3)
    inf2 = VelocityInfluence(P2, P3)
    inf4 = VelocityInfluence(P4, P3)
    final = inf1 .+ inf2 .+ inf4
    return final .* dt .+ P3
end

"""
Returns the new position of Vortex 2, following 
the induced velocities from the other vortexes
"""
function Vortex4()
    inf1 = VelocityInfluence(P1,P4)
    inf2 = VelocityInfluence(P2, P4)
    inf3 = VelocityInfluence(P3, P4)
    final = inf1 .+ inf2 .+ inf3
    return final .* dt .+ P4
end


"""
This is the main function called that calls the 
above functions to calculate the new position of 
each vortex
"""
function CalcFinalPosition()

    final1 = Vortex1()
    final2 = Vortex2()
    final3 = Vortex3()
    final4 = Vortex4()
    return final1, final2, final3, final4
end


"""
Plots the vortex paths with adjustable colors, linewidth, and
linestyle; Saves image as 'vortex_paths.png'
"""
function PlotVortexes(vortex1,vortex2,vortex3,vortex4)
    plot(vortex1[:,1],vortex1[:,2,],color=:orange,linewidth=1.5,linestyle=:solid,label="Toroidal Vortex 1")
    plot!(vortex2[:,1],vortex2[:,2],color=:orange,linewidth=1.5,linestyle=:solid,label=false)
    plot!(vortex3[:,1],vortex3[:,2],color=:purple,linewidth=1.5,linestyle=:solid,label="Toroidal Vortex 2")
    plot!(vortex4[:,1],vortex4[:,2],color=:purple,linewidth=1.5,linestyle=:solid,label=false)
    plot!(legend=:topright,grid=false, xticks=[0,2.5,5.0,7.5,10.0],yticks=[0])
    savefig("vortex_paths.png")
end


while time < 40 # Program runs until 40 seconds is reached
    P1, P2, P3, P4 = CalcFinalPosition()
    vortex1 = vcat(vortex1,P1') # Concatenate new position onto existing position record array
    vortex2 = vcat(vortex2,P2')
    vortex3 = vcat(vortex3,P3')
    vortex4 = vcat(vortex4,P4')
    time += 0.05 # Move time forward
end

PlotVortexes(vortex1,vortex2,vortex3,vortex4)

v1x = vortex1[:,1] #Separate x and y data for the animation below
v2x = vortex2[:,1]
v3x = vortex3[:,1]
v4x = vortex4[:,1]
v1y = vortex1[:,2]
v2y = vortex2[:,2]
v3y = vortex3[:,2]
v4y = vortex4[:,2]

n = 800
anim = @animate for i in 1:n
    plot(v1x[1:i], v1y[1:i], xlims=(0,12.5),xticks=[0.0,2.5,5.0,7.5,10.0],ylims=(-0.8,0.8),yticks=[0],color=:orange,legend=:topright,grid=false,label="Toroidal Vortex 1")
    plot!(v2x[1:i],v2y[1:i],color=:orange,label=false)
    plot!(v3x[1:i],v3y[1:i],color=:purple,label="Toroidal Vortex 2")
    plot!(v4x[1:i],v4y[1:i],color=:purple,label=false)
    plot!([v2x[i],v1x[i]],[v2y[i],v1y[i]],color=:orange,linestyle=:dash,alpha=0.65,label=false)
    plot!([v3x[i],v4x[i]],[v3y[i],v4y[i]],color=:purple,linestyle=:dash,alpha=0.65,label=false)
end

gif(anim, "vortex_anim.gif",fps=30)