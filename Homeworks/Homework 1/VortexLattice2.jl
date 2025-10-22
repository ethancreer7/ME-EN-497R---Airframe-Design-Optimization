using VortexLattice

xle = [0.0, 0.4] # leading edge x-position
yle = [0.0, 7.5] # leading edge y-position
zle = [0.0, 0.0] # leading edge z-position
chord = [2.2, 1.8] # chord length
theta = [2.0*pi/180, 2.0*pi/180] # twist (in radians)
phi = [0.0, 0.0] # section rotation about the x-axis
fc = fill((xc) -> 0, 2) # camberline function for each section (y/c = f(x/c))

ns = 12 # number of spanwise panels
nc = 6  # number of chordwise panels
spacing_s = Sine() # spanwise discretization scheme
spacing_c = Uniform() # chordwise discretization scheme

grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)

grids = [grid]
ratios = [ratio]
system = System(grids; ratios)

Sref = 30.0 # reference area
cref = 2.0  # reference chord
bref = 15.0 # reference span
rref = [0.50, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity (magnitude)
ref = Reference(Sref, cref, bref, rref, Vinf)

alpha = 1.0*pi/180 # angle of attack
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
fs = Freestream(Vinf, alpha, beta, Omega)

symmetric = false

steady_analysis!(system, ref, fs; symmetric)

CF, CM = body_forces(system; frame=Wind())

# extract aerodynamic forces
CD, CY, CL = CF
Cl, Cm, Cn = CM

CDiff = far_field_drag(system)

# calculate lifting line coefficients
cf, cm = lifting_line_coefficients(system; frame=Body())

dCFb, dCMb = body_derivatives(system)

# traditional names for each body derivative
CXu, CYu, CZu = dCFb.u
CXv, CYv, CZv = dCFb.v
CXw, CYw, CZw = dCFb.w
CXp, CYp, CZp = dCFb.p
CXq, CYq, CZq = dCFb.q
CXr, CYr, CZr = dCFb.r
Clu, Cmu, Cnu = dCMb.u
Clv, Cmv, Cnv = dCMb.v
Clw, Cmw, Cnw = dCMb.w
Clp, Cmp, Cnp = dCMb.p
Clq, Cmq, Cnq = dCMb.q
Clr, Cmr, Cnr = dCMb.r

dCFs, dCMs = stability_derivatives(system)

# traditional names for each stability derivative
CDa, CYa, CLa = dCFs.alpha
Cla, Cma, Cna = dCMs.alpha
CDb, CYb, CLb = dCFs.beta
Clb, Cmb, Cnb = dCMs.beta
CDp, CYp, CLp = dCFs.p
Clp, Cmp, Cnp = dCMs.p
CDq, CYq, CLq = dCFs.q
Clq, Cmq, Cnq = dCMs.q
CDr, CYr, CLr = dCFs.r
Clr, Cmr, Cnr = dCMs.r

write_vtk("simplewing", system)