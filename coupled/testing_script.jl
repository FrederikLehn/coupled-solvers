using JFVM
using PyPlot
include("CRT.jl")


# define the domain and constants
Lx = 1.0 # [m]
Nx = 100  # number of grids
ϕ  = 0.4 # porosity
u  = 1.0/(3600*24) # m/s, eq. 1 m/day
D  = [1e-9, 1e-9, 1e-9] # m^2/s diffusivity
K_eq = 10.0
k_f = 100.0
k_b = k_f/K_eq
nu_A = 1.0 # stoichiometric coefficient of component A
nu_B = 1.0 # stoichiometric coefficient of component B
nu_C = 1.0 # stoichiometric coefficient of component C
a = 1.0 #order of reaction
b = 1.0 #order of reaction
c = 1.0 #order of reaction
N_comp = 3 # number of components

# create the mesh
m = createMesh1D(Nx, Lx)  # mesh
# m = createMeshRadial2D(Nx, Nx, Lx, Lx)

# boundary conditions
C_bc = [1.0, 0.5, 0.0] # boundary concentration
C_0 = zeros(N_comp) # initial concentration

dt = Lx/(u/ϕ)/Nx
T = Nx/2*dt
tol = 1e-10 # convergence tolerance

#Calculating reference "true solution"
C_old = CRT.CRT_GIA_eps(m, Nx, ϕ, u, D, K_eq, nu_A, nu_B, nu_C, C_bc, C_0, N_comp, dt, T, tol)
#C_old = CRT.CRT_SNIA(m, Nx, ϕ, u, D, K_eq, nu_A, nu_B, nu_C, C_bc, C_0, N_comp, dt, T)
# Visualizing
fig, ax = subplots()
x = [C_old[1].domain.facecenters.x[1]; C_old[1].domain.cellcenters.x; C_old[1].domain.facecenters.x[end]]
  phi = [0.5*(C_old[1].value[1]+C_old[1].value[2]); C_old[1].value[2:end-1]; 0.5*(C_old[1].value[end-1]+C_old[1].value[end])]
  plot(x, phi, label=L"$C_{A},\ \nu_A=$""$nu_A"L",\ a=""$a")

  x = [C_old[2].domain.facecenters.x[1]; C_old[2].domain.cellcenters.x; C_old[2].domain.facecenters.x[end]]
    phi = [0.5*(C_old[2].value[1]+C_old[2].value[2]); C_old[2].value[2:end-1]; 0.5*(C_old[2].value[end-1]+C_old[2].value[end])]
    plot(x, phi, label=L"$C_{B},\ \nu_B=$""$nu_B"L",\ b=""$b")

    x = [C_old[3].domain.facecenters.x[1]; C_old[3].domain.cellcenters.x; C_old[3].domain.facecenters.x[end]]
      phi = [0.5*(C_old[3].value[1]+C_old[3].value[2]); C_old[3].value[2:end-1]; 0.5*(C_old[3].value[end-1]+C_old[3].value[end])]
      plot(x, phi, label=L"$C_{C},\ \nu_C=$""$nu_C"L",\ c=""$c")
ax[:legend]()
grid("on")
xlabel(L"$L_x\ [m]$")
ylabel(L"$C_{\alpha}\ [mol/m^3]$")
#ax[:set_title](L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(\mathcal{D}_{\alpha}\phi\nabla C_{\alpha})+\nu_{\alpha}\epsilon=0,\quad K_{eq}=$""$K_eq"" (GIA)")
ax[:set_title](L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}\epsilon=0,\quad K_{eq}=$""$K_eq"" (GIA)")
ax[:legend]()
savefig("C:/Users/Frederik/Dropbox/Backup/DTU/9. semester/Coupled Reactive Transport Modelling/Code/concentration_GIA_eps.png")
#savefig("C:/Users/Frederik/Dropbox/Backup/DTU/9. semester/Coupled Reactive Transport Modelling/Code/concentration_SNIA.png")


#CRT_SIA title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}\epsilon=0,\quad K_{eq}=$""$K_eq"" (SIA)"

#CRT_SNIA title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}\epsilon=0,\quad K_{eq}=$""$K_eq"" (SNIA)"

#CRT_GIA_eps title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}\epsilon=0,\quad K_{eq}=$""$K_eq"" (GIA)"

#CRT_GIA_kin_eq title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}(K_{eq}C_A^aC_B^b-C_C^c)=0,\quad K_{eq}=$""$K_eq"" (GIA)"

#CRT_GIA_kin_fb title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}(k_{f}C_A^aC_B^b-k_{b}C_C^c)=0,\quad K_{eq}=$""$K_eq"" (GIA)"
