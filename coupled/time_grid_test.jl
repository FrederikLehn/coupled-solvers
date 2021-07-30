using JFVM
using PyPlot
include("CRT.jl")

#Calculating error as a function of timesteps

# define the domain and constants
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

# Defining time-stepping parameters
tol = 1e-10 # convergence tolerance

# boundary conditions
C_bc = [1.0, 0.5, 0.0] # boundary concentration
C_0 = zeros(N_comp) # initial concentration

NX = 100:100:1000
time_array = zeros(5, length(NX))

for i=1:length(NX)
  Lx = 1.0 # [m]
  Nx = NX[i]  # number of grids
  m = createMesh1D(Nx, Lx)  # mesh

  #Time-stepping parameters
  dt = Lx/(u/ϕ)/Nx
  T = Nx/2*dt

  tic()
  C_SNIA = CRT.CRT_SNIA(m, Nx, ϕ, u, D, K_eq, nu_A, nu_B, nu_C, C_bc, C_0, N_comp, dt, T)
  time_array[1, i] = toc()

  tic()
  C_SIA = CRT.CRT_SIA(m, Nx, ϕ, u, D, K_eq, nu_A, nu_B, nu_C, C_bc, C_0, N_comp, dt, T, tol)
  time_array[2, i] = toc()

  tic()
  C_eps = CRT.CRT_GIA_eps(m, Nx, ϕ, u, D, K_eq, nu_A, nu_B, nu_C, C_bc, C_0, N_comp, dt, T, tol)
  time_array[3, i] = toc()

  tic()
  C_GIA_kin_eq = CRT.CRT_GIA_kin_eq(m, Nx, ϕ, u, D, K_eq, nu_A, nu_B, nu_C, a, b, c, C_bc, C_0, N_comp, dt, T, tol)
  time_array[4, i] = toc()

  tic()
  C_GIA_kin_fb = CRT.CRT_GIA_kin_fb(m, Nx, ϕ, u, D, k_f, k_b, nu_A, nu_B, nu_C, a, b, c, C_bc, C_0, N_comp, dt, T, tol)
  time_array[5, i] = toc()

  println(i)
end

figure()
plot(NX, time_array[1,:])
plot(NX, time_array[2,:])
plot(NX, time_array[3,:])
plot(NX, time_array[4,:])
plot(NX, time_array[5,:])
grid("on")
legend(["SNIA", "SIA", "GIA ("L"$\epsilon$"")", "GIA ("L"$K_{eq}$"")", "GIA ("L"$k_fk_b$"")"])
xlabel(L"$N_x$")
ylabel("Time [s]""")
savefig("time_grid_plot.png")

#CRT_SIA title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}\epsilon=0,\quad K_{eq}=$""$K_eq"" (SIA)"

#CRT_SNIA title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}\epsilon=0,\quad K_{eq}=$""$K_eq"" (SNIA)"

#CRT_GIA_eps title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}\epsilon=0,\quad K_{eq}=$""$K_eq"" (FIM)"

#CRT_GIA_kin_eq title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}(K_{eq}C_A^aC_B^b-C_C^c)=0,\quad K_{eq}=$""$K_eq"" (FIM)"

#CRT_GIA_kin_fb title:
#L"$\phi\frac{\partial C_{\alpha}}{\partial t}+\nabla\cdot(uC_{\alpha})-\nabla\cdot(D_{\alpha}\nabla C_{\alpha})+\nu_{\alpha}(k_{f}C_A^aC_B^b-k_{b}C_C^c)=0,\quad K_{eq}=$""$K_eq"" (FIM)"
