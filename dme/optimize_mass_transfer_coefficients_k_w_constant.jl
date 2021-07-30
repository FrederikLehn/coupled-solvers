using Roots, NLopt, PyPlot, LsqFit, JLD, Polynomials, JFVM, Dierckx
include("DME_solver.jl")

function rec_fact(param, G, K, MT, T, BC, IC, k_w, t_exp, R_exp)
  k_o = param[1]
  (t, R) = DME_solver(G, K, MT, T, BC, IC, k_o, k_w, false)
  R_int = Spline1D(t, R, k=1)
  return R_int(t_exp)
end

model = (t_exp, param) -> rec_fact(param, G, K, MT, T, BC, IC, k_w, t_exp, R_exp)

function error_calc(param, G, K, MT, T, BC, IC, k_w, t_exp, R_exp)
  return rec_fact(param, G, K, MT, T, BC, IC, k_w, t_exp, R_exp) - R_exp
end

function error_calculation(param, param_ind, G, K, MT, T, BC, IC, k_w, t_exp, R_exp; grad_calc=true)
  eps1 = 1e-8
  error_vals = error_calc(param, G, K, MT, T, BC, IC, k_w, t_exp, R_exp)
  if grad_calc
    jac  = zeros(length(error_vals), length(param_ind))
  else
    jac = zeros(0)
  end
  if grad_calc
    for j in eachindex(param_ind)
      param2 = copy(param)
      param2[param_ind[j]]+=eps1
      error_val2 = error_calc(param2, G, K, MT, T, BC, IC, k_w, t_exp, R_exp)
      jac[:, j] = (error_val2 - error_vals) / eps1
    end
  end
  return error_vals, jac
end

function objective_function(param::Vector{Float64}, grad::Vector{Float64}, param_all, param_ind, G, K, MT, T, BC, IC, k_w, t_exp, R_exp; w=1.0)
  param2 = copy(param_all)
  param2[param_ind] = param
  grad_calc = false
  if length(grad)>0
    grad_calc = true
  end
  error_val, jac = error_calculation(param2, param_ind, G, K, MT, T, BC, IC, k_w, t_exp, R_exp, grad_calc = grad_calc)
  obj_func_value = sum(abs2, w.*error_val)
  if length(grad)>0
      grad[:] = 2.0*sum(w.*error_val.*jac, 1)
  end
  return obj_func_value
end

obj_fun = (param, grad) -> objective_function(param::Vector{Float64}, grad::Vector{Float64}, param_all, param_ind, G, K, MT, T, BC, IC, k_w, t_sec_cor, R_oil, w=w)



  perm_val  = 0.01e-12 # [m^2] permeability
  poros_val = 0.40     # [-] porosity
  ρ_rock    = 2700     # [kg/m^3]

  u_inj = 1.0/(3600*24) # 1 m/day to m/s injection velocity
  c_inj1 = 0.0 # mol/m^3 DME in injection water at initial water flooding
  c_inj2 = 10.0 # mol/m^3 DME in injection water at DME flooding
  c_tracer1 = 0.0
  c_tracer2 = 1.0

  Kval = 2.0
  a = 2.0*2700.0*1000.0*(1.0-poros_val) # m^2/m^3 = 2 m^2/g * 2700 kg/m^3 * 1000 g/kg * (1-phi)
  D_oil = 0.5*10^(-9.0) # m^2/s
  D_water = 2*10^(-9.0) # m^2/s

  ## Creating DME parameters
  ρ_oil = c_o -> 748.0015057043045-0.0016431963703508173.*c_o+2.9364412314092913.*10.0.^(-7.0).*c_o.^2.0-4.107225978803957.*10.0.^(-11.0).*c_o.^3.0
  dρ_oil = c_o -> -0.0016431963703508173+2.0.*2.9364412314092913.*10.0.^(-7.0).*c_o-3.0.*4.107225978803957.*10.0.^(-11.0).*c_o.^2.0
  ρ_water = 983.6319389796079
  μ_oil = c_o -> 0.0010131130815849523-2.528969514668829.*10.0.^(-5.0).*c_o+2.3873417406201346.*10.0.^(-11.0).*c_o.^2.0
  dμ_oil = c_o -> -2.528969514668829.*10.0.^(-5.0)+2*2.3873417406201346.*10.0.^(-11.0).*c_o
  μ_water = c_w -> 0.00040706344808758737+1.47269029810352587.*10.0.^(-7.0).*c_w
  dμ_water = c_w -> 1.47269029810352587.*10.0.^(-7.0).*ones(size(c_w))

  MT = FluidParametersDME(ρ_oil, dρ_oil, ρ_water, μ_oil, dμ_oil, μ_water, dμ_water, Kval, a, D_oil, D_water)

  # Creating realitive permeability
  kro0 = 0.8
  krw0 = 0.2
  n_o  = 2.0
  n_w  = 2.0
  swc  = 0.08

  sor_min  = 0.05
  sor_max  = 0.3
  c_min = 0.0 # min DME mol/m^3 in oil
  c_max = Kval*c_inj2 # max DME mol/m^3 fraction in oil
  m_sor = -(sor_max-sor_min)/(c_max-c_min)
  sor_c = c -> m_sor*c+sor_max

  K = RelativePermeabilityDME(kro0, krw0, n_o, n_w, sor_c, swc)

  # =============================================================================
  # Preparing geology/mesh
  # =============================================================================
  dim = 1

  srand(123456543)
  λ_dp    = 0.8
  cor_x   = 0.1
  cor_y   = 0.1
  cor_z   = 0.1
  Lx      = 80 # [m]
  Ly      = 20 # [m]
  Lz      = 10 # [m]
  Nx      = 50  # number of grids in the x direction
  Ny      = 20  # number of grids in the y direction
  Nz      = 10 # number of grids in the z direction

  if dim == 1
      m           = createMesh1D(Nx, Lx)
      perm_mat    = permfieldlogrnde(Nx, perm_val, λ_dp, cor_x)
  elseif dim == 2
      m           = createMesh2D(Nx, Ny, Lx, Ly)
      perm_mat    = permfieldlogrnde(Nx, Ny, perm_val, λ_dp, cor_x, cor_y)
  else
      m           = createMesh3D(Nx, Ny, Nz, Lx, Ly, Lz)
      perm_mat    = permfieldlogrnde(Nx, Ny, Nz, perm_val, λ_dp, cor_x, cor_y, cor_z)
  end

  k       = createCellVariable(m, perm_mat)
  k_face  = harmonicMean(k)     # permeability on the cell faces
  ϕ       = createCellVariable(m, poros_val)

  G = Geology(m, perm_val, poros_val, k, ϕ)

  # =============================================================================
  # Preparing time-stepping parameters
  # =============================================================================

  ## Simulation Time
  n_pv_water1       = 0.8 # number of injected pore volumes
  t_final_water1    = n_pv_water1*Lx/(u_inj/poros_val) # [s] final time for water flood
  dt_water1         = t_final_water1/n_pv_water1/Nx/5 # [s] time step

  n_pv_DME          = 1.0
  t_final_DME       = t_final_water1+n_pv_DME*Lx/(u_inj/poros_val) # [s] final time for DME flood
  dt_DME            = t_final_DME/(n_pv_water1 + n_pv_DME)/Nx/5 # [s] time step

  n_pv_water2       = 0.5
  t_final_water2    = t_final_DME+n_pv_water2*Lx/(u_inj/poros_val) # [s] final time for DME flood
  dt_water2         = t_final_water2/(n_pv_water1 + n_pv_DME + n_pv_water2)/Nx/5 # [s] time step

  T_water1 = SimulationTime(t_final_water1, dt_water1)
  T_DME    = SimulationTime(t_final_DME, dt_DME)
  T_water2 = SimulationTime(t_final_water2, dt_water2)

  T = [T_water1, T_DME, T_water2]

  ## Boundary Conditions for water flooding
  #T0 = 70 + 273.15 # [K]

  c_oil0      = 0.0 # initial DME concentration in the oil phase
  c_water0    = 0.0 # initial DME concentration in the water phase
  p0          = 2000/14.7*1e5   # [Pa]
  sw0         = swc+0.05 # [vol frac] initial water saturation
  c_t0        = 0.0 # initial tracer concentration

  BCc_oil_water     = createBC(m) # concentration (DME in oil) boundary condition
  BCc_water_water   = createBC(m) # concentration (DME in water) boundary condition
  BCp_water         = createBC(m) # pressure boundary condition
  BCs_water         = createBC(m) # saturation boundary condition
  BCt_water         = createBC(m) # tracer concentration boundary condition

  BCp_water.right.a[:] = 0.0
  BCp_water.right.b[:] = 1.0 #Dirichlet
  BCp_water.right.c[:] = p0

  BCp_water.left.a[:]  = k_face.xvalue[1,:]./μ_water(c_inj1) #Neumann
  BCp_water.left.b[:]  = 0.0
  BCp_water.left.c[:]  = -u_inj

  BCs_water.left.a[:]  = 0.0
  BCs_water.left.b[:]  = 1.0 #Dirichlet
  BCs_water.left.c[:]  = 1.0 #Value S_w=1

  BCc_water_water.left.a[:]  = 0.0 #
  BCc_water_water.left.b[:]  = 1.0 #Dirichlet
  BCc_water_water.left.c[:]  = c_inj1 #injected concentration of DME in water

  BCc_oil_water.left.a[:]  = 0.0 #
  BCc_oil_water.left.b[:]  = 1.0 #Dirichlet
  BCc_oil_water.left.c[:]  = Kval*c_inj1 #injected concentration of DME in water

  BCt_water.left.a[:]  = 0.0
  BCt_water.left.b[:]  = 1.0 #Dirichlet
  BCt_water.left.c[:]  = c_tracer1

  BC_water = BoundaryConditionDME(BCc_oil_water, BCc_water_water, BCp_water, BCs_water, BCt_water)

  c_oil_init      = createCellVariable(m, c_oil0, BCc_oil_water)
  c_water_init    = createCellVariable(m, c_water0, BCc_water_water)
  p_init          = createCellVariable(m, p0, BCp_water)
  sw_init         = createCellVariable(m, sw0, BCs_water)
  c_t_init        = createCellVariable(m, c_t0, BCt_water)

  IC = InitialConditionDME(c_oil_init, c_water_init, p_init, sw_init, c_t_init)

  # Creating Boundary conditions for DME flooding
  BCc_oil_DME     = createBC(m) # concentration (DME in oil) boundary condition
  BCc_water_DME   = createBC(m) # concentration (DME in water) boundary condition
  BCp_DME         = createBC(m) # pressure boundary condition
  BCs_DME         = createBC(m) # saturation boundary condition
  BCt_DME         = createBC(m) # tracer concentration boundary condition


  BCc_water_DME.left.a[:] = 0.0
  BCc_water_DME.left.b[:] = 1.0 #Dirichlet
  BCc_water_DME.left.c[:] = c_inj2 #Injected concentration of DME

  BCc_oil_DME.left.a[:] = 0.0
  BCc_oil_DME.left.b[:] = 1.0 #Dirichlet
  BCc_oil_DME.left.c[:] = Kval*c_inj2 #Created from reaction at boundary

  BCp_DME.right.a[:] = 0.0
  BCp_DME.right.b[:] = 1.0 #Dirichlet
  BCp_DME.right.c[:] = p0

  BCp_DME.left.a[:]  = k_face.xvalue[1,:]/μ_water(c_inj2) #Neumann
  BCp_DME.left.b[:]  = 0.0
  BCp_DME.left.c[:]  = -u_inj

  BCs_DME.left.a[:]  = 0.0
  BCs_DME.left.b[:]  = 1.0 #Dirichlet
  BCs_DME.left.c[:]  = 1.0 #Value S_w=1

  BCt_DME.left.a[:]       = 0.0
  BCt_DME.left.b[:]       = 1.0 #Dirichlet
  BCt_DME.left.c[:]       = c_tracer2

  BC_DME = BoundaryConditionDME(BCc_oil_DME, BCc_water_DME, BCp_DME, BCs_DME, BCt_DME)

  BC = [BC_water, BC_DME, BC_water]


## Experimental results
t_sec_cor = [0, 16.0/6, 16.0/6*2, 16.0/6*3, 16.0/6*4, 16.0/6*5,  16.0, 82.0, 100.0, 112.0].*24.0.*60.0.*60.0
R_exp = [0.0, 10.0/100, 20.0/100, 30.0/100, 40.0/100, 50.0/100, 58.0/100, 62.0/100, 80.0/100, 88.0/100]
R_oil = R_exp

# =============================================================================
# Calculating heuristic for k_o and k_w
# =============================================================================
d = sqrt(perm_val) #Characteristic length of capillary in porous medium
c_o = Kval*c_inj2
c_w = c_inj2
v_o = u_inj/poros_val
v_w = u_inj/poros_val

k_o_ini = 1.1*10.0^(-0.0)*0.023*D_oil/d*((ρ_oil(c_o)*v_o*d)/μ_oil(c_o))^(0.83)*(μ_oil(c_o)/(ρ_oil(c_o)*D_oil))^(1.0/3.0)
k_w = 10.0^(-0.0)*0.023*D_water/d*((ρ_water*v_w*d)/μ_water(c_w))^(0.83)*(μ_water(c_w)/(ρ_water*D_water))^(1.0/3.0)


param_all = [k_o_ini]
param_ind = [1]
w = ones(length(R_exp))

x_init = copy(param_all)
x_lb = [0.0]
x_ub = [1e-5]

# =============================================================================
# Optimizing
# =============================================================================

opt_type = 2

if opt_type == 1
    fit = curve_fit(model, t_sec_cor, R_exp, w, x_init, lower = x_lb, upper = x_ub)
    print("k_o_ini=", k_o_ini)
    print("k_o_opt=", fit.param[1])

elseif opt_type == 2
    # algorithms
    # :LD_MMA
    # :LN_COBYLA
    # :LD_LBFGS
    # :GN_DIRECT
    # :GN_DIRECT_L
    opt_alg=:LD_MMA

    opt1 = Opt(opt_alg, length(x_init)) # choose the algorithm
    lower_bounds!(opt1, x_lb)
    upper_bounds!(opt1, x_ub)
    ftol_rel!(opt1, 1e-5)
    ftol_abs!(opt1, 1e-5)

    min_objective!(opt1, obj_fun)
    (fObjOpt, paramOpt, flag) = optimize(opt1, x_init)
    print("k_o_ini=", k_o_ini, "\n\n")
    print("k_o_opt=", paramOpt[1])
end
