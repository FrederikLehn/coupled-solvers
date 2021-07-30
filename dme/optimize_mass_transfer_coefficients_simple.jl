using Roots, NLopt, PyPlot, LsqFit, JLD, Polynomials, JFVM, Dierckx
include("DME_solver.jl")

function rec_fact(param, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME)
  (k_o, k_w) = param
  (t, R) = DME_solver_no_gravity_no_pc_incompressible(G, K, MT, T, BC, IC, k_o, k_w, false)

  PV = t.*((u_inj./poros_val)./Lx)
  R_int = Spline1D(PV, R, k=1)
  return R_int(PV_exp_DME)
end

model = (PV_exp_DME, param) -> rec_fact(param, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME)

function error_calc(param, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME, R_exp_DME)
  return rec_fact(param, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME) - R_exp_DME
end

function error_calculation(param, param_ind, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME, R_exp_DME; grad_calc=true)
  eps1 = 1e-8
  error_vals = error_calc(param, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME, R_exp_DME)
  if grad_calc
    jac  = zeros(length(error_vals), length(param_ind))
  else
    jac = zeros(0)
  end
  if grad_calc
    for j in eachindex(param_ind)
      param2 = copy(param)
      param2[param_ind[j]]+=eps1
      error_val2 = error_calc(param2, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME, R_exp_DME)
      jac[:, j] = (error_val2 - error_vals) / eps1
    end
  end
  return error_vals, jac
end

function objective_function(param::Vector{Float64}, grad::Vector{Float64}, param_all, param_ind, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME, R_exp_DME_DME; w=1.0)
  param2 = copy(param_all)
  param2[param_ind] = param
  grad_calc = false
  if length(grad)>0
    grad_calc = true
  end

  error_val, jac = error_calculation(param2, param_ind, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME, R_exp_DME, grad_calc = grad_calc)

  obj_func_value = sum(abs2, w.*error_val)

  if length(grad)>0
      grad[:] = 2.0*sum(w.*error_val.*jac, 1)
  end
  return obj_func_value
end

obj_fun = (param, grad) -> objective_function(param::Vector{Float64}, grad::Vector{Float64}, param_all, param_ind, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp_DME, R_exp_DME, w=w)


experiment = 1

if experiment == 1
    perm_val  = 1e-15 # [m^2] permeability
elseif experiment == 2
    perm_val  = 0.8e-15 # [m^2] permeability
end

poros_val = 0.256     # [-] porosity
ρ_rock    = 2700.0     # [kg/m^3]

D_core = 0.025
A_core  = π*D_core^2/4

if experiment == 1
    injection_rate = 4.5202702702702705e-10 # 0.00011/(3600*24) # m^3/day
    u_inj   = injection_rate/A_core # Darcy velocity
elseif experiment == 2
    injection_rate = 4.5202702702702705e-10 #u_inj = 0.00015/(3600*24) m^3/day
    u_inj   = injection_rate/A_core # Darcy velocity
end


c_inj1 = 0.0 # mol/m^3 DME in injection water at initial water flooding
c_inj2 = 10.0 # mol/m^3 DME in injection water at DME flooding
c_tracer1 = 0.0
c_tracer2 = 1.0


Kval = c -> 2.0
dKval_dc = c -> zeros(size(c))
a = 2.0*ρ_rock*1000.0*(1.0-poros_val) # m^2/m^3 = 2 m^2/g * 2700 kg/m^3 * 1000 g/kg * (1-phi)
D_oil = 0.5*10^(-9.0) # m^2/s
D_water = 2*10^(-9.0) # m^2/s

## Creating DME parameters
ρ_oil = c_o -> 748.0015057043045-0.0016431963703508173.*c_o+2.9364412314092913.*10.0.^(-7.0).*c_o.^2.0-4.107225978803957.*10.0.^(-11.0).*c_o.^3.0
dρ_dc_oil = c_o -> -0.0016431963703508173+2.0.*2.9364412314092913.*10.0.^(-7.0).*c_o-3.0.*4.107225978803957.*10.0.^(-11.0).*c_o.^2.0
dρ_dp_oil = c_o -> zeros(size(c_o))
ρ_water = c_w -> 983.6319389796079.*ones(size(c_w))
dρ_dp_water = c_w -> zeros(size(c_w))
μ_oil = c_o -> 0.0010131130815849523-2.528969514668829.*10.0.^(-5.0).*c_o+2.3873417406201346.*10.0.^(-11.0).*c_o.^2.0
dμ_dc_oil = c_o -> -2.528969514668829.*10.0.^(-5.0)+2*2.3873417406201346.*10.0.^(-11.0).*c_o
μ_water = c_w -> 0.00040706344808758737+1.47269029810352587.*10.0.^(-7.0).*c_w
dμ_dc_water = c_w -> 1.47269029810352587.*10.0.^(-7.0).*ones(size(c_w))

MT = FluidParametersDME(ρ_oil, dρ_dc_oil, dρ_dp_oil, ρ_water, dρ_dp_water, μ_oil, dμ_dc_oil, μ_water, dμ_dc_water, Kval, dKval_dc, a, D_oil, D_water)

# Creating realitive permeability

if experiment == 1
    #swc  = 0.392349
    #kro0 = 0.121566
    #krw0 = 0.269128
    #n_w  = 1.88088
    #n_o  = 1.91705
    #sor_min  = 0.03
    #sor_max  = 0.250158
    (sor_max, swc, kro0min, krw0, n_o, n_w) = [0.198781  0.205555  0.249794  0.149609  3.03704  3.22222]
    sor_min  = 0.01 # equivalent of equilibrium mass frac with 0.1 mol frac DME-water
    kro0max = 0.95
elseif experiment == 2
    swc  = 0.317254
    kro0 = 0.563868
    krw0 = 0.366759
    n_o  = 2.38513
    n_w  = 2.40567
    sor_min  = 0.1
    sor_max  = 0.3184591
end


m_kro = -(kro0max-kro0min)/(sor_max-sor_min)
kro0_c = c -> m_kro*c + kro0max


c_min = 0.0 # min DME mol/m^3 in oil
c_max = Kval(0)*c_inj2 # max DME mol/m^3 fraction in oil
c_max = c_max[1]
m_sor = -(sor_max-sor_min)/(c_max-c_min)
sor_c = c -> m_sor*c+sor_max
dsor_c = c -> m_sor + 0.0*c #irrelevant

K = RelativePermeabilityDME(kro0_c, krw0, n_o, n_w, sor_c, dsor_c, swc)

# =============================================================================
# Preparing geology/mesh
# =============================================================================
dim = 1

srand(123456543)
λ_dp    = 0.8
cor_x   = 0.1
cor_y   = 0.1
cor_z   = 0.1
Lx      = 0.104 # [m] length of core
Ly      = 0.025 # [m] diameter of core
Lz      = 0.025 # [m] diameter of core
Nx      = 50 # number of grids in the x direction
Ny      = 20 # number of grids in the y direction
Nz      = 10 # number of grids in the z direction

if dim == 1
    m           = createMesh1D(Nx, Lx)
    #perm_mat    = permfieldlogrnde(Nx, perm_val, λ_dp, cor_x)
elseif dim == 2
    m           = createMesh2D(Nx, Ny, Lx, Ly)
    #perm_mat    = permfieldlogrnde(Nx, Ny, perm_val, λ_dp, cor_x, cor_y)
else
    m           = createMesh3D(Nx, Ny, Nz, Lx, Ly, Lz)
    #perm_mat    = permfieldlogrnde(Nx, Ny, Nz, perm_val, λ_dp, cor_x, cor_y, cor_z)
end

k       = createCellVariable(m, perm_val)
k_face  = harmonicMean(k)     # permeability on the cell faces
ϕ       = createCellVariable(m, poros_val)

G = Geology(m, perm_val, poros_val, k, ϕ)


# =============================================================================
# Preparing time-stepping parameters
# =============================================================================

if experiment == 1
    n_pv_water1 = 9.0 # number of injected pore volumes
elseif experiment == 2
    n_pv_water1 = 13.3
end
t_final_water1  = n_pv_water1*Lx/(u_inj/poros_val) # [s] final time for water flood
dt_water1       = t_final_water1/n_pv_water1/Nx/5 # [s] time step

if experiment == 1
    n_pv_DME = 9.0
elseif experiment == 2
    n_pv_DME = 0.3
end
t_final_DME     = t_final_water1+n_pv_DME*Lx/(u_inj/poros_val) # [s] final time for DME flood
dt_DME          = t_final_DME/(n_pv_water1 + n_pv_DME)/Nx/5 # [s] time step

if experiment == 1
    n_pv_water2 = 3.8
elseif experiment == 2
    n_pv_water2 = 3.1
end

t_final_water2  = t_final_DME+n_pv_water2*Lx/(u_inj/poros_val) # [s] final time for DME flood
dt_water2       = t_final_water2/(n_pv_water1 + n_pv_DME + n_pv_water2)/Nx/5 # [s] time step

T_water1 = SimulationTime(t_final_water1, dt_water1)
T_DME    = SimulationTime(t_final_DME, dt_DME)
T_water2 = SimulationTime(t_final_water2, dt_water2)

T = [T_water1, T_DME, T_water2]

## Boundary Conditions for water flooding
#T0 = 70 + 273.15 # [K]

c_oil0      = 0.0 # initial DME concentration in the oil phase
c_water0    = 0.0 # initial DME concentration in the water phase
p0          = 2000/14.7*1e5   # [Pa]

if experiment == 1
    sw0 = 1.0-0.667 # [vol frac] initial water saturation
elseif experiment == 2
    sw0 = 1.0-0.738
end

c_t0 = 0.0 # initial tracer concentration

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
BCc_oil_water.left.c[:]  = Kval(0)*c_inj1 #injected concentration of DME in water

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
BCc_oil_DME.left.c[:] = Kval(0)*c_inj2 #Created from reaction at boundary

BCp_DME.right.a[:] = 0.0
BCp_DME.right.b[:] = 1.0 #Dirichlet
BCp_DME.right.c[:] = p0

BCp_DME.left.a[:]  = k_face.xvalue[1,:]/μ_water(c_inj2) #Neumann
BCp_DME.left.b[:]  = 0.0
BCp_DME.left.c[:]  = -u_inj

BCs_DME.left.a[:]  = 0.0
BCs_DME.left.b[:]  = 1.0 #Dirichlet
BCs_DME.left.c[:]  = 1.0 #Value S_w=1

BCt_DME.left.a[:] = 0.0
BCt_DME.left.b[:] = 1.0 #Dirichlet
BCt_DME.left.c[:] = c_tracer2

BC_DME = BoundaryConditionDME(BCc_oil_DME, BCc_water_DME, BCp_DME, BCs_DME, BCt_DME)

BC = [BC_water, BC_DME, BC_water]


## Experimental results
if experiment == 1
    data = [0, 0.03883495145631066,
            0.18404907975460127, 0.23300970873786409,
            0.3067484662576685, 0.4271844660194175,
            0.42944785276073616, 0.4951456310679613,
            0.8588957055214725, 0.5048543689320388,
            0.9815950920245402, 0.5242718446601942,
            1.288343558282209, 0.5339805825242718,
            1.5950920245398772, 0.5436893203883496,
            2.2085889570552153, 0.5631067961165048,
            2.883435582822086, 0.5825242718446602,
            3.619631901840491, 0.5825242718446602,
            4.294478527607362, 0.5922330097087378,
            4.969325153374233, 0.5922330097087378,
            5.644171779141105, 0.6019417475728156,
            6.380368098159509, 0.6019417475728156,
            7.0552147239263805, 0.6019417475728156,
            7.730061349693252, 0.6019417475728156,
            8.404907975460121, 0.6019417475728156,
            9.079754601226991, 0.6019417475728156,
            9.570552147239262, 0.6019417475728156,
            9.938650306748466, 0.6116504854368932,
            10.061349693251532, 0.6310679611650485,
            10.184049079754601, 0.6893203883495146,
            10.306748466257668, 0.7475728155339806,
            10.429447852760735, 0.7669902912621359,
            10.736196319018404, 0.8155339805825242,
            10.920245398773005, 0.8252427184466019,
            11.16564417177914, 0.8446601941747572,
            11.411042944785276, 0.854368932038835,
            11.717791411042944, 0.8640776699029126,
            12.085889570552146, 0.8737864077669903,
            12.760736196319018, 0.9029126213592233,
            13.43558282208589, 0.9223300970873787,
            14.478527607361961, 0.941747572815534,
            15.0920245398773, 0.9514563106796117,
            15.766871165644169, 0.9514563106796117,
            16.44171779141104, 0.9514563106796117,
            17.116564417177912, 0.9514563106796117,
            17.791411042944784, 0.9514563106796117,
            18.343558282208587, 0.9514563106796117,
            18.71165644171779, 0.9514563106796117,
            19.01840490797546, 0.9514563106796117,
            19.386503067484657, 0.9514563106796117,
            20.06134969325153, 0.9514563106796117,
            20.7361963190184, 0.9514563106796117,
            21.411042944785272, 0.9514563106796117]

elseif experiment == 2
    data = [0.046242774566473965, 0.08333333333333326,
            0.18497109826589608, 0.26041666666666663,
            0.323699421965318, 0.4375,
            0.601156069364162, 0.5,
            0.8786127167630056, 0.5,
            1.0635838150289019, 0.5,
            1.2023121387283238, 0.5312499999999999,
            1.38728323699422, 0.5312499999999999,
            1.6647398843930634, 0.5312499999999999,
            1.8497109826589597, 0.5312499999999999,
            2.0346820809248554, 0.5312499999999999,
            2.6820809248554918, 0.5520833333333333,
            3.421965317919075, 0.5520833333333333,
            4.069364161849711, 0.5520833333333333,
            4.809248554913295, 0.5520833333333333,
            5.456647398843931, 0.5520833333333333,
            6.1502890173410405, 0.5520833333333333,
            6.84393063583815, 0.5520833333333333,
            7.5375722543352595, 0.5520833333333333,
            8.23121387283237, 0.5520833333333333,
            8.92485549132948, 0.5520833333333333,
            9.61849710982659, 0.5520833333333333,
            10.312138728323701, 0.5520833333333333,
            11.00578034682081, 0.5520833333333333,
            11.653179190751446, 0.5520833333333333,
            12.346820809248555, 0.5520833333333333,
            13.086705202312139, 0.5520833333333333,
            13.271676300578035, 0.5520833333333333,
            13.45664739884393, 0.5520833333333333,
            13.641618497109828, 0.5520833333333333,
            13.826589595375722, 0.5520833333333333,
            14.011560693641618, 0.5520833333333333,
            14.150289017341041, 0.5520833333333333,
            14.335260115606937, 0.5729166666666666,
            14.47398843930636, 0.6041666666666665,
            14.612716763005782, 0.625,
            14.751445086705203, 0.6354166666666665,
            14.890173410404625, 0.6458333333333333,
            15.075144508670519, 0.6458333333333333,
            15.260115606936417, 0.65625,
            15.44508670520231, 0.65625,
            15.676300578034684, 0.65625,
            15.861271676300579, 0.65625,
            16, 0.65625,
            16.739884393063583, 0.65625]
end

N_d = length(data)
PV_exp = zeros(N_d/2)
R_exp = zeros(N_d/2)
count_PV = 1
count_R = 1
for i = 1:N_d
    if mod(i,2) == 1
        PV_exp[count_PV] = data[i]
        count_PV += 1
    elseif mod(i,2) == 0
        R_exp[count_R] = data[i]
        count_R += 1
    end
end

PV_exp_DME = zeros(0)
R_exp_DME = zeros(0)
N_obs = length(PV_exp)
for i = 1:N_obs
    if PV_exp[i] >= n_pv_water1
        PV_exp_DME = push!(PV_exp_DME, PV_exp[i])
        R_exp_DME = push!(R_exp_DME, R_exp[i])
    end
end

# =============================================================================
# Calculating heuristic for k_o and k_w
# =============================================================================
d = sqrt(perm_val) #Characteristic length, i.e. diameter of capillary in porous medium
c_o = Kval(0)*c_inj2
c_o = c_o[1]
c_w = c_inj2
v_o = u_inj/poros_val
v_w = u_inj/poros_val

k_o_ini = 0.023*D_oil/d*((ρ_oil(c_o)*v_o*d)/μ_oil(c_o))^(0.83)*(μ_oil(c_o)/(ρ_oil(c_o)*D_oil))^(1.0/3.0)
k_w_ini = 0.023*D_water/d*((ρ_water(c_w)[1]*v_w*d)/μ_water(c_w))^(0.83)*(μ_water(c_w)/(ρ_water(c_w)[1]*D_water))^(1.0/3.0)

R_heu = rec_fact([k_o_ini, k_w_ini], G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp)

param_all = [1e-3*k_o_ini, 1e-3*k_w_ini]
param_ind = [1, 2]
w = ones(length(R_exp_DME))



x_init = copy(param_all)
x_lb = [1e-13, 1e-13]
x_ub = [1e-10, 1e-10]

# =============================================================================
# Optimizing
# =============================================================================

opt_type = 2

if opt_type == 1
    fit = curve_fit(model, PV_exp_DME, R_exp_DME, w, x_init, lower = x_lb, upper = x_ub)
    print("k_o_ini=", k_o_ini, ", k_w_ini=", k_w_ini, "\n\n")
    print("k_o_opt=", fit.param[1], ", k_w_opt=", fit.param[2])

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
    print("k_o_ini=", k_o_ini, ", k_w_ini=", k_w_ini, "\n\n")
    print("k_o_opt=", paramOpt[1], ", k_w_opt=", paramOpt[2])
end

R_opt = rec_fact(paramOpt, G, K, MT, T, BC, IC, Lx, u_inj, poros_val, PV_exp)


#Visualizing
figure()
plot(PV_exp, R_exp.*100)
plot(PV_exp, R_heu.*100)
plot(PV_exp, R_opt.*100)
plot([n_pv_water1, n_pv_water1], [0.0, 100.0], color=:black, linestyle="--")
plot([n_pv_water1 + n_pv_DME, n_pv_water1 + n_pv_DME], [0.0, 100.0], color=:black, linestyle="--")
axis([0, n_pv_water1 + n_pv_DME + n_pv_water2, 0, 100])
grid("on")
title("Experimental vs. Simulated Recovery Factor")
xlabel(L"$PVI$")
ylabel(L"$RF$"" [%]")
legend(["Experimental", "Heuristic", "Optimized"])
savefig(string("Optimization Results/ExpVsSim", experiment, ".png"))
