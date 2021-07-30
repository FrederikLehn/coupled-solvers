using Polynomials, Roots, JLD, PyPlot
using JFVM
include("peteng-master/functions/rel_perms_real.jl")

struct RelativePermeabilityDME
    kro0::Function
    krw0::Real
    n_o::Real
    n_w::Real
    sor::Function
    dsor_dc::Function #used in DME_solver_gravity_pc_compressible
    swc::Real
end

struct CapillaryPressureDME
    pc::Function
    dpc_dS::Function
    d2pc_d2S::Function
    d2pc_dcdS::Function
    theta::Real #Contact angle
    lambda::Real
    b::Real
    pc0::Real #maximum allowed pressure at low water saturations
end

struct Geology
    m::MeshStructure
    perm_val::Real
    poros_val::Real
    k::CellValue
    ϕ::CellValue
end

struct SimulationTime
    t_final::Real # [s] final time for simulation
    dt::Real # [s] time step
end

struct BoundaryConditionDME
    BCc_oil::BoundaryCondition # concentration (DME in oil) boundary condition
    BCc_water::BoundaryCondition # concentration (DME in water) boundary condition
    BCp::BoundaryCondition # pressure boundary condition
    BCs::BoundaryCondition # saturation boundary condition
    BCt::BoundaryCondition # tracer concentration boundary condition
end

struct InitialConditionDME
    c_oil_init::CellValue
    c_water_init::CellValue
    p_init::CellValue
    sw_init::CellValue
    c_t_init::CellValue
end

struct FluidParametersDME
    ρ_oil::Function
    dρ_dc_oil::Function
    dρ_dp_oil::Function #used in DME_solver_gravity_pc_compressible
    ρ_water::Function
    dρ_dp_water::Function #used in DME_solver_gravity_pc_compressible
    μ_oil::Function
    dμ_dc_oil::Function
    μ_water::Function
    dμ_dc_water::Function
    K::Function
    dK_dc::Function #used in DME_solver_gravity_pc_compressible
    a::Real
    D_oil::Real
    D_water::Real
end


function DME_solver_no_gravity_no_pc_incompressible(G::Geology, K::RelativePermeabilityDME, MT::FluidParametersDME, T::Array{SimulationTime,1}, BC::Array{BoundaryConditionDME,1}, IC::InitialConditionDME, k_o::Real, k_w::Real, plotFigures::Bool)

    # =============================================================================
    # Unpacking the relative permeabillity parameters
    # =============================================================================
    krw0 = K.krw0
    kro0_c = K.kro0 #function
    n_w  = K.n_w
    n_o  = K.n_o
    swc  = K.swc
    sor_c = K.sor #function

    KRW  = (sw, sor) -> krw.(sw, krw0, sor, swc, n_w)
    dKRW = (sw, sor) -> dkrwdsw.(sw, krw0, sor, swc, n_w)
    KRO  = (sw, kro0, sor) -> kro.(sw, kro0, sor, swc, n_o)
    dKRO = (sw, kro0, sor) -> dkrodsw.(sw, kro0, sor, swc, n_o)

    # =============================================================================
    # Unpacking the boundary and initial conditions
    # =============================================================================

    BCc_water   = BC[1].BCc_water # concentration (DME in water) boundary condition
    BCc_oil     = BC[1].BCc_oil # concentration (DME in oil) boundary condition
    BCp         = BC[1].BCp # pressure boundary condition
    BCs         = BC[1].BCs # saturation boundary condition
    BCt         = BC[1].BCt # tracer concentration boundary condition

    # discretize boundary conditions
    M_bc_c_water, RHS_bc_c_water    = boundaryConditionTerm(BCc_water)
    M_bc_c_oil, RHS_bc_c_oil        = boundaryConditionTerm(BCc_oil)
    M_bc_p, RHS_bc_p                = boundaryConditionTerm(BCp)
    M_bc_s, RHS_bc_s                = boundaryConditionTerm(BCs)
    M_bc_t, RHS_bc_t                = boundaryConditionTerm(BCt)

    # Iinitial conditions
    c_water_init    = IC.c_water_init
    c_oil_init      = IC.c_oil_init
    p_init          = IC.p_init
    sw_init         = IC.sw_init
    c_t_init        = IC.c_t_init

    # new values of each variable
    c_water_val     = copyCell(c_water_init)
    c_oil_val       = copyCell(c_oil_init)
    p_val           = copyCell(p_init)
    c_t_val         = copyCell(c_t_init)
    sw_val          = copyCell(sw_init)
    c_water_face    = arithmeticMean(c_water_val)
    c_oil_face      = arithmeticMean(c_oil_val)
    sor_face        = faceEval(sor_c, c_oil_face)

    # =============================================================================
    # Unpacking the fluid parameters
    # =============================================================================
    D_water = MT.D_water# m^2/s
    D_oil   = MT.D_oil # m^2/s

    ρ_water         = MT.ρ_water
    ρ_oil           = MT.ρ_oil
    dρ_dc_oil       = MT.dρ_dc_oil
    μ_water         = MT.μ_water
    dμ_dc_water     = MT.dμ_dc_water
    μ_oil           = MT.μ_oil
    dμ_dc_oil       = MT.dμ_dc_oil


    Kval = MT.K #Equilibrium constant

    #surface area only active if DME flooding is executed, otherwise 0
    DME_injected = false
    DME_injection_ended = false
    if maximum(BCc_oil.left.c[:]) > 0
        a = MT.a # m^2/m^3
        DME_injected = true
    else
        a = 0
    end

    K_w = (Kval(0)*k_o*k_w)/(Kval(0)*k_o+k_w)
    K_o = (k_o*k_w)/(Kval(0)*k_o+k_w)

    # =============================================================================
    # Unpacking the geology/mesh
    # =============================================================================
    perm_val    = G.perm_val      # [m^2] permeability
    poros_val   = G.poros_val     # [-] porosity
    k           = G.k
    k_face      = harmonicMean(k)     # permeability on the cell faces
    poros_val   = G.poros_val
    ϕ           = G.ϕ

    m   = G.m
    Nx  = m.dims[1]  # number of grids in the x direction
    Lx  = sum(m.cellsize.x[2:end-1]) # [m]
    N   = (Nx+2)

    if m.dimension>=2
        Ny  = m.dims[2]  # number of grids in the y direction
        Ly  = sum(m.cellsize.y[2:end-1]) # [m]
        N   = (Nx+2)*(Ny+2)
    end

    if m.dimension>=3
        Nz  = m.dims[3]  # number of grids in the y direction
        Lz  = sum(m.cellsize.z[2:end-1]) # [m]
        N   = (Nx+2)*(Ny+2)*(Ny+2)
    end


    # =============================================================================
    # Unpacking the time related parameters
    # =============================================================================

    t_final1  = T[1].t_final
    dt0       = T[1].dt
    dt        = dt0

    t_final_end  = T[end].t_final
    no_of_switch = length(T)

    # =============================================================================
    # Preparing for time-stepping
    # =============================================================================

    # outside the loop: initialization
    uw       = gradientTerm(p_val) # only for initialization of the water velocity vector
    uo       = gradientTerm(p_val) # only for initialization of the oil velocity vector
    ut       = 0 # initialize total velocity

    # Diffusion terms in the DME equations will be constant throughout the simulation
    M_c_water_diff  = diffusionTerm(createFaceVariable(m, poros_val*D_water))
    M_c_oil_diff    = diffusionTerm(createFaceVariable(m, poros_val*D_oil))

    # this whole thing goes inside two loops (Newton and time)
    tol_s = 1e-7
    tol_c = 1e-7
    max_change_c = 0.1 # 10 % relative change
    max_change_s = 0.1 # 10 % relative change
    max_int_loop = 40
    t = 0.0
    oil_init=domainInt(cellEval(ρ_oil, c_oil_val).*(1-sw_init)) # initial oil volume in the core
    rec_fact=zeros(1)
    dp_hist = zeros(1)
    t_s=zeros(1)
    switch_count = 2
    FL = fluxLimiter("SUPERBEE") # flux limiter


    while t < t_final_end
        # Switching injected fluid
        if switch_count <= no_of_switch && t > T[switch_count-1].t_final
            dt0 = T[switch_count].dt

            BCc_water   = BC[switch_count].BCc_water
            BCc_oil     = BC[switch_count].BCc_oil
            BCp         = BC[switch_count].BCp
            BCt         = BC[switch_count].BCt

            c_water_init    = createCellVariable(m, internalCells(c_water_val), BCc_water)
            c_oil_init      = createCellVariable(m, internalCells(c_oil_val), BCc_oil)
            p_val           = createCellVariable(m, internalCells(p_val), BCp)
            c_t_init        = createCellVariable(m, internalCells(c_t_val), BCt)

            M_bc_c_water, RHS_bc_c_water    = boundaryConditionTerm(BCc_water)
            M_bc_c_oil, RHS_bc_c_oil        = boundaryConditionTerm(BCc_oil)
            M_bc_p, RHS_bc_p                = boundaryConditionTerm(BCp)
            M_bc_t, RHS_bc_t                = boundaryConditionTerm(BCt)

            #changing surface area if DME flooding is beginning
            if maximum(BCc_oil.left.c[:]) > 0
                a = MT.a # m^2/m^3
                DME_injected = true
            else
                if DME_injected
                    DME_injection_ended = true
                end
            end

            switch_count += 1
        end

        error_s = 2*tol_s
        error_c = 2*tol_c
        loop_countor = 0

        #Newton-loop
        while error_s>tol_s || error_c_oil>tol_c || error_c_water>tol_c
            loop_countor += 1
            if loop_countor > max_int_loop
              c_water_val   = copyCell(c_water_init)
              c_oil_val     = copyCell(c_oil_init)
              p_val         = copyCell(p_init)
              sw_val        = copyCell(sw_init)
              c_t_val       = copyCell(c_t_init)
              dt = dt/3.0
              break
            end

            ∇p0 = gradientTerm(p_val)

            c_water_face  = upwindMean(c_water_val, uw)
            c_oil_face    = upwindMean(c_oil_val, uo)
            sw_face       = upwindMean(sw_val, uw)
            sor_face      = faceEval(sor_c, c_oil_face)
            sor_cell      = cellEval(sor_c, c_oil_val)
            kro0_face     = faceEval(kro0_c, sor_face)
            kro0_cell     = cellEval(kro0_c, sor_cell)

            ρ_water_cell        = cellEval(ρ_water, c_water_val)
            ρ_water_face        = faceEval(ρ_water, c_water_face)
            ρ_oil_cell          = cellEval(ρ_oil, c_oil_val)
            ρ_oil_face          = faceEval(ρ_oil, c_oil_face)
            dρ_dc_oil_cell      = cellEval(dρ_dc_oil, c_oil_val)
            dρ_dc_oil_face      = faceEval(dρ_dc_oil, c_oil_face)
            μ_water_face        = faceEval(μ_water, c_water_face)
            dμ_dc_water_face    = faceEval(dμ_dc_water, c_water_face)
            μ_oil_face          = faceEval(μ_oil, c_oil_face)
            dμ_dc_oil_face      = faceEval(dμ_dc_oil,c_oil_face)

            dc_μ_dc_water   = (μ_water_face-c_water_face.*dμ_dc_water_face)./(μ_water_face.*μ_water_face)
            d1_μ_dc_water   = (-dμ_dc_water_face)./(μ_water_face.*μ_water_face)
            dc_μ_dc_oil     = (μ_oil_face-c_oil_face.*dμ_dc_oil_face)./((μ_oil_face).*(μ_oil_face))
            dρ_μ_dc_oil    = (dρ_dc_oil_face.*μ_oil_face-dμ_dc_oil_face.*ρ_oil_face)./(μ_oil_face.*μ_oil_face)

            krw_face        = faceEval(KRW, sw_face, sor_face)
            kro_face        = faceEval(KRO, sw_face, kro0_face, sor_face)
            dkrw_dc_face    = faceEval(dKRW, sw_face, sor_face)
            dkro_dc_face    = faceEval(dKRO, sw_face, kro0_face, sor_face)
            uw              = -k_face.*krw_face./μ_water_face.*∇p0
            uo              = -k_face.*kro_face./μ_oil_face.*∇p0
            ut              = uw + uo

            #System matrices
            M_c_water_1     = M_bc_c_water  - M_c_water_diff  + linearSourceTerm(createCellVariable(m, K_w.*a.*poros_val))
            M_c_oil_1       =                                   linearSourceTerm(createCellVariable(m, -K_o.*a.*poros_val))
            M_c_water_2     =                                   linearSourceTerm(createCellVariable(m, -K_w.*a.*poros_val))
            M_c_oil_2       = M_bc_c_oil    - M_c_oil_diff    + linearSourceTerm(createCellVariable(m, K_o.*a.*poros_val))
            M_p_3           = M_bc_p
            M_S_4           = M_bc_s

            #Jacobians related to coupled reactive transport of water (row 1)
            J_t_c_water_1, RHS_t_c_water_1  = transientTerm(c_water_init, dt, ϕ.*sw_val)
            J_adv_c_water_1                 = convectionUpwindTerm(-k_face.*dc_μ_dc_water.*krw_face.*∇p0, ut)
            J_c_water_1                     = J_t_c_water_1 + J_adv_c_water_1

            J_p_1 = diffusionTerm(-k_face.*(1.0./μ_water_face).*krw_face.*c_water_face)

            J_t_S_1, RHS_t_S_1  = transientTerm(sw_init, dt, ϕ.*c_water_val)
            J_adv_S_1           = convectionUpwindTerm(-k_face.*(c_water_face./μ_water_face).*dkrw_dc_face.*∇p0, ut)
            J_S_1               = J_t_S_1 + J_adv_S_1

            #Jacobians related to coupled reactive transport of oil (row 2)
            J_t_c_oil_2, RHS_t_c_oil_2  = transientTerm(c_oil_init, dt, ϕ.*(1.0-sw_val))
            J_adv_c_oil_2               = convectionUpwindTerm(-k_face.*dc_μ_dc_oil.*kro_face.*∇p0, ut)
            J_c_oil_2                   = J_t_c_oil_2 + J_adv_c_oil_2

            J_p_2 = diffusionTerm(-k_face.*(1.0./μ_oil_face).*kro_face.*c_oil_face)

            J_t_S_2, RHS_t_S_2  = transientTerm(sw_init, dt, -ϕ.*c_oil_val)
            J_adv_S_2           = convectionUpwindTerm(-k_face.*(c_oil_face./μ_oil_face).*dkro_dc_face.*∇p0, ut)
            J_S_2               = J_t_S_2 + J_adv_S_2

            #Jacobians related to two-phase flow of water (row 3)
            J_adv_c_water_3                 = convectionUpwindTerm(-k_face.*d1_μ_dc_water.*krw_face.*ρ_water_face.*∇p0, ut)
            J_c_water_3                     = J_adv_c_water_3

            J_p_3 = diffusionTerm(-k_face.*(1.0./μ_water_face).*krw_face.*ρ_water_face)

            J_t_S_3, RHS_t_S_3  = transientTerm(sw_init, dt, ϕ.*ρ_water_cell)
            J_adv_S_3           = convectionUpwindTerm(-k_face.*(1.0./μ_water_face).*dkrw_dc_face.*ρ_water_face.*∇p0, ut)
            J_S_3               = J_t_S_3 + J_adv_S_3

            #Jacobians related to two-phase flow of oil (row 4)
            J_t_c_oil_4, RHS_t_c_oil_4  = transientTerm(c_oil_init, dt, ϕ.*(1.0-sw_val).*dρ_dc_oil_cell)
            J_adv_c_oil_4               = convectionUpwindTerm(-k_face.*dρ_μ_dc_oil.*kro_face.*∇p0, ut)
            J_c_oil_4                   = J_t_c_oil_4 + J_adv_c_oil_4

            J_p_4 = diffusionTerm(-k_face.*(ρ_oil_face./μ_oil_face).*kro_face)

            J_t_S_4, RHS_t_S_4  = transientTerm(sw_init, dt, -ϕ.*ρ_oil_cell)
            J_adv_S_4           = convectionUpwindTerm(-k_face.*(ρ_oil_face./μ_oil_face).*dkro_dc_face.*∇p0, ut)
            J_S_4               = J_t_S_4 + J_adv_S_4

            # create the PDE system M [c_water,c_oil,p,S]=RHS
            # x_val = [c_water_val.value[:], c_oil_val.value[:], p_val.value[:]; sw_val.value[:]]
            M = [M_c_water_1+J_c_water_1 M_c_oil_1 J_p_1 J_S_1;
                 M_c_water_2 M_c_oil_2+J_c_oil_2 J_p_2 J_S_2;
                 J_c_water_3 zeros(N,N) M_p_3+J_p_3 J_S_3;
                 zeros(N,N) J_c_oil_4 J_p_4 M_S_4+J_S_4]

            RHS_1 = RHS_bc_c_water  + RHS_t_c_water_1   + RHS_t_S_1     + J_adv_c_water_1*c_water_val.value[:]  + J_adv_S_1*sw_val.value[:]
            RHS_2 = RHS_bc_c_oil    + RHS_t_c_oil_2     + RHS_t_S_2     + J_adv_c_oil_2*c_oil_val.value[:]      + J_adv_S_2*sw_val.value[:]
            RHS_3 = RHS_bc_p                            + RHS_t_S_3     + J_adv_c_water_3*c_water_val.value[:]  + J_adv_S_3*sw_val.value[:]
            RHS_4 = RHS_bc_s        + RHS_t_c_oil_4     + RHS_t_S_4     + J_adv_c_oil_4*c_oil_val.value[:]      + J_adv_S_4*sw_val.value[:]

            RHS = [RHS_1[:];
                   RHS_2[:];
                   RHS_3[:];
                   RHS_4[:]]

            x_sol = M\RHS

            c_water_new = x_sol[1:N]
            c_oil_new   = x_sol[(1:N)+N]
            p_new       = x_sol[(1:N)+(2*N)]
            s_new       = x_sol[(1:N)+(3*N)]

            error_c_water = sumabs(c_water_new[2:end-1]-c_water_val.value[2:end-1])
            error_c_oil = sumabs(c_oil_new[2:end-1]-c_oil_val.value[2:end-1])
            error_s = sumabs(s_new[2:end-1]-sw_val.value[2:end-1])

            c_water_val.value[:]    = c_water_new[:]
            c_oil_val.value[:]      = c_oil_new[:]
            p_val.value[:]          = p_new[:]

            #Appleyard approach to updating saturations
            if !DME_injection_ended
                dsw_apple = 0.2
                eps_apple = sqrt(eps()) # 1e-5
                for i in eachindex(s_new)
                    if s_new[i]>=(1-sor_cell.value[i])
                        if sw_val.value[i]<(1-sor_cell.value[i]-eps_apple)
                            sw_val.value[i] = 1-sor_cell.value[i]-eps_apple
                        else
                            sw_val.value[i] = 1-sor_cell.value[i]
                        end
                    elseif s_new[i]<=swc
                        if sw_val.value[i]> swc+eps_apple
                            sw_val.value[i] = swc+eps_apple
                        else
                            sw_val.value[i] = swc
                        end
                    elseif abs(s_new[i]-sw_val.value[i])>dsw_apple
                        sw_val.value[i] += dsw_apple*sign(s_new[i]-sw_val.value[i])
                    else
                        sw_val.value[i] = s_new[i]
                    end
                end
            else
                w_c_w = 0.6
                w_c_o = 0.6
                w_p = 0.9
                w_sw = 0.6
                c_water_val.value[:] = w_c_w*c_water_new[:]+(1-w_c_w)*c_water_val.value[:]
                c_oil_val.value[:] = w_c_o*c_oil_new[:]+(1-w_c_o)*c_oil_val.value[:]
                p_val.value[:] = w_p*p_new[:]+(1-w_p)*p_val.value[:]
                sw_val.value[:] = w_sw*s_new[:]+(1-w_sw)*sw_val.value[:]
            end

            sw_val = createCellVariable(m, internalCells(sw_val), BCs)
        end

        if loop_countor < max_int_loop
            for tracer_count = 1:4
                # tracer calculation
                M_t_t, RHS_t_t = transientTerm(c_t_init, dt, sw_val.*ϕ)
                M_s_t          = linearSourceTerm(ϕ.*(sw_val-sw_init)/dt)
                M_a_t, RHS_a_t = convectionTvdTerm(uw, c_t_val, FL, ut)
                c_t_val = solveLinearPDE(m, M_t_t+M_s_t+M_a_t+M_bc_t, RHS_a_t+RHS_bc_t+RHS_t_t)
            end
            c_water_init    = copyCell(c_water_val)
            c_oil_init      = copyCell(c_oil_val)
            p_init          = copyCell(p_val)
            sw_init         = copyCell(sw_val)
            c_t_init        = copyCell(c_t_val)

            sor_face        = faceEval(sor_c, c_oil_face)
            kro0_face       = faceEval(kro0_c, sor_face)
            t += dt
            dt = dt0
            #println("progress: $(t/t_final_end*100) [%]")
            rec_fact=push!(rec_fact, (oil_init-domainInt(cellEval(ρ_oil,c_oil_val).*(1-sw_val)))/oil_init)
            t_s=push!(t_s, t)
        end
    end

    #Rounding off solutions, to avoid Julia plotting issues
    c_water_init.value[:] = round(c_water_init.value[:],2)
    c_oil_init.value[:] = round(c_oil_init.value[:],2)
    c_t_init.value[:] = round(c_t_init.value[:],2)

    if m.dimension >= 1 && m.dimension < 2 && plotFigures
        plot_DME_1D(k, sor_c, kro0_c, KRO, KRW, BC_DME, c_water_init, c_oil_init, p_init, sw_init, c_t_init, t_s, rec_fact)
    elseif m.dimension >= 2 && m.dimension < 3 && plotFigures
        plot_DME_2D(k, sor_c, kro0_c, KRO, KRW, BC_DME, c_water_init, c_oil_init, p_init, sw_init, c_t_init, t_s, rec_fact)
    elseif plotFigures
        plot_DME_3D(k, KRO, KRW, BC_DME, c_water_init, c_oil_init, p_init, sw_init, c_t_init, t_s, rec_fact)
    end

    println("Finished Simulation")
    return t_s, rec_fact
end

function DME_solver_gravity_pc_compressible(G::Geology, K::RelativePermeabilityDME, PC::CapillaryPressureDME, MT::FluidParametersDME, T::Array{SimulationTime,1}, BC::Array{BoundaryConditionDME,1}, IC::InitialConditionDME, k_o::Real, k_w::Real, plotFigures::Bool)
    # =============================================================================
    # Unpacking the relative permeabillity parameters
    # =============================================================================
    krw0    = K.krw0
    kro0    = K.kro0
    n_w     = K.n_w
    n_o     = K.n_o
    swc     = K.swc
    sor_c   = K.sor
    dsor_c  = K.dsor_dc

    KRW             = (sw, sor_c) -> krw.(sw, krw0, sor_c, swc, n_w)
    dKRW_dc_oil     = (sw, sor_c, dsor_c) -> dkrw_dc_oil.(sw, krw0, sor_c, dsor_c, swc, n_w)
    dKRW_dS         = (sw, sor_c) -> dkrwdsw.(sw, krw0, sor_c, swc, n_w)
    KRO             = (sw, sor_c) -> kro.(sw, kro0, sor_c, swc, n_o)
    dKRO_dc_oil     = (sw, sor_c, dsor_c) -> dkro_dc_oil.(sw, kro0, sor_c, dsor_c, swc, n_o)
    dKRO_dS         = (sw, sor_c) -> dkrodsw.(sw, kro0, sor_c, swc, n_o)

    # =============================================================================
    # Unpacking the capillary pressure functions
    # =============================================================================

    pc = PC.pc
    dpc_dS = PC.dpc_dS
    d2pc_d2S = PC.d2pc_d2S
    d2pc_dcdS = PC.d2pc_dcdS

    # =============================================================================
    # Unpacking the boundary and initial conditions
    # =============================================================================

    BCc_water   = BC[1].BCc_water # concentration (DME in water) boundary condition
    BCc_oil     = BC[1].BCc_oil # concentration (DME in oil) boundary condition
    BCp         = BC[1].BCp # pressure boundary condition
    BCs         = BC[1].BCs # saturation boundary condition
    BCt         = BC[1].BCt # tracer concentration boundary condition

    # discretize boundary conditions
    M_bc_c_water, RHS_bc_c_water    = boundaryConditionTerm(BCc_water)
    M_bc_c_oil, RHS_bc_c_oil        = boundaryConditionTerm(BCc_oil)
    M_bc_p, RHS_bc_p                = boundaryConditionTerm(BCp)
    M_bc_s, RHS_bc_s                = boundaryConditionTerm(BCs)
    M_bc_t, RHS_bc_t                = boundaryConditionTerm(BCt)

    # Iinitial conditions
    c_water_init    = IC.c_water_init
    c_oil_init      = IC.c_oil_init
    p_init          = IC.p_init
    sw_init         = IC.sw_init
    c_t_init        = IC.c_t_init

    # new values of each variable
    c_water_val     = copyCell(c_water_init)
    c_oil_val       = copyCell(c_oil_init)
    p_val           = copyCell(p_init)
    c_t_val         = copyCell(c_t_init)
    sw_val          = copyCell(sw_init)
    c_water_face    = arithmeticMean(c_water_val)
    c_oil_face      = arithmeticMean(c_oil_val)
    sor_face        = faceEval(sor_c, c_oil_face)

    # =============================================================================
    # Unpacking the fluid parameters
    # =============================================================================
    D_water = MT.D_water# m^2/s
    D_oil   = MT.D_oil # m^2/s

    ρ_water         = MT.ρ_water
    dρ_dp_water     = MT.dρ_dp_water
    ρ_oil           = MT.ρ_oil
    dρ_dc_oil       = MT.dρ_dc_oil
    dρ_dp_oil       = MT.dρ_dp_oil
    μ_water         = MT.μ_water
    dμ_dc_water     = MT.dμ_dc_water
    μ_oil           = MT.μ_oil
    dμ_dc_oil       = MT.dμ_dc_oil

    Kval = MT.K #Equilibrium constant
    dKval_dc = MT.dK_dc

    #surface area only active if DME flooding is executed, otherwise 0
    DME_injected = false
    DME_injection_ended = false
    if maximum(BCc_oil.left.c[:]) > 0
        a = MT.a # m^2/m^3
        DME_injected = true
    else
        a = 0
    end

    K_w = c -> (Kval(c).*k_o.*k_w)./(Kval(c).*k_o+k_w)
    K_o = c -> (k_o.*k_w)./(Kval(c).*k_o+k_w)


    # =============================================================================
    # Unpacking the geology/mesh
    # =============================================================================
    perm_val    = G.perm_val      # [m^2] permeability
    poros_val   = G.poros_val     # [-] porosity
    k           = G.k
    k_face      = harmonicMean(k)     # permeability on the cell faces
    poros_val   = G.poros_val
    ϕ           = G.ϕ

    m       = G.m
    Nx      = m.dims[1]  # number of grids in the x direction
    Lx      = sum(m.cellsize.x[2:end-1]) # [m]
    N       = (Nx+2)
    g       = createFaceVariable(m, -9.82) #gravity in x-direction

    if m.dimension>=2
        Ny      = m.dims[2]  # number of grids in the y direction
        Ly      = sum(m.cellsize.y[2:end-1]) # [m]
        N       = (Nx+2)*(Ny+2)
        g       = createFaceVariable(m, [0.0, -9.82]) #gravity in y-direction
    end

    if m.dimension>=3
        Nz      = m.dims[3]  # number of grids in the y direction
        Lz      = sum(m.cellsize.z[2:end-1]) # [m]
        N       = (Nx+2)*(Ny+2)*(Ny+2)
        g       = createFaceVariable(m, [0.0, 0.0, -9.82]) #gravity in z-direction
    end


    # =============================================================================
    # Unpacking the time related parameters
    # =============================================================================

    t_final1  = T[1].t_final
    dt0       = T[1].dt
    dt        = dt0

    t_final_end  = T[end].t_final
    no_of_switch = length(T)

    # =============================================================================
    # Preparing for time-stepping
    # =============================================================================

    # outside the loop: initialization
    uw       = gradientTerm(p_val) # only for initialization of the water velocity vector
    uo       = gradientTerm(p_val) # only for initialization of the oil velocity vector
    ut       = 0 # initialize total velocity

    # Diffusion terms in the DME equations will be constant throughout the simulation
    M_c_water_diff  = diffusionTerm(createFaceVariable(m, poros_val*D_water))
    M_c_oil_diff    = diffusionTerm(createFaceVariable(m, poros_val*D_oil))

    # this whole thing goes inside two loops (Newton and time)
    tol_s = 1e-7
    tol_c = 1e-7
    max_change_c = 0.1 # 10 % relative change
    max_change_s = 0.1 # 10 % relative change
    max_int_loop = 15
    t = 0.0
    oil_init=domainInt(cellEval(ρ_oil, c_oil_val, p_val).*(1-sw_init)) # initial oil volume in the core
    rec_fact=zeros(1)
    dp_hist = zeros(1)
    t_s=zeros(1)
    switch_count = 2
    FL = fluxLimiter("SUPERBEE") # flux limiter

    while t < t_final_end
        # Switching injected fluid
        if switch_count <= no_of_switch && t > T[switch_count-1].t_final
            dt0 = T[switch_count].dt

            BCc_water   = BC[switch_count].BCc_water
            BCc_oil     = BC[switch_count].BCc_oil
            BCp         = BC[switch_count].BCp
            BCt         = BC[switch_count].BCt

            c_water_init    = createCellVariable(m, internalCells(c_water_val), BCc_water)
            c_oil_init      = createCellVariable(m, internalCells(c_oil_val), BCc_oil)
            p_val           = createCellVariable(m, internalCells(p_val), BCp)
            c_t_init        = createCellVariable(m, internalCells(c_t_val), BCt)

            M_bc_c_water, RHS_bc_c_water    = boundaryConditionTerm(BCc_water)
            M_bc_c_oil, RHS_bc_c_oil        = boundaryConditionTerm(BCc_oil)
            M_bc_p, RHS_bc_p                = boundaryConditionTerm(BCp)
            M_bc_t, RHS_bc_t                = boundaryConditionTerm(BCt)

            #changing surface area if DME flooding is executed, otherwise 0
            if maximum(BCc_oil.left.c[:]) > 0
                a = MT.a # m^2/m^3
                DME_injected = true
            else
                if DME_injected
                    DME_injection_ended = true
                end
            end

            switch_count += 1
        end

        error_s = 2*tol_s
        error_c = 2*tol_c
        loop_countor = 0

        #System matrices
        M_c_water_1     = M_bc_c_water  - M_c_water_diff
        M_c_oil_2       = M_bc_c_oil    - M_c_oil_diff
        M_p_3           = M_bc_p
        M_S_4           = M_bc_s

        #Newton-loop
        while error_s>tol_s || error_c_oil>tol_c || error_c_water>tol_c
            loop_countor += 1
            if loop_countor > max_int_loop
              c_water_val   = copyCell(c_water_init)
              c_oil_val     = copyCell(c_oil_init)
              p_val         = copyCell(p_init)
              sw_val        = copyCell(sw_init)
              c_t_val       = copyCell(c_t_init)
              dt = dt/3.0
              break
            end

            ∇p0 = gradientTerm(p_val)
            ∇S0 = gradientTerm(sw_val)

            c_water_face  = upwindMean(c_water_val, uw)
            c_oil_face    = upwindMean(c_oil_val, uo)
            p_face        = upwindMean(p_val, uw) #CORRECT TO USE uw?
            sw_face       = upwindMean(sw_val, uw) #CORRECT TO USE uw?
            sor_cell      = cellEval(sor_c, c_oil_val)
            sor_face      = faceEval(sor_c, c_oil_face)
            dsor_cell     = cellEval(dsor_c, c_oil_val)
            dsor_face     = faceEval(dsor_c, c_oil_face)

            ρ_water_cell        = cellEval(ρ_water, p_val)
            ρ_water_face        = faceEval(ρ_water, p_face)
            dρ_dp_water_cell    = cellEval(dρ_dp_water, p_val)
            dρ_dp_water_face    = faceEval(dρ_dp_water, p_face)
            ρ_oil_cell          = cellEval(ρ_oil, c_oil_val, p_val)
            ρ_oil_face          = faceEval(ρ_oil, c_oil_face, p_face)
            dρ_dc_oil_cell      = cellEval(dρ_dc_oil, c_oil_val, p_val)
            dρ_dc_oil_face      = faceEval(dρ_dc_oil, c_oil_face, p_face)
            dρ_dp_oil_cell      = cellEval(dρ_dp_oil, c_oil_val, p_val)
            dρ_dp_oil_face      = faceEval(dρ_dp_oil, c_oil_face, p_face)
            μ_water_cell        = cellEval(μ_water, c_water_val)
            μ_water_face        = faceEval(μ_water, c_water_face)
            dμ_dc_water_face    = faceEval(dμ_dc_water, c_water_face)
            μ_oil_cell          = cellEval(μ_oil, c_oil_val)
            μ_oil_face          = faceEval(μ_oil, c_oil_face)
            dμ_dc_oil_face      = faceEval(dμ_dc_oil,c_oil_face)

            krw_cell            = cellEval(KRW, sw_val, sor_cell)
            krw_face            = faceEval(KRW, sw_face, sor_face)
            kro_cell            = cellEval(KRO, sw_val, sor_cell)
            kro_face            = faceEval(KRO, sw_face, sor_face)
            dkrw_dc_face        = faceEval(dKRW_dc_oil, sw_face, sor_face, dsor_face)
            dkro_dc_face        = faceEval(dKRO_dc_oil, sw_face, sor_face, dsor_face)
            dkrw_dS_face        = faceEval(dKRW_dS, sw_face, sor_face)
            dkro_dS_face        = faceEval(dKRO_dS, sw_face, sor_face)

            #dpc_dS_face         = faceEval(dpc_dS, sw_face, sor_face)
            #d2pc_d2S_face       = faceEval(d2pc_d2S, sw_face, sor_face)
            #d2pc_dcdS_face      = faceEval(d2pc_dcdS, sw_face, sor_face, dsor_face)

            Kval_cell     = cellEval(Kval, c_oil_val)
            Kval_face     = faceEval(Kval, c_oil_face)
            dKval_dc_cell = cellEval(dKval_dc, c_oil_val)
            dKval_dc_face = faceEval(dKval_dc, c_oil_face)
            K_w_cell      = cellEval(K_w, c_oil_val)
            K_w_face      = faceEval(K_w, c_oil_face)
            K_o_cell      = cellEval(K_o, c_oil_val)
            K_o_face      = faceEval(K_o, c_oil_face)

            #Derivatives of c_water
            dc_μ_dc_water = (μ_water_face-c_water_face.*dμ_dc_water_face)./(μ_water_face.*μ_water_face)
            d1_μ_dc_water = -dμ_dc_water_face./(μ_water_face.*μ_water_face)

            #Derivatives of c_oil
            dc_μ_kro_dc_oil     = kro_face./μ_oil_face-c_oil_face.*kro_face./(μ_oil_face.*μ_oil_face).*dμ_dc_oil_face+c_oil_face./μ_oil_face.*dkro_dc_face
            dρ_μ_kro_dc_oil     = kro_face./μ_oil_face.*dρ_dc_oil_face-ρ_oil_face.*kro_face./(μ_oil_face.*μ_oil_face).*dμ_dc_oil_face+ρ_oil_face./μ_oil_face.*dkro_dc_face
            dρ2_μ_kro_dc_oil    = 2.0.*ρ_oil_face.*kro_face./μ_oil_face.*dρ_dc_oil_face-(ρ_oil_face.*ρ_oil_face).*kro_face./(μ_oil_face.*μ_oil_face).*dμ_dc_oil_face+(ρ_oil_face.*ρ_oil_face)./μ_oil_face.*dkro_dc_face
            dρ_μ_c_kro_dc_oil   = c_oil_face.*kro_face./μ_oil_face.*dρ_dc_oil_face+ρ_oil_face.*kro_face./μ_oil_face-ρ_oil_face.*c_oil_face.*kro_face./(μ_oil_face.*μ_oil_face).*dμ_dc_oil_face+ρ_oil_face.*c_oil_face./μ_oil_face.*dkro_dc_face
            #dkro_μ_c_pcS_dc_oil = c_oil_face.*dpc_dS_face./μ_oil_face.*dkro_dc_face-kro_face.*c_oil_face.*dpc_dS_face./(μ_oil_face.*μ_oil_face).*dμ_dc_oil_face+kro_face.*dpc_dS_face./μ_oil_face+kro_face.*c_oil_face./μ_oil_face.*d2pc_dcdS_face
            #dkro_μ_ρ_pcS_dc_oil = ρ_oil_face.*dpc_dS_face./μ_oil_face.*dkro_dc_face-kro_face.*ρ_oil_face.*dpc_dS_face./(μ_oil_face.*μ_oil_face).*dμ_dc_oil_face+kro_face.*dpc_dS_face./μ_oil_face.*dρ_dc_oil_face+kro_face.*ρ_oil_face./μ_oil_face.*d2pc_dcdS_face
            dKw_dc              = k_o.*k_w./(Kval_cell.*k_o+k_w).*dKval_dc_cell-Kval_cell.*(k_o.*k_o).*k_w./((Kval_cell.*k_o+k_w).*(Kval_cell.*k_o+k_w)).*dKval_dc_cell
            dKo_c_dc            = k_o.*k_w./(Kval_cell.*k_o+k_w)-(k_o.*k_o).*k_w.*c_oil_val./((Kval_cell.*k_o+k_w).*(Kval_cell.*k_o+k_w)).*dKval_dc_cell

            #Derivatives of p
            dρ2_dp_water    = 2.0.*ρ_water_face.*dρ_dp_water_face
            dρ2_dp_oil      = 2.0.*ρ_oil_face.*dρ_dp_oil_face

            #Derivatives of S
            #dkro_pcS_dS_oil     = dkro_dc_face.*dpc_dS_face+kro_face.*d2pc_d2S_face

            #Phase velocities
            uw = -k_face.*krw_face./μ_water_face.*(∇p0-ρ_water_face.*g)
            uo = -k_face.*kro_face./μ_oil_face.*(∇p0-ρ_oil_face.*g)#+dpc_dS_face.*∇S0
            ut = uw + uo

            #println("uw=", mean(uw.xvalue), " , uo=", mean(uo.xvalue))

            #Jacobians related to coupled reactive transport of water (row 1)
            J_t_c_water_1, RHS_t_c_water_1  = transientTerm(c_water_init, dt, ϕ.*sw_val)
            J_adv_c_water_p_1               = convectionUpwindTerm(-k_face.*dc_μ_dc_water.*krw_face.*∇p0, uw)
            J_adv_c_water_g_1               = convectionUpwindTerm(-k_face.*ρ_water_face.*dc_μ_dc_water.*krw_face.*g, uw)
            J_lin_c_water_cw_1              = linearSourceTerm(createCellVariable(m, internalCells(K_w_cell).*a.*internalCells(ϕ)))
            J_c_water_1                     = J_t_c_water_1 + J_adv_c_water_p_1 - J_adv_c_water_g_1 + J_lin_c_water_cw_1

            J_adv_c_oil_p_1                 = convectionUpwindTerm(-k_face.*c_water_face./μ_water_face.*dkrw_dc_face.*∇p0, uw)
            J_adv_c_oil_g_1                 = convectionUpwindTerm(-k_face.*ρ_water_face.*c_water_face./μ_water_face.*dkrw_dc_face.*g, uw)
            J_lin_c_oil_cw_1                = linearSourceTerm(createCellVariable(m, internalCells(dKw_dc).*a.*internalCells(ϕ).*internalCells(c_water_val)))
            J_lin_c_oil_co_1                = linearSourceTerm(createCellVariable(m, internalCells(dKo_c_dc).*a.*internalCells(ϕ)))
            J_c_oil_1                       = J_adv_c_oil_p_1 - J_adv_c_oil_g_1 + J_lin_c_oil_cw_1 - J_lin_c_oil_co_1

            J_adv_p_g_1                     = convectionUpwindTerm(-k_face.*dρ_dp_water_face.*c_water_face./μ_water_face.*krw_face.*g, uw)
            J_diff_p_1                      = diffusionTerm(-k_face.*c_water_face./μ_water_face.*krw_face)
            J_p_1                           = -J_adv_p_g_1 + J_diff_p_1

            J_t_S_1, RHS_t_S_1              = transientTerm(sw_init, dt, ϕ.*c_water_val)
            J_adv_S_p_1                     = convectionUpwindTerm(-k_face.*(c_water_face./μ_water_face).*dkrw_dS_face.*∇p0, uw)
            J_adv_S_g_1                     = convectionUpwindTerm(-k_face.*ρ_water_face.*c_water_face./μ_water_face.*dkrw_dS_face.*g, uw)
            J_S_1                           = J_t_S_1 + J_adv_S_p_1 - J_adv_S_g_1

            #Jacobians related to coupled reactive transport of oil (row 2)
            J_lin_c_water_cw_2              = linearSourceTerm(createCellVariable(m, internalCells(K_w_cell).*a.*internalCells(ϕ)))
            J_c_water_2                     = -J_lin_c_water_cw_2

            J_t_c_oil_2, RHS_t_c_oil_2      = transientTerm(c_oil_init, dt, ϕ.*(1.0-sw_val))
            J_adv_c_oil_p_2                 = convectionUpwindTerm(-k_face.*dc_μ_kro_dc_oil.*∇p0, uo)
            J_adv_c_oil_pc_2                = zeros(N,N)#convectionUpwindTerm(-k_face.*dkro_μ_c_pcS_dc_oil.*∇S0, uo)
            J_adv_c_oil_g_2                 = convectionUpwindTerm(-k_face.*dρ_μ_c_kro_dc_oil.*g, uo)
            J_lin_c_oil_cw_2                = linearSourceTerm(createCellVariable(m, internalCells(dKw_dc).*a.*internalCells(ϕ).*internalCells(c_water_val)))
            J_lin_c_oil_co_2                = linearSourceTerm(createCellVariable(m, internalCells(dKo_c_dc).*a.*internalCells(ϕ)))
            J_c_oil_2                       = J_t_c_oil_2 + J_adv_c_oil_p_2 + J_adv_c_oil_pc_2 - J_adv_c_oil_g_2 - J_lin_c_oil_cw_2 + J_lin_c_oil_co_2

            J_adv_p_g_2                     = convectionUpwindTerm(-k_face.*dρ_dp_oil_face.*c_oil_face./μ_oil_face.*kro_face.*g, uo)
            J_diff_p_2                      = diffusionTerm(-k_face.*c_oil_face./μ_oil_face.*kro_face)
            J_p_2                           = -J_adv_p_g_2 + J_diff_p_2

            J_t_S_2, RHS_t_S_2              = transientTerm(sw_init, dt, -ϕ.*c_oil_val)
            J_adv_S_p_2                     = convectionUpwindTerm(-k_face.*c_oil_face./μ_oil_face.*dkro_dS_face.*∇p0, uo)
            J_adv_S_pc_2                    = zeros(N,N)#convectionUpwindTerm(-k_face.*c_oil_face./μ_oil_face.*dkro_pcS_dS_oil.*∇S0, uo)
            J_adv_S_g_2                     = convectionUpwindTerm(-k_face.*ρ_oil_face.*c_oil_face./μ_oil_face.*dkro_dS_face.*g, uo)
            J_diff_S_pc_2                   = zeros(N,N)#diffusionTerm(-k_face.*c_oil_face./μ_oil_face.*kro_face.*dpc_dS_face)
            J_S_2                           = J_t_S_2 + J_adv_S_p_2 + J_adv_S_pc_2 - J_adv_S_g_2  + J_diff_S_pc_2

            #Jacobians related to two-phase flow of water (row 3)
            J_adv_c_water_p_3               = convectionUpwindTerm(-k_face.*ρ_water_face.*d1_μ_dc_water.*krw_face.*∇p0, uw)
            J_adv_c_water_g_3               = convectionUpwindTerm(-k_face.*(ρ_water_face.*ρ_water_face).*d1_μ_dc_water.*krw_face.*g, uw)
            J_c_water_3                     = J_adv_c_water_p_3 - J_adv_c_water_g_3

            J_adv_c_oil_p_3                 = convectionUpwindTerm(-k_face.*ρ_water_face.*(1.0./μ_water_face).*dkrw_dc_face.*∇p0, uw)
            J_adv_c_oil_g_3                 = convectionUpwindTerm(-k_face.*(ρ_water_face.*ρ_water_face).*(1.0./μ_water_face).*dkrw_dc_face.*g, uw)
            J_c_oil_3                       = J_adv_c_oil_p_3 - J_adv_c_oil_g_3

            J_t_p_3, RHS_t_p_3              = transientTerm(p_init, dt, ϕ.*sw_val.*dρ_dp_water_cell)
            J_adv_p_p_3                     = convectionUpwindTerm(-k_face.*dρ_dp_water_face.*1.0./μ_water_face.*krw_face.*∇p0, uw)
            J_adv_p_g_3                     = convectionUpwindTerm(-k_face.*dρ2_dp_water.*(1.0./μ_water_face).*krw_face.*g, uw)
            J_diff_p_3                      = diffusionTerm(-k_face.*ρ_water_face.*(1.0./μ_water_face).*krw_face)
            J_p_3                           = J_t_p_3 + J_adv_p_p_3 - J_adv_p_g_3 + J_diff_p_3

            J_t_S_3, RHS_t_S_3              = transientTerm(sw_init, dt, ϕ.*ρ_water_cell)
            J_adv_S_p_3                     = convectionUpwindTerm(-k_face.*ρ_water_face.*(1.0./μ_water_face).*dkrw_dS_face.*∇p0, uw)
            J_adv_S_g_3                     = convectionUpwindTerm(-k_face.*(ρ_water_face.*ρ_water_face).*(1.0./μ_water_face).*dkrw_dS_face.*g, uw)
            J_S_3                           = J_t_S_3 + J_adv_S_p_3 - J_adv_S_g_3

            #Jacobians related to two-phase flow of oil (row 4)
            J_t_c_oil_4, RHS_t_c_oil_4      = transientTerm(c_oil_init, dt, ϕ.*(1.0-sw_val).*dρ_dc_oil_cell)
            J_adv_c_oil_p_4                 = convectionUpwindTerm(-k_face.*dρ_μ_kro_dc_oil.*∇p0, uo)
            J_adv_c_oil_pc_4                = zeros(N,N)#convectionUpwindTerm(-k_face.*dkro_μ_ρ_pcS_dc_oil.*∇S0, uo)
            J_adv_c_oil_g_4                 = convectionUpwindTerm(-k_face.*dρ2_μ_kro_dc_oil.*g, uo)
            J_c_oil_4                       = J_t_c_oil_4 + J_adv_c_oil_p_4 + J_adv_c_oil_pc_4 - J_adv_c_oil_g_4

            J_t_p_4, RHS_t_p_4              = transientTerm(p_init, dt, ϕ.*(1.0-sw_val).*dρ_dp_oil_cell)
            J_adv_p_p_4                     = convectionUpwindTerm(-k_face.*dρ_dp_oil_face./μ_oil_face.*kro_face.*∇p0, uo)
            J_adv_p_pc_4                    = zeros(N,N)#convectionUpwindTerm(-k_face.*dρ_dp_oil_face./μ_oil_face.*kro_face.*dpc_dS_face.*∇S0, uo)
            J_adv_p_g_4                     = convectionUpwindTerm(-k_face.*dρ2_dp_oil.*(1.0./μ_oil_face).*kro_face.*g, uo)
            J_diff_p_4                      = diffusionTerm(-k_face.*ρ_oil_face./μ_oil_face.*kro_face)
            J_p_4                           = J_t_p_4 + J_adv_p_p_4 + J_adv_p_pc_4 - J_adv_p_g_4 + J_diff_p_4

            J_t_S_4, RHS_t_S_4              = transientTerm(sw_init, dt, -ϕ.*ρ_oil_cell)
            J_adv_S_p_4                     = convectionUpwindTerm(-k_face.*(ρ_oil_face./μ_oil_face).*dkro_dS_face.*∇p0, uo)
            J_adv_S_pc_4                    = zeros(N,N)#convectionUpwindTerm(-k_face.*ρ_oil_face./μ_oil_face.*dkro_pcS_dS_oil.*∇S0, uo)
            J_adv_S_g_4                     = convectionUpwindTerm(-k_face.*(ρ_oil_face.*ρ_oil_face)./μ_oil_face.*dkro_dS_face.*g, uo)
            J_diff_S_pc_4                   = zeros(N,N)#diffusionTerm(-k_face.*ρ_oil_face./μ_oil_face.*kro_face.*dpc_dS_face)
            J_S_4                           = J_t_S_4 + J_adv_S_p_4 + J_adv_S_pc_4 - J_adv_S_g_4 + J_diff_S_pc_4

            # Constant terms for all equations
            RHS_const_c_water_cw_1    = constantSourceTerm(createCellVariable(m, internalCells(K_w_cell).*a.*internalCells(ϕ).*internalCells(c_water_val)))
            RHS_const_c_oil_co_1      = constantSourceTerm(createCellVariable(m, internalCells(K_o_cell).*a.*internalCells(ϕ).*internalCells(c_oil_val)))
            RHS_div_g_1             = divergenceTerm(-k_face.*ρ_water_face.*c_water_face./μ_water_face.*krw_face.*g)

            RHS_const_c_water_cw_2    = constantSourceTerm(createCellVariable(m, internalCells(K_w_cell).*a.*internalCells(ϕ).*internalCells(c_water_val)))
            RHS_const_c_oil_co_2      = constantSourceTerm(createCellVariable(m, internalCells(K_o_cell).*a.*internalCells(ϕ).*internalCells(c_oil_val)))
            RHS_div_g_2             = divergenceTerm(-k_face.*ρ_oil_face.*c_oil_face./μ_oil_face.*kro_face.*g)

            RHS_div_g_3             = divergenceTerm(-k_face.*(ρ_water_face.*ρ_water_face)./μ_water_face.*krw_face.*g)

            RHS_div_g_4             = divergenceTerm(-k_face.*(ρ_oil_face.*ρ_oil_face)./μ_oil_face.*kro_face.*g)


            # create the PDE system M [c_water,c_oil,p,S]=RHS
            # x_val = [c_water_val.value[:], c_oil_val.value[:], p_val.value[:]; sw_val.value[:]]
            M = [M_c_water_1+J_c_water_1 J_c_oil_1 J_p_1 J_S_1;
                 J_c_water_2 M_c_oil_2+J_c_oil_2 J_p_2 J_S_2;
                 J_c_water_3 zeros(N,N) M_p_3+J_p_3 J_S_3;
                 zeros(N,N) J_c_oil_4 J_p_4 M_S_4+J_S_4]

            #RHS related to coupled reactive transport of water
            RHS_1 = RHS_bc_c_water + RHS_t_c_water_1 + RHS_t_S_1 + (J_adv_c_water_p_1 - J_adv_c_water_g_1 + J_lin_c_water_cw_1)*c_water_val.value[:] + (J_adv_c_oil_p_1 - J_adv_c_oil_g_1 + J_lin_c_oil_cw_1-J_lin_c_oil_co_1)*c_oil_val.value[:] - J_adv_p_g_1*p_val.value[:] + (J_adv_S_p_1 - J_adv_S_g_1)*sw_val.value[:] - RHS_const_c_water_cw_1 + RHS_const_c_oil_co_1 + RHS_div_g_1

            #RHS related to coupled reactive transport of oil
            RHS_2 = RHS_bc_c_oil + RHS_t_c_oil_2 + RHS_t_S_2 - J_lin_c_water_cw_2*c_water_val.value[:] + (J_adv_c_oil_p_2 + J_adv_c_oil_pc_2 - J_adv_c_oil_g_2 - J_lin_c_oil_cw_2 + J_lin_c_oil_co_2)*c_oil_val.value[:] - J_adv_p_g_2*p_val.value[:] + (J_adv_S_p_2 + J_adv_S_pc_2 - J_adv_S_g_2)*sw_val.value[:] + RHS_const_c_water_cw_2 - RHS_const_c_oil_co_2 + RHS_div_g_2

            #RHS related to two-phase flow of water
            RHS_3 = RHS_bc_p + RHS_t_p_3 + RHS_t_S_3 + (J_adv_c_water_p_3 - J_adv_c_water_g_3)*c_water_val.value[:] + (J_adv_c_oil_p_3 - J_adv_c_oil_g_3)*c_oil_val.value[:] + (J_adv_p_p_3 - J_adv_p_g_3)*p_val.value[:] + (J_adv_S_p_3 - J_adv_S_g_3)*sw_val.value[:] + RHS_div_g_3

            #RHS related to two-phase flow of oil
            RHS_4 = RHS_bc_s + RHS_t_c_oil_4 + RHS_t_p_4 + RHS_t_S_4 + (J_adv_c_oil_p_4 + J_adv_c_oil_pc_4 - J_adv_c_oil_g_4)*c_oil_val.value[:] + (J_adv_p_p_4 + J_adv_p_pc_4 - J_adv_p_g_4)*p_val.value[:] + (J_adv_S_p_4 + J_adv_S_pc_4 - J_adv_S_g_4)*sw_val.value[:] + RHS_div_g_4

            RHS = [RHS_1[:];
                   RHS_2[:];
                   RHS_3[:];
                   RHS_4[:]]

            x_sol = M\RHS

            c_water_new = x_sol[1:N]
            c_oil_new   = x_sol[(1:N)+N]
            p_new       = x_sol[(1:N)+(2*N)]
            s_new       = x_sol[(1:N)+(3*N)]

            error_c_water = sumabs(c_water_new[2:end-1]-c_water_val.value[2:end-1])
            error_c_oil = sumabs(c_oil_new[2:end-1]-c_oil_val.value[2:end-1])
            error_s = sumabs(s_new[2:end-1]-sw_val.value[2:end-1])

            c_water_val.value[:]    = c_water_new[:]
            c_oil_val.value[:]      = c_oil_new[:]
            p_val.value[:]          = p_new[:]

            #Appleyard approach to updating saturations
            if !DME_injection_ended
                dsw_apple = 0.2
                eps_apple = sqrt(eps()) # 1e-5
                for i in eachindex(s_new)
                    if s_new[i]>=(1-sor_cell.value[i])
                        if sw_val.value[i]<(1-sor_cell.value[i]-eps_apple)
                            sw_val.value[i] = 1-sor_cell.value[i]-eps_apple
                        else
                            sw_val.value[i] = 1-sor_cell.value[i]
                        end
                    elseif s_new[i]<=swc
                        if sw_val.value[i]> swc+eps_apple
                            sw_val.value[i] = swc+eps_apple
                        else
                            sw_val.value[i] = swc
                        end
                    elseif abs(s_new[i]-sw_val.value[i])>dsw_apple
                        sw_val.value[i] += dsw_apple*sign(s_new[i]-sw_val.value[i])
                    else
                        sw_val.value[i] = s_new[i]
                    end
                end
            else
                w_c_w = 0.6
                w_c_o = 0.6
                w_p = 0.9
                w_sw = 0.6
                c_water_val.value[:] = w_c_w*c_water_new[:]+(1-w_c_w)*c_water_val.value[:]
                c_oil_val.value[:] = w_c_o*c_oil_new[:]+(1-w_c_o)*c_oil_val.value[:]
                p_val.value[:] = w_p*p_new[:]+(1-w_p)*p_val.value[:]
                sw_val.value[:] = w_sw*s_new[:]+(1-w_sw)*sw_val.value[:]
            end

            sw_val = createCellVariable(m, internalCells(sw_val), BCs)
        end

        if loop_countor < max_int_loop
            for tracer_count = 1:4
                # tracer calculation
                M_t_t, RHS_t_t = transientTerm(c_t_init, dt, sw_val.*ϕ)
                M_s_t          = linearSourceTerm(ϕ.*(sw_val-sw_init)/dt)
                M_a_t, RHS_a_t = convectionTvdTerm(uw, c_t_val, FL, ut)
                c_t_val = solveLinearPDE(m, M_t_t+M_s_t+M_a_t+M_bc_t, RHS_a_t+RHS_bc_t+RHS_t_t)
            end
            c_water_init    = copyCell(c_water_val)
            c_oil_init      = copyCell(c_oil_val)
            p_init          = copyCell(p_val)
            sw_init         = copyCell(sw_val)
            c_t_init        = copyCell(c_t_val)

            sor_face        = faceEval(sor_c, c_oil_face)
            t += dt
            dt = dt0
            println("progress: $(t/t_final_end*100) [%]")
            rec_fact=push!(rec_fact, (oil_init-domainInt(cellEval(ρ_oil,c_oil_val, p_val).*(1-sw_val)))/oil_init)
            t_s=push!(t_s, t)
        end
    end

    if m.dimension >= 1 && m.dimension < 2 && plotFigures
        plot_DME_1D(k, sor_c, kro0_c, KRO, KRW, BC_DME, c_water_init, c_oil_init, p_init, sw_init, c_t_init, t_s, rec_fact)
    elseif m.dimension >= 2 && m.dimension < 3 && plotFigures
        plot_DME_2D(k, KRO, KRW, BC_DME, c_water_init, c_oil_init, p_init, sw_init, c_t_init, t_s, rec_fact)
    elseif plotFigures
        plot_DME_3D(k, KRO, KRW, BC_DME, c_water_init, c_oil_init, p_init, sw_init, c_t_init, t_s, rec_fact)
    end

    return t_s, rec_fact
end

function plot_DME_1D(k::CellValue, sor_c::Function, kro0_c::Function, KRO::Function, KRW::Function, BC_DME::BoundaryConditionDME, c_water_init::CellValue, c_oil_init::CellValue, p_init::CellValue, sw_init::CellValue, c_t_init::CellValue, t_s, rec_fact)
    figure()
    grid("on")
    visualizeCells(k)
    title("Permeability")
    xlabel(L"$L_x$"" [m]")
    ylabel("Permeability, "L"$k$"" "L"[m^2]")
    savefig("1D Results/perm_field.png")

    figure()
    grid("on")
    visualizeCells(cellEval(sor_c, c_oil_init))
    ylim([0, 1])
    title("Residual Oil Saturation")
    xlabel(L"$L_x$"" [m]")
    ylabel("Residual Oil Saturation, "L"$S_{or}$"" "L"[-]")
    savefig("1D Results/residual_oil.png")

    figure()
    grid("on")
    sw_plot = collect(linspace(0,1,1000))
    sor = sor_c(0)
    plot(sw_plot, KRO.(sw_plot, kro0_c(sor), sor), sw_plot, KRW.(sw_plot, sor))
    plot([swc, swc], [0, 1], color=:black, linestyle="--")
    plot([1-sor, 1-sor], [0, 1], color=:black, linestyle="--")
    axis([0,1,0,1])
    title("Relative permeability, "L"$k_{ro},\ k_{rw}$"" [-]")
    xlabel(L"$S_w$"" [-]")
    ylabel(L"$k_r$"" [-]")
    legend([L"$k_{ro}$", L"$k_{rw}$"])
    savefig("1D Results/rel_perm_start.png")

    figure()
    grid("on")
    sw_plot = collect(linspace(0,1,1000))
    sor = sor_c(maximum(BC_DME.BCc_oil.left.c[:]))
    plot(sw_plot, KRO.(sw_plot, kro0_c(sor), sor), sw_plot, KRW.(sw_plot, sor))
    plot([swc, swc], [0, 1], color=:black, linestyle="--")
    plot([1-sor, 1-sor], [0, 1], color=:black, linestyle="--")
    axis([0,1,0,1])
    title("Relative permeability, "L"$k_{ro},\ k_{rw}$"" [-]")
    xlabel(L"$S_w$"" [-]")
    ylabel(L"$k_r$"" [-]")
    legend([L"$k_{ro}$", L"$k_{rw}$"])
    savefig("1D Results/rel_perm_end.png")

    figure()
    grid("on")
    visualizeCells(p_init/1e6) #MPa
    title("Pressure")
    xlabel(L"$L_x$"" [m]")
    ylabel("Pressure, "L"$p$"" [MPa]")
    savefig("1D Results/p_init.png")

    figure()
    grid("on")
    visualizeCells(sw_init)
    title("Water saturation")
    xlabel(L"$L_x$"" [m]")
    ylabel("Water saturation, "L"$S_w$"" [-]")
    savefig("1D Results/sw_init.png")

    figure()
    grid("on")
    c_oil_init.value[:] = round.(c_oil_init.value[:],3)
    visualizeCells(c_oil_init)
    title("Concentration of DME in oil")
    xlabel(L"$L_x$"" [m]")
    ylabel("Concentration of DME in oil, "L"$c_o$"" "L"\left[\frac{mol}{m^3}\right]")
    savefig("1D Results/c_oil_init.png")

    figure()
    grid("on")
    c_water_init.value[:] = round.(c_water_init.value[:],3)
    visualizeCells(c_water_init)
    title("Concentration of DME in water")
    xlabel(L"$L_x$"" [m]")
    ylabel("Concentration of DME in water, "L"$c_w$"" "L"\left[\frac{mol}{m^3}\right]")
    savefig("1D Results/c_water_init.png")

    figure()
    grid("on")
    c_t_init.value[:] = round.(c_t_init.value[:],3)
    visualizeCells(c_t_init)
    title("Tracer")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$c_t$")
    savefig("1D Results/ct_init.png")

    figure()
    grid("on")
    plot(t_s/(60*60*24), rec_fact*100)
    ylim([0, 100])
    title("Recovery factor")
    xlabel(L"$t$"" [days]")
    ylabel(L"$RF$"" [%]")
    savefig("1D Results/recovery_factor.png")
end

function plot_DME_2D(k::CellValue, sor_c::Function, kro0_c::Function, KRO::Function, KRW::Function, BC_DME::BoundaryConditionDME, c_water_init::CellValue, c_oil_init::CellValue, p_init::CellValue, sw_init::CellValue, c_t_init::CellValue, t_s, rec_fact)
    figure()
    visualizeCells(k)
    colorbar()
    title("Permeability, "L"$k$"" "L"[m^2]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("2D Results/perm_field.png")

    figure()
    sor = cellEval(sor_c, c_oil_init)
    sor.value[:] = round(sor.value[:],2)
    visualizeCells(sor)
    colorbar()
    title("Residual Oil Saturation, "L"$S_{or}$"" "L"[-]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("2D Results/residual_oil.png")

    figure()
    grid("on")
    sw_plot = collect(linspace(0,1,1000))
    sor = sor_c(0)
    plot(sw_plot, KRO.(sw_plot, kro0_c(sor), sor), sw_plot, KRW.(sw_plot, sor))
    plot([swc, swc], [0, 1], color=:black, linestyle="--")
    plot([1-sor, 1-sor], [0, 1], color=:black, linestyle="--")
    axis([0,1,0,1])
    title("Relative permeability, "L"$k_{ro},\ k_{rw}$"" [-]")
    xlabel(L"$S_w$"" [-]")
    ylabel(L"$k_r$"" [-]")
    legend([L"$k_{ro}$", L"$k_{rw}$"])
    savefig("2D Results/rel_perm_start.png")

    figure()
    grid("on")
    sw_plot = collect(linspace(0,1,1000))
    sor = sor_c(maximum(BC_DME.BCc_oil.left.c[:]))
    plot(sw_plot, KRO.(sw_plot, kro0_c(sor), sor), sw_plot, KRW.(sw_plot, sor))
    plot([swc, swc], [0, 1], color=:black, linestyle="--")
    plot([1-sor, 1-sor], [0, 1], color=:black, linestyle="--")
    axis([0,1,0,1])
    title("Relative permeability, "L"$k_{ro},\ k_{rw}$"" [-]")
    xlabel(L"$S_w$"" [-]")
    ylabel(L"$k_r$"" [-]")
    legend([L"$k_{ro}$", L"$k_{rw}$"])
    savefig("2D Results/rel_perm_end.png")

    figure()
    visualizeCells(p_init/1e6) #MPa
    colorbar()
    title("Pressure, "L"$p$"" [MPa]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("2D Results/p_init.png")

    figure()
    visualizeCells(sw_init)
    colorbar()
    title("Water saturation, "L"$S_w$"" [-]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("2D Results/sw_init.png")

    figure()
    visualizeCells(c_oil_init)
    colorbar()
    title("Concentration of DME in oil, "L"$c_o$"" "L"\left[\frac{mol}{m^3}\right]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("2D Results/c_oil_init.png")

    figure()
    visualizeCells(c_water_init)
    colorbar()
    title("Concentration of DME in water, "L"$c_w$"" "L"\left[\frac{mol}{m^3}\right]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("2D Results/c_water_init.png")

    figure()
    visualizeCells(c_t_init)
    colorbar()
    title("Tracer")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("2D Results/ct_init.png")

    figure()
    plot(t_s/(60*60*24), rec_fact*100)
    ylim([0, 100])
    title("Recovery factor")
    xlabel(L"$t$"" [days]")
    ylabel(L"$RF$"" [%]")
    savefig("2D Results/recovery_factor.png")
end

function plot_DME_3D(k::CellValue, KRO::Function, KRW::Function, BC_DME::BoundaryConditionDME, c_water_init::CellValue, c_oil_init::CellValue, p_init::CellValue, sw_init::CellValue, c_t_init::CellValue, t_s, rec_fact)
    figure()
    visualizeCells(k)
    colorbar()
    title("Permeability, "L"$k$"" "L"[m^2]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("3D Results/perm_field.png")

    figure()
    sw_plot = collect(linspace(0,1,1000))
    plot(sw_plot, KRO.(sw_plot, sor_c(0)), sw_plot, KRW.(sw_plot, sor_c(0)))
    plot([swc, swc], [0, 1], color=:black)
    plot([1-sor_c(0),1-sor_c(0)], [0, 1], color=:black)
    axis([0,1,0,1])
    title("Relative permeability, "L"$k_{ro},\ k_{rw}$"" [-]")
    xlabel(L"$S_w$"" [-]")
    ylabel(L"$k_r$"" [-]")
    legend([L"$k_{ro}$", L"$k_{rw}$"])
    savefig("3D Results/rel_perm_start.png")

    figure()
    sw_plot = collect(linspace(0,1,1000))
    plot(sw_plot, KRO.(sw_plot, sor_c(maximum(BC_DME.BCc_oil.left.c[:]))), sw_plot, KRW.(sw_plot, sor_c(maximum(BC_DME.BCc_oil.left.c[:]))))
    plot([swc, swc], [0, 1], color=:black)
    plot([1-sor_c(maximum(BC_DME.BCc_oil.left.c[:])),1-sor_c(maximum(BC_DME.BCc_oil.left.c[:]))], [0, 1], color=:black)
    axis([0,1,0,1])
    title("Relative permeability, "L"$k_{ro},\ k_{rw}$"" [-]")
    xlabel(L"$S_w$"" [-]")
    ylabel(L"$k_r$"" [-]")
    legend([L"$k_{ro}$", L"$k_{rw}$"])
    savefig("3D Results/rel_perm_end.png")

    figure()
    visualizeCells(p_init/1e6) #MPa
    colorbar()
    title("Pressure, "L"$p$"" [MPa]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("3D Results/p_init.png")

    figure()
    visualizeCells(sw_init)
    colorbar()
    title("Water saturation, "L"$S_w$"" [-]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("3D Results/sw_init.png")

    figure()
    visualizeCells(c_oil_init)
    colorbar()
    title("Concentration of DME in oil, "L"$c_o$"" "L"\left[\frac{mol}{m^3}\right]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("3D Results/c_oil_init.png")

    figure()
    visualizeCells(c_water_init)
    colorbar()
    title("Concentration of DME in water, "L"$c_w$"" "L"\left[\frac{mol}{m^3}\right]")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("3D Results/c_water_init.png")

    figure()
    visualizeCells(c_t_init)
    colorbar()
    title("Tracer")
    xlabel(L"$L_x$"" [m]")
    ylabel(L"$L_y$"" [m]")
    savefig("3D Results/ct_init.png")

    figure()
    plot(t_s/(60*60*24), rec_fact*100)
    ylim([0, 100])
    title("Recovery factor")
    xlabel(L"$t$"" [days]")
    ylabel(L"$RF$"" [%]")
    savefig("3D Results/recovery_factor.png")
end
