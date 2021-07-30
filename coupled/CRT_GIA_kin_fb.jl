function CRT_GIA_kin_fb(m, Nx, ϕ, u, D, k_f, k_b, nu_A, nu_B, nu_C, a, b, c, C_bc, C_0, N_comp, dt, T, tol)
  BC = Array{Any}(N_comp)
  for i in 1:N_comp
      BC[i] = createBC(m)
      # Dirichlet BC on the left side
      BC[i].left.a[:] = 0.0
      BC[i].left.b[:] = 1.0
      BC[i].left.c[:] = C_bc[i]
  end

  # initial condition
  C_old = Array{Any}(N_comp)
  for i in 1:N_comp
      C_old[i] = createCellVariable(m, C_0[i], BC[i]) # initial value of concentration
  end

  # Advective term
  u_vec = createFaceVariable(m, u)
  M_adv = convectionUpwindTerm(u_vec) # matrix of coef for the advection term

  # Diffusive term
  M_diff = Array{Any}(N_comp)
  for i in 1:N_comp
    D_face = createFaceVariable(m, D[i]*ϕ)
    M_diff[i] = diffusionTerm(D_face)      # matrix of coef for the diffusion term
  end

  # Boundary terms
  M_bc = Array{Any}(N_comp)
  RHS_bc  = Array{Any}(N_comp)
  for i in 1:N_comp
      M_bc[i], RHS_bc[i] = boundaryConditionTerm(BC[i]) # matrix and RHS for the BC
  end

  # Constructing stationary part of each system matrix
  M_stat = Array{Any}(N_comp)
  for i in 1:N_comp
    M_stat[i] = M_adv - M_diff[i] + M_bc[i]
  end

  # Preparing transient terms
  M_t = Array{Any}(N_comp)
  RHS_t = Array{Any}(N_comp)

  # Preparing Jacobian terms
  J_eps = Array{Any}(N_comp)
  RHS_eps = Array{Any}(N_comp)

  # Preparing system matrices and RHS
  M_system = Array{Any}(N_comp)
  RHS_system = Array{Any}(N_comp)
  C_new = Array{Any}(N_comp)
  C_temp = Array{Any}(N_comp)

  # Defining time-stepping parameters
  maxit = 5000 # maximum number of iterations

  # Time-stepping
  for t in dt:dt:T
    iter = 0
    error = 1000 # forcing at least one Newton-iteration

    #Assigning C_old as the new start guess for the Newton-Iteration and creating the transport related system matrix
    for i in 1:N_comp
      C_temp[i] = copyCell(C_old[i])
      M_t[i], RHS_t[i] = transientTerm(C_old[i], dt, ϕ)
      M_system[i] = M_t[i] + M_stat[i]
    end

    #Newton iterations loop
    while (error > tol && iter < maxit)

      # Constructing Jacobian terms
      J_eps[1] = linearSourceTerm(createCellVariable(m, a*k_f.*(internalCells(C_temp[1]).^(a-1)).*(internalCells(C_temp[2]).^b), BC[1]))
      J_eps[2] = linearSourceTerm(createCellVariable(m, b*k_f.*(internalCells(C_temp[1]).^a).*(internalCells(C_temp[2]).^(b-1)), BC[2]))
      J_eps[3] = linearSourceTerm(createCellVariable(m, -c*k_b.*(internalCells(C_temp[3]).^(c-1)), BC[3]))

      RHS_reacA = constantSourceTerm(createCellVariable(m, -k_f.*(internalCells(C_temp[1]).^a).*(internalCells(C_temp[2]).^b).*(1-a-b)+k_b*(internalCells(C_temp[3]).^c)*(1-c), BC[1]))
      RHS_reacB = constantSourceTerm(createCellVariable(m, -k_f.*(internalCells(C_temp[1]).^a).*(internalCells(C_temp[2]).^b).*(1-a-b)+k_b*(internalCells(C_temp[3]).^c)*(1-c), BC[2]))
      RHS_reacC = constantSourceTerm(createCellVariable(m, -k_f.*(internalCells(C_temp[1]).^a).*(internalCells(C_temp[2]).^b).*(1-a-b)+k_b*(internalCells(C_temp[3]).^c)*(1-c), BC[3]))

      # Constructing the RHS's
      RHS_system[1] = RHS_t[1] + RHS_bc[1] - nu_A*RHS_reacA
      RHS_system[2] = RHS_t[2] + RHS_bc[2] - nu_B*RHS_reacB
      RHS_system[3] = RHS_t[3] + RHS_bc[3] + nu_C*RHS_reacC

      RHS = [RHS_system[1];
            RHS_system[2];
            RHS_system[3]]

      # Constructing block system matrix
      M = [M_system[1]-(nu_A*J_eps[1]) -(nu_A*J_eps[2]) -(nu_A*J_eps[3]);
          -(nu_B*J_eps[1]) M_system[2]-(nu_B*J_eps[2]) -(nu_B*J_eps[3]);
          (nu_C*J_eps[1]) (nu_C*J_eps[2]) M_system[3]+(nu_C*J_eps[3])]

      # Solving the system
      sol = M\RHS

      # Putting solution into JFVM structure
      for i in 1:N_comp
        C_new[i] = copyCell(C_temp[i])
        C_new[i].value[:] = sol[((i-1)*(Nx+2)+1):(i*(Nx+2))]
      end

      # Calculating convergence criteria
      error = 0
      for i in 1:N_comp
        error = error + sum(abs(C_new[i].value-C_temp[i].value))
        C_temp[i] = copyCell(C_new[i])
      end

      iter = iter + 1
    end

    #if (t == dt) break end

    for i in 1:N_comp
      C_old[i] = copyCell(C_new[i])
    end

    #print(iter, ", ", t, ", ", error, "\n")
    #print(sum(internalCells(C_old[3])-K_eq*internalCells(C_old[1]).*internalCells(C_old[2])), "\n") #mass balance
  end

  return C_old
end
