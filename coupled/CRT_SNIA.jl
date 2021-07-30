function CRT_SNIA(m, Nx, ϕ, u, D, K_eq, nu_A, nu_B, nu_C, C_bc, C_0, N_comp, dt, T)
  # boundary conditions
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

  # Time-stepping
  for t in dt:dt:T
      # transport equation
      for i in 1:N_comp
          M_t[i], RHS_t[i] = transientTerm(C_old[i], dt, ϕ)
          M = M_t[i]+M_stat[i]
          RHS = RHS_t[i]+RHS_bc[i]
          C_new = solveLinearPDE(m, M, RHS)
          C_old[i] = copyCell(C_new)
      end

      # reaction equation
      eps1 = epsilon.(internalCells(C_old[1]), internalCells(C_old[2]),
              internalCells(C_old[3]), nu_A, nu_B, nu_C, K_eq)
      C_old[1] = createCellVariable(m, internalCells(C_old[1])-nu_A*eps1, BC[1])
      C_old[2] = createCellVariable(m, internalCells(C_old[2])-nu_B*eps1, BC[2])
      C_old[3] = createCellVariable(m, internalCells(C_old[3])+nu_C*eps1, BC[3])

      #print(sum(internalCells(C_old[3])-K_eq*internalCells(C_old[1]).*internalCells(C_old[2])), "\n") #mass balance
  end

  return C_old
end
