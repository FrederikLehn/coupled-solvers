module CRT

using JFVM
using PyPlot

export CRT_SNIA, CRT_SIA, CRT_GIA_eps, CRT_GIA_kin_eq, CRT_GIA_kin_fb, plot_concentrations

#Define functions
function epsilon(C_A, C_B, C_C, nu_A, nu_B, nu_C, K_eq)
    # A + B <-> C
    return 1/(2*nu_A*nu_B*K_eq)*(nu_A*C_B*K_eq+nu_B*C_A*K_eq+nu_C-sqrt(nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2))
end

function epsilon_A(C_A, C_B, C_C, nu_A, nu_B, nu_C, K_eq)
  return 1/(2*nu_A*nu_B*K_eq)*(nu_B*K_eq-(C_A*K_eq^2*nu_B^2-C_B*K_eq^2*nu_A*nu_B+K_eq*nu_B*nu_C)/sqrt(nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2))
end

function epsilon_AA(C_A, C_B, C_C, nu_A, nu_B, nu_C, K_eq)
  return 1/(2*nu_A*nu_B*K_eq)*(-(K_eq^2*nu_B^2)/sqrt(nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2)+((C_A*K_eq^2*nu_B^2-C_B*K_eq^2*nu_A*nu_B+K_eq*nu_B*nu_C)*(2*C_A*K_eq^2*nu_B^2-2*C_B*K_eq^2*nu_A*nu_B+2*K_eq*nu_B*nu_C))/(2*(nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2)^(3/2)))
end

function epsilon_B(C_A, C_B, C_C, nu_A, nu_B, nu_C, K_eq)
  return 1/(2*nu_A*nu_B*K_eq)*(nu_A*K_eq-(C_B*K_eq^2*nu_A^2-C_A*K_eq^2*nu_A*nu_B+K_eq*nu_A*nu_C)/sqrt(nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2))
end

function epsilon_BB(C_A, C_B, C_C, nu_A, nu_B, nu_C, K_eq)
  return 1/(2*nu_A*nu_B*K_eq)*(-(K_eq^2*nu_A^2)/sqrt(nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2)+((-C_A*K_eq^2*nu_A*nu_B+C_B*K_eq^2*nu_A+K_eq*nu_A*nu_C)*(-2*C_A*K_eq^2*nu_A*nu_B+2*C_B*K_eq^2*nu_A^2+2*K_eq*nu_A*nu_C))/(2*(nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2)^(3/2)))
end

function epsilon_C(C_A, C_B, C_C, nu_A, nu_B, nu_C, K_eq)
  return -1/sqrt(nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2)
end

function epsilon_CC(C_A, C_B, C_C, nu_A, nu_B, nu_C, K_eq)
  return (2*K_eq*nu_A*nu_B)/((nu_A^2*C_B^2*K_eq^2-2*nu_A*nu_B*C_A*C_B*K_eq^2+nu_B^2*C_A^2*K_eq^2+4*nu_A*nu_B*C_C*K_eq+2*nu_A*nu_C*C_B*K_eq+2*nu_B*nu_C*C_A*K_eq+nu_C^2)^(3/2))
end

include("CRT_SNIA.jl")
include("CRT_SIA.jl")
include("CRT_GIA_eps.jl")
include("CRT_GIA_kin_eq.jl")
include("CRT_GIA_kin_fb.jl")
include("plot_concentrations.jl")



end
