using Polynomials, Roots, JLD, PyPlot
using JFVM
include("peteng-master/functions/rel_perms_real.jl")

pce = 1e6
swc = 0.12
sor = 0.08

figure()
sw_plot = collect(linspace(0,1,1000))
plot(sw_plot, pc_drain.(sw_plot, pce, swc), sw_plot, pc_imb.(sw_plot, pce, swc, sor))
title("Capillary Pressure, "L"$p_c$"" [Pa]")
xlabel(L"$S_w$"" [-]")
ylabel(L"$p_c$"" [Pa]")
legend(["Drainage", "Imbibition"])
savefig("Capillary Pressure/method1.png")


figure()
sw_plot = collect(linspace(0,1,1000))
plot(sw_plot, pc_drain2.(sw_plot, pce, swc), sw_plot, pc_imb2.(sw_plot, pce, swc, sor))
title("Capillary Pressure, "L"$p_c$"" [Pa]")
xlabel(L"$S_w$"" [-]")
ylabel(L"$p_c$"" [Pa]")
legend(["Drainage", "Imbibition"])
savefig("Capillary Pressure/method2.png")

figure()
sw_plot = collect(linspace(0,1,1000))
plot(sw_plot, pc_drain3.(sw_plot, pce, swc), sw_plot, pc_imb3.(sw_plot, pce, swc, sor))
title("Capillary Pressure, "L"$p_c$"" [Pa]")
xlabel(L"$S_w$"" [-]")
ylabel(L"$p_c$"" [Pa]")
legend(["Drainage", "Imbibition"])
savefig("Capillary Pressure/method3.png")

#Van-Genuchten
a_drain = 0.066
a_imb = 0.6
m_drain = 0.6
m_imb = 0.5
M_drain = 0.26
M_imb = 0.2

k = 1e-12
ϕ = 0.4
sigma = 25.0*1000 #N/m
theta = 0.785

J = (S, S_wc, S_or, a, m, M) -> a.*(((S-S_wc)./(1.0-S_or-S_wc)).^(-1.0./m)-1.0).^(1.0-m)-M
p_c = (S, S_wc, S_or, a, m, M) -> J.(S, S_wc, S_or, a, m, M).*sigma.*cos(theta).*sqrt(ϕ/k)

sw = collect(linspace(swc+0.02,1.0-sor,1000))

figure()
plot(sw, p_c.(sw, swc, sor, a_drain, m_drain, M_drain), sw, p_c.(sw, swc, sor, a_imb, m_imb, M_imb))
title("Capillary Pressure, "L"$p_c$"" [Pa]")
xlabel(L"$S_w$"" [-]")
ylabel(L"$p_c$"" [Pa]")
legend(["Drainage", "Imbibition"])
savefig("Capillary Pressure/Van_Genuchten.png")
