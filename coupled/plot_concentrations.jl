function plot_concentrations(C_old, savepath, title)
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
xlabel(L"$L_x$")
ylabel(L"$C_{\alpha}$")
ax[:set_title](title)
ax[:legend]()

savefig(savepath)
end
