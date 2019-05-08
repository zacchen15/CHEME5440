include("Q1Balance.jl")
using ODE

#Setup Time Vector
tStart = 0.0
tStep = 1.0
tStop = 30.0
tSim = collect(tStart:tStep:tStop)

#Setup Inhibitor Vectors
iStart = 0.01
iStep = 0.01
iStop = 1000.0
iSim = collect(iStart:iStep:iStop)

#Setup initial conditions
x0 = [100.0; #A
      0.0; #B
      0.0; #C
      0.0; #Inhibitor 1
      0.0; #Inhibitor 2
      ]

#Setup array for proteinA
proteinA = Array{Float64}(undef, length(tSim), length(iSim))
time = Array{Float64}(undef, length(tSim), length(iSim))
inhibitor1 = Array{Float64}(undef, length(tSim), length(iSim))
inhibitor2 = Array{Float64}(undef, length(tSim), length(iSim))

for j in 1:length(iSim)
      x0[4] = iSim[j]
      x0[5] = iSim[j]

      f(t,x) = Balances(t,x)
      t,X = ode23s(f,x0,tSim; points=:specified)

      time[:,j] = t
      proteinA[:,j] = [i[1] for i in X]
      proteinB = [i[2] for i in X]
      proteinC = [i[3] for i in X]
      inhibitor1[:,j] = [i[4] for i in X]
      inhibitor2[:,j] = [i[5] for i in X]
end

loginhibitor1 = Array{Float64}(undef, length(tSim), length(iSim))
loginhibitor2 = Array{Float64}(undef, length(tSim), length(iSim))
for i in 1:length(tSim)
      for j in 1:length(iSim)
            loginhibitor1[i,j] = log10(inhibitor1[i,j])
            loginhibitor2[i,j] = log10(inhibitor2[i,j])
      end
end

using PyPlot
close("all")
figure(1)
plot_surface(time, loginhibitor1, proteinA, color="blue")
xlabel("Time (seconds)")
ylabel("log([I1]) (units)")
zlabel("Protein A (units)")

figure(2)
plot_surface(time, loginhibitor2, proteinA, color="red")
xlabel("Time (seconds)")
ylabel("log([I2]) (units)")
zlabel("Protein A (units)")
