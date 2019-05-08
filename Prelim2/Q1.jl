include("Q1Balance.jl")
using ODE

#Setup Time Vector
tStart = 0.0
tStep = 1.0
tStop = 1000.0
tSim = collect(tStart:tStep:tStop)

#Setup initial conditions
x0 = [0.0; #A
      50.0; #B
      50.0; #C
      8.0; #Inhibitor 1
      9.0; #Inhibitor 2
      ]

f(t,x) = Balances(t,x)
t,X = ode23s(f,x0,tSim; points=:specified)

proteinA = [i[1] for i in X]
proteinB = [i[2] for i in X]
proteinC = [i[3] for i in X]
inhibitor1 = [i[4] for i in X]
inhibitor2 = [i[5] for i in X]

using PyPlot
close("all")
figure(1)
plot(t,proteinA,color="black",linestyle="solid")
plot(t,proteinB,color="blue",linestyle="dashed")
plot(t,proteinC,color="red",linestyle="dotted")
xlabel("Time (seconds)")
ylabel("Concentration (units)")
tight_layout()
