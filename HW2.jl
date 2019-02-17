include("Balances.jl")
using Pkg
Pkg.add("ODE")
using ODE
Pkg.add("PyPlot")


#Setup Time Vector
tStart = 0.0
tStep = 0.1
tStop = 10.0
tSim = collect(tStart:tStep:tStop)

#Setup initial conditions
x0 = [0.0; #m1
      0.0; #m2
      0.0; #m3
      0.0; #p1
      0.0; #p2
      0.0; #p3
      ]

f(t,x) = Balances(t,x)
t,X = ode23s(f,x0,tSim; points=:specified)

m1 = [i[1] for i in X]
m2 = [i[2] for i in X]
m3 = [i[3] for i in X]
p1 = [i[4] for i in X]
p2 = [i[5] for i in X]
p3 = [i[6] for i in X]

using PyPlot
#figure(figsize=(4,3))
figure(1)
plot(t,m1,color="black")
plot(t,m2,color="blue")
plot(t,m2,color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 10, 0, 1])
tight_layout()
