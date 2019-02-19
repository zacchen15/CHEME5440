include("Balances.jl")
using ODE

#Setup Time Vector
tStart = 0.0
tStep = 0.1
tStop = 1.0
tSim = collect(tStart:tStep:tStop)

#Setup initial conditions
x0 = [0.0; #m1
      0.0; #m2
      0.0; #m3
      0.0; #p1
      0.0; #p2
      0.0; #p3
      10.0; #I
      ]

f(t,x) = Balances(t,x)
t,X = ode23s(f,x0,tSim; points=:specified)

m1 = [i[1] for i in X]
m2 = [i[2] for i in X]
m3 = [i[3] for i in X]
p1 = [i[4] for i in X]
p2 = [i[5] for i in X]
p3 = [i[6] for i in X]
Inducer = [i[7] for i in X]

using PyPlot
#figure(figsize=(4,3))
figure(1)
plot(t,m1,color="black")
plot(t,m2,color="blue")
plot(t,m3,color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 5, 0, 1e-6])
tight_layout()

figure(2)
plot(t,p1,color="black")
plot(t,p2,color="blue")
plot(t,p3,color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 5, 0, 1e-3])
tight_layout()

#Setup Time Vector
tStart = 1
tStep = 0.1
tStop = 5
tSim = collect(tStart:tStep:tStop)

#Setup initial conditions
x0 = [m1[end];
      m2[end];
      m3[end];
      p1[end];
      p2[end];
      p3[end];
      0.0; #I
      ]


f(t,x) = Balances(t,x)
t,X = ode23s(f,x0,tSim; points=:specified)

m1 = [i[1] for i in X]
m2 = [i[2] for i in X]
m3 = [i[3] for i in X]
p1 = [i[4] for i in X]
p2 = [i[5] for i in X]
p3 = [i[6] for i in X]
Inducer = [i[7] for i in X]

using PyPlot
#figure(figsize=(4,3))
figure(1)
plot(t,m1,color="black")
plot(t,m2,color="blue")
plot(t,m3,color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 5, 0, 1e-6])

figure(2)
plot(t,p1,color="black")
plot(t,p2,color="blue")
plot(t,p3,color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 5, 0, 1e-3])
tight_layout()
