include("Balances.jl")
using ODE

#Setup Time Vector
tStart = 0.0
tStep = 10
tStop = 10000
tSim = collect(tStart:tStep:tStop)

#Setup initial conditions
x0 = [0.0; #E0
      0.0; #E1
      10.0; #E1*
      0.0; #B
      2.0; #Bp
      0.0; #E1*B
      0.0; #E1*Bp
      ]

f(t,x) = Balances(t,x)
t,X = ode23s(f,x0,tSim; points=:specified)

E0 = [i[1] for i in X]
E1 = [i[2] for i in X]
E1_star = [i[3] for i in X]
B = [i[4] for i in X]
Bp = [i[5] for i in X]
E1_starB = [i[6] for i in X]
E1_starBp = [i[7] for i in X]

#Protein concetration
Etot = 10
R = 0.2
#Constants
kr = 0.1
kbp = 1.0
dbp = 0.01
abp = 100.0
#Calculating parameters for activity
Kb = (kbp + dbp) / abp
Vr_max = kr*R
Ast = similar(Bp)
y = similar(Ast)

#Calculate activity
for i = 1:length(Bp)
      Vb_max = kbp*Bp[i]
      Ast[i] = Kb * (Vr_max / (Vr_max + Vb_max))
      #y[i] = E1_star[i] / Ast[i]
end

using PyPlot
figure(1)
plot(t,Ast,color="black")
xlabel("Time (s)")
ylabel("Activity")
