include("RateTranscription.jl")
include("RateTranslation.jl")

#Setup Time Vector
tStart = 0.0
tStep = 0.01
tStop = 5.0
tSim = collect(tStart:tStep:tStop)

#Given Parameters
copies = 200
Lx1 = 1200
Lx2 = 2400
Lx3 = 600

#Setup initial conditions
x0 = [0.0; #m1
      0.0; #m2
      0.0; #m3
      0.0; #p1
      0.0; #p2
      0.0; #p3
      ] #I
Inducer = 10

#Stoichiometric matrix parameters
rx1 = RatemRNA(Lx1, copies)
rx2 = RatemRNA(Lx2, copies)  #M/s per gDCW
rx3 = RatemRNA(Lx3, copies)
rL1 = RateAA(Lx1, x0[1])
rL2 = RateAA(Lx2, x0[2])
rL3 = RateAA(Lx3, x0[3])
#Binding Constants for RNAP/transcription
n = 1.5
Kd1 = 1
Kd = 1e-5
Kx = 0
Ki = 1
K12 = 1
K32 = 1
K13 = 1
K23 = 1
#Binding Constants for Ribosome/translation
TL1 = 1
TL2 = 1
TL3 = 1
KL1 = 1
KL2 = 1
KL3 = 1
#First Order rxn parameters
k_xd = 1;   #mRNA degradation 1/s
k_ld = 1;   #protein degradation 1/s
mu_xd = log(2) / 0.5;  #mRNA dilution
mu_ld = log(2) / 0.5;  #protein dilution
#Setup fraction bound equations
fi = (Inducer)^n / ( (Kd)^n + (Inducer)^n )
f12 = (x0[4])^n / ( (Kd)^n + (x0[4])^n )
f32 = (x0[6])^n / ( (Kd)^n + (x0[6])^n )
f13 = (x0[4])^n / ( (Kd)^n + (x0[4])^n )
f23 = (x0[5])^n / ( (Kd)^n + (x0[5])^n )
#Setup Control Functions
ux1 = (Kx + Ki*fi) / (1 + Kx + Ki*fi)
ux2 = (Kx + K12*f12 + K32*f32) / (1 + Kx + K12*f12 + K32*f32)
ux3 = (Kx + K13*f13 + K23*f23) / (1 + Kx + K13*f13 + K23*f23)
uL1 = 1 #(x0[4]) / (TL1*KL1 + x0[4])
uL2 = 1 #(x0[5]) / (TL2*KL2 + x0[5])
uL3 = 1 #(x0[6]) / (TL3*KL3 + x0[6])

#Setup initial conditions for r
r = [rx1*ux1;
      rx2*ux2;
      rx3*ux3;
      rL1*uL1;
      rL2*uL2;
      rL3*uL3;
      ]

#Setup the species and stoichmetric matrices
using LinearAlgebra
A = [-k_xd 0.0 0.0 0.0 0.0 0.0;         #m1
      0.0 -k_xd-mu_xd 0.0 0.0 0.0 0.0;  #m2
      0.0 0.0 -k_xd-mu_xd 0.0 0.0 0.0;  #m3
      0.0 0.0 0.0 -k_ld 0.0 0.0;        #p1
      0.0 0.0 0.0 0.0 -k_ld-mu_ld 0.0;  #p2
      0.0 0.0 0.0 0.0 0.0 -k_ld-mu_ld;  #p3
      ]
S = Matrix{Float64}(I, 6, 6)
A_hat = exp(A*tStep)
Iden = Matrix{Float64}(I, 6, 6)
S_hat = inv(A) * (A_hat - Iden) * S

xk = Vector(undef, length(tSim))
m1 = Vector{Float64}(undef, length(xk))
m2 = Vector{Float64}(undef, length(xk))
m3 = Vector{Float64}(undef, length(xk))
p1 = Vector{Float64}(undef, length(xk))
p2 = Vector{Float64}(undef, length(xk))
p3 = Vector{Float64}(undef, length(xk))

xk[1] = x0
m1[1] = x0[1]
m2[1] = x0[2]
m3[1] = x0[3]
p1[1] = x0[4]
p2[1] = x0[5]
p3[1] = x0[6]

i = 1;
while i < (length(tSim)-1)
      if i < 100
            Inducer = 10
      else
            Inducer = 0
      end
      #Do the matrix calculations
      xk[i+1] = A_hat*xk[i] + S_hat*r
      #Update the amount
      m1[i+1] = xk[i+1][1]
      m2[i+1] = xk[i+1][2]
      m3[i+1] = xk[i+1][3]
      p1[i+1] = xk[i+1][4]
      p2[i+1] = xk[i+1][5]
      p3[i+1] = xk[i+1][6]
      #Update the translation rate based on mRNA conc
      r[4] = RateAA(Lx1, m1[i+1])
      r[5] = RateAA(Lx2, m2[i+1])
      r[6] = RateAA(Lx3, m3[i+1])
      #Setup fraction bound equations
      fi = (Inducer)^n / ( (Kd)^n + (Inducer)^n )
      f12 = (p1[i+1])^n / ( (Kd)^n + (p1[i+1])^n )
      f32 = (p3[i+1])^n / ( (Kd)^n + (p3[i+1])^n )
      f13 = (p1[i+1])^n / ( (Kd)^n + (p1[i+1])^n )
      f23 = (p2[i+1])^n / ( (Kd)^n + (p2[i+1])^n )
      #Setup Control Functions
      ux1 = (Kx + Ki*fi) / (1 + Kx + Ki*fi)
      ux2 = (Kx + K12*f12 + K32*f32) / (1 + Kx + K12*f12 + K32*f32)
      ux3 = (Kx + K13*f13 + K23*f23) / (1 + Kx + K13*f13 + K23*f23)
      #Update the transcription rate based on protein conc
      r[1] = rx1*ux1;
      r[2] = rx2*ux2;
      r[3] = rx3*ux3;
      global i += 1
end

using PyPlot
figure(1)
plot(tSim[1:end-1],m1[1:end-1],color="black")
plot(tSim[1:end-1],m2[1:end-1],color="blue")
plot(tSim[1:end-1],m3[1:end-1],color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 5, 0, 3e-7])
tight_layout()

figure(2)
plot(tSim[1:end-1],p1[1:end-1],color="black")
plot(tSim[1:end-1],p2[1:end-1],color="blue")
plot(tSim[1:end-1],p3[1:end-1],color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 5, 0, 3e-4])
tight_layout()

#----------------------------------------------------------#
#Simulation 2 with Inducer concentration zero

#Setup Time Vector
#=
tStart = 1.0
tStep = 0.01
tStop = 5
tSim = collect(tStart:tStep:tStop)

#Given Parameters
copies = 200
Lx1 = 1200
Lx2 = 2400
Lx3 = 600

#Setup initial conditions
x0 = [m1[end-1]; #m1
      m2[end-1]; #m2
      m3[end-1]; #m3
      p1[end-1]; #p1
      p2[end-1]; #p2
      p3[end-1]; #p3
      ] #I
Inducer = 0

#Stoichiometric matrix parameters
rx1 = RatemRNA(Lx1, copies)
rx2 = RatemRNA(Lx2, copies)  #M/s per gDCW
rx3 = RatemRNA(Lx3, copies)
rL1 = RateAA(Lx1, x0[1])
rL2 = RateAA(Lx2, x0[2])
rL3 = RateAA(Lx3, x0[3])

#Setup fraction bound equations
fi = (Inducer)^n / ( (Kd)^n + (Inducer)^n )
f12 = 1
f32 = 1
f13 = 1
f23 = 1
#Setup Control Functions
ux1 = (Kx + Ki*fi) / (1 + Kx + Ki*fi)
ux2 = (Kx + K12*f12 + K32*f32) / (1 + Kx + K12*f12 + K32*f32)
ux3 = (Kx + K13*f13 + K23*f23) / (1 + Kx + K13*f13 + K23*f23)
uL1 = (x0[4]) / (TL1*KL1 + x0[4])
uL2 = (x0[5]) / (TL2*KL2 + x0[5])
uL3 = (x0[6]) / (TL3*KL3 + x0[6])
#First Order rxn parameters
k_xd = 1;   #mRNA degradation 1/s
k_ld = 1;   #protein degradation 1/s
mu_xd = log(2) / 0.5;  #mRNA dilution
mu_ld = log(2) / 0.5;  #protein dilution

#Setup initial conditions for r
r = [rx1*ux1;
      rx2*ux2;
      rx3*ux3;
      rL1*uL1;
      rL2*uL2;
      rL3*uL3;
      ]

#Setup the species and stoichmetric matrices
using LinearAlgebra
A = [-k_xd 0.0 0.0 0.0 0.0 0.0;         #m1
      0.0 -k_xd-mu_xd 0.0 0.0 0.0 0.0;  #m2
      0.0 0.0 -k_xd-mu_xd 0.0 0.0 0.0;  #m3
      0.0 0.0 0.0 -k_ld 0.0 0.0;        #p1
      0.0 0.0 0.0 0.0 -k_ld-mu_ld 0.0;  #p2
      0.0 0.0 0.0 0.0 0.0 -k_ld-mu_ld;  #p3
      ]
S = Matrix{Float64}(I, 6, 6)
A_hat = exp(A*tStep)
Iden = Matrix{Float64}(I, 6, 6)
S_hat = inv(A) * (A_hat - Iden) * S

xk = Vector(undef, length(tSim))
m1 = Vector{Float64}(undef, length(xk))
m2 = Vector{Float64}(undef, length(xk))
m3 = Vector{Float64}(undef, length(xk))
p1 = Vector{Float64}(undef, length(xk))
p2 = Vector{Float64}(undef, length(xk))
p3 = Vector{Float64}(undef, length(xk))

xk[1] = x0
m1[1] = x0[1]
m2[1] = x0[2]
m3[1] = x0[3]
p1[1] = x0[4]
p2[1] = x0[5]
p3[1] = x0[6]

i = 1;
while i < (length(tSim)-1)
      #Do the matrix calculations
      xk[i+1] = A_hat*xk[i] + S_hat*r
      #Update the amount
      m1[i+1] = xk[i+1][1]
      m2[i+1] = xk[i+1][2]
      m3[i+1] = xk[i+1][3]
      p1[i+1] = xk[i+1][4]
      p2[i+1] = xk[i+1][5]
      p3[i+1] = xk[i+1][6]
      #Update the translation rate based on mRNA conc
      r[4] = RateAA(Lx1, xk[i+1][1])
      r[5] = RateAA(Lx2, xk[i+1][2])
      r[6] = RateAA(Lx3, xk[i+1][3])
      global i += 1
end

using PyPlot
figure(1)
plot(tSim[1:end-1],m1[1:end-1],color="black")
plot(tSim[1:end-1],m2[1:end-1],color="blue")
plot(tSim[1:end-1],m3[1:end-1],color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 5, 0, 3e-7])
tight_layout()

figure(2)
plot(tSim[1:end-1],p1[1:end-1],color="black")
plot(tSim[1:end-1],p2[1:end-1],color="blue")
plot(tSim[1:end-1],p3[1:end-1],color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 5, 0, 2e-4])
tight_layout()

=#
