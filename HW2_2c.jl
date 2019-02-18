include("RateTranscription.jl")
#Setup Time Vector
tStart = 0.0
tStep = 0.01
tStop = 60.0
tSim = collect(tStart:tStep:tStop)

#Setup initial conditions
x0 = [0.0; #m1
      0.0; #m2
      0.0; #m3
      0.0; #p1
      0.0; #p2
      0.0; #p3
      10.0] #I

#Given Parameters
copies = 200
Lx1 = 1200
Lx2 = 2400
Lx3 = 600
#Stoichiometric matrix parameters
rx1 = RatemRNA(Lx1, copies)
rx2 = RatemRNA(Lx2, copies)  #M/s per gDCW
rx3 = RatemRNA(Lx3, copies)
rL1 = 1
rL2 = 1
rL3 = 1
#Binding Constants for RNAP/transcription
n = 1.5
Kd = 0.3
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
#Setup fraction bound equations
fi = (x0[7])^n / ( (Kd)^n + (x0[7])^n )
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
k_ld = 10;   #protein degradation 1/s
mu_xd = log(2) / 30*60;  #mRNA dilution
mu_ld = log(2) / 30*60;  #protein dilution

#Setup initial conditions for r
r = [rx1*ux1;
      rx2*ux2;
      rx3*ux3;
      rL1*uL1;
      rL2*uL2;
      rL3*uL3;
      0]

#Setup the species and stoichmetric matrices
using LinearAlgebra
A = [-k_xd 0.0 0.0 0.0 0.0 0.0 0.0;         #m1
      0.0 -k_xd-mu_xd 0.0 0.0 0.0 0.0 0.0;  #m2
      0.0 0.0 -k_xd-mu_xd 0.0 0.0 0.0 0.0;  #m3
      0.0 0.0 0.0 -k_ld 0.0 0.0 0.0;        #p1
      0.0 0.0 0.0 0.0 -k_ld-mu_ld 0.0 0.0;  #p2
      0.0 0.0 0.0 0.0 0.0 -k_ld-mu_ld 0.0;  #p3
      0.0 0.0 0.0 0.0 0.0 0.0 1.0]          #I
S = Matrix{Float64}(I, 7, 7)
A_hat = exp(A*tStep)
Iden = Matrix{Float64}(I, 7, 7)
S_hat = inv(A) * (A_hat - Iden) * S

xk = Vector(undef, length(tSim)+1)
xk[1] = x0
for i = 1:length(tSim)
      xk[i+1] = A_hat*xk[i] + S_hat*r
end

m1 = Vector(undef, length(xk))
m2 = Vector(undef, length(xk))
m3 = Vector(undef, length(xk))
p1 = Vector(undef, length(xk))
p2 = Vector(undef, length(xk))
p3 = Vector(undef, length(xk))

for i in length(xk)
      m1[i] = xk[i][1]
      m2[i] = xk[i][2]
      m3[i] = xk[i][3]
      p1[i] = xk[i][4]
      p2[i] = xk[i][5]
      p3[i] = xk[i][6]
end
