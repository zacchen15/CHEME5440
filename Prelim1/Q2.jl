include("RateTranscription.jl")
include("RateTranslation.jl")
include("BindingFunction.jl")

#------------------------------------------------------------------------------#
#-----------------------------Question 2a--------------------------------------#

# Saturation constants
transcription_saturation_constant = 0.24 # nmol/gDW
Kx = transcription_saturation_constant
translation_saturation_constant = 465.64 # nmol/gDW
Kl = translation_saturation_constant

# Setup Time Vector (this runs in hours)
tStart = 0.0
tStep = 1/60
tStop = 6
tSim = collect(tStart:tStep:tStop)

# Gene and protein lengths
copies = 200
Lx = 1000
Lx1 = 1200
Lx2 = 2400
Lx3 = 600
Lt = 333
Lt1 = Lx1 / 3
Lt2 = Lx2 / 3
Lt3 = Lx3 / 3

# Dilution rate
cell_doubling_time = 40 * (1/60)        # hr
mu_xd = log(2) / cell_doubling_time     # mRNA dilution 1/hr
mu_ld = log(2) / cell_doubling_time     # protein dilution 1/hr

# Degradation rate
mrna_half_life = 2.1 * (1/60)           # hr
protein_half_life = 24                  # hr
k_xd = log(2) / mrna_half_life          # mRNA degradation 1/hr
k_ld = log(2) / protein_half_life       # protein degradation 1/hr

# Setup initial conditions
x0 = [0.0; #m1
      0.0; #m2
      0.0; #m3
      0.0; #p1
      0.0; #p2
      0.0] #p3

inducer_concentration = 0 # mM

# Stoichiometric matrix parameters
rx1 = transcription_rate(Lx1, copies)   # nmol/gDW-hour
rx2 = transcription_rate(Lx2, copies)   # nmol/gDW-hour
rx3 = transcription_rate(Lx3, copies)   # nmol/gDW-hour
rL1 = translation_rate(Lt1, x0[1])      # nmol/gDW-hour
rL2 = translation_rate(Lt2, x0[2])      # nmol/gDW-hour
rL3 = translation_rate(Lt3, x0[3])      # nmol/gDW-hour

# Weights for setting up control function
WI1 = 100
W11 = 0.000001
W12 = 10.0
W13 = 5.0
W22 = 0.000001
W23 = 100.0
W33 = 0.000001

# Promoter binding parameters
kI1 = 0.30 # mM
nI1 = 1.5
k12 = 1000.0 # nmol/gDW
n12 = 1.5
k13 = 1000.0 # nmol/gDW
n13 = 1.5
k23 = 1000.0 # nmol/gDW
n23 = 10.0

#Setup fraction bound equations
fI1 = binding_function(kI1, inducer_concentration, nI1)
f12 = binding_function(k12, x0[4], n12)
f13 = binding_function(k13, x0[4], n13)
f23 = binding_function(k23, x0[5], n23)

#Setup control functions for transcription (u) and translation (w)
u1 = (W11 + WI1*fI1) / (1 + W11 + WI1*fI1)
u2 = (W22 + W12*f12) / (1 + W22 + W12*f12)
u3 = (W33 + W13*f13) / (1 + W33 + W13*f13 + W23*f23)

#Setup the species and stoichmetric matrices
using LinearAlgebra
A = [-k_xd-mu_xd 0.0 0.0 0.0 0.0 0.0;   # m1
      0.0 -k_xd-mu_xd 0.0 0.0 0.0 0.0;  # m2
      0.0 0.0 -k_xd-mu_xd 0.0 0.0 0.0;  # m3
      0.0 0.0 0.0 -k_ld-mu_ld 0.0 0.0;  # p1
      0.0 0.0 0.0 0.0 -k_ld-mu_ld 0.0;  # p2
      0.0 0.0 0.0 0.0 0.0 -k_ld-mu_ld]  # p3
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

# Setup initial conditions for r
r = Array{Float64}(undef, 6, length(xk))
r[:,1] = [rx1*u1;
          rx2*u2;
          rx3*u3;
          rL1;
          rL2;
          rL3]

# Record all u values
u = Array{Float64}(undef, 3, length(xk))
u[:,1] = [u1;
          u2;
          u3]

# Record all f values
f = Array{Float64}(undef, 4, length(xk))
f[:,1] = [fI1;
          f12;
          f13;
          f23]

i = 1;
while i < (length(tSim))
      # Figure out the time
      time = tSim[i]
      if time <= 1
            inducer_concentration = 0 # mM
      elseif (time > 1) & (time <= 6)
            inducer_concentration = 10
      else
            inducer_concentration = 10 # mM
      end

      # Do the matrix calculations
      xk[i+1] = A_hat*xk[i] + S_hat*r[:,i]

      # Update the concentration of mRNA and protein
      m1[i+1] = xk[i+1][1]
      m2[i+1] = xk[i+1][2]
      m3[i+1] = xk[i+1][3]
      p1[i+1] = xk[i+1][4]
      p2[i+1] = xk[i+1][5]
      p3[i+1] = xk[i+1][6]

      # Update the binding function
      f[1,i+1] = binding_function(kI1, inducer_concentration, nI1)
      f[2,i+1] = binding_function(k12, p1[i+1], n12)
      f[3,i+1] = binding_function(k13, p1[i+1], n13)
      f[4,i+1] = binding_function(k23, p2[i+1], n23)

      # Update the control function
      u[1,i+1] = (W11 + WI1*f[1,i+1]) / (1 + W11 + WI1*f[1,i+1])
      u[2,i+1] = (W22 + W12*f[2,i+1]) / (1 + W22 + W12*f[2,i+1])
      u[3,i+1] = (W33 + W13*f[3,i+1]) / (1 + W33 + W13*f[3,i+1] + W23*f[4,i+1])

      #Update the transcription rate based on protein concentration
      r[1,i+1] = rx1*u[1,i+1]
      r[2,i+1] = rx2*u[2,i+1]
      r[3,i+1] = rx3*u[3,i+1]
      r[4,i+1] = translation_rate(Lt1, m1[i+1])      # nmol/gDW-hour
      r[5,i+1] = translation_rate(Lt2, m2[i+1])
      r[6,i+1] = translation_rate(Lt3, m3[i+1])

      global i += 1
end

using PyPlot
figure(1)
plot(tSim[1:end-1],m1[1:end-1],color="black",label="m1")
plot(tSim[1:end-1],m2[1:end-1],color="blue",label="m2")
plot(tSim[1:end-1],m3[1:end-1],color="red",label="m3")
legend(loc="upper right")
xlabel("Time (hours)")
ylabel("mRNA Concentration (nmol/gDW)")
# axis([0, 15, 0, 3e-7])
tight_layout()

figure(2)
plot(tSim[1:end-1],p1[1:end-1],color="black",label="p1")
plot(tSim[1:end-1],p2[1:end-1],color="blue",label="p2")
plot(tSim[1:end-1],p3[1:end-1],color="red",label="p3")
legend(loc="upper right")
xlabel("Time (hours)")
ylabel("Protein Concentration (nmol/gDW)")
# axis([0, 15, 0, 3e-4])
tight_layout()

#------------------------------------------------------------------------------#
#-----------------------------Question 2b--------------------------------------#

# Set up sensistivty for p1, p2, and p3 for Phase I
# p1 has parameters m1, fiI, u1, r1, r4
sp11 = Array{Float64}(undef,5,20)
# p2 has parameters m2, f12, u2, r2, r5
sp21 = Array{Float64}(undef,5,20)
# p2 has parameters m2, f13, f23, u3, r3, r6
sp31 = Array{Float64}(undef,6,20)

# looking at the sensitivty in the middle 20 minutes Phase 1
time_interval1 = collect(20:1:40)
for index1 in 1:length(time_interval1)-1
      time = time_interval1[index1]
      # p1 sensetivity
      sp11[1,index1] = (m1[time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (m1[time+1] - m1[time-1]) )
      sp11[2,index1] = (f[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[1,time+1] - f[1,time-1]) )
      sp11[3,index1] = (u[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[1,time+1] - u[1,time-1]) )
      sp11[4,index1] = (r[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[1,time+1] - r[1,time-1]) )
      sp11[5,index1] = (r[4,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[4,time+1] - r[4,time-1]) )
      # p2 sensetivity
      sp21[1,index1] = (m2[time] / p2[time]) * ( (p2[time+1] - p2[time-1]) / (m2[time+1] - m2[time-1]) )
      sp21[2,index1] = (f[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[2,time+1] - f[2,time-1]) )
      sp21[3,index1] = (u[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[2,time+1] - u[2,time-1]) )
      sp21[4,index1] = (r[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[2,time+1] - r[2,time-1]) )
      sp11[5,index1] = (r[5,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[5,time+1] - r[5,time-1]) )
      # p3 sensetivity
      sp31[1,index1] = (m2[time] / p2[time]) * ( (p2[time+1] - p2[time-1]) / (m2[time+1] - m2[time-1]) )
      sp31[2,index1] = (f[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[3,time+1] - f[3,time-1]) )
      sp31[3,index1] = (f[4,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[4,time+1] - f[4,time-1]) )
      sp31[4,index1] = (u[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[3,time+1] - u[3,time-1]) )
      sp31[5,index1] = (r[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[3,time+1] - r[3,time-1]) )
      sp31[6,index1] = (r[6,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[6,time+1] - r[6,time-1]) )
end

# Set up sensistivty for p1, p2, and p3 for Early Phase II (between hour 1 and 2)
# p1 has parameters m1, fiI, u1, r1, r4
sp12 = Array{Float64}(undef,5,20)
# p2 has parameters m2, f12, u2, r2, r5
sp22 = Array{Float64}(undef,5,20)
# p2 has parameters m2, f13, f23, u3, r3, r6
sp32 = Array{Float64}(undef,6,20)

# looking at the sensitivty in the middle 20 minutes Early Phase II (between hour 1 and 2)
time_interval2 = collect(80:1:100)
for index2 in 1:length(time_interval2)-1
      time = time_interval2[index2]
      # p1 sensetivity
      sp12[1,index2] = (m1[time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (m1[time+1] - m1[time-1]) )
      sp12[2,index2] = (f[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[1,time+1] - f[1,time-1]) )
      sp12[3,index2] = (u[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[1,time+1] - u[1,time-1]) )
      sp12[4,index2] = (r[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[1,time+1] - r[1,time-1]) )
      sp12[5,index2] = (r[4,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[4,time+1] - r[4,time-1]) )
      # p2 sensetivity
      sp22[1,index2] = (m2[time] / p2[time]) * ( (p2[time+1] - p2[time-1]) / (m2[time+1] - m2[time-1]) )
      sp22[2,index2] = (f[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[2,time+1] - f[2,time-1]) )
      sp22[3,index2] = (u[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[2,time+1] - u[2,time-1]) )
      sp22[4,index2] = (r[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[2,time+1] - r[2,time-1]) )
      sp12[5,index2] = (r[5,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[5,time+1] - r[5,time-1]) )
      # p3 sensetivity
      sp32[1,index2] = (m2[time] / p2[time]) * ( (p2[time+1] - p2[time-1]) / (m2[time+1] - m2[time-1]) )
      sp32[2,index2] = (f[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[3,time+1] - f[3,time-1]) )
      sp32[3,index2] = (f[4,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[4,time+1] - f[4,time-1]) )
      sp32[4,index2] = (u[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[3,time+1] - u[3,time-1]) )
      sp32[5,index2] = (r[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[3,time+1] - r[3,time-1]) )
      sp32[6,index2] = (r[6,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[6,time+1] - r[6,time-1]) )
end

# Set up sensistivty for p1, p2, and p3 for Late Phase II (between hour 5 and 6)
# p1 has parameters m1, fiI, u1, r1, r4
sp13 = Array{Float64}(undef,5,20)
# p2 has parameters m2, f12, u2, r2, r5
sp23 = Array{Float64}(undef,5,20)
# p2 has parameters m2, f13, f23, u3, r3, r6
sp33 = Array{Float64}(undef,6,20)

# looking at the sensitivty in the middle 20 minutes Late Phase II (between hour 5 and 6)
time_interval3 = collect(320:1:340)
for index3 in 1:length(time_interval3)-1
      time = time_interval3[index3]
      # p1 sensetivity
      sp13[1,index3] = (m1[time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (m1[time+1] - m1[time-1]) )
      sp13[2,index3] = (f[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[1,time+1] - f[1,time-1]) )
      sp13[3,index3] = (u[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[1,time+1] - u[1,time-1]) )
      sp13[4,index3] = (r[1,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[1,time+1] - r[1,time-1]) )
      sp13[5,index3] = (r[4,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[4,time+1] - r[4,time-1]) )
      # p2 sensetivity
      sp23[1,index3] = (m2[time] / p2[time]) * ( (p2[time+1] - p2[time-1]) / (m2[time+1] - m2[time-1]) )
      sp23[2,index3] = (f[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[2,time+1] - f[2,time-1]) )
      sp23[3,index3] = (u[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[2,time+1] - u[2,time-1]) )
      sp23[4,index3] = (r[2,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[2,time+1] - r[2,time-1]) )
      sp13[5,index3] = (r[5,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[5,time+1] - r[5,time-1]) )
      # p3 sensetivity
      sp33[1,index3] = (m2[time] / p2[time]) * ( (p2[time+1] - p2[time-1]) / (m2[time+1] - m2[time-1]) )
      sp33[2,index3] = (f[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[3,time+1] - f[3,time-1]) )
      sp33[3,index3] = (f[4,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (f[4,time+1] - f[4,time-1]) )
      sp33[4,index3] = (u[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (u[3,time+1] - u[3,time-1]) )
      sp33[5,index3] = (r[3,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[3,time+1] - r[3,time-1]) )
      sp33[6,index3] = (r[6,time] / p1[time]) * ( (p1[time+1] - p1[time-1]) / (r[6,time+1] - r[6,time-1]) )
end

#------------------------------------------------------------------------------#
#-----------------------------Question 2c--------------------------------------#

# Create Single Value Decomposition (svd) for each sensitivity array
# Phase I
Fp11 = svd(sp11)
Fp21 = svd(sp21)
Fp31 = svd(sp31)

# Early Phase II
Fp12 = svd(sp12)
Fp22 = svd(sp22)
Fp32 = svd(sp32)

# Early Phase II
Fp13 = svd(sp31)
Fp23 = svd(sp23)
Fp33 = svd(sp33)

# Find the U matrix
# Phase I
Up11 = Fp11.U
Up21 = Fp21.U
Up31 = Fp31.U

# Early Phase II
Up12 = Fp12.U
Up22 = Fp22.U
Up32 = Fp32.U

# Early Phase II
Up13 = Fp13.U
Up23 = Fp23.U
Up33 = Fp33.U
