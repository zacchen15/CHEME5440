# given data
# col 1: AMP concentration (mM)
# col 2: overall rate (uM/h)
# col 3: 95% confidence
data = [0.000 3.003 0.59;
        0.055 6.302 1.20;
        0.093 29.761 5.7;
        0.181 52.002 10.2;
        0.405 60.306 11.8;
        0.990 68.653 13.3]

# substrate concentration
# note that PFK is in uM
ATP_concentration = 2.3 #mM
F6P_concentration = 0.1 #mM
PFK_concentration = 0.12 #uM

# rate constants associated with each compound
ATP_binding_constant = 0.42 #mM
F6P_binding_constant = 0.11 #mM
k_cat = 0.4*3600 #1/h

# concentration range to be studied
AMPstart = 0.0
AMPstep = 0.01
AMPstop = 1.0
AMP_concentration = collect(AMPstart:AMPstep:AMPstop)

# things that need estimating
AMP_binding_constant = 1.0
n = 1.0
W1 = 0.01
W2 = 15.0

# calculate things that are constant
# rates and kinetics
saturation_ATP = (ATP_concentration) / (ATP_binding_constant + ATP_concentration)
saturation_F6P = (F6P_concentration) / (F6P_binding_constant + F6P_concentration)
kinetic_rate_limit = k_cat * PFK_concentration * saturation_ATP * saturation_F6P

u = similar(AMP_concentration)
f_AMP = similar(AMP_concentration)
for i in 1:length(AMP_concentration)
        f_AMP[i] = (AMP_concentration[i]/AMP_binding_constant)^n / (1.0 + (AMP_concentration[i]/AMP_binding_constant)^n)
        u[i] = (W1 + W2*f_AMP[i]) / (1.0 + W1 + W2*f_AMP[i])
end

overall_rate = kinetic_rate_limit * u

# plot stuff
using PyPlot
close("all")
figure(1)
errorbar(data[:,1],data[:,2],data[:,3],color="black",marker=".")
plot(AMP_concentration,overall_rate,color="red")
xlabel("AMP Concentration (mM)")
ylabel("Overall Rate (uM/h)")
