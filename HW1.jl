#Part A

#Find the elongation rate constant
e_x = 27 # nt/s from BioNumbers (Proshkin S., et al.)
L = 924 # average gene length in prokaryotes in bp or nt from BioNumbers (Xu L., et al)
L_j = 3075 # nt given in HW
k_Ej = (e_x/L)*(L/L_j) # 1/s

#Finding K_xj from McClure
intercept = 1.04e-6 # M-s
slope = 42 # s
K_xj = intercept/slope # Molar

#Parameters for Txj
k_I = 0.09 # 1/s at 37 degC given by BioNumbers (Buc, et al.)
k_A = 0 # assuming abortive initiation is low
T_xj = (k_Ej+k_A)/k_I # dimensionless

Vol_cell = (pi*(.308e-6)^2*(2.62e-6))/1000 # L
G_j = 2500 # copies/cell
G_j = G_j*(1/Vol_cell)*(1/6.02e23) # M per cell

R_xt = 1000 # RNAP/cell from BioNumbers
R_xt = R_xt*(1/Vol_cell)*(1/6.02e23) # M per cell

r_xj = k_Ej*R_xt*G_j*(1/((K_xj*T_xj)+(T_xj+1)*G_j)) # M/s

#-----------------------------------------------------------------------------#

# Part B
T_xj = (k_Ej+k_A)/k_I # dimensionless

#-----------------------------------------------------------------------------#

# Part C
I = [0.0001, 0.001, 0.01, 0.1, 1, 10]
n = 1.5
K = 0.30
W1 = 0.26
W2 = 300
Kd = log(2)/(4*60) # 1/s
M = log(2)/(30*60) #1/s

f_I = zeros(Float64,length(I))
u = zeros(Float64,length(I))

for i = 1:length(I)
    f_I[i] = (I[i]^n)/(K+I[i]^n)
    u[i] = (W1 + W2*f_I[i])/(1 + W1 + W2*f_I[i])
end

r_hat = r_xj * u
m_j = r_hat / (Kd + M)

using Pkg
Pkg.add("Plots")
using Plots
plot(I,m_j,linewidth=2,xaxis=:log, label="Steady-State mRNA concentration")
xlabel!("Inducer Concentration (mM)")
ylabel!("mRNA Concentration (M/cell)")
