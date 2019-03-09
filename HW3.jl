include("Flux.jl")

#constants
mass_of_cell = 2.3e-9 #g
vol_of_cell = 3.7e-12 #l
double_time = 20 #hr-1

#note that b is per hour while kcat is per second and E is in micromole
SS_enzyme_conc = 0.01 #umol/gDW
SS_enzyme_conc = SS_enzyme_conc * 10e-3 #mmol/gDW
b_upper = 10 #mmol/gDW-hr
b_upper = b_upper * (1/3600) #mmol/gDW-s

#Michaelis constants
Km_v1_asp = 0.15 #mM
Km_v1_cit = 0.056 #mM
Km_v1_ATP = 0.051 #mM
Km_v2 = 3  #mM
Km_v3 = 1.4 #mM
Km_v4_cp = 0.13 #mM
Km_v4_orn = 0.36 #mM
Km_v5f_arg = 0.0044 #mM
Km_v5f_NADPH = 0.0003 #mM
Km_v5r = 0 #mM (assuming saturation)

#kcat values
Kcat_v1 = 203 #s-1
Kcat_v2 = 34.5 #s-1
Kcat_v3 = 249 #s-1
Kcat_v4 = 88.1 #s-1
Kcat_v5f = 13.7 #s-1
Kcat_v5r = 13.7 #s-1 (is this reverse the same?)

#metabolite concentrations (note that these values were converted from M)
# values for argininosccuinate, urea, and CP not found
# these were assumed to be >> Km so saturation constant at 1
# parameters were given a dummy value of 1
citrulline = 2.7e1 #mM
aspartate = 1.49e1 #mM
argininosuccinate = 1 #mM
fumarate = 4.85e-1 #mM
arginine = 2.18e1 #mM
urea = 1 #mM
ornithine = 4.49 #mM
carbamoyl_phosphate = 1 #mM
ATP = 2.25 #mM
NADPH = 100 #mM

#flux reaction velocity
v1 = (Kcat_v1*SS_enzyme_conc) * (aspartate / (Km_v1_asp + aspartate)) *(ATP / (Km_v1_ATP + ATP)) #* (citrulline / (Km_v1_cit + citrulline))
v2 = (Kcat_v2*SS_enzyme_conc) * 1 #(argininsuccinate / (Km_v2 + argininsuccinate))
v3 = (Kcat_v3*SS_enzyme_conc) * (arginine / (Km_v3 + arginine))
v4 = (Kcat_v4*SS_enzyme_conc) * (ornithine / (Km_v4_orn + ornithine)) * 1 #(carbamoyl_phosphate / (Km_v4_cp + carbamoyl_phosphate))
v5f = (Kcat_v5f*SS_enzyme_conc) * (arginine / (Km_v5f_arg + arginine)) * (NADPH / (Km_v5f_NADPH + NADPH))
v5r = (Kcat_v5r*SS_enzyme_conc) * (citrulline / (Km_v5r + citrulline))

#stoichiometic matrix, should be 18x20
#there are 20 fluxes for reactions or simple metabolite transport
#   and 18 metabolites
S = Array{Float64}(undef,18,20)
S[1,:] = [-1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0] #aspartate
S[2,:] = [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] #argininsuccinate
S[3,:] = [0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0] #fumarate
S[4,:] = [0 1 -1 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0] #arginine
S[5,:] = [0 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0] #urea
S[6,:] = [0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] #ornithine
S[7,:] = [0 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0] #carbamoyl_phosphate
S[8,:] = [-1 0 0 1 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0] #citrulline
S[9,:] = [-1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0] #ATP
S[10,:] = [0 0 -1 0 2 -2 0 0 0 0 0 1 0 0 0 0 0 0 0 0] #H20
S[11,:] = [0 0 0 1 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0] #Pi
S[12,:] = [1 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0] #AMP
S[13,:] = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0] #PPi
S[14,:] = [0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 -1 0 0 0 0] #NO
S[15,:] = [0 0 0 0 -2 2 0 0 0 0 0 0 0 0 0 0 1 0 0 0] #O2
S[16,:] = [0 0 0 0 -1.5 1.5 0 0 0 0 0 0 0 0 0 0 0 1 0 0] #H+
S[17,:] = [0 0 0 0 -1.5 1.5 0 0 0 0 0 0 0 0 0 0 0 0 1 0] #NADPH
S[18,:] = [0 0 0 0 1.5 -1.5 0 0 0 0 0 0 0 0 0 0 0 0 0 -1] #NADP+

A = Array{Float64}(undef,6,18)
A[1,:] = [4 10 4 6 1 5 1 6 10 0 0 10 0 0 0 0 21 21] #C
A[2,:] = [7 18 4 14 4 12 4 13 16 2 3 14 4 0 0 1 30 29] #H
A[3,:] = [1 4 0 4 2 2 1 3 5 0 0 5 0 1 0 0 7 7] #N
A[4,:] = [4 6 4 2 1 2 5 3 13 1 4 7 7 1 2 0 17 17] #O
A[5,:] = [0 0 0 0 0 0 1 0 3 0 1 1 2 0 0 0 3 3] #P
A[6,:] = [0 0 0 0 0 0 1 0 3 0 1 1 2 0 0 0 0 0] #S

default_bounds_array = [0.0 v1; #v1
                        0.0 v2; #v2
                        0.0 v3; #v3
                        0.0 v4; #v4
                        0.0 v5f; #v5f
                        0.0 v5r; #v5r
                        0.0 b_upper; #b1 carbamoyl phosphate
                        0.0 b_upper; #b2 aspartate
                        0.0 b_upper; #b3 fumarate
                        0.0 b_upper; #b4 urea
                        0.0 b_upper; #b5 ATP
                        0.0 b_upper; #b6 H2O
                        0.0 b_upper; #b7 Pi
                        0.0 b_upper; #b8 AMP
                        0.0 b_upper; #b9 PPi
                        0.0 b_upper; #b10 NO
                        0.0 b_upper; #b11 O2
                        0.0 b_upper; #b12 H+
                        0.0 b_upper; #b13 NADPH
                        0.0 b_upper; #b14 NADP+
                        ]


species_bounds = zeros(length(default_bounds_array), 2)

obj_func_array = [0.0; #v1
                  0.0; #v2
                  0.0; #v3
                  0.0; #v4
                  0.0; #v5f
                  0.0; #v5r
                  0.0; #b1 carbamoyl phosphate
                  0.0; #b2 aspartate
                  0.0; #b3 fumarate
                  -1.0; #b4 urea
                  0.0; #b5 ATP
                  0.0; #b6 H2O
                  0.0; #b7 Pi
                  0.0; #b8 AMP
                  0.0; #b9 PPi
                  0.0; #b10 NO
                  0.0; #b11 O2
                  0.0; #b12 H+
                  0.0; #b13 NADPH
                  0.0; #b14 NADP+
                  ]

#element matrix to check if the reactions in the urea cycle reconstruction are balanced
E = A*S

answer = calculate_optimal_flux_distribution(S, default_bounds_array, species_bounds, obj_func_array)
obj_func = answer[1]
flux_array = answer[2]
dual_value_array = answer[3]
uptake_array = answer[4]

print("Maximimum flux of Urea is ", -1*obj_func, " mmol/gDW-s")
