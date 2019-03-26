# For flux balance analysis on Prelim 1 Q3
include("Flux.jl")

# Inducer concentration
I = collect(0.0001:0.0001:10) # mM
obj_func = zeros(Float64, length(I))
flux_array = zeros(15, length(I))
dual_value_array = zeros(15, length(I))
uptake_array = zeros(17, length(I))

for i = 1:length(I)
    # Given information on the gene
    n = 924 # nt
    a = 308 # aa
    G = 5e-3 # uM
    Vol = 15 #uL

    S = zeros(17,15)
    S[1,:] = [-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0] #G
    S[2,:] = [-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0] #RNAP
    S[3,:] = [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0] #G*
    S[4,:] = [0 -n 0 0 0 0 0 1 0 0 0 0 0 0 0] #NTP
    S[5,:] = [0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0] #mRNA
    S[6,:] = [0 2n 0 0 2a 2 0 0 0 0 0 0 0 0 -1] #Pi
    S[7,:] = [0 0 n 0 0 0 0 0 0 -1 0 0 0 0 0] #NMP
    S[8,:] = [0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0] #rib
    S[9,:] = [0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0] #rib*
    S[10,:] = [0 0 0 0 -1a 1 0 0 0 0 0 0 0 0 0] #AAtRNA
    S[11,:] = [0 0 0 0 -2a 0 0 0 0 0 0 0 1 0 0] #GTP
    S[12,:] = [0 0 0 0 1a -1 0 0 0 0 0 0 0 0 0] #tRNA
    S[13,:] = [0 0 0 0 2a 0 0 0 0 0 0 0 0 -1 0] #GDP
    S[14,:] = [0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0] #protein
    S[15,:] = [0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0] #AA
    S[16,:] = [0 0 0 0 0 -1 0 0 0 0 1 0 0 0 0] #ATP
    S[17,:] = [0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0] #AMP

    # Determine transciption and translation rate
    Rx = 0.15 # uM
    Rl = 1.6 # uM
    vx = 60 # nt/s
    vx = vx * 3600 # nt/hr
    vl = 16.5 # aa/s
    vl = vl * 3600 # aa/hr
    Kx = 0.3 # uM
    Kl = 57.0 # uM
    Tx = 2.7
    Tl = 0.8
    kd_x = 8.35 # 1/h
    kd_l = 0.0099 # 1/h

    # Control function constants
    W1 = 0.26
    W2 = 300.0
    K = 0.30 #mM
    n_coop = 1.5

    # Elongation rate constant
    kE_x = vx / n
    kE_l = vl / a

    # Control function
    fi = ( I[i]^n_coop ) / ( K + I[i]^n_coop )
    ui = ( W1 + W2*fi ) / ( 1 + W1 + W2*fi )

    # Rates and fluxes
    rx = kE_x*Rx*( G / (Kx*Tx + (Tx+1)*G )) # uM / hr
    v1 = rx*ui
    v2 = v1
    v3 = v1

    mRNA = v3 / kd_x
    rl = kE_l*Rl*( mRNA / (Kl*Tl + (Tl+1)*mRNA )) # uM / hr

    # Setting the bounds on each flux
    b_upper = 100000.0 # uM / hr
    default_bounds_array = [v1 v1; #v1 transcription initiation
                            v2 v2; #v2 transcription elongation
                            v3 v3; #v3 mRNA decay
                            0.0 rl; #v4 translation initiation
                            0.0 rl; #v5 translation elongation
                            0.0 Inf; #v6 tRNA charging
                            -b_upper b_upper; #b1
                            -b_upper b_upper; #b2
                            -b_upper b_upper; #b3
                            -b_upper b_upper; #b4
                            -b_upper b_upper; #b5
                            -b_upper b_upper; #b6
                            -b_upper b_upper; #b7
                            -b_upper b_upper; #b8
                            -b_upper b_upper; #b9
                            ]

    # Bounds on each ofo the species
    species_bounds = zeros(length(default_bounds_array), 2)

    # Setting coefficients of the objective function
    obj_func_array = [0.0; #v1
                      0.0; #v2
                      0.0; #v3
                      0.0; #v4
                      -1.0; #v5
                      0.0; #v6
                      0.0; #b1
                      0.0; #b2
                      0.0; #b3
                      0.0; #b4
                      0.0; #b5
                      0.0; #b6
                      0.0; #b7
                      0.0; #b8
                      0.0; #b9
                      ]

    answer = calculate_optimal_flux_distribution(S, default_bounds_array, species_bounds, obj_func_array)
    obj_func[i] = answer[1]
    flux_array[:,i] = answer[2]
    dual_value_array[:,i] = answer[3]
    uptake_array[:,i] = answer[4]
    println(i)
end

kd_l = 0.0099 # 1/h
translation_rate = -obj_func
protein = translation_rate / kd_l

using PyPlot
figure(1)
semilogx(I,protein,color="blue")
xlabel("Inducer Concentration (mM)")
ylabel("Protein (uM)")
