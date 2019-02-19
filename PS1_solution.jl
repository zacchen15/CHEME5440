# Script to simulate the steady-state mRNA balance for PS1 CHEME 7770 S19
# Author: J.Varner
# version: 1.0

# Include some stuff -
using PyPlot
using PyCall
using DelimitedFiles
@pyimport numpy as np

# setup some constants -
fraction_water_cell = 0.70                      # dimensionless
fraction_dry_cell = 1 - fraction_water_cell     # dimensionless
volume_of_single_cell = 9e-17                   # L/cell BIND: 114922
mass_of_single_cell = 2.8e-13                   # g/cell BIND:103904
number_of_rnapII = 4600            	            # copies/cells
fraction_RNAP_available = 0.25                  # dimensionless
number_of_ribosome = 50000         	            # copies/cells
mRNA_half_life_TF = 0.022                       # hrs
protein_half_life = 70                          # hrs
doubling_time_cell = 0.5                        # hrs
max_translation_rate = 16.5                     # aa/sec
max_transcription_rate = 48.0                   # nt/sec
transcription_initiation_time_contstant = 42    # sec
average_transcript_length = 1000   	            # nt
av_number = 6.02e23                             # number/mol
avg_gene_number = 2500                          # number of copies of a gene
sample_size = 1.0                               # ml
population_size = 1e8                           # cells/ml
lacZ_gene_length = 3075                         # nt

# What is the saturation constant, gene concentration et al?
length_factor = (average_transcript_length/lacZ_gene_length)                                    # dimensionless
rnapII_concentration = (fraction_RNAP_available*number_of_rnapII)*(1/av_number)*(1/mass_of_single_cell)*1e9               # nmol/gdw
avg_gene_concentration = avg_gene_number*(1/mass_of_single_cell)*(1/av_number)*1e9              # nmol/gdw
kcat_transcription = max_transcription_rate*(3600/average_transcript_length)                    # hr^-1
kcat_transcription_initiation = (1/transcription_initiation_time_contstant)*(3600)              # hr^-1
tau_factor = (kcat_transcription)/(kcat_transcription_initiation)                               # dimensionless
m = 0.03                                                                                        # Slope muM, McClure 1980

# Compute the mu and degradation time from
mugmax = (1/doubling_time_cell)*log(2)                                                          # hr^-1
kdT = -(1/mRNA_half_life_TF)*log(0.5)                                                           # hr^-1

# ---------------------------------------------------------------------------------------------- #
# 1c - compute the mRNA predicted by the model -
saturation_transcription = m*(volume_of_single_cell)*(1/mass_of_single_cell)*(1e9/1e6)*(1/fraction_dry_cell)    # nmol/gDW
kinetic_limit = (kcat_transcription*length_factor)*(rnapII_concentration)*((avg_gene_concentration)/(saturation_transcription*tau_factor+(1+tau_factor)*avg_gene_concentration))
delta_term = (mugmax + kdT)

# Golding parameters -
K1 = 0.26
K2 = 300.0
k_binding_function = 1.5*0.30
n = 1.5

# initialize simulated data -
number_of_simulation_points = 1000
inducer_array = collect(exp10.(range(-5,stop=2,length=number_of_simulation_points)))
simulated_mRNA_concentration = zeros(length(inducer_array),3)
for step_index = 1:number_of_simulation_points

    # compute the u function -
    inducer = inducer_array[step_index]
    f_binding = ((inducer)^n)/(k_binding_function^n+(inducer)^n)
    u_value = (K1+K2*f_binding)/(1+K1+K2*f_binding)

    # compute the mRNA level -
    mRNA_level = (kinetic_limit/delta_term)*u_value

    # cache -
    simulated_mRNA_concentration[step_index,1] = inducer
    simulated_mRNA_concentration[step_index,2] = mRNA_level
    simulated_mRNA_concentration[step_index,3] = u_value
end
# ---------------------------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------------------------- #
# 1d Make a plot - mRNA (y-axis) versus I (x-axis) -
semilogx(simulated_mRNA_concentration[:,1],simulated_mRNA_concentration[:,2],color="black",lw=2)

# label the axes -
xlabel("Inducer [mM]",fontsize=16)
ylabel("mRNA [nmol/gDW]",fontsize=16)

# dump to disk -
savefig("Solution-mRNA-PS1-CHEME-7770-S19.pdf")

# Make a plot - u (y-axis) versus I (x-axis) -
figure()
semilogx(simulated_mRNA_concentration[:,1],simulated_mRNA_concentration[:,3],color="black",lw=2)

# label the axes -
xlabel("Inducer [mM]",fontsize=16)
ylabel("u (dimensionless)",fontsize=16)
savefig("Solution-u-PS1-CHEME-7770-S19.pdf")
# ---------------------------------------------------------------------------------------------- #
