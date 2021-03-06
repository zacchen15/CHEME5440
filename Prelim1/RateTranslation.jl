function translation_rate(protein_length, mRNA_concentration)
    # Parameters given
    characteristic_transcript_length = 1000 # nt
    characteristic_protein_length = 333 # aa
    RNAPII_copy_number = 1150 # copies/cell
    ribosome_copy_number = 45000 # copies/cell
    mass_of_single_cell = 2e-13 # g/cell
    fraction_of_water_per_cell = 0.7 #dimensionless
    transcription_saturation_constant = 0.24 # nmol/gDW
    translation_saturation_constant = 465.64 # nmol/gDW
    characteristic_initiation_time_transcription = 42 # seconds
    characteristic_initiation_time_translation = 15 # seconds
    transcription_elongation_rate = 60 # nt/s
    translation_elongation_rate = 16.5 # aa/s

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the RNAP concentration -
    ribosome_concentration = ribosome_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # Calcualte the specific elongation rate
    kE = translation_elongation_rate / protein_length

    # Calculate kE, kI, time constant
    kI = 1 / characteristic_initiation_time_translation
    tau_factor = kI / kE

    # Saturation constant
    KX = translation_saturation_constant

    # mRNA concentration
    mRNA = mRNA_concentration

    # Calcualte the maximum transcription rate (nmol/gDW-hour)
    rate = kE*ribosome_concentration*(mRNA/(KX*tau_factor+(1+tau_factor)*mRNA))*(3600) # nmol/gDW-hour
    return rate
end
