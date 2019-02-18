function RatemRNA(gene_length, copies)
    #Find the elongation rate constant
    e_x = 27 # nt/s from BioNumbers (Proshkin S., et al.)
    L = 1000 # average gene length in prokaryotes in bp or nt from BioNumbers (Xu L., et al)
    L_j = gene_length # nt given in HW
    k_Ej = (e_x/L)*(L/L_j) # 1/s

    #Finding K_xj from McClure
    intercept = 1.04e-6 # M-s
    slope = 42 # s
    K_xj = intercept/slope # Molar

    #Parameters for Txj
    k_I = 0.04 # 1/s at 37 degC given by BioNumbers (McClure, et al.)
    k_A = 0 # assuming abortive initiation is low
    T_xj = (k_Ej+k_A)/k_I # dimensionless

    Vol_cell = 6.7e-16 # L/cell
    G_j = copies # copies/cell
    G_j = G_j*(1/Vol_cell)*(1/6.02e23) # M/cell
    R_xt = 30e-9 # M/cell from BioNumbers
    #R_xt = R_xt*(1/Vol_cell)*(1000)*(1/0.3) # M/gDCW

    r_xj = k_Ej*R_xt*G_j*(1/((K_xj*T_xj)+(T_xj+1)*G_j)) # M/s per cell
    r_xj = r_xj*(1/6.02e23)*(1/Vol_cell)*(1/1000)*(1/0.3) # M/s per gCDW
    r_xj
end
