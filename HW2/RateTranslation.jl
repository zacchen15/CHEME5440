function RateAA(mRNA_length, mRNA_conc)
    #mRNA concentration should be input in units of mM
    #Find the elongation rate constant
    e_x = 14.5 # amino acids/s from BioNumbers (Dalbow DG et al.)
    L = 1000/3 # average gene length in prokaryotes in bp or nt from BioNumbers (Xu L., et al)
    L_j = mRNA_length/3 # amino acid given in HW
    k_Ej = (e_x/L)*(L/L_j) # 1/s

    #Finding K_xj from McClure
    intercept = 1.04e-3 # mM-s
    slope = 42 # s
    K_xj = intercept/slope # mM
    #print(K_xj)

    #Parameters for Txj
    k_I = 0.04 # 1/s at 37 degC given by BioNumbers (McClure, et al.)
    k_A = 0 # assuming abortive initiation is low
    T_xj = (k_Ej+k_A)/k_I # dimensionless
    #print(T_xj)

    Vol_cell = 9e-17 # L/cell
    G_j = mRNA_conc # mM
    #G_j = G_j*(1/Vol_cell)*(1/6.02e23) # M/cell
    R_xt = 50000 # number of ribosome per cell from HW1 solutions
    R_xt = R_xt*(1/Vol_cell)*(1/6.02e23)*1000 # mM
    #print(R_xt)

    r_xj = k_Ej*R_xt*G_j*(1/((K_xj*T_xj)+(T_xj+1)*G_j)) # M/s per cell
    r_xj
end
