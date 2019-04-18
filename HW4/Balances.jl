function Balances(t,x)
    #Define ligand concentration
    ligand = 0


    #Define x species vector
    E0 = x[1]
    E1 = x[2]
    E1_star = x[3]
    B = x[4]
    Bp = x[5]
    E1_starB = x[6]
    E1_starBp = x[7]

    #Constants from Barkai (2001) page 877
    ar = 0.2 #s^-1 * uM^-1
    dr = 0.1
    kr = 0.1
    ab = 1000
    db = 1
    kb = 0
    abp = 100
    dbp = 0.01
    kbp = 1
    k_plus = 1
    k_minus = 1
    a0_plus = 10
    a1_plus = 1 / (1+ligand)
    a0_minus = 0
    a1_minus = ligand / (1+ligand)
    B0 = 0
    B1 = 2.5*ligand / (1+ligand)

    ab = 100
    abp = 10
    k_plus = 10
    k_minus = 1000

    Etot = 10
    R = 0.2
    B = 2
    A = E1_star

    #Setup Mass Balances
    dxdt = similar(x)
    dxdt[1] = -kr*R + kb*E1_starB + kbp*E1_starBp
    dxdt[2] = a1_minus*E1_star - a1_plus*E1 + B1*E1_starB + B1*E1_starBp + kr*R
    dxdt[3] = a1_plus*E1 - a1_minus*E1_star - ab*E1_star*B - abp*E1_star*Bp + db*E1_starB + dbp*E1_starBp
    dxdt[4] = k_minus*Bp - k_plus*A*B + (B1 + db + kb)*E1_starB - ab*E1_star*B
    dxdt[5] = k_plus*A*B - k_minus*Bp + (B1 + dbp + kbp)*E1_starBp - abp*E1_star*Bp
    dxdt[6] = ab*E1_star*B - (B1 + db + kb)*E1_starB
    dxdt[7] = abp*E1_star*Bp - (B1 + dbp + kbp)*E1_starBp
    dxdt

end
