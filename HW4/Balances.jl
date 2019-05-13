function Balances(t,x)
    #Define ligand concentration
    if t < 100
        ligand = 0.0
    else
        ligand = 1.0
    end

    #Define x species vector
    E0 = x[1]
    E1 = x[2]
    E1_star = x[3]
    B = x[4]
    Bp = x[5]
    E1_starB = x[6]
    E1_starBp = x[7]

    #Constants from Barkai (2001) page 877
    kr = 0.1
    ab = 1000
    db = 1
    kb = 0
    abp = 100
    dbp = 0.01
    kbp = 1
    k_plus = 0#1
    k_minus = 0#1
    a1_plus = 1 / (1+ligand)
    a1_minus = ligand / (1+ligand)
    b1 = 2.5*ligand / (1+ligand)

    Etot = 10
    R = 0.2
    A = E1_star
    Vrmax = kr*R

    #Setup Mass Balances
    dxdt = similar(x)
    dxdt[1] = kb*E1_starB + kbp*E1_starBp - Vrmax
    dxdt[2] = a1_minus*E1_star - a1_plus*E1 + b1*E1_starB + b1*E1_starBp + Vrmax
    dxdt[3] = a1_plus*E1 - a1_minus*E1_star - abp*E1_star*Bp + dbp*E1_starBp - ab*E1_star*B + db*E1_star*B
    dxdt[4] = k_minus*Bp - k_plus*A*B - ab*E1_star*B + (b1+db+kb)*E1_starB
    dxdt[5] = k_plus*A*B - k_minus*Bp - abp*E1_star*Bp + (b1+dbp+kbp)*E1_starBp
    dxdt[6] = ab*E1_star*B - (b1+kb+db)*E1_starB
    dxdt[7] = abp*E1_star*Bp - (b1+kbp+dbp)*E1_starBp
    dxdt

    # dxdt[1]= kbp*E1_starBp + kb*E1_starB - Vrmax
    # dxdt[2]= a1_minus*E1_star - a1_plus*E1 + E1_starB*b1 + E1_starBp*b1 + Vrmax
    # dxdt[3]= a1_plus*E1 - a1_minus*E1_star - E1_star*Bp*abp + E1_starBp*dbp - E1_star*B*ab + E1_starB*db
    # dxdt[4]= -ab*E1_star*B + kb*E1_starB + E1_starB*b1 + db*E1_starB + k_minus*Bp - k_plus*E1_star*B
    # dxdt[5]= -abp*E1_star*Bp + (kbp+b1+dbp)*E1_starBp - k_minus*Bp + k_plus*E1_star*B
    # dxdt[6]= -db*E1_starB +ab*E1_star*B - b1*E1_starB - kb*E1_starB
    # dxdt[7]= -dbp*E1_starBp + abp*E1_star*Bp - b1*E1_starBp - kbp*E1_starBp
    # dxdt
end
