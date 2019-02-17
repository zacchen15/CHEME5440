function Balances(t,x)

    #Define x species vector
    m1 = x[1]
    m2 = x[2]
    m3 = x[3]
    p1 = x[4]
    p2 = x[5]
    p3 = x[6]

    #First Order rxn parameters
    k_xd = 1;
    k_ld = 1;
    mu_xd = 1;
    mu_ld = 1;

    #Stoichiometric matrix parameters
    rx1 = 1
    rx2 = 1
    rx3 = 1
    rL1 = 1
    rL2 = 1
    rL3 = 1

    #Setup Control Functions
    ux1 = 1
    ux2 = 1
    ux3 = 1
    uL1 = 1
    uL2 = 1
    uL3 = 1

    #Setup Mass Balances
    dxdt = similar(x)
    dxdt[1] = -k_xd*m1 + rx1*ux1 #m1
    dxdt[2] = (-k_xd-mu_xd)*m2 + rx2*ux2 #m2
    dxdt[3] = (-k_xd-mu_xd)*m3 + rx3*ux3 #m3
    dxdt[4] = -k_xd*p1 + rL1*uL1 #p1
    dxdt[5] = (-k_xd-mu_xd)*p2 + rL2*uL2 #p2
    dxdt[6] = (-k_xd-mu_xd)*p3 + rL3*uL3 #p3
    dxdt
end
