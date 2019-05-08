function Balances(t,x)
    # Constants
    vmax1 = 5.0
    vmax2 = 5.0
    vmax3 = 1.0
    vmax4 = 1.0
    Ks1 = 5.0
    Ks2 = 5.0
    Ks3 = 5.0
    Ks4 = 5.0
    Ki1 = 1.0
    Ki2 = 1.0
    Stot = 100.0

    # Calculate the reaction rates
    v1 = (vmax1*x[1]) / ( (1 + (x[4]/Ki1)) * (Ks1 + x[1]) )
    v2 = (vmax2*x[1]) / ( (1 + (x[5]/Ki2)) * (Ks2 + x[1]) )
    v3 = (vmax3*x[2]) / ( (Ks3 + x[2]) )
    v4 = (vmax4*x[3]) / ( (Ks4 + x[3]) )

    #Setup Mass Balances
    dxdt = similar(x)
    dxdt[1] = -(v1 - v3) - (v2 - v4)    #A
    dxdt[2] = v1 - v3                   #B
    dxdt[3] = v2 - v4                   #C
    dxdt[4] = 0.0                       #Inhibitor 1
    dxdt[5] = 0.0                       #Inhibitor 2
    dxdt
end
