function Balances(t,x)

    include("RateTranscription.jl")
    #Define x species vector
    m1 = x[1]
    m2 = x[2]
    m3 = x[3]
    p1 = x[4]
    p2 = x[5]
    p3 = x[6]
    I = x[7]

    #First Order rxn parameters
    k_xd = 1;   #mRNA degradation 1/s
    k_ld = 10;   #protein degradation 1/s
    mu_xd = log(2) / 30*60;  #mRNA dilution
    mu_ld = log(2) / 30*60;  #protein dilution

    #Given Parameters
    copies = 200
    Lx1 = 1200
    Lx2 = 2400
    Lx3 = 600

    #Stoichiometric matrix parameters
    rx1 = RatemRNA(Lx1, copies)
    rx2 = RatemRNA(Lx2, copies)  #M/s per gDCW
    rx3 = RatemRNA(Lx3, copies)
    rL1 = 1
    rL2 = 1
    rL3 = 1

    #Binding Constants for RNAP/transcription
    n = 1.5
    Kd = 0.3
    Kx = 0
    Ki = 1
    K12 = 1
    K32 = 1
    K13 = 1
    K23 = 1

    #Binding Constants for Ribosome/translation
    TL1 = 1
    TL2 = 1
    TL3 = 1
    KL1 = 1
    KL2 = 1
    KL3 = 1

    #Setup fraction bound equations
    fi = (I)^n / ( (Kd)^n + (I)^n )
    f12 = 1
    f32 = 1
    f13 = 1
    f23 = 1

    #Setup Control Functions
    ux1 = (Kx + Ki*fi) / (1 + Kx + Ki*fi)
    ux2 = (Kx + K12*f12 + K32*f32) / (1 + Kx + K12*f12 + K32*f32)
    ux3 = (Kx + K13*f13 + K23*f23) / (1 + Kx + K13*f13 + K23*f23)
    uL1 = (m1) / (TL1*KL1 + m1)
    uL2 = (m2) / (TL2*KL2 + m2)
    uL3 = (m3) / (TL3*KL3 + m3)

    #Setup Mass Balances
    dxdt = similar(x)
    dxdt[1] = -k_xd*m1 + rx1*ux1 #m1
    dxdt[2] = (-k_xd-mu_xd)*m2 + rx2*ux2 #m2
    dxdt[3] = (-k_xd-mu_xd)*m3 + rx3*ux3 #m3
    dxdt[4] = -k_ld*p1 + rL1*uL1 #p1
    dxdt[5] = (-k_ld-mu_ld)*p2 + rL2*uL2 #p2
    dxdt[6] = (-k_ld-mu_ld)*p3 + rL3*uL3 #p3
    dxdt[7] = 0 #I
    dxdt
end
