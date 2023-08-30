# F = Function
# FD = dF/dP
# P = Current guess of pressure
# PK = Pressure data 
# CK = Sound speed
from numpy import sqrt

def PressureFunction(P, DK, PK, CK, gamma):

    G0 = gamma + 1.0
    G1 = gamma - 1.0
    G2 = G1/(2*gamma)
    G3 = G0/(2*gamma)
    G4 = 2.0/G0
    G5 = G1/G0

    F = float('nan')
    FD = float('nan')

    if (P <= PK):
        
        # Rarefaction wave
        P_ratio = P/PK
        F = 2.0*CK/G1*(P_ratio**G2 - 1)
        FD = (1.0/(DK*CK))*P_ratio**(-G3)
        

    else:

        # Shock wave
        AK = G4/DK
        BK = G5*PK
        Q = sqrt(AK/(P+BK))
        F = (P - PK)*Q
        FD = (1.0 - 0.5*(P - PK)/(BK + P))*Q

    return F, FD   
