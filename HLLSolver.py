import numpy as np
from cVars import v2matrix

def F(U): # U is a matrix
    
    u = U[1]/U[0]
    p = U[0]
    e = U[2]/U[0] - 0.5*u**2
    P = abs(e)*U[0]*(1.4 - 1)

    F = np.array([U[1], 
        P + p*u**2,
        u*U[2]])
    
    return F

def UHLL(UL, UR, SL, SR): # UR/UL is a matrix 
    FL = F(UL)
    FR = F(UR)
    UHLL = (SR*UR - SL*UL + FL - FR)/(SR - SL)
    return UHLL


def HLL(data_states, U, domain_len, time, x, cells):
    # UL = n
    # UR = n + 1
    x0 = domain_len/2
    mid_cells = len(U.mass)//2
    
    xqt = (x-x0)/time
    uL = U.momentum[:mid_cells]/U.mass[:mid_cells]
    uR = U.momentum[mid_cells:cells]/U.mass[mid_cells:cells]

    SL = max(uL)
    SR = max(uR)
    
    if xqt < SL:
        u_sol = data_states.UL
        p_sol = data_states.DL
        P_sol = data_states.PL

    elif xqt > SR:
        u_sol = data_states.UR
        p_sol = data_states.DR
        P_sol = data_states.PR
    else:
        UL = v2matrix(data_states.UL,data_states.DL,data_states.PL)
        UR = v2matrix(data_states.UR,data_states.DR,data_states.PR)
        UH = UHLL(UL, UR, SL, SR)
        u_sol = UH[1]/UH[0]
        p_sol = UH[0]

        e = UH[2]/UH[0] - 0.5*u_sol**2
        P_sol = abs(e)*UH[0]*(1.4 - 1)

    return u_sol, p_sol, P_sol
    

            