
import numpy as np
import matplotlib.pyplot as plt
from Exact_Reimann.States import States
from Exact_Reimann.StarPU import StarPU
from PVRS import PVRS
from TRRS import TRRS
from TSRS import TSRS
from AdaptiveReimann import AdaptiveReimann
from HLLSolver import HLL
from pVar import pVariables
from streamLoop import streamLoop



def ReimannSolver(W, U_old, time, cells,scheme):#pL, uL, PL, pR, uR, PR, time, cells, iteration): # L = ii; R= ii+1

    # Parameters of the problem
    domain_len = 1
    gamma = 1.4
    dx = domain_len/float(cells) 

    u_sol = np.zeros(cells,dtype=float)
    p_sol = np.zeros(cells,dtype=float)
    P_sol = np.zeros(cells,dtype=float)

    for ii in range(cells):

        x = (ii + 0.5)*dx

        # Left and Right data states
        if ii == 0 or ii == cells-1:
            # print('The ends')
            DL = W.p[ii]
            UL = W.u[ii]
            PL = W.P[ii]
            DR = W.p[ii]
            UR = W.u[ii]
            PR = W.P[ii]
        else:
            # print('not the ends')
            DL = W.p[ii]
            UL = W.u[ii]
            PL = W.P[ii]
            DR = W.p[ii+1]
            UR = W.u[ii+1]
            PR = W.P[ii+1]

        # Sound speed
        CL = np.sqrt(gamma*PL/DL)
        CR = np.sqrt(gamma*PR/DR)
        
        data_states = States(gamma, DL, UL, PL, CL, DR, UR, PR, CR)

        if ( (2/(gamma - 1.0))*(CL + CR) <= UR - UL):
            print('Vaccuum is generated by initial data')
            exit()

        # Find star solution
        if scheme == 1:
            PM, UM = StarPU(data_states)

        elif scheme == 2: #Primitive Variable 
            PM, UM, pML, pMR = PVRS(data_states)

        elif scheme == 3: #TRRS
            PM, UM, pML, pMR = TRRS(data_states)

        elif scheme == 4: # TSRS
            Ppvrs, _, _, _ = PVRS(data_states)
            PM, UM, pML, pMR = TSRS(data_states, Ppvrs, W)

        elif scheme == 5: # Adaptive Scheme, switches between TRRS & TSRS
            PM, UM, pML, pMR = AdaptiveReimann(data_states, W)

        elif scheme == 6: # UHLL
            # print('HLL')
            U, D, P = HLL(data_states, U_old, domain_len, time, x, cells)


        if scheme == 6:
            pass
        else:
            D, U, P = streamLoop(data_states,time,x, PM, UM, domain_len)

        # E = P/D
        # E /= gamma - 1.0  # Internal energy from EOS
        u_sol[ii] = U 
        p_sol[ii] = D
        P_sol[ii] = P
        
        # x_vec = xpos
        # sol_vec = D, U, P, E
            
    return u_sol, p_sol, P_sol

 
    