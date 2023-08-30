########################################
# Author: Anthony Flores
# Class: Gas Dynamics
# Description: To gain a basic understanding of Finite Volume method, solve for 1D shocks and/or rarefractions travelling across two different
# fluids using seven different solvers. 
########################################


from ReimannSolvers import ReimannSolver
import numpy as np
from Gudonov import Guddy 
import matplotlib.pyplot as plt
from pVar import pVariables
from cVars import conVar

def FVM():
     # ccase (differnt solvers)
    # 1: Exact Reimann
    # 2: PVRS
    # 3: TRRS
    # 4: TSRS
    # 5: Adaptive
    # 6: HLL
    # 7: HLLC
    
    # Initialization
    ccase = 6
    test = 0
    cells = 100
    gamma = 1.4
    u = np.zeros(cells,dtype=float)
    p = np.zeros(cells,dtype=float)
    P = np.zeros(cells,dtype=float)

    # Inital value given
    pL = (1, 1, 1, 1, 5.99924)
    uL = (0.75, -2, 0, 0, 19.5975)
    PL = (1, 0.4, 1000, 0.01, 460.894)
    pR = (0.125, 1, 1, 1, 5.99242)
    uR = (0, 2, 0, 0, -6.19633)
    PR = (0.1, 0.4, 0.01, 100, 46.0950)

    x = np.linspace(0,1, num = cells)

    tf = 0.15
    t0 = 0.001
    t = np.linspace(t0,tf,num=10)
    
    # Get initial primitive values across the domain
    u[:(len(u)//2)] = uL[test]
    u[(len(u)//2):len(u)] = uR[test]
    
    p[:(len(u)//2)] = pL[test]
    p[(len(u)//2):len(u)] = pR[test]

    P[:(len(u)//2)] = PL[test]
    P[(len(u)//2):len(u)] = PR[test]
    
    W_old = pVariables(u,p,P)
    
    # Plot velocity Initial Condition
    plt.plot(x,W_old.u)

    # Get the Conservative values (p, pu, E) across the domain
    U = conVar(W_old)

    # Time loop, solves for values across a time interval
    for tt in range(len(t)):
        
        time = t[0] + t[tt] 
        u_sol, p_sol, P_sol = ReimannSolver(W_old,U, time, cells, ccase) # ccase
        
        W_reimann = pVariables(u_sol,p_sol, P_sol)
        U_new = Guddy(W_reimann, U, cells, ccase)
        
        # re-initialize primitive and conservative variables
        U.mass = U_new[0,:]
        U.momentum = U_new[1,:]
        U.E = abs(U_new[2,:])

        W_old.u = U.momentum/U.mass
        W_old.p = U.mass
        e = U.E/U.mass - 0.5*W_old.u**2
        W_old.P = abs(e)*U.mass*(gamma - 1)
        
    plt.plot(x,W_old.u)
    plt.legend(['Initial Condition','Solution after 0.15s'])
    plt.xlabel('x')
    plt.ylabel('u')
    plt.show()


if __name__ == "__main__":
    FVM()
   
