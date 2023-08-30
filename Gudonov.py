
import numpy as np
import matplotlib.pyplot as plt
from flux import flux
from HLLSolver import F
from cVars import conVar
from pVar import pVariables


def Guddy(W_reimann, U,  cells, ccase):
    
    # Intialize variables
    L = 1
    N = cells          # Cells
    dx = L/cells     
    CFL = 0.8           # Stability condition
    
    U_old = np.array([U.mass, U.momentum, U.E]) # Density Momentum Energy
    W_reim = np.array([W_reimann.u, W_reimann.p, W_reimann.P]) # Velocity Density Pressure
    
    U_new = np.zeros(np.shape(U_old))
    
    # ONLY IF HLL
    if ccase == 6:
        W2 = pVariables(W_reimann.u, W_reimann.p, W_reimann.P)
        U_reim = conVar(W2)
        U_reim = np.array([U_reim.mass, U_reim.momentum, U_reim.E])
        # print(U_reim)
        mid_cells = len(U.mass)//2
        uL = U.momentum[:mid_cells]/U.mass[:mid_cells]
        uR = U.momentum[mid_cells:cells]/U.mass[mid_cells:cells]

        SL = max(uL)
        SR = max(uR)
        
    else:
        pass

    u = W_reimann.u
    
    u_avg = 0.5*(u[0:-2] + u[1:-1])
    u_avg = max(abs(u_avg))
    u_max = max(abs(u))
    u_array = [u_avg, u_max]
    S_max = max(u_array)
    
    dt = CFL*dx/(S_max*10)
    
    for jj in range(cells):
        
        # Finding the Fluxes
        if jj == 0: # Stationary Ghost boundary at the start
                fp1 = flux(W_reim[:,jj],W_reim[:,jj+1])
                fm1 = flux(W_reim[:,jj],W_reim[:,jj])
            
                df = fp1 - fm1
                

        elif jj == cells - 1: # Stationary Ghost Boundary at the end
                fp1 = flux(W_reim[:,jj],W_reim[:,jj])
                fm1 = flux(W_reim[:,jj-1],W_reim[:,jj])
                df = fp1 - fm1

        else:
            if ccase == 6: # For HLL Solver
                if SL <= 0:
                    # print('SL')
                    df = F(U_reim[:,jj]) # Gets flux from UL
                elif SR <= 0:
                    # print('SR')
                    df = F(U_reim[:,jj+1]) # Gets flux from UR
                else:
                    # Gets flux from both sides
                    df = (SR*F(U_reim[:,jj]) - SL*F(U_reim[:,jj+1] + SL*SR*(U_reim[:,jj+1]-U_reim[:,jj])) ) / (SR -SL)      
            
            else:
                fp1 = flux(W_reim[:,jj+1],W_reim[:,jj])
                fm1 = flux(W_reim[:,jj-1],W_reim[:,jj])
                df = fp1 - fm1
        
        U_new[:,jj] = U_old[:,jj] - (dt / dx)*(df)

    U_new[:,0] = U_old[:,0] 
    U_new[:,-1] = U_old[:,-1]
    
    return U_new
    # End

    
