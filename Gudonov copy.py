
import numpy as np
import matplotlib.pyplot as plt
from flux import flux

def Guddy(W, cells):
    
    # Intialize variables
    L = 1
    N = cells          # Cells
    dx = L/cells     
    CFL = 0.9      
    
    U_old = np.array([W[0,:], W[0,:]*W[1,:], W[3,:]]) # Density Momentum Energy

    U_new = np.zeros(np.shape(U_old))

    xpos = np.zeros((N))
    # u_avg = np.zeros((N))
    # u_max = np.zeros((N))

    for kk in range(N):
        xpos[kk] = (kk + 0.5)*dx

    # End
    
    # print( U_old[1,0:-2] )
    for tt in range(10):
        
        u = U_old[1]/U_old[0]
        
        u_avg = 0.5*(u[0:-2] + u[1:-1])
        u_avg = max(abs(u_avg))
        u_max = max(abs(u))
        u_array = [u_avg, u_max]
        S_max = max(u_array)
        
        dt = CFL*dx/(S_max)
        print(dt)
        
        for jj in range(cells-1):
            print(jj)
            if jj == 1: # Stationary Ghost boundary at the start
                fp1 = flux(W[:,jj],W[:,jj+1])
                fm1 = flux(W[:,jj],W[:,jj])
            
            elif jj == cells: # Stationary Ghost Boundary at the end
                fp1 = flux(W[:,jj],W[:,jj])
                fm1 = flux(W[:,jj-1],W[:,jj])
            
            else:
                fp1 = flux(W[:,jj],W[:,jj+1])
                fm1 = flux(W[:,jj-1],W[:,jj])

            U_new[:,jj] = U_old[:,jj] - (dt / dx)*( fp1- fm1)
            
            U_new[:,-1] = U_new[:,-2]

        U_old = U_new
        

    # End
    
    plt.plot(xpos,U_new[0,:]) # Density
    plt.show()
