
from .PressureFunction import PressureFunction

tol_pressure = 1e-6
n_iter = 100

def StarPU(data_states):

    P, U = float('nan'), float('nan')
    
    # First find P*
    Pstart = 0.5*(data_states.PR + data_states.PL)
    Pold = Pstart 
    Udiff = data_states.UR - data_states.UL
    
    FL  = float('nan')
    FLD = float('nan')
    FR  = float('nan')
    FRD = float('nan')  
    delta_P = float('nan')
    converged_its = -100 

    counter = 0
    for i in range(n_iter):

        counter += 1

        FL, FLD = PressureFunction(Pold, data_states.DL, data_states.PL, data_states.CL, data_states.gamma)
        FR, FRD = PressureFunction(Pold, data_states.DR, data_states.PR, data_states.CR, data_states.gamma)

        P = Pold - (FL + FR + Udiff)/(FLD + FRD)
        
        delta_P = 2.0*abs( (P - Pold)/(P + Pold)) 

        if delta_P <= tol_pressure:
            converged_its = i+1
            break

        if P <= 0.0:
            P = tol_pressure

        Pold = P    
    
    # print(f'Newton-Raphson method for finding pressure converged in {converged_its} iterations.')

    if delta_P > tol_pressure or counter > n_iter:
        print('Newton-Raphson method for pressure did not converge.')
        print('Terminating...')
        exit()

    
    # Find U* next
    U = 0.5*(data_states.UL + data_states.UR + FR - FL)

    return P, U
    


