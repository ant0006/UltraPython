from numpy import array

def flux(W1, W2):
# Returns flux
# For our case it will return flux (F) of W_i+.5 and W_i-.5
# W_i+.5 = average(W_i , W_i+1)
# W_i-.5 = average(W_i-1 , W_i)
# U = [p  u  E  P]  -> (Density, velocity, internal Energy, Pressure)
# F = [pu  pu^2 + P  u(E+P)]
     
     W = 0.5*(W1+W2)
     
     u = W[0]    # velocity
     p = W[1]    # density
     P = W[2]    # Pressure

     E1 = P*(1.4/(1.4 - 1.0))
     E2 = 0.5*p*u**2
     ET = E1+E2
     # E =  P/p
     # E /= 1.4 - 1.0
     
     F = array([p*u, 
          P + p*u**2,
          u*ET])
     # u*(E+P)
     return F
