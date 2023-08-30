import numpy as np
from Exact_Reimann.States import States
from streamLoop import streamLoop
from PVRS import PVRS

# Description: Two Shocks exist

def Ak(pK, gamma):
    A = 2/(pK*(gamma + 1))
    
    return A

def Bk(PK, gamma):
    B = PK*(gamma - 1)/(gamma + 1)
    
    return B
    
def gk(Po, A, B):
    
    g = (A/(Po + B))**0.5
    
    return g

def TSRS(data_states, Ppvrs, W):

    tol_pressure = 1e-6
    n = 20
    counter = 0
    Udiff = W.u[-2] - W.u[1]
    Usum = W.u[-2] + W.u[1]
    # Udiff = data_states.UR - data_states.UL
    # Usum =  data_states.UR + data_states.UL
    
    gamma = data_states.gamma

    Po = max(0, Ppvrs)  # This needs to change to Po = max(0,Ppvrs)
    
    for ii in range(n):

        counter += 1
        converged_its = 0
        AL = Ak(data_states.DL, gamma)
        AR = Ak(data_states.DR, gamma)
        
        BL = Bk(data_states.PL, gamma)
        BR = Bk(data_states.PR, gamma)
        
        gL = gk(Po, AL, BL)
        gR = gk(Po, AR, BR)
        N = Udiff - data_states.PL*gL + data_states.PR*gR
        # print(data_states.PL*gL, data_states.PR*gR)
        D = gR + gL
        PM = N/D

        PLdiff = PM - data_states.PL
        PRdiff = PM - data_states.PR
        UM = 0.5*Usum + 0.5*(PRdiff*gR - PLdiff*gL)
    
        delta_P = 2.0*abs( (PM - Po)/(PM + Po))
        
        if delta_P <= tol_pressure:
            converged_its = ii+1
            break

        if PM <= 0.0:
            PM = tol_pressure

        Po = PM  

    # if delta_P > tol_pressure or counter > n:
    #     print('Your iteration method for pressure did not converge.')
    #     print('Terminating...')
    #     exit()

    # PLdiff = PM - data_states.PL
    # PRdiff = PM - data_states.PR
    # UM = 0.5*(data_states.UL + data_states.UR) + 0.5*()
    TopL = (PM/data_states.PL) +((gamma - 1)/(gamma + 1))
    BotL = ((gamma - 1)/(gamma + 1))*(PM/data_states.PL) + 1
    TopR = (PM/data_states.PR) +((gamma - 1)/(gamma + 1))
    BotR = ((gamma - 1)/(gamma + 1))*(PM/data_states.PR) + 1
    
    pML = data_states.DL *(TopL/ BotL)
    pMR = data_states.DR * (TopR/ BotR)

    return PM, UM, pML, pMR


