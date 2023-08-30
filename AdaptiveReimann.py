
import numpy as np
from Exact_Reimann.States import States

from PVRS import PVRS
from TRRS import TRRS
from TSRS import TSRS 


def AdaptiveReimann(data_states, W):
    
    # Find Q = pmax / pmin < Quser = 2
    Quser = 2
    Pmin = min(data_states.PL, data_states.PR)
    Pmax = max(data_states.PL, data_states.PR)
    Q = Pmax/Pmin
    
    Ppvrs, Upvrs, ppvrsL, ppvrsR = PVRS(data_states)
    
    # Find star solution
    if Q < Quser and Pmin < Ppvrs and Ppvrs < Pmax:
        PM = Ppvrs
        UM = Upvrs
        pML = ppvrsL
        pMR = ppvrsR
    else:
        if Ppvrs > Pmax: #Two Shocks
            PM, UM, pML, pMR = TSRS(data_states,Ppvrs,W)
            
        else: #Two Rarefractions
            PM, UM, pML, pMR = TRRS(data_states)


    return PM, UM, pML, pMR


    


