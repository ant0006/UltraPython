
# PM = P*
# Uses ONLY primitive variables to sove the problem
def PVRS(data_states):
    p_ = 0.5*(data_states.DL + data_states.DR)
    a_ = 0.5*(data_states.CR + data_states.CL)
    
    UM = (data_states.PL - data_states.PR)/(2*a_*p_) + (data_states.UR + data_states.UL)/(2)
    PM = (data_states.PL + data_states.PR)/(2) + (data_states.UL - data_states.UR)/(2*a_*p_)
    pML = data_states.DL + (data_states.UL - UM)*p_/a_
    pMR = data_states.DR + (UM - data_states.UR)*p_/a_
    
    return PM, UM, pML, pMR



    


