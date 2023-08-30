
# Description: Two Rarefractions exist
def TRRS(data_states):
    z = (data_states.gamma-1)/(2*data_states.gamma)
    PM_top = data_states.CL + data_states.CR - data_states.gamma*z*(data_states.UR - data_states.UL)
    PM_low = (data_states.CL/data_states.PL**z) + (data_states.CR/data_states.PR**z)
    PM = (PM_top/PM_low)**(1/z)
    fR = (2*data_states.CR/(data_states.gamma - 1))*( -1 + (PM/data_states.PR)**z )
    fL = (2*data_states.CL/(data_states.gamma - 1))*( -1 + (PM/data_states.PL)**z )
    
    UM =  0.5*(data_states.UL + data_states.UR) + 0.5*(fR - fL)
    
    pML = data_states.DL * (PM/data_states.PL)**(1/data_states.gamma)
    pMR = data_states.DR * (PM/data_states.PR)**(1/data_states.gamma)
    
    return PM, UM, pML, pMR



    


