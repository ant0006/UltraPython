from Exact_Reimann.Sample import Sample

def streamLoop(data_states,time,x, PM, UM, domain_len):

    x0 = domain_len/2
    Pratio = data_states.PR/ data_states.PL
    # S = max(data_states.UR -data_states.CR,data_states.UL+data_states.CL)#(x-x0)/time 
    # S = data_states.UL +  data_states.CL* np.sqrt(1 + (Pratio - 1)*( data_states.gamma+1)/(2* data_states.gamma))
    S1 = 0.5*(data_states.UR + data_states.UL)
    S = abs(max(data_states.UR,data_states.UL,S1))
    D, U, P = Sample(PM, UM, S, data_states)
    
    E = P/D
    E /=  data_states.gamma - 1.0  # Internal energy from EOS

    return D, U, P
    