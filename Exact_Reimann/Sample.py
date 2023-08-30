
# PM : Presure in * region
# UM : Velocity in * region
from numpy import sqrt

def Sample(PM, UM, S, data_states):

    D, U, P = float('nan'), float('nan'), float('nan')

    G1 = (data_states.gamma - 1.0)/(2.0*data_states.gamma)
    G2 = 1.0/data_states.gamma
    G3 = 2.0/(data_states.gamma + 1.0)
    G4 = (data_states.gamma - 1.0)/2.0
    G5 = (data_states.gamma - 1.0)/(data_states.gamma + 1.0)
    G6 = 2.0/(data_states.gamma - 1.0)
    G7 = G6*data_states.gamma
    G8 = (data_states.gamma + 1)/(2*data_states.gamma)
    G9 = (data_states.gamma - 1.0)/(2*data_states.gamma)

    if S <= UM :
        # Sampling point lies to the left of the contact discontinuity

        if PM <= data_states.PL :

            # Left rarefaction

            SHL = data_states.UL - data_states.CL # Head characteristic velocity

            if (S <= SHL):
                
                # Sampled point is the left data state

                D = data_states.DL
                U = data_states.UL
                P = data_states.PL
            else:

                CML = data_states.CL * (PM/data_states.PL)**G1
                STL = UM - CML    # Tail characteristic velocity

                if (S >= STL):

                    # Sampled point is * left state

                    D = data_states.DL * (PM/data_states.PL)**(G2)
                    U = UM
                    P = PM
                else :

                    # Sampled point is inside the left fan

                    U =  G3 * (data_states.CL  + G4* data_states.UL + S)
                    C = G3 + G5*(data_states.UL - S)/data_states.CL
                    D = data_states.DL * C**G6
                    P = data_states.PL * C**G7
        else:

            # Left shock

            PML = PM/data_states.PL
            SL = data_states.UL - data_states.CL*sqrt(PML*G8 + G9)  

            if (S <= SL):

                # Sampled point is the left data state

                D = data_states.DL
                U = data_states.UL
                P = data_states.PL
            else:

                # Sampled point is the * left state

                D = data_states.DL*((PML + G5)/(G5*PML + 1.0))
                U = UM
                P = PM
    else:

        # Sampling point lies to the right of the contact discontinuity
        
        if (PM >= data_states.PR):

            # Right shock
            PMR = PM/data_states.PR
            SR = data_states.UR + data_states.CR*sqrt(PMR*G8 + G9) 

            if (S >= SR):

                # Sampled point is the right data state

                D = data_states.DR
                U = data_states.UR
                P = data_states.PR
            else:

                # Sampled point is the * right state

                D = data_states.DR*((PMR + G5)/(G5*PMR + 1.0))
                U = UM
                P = PM
        else:

            # Right rarefaction        
            SHR = data_states.UR + data_states.CR # Head characteristic velocity

            if (S >= SHR):
                
                # Sampled point is the right data state

                D = data_states.DR
                U = data_states.UR
                P = data_states.PR
            else:

                CMR = data_states.CR * (PM/data_states.PR)**G1
                STR = UM + CMR    # Tail characteristic velocity

                if (S <= STR):

                    # Sampled point is * right state

                    D = data_states.DR * (PM/data_states.PR)**(G2)
                    U = UM
                    P = PM
                else :

                    # Sampled point is inside the right fan

                    U =  G3 * (-data_states.CR  + G4* data_states.UR + S)
                    C = G3 - G5*(data_states.UR - S)/data_states.CR
                    D = data_states.DR * C**G6
                    P = data_states.PR * C**G7

    return  D, U, P              
