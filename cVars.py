from numpy import array


class conVar():
    #Conservative values mass, momentum & energy
    def __init__(self, W):
        self.mass = W.p
        self.momentum = W.p*W.u
        self.E =  W.P/W.p
        self.E /= (1.4 - 1.0)

def v2matrix( u, p, P):
    # Converts a vector to matix
    mass = p
    momentum = p*u
    E =  P/p
    E /= (1.4 - 1.0)
    U = array([mass, momentum, E])
    return U