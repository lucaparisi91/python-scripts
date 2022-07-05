from pimc import gaussianPotentialImpl

class gaussianPotential:

    def __init__(self,a,R):
        '''
        a: scattering length
        R: range of the interaction
        '''
        self.a=a
        self.R=R
        
        self.V0=gaussianPotentialImpl.V0_R0(self.a)/self.R**2
        self.alpha=1/( 2*self.R**2 )
    
    def toJson(self):
        return {
                "kind": "gaussian",
                "V0": self.V0,
                "alpha": self.alpha
            }

class harmonicPotential:
    def __init__(self, omega=1):
        '''
        omega : frequency
        V(r) = 0.5 * omega**2 * r**2
        '''
        self.omega=omega
    def toJson(self):
        return {
            "kind" : "harmonic",
            "omega" : self.omega
        }
    