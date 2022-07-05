
class magnetization:
    def __init__(self, groups , label ):
        self.setA,self.setB=groups
        self.label=label

    def toJson(self):
        return   {
            "kind": "magnetization",
            "label": self.label,
            "groupA" : self.setA,
            "groupB" : self.setB
        }


class magnetizationSquared:

    def __init__(self, groups , label ):
        self.setA,self.setB=groups
        self.label=label

    def toJson(self):
        return   {
            "kind": "magnetizationSquared",
            "label": self.label,
            "groupA" : self.setA,
            "groupB" : self.setB
        }

class pairCorrelation:

    def __init__(self,groups,label,xRange,bins=100 ):

        self.setA,self.setB=groups
        self.label=label
        self.bins=bins
        self.xRange=xRange

    def toJson(self):

        return {
            "kind": "pairCorrelation",
            "label": self.label,
            "setA" : self.setA,
            "setB" : self.setB,
            "bins" : self.bins,
            "minx" : self.xRange[0],
            "maxx" : self.xRange[1]
        }

class oneBody:
    def __init__(self,group,label,xRange,bins=100 ):
        self.group=group
        self.label=label
        self.bins=bins
        self.xRange=xRange
    
    def toJson(self):
        return {
            "kind": "oneBody",
            "label": self.label,
            "set" : self.group,
            "bins" : self.bins,
            "minx" : self.xRange[0],
            "maxx" : self.xRange[1]
        }


class superfluidFraction:
    def __init__(self,groups,label ):
        self.groups=groups
        self.label=label
        
    def toJson(self):
        return {
            "kind": "superfluidFraction",
            "label": self.label,
            "size" : len(self.groups)*(len(self.groups) + 1  )/2
        }

class angleEstimator:
    def __init__(self,groups,label,bins=100 ):
        self.setA,self.setB=groups
        self.label=label
        self.bins=bins

    def toJson(self):

        return {
            "kind": "angleEstimator",
            "label": self.label,
            "setA" : self.setA,
            "setB" : self.setB,
            "bins" : self.bins,
            "minx" : -1.0,
            "maxx" : 1.0
        }


class magnetizationDistribution:

    def __init__(self,groups,label,magRange):

        self.groups=groups
        self.label=label
        self.range=magRange

    def toJson(self):
        return {
            "kind" : "magnetizationDistribution",
            "label" : self.label,
            "sets" : self.groups,
            "min" : self.range[0],
            "max" : self.range[1]
        }


class thermodynamicEnergy:

    def __init__(self,label,N=None, magnetizationDistribution=False):
        self.label=label
        self.magnetizationDistribution=magnetizationDistribution
        self.N=N

    def toJson(self):

        if (self.magnetizationDistribution):

            if self.N is None:
                raise RuntimeError("N needs to be specified")
            return {
            "kind" : "thermodynamicEnergyMagnetization",
            "label" : self.label,
            "size" : self.N + 1
        }
        else:
            return {
            "kind" : "thermalEnergy",
            "label" : self.label
        }        


class virialEnergy:
    def __init__(self,label,N=None, magnetizationDistribution=False):
        self.label=label
        self.magnetizationDistribution=magnetizationDistribution
        self.N=N

    def toJson(self):

        if (self.magnetizationDistribution):

            if self.N is None:
                raise RuntimeError("N needs to be specified")
            return {
            "kind" : "virialEnergyMagnetization",
            "label" : self.label,
            "size" : self.N + 1
        }
        else:
            return {
            "kind" : "virialEnergy",
            "label" : self.label
        }        
    