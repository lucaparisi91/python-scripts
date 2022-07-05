import copy


class restriction:
    def __init__(self , sets, particleRanges):

        self.particleRanges=particleRanges
        self.sets=sets
    def toJson(self):
        return {
            "sets" : self.sets,
            "minParticleNumber" : [  setRange[0]    for setRange in self.particleRanges ],
            "maxParticleNumber" : [  setRange[1]    for setRange in self.particleRanges ]
        }


class semiOpenClose:
    def __init__( self, setA, setB , C, l=1 ):
        self._jOpen={"setA" : setA, "setB" : setB , "C" : C ,  "reconstructionMaxLength" : l,"kind" : "semiOpen"}
        
        self.name="semiOpenClose"

    def toJson(self):
        jClose=copy.deepcopy(self._jOpen)
        jClose["kind"]="semiClose"
        return [self._jOpen,jClose ] 

class levy:
    def __init__( self, l = 1 ):
        self._j={ "reconstructionMaxLength" : l,"kind" : "levy"}
        self.name="levy"

    def toJson(self):
        return [self._j]

class swap:
    def __init__( self, l = 1 ):
        self._j={ "reconstructionMaxLength" : l,"kind" : "swap"}
        self.name="swap"

    def toJson(self):
        return [self._j]

class advanceRecedeHead:
    def __init__( self, l = 1 , restriction=None):
        self._jAdvance={ "reconstructionMaxLength" : l,"kind" : "advanceHead"}
        self._jRecede={ "reconstructionMaxLength" : l,"kind" : "recedeHead"}
        self.name="advanceRecedeHead"

        

    def toJson(self):
        return [self._jAdvance,self._jRecede]

    
class advanceRecedeHeadTail:
    def __init__( self, setA , setB, l = 1 , restriction=None ):
        self._jAdvance={ "reconstructionMaxLength" : l,"kind" : "advanceHeadTail","setA" : setA, "setB" : setB}
        self._jRecede={ "reconstructionMaxLength" : l,"kind" : "recedeHeadTail", "setA" : setA, "setB" : setB}
        self.name="advanceRecedeHeadTail"
        if restriction is not None:
            self._jAdvance["restriction"]=restriction.toJson()
            self._jRecede["restriction"]=restriction.toJson()
            

    def toJson(self):
        return [ self._jAdvance, self._jRecede ]


class fullOpenClose:
    def __init__( self, setA , setB, C,  l = 1 , restriction=None ):
        self._jOpen={ "reconstructionMaxLength" : l,"kind" : "fullOpen","setA" : setA, "setB" : setB,"C" : C }
        self._jClose={ "reconstructionMaxLength" : l,"kind" : "fullClose","setA" : setA, "setB" : setB,"C" : C}
        self.name="fullOpenClose"
        if restriction is not None:
            self._jClose["restriction"]=restriction.toJson()
            
    def toJson(self):
        return [ self._jOpen, self._jClose ]

class createDeleteTwoWorms:
    def __init__( self, setA , setB, C , l = 1 ,restriction=None):
        self._jCreate={ "reconstructionMaxLength" : l,"kind" : "createTwoWorms","setA" : setA, "setB" : setB , "CA" : C[0],
            "CB" : C[1], "initialBeadSampling" : { "kind" : "uniform"} }
        self._jDelete=copy.copy(self._jCreate)
        self._jDelete["kind"]="deleteTwoWorms"
        self.name="createDeleteTwoWorms"
        if restriction is not None:
            self._jCreate["restriction"]=restriction.toJson()
            self._jDelete["restriction"]=restriction.toJson()


    def toJson(self):
        return [ self._jCreate, self._jDelete ]



class translate:
    def __init__( self, delta , l= 1 ):
        self._j={ "delta" : delta,"kind" : "translate"}
        self.name="translate"
    def toJson(self):
        return [self._j]
    

class moveHead:
    def __init__( self, l = 1 ):
        self._j={ "reconstructionMaxLength" : l,"kind" : "moveHead"}
        self.name="moveHead"
    
    def toJson(self):
        return [self._j]

    
class moveTail:
    def __init__( self, l = 1 ):
        self._j={ "reconstructionMaxLength" : l,"kind" : "moveTail"}
        self.name="moveTail"
    def toJson(self):
        return [self._j]


class tableMove:
    
    def __init__(self):
        self.moves=[]

    def addMove(self,move, sector, p, group=0, *args , **kwds):

        self.moves.append( {"move" : move ,"p": p , "sector" : sector , "set" :  group } )
    
    def table( self):

        rows= [ { **move, "move" : move["move"].name  } for move in self.moves ]

        
        return pd.DataFrame( rows )


        return data
    def closedSectorMoves(self,group):

        return [ move for move in self.moves if ( move["set"]==group and move["sector"]=="closed"  ) ]

    def openSectorMoves(self,group):

        return [ move for move in self.moves if ( move["set"]==group and move["sector"]=="open"  ) ]

    def totProbability(self,group,sector):

        data=self.table()

        data=data[ ( (data["sector"] == sector ) | (data["sector"]=="closed/open") ) & (data["set"]== group) ]

        return np.sum( data["p"])

    def __str__(self):

        return str(self.table() )
    
    def __repr__(self):
        return repr( self.table() )
    
    
    def toJson(self):

        j=[]

        for move in self.moves:
                jMoves=move["move"].toJson()
                n=len(jMoves)
                sectors=[move["sector"] for w in range(n) ]
                p=move["p"]/n
                if (move["sector"] == "closed/open"):
                    sectors[0]="closed"
                    sectors[1]="open"
                    p*=2

                for sector,jMove in zip(sectors,jMoves ):
                    j.append( { "move" : jMove , "weight" : p, "sectors" : [sector], "sets" : [ move["set"] ] } )
        
        return j


def createTableSemiCanonical( C,l,lShort,groups=None,uniform=True,delta=1, restriction=None ):
    if groups is None:
        groups=[0,1]
    
    tab=tableMove( )
    setA,setB=groups

    if not uniform:
        raise NotImplementedError("not uniform sampling of initial bead")
    CAB=C[2]
    for group in [setA,setB]:
        tab.addMove(   levy(l=l),p=0.8,sector="closed",group=group)
        tab.addMove(   translate(delta=delta),p=0.05,sector="closed",group=group)

        tab.addMove(   levy(l=lShort),p=0.1,sector="open",group=group)
        tab.addMove(   translate(delta=delta),p=0.05,sector="open",group=group)
        tab.addMove(   moveHead( l=lShort ),p=0.05,sector="open",group=group)
        tab.addMove(   moveTail( l=lShort ),p=0.05,sector="open",group=group)
        tab.addMove(   swap( l=l ),p=0.1,sector="open",group=group)
        tab.addMove(   advanceRecedeHeadTail( l=lShort, setA=group,setB=(group + 1)%2 ,restriction=restriction ) ,p=0.5,sector="open",group=group)

        tab.addMove(   semiOpenClose(l=lShort,setA=group,setB=(group + 1)%2,C=C[group]),p=0.05 ,sector="closed/open" , group=group),
        tab.addMove(   fullOpenClose(l=lShort,setA=group,setB=(group + 1)%2,C=CAB/C[(group + 1)%2] , restriction=restriction),p=0.05 ,sector="closed/open" , group=group),
        tab.addMove(   createDeleteTwoWorms(l=lShort,setA=group,setB=(group + 1)%2,C=[CAB,1] , restriction=restriction ),p=0.05 ,sector="closed/open" , group=group)

    return tab
        

def createTableCanonical(C,l,lShort,groups=None,uniform=True,delta=1):

    if groups is None:
        groups=[ 0 ]
    
    
    tab=tableMove( )
    
    for group in groups:
        tab.addMove( levy(l=l),p=0.8,sector="closed",group=group)
        tab.addMove(  translate(delta=delta),p=0.05,sector="closed",group=group)
        tab.addMove(  levy(l=lShort),p=0.6 ,sector="open",group=group)
        tab.addMove(  translate(delta=delta),p=0.05,sector="open",group=group)
        tab.addMove(  moveHead( l=lShort ),p=0.05,sector="open",group=group)
        tab.addMove(  moveTail( l=lShort ),p=0.05,sector="open",group=group)
        tab.addMove(  swap( l=l ),p=0.1,sector="open",group=group)
        tab.addMove(  semiOpenClose(l=lShort,setA=group,setB=-1,C=C[group]),p=0.15 ,sector="closed/open" , group=group)
    
    return tab



def createTable(C,l,lShort,ensamble="canonical",groups=None,uniform=True,delta=1):

    if ensamble == "semiCanonical":
        return createTableSemiCanonical(C,l,lShort,groups=groups,uniform=uniform,delta=delta)

    tab=tableMove()

    for group in groups:
        tab.addMove("levy","closed",l=l,weight=2,group=group)
        tab.addMove("translate","closed",delta=delta,group=group)

        tab.addMove("levy","open",l=l,group=group)
        tab.addMove("moveTail","open",l=lShort,group=group)
        tab.addMove("moveHead","open",l=lShort,group=group)
        tab.addMove("swap","open",l=lShort,group=group)
        tab.addMove("open/close","closed/open",l=lShort,C=C,group=group)
        tab.addMove("translate","open",delta=delta,group=group)

        if (ensamble == "grandCanonical"):
        
            tab.addMove("advanceHead","open",l=lShort,group=group)
            tab.addMove("recedeHead","open",l=lShort,group=group)

            if uniform:
                firstParticleDistribution="uniform"
            else:
                firstParticleDistribution="gaussian"
            
            tab.addMove("createWorm/deleteWorm","closed/open",l=lShort,C=C,alpha=1,group=group,firstParticleDistribution=firstParticleDistribution)

    
    return tab