from pytest import param
from . import simulation
from . import action
from . import observable
from . import moves
import argparse
import pandas as pd
import numpy as np
from .theory import free
from . import inputFileTools
import os


def createSim( a,ratio, boxSize, N, T, C , nBeads,nBlocks=10000,maxTime=None):
    ensamble=simulation.canonicalEnsamble(N=N,boxSize=boxSize,T=T)
    run=simulation.run(nBlocks=nBlocks,stepsPerBlock=1000,correlationSteps=np.sum(N),maxTime=maxTime )

    run.saveConfigurations=False
    a12=a*ratio
    actions=[]
    if a!=0:
        actions=[action.caoBerneAction(a=a,groups=[0,0] ),action.caoBerneAction(a=a,groups=[1,1]),action.caoBerneAction(a=a12,groups=[0,1] ) ]
    

    observables=[observable.thermodynamicEnergy(label="energy") , observable.virialEnergy(label="eV") , observable.superfluidFraction(groups=[0,1],label="rho") ]


    model=simulation.model( ensamble=ensamble,actions=actions,nBeads=nBeads)

    tab= moves.createTableCanonical(C=C,l=int(0.6*nBeads),lShort=int(0.3*nBeads),groups=[0,1],delta=0.3*boxSize[0])

    sim=simulation.simulation(model=model,run=run, observables=observables,moves=tab,N0=N )

    return(sim)


def generateInputFiles(data):
    js=[]
    for i,row in data.iterrows():
        N=int(row["N"])
        P=float(row["P"])
        T=float(row["T"])
        na3=float(row["na3"])
        C=[ float(row["CA"]),float(row["CB"]) ]
        nBeads=int(row["nBeads"])
        ratio=int(row["ratio"])
        nBlocks=int(row["nBlocks"])
        
        if "maxTime" in data.columns:
            maxTime=data["maxTime"]
        

        Na= int(N/2 * ( 1 + P ) )
        Nb= N - Na

        landaC=np.sqrt(2*np.pi)
        L=( (N/2) * landaC**3/free.G(3/2,1) )**(1/3)
        n=N/L**3
        a=(na3/n)**(1/3)

        sim=createSim(a=a,N=[Na,Nb],T=T,boxSize=[L,L,L] ,C=C,nBeads=nBeads,ratio=ratio,nBlocks=nBlocks,maxTime=maxTime)

        js.append(sim.toJson())

    return(js)



def generateLabels(data):
    labels=[]
    for i,row in data.iterrows():
        values=[ value for value in row.values ]

        formatStrings=["{}{{}}".format( key) for key in row.keys() ]   
        formatString="_".join(formatStrings)
        #return values
        labels.append(formatString.format(*values) )
    return (labels)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Generate list of json files')
    parser.add_argument('parameters_file', type=str )

    parser.add_argument('--folder', type=str , default="." )

    args = parser.parse_args()
        
    parameters=args.parameters_file

    data=pd.read_csv(parameters,delim_whitespace=True)

    print(data)

    js=generateInputFiles(data)

    labels=generateLabels(data)


    settings=[ {"folder" : os.path.join(args.folder,label) , "jSon" : [ ["input.json", j  ] ] } for label,j in zip(labels,js) ]

    inputFileTools.createSimFolders(settings)






