import numpy as np 
import scipy as sp
import pandas as pd
import os 
import json
from pathlib import Path
import re
import tqdm
import re

def loadScalarData( inputFile , label , minIteration=0 , maxIteration=None ):
    with open (inputFile) as f:
        j=json.load(f)
    path=Path(inputFile)
    simDir=path.parent
    dataFile=os.path.join(simDir,label + ".dat")
    data=pd.read_csv(dataFile,delim_whitespace=True,names=["iteration",label])
    data["iteration"]=data.index
    data["folder"]=simDir
    dataFilter=data["iteration"]>= minIteration
    if (maxIteration is not None):
        dataFilter=dataFilter &  (data["iteration"]<= maxIteration )
    
    data=data[dataFilter]
    return data



def getJson(j,queryString):
    '''
    Use the queryString to obtain an unique value from the json object j
    queryString format : name[index]/.........
    '''    
    if queryString=="":
        return j
    pattern="([a-zA-Z0-9]+)(?:\[(\d+)\])?(?:/(.+))?"
    m=re.match(pattern,queryString)
    if m is not None:
        key=m[1]
        index=m[2]
        remainder=m[3]
        
        result=j[key]
        
        if index is not None:
            result=result[int(index)]
        if remainder is not None:    
            return getJson(result,remainder)
        else:
            return result


def setJson(j,queryString,value):
    '''
    Use the queryString to set an unique value from the json object j
    queryString format : name[index]/.........
    TODO : use mangoDB sintax ?
    '''
    if queryString=="":
        return j
    pattern="([a-zA-Z0-9]+)(?:\[(\d+)\])?(?:/(.+))?"
    m=re.match(pattern,queryString)
    key=m[1]
    index=m[2]
    remainder=m[3]
    
    
    if m is not None:
        
        if remainder is None:
            
            if index is None:
                j[key]=value
            else:
                j[key][int(index)]=value
        else:
                result=j[key]
        
                if index is not None:
                    result=result[int(index)]
                setJson(result,remainder,value)






class jSonParameter:
    def __init__(self,query,label=None,dtype=None):
        ''''
        accept a string wich determines the query .
        Format: "name[/name...]"
        '''
        if label is None:
            label=query
        
        self.query=query
        self.label=label
        self.dtype=dtype
    
    def __getitem__(self, j): # j : json object file, not index
      return getJson(j,self.query)
    
    def __setitem__(self,j,value): 

        setJson(j,self.query,value)
    
    


class systemParameters:
    def __init__(self,parameters=None):
        self._recordedParameters={}

        if parameters is not None:
            for p in parameters:
                self.register(p)
    def register(self,p):
        self._recordedParameters[p.label]=p
    def __getitem__(self,label):
        return self._recordedParameters[label]



class volumeParameter(jSonParameter):

    def __init__(self):
        super(volumeParameter,self).__init__( "lBox" ,"volume")
    def __getitem__(self,j):
        lBox= super(volumeParameter,self).__getitem__(j)
        V=1.0
        for l in lBox:
            V*=l
        return V
    def __setitem__(self,j,V):
        raise NotImplementedError("Cannot set the Volume. Is a derived parameter")

    


class nParticlesParameter(jSonParameter):
    def __init__(self):
        super(nParticlesParameter,self).__init__("particles" ,label="N")
    def __getitem__(self,j):
        ns= super(nParticlesParameter,self).__getitem__(j)
        return np.sum(np.array(ns).astype(int))
    def __setitem__(self,j,N):
        ns= super(nParticlesParameter,self).__getitem__(j)
        if len(ns)==1:
            super(nParticlesParameter,self).__setitem__(j,[int(N)] )
        else:
            raise NotImplementedError("Cannot set the total number of particles. Is a derived parameter")



class densityParameter:
    def __init__(self):
        self.label="density"

    def __getitem__(self,j):
        V=volumeParameter()
        N=nParticlesParameter()

        return N[j]/V[j]
    
    def __setitem__(self,j,n):
        raise NotImplementedError("Cannot set the total number of particles. Is a derived parameter")






def loadParametersFromFile(filename,parameters,recordFolder=True):
    with open(filename) as f:
        j=json.load(f)
    data=loadParameters(j,parameters)
    
    path=Path(filename)
    data["folder"]=path.parent
    return data






def loadParameters(jData, parameters):

    data={}

    for p in parameters:
        data[p.label]=p[jData]
        
    return pd.DataFrame(data,index=[0])


def queryTableFromFolder(jsonFile, label ,parameters=None,minIteration=0,maxIteration=None,systemParameters=None):
    '''
    Construct a trace table for an observable indecated by label from a jSon file
    Returns:
    pandas dataframe with the data
    '''
    if parameters is None:
        parameters=[]
    # parameters which are string are convert to the known supplied parameters
    if systemParameters is not None:
        parameters=[  p if not isinstance(p, str) else systemParameters[p]   for p in parameters      ]


    joinColumn="folder"
    data=loadScalarData(jsonFile,label=label,minIteration=minIteration,maxIteration=maxIteration)
    parameters=loadParametersFromFile(jsonFile,parameters)
    data=pd.merge(data,parameters,on=joinColumn,how="inner")
    data=data.drop(joinColumn,axis=1)

    return data

def queryTable(jsonFiles, *args,**kwds):
    '''
    query the data containted in a list of json files
    Returns a pandas dataframe
    '''
    datas=[]
    for file in tqdm.tqdm(jsonFiles):
        datas.append(queryTableFromFolder(file,*args,**kwds))

    data=pd.concat(datas)
    return data



def scanJsonFiles(folder,maxLevel=1,_currentLevel=0):
    '''
    Iterator on all json file in a certain folder up to maxLevel.
    _currentLevel internal argument used in recursive calls , should not be set by the user
    '''
    if not (  (maxLevel is not None) and (_currentLevel > maxLevel) ):        
        for entry in os.scandir(folder):
            if os.path.isdir(entry):
                yield from scanJsonFiles(entry,maxLevel=maxLevel,_currentLevel=_currentLevel+1)
            if os.path.isfile(entry):
                
                m=re.match(".*\.json",entry.name)
                if m is not None:
                    yield entry.path


def defaultPimcParameters():
    nBeads=jSonParameter("nBeads","nBeads")
    volume=volumeParameter()

    pimcParameters=systemParameters([nBeads,volume,nParticlesParameter(),densityParameter()])

    return pimcParameters




