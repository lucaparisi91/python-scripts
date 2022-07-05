from functools import reduce
import re
import pandas as pd
import numpy as np
import re
import copy
import json
import itertools
import os
from . import moves



    

        

def getJson(j,queryString):
    
    if queryString=="":
        return j
    pattern="([a-zA-Z0-9]+)(?:\[(\d+)\])?(?:/(.+))?"
    m=re.match(pattern,queryString)
    print(m[3])
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



def expandJson(jTtemplate, **kwds):
    '''
    jTemplate : template json input file
    Takes in input a series of kwd , valude with value a vector and kwd the value to update
    TODO : mondodb queries 
    '''
    js=[]
    #alues=[ v.flatten() for v in np.meshgrid(*kwds.values()) ]
    for currentValues in itertools.product(*kwds.values() ):
        j=copy.deepcopy(jTemplate)
        for kwd,value in zip( kwds.keys(), currentValues):
            
            setJson(j,kwd,value)
        js.append(j)
    return js

def createInputFilesFromTemplate(template,params):
    keys=list(params.keys())
    N=len(params[keys[0]])
    js=[]
    
    for i in range(N):
        
        j=copy.deepcopy(template)
        for key in keys:
            value=params[key][i]
            setJson(j,str(key),value)
        js.append(j)
    return js

class numpyCompatibleEncoder(json.JSONEncoder):
    # Handles the default behavior of
    # the encoder when it parses an object 'obj'
    def default(self, obj):
        # If the object is a numpy array
        if isinstance(obj, np.ndarray):
            # Convert to Python List
            return obj.tolist()
        else:
            if isinstance(obj, np.int64) or isinstance(obj, np.int32):
                return int(obj)
            
            # Let the base class Encoder handle the object
            return json.JSONEncoder.default(self, obj)


        
def createSimFolders(settings):
    '''
    vector of json file of format
    "folder" : name of the folder to create
    "jSon" : vector of {name : content}
    '''
    
    for j in settings:
        folder=j["folder"]
        if not os.path.exists(folder):
            os.makedirs(folder)
        for jSonInput in j["jSon"]:
            filename=os.path.join( folder, jSonInput[0] )
            content=jSonInput[1]
            with open(filename,"w") as f:
                json.dump(content,f,indent=4,cls=numpyCompatibleEncoder)

def cartesianProduct(kwds):
    rows=[]
    for currentValues in itertools.product(*kwds.values() ):
          rows.append(np.array(currentValues)) 
    return pd.DataFrame(rows,columns=kwds.keys())
def setMoveParameters(j):
    for move in j["movesTable"]:
        if "reconstructionMaxLength" in move["move"].keys():
            
            move["move"]["reconstructionMaxLength"]=max( int(0.2 * int(j["nBeads"]) ) , 1 )
            
            
            if move["move"]["kind"]=="open" or move["move"]["kind"]=="close":
                move["move"]["C"]=1
                move["move"]["reconstructionMaxLength"]=1
            if move["move"]["kind"] != "levy":
                move["move"]["reconstructionMaxLength"]= max(1, move["move"]["reconstructionMaxLength"]//3) 
    return j



