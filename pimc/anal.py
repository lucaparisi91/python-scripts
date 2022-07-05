import numpy as np
from numpy.lib.function_base import average
import scipy as sp
import pandas as pd
import copy
import arviz as az
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import colorsys
from scipy.optimize import curve_fit
import inspect
import colorsys
import tqdm


defaultTemplate="plotly_dark"


def stringHexToRGB(color):
    h = color.lstrip("#")
    return tuple(   int(h[i:i+2], 16)/255 for i in (0, 2, 4))


def meanLineColor(color):
    color=stringHexToRGB(color)
    hlsColor=list(colorsys.rgb_to_hls(*np.array(color)))
    hlsColor[1]=hlsColor[1]**2

    rgbColor=tuple(np.ceil(np.array(colorsys.hls_to_rgb(*hlsColor))*255).astype(int))
    return '#{:02X}{:02X}{:02X}'.format(*rgbColor)

def flattenMultiIndices(agg):
    df=agg.reset_index()
    df.columns=["_".join(filter(lambda x: x!="" ,list(index))) for index in agg.reset_index().columns ]
    return df

def aggregate(data,columns,parameters=None,estimates=["mean","error","correlationLength"]):

    dropIndex=False
    if parameters is None:
        data=data.assign(__dummy_key=1)
        parameters=["__dummy_key"]
        dropIndex=True
    def ess(x) :
        return az.ess(np.array(x))
    groups=data.groupby(parameters)
    
    def sampleSize(x):
        return len(x)
    def var(x):
        return np.var(x)
    
    agg=groups.agg({column : [np.mean,ess,sampleSize,var] for column in columns} )
    for level in agg.columns.levels[0]:
        agg1=agg[level].assign(error=lambda x : np.sqrt(x["var"]/x["ess"]))
        agg.loc[:,(level,"error")]=agg1["error"]
    
    if dropIndex:
        agg=agg.reset_index(drop=True)


    
    return flattenMultiIndices(agg)





def addHLine(fig,rangeX,a , col=1,row=1,name=None,*args,**kwds):
    
    fig.add_trace( go.Scattergl(name=name,
                x=np.array(rangeX),y=np.zeros(2)+ a  ,
                mode = "lines", *args , **kwds
                ) ,row=row, col=col )


def addVLine(fig,x,yRange , col=1,row=1,name=None,*args,**kwds):
    
    addLine(fig,np.array([x,x]),np.array(yRange),name=name,*args,**kwds)


def addLine(fig,x,y , col=1,row=1,name=None,*args ,**kwds):
    
    fig.add_trace( go.Scattergl(name=name,
                x=x,y=y  ,
                mode = "lines",
                *args, **kwds
                ) ,row=row, col=col )




def runningMean(x):
    return [ np.mean(x[0:i+1])    for i in range( len(x)) ]


def runningError(x,corrLength=None):
    if corrLength is None:
        corrLength=len(x)/az.ess(np.array(x))
    
    return [ np.sqrt(np.var(x[0:i+1])/( ( i+1 )/corrLength))    for i in range( len(x)) ]


def runningStatistics(data,columns,corrLength=None):
    runningStats={}
    for column in columns:
        runningStats[(column,"runningError")]=runningError(data[column],corrLength=corrLength)
        runningStats[(column,"runningMean")]=runningMean(data[column])
        
    return pd.DataFrame(runningStats)


def plotTimeSeries(data,observables,x=None,parameters=None,width=1000,height=500,showMean=True,runningMean=True,alpha=0.8,overlap=False,corrLength=None,sampleFraction=None):
    noKeys=False
    if parameters is None:
        parameters=["__dummy_key"]
        data=data.assign(__dummy_key=0)
        noKeys=True        

    groups=data.groupby(parameters)
    fig=go.Figure()

    nPlots=len(groups)

    if overlap:
        nPlots=1
    else:
        nPlots=len(groups)


    fig = make_subplots(rows=nPlots, cols=1)

    def sampleData(df):
        if sampleFraction is not None:
            return df.sample(frac=sampleFraction).sort_index()
        else:
            return df

    for i, (key,df) in tqdm.tqdm(enumerate(groups)):
        
        df=sampleData(df)
        agg=aggregate(df,observables)
        color="orange"
        runStats=runningStatistics(df,observables,corrLength=corrLength)

        
        iPlot = i % nPlots

        if x is None:
            iterations=df.index
        else:
            iterations=df[x]
        
        
        keys=[key]
        for iOb,observable in enumerate(observables):
            color=px.colors.qualitative.Plotly[ ( iOb  + i*len(observables) )% len(px.colors.qualitative.Plotly)]
            label=observable + " "

            if not noKeys:
                label += ",".join( [ "{}={}".format(parameter,key) for parameter,key in zip(parameters,keys) ] )


            fig.add_trace( go.Scattergl(name=label,
                x=iterations,y=df[observable]  ,
                mode = "markers",
                opacity=alpha,
                line = {"color":color,"dash":"dash"}) ,row=iPlot+1, col=1 )
            
            if showMean and not runningMean:
                fig.add_trace( go.Scattergl(
                x=iterations,y=df.index*0 + float(agg[observable]["mean"])  , mode = "lines",showlegend=False,line_color=color,fill=None) ,row=iPlot+1, col=1 )
                
            if runningMean:
               
                
                

                xShaded=np.concatenate([iterations,np.flip(iterations)])
                yShaded=np.concatenate( [ 
                    runStats[observable]["runningMean"] + runStats[observable]["runningError"],
                    np.flip(runStats[observable]["runningMean"] - runStats[observable]["runningError"])
                ]
                )


                fig.add_trace( go.Scattergl(
                                x=xShaded,y=yShaded   , mode = "lines",line_color=meanLineColor(color),opacity=0.8,showlegend=False,fill="toself") ,
                                row=iPlot+1, col=1 )
                
            
                fig.add_trace( go.Scattergl(
                x=iterations,y=runStats[observable]["runningMean"]  ,
                    mode = "lines",fill=None,line_color=meanLineColor(color),opacity=1,showlegend=False) ,row=iPlot+1, col=1 )
                
            
                
            
                
            
        
    fig.update_layout(
        autosize=False,
        template=defaultTemplate,
        height=height,
        width=width)
    return fig



def scatterPlot(data,x,y,delta=None,parameters=None,width=1000,height=500,alpha=0.8,overlap=False): 

    noKeys=False
    if parameters is None:
        parameters=["__dummy_key"]
        data=data.assign(__dummy_key=0)
        noKeys=True        

    groups=data.groupby(parameters)
    fig=go.Figure()

    nPlots=len(groups)

    if overlap:
        nPlots=1
    else:
        nPlots=len(groups)


    fig = make_subplots(rows=nPlots, cols=1)


    for i, (key,df) in enumerate(groups):
        keys=[key]
        label=y
        if not noKeys:
            label += ",".join( [ "{}={}".format(parameter,key) for parameter,key in zip(parameters,keys) ] )

        color=px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)]
        iPlot = i % nPlots

        if delta is not None:
            error_y={
            "type" : 'data', # value of error bar given in data coordinates
            "array" : df[delta],
            "visible" : True
            }
        else:
            error_y={}

        fig.add_trace( go.Scattergl(name=label,
                x=df[x],y=df[y]  ,
                mode = "markers",
                opacity=alpha,
                line = {"color":color,"dash":"dash"},
                error_y= error_y,
                ) ,row=iPlot+1, col=1 )


    fig.update_layout(
        autosize=False,
        template=defaultTemplate,
        height=height,
        width=width)

    return fig



def blockedError(data, nBlocks):
    blockedPopulations=list(map( np.mean, np.array_split( data  ,nBlocks) ))
    error=np.sqrt(np.var(blockedPopulations)/len(blockedPopulations) )
    return error

def blockedErrors(data,nBlockings=10, minSampleSize=10, maxSampleSize=None):
    if maxSampleSize is None:
        maxSampleSize=len(data)//10
    else:
        maxSampleSize=min(maxSampleSize,len(data)//2)
    
    blockedSizes=range(minSampleSize,maxSampleSize , (maxSampleSize - minSampleSize)//nBlockings )
    
    blockings=len(data)/np.array(blockedSizes) 
    errors= [ blockedError(data,nBlocks)  for nBlocks in blockings  ]
    return len(data)/np.array(blockings),errors



def concatenateData(f):
    '''
    wrap a function f(data, ...) returning a dataset.
    '''

    def wrapped(data,*args,parameters=None,**kwds):
        noKeys=False
        if parameters is None:
            parameters=["__dummy_key"]
            data=data.assign(__dummy_key=0)
            noKeys=True        
        
        groups=data.groupby(parameters)

        results=[]
        for i, (key,df) in tqdm.tqdm(enumerate(groups)):
            try:
                result=f(df,*args,**kwds)
            except TypeError as e:
                
                print ("Warning. Error for key=  " + str(key) )
                print(str(e) )
                continue
        
            if (not noKeys) :
                
                if not ( hasattr(key, '__iter__') ) or isinstance( key, str)  :
                    keys=[key]

                else:
                    keys=key
                
                for value,name in zip(keys,parameters):
                    result[name]=value
            results.append(result)
        return pd.concat(results).reset_index(drop=True)
    return wrapped




def concatenateFigures(f):
    '''
    wrap a function f(fig,data, ...) returning a figure.
    '''

    def wrapped(fig,data,*args,parameters=None,**kwds):
        noKeys=False
        if parameters is None:
            parameters=["__dummy_key"]
            data=data.assign(__dummy_key=0)
            noKeys=True        
        
        groups=data.groupby(parameters)

        results=[]
        for i, (key,df) in enumerate(groups):
            try:
                newFig=f(fig,df,*args,**kwds)
                if newFig is not None:
                    fig=newFig
            except TypeError as e:
                
                print ("Warning. Error for key=  " + str(key) )
                print(str(e) )
                continue
        
                
            if not ( hasattr(key, '__iter__') ) or isinstance( key, str)  :
                keys=[key]
                
        
        return fig
    return wrapped



@concatenateData
def fitModel(data , f , x , y, delta=None  ):
    if delta is None:
        parameters,cov=curve_fit(f,np.array(data[x]),np.array(data[y])  )
    else:
        parameters,cov=curve_fit(f,np.array(data[x]),np.array(data[y]) ,sigma=np.array(data[delta]) )
    
    errors=np.sqrt(np.diag(cov))
    names=inspect.getfullargspec(f).args[1:]
    
    result={ name : p for p,name in zip(parameters,names) }
    result.update( { name + "_error" : err for err,name in zip(errors,names)   })
    
    
    return pd.DataFrame(result,index=[0])

@concatenateData
def rebatch(data,nSamples):
    
    
    dfs=np.array_split(data.reset_index(drop=True),nSamples)
    for df in dfs:
        df["iteration"]=np.mean(df["iteration"])
    averagedDfs=[ df.mean(axis=0).to_frame().T for df in dfs]
    if len(averagedDfs) !=0:
        result=pd.concat(averagedDfs)
        
        return result
    



@concatenateFigures
def addFittedModel(fig , data , f, name, row =1 , col=1  ,nSamples=1000,*args, **kwds ):
    
    full_fig = fig.full_figure_for_development(warn=False)
    xRange=full_fig.layout.xaxis.range
    
    names=inspect.getfullargspec(f).args[1:]
    
    parameters= np.array([ data[name] for name in names]).astype(float)
    errors= np.array([ data[name+"_error"] for name in names]).astype(float)
    
    
    x=np.linspace(xRange[0],xRange[1],num=nSamples)
    y=np.array(f(x,*parameters))
    
    yp=np.array(f(x,*(parameters + errors)))
    ym=np.array(f(x,*(parameters-errors)))
    
    fig.add_trace( go.Scattergl(name=name,
                x=x,y=y  ,
                mode = "lines" , *args , **kwds
                ) ,row=row, col=col )
    #fig.add_trace( go.Scattergl(name=name,showlegend=False,
    #            x=x,y=ym  ,
    #            mode = "lines", *args,**kwds
    #            ) ,row=row, col=col )
    #fig.add_trace( go.Scattergl(name=name,showlegend=False,
    #            x=x,y=yp  ,
    #            mode = "lines", *args, **kwds
    #            ) ,row=row, col=col )
    
    xShaded=np.concatenate([x,np.flip(x)])
        
    yShaded=np.concatenate( [ ym,np.flip(yp) ] )
    
    fig.add_trace( go.Scattergl( 
                                x=xShaded,y=yShaded   , mode = "lines",opacity=0.4,
                                showlegend=False,fill="toself",*args,**kwds) ,
                                row=row, col=col 
                )



def setLogScale(fig,axis):
    if axis == "x":
        fig.update_xaxes(type="log")
    else:
        fig.update_yaxes(type="log")
    return fig
