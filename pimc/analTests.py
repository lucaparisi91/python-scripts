# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
from dash_core_components.Dropdown import Dropdown
import dash_html_components as html
import plotly.express as px
import pandas as pd
import os
import anal
import numpy as np
import json
from dash.dependencies import Input, Output


knownRuns = { "semiOpenCloseMove-1" : 
                {"observables" : {"lengthOpen" : 2.95, "lengthClosed":2.6756971068005466, "openRatio":0.5826452658869408 } ,
                    "burnIn" : 200
                } 
            
            }
knownRuns["semiOpenCloseMove-2"]=knownRuns["semiOpenCloseMove-1"]
knownRuns["semiOpenCloseMove-3"]=knownRuns["semiOpenCloseMove-1"]
knownRuns["semiOpenCloseMove-4"]=knownRuns["semiOpenCloseMove-1"]



def findDirectories(folder):
    dirs=[]
    for file in os.listdir(folder):
        d = os.path.join(folder, file)
        print(d)
        if os.path.isdir(d) and file in knownRuns.keys():
            dirs.append(d)
    return dirs





def createObservableReports(dir):
    baseDirName=os.path.basename(dir)
    run=knownRuns[ baseDirName ]
    observables=run["observables"]
    observableReports=[]
    for name,value in observables.items():
        source=os.path.join(dir,name+".dat")

        data=pd.read_csv(source,
                    delim_whitespace=True,
                    names=["iteration","l2"]
                )
        data=data[data["iteration"] > run["burnIn"] ] 

        fig=anal.plotTimeSeries(data.query("iteration>=0") ,["l2"],parameters=None,height=500,width=900,overlap=False,x="iteration",corrLength=None)
        xMax=np.max(data["iteration"])
        anal.addHLine(fig,[0,xMax], value)

        observableReports.append(html.Div(children=[
                html.H3(children=name),
            dcc.Graph(
            figure=fig
            ) ]
        ) )
    return observableReports

def createFolderReport(dir):
    baseDirName=os.path.basename(dir)
    observableReports=createObservableReports(dir)
    return html.Div(children=[
    html.H2(children=baseDirName),
    html.Div(children='''
    Test in the semi canonical opening and closing a chain
'''),
    *observableReports

    ] )
    


app = dash.Dash(__name__)

dirs=findDirectories("/home/luca/source/qmc4/PPA/run/semiOpenCloseMove")

reports=[]

defaultDir=dirs[0]

report=html.Div(createFolderReport(defaultDir) , id="report" )


selectFolder=dcc.Dropdown(
            options=[
            {"label" : dir , "value": dir} for dir in dirs ],
            value=defaultDir,id="folder-select"
        )


@app.callback(
    Output("report",component_property="children"),
    Input('folder-select', component_property='value')
)
def updateFolderReport(dir):
    return createFolderReport(dir)

app.layout=html.Div(children=[selectFolder,report])

if __name__ == '__main__':
    app.run_server(debug=True)


    