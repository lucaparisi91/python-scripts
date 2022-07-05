import argparse
from pimc import *




if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Gather observables from a folder')
    parser.add_argument('folder', type=str )
    parser.add_argument('label', type=str )
    parser.add_argument('--parameters',action="extend", nargs="+", type=str )



    args = parser.parse_args()
        
    parameters=args.parameters
    

    label = args.label

    folder= args.folder

    
    files=list(scanJsonFiles(folder))


    systemParameters = defaultPimcParameters()

    systemParameters.register( jSonParameter("action[0]/potential/V0",label="V0"))

    data=queryTable(files,label=label,maxIteration=None,minIteration=0,parameters=parameters,systemParameters=systemParameters )


    print( data.to_csv(sep="\t") )