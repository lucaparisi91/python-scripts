import pandas as pd
import numpy as np
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate list of json files')
    parser.add_argument('parameters_file', type=str )
    parser.add_argument('out', type=str )
    parser.add_argument('minC', type=float , default="." )
    parser.add_argument('maxC', type=float , default="." )
    parser.add_argument('--scale', type=str , default="log" )
    parser.add_argument('--num', type=int , default=10 )

    args = parser.parse_args()

    CA=[]
    if args.scale == "log":
        CA=np.logspace( args.minC,args.maxC,num=args.num)
    else:
        raise( RuntimeError("Scale not known"))
    
    data=pd.read_csv(args.parameters_file,delim_whitespace=True)
    data=data.drop("CA",axis=1)
    data=pd.merge(data,pd.DataFrame({"CA":CA}),how="cross")
    data.to_csv(args.out,sep=" ")
    