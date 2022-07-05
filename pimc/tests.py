import argparse
from pimc import *
import os
import subprocess
import json
import tqdm

class testSemiOpenCloseMove:

    def __init__(self,pimc_test):
        self.testEx=pimc_test
        self.subDirectoryBase="."    
        self.inputFileName="input.json"
        
    def runQMC(self,settings, name):
        if not os.path.exists(self.subDirectoryBase):
            os.makedirs(self.subDirectoryBase)

        subDirectory=os.path.join(self.subDirectoryBase, name)

        if not os.path.exists(subDirectory):
            os.makedirs(subDirectory)

        # create the input file
        inputFilePath=os.path.join(subDirectory,self.inputFileName)
        
        with open(inputFilePath, "w") as f:
            json.dump(settings,f)
        
        logFileName=os.path.join(subDirectory,"run.log")
        with open(logFileName,"w") as logFile:
            result = subprocess.run([ self.testEx, "--gtest_filter=configurationsTest.semiOpenCloseMove", self.inputFileName],cwd=subDirectory,stdout=logFile)
        return result


    def run(self):

        def generator():


            yield self.runQMC({"C": 1e-1,"l":4,"t0":8,"nBlocks": 3000 , "stepsPerBlock" : 10000 } , name="semiOpenCloseMove-1")
            yield self.runQMC({"C": 1e-1,"l":6,"t0":2,"nBlocks": 3000 , "stepsPerBlock" : 10000 } , name="semiOpenCloseMove-2")
            yield self.runQMC({"C": 1e-1,"l":2,"t0":2,"nBlocks": 3000 , "stepsPerBlock" : 10000 } , name="semiOpenCloseMove-3")
            yield self.runQMC({"C": 1e-1,"l":9,"t0":9,"nBlocks": 10000 , "stepsPerBlock" : 10000 } , name="semiOpenCloseMove-4")
        
        for result in tqdm.tqdm( generator() , total=4):
            pass





if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run some example test for pimc')
    parser.add_argument('pimc_test', type=str )

    args = parser.parse_args()

    currentTest=testSemiOpenCloseMove(args.pimc_test)
    currentTest.runQMC({"C": 1e-1,"l":4,"t0":8,"nBlocks": 10000 , "stepsPerBlock" : 100000 } , name="semiOpenCloseMove-1")
