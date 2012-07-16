import sys
import imp
random_sampling = imp.load_source("random_sampling", "/home/joosep/susy/random_sampling.py")

import subprocess

N = 100
modelname = "NMSSM"

def call(varlist):
    cmd = ["./main"] + varlist
    #print "calling %s" % cmd
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = proc.stdout.readlines()
    return out

if __name__=="__main__":

    conf = random_sampling.read_parameters("./sampling_NMSSM.cfg")

    for i in range(100):
        s = random_sampling.sample_nmssm(conf)
        out = call( map(lambda x: str(x), [s["m0"], s["mhf"], s["a0"], s["tb"], s["Lambda"], s["aKappa"] ]))
        if out[-1][:-1] == "success":
            print s

