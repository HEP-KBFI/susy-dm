import sys
import imp
random_sampling = imp.load_source("random_sampling", "/home/joosep/susy/random_sampling.py")

import subprocess

N = 100
modelname = "NMSSM"

def call(infile):
    cmd = ["./main", infile]
    #print "calling %s" % cmd
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    ret = proc.wait()
    out = proc.stdout.readlines()
    return (out, ret)

def prepare_input_file(parameters, of="par.dat"):
    f = open(of, "w")
    for (name, val) in parameters.items():
        f.write("%s %f\n" % (name, val))
    f.close()

if __name__=="__main__":

    conf = random_sampling.read_parameters("./sampling_NMSSM_EWSB.cfg")

    good_points = []
    for i in range(100):
        s = random_sampling.sample(conf)
        prepare_input_file(s)
        (out, retcode) = call("par.dat")
        print "%d %d %s" % (i, retcode, out[-1][:-1])
        if retcode==0:
            good_points.append(s)
    if len(good_points)>0:
        keys = good_points[0].keys()
        for p in good_points:
            s = ""
            for k in keys:
                s += "%s: %.2f " % (k, p[k])
            print s
