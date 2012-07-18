import sys
import imp
#random_sampling = imp.load_source("random_sampling", "random_sampling.py")
#from . import random_sampling
import random_sampling

import subprocess
import os
import shutil

N = 1000000
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

def save_output(d, output, ofdir="out"):
    if not os.path.exists(ofdir):
        os.makedirs(ofdir)
    shutil.move("par.dat", "%s/par_%d.dat" % (ofdir, d))
    shutil.move("spectr.dat", "%s/spectr_%d.dat" % (ofdir, d))
    shutil.move("decay.dat", "%s/decay_%d.dat" % (ofdir, d))
    of = open("%s/out_%d.dat" % (ofdir, d), "w")
    for l in output:
        of.write(l)
    of.close()
    return

if __name__=="__main__":

    fromFile = False
    saveGood = True

    if not fromFile:
        conf = random_sampling.read_parameters("./sampling_NMSSM_EWSB.cfg")
    else:
        pointFile = open("good.txt")

    good_points = []

    if not fromFile:
        for i in range(N):
            s = random_sampling.sample(conf)
            prepare_input_file(s)
            (out, retcode) = call("par.dat")
            print "%d %d %s" % (i, retcode, out[-1][:-1])
            if retcode==0:
                save_output(i, out)
                good_points.append(s)

    if fromFile:
        i = 0
        for line in pointFile:
            s = {}
            coords = line.split(" ")
            coords = filter(lambda x: "=" in x, coords)
            print coords
            for c in coords:
                coord = c.split("=")
                s[coord[0]] = float(coord[1])
            prepare_input_file(s)

            (out, retcode) = call("par.dat")
            print "%d %d %s" % (i, retcode, out[-1][:-1])
            if retcode==0:
                save_output(i, out)
                good_points.append(s)
            i += 1

    if len(good_points)>0:
        keys = good_points[0].keys()
        if saveGood:
            goodFile = open("good.txt", "w")
        for p in good_points:
            s = ""
            for k in keys:
                if saveGood:
                    goodFile.write("%s=%f " % (k, p[k]))
                s += "%s: %.2f " % (k, p[k])
            if saveGood: goodFile.write("\n")
            print s
        if saveGood:
            goodFile.close()
