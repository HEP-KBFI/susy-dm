import tables
import glob
import sys
import bz2
import time
from multiprocessing import Pool


outdir = "/Users/joosep/Desktop/test/"

NMSSMPoint = {
                "Lambda": tables.Float32Col(), #0
                "Kappa": tables.Float32Col(), #1
                "tanbeta": tables.Float32Col(), #2
                "mu": tables.Float32Col(), #3
                "Alambda": tables.Float32Col(), #4
                "Akappa": tables.Float32Col(), #5
                "M1": tables.Float32Col(), #6
                "M2": tables.Float32Col(), #7
                "M3": tables.Float32Col(), #8
                "MA": tables.Float32Col(), #9
                "MP": tables.Float32Col(), #10
                "h1_mass": tables.Float32Col(), #11
                "h2_mass": tables.Float32Col(), #13
                "h3_mass": tables.Float32Col(), #15
                "chi1_mass": tables.Float32Col(), #22
                "omg": tables.Float32Col(), #26
                #"PROB": tables.BoolCol(),
                "IFAIL": tables.Int32Col(),
                "PROB": tables.Int64Col(),
}
NProb = 54

format_nmhdecay = dict()
format_nmhdecay["Lambda"] = 0
format_nmhdecay["Kappa"] = 1
format_nmhdecay["tanbeta"] = 2
format_nmhdecay["mu"] = 3
format_nmhdecay["Alambda"] = 4
format_nmhdecay["Akappa"] = 5
format_nmhdecay["M1"] = 6
format_nmhdecay["M2"] = 7
format_nmhdecay["M3"] = 8
format_nmhdecay["MA"] = 9
format_nmhdecay["MP"] = 10
format_nmhdecay["h1_mass"] = 11
format_nmhdecay["h2_mass"] = 13
format_nmhdecay["h3_mass"] = 15
format_nmhdecay["chi1_mass"] = 22
format_nmhdecay["omg"] = 26
for i in range(1,NProb):
    format_nmhdecay["PROB%d"%i] = 26+i

format_nmspec = dict()
format_nmspec["IFAIL"]=78

#for i in range(1,NProb):
#    NMSSMPoint["PROB%d"%i] = tables.BoolCol()

def bool2int(x):
    y = 0L
    for i,j in enumerate(x):
        if j: y += 1<<i
    return y

def process_line(point, line, form=format_nmhdecay):
    try:
        line_data = map(float, line.split())
    except ValueError:
        print "malformed line: %s" % line
        return

#    if not len(line_data)==81:
#        return

    try:
        for (k, v) in form.items():
            point[k] = line_data[v]
    except:
        print "malformed line: %s" % line
        return

    point["IFAIL"] = int(point["IFAIL"])
    #point["Lambda"] = line_data[0]
    #point["Kappa"] = line_data[1]
    #point["tanbeta"] = line_data[2]
    #point["mu"] = line_data[3]
    #point["Alambda"] = line_data[4]
    #point["Akappa"] = line_data[5]
    #point["M1"] = line_data[6]
    #point["M2"] = line_data[7]
    #point["M3"] = line_data[8]
    #point["MA"] = line_data[9]
    #point["MP"] = line_data[10]
    #point["h1_mass"] = line_data[11]
    #point["h2_mass"] = line_data[13]
    #point["h3_mass"] = line_data[15]
    #point["chi1_mass"] = line_data[22]
    #point["omg"] = line_data[26]
    #point["IFAIL"] = int(line_data[80])
    #probs = [bool(line_data[26+i]) for i in range(1,NProb)]
    probs = list()
    for i in range(1,NProb):
        p = "PROB%d" % i
        if p in point:
            probs.append(point["PROB%d" % i])
    if len(probs)==NProb-1:
        point["PROB"] = bool2int(probs)

    point.append()
    return

def process_file(infn):
    ofn = outdir+infn.replace(".dat.bz2", ".h5").split("/")[-1]
    print "begin processing %s -> %s" % (infn, ofn)
    h5file = tables.openFile(ofn, mode = "w")
    filters = tables.Filters(complevel=9, complib='blosc', fletcher32=False)
    table = h5file.createTable("/", 'parspace', NMSSMPoint, "parameter space", expectedrows=2000000, filters=filters)
    point = table.row
    t0 = time.time()
    f = bz2.BZ2File(infn, "r")
    lines = True
    while lines:
        lines = f.readlines(100000000)
        map(lambda x: process_line(point, x, form=format_nmspec), lines)
    f.close()
    t1 = time.time()
    elapsed = t1-t0
    h5file.flush()
    h5file.close()
    return [ofn, elapsed]

if __name__=="__main__":
    files = glob.glob(sys.argv[1])
    outdir = sys.argv[2]
    if outdir[-1]!="/": outdir += "/"

    print "Input files: %d" % len(files)
    for f in files:
        print f
    print 80*"-"

    ncores = int(sys.argv[3])
    if ncores>1:
        p = Pool(int(sys.argv[3]))
        res = p.map(process_file, files)
    else:
        res = map(process_file, files)
    print res
    print "All done"
