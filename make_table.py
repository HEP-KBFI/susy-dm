import tables
import glob
import sys
import bz2
import time

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
                "PROB": tables.BoolCol(),
                "IFAIL": tables.Int32Col(),
}
NProb = 54

for i in range(1,NProb):
    NMSSMPoint["PROB%d"%i] = tables.BoolCol()


def process_line(point, line):
    try:
        line_data = map(float, line.split())
    except ValueError:
        print "bad line in %s" % fn
        return

    if not len(line_data)==81:
        return

    point["Lambda"] = line_data[0]
    point["Kappa"] = line_data[1]
    point["tanbeta"] = line_data[2]
    point["mu"] = line_data[3]
    point["Alambda"] = line_data[4]
    point["Akappa"] = line_data[5]
    point["M1"] = line_data[6]
    point["M2"] = line_data[7]
    point["M3"] = line_data[8]
    point["MA"] = line_data[9]
    point["MP"] = line_data[10]
    point["h1_mass"] = line_data[11]
    point["h2_mass"] = line_data[13]
    point["h3_mass"] = line_data[15]
    point["chi1_mass"] = line_data[22]
    point["omg"] = line_data[26]
    point["IFAIL"] = int(line_data[80])

    prob = False
    for i in range(1,NProb):
        p = bool(line_data[26+i])
        point["PROB%d"%i] = p
        if p: prob=True
    point["PROB"] = prob #Set to True if any of PROB1-53 is True
    point.append()
    return

if __name__=="__main__":
    h5file = tables.openFile("nmssm1.h5", mode = "w", title = "SUSY parameter space points")
    group = h5file.createGroup("/", 'NMSSM1', 'NMSSM points')
    filters = tables.Filters(complevel=9, complib='blosc', fletcher32=False)
    table = h5file.createTable(group, 'parspace', NMSSMPoint, "parameter space", expectedrows=int(sys.argv[2]))
    point = table.row
    files = glob.glob(sys.argv[1])
    print files
    for fn in files:
        t0 = time.time()
        print "Processing file %s" % fn
        f = bz2.BZ2File(fn, "r")
        lines = True
        while lines:
            lines = f.readlines(1000000)
            for line in lines:
                process_line(point, line)
        f.close()
        t1 = time.time()
        elapsed = t1-t0
        print "Processing took %d seconds" % elapsed
    h5file.close()
    print "All done"
