import tables
import glob
import sys

NMSSMPoint = {
                "lambda": tables.Float32Col(), #0
                "kappa": tables.Float32Col(), #1
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
                "chi1_mass": tables.Float32Col(), #22
                "omg": tables.Float32Col(), #26
                "PROB": tables.BoolCol(),
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

    if not len(line_data)==80:
        return

    point["lambda"] = line_data[0]
    point["kappa"] = line_data[1]
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
    point["chi1_mass"] = line_data[22]
    point["omg"] = line_data[26]
    
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
        print "Processing file %s" % fn
        f = open(fn, "r")
        lines = True
        while lines:
            lines = f.readlines(sizehint=10000000)
            for line in lines:
                process_line(point, line)
        f.close()
    h5file.close()
    print "All done"
