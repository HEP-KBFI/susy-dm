import tables
import glob
import sys

NMSSMPoint = {
                "h1_mass": tables.Float32Col(),
                "chi1_mass": tables.Float32Col(),
                "PROB": tables.BoolCol(),
}
NProb = 54

for i in range(1,NProb):
    NMSSMPoint["PROB%d"%i] = tables.BoolCol()

if __name__=="__main__":
    h5file = tables.openFile("nmssm1.h5", mode = "w", title = "SUSY parameter space points")
    group = h5file.createGroup("/", 'NMSSM1', 'NMSSM points')
    table = h5file.createTable(group, 'parspace', NMSSMPoint, "parameter space")
    point = table.row
    files = glob.glob(sys.argv[1])
    print files
    for fn in files:
        print "Processing file %s" % fn
        f = open(fn, "r")

        for line in f:
            try:
                line_data = map(float, line.split())
            except ValueError, e:
                print "bad line in %s" % fn
                continue

            if not len(line_data)==80:
                continue

            point["h1_mass"] = line_data[11]
            point["chi1_mass"] = line_data[22]
            prob = False
            for i in range(1,NProb):
                p = bool(line_data[26+i])
                point["PROB%d"%i] = p
                if p: prob=True
            point["PROB"] = prob #Set to True if any of PROB1-53 is True


            point.append()
        f.close()
    h5file.close()
    print "All done"
