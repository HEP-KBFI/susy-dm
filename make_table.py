import tables
import glob
import sys

NMSSMPoint = {
				"h1_mass": tables.Float32Col(),
				"chi1_mass": tables.Float32Col(),
}
NProb = 54

for i in range(1,NProb):
	NMSSMPoint["PROB%d"%i] = tables.BoolCol()

if __name__=="__main__":
	h5file = tables.openFile("test.h5", mode = "w", title = "SUSY parameter space points")
	group = h5file.createGroup("/", 'NMSSM1', 'NMSSM points')
	table = h5file.createTable(group, 'parspace', NMSSMPoint, "parameter space")
	point = table.row
	files = glob.glob(sys.argv[1])
	for fn in files:
		print "Processing file %s" % fn
		f = open(fn, "r")

		for line in f:
			line_data = map(float, line.split())
			if not len(line_data)==80:
				continue

			point["h1_mass"] = line_data[11]
			point["chi1_mass"] = line_data[22]
			for i in range(1,NProb):
				point["PROB%d"%i] = bool(line_data[26+i])
			point.append()
		h5file.close()
		f.close()
	print "All done"