import tables
import glob
import sys

if __name__=="__main__":
	of = tables.openFile("out.h5", "w")
	infiles = map(tables.openFile, glob.glob(sys.argv[1]))
	print "Copying first table"
	infiles[0].root.parspace.copy(of.root)
	ot = of.root.parspace
	print "Output table has %d rows" % ot.nrows
	colname = of.root.parspace.colnames[0]
	dummySel = colname +"=="+colname

	n = 1
	for inf in infiles[1:]:
		t = inf.root.parspace
		print "Putting %d rows from %s (%d/%d)" % (t.nrows, str(inf.filename), n, len(infiles))
		t.whereAppend(ot, dummySel)
		inf.close()
		of.flush()
		print "Output table has %d rows" % ot.nrows
		n += 1
	of.close()