import tables
import glob
import sys

if __name__=="__main__":
	of = tables.openFile("out.h5", "w")
	infiles = map(tables.openFile, glob.glob(sys.argv[1]))
	infiles[0].root.parspace.copy(of.root)
	ot = of.root.parspace
	colname = of.root.parspace.colnames[0]
	dummySel = colname +"=="+colname
	for inf in infiles[1:]:
		t = inf.root.parspace
		print "Putting %d rows from %s" % (t.nrows, str(inf))
		t.whereAppend(ot, dummySel)
		inf.close()
		of.flush()
		print "Output table has %d rows" % ot.nrows
	of.close()