import tables
import glob
import sys
import time

if __name__=="__main__":
	of = tables.openFile("out.h5", "w")
	infiles = map(tables.openFile, glob.glob(sys.argv[1]))
	print "Copying first table from file %s" % infiles[0].filename
	infiles[0].root.parspace.copy(of.root)
	ot = of.root.parspace
	print "Output table has %d rows" % ot.nrows
	colname = of.root.parspace.colnames[0]
	dummySel = colname +"=="+colname
	inrows = sum(map(lambda x: x.root.parspace.nrows, infiles))
	n = 1
	for inf in infiles[1:]:
		t0 = time.time()
		t = inf.root.parspace
		print "Putting %d rows from %s (%d/%d)" % (t.nrows, str(inf.filename), n, len(infiles))
		t.whereAppend(ot, dummySel)
		inf.close()
		of.flush()
		print "Output table has %d rows" % ot.nrows
		n += 1
		t1 = time.time()
		print "Copying took %d seconds" % (t1-t0)

	outrows = ot.nrows
	print "Copied %d rows" % outrows
	if inrows != outrows:
		raise Exception("Inrows =/= outrows")
	of.close()