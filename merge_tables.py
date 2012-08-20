import tables
import glob
import sys

if __name__=="__main__":
	of = tables.openFile("out.h5", "w")
	infiles = map(tables.openFile, glob.glob(sys.argv[1]))
	infiles[0].root.parspace.copy(of.root)
	for inf in infiles[1:]:
		inf.root.parspace.whereAppend(of.root.parspace, of.parspace.colnames[0] + "==" + of.parspace.colnames[0])
		inf.close()
		of.flush()
	of.close()