import tables
import glob
import sys

if __name__=="__main__":
	of = tables.openFile("out.h5", "w")
	infiles = glob.glob(sys.argv[1])
	infiles[0].parspace.copy(of.root)
	for infn in infiles[1:]:
		inf = tables.openFile(infn)
		inf.NMSSM1.parspace.whereAppend(of.root.parspace, of.parspace.colnames[0] + "==" + of.parspace.colnames[0])
		inf.close()
		of.flush()
	of.close()