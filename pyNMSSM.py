import string
import subprocess
import os
import multiprocessing
import sys

dev_null = open(os.devnull, "w")
class Input:
	def __init__(self, templateFileName):
		f = open(templateFileName)
		inp_template_txt = f.read()
		f.close()
		self.template = string.Template(inp_template_txt)

	def writeInput(self, values, ofName):
		for (k,v) in values.items():
			values[k] = Input.toFortranString(v)

		t = self.template.substitute(values)
		f = open(ofName, "w")
		f.write(t)
		f.close()

	@staticmethod
	def toFortranString(d):
		if type(d)==float:
			return str(d) + "D0"
		elif type(d)==long:
			return str(d)

class Point:
	def __init__(self, **kwargs):
		self.values = {}
		for (k, v) in kwargs.items():
			self.values[k] = v

	def dict(self):
		return self.values

class Runner:
	def __init__(self):
		pass

	def call(self, inFileName):
		p = subprocess.Popen("./run %s" % inFileName, shell=True, stdout=dev_null)
		ret = p.wait()
		return ret

outDir = "out/"
def runWithSeed(s):
	sys.stdout.write(".")
	sys.stdout.flush()
	tag = str(s)
	#print "running with seed %s" % s
	f = outDir+"randinp%s.dat"%tag
	i = Input("randinp.template")
	i.writeInput({"SEED": s}, f)
	r = Runner()
	ret = r.call(f)
	return ret
	#print "NMSSMtools reports %d" % ret

if __name__=="__main__":

	seeds = [long(os.urandom(8).encode("hex"), 16) for i in range(100)]
	pool = multiprocessing.Pool(12)
	ret = pool.map(runWithSeed, seeds)
	print ""
	print ret
	print "All done!"
	dev_null.close()
	