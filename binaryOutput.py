import numpy
import sys

def appendArray(arr):
    f = open(fn, "a+b")
    numpy.save(f, arr)
    f.close()

if __name__=="__main__":
    fn = "test.npy"


    for l in sys.stdin.readlines():
        try:
            line_data = map(float, l.split())
            appendArray(numpy.array(line_data))
        except ValueError:
            print "could not convert to array: %s" % l
            pass

f = open(fn, "rb")
while True:
    try:
        A = numpy.load(f)
        print A
    except IOError:
        print "EOF"
        break
f.close()
