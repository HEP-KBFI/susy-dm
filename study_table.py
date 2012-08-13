from __future__ import division
import tables
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy

h5file = tables.openFile("nmssm1.h5", mode = "r")
t = h5file.root.NMSSM1.parspace

m_low = 123
m_high = 129

def get_points_with_sel(sel, step=1):
    s = t.where(sel, step=step)
    p = [(x["h1_mass"], x["chi1_mass"], exclusions(x)) for x in s]
    return p
    #return (map(lambda x: x[0], p), map(lambda x: x[1], p), map(lambda x: x[2], p))

def not_excluded(i, step=1):
	r = ""
	for a in i:
		r += "(PROB%d==0)&"%a
	r = r[:-1]
	return r

def exclusions(x):
	excl = []
	for i in range(1,54):
		if x["PROB%d"%i]!=0:
			excl.append(i)
	return excl

def get_primary_exclusions(p):
	return sorted([(y, len(filter(lambda x: y in x[2], p))/len(p)) for y in range(1,54)], key=lambda x: x[1], reverse=True)[0:3]


def exclusion_chain():
	A = get_points_with_sel("h1_mass>0", step=1000)
	excls = []
	for x in range(10):
		c_excl = get_primary_exclusions(A)[0][0]
		excls.append(c_excl)
		A = get_points_with_sel("(h1_mass>0)&" + not_excluded(excls))
	return excls

def draw_with_excl(excl=None, tag=None):
	fig = plt.figure(figsize=(20,20), dpi=1000)
	ax1 = fig.add_subplot(111)
	plt.xlabel("Higgs mass GeV/c**2")
	plt.ylabel("chi0 mass GeV/c**2")
	sel = "(h1_mass>0)&"+not_excluded(excl)
	s = max(int(len(t.getWhereList(sel))/50000),1)
	A = get_points_with_sel("(h1_mass>0)&"+not_excluded(excl), step=s)
	ax1.scatter(map(lambda x: x[0], A), map(lambda x: x[1], A), s=10.0, marker="o", c="b", alpha=0.05)
	plt.show()
	plt.savefig("/home/joosep/web/nmssm_%s.png"%tag)

draw_with_excl(excl=[30], tag="allowed_wmap")
r = range(1,54)
draw_with_excl(excl=r, tag="")
#def probs(x):
#    ret = []
#    for i in range(1,54):
#        if not x["PROB%d" % i]==0:
#            ret.append(i)
#    return ret
#
#sel = t.iterrows()
#points = [(x["h1_mass"], x["chi1_mass"]) for x in sel]
#h1_mass = map(lambda x: x[0], points)
#chi1_mass = map(lambda x: x[1], points)
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#plt.xlabel("Higgs mass GeV/c**2")
#plt.ylabel("chi0 mass GeV/c**2")
#ax1.scatter(h1_mass, chi1_mass, s=2.0, marker=".")
#sel = t.where("PROB==0")
#points = [(x["h1_mass"], x["chi1_mass"]) for x in sel]
#h1_mass = map(lambda x: x[0], points)
#chi1_mass = map(lambda x: x[1], points)
#ax1.scatter(h1_mass, chi1_mass, s=2.0, marker=".", c="r")
#plt.show()
#plt.savefig("/home/joosep/web/nmssm1.png")