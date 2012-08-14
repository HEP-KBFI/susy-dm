from __future__ import division
import tables
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy
import sys
import scipy
import scipy.stats
import logging

h5file = tables.openFile("/Users/joosep/Desktop/nmssm1.h5", mode = "r")
t = h5file.root.NMSSM1.parspace
logging.basicConfig(level=logging.DEBUG)

def chi_h_points(selection, recreate=True):
	logger = logging
	logger.debug("Starting to get chi1/h1 points")
	li = t.getWhereList(selection)
	sel = t.readCoordinates(li)
	l = len(li)

	if recreate:
		logger.debug("Recreating file")
		tempfile = tables.openFile("temp.h5", mode="w")
		filters = tables.Filters(complevel=9, complib='zlib', fletcher32=True)
		chi1_h1_arr = tempfile.createCArray(tempfile.root,'chi1_h1',tables.Float32Atom(),shape=(l,2), filters=filters)
		logger.debug("Putting chi1_mass to file")
		chi1_h1_arr[:,0] = sel[:]["chi1_mass"]
		logger.debug("Putting h1_mass to file")
		chi1_h1_arr[:,1] = sel[:]["h1_mass"]
		tempfile.flush()
	else:
		logger.debug("Opening existing file")
		tempfile = tables.openFile("temp.h5", mode="r")
		chi1_h1_arr = tempfile.root.chi1_h1
	logger.debug("Got chi1/h1 points: %d" % len(chi1_h1_arr[:,0]))
	return ((chi1_h1_arr[:,0], chi1_h1_arr[:,1]), tempfile)

m_low = 123
m_high = 129

def eval_density(points_x, points_y):
	logging.debug("Starting to evaluate density with KDE")
	xlow, xhigh = 0,500
	ylow, yhigh = 100,130
	x,y = numpy.mgrid[xlow:xhigh:1, ylow:yhigh:1]

	logging.debug("Preparing grid points for evaluation")
	evalpoints = numpy.vstack([x.ravel(), y.ravel()])
	logging.debug("Preparing data points for evaluation")
	density_points = numpy.vstack([points_x.T, points_y.T])
	density_points_for_calc = density_points[:,0:1000000]
	logging.debug("Preparing KDE")
	gkde = scipy.stats.kde.gaussian_kde(density_points_for_calc)
	logging.debug("Evaluating KDE over grid on %d points with %d data points"% (len(evalpoints.T), len(density_points_for_calc.T))
	densities = gkde.evaluate(evalpoints)
	densities = numpy.reshape(densities.T, x.shape)

	densities = densities*10000.0+1.0
	densities = numpy.log(densities)
	logging.debug("Done evaluating density")
	return (x,y,density_points_for_calc, densities)

if __name__=="__main__":
	logging.debug("Getting data points")
	((chi1, h1), tempfile) = chi_h_points("h1_mass>0")
	logging.debug("Calculating density")
	(x,y, density_points_for_calc, densities) = eval_density(chi1, h1)
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	xlow, xhigh = min(chi1)-10,max(chi1)+10
	ylow, yhigh = 123,129
	plt.xlim(xlow, xhigh)
	plt.ylim(ylow, yhigh)
	logging.debug("Plotting density contours")
	ax1.contour(x,y, densities, 5)
	logging.debug("Plotting scatterpoints")
	ax1.scatter(density_points_for_calc[0], density_points_for_calc[1], c="k", alpha=0.1, s=10.0)
	plt.ylabel("Higgs mass GeV/c**2")
	plt.xlabel("chi0 mass GeV/c**2")
	plt.suptitle("NMSSM parameter scan")
	plt.show()

def points(sel):
	#A = [(x["h1_mass"], x["chi1_mass"]) for x in t.where(sel)]
	
	l = len(t.getWhereList(sel))
	arr = numpy.empty((2,l))

	i=0
	for x in t.where(sel):
		arr[0,i] = x["h1_mass"]
		arr[1,i] = x["chi1_mass"]
		i += 1
	return arr

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
	#plt.savefig("/home/joosep/web/nmssm_%s.png"%tag)

#draw_with_excl(excl=[30], tag="allowed_wmap")
#r = range(1,54)
#draw_with_excl(excl=r, tag="")


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