"""
*      PROB(1) =/= 0   chargino too light
*      PROB(2) =/= 0   excluded by Z -> neutralinos
*      PROB(3) =/= 0   charged Higgs too light
*      PROB(4) =/= 0   excluded by ee -> hZ 
*      PROB(5) =/= 0   excluded by ee -> hZ, h -> bb
*      PROB(6) =/= 0   excluded by ee -> hZ, h -> tautau
*      PROB(7) =/= 0   excluded by ee -> hZ, h -> invisible 
*      PROB(8) =/= 0   excluded by ee -> hZ, h -> 2jets
*      PROB(9) =/= 0   excluded by ee -> hZ, h -> 2photons
*      PROB(10) =/= 0  excluded by ee -> hZ, h -> AA -> 4bs
*      PROB(11) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus
*      PROB(12) =/= 0  excluded by ee -> hZ, h -> AA -> 2bs 2taus
*      PROB(13) =/= 0  excluded by Z -> hA (Z width)
*      PROB(14) =/= 0  excluded by ee -> hA -> 4bs
*      PROB(15) =/= 0  excluded by ee -> hA -> 4taus
*      PROB(16) =/= 0  excluded by ee -> hA -> 2bs 2taus
*      PROB(17) =/= 0  excluded by ee -> hA -> AAA -> 6bs
*      PROB(18) =/= 0  excluded by ee -> hA -> AAA -> 6taus
*      PROB(19) =/= 0  excluded by ee -> Zh -> ZAA -> Z + light pairs
*      PROB(20) =/= 0  excluded by stop -> b l sneutrino
*      PROB(21) =/= 0  excluded by stop -> neutralino c
*      PROB(22) =/= 0  excluded by sbottom -> neutralino b
*      PROB(23) =/= 0  squark/gluino too light
*      PROB(24) =/= 0  selectron/smuon too light
*      PROB(25) =/= 0  stau too light
*      PROB(26) =/= 0  lightest neutralino is not LSP
*      PROB(27) =/= 0  Landau Pole in l, k, ht, hb below MGUT
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(30) =/= 0  excluded by WMAP (checked only if OMGFLAG=1)
*      PROB(31) =/= 0  eff. Higgs self-couplings in Micromegas > 1
*      PROB(32) =/= 0  b->s gamma more than 2 sigma away
*      PROB(33) =/= 0  Delta M_s more than 2 sigma away
*      PROB(34) =/= 0  Delta M_d more than 2 sigma away
*      PROB(35) =/= 0  B_s->mu+mu- more than 2 sigma away
*      PROB(36) =/= 0  B+-> tau+nu_tau more than 2 sigma away
*      PROB(37) =/= 0  (g-2)_muon more than 2 sigma away
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> A gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
*      PROB(40) =/= 0  BR(B-->X_s mu+ mu-) more than 2 sigma away
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
*      PROB(45) =/= 0  excluded by A/H -> tautau MSSM (LHC)
*      PROB(46) =/= 0  excluded by H -> tautau SM (LHC)
*      PROB(47) =/= 0  excluded by WH -> bb SM (LHC)
*      PROB(48) =/= 0  excluded by H -> ZZ SM (LHC)
*      PROB(49) =/= 0  excluded by H -> WW SM (LHC)
*      PROB(50) =/= 0  excluded by H -> gammagamma SM (LHC)
*      PROB(51) =/= 0  excluded by t -> bH+ (LHC)
*      PROB(52) =/= 0  excluded by squark/gluino searches (LHC)
*      PROB(53) =/= 0  excluded by Xenon100 (checked only if OMGFLAG=2 or 4)
"""

import numpy
import matplotlib.pyplot as plt


# from tables import *
# class Point(IsDescription):
#     h1_mass      = Float32Col() #Higgs mass
#     chi1_mass      = Float32Col()   #lightest neutralino mass
#     prob1 = BoolCol()


def get_prob_indices(row):
	"""
	Get's the indices of the PROB(I) elements of a row that are nonzero.
	1. zip (I, PROB(I)) together
	2. select PROB(I)!=0
	3. get indices
	"""
	return map(lambda x: x[0], filter(lambda x: x[1]!=0.0, zip(range(1,54), row[27:80])))

f = open("/Users/joosep/Desktop/NMSSM1.dat")
l = f.readlines()
f.close()
l = filter(lambda x: len(x.split())==80, l)
l = map(lambda x: map(lambda y: float(y), x.split()), l)
arr = numpy.array(l)


def not_excluded(excl, a=arr):
	t = filter(lambda x: len(set(get_prob_indices(x)).intersection(set(excl)))==0, a)
	t = numpy.array(t)
	if len(t)==0:
		raise Exception("All points are excluded with %s" % excl)
	return t

arr = filter(lambda x: x[11]>122 and x[11]<128, arr) #Remove unphysical Higgs
arr = not_excluded([1,3], arr) #Remove chargino, charged higgs too light
arr = not_excluded([26], arr) #Lightest neutralino is not LSP
arr = not_excluded(range(4,26), arr) #Various exclusions


def draw(excl, axes, a=arr, c='b'):
	t = not_excluded(excl, a=arr)
	h_mass = t[:,11]
	chi0_mass = t[:,22]
	axes.scatter(h_mass, chi0_mass, c=c)

def get_primary_exclusions(a=arr):
	s = filter(lambda x: x[1]!=0, [(i, len(filter(lambda x: i in get_prob_indices(x), a))) for i in range(54)])
	return sorted(s, key=lambda x: x[1], reverse=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)
draw([], ax1, c="b",)
draw([30], ax1, c="r")
draw([30,31], ax1, c="g")
draw([30,31,37], ax1, c="y")
draw([30,31,37,32], ax1, c="k")
draw(range(1,55), ax1, c="m")
plt.show()
