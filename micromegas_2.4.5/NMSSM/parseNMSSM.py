import slhautil
import glob
import numpy
from tables import *

def getParticle(slha, id_n, block):
    f = filter(lambda x: int(x[0])==id_n, slha[block])
    if len(f)==1:
        return float(f[0][1])
    else:
        raise Exception("Error parsing SLHA block %s for PDG id %d" % (block, id_n))

def getWarnings(slha):
    warns = filter(lambda x: int(x[0])==3, slha["SPINFO"])
    os = ""
    for w in warns:
        os += ":".join(w) + ";"
    return os[:-1]

class ModelPoint(IsDescription):
    warnings = StringCol(1000)
    h1_mass = Float32Col()
    h2_mass = Float32Col()
    n1_mass = Float32Col()
    run_id = Int32Col()
    n_warnings = Int32Col()

h5file = openFile("tutorial1.h5", mode = "w", title = "Test file")
table = h5file.createTable("/", "modelPoint", ModelPoint)

point = table.row

files = glob.glob("out/spectr*.dat")

masses = []
for f in files:
    s = slhautil.readSLHA(f)

    run_id = int(f[f.index("_")+1:f.index(".dat")])

    point["h1_mass"] = getParticle(s,25, "MASS")
    point["h2_mass"] = getParticle(s,35, "MASS")
    point["n1_mass"] = getParticle(s, 1000022, "MASS")
    warns = getWarnings(s)
    point["warnings"] = warns
    point["run_id"] = run_id
    point["n_warnings"] = warns.count("3")

    point.append()

table.flush()

out = [(x["run_id"], x["h1_mass"], x["n1_mass"], x["warnings"].split(";")) for x in table.where("(h1_mass>124.0) & (h1_mass<126.0)")]

no_warns = numpy.array([(x["run_id"], x["h1_mass"], x["n1_mass"]) for x in table.where("n_warnings==0")])

h1_masses = no_warns[:,1]
n1_masses = no_warns[:,2]

import matplotlib.pyplot as plt

plt.scatter(h1_masses, n1_masses)
plt.ylabel("neutralino(1) mass (GeV)")
plt.xlabel("higgs(1) mass (GeV)")
plt.xlim(100,130)
plt.ylim(0,2000)
plt.title("NMSSM output points")
plt.show()
plt.savefig("/home/joosep/web/nmssm.png")

#[['25', '1.13418412E+02', ' lightest neutral scalar'],
# ['35', '5.43710780E+03', ' second neutral scalar'],
# ['45', '1.74086057E+04', ' third neutral scalar'],
# ['36', '3.34443744E+03', ' lightest pseudoscalar'],
# ['46', '1.74752434E+04', ' second pseudoscalar'],
# ['37', '1.74752257E+04', ' charged Higgs'],
# ['1000001', '1.57949833E+03', ' ~d_L'],
# ['2000001', '1.34530732E+03', ' ~d_R'],
# ['1000002', '1.57761248E+03', ' ~u_L'],
# ['2000002', '8.46793897E+02', ' ~u_R'],
# ['1000003', '1.57949833E+03', ' ~s_L'],
# ['2000003', '1.34530732E+03', ' ~s_R'],
# ['1000004', '1.57761248E+03', ' ~c_L'],
# ['2000004', '8.46793897E+02', ' ~c_R'],
# ['1000005', '1.37465507E+03', ' ~b_1'],
# ['2000005', '1.82647816E+03', ' ~b_2'],
# ['1000006', '1.38872780E+03', ' ~t_1'],
# ['2000006', '1.95399677E+03', ' ~t_2'],
# ['1000011', '2.27729024E+02', ' ~e_L'],
# ['2000011', '4.95383385E+02', ' ~e_R'],
# ['1000012', '2.14203337E+02', ' ~nue_L'],
# ['1000013', '2.26446142E+02', ' ~mu_1'],
# ['2000013', '4.95971142E+02', ' ~mu_2'],
# ['1000014', '2.14203337E+02', ' ~numu_L'],
# ['1000015', '5.07450368E+02', ' ~tau_1'],
# ['2000015', '1.50393015E+03', ' ~tau_2'],
# ['1000016', '1.49657403E+03', ' ~nutau_L'],
# ['1000021', '1.22345813E+03', ' ~g'],
# ['1000022', '4.89273286E+02', ' neutralino(1)'],
# ['1000023', '1.60318654E+03', ' neutralino(2)'],
# ['1000025', '-1.66477122E+03', ' neutralino(3)'],
# ['1000035', '1.71543361E+03', ' neutralino(4)'],
# ['1000045', '-5.88681368E+03', ' neutralino(5)'],
# ['1000024', '1.60320122E+03', ' chargino(1)'],
# ['1000037', '-1.71548781E+03', ' chargino(2)']]

