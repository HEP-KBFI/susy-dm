import ConfigParser
import sys
import numpy

numpy.random.seed(None) #None should read from urandom according to docs


def read_parameters(filename):
  """ Reads the model parameter configuration file 'filename' and return a list with the parameter information dictionaries. """
  conf = ConfigParser.RawConfigParser()
  conf.read(filename)
  vars = []
  for s in conf.sections():
    varname = s
    var = {"name": varname}
    for (item_name, val) in conf.items(s):
      var[item_name] = val
    vars.append(var)
  return vars

def process_var(v):
  """ Returns a random value based on the configuration dictionary in 'var' """
  ran = None

  #If the distribution is set as uniform, use 'min' and 'max'
  if v["dist"] == "uniform":
    ran = numpy.random.uniform(v["min"], v["max"])

  #If the distribution is set as 'gauss', use 'mean' and 'sigma'
  elif v["dist"] == "gauss":
    ran = numpy.random.normal(v["mean"], v["sigma"])

  #If the distribution is set as 'sign', create +1/-1 with equal probability
  elif v["dist"] == "sign":
    r = numpy.random.random()
    if r>0.5: ran = 1
    else: ran = -1

  return ran

#def sample_gut_mssm(vars):
#  s = {}
#  for v in vars:
#
#    #Add the possible dependence on previous vars here
#    pass
#
#    #Sample the variable
#    ran = process_var(v)
#    if ran != None: s[v["name"]] = ran
#
#  #Set all the parameters for suspectSUGRA. This effectively defines the model
#  #                  &tb, &Mtp, &MbMb, &alfSMZ, &sgn,
#  #                  &gMG1, &gMG2, &gMG3,
#  #                  &gAl, &gAt, &gAb,
#  #                  &gMHu, &gMHd,
#  #                  &gMl2, &gMl3, &gMr2 ,&gMr3,
#  #                  &gMq2, &gMq3, &gMu2, &gMu3, &gMd2, &gMd3); //23 parameters in total
#  o = {}
#
#  o["tb"] = s["tb"]
#  o["Mtp"] = s["Mtp"]
#  o["MbMb"] = s["MbMb"]
#  o["alfSMZ"] = s["alfSMZ"]
#  o["sgn"] = s["sgn"]
#
#  o["gMG1"] = s["M1"]
#  o["gMG2"] = s["M2"]
#  o["gMG3"] = s["M3"]
#
#  o["gAl"] = s["Asl"]
#
#  o["gAt"] = s["Asq"]
#  o["gAb"] = s["Asq"]
#
#  o["gMHu"] = s["Mhu"]
#  o["gMHd"] = s["Mhd"]
#
#  o["gMl2"] = s["msl"]
#  o["gMl3"] = s["msl"]
#  o["gMr2"] = s["msl"]
#  o["gMr3"] = s["msl"]
#
#  o["gMq2"] = s["msq"]
#  o["gMq3"] = s["msq"]
#  o["gMu2"] = s["msq"]
#  o["gMu3"] = s["msq"]
#  o["gMd2"] = s["msq"]
#  o["gMd3"] = s["msq"]
#
#  return o
#
#def sample_3g_mssm(vars):
#  s = {}
#  for v in vars:
#
#    #Add the possible dependence on previous vars here
#    pass
#
#    #Sample the variable
#    ran = process_var(v)
#    if ran != None: s[v["name"]] = ran
#
#  #Set all the parameters for suspectSUGRA. This effectively defines the model
#  #                  &tb, &Mtp, &MbMb, &alfSMZ, &sgn,
#  #                  &gMG1, &gMG2, &gMG3,
#  #                  &gAl, &gAt, &gAb,
#  #                  &gMHu, &gMHd,
#  #                  &gMl2, &gMl3, &gMr2 ,&gMr3,
#  #                  &gMq2, &gMq3, &gMu2, &gMu3, &gMd2, &gMd3); //23 parameters in total
#  o = {}
#
#  o["tb"] = s["tb"]
#  o["Mtp"] = s["Mtp"]
#  o["MbMb"] = s["MbMb"]
#  o["alfSMZ"] = s["alfSMZ"]
#  o["sgn"] = s["sgn"]
#
#  o["gMG1"] = s["m1/2"]
#  o["gMG2"] = s["m1/2"]
#  o["gMG3"] = s["m1/2"]
#
#  o["gAl"] = s["A0"]
#
#  o["gAt"] = s["A0"]
#  o["gAb"] = s["A0"]
#
#  o["gMHu"] = s["mhu"]
#  o["gMHd"] = s["mhd"]
#
#  o["gMl2"] = s["msl2"]
#  o["gMl3"] = s["msl3"]
#  o["gMr2"] = s["msl2"]
#  o["gMr3"] = s["msl3"]
#
#  o["gMq2"] = s["msq2"]
#  o["gMq3"] = s["msq3"]
#  o["gMu2"] = s["msq2"]
#  o["gMu3"] = s["msq3"]
#  o["gMd2"] = s["msq2"]
#  o["gMd3"] = s["msq3"]
#
#  return o

def sample(var_defs):
  s = {}
  for v in var_defs:
    ran = process_var(v)
    if ran!=None:
      s[v["name"]] = ran

  return s

#def main():
#  #Initialize with random seed
#  numpy.random.seed(None) #None should read from urandom according to docs
#
#  v = read_parameters(sys.argv[1])
#  modelname = sys.argv[2]
#  N = int(sys.argv[3])
#
#  f = None
#  if modelname == "3GMSSM":
#    f = sample_3g_mssm
#  elif modelname == "GUTMSSM":
#    f = sample_gut_mssm
#  elif modelname == "NMSSM":
#    f = sample_nmssm
#  if f == None:
#    sys.exit(1)
#
#  for n in range(N):
#    s = f(v)
#    if modelname == "NMSSM":
#      print s["m0"], s["mhf"], s["a0"], s["tb"], s["Lambda"], s["aKappa"]#, s["sgn"], s["aLambda"]
#
#if __name__ == "__main__":
#  main()
