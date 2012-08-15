# RELIC DENSITY CALCULATED BY MICROMEGAS
#
BLOCK RDINFO   # Program information
     1   MicrOmegas # Dark matter package
     2   3.2.0      # Version number
     3   # Relic density too large (WMAP)
#
BLOCK RELDEN
  0.14560865E+00   # omega h^2
  0.23845793E+02   # Xf
#
BLOCK LSP
  0.61176104E+02   # LSP mass
 ~o1 = -0.035*bino+0.037*wino+0.038*higgsino1+0.344*higgsino2+0.937*singlino
#
BLOCK CHANNELS
31.1% ~o1 ~o1 -> h1 h1
16.4% ~o1 ~o1 -> b B
 9.2% ~o1 ~o1 -> s S
 9.2% ~o1 ~o1 -> d D
 7.4% ~o1 ~o1 -> c C
 7.1% ~o1 ~o1 -> u U
 4.2% ~o1 ~o1 -> nm Nm
 4.2% ~o1 ~o1 -> ne Ne
 4.2% ~o1 ~o1 -> nl Nl
 2.8% ~o1 ~o1 -> l L
 2.1% ~o1 ~o1 -> m M
 2.1% ~o1 ~o1 -> e E
 0.1% ~o1 ~o1 -> Z h1
