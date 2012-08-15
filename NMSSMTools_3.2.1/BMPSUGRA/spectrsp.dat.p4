# NMSSMTools OUTPUT IN SLHA FORMAT
# Info about spectrum calculator
BLOCK SPINFO   # Program information
     1   NMSSMTools # Spectrum calculator
     2   3.2.0      # Version number
     8   0          # Higgs mass precision
     3   # Relic density too large (WMAP)
     3   # Muon magn. mom. more than 2 sigma away
# Input parameters
BLOCK MODSEL
    3     1         # NMSSM PARTICLE CONTENT
BLOCK SMINPUTS
     1     1.27920000E+02   # ALPHA_EM^-1(MZ)
     2     1.16639000E-05   # GF
     3     1.17200000E-01   # ALPHA_S(MZ)
     4     9.11870000E+01   # MZ
     5     4.21400000E+00   # MB(MB)
     6     1.71400000E+02   # MTOP (POLE MASS)
     7     1.77700000E+00   # MTAU
# SMINPUTS Beyond SLHA:
# MW:     0.80420000E+02
# MS:     0.19000000E+00
# MC:     0.14000000E+01
# VUS:     0.22000000E+00
# VCB:     0.40000000E-01
# VUB:     0.40000000E-02
BLOCK MINPAR
     4    -1.00000000E+00   # SIGMU
     3     3.05000000E+00   # TANBETA(MZ)
     1     7.80000000E+02   # M0(MGUT)
     2     7.75000000E+02   # M12(MGUT)
     5    -2.25000000E+03   # A0(MGUT)
BLOCK EXTPAR
    21     7.25000000E+05   # MHD^2 AT THE GUT SCALE
    22     4.65000000E+06   # MHU^2 AT THE GUT SCALE
    61     4.90000000E-01   # LAMBDA AT THE SUSY SCALE
    64    -9.65000000E+02   # AKAPPA AT THE GUT SCALE
# 
BLOCK MASS   # Mass spectrum 
#  PDG Code     mass             particle 
        25     4.05984124E+01   # lightest neutral scalar
        35     1.20776865E+02   # second neutral scalar
        45     6.83326214E+02   # third neutral scalar
        36     1.15328840E+02   # lightest pseudoscalar
        46     6.82970046E+02   # second pseudoscalar
        37     6.78443288E+02   # charged Higgs
   1000001     1.74230770E+03   #  ~d_L
   2000001     1.67388922E+03   #  ~d_R
   1000002     1.74089788E+03   #  ~u_L
   2000002     1.73674858E+03   #  ~u_R
   1000003     1.74230770E+03   #  ~s_L
   2000003     1.67388922E+03   #  ~s_R
   1000004     1.74089788E+03   #  ~c_L
   2000004     1.73674858E+03   #  ~c_R
   1000005     1.26382019E+03   #  ~b_1
   2000005     1.62512437E+03   #  ~b_2
   1000006     4.13282878E+02   #  ~t_1
   2000006     1.28743584E+03   #  ~t_2
   1000011     9.80803399E+02   #  ~e_L
   2000011     7.03205071E+02   #  ~e_R
   1000012     9.78292169E+02   #  ~nue_L
   1000013     9.80803399E+02   #  ~mu_L
   2000013     7.03205071E+02   #  ~mu_R
   1000014     9.78292169E+02   #  ~numu_L
   1000015     6.99984262E+02   #  ~tau_1
   2000015     9.79679755E+02   #  ~tau_2
   1000016     9.77151632E+02   #  ~nutau_L
   1000021     1.75889170E+03   #  ~g
   1000022    -6.11761043E+01   # neutralino(1)
   1000023     2.25747587E+02   # neutralino(2)
   1000025    -2.27459591E+02   # neutralino(3)
   1000035     3.34835191E+02   # neutralino(4)
   1000045     6.40538028E+02   # neutralino(5)
   1000024     2.13155820E+02   # chargino(1)
   1000037    -6.40512337E+02   # chargino(2)
# 
# Low energy observables
BLOCK LOWEN
# Exp. 2 Sigma: 3.03E-04 < BR(b -> s gamma) < 4.01E-04:
         1     3.89763468E-04   # BR(b -> s gamma)
        11     4.35529877E-04   # (BR(b -> s gamma)+Theor.Err.)
        12     3.24668457E-04   # (BR(b -> s gamma)-Theor.Err.)
# Exp. 2 Sigma: 4.99E-01 < Delta M_d < 5.15E-01:
         2     6.29357630E-01   # Delta M_d in ps^-1
        21     1.10722300E+00   # Delta M_d +Theor.Err.
        22     1.67559101E-01   # Delta M_d -Theor.Err.
# Exp. 2 Sigma: 1.753E+01 < Delta Ms < 1.801E+01:
         3     2.18044316E+01   # Delta M_s in ps^-1
        31     2.89868186E+01   # Delta M_s +Theor.Err.
        32     1.48560946E+01   # Delta M_s -Theor.Err.
# Exp. 95% C.L.: BR(Bs->mu+mu-) < 4.5E-08:
         4     3.54042790E-09   # BR(Bs -> mu+mu-)
        41     6.01335015E-09   # BR(Bs -> mu+mu-)+Theor.Err.
        42     1.71864129E-09   # BR(Bs -> mu+mu-)-Theor.Err.
# Exp. 2 Sigma: 8.90E-05 < BR(B+ > tau+ + nu_tau) < 2.45E-04:
         5     1.31652146E-04   # BR(B+ -> tau+ + nu_tau)
        51     2.63355790E-04   # BR(B+ -> tau+ + nu_tau) + Theor.Err.
        52     5.68179071E-05   # BR(B+ -> tau+ + nu_tau) - Theor.Err.
# BSM contr. to the muon anomalous magn. moment:
         6    -8.45333320E-11   # Del_a_mu
        61     1.97313617E-10   # Del_a_mu + Theor.Err.
        62    -3.66380281E-10   # Del_a_mu - Theor.Err.
# Minimal Exp.-SM (2 sigma):  8.77306222E-10
# Maximal Exp.-SM (2 sigma):  4.61144414E-09
# Omega h^2 (allowed: 0.094 < Omega h^2 < 0.136):
    10     1.45608646E-01   # Omega h^2
# 
BLOCK HMIX Q=  8.14769482E+02 # (STOP/SBOTTOM MASSES)
     1    -2.06112690E+02   # MUEFF
     2     3.04128408E+00   # TAN(BETA)
     3     2.42973411E+02   # V(Q)
     4     4.65393417E+05   # MA^2
     5     1.88249763E+04   # MP^2
# 
# 3*3 Higgs mixing
BLOCK NMHMIX
  1  1    -1.00428844E-01   # S_(1,1)
  1  2     1.58506189E-02   # S_(1,2)
  1  3     9.94817976E-01   # S_(1,3)
  2  1     3.16199508E-01   # S_(2,1)
  2  2     9.48543818E-01   # S_(2,2)
  2  3     1.68076421E-02   # S_(2,3)
  3  1     9.43362029E-01   # S_(3,1)
  3  2    -3.16248926E-01   # S_(3,2)
  3  3     1.00273117E-01   # S_(3,3)
# 
# 3*3 Pseudoscalar Higgs mixing
BLOCK NMAMIX
  1  1     1.05684897E-01   # P_(1,1)
  1  2     3.46507860E-02   # P_(1,2)
  1  3     9.93795767E-01   # P_(1,3)
  2  1     9.44334095E-01   # P_(2,1)
  2  2     3.09617736E-01   # P_(2,2)
  2  3    -1.11220387E-01   # P_(2,3)
# 
# 3rd generation sfermion mixing
BLOCK STOPMIX  # Stop mixing matrix
  1  1     1.71797952E-01   # Rst_(1,1)
  1  2     9.85132206E-01   # Rst_(1,2)
  2  1    -9.85132206E-01   # Rst_(2,1)
  2  2     1.71797952E-01   # Rst_(2,2)
BLOCK SBOTMIX  # Sbottom mixing matrix
  1  1     9.99972303E-01   # Rsb_(1,1)
  1  2     7.44267280E-03   # Rsb_(1,2)
  2  1    -7.44267280E-03   # Rsb_(2,1)
  2  2     9.99972303E-01   # Rsb_(2,2)
BLOCK STAUMIX  # Stau mixing matrix
  1  1     7.17904090E-03   # Rsl_(1,1)
  1  2     9.99974230E-01   # Rsl_(1,2)
  2  1    -9.99974230E-01   # Rsl_(2,1)
  2  2     7.17904090E-03   # Rsl_(2,2)
# 
# Gaugino-Higgsino mixing
BLOCK NMNMIX  # 5*5 Neutralino Mixing Matrix
  1  1    -3.46264847E-02   # N_(1,1)
  1  2     3.65241481E-02   # N_(1,2)
  1  3     3.84815711E-02   # N_(1,3)
  1  4     3.43603621E-01   # N_(1,4)
  1  5     9.36975300E-01   # N_(1,5)
  2  1    -1.58673296E-01   # N_(2,1)
  2  2     7.57618533E-02   # N_(2,2)
  2  3     7.04035305E-01   # N_(2,3)
  2  4     6.32933195E-01   # N_(2,4)
  2  5    -2.69838445E-01   # N_(2,5)
  3  1    -6.57813775E-02   # N_(3,1)
  3  2     7.97809569E-02   # N_(3,2)
  3  3    -7.05886640E-01   # N_(3,3)
  3  4     6.65149870E-01   # N_(3,4)
  3  5    -2.20471112E-01   # N_(3,5)
  4  1     9.84403486E-01   # N_(4,1)
  4  2     3.46558430E-02   # N_(4,2)
  4  3     6.76528168E-02   # N_(4,3)
  4  4     1.56643290E-01   # N_(4,4)
  4  5    -2.51937605E-02   # N_(4,5)
  5  1    -1.56964666E-02   # N_(5,1)
  5  2     9.92653061E-01   # N_(5,2)
  5  3    -7.78499521E-04   # N_(5,3)
  5  4    -1.19877654E-01   # N_(5,4)
  5  5     4.71840918E-03   # N_(5,5)
# 
BLOCK UMIX  # Chargino U Mixing Matrix
  1  1     9.47279966E-04   # U_(1,1)
  1  2     9.99999551E-01   # U_(1,2)
  2  1    -9.99999551E-01   # U_(2,1)
  2  2     9.47279966E-04   # U_(2,2)
# 
BLOCK VMIX  # Chargino V Mixing Matrix
  1  1     1.69036040E-01   # V_(1,1)
  1  2    -9.85609871E-01   # V_(1,2)
  2  1     9.85609871E-01   # V_(2,1)
  2  2     1.69036040E-01   # V_(2,2)
# 
# Higgs reduced couplings
# (as compared to a SM Higgs with same mass)
BLOCK REDCOUP
# H1
  1  1     1.66808315E-02   # U-type fermions
  1  2    -3.22351560E-01   # D-type fermions
  1  3    -1.62269473E-02   # W,Z bosons
  1  4     1.68023611E-01   # Gluons
  1  5     6.99803585E-02   # Photons
# H2
  2  1     9.98225983E-01   # U-type fermions
  2  2     1.01492161E+00   # D-type fermions
  2  3     9.99846525E-01   # W,Z bosons
  2  4     9.73353180E-01   # Gluons
  2  5     9.95465340E-01   # Photons
# H3
  3  1    -3.32813192E-01   # U-type fermions
  3  2     3.02795700E+00   # D-type fermions
  3  3    -6.60400048E-03   # W,Z bosons
  3  4     3.18951773E-01   # Gluons
  3  5     1.36546515E+00   # Photons
# A1
  4  1     3.64657006E-02   # U-type fermions
  4  2     3.39222180E-01   # D-type fermions
  4  3     0.00000000E+00   # W,Z bosons
  4  4     3.23086961E-02   # Gluons
  4  5     2.75940345E-01   # Photons
# A2
  5  1     3.25834678E-01   # U-type fermions
  5  2     3.03107709E+00   # D-type fermions
  5  3     0.00000000E+00   # W,Z bosons
  5  4     3.32384819E-01   # Gluons
  5  5     2.87544317E-01   # Photons
# 
# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE
BLOCK GAUGE Q=  1.66075080E+03 # (SUSY SCALE)
         1     3.64129120E-01   # g1(Q,DR_bar)
         2     6.42230470E-01   # g2(Q,DR_bar)
         3     1.03597647E+00   # g3(Q,DR_bar)
BLOCK YU Q=  1.66075080E+03 # (SUSY SCALE)
  3  3     8.82706505E-01   # HTOP(Q,DR_bar)
BLOCK YD Q=  1.66075080E+03 # (SUSY SCALE)
  3  3     4.38586620E-02   # HBOT(Q,DR_bar)
BLOCK YE Q=  1.66075080E+03 # (SUSY SCALE)
  3  3     3.18976708E-02   # HTAU(Q,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE
BLOCK AU Q=  1.66075080E+03 # (SUSY SCALE)
  3  3    -1.76576688E+03   # ATOP
BLOCK AD Q=  1.66075080E+03 # (SUSY SCALE)
  3  3    -3.73049498E+03   # ABOT
BLOCK AE Q=  1.66075080E+03 # (SUSY SCALE)
  2  2    -2.53124214E+03   # AMUON
  3  3    -2.52475759E+03   # ATAU
# 
# SOFT MASSES AT THE SUSY SCALE
BLOCK MSOFT Q=  1.66075080E+03 # (SUSY SCALE)
         1     3.33802880E+02   # M1
         2     6.13911717E+02   # M2
         3     1.67034352E+03   # M3
        21     3.75608745E+05   # M_HD^2
        22     8.11433609E+04   # M_HU^2
        31     9.79927280E+02   # M_eL
        32     9.79927280E+02   # M_muL
        33     9.78788649E+02   # M_tauL
        34     7.02149057E+02   # M_eR
        35     7.02149057E+02   # M_muR
        36     6.98938444E+02   # M_tauR
        41     1.67918666E+03   # M_q1L
        42     1.67918666E+03   # M_q2L
        43     1.25768079E+03   # M_q3L
        44     1.67451349E+03   # M_uR
        45     1.67451349E+03   # M_cR
        46     5.27836090E+02   # M_tR
        47     1.60904997E+03   # M_dR
        48     1.60904997E+03   # M_sR
        49     1.60663268E+03   # M_bR
# 
# NMSSM SPECIFIC PARAMETERS THE SUSY SCALE
BLOCK NMSSMRUN Q=  1.66075080E+03 # (SUSY SCALE)
     1     4.90000000E-01   # LAMBDA(Q,DR_bar)
     2     5.60657298E-02   # KAPPA(Q,DR_bar)
     3    -6.84452810E+02   # ALAMBDA
     4     1.51865898E+02   # AKAPPA
     5    -2.07574862E+02   # MUEFF
     6     0.00000000E+00   # XIF
     7     0.00000000E+00   # XIS
     8     0.00000000E+00   # MU'
     9     0.00000000E+00   # MS'^2
    10     6.10068505E+03   # MS^2
    12     0.00000000E+00   # M3H^2
# 
# GAUGE AND YUKAWA COUPLINGS AT THE GUT SCALE
BLOCK GAUGE Q=  1.74530370E+16 # (GUT SCALE)
         1     7.09781570E-01   # g1(MGUT,DR_bar), GUT normalization
         2     7.09781568E-01   # g2(MGUT,DR_bar)
         3     6.99955440E-01   # g3(MGUT,DR_bar)
BLOCK YU Q=  1.74530370E+16 # (GUT SCALE)
  3  3     6.46598753E-01   # HTOP(MGUT,DR_bar)
BLOCK YD Q=  1.74530370E+16 # (GUT SCALE)
  3  3     1.91123010E-02   # HBOT(MGUT,DR_bar)
BLOCK YE Q=  1.74530370E+16 # (GUT SCALE)
  3  3     2.35546217E-02   # HTAU(MGUT,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE GUT SCALE
BLOCK AU Q=  1.74530370E+16 # (GUT SCALE)
  3  3    -2.24999930E+03   # ATOP
BLOCK AD Q=  1.74530370E+16 # (GUT SCALE)
  3  3    -2.24999902E+03   # ABOT
BLOCK AE Q=  1.74530370E+16 # (GUT SCALE)
  2  2    -2.24999906E+03   # AMUON
  3  3    -2.24999904E+03   # ATAU
# 
# SOFT MASSES SQUARED AT THE GUT SCALE
BLOCK MSOFT Q=  1.74530370E+16 # (GUT SCALE)
         1     7.75000097E+02   # M1
         2     7.75000044E+02   # M2
         3     7.74999942E+02   # M3
        21     7.25001434E+05   # M_HD^2
        22     4.65000396E+06   # M_HU^2
        31     6.08399918E+05   # M_eL^2
        32     6.08399918E+05   # M_muL^2
        33     6.08399918E+05   # M_tauL^2
        34     6.08400026E+05   # M_eR^2
        35     6.08400026E+05   # M_muR^2
        36     6.08400026E+05   # M_tauR^2
        41     6.08399868E+05   # M_q1L^2
        42     6.08399868E+05   # M_q2L^2
        43     6.08400677E+05   # M_q3L^2
        44     6.08399867E+05   # M_uR^2
        45     6.08399867E+05   # M_cR^2
        46     6.08401519E+05   # M_tR^2
        47     6.08399921E+05   # M_dR^2
        48     6.08399921E+05   # M_sR^2
        49     6.08399921E+05   # M_bR^2
# 
# NMSSM SPECIFIC PARAMETERS AT THE GUT SCALE
BLOCK NMSSMRUN Q=  1.74530370E+16 # (GUT SCALE)
     1     6.53544839E-01   # LAMBDA(MGUT,DR_bar)
     2     8.21357379E-02   # KAPPA(MGUT,DR_bar)
     3    -2.24999465E+03   # ALAMBDA
     4    -9.64988758E+02   # AKAPPA
     6     0.00000000E+00   # XIF
     7     0.00000000E+00   # XIS
     8     0.00000000E+00   # MU'
     9     0.00000000E+00   # MS'^2
    10     1.41813436E+06   # MS^2
    12     0.00000000E+00   # M3H^2
# 
# FINE-TUNING parameter d(ln Mz^2)/d(ln PG^2)
         1     3.08687258E+02   # PG=MHU
         2    -8.52000265E+00   # PG=MHD
         3     9.66570525E+01   # PG=MS
         4    -5.81241008E+01   # PG=M0
         5     0.00000000E+00   # PG=M1
         6     0.00000000E+00   # PG=M2
         7     0.00000000E+00   # PG=M3
         8    -1.64429079E+02   # PG=M12
         9     0.00000000E+00   # PG=ALAMBDA
        10     1.16814046E+00   # PG=AKAPPA
        11    -1.63779061E+02   # PG=A0
        12     0.00000000E+00   # PG=XIF
        13     0.00000000E+00   # PG=XIS
        14     0.00000000E+00   # PG=MUP
        15     0.00000000E+00   # PG=MSP
        16     0.00000000E+00   # PG=M3H
        17    -9.02991457E+01   # PG=LAMBDA
        18    -5.78457201E-01   # PG=KAPPA
        19     7.60323988E+00   # PG=HTOP
        20    -6.30920683E+02   # PG=G0
        21    -1.16688440E+01   # PG=MGUT
        22     3.08687258E+02   # MAX
        23                  1   # IMAX
# 
# REDUCED CROSS SECTIONS AT LHC
        11     2.78887576E-04   # VBF -> H1 -> tautau
        12     0.00000000E+00   # ggF -> H1 -> ZZ
        13     0.00000000E+00   # ggF -> H1 -> WW
        14     1.40933512E-03   # ggF -> H1 -> gammagamma
        15     1.31445976E-05   # VBF -> H1 -> gammagamma
        21     1.92939637E-01   # VBF -> H2 -> tautau
        22     1.77458726E-01   # ggF -> H2 -> ZZ
        23     1.77458726E-01   # ggF -> H2 -> WW
        24     1.75917205E-01   # ggF -> H2 -> gammagamma
        25     1.85623987E-01   # VBF -> H2 -> gammagamma
        31     8.81993419E-03   # VBF -> H3 -> tautau
        32     9.78624263E-05   # ggF -> H3 -> ZZ
        33     9.78624263E-05   # ggF -> H3 -> WW
        34     4.18396825E+00   # ggF -> H3 -> gammagamma
        35     1.79371125E-03   # VBF -> H3 -> gammagamma
