      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*
*   On input:      
*
*      tan(beta) at the scale MZ, lambda at the scale Q2
*      m0, M1/2, A0 at the scale MGUT      
*      non-universal parameters are allowed in the Higgs/gaugino sector
*
*      The input file contains lower and upper bounds for these parameters
*      NTOT: Total number of points scanned
*
*   On output:
*
*      PAR(1) = lambda
*      PAR(2) = kappa
*      PAR(3) = tan(beta)
*      PAR(4) = mu (effective mu term = lambda*s)
*      PAR(5) = Alambda
*      PAR(6) = Akappa
*      PAR(7) = mQ3**2
*      PAR(8) = mU3**2
*      PAR(9) = mD3**2
*      PAR(10) = mL3**2
*      PAR(11) = mE3**2
*      PAR(12) = AU3
*      PAR(13) = AD3
*      PAR(14) = AE3
*      PAR(15) = mQ2**2
*      PAR(16) = mU2**2
*      PAR(17) = mD2**2
*      PAR(18) = mL2**2
*      PAR(19) = mE2**2
*      PAR(20) = M1
*      PAR(21) = M2
*      PAR(22) = M3
*      PAR(23) = MA (diagonal doublet CP-odd mass matrix element)
*      PAR(24) = MP (diagonal singlet CP-odd mass matrix element)
*      PAR(25) = AE2
*
*      All these parameters are assumed to be defined in DRbar at the scale
*      Q2 which is either user defined or computed as (2*mQ2+mU2+mD2)/4
*
*      SMASS(1-3): CP-even masses (ordered)
*
*      SCOMP(1-3,1-3): Mixing angles: if HB(I) are the bare states,
*        HB(I) = Re(H1), Re(H2), Re(S), and HM(I) are the mass eigenstates, 
*        the convention is HB(I) = SUM_(J=1,3) SCOMP(J,I)*HM(J)
*        which is equivalent to HM(I) = SUM_(J=1,3) SCOMP(I,J)*HB(J)
*
*      PMASS(1-2): CP-odd masses (ordered)
*
*      PCOMP(1-2,1-2): Mixing angles: if AB(I) are the bare states,
*        AB(I) = Im(H1), Im(H2), Im(S), and AM(I) are the mass eigenstates, 
*        the convention is 
*        AM(I) = PCOMP(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
*              + PCOMP(I,2)*AB(3)
*
*      CMASS: Charged Higgs mass
*
*      CU,CD,CV,CJ,CG(i) Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                        a1,a2 (i=4,5) to up type fermions, down type
*                        fermions, gauge bosons, gluons and photons
*                        Note: CV(4)=CV(3)=0
*
*      WIDTH(i) Total decay width of h1,h2,h3,a1,a2 (i=1..5)
*               with the following branching ratios:
*      BRJJ(i) h1,h2,h3,a1,a2 -> gluon gluon
*      BRMM(i)        "       -> mu mu
*      BRLL(i)        "       -> tau tau
*      BRSS(i)        "       -> ss
*      BRCC(i)        "       -> cc
*      BRBB(i)        "       -> bb
*      BRTT(i)        "       -> tt
*      BRWW(i)        "       -> WW (BRWW(4)=BRWW(5)=0)
*      BRZZ(i)        "       -> ZZ (BRZZ(4)=BRZZ(5)=0)
*      BRGG(i)        "       -> gamma gamma
*      BRZG(i)        "       -> Z gamma
*      BRHIGGS(i)   (i=1..5)  -> other Higgses, including:
*        BRHAA(i,j)   hi      -> a1a1, a1a2, a2a2 (i=1..3, j=1..3)
*        BRHCHC(i)    hi      -> h+h- (i=1..3)
*        BRHAZ(i,j)   hi      -> Zaj  (i=1..3)
*        BRHCW(i)  h1,h2,h3   -> h+W- (i=1..3), a1,a2 -> h+W- (i=4,5)
*        BRHHH(i)     h2      -> h1h1, h3-> h1h1, h1h2, h2h2 (i=1..4)
*        BRAHA(i)     a2      -> a1hi (i=1..3)
*        BRAHZ(i,j)   ai      -> Zhj  (i=1,2, j=1..3)
*      BRSUSY(i)    (i=1..5)  -> susy particles, including:
*        BRNEU(i,j,k)         -> neutralinos j,k (i=1..5, j,k=1..5)
*        BRCHA(i,j)           -> charginos 11, 12, 22 (i=1..5, j=1..3)
*        BRHSQ(i,j)   hi      -> uLuL, uRuR, dLdL, dRdR, t1t1, t2t2,
*                                t1t2, b1b1, b2b2, b1b2 (i=1..3, j=1..10)
*        BRASQ(i,j)   ai      -> t1t2, b1b2 (i=1,2, j=1,2)
*        BRHSL(i,j)   hi      -> lLlL, lRlR, nLnL, l1l1, l2l2, l1l2,
*                                ntnt (i=1..3, j=1..7)
*        BRASL(i)     ai      -> l1l2 (i=1,2)
*
*      HCWIDTH  Total decay width of the charged Higgs
*               with the following branching ratios:
*      HCBRM         h+ -> mu nu_mu
*      HCBRL         "  -> tau nu_tau
*      HCBRSU        "  -> s u
*      HCBRBU        "  -> b u
*      HCBRSC        "  -> s c
*      HCBRBC        "  -> b c
*      HCBRBT        "  -> b t
*      HCBRWHT       "  -> neutral Higgs W+, including:
*        HCBRWH(i)   "  -> H1W+, H2W+, h3W+, a1W+, a2W+ (i=1..5)
*      HCBRSUSY      "  -> susy particles,including
*        HCBRNC(i,j) "  -> neutralino i chargino j (i=1..5, j=1,2)
*        HCBRSQ(i)   "  -> uLdL, t1b1, t1b2, t2b1, t2b2 (i=1..5)
*        HCBRSL(i)   "  -> lLnL, t1nt, t2nt (i=1..3)
*
*      MNEU(i)   Mass of neutralino chi_i (i=1,5, ordered in mass)
*      NEU(i,j)  chi_i components of bino, wino, higgsino u&d, singlino 
*                (i,j=1..5)
*
*      MCHA(i)       Chargino masses
*      U(i,j),V(i,j) Chargino mixing matrices
*
*  ERRORS: IFAIL = 0..16
*
*  IFAIL = 0         OK
*          1,3,5,7   m_h1^2 < 0
*          2,3,6,7   m_a1^2 < 0
*          4,5,6,7   m_h+^2 < 0
*          8         m_sfermion^2 < 0
*          9         l, k, tan(beta) or mu = 0
*          10        Violation of phenomenological constraint(s)
*          11,12,13  Problem in integration of RGEs
*          14,15     Convergence problem
*          16        No electroweak symmetry breaking
*          17        Scenario
*
*  Phenomenological constraints:
*
*      PROB(I)  = 0, I = 1..53: OK
*            
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
*
************************************************************************

      IMPLICIT NONE

      INTEGER NFL,NPROB,NPAR
      PARAMETER (NFL=17,NPROB=53,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,TOT,ITOT,NTOT,IDUM
      INTEGER ITER,ITERMU,MUFLAG,SFFLAG,Q2FIX,IM
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER AKFLAG,ALFLAG,OMGFLAG,MAFLAG,IHIGGS,IGMU

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),CHECK,GUTEST,DELMB
      DOUBLE PRECISION DETM,R,K,MUM,MDM,MSM,MUT,MDT,MST,Q2,Q2MIN
      DOUBLE PRECISION M0,M12,A0,SIGMU,MUFAIL,MUSTEP,MUINIT
      DOUBLE PRECISION M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN,LN,LNN
      DOUBLE PRECISION AKN,AKNN,ALN,ALNN,MHDN,MHDNN,MHUN,MHUNN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN,KN,KNN,MUN,MUNN
      DOUBLE PRECISION M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT
      DOUBLE PRECISION ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT
      DOUBLE PRECISION MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      DOUBLE PRECISION MU3GUT,MD3GUT,GAU,CMS_SUSY
      DOUBLE PRECISION M0CEN,M12CEN,TBCEN,A0CEN,M1CEN,M2CEN
      DOUBLE PRECISION M3CEN,MHDCEN,MHUCEN,LCEN,ALCEN,AKCEN
      DOUBLE PRECISION KCEN,MUCEN,XCEN,XDEV,X,M0DEV,M12DEV
      DOUBLE PRECISION TBDEV,A0DEV,M1DEV,M2DEV,M3DEV,MHDDEV
      DOUBLE PRECISION MHUDEV,LDEV,ALDEV,AKDEV,KDEV,MUDEV

      COMMON/SIGMU/SIGMU
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,
     . MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/DELMB/DELMB
      COMMON/DETM/DETM
      COMMON/MUFAIL/MUFAIL
      COMMON/BOUNDS/M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN,MHUN,MHUNN,
     . LN,LNN,ALN,ALNN,AKN,AKNN,KN,KNN,MUN,MUNN
      COMMON/INPPAR/M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP,
     . XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      COMMON/MCMCPAR/M0CEN,M12CEN,TBCEN,A0CEN,M1CEN,M2CEN,M3CEN,
     . MHDCEN,MHUCEN,LCEN,ALCEN,AKCEN,KCEN,MUCEN,XCEN,XDEV,X,
     . M0DEV,M12DEV,TBDEV,A0DEV,M1DEV,M2DEV,M3DEV,
     . MHDDEV,MHUDEV,LDEV,ALDEV,AKDEV,KDEV,MUDEV,IHIGGS,IGMU
      COMMON/STEPS/NTOT,IDUM
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/MSAVE/MUM,MDM,MSM,MUT,MDT,MST,IM

*   Initialization

      CALL INITIALIZE()
       DO I=1,NFL
       NFAIL(I)=0
      ENDDO
      TOT=0
      IFAIL=20

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

*   Initialization of the range of parameters that has passed all tests

      M0N=1d99
      M0NN=-1d99
      M12N=1d99
      M12NN=-1d99
      TBN=1d99
      TBNN=-1d99
      A0N=1d99
      A0NN=-1d99
      LN=1d99
      LNN=-1d99
      IF(M1FLAG.EQ.1)THEN
       M1N=1d99
       M1NN=-1d99
      ENDIF
      IF(M2FLAG.EQ.1)THEN
       M2N=1d99
       M2NN=-1d99
      ENDIF
      IF(M3FLAG.EQ.1)THEN
       M3N=1d99
       M3NN=-1d99
      ENDIF
      IF(MHUFLAG.EQ.1)THEN
       MHUN=1d99
       MHUNN=-1d99
      ENDIF
      IF(MHDFLAG.EQ.1)THEN
       MHDN=1d99
       MHDNN=-1d99
      ENDIF
      IF(ALFLAG.EQ.1)THEN
       ALN=1d99
       ALNN=-1d99
      ENDIF
      IF(AKFLAG.EQ.1)THEN
       AKN=1d99
       AKNN=-1d99
      ENDIF
      IF(MAFLAG.EQ.-5)THEN
       KN=1d99
       KNN=-1d99
       MUN=1d99
       MUNN=-1d99
      ENDIF

*   Beginning of the random scan

      DO ITOT=1,NTOT

      IF(ITOT.EQ.1)THEN

       M0=M0CEN
       M12=M12CEN
       PAR(3)=TBCEN
       A0=A0CEN
       PAR(1)=LCEN
       IF(ALFLAG.EQ.0)THEN
        ALINP=A0
       ELSE
        ALINP=ALCEN
       ENDIF
       IF(AKFLAG.EQ.0)THEN
        AKINP=A0
       ELSE
        AKINP=AKCEN
       ENDIF
       IF(M1FLAG.EQ.0)THEN
        M1INP=M12
       ELSE
        M1INP=M1CEN
       ENDIF
       IF(M2FLAG.EQ.0)THEN
        M2INP=M12
       ELSE
        M2INP=M2CEN
       ENDIF
       IF(M3FLAG.EQ.0)THEN
        M3INP=M12
       ELSE
        M3INP=M3CEN
       ENDIF
       IF(MHDFLAG.EQ.0)THEN
        MHDINP=M0**2
       ELSE
        MHDINP=MHDCEN
       ENDIF
       IF(MHUFLAG.EQ.0)THEN
        MHUINP=M0**2
       ELSE
        MHUINP=MHUCEN
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        PAR(2)=KCEN
        PAR(4)=MUCEN
       ENDIF

      ELSE

 111   M0=M0CEN*(1d0+M0DEV*GAU(IDUM))
       IF(M0.LE.0d0)GOTO 111
       IF(M0.GE.2000d0)GOTO 111
 222   M12=M12CEN*(1d0+M12DEV*GAU(IDUM))
       IF(M12.LE.0d0)GOTO 222
       IF(M12.GE.1000d0)GOTO 222
 333   PAR(1)=LCEN*(1d0+LDEV*GAU(IDUM))
       IF(PAR(1).LE.0d0 .OR. PAR(1).GT.7d-1)GOTO 333
 444   PAR(3)=TBCEN*(1d0+TBDEV*GAU(IDUM))
       IF(PAR(3).LT.1.5d0 .OR. PAR(3).GT.60d0)GOTO 444
       A0=A0CEN*(1d0+A0DEV*GAU(IDUM))
       IF(ALFLAG.EQ.0)THEN
        ALINP=A0
       ELSE
        ALINP=ALCEN*(1d0+ALDEV*GAU(IDUM))
       ENDIF
       IF(AKFLAG.EQ.0)THEN
        AKINP=A0
       ELSE
        AKINP=AKCEN*(1d0+AKDEV*GAU(IDUM))
       ENDIF
       IF(M1FLAG.EQ.0)THEN
        M1INP=M12
       ELSE
 555    M1INP=M1CEN*(1d0+M1DEV*GAU(IDUM))
        IF(M1INP.LT.0d0)GOTO 555
       ENDIF
       IF(M2FLAG.EQ.0)THEN
        M2INP=M12
       ELSE
 666    M2INP=M2CEN*(1d0+M2DEV*GAU(IDUM))
        IF(M2INP.LT.0d0)GOTO 666
       ENDIF
       IF(M3FLAG.EQ.0)THEN
        M3INP=M12
       ELSE
 777    M3INP=M3CEN*(1d0+M3DEV*GAU(IDUM))
        IF(M3INP.LT.0d0)GOTO 777
       ENDIF
       IF(MHDFLAG.EQ.0)THEN
        MHDINP=M0**2
       ELSE
        MHDINP=MHDCEN*(1d0+MHDDEV*GAU(IDUM))
       ENDIF
       IF(MHUFLAG.EQ.0)THEN
        MHUINP=M0**2
       ELSE
        MHUINP=MHUCEN*(1d0+MHUDEV*GAU(IDUM))
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        PAR(2)=KCEN*(1d0+KDEV*GAU(IDUM))
        PAR(4)=MUCEN*(1d0+MUDEV*GAU(IDUM))
       ENDIF
      ENDIF

      XIFINP=0d0
      XISINP=0d0
      MUPINP=0d0
      MSPINP=0d0
      MSINP=0d0
      M3HINP=0d0

      !WRITE(0,*)""
      !WRITE(0,*)"------------------------------------------------------"
      !WRITE(0,*)""
      !WRITE(0,*)"Point ",ITOT
      !WRITE(0,*)""
      !WRITE(0,*)M0,M12,PAR(3),A0,PAR(1),AKINP,ALINP,MHDINP,MHUINP
      !WRITE(0,*)""

*   LHC constraint from squark/gluino searches

c      PROB(52)=DDIM(1d0,M12/CMS_SUSY(M0))
c      IF(PROB(52).NE.0d0)THEN
c       IFAIL=10
c       GOTO 11
c      ENDIF

*   Check for singular parameters l, tan(beta)

      IF(PAR(1)*PAR(3).EQ.0d0)THEN
       IFAIL=9
       GOTO 11
      ENDIF

*   Initialization of algorithm parameters

      SFFLAG=0
      MUSTEP=SIGMU*10d0
      IF(IFAIL.NE.0 .AND. IFAIL.NE.10 .AND. IFAIL.NE.17)THEN
       MUINIT=SIGMU*DSQRT(PAR(3)*Q2MIN)-MUSTEP
      ELSE
       MUINIT=PAR(4)-MUSTEP
      ENDIF
 1    MUINIT=MUINIT+MUSTEP
      MUFAIL=MUINIT-MUSTEP
      MUFLAG=0
      MUT=0d0
      MDT=0d0
      MST=0d0
      IM=0

*   Guess parameters at Q2

      IF(IFAIL.NE.0 .AND. IFAIL.NE.10 .AND. IFAIL.NE.17)THEN
       IF(Q2FIX.EQ.0)THEN
        Q2=MAX(M0**2+4d0*M12**2,Q2MIN)
       ENDIF
       R=(1d0+PAR(3)**2)/(1.29d0*PAR(3)**2)
       K=(1d0-R)*(A0-2.24d0*M12)**2+7.84d0*M12**2
       PAR(5)=ALINP-R/2d0*(A0-2.24d0*M12)-.59d0*M12
       PAR(6)=AKINP
       PAR(7)=(1d0-R/2d0)*M0**2+7.02d0*M12**2-R/6d0*K
       PAR(8)=(1d0-R)*M0**2+6.6d0*M12**2-R/3d0*K
       PAR(9)=M0**2+6.55d0*M12**2
       PAR(10)=M0**2+.52d0*M12**2
       PAR(11)=M0**2+.15d0*M12**2
       PAR(12)=A0-R*(A0-2.24d0*M12)-3.97d0*M12
       PAR(13)=A0-R/6d0*(A0-2.24d0*M12)-3.93d0*M12
       PAR(14)=A0-.69d0*M12
       PAR(15)=M0**2+7.02d0*M12**2
       PAR(16)=M0**2+6.6d0*M12**2
       PAR(17)=M0**2+6.55d0*M12**2
       PAR(18)=M0**2+.52d0*M12**2
       PAR(19)=M0**2+.15d0*M12**2
       PAR(20)=.4d0*M1INP
       PAR(21)=.8d0*M2INP
       PAR(22)=2.4d0*M3INP
       PAR(23)=DSQRT(Q2)
       PAR(24)=DSQRT(Q2)
       PAR(25)=A0-.69d0*M12
       DELMB=.1d0
       IF(MAFLAG.NE.-5)THEN
        PAR(2)=PAR(1)/5d0
        PAR(4)=MUINIT
       ENDIF
       MUM=10D0*MHUINP
       MDM=MHDINP
       MSM=Q2
      ENDIF

      M1GUT=M1INP
      M2GUT=M2INP
      M3GUT=M3INP
      ALGUT=ALINP
      AKGUT=AKINP
      XIFGUT=XIFINP
      XISGUT=XISINP
      MUPGUT=MUPINP
      MSPGUT=MSPINP
      M3HGUT=M3HINP

!      WRITE(0,*)"New try"
!      WRITE(0,*)""
!      WRITE(0,*)"MUINIT =",MUINIT
!      WRITE(0,*)""

!      WRITE(0,*)"M1 =",PAR(20)
!      WRITE(0,*)"M2 =",PAR(21)
!      WRITE(0,*)"M3 =",PAR(22)
!      WRITE(0,*)"AL =",PAR(5)
!      WRITE(0,*)"AK =",PAR(6)
!      WRITE(0,*)"ATOP =",PAR(12)
!      WRITE(0,*)"ABOT =",PAR(13)
!      WRITE(0,*)"ATAU =",PAR(14)
!      WRITE(0,*)"AMUON =",PAR(25)
!      WRITE(0,*)"MQ3 =",PAR(7)
!      WRITE(0,*)"MU3 =",PAR(8)
!      WRITE(0,*)"MD3 =",PAR(9)
!      WRITE(0,*)"MQ =",PAR(15)
!      WRITE(0,*)"MU =",PAR(16)
!      WRITE(0,*)"MD =",PAR(17)
!      WRITE(0,*)"ML3 =",PAR(10)
!      WRITE(0,*)"ME3 =",PAR(11)
!      WRITE(0,*)"ML =",PAR(18)
!      WRITE(0,*)"ME =",PAR(19)
!      IF(MAFLAG.NE.-5)THEN
!       WRITE(0,*)"KAPPA =",PAR(2)
!      ENDIF
!      WRITE(0,*)""

*   Initialization of PROB and IFAIL

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO      
      IFAIL=0

*   Guess for GUT scale and GUT couplings

      CALL RGES(PAR,PROB,IFAIL)
      IF(IFAIL.NE.0.OR.PROB(27).NE.0d0)THEN
!       WRITE(0,*)"RGE integration problem 1"
       IF(DABS(MUINIT).LT.100d0 .AND. DABS(MUSTEP).GT.1.D-1
     .    .AND. MAFLAG.NE.-5)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=11)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=11
        GOTO 11
       ENDIF
      ENDIF

*   External loop to compute the soft parameters at Q2

      ITER=0
 21   ITER=ITER+1
!      WRITE(0,*)"ITER =",ITER
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RGESINV(PAR,IFAIL)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 2"
       IF(DABS(MUINIT).LT.100d0 .AND. DABS(MUSTEP).GT.1.D-1
     .    .AND. MAFLAG.NE.-5)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=13)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=13
        GOTO 11
       ENDIF
      ENDIF

*   Internal loop to compute mu, k and ms

      ITERMU=0
 22   ITERMU=ITERMU+1
!      WRITE(0,*)"ITERMU =",ITERMU
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RUNPAR(PAR)

      CALL MSFERM(PAR,IFAIL,0)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Negative sfermion mass squared"
       IF(MAFLAG.NE.-5)THEN
        IF(DABS(MUINIT).GT.30d0)THEN
         SFFLAG=1
         MUINIT=MUINIT-2d0*MUSTEP
         GOTO 1
        ELSE
!         WRITE(0,*)"Exit (IFAIL=8)"
!         WRITE(0,*)""
!         WRITE(0,*)""
         IFAIL=8
         GOTO 11
        ENDIF
       ENDIF
      ENDIF

      CALL MINIMIZE(PAR,CHECK)
      IF(ITER.GT.10.AND.DETM.GE.0d0)MUFLAG=1

      IF(CHECK.GT.1.D-12.AND.ITERMU.LT.10) GOTO 22
      IF(CHECK.GT.1.D-9.AND.ITERMU.LT.100) GOTO 22
      IF(CHECK.GT.1.D-6)THEN
       IFAIL=14
       GOTO 11
      ENDIF
      
      CALL RGES(PAR,PROB,IFAIL)
      IF(IFAIL.NE.0.OR.PROB(27).NE.0d0)THEN
!       WRITE(0,*)"RGE integration problem 3"
       IF(DABS(MUINIT).LT.100d0 .AND. DABS(MUSTEP).GT.1.D-1
     .    .AND. MAFLAG.NE.-5)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=11)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=11
        GOTO 11
       ENDIF
      ENDIF

      CALL RGESUNI(PAR,IFAIL,GUTEST)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 4"
       IF(DABS(MUINIT).LT.100d0 .AND. DABS(MUSTEP).GT.1.D-1
     .    .AND. MAFLAG.NE.-5)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=12)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=12
        GOTO 11
       ENDIF
      ENDIF

      IF(GUTEST.GT.1.D-12.AND.ITER.LT.10) GOTO 21
      IF(GUTEST.GT.1.D-9.AND.ITER.LT.100) GOTO 21
      IF(GUTEST.GT.1.D-6.AND.ITER.LT.1000) GOTO 21
      IF(GUTEST.GT.1.D-3)THEN
!       WRITE(0,*)"No convergence"
       IF(DABS(MUSTEP).GT.1.D-1 .AND. MAFLAG.NE.-5)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL+MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
        ITER=0
!        WRITE(0,*)"New try"
!        WRITE(0,*)"MUFAIL =",MUFAIL
!        WRITE(0,*)""
!        WRITE(0,*)""
        GOTO 21
       ELSE
!        WRITE(0,*)"Exit (IFAIL=15)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=15
        GOTO 11
       ENDIF
      ENDIF

*   Check if correct EWSB

      IF(DETM.LE.0d0 .AND. MAFLAG.NE.-5)THEN
!       WRITE(0,*)"Convergence in a false minimum"
       IF(MUFLAG.EQ.1 .AND. DABS(MUSTEP).GT.1.D-1)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL+MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
        ITER=0
!        WRITE(0,*)"New try"
!        WRITE(0,*)"MUFAIL =",MUFAIL
!        WRITE(0,*)""
!        WRITE(0,*)""
        GOTO 21
       ELSEIF(MUFLAG.EQ.0 .AND. DABS(MUSTEP).LT.DABS(MUFAIL))THEN
        MUFAIL=MUFAIL-MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
        ITER=0
!        WRITE(0,*)"New try"
!        WRITE(0,*)"MUFAIL =",MUFAIL
!        WRITE(0,*)""
!        WRITE(0,*)""
        GOTO 21
       ELSEIF(DETM.LT.0d0)THEN
!        WRITE(0,*)"Exit (IFAIL=16)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=16
        GOTO 11
       ENDIF
      ENDIF

*   Computation of sfermion masses:      

      CALL MSFERM(PAR,IFAIL,1)
      IF(IFAIL.NE.0)GOTO 11
      
*   Computation of Higgs masses

      CALL MHIGGS(PAR,IFAIL)
      IF(IFAIL.NE.0)GOTO 11

*   Computation of gluino mass

      CALL GLUINO(PAR)      
      
*   Computation of chargino masses

      CALL CHARGINO(PAR)

*   Computation of neutralino masses

      CALL NEUTRALINO(PAR)
      
*   Computation of Higgs + top branching ratios

      CALL DECAY(PAR)
      CALL TDECAY(PAR)

*   Exp. constraints

      CALL SUBEXP(PAR,PROB)

*   b -> s gamma + B physics

      CALL BSG(PAR,PROB)

*   Anom. magn. moment of the Muon

      CALL MAGNMU(PAR,PROB)

*   Global minimum?

      CALL CHECKMIN(PROB)

*   Relic density

      CALL RELDEN(PAR,PROB)

*   Check for problems

      DO I=1,45
       IF(PROB(I).NE.0d0 .AND. I.NE.31 .AND. I.NE.37)IFAIL=10
      ENDDO
      IF(PROB(37).NE.0d0 .AND. IGMU.EQ.1)IFAIL=10

      DO I=51,NPROB
       IF(PROB(I).NE.0d0)IFAIL=10
      ENDDO

*   Computation of the fine-tuning

      CALL FTPAR(PAR,1)

*   Sparticle decays

      CALL NMSDECAY(PAR)

*   Recording of the results

11    CALL MCMCSTEP(PAR,PROB,IFAIL)
      CALL OUTPUT(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       M0N=MIN(M0,M0N)
       M0NN=MAX(M0,M0NN)
       M12N=MIN(M12,M12N)
       M12NN=MAX(M12,M12NN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       A0N=MIN(A0,A0N)
       A0NN=MAX(A0,A0NN)
       LN=MIN(PAR(1),LN)
       LNN=MAX(PAR(1),LNN)
       ALN=MIN(ALINP,ALN)
       ALNN=MAX(ALINP,ALNN)
       AKN=MIN(AKINP,AKN)
       AKNN=MAX(AKINP,AKNN)
       KN=MIN(PAR(2),KN)
       KNN=MAX(PAR(2),KNN)
       MUN=MIN(PAR(4),MUN)
       MUNN=MAX(PAR(4),MUNN)
       MHDN=MIN(MHDINP,MHDN)
       MHDNN=MAX(MHDINP,MHDNN)
       MHUN=MIN(MHUINP,MHUN)
       MHUNN=MAX(MHUINP,MHUNN)
       M1N=MIN(M1INP,M1N)
       M1NN=MAX(M1INP,M1NN)
       M2N=MIN(M2INP,M2N)
       M2NN=MAX(M2INP,M2NN)
       M3N=MIN(M3INP,M3N)
       M3NN=MAX(M3INP,M3NN)
      ELSE
       NFAIL(IFAIL)=NFAIL(IFAIL)+1
      ENDIF

      ENDDO

*   Summary of the scanning:
*   Number of points that passed/failed the tests
*   and range for scanned parameters

      CALL ERROR(TOT,NTOT,NFAIL)

      END


      SUBROUTINE INPUT(PAR,NPAR)
      
*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX,ISEED
      INTEGER NTOT,N0,NLOOP,NBER,NPAR,ERR,NMSFLAG
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER AKFLAG,ALFLAG,OMGFLAG,MAFLAG,PFLAG,IHIGGS,IGMU

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,SIGMU,Q2MIN
      DOUBLE PRECISION M0CEN,M12CEN,TBCEN,A0CEN,M1CEN,M2CEN
      DOUBLE PRECISION M3CEN,MHDCEN,MHUCEN,LCEN,ALCEN,AKCEN
      DOUBLE PRECISION KCEN,MUCEN,XCEN,XDEV,X,M0DEV,M12DEV
      DOUBLE PRECISION TBDEV,A0DEV,M1DEV,M2DEV,M3DEV,MHDDEV
      DOUBLE PRECISION MHUDEV,LDEV,ALDEV,AKDEV,KDEV,MUDEV

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/SIGMU/SIGMU
      COMMON/MCMCPAR/M0CEN,M12CEN,TBCEN,A0CEN,M1CEN,M2CEN,M3CEN,
     . MHDCEN,MHUCEN,LCEN,ALCEN,AKCEN,KCEN,MUCEN,XCEN,XDEV,X,
     . M0DEV,M12DEV,TBDEV,A0DEV,M1DEV,M2DEV,M3DEV,
     . MHDDEV,MHUDEV,LDEV,ALDEV,AKDEV,KDEV,MUDEV,IHIGGS,IGMU
      COMMON/STEPS/NTOT,ISEED
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=0d0
      ENDDO
      SIGMU=1d99
      M0CEN=1d99
      M12CEN=1d99
      TBCEN=1d99
      A0CEN=1d99
      M1CEN=1d99
      M2CEN=1d99
      M3CEN=1d99
      MHDCEN=1d99
      MHUCEN=1d99
      LCEN=1d99
      ALCEN=1d99
      AKCEN=1d99
      KCEN=1d99
      MUCEN=1d99
      XCEN=1d99
      XDEV=1d99
      M0DEV=1d99
      M12DEV=1d99
      TBDEV=1d99
      A0DEV=1d99
      M1DEV=1d99
      M2DEV=1d99
      M3DEV=1d99
      MHDDEV=1d99
      MHUDEV=1d99
      LDEV=1d99
      ALDEV=1d99
      AKDEV=1d99
      KDEV=1d99
      MUDEV=1d99
      NTOT=0
      IHIGGS=2
      IGMU=0

*   DEFAULT VALUES FOR FLAGS
      OMGFLAG=0
      PFLAG=0
      NMSFLAG=0
      M1FLAG=0
      M2FLAG=0
      M3FLAG=0
      MHDFLAG=0
      MHUFLAG=0
      ALFLAG=0
      AKFLAG=0

*   DEFAULT VALUE FOR THE RANDOM SEED (DEFAULT: -1)
      ISEED=-1

*   DEFAULT VALUE FOR THE RENSCALE Q2 (DEFAULT: 0d0)
      Q2=0d0

*   INITIALIZE READ LOOP
      NLINE=0
      CHBLCK=' '

*   START TO READ NEW LINE INTO CHINL
 21   CHINL=' '

*   LINE NUMBER
      NLINE=NLINE+1

      READ(5,'(A120)',END=29,ERR=999) CHINL
      
*   CHECK FOR EMPTY OR COMMENT LINES
      IF(CHINL.EQ.' '.OR.CHINL(1:1).EQ.'#'
     .  .OR.CHINL(1:1).EQ.'*') GOTO 21

*   FORCE UPPER CASE LETTERS IN CHINL (AS REQUIRED BELOW)
      INL=0
 22   INL=INL+1
      IF(CHINL(INL:INL).NE.'#')THEN
       DO ICH=97,122
        IF(CHINL(INL:INL).EQ.CHAR(ICH)) CHINL(INL:INL)=CHAR(ICH-32)
       ENDDO
       IF(INL.LT.120) GOTO 22
      ENDIF

*   CHECK FOR BLOCK STATEMENT
      IF(CHINL(1:1).EQ.'B')THEN
       READ(CHINL,'(A6,A)',ERR=999) CHDUM,CHBLCK
       GOTO 21
      ENDIF

*   CHECK FOR NMSSM MODEL IN MODSEL
*   IF THE RELIC DENSITY SHOULD BE COMPUTED
*   THE BLOCK MODSEL MUST CONTAIN THE LINE "  9     1    "
      IF(CHBLCK(1:6).EQ.'MODSEL')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.8) PFLAG=IVAL
       IF(IX.EQ.9) OMGFLAG=IVAL
       IF(IX.EQ.13) NMSFLAG=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL
 
*   READ GUT PARAMETERS, SIGMU, Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.1) M0CEN=VAL
       IF(IX.EQ.16) M0DEV=VAL
       IF(IX.EQ.2) M12CEN=VAL
       IF(IX.EQ.26) M12DEV=VAL
       IF(IX.EQ.3) TBCEN=VAL
       IF(IX.EQ.36) TBDEV=VAL
       IF(IX.EQ.4) SIGMU=VAL
       IF(IX.EQ.5) A0CEN=VAL
       IF(IX.EQ.56) A0DEV=VAL
 
*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.-1) XDEV=VAL
       IF(IX.EQ.1) M1CEN=VAL
       IF(IX.EQ.106) M1DEV=VAL
       IF(IX.EQ.2) M2CEN=VAL
       IF(IX.EQ.206) M2DEV=VAL
       IF(IX.EQ.3) M3CEN=VAL
       IF(IX.EQ.306) M3DEV=VAL
       IF(IX.EQ.21) MHDCEN=VAL
       IF(IX.EQ.216) MHDDEV=VAL
       IF(IX.EQ.22) MHUCEN=VAL
       IF(IX.EQ.226) MHUDEV=VAL
       IF(IX.EQ.61) LCEN=VAL
       IF(IX.EQ.616) LDEV=VAL
       IF(IX.EQ.62) KCEN=VAL
       IF(IX.EQ.626) KDEV=VAL
       IF(IX.EQ.63) ALCEN=VAL
       IF(IX.EQ.636) ALDEV=VAL
       IF(IX.EQ.64) AKCEN=VAL
       IF(IX.EQ.646) AKDEV=VAL
       IF(IX.EQ.65) MUCEN=VAL
       IF(IX.EQ.656) MUDEV=VAL

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.0) NTOT=IVAL
       IF(IX.EQ.1) ISEED=IVAL
       IF(IX.EQ.2) IHIGGS=IVAL
       IF(IX.EQ.3) IGMU=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   ERR=0
      IF(XDEV.EQ.1d99)THEN
       WRITE(0,1)"XDEV MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(M0CEN.EQ.1d99)THEN
       WRITE(0,1)"M0CEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(A0CEN.EQ.1d99)THEN
       WRITE(0,1)"A0CEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(TBCEN.EQ.1d99)THEN
       WRITE(0,1)"TBCEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(DABS(SIGMU).NE.1d0 .AND. MUCEN.EQ.1d99)THEN
       WRITE(0,1)"SIGMU IS EITHER 1 OR -1 IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(SIGMU.NE.1d99 .AND. MUCEN.NE.1d99)THEN
       WRITE(0,1)"BOTH MUCEN AND SIGMU CANNOT BE GIVEN"
       ERR=1
      ENDIF
      IF(LCEN.EQ.1d99)THEN
       WRITE(0,1)"LCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KCEN.NE.1d99 .AND. MUCEN.EQ.1d99)THEN
       WRITE(0,1)"IF KCEN IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MUCEN MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUCEN.NE.1d99 .AND. KCEN.EQ.1d99)THEN
       WRITE(0,1)"IF MUCEN IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN KCEN MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUCEN.NE.1d99 .AND. MHDCEN.NE.1d99)THEN
       WRITE(0,1)"BOTH MUCEN AND MHDCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUCEN.NE.1d99 .AND. MHUCEN.NE.1d99)THEN
       WRITE(0,1)"BOTH MUCEN AND MHUCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M0DEV.EQ.1d99)THEN
       WRITE(0,1)"M0DEV MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(A0DEV.EQ.1d99)THEN
       WRITE(0,1)"A0DEV MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(TBDEV.EQ.1d99)THEN
       WRITE(0,1)"TBDEV MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(SIGMU.NE.1d99 .AND. MUDEV.NE.1d99)THEN
       WRITE(0,1)"BOTH MUDEV AND SIGMU CANNOT BE GIVEN"
       ERR=1
      ENDIF
      IF(LDEV.EQ.1d99)THEN
       WRITE(0,1)"LDEV MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KDEV.NE.1d99 .AND. MUDEV.EQ.1d99)THEN
       WRITE(0,1)"IF KDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MUDEV MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUDEV.NE.1d99 .AND. KDEV.EQ.1d99)THEN
       WRITE(0,1)"IF MUDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN KDEV MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUDEV.NE.1d99 .AND. MHDDEV.NE.1d99)THEN
       WRITE(0,1)"BOTH MUDEV AND MHDDEV CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUDEV.NE.1d99 .AND. MHUDEV.NE.1d99)THEN
       WRITE(0,1)"BOTH MUDEV AND MHUDEV CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KDEV.NE.1d99 .AND. MUDEV.EQ.1d99)THEN
       WRITE(0,1)"IF KDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MUDEV MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUCEN.NE.1d99 .AND. KDEV.EQ.1d99)THEN
       WRITE(0,1)"IF MUCEN IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN KDEV MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUCEN.NE.1d99 .AND. MHDDEV.NE.1d99)THEN
       WRITE(0,1)"BOTH MUCEN AND MHDDEV CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUCEN.NE.1d99 .AND. MHUDEV.NE.1d99)THEN
       WRITE(0,1)"BOTH MUCEN AND MHUDEV CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KDEV.NE.1d99 .AND. MUCEN.EQ.1d99)THEN
       WRITE(0,1)"IF KDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MUCEN MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUDEV.NE.1d99 .AND. KCEN.EQ.1d99)THEN
       WRITE(0,1)"IF MUDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN KCEN MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUDEV.NE.1d99 .AND. MHDDEV.NE.1d99)THEN
       WRITE(0,1)"BOTH MUDEV AND MHDDEV CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUDEV.NE.1d99 .AND. MHUCEN.NE.1d99)THEN
       WRITE(0,1)"BOTH MUDEV AND MHUCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KCEN.NE.1d99 .AND. KDEV.EQ.1d99)THEN
       WRITE(0,1)"IF KCEN IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN KDEV MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUCEN.NE.1d99 .AND. MUDEV.EQ.1d99)THEN
       WRITE(0,1)"IF MUCEN IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MUDEV MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MHUCEN.NE.1d99 .AND. MHUDEV.EQ.1d99)THEN
       WRITE(0,1)"IF MHUCEN IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MHUDEV MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MHDCEN.NE.1d99 .AND. MHDDEV.EQ.1d99)THEN
       WRITE(0,1)"IF MHDCEN IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MHDDEV MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(KDEV.NE.1d99 .AND. KCEN.EQ.1d99)THEN
       WRITE(0,1)"IF KDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN KCEN MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MUDEV.NE.1d99 .AND. MUCEN.EQ.1d99)THEN
       WRITE(0,1)"IF MUDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MUCEN MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MHUDEV.NE.1d99 .AND. MHUCEN.EQ.1d99)THEN
       WRITE(0,1)"IF MHUDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MHUCEN MUST BE AS WELL"
       ERR=1
      ENDIF
      IF(MHDDEV.NE.1d99 .AND. MHDCEN.EQ.1d99)THEN
       WRITE(0,1)"IF MHDDEV IS GIVEN IN BLOCK EXTPAR"
       WRITE(0,1)"THEN MHDCEN MUST BE AS WELL"
       ERR=1
      ENDIF

*   Total number of points

      IF(NTOT.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF

*   Set MAFLAG, SCANFLAGS

      IF(MUCEN.NE.1d99)THEN
       MAFLAG=-5
      ELSE
       MAFLAG=-1
      ENDIF
      IF(M1CEN.NE.1d99)M1FLAG=1
      IF(M2CEN.NE.1d99)M2FLAG=1
      IF(M3CEN.NE.1d99)M3FLAG=1
      IF(MHDCEN.NE.1d99)MHDFLAG=1
      IF(MHUCEN.NE.1d99)MHUFLAG=1
      IF(ALCEN.NE.1d99)ALFLAG=1
      IF(AKCEN.NE.1d99)AKFLAG=1

      IF(M1FLAG*M2FLAG*M3FLAG.EQ.0)THEN
       IF(M12CEN.EQ.1d99)THEN
        WRITE(0,1)"M12 MUST BE GIVEN IN BLOCK CENPAR"
        ERR=1
       ENDIF
      ELSE
       IF(M12CEN.NE.1d99)THEN
        WRITE(0,1)"WARNING: M12 IS NOT USED"
       ENDIF
      ENDIF

*   Set default values and guesses

      IF(M1CEN.EQ.1d99)M1CEN=M12CEN
      IF(M2CEN.EQ.1d99)M2CEN=M12CEN
      IF(M3CEN.EQ.1d99)M3CEN=M12CEN
      IF(MHDCEN.EQ.1d99)MHDCEN=M0CEN**2
      IF(MHUCEN.EQ.1d99)MHUCEN=M0CEN**2
      IF(ALCEN.EQ.1d99)ALCEN=A0CEN
      IF(AKCEN.EQ.1d99)AKCEN=A0CEN

*   Stop if error

      IF(ERR.EQ.1)THEN
       WRITE(0,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF

*   Set Q2MIN, Q2FIX:
      Q2MIN=100d0**2
      Q2FIX=1
      IF(Q2.LE.Q2MIN)THEN
       Q2FIX=0
      ENDIF

*   Initialization for ALPHAS and RUNM (as in hdecay)
*   The bottom quark pole mass MBP is set in INIT and can be changed
*   only there (changing its running mass MB above has no effect
*   on MBP, since one would have to compute alpha_s(MB) first)

      MC0=MC
      MB0=MBP
      MT0=MT
      N0=5
      NLOOP=2
      NBER=18
      ACC=1.D-8
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      CALL ALSINI(ACC)
      CALL BERNINI(NBER)

*    g1,g2  and sin(theta)^2 in the on-shell scheme in terms of 
*    GF, MZ(pole) and MW(pole)

      g2=4d0*DSQRT(2d0)*GF*MW**2
      g1=4d0*DSQRT(2d0)*GF*(MZ**2-MW**2)
      S2TW=1d0-(MW/MZ)**2

      RETURN

 999  WRITE(0,1)"READ ERROR ON LINE:", NLINE
      WRITE(0,*)CHINL(1:80)
      STOP 1

 1    FORMAT(A)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL) 

*********************************************************************      
*   Subroutine writing all the results in the the output file.
*********************************************************************      
 
      IMPLICIT NONE 
 
      INTEGER NBIN,IFAIL,I,J,IHIGGS,IGMU
      INTEGER NSUSY,NGUT,NMES,IMAX
      PARAMETER (NSUSY=14,NGUT=21,NMES=15,IMAX=100)
 
      DOUBLE PRECISION RES(IMAX),PAR(*),PROB(*)
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,Q2
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
      DOUBLE PRECISION HBOTGUT,HTAUGUT,M1GUT,M2GUT,M3GUT,ALGUT,AKGUT
      DOUBLE PRECISION ATGUT,ABGUT,ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT
      DOUBLE PRECISION MQ3GUT,MU3GUT,MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT
      DOUBLE PRECISION ME3GUT,MLGUT,MEGUT,M0,M12,A0,SIGMU
      DOUBLE PRECISION M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX
      DOUBLE PRECISION Xf,sigmaV,xx(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigma0,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION MHUQ,MHDQ,MSQ,LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION SIG(5,11),SIGT(5,6)
      DOUBLE PRECISION M0CEN,M12CEN,TBCEN,A0CEN,M1CEN,M2CEN
      DOUBLE PRECISION M3CEN,MHDCEN,MHUCEN,LCEN,ALCEN,AKCEN
      DOUBLE PRECISION KCEN,MUCEN,XCEN,XDEV,X,M0DEV,M12DEV
      DOUBLE PRECISION TBDEV,A0DEV,M1DEV,M2DEV,M3DEV,MHDDEV
      DOUBLE PRECISION MHUDEV,LDEV,ALDEV,AKDEV,KDEV,MUDEV
*
      DOUBLE PRECISION chartot2(2),chartot(2),chartotmulti(2)
      DOUBLE PRECISION brcharst1(2),brcharst2(2),brcharsb1(2),
     .          brcharsb2(2),brcharsupl(2),brcharsupr(2),
     .          brcharsdownl(2),
     .          brcharsdownr(2),brcharsnel(2),brcharsn1(2),
     .          brcharsn2(2),brcharsell(2),brcharstau1(2),brcharselr(2),
     .          brcharstau2(2),brcharhcneut(2,5),brcharwneut(2,5),
     .          brcharzchic,brcharHchic(3),brcharAchic(2)
      DOUBLE PRECISION brntaunut(2,5),brnelnue(2,5),brnmunumu(2,5),
     .          brnupdb(2,5),brnchsb(2,5),brntopbb(2,5),
     .          brglupdb(2),brglchsb(2),brgltopbb(2),
     .          brchee,brchmumu,brchtautau,brchnene,
     .          brchnmunmu,brchntauntau,brchupup,brchdodo,
     .          brchchch,brchstst,brchtoptop,brchbotbot
*
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttot3(5),neuttotrad(5)
      DOUBLE PRECISION brneutst1(5),brneutst2(5),brneutsb1(5),
     .          brneutsb2(5),
     .          brneutsupl(5),brneutsupr(5),brneutsdownl(5),
     .          brneutsdownr(5),brneutsnel(5),brneutsn1(5),
     .          brneutsn2(5),brneutsell(5),brneutselr(5),
     .          brneutstau1(5),brneutstau2(5),brneutwchar(5,2),
     .          brneuthcchar(5,2),brneutzneut(5,5),
     .          brneutHneut(5,5,3),brneutAneut(5,5,2),brnraddec(5,5)
      DOUBLE PRECISION brneutup(5,5),brneutdow(5,5),brneutch(5,5),
     .          brneutst(5,5),brneutbot(5,5),brneuttop(5,5),
     .          brneutel(5,5),brneutmu(5,5),brneuttau(5,5),
     .          brneutnue(5,5),brneutnumu(5,5),brneutnutau(5,5),
     .          brchubd(5,2),brchcbs(5,2),brchtbb(5,2),brchelne(5,2),
     .          brchmunmu(5,2),brchtauntau(5,2),brglup(5),brgldo(5),
     .          brglch(5),brglst(5),brgltop(5),brglbot(5)
*
      DOUBLE PRECISION selltot,selltot2,selltot3,
     . selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     . sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .          brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .          brsntauneut(5),brsntauchar(2),brsntau1hcstau(2),
     .          brsntau1wstau(2),brstau1neut(5),brstau2neut(5),
     .          brstau1char(2),
     .          brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .          brstau1wsn(2),brstau2wsn(2),brstau2H(3),brstau2A(2),
     .          brstau2ztau
      DOUBLE PRECISION brselrstau,brselrstaustar,brstau2stau1star,
     .    brstau2stau1,brstau2stau1nn,
     .    brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .    brsellstau1star,brsellstau1,brsellstau1nutau,
     .    brsnestau1star,brsnestau1,brsnestau1nutau
*
      DOUBLE PRECISION supltot2,suprtot2,sdowltot2,sdowrtot2
      DOUBLE PRECISION brsuplnup(5),brsuplcdow(2),brsuplglui,
     .          brsuprnup(5),brsuprcdow(2),brsuprglui,
     .          brsdowlndow(5),brsdowlchup(2),brsdowlglui,
     .          brsdowrndow(5),brsdowrchup(2),brsdowrglui
*
      DOUBLE PRECISION 
     .          stoptot(2),stoptot2(2),stoptotmulti(2),stoptotrad(2)
      DOUBLE PRECISION brst1neutt(5),brst2neutt(5),brst1charb(2),
     .          brst2charb(2),brst1hcsb(2),brst2hcsb(2),brst1wsb(2),
     .          brst2wsb(2),brst1glui,brst2glui,brst2H(3),brst2A(2),
     .          brst2ztop,brgamma,brgammaup,brgammagluino
      DOUBLE PRECISION brstopw(2,5),brstoph(2,5),brststau(2,2),
     .          brstsntau(2,2),brstsel(2,2),brstsnel(2),
     .          brstbsbst(2,2),brstbbsbt(2,2),brsttausbnu(2,2),
     .          brstelsbnu(2,2),brstupsbdow(2,2),
     .          brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .          brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      DOUBLE PRECISION sbottot(2),sbottot2(2),sbottotmulti(2)
      DOUBLE PRECISION brsb1neutt(5),brsb2neutt(5),brsb1chart(2),
     .          brsb2chart(2),brsb1hcst(2),brsb2hcst(2),
     .          brsb1glui,brsb2glui,brsb1wst(2),
     .          brsb2wst(2),brsb2H(3),brsb2A(2),brsb2zbot
      DOUBLE PRECISION  brsbstau(2,2),brsbsntau(2,2),brsbsel(2,2),
     .          brsbtstsb(2,2),brsbtbstb(2,2),brsbtaustnu(2,2),
     .          brsbelstnu(2,2),brsbupstdow(2,2),brsbsnel(2),
     .          brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .          brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .          brsb2sb1tautau
      DOUBLE PRECISION gluitot,gluitot2,gluitotmulti,gluitotrad
      DOUBLE PRECISION brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon(5)
      DOUBLE PRECISION brgoup(5),brgoch(5),brgodn(5),brgost(5),
     .         brgotp(5),
     .         brgobt(5),brgoud(2),brgocs(2),brgotb(2),brhcst1b,brwst1b
*************************************************************
*
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartotmulti
      COMMON/CHARGINO_BR_2BD/brcharst1,brcharst2,brcharsb1,
     .          brcharsb2,brcharsupl,brcharsupr,brcharsdownl,
     .          brcharsdownr,brcharsnel,brcharsn1,
     .          brcharsn2,brcharsell,brcharstau1,brcharselr,
     .          brcharstau2,brcharhcneut,brcharwneut,
     .          brcharzchic,brcharHchic,brcharAchic
      COMMON/CHARGINO_BR_3BD/brntaunut,brnelnue,brnmunumu,
     .          brnupdb,brnchsb,brntopbb,
     .          brglupdb,brglchsb,brgltopbb,
     .          brchee,brchmumu,brchtautau,brchnene,
     .          brchnmunmu,brchntauntau,brchupup,brchdodo,
     .          brchchch,brchstst,brchtoptop,brchbotbot
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttot3,neuttotrad
      COMMON/NEUTRALINO_BR_2BD/brneutst1,brneutst2,brneutsb1,brneutsb2,
     .         brneutsupl,brneutsupr,brneutsdownl,brneutsdownr,
     .         brneutsnel,brneutsn1,brneutsn2,brneutsell,brneutselr,
     .         brneutstau1,brneutstau2,brneutwchar,brneuthcchar,
     .         brneutzneut,brneutHneut,brneutAneut,brnraddec
      COMMON/NEUTRALINO_BR_3BD/brneutup,brneutdow,brneutch,brneutst,
     .         brneutbot,brneuttop,brneutel,brneutmu,brneuttau,
     .         brneutnue,brneutnumu,brneutnutau,brchubd,brchcbs, 
     .         brchtbb,brchelne,brchmunmu,brchtauntau,brglup,brgldo,
     .         brglch,brglst,brgltop,brglbot
*
      COMMON/SLEPTON_WIDTH/selltot,selltot2,selltot3,
     .selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     .sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
      COMMON/SLEPTON_BR_2BD/brsellneute,brselrneute,brsellcharnue,
     .          brselrcharnue,brsnellneut,brsnellchar,
     .          brsntauneut,brsntauchar,brsntau1hcstau,
     .          brsntau1wstau,brstau1neut,brstau2neut,brstau1char,
     .          brstau2char,brstau1hcsn,brstau2hcsn,
     .          brstau1wsn,brstau2wsn,brstau2H,brstau2A,brstau2ztau
      COMMON/SLEPTON_BR_3BD/brselrstau,brselrstaustar,brstau2stau1star,
     .   brstau2stau1,brstau2stau1nn,
     .   brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .   brsellstau1star,brsellstau1,brsellstau1nutau,
     .   brsnestau1star,brsnestau1,brsnestau1nutau
*
      COMMON/SQUARK_WIDTH/supltot2,suprtot2,sdowltot2,sdowrtot2
      COMMON/SQUARK_BR_2BD/brsuplnup,brsuplcdow,brsuplglui,
     .          brsuprnup,brsuprcdow,brsuprglui,
     .          brsdowlndow,brsdowlchup,brsdowlglui,
     .          brsdowrndow,brsdowrchup,brsdowrglui
*
      COMMON/STOP_WIDTH/stoptot,stoptot2,stoptotmulti,stoptotrad
      COMMON/STOP_BR_2BD/brst1neutt,brst2neutt,brst1charb,
     .          brst2charb,brst1hcsb,brst2hcsb,brst1wsb,
     .          brst2wsb,brst1glui,brst2glui,brst2H,brst2A,
     .          brst2ztop,brgamma,brgammaup,brgammagluino
      COMMON/STOP_BR_3BD/brstopw,brstoph,brststau,
     .          brstsntau,brstsel,brstsnel,
     .          brstbsbst,brstbbsbt,brsttausbnu,
     .          brstelsbnu,brstupsbdow,
     .          brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .          brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      COMMON/SBOTTOM_WIDTH/sbottot,sbottot2,sbottotmulti
      COMMON/SBOTTOM_BR_2BD/brsb1neutt,brsb2neutt,brsb1chart,
     .          brsb2chart,brsb1hcst,brsb2hcst,
     .          brsb1glui,brsb2glui,brsb1wst,
     .          brsb2wst,brsb2H,brsb2A,brsb2zbot
      COMMON/SBOTTOM_BR_3BD/brsbstau,brsbsntau,brsbsel,
     .          brsbtstsb,brsbtbstb,brsbtaustnu,
     .          brsbelstnu,brsbupstdow,brsbsnel,
     .          brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .          brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .          brsb2sb1tautau
*
      COMMON/GLUINO_WIDTH/gluitot,gluitot2,gluitotmulti,gluitotrad
      COMMON/GLUINO_BR_2BD/brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon
      COMMON/GLUINO_BR_3BD/brgoup,brgoch,brgodn,brgost,brgotp,
     .         brgobt,brgoud,brgocs,brgotb,brhcst1b,brwst1b
*
      COMMON/SIGMU/SIGMU
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/INPPAR/M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP,
     . XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     . MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,xx,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigma0,csPsi,csNsi,csPsd,csNsd
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/LHCSIG/SIG,SIGT
      COMMON/MCMCPAR/M0CEN,M12CEN,TBCEN,A0CEN,M1CEN,M2CEN,M3CEN,
     . MHDCEN,MHUCEN,LCEN,ALCEN,AKCEN,KCEN,MUCEN,XCEN,XDEV,X,
     . M0DEV,M12DEV,TBDEV,A0DEV,M1DEV,M2DEV,M3DEV,
     . MHDDEV,MHUDEV,LDEV,ALDEV,AKDEV,KDEV,MUDEV,IHIGGS,IGMU

      RES(1)=PAR(1)
      RES(2)=PAR(2)
      RES(3)=PAR(3)
      RES(4)=PAR(4)
      RES(5)=M0
      RES(6)=M12
      RES(7)=A0
      RES(8)=ALINP
      RES(9)=AKINP
      RES(10)=MHUGUT
      RES(11)=MHDGUT
      RES(12)=MSGUT

      DO I=1,3
       RES(4+9*I)=SMASS(I)
       DO J=1,3
        RES(4+9*I+J)=SCOMP(I,J)**2
       ENDDO
       RES(8+9*I)=SIG(I,7)
       RES(9+9*I)=SIG(I,8)
       RES(10+9*I)=SIG(I,9)
       RES(11+9*I)=SIG(I,10)
       RES(12+9*I)=SIG(I,11)
      ENDDO
      DO I=1,2
       RES(38+2*I)=PMASS(I)
       RES(39+2*I)=PCOMP(I,2)**2
      ENDDO
      RES(44)=CMASS
      RES(45)=DABS(MNEU(1))
      RES(46)=NEU(1,1)**2
      RES(47)=NEU(1,2)**2
      RES(48)=NEU(1,3)**2+NEU(1,4)**2
      RES(49)=NEU(1,5)**2
      RES(50)=OMG
      RES(51)=CSPSI
      RES(52)=MGL
      RES(53)=MIN(MUR,MDR,MUL,MDL)
      RES(54)=MST1
      RES(55)=MSL1
      RES(56)=PROB(37)
      RES(57)=FTGUT(NGUT+1)

      IF(IFAIL.EQ.0)WRITE(6,11)(RES(I),I=1,57)
      
 11   FORMAT(57E16.8)

      END


      SUBROUTINE ERROR(TOT,NTOT,NFAIL)

*********************************************************************      
*   Subroutine for the error file. It contains a summary of the scan:
*   Number of points that passed/failed the tests
*   and ranges for scanned parameters that passed the tests
*********************************************************************      
 
      IMPLICIT NONE

      INTEGER I,S,TOT,NTOT,NFAIL(*)
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER AKFLAG,ALFLAG,OMGFLAG,MAFLAG

      DOUBLE PRECISION M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN,LN,LNN
      DOUBLE PRECISION AKN,AKNN,ALN,ALNN,MHDN,MHDNN,MHUN,MHUNN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN,KN,KNN,MUN,MUNN

      COMMON/BOUNDS/M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN,MHUN,MHUNN,
     . LN,LNN,ALN,ALNN,AKN,AKNN,KN,KNN,MUN,MUNN
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG

      WRITE(0,*)
      WRITE(0,20)"Number of points:                       "
      WRITE(0,*)
      WRITE(0,20)"  scanned                               ",NTOT
      WRITE(0,20)"  l, k, tan(beta) or mu=0               ",NFAIL(9)
      WRITE(0,20)"  no electroweak symmetry breaking      ",NFAIL(16)
      S=0
      DO I=1,7
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,20)"  with mh1^2 or ma1^2 or mhc^2 < 0      ",S
      WRITE(0,20)"  with m_sfermion^2 < 0                 ",NFAIL(8)
      WRITE(0,20)"  violating phenomenological constraints",NFAIL(10)
      S=NFAIL(11)+NFAIL(12)+NFAIL(13)
      WRITE(0,20)"  RGE integration problem               ",S
      S=NFAIL(14)+NFAIL(15)
      WRITE(0,20)"  convergence problem                   ",S
      WRITE(0,20)"  violating scenario                    ",NFAIL(17)
      WRITE(0,*)
      WRITE(0,20)"Remaining good points                   ",TOT
      IF(TOT.GT.0)THEN
       WRITE(0,*)
       WRITE(0,20)"Parameter ranges for good points:       "
       WRITE(0,*)
       WRITE(0,30)" M0: ",M0N,M0NN
       WRITE(0,30)" M12: ",M12N,M12NN
       WRITE(0,30)" TANB: ",TBN,TBNN
       WRITE(0,30)" A0: ",A0N,A0NN
       IF(M1FLAG.EQ.1)THEN
        WRITE(0,30)"M1: ",M1N,M1NN
       ENDIF
       IF(M2FLAG.EQ.1)THEN
        WRITE(0,30)"M2: ",M2N,M2NN
       ENDIF
       IF(M3FLAG.EQ.1)THEN
        WRITE(0,30)"M3: ",M3N,M3NN
       ENDIF
       IF(MHDFLAG.EQ.1)THEN
        WRITE(0,30)"MHD^2: ",MHDN,MHDNN
       ENDIF
       IF(MHUFLAG.EQ.1)THEN
        WRITE(0,30)"MHU^2: ",MHUN,MHUNN
       ENDIF
       WRITE(0,30)" LAMBDA: ",LN,LNN
       IF(ALFLAG.EQ.1)THEN
        WRITE(0,30)"ALAMBDA: ",ALN,ALNN
       ENDIF
       IF(AKFLAG.EQ.1)THEN
        WRITE(0,30)"AKAPPA: ",AKN,AKNN
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        WRITE(0,30)"KAPPA: ",KN,KNN
        WRITE(0,30)"MUEFF: ",MUN,MUNN
       ENDIF
      ENDIF

 20   FORMAT(A40,I10)
 30   FORMAT(A15,2E15.4)

      END
