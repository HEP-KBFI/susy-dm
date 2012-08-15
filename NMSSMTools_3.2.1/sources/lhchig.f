      SUBROUTINE LHCHIG(PAR,PROB)

*   Subroutine to check LHC constraints

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION PAR(*),PROB(*),M(5),D(6),SIG(5,11),SIGT(5,6),S,SM
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION SIGMA_TTH,SIGMA_BBH,SIGMA_GGF,SIGMA_VBF,SIGMA_WH
      DOUBLE PRECISION LHC_TAUTAU,LHC_LL,LHC_BB,LHC_ZZ,LHC_WW,LHC_GG
      DOUBLE PRECISION ATL_GG,LHC_TBH
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH

      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     .      BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     .       HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     .       HCBRSUSY,HCWIDTH
      COMMON/LHCSIG/SIG,SIGT

* Loop over H1, H2, H3

      DO I=1,3

       M(I)=SMASS(I)

       DO J=1,11
        SIG(I,J)=0d0
       ENDDO
       DO J=1,6
        SIGT(I,J)=0d0
       ENDDO

       S=CU(I)**2*SIGMA_TTH(M(I))+CD(I)**2*SIGMA_BBH(M(I))
     .  +CJ(I)**2*SIGMA_GGF(M(I))+CV(I)**2*SIGMA_VBF(M(I))
     .  +CV(I)**2*SIGMA_WH(M(I))

       SM=SIGMA_TTH(M(I))+SIGMA_BBH(M(I))
     .  +SIGMA_GGF(M(I))+SIGMA_VBF(M(I))
     .  +SIGMA_WH(M(I))

       CALL HDECAY(M(I),BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .      BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   H -> tautau

       SIG(I,1)=S*BRLL(I)
       IF(SM*BRLLSM.NE.0d0)SIG(I,2)=S/SM*BRLL(I)/BRLLSM
       IF(BRLLSM.NE.0d0)SIG(I,7)=CV(I)**2*BRLL(I)/BRLLSM

*   H -> bb

       IF(BRBBSM.NE.0d0)SIG(I,3)=CV(I)**2*BRBB(I)/BRBBSM

*   H -> ZZ

       IF(SM*BRZZSM.NE.0d0)SIG(I,4)=S/SM*BRZZ(I)/BRZZSM
       IF(BRZZSM.NE.0d0)SIG(I,8)=CJ(I)**2*BRZZ(I)/BRZZSM

*   H -> WW

       IF(SM*BRWWSM.NE.0d0)SIG(I,5)=S/SM*BRWW(I)/BRWWSM
       IF(BRWWSM.NE.0d0)SIG(I,9)=CJ(I)**2*BRWW(I)/BRWWSM

*   H -> gammagamma

       IF(SM*BRGGSM.NE.0d0)SIG(I,6)=S/SM*BRGG(I)/BRGGSM
       IF(BRGGSM.NE.0d0)SIG(I,10)=CJ(I)**2*BRGG(I)/BRGGSM
       IF(BRGGSM.NE.0d0)SIG(I,11)=CV(I)**2*BRGG(I)/BRGGSM

      ENDDO

* Loop over A1, A2

      DO I=4,5

       M(I)=PMASS(I-3)

       DO J=1,11
        SIG(I,J)=0d0
       ENDDO
       DO J=1,6
        SIGT(I,J)=0d0
       ENDDO

       S=CU(I)**2*SIGMA_TTH(M(I))+CD(I)**2*SIGMA_BBH(M(I))
     .  +CJ(I)**2*SIGMA_GGF(M(I))

       SM=SIGMA_TTH(M(I))+SIGMA_BBH(M(I))
     .  +SIGMA_GGF(M(I))

       CALL ADECAY(M(I),BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .      BRBBSM,BRTTSM,BRGGSM,BRZGSM)

*   A -> tautau

       SIG(I,1)=S*BRLL(I)
       IF(SM*BRLLSM.NE.0d0)SIG(I,2)=SIG(I,1)/SM/BRLLSM

*   A -> gammagamma

       IF(SM*BRGGSM.NE.0d0)SIG(I,6)=S/SM*BRGG(I)/BRGGSM

      ENDDO

*   Modify SIG(IJ,K) if |M(I)-M(J)| < DM(K)*(M(I)+M(J))/2 :

      D(1)=21d-2
      D(2)=20d-2
      D(3)=10d-2
      D(4)=2d-2
      D(5)=20d-2
      D(6)=2d-2
c      D(2)=1d-6
c      D(3)=1d-6
c      D(4)=1d-6
c      D(5)=1d-6
c      D(6)=1d-6

      CALL COMBINE(M,D,SIG,SIGT,6)

* Compare to LHC bounds

      PROB(45)=DDIM(SIGT(1,1)/LHC_TAUTAU(M(1)),1d0)
      PROB(46)=DDIM(SIGT(1,2)/LHC_LL(M(1)),1d0)
      PROB(47)=DDIM(SIGT(1,3)/LHC_BB(M(1)),1d0)
      PROB(48)=DDIM(SIGT(1,4)/LHC_ZZ(M(1)),1d0)
      PROB(49)=DDIM(SIGT(1,5)/LHC_WW(M(1)),1d0)
      PROB(50)=DDIM(SIGT(1,6)/LHC_GG(M(1)),1d0)
     .        +DDIM(SIGT(1,6)/ATL_GG(M(1)),1d0)
      DO I=2,5
        IF(SIGT(I,1).NE.SIGT(I-1,1))
     .    PROB(45)=PROB(45)+DDIM(SIGT(I,1)/LHC_TAUTAU(M(I)),1d0)
        IF(SIGT(I,2).NE.SIGT(I-1,2))
     .    PROB(46)=PROB(46)+DDIM(SIGT(I,2)/LHC_LL(M(I)),1d0)
        IF(SIGT(I,3).NE.SIGT(I-1,3))
     .    PROB(47)=PROB(47)+DDIM(SIGT(I,3)/LHC_BB(M(I)),1d0)
        IF(SIGT(I,4).NE.SIGT(I-1,4))
     .    PROB(48)=PROB(48)+DDIM(SIGT(I,4)/LHC_ZZ(M(I)),1d0)
        IF(SIGT(I,5).NE.SIGT(I-1,5))
     .    PROB(49)=PROB(49)+DDIM(SIGT(I,5)/LHC_WW(M(I)),1d0)
        IF(SIGT(I,6).NE.SIGT(I-1,6))
     .    PROB(50)=PROB(50)+DDIM(SIGT(I,6)/LHC_GG(M(I)),1d0)
     .            +DDIM(SIGT(I,6)/ATL_GG(M(I)),1d0)
      ENDDO

* Bound on Br(t->bH+)*BR(H+->tau nu)

      PROB(51)=DDIM(brtopbh*HCBRL/LHC_TBH(CMASS),1d0)

      END


      DOUBLE PRECISION FUNCTION SIGMA_GGF(M)

*   sigma(pp(gg)->H) from 1101.0593, table 5

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=50)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,95d0,100d0,105d0,110d0,115d0,120d0,125d0,130d0,135d0,
     .140d0,145d0,150d0,155d0,160d0,165d0,170d0,175d0,180d0,185d0,190d0,
     .195d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,290d0,
     .300d0,320d0,340d0,360d0,380d0,400d0,450d0,500d0,550d0,600d0,650d0,
     .700d0,750d0,800d0,850d0,900d0,950d0,1000d0/
      DATA Y/29.79d0,26.77d0,24.25d0,22.01d0,20.06d0,18.35d0,16.84d0,
     .15.51d0,14.32d0,13.26d0,12.31d0,11.45d0,10.67d0,9.94d0,9.21d0,
     .8.47d0,7.87d0,7.35d0,6.86d0,6.42d0,6.01d0,5.65d0,5.34d0,4.81d0,
     .4.36d0,3.97d0,3.65d0,3.37d0,3.11d0,2.89d0,2.71d0,2.55d0,2.42d0,
     .2.23d0,2.19d0,2.31d0,2.18d0,1.93d0,1.27d0,0.79d0,0.49d0,0.31d0,
     .0.20d0,0.13d0,0.08d0,0.06d0,0.04d0,0.03d0,0.02d0,0.01d0/

      IF(M.LE.X(2))THEN
       SIGMA_GGF=1d3*(Y(1)+(Y(2)-Y(1))*(M-X(1))/(X(2)-X(1)))
       RETURN
      ENDIF
      DO I=2,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        SIGMA_GGF=1d3*(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO
      SIGMA_GGF=0d0

      END


      DOUBLE PRECISION FUNCTION SIGMA_VBF(M)

*   sigma(pp(WW)->qqH) from 1101.0593, table 7

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=50)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,95d0,100d0,105d0,110d0,115d0,120d0,125d0,130d0,135d0,
     .140d0,145d0,150d0,155d0,160d0,165d0,170d0,175d0,180d0,185d0,190d0,
     .195d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,290d0,
     .300d0,320d0,340d0,360d0,380d0,400d0,450d0,500d0,550d0,600d0,650d0,
     .700d0,750d0,800d0,850d0,900d0,950d0,1000d0/
      DATA Y/1682d0,1598d0,1530d0,1445d0,1385d0,1312d0,1257d0,1193d0,
     .1144d0,1087d0,1042d0,992d0,951d0,907d0,869d0,842d0,808d0,772d0,
     .738d0,713d0,684d0,658d0,630d0,580d0,535d0,495d0,458d0,425d0,
     .395d0,368d0,343d0,320d0,298d0,260d0,227d0,200d0,180d0,161d0,
     .125d0,94.6d0,74.8d0,57.6d0,46.6d0,36.4d0,30.0d0,23.7d0,19.9d0,
     .15.9d0,13.6d0,11.0d0/

      IF(M.LE.X(2))THEN
       SIGMA_VBF=Y(1)+(Y(2)-Y(1))*(M-X(1))/(X(2)-X(1))
       RETURN
      ENDIF
      DO I=2,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        SIGMA_VBF=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO
      SIGMA_VBF=0d0

      END


      DOUBLE PRECISION FUNCTION SIGMA_WH(M)

*   sigma(pp->WH) from 1101.0593, table 13

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=33)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,95d0,100d0,105d0,110d0,115d0,120d0,125d0,130d0,135d0,
     .140d0,145d0,150d0,155d0,160d0,165d0,170d0,175d0,180d0,185d0,190d0,
     .195d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,290d0,
     .300d0/
      DATA Y/1.640d0,1.392d0,1.186d0,1.018d0,0.8754d0,0.7546d0,0.6561d0,
     .0.5729d0,0.5008d0,0.4390d0,0.3857d0,0.3406d0,0.3001d0,0.2646d0,
     .0.2291d0,0.2107d0,0.1883d0,0.1689d0,0.1521d0,0.1387d0,0.1253d0,
     .0.1138d0,0.1032d0,0.08557d0,0.07142d0,0.06006d0,0.05075d0,
     .0.04308d0,0.03674d0,0.03146d0,0.02700d0,0.02333d0,0.02018d0/

      IF(M.LE.X(2))THEN
       SIGMA_WH=1d3*(Y(1)+(Y(2)-Y(1))*(M-X(1))/(X(2)-X(1)))
       RETURN
      ENDIF
      DO I=2,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        SIGMA_WH=1d3*(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO
      SIGMA_WH=0d0

      END


      DOUBLE PRECISION FUNCTION SIGMA_TTH(M)

*   sigma(pp->ttH) from 1101.0593, table 13

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=33)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,95d0,100d0,105d0,110d0,115d0,120d0,125d0,130d0,135d0,
     .140d0,145d0,150d0,155d0,160d0,165d0,170d0,175d0,180d0,185d0,190d0,
     .195d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,290d0,
     .300d0/
      DATA Y/216.2d0,188.0d0,163.8d0,143.3d0,125.7d0,110.6d0,97.56d0,
     .86.34d0,76.58d0,68.10d0,60.72d0,54.35d0,48.69d0,43.74d0,39.42d0,
     .35.59d0,32.19d0,29.18d0,26.52d0,24.14d0,22.06d0,20.16d0,18.49d0,
     .15.62d0,13.30d0,11.43d0,9.873d0,8.593d0,7.524d0,6.636d0,5.889d0,
     .5.256d0,4.719d0/

      IF(M.LE.X(2))THEN
       SIGMA_TTH=Y(1)+(Y(2)-Y(1))*(M-X(1))/(X(2)-X(1))
       RETURN
      ENDIF
      DO I=2,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        SIGMA_TTH=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO
      SIGMA_TTH=0d0

      END


      DOUBLE PRECISION FUNCTION SIGMA_BBH(M)

*   sigma(pp->ttH) from 1101.0593, figure 22

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=52)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/70d0,80d0,90d0,95d0,100d0,105d0,110d0,115d0,120d0,125d0,
     .130d0,135d0,140d0,145d0,150d0,155d0,160d0,165d0,170d0,175d0,180d0,
     .185d0,190d0,195d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,
     .280d0,290d0,300d0,320d0,340d0,360d0,380d0,400d0,450d0,500d0,550d0,
     .600d0,650d0,700d0,750d0,800d0,850d0,900d0,950d0,1000d0/
      DATA Y/.878d3,.758d3,.533d3,.446d3,.375d3,.324d3,.276d3,.235d3,
     ..203d3,.177d3,.153d3,.133d3,.117d3,.103d3,.911d2,.811d2,.723d2,
     ..645d2,.571d2,.512d2,.456d2,.413d2,.370d2,.333d2,.299d2,.245d2,
     ..204d2,.170d2,.141d2,.119d2,.101d2,.853d1,.722d1,.622d1,.526d1,
     ..389d1,.297d1,.227d1,.174d1,.137d1,.748d0,.437d0,.262d0,.161d0,
     ..103d0,.660d-1,.439d-1,.293d-1,.200d-1,.138d-1,.961d-2,.703d-2/

      IF(M.LE.X(2))THEN
       SIGMA_BBH=Y(1)+(Y(2)-Y(1))*(M-X(1))/(X(2)-X(1))
       RETURN
      ENDIF
      DO I=2,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        SIGMA_BBH=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO
      SIGMA_BBH=0d0

      END


      DOUBLE PRECISION FUNCTION LHC_TBH(M)

* ATLAS constraints on BR(t->bH+)*BR(H+->taunu), ATLAS-CONF-2011-151 tab.5

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=8)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,100d0,110d0,120d0,130d0,140d0,150d0,160d0/ 
      DATA Y/.104d0,.098d0,.095d0,.077d0,.066d0,.071d0,.052d0,.141d0/ 

      LHC_TBH=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_TBH=1d3*(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LHC_TAUTAU(M)

* Parametrisation of constraints on sigma(pp->H/A->tautau)
* such that the constraints on tanbeta vs mA constraint in the MSSM
* limit from CMS PAS HIG-11-029 tab. 4 are reproduced

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=14)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,100d0,120d0,130d0,140d0,160d0,180d0,200d0,
     .250d0,300d0,350d0,400d0,450d0,500d0/ 
      DATA Y/31.08d0,20.27d0,8.64d0,3.91d0,4.68d0,1.72d0,.859d0,
     ..5415d0,.529d0,.5355d0,.516d0,.4725d0,.3855d0,.317d0/

      LHC_TAUTAU=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_TAUTAU=1d3*(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LHC_LL(M)

* CMS constraints on sigma(pp->H->tautau), CMS PAS HIG-11-029 tab. 5

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=8)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/110d0,115d0,120d0,125d0,130d0,135d0,140d0,145d0/ 
      DATA Y/3.48d0,2.86d0,3.15d0,3.55d0,4.03d0,4.55d0,4.89d0,6.28d0/ 

      LHC_LL=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_LL=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LHC_BB(M)

* CMS constraints on sigma(pp->WH->Wbb), CMS PAS HIG-11-031 tab. 14

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=6)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/110d0,115d0,120d0,125d0,130d0,135d0/ 
      DATA Y/3.14d0,5.18d0,4.38d0,5.72d0,9.00d0,7.53d0/ 

      LHC_BB=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_BB=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LHC_ZZ(M)

* LHC constraints on sigma(pp->H->ZZ), CMS PAS HIG-11-032 fig. 3

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=131)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/1.11d2,1.12d2,1.13d2,1.14d2,1.15d2,1.16d2,1.17d2,1.18d2,
     .1.20d2,1.21d2,1.22d2,1.23d2,1.24d2,1.26d2,1.28d2,1.30d2,1.31d2,
     .1.32d2,1.33d2,1.34d2,1.35d2,1.37d2,1.41d2,1.43d2,1.46d2,1.48d2,
     .1.51d2,1.53d2,1.55d2,1.57d2,1.58d2,1.59d2,1.60d2,1.60d2,1.61d2,
     .1.62d2,1.63d2,1.65d2,1.66d2,1.66d2,1.68d2,1.70d2,1.71d2,1.73d2,
     .1.74d2,1.75d2,1.76d2,1.78d2,1.79d2,1.80d2,1.80d2,1.81d2,1.82d2,
     .1.83d2,1.83d2,1.84d2,1.85d2,1.86d2,1.89d2,1.90d2,1.91d2,1.91d2,
     .1.92d2,1.94d2,1.96d2,1.98d2,2.01d2,2.02d2,2.05d2,2.07d2,2.09d2,
     .2.10d2,2.13d2,2.15d2,2.16d2,2.18d2,2.23d2,2.25d2,2.28d2,2.30d2,
     .2.35d2,2.38d2,2.42d2,2.47d2,2.49d2,2.50d2,2.53d2,2.55d2,2.58d2,
     .2.60d2,2.62d2,2.66d2,2.70d2,2.72d2,2.75d2,2.76d2,2.78d2,2.80d2,
     .2.86d2,2.90d2,2.93d2,2.96d2,3.00d2,3.03d2,3.07d2,3.10d2,3.17d2,
     .3.21d2,3.29d2,3.33d2,3.36d2,3.40d2,3.42d2,3.46d2,3.49d2,3.51d2,
     .3.60d2,3.69d2,3.81d2,3.91d2,4.01d2,4.20d2,4.41d2,4.60d2,4.80d2,
     .5.01d2,5.21d2,5.41d2,5.59d2,5.82d2,6.00d2/
      DATA Y/5.53d0,5.13d0,5.58d0,5.63d0,6.23d0,6.96d0,7.02d0,6.18d0,
     .5.30d0,4.36d0,3.74d0,3.13d0,2.57d0,2.36d0,2.15d0,1.89d0,1.70d0,
     .1.47d0,1.26d0,1.05d0,9.03d-1,7.95d-1,7.49d-1,8.08d-1,7.18d-1,
     .6.42d-1,5.75d-1,5.24d-1,5.19d-1,6.10d-1,7.62d-1,9.34d-1,1.14d0,
     .1.47d0,1.77d0,2.17d0,2.55d0,3.08d0,3.49d0,3.87d0,3.87d0,3.41d0,
     .2.87d0,2.44d0,2.13d0,1.86d0,1.61d0,1.41d0,1.23d0,1.06d0,8.95d-1,
     .7.75d-1,6.70d-1,5.90d-1,5.10d-1,4.42d-1,3.82d-1,3.36d-1,3.31d-1,
     .4.06d-1,4.77d-1,5.80d-1,6.37d-1,6.99d-1,7.24d-1,6.65d-1,5.70d-1,
     .4.89d-1,4.23d-1,4.89d-1,5.28d-1,5.61d-1,5.19d-1,4.61d-1,3.92d-1,
     .3.39d-1,3.01d-1,3.22d-1,3.79d-1,4.38d-1,4.89d-1,5.56d-1,6.53d-1,
     .5.61d-1,5.15d-1,4.61d-1,4.85d-1,5.51d-1,6.48d-1,6.76d-1,5.75d-1,
     .5.15d-1,5.80d-1,5.42d-1,5.46d-1,6.00d-1,6.21d-1,6.82d-1,6.16d-1,
     .6.99d-1,6.10d-1,5.42d-1,5.61d-1,6.42d-1,7.49d-1,8.73d-1,8.88d-1,
     .9.34d-1,8.95d-1,8.22d-1,7.42d-1,6.48d-1,5.80d-1,5.19d-1,4.57d-1,
     .4.13d-1,4.16d-1,3.82d-1,3.76d-1,3.36d-1,3.92d-1,4.09d-1,5.75d-1,
     .6.76d-1,6.42d-1,8.51d-1,7.68d-1,9.83d-1,9.03d-1,9.34d-1,1.15d0/

      LHC_ZZ=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_ZZ=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LHC_WW(M)

* LHC constraints on sigma(pp->H->WW), CMS PAS HIG-11-032 fig. 3

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=94)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/1.11d2,1.12d2,1.13d2,1.15d2,1.16d2,1.17d2,1.18d2,1.19d2,
     .1.20d2,1.22d2,1.23d2,1.25d2,1.25d2,1.27d2,1.28d2,1.29d2,1.31d2,
     .1.33d2,1.35d2,1.37d2,1.39d2,1.43d2,1.45d2,1.46d2,1.50d2,1.53d2,
     .1.55d2,1.55d2,1.55d2,1.58d2,1.62d2,1.66d2,1.71d2,1.74d2,1.78d2,
     .1.81d2,1.85d2,1.89d2,1.93d2,1.96d2,1.99d2,2.02d2,2.04d2,2.06d2,
     .2.09d2,2.15d2,2.19d2,2.25d2,2.31d2,2.37d2,2.42d2,2.48d2,2.55d2,
     .2.62d2,2.69d2,2.76d2,2.83d2,2.90d2,2.98d2,3.02d2,3.05d2,3.12d2,
     .3.20d2,3.29d2,3.37d2,3.46d2,3.56d2,3.65d2,3.74d2,3.83d2,3.92d2,
     .4.02d2,4.11d2,4.20d2,4.25d2,4.28d2,4.34d2,4.37d2,4.42d2,4.51d2,
     .4.59d2,4.73d2,4.84d2,4.95d2,5.05d2,5.16d2,5.25d2,5.34d2,5.45d2,
     .5.55d2,5.65d2,5.75d2,5.86d2,5.99d2/
      DATA Y/7.46d0,6.51d0,5.68d0,4.91d0,4.36d0,3.77d0,3.24d0,2.80d0,
     .2.47d0,2.23d0,1.98d0,1.70d0,1.47d0,1.28d0,1.11d0,9.58d-1,8.51d-1,
     .8.15d-1,7.81d-1,6.99d-1,6.16d-1,5.70d-1,5.10d-1,4.49d-1,4.23d-1,
     .4.06d-1,3.54d-1,3.06d-1,2.63d-1,2.41d-1,2.31d-1,2.45d-1,2.61d-1,
     .2.74d-1,2.54d-1,2.79d-1,3.14d-1,3.48d-1,3.79d-1,4.27d-1,4.73d-1,
     .4.85d-1,4.16d-1,4.98d-1,5.37d-1,5.70d-1,6.05d-1,6.10d-1,6.31d-1,
     .6.70d-1,7.05d-1,7.42d-1,7.75d-1,8.22d-1,8.58d-1,9.03d-1,9.50d-1,
     .9.75d-1,1.04d0,9.67d-1,8.36d-1,8.22d-1,8.88d-1,8.95d-1,9.03d-1,
     .9.18d-1,9.18d-1,9.58d-1,1.01d0,1.08d0,1.17d0,1.26d0,1.36d0,1.49d0,
     .1.29d0,1.13d0,9.67d-1,8.43d-1,7.55d-1,8.36d-1,9.34d-1,9.11d-1,
     .9.26d-1,1.02d0,1.13d0,1.24d0,1.37d0,1.53d0,1.72d0,1.91d0,1.88d0,
     .1.65d0,1.64d0,1.86d0/

      LHC_WW=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_WW=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LHC_GG(M)

* LHC constraints on sigma(pp->H->gammagamma), CMS 1202.1487, fig. 2

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=391)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/1.1000d2,1.1030d2,1.1057d2,1.1072d2,1.1087d2,1.1108d2,
     .1.1118d2,1.1125d2,1.1134d2,1.1139d2,1.1144d2,1.1151d2,1.1154d2,
     .1.1159d2,1.1164d2,1.1172d2,1.1179d2,1.1188d2,1.1197d2,1.1209d2,
     .1.1221d2,1.1236d2,1.1257d2,1.1276d2,1.1292d2,1.1306d2,1.1319d2,
     .1.1328d2,1.1341d2,1.1352d2,1.1368d2,1.1387d2,1.1402d2,1.1413d2,
     .1.1425d2,1.1436d2,1.1447d2,1.1448d2,1.1463d2,1.1485d2,1.1506d2,
     .1.1525d2,1.1548d2,1.1565d2,1.1580d2,1.1595d2,1.1610d2,1.1628d2,
     .1.1645d2,1.1665d2,1.1684d2,1.1707d2,1.1722d2,1.1732d2,1.1755d2,
     .1.1782d2,1.1804d2,1.1820d2,1.1835d2,1.1849d2,1.1869d2,1.1892d2,
     .1.1914d2,1.1939d2,1.1959d2,1.1983d2,1.2000d2,1.2013d2,1.2025d2,
     .1.2026d2,1.2034d2,1.2043d2,1.2054d2,1.2063d2,1.2071d2,1.2079d2,
     .1.2086d2,1.2092d2,1.2099d2,1.2106d2,1.2111d2,1.2116d2,1.2121d2,
     .1.2121d2,1.2126d2,1.2131d2,1.2136d2,1.2141d2,1.2146d2,1.2151d2,
     .1.2156d2,1.2161d2,1.2166d2,1.2171d2,1.2176d2,1.2181d2,1.2186d2,
     .1.2186d2,1.2191d2,1.2196d2,1.2200d2,1.2205d2,1.2210d2,1.2215d2,
     .1.2220d2,1.2225d2,1.2230d2,1.2235d2,1.2240d2,1.2245d2,1.2250d2,
     .1.2253d2,1.2258d2,1.2265d2,1.2275d2,1.2284d2,1.2291d2,1.2301d2,
     .1.2315d2,1.2326d2,1.2342d2,1.2360d2,1.2377d2,1.2402d2,1.2422d2,
     .1.2437d2,1.2450d2,1.2460d2,1.2469d2,1.2476d2,1.2484d2,1.2489d2,
     .1.2496d2,1.2504d2,1.2509d2,1.2514d2,1.2519d2,1.2524d2,1.2529d2,
     .1.2529d2,1.2534d2,1.2539d2,1.2544d2,1.2549d2,1.2554d2,1.2559d2,
     .1.2564d2,1.2564d2,1.2569d2,1.2574d2,1.2576d2,1.2579d2,1.2584d2,
     .1.2584d2,1.2587d2,1.2589d2,1.2594d2,1.2595d2,1.2599d2,1.2604d2,
     .1.2607d2,1.2609d2,1.2614d2,1.2619d2,1.2621d2,1.2624d2,1.2629d2,
     .1.2631d2,1.2634d2,1.2634d2,1.2639d2,1.2644d2,1.2648d2,1.2650d2,
     .1.2656d2,1.2661d2,1.2665d2,1.2672d2,1.2676d2,1.2681d2,1.2686d2,
     .1.2692d2,1.2699d2,1.2699d2,1.2704d2,1.2710d2,1.2719d2,1.2725d2,
     .1.2733d2,1.2743d2,1.2753d2,1.2771d2,1.2791d2,1.2811d2,1.2830d2,
     .1.2846d2,1.2866d2,1.2880d2,1.2885d2,1.2908d2,1.2924d2,1.2938d2,
     .1.2958d2,1.2971d2,1.2985d2,1.2999d2,1.3014d2,1.3037d2,1.3060d2,
     .1.3080d2,1.3100d2,1.3120d2,1.3139d2,1.3160d2,1.3173d2,1.3182d2,
     .1.3187d2,1.3197d2,1.3209d2,1.3220d2,1.3232d2,1.3242d2,1.3254d2,
     .1.3264d2,1.3273d2,1.3285d2,1.3294d2,1.3304d2,1.3315d2,1.3321d2,
     .1.3326d2,1.3338d2,1.3354d2,1.3374d2,1.3393d2,1.3408d2,1.3422d2,
     .1.3433d2,1.3444d2,1.3453d2,1.3464d2,1.3473d2,1.3482d2,1.3488d2,
     .1.3491d2,1.3500d2,1.3508d2,1.3519d2,1.3528d2,1.3538d2,1.3548d2,
     .1.3555d2,1.3562d2,1.3568d2,1.3575d2,1.3580d2,1.3585d2,1.3585d2,
     .1.3590d2,1.3599d2,1.3606d2,1.3615d2,1.3625d2,1.3635d2,1.3650d2,
     .1.3670d2,1.3692d2,1.3715d2,1.3737d2,1.3752d2,1.3765d2,1.3775d2,
     .1.3785d2,1.3795d2,1.3803d2,1.3813d2,1.3822d2,1.3832d2,1.3842d2,
     .1.3852d2,1.3862d2,1.3862d2,1.3871d2,1.3881d2,1.3891d2,1.3902d2,
     .1.3912d2,1.3924d2,1.3936d2,1.3952d2,1.3974d2,1.3991d2,1.4019d2,
     .1.4044d2,1.4061d2,1.4064d2,1.4083d2,1.4108d2,1.4128d2,1.4148d2,
     .1.4171d2,1.4196d2,1.4215d2,1.4230d2,1.4245d2,1.4258d2,1.4267d2,
     .1.4280d2,1.4289d2,1.4300d2,1.4313d2,1.4318d2,1.4330d2,1.4350d2,
     .1.4372d2,1.4392d2,1.4415d2,1.4434d2,1.4451d2,1.4459d2,1.4469d2,
     .1.4477d2,1.4487d2,1.4495d2,1.4504d2,1.4514d2,1.4524d2,1.4532d2,
     .1.4534d2,1.4543d2,1.4554d2,1.4564d2,1.4575d2,1.4587d2,1.4597d2,
     .1.4609d2,1.4618d2,1.4629d2,1.4638d2,1.4653d2,1.4674d2,1.4691d2,
     .1.4709d2,1.4723d2,1.4737d2,1.4751d2,1.4763d2,1.4776d2,1.4787d2,
     .1.4797d2,1.4808d2,1.4818d2,1.4826d2,1.4831d2,1.4836d2,1.4841d2,
     .1.4848d2,1.4855d2,1.4861d2,1.4866d2,1.4870d2,1.4875d2,1.4880d2,
     .1.4885d2,1.4890d2,1.4895d2,1.4900d2,1.4905d2,1.4905d2,1.4910d2,
     .1.4915d2,1.4919d2,1.4923d2,1.4927d2,1.4930d2,1.4935d2,1.4940d2,
     .1.4943d2,1.4946d2,1.4950d2,1.4955d2,1.4960d2,1.4960d2,1.4965d2,
     .1.4970d2,1.4975d2,1.4979d2,1.4985d2,1.4987d2,1.4993d2,1.5000d2/
      DATA Y/8.5000d-1,8.6000d-1,8.9000d-1,9.2330d-1,9.6117d-1,1.0292d0,
     .1.0748d0,1.1127d0,1.1600d0,1.2027d0,1.2405d0,1.2863d0,1.3021d0,
     .1.3352d0,1.3731d0,1.4157d0,1.4583d0,1.5007d0,1.5388d0,1.5720d0,
     .1.6059d0,1.6349d0,1.6477d0,1.6335d0,1.6054d0,1.5767d0,1.5436d0,
     .1.5057d0,1.4725d0,1.4360d0,1.4084d0,1.3793d0,1.3496d0,1.3163d0,
     .1.2736d0,1.2358d0,1.2027d0,1.1932d0,1.1742d0,1.1906d0,1.1983d0,
     .1.1864d0,1.1742d0,1.1979d0,1.2358d0,1.2737d0,1.3026d0,1.3274d0,
     .1.3511d0,1.3766d0,1.3990d0,1.4205d0,1.4299d0,1.4392d0,1.4530d0,
     .1.4692d0,1.4716d0,1.4463d0,1.4162d0,1.3875d0,1.3669d0,1.3542d0,
     .1.3579d0,1.3782d0,1.3959d0,1.4181d0,1.4394d0,1.4732d0,1.5009d0,
     .1.5048d0,1.5483d0,1.5862d0,1.6335d0,1.6809d0,1.7235d0,1.7661d0,
     .1.8087d0,1.8521d0,1.8892d0,1.9318d0,1.9697d0,2.0170d0,2.0549d0,
     .2.0407d0,2.0928d0,2.1354d0,2.1733d0,2.2112d0,2.2491d0,2.2917d0,
     .2.3295d0,2.3769d0,2.4148d0,2.4621d0,2.5000d0,2.5473d0,2.5947d0,
     .2.5994d0,2.6420d0,2.6799d0,2.7273d0,2.7652d0,2.8125d0,2.8504d0,
     .2.8977d0,2.9356d0,2.9830d0,3.0208d0,3.0634d0,3.1013d0,3.1392d0,
     .3.1581d0,3.1866d0,3.2292d0,3.2718d0,3.3191d0,3.3570d0,3.3949d0,
     .3.4280d0,3.4637d0,3.4908d0,3.5133d0,3.5322d0,3.5288d0,3.5083d0,
     .3.4801d0,3.4422d0,3.4091d0,3.3665d0,3.3191d0,3.2812d0,3.2386d0,
     .3.1939d0,3.1439d0,3.1061d0,3.0682d0,3.0303d0,2.9924d0,2.9545d0,
     .2.9403d0,2.9072d0,2.8693d0,2.8220d0,2.7746d0,2.7273d0,2.6799d0,
     .2.6326d0,2.5947d0,2.5521d0,2.5095d0,2.4669d0,2.4290d0,2.3864d0,
     .2.3816d0,2.3438d0,2.3059d0,2.2680d0,2.2254d0,2.1828d0,2.1449d0,
     .2.1070d0,2.0644d0,2.0218d0,1.9792d0,1.9365d0,1.8987d0,1.8608d0,
     .1.8225d0,1.7803d0,1.8087d0,1.7377d0,1.6951d0,1.6477d0,1.6098d0,
     .1.5668d0,1.5199d0,1.4773d0,1.4347d0,1.3968d0,1.3589d0,1.3210d0,
     .1.2832d0,1.2405d0,1.2405d0,1.2027d0,1.1600d0,1.1174d0,1.0736d0,
     .1.0322d0,9.8011d-1,9.4697d-1,9.1502d-1,8.9531d-1,8.6687d-1,
     .8.3619d-1,8.0966d-1,7.9077d-1,7.8598d-1,7.8598d-1,7.9545d-1,
     .8.2962d-1,8.6174d-1,8.6174d-1,8.1913d-1,7.8125d-1,7.4272d-1,
     .7.1153d-1,6.9129d-1,6.8979d-1,7.0757d-1,7.3345d-1,7.5117d-1,
     .7.6705d-1,7.9946d-1,8.2883d-1,8.4754d-1,8.6174d-1,8.9962d-1,
     .9.3277d-1,9.7064d-1,1.0038d0,1.0511d0,1.0938d0,1.1316d0,1.1695d0,
     .1.2169d0,1.2547d0,1.2926d0,1.3301d0,1.3589d0,1.3731d0,1.4062d0,
     .1.4441d0,1.4730d0,1.4939d0,1.5199d0,1.5611d0,1.5956d0,1.6335d0,
     .1.6714d0,1.7188d0,1.7661d0,1.8134d0,1.8419d0,1.8561d0,1.9034d0,
     .1.9460d0,1.9934d0,2.0407d0,2.0881d0,2.1354d0,2.1780d0,2.2206d0,
     .2.2585d0,2.3011d0,2.3390d0,2.3769d0,2.3816d0,2.4148d0,2.4574d0,
     .2.5047d0,2.5473d0,2.5852d0,2.6326d0,2.6745d0,2.7083d0,2.7178d0,
     .2.7083d0,2.6851d0,2.6562d0,2.6231d0,2.5758d0,2.5284d0,2.4905d0,
     .2.4479d0,2.4006d0,2.3532d0,2.3059d0,2.2585d0,2.2112d0,2.1638d0,
     .2.1638d0,2.1165d0,2.0691d0,2.0218d0,1.9744d0,1.9366d0,1.8939d0,
     .1.8608d0,1.8230d0,1.7992d0,1.7803d0,1.7751d0,1.7519d0,1.7330d0,
     .1.7330d0,1.7140d0,1.6951d0,1.6951d0,1.7045d0,1.7140d0,1.7277d0,
     .1.7566d0,1.7865d0,1.8229d0,1.8561d0,1.8939d0,1.9348d0,1.9697d0,
     .2.0028d0,2.0360d0,2.0455d0,2.0687d0,2.0964d0,2.1117d0,2.1117d0,
     .2.0898d0,2.0595d0,2.0214d0,1.9839d0,1.9366d0,1.8939d0,1.8466d0,
     .1.8040d0,1.7566d0,1.7093d0,1.6619d0,1.6288d0,1.6146d0,1.5720d0,
     .1.5294d0,1.4915d0,1.4536d0,1.4110d0,1.3731d0,1.3305d0,1.2926d0,
     .1.2453d0,1.2074d0,1.1683d0,1.1553d0,1.1742d0,1.1979d0,1.2263d0,
     .1.2643d0,1.3068d0,1.3400d0,1.3826d0,1.4157d0,1.4536d0,1.5009d0,
     .1.5483d0,1.5909d0,1.6051d0,1.6383d0,1.6761d0,1.7188d0,1.7614d0,
     .1.7992d0,1.8371d0,1.8750d0,1.9129d0,1.9508d0,1.9886d0,2.0360d0,
     .2.0739d0,2.1212d0,2.1591d0,2.1638d0,2.2064d0,2.2538d0,2.2948d0,
     .2.3520d0,2.3958d0,2.4337d0,2.4716d0,2.5189d0,2.5606d0,2.6042d0,
     .2.6420d0,2.6799d0,2.7273d0,2.7225d0,2.7746d0,2.8220d0,2.8693d0,
     .2.9167d0,2.9640d0,3.0019d0,3.0445d0,3.0934d0/

      LHC_GG=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        LHC_GG=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION ATL_GG(M)

* LHC constraints on sigma(pp->H->gammagamma), ATLAS-CONF-2011-161 tab.7

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=41)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/110d0,111d0,112d0,113d0,114d0,115d0,116d0,117d0,118d0,
     .119d0,120d0,121d0,122d0,123d0,124d0,125d0,126d0,127d0,128d0,
     .129d0,130d0,131d0,132d0,133d0,134d0,135d0,136d0,137d0,138d0,
     .139d0,140d0,141d0,142d0,143d0,144d0,145d0,146d0,147d0,148d0,
     .149d0,150d0/
      DATA Y/1.94d0,1.67d0,1.32d0,1.01d0,0.86d0,0.93d0,1.28d0,1.83d0,
     .2.12d0,1.85d0,1.41d0,1.15d0,1.23d0,1.73d0,2.62d0,3.55d0,4.04d0,
     .3.82d0,3.06d0,2.19d0,1.65d0,1.43d0,1.34d0,1.28d0,1.14d0,0.98d0,
     .0.95d0,1.21d0,1.68d0,1.96d0,1.76d0,1.46d0,1.46d0,1.87d0,2.47d0,
     .2.88d0,2.85d0,2.54d0,2.25d0,2.02d0,1.92d0/

      ATL_GG=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        ATL_GG=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      SUBROUTINE COMBINE(M,D,S,ST,N)

      IMPLICIT NONE
      INTEGER N,I,J,K
      DOUBLE PRECISION M(5),D(*),ST(5,*),S(5,*),W

      DO J=1,5
       ST(J,1)=0d0
       DO K=1,5
        W=D(1)*M(K)*(1d0+(M(K)-130d0)/13d0)
        IF(M(K).LE.120d0) W=D(1)*M(K)*(1d0-10d0/13d0)
        IF(M(K).GE.160d0) W=D(1)*M(K)*(1d0+30d0/13d0)
        ST(J,1)=ST(J,1)+S(K,1)*DEXP(-(M(J)-M(K))**2/(2d0*W**2))
       ENDDO
      ENDDO

      DO I=2,N
       DO J=1,5
        IF(S(J,I).NE.0d0)THEN
         ST(J,I)=0d0
         DO K=1,5
          W=D(I)*M(K)
          ST(J,I)=ST(J,I)+S(K,I)*DEXP(-(M(J)-M(K))**2/(2d0*W**2))
         ENDDO
        ENDIF
       ENDDO
      ENDDO


      END
