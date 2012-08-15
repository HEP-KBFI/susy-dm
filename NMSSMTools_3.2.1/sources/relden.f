      SUBROUTINE RELDEN(PAR,PROB)

**********************************************************************
*   Subroutine for the computation of the dark matter relic density
*   PROB(30) =/= 0 excluded by WMAP,
*   PROB(30) = -1  LSP is not the lightest neutralino,
*   PROB(31) =/= 0 Higgs eff. self-couplings in Micromegas > 1.
*
**********************************************************************

      IMPLICIT NONE

      CHARACTER name*10,mess*20

      INTEGER NORD(5),HORD(3),NBIN,OMGFLAG,MAFLAG
      INTEGER sortOddParticles,err,i,j
      INTEGER nucleonAmplitudes

      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,PI
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MH(3),MA(2),MHC,XENON100
      DOUBLE PRECISION MGL,MCHA(2),UU(2,2),VV(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION SST,SSB,SSL,COSB,SINB,TANB
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION tab(250),OMG,OMGMIN,OMGMAX,Xf
      DOUBLE PRECISION sigmaV,x(100),dNdx(100),EMIN,LAM
      DOUBLE PRECISION sigmaPiN,sigma0,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION higgsPotent,darkOmega,calcSpectrum,zInterp
      DOUBLE PRECISION FeScLoop,LOPmass,NOFF,Nmass,SCcoeff
      DOUBLE PRECISION pA0(2),pA5(2),nA0(2),nA5(2),ffS0P(3),ffS0N(3)

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/HMO/MH,MA,MHC
      COMMON/SUSYSPEC/MGL,MCHA,UU,VV,MNEU,NEU
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigma0,csPsi,csNsi,csPsd,csNsd
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/LAM/LAM

      EXTERNAL FeScLoop,LOPmass

      DATA NOFF/-12345d0/
      DATA NORD/1,2,4,3,5/
      DATA HORD/2,1,3/

      IF (OMGFLAG.EQ.0) RETURN

*   Input parameters:

      PI=4d0*DATAN(1d0)
      TANB=PAR(3)
      COSB=1d0/DSQRT(1d0+TANB**2)
      SINB=TANB*COSB

      SST=DSQRT(1-CST**2)
      SSB=DSQRT(1-CSB**2)
      SSL=DSQRT(1-CSL**2)

      CALL assignValW('alfEMZ',ALEMMZ)
      CALL assignValW('alfSMZ',ALSMZ)
      CALL assignValW('MbMb',MB)
      CALL assignValW('Mtp',MT)

      CALL assignValW('At',PAR(12))
      CALL assignValW('Ab',PAR(13))
      CALL assignValW('Al',PAR(14))

      CALL assignValW('Lambda',LQ/DSQRT(ZHU*ZHD*ZS))
      CALL assignValW('Kappa',KQ/DSQRT(ZS**3))
      CALL assignValW('tb',PAR(3))
      CALL assignValW('aLmbd0',PAR(5))

      CALL assignValW('Mha',MA(1))
      CALL assignValW('Mhb',MA(2))
      CALL assignValW('MHc',MHC)
      DO i=1,3
       WRITE(name,fmt='(A2,I1)') 'Mh',i
       CALL assignValW(name,MH(i))
      DO j=1,3
       WRITE(name,fmt='(A2,I1,I1)') 'Zh',i,j
       CALL assignValW(name,SCOMP(i,HORD(j)))
      ENDDO
      ENDDO
      CALL assignValW('Pa11',PCOMP(1,1))
      CALL assignValW('Pa12',PCOMP(1,2))
      CALL assignValW('Pa21',PCOMP(2,1))
      CALL assignValW('Pa22',PCOMP(2,2))

      DO i=1,5
       WRITE(name,fmt='(A3,I1)') 'MNE',i
       CALL assignValW(name,MNEU(i))
       DO j=1,5
         WRITE(name,fmt='(A2,I1,I1)') 'Zn',i,j
         CALL assignValW(name,NEU(i,NORD(j)))
        ENDDO
      ENDDO

      CALL assignValW('MSl1',MSL1)
      CALL assignValW('MSl2',MSL2)
      CALL assignValW('Zl11',CSL)
      CALL assignValW('Zl12',SSL)
      CALL assignValW('Zl21',-SSL)
      CALL assignValW('Zl22',CSL)

      CALL assignValW('MSb1',MSB1)
      CALL assignValW('MSb2',MSB2)
      CALL assignValW('Zb11',CSB)
      CALL assignValW('Zb12',SSB)
      CALL assignValW('Zb21',-SSB)
      CALL assignValW('Zb22',CSB)

      CALL assignValW('MSt1',MST1)
      CALL assignValW('MSt2',MST2)
      CALL assignValW('Zt11',CST)
      CALL assignValW('Zt12',SST)
      CALL assignValW('Zt21',-SST)
      CALL assignValW('Zt22',CST)

      CALL assignValW('MSeL',MLL)
      CALL assignValW('MSeR',MLR)
      CALL assignValW('MSmL',MLL)
      CALL assignValW('MSmR',MLR)
      CALL assignValW('MSne',MNL)
      CALL assignValW('MSnm',MNL)
      CALL assignValW('MSnl',MSNT)
      CALL assignValW('MSuL',MUL)
      CALL assignValW('MSuR',MUR)
      CALL assignValW('MSdL',MDL)
      CALL assignValW('MSdR',MDR)
      CALL assignValW('MScL',MUL)
      CALL assignValW('MScR',MUR)
      CALL assignValW('MSsL',MDL)
      CALL assignValW('MSsR',MDR)
      CALL assignValW('MSG',MGL)

*   Improved Higgs potential

      LAM=higgsPotent()
      IF(LAM.EQ.-1d0 .OR. LAM.GE.1d0) THEN
        PROB(31)=LAM
      ENDIF

*   Sorting sparticles

      err=sortOddParticles(mess)
      IF(err.ne.0 .OR. mess.ne.'~o1') THEN
        OMG=-1d0
        PROB(30)=-1d0
        RETURN
      ENDIF

*   Computing relic density

      OMG=darkOmega(Xf,1,1.D-6)
      PROB(30)=DDIM(OMG/OMGMAX,1d0)-DDIM(1d0,OMG/OMGMIN)

      IF (OMGFLAG.EQ.1) RETURN
      IF (OMGFLAG.EQ.3) GOTO 1

*  Computing WIMP-Nucleon cross sections
*  Muq/Mdq=0.553d0, Msq/Mdq=18.9d0

      CALL getScalarFF(0.553d0,18.9d0,sigmaPiN,sigma0,ffS0P,ffS0N)
      CALL setProtonFF(ffS0P,NOFF,NOFF)
      CALL setNeutronFF(ffS0N,NOFF,NOFF)
      err=nucleonAmplitudes(FeScLoop,pA0,pA5,nA0,nA5)
      Nmass=0.939d0

      SCcoeff=4/PI*3.8937966D8*(Nmass*lopmass()/(Nmass+ lopmass()))**2
      csPsi=SCcoeff*pA0(1)**2
      csNsi=SCcoeff*nA0(1)**2
      csPsd=3*SCcoeff*pA5(1)**2
      csNsd=3*SCcoeff*nA5(1)**2
      IF( pA0(1)*nA0(1) .lt. 0) csNsi=-csNsi
      IF( pA5(1)*nA5(1) .lt. 0) csNsd=-csNsd
      PROB(53)=DDIM(csPsi/XENON100(DABS(MNEU(1))),1d0)

      IF (OMGFLAG.EQ.2) RETURN

*  Computing indirect detection rate

 1      sigmaV=calcSpectrum(1.D-3,0,tab,err)
      IF (err.EQ.0) sigmaV=0d0
      IF (sigmaV.NE.0d0) THEN
       DO I=1,NBIN
        dNdx(I)=zInterp(DLOG(10d0)*x(I),tab)*DLOG(10d0)
       ENDDO
      ENDIF

      END


      DOUBLE PRECISION FUNCTION XENON100(M)

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=117)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/6.6499d0,6.7050d0,6.7387d0,6.7789d0,6.8112d0,6.8530d0,
     .6.8844d0,6.9278d0,6.9589d0,7.0035d0,7.0371d0,7.1141d0,7.2082d0,
     .7.3145d0,7.4227d0,7.5285d0,7.6749d0,7.7572d0,7.8661d0,8.0189d0,
     .8.1547d0,8.2812d0,8.4380d0,8.5784d0,8.7594d0,8.8924d0,9.0539d0,
     .9.3423d0,9.5459d0,9.7535d0,9.9660d0,1.0209d1,1.0546d1,1.0895d1,
     .1.1129d1,1.1464d1,1.1887d1,1.2255d1,1.2684d1,1.3075d1,1.3531d1,
     .1.3979d1,1.4580d1,1.5065d1,1.5729d1,1.6423d1,1.7166d1,1.7894d1,
     .1.8733d1,1.9771d1,2.0874d1,2.1912d1,2.3103d1,2.4274d1,2.5091d1,
     .2.6356d1,2.8131d1,3.0026d1,3.2048d1,3.4021d1,3.5957d1,3.8345d1,
     .4.0479d1,4.2971d1,4.5370d1,4.7902d1,5.1128d1,5.3982d1,5.7305d1,
     .6.1245d1,6.4931d1,6.8555d1,7.2381d1,7.6900d1,8.1568d1,8.6121d1,
     .9.0928d1,9.6548d1,1.0359d2,1.1005d2,1.1785d2,1.2568d2,1.3299d2,
     .1.4195d2,1.4987d2,1.5824d2,1.6707d2,1.7639d2,1.8624d2,1.9663d2,
     .2.0761d2,2.1920d2,2.3213d2,2.4972d2,2.6653d2,2.8245d2,2.9875d2,
     .3.2059d2,3.3849d2,3.5738d2,3.7733d2,3.9839d2,4.2063d2,4.4411d2,
     .4.7402d2,5.0047d2,5.2841d2,5.6400d2,5.9548d2,6.2872d2,6.6402d2,
     .7.0698d2,7.5623d2,7.9844d2,8.5222d2,9.0961d2,9.6514d2/
      DATA Y/9.3189d-4,7.5068d-4,6.0470d-4,4.8711d-4,3.9238d-4,
     .3.1608d-4,2.5461d-4,2.0510d-4,1.6522d-4,1.3309d-4,1.0721d-4,
     .8.6360d-5,6.9567d-5,5.4755d-5,4.1401d-5,3.1938d-5,2.4638d-5,
     .1.9847d-5,1.5311d-5,1.1812d-5,9.4389d-6,7.6644d-6,6.1740d-6,
     .4.9734d-6,3.8367d-6,3.0907d-6,2.3842d-6,1.8401d-6,1.4822d-6,
     .1.1935d-6,9.6141d-7,7.7948d-7,6.0169d-7,4.6444d-7,3.7127d-7,
     .3.0797d-7,2.4091d-7,1.9406d-7,1.5743d-7,1.2868d-7,1.0738d-7,
     .8.9095d-8,7.0234d-8,5.9077d-8,4.9615d-8,4.1804d-8,3.5189d-8,
     .2.9183d-8,2.5898d-8,2.2742d-8,1.9141d-8,1.6754d-8,1.5223d-8,
     .1.3064d-8,1.0939d-8,9.9572d-9,8.7390d-9,8.0769d-9,8.0459d-9,
     .7.9072d-9,7.7230d-9,7.2490d-9,6.9450d-9,7.1034d-9,7.1837d-9,
     .7.0974d-9,7.0963d-9,7.2382d-9,7.2520d-9,7.4232d-9,7.7400d-9,
     .7.9072d-9,8.0962d-9,8.4464d-9,8.7796d-9,9.0026d-9,9.1569d-9,
     .9.5793d-9,1.0470d-8,1.0949d-8,1.1404d-8,1.1913d-8,1.2403d-8,
     .1.3031d-8,1.3595d-8,1.4108d-8,1.4773d-8,1.5396d-8,1.6147d-8,
     .1.6829d-8,1.7567d-8,1.8354d-8,1.9611d-8,2.0916d-8,2.1997d-8,
     .2.4353d-8,2.5961d-8,2.7401d-8,2.8492d-8,3.0121d-8,3.1568d-8,
     .3.2963d-8,3.4508d-8,3.6283d-8,3.8186d-8,3.9928d-8,4.1873d-8,
     .4.4618d-8,4.6937d-8,4.9574d-8,5.3023d-8,5.6418d-8,5.9839d-8,
     .6.3063d-8,6.7203d-8,7.1799d-8,7.4469d-8/

      IF(M.LE.X(1))THEN
       XENON100=Y(1)
       RETURN
      ENDIF
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        XENON100=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO
      XENON100=Y(N)

      END
