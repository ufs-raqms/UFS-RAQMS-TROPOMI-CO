C----jv_cmn.h---COMMON BLOCKS for new FAST-J code (wild/prather 7/99)
c
c  Parameters
c  ----------
c
c     NB    Number of levels in CTM plus one for above model top
c     NC    Number of levels in the fundamental Fast-J grid
c           (now defined in individual subroutines, hyl, 04/14/03)
c     NS    Maximum number of species which require J-values calculating
c     NW    Maximum number of wavelength bins that can be used
c     NP    Maximum number of aerosol/cloud types that can be used
c     NH    Maximum number of Herzberg X-sections that can be used
c     MX    Number of aerosol/cloud types supplied from CTM
CCCMNN0308
C     NDUST Number of dust size bins, treated as separate aerosols
C     NAER  Number of aerosol types (5) X Number of relative humidity
C           bins (5)
c
c                                        Note: THETA(NL) no longer used
c  =====================================================================
c hyl
      INTEGER        NH
      PARAMETER   (           NH=7)
      CHARACTER*20 TITLEA(NP)
      CHARACTER*78 TITLE0
      CHARACTER*7  TITLEJ(3,NS), jlabel(jppj), hzlab(nh)
      INTEGER jind(jppj),nhz,hzind(nh)
      INTEGER NJVAL,NW1,NW2,MIEDX,NAA,NLBATM,npdep,jpdep(NS)
c hyl
      REAL*8 ZZHT,RAD
      REAL*8 WBIN,WL,FL,QO2,QO3,Q1D,QQQ,QRAYL,TQQ,WAA,QAA,PAA
      REAL*8 RAA,SSA,TREF,OREF,BREF,QBC,DBC,zpdep(NW,3)
      REAL*8 dtaumax,szamax,jfacta(jppj)
      REAL*8 dtausub,dsubdiv
C----------------------------------------------------------------------
cJAA  For parallelization, separated common blocks into those which
cJAA   need to be threadprivate and those which don't.
cJAA Global blocks:
      COMMON /TITLS/TITLE0,TITLEJ,TITLEA
      COMMON /ATMOS/RAD,ZZHT
      COMMON /CCWVL/WBIN(NW+1),WL(NW),FL(NW),QO2(NW,3),QO3(NW,3),
     .              Q1D(NW,3),QQQ(NW,2,NS-3),QRAYL(NW+1),TQQ(3,NS),
     .              WAA(4,NP),QAA(4,NP),
     .              PAA(8,4,NP),RAA(4,NP),SSA(4,NP),QBC(NW),
     .              NJVAL,NW1,NW2,MIEDX(MX),NAA,NLBATM
      COMMON /CLIM/ TREF(51,18,12),OREF(51,18,12),BREF(51)
      COMMON /JCNTR/dtaumax,szamax
      COMMON /JVALS/jfacta,zpdep,npdep,jpdep,jind,jlabel
      COMMON /JVSUB/dtausub,dsubdiv
C-----------------------------------------------------------------------
