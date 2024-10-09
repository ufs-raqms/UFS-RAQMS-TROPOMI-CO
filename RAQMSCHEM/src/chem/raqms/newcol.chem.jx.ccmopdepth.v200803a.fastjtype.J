#include <options.h>
#ifdef JVALUEJ
      module iatjat
      logical doiat
      end module iatjat
#endif
!    #define DIAGCHEM
c---(pphot.f)-------generic CTM shell from UCIrvine (p-code 4.0, 7/99)
c---------Oliver Wild (7/99,troposphere), Huisheng Bian (8/02,stratosphere)
c---------PPHOT calculates photolysis rates with the Fast-J2 scheme
c---------subroutines:  inphot, photoj, Fast-J2 schemes...
c-----------------------------------------------------------------------
cMNN 083007   This version includes modification adopted from FAST_JX V 5.7
cMNN          (Algorithm for MIE code levels considers all the aerosols and
cMNN           cloud).  Both BC and OC are included and are assumed to have
cMNN           same absorption coefficients. The optical depth for BC passed
cMNN           from the driver routine includes the OD due to OC.
c
c Revised for interfacing with RAQMSregional (hyl, 11/22/02, 03/31/03)
c
      subroutine inphot(nlayer,nreacs)
      use raqmschem_pmgrid_mod, only : prefi
C     subroutine inphot
C-----------------------------------------------------------------------
C  Routine to initialise photolysis rate data, called directly from the
C  cinit routine in ASAD. Currently use it to read the JPL spectral data
C  and standard O3 and T profiles and to set the appropriate reaction index.
C-----------------------------------------------------------------------
C  Add the following input variables (hyl, 11/13/02)
C
C  Variable  Type    Dimensn Units   Description
C  --------  ----    ------- -----   -----------
C  nlayer    int        -      -     Top layer where J-values required
C  nreacs    int        -      -     Total # of photolysis reactions
C-----------------------------------------------------------------------
c
c     iph       Channel number for reading all data files
c     rad       Radius of Earth (cm)
c     zzht      Effective scale height above top of atmosphere (cm)
c     dtaumax   Maximum opt.depth above which a new level should be inserted
c     dtausub   No. of opt.depths at top of cloud requiring subdivision
c     dsubdiv   Number of additional levels to add at top of cloud
c     szamax    Solar zenith angle cut-off, above which to skip calculation
c
C-----------------------------------------------------------------------
c
      use raqmschem_pmgrid_mod, only : nlev
      use raqmschemcomm_mod, only : np,ns,nw,mx,lpar,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

C=============== INPUT PARAMETERS ======================================
      integer nlayer, nreacs

C=============== LOCAL VARIABLES =======================================
      integer iph
!      integer lpar
!      lpar=nlev

      IF ( nlayer .ne. lpar ) THEN
         PRINT*, 'nlayer = ', nlayer, ' ne lpar! ',lpar
         PRINT*, 'STOP in subroutine inphot.f'
         call killit('layer')
!         STOP
      ENDIF

      IF ( nreacs .ne. jppj ) THEN
         PRINT*, 'nreacs = ', nreacs, ' ne JPPJ! ',jppj
         PRINT*, 'STOP in subroutine inphot.f'
         call killit('nread')
!         STOP
      ENDIF

c Use channel 8 to read files at the moment
      iph=8

c Defaults & constants
      RAD = 6375.d5
      ZZHT = 5.d5

c Calculate new levels if tau exceeds DTAUMAX
      dtaumax = 1.d0
c Add dsubdiv additional levels in first dtausub of cloud
      dtausub = 1.0d0
      dsubdiv = 10.d0

c Maximum Zenith Angle  (98 degrees at 63 km; 99 degrees at 80 km)
      szamax = 98.0d0

c You may need to change the directory paths for data files 'ratj.d', 
c 'jv_spec.dat' and 'jv_atms.dat'. Also change the string lengths in 
c the subroutines using them. (hyl, 11/22/02)

c Read in labels of photolysis rates required
!      CALL RD_JS(iph, './Input/Fastjx/ratj.d')
      CALL RD_JS(iph, trim(prefi)//'/Fastjx/ratj.d')

c Read in JPL spectral data set
!      CALL RD_TJPL(iph, './Input/Fastjx/jx_spec.dat')
      CALL RD_TJPL(iph,trim(prefi)//'/Fastjx/jx_spec.dat')

c Read in T & O3 climatology
!      CALL RD_PROF(iph, './Input/Fastjx/jv_atms.dat')
      CALL RD_PROF(iph,trim(prefi)//'/Fastjx/jv_atms.dat')

c Select Aerosol/Cloud types to be used
      CALL SET_AER

      return
      end

      subroutine RD_JS(nj1,namfil)
C-----------------------------------------------------------------------
c  Reread the ratj.d file to map photolysis rate to reaction
c  Read in quantum yield 'jfacta' and fastj2 label 'jlabel'
C-----------------------------------------------------------------------
c
c     jfacta    Quantum yield (or multiplication factor) for photolysis
c     jlabel    Reference label identifying appropriate J-value to use
c     ipr       Photolysis reaction counter - should total 'jppj'
c
C-----------------------------------------------------------------------
      use raqmschem_pmgrid_mod, only : masterproc,iam
#ifdef DOMPI
      use mpishorthand, only  : mpireal,mpicom,mpiint,mpichar
#endif
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
      use raqmschemcomm_mod, only : np,ns,nw,mx,jppj
      implicit none
 
!      include 'cmn_fj.h'
#include "jv_cmn.h"
c
      integer nj1, ipr, i
c     character*6  namfil    !original
!      character*21 namfil    !hyl
      character*(*) namfil    !hyl
      character*120 cline
      integer localrc
c
c Reread the ratj.d file to map photolysis rate to reaction
c Read in quantum yield jfacta and fastj2 label jlabel
      if(masterproc)then
        ipr=0
        open(nj1,file=namfil,status='old',form='formatted')
 10     read(nj1,'(a)',err=20) cline
        if(cline(2:5).eq.'9999') then
          go to 20
        elseif(cline(1:1).eq.'#') then
          go to 10
        elseif(cline(5:5).eq.'$') then
          go to 10
        else
          ipr=ipr+1
          read(cline(79:83),'(f5.1)') jfacta(ipr)
          read(cline(86:92),'(a7)')   jlabel(ipr)
          jfacta(ipr)=jfacta(ipr)/100.d0
c -----------------------------------------------------------------------
c Later, if necessary, add code to read reaction names and branch numbers
c (hyl, 11/12/02)
c         read (cline(7:10),"(a4)") rnames(ipr)
c         rnames(ipr) = trim(rnames(ipr))
c         branch(ipr) = 1
c         do i=1,ipr-1
c            if (rnames(ipr) == rnames(i)) branch(ipr) = branch(i) + 1
c         enddo
c ----------------------------------------------------------------------

          go to 10
        endif
 20     close(nj1)
        if(ipr.ne.jppj) then
          write(6,1000) ipr,jppj
          call flush(6)
          write(0,*)'error aa'
          call flush(0)
          call killit('error aa')
!          stop
        endif
c
c       Print details to standard output
        write(6,1100) ipr
        write(6,1200) (i, jlabel(i), jfacta(i),i=1,ipr)
      endif
#ifdef DOMPI
!      call mpibcast(ipr,1,mpiint,0,mpicom)
      call raqmschem_comm_all_bcast(ipr,rc=localrc)
!      call mpibcast(jlabel,7*jppj,mpichar,0,mpicom)
      call raqmschem_comm_all_bcast(jlabel,rc=localrc)
!      call mpibcast(jfacta,jppj,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(jfacta,rc=localrc)
#ifdef DIAGCHEM
      if(iam.eq.1)then
         write(6,*)'ipr',ipr
         write(6,*)'jlabel',jlabel
         write(6,*)'jfacta',jfacta
      endif
#endif
#endif
c
      return
 1000 format(' Error: ',i3,' photolysis labels but ',i3,' reactions')
 1100 format(' Fast-J2 Photolysis Scheme: considering ',i2,' reactions')
 1200 format(3x,10(3(i2,': ',a7,' (Q.Y. ',f6.3,') '),/,3x))
      end


      subroutine photoj(nslon, nslat, lbot, rlat, gmt, julday,imon, sza,
     &                  sp, zsurf, pgrd, zgrd, t, dbox, dbcolo3v, 
cjaa     &                  od, sa, zpj )
cjaa 12/7/2004 adapt for use with ccm3 radiative code
c    ccm3 provides optical depth for liquid and ice separately so pass 
c    in to photoj as odl (liquid) and odi (ice).
     &                  odl, odi, aer_opd,dust_opd, sa, zpj )

C-----------------------------------------------------------------------
C----jv_trop.f:  new FAST J-Value code, troposphere only (mjprather 6/96)
C----     uses special wavelength quadrature spectral data (jv_spec.dat)
C---      that includes only 289 nm - 800 nm  (later a single 205 nm add-on)
C---      uses special compact Mie code based on Feautrier/Auer/Prather vers.
C-----------------------------------------------------------------------
c
c     zpj      External array providing J-values to main CTM code
c     timej    Offset in hours from start of timestep to time J-values
c              required for - take as half timestep for mid-step Js.
c              (timej is now obsolete, hyl)
c     solf     Solar distance factor, for scaling; normally given by:
c                      1.0-(0.034*cos(real(iday-186)*2.0*pi/365.))
c
CCCMNN0308  
CCC   aer_opd  optical depth due to aerosols (2D array)
CCC   dust_opd optical depth due to dust (5 size bins)
C----------basic common blocks:-----------------------------------------
#ifdef JVALUEJ
!      use ijprn
      use raqmschem_pmgrid_mod, only : iprn,jprn,kprn
#endif
      use raqmschem_pmgrid_mod, only : ibeg,tile,nlev
      use raqmschemcomm_mod, only : ndust,naer
      use raqmschemcomm_mod, only : np,ns,nw,mx,fastj,jppj,jpnl,lpar
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

chyl  --- Input parameters ---
      integer imon, nslat, nslon, julday, lbot
      real*8  gmt, rlat, sp, zsurf
cjaa      real*8  pgrd(lbot),zgrd(lbot),t(lbot), od(lbot)
      real*8  pgrd(lbot),zgrd(lbot),t(lbot), odl(lbot), odi(lbot)
CC      real*8  odbc(lbot)
      real*8  aer_opd(NAER,lbot),dust_opd(NDUST,lbot)

      !dbox is o3 mixing ratio, dbcolo3v is vert o3 column above 
      !the top of each grid box (see subroutine calc_colo3_col)
      real*8  dbox(lbot), dbcolo3v(lbot)
      real*8  sza, sa
      real*8  timej, solf

chyl  --- local variables ---
      integer i,j

C save J-values in a 2D array, instead of a 4D array which requires 
C too large a memory (hyl, 02/20/03) 
      real*8 zpj(jpnl,jppj)
!      integer lpar
!      lpar=nlev

C-----------------------------------------------------------------------
C---Initialize ZJ and ZPJ
#ifdef JRATETJ
      if(nslon.eq.iprn.and.nslat.eq.jprn+1)then
        write(6,*)'top of photoj jpnl',jpnl,'jppj',jppj,'lbot',lbot,
     * 'sza',sza,'szamax',szamax
      write(6,*)'gmt',  gmt, rlat, sp, zsurf
      endif
#endif
      do i=1,jpnl
        do j=1,jppj
!          zj(i,j)=0.D0
          fastj%zj(i,j)=0.D0
          zpj(i,j)=0.D0
        enddo
      enddo

C---Calculate new solar zenith angle
C     CALL SOLAR2(timej)
C now import sza (deg) from CTM (hyl,02/20/03)
!      U0 = cos(sza*3.141592653589793D0/180.D0)
      fastj%U0 = cos(sza*3.141592653589793D0/180.D0)

C-----------------------------------------------------------------------
C---If you want to set sza = 0 degrees for testing,
C---then uncomment the following lines (hyl, 11/13/02) 
C     sza = 0.0d0
C-----------------------------------------------------------------------
!      if(nslon.eq.1.and.nslat.eq.1.and.tile.eq.1)then
!         write(6,*)'sza',sza,'szamax',szamax
!         call flush(6)
!      endif
      if(sza.gt.szamax) return         !hyl
!       aer_opd=0.0
!       dust_opd=0.0
!       odl=0.0
!       odi=0.0

C---Set up profiles on model levels
       call set_prof (lbot, sp, zsurf, pgrd, zgrd, t, dbox, dbcolo3v, 
cjaa     &                od, sa, imon, rlat )
CMNN 0308     &                odl, odi, odbc, sa, imon, rlat )
     &    odl, odi,aer_opd,dust_opd,sa,imon,rlat)

C---Print out atmosphere (here sza is in degree, hyl,11/13/02)
#ifdef DOPRTATM
       call PRTATM (3, lbot, gmt, imon, rlat, sza)
#endif

#if (defined JVALUEJ || defined DIAGOPMIE)
       call JVALUE (lbot, sza, sa,nslon,nslat)
#else
C-----------------------------------------------------------------------
       call JVALUE (lbot, sza, sa)
#endif
C-----------------------------------------------------------------------

C---Print solar flux terms
c      WRITE(6,'(A16,I5,20I9)') '  wave (beg/end)',(i,i=1,jpnl)  
c      DO j=NW1,NW2
c        WRITE(6,'(2F8.2,20F9.6)') WBIN(j),WBIN(j+1),
c     $                            (FFF(j,i)/FL(j),i=1,jpnl)
c      ENDDO
c
C---Include variation in distance to sun
cc     solf=1.d0-(0.034d0*cos(dble(iday-186)*2.d0*pi/365.d0))
c      solf=1.d0
C-----------------------------------------------------------------------
c now calculate solf within JRATET (hyl, 11/12/02)
#ifdef JRATETJ
       CALL JRATET(lbot, t,julday,nslon,nslat)
#else
       CALL JRATET(lbot, t,julday)
#endif
C-----------------------------------------------------------------------
C
C  "zj" updated in JRATET - pass this back to chem_driver as "zpj" (hyl)
      do i=1,lbot
        do j=1,jppj
!           zpj(i+lpar-lbot,j)= zj(i,j)
           zpj(i+lpar-lbot,j)= fastj%zj(i,j)
!           if(isnan(zpj(i+lpar-lbot,j)))then
!            write(6,*)'photoj zpj',i,lpar,lbot,j,zpj(i+lpar-lbot,j)
!            call flush(6)
!            stop
!          endif
#ifdef DOCHEMJJJ
           if(nslon.eq.iprn.and.nslat.eq.jprn+1)then
             write(6,*)'return zpj',i+lpar-lbot,j,fastj%zj(i,j),'i',i,j
           endif
#endif
        enddo
      enddo

      return
      end


       subroutine set_prof(lbot,sp,zsurf,pgrd,zgrd,t,dbox,dbcolo3v,
cjaa     &                     odcol,sa,imon, rlat)
c jaa 12/7/2004 adapt for use with ccm3 radiative code
c    ccm3 provides optical depth for liquid and ice separately so pass 
c    in to fastj as odcoll (liquid) and odcoli (ice).
CCCMNN 0308     &                     odcoll,odcoli,odcolbc,sa,imon, rlat)
     &  odcoll,odcoli,aer_opd,dust_opd,sa,imon,rlat)
c
C-----------------------------------------------------------------------
c  Routine to set up atmospheric profiles required by Fast-J2 using a
c  doubled version of the level scheme used in the CTM. First pressure
c  and z* altitude are defined, then O3 and T are taken from the supplied
c  climatology and integrated to the CTM levels (may be overwritten with
c  values directly from the CTM, if desired) and then black carbon and
c  aerosol profiles are constructed.
c                                                     Oliver (04/07/99)
C---------------------------------------------------------------------------------
C  Add the following input variables for CTM interface (hyl, 11/13/02)
C
C  Variable  Type    Dimension   Units   Description
C  --------  ----    ---------   -----   -----------
C  SP        dble     -        [mb]     Surface Pressure
C  T         dble   [lbot]      [K]     Vertical temperature profile
C  DBOX      dble   [lbot]       -      Ozone mixing ratio profile
C  DBCOLO3V  dble   [lbot] [molec cm-2] O3 column above top of each model level 
C  ODCOL     dble   [lbot]       -      Vertical optical depth profile
C  SA        dble     -          -      Surface Albedo
C--------------------------------------------------------------------------------
c
c     pj       Pressure at boundaries of model levels (hPa)
c     z        Altitude of boundaries of model levels (cm)
c     odcol    Optical depth at each model level
c     masfac   Conversion factor for pressure to column density
c
c     TJ       Temperature profile on model grid
c     DM       Air column for each model level (molecules.cm-2)
c     DO3      Ozone column for each model level (molecules.cm-2)
c     DBC      Mass of Black Carbon at each model level (g.cm-3)  !  .....!
c     PSTD     Approximate pressures of levels for supplied climatology
c
C-----------------------------------------------------------------------
      use raqmschemcomm_mod, only : ndust,naer
      use raqmschemcomm_mod, only : np,ns,nw,mx
      use raqmschemcomm_mod, only : fastj,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"
cMNN 083007
#include "jv_mie.h"

Chyl
C=============== INPUT PARAMETERS ======================================
      integer imon, lbot 
      real*8  sp, t(lbot), sa, rlat 
cjaa      real*8  odcol(lbot)
CCCMNN 0308      real*8  odcoll(lbot),odcoli(lbot),odcolbc(lbot)
      real*8  odcoll(lbot),odcoli(lbot),aer_opd(NAER,lbot)
      real*8  dust_opd(NDUST,lbot)
      !dbox is o3 mixing ratio, dbcolo3v is vert o3 column above
      !the top of each grid box (see subroutine calc_colo3_col)
      real*8  dbox(lbot), dbcolo3v(lbot), pgrd(lbot), zgrd(lbot), zsurf
#ifdef DOINT
      real r1000,zero8
#endif

C=============== LOCAL VARIABLES =======================================
      integer i, k, l, mn, ndum

      real*8  dlogp,f0,t0,b0,pb,pc,xc,masfac,scaleh
      real*8  pstd(52),oref2(51),tref2(51),bref2(51)
      real*8  ATAU
      integer JTAUMX

c  First calculate pressure (pj) at boundaries of CTM levels (lbot+1)

c  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      masfac=100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

c sp is air pressure at Earth's surface and topography has been accounted for (as of 12/16/02).
!      pj(1) = sp
      fastj%pj(1) = sp

c valid for constant grid spacing (hyl,04/01/03)
!      z(1) = zsurf*100.d0           !m-->cm
      fastj%z(1) = zsurf*100.d0           !m-->cm
      scaleh=1.3806d-19*masfac*t(1)
!      pj(2) = pgrd(1)*dexp(-(zgrd(2)-zgrd(1))*100.d0/2./scaleh) 
      fastj%pj(2) = pgrd(1)*dexp(-(zgrd(2)-zgrd(1))*100.d0/2./scaleh) 
!      z(2) = zgrd(1)*100.d0 + (zgrd(2)-zgrd(1)) *100.d0 / 2.
      fastj%z(2) = zgrd(1)*100.d0 + (zgrd(2)-zgrd(1)) *100.d0 / 2.
      do i=2,lbot 
        scaleh=1.3806d-19*masfac*t(i)
!        pj(i+1) = pj(i)*dexp(-(zgrd(i)-zgrd(i-1))*100.d0/scaleh) 
        fastj%pj(i+1) = fastj%pj(i)*dexp(-(zgrd(i)-zgrd(i-1))*100.d0/scaleh) 
!        z(i+1) = z(i) + (zgrd(i)-zgrd(i-1))*100.d0     !deltaZ=400*100cm
        fastj%z(i+1) = fastj%z(i) + (zgrd(i)-zgrd(i-1))*100.d0     !deltaZ=400*100cm
      enddo
!      pj(lbot+1+1) = 0.d0     !pj(NB+1) = 0.d0
      fastj%pj(lbot+1+1) = 0.d0     !pj(NB+1) = 0.d0

c  Set up cloud and surface properties
cjaa       call cldsrf(lbot,odcol,sa)
CCCMNN 0308       call cldsrf(lbot,odcoll,odcoli,odcolbc,sa)
!      if(isnan(aer_opd(1,1)))then
         !write(6,*)'call cldsrf',aer_opd(1:4,1)
!      endif
      call cldsrf(lbot,odcoll,odcoli,aer_opd,dust_opd,sa)
c
c  Set up pressure levels for O3/T climatology - assume that value
c  given for each 2 km z* level applies from 1 km below to 1 km above,
c  so select pressures at these boundaries. Surface level values at
c  1000 mb are assumed to extend down to the actual P(nslon,nslat).
c
#ifdef DOINT
      r1000=1000.
!      pstd(1) = max(pj(1),r1000)
      pstd(1) = max(fastj%pj(1),r1000)
#else
!      pstd(1) = max(pj(1),1000.d0)
      pstd(1) = max(fastj%pj(1),1000.d0)
#endif
      pstd(2) = 1000.d0*10.d0**(-1.d0/16.d0)
      dlogp = 10.d0**(-2.d0/16.d0)
      do i=3,51
        pstd(i) = pstd(i-1)*dlogp
      enddo
      pstd(52) = 0.d0
c
c  Select appropriate monthly and latitudinal profiles
c  Now use RLAT instead of Oliver's YDGRD(NSLAT) (hyl, 11/12/00)
       mn = max(1,min(12,imon))
       l = max(1,min(18,(int(rlat)+99)/10))
c
c  Temporary arrays for climatology data
      do i=1,51
        oref2(i)=oref(i,l,mn)
        tref2(i)=tref(i,l,mn)
        bref2(i)=bref(i)
      enddo
c
c  Apportion O3 and T on supplied climatology z* levels onto CTM levels 
c  with mass (pressure) weighting, assuming constant mixing ratio and
c  temperature half a layer on either side of the point supplied.
c
      do i = 1,lbot+1        !do i = 1,NB
        F0 = 0.d0
        T0 = 0.d0
        B0 = 0.d0
        do k = 1,51
!          PC = min(pj(i),pstd(k))
          PC = min(fastj%pj(i),pstd(k))
!          PB = max(pj(i+1),pstd(k+1))
          PB = max(fastj%pj(i+1),pstd(k+1))
          if(PC.gt.PB) then
!            XC = (PC-PB)/(pj(i)-pj(i+1))
            XC = (PC-PB)/(fastj%pj(i)-fastj%pj(i+1))
            F0 = F0 + oref2(k)*XC
            T0 = T0 + tref2(k)*XC
            B0 = B0 + bref2(k)*XC
          endif
        enddo
!        TJ(i) = T0
!        DO3(i)= F0*1.d-6
!        DBC(i)= B0
        fastj%TJ(i) = T0
        fastj%DO3(i)= F0*1.d-6
        fastj%DBC(i)= B0
      enddo
c --------------------------------------------------------------------------
c from original Fast-J2 code distribution:
c Insert model values here to replace or supplement climatology.
c Note that CTM temperature is always used in x-section calculations
c (see JRATET); TJ is used in actinic flux calculation only.
c     do i=1,lpar
c       DO3(i) = my_ozone(i)        ! Volume Mixing Ratio
c       TJ(i)  = T(nslon,nslat,I)   ! Kelvin
c     enddo
c     DO3(lpar+1) = my_ozone*exp()  ! Above top of model (or use climatology)
c     TJ(lpar+1)  = my_temp(lpar)   ! Above top of model (or use climatology)
c ---------------------------------------------------------------------------
c
c  Add Aerosol Column - include aerosol types here. Currently use soot
c  water and ice; assume black carbon x-section of 10 m2/g, independent
c  of wavelength; assume limiting temperature for ice of -40 deg C.
c
      do i=1,lbot
cjaa        AER(1,i) = DBC(i)*10.d0*(z(i+1)-z(i))
cjaa        if(T(I).gt.233.d0) then    !hyl
cjaa          AER(2,i) = odcol(i)
cjaa          AER(3,i) = 0.d0
cjaa        else
cjaa          AER(2,i) = 0.d0
cjaa          AER(3,i) = odcol(i)
cjaa        endif          
CCCMNN 0308        AER(1,i) =  odcolbc(i)
!        AER(2,i) = odcoll(i)
!        AER(3,i) = odcoli(i)
        fastj%AER(2,i) = odcoll(i)
        fastj%AER(3,i) = odcoli(i)
c Later, we can also add in aerosol optical depth columns (hyl, 11/12/00)
!        if(isnan(odcoll(i)).or.isnan(odcoli(i)))then
!          write(6,*)i,'odcoll',odcoll(i),odcoli(i)
!         endif
CCCMNN 0308
       do ndum=1,NDUST
!         AER(3+ndum,i) = dust_opd(ndum,i)     
         fastj%AER(3+ndum,i) = dust_opd(ndum,i)     
!          if(isnan(dust_opd(ndum,i)))then
            !write(6,*)'dustopd',ndum,i
!           endif
         enddo
CCCMNN 0308    add aerosol optical depth columns
       do ndum=1,NAER
!        AER(3+ndust+ndum,i)=aer_opd(ndum,i)
        fastj%AER(3+ndust+ndum,i)=aer_opd(ndum,i)
!          if(isnan(aer_opd(ndum,i)))then
!            !write(6,*)'aeropd',ndum,i,aer_opd(ndum,i)
!          endif
        end do
      enddo

      do k=1,MX
!        AER(k,lbot+1) = 0.d0      
        fastj%AER(k,lbot+1) = 0.d0      
      enddo
c
CMNN 083007  Calculate extra levels for scattering code (Adopted from V 5.7)
CCC
      JTAUMX = 250
      ATAU = 1.12D0
      ATAU0 = 0.01D0
      CALL EXTRAL(JTAUMX,ATAU,ATAU0)
CCC

c  Calculate column quantities for Fast-J2
      do i=1,lbot+1      !do i=1,NB
!        DM(i)  = (PJ(i)-PJ(i+1))*masfac
        fastj%DM(i)  = (fastj%PJ(i)-fastj%PJ(i+1))*masfac
!        DO3(i) = DO3(i)*DM(i) 
        fastj%DO3(i) = fastj%DO3(i)*fastj%DM(i) 
cc      DO3(i) = DO3(i)*DM(i) / (2.69*1e16)      !Dobson Unit 
      enddo

c Here replace O3 climatology with RAQMS model level column o3 (density)
c (hyl, 01/23/03)
!      DO3(lbot+1) = dbcolo3v(lbot)      
      fastj%DO3(lbot+1) = dbcolo3v(lbot)      
      do i=lbot,2,-1                   
!         DO3(i) = dbcolo3v(i-1) - dbcolo3v(i)
         fastj%DO3(i) = dbcolo3v(i-1) - dbcolo3v(i)
      enddo
#ifdef DOINT
      zero8=0.0
!      DO3(1) = max(dbox(1),zero8) * DM(1) 
      fastj%DO3(1) = max(dbox(1),zero8) * fastj%DM(1) 
#else
!      DO3(1) = max(dbox(1),0.D0) * DM(1) 
      fastj%DO3(1) = max(dbox(1),0.D0) * fastj%DM(1) 
#endif
c     print *, dbcolo3v(2),dbcolo3v(1), DO3(2), DO3(1)

c Replace T climatology with RAQMS T (kelvin), and use climatology above 
c top of model
      do i=1,lbot
!        TJ(i)  = t(i)              
        fastj%TJ(i)  = t(i)              
      enddo

      return
      end

cjaa       subroutine cldsrf(lbot,odcol,sa)
cjaa 12/7/2004 adapt for use with ccm3 radiative code
c    ccm3 provides optical depth for liquid and ice separately so pass 
c    in to fastj as odcoll (liquid) and odcoli (ice).
CCCMNN 0308       subroutine cldsrf(lbot,odcoll,odcoli,odcolbc,sa)
      subroutine cldsrf(lbot,odcoll,odcoli,aer_opd,dust_opd,sa)
c
C-----------------------------------------------------------------------
c  Routine to set cloud and surface properties
C-----------------------------------------------------------------------
c     rflect   Surface albedo (Lambertian)
c     odmax    Maximum allowed optical depth, above which they are scaled
c     odcol    Optical depth at each model level
c     odsum    Column optical depth
c     nlbatm   Level of lower photolysis boundary - usually surface
C-----------------------------------------------------------------------
      use raqmschemcomm_mod, only : ndust,naer
      use raqmschemcomm_mod, only : np,ns,nw,mx,fastj,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

Chyl
C=============== INPUT PARAMETERS ======================================
      integer lbot
cjaa      real*8  odcol(lbot), sa       
CCCMNN 0308      real*8  odcoll(lbot), odcoli(lbot), odcolbc(lbot), sa       
      real*8  odcoll(lbot), odcoli(lbot), aer_opd(NAER,lbot)
      real*8  dust_opd(NDUST,lbot), sa
C=============== LOCAL VARIABLES =======================================
      integer i, j, k, ia
      real*8  odsum, odmax, odtot
#ifdef DOINT
      real zero,one
#endif
c
c Default lower photolysis boundary as bottom of level 1
      nlbatm = 1
c
c Set surface albedo
!      RFLECT = dble(sa)     !hyl
      fastj%RFLECT = dble(sa)     !hyl
#ifdef DOINT
      zero=0.0
      one=1.0
!      RFLECT = max(zero,min(one,RFLECT))
      fastj%RFLECT = max(zero,min(one,fastj%RFLECT))
#else
!      RFLECT = max(0.d0,min(1.d0,RFLECT))
      fastj%RFLECT = max(0.d0,min(1.d0,fastj%RFLECT))
#endif
c
c Zero aerosol column
      do k=1,MX
        do i=1,lbot+1            
!          AER(k,i) = 0.d0
          fastj%AER(k,i) = 0.d0
        enddo
      enddo
c
c Scale optical depths as appropriate - limit column to 'odmax'
      odmax = 200.d0
      odsum =   0.d0
      do i=1,lbot    
cjaa        odcol(i) = dble(odcol(i))       !hyl
cjaa        odsum    = odsum + odcol(i)
CCCMNN 0308        odsum    = odsum + odcoll(i)+odcoli(i)+odcolbc(i)
      odsum    = odsum + odcoll(i)+odcoli(i)
       do ia=1,NAER
        odsum = odsum + aer_opd(ia,i)
        end do
       do ia=1,NDUST
        odsum = odsum + dust_opd(ia,i)
        end do
       enddo
      if(odsum.gt.odmax) then
        odsum = odmax/odsum
        do i=1,lbot       
cjaa          odcol(i) = odcol(i)*odsum
          odcoll(i) = odcoll(i)*odsum
          odcoli(i) = odcoli(i)*odsum
CCCMNN 0308 odcolbc(i) = odcolbc(i)*odsum
          do ia = 1,NAER
           aer_opd(ia,i)= aer_opd(ia,i)*odsum
           end do
          do ia = 1,NDUST
           dust_opd(ia,i) = dust_opd(ia,i)*odsum
           end do
        enddo
        odsum = odmax
      endif
c
c  Use clear-sky conditions
c     do i=1,jpnl
c       odcol(i)=0.d0
c     enddo
c
c Set sub-division switch if appropriate
cMNN 083007   Number of extra levels for scattering code is
cMNN          evaluated in a separate subroutine called from SET_PROF
c     odtot=0.d0
c     jadsub(lbot+1)=0       
c     jadsub(lbot)=0       
c     do i=lbot,1,-1   
c       k=2*i
c       jadsub(k)=0
c       jadsub(k-1)=0
cjaac       odtot=odtot+odcol(i)
cjaac       if(odcol(i).gt.0.d0.and.dtausub.gt.0.d0) then
c       odtot=odtot+odcoll(i)+odcoli(i)+odcolbc(i)
c       if((odcoll(i)+odcoli(i)+odcolbc(i)).gt.0.d0.and.
c    &      dtausub.gt.0.d0) then
c         if(odtot.le.dtausub) then
c           jadsub(k)=1
c           jadsub(k-1)=1
c         elseif(odtot.gt.dtausub) then
c           jadsub(k)=1
c           jadsub(k-1)=0
c           do j=1,2*(i-1)
c             jadsub(j)=0
c           enddo
c           go to 20
c         endif
c       endif
c     enddo
c 20   continue

      return
      end


      subroutine rd_prof(nj2,namfil)
C-----------------------------------------------------------------------
c  Routine to input T and O3 reference profiles
C-----------------------------------------------------------------------
      use raqmschem_pmgrid_mod, only : masterproc,iam
#ifdef DOMPI
      use mpishorthand, only : mpireal,mpicom
#endif
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
      use raqmschemcomm_mod, only : np,ns,nw,mx,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

      integer ia, i, m, l, lat, mon, ntlats, ntmons, n216, nj2
      real*8  ofac,ofak   
c     character*11 namfil   !original
!      character*26 namfil   !hyl
      character*(*) namfil   !hyl
      integer localrc
c
      if(masterproc)then
        open(NJ2,file=namfil,status='unknown')
        read(NJ2,'(A)') TITLE0
        read(NJ2,'(2I5)') NTLATS,NTMONS
        write(6,'(1X,A)') TITLE0
        write(6,1000) NTLATS,NTMONS
        N216 = MIN0(216, NTLATS*NTMONS)
        do IA=1,N216
          read(NJ2,'(1X,I3,3X,I2)') LAT, MON
          M = MIN(12, MAX(1, MON))
          L = MIN(18, MAX(1, (LAT+95)/10))
          read(NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
          read(NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
        enddo
        close(NJ2)
c
c       Extend climatology to 100 km
        ofac=exp(-2.d5/ZZHT)
        do i=32,51
          ofak=ofac**(i-31)
          do m=1,ntmons
            do l=1,ntlats
              oref(i,l,m)=oref(31,l,m)*ofak
            enddo
          enddo
        enddo
        do l=1,ntlats
          do m=1,ntmons
            do i=42,51
              tref(i,l,m)=tref(41,l,m)
            enddo
          enddo
        enddo
      endif
#ifdef DOMPI
!      call mpibcast(tref,51*18*12,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(tref,rc=localrc)
!      call mpibcast(oref,51*18*12,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(oref,rc=localrc)
#ifdef DIAGCHEM
      if(iam.eq.1)then
        write(6,*)'tref',tref
      endif
#endif
#endif
c
      return
 1000 format(1x,'Data: ',i3,' Lats x ',i2,' Months')
      end


      subroutine set_aer
C-----------------------------------------------------------------------
c  Set aerosol/cloud types and define black carbon profile
C-----------------------------------------------------------------------
c     MX       Number of different types of aerosol to be considered
c     MIEDX    Index of aerosol types in jv_spec.dat - hardwire in here
C-----------------------------------------------------------------------
      use raqmschem_pmgrid_mod, only : iam
      use raqmschemcomm_mod, only : np,ns,nw,mx,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

      integer i
c
c  Initialise aerosol index
      do i=1,MX
        MIEDX(i) = 0
      enddo
c
c  Select Aerosol/Cloud types to be used - define types here
      MIEDX(1) =  3    !  Black carbon absorber
      MIEDX(2) = 10    !  Water Cloud (Deirmenjian 8 micron)
      MIEDX(3) = 14    !  Irregular Ice Cloud (Mishchenko)
CCCMNN 0308  Dust size bins (5 bins, 18,19,20,21,21)
      do i=1,4
       MIEDX(3+i) = 17+i
       end do
       MIEDX(8) = 21   ! Dust bin 5 parameters same as for bin 4
CCCMNN 0308  Aerosols  Sulfate, BC, OC, Seasalt(a),and Seasalt (c)
CCCMNN       highest relative humidity is 0.9
CCC
      do i=1,5
CCCCCC  Sulfate
       MIEDX(8+i) = 21+i
CCCCCC  BC
       MIEDX(13+i) = 28+i
CCCCCC  OC
       MIEDX(18+i) = 35+i
CCCCCC  Seasalt(A)
       MIEDX(23+i) = 42+i
CCCCCC  Seasalt(C)
       MIEDX(28+i) = 49+i
       end do
c
c  Ensure all 'MX' types are valid selections

      do i=1,MX

c comment out this printout for now (hyl,02/19/03)
c "Using Aerosol type:  3 ABSRB  = fully absor
c  Using Aerosol type: 10 W_C08 = water cloud
c  Using Aerosol type: 14 Ice-I = irregular ic"
c       write(6,1000) MIEDX(i),TITLEA(MIEDX(i))

        if(MIEDX(i).gt.NAA.or.MIEDX(i).le.0) then
          write(6,*)'i',i,'mx',mx
          write(6,1200)iam, MIEDX(i),NAA
          call flush(6)
          write(0,*)'error bb'
          call flush(0)
          call killit('error bb')
!          stop
        endif
      enddo
c
c Approximate Black Carbon up to 10 km; surface 200 ng/m3  (Liousse et al)
c Scale: 1 ng/m3 = 1.0d-15 g/cm3 (1.0d-11 g/m2/cm as BREF is in cm))
c
c Simple place-holder profile
      do i=1,51
        BREF(i)=10.d0*1.0d-11
        if(i.gt.6) BREF(i)=0.d0
      enddo
c
      return
 1000 format('Using Aerosol type: ',i2,1x,a)
 1200 format(
     *i3,'Aerosol type ',i10,' unsuitable; supplied values must be ',
     $       'between 1 and ',i10)
      end


#ifdef JRATETJ
      SUBROUTINE JRATET(lbot, t, iday,iat,jat)       
#else
      SUBROUTINE JRATET(lbot, t, iday)       
#endif
C-----------------------------------------------------------------------
c  Calculate and print J-values. Note that the loop in this routine 
c  only covers the jpnl levels actually needed by the CTM.
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (hyl, 11/13/02)
C
C  Variable  Type    Dimensn Units   Description
C  --------  ----    ------- -----   -----------
C  T         dble     [lpar]  [K]    Vertical temperature profile
C  IDAY      int        -      -     Day of Year (0-365 or 0-366)
C-----------------------------------------------------------------------
c
c     FFF    Actinic flux at each level for each wavelength bin
c     QQQ    Cross sections for species (read in in RD_TJPL)
c     SOLF   Solar distance factor, for scaling; normally given by:
c                      1.0-(0.034*cos(real(iday-186)*2.0*pi/365.))
c                      Assumes aphelion day 186, perihelion day 3.
c     TQQ    Temperatures at which QQQ cross sections supplied
c
C-----------------------------------------------------------------------
#ifdef JVALUEJ
!      use ijprn
      use raqmschem_pmgrid_mod, only : iprn,jprn,kprn
#endif
      use raqmschemcomm_mod, only : np,ns,nw,mx,fastj,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

Chyl
C=============== INPUT PARAMETERS ======================================
      integer   iday, lbot
      real*8  t(lbot)      

C=============== LOCAL VARIABLES =======================================
#ifdef JRATETJ
      integer iat,jat
#endif
      integer i, j, k
      real*8 qo2tot, qo3tot, qo31d, qo33p, qqqt
      real*8 xseco2, xseco3, xsec1d, solf, tfact

c hyl, 11/12/02
      integer L
C     Parameters for Solar distance compensation
      real*8  PI
      PARAMETER (PI=3.14159265358979324D0)

C     Physical constants
      REAL*8  Na, R
      PARAMETER (Na=6.02217d23, R=8.3143d0)
#ifdef DOINT
      real zero8,one8
      parameter (zero8=0.0,one8=1.0)
#endif

C     Scale actinic flux (FFF) by Solar distance factor (SOLF)
C     solf=1.d0-(0.034d0*cos(dble(iday-172)*2.d0*pi/365.d0))   !GEOS-CHEM
      solf=1.d0-(0.034d0*cos(dble(iday-186)*2.d0*pi/365.d0))   !UCI Fast-J2
C
C----------------------------------------------------------------------
C If you want to set SOLF = 1.0 for testing, uncomment the next line
C      SOLF = 1d0
C----------------------------------------------------------------------
C
#ifdef JRATETJ
       if(iat.eq.iprn.and.jat.eq.jprn+1)then
         print *,'top jratet lbot',lbot
       endif
#endif
      do I=1,lbot                         !do I=1,jpnl
#ifdef JRATETJ
       if(iat.eq.iprn.and.jat.eq.jprn+1.and.i.eq.1)then
          write(6,*)'k',i,'t',t(i),'nw1',nw1,nw2,'njval',njval
       endif
#endif
!       VALJ(1) = 0.d0
!       VALJ(2) = 0.d0
!       VALJ(3) = 0.d0
       fastj%valj(1:3)=0.d0
       do K=NW1,NW2                       ! Using model 'T's here
#ifdef DOINT
         QO2TOT= XSECO2(K,T(I))          !hyl
#else
         QO2TOT= XSECO2(K,dble(T(I)))          !hyl
#endif
!         VALJ(1) = VALJ(1) + QO2TOT*FFF(K,I)
!         VALJ(1) = VALJ(1) + QO2TOT*fastj%FFF(K,I)
         fastj%VALJ(1) = fastj%VALJ(1) + QO2TOT*fastj%FFF(K,I)
!         if(isnan(fastj%valj(1)))then
!            write(6,*)'k',k,nw1,nw2,'T',i,t(i),'qo2tot',QO2TOT,'FFF',fastj%FFF(K,I)
         !endif
#ifdef JRATETJ
       if(iat.eq.iprn.and.jat.eq.jprn+1.and.i.eq.1)then
!          write(6,*)'valj1',valj(1),'qo2tot',qo2tot,'fff',fff(k,i),
!          write(6,*)'valj1',valj(1),'qo2tot',qo2tot,'fff',fastj%fff(k,i),
          write(6,*)'valj1',fastj%valj(1),'qo2tot',qo2tot,'fff',fastj%fff(k,i),
     * 'k',k,'t',t(i)
       endif
#endif

#ifdef DOINT
         QO3TOT= XSECO3(K,T(I))          !hyl
         QO31D = XSEC1D(K,T(I))*QO3TOT   !hyl
#else
         QO3TOT= XSECO3(K,dble(T(I)))          !hyl
         QO31D = XSEC1D(K,dble(T(I)))*QO3TOT   !hyl
#endif
         QO33P = QO3TOT - QO31D
!         VALJ(2) = VALJ(2) + QO33P*FFF(K,I)
!         VALJ(3) = VALJ(3) + QO31D*FFF(K,I)
!         VALJ(2) = VALJ(2) + QO33P*fastj%FFF(K,I)
!         VALJ(3) = VALJ(3) + QO31D*fastj%FFF(K,I)
         fastj%VALJ(2) = fastj%VALJ(2) + QO33P*fastj%FFF(K,I)
!         VALJ(3) = VALJ(3) + QO31D*fastj%FFF(K,I)
         fastj%VALJ(3) = fastj%VALJ(3) + QO31D*fastj%FFF(K,I)
       enddo
C------Calculate remaining J-values with T-dep X-sections 
       do J=4,NJVAL
!         VALJ(J) = 0.d0
         fastj%VALJ(J) = 0.d0
         TFACT = 0.d0
         L = jpdep(J)    !for pressure dependencies (hyl)
#ifdef JRATETJ
       if(iat.eq.iprn.and.jat.eq.jprn+1.and.i.eq.1)then
        print *,'j',j,'l',l,'tqq2',tqq(2,j),tqq(1,j)
       endif
#endif
#ifdef DOINT
         if(TQQ(2,J).gt.TQQ(1,J)) TFACT = max(zero8,min(one8,
     $        (T(I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))            !hyl
#else
         if(TQQ(2,J).gt.TQQ(1,J)) TFACT = max(0.d0,min(1.d0,
     $        (T(I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))            !hyl
#endif
#ifdef JRATETJ
       if(iat.eq.iprn.and.jat.eq.jprn+1.and.i.eq.1)then
         print *,'tfact',tfact
       endif
#endif
         do K=NW1,NW2
           QQQT = QQQ(K,1,J-3) + (QQQ(K,2,J-3) - QQQ(K,1,J-3))*TFACT
           if (L.eq.0) then   !for pressure dependencies (hyl)
!             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)    
!             VALJ(J) = VALJ(J) + QQQT*fastj%FFF(K,I)    
             fastj%VALJ(J) = fastj%VALJ(J) + QQQT*fastj%FFF(K,I)    
           else               !for pressure dependencies (hyl)

C----------------------------------------------------------------------
C Prior to 9/17/99
C Original form for acetaldehyde P-dep -- believed to be incorrect (pjc)
C             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)*
C     $                   (1.d0+zpdep(K,L)*(pj(i)+pj(i+1))*0.5d0)
C----------------------------------------------------------------------
C Essentially the change is the replacement of the factor
C
C   (1 + a P)     with               1
C                           ---------------------
C                             (1 + b density)
C
C where a and b are constants, P is pressure, and density is the 
C density of air in molec-cm(-3)   (pjc, 9/17/99)
C----------------------------------------------------------------------
cfor pressure dependencies (hyl)
c              VALJ(J)=VALJ(J)+QQQT*FFF(K,I)/(1 +
c     $                 (zpdep(K,L)*Na*1d-6 /(R*T(I))) *
c     $                 (pj(i)+pj(i+1))*0.5d0*1d2)
cccccccc
CCCMNN   TOTAL DENSITY DEPENDENCE FOR ACETONE
C
!            IF (L .EQ. 1) VALJ(J)=VALJ(J)+QQQT*FFF(K,I)/(1.0+
!            IF (L .EQ. 1) VALJ(J)=VALJ(J)+QQQT*fastj%FFF(K,I)/(1.0+
            IF (L .EQ. 1) fastj%VALJ(J)=fastj%VALJ(J)+QQQT*fastj%FFF(K,I)/(1.0+
     $              (ZPDEP(K,L)* 7.2435423D18/T(I))*
!     $              (pj(i)+pj(i+1))*0.5d0)
     $              (fastj%pj(i)+fastj%pj(i+1))*0.5d0)
C
CCCMNN   PRESURE DEPENDENCE FOR MGLY
C
!            IF (L .EQ. 2) VALJ(J) = VALJ(J)+QQQT*FFF(K,I)/(1.0+
!            IF (L .EQ. 2) VALJ(J) = VALJ(J)+QQQT*fastj%FFF(K,I)/(1.0+
            IF (L .EQ. 2) fastj%VALJ(J) = fastj%VALJ(J)+QQQT*fastj%FFF(K,I)/(1.0+
!     $                   zpdep(K,L)*(pj(i)+pj(i+1))*0.5d0)
     $                   zpdep(K,L)*(fastj%pj(i)+fastj%pj(i+1))*0.5d0)
C
CCCMNN   TOTAL NUMBER DENSITY DEPENDENCE FOR MVK
C
!            IF (L .EQ. 3) VALJ(J) = VALJ(J)+QQQT*FFF(K,I)/(1.0+
!            IF (L .EQ. 3) VALJ(J) = VALJ(J)+QQQT*fastj%FFF(K,I)/(1.0+
            IF (L .EQ. 3) fastj%VALJ(J) = fastj%VALJ(J)+QQQT*fastj%FFF(K,I)/(1.0+
     $              (ZPDEP(K,L)* 7.2435423D18/T(I))*
!     $              (pj(i)+pj(i+1))*0.5d0)
     $              (fastj%pj(i)+fastj%pj(i+1))*0.5d0)

           endif

c------Additional code for pressure dependencies
c           if(jpdep(J).ne.0) then
c             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)*
c     $                   (zpdep(K,L)*(pj(i)+pj(i+1))*0.5d0)
c           endif
#ifdef JRATETJ
       if(iat.eq.iprn.and.jat.eq.jprn+1.and.i.eq.1)then
           write(6,*)'k',k,'j',j,'valj',valj(j),'l',l
       endif
#endif

         enddo
       enddo
       do j=1,jppj
!         zj(i,j)=VALJ(jind(j))*jfacta(j)*SOLF
         fastj%zj(i,j)=fastj%VALJ(jind(j))*jfacta(j)*SOLF
!         if(isnan(fastj%zj(i,j)))then
!            write(6,*)'i',i,j,'valj',jind(j),'solf',solf,'jfacta',jfacta(j),'valj',fastj%VALJ(jind(j))
!          endif
#ifdef JRATETJ
       if(iat.eq.iprn.and.jat.eq.jprn+1.and.i.eq.1.and.j.eq.1)then
           print *,'set zj',i,j,zj(i,j),'jind',jind(j),'valj',
     *  valj(jind(j)),'jfacta',jfacta(j),'solf',solf
         endif
#endif
       enddo
cc     write(6,'(I5,1P,7E10.3/(5X,7E10.3))') I, (VALJ(J), J=1,NJVAL)
c      write(6,'(I5,1P,4E10.3)') I, VALJ(2),VALJ(3),VALJ(25),VALJ(27)
      enddo
      return
      end subroutine jratet


      SUBROUTINE RD_TJPL(NJ1,NAMFIL)
C-----------------------------------------------------------------------
c  Read in wavelength bins, solar fluxes, Rayleigh parameters, temperature-
c  dependent cross sections and Rayleigh/aerosol scattering phase functions
c  with temperature dependences. Current data originates from JPL'97
C-----------------------------------------------------------------------
c
c     NAMFIL   Name of spectral data file (jv_spec.dat)
c     NJ1      Channel number for reading data file
c     NJVAL    Number of species to calculate J-values for
c     NWWW     Number of wavelength bins, from NW1:NW2
c     WBIN     Boundaries of wavelength bins
c     WL       Centres of wavelength bins - 'effective wavelength'
c     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
c     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
c     QBC      Black Carbon absorption extinct. (specific cross-sect.) (m2/g)
c     QO2      O2 cross-sections
c     QO3      O3 cross-sections
c     Q1D      O3 => O(1D) quantum yield
c     TQQ      Temperature for supplied cross sections
c     QQQ      Supplied cross sections in each wavelength bin (cm2)
c     NAA      Number of categories for scattering phase functions
c     QAA      Aerosol scattering phase functions
c     NK       Number of wavelengths at which functions supplied (set as 4)
c     WAA      Wavelengths for the NK supplied phase functions
c     PAA      Phase function: first 8 terms of expansion
c     RAA      Effective radius associated with aerosol type
c     SSA      Single scattering albedo
c
c     npdep    Number of pressure dependencies
c     zpdep    Pressure dependencies by wavelength bin
c     jpdep    Index of cross sections requiring pressure dependence
c     lpdep    Label for pressure dependence
c
C-----------------------------------------------------------------------
      use raqmschem_pmgrid_mod, only : masterproc,iam
#ifdef DOMPI
      use mpishorthand, only : mpireal,mpicom,mpichar,mpiint
#endif
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
      use raqmschemcomm_mod, only : ns,nw,np,mx
      use raqmschemcomm_mod, only : fastj,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"
      integer i, j, k, iw, nk, nqqq, nwww, nj1
      character*7  lpdep(3)
c     character*11 NAMFIL   !original
!      character*26 NAMFIL   !hyl
      character*(*) NAMFIL   !hyl
      integer localrc
      logical exist
      do J=1,NS
        do K=1,3
          TQQ(K,J) = 0.d0
        enddo
      enddo
C-------------spectral data---------------------------------------------
      if(masterproc)then
        inquire(file=namfil,exist=exist)
        if(.not.exist)then
            write(6,*)'error namfil does not exist ',trim(namfil)
            call flush(6)
            call killit('rd_tjpl')
        endif
        open(NJ1, FILE=NAMFIL,status='unknown')
        read(NJ1,'(A)') TITLE0
        write(6,'(1X,A)') TITLE0
        read(NJ1,'(10X,14I5)') NJVAL,NWWW,NW1,NW2
        if(NJVAL.gt.NS) then
          write(6,300) NJVAL,NS
          call killit('rd_tjpl')
        endif
C------------NQQQ = no. additional J-values from X-sects (O2,O3P,O3D+NQQQ)
        NQQQ = NJVAL-3
c      read(NJ1,102) (WBIN(IW),IW=1,NWWW)
c      read(NJ1,102) (WBIN(IW+1),IW=1,NWWW)
        read(NJ1,102) (WL(IW),IW=1,NWWW)
        read(NJ1,102) (FL(IW),IW=1,NWWW)
        read(NJ1,102) (QRAYL(IW),IW=1,NWWW)
        read(NJ1,102) (QBC(IW),IW=1,NWWW)   !  From Loiusse et al. [JGR, 1996]
c
C---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
        do K=1,3
          read(NJ1,103) TITLEJ(K,1),TQQ(K,1), (QO2(IW,K),IW=1,NWWW)
        enddo
        do K=1,3
          read(NJ1,103) TITLEJ(K,2),TQQ(K,2), (QO3(IW,K),IW=1,NWWW)
        enddo
        do K=1,3
          read(NJ1,103) TITLEJ(K,3),TQQ(K,3), (Q1D(IW,K),IW=1,NWWW)
        enddo
chyl
c       do K=1,3
c         write(6,200) titlej(1,k),(tqq(i,k),i=1,3)
c       enddo
c
C---Read remaining species:  X-sections at 2 T's
        do J=1,NQQQ
          read(NJ1,103) TITLEJ(1,J+3),TQQ(1,J+3),(QQQ(IW,1,J),IW=1,NWWW)
          read(NJ1,103) TITLEJ(2,J+3),TQQ(2,J+3),(QQQ(IW,2,J),IW=1,NWWW)
#ifdef PRINTTITLE
          write(6,200) titlej(1,j+3),(tqq(i,j+3),i=1,2)
#endif
        enddo
        read(NJ1,'(A)') TITLE0

        npdep=0
cc---Pressure dependencies - see utils/pdep if required
        read(NJ1,104) npdep
        do k=1,npdep
          read(NJ1,105) lpdep(k),(zpdep(iw,k),iw=1,nwww)
#ifdef PRINTTITLE
          if(masterproc)write(6,201)  lpdep(k),(zpdep(iw,k),iw=1,nwww)
#endif
        enddo
        read(NJ1,'(A)') TITLE0

c Comment out for now. zpdep for 18 wavelengths are required,
c but only 7 of them are available. Uncomment these lines if
c pressure-dependencies are included in jv_spec.dat. (hyl,11/20/02)
c---Pressure dependencies
c      read(NJ1,104) npdep
c      do k=1,npdep
c        read(NJ1,105) lpdep(k),(zpdep(iw,k),iw=1,nwww)
c        write(6,201)  lpdep(k),(zpdep(iw,k),iw=1,nwww)
c      enddo
c      read(NJ1,'(A)') TITLE0

c
C---Read aerosol phase functions:
        read(NJ1,'(A10,I5,/)') TITLE0,NAA
!        write(6,*)'title0',title0
!        write(6,*)'NAA=',NAA
!        call flush(6)
        if(NAA.gt.NP) then
          write(6,350) NAA
          call flush(6)
          write(0,*)'error cd'
          call flush(0)
          call killit('error cd')
!          stop
        endif
      endif
      NK=4        ! Fix number of wavelengths at 4
      if(masterproc)then
        do j=1,NAA
          read(NJ1,110) TITLEA(j)
          do k=1,NK
            read(NJ1,*) WAA(k,j),QAA(k,j),RAA(k,j),SSA(k,j),
     &                                             (PAA(i,k,j),i=1,8)
          enddo
        enddo
c
c comment out for now (hyl,02/19/03)
c     write(6,*) 'Aerosol phase functions & wavelengths'
c     do J=1,NAA
c       write(6,'(1x,A8,I2,A,9F8.1)')
c    $                   TITLEA(J),J,'  wavel=',(WAA(K,J),K=1,NK)
c       write(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
c      enddo
c
c---Zero index arrays
        do j=1,jppj
          jind(j)=0
        enddo
        do j=1,NJVAL
          jpdep(j)=0
        enddo
c
C---Set mapping index
        do j=1,NJVAL
          do k=1,jppj
            if(jlabel(k).eq.titlej(1,j)) jind(k)=j
          enddo
          do k=1,npdep
            if(lpdep(k).eq.titlej(1,j)) jpdep(j)=k
          enddo
        enddo
        do k=1,jppj
c comment out for now (hyl, 02/19/03)
c       if(jfacta(k).eq.0.d0)
c    &             write(6,*) 'Not using photolysis reaction ',k
          if(jind(k).eq.0) then
            if(jfacta(k).eq.0.d0) then
              jind(k)=1
            else
              write(6,*) 'Which J-rate for photolysis reaction ',k,' ?'
              call flush(6)
              write(0,*)'error ee'
              call flush(0)
              call killit('error ee')
!              stop
            endif
          endif
        enddo
      endif
#ifdef DIAGCHEM
      if(masterproc)then
       write(6,*)'titlej master'
       call flush(6)
         do j=1,ns
       do i=1,3
       write(6,*)'titlej',titlej(i,j),' i ',i,j
       call flush(6)
       end do
       end do
      endif
#endif
      call raqmschem_comm_all_bcast(titlej,rc=localrc)
      call raqmschem_comm_all_bcast(wl,rc=localrc)
      call raqmschem_comm_all_bcast(fl,rc=localrc)
      call raqmschem_comm_all_bcast(qrayl,rc=localrc)
      call raqmschem_comm_all_bcast(qbc,rc=localrc)
      call raqmschem_comm_all_bcast(jfacta,rc=localrc)
      call raqmschem_comm_all_bcast(waa,rc=localrc)
      call raqmschem_comm_all_bcast(qaa,rc=localrc)
      call raqmschem_comm_all_bcast(raa,rc=localrc)
      call raqmschem_comm_all_bcast(ssa,rc=localrc)
      call raqmschem_comm_all_bcast(paa,rc=localrc)
      call raqmschem_comm_all_bcast(tqq,rc=localrc)
      call raqmschem_comm_all_bcast(qqq,rc=localrc)
      call raqmschem_comm_all_bcast(zpdep,rc=localrc)
      call raqmschem_comm_all_bcast(jpdep,rc=localrc)
      call raqmschem_comm_all_bcast(jind,rc=localrc)
      call raqmschem_comm_all_bcast(npdep,rc=localrc)
      call raqmschem_comm_all_bcast(njval,rc=localrc)
      call raqmschem_comm_all_bcast(nw1,rc=localrc)
      call raqmschem_comm_all_bcast(nw2,rc=localrc)
      call raqmschem_comm_all_bcast(naa,rc=localrc)
      call raqmschem_comm_all_bcast(qo2,rc=localrc)
      call raqmschem_comm_all_bcast(qo3,rc=localrc)
      call raqmschem_comm_all_bcast(q1d,rc=localrc)
!      write(6,*)'naa',naa,'after bcast'
!      call flush(6)
#ifdef DIAGCHEM
      if(iam.eq.1)then
!         do j=1,nw
         do j=1,3
           do i=1,nw
             write(6,*)'qo2',qo2(i,j),qo3(i,j),q1d(i,j)
           end do
         end do
         do j=1,ns
       do i=1,3
       write(6,*)'titlej',titlej(i,j),' i ',i,j
       end do
       end do
       write(6,*)'wl',wl
       write(6,*)'fl',fl
       write(6,*)'qrayl',qrayl
       write(6,*)'qbc',qbc
       write(6,*)'jfacta',jfacta
       write(6,*)'waa',waa
       write(6,*)'qaa',qaa
       write(6,*)'raa',raa
       write(6,*)'ssa',ssa
       write(6,*)'paa',paa
       write(6,*)'tqq',tqq
       write(6,*)'qqq',qqq
       write(6,*)'zpdep',zpdep
       write(6,*)'jpdep',jpdep
       write(6,*)'jind',jind
       write(6,*)'npdep',npdep,njval,nw1,nw2,naa
      endif
#endif
c      
C--------
  101 FORMAT(8E10.3)
  102 FORMAT((10X,6E10.3)/(10X,6E10.3)/(10X,6E10.3))
  103 FORMAT(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))  
c 103 FORMAT(A7,F3.0,7E10.3/(10X,7E10.3))
  104 FORMAT(13x,i2)
  105 FORMAT(A7,3x,6E10.3/(10X,6E10.3))
  110 format(3x,a20)
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(1x,' pr.dep:',a10,6(1pE10.3))
  300 format(' Number of x-sections supplied to Fast-J2: ',i3,/,
     &       ' Maximum number allowed (NS) only set to: ',i3,
     &       ' - increase in jv_cmn.h')
  350 format(' Too many phase functions supplied; increase NP to ',i2)
  400 format(' Not using Herzberg bin')
  420 format(' Using Herzberg bin for: ',10a7)
      close(NJ1)
      return
      end

      subroutine PRTATM(N, lbot, gmt, imon, rlat, sza_deg)
C-----------------------------------------------------------------------
C  Print out the atmosphere and calculate appropriate columns
C     N=1    Print out column totals only
C     N=2    Print out full columns
C     N=3    Print out full columns and climatology
C-----------------------------------------------------------------------
      use raqmschemcomm_mod, only : np,ns,nw,mx
      use raqmschemcomm_mod, only : fastj,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"
chyl
      real*8 gmt, rlat, sza_deg
      integer imon, lbot

      integer n, i, k, l, m
c     real*8 COLO3(NB),COLO2(NB),COLAX(MX,NB),ZKM,ZSTAR,PJC
c     real*8 climat(9),masfac,dlogp
      real*8 COLO3(lbot+1),COLO2(lbot+1),COLAX(MX,lbot+1),ZKM,ZSTAR,PJC
      real*8 climat(9),masfac,dlogp

      if(N.eq.0) return

C---Calculate columns, for diagnostic output only:
c      COLO3(NB) = DO3(NB)
c      COLO2(NB) = DM(NB)*0.20948d0
c      do K=1,MX
c        COLAX(K,NB) = AER(K,NB)
c      enddo
!      COLO3(lbot+1) = DO3(lbot+1) 
      COLO3(lbot+1) = fastj%DO3(lbot+1) 
!      COLO2(lbot+1) = DM(lbot+1)*0.20948d0
      COLO2(lbot+1) = fastj%DM(lbot+1)*0.20948d0
      do K=1,MX
!        COLAX(K,lbot+1) = AER(K,lbot+1)
        COLAX(K,lbot+1) = fastj%AER(K,lbot+1)
      enddo  

      do I=lbot,1,-1      !do I=NB-1,1,-1
!        COLO3(i) = COLO3(i+1)+DO3(i)
        COLO3(i) = COLO3(i+1)+fastj%DO3(i)
!        COLO2(i) = COLO2(i+1)+DM(i)*0.20948d0
        COLO2(i) = COLO2(i+1)+fastj%DM(i)*0.20948d0
        do K=1,MX
!          COLAX(k,i) = COLAX(k,i+1)+AER(k,i)
          COLAX(k,i) = COLAX(k,i+1)+fastj%AER(k,i)
        enddo
      enddo

chyl, 11/13/02
c     write(6,1200) ' GMT=',gmt,'  SZA=',sza_deg
c     write(6,1200) ' RAQMS O3-column(DU)=',COLO3(1)/2.687d16,
c    $              '  column aerosol @1000nm=',(COLAX(K,1),K=1,MX)

C---Print out atmosphere
      if(N.gt.1) then
chyl    write(6,1000) (' AER-X ','col-AER',k=1,mx)
        do I=lbot+1,1,-1         !do I=NB,1,-1
!          PJC = PJ(I)
          PJC = fastj%PJ(I)
!          ZKM =1.d-5*Z(I)
          ZKM =1.d-5*fastj%Z(I)
          ZSTAR = 16.d0*DLOG10(1000.d0/PJC)
chyl      write(6,1100) I,ZKM,ZSTAR,DM(I),DO3(I),1.d6*DO3(I)/DM(I),
c    $         TJ(I),PJC,COLO3(I),COLO2(I),(AER(K,I),COLAX(K,I),K=1,MX)
        enddo
      endif
c
C---Print out climatology
      if(N.gt.2) then
        do i=1,9
          climat(i)=0.d0
        enddo
         m = max(1,min(12,imon))
         l = max(1,min(18,(int(rlat)+99)/10))
        masfac=100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)
c       write(6,*) 'UCI Specified Climatology'
c       write(6,1000)
        do i=51,1,-1
          dlogp=10.d0**(-1.d0/16.d0)
          PJC = 1000.d0*dlogp**(2*i-2)
          climat(1) = 16.d0*DLOG10(1000.D0/PJC)
          climat(2) = climat(1)
          climat(3) = PJC*(1.d0/dlogp-dlogp)*masfac
          if(i.eq.1) climat(3)=PJC*(1.d0-dlogp)*masfac
          climat(4)=climat(3)*oref(i,l,m)*1.d-6
          climat(5)=oref(i,l,m)
          climat(6)=tref(i,l,m)
          climat(7)=PJC
          climat(8)=climat(8)+climat(4)
          climat(9)=climat(9)+climat(3)*0.20948d0
Chyl      write(6,1100) I,(climat(k),k=1,9)
        enddo
c       write(6,1200) 'UCI O3-column(DU)=',climat(8)/2.687d16
      endif
      return
 1000 format(5X,'Zkm',3X,'Z*',8X,'M',8X,'O3',6X,'f-O3',5X,'T',7X,'P',6x,
     $    'col-O3',3X,'col-O2',2X,10(a7,2x))
 1100 format(1X,I2,0P,2F6.2,1P,2E10.3,0P,F7.3,F8.2,F10.4,1P,10E9.2)
 1200 format(A,F8.1,A,10(1pE10.3))
      end

#if (defined JVALUEJ || defined DIAGOPMIE)
       subroutine jvalue(lbot, sza, sa,iat,jat)
#else
       subroutine jvalue(lbot, sza, sa)
#endif
c      SUBROUTINE JVALUE
C-----------------------------------------------------------------------
c  Calculate the actinic flux at each level for the current SZA value.
C        quit when SZA > 98.0 deg ==> tangent height = 63 km
C             or         99.                           80 km
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (hyl, 11/13/02)
C
C  Variable  Type    Dimensn Units   Description
C  --------  ----    ------- -----   -----------
C  SA        dble    -       -       Surface Albedo
C-----------------------------------------------------------------------
c
c     AVGF   Attenuation of beam at each level for each wavelength
c     FFF    Actinic flux at each desired level
c     WAVE   Effective wavelength of each wavelength bin
c     XQO2   Absorption cross-section of O2
c     XQO3   Absorption cross-section of O3
c
C-----------------------------------------------------------------------
#ifdef JVALUEJ
!      use ijprn
      use raqmschem_pmgrid_mod, only : iprn,jprn,kprn
      use iatjat
#endif
      use raqmschemcomm_mod, only : np,ns,nw,mx,fastj,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

chyl
C=============== INPUT PARAMETERS ======================================
      integer lbot
      real*8  sza, sa 

C=============== LOCAL VARIABLES =======================================
      integer j, k
#if (defined JVALUEJ || defined DIAGOPMIE)
      integer iat,jat
#endif
      real*8  wave, xseco3, xseco2
      real*8  AVGF(lbot),XQO3(lbot+1),XQO2(lbot+1)  
#ifdef JVALUEJ
      if(iat.eq.iprn.and.jat.eq.jprn+1)then
        doiat=.true.
      else
        doiat=.false.
      endif
#endif
C
      do J=1,lbot         !do J=1,jpnl
        do K=NW1,NW2
          !FFF(K,J) = 0.d0
          fastj%FFF(K,J) = 0.d0
        enddo
      enddo
c
c---SZA check
c      write(6,1000) SZA, RFLECT, (OD(nslon,nslat,j),j=1,lpar)
      if(sza.gt.szamax) GOTO 99
c
C---Calculate spherical weighting functions
cSPHERE was a common block name in vut_misc.f, so use a different name here.
cOtherwise, you got 'segmentation fault' error message (hyl, 11/20/02)
       CALL SPHERE_UCI (lbot)
c
C---Loop over all wavelength bins
      do K=NW1,NW2
        WAVE = WL(K)
#ifdef DOINT
        do J=1,lbot+1       
!          XQO3(J) = XSECO3(K,TJ(J))
          XQO3(J) = XSECO3(K,fastj%TJ(J))
        enddo
        do J=1,lbot+1     
!          XQO2(J) = XSECO2(K,TJ(J))
          XQO2(J) = XSECO2(K,fastj%TJ(J))
        enddo
#else
        do J=1,lbot+1       
!          XQO3(J) = XSECO3(K,dble(TJ(J)))
          XQO3(J) = XSECO3(K,dble(fastj%TJ(J)))
        enddo
        do J=1,lbot+1     
!          XQO2(J) = XSECO2(K,dble(TJ(J)))
          XQO2(J) = XSECO2(K,dble(fastj%TJ(J)))
        enddo
#endif
#ifdef JVALUEJ
        if(doiat)then
          write(6,*)'call opmie lbot',lbot,'k',k,'wave',wave
          do j=1,lbot+1
!            write(6,*)'xqo3',xqo3(j),'tj',tj(j)
            write(6,*)'xqo3',xqo3(j),'tj',fastj%tj(j)
            write(6,*)'xqo2',xqo2(j)
          end do
        endif
#endif
C----------------------------------------------------------
#if (defined JVALUEJ || defined DIAGOPMIE)
        CALL OPMIE(lbot, 2*(lbot+1), K,WAVE,XQO2,XQO3,AVGF,iat,jat)
#else
        CALL OPMIE(lbot, 2*(lbot+1), K,WAVE,XQO2,XQO3,AVGF)
#endif
C----------------------------------------------------------
        do J=1,lbot              !do J=1,jpnl
!          FFF(K,J) = FFF(K,J) + FL(K)*AVGF(J)
          fastj%FFF(K,J) = fastj%FFF(K,J) + FL(K)*AVGF(J)
!          if(isnan(fl(k)).or.isnan(avgf(j)))then
!              !!write(6,*)'fl',k,fl(k),'j',j,'avgf',avgf(j)
!          endif
#ifdef JVALUEJ
          if(iat.eq.iprn.and.jat.eq.jprn+1.and.k.eq.18.and.j.eq.1)then
!            write(6,*)'fff',fff(k,j),'fl',fl(k),'avgf',avgf(j)
            write(6,*)'fff',fastj%fff(k,j),'fl',fl(k),'avgf',avgf(j)
          endif
#endif
        enddo
      enddo
c
   99 continue
 1000 format('  SZA=',f6.1,' Reflectvty=',f6.3,' OD=',10(1pe10.3))
      return
      end subroutine jvalue

      FUNCTION XSECO3(K,TTT)
C-----------------------------------------------------------------------
c  Cross-sections for O3 for all processes interpolated across 3 temps
C-----------------------------------------------------------------------
      use raqmschemcomm_mod, only : np,ns,nw,mx,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

      integer k
      real*8 ttt, flint, xseco3
      XSECO3  = 
     F  FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))
      return
      end

      FUNCTION XSEC1D(K,TTT)
C-----------------------------------------------------------------------
c  Quantum yields for O3 --> O2 + O(1D) interpolated across 3 temps
C-----------------------------------------------------------------------
      use raqmschemcomm_mod, only : np,ns,nw,mx,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

      integer k
      real*8 ttt, flint, xsec1d
      XSEC1D =
     F  FLINT(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))
      return
      end

      FUNCTION XSECO2(K,TTT)
C-----------------------------------------------------------------------
c  Cross-sections for O2 interpolated across 3 temps; No S_R Bands yet!
C-----------------------------------------------------------------------
      use raqmschemcomm_mod, only : np,ns,nw,mx,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"

      integer k
      real*8 ttt, flint, xseco2
      XSECO2 =
     F  FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))
      return
      end


       REAL*8 FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
C-----------------------------------------------------------------------
c  Three-point linear interpolation function
C-----------------------------------------------------------------------
      implicit none

      real*8 TINT,T1,T2,T3,F1,F2,F3

      IF (TINT .LE. T2)  THEN
        IF (TINT .LE. T1)  THEN
          FLINT  = F1
        ELSE
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT .GE. T3)  THEN
          FLINT  = F3
        ELSE
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      return
      end


       SUBROUTINE SPHERE_UCI (lbot)
C-----------------------------------------------------------------------
c  Calculation of spherical geometry; derive tangent heights, slant path
c  lengths and air mass factor for each layer. Not called when
c  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
c  beam (where tangent height is below altitude J-value desired at).
C-----------------------------------------------------------------------
c
c     GMU     MU, cos(solar zenith angle)
c     RZ      Distance from centre of Earth to each point (cm)
c     RQ      Square of radius ratios
c     TANHT   Tangent height for the current SZA
c     XL      Slant path between points
c     AMF     Air mass factor for slab between level and level above
c
C-----------------------------------------------------------------------
      use raqmschemcomm_mod, only : np,ns,nw,mx,fastj,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"
chyl
      integer lbot

      integer i, j, k, ii
      real*8 airmas, gmu, xmu1, xmu2, xl, diff
      REAL*8 Ux,H,RZ(lbot+1),RQ(lbot+1),ZBYR        
#ifdef DOINT
      real zero8
      parameter (zero8=0.0)
#endif
c
c  Inlined air mass factor function for top of atmosphere
      AIRMAS(Ux,H) = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0-
     $         0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*H))/
     $                                             (1.0d0+0.625d0*H)))
c
!      GMU = U0
      GMU = fastj%U0

!      RZ(1)=RAD+Z(1)
      RZ(1)=RAD+fastj%Z(1)
      ZBYR = ZZHT/RAD
      DO 2 II=2,lbot+1           
!        RZ(II) = RAD + Z(II)
        RZ(II) = RAD + fastj%Z(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2
    2 CONTINUE
      IF (GMU.LT.0.0D0) THEN
!        TANHT = RZ(nlbatm)/DSQRT(1.0D0-GMU**2)
        fastj%TANHT = RZ(nlbatm)/DSQRT(1.0D0-GMU**2)
      ELSE
!        TANHT = RZ(nlbatm)
        fastj%TANHT = RZ(nlbatm)
      ENDIF
c
c  Go up from the surface calculating the slant paths between each level
c  and the level above, and deriving the appropriate Air Mass Factor
      DO 16 J=1,lbot+1           
        DO K=1,lbot+1         
!          AMF(K,J)=0.D0
          fastj%AMF(K,J)=0.D0
        ENDDO
c
c  Air Mass Factors all zero if below the tangent height
!        IF (RZ(J).LT.TANHT) GOTO 16
        IF (RZ(J).LT.fastj%TANHT) GOTO 16
c  Ascend from layer J calculating AMFs
        XMU1=ABS(GMU)
        DO 12 I=J,lbot        
          XMU2=DSQRT(1.0D0-RQ(I)*(1.0D0-XMU1**2))
          XL=RZ(I+1)*XMU2-RZ(I)*XMU1
!          AMF(I,J)=XL/(RZ(I+1)-RZ(I))
          fastj%AMF(I,J)=XL/(RZ(I+1)-RZ(I))
          XMU1=XMU2
   12   CONTINUE
c  Use function and scale height to provide AMF above top of model
!        AMF(lbot+1,J)=AIRMAS(XMU1,ZBYR)           !AMF(NB,J)=AIRMAS(XMU1,ZBYR)
        fastj%AMF(lbot+1,J)=AIRMAS(XMU1,ZBYR)           !AMF(NB,J)=AIRMAS(XMU1,ZBYR)
c
c  Twilight case - Emergent Beam
        IF (GMU.GE.0.0D0) GOTO 16
        XMU1=ABS(GMU)
c  Descend from layer J
        DO 14 II=J-1,1,-1
          DIFF=RZ(II+1)*DSQRT(1.0D0-XMU1**2)-RZ(II)
#ifdef DOINT
          if(II.eq.1) DIFF=max(DIFF,zero8)   ! filter
#else
          if(II.eq.1) DIFF=max(DIFF,0.d0)   ! filter
#endif
c  Tangent height below current level - beam passes through twice
          IF (DIFF.LT.0.0D0) THEN
            XMU2=DSQRT(1.0D0-(1.0D0-XMU1**2)/RQ(II))
            XL=ABS(RZ(II+1)*XMU1-RZ(II)*XMU2)
!            AMF(II,J)=2.d0*XL/(RZ(II+1)-RZ(II))
            fastj%AMF(II,J)=2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1=XMU2
c  Lowest level intersected by emergent beam
          ELSE
            XL=RZ(II+1)*XMU1*2.0D0
c            WTING=DIFF/(RZ(II+1)-RZ(II))
c            AMF(II,J)=(1.0D0-WTING)*2.D0**XL/(RZ(II+1)-RZ(II))
!            AMF(II,J)=XL/(RZ(II+1)-RZ(II))
            fastj%AMF(II,J)=XL/(RZ(II+1)-RZ(II))
            GOTO 16
          ENDIF
   14   CONTINUE
c
   16 CONTINUE
      RETURN
      END


#if (defined JVALUEJ || defined DIAGOPMIE)
      SUBROUTINE OPMIE(lbot, NC, KW,WAVEL,XQO2,XQO3,FMEAN,iat,jat)
#else
      SUBROUTINE OPMIE(lbot, NC, KW,WAVEL,XQO2,XQO3,FMEAN)
#endif
C-----------------------------------------------------------------------
C  NEW Mie code for J's, only uses 8-term expansion, 4-Gauss pts
C  Currently allow up to NP aerosol phase functions (at all altitudes) to
C  be associated with optical depth AER(1:NC) = aerosol opt.depth @ 1000 nm
C  
C  Pick Mie-wavelength with phase function and Qext:
C
C  01 RAYLE = Rayleigh phase
C  02 ISOTR = isotropic
C  03 ABSRB = fully absorbing 'soot', wavelength indep.
C  04 S_Bkg = backgrnd stratospheric sulfate (n=1.46,log-norm:r=.09um/sigma=.6)
C  05 S_Vol = volcanic stratospheric sulfate (n=1.46,log-norm:r=.08um/sigma=.8)
C  06 W_H01 = water haze (H1/Deirm.) (n=1.335, gamma:  r-mode=0.1um /alpha=2)
C  07 W_H04 = water haze (H1/Deirm.) (n=1.335, gamma:  r-mode=0.4um /alpha=2)
C  08 W_C02 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=2.0um /alpha=6)
C  09 W_C04 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=4.0um /alpha=6)
C  10 W_C08 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=8.0um /alpha=6)
C  11 W_C13 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=13.3um /alpha=6)
C  12 W_L06 = water cloud (Lacis) (n=1.335, r-mode=5.5um / alpha=11/3)
C  13 Ice-H = hexagonal ice cloud (Mishchenko)
C  14 Ice-I = irregular ice cloud (Mishchenko)
C
C  Choice of aerosol index MIEDX is made in SET_AER; optical depths are
C  apportioned to the AER array in SET_PROF
C
C-----------------------------------------------------------------------
C  FUNCTION RAYLAY(WAVE)---RAYLEIGH CROSS-SECTION for wave > 170 nm
C       WSQI = 1.E6/(WAVE*WAVE)
C       REFRM1 = 1.0E-6*(64.328+29498.1/(146.-WSQI)+255.4/(41.-WSQI))
C       RAYLAY = 5.40E-21*(REFRM1*WSQI)**2
C-----------------------------------------------------------------------
c
c     DTAUX    Local optical depth of each CTM level
c     PIRAY    Contribution of Rayleigh scattering to extinction
c     PIAER    Contribution of Aerosol scattering to extinction
c     TTAU     Optical depth of air vertically above each point (to top of atm)
c     FTAU     Attenuation of solar beam
c     POMEGA   Scattering phase function
c     FMEAN    Mean actinic flux at desired levels
c
c-----------------------------------------------------------------------
#ifdef JVALUEJ
      use iatjat
#endif
      use raqmschemcomm_mod, only : np,ns,nw,mx,fastj,nb,jppj
      implicit none

!      include 'cmn_fj.h'
#include "jv_cmn.h"
#include "jv_mie.h"

c NC=2*lbot+2 is number of levels in the fundamental Fast-J grid, and
c was previously defined in jv_cmn.h (hyl, 04/14/03)
#if (defined JVALUEJ || defined DIAGOPMIE)
      integer iat,jat
#endif
      integer lbot, NC     

      integer jndlev(lbot),jaddlv(NC),jaddto(NC+1)
      integer KW,km,i,j,k,l,ix,j1
      integer L2,L2L,L22,Lz,LZZ,NDZ,l2lev(nc+1)
      real*8 QXMIE(MX),XLAER(MX),SSALB(MX)
      real*8 xlo2,xlo3,xlray,xltau,zk,taudn,tauup,zk2
      real*8 WAVEL,XQO2(lbot+1),XQO3(lbot+1),FMEAN(lbot+1),
     &       POMEGAJ(2*M__,NC+1)
      real*8 DTAUX(lbot+1),PIRAY(lbot+1),PIAER(MX,lbot+1),
     &       TTAU(NC+1),FTAU(NC+1) 
      real*8 ftaulog,dttau,dpomega(2*M__)
      real*8 ftaulog2,dttau2,dpomega2(2*M__)
      real*8 TAUBTM,TAUTOP,FBTM,FTOP,POMEGAB(2*M__),ATAUA,ATAUZ
c
C---Pick nearest Mie wavelength, no interpolation--------------
                              KM=1
      if( WAVEL .gt. 355.d0 ) KM=2
      if( WAVEL .gt. 500.d0 ) KM=3
C     if( WAVEL .gt. 800.d0 ) KM=4  !drop the 1000 nm wavelength
c
C---For Mie code scale extinction at 1000 nm to wavelength WAVEL (QXMIE)
      do I=1,MX
        QXMIE(I) = QAA(KM,MIEDX(I))/QAA(4,MIEDX(I))
!        if(i.eq.2.or.i.eq.3)then
!          write(6,*)'qxmie',i,qxmie(i),'miedx',miedx(i),'km',km
!          call flush(6)
!        endif
        SSALB(I) = SSA(KM,MIEDX(I))
#ifdef JVALUEJ
        if(doiat)then
          write(6,*)'qxmie',i,qxmie(i),'km',km,'miedx',miedx(i),
     *   'qaa',qaa(km,miedx(i)),qaa(4,miedx(i))
          write(6,*)'ssalb',ssalb(i)
        endif
#endif
      enddo
c
C---Reinitialize arrays
      do j=1,NC+1
        ttau(j)=0.d0
        ftau(j)=0.d0
      enddo
c
C---Set up total optical depth over each CTM level, DTAUX
      J1 = NLBATM
      do J=J1,lbot+1                
!        XLO3=DO3(J)*XQO3(J)
!        XLO2=DM(J)*XQO2(J)*0.20948d0
!        XLRAY=DM(J)*QRAYL(KW)
        XLO3=fastj%DO3(J)*XQO3(J)
        XLO2=fastj%DM(J)*XQO2(J)*0.20948d0
        XLRAY=fastj%DM(J)*QRAYL(KW)
CCCMNN_JX  HAVE DROPPED FAST_J2 FORMULATION FOR FULL SCATTERING
CCC	if (WAVEL.le.291.d0) then
CCC        XLRAY=XLRAY * 0.57
CCC	endif
CCCCCC
c  Zero absorption for testing purposes
c        call NOABS(XLO3,XLO2,XLRAY,AER(1,j),RFLECT)
        do I=1,MX
!          XLAER(I)=AER(I,J)*QXMIE(I)
          XLAER(I)=fastj%AER(I,J)*QXMIE(I)
!          if(isnan(xlaer(i)))then
!            write(6,*)'aer',i,j,fastj%aer(i,j),'qxmie',i,qxmie(i)
!          endif
        enddo
c  Total optical depth from all elements
        DTAUX(J)=XLO3+XLO2+XLRAY
!        if(isnan(dtaux(j)))then
!           write(6,*)j,'xlo',xlo3,xlo2,xlray
!        endif
        do I=1,MX
          DTAUX(J)=DTAUX(J)+XLAER(I)
        enddo
c  Fractional extinction for Rayleigh scattering and each aerosol type
        PIRAY(J)=XLRAY/DTAUX(J)
        do I=1,MX
          PIAER(I,J)=SSALB(I)*XLAER(I)/DTAUX(J)
        enddo
      enddo
c
C---Calculate attenuated incident beam EXP(-TTAU/U0) and flux on surface
      do J=J1,lbot+1              
!        if(AMF(J,J).gt.0.0D0) then
        if(fastj%AMF(J,J).gt.0.0D0) then
          XLTAU=0.0D0
          do I=1,lbot+1            
!            XLTAU=XLTAU + DTAUX(I)*AMF(I,J)
            XLTAU=XLTAU + DTAUX(I)*fastj%AMF(I,J)
!            if(isnan(xltau))then
!              write(6,*)'dtaux',i,dtaux(i),'amf',fastj%AMF(I,J)
!            endif
          enddo
!          if(isnan(xltau))then
!             write(6,*)'xltau',xltau,'i',i,'j',j
!          endif
          if(XLTAU.gt. 76.d0) then   ! zero out flux at 1.0e-33
            FTAU(j)=0.d0
          else
#ifdef DOINT
            FTAU(J)=EXP(-XLTAU)
#else
            FTAU(J)=DEXP(-XLTAU)
#endif
          endif
        else
          FTAU(J)=0.0D0
        endif
      enddo
c
CCCMNN JX  HAVE DROPPED FAST_J2 FORMULATION FOR FULL SCATTERING
C---in UV region, use pseudo-Rayleigh absorption instead of scattering
CCC      if (WAVEL.le.291.d0) then   !<<< NEED TO FIND ENDIF
c
C---Accumulate attenuation for level centers
c
CCC	 do j=1,lbot           
CCC	  if (j.lt.J1) then
CCC	  FMEAN(J) = 0.d0
CCC	  else
CCC	  FMEAN(J) = sqrt(FTAU(J)*FTAU(J+1))
CCC	  endif
CCC	 enddo
CCC         return
c
C---in visible region, consider scattering
C---Define the scattering phase fn. with mix of Rayleigh(1) & Mie(MIEDX)
C   No. of quadrature pts fixed at 4 (M__), expansion of phase fn @ 8
CCC      else
      N = M__
      MFIT = 2*M__
      do j=j1,lbot+1         
        do i=1,MFIT
          pomegaj(i,j) = PIRAY(J)*PAA(I,KM,1)
#ifdef JVALUEJ
        if(doiat.and.i.eq.2.and.j<3)then
       write(6,*)j,'pomegaja',pomegaj(i,j),'piray',piray(j),paa(i,km,1)
        endif
#endif
          do k=1,MX
            pomegaj(i,j) = pomegaj(i,j) + PIAER(K,J)*PAA(I,KM,MIEDX(K))
#ifdef JVALUEJ
        if(doiat.and.i.eq.2.and.j<3)then
       write(6,*)j,'pomegajk',k,pomegaj(i,j),'piaer',piaer(k,j),'paa',
     *  paa(i,km,miedx(k)),km,miedx(k)
        endif
#endif
          enddo
        enddo
      enddo
c
C------------------------------------------------------------------------
c  Take optical properties on CTM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of CTM
c  layers; the index of these layers is kept in the jndlev array.
C------------------------------------------------------------------------
c
c  Set lower boundary and levels to calculate J-values at 
      J1=2*J1-1
      do j=1,lbot        
        jndlev(j)=2*j
      enddo
c
c  Calculate column optical depths above each level, TTAU
      TTAU(NC+1)=0.0D0
      do J=NC,J1,-1
        I=(J+1)/2
        TTAU(J)=TTAU(J+1) + 0.5d0*DTAUX(I)
CCCMNN 083007
CCCMNN JX-- NEW, FAST-JX FIX, ADD SCATTERING AND SPHERICAL TO FAST-J2
C---HOW TO ALLOW FOR LARGE TOTAL TAUS IN UV (<291 NM) AND
C---  AVOID ADDING MORE LAYERS:
C        IF (TTAU(J) .GT. 50.D0 .AND. WAVEL .LT. 292.D0) THEN
cjaa Test: include 1 more wavelength band for 10x fires 
cjaa        IF (TTAU(J) .GT. 50.D0 .AND. WAVEL .LT. 300.D0) THEN
C         JADDLV(J) = 0
C         ELSE
C         jaddlv(j)=int(0.5d0*DTAUX(I)/dtaumax)
C         ENDIF
C
c  Subdivide cloud-top levels if required
c       if(jadsub(j).gt.0) then
c         jadsub(j)=min(jaddlv(j)+1,nint(dtausub))*(nint(dsubdiv)-1)
c         jaddlv(j)=jaddlv(j)+jadsub(j)
c       endif
CCCCCCCCCCCCCCCCC
      enddo
c
C---reflect flux from surface
c
!      if(U0.gt.0.D0) then
      if(fastj%U0.gt.0.D0) then
!        ZFLUX = U0*FTAU(J1)*RFLECT/(1.d0+RFLECT)
!        ZFLUX = U0*FTAU(J1)*fastj%RFLECT/(1.d0+fastj%RFLECT)
        ZFLUX = fastj%U0*FTAU(J1)*fastj%RFLECT/(1.d0+fastj%RFLECT)
      else
        ZFLUX = 0.d0
      endif
c
c  Calculate attenuated beam, FTAU, level boundaries then level centres
cMNN 083007   MOD V 5.7  FAST_JX
c---version 5.6 fix   (mp, 3/2007)
c---solar flux at interp-levels: use exp(-TAU/U0) if U0>0
!      if (U0 .gt. 0.0d0) then
      if (fastj%U0 .gt. 0.0d0) then
       FTAU(NC+1) = 1.0D0
       DO J = NC,1,-2
        I = J/2
!        FTAU(J) =  FTAU(J+1)*EXP((TTAU(J+1)-TTAU(J))/U0)
        FTAU(J) =  FTAU(J+1)*EXP((TTAU(J+1)-TTAU(J))/fastj%U0)
#ifdef CHECKNANOPT
       if(.not.ftau(j)<=0.0.and..not.ftau(j)>=0.0.and.
     *  abs(ftau(j))>1.e30)then
         write(6,*)'ftau nan ',j,ftau(j),'plus ',ftau(j+1),
!     * 'ttau',j,ttau(j:j+1),'u0',u0,'i',i,ftau(i)
     * 'ttau',j,ttau(j:j+1),'u0',fastj%u0,'i',i,ftau(i)
         write(6,*)'iat',iat,jat
         call flush(6)
       endif
#endif
        FTAU(J-1) = FTAU(I)
        ENDDO
c---old, pre-5.6 version dropped flux in deep clouds to zero (FBTM)
       ELSE
       FTAU(NC+1)=1.0D0
       FTAU(NC) = SQRT(1.0D0*FTAU(NB))
       FTAU(NC-1) = FTAU(NB)
       DO J=NC-3,1,-2
        I = (J+1)/2
        FTAU(J) = FTAU(I)
        FTAU(J+1) = SQRT(FTAU(I+1)*FTAU(I))
        ENDDO
       ENDIF

c     FTAU(NC+1)=1.0D0
c     do J=NC-1,J1,-2
c       I=(J+1)/2
c       FTAU(J)=FTAU(I)
c     enddo
c     do J=NC,J1,-2
c       FTAU(J)=sqrt(FTAU(J+1)*FTAU(J-1))
c     enddo
c
c  Calculate scattering properties, level centres then level boundaries
c  using an inverse interpolation to give correctly-weighted values
      do j=NC,2,-2
        do i=1,MFIT
          pomegaj(i,j) = pomegaj(i,j/2)
#ifdef JVALUEJ
        if(doiat.and.i.eq.2.and.j.eq.1)then
       write(6,*)'pomegajs',pomegaj(i,j),'j/2',j/2
        endif
#endif
        enddo
      enddo
      do j=J1+2,NC-1,2
        taudn = ttau(j-1)-ttau(j)
        tauup = ttau(j)-ttau(j+1)
        do i=1,MFIT
          pomegaj(i,j) = (pomegaj(i,j-1)*taudn + 
     $                    pomegaj(i,j+1)*tauup) / (taudn+tauup)
        enddo
      enddo
c  Define lower and upper boundaries
      do i=1,MFIT
        pomegaj(i,J1)   = pomegaj(i,J1+1)
        pomegaj(i,NC+1) = pomegaj(i,NC)
#ifdef JVALUEJ
        if(doiat.and.i.eq.2.and.j1.eq.1)then
       write(6,*)'pomegajl',pomega(i,j1)
        endif
#endif
      enddo
c
C------------------------------------------------------------------------
c  Calculate cumulative total and define levels we want J-values at.
c  Sum upwards for levels, and then downwards for Mie code readjustments.
c
c     jaddlv(i)   Number of new levels to add between (i) and (i+1)
c     jaddto(i)   Total number of new levels to add to and above level (i)
c     jndlev(j)   Level needed for J-value for CTM layer (j)
c
C------------------------------------------------------------------------
c
CCCMNN081007
C
      DO J=1,NC
!       JADDLV(J)=JADSUB(J)
       JADDLV(J)=fastj%JADSUB(J)
       ENDDO
C
c  Reinitialize level arrays
      do j=1,NC+1
        jaddto(j)=0
      enddo
c
      jaddto(J1)=jaddlv(J1)
      do j=J1+1,NC
        jaddto(j)=jaddto(j-1)+jaddlv(j)
      enddo
      if((jaddto(NC)+NC).gt.nl) then
         write(6,1500)  jaddto(NC)+NC, 'NL',NL
         print*,lbot,nc,wavel
         print*,jaddto
         print*,jaddlv
         print*,fastj%do3
         print*,fastj%pj  !pressure
         print*,dtaux
         print*,ttau
         call flush(6)
c  Recalculate added levels but limit the large refinements   jaa 20050325
         do j=1,NC+1
           jaddto(j)=0
           jaddlv(j)=min(jaddlv(j),100)
         enddo
         jaddto(J1)=jaddlv(J1)
         do j=J1+1,NC
           jaddto(j)=jaddto(j-1)+jaddlv(j)
         enddo
         if((jaddto(NC)+NC).gt.nl) then
           print*,'Still too many levels in photolysis code...'
           print*,jaddlv
           call flush(6)
           do j=1,NC+1
             jaddto(j)=0
             jaddlv(j)=min(jaddlv(j),20)
           enddo
           jaddto(J1)=jaddlv(J1)
           do j=J1+1,NC
             jaddto(j)=jaddto(j-1)+jaddlv(j)
           enddo
           if((jaddto(NC)+NC).gt.nl) then
             print*,'Really a problem in photolysis code...'
             print*,jaddlv
             call flush(6)
             write(0,*)'error ff'
             call flush(0)
             call killit('error ff')
!             stop
           endif
         endif
      endif

c      write(6,1300) jndlev
c      write(6,1300) jaddto
      do i=1,lbot      
        jndlev(i)=jndlev(i)+jaddto(jndlev(i)-1)
      enddo
cMNN 083007   (MOD ADOPTED FROM  V 5.7)
c
      L2LEV(1) = 1
      DO J=2,NC+1
       L2LEV(J) = L2LEV(J-1)+1+JADDLV(J-1)
       ENDDO
C
      jaddto(NC)=jaddlv(NC)
      do j=NC-1,J1,-1
        jaddto(j)=jaddto(j+1)+jaddlv(j)
      enddo
c      write(6,1300) jndlev
c      write(6,1300) jaddto
c
C---------------------SET UP FOR MIE CODE-------------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd point added at top of grid for unattenuated beam   (Z='inf')
c  
c        Surface:   TTAU(1)   now use ZTAU(2*NC+1)
c        Top:       TTAU(NC)  now use ZTAU(3)
c        Infinity:            now use ZTAU(1)
c
c  Mie scattering code only used from surface to level NC
C------------------------------------------------------------------------
C
c  Initialise all Fast-J2 optical property arrays
      do k=1,N__
        do i=1,MFIT
          pomega(i,k) = 0.d0
        enddo
        ztau(k) = 0.d0
        fz(k)   = 0.d0
      enddo
c
cMNN  083007
      NDZ = 2*NC + 2*JADDTO(1) + 1

c   Note that the successive sub-layers have the ratio in OD of ATAU
c      ATAUA = (ATAU - 1.d0)/ATAU     ! this is the limit for L22=>inf

      do L2 = 1,NC+1          ! L2 = index of CTM edge- and mid-layers
       L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = NDZ + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ) = TTAU(L2)
          FZ(LZ)   = FTAU(L2)
!          if(isnan(fz(lz)))then
!           write(6,*)'ftau',l2,ftau(l2),'lz',lz
!          endif
        do I=1,MFIT
          POMEGA(I,LZ) = POMEGAJ(I,L2)
        enddo
      enddo

c   Now go thru the pairs of L2 levels to see if we need JADD levels
      do L2 = 1,NC             ! L2 = index of CTM edge- and mid-layers
       L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)
        LZ  = NDZ + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
        L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels
      if (L22 .gt. 0) then
          TAUBTM = TTAU(L2)
          TAUTOP = TTAU(L2+1)
          FBTM   = FTAU(L2)
          FTOP   = FTAU(L2+1)
         do I = 1,MFIT
          POMEGAB(I) = POMEGAJ(I,L2)
         enddo
c---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
c---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM

          ATAUZ = exp(-log(TAUBTM/max(TAUTOP,ATAU0))/float(L22+1))

        do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt array
          ZTAU(LZZ) = TAUBTM * ATAUZ

          ATAUA=(TAUBTM-ZTAU(LZZ))/(TAUBTM-TAUTOP) !fraction from TAUBTM=>TAUTO


c---version 5.6 fix   (mp, 3/2007)
c---solar flux at interp-levels: use exp(TAU/U0) if U0>0, else scale by TAU
!          if (U0 .gt. 0.d0) then
          if (fastj%U0 .gt. 0.d0) then
!            FZ(LZZ) = FTOP * exp((TAUTOP-ZTAU(LZZ))/U0)
            FZ(LZZ) = FTOP * exp((TAUTOP-ZTAU(LZZ))/fastj%U0)
          else
            if (FBTM .lt. 1.d-32) then
              FZ(LZZ) = 0.d0
            else
              FZ(LZZ) = FBTM * (FTOP/FBTM)**ATAUA
            endif
          endif


          do I = 1,MFIT
            POMEGA(I,LZZ) = POMEGAB(I) +
     &               ATAUA*(POMEGAJ(I,L2+1)-POMEGAB(I))
          enddo
            TAUBTM       = ZTAU(LZZ)
            FBTM         = FZ(LZZ)
          do I = 1,MFIT
            POMEGAB(I) = POMEGA(I,LZZ)
          enddo

        enddo
       endif
      enddo

c   Now fill in the even points with simple interpolation in scatter arrays:
      do LZ = 2,NDZ-1,2
        ZTAU(LZ) = 0.5d0*(ZTAU(LZ-1)+ZTAU(LZ+1))
        FZ(LZ)   = sqrt(FZ(LZ-1)*FZ(LZ+1))
       do I=1,MFIT
        POMEGA(I,LZ) = 0.5d0*(POMEGA(I,LZ-1)+POMEGA(I,LZ+1))
       enddo
      enddo

      ND = NDZ
!      ZU0 = U0
      ZU0 = fastj%U0
!      ZREFL = RFLECT
      ZREFL = fastj%RFLECT
c
      if(ND .gt. N__) then
        print *, ' overflow of scatter arrays:',ND,N__
        call flush(6)
        write(0,*)'error gg'
        call flush(0)
        call killit('error gg')
!        stop
      endif
#ifdef JVALUEJ
        if(doiat)then
          do j=1,nc+1
          do i=1,2*M__
          write(6,*)'pomega',i,j,pomega(i,j),pomegaj(i,j)
          enddo
          end do
          do i=1,2*M__
            write(6,*)'pomegab',i,pomegab(i)
          end do
        endif
#endif


C-----------------------------------------
      CALL MIESCT
C-----------------------------------------
c  Accumulate attenuation for selected levels
      l=2*(NC+jaddto(J1))+3
      do j=1,lbot
        k=l-(2*jndlev(j))
        if(k.gt.ND-2) then
          FMEAN(j) = 0.d0
        else
          FMEAN(j) = FJ(k)
!          if(isnan(fmean(j)))then
!            write(6,*)'error opmie fmean j',j,'k',k,'fj',fj(k)
!            call flush(6)
!          endif
#ifdef JVALUEJ
          if(doiat)then
            write(6,*)'fj',k,fj(k),'j',j,'l',l,'j1',j1
          endif
#endif
        endif
      enddo
c
CCCMNN JX      endif                    ! FAST-J2 ENDIF 
c
      return
 1000 format(1x,i3,3(2x,1pe11.4),1x,i3)
 1300 format(1x,50(i3))
 1500 format(' Too many levels in photolysis code: need ',i5,' but ',a,
     $       ' dimensioned as ',i5)
      END


      subroutine EFOLD (F0, F1, N, F)
C-----------------------------------------------------------------------
C---  calculate the e-fold between two boundaries, given the value
C---     at both boundaries F0(x=0) = top, F1(x=1) = bottom.
C---  presume that F(x) proportional to exp[-A*x] for x=0 to x=1
C---          d2F/dx2 = A*A*F  and thus expect F1 = F0 * exp[-A]
C---           alternatively, could define A = ln[F0/F1]
C---  let X = A*x, d2F/dX2 = F
C---  assume equal spacing (not necessary, but makes this easier)
C---      with N-1 intermediate points (and N layers of thickness dX = A/N)
C---
C---  2nd-order finite difference:  (F(i-1) - 2F(i) + F(i+1)) / dX*dX = F(i)
C---      let D = 1 / dX*dX:
C
C  1  |   1        0        0        0        0        0   |    | F0 |
C     |                                                    |    | 0  |
C  2  |  -D      2D+1      -D        0        0        0   |    | 0  |
C     |                                                    |    | 0  |
C  3  |   0       -D      2D+1      -D        0        0   |    | 0  |
C     |                                                    |    | 0  |
C     |   0        0       -D      2D+1      -D        0   |    | 0  |
C     |                                                    |    | 0  |
C  N  |   0        0        0       -D      2D+1      -D   |    | 0  |
C     |                                                    |    | 0  |
C N+1 |   0        0        0        0        0        1   |    | F1 |
C      
C-----------------------------------------------------------------------
C  Advantage of scheme over simple attenuation factor: conserves total
C  number of photons - very useful when using scheme for heating rates.
C  Disadvantage: although reproduces e-folds very well for small flux
C  differences, starts to drift off when many orders of magnitude are
C  involved.
C-----------------------------------------------------------------------
      implicit none
      real*8 F0,F1,F(250)  !F(N+1)
      integer N
      integer I
      real*8 A,DX,D,DSQ,DDP1, B(101),R(101)
C
      if(F0.eq.0.d0) then
        do I=1,N
          F(I)=0.d0
        enddo
        return
      elseif(F1.eq.0.d0) then
        A = DLOG(F0/1.d-250)
      else
#ifdef DOINT
        A = LOG(F0/F1)
#else
        A = DLOG(F0/F1)
#endif
      endif
C
      DX = float(N)/A
      D = DX*DX
      DSQ = D*D
      DDP1 = D+D+1.d0
C
      B(2) = DDP1
      R(2) = +D*F0
      do I=3,N
        B(I) = DDP1 - DSQ/B(I-1)
        R(I) = +D*R(I-1)/B(I-1)
      enddo
      F(N+1) = F1
      do I=N,2,-1
        F(I) = (R(I) + D*F(I+1))/B(I)
      enddo
      F(1) = F0
      return
      end


      subroutine CH_PROF
C-----------------------------------------------------------------------
C  Check profiles to be passed to MIESCT
C-----------------------------------------------------------------------
      implicit none
#include "jv_mie.h"
      integer i,j
      write(6,1100) 'lev','ztau','fz  ','pomega( )'
      do i=1,ND
        if(ztau(i).ne.0.d0) then
          write(6,1200) i,ztau(i),fz(i),(pomega(j,i),j=1,8)
        endif
      enddo
      return
 1100 format(1x,a3,4(a9,2x))
 1200 format(1x,i3,11(1x,1pe10.3))
      end


      SUBROUTINE MIESCT
C-----------------------------------------------------------------------
C   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
C     Prather, 1974, Astrophys. J. 192, 787-792.
C         Sol'n of inhomogeneous Rayleigh scattering atmosphere. 
C         (original Rayleigh w/ polarization)
C     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
C         Raman scattering in the atmospheres of the major planets.
C         (first use of anisotropic code)
C     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
C         Chemistry of a polluted cloudy boundary layer,
C         (documentation of extension to anisotropic scattering)
C
C    takes atmospheric structure and source terms from std J-code
C    ALSO limited to 4 Gauss points, only calculates mean field!
C
C   mean rad. field ONLY (M=1)
C   initialize variables FIXED/UNUSED in this special version:
C   FTOP = 1.0 = astrophysical flux (unit of pi) at SZA, -ZU0, use for scaling
C   FBOT = 0.0 = external isotropic flux on lower boundary 
C   SISOTP = 0.0 = Specific Intensity of isotropic radiation incident from top
C
C   SUBROUTINES:  MIESCT              needs 'jv_mie.cmn'
C                 BLKSLV              needs 'jv_mie.cmn'
C                 GEN (ID)            needs 'jv_mie.cmn'
C                 LEGND0 (X,PL,N)
C                 MATIN4 (A)
C                 GAUSSP (N,XPT,XWT)
C-----------------------------------------------------------------------
#ifdef JVALUEJ
      use iatjat
#endif
      IMPLICIT NONE
#include "jv_mie.h"
      integer i, id, im
      real*8  cmeq1
C-----------------------------------------------------------------------
C---fix scattering to 4 Gauss pts = 8-stream
      CALL GAUSSP (N,EMU,WT)
C---solve eqn of R.T. only for first-order M=1
C      ZFLUX = (ZU0*FZ(ND)*ZREFL+FBOT)/(1.0d0+ZREFL)
      ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M=1
      DO I=1,N
        CALL LEGND0 (EMU(I),PM0,MFIT)
        DO IM=M,MFIT
          PM(I,IM) = PM0(IM)
        ENDDO
      ENDDO
C
      CMEQ1 = 0.25D0
      CALL LEGND0 (-ZU0,PM0,MFIT)
      DO IM=M,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      ENDDO
C
      CALL BLKSLV
C
      DO ID=1,ND,2
        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
!        if(isnan(fj(id)))then
!          write(6,*)'at  2593 ',id,'fz',fz(id)
!        endif

      ENDDO

      RETURN
      END

      SUBROUTINE BLKSLV
C-----------------------------------------------------------------------
C  Solves the block tri-diagonal system:
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
C-----------------------------------------------------------------------
      IMPLICIT NONE
#include "jv_mie.h"
      integer i, j, k, id
      real*8  sum
C-----------UPPER BOUNDARY ID=1
      CALL GEN(1)
      CALL MATIN4 (B)
      DO I=1,N
         RR(I,1) = 0.0d0
        DO J=1,N
          SUM = 0.0d0
         DO K=1,N
          SUM = SUM - B(I,K)*CC(K,J)
         ENDDO
         DD(I,J,1) = SUM
         RR(I,1) = RR(I,1) + B(I,J)*H(J)
        ENDDO
      ENDDO
C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      DO ID=2,ND-1
        CALL GEN(ID)
        DO I=1,N
          DO J=1,N
          B(I,J) = B(I,J) + A(I)*DD(I,J,ID-1)
          ENDDO
          H(I) = H(I) - A(I)*RR(I,ID-1)
        ENDDO
        CALL MATIN4 (B)
        DO I=1,N
          RR(I,ID) = 0.0d0
          DO J=1,N
          RR(I,ID) = RR(I,ID) + B(I,J)*H(J)
          DD(I,J,ID) = - B(I,J)*C1(J)
          ENDDO
        ENDDO
      ENDDO
C---------FINAL DEPTH POINT: ND
      CALL GEN(ND)
      DO I=1,N
        DO J=1,N
          SUM = 0.0d0
          DO K=1,N
          SUM = SUM + AA(I,K)*DD(K,J,ND-1)
          ENDDO
        B(I,J) = B(I,J) + SUM
        H(I) = H(I) - AA(I,J)*RR(J,ND-1)
        ENDDO
      ENDDO
      CALL MATIN4 (B)
      DO I=1,N
        RR(I,ND) = 0.0d0
        DO J=1,N
        RR(I,ND) = RR(I,ND) + B(I,J)*H(J)
        ENDDO
      ENDDO
C-----------BACK SOLUTION
      DO ID=ND-1,1,-1
       DO I=1,N
        DO J=1,N
         RR(I,ID) = RR(I,ID) + DD(I,J,ID)*RR(J,ID+1)
        ENDDO
       ENDDO
      ENDDO
C----------MEAN J & H
      DO ID=1,ND,2
        FJ(ID) = 0.0d0
       DO I=1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)
       ENDDO
      ENDDO
      DO ID=2,ND,2
        FJ(ID) = 0.0d0
       DO I=1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)*EMU(I)
       ENDDO
      ENDDO
C Output fluxes for testing purposes
c      CALL CH_FLUX
c
      RETURN
      END


      SUBROUTINE CH_FLUX
C-----------------------------------------------------------------------
C  Diagnostic routine to check fluxes at each level - makes most sense
C  when running a conservative atmosphere (zero out absorption in
C  OPMIE by calling the NOABS routine below)
C-----------------------------------------------------------------------
      IMPLICIT NONE
#include "jv_mie.h"
      integer I,ID
      real*8 FJCHEK(N__),FZMEAN
c
c  Odd (h) levels held as actinic flux, so recalculate irradiances
      DO ID=1,ND,2
        FJCHEK(ID) = 0.0d0
       DO I=1,N
        FJCHEK(ID) = FJCHEK(ID) + RR(I,ID)*WT(I)*EMU(i)
       ENDDO
      ENDDO
c
c  Even (j) levels are already held as irradiances
      DO ID=2,ND,2
       DO I=1,N
        FJCHEK(ID) = FJ(ID)
       ENDDO
      ENDDO
c
c  Output Downward and Upward fluxes down through atmosphere
      WRITE(6,1200)
      DO ID=2,ND,2
        FZMEAN=sqrt(FZ(ID)*FZ(ID-1))
        WRITE(6,1000) ID, ZU0*FZMEAN-2.0*(FJCHEK(id)-FJCHEK(id-1)),
     $                               2.0*(FJCHEK(id)+FJCHEK(id-1)),
     $                               2.0*(FJCHEK(id)+FJCHEK(id-1))/
     $                   (ZU0*FZMEAN-2.0*(FJCHEK(id)-FJCHEK(id-1)))
      ENDDO
      RETURN
 1000 FORMAT(1x,i3,1p,2E12.4,1x,0p,f9.4)
 1200 FORMAT(1x,'Lev',3x,'Downward',4x,'Upward',7x,'Ratio')
      END

      SUBROUTINE NOABS(XLO3,XLO2,XLRAY,BCAER,RFLECT)
C-----------------------------------------------------------------------
C  Zero out absorption terms to check scattering code. Leave a little
C  Rayleigh to provide a minimal optical depth, and set surface albedo
C  to unity.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      real*8 XLO3,XLO2,XLRAY,BCAER,RFLECT
      XLO3=0.d0
      XLO2=0.d0
      XLRAY=XLRAY*1.d-10
      BCAER=0.d0      
      RFLECT=1.d0
      RETURN
      END


      SUBROUTINE GEN(ID)
C-----------------------------------------------------------------------
C  Generates coefficient matrices for the block tri-diagonal system:
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
C-----------------------------------------------------------------------
      IMPLICIT NONE
#include "jv_mie.h"
      integer id, id0, id1, im, i, j, k, mstart
      real*8  sum0, sum1, sum2, sum3
      real*8  deltau, d1, d2, surfac
C---------------------------------------------
      IF(ID.EQ.1 .OR. ID.EQ.ND) THEN
C---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       IF(ID.GE.ND) ID1 = ID-1
       DO 10 I=1,N
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
        DO IM=M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
        DO IM=M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
         H(I) = 0.5d0*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         A(I) = 0.5d0*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
        DO J=1,I
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
         DO IM=M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         DO IM=M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         W(I,J) = - SUM1*WT(J)
         W(J,I) = - SUM1*WT(I)
         U1(I,J) = - SUM3*WT(J)
         U1(J,I) = - SUM3*WT(I)
          SUM0 = 0.5d0*(SUM0 + SUM2)
         B(I,J) = - SUM0*WT(J)
         B(J,I) = - SUM0*WT(I)
        ENDDO
         S(I,I) = S(I,I) + 1.0d0
         W(I,I) = W(I,I) + 1.0d0
         U1(I,I) = U1(I,I) + 1.0d0
         B(I,I) = B(I,I) + 1.0d0
   10  CONTINUE
       DO I=1,N
         SUM0 = 0.0d0
        DO J=1,N
         SUM0 = SUM0 + S(I,J)*A(J)/EMU(J)
        ENDDO
        C1(I) = SUM0
       ENDDO
       DO I=1,N
        DO J=1,N
          SUM0 = 0.0d0
          SUM2 = 0.0d0
         DO K=1,N
          SUM0 = SUM0 + S(J,K)*W(K,I)/EMU(K)
          SUM2 = SUM2 + S(J,K)*U1(K,I)/EMU(K)
         ENDDO
         A(J) = SUM0
         V1(J) = SUM2
        ENDDO
        DO J=1,N
         W(J,I) = A(J)
         U1(J,I) = V1(J)
        ENDDO
       ENDDO
       IF (ID.EQ.1) THEN
C-------------upper boundary, 2nd-order, C-matrix is full (CC)
        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25d0*DELTAU
        DO I=1,N
          D1 = EMU(I)/DELTAU
          DO J=1,N
           B(I,J) = B(I,J) + D2*W(I,J)
           CC(I,J) = D2*U1(I,J)
          ENDDO
          B(I,I) = B(I,I) + D1
          CC(I,I) = CC(I,I) - D1
C         H(I) = H(I) + 2.0d0*D2*C1(I) + D1*SISOTP
          H(I) = H(I) + 2.0d0*D2*C1(I)
          A(I) = 0.0d0
        ENDDO
       ELSE
C-------------lower boundary, 2nd-order, A-matrix is full (AA)
        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25d0*DELTAU
        SURFAC = 4.0d0*ZREFL/(1.0d0 + ZREFL)
        DO I=1,N
          D1 = EMU(I)/DELTAU
          H(I) = H(I) - 2.0d0*D2*C1(I)
           SUM0 = 0.0d0
          DO J=1,N
           SUM0 = SUM0 + W(I,J)
          ENDDO
           SUM0 = D1 + D2*SUM0
           SUM1 = SURFAC*SUM0
          DO J=1,N
           B(I,J) = B(I,J) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
          ENDDO
          B(I,I) = B(I,I) + D1
          H(I) = H(I) + SUM0*ZFLUX
          DO J=1,N
           AA(I,J) = - D2*U1(I,J)
          ENDDO
           AA(I,I) = AA(I,I) + D1
           C1(I) = 0.0d0
        ENDDO
       ENDIF
C------------intermediate points:  can be even or odd, A & C diagonal
      ELSE
        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = M + MOD(ID+1,2)
        DO I=1,N
          A(I) = EMU(I)/DELTAU
          C1(I) = -A(I)
           SUM0 = 0.0d0
          DO IM=MSTART,MFIT,2
           SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM0(IM)
          ENDDO
          H(I) = SUM0*FZ(ID)
          DO J=1,I
            SUM0 = 0.0d0
           DO IM=MSTART,MFIT,2
            SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM(J,IM)
           ENDDO
            B(I,J) =  - SUM0*WT(J)
            B(J,I) =  - SUM0*WT(I)
          ENDDO
          B(I,I) = B(I,I) + 1.0d0
        ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE LEGND0 (X,PL,N)
C---Calculates ORDINARY LEGENDRE fns of X (real) 
C---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      IMPLICIT NONE
      INTEGER N,I
      REAL*8 X,PL(N),DEN
C---Always does PL(2) = P[1]
        PL(1) = 1.D0
        PL(2) = X
        DO I=3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.D0/DEN) - PL(I-2)*(1.d0-1.D0/DEN)
        ENDDO
      RETURN
      END

      SUBROUTINE MATIN4 (A)
C-----------------------------------------------------------------------
C  invert 4x4 matrix A(4,4) in place with L-U decomposition (mjp, old...)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 A(4,4)
C---SETUP L AND U
      A(2,1) = A(2,1)/A(1,1)
      A(2,2) = A(2,2)-A(2,1)*A(1,2)
      A(2,3) = A(2,3)-A(2,1)*A(1,3)
      A(2,4) = A(2,4)-A(2,1)*A(1,4)
      A(3,1) = A(3,1)/A(1,1)
      A(3,2) = (A(3,2)-A(3,1)*A(1,2))/A(2,2)
      A(3,3) = A(3,3)-A(3,1)*A(1,3)-A(3,2)*A(2,3)
      A(3,4) = A(3,4)-A(3,1)*A(1,4)-A(3,2)*A(2,4)
      A(4,1) = A(4,1)/A(1,1)
      A(4,2) = (A(4,2)-A(4,1)*A(1,2))/A(2,2)
      A(4,3) = (A(4,3)-A(4,1)*A(1,3)-A(4,2)*A(2,3))/A(3,3)
      A(4,4) = A(4,4)-A(4,1)*A(1,4)-A(4,2)*A(2,4)-A(4,3)*A(3,4)
C---INVERT L
      A(4,3) = -A(4,3)
      A(4,2) = -A(4,2)-A(4,3)*A(3,2)
      A(4,1) = -A(4,1)-A(4,2)*A(2,1)-A(4,3)*A(3,1)
      A(3,2) = -A(3,2)
      A(3,1) = -A(3,1)-A(3,2)*A(2,1)
      A(2,1) = -A(2,1)
C---INVERT U
      A(4,4) = 1.D0/A(4,4)
      A(3,4) = -A(3,4)*A(4,4)/A(3,3)
      A(3,3) = 1.D0/A(3,3)
      A(2,4) = -(A(2,3)*A(3,4)+A(2,4)*A(4,4))/A(2,2)
      A(2,3) = -A(2,3)*A(3,3)/A(2,2)
      A(2,2) = 1.D0/A(2,2)
      A(1,4) = -(A(1,2)*A(2,4)+A(1,3)*A(3,4)+A(1,4)*A(4,4))/A(1,1)
      A(1,3) = -(A(1,2)*A(2,3)+A(1,3)*A(3,3))/A(1,1)
      A(1,2) = -A(1,2)*A(2,2)/A(1,1)
      A(1,1) = 1.D0/A(1,1)
C---MULTIPLY (U-INVERSE)*(L-INVERSE)
      A(1,1) = A(1,1)+A(1,2)*A(2,1)+A(1,3)*A(3,1)+A(1,4)*A(4,1)
      A(1,2) = A(1,2)+A(1,3)*A(3,2)+A(1,4)*A(4,2)
      A(1,3) = A(1,3)+A(1,4)*A(4,3)
      A(2,1) = A(2,2)*A(2,1)+A(2,3)*A(3,1)+A(2,4)*A(4,1)
      A(2,2) = A(2,2)+A(2,3)*A(3,2)+A(2,4)*A(4,2)
      A(2,3) = A(2,3)+A(2,4)*A(4,3)
      A(3,1) = A(3,3)*A(3,1)+A(3,4)*A(4,1)
      A(3,2) = A(3,3)*A(3,2)+A(3,4)*A(4,2)
      A(3,3) = A(3,3)+A(3,4)*A(4,3)
      A(4,1) = A(4,4)*A(4,1)
      A(4,2) = A(4,4)*A(4,2)
      A(4,3) = A(4,4)*A(4,3)
      RETURN
      END

      SUBROUTINE GAUSSP (N,XPT,XWT)
C-----------------------------------------------------------------------
C  Loads in pre-set Gauss points for 4 angles from 0 to +1 in cos(theta)=mu
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      REAL*8  XPT(N),XWT(N)
      REAL*8 GPT4(4),GWT4(4)
      DATA GPT4/.06943184420297D0,.33000947820757D0,.66999052179243D0,
     G          .93056815579703D0/
      DATA GWT4/.17392742256873D0,.32607257743127D0,.32607257743127D0,
     W          .17392742256873D0/
      N = 4
      DO I=1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      ENDDO
      RETURN
      END

      SUBROUTINE EXTRAL (JTAUMX,ATAU,ATAU0)
      use raqmschem_pmgrid_mod, only : nlev
      use raqmschemcomm_mod, only : np,ns,nw,mx,fastj,nb,lpar,jppj
C
      IMPLICIT NONE

!      include 'cmn_fj.h'
#include "jv_cmn.h"
c
      INTEGER NC
!      PARAMETER(NC=2*(LPAR+1))
c
      real*8 ttau(2*(nlev+1)+1)
!      REAL*8 DTAUX(NB),TTAU(NC+1),ATAU,ATAU0
      REAL*8 DTAUX(NB),ATAU,ATAU0
      real*8 ATAU1,ATAULN,ATAUM,ATAUN1
      INTEGER I,J,JTAUMX,JTOTL,NX,L
!      lpar=nlev
      nc=2*(lpar+1)
c
      NX=501
      DO J=1,NC+1
       TTAU(J) = 0.0D0
       ENDDO
C
      DO J=1,NB
       DTAUX(J) = 0.0D0
       DO I=1,MX
!        DTAUX(J) =DTAUX(J)+AER(I,J)
        DTAUX(J) =DTAUX(J)+fastj%AER(I,J)
        ENDDO
       ENDDO
C
CCCMNN081007   SEE V5.7 OF FAST_JX FOR COMMENTS ABOUT THE EVALUATION
CCC            OF EXTRA LEVELS.
      ATAU1 = ATAU - 1.0D0
      ATAULN = LOG(ATAU)
      TTAU(NC+1) = 0.0D0
      DO J=NC,1,-1
       I = (J+1)/2
       TTAU(J) = TTAU(J+1) + 0.5D0*DTAUX(I)
c---Now compute the number of log-spaced sub-layers to be added in

c---The objective is to have successive TAU-layers increasing by factor
cc     ATAU >1
c---the number of sub-layers + 1
       IF (TTAU(J) .LT. ATAU0) THEN
!        JADSUB(J) = 0
         fastj%JADSUB(J) = 0
        ELSE
          ATAUM = MAX(ATAU0,TTAU(J+1))
          ATAUN1 = LOG(TTAU(J)/ATAUM) / ATAULN
!         JADSUB(J) = MIN(JTAUMX,MAX(0,INT(ATAUN1 - 0.5D0)))
          fastj%JADSUB(J) = MIN(JTAUMX,MAX(0,INT(ATAUN1 - 0.5D0)))
        ENDIF
       ENDDO
c
C      CHECK ON OVERFLOW OF ARRAYS, CUT OFF JADSUB AT LOWER L IF THERE
C      ARE TOO MANY LEVELS
C
       JTOTL = NC+2
       DO J=NC,1,-1
!        JTOTL = JTOTL + JADSUB(J)
        JTOTL = JTOTL + fastj%JADSUB(J)
        IF (JTOTL .GT. NX/2) THEN
         WRITE (6,*) NX,NC,J
         DO L = J,1,-1
!          JADSUB(L) = 0
          fastj%JADSUB(L) = 0
          ENDDO
         GO TO 10
         ENDIF
        ENDDO
10     CONTINUE
       RETURN
       END
