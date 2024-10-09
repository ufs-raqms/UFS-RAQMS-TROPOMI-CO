#include <options.h>
module computeaod
  contains
      subroutine compute_aod_raqms(n1,n2,nchem,den,rh,dz, &
                    chemp,aod,nlay,ext_3d, &
                    mass_ext_3d,mole_frac_3d,aodpart)
!     ratioaod version
! ajl iac order
! 1 - s04
! 2 -bc1
! 3 -bc2
! 4 -oc1
! 5 -oc2
! 6- no3
! 7- du1 -4
! 11 ss1 -4 for ext_3d and mass_ext_3d, mole_frac_3d
! ajl aodpart ordef
! 1 - s04 2- bc1, bc2 4-oc1, oc2, 6 du1-4, 10 ss1-4, 14 no3
! new 1 - s04 2- bc1, bc2 4-oc1, oc2, 6-no3 7 du1-4, 11 ss1-4

! tks
! note subscripts on some arrays were changed.
!!!!!! arrays have different indices than in RAQMS
! rv is mixing ratio
! theta is theta here; in original code it was virtual theta

! This is based on 'compute_aod_RAQMS_ck.pro'

  use funcphys
  use raqmschem_pmgrid_mod, only : iam,beglat,endlat,masterproc,masterproct
!  use raqmschem_pmgrid_mod, only : aerosol_ugpkg
  use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
  use raqmschem_comm_mod, only : raqmschem_comm_all_bcast_tile
  use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
  use raqmschem_const_mod, only : mw_so4_aer

  implicit none
  integer n1,n2,nchem,nlay,localrc
! n1 is number of layers in atmosphere
! nlay is number of layers (from bottom) to do this calculation
  logical fexist
  real(CHEM_KIND_R4),optional :: aodpart(n2,beglat:endlat,14)
  real(CHEM_KIND_R4) :: aod(n2,beglat:endlat)
  real aod_so4(n2,beglat:endlat),aod_bc1(n2,beglat:endlat)
  real aod_bc2(n2,beglat:endlat),aod_oc1(n2,beglat:endlat)
  real aod_oc2(n2,beglat:endlat),aod_no3(n2,beglat:endlat)
  real aod_so4_hlut(n2,beglat:endlat),aod_bc1_hlut(n2,beglat:endlat)
  real aod_bc2_hlut(n2,beglat:endlat),aod_oc1_hlut(n2,beglat:endlat)
  real aod_oc2_hlut(n2,beglat:endlat),aod_no3_hlut(n2,beglat:endlat)
  real aod_du1_hlut(n2,beglat:endlat),aod_du2_hlut(n2,beglat:endlat)
  real aod_du3_hlut(n2,beglat:endlat),aod_du4_hlut(n2,beglat:endlat)
  real aod_ss1_hlut(n2,beglat:endlat),aod_ss2_hlut(n2,beglat:endlat)
  real aod_ss3_hlut(n2,beglat:endlat),aod_ss4_hlut(n2,beglat:endlat)
  real time,aod_cc(n2,beglat:endlat)
  real(CHEM_KIND_R4) ,dimension(n2,beglat:endlat,n1) :: den,dz
  real(CHEM_KIND_R4) ,dimension(n1,n2,beglat:endlat) :: rh
  real(CHEM_KIND_R4) ,dimension(n2,beglat:endlat,n1,nchem) :: chemp 


  real rmin,rms,rmax

  integer nac,iiac,iz,kk
!  parameter(nac=6)
  parameter(nac=5)
  real(chem_kind_r4), dimension(n1,nac) :: chemp_ppv

  real cp,R,cpr, cappa
!  parameter(cp=1004.,R=287.,cpr=cp/R, cappa= R/cp)
!match values in UW model (CCM3 physics)
  parameter(cp=1004.64,R=287.04,cpr=cp/R, cappa= R/cp)
  integer nest,iest,es2,i,j,k,idat
  parameter(nest=150)
  real(CHEM_KIND_R4) :: estbar_tmp,estbar(nest),es,es1,es3,es4,es5,qs


 integer nac_hlut
 parameter (nac_hlut=14)

  real mw_aero(nac_hlut),mw_air,mw_w
  parameter(mw_w=18.)
!tks
  real aero_min
!tks
  real den_aero(nac_hlut)
  real aero_conc(n1,nac_hlut)
  integer iac,itau

  character*100 infl,indir,date*10
  integer irh,nrh,irh1,irh2,irh1save(n2,nlay),irh2save(n2,nlay)
  parameter(nrh=99)
  integer rh_lut(nrh),rh_tmp
  real mole_frac_so4(nrh),den_mix_so4(nrh)
  real mole_frac_bc(nrh),den_mix_bc(nrh)
  real mole_frac_oc(nrh),den_mix_oc(nrh)
  real mole_frac_no3(nrh),den_mix_no3(nrh)
  real ratio_tmp,mole_frac_tmp,den_mix_tmp
  integer im,nmole_frac,im1,im2
  parameter(nmole_frac=100)
  real(CHEM_KIND_R4) :: moledensend(2,4,nrh)
  real(CHEM_KIND_R4) :: qrsend(2,4,nmole_frac)
  real qext_tmp,reff_tmp
  real(CHEM_KIND_R4) :: mole_frac_lut(nmole_frac)
  real qext_so4(nmole_frac)
  real qext_bc(nmole_frac)
  real qext_oc(nmole_frac)
  real qext_no3(nmole_frac)
  real reff_so4(nmole_frac)
  real reff_bc(nmole_frac)
  real reff_oc(nmole_frac)
  real reff_no3(nmole_frac)
  real molef_rh(n1,nac),vmxr_w(n1,nac)
  real mass_w(n1,nac),ext(n1,nac)
  real reff_molef(n1,nac),mole_frac(n1,nac)
  real den_mix_rh(n1,nac),mass_ext(n1,nac)
  real qext_molef(n1,nac)
!  real pi,dz_in
  real pi
!  real q_column(nac)
  real q_column(n2,nac)




!  parameter(pi=3.14159265,dz_in=400.)
  parameter(pi=3.14159265)

! Wavelength for which AOD is calculated
  integer n_atau
  parameter(n_atau=1)
  real a_wli
  parameter(a_wli=0.55)

! Output AOD array
  real a_taus(n_atau,nac)

  integer nc
  character*10 ct,ct_new
  character*200 fl_aod

  real vmax,vmin
  real deltax,deltay,centrallat,centrallon
  real(CHEM_KIND_R4) ext_3d(n2,beglat:endlat,n1,nac_hlut)
  real(CHEM_KIND_R4),optional :: mass_ext_3d(n2,beglat:endlat,n1,nac_hlut)
  
  real(CHEM_KIND_R4),optional :: mole_frac_3d(n2,beglat:endlat,n1,nac_hlut)

! Variables for Harvard Lookup table

 character(len=50)  infl_qext
 character(len=30)  spec_name(nac_hlut)
 integer ispec,ilamda,iterm
 integer nspec_hlut,nlamda_hlut,nterm_hlut
 parameter (nspec_hlut=42,nlamda_hlut=7,nterm_hlut=8)
! ajl speed up since ispecac doesn't change with time or i,j
 integer ispecac(nac_hlut),ilamdaac

 real rh_hlut(7)

 real(CHEM_KIND_R4) :: q_hlut(nlamda_hlut,nspec_hlut)
 real(CHEM_KIND_R4) :: reff_hlut(nlamda_hlut,nspec_hlut)
 real(CHEM_KIND_R4) :: ssa_hlut(nlamda_hlut,nspec_hlut)
 real(CHEM_KIND_R4) :: pf_hlut(nterm_hlut,nlamda_hlut,nspec_hlut)
 real(CHEM_KIND_R4) :: lamda_hlut(nlamda_hlut)
 real(CHEM_KIND_R4) :: tmp_hlut(12)

  real q_tmp_hlut(7),reff_tmp_hlut(7)
  real ssa_tmp_hlut(7),pf_tmp_hlut(nterm_hlut,7)
  real q_rh_hlut(n1,nac_hlut),reff_rh_hlut(n1,nac_hlut)
  real ssa_rh_hlut(n1,nac_hlut),pf_rh_hlut(n1,nterm_hlut,nac_hlut)
  real beta_aero_hlut(n1,nac_hlut),q_column_hlut(n2,nac_hlut)
  real f_w,den_aerow,mass_aerow,mass_aerow_min
  real Hterm,KEaod,Gterm

! Wavelength for which AOD is calculated
  integer n_atau_hlut
  parameter(n_atau_hlut=1)
  real a_wli_hlut
  parameter(a_wli_hlut=0.55)

! Output AOD array
  real a_taus_hlut(n_atau_hlut,nac_hlut)

 character(len=30)  spec_line(nspec_hlut)
 character(len=30)  line

 save rh_hlut
 data rh_hlut/0.,50.,70.,80.,90.,95.,99./
 logical first
 save first
 data first/.true./
!save and read once
 save mole_frac_so4,den_mix_so4,mole_frac_bc,den_mix_bc
 save mole_frac_no3,den_mix_no3,mole_frac_oc,den_mix_oc
 save estbar,qext_so4,reff_so4,qext_bc,reff_bc
 save qext_no3,reff_no3,qext_oc,reff_oc
 save lamda_hlut,q_hlut,reff_hlut,ssa_hlut,pf_hlut
 save ispecac,ilamdaac




!this is for GFS physics
!#ifdef USEGFUNCPHYS
!      call gfuncphys
!#endif


!tks 5-24-07 can be Nan otherwise
reff_rh_hlut=0.
q_rh_hlut=0


! aerosol physical/chemical properties
!20060223 den_aero(1)=1.7e+3   ! sulfate
den_aero(1)=1.769e+3   ! sulfate
den_aero(2)=1.0e+3   ! bc
den_aero(3)=1.0e+3   ! bc
den_aero(4)=1.4e+3   ! oc
den_aero(5)=1.4e+3   ! oc
den_aero(6)=1.725e+3   ! no3
den_aero(7)=2500.
den_aero(8)=2650.
den_aero(9)=2650.
den_aero(10)=2650.
den_aero(11)=2200.
den_aero(12)=2200.
den_aero(13)=2200.
den_aero(14)=2200.

!mw_aero(1)=96.   ! sulfate
mw_aero(1)=132.   ! sulfate assume ammonium sulfate here
mw_aero(2)=12.   ! bc
mw_aero(3)=12.   ! bc
mw_aero(4)=16.8   ! oc
mw_aero(5)=16.8   ! oc
mw_aero(6)=79.    ! no3
mw_aero(7)=28.97    !dust
mw_aero(8)=28.97    !dust
mw_aero(9)=28.97    !dust
mw_aero(10)=28.97   !dust
!Chieko said to use 29 for dust, now consistent with calling program
!mw_aero(7)=60.1
!mw_aero(8)=60.1
!mw_aero(9)=60.1
!mw_aero(10)=60.1
!J said to use 28.97 (air) for seasalt
mw_aero(11)=28.97   !seasalt
mw_aero(12)=28.97   !seasalt
mw_aero(13)=28.97   !seasalt
mw_aero(14)=28.97   !seasalt

spec_name(1)='sulfate'
spec_name(2)='Black C'
spec_name(3)='Black C'
spec_name(4)='Organic C'
spec_name(5)='Organic C'
spec_name(6)='sulfate' ! use sulfate extinction efficiency for now
spec_name(7)='Mdust 0.8'
spec_name(8)='Mdust 1.5'
spec_name(9)='Mdust 2.5'
spec_name(10)='Mdust 4.0'
spec_name(11)='SSa00 Sea Salt (accum)'
!spec_name(12)='SSa00 Sea Salt (accum)'
!changed 3-19-08
spec_name(12)='SSc00 Sea Salt (coarse)'
spec_name(13)='SSc00 Sea Salt (coarse)'
spec_name(14)='SSc00 Sea Salt (coarse)'

! Initialize mole_frac_3d as 1.
if(present(mole_frac_3d))then
 do iac=1,nac_hlut
! do k=1,n1

  do k=1,nlay
   do j= beglat,endlat
    do i=1,n2
      mole_frac_3d(i,j,k,iac)=1.
    enddo
  enddo
  enddo
  enddo
endif

! others
!tks 20080328mw_air=28.9
mw_air=28.97

if(first)then
  if(masterproc.or.masterproct)then
     if(masterproc)then
     write(6,*)'AOD DUST fix on 1/28/2011'
#ifndef NOBRADSSSFIX
     write(6,*)'BRAD SS fix on 2/2/2011 faster'
#endif
     endif

! Read in values from LUP

    infl='hygroscopic_growth_factors_nh4so4_0.1um_273K.txt'
    inquire(file=infl,exist=fexist)
    if(.not.fexist)then
      inquire(file=infl,exist=fexist)
      if(.not.fexist)then
         write(6,*)'INFL',infl,'does not exist'
         call killit('COMP')
      endif
    endif
    open (1,file=trim(infl),status='old', &
    form='formatted')

    do irh=1,nrh
      read (1,"(i12,f13.5,f13.7,f13.2)",advance="yes") &
      rh_tmp,ratio_tmp,mole_frac_tmp,den_mix_tmp
      rh_lut(irh)=rh_tmp
      mole_frac_so4(irh)=mole_frac_tmp
     den_mix_so4(irh)=den_mix_tmp
    enddo
    close(1)

    infl='hygroscopic_growth_factors_bc_0.1um_273K.txt'
    inquire(file=infl,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl),status='old', &
        form='formatted')
    else
      call killit('groth')
    endif
    do irh=1,nrh
      read (1,"(i12,f13.5,f13.7,f13.2)",advance="yes") &
      rh_tmp,ratio_tmp,mole_frac_tmp,den_mix_tmp
      mole_frac_bc(irh)=mole_frac_tmp
      den_mix_bc(irh)=den_mix_tmp
    enddo
    close(1)

    infl='hygroscopic_growth_factors_oc_0.1um_273K.txt'
    inquire(file=infl,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl),status='old', &
        form='formatted')
    else
        call killit('goc')
    open (1,file=trim(infl),status='old', &
        form='formatted')
    endif
    do irh=1,nrh
      read (1,"(i12,f13.5,f13.7,f13.2)",advance="yes") &
      rh_tmp,ratio_tmp,mole_frac_tmp,den_mix_tmp
      mole_frac_oc(irh)=mole_frac_tmp
      den_mix_oc(irh)=den_mix_tmp
    enddo
    close(1)

    infl='hygroscopic_growth_factors_nh4no3_0.1um_273K.txt'
    inquire(file=infl,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl),status='old', &
        form='formatted')
    else
       call killit('nh4')
    endif
    do irh=1,nrh
      read (1,"(i12,f13.5,f13.7,f13.2)",advance="yes") &
      rh_tmp,ratio_tmp,mole_frac_tmp,den_mix_tmp
      mole_frac_no3(irh)=mole_frac_tmp
      den_mix_no3(irh)=den_mix_tmp
    enddo
    close(1)
    infl='estbar.dat'
    inquire(file=infl,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl),status='old', &
        form='formatted')
    else
       call killit('est')
    endif
    estbar_tmp=0.
    do iest=1,nest
      read(1,*) estbar_tmp
      estbar(iest)=estbar_tmp
    enddo
    close(1)

! Read in values from LUT
    infl='scat_effic_nh4so4_r0.075_s2.txt'
    inquire(file=infl,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl),status='old', &
        form='formatted')
    else
       call killit('so4')
    endif
    do im=1,nmole_frac
      read (1,"(3f13.6)",advance="yes") &
      mole_frac_tmp,qext_tmp,reff_tmp
      mole_frac_lut(im)=mole_frac_tmp
      qext_so4(im)=qext_tmp
      reff_so4(im)=reff_tmp
    enddo
    close(1)

!   Read in values from LUT
    infl='scat_effic_bc_r0.075_s2.txt'
    inquire(file=infl,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl),status='old', &
        form='formatted')
    else
        call killit('bc')
    endif
    do im=1,nmole_frac
      read (1,"(3f13.6)",advance="yes") &
      mole_frac_tmp,qext_tmp,reff_tmp
      qext_bc(im)=qext_tmp
      reff_bc(im)=reff_tmp
    enddo
    close(1)

!   Read in values from LUT
    infl='scat_effic_oc_r0.075_s2.txt'
    inquire(file=infl,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl),status='old', &
        form='formatted')
    else
       call killit('oc')
    endif
    do im=1,nmole_frac
      read (1,"(3f13.6)",advance="yes") &
      mole_frac_tmp,qext_tmp,reff_tmp
      qext_oc(im)=qext_tmp
      reff_oc(im)=reff_tmp
    enddo
    close(1)

!   Read in values from LUT
   
    infl='scat_effic_nh4no3_r0.075_s2.txt'
    inquire(file=infl,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl),status='old', &
        form='formatted')
    else
        call killit('no3')
    endif
    do im=1,nmole_frac
      read (1,"(3f13.6)",advance="yes") &
     mole_frac_tmp,qext_tmp,reff_tmp
     qext_no3(im)=qext_tmp
     reff_no3(im)=reff_tmp
    enddo
    close(1)
  endif
! bcast parameters
  moledensend(1,1,:)=mole_frac_so4(:)
  moledensend(1,2,:)=mole_frac_bc(:)
  moledensend(1,3,:)=mole_frac_oc(:)
  moledensend(1,4,:)=mole_frac_no3(:)
  moledensend(2,1,:)=den_mix_so4(:)
  moledensend(2,2,:)=den_mix_bc(:)
  moledensend(2,3,:)=den_mix_oc(:)
  moledensend(2,4,:)=den_mix_no3(:)
  call raqmschem_comm_all_bcast_tile(moledensend,rc=localrc)
  call raqmschem_comm_all_bcast_tile(rh_lut,rc=localrc)
  mole_frac_so4(:)=moledensend(1,1,:)
  mole_frac_bc(:)=moledensend(1,2,:)
  mole_frac_oc(:)=moledensend(1,3,:)
  mole_frac_no3(:)=moledensend(1,4,:)
  den_mix_so4(:)=moledensend(2,1,:)
  den_mix_bc(:)=moledensend(2,2,:)
  den_mix_oc(:)=moledensend(2,3,:)
  den_mix_no3(:)=moledensend(2,4,:)
  call raqmschem_comm_all_bcast_tile(estbar,rc=localrc)
  call raqmschem_comm_all_bcast_tile(mole_frac_lut,rc=localrc)
  qrsend(1,1,:)=qext_so4(:)
  qrsend(1,2,:)=qext_bc(:)
  qrsend(1,3,:)=qext_oc(:)
  qrsend(1,4,:)=qext_no3(:)
  qrsend(2,1,:)=reff_so4(:)
  qrsend(2,2,:)=reff_bc(:)
  qrsend(2,3,:)=reff_oc(:)
  qrsend(2,4,:)=reff_no3(:)
  call raqmschem_comm_all_bcast_tile(qrsend,rc=localrc)
  qext_so4(:)=qrsend(1,1,:)
  qext_bc(:)=qrsend(1,2,:)
  qext_oc(:)=qrsend(1,3,:)
  qext_no3(:)=qrsend(1,4,:)
  reff_so4(:)=qrsend(2,1,:)
  reff_bc(:)=qrsend(2,2,:)
  reff_oc(:)=qrsend(2,3,:)
  reff_no3(:)=qrsend(2,4,:)

! um to m
  reff_so4=reff_so4*1.e-6
  reff_bc=reff_bc*1.e-6
  reff_oc=reff_oc*1.e-6
  reff_no3=reff_no3*1.e-6
! end of first
endif

aod=0.
aod_so4=0.
aod_bc1=0.
aod_bc2=0.
aod_oc1=0.
aod_oc2=0.
aod_no3=0.
aod_so4_hlut=0.
aod_bc1_hlut=0.
aod_bc2_hlut=0.
aod_oc1_hlut=0.
aod_oc2_hlut=0.
aod_no3_hlut=0.
aod_du1_hlut=0.
aod_du2_hlut=0.
aod_du3_hlut=0.
aod_du4_hlut=0.
aod_ss1_hlut=0.
aod_ss2_hlut=0.
aod_ss3_hlut=0.
aod_ss4_hlut=0.


!rh brought in now
do j=beglat,endlat
  do k=1,nlay
    do i=1,n2
      do irh=1,nrh
       if(rh(k,i,j).le.rh_lut(irh)) goto 18
      enddo
!     safer to say irh is nrh else could be nrh+1 which is out of bounds ajl
      irh=nrh
18    continue

      if(irh.eq.1) irh1=1
      if(irh.ge.2) irh1=irh-1
      if(irh.eq.1) irh2=1
      if(irh.ge.2) irh2=irh
      irh1save(i,k)=irh1
      irh2save(i,k)=irh2
    end do
  end do
!do i=1,n2

! in RAQMS theta is theta tks
! temp is T in K
! rv is specific humidity for UWCCM3

!  do k=1,n1
!  do k=1,nlay
!    p3d=pp(i,j,k)                                   ! tks pressure (hPa)
!    temp=theta(i,j,k)*(p3d/1000.)**cappa            ! tks T (K)

!tks    p3d=pp(k,i,j)                                   ! tks pressure (hPa)
!tks    temp=theta(k,i,j)*(p3d/1000.)**cappa            ! tks T (K)

! CCM3 Physics
!    ES1=ESTBLF(TEMP)*.01

!NCEP physics for now
!    ES1=fpvs(TEMP)*.01      !mb 
!    epsqs=0.622

!    ES1=MIN(ES1,p3d)
!    QS=EPSQS*ES1/(p3d - (1. - EPSQS)*ES1)
!tks    rh(k,i,j)=rv(k,i,j)/QS
! leaving rh subscipted this way since it is used later
! I did change rv since it comes from the model
!    rh(k,i,j)=100.*rv(i,j,k)/QS
!    if(rh(k,i,j).lt.0.0)print*,'neg rh ',i,j,k,rh(k,i,j)
!    if(rh(k,i,j).gt.100.0)print*,'>100 rh ',i,j,k,rh(k,i,j)

!    rh(k,i,j)=max(rh(k,i,j),0.)
!    rh(k,i,j)=min(rh(k,i,j),99.)

112  format(1x,i3,2x,3(f7.1,2x),e8.1,2x,f7.1)

!    es=0.
!    es1=min(temp-163.,148.)
!    es2=max(int(es1),0)
!    es3=min((temp-(es2+163.)*1.),0.)
!    es4=max(es3,0.)
!    es=estbar(es2)*(1.-es4)+estbar(es2+1)*es4
!    es5=0.
!    es5=max(1.e-10,(p3d*100.-es))
!    qs=.622*es/es5
!    rh(k,i,j)=rv(k,i,j)/qs*100.
!    rh(k,i,j)=max(rh(k,i,j),0.)
!    rh(k,i,j)=min(rh(k,i,j),99.)
!  enddo


!============================================================
do itau=1,n_atau
do iac=1,nac
!Chieko's
!if(iac.eq.1) iiac=4
!if(iac.eq.2) iiac=17
!if(iac.eq.3) iiac=18
!if(iac.eq.4) iiac=19
!if(iac.eq.5) iiac=20
!if(iac.eq.6) iiac=44

!RAQMS Global
 !if(iac.eq.1) iiac=1
 !if(iac.eq.2) iiac=2
 !if(iac.eq.3) iiac=3
 !if(iac.eq.4) iiac=4
 !if(iac.eq.5) iiac=5
 !if(iac.eq.6) iiac=6
! since equal just set
  iiac=iac
 q_column(:,iac)=0.

! Convert mixing ratio to mass concentration

!do k=1,n1
do k=1,nlay
!  kk=n1-k+1
 kk=k ! for fv3gfs
 do i=1,n2
!tks  chemp(k,i,j,1,iiac)=max(chemp(k,i,j,1,iiac),0.)
!tks aero_conc(k,iac)=chemp(k,i,j,1,iiac)*mw_aero(iac)/ &
!                 mw_air*den(k,i,j) ! in kg/m3
!tks
 chemp(i,j,k,iiac)=max(chemp(i,j,k,iiac),0.)
! if(aerosol_ugpkg)then
   if(iac.eq.1)then
!    for now tracer is so4 use its mw to convert to ppv
     aero_conc(k,iac)=chemp(i,j,k,iiac)*mw_air/mw_so4_aer*1.e-9
     aero_conc(k,iac)=aero_conc(k,iac)*mw_aero(iac)/mw_air*den(i,j,k)
     chemp_ppv(k,iac)=chemp(i,j,k,iiac)*mw_air/mw_so4_aer*1.e-9
   else
     aero_conc(k,iac)=chemp(i,j,k,iiac)*1.e-9*den(i,j,k) 
     chemp_ppv(k,iac)=chemp(i,j,k,iiac)*mw_air/mw_aero(iac)*1.e-9
   endif
! else
!   aero_conc(k,iac)=chemp(i,j,k,iiac)*mw_aero(iac)/ &
!                  mw_air*den(i,j,k) ! in kg/m3
! endif
!tks
! enddo


! Relative humidity (rh) dependency
!do k=1,n1

!do k=1,nlay
  ! Find two closest relative humidities from rh_lut
  
  if(iac.eq.2) then
     molef_rh(k,iac)=1.
     den_mix_rh(k,iac)=den_mix_bc(1)
  endif
  if(iac.eq.4) then
     molef_rh(k,iac)=1.
     den_mix_rh(k,iac)=den_mix_oc(1)
  endif
  irh1=irh1save(i,k)
  irh2=irh2save(i,k)

!  if((rh_lut(irh2)-rh_lut(irh1))<1.) then
  if((rh_lut(irh2)-rh_lut(irh1))<.9) then
    if(iac.eq.1) then
      molef_rh(k,iac)=mole_frac_so4(irh1)
      den_mix_rh(k,iac)=den_mix_so4(irh1)
    endif
    if(iac.eq.3) then
      molef_rh(k,iac)=mole_frac_bc(irh1)
      den_mix_rh(k,iac)=den_mix_bc(irh1)
    endif
    if(iac.eq.5) then
      molef_rh(k,iac)=mole_frac_oc(irh1)
      den_mix_rh(k,iac)=den_mix_oc(irh1)
    endif
    if(iac.eq.6) then
      molef_rh(k,iac)=mole_frac_no3(irh1)
     den_mix_rh(k,iac)=den_mix_no3(irh1)
    endif

  endif

!  if((rh_lut(irh2)-rh_lut(irh1))>=1.) then
  if((rh_lut(irh2)-rh_lut(irh1))>=.9) then
    if(iac.eq.1) then 
      molef_rh(k,iac)=mole_frac_so4(irh1)+ &
      (mole_frac_so4(irh2)-mole_frac_so4(irh1))/ &
         (rh_lut(irh2)-rh_lut(irh1))* &
         (rh(k,i,j)-rh_lut(irh1))
    endif
    if(iac.eq.3) then
      molef_rh(k,iac)=mole_frac_bc(irh1)+ &
      (mole_frac_bc(irh2)-mole_frac_bc(irh1))/ &
         (rh_lut(irh2)-rh_lut(irh1))* &
         (rh(k,i,j)-rh_lut(irh1))
    endif
    if(iac.eq.5) then
      molef_rh(k,iac)=mole_frac_oc(irh1)+ &
     (mole_frac_oc(irh2)-mole_frac_oc(irh1))/ &
         (rh_lut(irh2)-rh_lut(irh1))* &
         (rh(k,i,j)-rh_lut(irh1))
    endif
    if(iac.eq.6) then
      molef_rh(k,iac)=mole_frac_no3(irh1)+ &
      (mole_frac_no3(irh2)-mole_frac_no3(irh1))/ &
         (rh_lut(irh2)-rh_lut(irh1))* &
         (rh(k,i,j)-rh_lut(irh1))
    endif
    if(iac.eq.1) then
      den_mix_rh(k,iac)=den_mix_so4(irh1)+ &
            (den_mix_so4(irh2)-den_mix_so4(irh1))/ &
            (rh_lut(irh2)-rh_lut(irh1))* &
            (rh(k,i,j)-rh_lut(irh1))
    endif
    if(iac.eq.3) then
      den_mix_rh(k,iac)=den_mix_bc(irh1)+ &
            (den_mix_bc(irh2)-den_mix_bc(irh1))/ &
            (rh_lut(irh2)-rh_lut(irh1))* &
            (rh(k,i,j)-rh_lut(irh1))
    endif
    if(iac.eq.5) then
      den_mix_rh(k,iac)=den_mix_oc(irh1)+ &
            (den_mix_oc(irh2)-den_mix_oc(irh1))/ &
            (rh_lut(irh2)-rh_lut(irh1))* &
            (rh(k,i,j)-rh_lut(irh1))
    endif
    if(iac.eq.6) then
      den_mix_rh(k,iac)=den_mix_no3(irh1)+ &
            (den_mix_no3(irh2)-den_mix_no3(irh1))/ &
            (rh_lut(irh2)-rh_lut(irh1))* &
            (rh(k,i,j)-rh_lut(irh1))
    endif

  endif ! end if for rh

  den_mix_rh(k,iac)=min(den_mix_rh(k,iac),1800.)
  den_mix_rh(k,iac)=max(den_mix_rh(k,iac),1000.)

!  if(aerosol_ugpkg)then
    vmxr_w(k,iac)=(1.-molef_rh(k,iac))/molef_rh(k,iac)* &
                chemp_ppv(k,iiac)
!  else
!    vmxr_w(k,iac)=(1.-molef_rh(k,iac))/molef_rh(k,iac)* &
!                chemp(i,j,k,iiac)
!  endif
!tks                chemp(k,i,j,1,iiac)

!tks   mass_w(k,iac)=vmxr_w(k,iac)*mw_w/mw_air*den(k,i,j)

  mass_w(k,iac)=vmxr_w(k,iac)*mw_w/mw_air*den(i,j,k)
  Hterm=mw_w/mw_aero(iac)*(1.-molef_rh(k,iac))/molef_rh(k,iac)

  mole_frac(k,iac)=1.

!tks
! if((chemp(k,i,j,1,iiac)+vmxr_w(k,iac)).gt.1e-10) &
! mole_frac(k,iac)=chemp(k,i,j,1,iiac)/ &
!   (chemp(k,i,j,1,iiac)+vmxr_w(k,iac))
! ajl it turnout that since vmxr is prop to chemp that mole_frac = molef_rh
  mole_frac(k,iac)=molef_rh(k,iac)
!  if(doratioaod)then
!    if(chemp(i,j,k,iiac).ne.0.0)then
!      mole_frac(k,iac)=chemp(i,j,k,iiac)/ &
!      (chemp(i,j,k,iiac)+vmxr_w(k,iac))
!      if(abs(molef_rh(k,iac)-mole_frac(k,iac))>1.e-10)then
!        print *,'molef error ',i,j,k,iac,molef_rh(k,iac),mole_frac(k,iac)
!      endif
!    endif
!  else
!
!    if((chemp(i,j,k,iiac)+vmxr_w(k,iac)).gt.1e-10) &
!    mole_frac(k,iac)=chemp(i,j,k,iiac)/ &
!    (chemp(i,j,k,iiac)+vmxr_w(k,iac))
!  endif



! Find two closest mole fractions from lut
  do im=1,nmole_frac
   if(mole_frac(k,iac).le.mole_frac_lut(im)) goto 19
  enddo
  
19 continue
  if(im.eq.1) im1=1
  if(im.ge.2) im1=im-1
  if(im.eq.1) im2=1
  if(im.ge.2) im2=im

  if(iac.eq.2) then
   qext_molef(k,iac)=qext_bc(nmole_frac)
   reff_molef(k,iac)=reff_bc(nmole_frac)
  endif
  if(iac.eq.4) then
   qext_molef(k,iac)=qext_oc(nmole_frac)
   reff_molef(k,iac)=reff_oc(nmole_frac)
  endif
  
!  if((mole_frac_lut(im2)-mole_frac_lut(im1))<0.01) then
  if((mole_frac_lut(im2)-mole_frac_lut(im1))<0.005) then
    if(iac.eq.1) then
      qext_molef(k,iac)=qext_so4(im1)
      reff_molef(k,iac)=reff_so4(im1)
    endif
    if(iac.eq.3) then
      qext_molef(k,iac)=qext_bc(im1)
      reff_molef(k,iac)=reff_bc(im1)
    endif
    if(iac.eq.5) then
      qext_molef(k,iac)=qext_oc(im1)
      reff_molef(k,iac)=reff_oc(im1)
    endif
    if(iac.eq.6) then
      qext_molef(k,iac)=qext_no3(im1)
      reff_molef(k,iac)=reff_no3(im1)
    endif
  endif


  if((mole_frac_lut(im2)-mole_frac_lut(im1))>=0.005) then

    if(iac.eq.1) then
      qext_molef(k,iac)=qext_so4(im1)+ &
            (qext_so4(im2)-qext_so4(im1))/ &
            (mole_frac_lut(im2)-mole_frac_lut(im1))* &
            (mole_frac(k,iac)-mole_frac_lut(im1))

     reff_molef(k,iac)=reff_so4(im1)+ &
            (reff_so4(im2)-reff_so4(im1))/ &
            (mole_frac_lut(im2)-mole_frac_lut(im1))* &
            (mole_frac(k,iac)-mole_frac_lut(im1))
    endif

    if(iac.eq.3) then
      qext_molef(k,iac)=qext_bc(im1)+ &
            (qext_bc(im2)-qext_bc(im1))/ &
            (mole_frac_lut(im2)-mole_frac_lut(im1))* &
            (mole_frac(k,iac)-mole_frac_lut(im1))
      reff_molef(k,iac)=reff_bc(im1)+ &
            (reff_bc(im2)-reff_bc(im1))/ &
            (mole_frac_lut(im2)-mole_frac_lut(im1))* &
            (mole_frac(k,iac)-mole_frac_lut(im1))
    endif

    if(iac.eq.5) then
      qext_molef(k,iac)=qext_oc(im1)+ &
            (qext_oc(im2)-qext_oc(im1))/ &
            (mole_frac_lut(im2)-mole_frac_lut(im1))* &
            (mole_frac(k,iac)-mole_frac_lut(im1))
      reff_molef(k,iac)=reff_oc(im1)+ &
            (reff_oc(im2)-reff_oc(im1))/ &
            (mole_frac_lut(im2)-mole_frac_lut(im1))* &
            (mole_frac(k,iac)-mole_frac_lut(im1))
    endif

    if(iac.eq.6) then
      qext_molef(k,iac)=qext_no3(im1)+ &
            (qext_no3(im2)-qext_no3(im1))/ &
            (mole_frac_lut(im2)-mole_frac_lut(im1))* &
            (mole_frac(k,iac)-mole_frac_lut(im1))
     reff_molef(k,iac)=reff_no3(im1)+ &
            (reff_no3(im2)-reff_no3(im1))/ &
            (mole_frac_lut(im2)-mole_frac_lut(im1))* &
            (mole_frac(k,iac)-mole_frac_lut(im1))
    endif

  endif ! end of mole_frac


  mass_ext(k,iac)=3./4.*qext_molef(k,iac)/ &
                  (den_mix_rh(k,iac)*reff_molef(k,iac))
  KEaod=mass_ext(k,iac)*(1.+Hterm)

  ext(k,iac)=mass_ext(k,iac)*(aero_conc(k,iac)+mass_w(k,iac))
!  q_column(iac)=q_column(iac)+ext(k,iac)*dz_in
  q_column(i,iac)=q_column(i,iac)+ext(k,iac)*dz(i,j,k)
 
  ext_3d(i,j,kk,iac)=ext(k,iac)
  if(present(mass_ext_3d))then
    mass_ext_3d(i,j,k,iac)=mass_ext(k,iac)
  endif
  if(present(mole_frac_3d))then
    mole_frac_3d(i,j,k,iac)=mole_frac(k,iac)
  endif
 enddo ! new i

enddo   ! end loop over layers
  do i=1,n2 

    a_taus(itau,iac)=q_column(i,iac)


!enddo  ! end loop over species
!enddo  ! end loop over wavelength

!do itau=1,n_atau
! do iac=1,nac
    aod(i,j)=aod(i,j)+a_taus(itau,iac)
    if(iac.eq.1) aod_so4(i,j)=aod_so4(i,j)+a_taus(itau,iac)
    if(iac.eq.2) aod_bc1(i,j)=aod_bc1(i,j)+a_taus(itau,iac)
    if(iac.eq.3) aod_bc2(i,j)=aod_bc2(i,j)+a_taus(itau,iac)
    if(iac.eq.4) aod_oc1(i,j)=aod_oc1(i,j)+a_taus(itau,iac)
    if(iac.eq.5) aod_oc2(i,j)=aod_oc2(i,j)+a_taus(itau,iac)
    if(iac.eq.6) aod_no3(i,j)=aod_no3(i,j)+a_taus(itau,iac)
  enddo ! i
 enddo ! iac
enddo ! itau
do i=1,n2 

  aod_cc(i,j)=aod_bc1(i,j)+aod_bc2(i,j)+aod_oc1(i,j)+aod_oc2(i,j)


enddo  ! end loop over i
enddo  ! end loop over j



!tks jump around dust contributions
!!!!!!!!!!!!!!!!!!!!!      go to 1147


! Write AOD into a file

!      write(ct,'(I10)')INT(TIME+.1)
!      DO NC=1,10
!       if(ct(nc:nc).ne.' ') then
!        ct_new=ct(nc:10)
!        goto 999
!       endif
!      ENDDO
!999    continue

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
if (first)then
! Compute AOD using Harvard Lookup table

! Read in values from LUP

!tks
  if(masterproc.or.masterproct)then
    infl_qext='Qext_Reff_GEOS-CHEM.dat'
    inquire(file=infl_qext,exist=fexist)
    if(fexist)then
      open (1,file=trim(infl_qext),status='old', &
        form='formatted')
    else
      call killit('qext')

    endif
    read (1,"(a20)",advance="yes") line
    do ispec=1,nspec_hlut
      read (1,"(a20)",advance="yes") line
      spec_line(ispec)=line
      do ilamda=1,nlamda_hlut
        read (1,*) tmp_hlut
        lamda_hlut(ilamda)=tmp_hlut(1)
        q_hlut(ilamda,ispec)=tmp_hlut(2)
        reff_hlut(ilamda,ispec)=tmp_hlut(3)*1.e-6
        ssa_hlut(ilamda,ispec)=tmp_hlut(4)
        pf_hlut(1,ilamda,ispec)=tmp_hlut(5)
        pf_hlut(2,ilamda,ispec)=tmp_hlut(6)
        pf_hlut(3,ilamda,ispec)=tmp_hlut(7)
        pf_hlut(4,ilamda,ispec)=tmp_hlut(8)
        pf_hlut(5,ilamda,ispec)=tmp_hlut(9)
        pf_hlut(6,ilamda,ispec)=tmp_hlut(10)
        pf_hlut(7,ilamda,ispec)=tmp_hlut(11)
        pf_hlut(8,ilamda,ispec)=tmp_hlut(12)
     enddo
    enddo
    close(1)
  endif
  call raqmschem_comm_all_bcast_tile(spec_line,rc=localrc)
  call raqmschem_comm_all_bcast_tile(pf_hlut,rc=localrc)
  call raqmschem_comm_all_bcast_tile(lamda_hlut,rc=localrc)
  call raqmschem_comm_all_bcast_tile(q_hlut,rc=localrc)
  call raqmschem_comm_all_bcast_tile(reff_hlut,rc=localrc)
  call raqmschem_comm_all_bcast_tile(ssa_hlut,rc=localrc)
!  moved here since doesn't change
   do ilamda=1,nlamda_hlut
     if(abs(lamda_hlut(ilamda)-a_wli_hlut*1000.).lt.1.) goto 9
   enddo
   write(6,*)'lamda is not found'
   call flush(6)
   call killit('no lamda')
9 continue
  ilamdaac=ilamda
  if(masterproc)print *,'ilamdaac=',ilamdaac
  do iac=1,nac_hlut
    do ispec=1,nspec_hlut
     if(index(spec_line(ispec),trim(spec_name(iac)(1:17))).gt.0) goto 99
    enddo
    do ispec=1,nspec_hlut
      write(100+iam,*)'spec is not found',ispec,spec_line(ispec)
      call flush(100+iam)
    end do
     call killit('spec not found')
99 continue
  ispecac(iac)=ispec
end do ! iac
! end of first
  first=.false.
endif

ilamda=ilamdaac
do j=beglat,endlat 
  do k=1,nlay
    do i=1,n2
      do irh=1,7
        if(rh(k,i,j).le.rh_hlut(irh)) goto 8
      enddo 
!     ajl safer to set to 7 else is 8 out of bounds
      irh=7
8     continue

      if(irh.eq.1) irh1=1
      if(irh.ge.2) irh1=irh-1
      if(irh.eq.1) irh2=1
      if(irh.ge.2) irh2=irh
      irh1save(i,k)=irh1
      irh2save(i,k)=irh2
    end do
  end do
!do i=1,n2

do itau=1,n_atau_hlut  ! loop ever wavelengths



do iac=1,nac_hlut
  if(iac.eq.6)cycle ! no no3aer
!Chieko's
!if(iac.eq.1) iiac=4
!if(iac.eq.2) iiac=17
!if(iac.eq.3) iiac=18
!if(iac.eq.4) iiac=19
!if(iac.eq.5) iiac=20
!if(iac.eq.6) iiac=44
!if(iac.eq.7) iiac=21
!if(iac.eq.8) iiac=22
!if(iac.eq.9) iiac=23
!if(iac.eq.10) iiac=24
!if(iac.eq.11) iiac=26
!if(iac.eq.12) iiac=27
!if(iac.eq.13) iiac=28
!if(iac.eq.14) iiac=29

!RAQMS Global (index of chemp)
! since equal just set don't test

! if(iac.eq.1) iiac=1
! if(iac.eq.2) iiac=2
! if(iac.eq.3) iiac=3
! if(iac.eq.4) iiac=4
! if(iac.eq.5) iiac=5
! if(iac.eq.6) iiac=6
! if(iac.eq.7) iiac=7
! if(iac.eq.8) iiac=8
! if(iac.eq.9) iiac=9
! if(iac.eq.10) iiac=10
! if(iac.eq.11) iiac=11
! if(iac.eq.12) iiac=12
! if(iac.eq.13) iiac=13
! if(iac.eq.14) iiac=14
 q_column_hlut(:,iac)=0.
! ajl new faster
 iiac=iac
 ispec=ispecac(iac)

!============================================================



! Convert mixing ratio to mass concentration
! do k=1,n1
! do k=1,nlay
!tks aero_conc(k,iac)=chemp(k,i,j,1,iiac)*mw_aero(iac)/ &
!tks                  mw_air*den(k,i,j) ! in kg/m3

! aero_conc(k,iac)=chemp(i,j,k,iiac)*mw_aero(iac)/ &
!                  mw_air*den(i,j,k) ! in kg/m3
! enddo

!-------------------------------------------!
!  q(nlamda,nspec),  reff(nlamda,nspec)
!  ssa(nlamda,nspec),pf(nterm,nlamda,nspec)
!-------------------------------------------!


! Exitinction coefficient at lamda

! ispec:ispec+6 code was put in place of original
! ispec:ispec+7 from Chieko.  Chieko and Brad 
! reviewed this code and recommended the change 4-19-07

if(iac.le.6.or.iac.ge.11) & 
 q_tmp_hlut(1:7)=q_hlut(ilamda,ispec:ispec+6)
if(iac.ge.7.and.iac.le.10) &
 q_tmp_hlut(1)=q_hlut(ilamda,ispec)

! Effective radius
if(iac.le.6.or.iac.ge.11) &
 reff_tmp_hlut(1:7)=reff_hlut(ilamda,ispec:ispec+6)
if(iac.ge.7.and.iac.le.10) &
 reff_tmp_hlut(1)=reff_hlut(ilamda,ispec)

! Single scattering albedo
if(iac.le.6.or.iac.ge.11) &
 ssa_tmp_hlut(1:7)=ssa_hlut(ilamda,ispec:ispec+6)
if(iac.ge.7.and.iac.le.10) &
 ssa_tmp_hlut(1)=ssa_hlut(ilamda,ispec)

! Phase function terms
if(iac.le.6.or.iac.ge.11) &
 pf_tmp_hlut(1:nterm_hlut,1:7)=pf_hlut(1:nterm_hlut,ilamda,ispec:ispec+6)

if(iac.ge.7.and.iac.le.10) &
 pf_tmp_hlut(1:nterm_hlut,1)=pf_hlut(1:nterm_hlut,ilamda,ispec)

! Relative humidity (rh) dependency
do k=1,nlay
! kk=n1-k+1
 kk=k ! for fv3gfs
 do i=1,n2

  ! Find two closest relative humidities from rh_lut
!  if(aerosol_ugpkg)then
   if(iac.eq.1)then
!    for now tracer is so4 use its mw to convert to ppv
     aero_conc(k,iac)=chemp(i,j,k,iiac)*mw_air/mw_so4_aer*1.e-9
     aero_conc(k,iac)=aero_conc(k,iac)*mw_aero(iac)/mw_air*den(i,j,k)
   else
     aero_conc(k,iac)=chemp(i,j,k,iiac)*1.e-9*den(i,j,k) 
   endif
!  else
!    aero_conc(k,iac)=chemp(i,j,k,iiac)*mw_aero(iac)/ &
!                  mw_air*den(i,j,k) ! in kg/m3
!  endif
  irh1=irh1save(i,k)
  irh2=irh2save(i,k)
  if((rh_hlut(irh2)-rh_hlut(irh1))>=1.) &
  q_rh_hlut(k,iac)=q_tmp_hlut(irh1)+ &
         (q_tmp_hlut(irh2)-q_tmp_hlut(irh1))/ &
         (rh_hlut(irh2)-rh_hlut(irh1))* &
         (rh(k,i,j)-rh_hlut(irh1))
  if(iac.eq.2.or.iac.eq.4) q_rh_hlut(k,iac)=q_tmp_hlut(1)

  if(iac.ge.7.and.iac.le.10) q_rh_hlut(k,iac)=q_tmp_hlut(1)


! if((rh_hlut(irh2)-rh_hlut(irh1))>=1.) &
! reff_rh_hlut(k,iac)=reff_tmp_hlut(irh1)+ &
!           (reff_tmp_hlut(irh2)-reff_tmp_hlut(irh1))/ &
!           (rh_hlut(irh2)-rh_hlut(irh1))* &
!           (rh(k,i,j)-rh_hlut(irh1))

  if((rh_hlut(irh2)-rh_hlut(irh1))>=1.) then
    reff_rh_hlut(k,iac)=reff_tmp_hlut(irh1)+ &
            (reff_tmp_hlut(irh2)-reff_tmp_hlut(irh1))/ &
            (rh_hlut(irh2)-rh_hlut(irh1))* &
            (rh(k,i,j)-rh_hlut(irh1))
!
! rbp fix hydrophillic bc, oc and dust reff bug (no re-initialize)
!
!   may not want inside of this loop
    if(iac.eq.2.or.iac.eq.4) reff_rh_hlut(k,iac)=reff_tmp_hlut(1)

    if(iac.ge.7.and.iac.le.10) reff_rh_hlut(k,iac)=reff_tmp_hlut(1)
! else ajl
  else
   write(6,*)'i',i,j,k,'at zero diff irh1',irh1,irh2,'reff',reff_rh_hlut(k,iac),' rh ',rh(k,i,j)
   call flush(6)

  endif

  reff_rh_hlut(k,iac)=max(reff_rh_hlut(k,iac),1.e-10)
  reff_rh_hlut(k,iac)=min(reff_rh_hlut(k,iac),10.e-6)
!
! rbp fix hydrophillic bc, oc and dust ssa bug (wrong place)
!
!  if(iac.eq.2.or.iac.eq.4) ssa_rh_hlut(k,iac)=ssa_tmp_hlut(1)
!
!  if(iac.ge.7.and.iac.le.10) ssa_rh_hlut(k,iac)=ssa_tmp_hlut(1)


  if((rh_hlut(irh2)-rh_hlut(irh1))>=1.) &
  ssa_rh_hlut(k,iac)=ssa_tmp_hlut(irh1)+ &
           (ssa_tmp_hlut(irh2)-ssa_tmp_hlut(irh1))/ &
           (rh_hlut(irh2)-rh_hlut(irh1))* &
           (rh(k,i,j)-rh_hlut(irh1))
  if((rh_hlut(irh2)-rh_hlut(irh1))<1.) then
    write(6,*)'rhdiff too small',i,j,k,iac,' irh1 ',irh1,irh2
    call flush(6)
    call killit('rhdiff')
  endif
!
! rbp fix hydrophillic bc, oc and dust ssa bug (no re-initialize)
!

  if(iac.eq.2.or.iac.eq.4) ssa_rh_hlut(k,iac)=ssa_tmp_hlut(1)

  if(iac.ge.7.and.iac.le.10) ssa_rh_hlut(k,iac)=ssa_tmp_hlut(1)



  do iterm=1,nterm_hlut
   if((rh_hlut(irh2)-rh_hlut(irh1))>=1.) &
   pf_rh_hlut(k,iterm,iac)=pf_tmp_hlut(iterm,irh1)+ &
           (pf_tmp_hlut(iterm,irh2)-pf_tmp_hlut(iterm,irh1))/ &
            (rh_hlut(irh2)-rh_hlut(irh1))* &
            (rh(k,i,j)-rh_hlut(irh1))
   if(iac.eq.2.or.iac.eq.4) &
    pf_rh_hlut(k,iterm,iac)=pf_tmp_hlut(iterm,1)
   if(iac.ge.7.and.iac.le.10) &
    pf_rh_hlut(k,iterm,iac)=pf_tmp_hlut(iterm,1)
  enddo

!tks    reff_rh_hlut can't be less that 1e-10 because of 
!of the max(), min() 20 lines above

  if(reff_rh_hlut(k,iac)>=1.e-15) &
   f_w=reff_tmp_hlut(1)**3/reff_rh_hlut(k,iac)**3
  if(reff_rh_hlut(k,iac)<1.e-15) &
   f_w=1.
  den_aerow=f_w*den_aero(iac)+(1.-f_w)*1.e+3

!  mass_aerow=aero_conc(k,iac)
!   +4./3.*pi*(reff_rh(k,iac)**3)*den_aerow

  mass_aerow=1./f_w*den_aerow/den_aero(iac)*aero_conc(k,iac)

  ! beta: specific extinction (m2/kg)
!tks 4-24-07

  !aero_min=1.e-20
  aero_min=1.e-40

!    aero_min=1.e-12

!tks 3-27-08 added after discussion with Brad, avoids possible division
! by zero

  mass_aerow_min=1./f_w*den_aerow/den_aero(iac)*aero_min
  if(mass_aerow.lt.mass_aerow_min)mass_aerow=mass_aerow_min
! tks 3-27-08


!tks  if(aero_conc(k,iac)>=1.e-12) &
!tks 3-27-08  if(aero_conc(k,iac)>=aero_min) &

  beta_aero_hlut(k,iac)=3.*q_rh_hlut(k,iac) &
   /(4.*den_aerow*reff_rh_hlut(k,iac))


!tks  if(aero_conc(k,iac) < 1.e-12) beta_aero_hlut(k,iac)=0.
!tks 3-27-08  if(aero_conc(k,iac) < aero_min) beta_aero_hlut(k,iac)=0.
#ifndef NOBRADSSSFIX
  q_column_hlut(i,iac)=q_column_hlut(i,iac)+ &
    beta_aero_hlut(k,iac)*(aero_conc(k,iac)+mass_aerow)*dz(i,j,k)
  if(iac.ge.7.and.iac.le.14) then 
    ext_3d(i,j,kk,iac)=beta_aero_hlut(k,iac)*(aero_conc(k,iac)+mass_aerow)
    if(present(mass_ext_3d))then
      mass_ext_3d(i,j,k,iac)=beta_aero_hlut(k,iac)
     endif
  endif
#else

  q_column_hlut(i,iac)=q_column_hlut(i,iac)+ &
!   beta_aero_hlut(k,iac)*mass_aerow*dz_in
   beta_aero_hlut(k,iac)*mass_aerow*dz(i,j,k)
!  if(i.eq.1.and.j.eq.91)then
!     write(6,*)'q_column',k,q_column_hlut(i,iac),'beta',beta_aero_hlut(k,iac),'mass',mass_aerow,'dz',dz(i,j,k)
  !endif
!
! rbp update (add aerosol and water mass, will effect ssalt)
!
! to be tested   beta_aero_hlut(k,iac)*(aero_conc(k,iac)+mass_aerow)*dz(i,j,k)


  if(iac.ge.7.and.iac.le.14) then 
  ext_3d(i,j,k,iac)=beta_aero_hlut(k,iac)*aero_conc(k,iac)
!
! rbp update (add aerosol and water mass, will effect ssalt)
!
! to be tested   ext_3d(i,j,k,iac)=beta_aero_hlut(k,iac)*(aero_conc(k,iac)+mass_aerow)

  if(present(mass_ext_3d))then
    mass_ext_3d(i,j,k,iac)=beta_aero_hlut(k,iac)
  endif

  endif
#endif
 enddo ! i new i

enddo   ! end loop over layers
do i=1,n2 ! new i

 a_taus_hlut(itau,iac)=q_column_hlut(i,iac)


!enddo  ! end loop over species
!enddo  ! end loop over wavelengths

!do itau=1,n_atau_hlut
!  do iac=1,nac_hlut
!20880327  if(iac.ge.7.and.iac.le.10) &
  if(iac.ge.7.and.iac.le.14) &
   aod(i,j)=aod(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.1) aod_so4_hlut(i,j)=aod_so4_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.2) aod_bc1_hlut(i,j)=aod_bc1_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.3) aod_bc2_hlut(i,j)=aod_bc2_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.4) aod_oc1_hlut(i,j)=aod_oc1_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.5) aod_oc2_hlut(i,j)=aod_oc2_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.6) aod_no3_hlut(i,j)=aod_no3_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.7) aod_du1_hlut(i,j)=aod_du1_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.8) aod_du2_hlut(i,j)=aod_du2_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.9) aod_du3_hlut(i,j)=aod_du3_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.10) aod_du4_hlut(i,j)=aod_du4_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.11) aod_ss1_hlut(i,j)=aod_ss1_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.12) aod_ss2_hlut(i,j)=aod_ss2_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.13) aod_ss3_hlut(i,j)=aod_ss3_hlut(i,j)+a_taus_hlut(itau,iac)
  if(iac.eq.14) aod_ss4_hlut(i,j)=aod_ss4_hlut(i,j)+a_taus_hlut(itau,iac)
 end do ! new i
 enddo
enddo

!enddo ! end loop over i
enddo ! end loop over j


!!!!!!!!!!!!!!!!!!!!!!!!1147  continue

!============================================================
!!!tks change
!!!fl_aod='/idea_data1/Aero_fx/Output/AOD/uwnms_aod_'// &
!!!trim(ct_new)//'.dat'


 if(present(aodpart))then
  do j= beglat,endlat
   do i=1,n2
    aodpart(i,j,1)=aod_so4(i,j)
    aodpart(i,j,2)=aod_bc1(i,j)
    aodpart(i,j,3)=aod_bc2(i,j)
    aodpart(i,j,4)=aod_oc1(i,j)
    aodpart(i,j,5)=aod_oc2(i,j)
    aodpart(i,j,6)=aod_no3(i,j)
    aodpart(i,j,7)=aod_du1_hlut(i,j)
    aodpart(i,j,8)=aod_du2_hlut(i,j)
    aodpart(i,j,9)=aod_du3_hlut(i,j)
    aodpart(i,j,10)=aod_du4_hlut(i,j)
    aodpart(i,j,11)=aod_ss1_hlut(i,j)
    aodpart(i,j,12)=aod_ss2_hlut(i,j)
    aodpart(i,j,13)=aod_ss3_hlut(i,j)
    aodpart(i,j,14)=aod_ss4_hlut(i,j)
   end do
  end do
! call mxmn2(aodpart(1,beglat,10),'ss1',0,1,rmax,rmin,rms)
 endif
#ifdef SAVEAODDAT
!write(1) aod_so4
!write(1) aod_bc1
!write(1) aod_bc2
!write(1) aod_oc1
!write(1) aod_oc2
!write(1) aod_no3
!write(1) aod_du1_hlut
!write(1) aod_du2_hlut
!write(1) aod_du3_hlut
!write(1) aod_du4_hlut
!write(1) aod_ss1_hlut
!write(1) aod_ss2_hlut
!write(1) aod_ss3_hlut
!write(1) aod_ss4_hlut
!write(1) aod_so4_hlut
!write(1) aod_bc1_hlut
!write(1) aod_bc2_hlut
!write(1) aod_oc1_hlut
!write(1) aod_oc2_hlut
!write(1) aod_no3_hlut
!write(1) aod_cc
!write(1) chemp
 if(masterproc.or.masterproct)then
 close(1)
 endif
#endif

return
end subroutine compute_aod_raqms 
end module computeaod
