#include <choosevan.h>
#include <options.h>
module raqmschemcomm_mod
  use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
  implicit none
  logical,parameter :: use_fv3_aero=.true.
  integer julday
  integer, allocatable :: i1grd(:,:),j1grd(:,:),ltb(:)
! for wetdep
  real(CHEM_KIND_R8) :: dy_jdgmt
#ifdef WETDEPDIAG
  real(CHEM_KIND_R4), allocatable, dimension(:,:,:,:) :: wetc,wetl,wetcf,wetlf
  real(CHEM_KIND_R4), allocatable, dimension(:,:) :: lpaveij,cpaveij
  real(CHEM_KIND_R4), allocatable, dimension(:,:,:) :: wetccol,wetlcol
  real(CHEM_KIND_R4), allocatable, dimension(:,:,:) :: gsdcol
  real(CHEM_KIND_R4), allocatable, dimension(:,:,:,:) :: gsd3d
  integer nwetdepdiag
#endif
  real(CHEM_KIND_R4), allocatable, dimension(:,:) :: cprate,lprate
  real(CHEM_KIND_R4), allocatable, dimension(:,:) :: ptropraq,ztropraq
  real(CHEM_KIND_R4), allocatable, dimension(:,:) :: ptropgsi,ztropgsi
  integer,allocatable,dimension(:,:) :: ktrop
  integer ncplprate
! ajl add angles for vort calc
!  real(CHEM_KIND_R4),allocatable, dimension(:,:) :: angledx,angledy
  real(CHEM_KIND_R4),allocatable, dimension(:,:) :: usina,ucosa,sinb,cosb
  real(CHEM_KIND_R4),allocatable, dimension(:,:) :: rcav_save,rnav_save ! these are cumulative so keep
  real(CHEM_KIND_R4),allocatable, dimension(:,:) :: rcav_write,rnav_write ! these are cumulative so keep

  real, allocatable :: dist1grd(:,:)
  real(CHEM_KIND_R8) :: gmt
  real*4, allocatable :: dpo3(:,:,:),dphno3(:,:,:),dph2o2(:,:,:),dpno2(:,:,:)
  real*4, allocatable :: dpco(:,:,:)
  real*4, allocatable :: srcnbio(:,:,:),srcnsoil(:,:,:)
  real*4, allocatable :: srcnair(:,:,:),slghtng(:,:,:),srcco_totl(:,:,:)
  real*4, allocatable :: sourceno(:,:)
  real, allocatable :: pksrcnair(:)
! ajl new aircraft sources
  real*4, allocatable :: srccoair(:,:,:)
#ifdef AERCO
  real*4, allocatable :: srcbcair(:,:,:),
  real*4, allocatable :: srcocair(:,:,:)
#endif
  real*4, allocatable :: srcso2air(:,:,:)
  real*4, allocatable :: bbefac(:,:,:),scale_factor(:,:)
  real*4, allocatable :: srcisop(:,:,:),srcterp(:,:,:)
  real*4, allocatable :: ch4clim(:,:,:)
  real*4, allocatable :: szagrd(:,:),spgrd(:,:),albgrd(:,:),zlwigrd(:,:)
  real*4, allocatable :: latgrd(:,:),longrd(:,:)
  real*4, allocatable :: tgrd(:,:,:),pgrd(:,:,:),zgrd(:,:,:),zeupgrd(:,:,:)
!  real*4, allocatable :: dpmgrd(:,:,:),szgrd(:,:)
  real*4, allocatable :: dpmgrd(:,:,:)
  real*4, allocatable :: qgrd(:,:,:),thgrd(:,:,:),ugrd(:,:,:),vgrd(:,:,:)
  real*4, allocatable :: rcraingrd(:,:,:),zsurf(:,:),slp(:,:),thbl(:,:)
  real*4, allocatable :: diurngrd(:,:),tskin(:,:),sfcwind(:,:)
#ifdef AERO
  real*4, allocatable :: seaicethk(:,:)
#endif
  real*4, allocatable :: srcnind(:,:,:)
  real*4, allocatable :: srcethane(:,:,:),srcpropane(:,:,:),srcbutane(:,:,:)
  real*4, allocatable :: srcpentane(:,:,:),srchexane(:,:,:),srcethene(:,:,:)
  real*4, allocatable :: srcpropene(:,:,:),srcch2o(:,:,:)
#ifdef SRCBC
  real*4, allocatable :: srcbc(:,:,:),srcoc(:,:,:)
#endif
  real*4, allocatable :: srcald(:,:,:),srcalkanone(:,:,:)
!  real*4, allocatable :: cod50grd(:,:),cod25grd(:,:)
#ifdef CO25D
  real*4, allocatable :: bbcod50grd(:,:),bbcod25grd(:,:)
#else
!  real*4, allocatable :: coanth25grd(:,:),bbcod25grd(:,:)
#endif
  real sorccod5025,xcod50,xcod25,kd50,kd25,sorcbbcod5025,xbbcod50
  real xcoanth25,xbbcod25,sorccoanth25,sorcbbcod25
#ifdef LANDMASKONCE
  integer , allocatable :: landlake(:,:)
#endif
! module savetendmod

  real*4, allocatable :: o3ften_save(:,:,:),o3dten_save(:,:,:)
  real*4, allocatable :: phytendco(:,:,:),phytendn2o(:,:,:),phytendch4(:,:,:)
#ifdef WETDEPDIAG
   real*4,allocatable :: ch4ften_save(:,:,:),ch4dten_save(:,:,:)
   real*4,allocatable :: n2often_save(:,:,:),n2odten_save(:,:,:)
#endif
#ifdef NEWTERM
!  real, allocatable :: o3ften_save_new(:,:,:),o3dten_save_new(:,:,:)
!  real, allocatable :: prodlossnsave(:,:,:,:),prodlossfsave(:,:,:,:)
   real*4,allocatable :: ch4ften_save(:,:,:),ch4dten_save(:,:,:)
   real*4,allocatable :: n2often_save(:,:,:),n2odten_save(:,:,:)
   real*4,allocatable :: ch4sfcten_save(:,:),n2osfcten_save(:,:)
   real*4,allocatable :: o3ten_part(:,:,:,:)
#endif
  real*4, allocatable :: o3dep_save(:,:)
  real*4, allocatable :: noyften_save(:,:,:),noydten_save(:,:,:)
  real*4, allocatable :: noydep_save(:,:),coften_save(:,:,:)
  real*4, allocatable :: codten_save(:,:,:),codep_save(:,:)
  real*4, allocatable :: lnox_save(:,:,:),no2ten(:,:,:)
  real*4, allocatable :: tau6hr(:,:,:)
#ifdef CLD6HR
  real*4, allocatable :: cld6hr(:,:,:)
#endif

  real*4, allocatable :: o3vmr_inst(:,:,:),o3mr(:,:,:)
  real*4, allocatable :: oh_inst(:,:,:),ho2_inst(:,:,:)
  real*4, allocatable :: no_inst(:,:,:),bro_inst(:,:,:)
  real*4, allocatable :: jo1d_inst(:,:,:),jno_inst(:,:,:)
  real*4, allocatable :: colnox_save(:,:)
  real, allocatable :: pcolnox(:,:)
  integer navg,navgfire
! end module savetendmod
! module cldmod
!  real*4, allocatable :: tauxcl3d(:,:,:),tauxci3d(:,:,:)
  real*4, allocatable :: clw3dconv(:,:,:),clw3dlay(:,:,:)
!  real*4, allocatable :: cloud3d(:,:,:)
  !real*4, allocatable,dimension(:,:,:) :: taucldfrc_liq,taucldfrc_ice
  real*4, allocatable,dimension(:,:,:) :: taucldfrc
#ifdef OLDTAUWATER
  real*4, allocatable :: tau3dl(:,:,:),tau3dc(:,:,:)
#endif
  real*4, allocatable :: xgrid(:,:),ygrid(:,:)
  real(CHEM_KIND_R4), allocatable :: bbco_d(:,:),firehtkm_d(:,:)
  real(CHEM_KIND_R4), allocatable :: bbco_hr(:,:,:)
  real(CHEM_KIND_R4), allocatable :: colemisco(:,:),colemisco_chem(:,:)
  real(chem_kind_r4),allocatable :: colemisoc(:,:)
  real(CHEM_KIND_R4), allocatable :: emisco3d(:,:,:),co3dtend(:,:,:)
  real(CHEM_KIND_R4), allocatable :: emisco3dave(:,:,:)
  integer,allocatable :: ktopco(:,:),kbotco(:,:)
  logical doplumerise
  logical :: do85percentvoc=.false.
  integer plumerisefire_frq,wetdep_ls_opt,chem_conv_tr
  real*4, allocatable :: bbco_n(:,:),firehtkm_n(:,:)
  real*4, allocatable :: pblht(:,:)
  real*4, allocatable :: slmsk2d(:,:),ustar(:,:),cmdrag(:,:),areagrd(:,:)
  real*4, allocatable :: w10(:,:)
  real*4, allocatable :: gsiinc(:,:,:,:)
  real*4, allocatable :: incMLS3d(:,:,:),incOMI3d(:,:,:)
  real*4, allocatable :: incnucapsco3d(:,:,:),inctropomico3d(:,:,:)
!  real*4, allocatable :: percov(:,:,:),intco(:,:,:)
  real*4, allocatable :: pergsiv(:,:,:,:),colgsiinc(:,:,:,:)
   real*4 :: brintgsmall
  integer,parameter :: idbc1=1
! error found 9/13/2021
!  integer,parameter :: idbc2=2
!  integer,parameter :: idoc1=3
  integer,parameter :: idbc2=3
  integer,parameter :: idoc1=2
  integer,parameter :: idoc2=4
  integer,parameter :: idso4=5
  integer nlnair,nlnlgt,nl_clim,nl_clim_dyn
  parameter (nlnair=19,nlnlgt=16)
  parameter (nl_clim=19,nl_clim_dyn=24)
!  real*4, parameter,dimension(5) :: tcmw=[12.01,16.80,12.01,16.80,96.0]
  real*4, parameter,dimension(5) :: tcmw=[12.01,16.80,12.01,16.80,96.06] ! ajl use more accurate value
  integer, allocatable :: lconvcld(:,:)
  real*4, allocatable :: ccthk(:,:)
  real*4, allocatable :: covermx(:,:,:),oxvermx(:,:,:)
! move these from jv_
  integer,parameter :: ndust=5,naer=25
! for fastj
  integer,parameter :: ns=55,nw=18,np=56,mx=33,jppj=53,jpmax=55
  integer lpar,jpnl,nb
  type fasttype
    real,allocatable,dimension(:) :: tj,pj,dm,do3,dbc,z,valj
    real rflect,u0,tanht
    real,allocatable,dimension(:,:) :: fff,zj,aer,amf
    integer,allocatable,dimension(:) :: jadsub
  end type
  type(fasttype),allocatable :: fastj
  logical,allocatable :: firstmap(:)
  logical :: docofire85percent=.false.,donofire85percent=.false.
! for nox fire source correction for FIREX-AQ TROPOMI-EXP
  real*4,allocatable :: bb_nox_scale_factor(:,:),bb_nox_adjust(:,:)
! contains 
!end module cldmod
! set up variables for atmos_loc thread private
  
  contains
    subroutine allocatechemcommperm(its,ite,jts,jte)
    use raqmschem_pmgrid_mod, only : nc,nr=>plat,nl=>plev,begj,endj,iam,ibeg,iend
    use raqmschem_pmgrid_mod, only : nmonHTAP,tile
!    use raqmschem_species_mod, only : nsol
    implicit none
    integer its,ite,jts,jte
    character *10 cenv
    cenv=' '
    call getenv('COFIRE85PERCENT',cenv)
    if(cenv.eq.'YES')then
      docofire85percent=.true.
      if(iam.eq.0.and.tile.eq.1)then
        write(6,*)'DO COFIRE 85PERCENT'
      endif
    endif
    cenv=' '
    call getenv('NOFIRE85PERCENT',cenv)
    if(cenv.eq.'YES')then
      donofire85percent=.true.
      if(iam.eq.0.and.tile.eq.1)then
        write(6,*)'DO NOFIRE 85PERCENT'
      endif
    endif

!#include <comm_defs_24t.chem>
    brintgsmall=1.e-15
    kd50=1./(50.*24.*3600.)
    kd25=1./(25.*24.*3600.)
    allocate (dpo3(nc,begj:endj,12),dphno3(nc,begj:endj,12))
    allocate (dph2o2(nc,begj:endj,12),dpno2(nc,begj:endj,12))
    allocate (dpco(nc,begj:endj,12))
!   write(6,*)'allocate srcnind',iam,'nc',nc,begj,endj
!   call flush(6)
    dpo3=0.0
    dphno3=0.0
    dph2o2=0.0
    dpno2=0.0
    dpco=0.0
    allocate (firstmap(nr))
    firstmap=.true.
!   for wetdep
    allocate (rcav_save(ibeg:iend,begj:endj),rnav_save(ibeg:iend,begj:endj))
!    allocate (angledx(its:ite,jts:jte),angledy(its:ite,jts:jte))
    allocate (usina(its:ite,jts:jte),ucosa(its:ite,jts:jte), &
      sinb(its:ite,jts:jte),cosb(its:ite,jts:jte))
    allocate (lprate(ibeg:iend,begj:endj),cprate(ibeg:iend,begj:endj))
    ncplprate=0
    lprate=0.0
    cprate=0.0
#ifdef WETDEPDIAG
    nwetdepdiag=0
    allocate (wetc(nc,begj:endj,nl,14),wetl(nc,begj:endj,nl,14))
    allocate (wetccol(nc,begj:endj,14),wetlcol(nc,begj:endj,14))
    allocate (wetcf(nc,begj:endj,nl,4),wetlf(nc,begj:endj,nl,4))
    allocate (lpaveij(ibeg:iend,begj:endj))
    allocate (cpaveij(ibeg:iend,begj:endj))
    wetc=0.0
    wetccol=0.0
    wetlcol=0.0
    wetl=0.0
    wetcf=0.0
    wetlf=0.0
    lpaveij=0.0
    cpaveij=0.0
    allocate(gsdcol(nc,begj:endj,4))
    allocate(gsd3d(nc,begj:endj,nl,4))
    gsdcol=0.0
    gsd3d=0.0
#endif
    allocate (ltb(nl))
    allocate (i1grd(nc,begj:endj),j1grd(nc,begj:endj),dist1grd(nc,begj:endj))
    allocate (srcnbio(nc,begj:endj,12))
    allocate (srcnsoil(nc,begj:endj,12),srcnair(nc,begj:endj,nlnair))
    allocate (pksrcnair(0:nlnair))
    allocate (sourceno(nc,begj:endj))
    srcnbio=0.0
    srcnsoil=0.0
    srcnair=0.0
!   new aircraft sources
    allocate (srccoair(nc,begj:endj,nlnair))
    srccoair=0.0
#ifdef AERCO
    allocate (srcbcair(nc,begj:endj,nlnair),srcocair(nc,begj:endj,nlnair))
    srcbcair=0.0
    srcocair=0.0
    allocate (srcso2air(nc,begj:endj,nlnair))
    srcso2air=0.0
#endif
    allocate (slghtng(nc,begj:endj,12),srcco_totl(nc,begj:endj,12))
    slghtng=0.0
    srcco_totl=0.0
    allocate (bbefac(nc,begj:endj,2),scale_factor(nc,begj:endj))
!  real*4,allocatable :: bb_nox_scale_factor
    bbefac=0.0
    scale_factor=0.0
    allocate (srcisop(nc,begj:endj,12),srcterp(nc,begj:endj,12))
    srcisop=0.0
    srcterp=0.0
    allocate (ch4clim(nc,begj:endj,12))
    ch4clim=0.0
!    write(6,*)'allocate szagrd',nc,begj,endj
!    call flush(6)
    allocate (szagrd(nc,begj:endj))
!    write(6,*)'spgrd'
!    call flush(6)
    allocate (spgrd(nc,begj:endj))
!    write(6,*)'szgrd'
!    call flush(6)
!    allocate (szgrd(nc,begj:endj))
!    write(6,*)'zero'
!    call flush(6)
    allocate (zsurf(nc,begj:endj))
    allocate (latgrd(nc,begj:endj),longrd(nc,begj:endj))
    szagrd=0.0
    spgrd=0.0
    zsurf=0.0
!    szgrd=0.0
!    write(6,*)'allocate albgrd',nc,begj,endj
!    call flush(6)
    allocate (albgrd(nc,begj:endj),zlwigrd(nc,begj:endj))
    albgrd=0.0
    zlwigrd=0.0
!    write(6,*)'allcoate dirungrd'
!    call flush(6)
    allocate (diurngrd(nc,begj:endj))
    allocate (tskin(nc,begj:endj),sfcwind(nc,begj:endj))
    tskin=0.0
#ifdef AERO
    allocate (seaicethk(nc,begj:endj))
    seaicethk=0.0
#endif
    diurngrd=0.0
!   new HTAP dimsions
    allocate (srcnind(nc,begj:endj,nmonHTAP))
    srcnind=0.0
    allocate (srcethane(nc,begj:endj,nmonHTAP),srcpropane(nc,begj:endj,nmonHTAP))
    srcethane=0.0
    srcpropane=0.0
    allocate (srcbutane(nc,begj:endj,nmonHTAP),srcpentane(nc,begj:endj,nmonHTAP))
    srcbutane=0.0
    srcpentane=0.0
    allocate (srchexane(nc,begj:endj,nmonHTAP),srcethene(nc,begj:endj,nmonHTAP))
    srchexane=0.0
    srcethene=0.0
    allocate (srcpropene(nc,begj:endj,nmonHTAP),srcch2o(nc,begj:endj,nmonHTAP))
    srcpropene=0.0
    srcch2o=0.0
    allocate (srcald(nc,begj:endj,nmonHTAP),srcalkanone(nc,begj:endj,nmonHTAP))
    srcald=0.0
    srcalkanone=0.0
#ifdef SRCBC
    allocate (srcbc(nc,begj:endj,nmonHTAP),srcoc(nc,begj:endj,nmonHTAP))
    srcbc=0.0
    srcoc=0.0
#endif
    allocate (pblht(nc,begj:endj))

!   subroutine allocatesavetend
    allocate (o3ften_save(nc,begj:endj,nl),o3dten_save(nc,begj:endj,nl))
    allocate (phytendco(nc,nl,begj:endj),phytendn2o(nc,nl,begj:endj),phytendch4(nc,nl,begj:endj))
#ifdef WETDEPDIAG
     allocate(ch4ften_save(nc,begj:endj,nl),ch4dten_save(nc,begj:endj,nl))
     allocate(n2often_save(nc,begj:endj,nl),n2odten_save(nc,begj:endj,nl))
#endif
#ifdef NEWTERM
!    allocate (o3ften_save_new(nc,begj:endj,nl),o3dten_save_new(nc,begj:endj,nl))
!    allocate (prodlossnsave(nc,begj:endj,nl,4),prodlossfsave(nc,begj:endj,nl,34))
!    allocate (prodlossnsave(nc,begj:endj,nl,4),prodlossfsave(nc,begj:endj,nl,2))
     allocate(ch4ften_save(nc,begj:endj,nl),ch4dten_save(nc,begj:endj,nl))
     allocate(n2often_save(nc,begj:endj,nl),n2odten_save(nc,begj:endj,nl))
     allocate(ch4sfcten_save(nc,begj:endj),n2osfcten_save(nc,begj:endj))
     allocate(o3ten_part(nc,begj:endj,nl,4))
#endif
    allocate (o3dep_save(nc,begj:endj))
    allocate (noyften_save(nc,begj:endj,nl),noydten_save(nc,begj:endj,nl))
    allocate (noydep_save(nc,begj:endj),coften_save(nc,begj:endj,nl))
    allocate (codten_save(nc,begj:endj,nl),codep_save(nc,begj:endj))
    allocate (lnox_save(nc,begj:endj,nl),colnox_save(nc,begj:endj))
    allocate (no2ten(nc,begj:endj,nl))
    allocate (pcolnox(nc,begj:endj))
!    write(6,*)'allocate xgrid nc,nr',nc,nr,'begj',begj,endj
!    call flush(6)
    allocate (xgrid(nc,begj:endj),ygrid(nc,begj:endj))
    allocate (bbco_d(nc,begj:endj),firehtkm_d(nc,begj:endj))
!    allocate (bbco_hr(nc,begj:endj,24))
!    allocate (colemisco(nc,begj:endj),colemisco_chem(nc,begj:endj))
    allocate (bbco_n(nc,begj:endj),firehtkm_n(nc,begj:endj))
    allocate (slmsk2d(nc,begj:endj))
    allocate (w10(nc,begj:endj))
    allocate (ustar(nc,begj:endj),cmdrag(nc,begj:endj),areagrd(nc,begj:endj))
!    write(200+iam,*)'allocate cmdrag',nc,begj,":",endj
!    call flush(200+iam)
!   return
!   end subroutine allocatesavetend
!   subroutine allocatecldmod
!    allocate (tauxcl3d(nc,begj:endj,nl),tauxci3d(nc,begj:endj,nl))
    allocate (clw3dconv(nc,begj:endj,nl),clw3dlay(nc,begj:endj,nl))
#ifdef DIAGCLDFRC
    write(300+iam,*)'did allocate clw3dconv',nc,begj,":",endj,'nl',nl
    call flush(300+iam)
#endif
!    allocate (cloud3d(nc,begj:endj,nl))
!    allocate (taucldfrc_liq(nc,begj:endj,nl),taucldfrc_ice(nc,begj:endj,nl))
    allocate (taucldfrc(nc,begj:endj,nl))
    taucldfrc=0.0
#ifdef OLDTAUWATER
    allocate (tau3dl(nc,begj:endj,nl),tau3dc(nc,begj:endj,nl))
#endif
#ifdef CLD6HR
    allocate (cld6hr(nc,begj:endj,nl))
#endif
    allocate (tau6hr(nc,begj:endj,nl))
    tau6hr=0.0
    allocate (covermx(nc,begj:endj,nl),oxvermx(nc,begj:endj,nl))
!   return
!   end subroutine allocatecldmod

    return
    end subroutine allocatechemcommperm
    subroutine deallocatechemcommperm
    use raqmschem_pmgrid_mod, only : gsicem,loccem
    write(6,*)'call deallocatechemcomperm'
    call flush(6)
    if(allocated(gsicem))then
      deallocate(gsicem,loccem)
    endif
    if(allocated(dpo3))then
    deallocate (dpo3,dphno3,dph2o2,dpno2,dpco)
    endif
    if(allocated(srcnbio))then
    deallocate (srcnbio,srcnsoil,srcnair)
    endif
    deallocate (pksrcnair)
!   new aircraft sources
    deallocate (srccoair)
#ifdef AERCO
    deallocate(srcbcair)
    deallocate(srcocair)
    deallocate(seaicethk)
#endif
    deallocate (tskin,sfcwind)
    deallocate (srcso2air,slghtng,srcco_totl,bbefac,scale_factor,srcisop,srcterp, &
    ch4clim,szagrd,spgrd,albgrd,zlwigrd,diurngrd)
    deallocate (latgrd,longrd)
!   new HTAP dimsions
    deallocate (srcnind,srcethane,srcpropane,srcbutane,srcpentane)
    deallocate (srchexane,srcethene,srcpropene,srcch2o,srcald,srcalkanone)
    deallocate (i1grd,j1grd,dist1grd)
#ifdef SRCBC
    ,srcbc,srcoc)
#endif
    if(allocated(covermx))then
      deallocate(covermx,oxvermx)
    endif

!   subroutine allocatesavetend
    deallocate (o3ften_save,o3dten_save,phytendco,phytendn2o,phytendch4)
#ifdef WETDEPDIAG
     deallocate(ch4ften_save,ch4dten_save,n2often_save,n2odten_save) 
#endif
#ifdef NEWTERM
     deallocate(ch4ften_save,ch4dten_save,n2often_save,n2odten_save,ch4sfcten_save,n2osfcten_save), &
     o3ten_part)
#endif
    deallocate (o3dep_save,noyften_save,noydten_save,noydep_save,coften_save,codten_save,codep_save)
    deallocate (lnox_save,colnox_save,pcolnox,xgrid,ygrid)
    deallocate (bbco_d,firehtkm_d,bbco_n,firehtkm_n)
    if(allocated (bbco_hr))then
      deallocate (bbco_hr)
    endif
    if(allocated(ktopco))then
      deallocate(ktopco)
    endif
    if(allocated(kbotco))then
      deallocate(kbotco)
    endif
    if(allocated(emisco3d))then
      deallocate(emisco3d)
    endif
    if(allocated(emisco3dave))then
      deallocate(emisco3dave)
    endif
    if(allocated(co3dtend))then
      deallocate(co3dtend)
    endif
    if(allocated(colemisco))then
      deallocate(colemisco)
    endif
    if(allocated(colemisoc))then
      deallocate(colemisoc)
    endif
    if(allocated(colemisco_chem))then
      deallocate(colemisco_chem)
    endif
    deallocate (slmsk2d,w10)
!   return
!   end subroutine allocatesavetend
!   subroutine allocatecldmod
!    deallocate (tauxcl3d,tauxci3d,clw3dconv,clw3dlay)
    deallocate (clw3dconv,clw3dlay)
#ifdef DIAGCLDFRC
    write(300+iam,*)'did allocate clw3dconv',nc,begj,":",endj,'nl',nl
    call flush(300+iam)
#endif
!    deallocate (cloud3d)
!    deallocate (taucldfrc_liq,taucldfrc_ice)
    deallocate (taucldfrc)
#ifdef OLDTAUWATER
    deallocate (tau3dl,tau3dc)
#endif
#ifdef CLD6HR
    deallocate (cld6hr)
#endif
    deallocate (tau6hr)
    return
    end subroutine deallocatechemcommperm
    subroutine allocatechemcommtemp
!    use pmgrid, only : nc,nr=>plat,nl=>plev,begj,endj
    use raqmschem_pmgrid_mod, only : nc,nr=>plat,nl=>plev,begj,endj
#ifdef ALLOCATEONCECHEM
    logical firstalloc
    save firstalloc
    data firstalloc/.true./
    if(.not.firstalloc)return
    firstalloc=.false.
#endif
    allocate (tgrd(nc,begj:endj,nl),pgrd(nc,begj:endj,nl))
    allocate (dpmgrd(nc,begj:endj,nl))
    allocate (zgrd(nc,begj:endj,nl),zeupgrd(nc,begj:endj,nl))
    allocate (ugrd(nc,begj:endj,nl),vgrd(nc,begj:endj,nl))
    allocate (thgrd(nc,begj:endj,nl),qgrd(nc,begj:endj,nl))
    allocate (slp(nc,begj:endj),thbl(nc,begj:endj))
    allocate (rcraingrd(nc,begj:endj,nl))
    allocate (o3vmr_inst(nc,begj:endj,nl),o3mr(nc,begj:endj,nl))
    allocate (oh_inst(nc,begj:endj,nl),ho2_inst(nc,begj:endj,nl))
    allocate (no_inst(nc,begj:endj,nl),bro_inst(nc,begj:endj,nl))
    allocate (jo1d_inst(nc,begj:endj,nl),jno_inst(nc,begj:endj,nl))
#ifdef CO25D
#ifdef COTRACERS
!    allocate (cod25grd(nc,nl))
    allocate (bbcod25grd(nc,nl))
#endif
#ifdef COTRACERS50
!    allocate (cod50grd(nc,nl))
    allocate (bbcod50grd(nc,nl))
#endif
#else
!    allocate (coanth25grd(nc,nl))
!    allocate (bbcod25grd(nc,nl))
#endif
    allocate (lconvcld(nc,begj:endj),ccthk(nc,begj:endj))
    allocate (fastj)
    nb=nl+1
    jpnl=nl
    lpar=nl
    
    allocate (fastj%tj(nb),fastj%pj(nb+1),fastj%dm(nb),fastj%do3(nb), &
              fastj%dbc(nb),fastj%z(nb),fastj%aer(mx,nb),fastj%amf(nb,nb))
    allocate (fastj%fff(nw,jpnl),fastj%valj(ns),fastj%jadsub(2*nb),fastj%zj(jpnl,jppj))
    return
    end subroutine allocatechemcommtemp
    subroutine deallocatechemcommtemp
#ifdef ALLOCATEONCECHEM
    return
#endif
    deallocate (tgrd,pgrd,zgrd,zeupgrd,ugrd,vgrd,thgrd,slp)
    deallocate (rcraingrd,qgrd,thbl,o3vmr_inst,o3mr)
    deallocate (oh_inst,ho2_inst,no_inst,bro_inst,jo1d_inst,jno_inst)
    deallocate (dpmgrd)
#ifdef CO25D
#ifdef COTRACERS
!    deallocate (cod25grd)
    deallocate (bbcod25grd)
#endif
#ifdef COTRACERS50
!    deallocate (cod50grd)
    deallocate (bbcod50grd)
#endif
#else
!    deallocate (coanth25grd,bbcod25grd)
#endif
    deallocate (lconvcld,ccthk)
    deallocate (fastj%tj,fastj%pj,fastj%dm,fastj%do3,fastj%dbc,fastj%z, &
                fastj%aer,fastj%amf)
    deallocate (fastj)
    return
    end subroutine deallocatechemcommtemp
end module raqmschemcomm_mod
