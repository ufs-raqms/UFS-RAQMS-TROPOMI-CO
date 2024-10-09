module raqmschem_map_mod
use raqmschemlocaltype_mod
contains
subroutine maptochem(jin,chemlocal,t3d_in,nchem,pr3d,prl3d,tk3d,u3d,v3d,ph3d,phl3d)
use raqmschem_pmgrid_mod, only : nc,nr,begj,endj,nlev,iam,ibeg,tile
use chem_types_mod, only : CHEM_KIND_R8
use raqmschem_species_mod
use raqmschemcomm_mod
!real(CHEM_KIND_R8),intent(in) :: t3d_in(nc,begj:endj,nlev,nchem)
!real(CHEM_KIND_R8),intent(in) :: t3d_in(:,:,:,:)
real(CHEM_KIND_R4),intent(in) :: t3d_in(:,:,:,:)
real(CHEM_KIND_R8),intent(in),dimension(:,:,:) :: pr3d,prl3d,tk3d,u3d,v3d,ph3d,phl3d
integer i,j,k,nchem,jin,lb(3),iat
type(chemlocaltype),target :: chemlocal
!data first/96*.true./
!save first
#include <choosechem.h>
#include <chemlocaldefinepointer.h>
#include <chemlocaldefinepointer2.h>
!  write(6,*)'top maptochem jin',jin,'begj',begj,'endj',endj
!  call flush(6)
!  write(6,*)'lbound pgrd',lbound(pgrd)
!  write(6,*)'ubound pgrd',ubound(pgrd)
!  lb=lbound(pgrd)
!  call flush(6)
!  !write(6,*)'shape tr3d_in',shape(tr3d_in) 
!write(6,*)'at 4 nc',nc,nlev,'nchem',nchem,'begj',begj,endj
!  write(100+iam,*)'at 4 nc',nc
!  call flush(100+iam)
  !call flush(6)
  call allocatechemlocal(chemlocal,nc,nlev)
!  write(100+iam,*)'call allocatechem2',nc,nl
!  call flush(100+iam)
  call allocatechemlocal2(chemlocal,nc,nlev)
!        write(100+iam,*)'did allocate'
!        call flush(100+iam)
#include <chemlocalsetpointer.h>
#include <chemlocalsetpointer2.h>
!if(jin<lb(2))then
!  write(6,*)'error map jin ',jin,'lb',lb
!  call flush(6)
!endif
!write(200+iam,*)'jin',jin,'begj',begj
!call flush(200+iam)
j=jin-begj+1
!spgrd(:,jin)=prl3d(:,j,1)
!zsurf(:,jin)=phl3d(:,j,1)/9.80616
!write(200+iam,*)'j',j,'lbound',lbound(prl3d),ubound(prl3d)
!call flush(200+iam)
!pgrd(:,jin,:)=prl3d(:,j,:)*.01 ! pascals in need mb for our chemistry
!write(200+iam,*)'pgrd',jin,maxval(pgrd(:,jin,:)),'sp',maxval(spgrd(:,jin))
!call flush(200+iam)
!ugrd(:,jin,:)=u3d(:,j,:)
!vgrd(:,jin,:)=v3d(:,j,:)
!tgrd(:,jin,:)=tk3d(:,j,:)
!thgrd(:,jin,:)=tk3d(:,j,:)*(1000./pgrd(:,jin,:))**.286
!zgrd(:,jin,:)=phl3d(:,j,:)/9.80616
!  if(jin.eq.49.and.iam.eq.7)then
!    do k=1,nlev
!      write(6,*)'pgrd',pgrd(33,jin,k),'zgrd',zgrd(33,jin,k)
!    end do
!  endif
!if(jin.eq.80.and.ibeg.eq.1)then
!  do k=1,nlev
!    write(6,*)'thgrd',k,thgrd(6,jin,k),' T ',tk3d(6,j,k),' P ',pgrd(6,jin,k),' z ',zgrd(6,jin,k)
!  end do
!!  call flush(6)
! endif
do k=2,nlev+1
! zeupgrd is solids starting with one above surface
  zeupgrd(:,jin,k-1)=ph3d(:,j,k)/9.80616+zsurf(:,jin)
!  if(iam.eq.9.and.jin.eq.18)then
!    write(101,*)k-1,zeupgrd(18,jin,k-1)
!    write(6,*)'zeupgrd 18,18 ',k-1,zeupgrd(18,jin,k-1)
!  endif
!  if(iam.eq.9.and.jin.eq.20)then
!    write(102,*)k-1,zeupgrd(18,jin,k-1)
!    write(6,*)'zeupgrd 15,20 ',k-1,zeupgrd(15,jin,k-1)
!  endif
end do
qgrd(:,jin,:)=t3d_in(:,j,:,1)
! ajl add some arrays

do k=1,nlev
 do i=1,nc
   iat=ibeg+i-1
!c  fill rainout array
!c
   if (pgrd(i,jin,k) .gt. 315.) then
!cjjj            rcraingrd(i,j,l)=1./(6.912e+06-5184.*pgrd(i,j,l))
!cjjj rainout rate has been 20-day lifetime @1000mb and 60-day @300mb
!cjjj increase to 10-day lifetime @1000mb (still 60-day @300mb) JAA 11/18/2003
      rcraingrd(i,jin,k)=1./(7.309e+06-6445.*pgrd(i,jin,k))
   else
     rcraingrd(i,jin,k)=0.
   endif
!  backwards ? done in raqms_model_mod
!  dpmgrd(i,jin,k)=(pr3d(i,j,k+1)-pr3d(i,j,k))*.01 ! make mb
  brclgrd(i,k)=t3d_in(i,j,k,p_brcl)
  brno3grd(i,k)=t3d_in(i,j,k,p_brno3)
  brygrd(i,k)=t3d_in(i,j,k,p_bry)
  ccl4grd(i,k)=t3d_in(i,j,k,p_ccl4)
  xmrpgrd(i,k)=t3d_in(i,j,k,p_xmrp)
  prdpgrd(i,k)=t3d_in(i,j,k,p_prdp)
  ch2ogrd(i,k)=t3d_in(i,j,k,p_ch2o)
  ch3brgrd(i,k)=t3d_in(i,j,k,p_ch3br)
  ch3clgrd(i,k)=t3d_in(i,j,k,p_ch3cl)
  ch3oohgrd(i,k)=t3d_in(i,j,k,p_ch3ooh)
  ch4grd(i,k)=t3d_in(i,j,k,p_ch4)
  cl2grd(i,k)=t3d_in(i,j,k,p_cl2)
  clno3grd(i,k)=t3d_in(i,j,k,p_clno3)
  clygrd(i,k)=t3d_in(i,j,k,p_cly)
  cogrd(i,k)=t3d_in(i,j,k,p_co)
  f11grd(i,k)=t3d_in(i,j,k,p_f11)
  f1211grd(i,k)=t3d_in(i,j,k,p_f1211)
  f12grd(i,k)=t3d_in(i,j,k,p_f12)
  f1301grd(i,k)=t3d_in(i,j,k,p_f1301)
  h2o2grd(i,k)=t3d_in(i,j,k,p_h2o2)
  vrpgrd(i,k)=t3d_in(i,j,k,p_vrp)
  hbrgrd(i,k)=t3d_in(i,j,k,p_hbr)
  hclgrd(i,k)=t3d_in(i,j,k,p_hcl)
  ripgrd(i,k)=t3d_in(i,j,k,p_rip)
  hno3tgrd(i,k)=t3d_in(i,j,k,p_hno3t)
  hno4grd(i,k)=t3d_in(i,j,k,p_hno4)
  hobrgrd(i,k)=t3d_in(i,j,k,p_hobr)
  hoclgrd(i,k)=t3d_in(i,j,k,p_hocl)
  xmtcfmgrd(i,k)=t3d_in(i,j,k,p_mtcfm)
  xn2o5grd(i,k)=t3d_in(i,j,k,p_n2o5)
  xn2ogrd(i,k)=t3d_in(i,j,k,p_n2o)
  xno2grd(i,k)=t3d_in(i,j,k,p_no2)
  xno3grd(i,k)=t3d_in(i,j,k,p_no3)
  xnoygrd(i,k)=t3d_in(i,j,k,p_noy)
#if 0
  if(tile.eq.1.and.iin.eq.139.and.jin.eq.101.and.k.eq.1)then
    write(200+iam,*)'xo2in',xno2grd(i,k),' noy ',xnoygrd(i,k)
  endif
#endif
  oclogrd(i,k)=t3d_in(i,j,k,p_oclo)
  oxgrd(i,k)=t3d_in(i,j,k,p_ox)
!  if(first(jin))then
!  if(iam.eq.9.and.i.eq.18.and.j.eq.18)then
!    write(101,*)i,j,k,pgrd(i,jin,k),oxgrd(i,k),xno2grd(i,k)
!    call flush(101)
!  endif
!  if(iam.eq.9.and.i.eq.15.and.j.eq.20)then
!    write(102,*)i,j,k,pgrd(i,jin,k),oxgrd(i,k),xno2grd(i,k)
!    call flush(102)
!  endif
  !endif
  c2h6grd(i,k)=t3d_in(i,j,k,p_c2h6)
  ald2grd(i,k)=t3d_in(i,j,k,p_ald2)
  ethoohgrd(i,k)=t3d_in(i,j,k,p_ethooh)
  pangrd(i,k)=t3d_in(i,j,k,p_pan)
  pargrd(i,k)=t3d_in(i,j,k,p_par)
  xonitgrd(i,k)=t3d_in(i,j,k,p_onit)
  aonegrd(i,k)=t3d_in(i,j,k,p_aone)
  roohgrd(i,k)=t3d_in(i,j,k,p_rooh)
  xmglygrd(i,k)=t3d_in(i,j,k,p_mgly)
  ethgrd(i,k)=t3d_in(i,j,k,p_eth)
  xoletgrd(i,k)=t3d_in(i,j,k,p_olet)
!  if (firstmap(jin))then
  if(p_olei>0)then
    xoleigrd(i,k)=0.0
  else
    xoleigrd(i,k)=0.0
  endif
!  !else
!    xoleigrd(i,k)=t3d_in(i,j,k,p_olei)
!  endif
  xisopgrd(i,k)=t3d_in(i,j,k,p_isop)
  xisoprdgrd(i,k)=t3d_in(i,j,k,p_isoprd)
  prop_pargrd(i,k)=t3d_in(i,j,k,p_propar)
  ch3ohgrd(i,k)=t3d_in(i,j,k,p_ch3oh)
  xmvkgrd(i,k)=t3d_in(i,j,k,p_mvk)
  xmacrgrd(i,k)=t3d_in(i,j,k,p_macr)
  xmpangrd(i,k)=t3d_in(i,j,k,p_mpan)
  hcngrd(i,k)=t3d_in(i,j,k,p_hcn)
  ch3cngrd(i,k)=t3d_in(i,j,k,p_ch3cn)
  so2grd(i,k)=t3d_in(i,j,k,p_so2)
  so4aergrd(i,k)=t3d_in(i,j,k,p_sulf)
  dmsgrd(i,k)=t3d_in(i,j,k,p_dms)
  msagrd(i,k)=t3d_in(i,j,k,p_msa)
  if(p_no3aer>0)then
    no3aergrd(i,k)=t3d_in(i,j,k,p_no3aer)
  else
    no3aergrd(i,k)=0.0
  endif
  if(p_nh3>0)then
    nh3grd(i,k)=t3d_in(i,j,k,p_nh3)
  else
    nh3grd(i,k)=0.0
  endif
  if(p_nh4aer>0)then
    nh4aergrd(i,k)=t3d_in(i,j,k,p_nh4aer)
  else
    nh4aergrd(i,k)=0.0
  endif
  bc1grd(i,k)=t3d_in(i,j,k,p_bc1)
  bc2grd(i,k)=t3d_in(i,j,k,p_bc2)
  oc1grd(i,k)=t3d_in(i,j,k,p_oc1)
  oc2grd(i,k)=t3d_in(i,j,k,p_oc2)
  du1grd(i,k)=t3d_in(i,j,k,p_dust1)
  du2grd(i,k)=t3d_in(i,j,k,p_dust2)
  du3grd(i,k)=t3d_in(i,j,k,p_dust3)
  du4grd(i,k)=t3d_in(i,j,k,p_dust4)
  du5grd(i,k)=t3d_in(i,j,k,p_dust5)
  ss1grd(i,k)=t3d_in(i,j,k,p_seas1)
  ss2grd(i,k)=t3d_in(i,j,k,p_seas2)
  ss3grd(i,k)=t3d_in(i,j,k,p_seas3)
  ss4grd(i,k)=t3d_in(i,j,k,p_seas4)
  ss5grd(i,k)=t3d_in(i,j,k,p_seas5)
#ifdef CO25D
  if(p_co50d>0)then
    cod50grd(i,k)=t3d_in(i,j,k,p_co50d)
  else
    cod50grd(i,k)=0.0
  endif
  if(p_co25d>0)then
    cod25grd(i,k)=t3d_in(i,j,k,p_co25d)
  else
    cod25grd(i,k)=0.0
  endif
#else
  if(p_coanth25>0)then
    coanth25grd(i,k)=t3d_in(i,j,k,p_coanth25)
  endif
  if(p_bbcod25>0)then
    bbcod25grd(i,k)=t3d_in(i,j,k,p_bbcod25)
  endif
#endif
  br2grd(i,k)=t3d_in(i,j,k,p_br2)
!  if(iam.eq.11)write(6,*)'br2grd in',br2grd(i,k),'i,',i,j,k,p_br2
 end do
end do
#ifdef DIAGPROF
if(firstmap(jin))then
  if(iam.eq.9.and.jin.eq.18)then
    do k=1,nlev
    write(101,*)k,pgrd(18,jin,k),zgrd(18,jin,k),zeupgrd(18,jin,k),tgrd(18,jin,k)
    end do
    do k=1,nlev
      write(101,*)k,pgrd(18,jin,k),oxgrd(18,k),xno2grd(18,k),qgrd(18,jin,k)
    end do
    call flush(101)
  endif
  if(iam.eq.9.and.jin.eq.20)then
    do k=1,nlev
    write(102,*)k,pgrd(15,jin,k),zgrd(15,jin,k),zeupgrd(15,jin,k),tgrd(15,jin,k)
    end do
    do k=1,nlev
      write(102,*)k,pgrd(15,jin,k),oxgrd(15,k),xno2grd(15,k),qgrd(15,jin,k)
    end do
    call flush(102)
  endif
endif
#endif
firstmap(jin)=.false.
return
end subroutine maptochem
subroutine mapfromchem(j,chemlocal,t3d_out,nchem)
use raqmschem_pmgrid_mod, only : nc,nr,begj,endj,nlev,tile
use chem_types_mod, only : CHEM_KIND_R8,CHEM_KIND_R4
use raqmschem_species_mod
!real(CHEM_KIND_R4),intent(out) :: t3d_out(nc,begj:endj,nlev,nchem)
real(CHEM_KIND_R8),intent(out) :: t3d_out(nc,begj:endj,nlev,nchem)
integer i,j,k,nchem
type(chemlocaltype),target :: chemlocal
#include <chemlocaldefinepointer.h>
#include <chemlocaldefinepointer2.h>
#include <chemlocalsetpointer.h>
#include <chemlocalsetpointer2.h>
do k=1,nlev
 do i=1,nc
  t3d_out(i,j,k,p_brcl)=brclgrd(i,k)
  t3d_out(i,j,k,p_brno3)=brno3grd(i,k)
  t3d_out(i,j,k,p_bry)=brygrd(i,k)
  t3d_out(i,j,k,p_ccl4)=ccl4grd(i,k)
  t3d_out(i,j,k,p_xmrp)=xmrpgrd(i,k)
  t3d_out(i,j,k,p_prdp)=prdpgrd(i,k)
  t3d_out(i,j,k,p_ch2o)=ch2ogrd(i,k)
  t3d_out(i,j,k,p_ch3br)=ch3brgrd(i,k)
  t3d_out(i,j,k,p_ch3cl)=ch3clgrd(i,k)
  t3d_out(i,j,k,p_ch3ooh)=ch3oohgrd(i,k)
  t3d_out(i,j,k,p_ch4)=ch4grd(i,k)
  t3d_out(i,j,k,p_cl2)=cl2grd(i,k)
#if 0
  if(isnan(cl2grd(i,k)))then
   write(6,*)'cl2grd',i,j,k,p_cl2
  endif
#endif
  t3d_out(i,j,k,p_clno3)=clno3grd(i,k)
  t3d_out(i,j,k,p_cly)=clygrd(i,k)
  t3d_out(i,j,k,p_co)=cogrd(i,k)
  t3d_out(i,j,k,p_f11)=f11grd(i,k)
  t3d_out(i,j,k,p_f1211)=f1211grd(i,k)
  t3d_out(i,j,k,p_f12)=f12grd(i,k)
  t3d_out(i,j,k,p_f1301)=f1301grd(i,k)
  t3d_out(i,j,k,p_h2o2)=h2o2grd(i,k)
  t3d_out(i,j,k,p_vrp)=vrpgrd(i,k)
  t3d_out(i,j,k,p_hbr)=hbrgrd(i,k)
  t3d_out(i,j,k,p_hcl)=hclgrd(i,k)
  t3d_out(i,j,k,p_rip)=ripgrd(i,k)
  t3d_out(i,j,k,p_hno3t)=hno3tgrd(i,k)
  t3d_out(i,j,k,p_hno4)=hno4grd(i,k)
  t3d_out(i,j,k,p_hobr)=hobrgrd(i,k)
  t3d_out(i,j,k,p_hocl)=hoclgrd(i,k)
  t3d_out(i,j,k,p_mtcfm)=xmtcfmgrd(i,k)
#if 0
  if(isnan(t3d_out(i,j,k,p_mtcfm)))then
     write(6,*)'P-mtcfm',p_mtcfm,i,j,k
  endif
#endif
  t3d_out(i,j,k,p_n2o5)=xn2o5grd(i,k)
#if 0
  if(isnan(t3d_out(i,j,k,p_n2o5)))then
     write(6,*)'p-n2o5',p_n2o5,i,j,k
  endif
#endif
  t3d_out(i,j,k,p_n2o)=xn2ogrd(i,k)
  t3d_out(i,j,k,p_no2)=xno2grd(i,k)
  t3d_out(i,j,k,p_no3)=xno3grd(i,k)
  t3d_out(i,j,k,p_noy)=xnoygrd(i,k)
#if 0
  if(tile.eq.1.and.iin.eq.139.and.jin.eq.101.and.k.eq.1)then
    write(200+iam,*)'xo2out',xno2grd(i,k),' noy ',xnoygrd(i,k)
  endif
#endif
  t3d_out(i,j,k,p_oclo)=oclogrd(i,k)
  t3d_out(i,j,k,p_ox)=oxgrd(i,k)
  t3d_out(i,j,k,p_c2h6)=c2h6grd(i,k)
  t3d_out(i,j,k,p_ald2)=ald2grd(i,k)
  t3d_out(i,j,k,p_ethooh)=ethoohgrd(i,k)
  t3d_out(i,j,k,p_pan)=pangrd(i,k)
  t3d_out(i,j,k,p_par)=pargrd(i,k)
  t3d_out(i,j,k,p_onit)=xonitgrd(i,k)
  t3d_out(i,j,k,p_aone)=aonegrd(i,k)
  t3d_out(i,j,k,p_rooh)=roohgrd(i,k)
  t3d_out(i,j,k,p_mgly)=xmglygrd(i,k)
  t3d_out(i,j,k,p_eth)=ethgrd(i,k)
  t3d_out(i,j,k,p_olet)=xoletgrd(i,k)
  if(p_olei>0)then
    t3d_out(i,j,k,p_olei)=xoleigrd(i,k)
  endif
  t3d_out(i,j,k,p_isop)=xisopgrd(i,k)
  t3d_out(i,j,k,p_isoprd)=xisoprdgrd(i,k)
  t3d_out(i,j,k,p_propar)=prop_pargrd(i,k)
  t3d_out(i,j,k,p_ch3oh)=ch3ohgrd(i,k)
  t3d_out(i,j,k,p_mvk)=xmvkgrd(i,k)
  t3d_out(i,j,k,p_macr)=xmacrgrd(i,k)
  t3d_out(i,j,k,p_mpan)=xmpangrd(i,k)
  t3d_out(i,j,k,p_hcn)=hcngrd(i,k)
  t3d_out(i,j,k,p_ch3cn)=ch3cngrd(i,k)
  t3d_out(i,j,k,p_so2)=so2grd(i,k)
#ifdef DO_AEROSOL
  t3d_out(i,j,k,p_sulf)=so4aergrd(i,k)
  t3d_out(i,j,k,p_dms)=dmsgrd(i,k)
  t3d_out(i,j,k,p_msa)=msagrd(i,k)
  if(p_no3aer>0)then
    t3d_out(i,j,k,p_no3aer)=no3aergrd(i,k)
  endif
  if(p_nh3>0)then
    t3d_out(i,j,k,p_nh3)=nh3grd(i,k)
  endif
  if(p_nh4aer>0)then
    t3d_out(i,j,k,p_nh4aer)=nh4aergrd(i,k)
  endif
  t3d_out(i,j,k,p_bc1)=bc1grd(i,k)
  t3d_out(i,j,k,p_bc2)=bc2grd(i,k)
  t3d_out(i,j,k,p_oc1)=oc1grd(i,k)
  t3d_out(i,j,k,p_oc2)=oc2grd(i,k)
  t3d_out(i,j,k,p_dust1)=du1grd(i,k)
  t3d_out(i,j,k,p_dust2)=du2grd(i,k)
  t3d_out(i,j,k,p_dust3)=du3grd(i,k)
  t3d_out(i,j,k,p_dust4)=du4grd(i,k)
  t3d_out(i,j,k,p_dust5)=du5grd(i,k)
  t3d_out(i,j,k,p_seas1)=ss1grd(i,k)
  t3d_out(i,j,k,p_seas2)=ss2grd(i,k)
  t3d_out(i,j,k,p_seas3)=ss3grd(i,k)
  t3d_out(i,j,k,p_seas4)=ss4grd(i,k)
  t3d_out(i,j,k,p_seas5)=ss5grd(i,k)
  t3d_out(i,j,k,p_seas5)=ss5grd(i,k)
#else
! not trully an aerosol so need to return back
  t3d_out(i,j,k,p_seas5)=ss5grd(i,k)
#endif
#ifdef CO25D
  if(p_co50d>0)then
    t3d_out(i,j,k,p_co50d)=cod50grd(i,k)
  endif
  if(p_co25d>0)then
    t3d_out(i,j,k,p_co25d)=cod25grd(i,k)
  endif
#else
  if(p_coanth25>0)then
    t3d_out(i,j,k,p_coanth25)=coanth25grd(i,k)
  endif
  if(p_bbcod25>0)then
    t3d_out(i,j,k,p_bbcod25)=bbcod25grd(i,k)
  endif
#endif
  t3d_out(i,j,k,p_br2)=br2grd(i,k)
!  if(iam.eq.11)write(6,*)'br2 out',br2grd(i,k),i,j,k
 end do
end do
return
end subroutine mapfromchem
end module raqmschem_map_mod
