#include <choosechem.h>
module raqmschemlocaltype_mod
  type chemlocaltype
    real*4,allocatable,dimension(:,:) :: xno2grd,ch4grd,hclgrd,xnoygrd,hno3tgrd,clygrd,xn2o5grd,h2o2grd,clno3grd, &
    xn2ogrd,f11grd,f12grd,ccl4grd,ch3clgrd,xmtcfmgrd,brygrd,ch3brgrd,f1301grd,f1211grd,hno4grd,hoclgrd,cogrd,oclogrd
    real*4,allocatable,dimension(:,:) :: oxgrd,xno3grd,ch2ogrd,ch3oohgrd,hbrgrd,brno3grd,hobrgrd,brclgrd,cl2grd,c2h6grd,&
    ald2grd,ethoohgrd,pangrd,pargrd,xonitgrd,aonegrd,roohgrd,xmglygrd,ethgrd,xoletgrd,xoleigrd,xisopgrd,xisoprdgrd &
    ,prop_pargrd,ch3ohgrd,xmvkgrd,xmacrgrd,xmpangrd
    real*4,allocatable,dimension(:,:) ::  &
#ifdef ISOPRENE_PEROX
        vrpgrd,ripgrd,xmrpgrd,prdpgrd
#else
        h2ogrd,hfgrd,cf2ogrd,cfclogrd
#endif
    real*4,allocatable,dimension(:,:) :: &
    hcngrd,ch3cngrd,so2grd,so4aergrd,dmsgrd,msagrd,no3aergrd,nh3grd,nh4aergrd,bc1grd,bc2grd,oc1grd,oc2grd, &
    du1grd,du2grd,du3grd,du4grd,du5grd,ss1grd,ss2grd,ss3grd,ss4grd,ss5grd
#ifdef CO25D
    real*4,allocatable,dimension(:,:) :: cod50grd,cod25grd,br2grd
#else
    real*4,allocatable,dimension(:,:) :: coanth25grd,bbcod25grd,br2grd
#endif
  end type chemlocaltype
contains
  subroutine allocatechemlocal(localarea,nc,nl)
    integer nc,nl
    type(chemlocaltype) :: localarea
!    write(6,*)'nc',nc,'nl',nl,allocated(localarea%xno2grd),allocated(localarea%cod25grd),allocated(localarea%ch4grd), &
!       allocated(localarea%ccl4grd),allocated(localarea%oxgrd)
!    call flush(6)
    allocate(localarea%xno2grd(nc,nl),localarea%ch4grd(nc,nl),localarea%hclgrd(nc,nl),localarea%xnoygrd(nc,nl),&
    localarea%hno3tgrd(nc,nl),localarea%clygrd(nc,nl),localarea%xn2o5grd(nc,nl),localarea%h2o2grd(nc,nl),&
    localarea%clno3grd(nc,nl),localarea%xn2ogrd(nc,nl),localarea%f11grd(nc,nl),localarea%f12grd(nc,nl),&
    localarea%ccl4grd(nc,nl),localarea%ch3clgrd(nc,nl),localarea%xmtcfmgrd(nc,nl),localarea%brygrd(nc,nl),&
    localarea%ch3brgrd(nc,nl),localarea%f1301grd(nc,nl),localarea%f1211grd(nc,nl),localarea%hno4grd(nc,nl),&
    localarea%hoclgrd(nc,nl),localarea%cogrd(nc,nl),localarea%oclogrd(nc,nl),localarea%xno3grd(nc,nl),&
    localarea%ch2ogrd(nc,nl),localarea%ch3oohgrd(nc,nl),localarea%hbrgrd(nc,nl),localarea%brno3grd(nc,nl),&
    localarea%hobrgrd(nc,nl),localarea%brclgrd(nc,nl),localarea%cl2grd(nc,nl),localarea%c2h6grd(nc,nl),&
    localarea%ald2grd(nc,nl),localarea%ethoohgrd(nc,nl),localarea%pangrd(nc,nl),localarea%pargrd(nc,nl),&
    localarea%xonitgrd(nc,nl),localarea%aonegrd(nc,nl),localarea%roohgrd(nc,nl),localarea%xmglygrd(nc,nl),&
    localarea%ethgrd(nc,nl),localarea%xoletgrd(nc,nl),localarea%xoleigrd(nc,nl),localarea%xisopgrd(nc,nl),&
    localarea%xisoprdgrd(nc,nl),localarea%prop_pargrd(nc,nl),localarea%ch3ohgrd(nc,nl),localarea%xmvkgrd(nc,nl))
    allocate(localarea%oxgrd(nc,nl),localarea%xmpangrd(nc,nl),localarea%xmacrgrd(nc,nl))
#ifdef ISOPRENE_PEROX
    allocate(localarea%vrpgrd(nc,nl),localarea%ripgrd(nc,nl),localarea%xmrpgrd(nc,nl),localarea%prdpgrd(nc,nl))
#else
    allocate(localarea%h2ogrd(nc,nl),localarea%hfgrd(nc,nl),localarea%cf2ogrd(nc,nl),localarea%cfclogrd(nc,nl))
#endif
#ifdef CO25D
    allocate(localarea%cod50grd(nc,nl),localarea%cod25grd(nc,nl),localarea%br2grd(nc,nl))
#else
    allocate(localarea%coanth25grd(nc,nl),localarea%bbcod25grd(nc,nl),localarea%br2grd(nc,nl))
#endif

  end subroutine allocatechemlocal
  subroutine deallocatechemlocal(localarea)
    type(chemlocaltype) :: localarea
    deallocate(localarea%xno2grd,localarea%ch4grd,localarea%hclgrd,localarea%xnoygrd,&
    localarea%hno3tgrd,localarea%clygrd,localarea%xn2o5grd,localarea%h2o2grd,&
    localarea%clno3grd,localarea%xn2ogrd,localarea%f11grd,localarea%f12grd,&
    localarea%ccl4grd,localarea%ch3clgrd,localarea%xmtcfmgrd,localarea%brygrd,&
    localarea%ch3brgrd,localarea%f1301grd,localarea%f1211grd,localarea%hno4grd,&
    localarea%hoclgrd,localarea%cogrd,localarea%oclogrd,localarea%xno3grd,&
    localarea%ch2ogrd,localarea%ch3oohgrd,localarea%hbrgrd,localarea%brno3grd,&
    localarea%hobrgrd,localarea%brclgrd,localarea%cl2grd,localarea%c2h6grd,&
    localarea%ald2grd,localarea%ethoohgrd,localarea%pangrd,localarea%pargrd,&
    localarea%xonitgrd,localarea%aonegrd,localarea%roohgrd,localarea%xmglygrd,&
    localarea%ethgrd,localarea%xoletgrd,localarea%xoleigrd,localarea%xisopgrd,&
    localarea%xisoprdgrd,localarea%prop_pargrd,localarea%ch3ohgrd,localarea%xmvkgrd,&
    localarea%xmacrgrd,localarea%oxgrd)
    deallocate (localarea%xmpangrd)
#ifdef CO25D
    deallocate (localarea%cod50grd,localarea%cod25grd,localarea%br2grd)
#else
    deallocate (localarea%coanth25grd,localarea%bbcod25grd,localarea%br2grd)
#endif
#ifdef ISOPRENE_PEROX
    deallocate (localarea%vrpgrd,localarea%ripgrd,localarea%xmrpgrd,localarea%prdpgrd)
#else
    deallocate (localarea%h2ogrd,localarea%hfgrd,localarea%cf2ogrd,localarea%cfclogrd)
#endif
  end subroutine deallocatechemlocal
  subroutine allocatechemlocal2(localarea,nc,nl)
    integer nc,nl
    type(chemlocaltype) :: localarea
    allocate(localarea%hcngrd(nc,nl),localarea%ch3cngrd(nc,nl),localarea%so2grd(nc,nl),localarea%so4aergrd(nc,nl),&
    localarea%dmsgrd(nc,nl),localarea%msagrd(nc,nl),localarea%no3aergrd(nc,nl),localarea%nh3grd(nc,nl),&
    localarea%nh4aergrd(nc,nl),localarea%bc1grd(nc,nl),localarea%bc2grd(nc,nl),localarea%oc1grd(nc,nl),&
    localarea%oc2grd(nc,nl),localarea%du1grd(nc,nl),localarea%du2grd(nc,nl),localarea%du3grd(nc,nl),&
    localarea%du4grd(nc,nl),localarea%du5grd(nc,nl),localarea%ss1grd(nc,nl),localarea%ss2grd(nc,nl),&
    localarea%ss3grd(nc,nl),localarea%ss4grd(nc,nl),localarea%ss5grd(nc,nl))
  end subroutine allocatechemlocal2
  subroutine deallocatechemlocal2(localarea)
    type(chemlocaltype) :: localarea
    deallocate( localarea%hcngrd,localarea%ch3cngrd,localarea%so2grd,localarea%so4aergrd,localarea%dmsgrd, &
    localarea%msagrd,localarea%no3aergrd,localarea%nh3grd,localarea%nh4aergrd,localarea%bc1grd,localarea%bc2grd, &
    localarea%oc1grd,localarea%oc2grd,localarea%du1grd,localarea%du2grd,localarea%du3grd,localarea%du4grd, &
    localarea%du5grd,localarea%ss1grd,localarea%ss2grd,localarea%ss3grd,localarea%ss4grd,localarea%ss5grd) 
  end subroutine deallocatechemlocal2 
end module raqmschemlocaltype_mod
