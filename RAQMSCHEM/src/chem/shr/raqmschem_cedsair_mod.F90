#include <choosechem.h>
module raqmschem_cedsair_mod
  use chem_types_mod, only : chem_kind_r4,chem_kind_r8
  use raqmschemcomm_mod, only : pgrd,tgrd,qgrd,zeupgrd,zsurf,nlnair,tcmw,airmw
  character*10  ccedsair(3)
  integer :: p_ceds_co,p_ceds_NOx
  integer :: p_ceds_n2o,numcedsair
  integer*4,allocatable  :: ktopcedsair(:,:,:),kbotcedsair(:,:,:)
  real*4, allocatable :: intaltair(:)
  real*4, allocatable :: intcedsair(:,:,:,:)
  real(chem_kind_r4),allocatable,dimension(:,:,:,:) :: cedsairemis
  logical :: lcedsair=.false.
  logical :: cedsairstep=.false.
contains
  subroutine docedsairm(j,dt)
  use raqmschem_pmgrid_mod,only : iam,masterproc,nlev,nc,begj,endj,tile,ibeg,iend
  use chem_raqms_mod, only : cedsairinc
  implicit none
  integer i,j,k,m,kbot,ktop,klo,kup,kk,kkb,ii,lbsrc(3),ubsrc(3)
  real(chem_kind_r4) :: zgeoh(nlev),zztop,zzbot,zup,zlo,dens,dt
  real(chem_kind_r8) :: aflo,afup,intvallo(numcedsair),intvalup,sac
  logical zloset,zupset
  if(.not.allocated(cedsairemis))then
     allocate (cedsairemis(nc,nlev,begj:endj,numcedsair))
  endif
  cedsairemis(:,:,j,:)=0.0
  do i=1,nc
    ii=ibeg+i-1
    ktop=maxval(ktopcedsair(i,j,:))
    if(ktop<1)cycle
    kbot=max(minval(kbotcedsair(i,j,:))-1,0)
    zztop=intaltair(ktop)
    zzbot=intaltair(kbot)
    zgeoh(1)=zsurf(i,j)*.001 ! make km
    do k=2,nlev
      zgeoh(k)=zeupgrd(i,j,k-1)*.001 ! make km
    end do
    zloset=.false.
    klo=-1
    aflo=0.0
    zlo=0.0
    if(zlo>=zztop)cycle ! skip this i,j no emissions
    if(kbot==0)then
      intvallo=0.0
      kkb=1
      zloset=.true.
      aflo=0.0
    else
      do k=1,nlev
        zup=zgeoh(k+1)-zgeoh(1)
        if(zup<=zzbot)cycle
        zlo=zgeoh(k)-zgeoh(1)
        if(zup>=zzbot.or.(zlo<intaltair(kbot+1).and.zlo>=intaltair(kbot)))then
          klo=kbot
          aflo=0.0
          intvallo=0.0
          zloset=.true.
          kkb=k
          exit
        endif
      end do
    endif
      
DOK:  do k=kkb,nlev
      dens=7.2431122e+18*pgrd(i,j,k)/(tgrd(i,j,k)*(1.+.608*qgrd(i,j,k)))
      zupset=.false.
      kup=-1
      afup=0.0
      zlo=zgeoh(k)-zgeoh(1)
      if(zlo>=zztop)exit DOK
      zup=zgeoh(k+1)-zgeoh(1)
      if(zup<=zzbot)cycle
      if(zup>=zztop-.000001)then
        zupset=.true.
        afup=1.
        kup=ktop
      endif
      if(.not.zloset.or..not.zupset)then
      do kk=kbot,ktop
        if(.not.zupset)then
          if(zup<intaltair(kk+1).and.zup>=intaltair(kk))then
            kup=kk
            afup=(zup-intaltair(kk))/(intaltair(kk+1)-intaltair(kk))
            zupset=.true.
            if(kk==ktop)then
                write(6,*)'hit ktop zupset',ktop,'i',i,j,k
            endif
          endif
        endif
        if(zupset)then
           exit
         endif
      end do ! kk
      endif
      if(.not.zupset.or.kup<0)then
         write(6,*)'error ',i,j,k,'tile',tile
         write(6,*)'klo',klo,kup,zupset,zloset
         write(6,*)'zlo',zlo,'zup',zup
         flush(6)
         write(6,*)'kbot',kbot,ktop
         write(6,*)'zztop',zztop,zzbot
         flush(6)
         do kk=kbot,ktop
           write(6,*)i,j,k,'intaltair',kk,intaltair(kk)
           flush(6)
         end do
         call killit('aa error')
      endif
      do m=1,numcedsair
        intvalup=intcedsair(i,j,kup+1,m)*afup+intcedsair(i,j,kup,m)*(1.-afup)
!       here zup, zlo are in km
        sac=(intvalup-intvallo(m))/(zup-zlo)*1.e-5 ! km to cm
          if(sac<-1.e-9)then
             write(6,*)i,j,k,m,'tile',tile
             write(6,*)'sacneg',sac,'intvallo',intvallo(m),intvalup,'klo',klo,kup,'aflo',aflo,afup
             write(6,*)'zztop',zztop,zzbot
             write(6,*)'klo',klo,kup,'zupset',zupset,zupset,zloset
             write(6,*)'intaltairklo',klo,intaltair(klo:klo+1)
             write(6,*)'intaltairkup',kup,intaltair(kup:kup+1)
             write(6,*)m,'intcedsairklo',klo,intcedsair(i,j,klo:klo+1,m)
             write(6,*)m,'intcedsairkup',kup,intcedsair(i,j,kup:kup+1,m)
             write(6,*)m,'zup',zup,'zlo',zlo
             write(6,*)'slope',(intvalup-intvallo(m))/(zup-zlo)
             flush(6)
         call killit('bb error')
        endif
        sac=max(0.0,sac)
        cedsairemis(i,k,j,m)=sac ! ajl
        intvallo(m)=intvalup
      end do ! m
      aflo=afup
      klo=kup
      kbot=kup
      if(kup<0)then
         write(6,*)'kup error ',kup,i,j,k
         call killit('kup')
      endif
    end do DOK ! k
  end do !i
  return
  end subroutine docedsairm
  subroutine addcedsairemis(i,j,k,sacco,sacnox,sacn2o)
  use raqmschemlocaltype_mod
  use raqmschem_pmgrid_mod, only : tile
  implicit none
  real sacco,sacnox,sacn2o
  integer i,j,k,m
  do m=1,numcedsair
    select case (ccedsair(m))
    case ('CO')
      sacco=cedsairemis(i,k,j,m)
    case ('NOx')
      sacnox=cedsairemis(i,k,j,m)
    case ('N2O')
       sacn2o=cedsairemis(i,k,j,m)
    end select
  end do
  return
  end subroutine addcedsairemis
end module raqmschem_cedsair_mod
