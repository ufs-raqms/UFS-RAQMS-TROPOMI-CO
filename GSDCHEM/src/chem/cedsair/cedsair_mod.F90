module cedsair_mod
  use  chem_types_mod
  real(chem_kind_r4),allocatable :: cedsairemis(:,:,:,:)
contains
!  subroutine cedsair_driver(dt,rri,chem,dz,zatw,num_chem, &
  subroutine cedsair_driver(dt,rri,dz,zatw, &
  is,ie,js,je,ks,ke,kte,tile)
  use  chem_types_mod
  use cedsair_data_mod,only : intcedsair,intaltair,ktopcedsair,kbotcedsair,num_cedsair 
  use cedsair_data_mod,only : cem_cedsair,num_cedsair
  use cedsair_data_mod,only : p_cedsair_bc,p_cedsair_oc,p_cedsair_so2
  use chem_tracers_mod, only : p_bc1,p_bc2,p_oc1,p_oc2,p_so2
  use chem_raqms_mod, only : cedsairinc
  implicit none
  integer hi,is,ie,js,je,i,j,k,n,kup,klo,kk,tile,kte
  integer ks,ke,m,ktop,kbot,kkb
  real(chem_kind_r4) :: dt,zlo,zup,sac,zztop,zzbot
  real(chem_kind_r8) :: aflo,afup,intvallo(num_cedsair),intvalup
  real(chem_kind_r4),dimension(is:ie,ks:ke,js:je) :: rri,dz,zatw
!  real(chem_kind_r4),dimension(is:ie,ks:ke,js:je,num_chem) :: chem
  logical zloset,zupset

! bc1,2 and oc1,oc2 are in ug/kg
! rho is in Kg/m2
! dz is in m
! best if intairbc, intairoc are in ug/m2/sec
! then dbc or doc in ug/Kg is delta(intairceds)/rho/dz*dt = delta(intairceds)*rri/dz*dt
! so2 is in ppmv
! best if intairso2 is in umole/m2/sec /mwair
! then dso2 in ppmv is delta(intairceds)*rri/dz/*dt/mwair
  if(.not.allocated(cedsairemis))then
    allocate(cedsairemis(is:ie,ks:kte,js:je,3))
  endif
  if(.not.allocated(cedsairinc))then
    allocate(cedsairinc(is:ie,js:je,ks:kte,3))
  endif
  cedsairemis=0.0
  cedsairinc=0.0
  do j=js,je
    do i=is,ie
      if(ktopcedsair(i,j)<1)cycle
      zztop=intaltair(ktopcedsair(i,j))
      zzbot=intaltair(kbotcedsair(i,j)-1)
      kbot=max(kbotcedsair(i,j)-1,0)
      ktop=ktopcedsair(i,j)
      zloset=.false. ! can reuse since zlo=zup-1
      klo=-1
      aflo=0.
      if(kbot==0)then
        intvallo=0.0
        kkb=1
        zloset=.true.
        aflo=0.0
      else
        do k=1,kte
          zup=zatw(i,k+1,j)
          if(zup<=zzbot)cycle
          zlo=zatw(i,k,j)
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
DOK:  do k=kkb,kte-1
        zupset=.false.
        kup=-1
        afup=0.
        zlo=zatw(i,k,j)
        if(zlo>=zztop)exit DOK
        zup=zatw(i,k+1,j)
        if(zup<=zzbot)cycle
        if(zup>=zztop-.000001)then
          zupset=.true.
          afup=1.
          kup=ktop
        endif
        if(.not.zupset)then
        do kk=kbot,ktop
          if(.not.zupset)then
            if(zup<intaltair(kk+1).and.zup>=intaltair(kk))then
              kup=kk
              afup=(zup-intaltair(kk))/(intaltair(kk+1)-intaltair(kk))
              zupset=.true.
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
           do kk=kbot-1,ktop
             write(6,*)i,j,k,'intaltair',kk,intaltair(kk)
             flush(6)
           end do
!           call killit('error')
        endif
        do m=1,num_cedsair
          intvalup=intcedsair(i,j,kup+1,m)*afup+intcedsair(i,j,kup,m)*(1.-afup)
!         here zup, zlo are in meters
          sac=(intvalup-intvallo(m))/(zup-zlo)*rri(i,k,j)*dt
          sac=max(0.0,sac)
          cedsairemis(i,k,j,m)=sac ! ajl
!          if(isnan(cedsairemis(i,k,j,m)))then
!             write(6,*)'cedsairemis nan',i,k,j,m
!             flush(6)
!             call killit('ceds nan')
!          endif
          cedsairinc(i,j,k,m)=sac*1.e-6 ! make ppv
             intvallo(m)=intvalup
          end do ! m
          aflo=afup
          klo=kup ! for next level since zlo=zup-1
          kbot=kup
      end do DOK
    end do
  end do
  end subroutine cedsair_driver
  subroutine addcedsairemis(chem,num_chem,is,ie,js,je,ks,ke,kte,tile)
  use cedsair_data_mod,only : cem_cedsair,num_cedsair
  use cedsair_data_mod,only : p_cedsair_bc,p_cedsair_oc,p_cedsair_so2
  use chem_tracers_mod, only : p_bc1,p_bc2,p_oc1,p_oc2,p_so2
  implicit none
  integer is,ie,js,je,ks,ke,num_chem,i,j,k,m,tile,kte
  real(chem_kind_r4),dimension(is:ie,ks:ke,js:je,num_chem) :: chem
  do j=js,je
    do k=ks,kte
      do i=is,ie
        do m=1,num_cedsair
          select case (cem_cedsair(m))
          case ('BC')
            chem(i,k,j,p_bc1)=chem(i,k,j,p_bc1)+.8*cedsairemis(i,k,j,p_cedsair_bc)
            chem(i,k,j,p_bc2)=chem(i,k,j,p_bc2)+.2*cedsairemis(i,k,j,p_cedsair_bc)
          case ('OC')
            chem(i,k,j,p_oc1)=chem(i,k,j,p_oc1)+.5*cedsairemis(i,k,j,p_cedsair_oc)
            chem(i,k,j,p_oc2)=chem(i,k,j,p_oc2)+.5*cedsairemis(i,k,j,p_cedsair_oc)
          case ('SO2')
            chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+cedsairemis(i,k,j,p_cedsair_so2)
          case default
          end select
        end do
      end do
    end do
  end do
  return
  end subroutine addcedsairemis
end module cedsair_mod
