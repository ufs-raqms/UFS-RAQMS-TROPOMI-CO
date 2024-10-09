module calcvort_dgrid_mod
use raqmschem_data_mod
use chem_types_mod
!use raqmschem_state_mod
!use raqmschem_config_mod
use raqmschem_model_mod
use chem_const_mod, only : omegx2,degrad,g=>grvity
use raqmschem_const_mod, only : kappa,p00
contains
  subroutine calcabsvort(vort,p3d,t3d,lat,ids,ide,jds,jde,nl,de,mype)
  use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
  implicit none
  integer ids,ide,jds,jde,nl,ihs,ihe,jhs,jhe,mype,de
  integer j,localrc,i,k,kb,kt
  real(CHEM_KIND_R8),dimension(ids:ide,jds:jde,nl),intent(in) :: p3d,t3d,vort
  real(CHEM_KIND_R8),dimension(ids:ide,jds:jde),intent(in) :: lat
  real(CHEM_KIND_R4),dimension(:,:,:),pointer :: absvort
  real(CHEM_KIND_R4) :: f
  real(CHEM_KIND_R4),dimension(ids:ide,jds:jde,nl) :: pottemp,dthdp
  type(chem_data_type),   pointer :: data
  type(raqmschem_config_type), pointer :: config
  call raqmschem_model_get(de=de, config=config, data=data, rc=localrc)
!  write(6,*)'p3d',maxval(p3d),'t3d',maxval(t3d)
!  if(.not.allocated(data%vort))then
!    allocate (data%vort(ids:ide,jds:jde,nl))
!  endif
  if(.not.allocated(data%absvort))then
    allocate (data%absvort(ids:ide,jds:jde,nl))
  endif
  if(.not.allocated(data%potvort))then
    allocate (data%potvort(ids:ide,jds:jde,nl))
  endif
!  vort=>data%vort
  absvort=>data%absvort
  do j=jds,jde
    do i=ids,ide
      f=omegx2*sin(lat(i,j)*degrad)
      do k=1,nl
        absvort(i,j,k)=vort(i,j,k)+f
      end do
    end do
  end do
  do k=1,nl
    do j=jds,jde
      do i=ids,ide
        pottemp(i,j,k)=t3d(i,j,k)*(p00/p3d(i,j,k))**kappa
      end do
    end do
!    write(6,*)'pottemp',k,maxval(pottemp(:,:,k)),minval(pottemp(:,:,k))
  end do
  do k=1,nl
    kb=max(k-1,1)
    kt=min(k+1,nl)
    do j=jds,jde
      do i=ids,ide
        dthdp(i,j,k)=g*(pottemp(i,j,kt)-pottemp(i,j,kb))/(p3d(i,j,kb)-p3d(i,j,kt))
        data%potvort(i,j,k)=absvort(i,j,k)*dthdp(i,j,k)
      end do
    end do
!    write(6,*)'dthdp',k,maxval(dthdp(:,:,k)),minval(dthdp(:,:,k))
!    write(6,*)'potvort',k,maxval(data%potvort(:,:,K)),minval(data%potvort(:,:,k))
  end do
  end subroutine calcabsvort
#ifdef OLDCALC
  subroutine calcvortdgrid(u3d,v3d,p3d,t3d,area,lat,ids,ide,jds,jde,nl,ihs,ihe,jhs,jhe,de,mype)
  use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
  use raqmschem_comm_mod, only : chem_reducetile_pushwithhalo
  implicit none
  integer ids,ide,jds,jde,nl,ihs,ihe,jhs,jhe,mype,de
  integer i,j,rc,k,ncu,nru,ncv,nrv
  real(CHEM_KIND_R8),dimension(ids:ide,jds:jde,nl),intent(in) :: u3d,v3d,p3d,t3d
  real(chEM_KIND_R8),dimension(ids:ide,jds:jde),intent(in) :: area
  real(CHEM_KIND_R4),dimension(ids:ide,jds:jde) :: rarea
!  real(CHEM_KIND_R4),dimension(ids:ide,jds:jde),intent(in) :: griddx4,griddy4
!  real(CHEM_KIND_R8),dimension(ids:ide,jds:jde) :: griddx,griddy
  real(CHEM_KIND_R4),dimension(ihs:ihe,jhs:jhe) :: dx,dy
  real(CHEM_KIND_R8),dimension(ids:ide,jds:jde),intent(in) :: lat
  real(CHEM_KIND_R8),dimension(ihs:ihe,jhs:jhe,nl) :: u3dearthwithhalo,v3dearthwithhalo
  real(CHEM_KIND_R8),dimension(ihs:ihe,jhs:jhe,nl) :: u3dwithhalo,v3dwithhalo
!  real(CHEM_KIND_R4),dimension(ihs:ihe,jhs:jhe) :: dxwithhalo,dywithhalo
  real(CHEM_KIND_R8),dimension(ids:ide,jhs:jhe+1) :: ud,dxd,utmp,ugrd
  real(CHEM_KIND_R8),dimension(ihs:ihe+1,jds:jde) :: vd,dyd,vtmp,vgrd
  real(CHEM_KIND_R8) :: first,second 
  real(CHEM_KIND_R4),dimension(:,:,:),pointer :: vort
  real(CHEM_KIND_R4),dimension(:,:,:),pointer :: absvort
  real(CHEM_KIND_R4) :: f
  real(CHEM_KIND_R4),dimension(ids:ide,jds:jde,nl) :: pottemp,dthdp
  type(chem_data_type),   pointer :: data
  integer imax,imin,jmax,jmin,jbar,ibar,localrc,kb,kt
  type(raqmschem_config_type), pointer :: config
  logical keepgrad
  
  call raqmschem_model_get(de=de, config=config, data=data, rc=localrc)
!  write(6,*)'p3d',maxval(p3d),'t3d',maxval(t3d)
  if(.not.allocated(data%vort))then
    allocate (data%vort(ids:ide,jds:jde,nl))
  endif
  if(.not.allocated(data%absvort))then
    allocate (data%absvort(ids:ide,jds:jde,nl))
  endif
  if(.not.allocated(data%potvort))then
    allocate (data%potvort(ids:ide,jds:jde,nl))
  endif
  vort=>data%vort
  absvort=>data%absvort
!  write(6,*)'vort bounds',lbound(vort),ubound(vort),'shape',shape(vort)
  keepgrad=.true.
  do j=jds,jde
    do i=ids,ide
      if(area(i,j)>0.0)then
        rarea(i,j)=1./area(i,j)
      else
         rarea(i,j)=0.0
      endif
    end do
  end do
!  griddx=griddx4
!  !griddy=griddy4
  call chem_reducetile_pushwithhalo(u3d,ids,ide,jds,jde,nl,u3dearthwithhalo, &
  ihs,ihe,jhs,jhe,de,rc)
  call chem_reducetile_pushwithhalo(v3d,ids,ide,jds,jde,nl,v3dearthwithhalo, &
  ihs,ihe,jhs,jhe,de,rc)
  !call chem_reducetile_pushwithhalo(data%griddx,ids,ide,jds,jde,dxwithhalo, &
  call chem_reducetile_pushwithhalo(data%griddx,ids,ide,jds,jde,dx, &
  ihs,ihe,jhs,jhe,de,rc)
!  call chem_reducetile_pushwithhalo(data%griddy,ids,ide,jds,jde,dywithhalo, &
  call chem_reducetile_pushwithhalo(data%griddy,ids,ide,jds,jde,dy, &
  ihs,ihe,jhs,jhe,de,rc)
! calculate vorticity
! determine where have 2dx and 2dy
!  dx=dxwithhalo
!  dy=dywithhalo
  do j=jhs+1,jhe
    do i=ids,ide
      dxd(i,j)=.5*(dx(i,j-1)+dx(i,j))
    end do
  end do
  if(jhs.eq.jds)then
    do i=ids,ide
       dxd(i,jhs)=dxd(i,jhs+1)
    end do
  endif
  if(jhe.eq.jde)then
    do i=ids,ide
       dxd(i,jde+1)=dxd(i,jde)
    end do
  endif
  do j=jds,jde
    do i=ihs+1,ihe
      dyd(i,j)=.5*(dy(i-1,j)+dy(i,j))
    end do
    if(ihs.eq.ids)then
       dyd(ihs,j)=dyd(ihs+1,j)
    endif
    if(ihe.eq.ide)then
      !dyd(ide+1,j)=1.5*dy(ide,j)-.5*dy(ide-1,j)
      dyd(ide+1,j)=dyd(ide,j)
    endif
  end do
  call uvtouvgrid(u3dearthwithhalo,v3dearthwithhalo,u3dwithhalo,v3dwithhalo,ihs,ihe,jhs,jhe,nl)
  do k=1,nl
!   convert from earth to grid relative

!   inside utmp
!    do j=jds,jde
!      do i=ids,ide
!        data%ugrid(i,j,k)=u3dwithhalo(i,j,k)
!        data%vgrid(i,j,k)=v3dwithhalo(i,j,k)
!      end do
    !end do
    do j=jhs+1,jhe
      do i=ids,ide
        ud(i,j)=.5*(u3dwithhalo(i,j-1,k)+u3dwithhalo(i,j,k))
        utmp(i,j)=ud(i,j)*dxd(i,j)
          
      end do
    end do
    if(jhs.eq.jds)then
      do i=ids,ide
!       now try ud on actual boundary with a half box
        ud(i,jhs)=u3dwithhalo(i,jhs,k)
        utmp(i,jhs)=ud(i,jhs)*dxd(i,jhs)

      end do
    endif
    if(jhe.eq.jde)then
      do i=ids,ide
        ud(i,jde+1)=1.5*u3dwithhalo(i,jde,k)-.5*u3dwithhalo(i,jde-1,k)
        utmp(i,jde+1)=ud(i,jde+1)*dxd(i,jde+1)
      end do
    endif
!   now vtmp
    do j=jds,jde
      do i=ihs+1,ihe
        vd(i,j)=.5*(v3dwithhalo(i-1,j,k)+v3dwithhalo(i,j,k))
        vtmp(i,j)=vd(i,j)*dyd(i,j)
      end do
      if(ihs.eq.ids)then
        vd(ihs,j)=1.5*v3dwithhalo(ihs,j,k)-.5*v3dwithhalo(ihs+1,j,k)
        vtmp(ihs,j)=vd(ihs,j)*dyd(ihs,j)
      endif
      if(ihe.eq.ide)then
        vd(ide+1,j)=1.5*v3dwithhalo(ide,j,k)-.5*v3dwithhalo(ide-1,j,k)
        vtmp(ide+1,j)=vd(ide+1,j)*dyd(ide+1,j)
      endif
      
    end do
!   redefine j=jds=jhs
    if(jhs.eq.jds)then
      do i=ihs+1,ihe
        vd(i,jhs)=.75*vd(i,jhs)+.26*vd(i,jhs+1)
        vtmp(i,jhs)=vd(i,jhs)*dyd(i,jds)*.5
      end do
    endif
    do j=jds,jde
      do i=ids,ide
        vort(i,j,k)=rarea(i,j)*(utmp(i,j)-utmp(i,j+1)-vtmp(i,j)+vtmp(i+1,j))
!        vort(i,j,k)=rarea(i,j)*(ugrd(i,j)-ugrd(i,j+1)-vgrd(i,j)+vgrd(i+1,j))
      end do
    end do
    if(jhs.eq.jds)then
      do i=ids,ide
!       since only half the areaa
        vort(i,jhs,k)=2.*vort(i,jhs,k)
      end do
    endif
!   now try justs extrapolation with 2nd order taylor series
    if(jhs.eq.jds)then
      do i=ids,ide
        first=vort(i,jds+2,k)-vort(i,jds+1,k)
        second=.5*(vort(i,jds+3,k)-2.*vort(i,jds+2,k)+vort(i,jds+1,k))
        vort(i,jds,k)=vort(i,jds+1,k)-first+second
      end do
    endif
    if(jhe.eq.jde)then
      do i=ids,ide
        first=vort(i,jde-1,k)-vort(i,jde-2,k)
        second=.5*(vort(i,jde-1,k)-2.*vort(i,jde-2,k)+vort(i,jde-3,k))
        vort(i,jde,k)=vort(i,jde-1,k)+first+second
      end do
    endif
!   now handle west edge
    if(ihs.eq.ids)then
      do j=jhs+1,jhe-1 ! skip corner for now
        first=vort(ids+2,j,k)-vort(ids+1,j,k)
        second=.5*(vort(ids+3,j,k)-2.*vort(ids+2,j,k)+vort(ids+1,j,k))
        vort(ids,j,k)=vort(ids+1,j,k)-first+second
      end do
      if(jhs.eq.jds)then
!       lower left corner
        vort(ids,jds,k)=(vort(ids+1,jds,k)+vort(ids+1,jds+1,k)+vort(ids,jds+1,k))/3.
      endif
      if(jhe.eq.jde)then
!       upper left corner
        vort(ids,jde,k)=(vort(ids+1,jde,k)+vort(ids+1,jde-1,k)+vort(ids,jde-1,k))/3.
      endif
    endif
    if(ihe.eq.ide)then
      do j=jhs+1,jhe-1 ! skip corner for now
        first=vort(ide-1,j,k)-vort(ide-2,j,k)
        second=.5*(vort(ide-1,j,k)-2.*vort(ide-2,j,k)+vort(ide-3,j,k))
        vort(ide,j,k)=vort(ide-1,j,k)+first+second
      end do
      if(jhs.eq.jds)then
!       loser left corner
        vort(ide,jds,k)=(vort(ide-1,jds,k)+ vort(ide-1,jds+1,k)+vort(ide,jds+1,k))/3.
      endif
      if(jhe.eq.jde)then
!       upper right corner
        vort(ide,jde,k)=(vort(ide-1,jde,k)+vort(ide-1,jde-1,k)+vort(ide,jde-1,k))/3.
      endif
    endif
  end do
!  vortout=vort
    do j=jds,jde
      do i=ids,ide
        f=omegx2*sin(lat(i,j)*degrad)
        do k=1,nl
          absvort(i,j,k)=vort(i,j,k)+f
        end do
      end do
    end do
  do k=1,nl
    do j=jds,jde
      do i=ids,ide
        pottemp(i,j,k)=t3d(i,j,k)*(p00/p3d(i,j,k))**kappa
      end do
    end do
!    write(6,*)'pottemp',k,maxval(pottemp(:,:,k)),minval(pottemp(:,:,k))
  end do
  do k=1,nl
    kb=max(k-1,1)
    kt=min(k+1,nl)
    do j=jds,jde
      do i=ids,ide
        dthdp(i,j,k)=g*(pottemp(i,j,kt)-pottemp(i,j,kb))/(p3d(i,j,kb)-p3d(i,j,kt))
        data%potvort(i,j,k)=absvort(i,j,k)*dthdp(i,j,k)
      end do
    end do
!    write(6,*)'dthdp',k,maxval(dthdp(:,:,k)),minval(dthdp(:,:,k))
!    write(6,*)'potvort',k,maxval(data%potvort(:,:,K)),minval(data%potvort(:,:,k))
  end do

  return
  end subroutine calcvortdgrid


subroutine uvtouvgrid(u,v,ugrd,vgrd,ihs,ihe,jhs,jhe,nl)
use raqmschem_pmgrid_mod,only : tile,iam
use raqmschemcomm_mod,only : ucosa,usina,cosb,sinb
!real u(nc,nr),v(nc,nr),ugrd(nc,nr),vgrd(nc,nr)
integer ihs,ihe,jhs,jhe,k,nl
real u(ihs:ihe,jhs:jhe,nl),v(ihs:ihe,jhs:jhe,nl),ugrd(ihs:ihe,jhs:jhe,nl),vgrd(ihs:ihe,jhs:jhe,nl)
!real cosa(nc,nr),sina(nc,nr),cosb(nc,nr),sinb(nc,nr)
write(6,*)'uvtouvgrid ',tile,ihs,ihe,jhs,jhe,nl
write(200+iam,*)'uvtouvgrid tile ',tile,ihs,ihe,jhs,jhe,nl
call flush(200+iam)
call flush(6)
select case(tile)
  case(1,2,6)
    if(tile==6)then
       write(200+iam,*)'ihs',ihs,ihe,jhs,jhe,'nl',nl
       call flush(6)
       do j=jhs,jhe
         do i=ihs,ihe
           write(200+iam,'(2i3," cosa ",2f9.2," cosb ",2f9.2)')i,j,ucosa(i,j),usina(i,j),cosb(i,j),sinb(i,j)
         end do
      end do
    endif
    do k=1,nl
    do j=jhs,jhe
      do i=ihs,ihe
        ugrd(i,j,k)=u(i,j,k)*ucosa(i,j)-v(i,j,k)*usina(i,j)
        vgrd(i,j,k)=u(i,j,k)*cosb(i,j)+v(i,j,k)*sinb(i,j)
    if(tile==6.and.k.eq.nl)then
      write(200+iam,'(2i3," u ",2f9.2," ugrd ",2f9.2)')i,j,u(i,j,k),v(i,j,k),ugrd(i,j,k),vgrd(i,j,k)
    endif
      end do
    end do
    end do
  caSe(3,4,5)
    do k=1,nl
    do j=jhs,jhe
      do i=ihs,ihe
        ugrd(i,j,k)=u(i,j,k)*ucosa(i,j)+v(i,j,k)*usina(i,j)
        vgrd(i,j,k)=u(i,j,k)*cosb(i,j)+v(i,j,k)*sinb(i,j)
      end do
    end do
    end do
  case default
    write(6,*)'illegal case'
    call flush(6)
end select

end subroutine uvtouvgrid
#endif
end module calcvort_dgrid_mod
