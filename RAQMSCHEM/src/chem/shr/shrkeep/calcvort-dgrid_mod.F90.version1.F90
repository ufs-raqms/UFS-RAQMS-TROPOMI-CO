module calcvort_dgrid_mod
contains
  subroutine calcvortdgrid(u3d,v3d,area,griddx4,griddy4,lat,ids,ide,jds,jde,nl,ihs,ihe,jhs,jhe,vort,de,mype)
  use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
  use chem_io_mod, only : chem_reducetile_pushwithhalo_3dr8,chem_reducetile_pushwithhalo_2dr8
  implicit none
  integer ids,ide,jds,jde,nl,ihs,ihe,jhs,jhe,mype,de
  integer i,j,rc,k
  real(CHEM_KIND_R8),dimension(ids:ide,jds:jde,nl),intent(in) :: u3d,v3d
  real(chEM_KIND_R8),dimension(ids:ide,jds:jde),intent(in) :: area
  real(CHEM_KIND_R4),dimension(ids:ide,jds:jde) :: rarea
  real(CHEM_KIND_R4),dimension(ids:ide,jds:jde),intent(in) :: griddx4,griddy4
  real(CHEM_KIND_R8),dimension(ids:ide,jds:jde) :: griddx,griddy
  real(CHEM_KIND_R4),dimension(ihs:ihe,jhs:jhe) :: dx,dy
  real(CHEM_KIND_R8),dimension(jds:jde),intent(in) :: lat
  real(CHEM_KIND_R4),allocatable,dimension(:,:,:),intent(out) :: vort
  real(CHEM_KIND_R8),dimension(ihs:ihe,jhs:jhe,nl) :: u3dwithhalo,v3dwithhalo
  real(CHEM_KIND_R8),dimension(ihs:ihe,jhs:jhe) :: dxwithhalo,dywithhalo
  real(CHEM_KIND_R8),dimension(ids:ide,jhs:jhe+1) :: ud,dxd,utmp
  real(CHEM_KIND_R8),dimension(ihs:ihe+1,jds:jde) :: vd,dyd,vtmp
  integer imax,imin,jmax,jmin,jbar,ibar
  logical keepgrad
  keepgrad=.true.
  if(mype.eq.1)then
    write(6,*)'u3d',maxval(u3d),minval(u3d)
    write(6,*)'v3d',maxval(v3d),minval(v3d)
    write(6,*)'area',maxval(area),minval(area)
    write(6,*)'griddx',maxval(griddx),minval(griddx)
    write(6,*)'griddy',maxval(griddy),minval(griddy)
    write(6,*)'lat',maxval(lat),minval(lat)
    write(6,*)'ids',ids,ide,'jds',jds,jde
    call flush(6)
  endif
  do j=jds,jde
    do i=ids,ide
      if(area(i,j)>0.0)then
        rarea(i,j)=1./area(i,j)
      else
         write(6,*)'area zero',i,j,area(i,j)
         call flush(6)
         rarea(i,j)=0.0
      endif
!     need to be half till do on average box instead
      rarea(i,j)=rarea(i,j)
#if 1
      if(j.eq.jbar.and.abs(i-ibar)<3)then
        if(mype<=1)then
          write(6,*)'rarea',i,j,rarea(i,j)
          write(80+mype,*)'rarea',i,j,rarea(i,j)
        endif
      endif
#endif
    end do
  end do
!  write(6,*)'rarea',maxval(rarea),minval(rarea)
  if(.not.allocated(vort))then
    allocate(vort(ids:ide,jds:jde,nl))
  endif
  griddx=griddx4
  griddy=griddy4
  call chem_reducetile_pushwithhalo_3dr8(u3d,ids,ide,jds,jde,nl,u3dwithhalo, &
  ihs,ihe,jhs,jhe,de,rc)
  call chem_reducetile_pushwithhalo_3dr8(v3d,ids,ide,jds,jde,nl,v3dwithhalo, &
  ihs,ihe,jhs,jhe,de,rc)
  call chem_reducetile_pushwithhalo_2dr8(griddx,ids,ide,jds,jde,dxwithhalo, &
  ihs,ihe,jhs,jhe,de,rc)
  call chem_reducetile_pushwithhalo_2dr8(griddy,ids,ide,jds,jde,dywithhalo, &
  ihs,ihe,jhs,jhe,de,rc)
  write(6,*)'after reduce'
  call flush(6)
#if 0
  if(mype<=1)then
    write(6,*)'shape v3dwithhalo',shape(v3dwithhalo)
    do i=ihs,ihe
      write(6,*)'calcvort vtmp63 ',i,v3dwithhalo(i,48,63),'u',u3dwithhalo(i,48,65)
    end do
  endif
  if(mype.eq.1)then
  
    write(6,*)'u3dwithalo',maxval(u3dwithhalo),minval(u3dwithhalo)
    write(6,*)'v3dwithalo',maxval(u3dwithhalo),minval(u3dwithhalo)
  endif
#endif
!  vort=1000.
! calculate vorticity
! determine where have 2dx and 2dy
  dx=dxwithhalo
  dy=dywithhalo
#if 0
  imin=min(ihs+1,ids)
  jmin=min(jhs+1,jds)
  imax=max(ihe-1,ide)
  jmax=max(jhe-1,jde)
  if(mype.eq.1)then
     write(6,*)'imin',imin,imax,'jmin',jmin,jmax,'ids',ids,ide,'j',jds,jde
  endif
! move to calc
!  do j=jmin,jmax
!    do i=imin,imax
!      dx(i,j)=.5*dx(i,j)
!      dy(i,j)=.5*dy(i,j)
!    end do
!  end do
  if(mype.eq.1)then
      write(6,*)'dx',maxval(dx),minval(dx)
      write(6,*)'dy',maxval(dy),minval(dy)
  endif
#endif
  if(mype<=1)then
    write(80+mype,*)'top vort ids',ids,ide,' jds ',jds,jde
    write(80+mype,*)'top vort ihs',ihs,ihe,' jhs ',jhs,jhe
      call flush(80+mype)
  endif
  do j=jhs+1,jde
    do i=ids,ide
      dxd(i,j)=.5*(dx(i,j-1)+dx(i,j))
    end do
  end do
  if(jhs.eq.jds)then
    do i=ids,ide
!      dxd(i,jhs)=1.5*dx(i,jhs)-.5*dx(i,jhs+1)
       dxd(i,jhs)=dxd(i,jhs+1)
    end do
  endif
  if(jhe.eq.jde)then
    do i=ids,ide
!      dxd(i,jde+1)=1.5*dx(i,jde)-.5*dx(i,jde-1)
       dxd(i,jde+1)=dxd(i,jde)
    end do
  endif
  do j=jds,jde
    do i=ihs+1,ihe
      dyd(i,j)=.5*(dy(i-1,j)+dy(i,j))
    end do
    if(ihs.eq.ids)then
!      dyd(ihs,j)=1.5*dy(ihs,j)-.5*dy(ihs+1,j)
       dyd(ihs,j)=dyd(ihs+1,j)
    endif
    if(ihe.eq.ide)then
      !dyd(ide+1,j)=1.5*dy(ide,j)-.5*dy(ide-1,j)
      dyd(ide+1,j)=dyd(ide,j)
    endif
  end do
  do k=1,nl
!   inside utmp
    do j=jhs+1,jde
      do i=ids,ide
        ud(i,j)=.5*(u3dwithhalo(i,j-1,k)+u3dwithhalo(i,j,k))
        utmp(i,j)=ud(i,j)*dxd(i,j)
        if(k.eq.63.and.i.eq.40.and.mype<=1)then
          write(70+mype,*)i,j,'ud',ud(i,j),'dxd',dxd(i,j)
          write(70+mype,*)'u3with',i,j-1,u3dwithhalo(i,j-1,k),i,u3dwithhalo(i,j,k)
          write(70+mype,*)'utmp',j,utmp(i,j)
          call flush(70+mype)
        endif
          
      end do
    end do
    if(jhs.eq.jds)then
      do i=ids,ide
        ud(i,jhs)=1.5*u3dwithhalo(i,jhs,k)-.5*u3dwithhalo(i,jhs+1,k)
!        utmp(i,jhs)=ud(i,jhs)*dxd(i,jhs)
!       keep ud gradient
        ud(i,jhs)=2.*ud(i,jhs+1)-ud(i,jhs+2)
        utmp(i,jhs)=ud(i,jhs)*dxd(i,jhs)

        if(k.eq.63.and.i.eq.40.and.mype<=1)then
          write(70+mype,*)i,'jhs',jhs,'ud',ud(i,jhs),'dxd',dxd(i,jhs)
          write(70+mype,*)'u3withalo',jhs,u3dwithhalo(i,jhs,k),jhs+1,u3dwithhalo(i,jhs+1,k)
          write(70+mype,*)'utmp',ihs,utmp(i,jhs)
          call flush(70+mype)
        endif
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
        if(k.eq.63.and.j.eq.10.and.mype<=1)then
          write(70+mype,*)i,' cen vd ',vd(i,j),'dyd',dyd(i,j)
          call flush(70+mype)
        endif
      end do
      if(ihs.eq.ids)then
        vd(ihs,j)=1.5*v3dwithhalo(ihs,j,k)-.5*v3dwithhalo(ihs+1,j,k)
        vtmp(ihs,j)=vd(ihs,j)*dyd(ihs,j)
        if(k.eq.63.and.j.eq.10.and.mype<=1)then
          write(70+mype,*)i,' west vd ',vd(ihs,j),'dyd',dyd(ihs,j)
          call flush(70+mype)
        endif
      endif
      if(ihe.eq.ide)then
        vd(ide+1,j)=1.5*v3dwithhalo(ide,j,k)-.5*v3dwithhalo(ide-1,j,k)
        vtmp(ide+1,j)=vd(ide+1,j)*dyd(ide+1,j)
        if(k.eq.63.and.j.eq.10.and.mype<=1)then
           write(70+mype,*)i,' east vd ',vd(ide+1,j),'dyd',dyd(ide+1,j)
          call flush(70+mype)
        endif 
      endif
    end do
    do j=jds,jde
      do i=ids,ide
        vort(i,j,k)=rarea(i,j)*(utmp(i,j)-utmp(i,j+1)-vtmp(i,j)+vtmp(i+1,j))
!        if(mype<=1.and.abs(j-48)<=5.and.k.eq.63.and.i.eq.96)then
        if(mype<=1.and.i.eq.40.and.k.eq.63)then
          write(80+mype,*)'vort ',vort(i,j,k),'utmp',j,utmp(i,j),j+1,utmp(i,j+1)
          write(80+mype,*)'vtmp',i,j,vtmp(i,j),i+1,j,vtmp(i+1,j),'rarea',rarea(i,j)
          write(80+mype,*)'dutmp',utmp(i,j)-utmp(i,j+1)
          write(80+mype,*)'dvtmp',-vtmp(i,j)+vtmp(i+1,j)
      call flush(80+mype)
        endif
#if 0
        if(j.eq.jbar.and.abs(i-ibar)<3.and.k.eq.63)then
          if(mype<=1)then
            write(6,*)'vort',i,vort(i,j,k),'vtmp',i-1,vtmp(i-1,j),i+1,vtmp(i+1,j)
            write(6,*)'utmp',i,j-1,utmp(i,j-1),utmp(i,j+1)
            write(80+mype,*)'vort',i,vort(i,j,k),'vtmp',i-1,vtmp(i-1,j),i+1,vtmp(i+1,j)
            write(80+mype,*)'utmp',i,j-1,utmp(i,j-1),utmp(i,j+1)
          endif
        endif
#endif
      end do
    end do
  end do
  return
  end subroutine calcvortdgrid
end module calcvort_dgrid_mod
