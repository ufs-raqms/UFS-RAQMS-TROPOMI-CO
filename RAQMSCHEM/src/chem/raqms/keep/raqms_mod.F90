module raqms_mod
  use chem_rc_mod
  use chem_types_mod,only : CHEM_KIND_R8,CHEM_KIND_R4
contains
  subroutine raqms_advance(de,dts,advancecount, &
    pr3d, &
    prl3d, &
    tk3d, &
    u3d, &
    v3d, &
    ws3d, &
    ph3d, &
    phl3d, &
    rn2d, &
    rc2d, &
    tr3d_in8, &
    tr3d_out, &
    ntra, &
    nbegin, &
    numgas, &
    num_moist, &
    nchem, &
    deg_lon,deg_lat, &
    yy,mm,dd,h,m,s,tz,juldayin, &
    is, ie, js, je, nb, nl, &
    iso, ieo, jso, jeo, nbo, nlo, &
    data, dz,&
    rc)
  use chem_comm_mod, only : chem_comm_get
  use raqmschem_data_mod, only : chem_data_type
  use field_manager_mod, only : find_field_index
  use raqmschem_species_mod
  use raqmschem_pmgrid_mod, only : initraqmschem_pmgrid,iam,masterproc,begj,endj,nc,nstepat,iprn,jprn,iprnin,jprnin
!  use raqmschemcomm_mod, only : allocatechemcommperm,allocatechemcommtemp,deallocatechemcommtemp
!  use raqmschemcomm_mod, only : deallocatechemcommperm
  use raqmschemcomm_mod
  use raqmschemlocaltype_mod
  use raqmschem_map_mod
  use wetdep_mod
  implicit none
  include 'comtim.H'
#include <comcdate.h>
  integer localrec,is,ie,js,je,nl,yy,mm,dd,h,m,ntra,numgas,nb,s,nbegin
  integer iso,ieo,jso,jeo,nbo,nlo,tz,npts,de,numi,numj,num_moist
  integer,save :: idate(4)
  integer advancecount,mype,iprint
  integer,optional :: rc
  real(CHEM_KIND_R8), intent(in) :: dts
  real(CHEM_KIND_R8), dimension(:, :, :), intent(in)  :: pr3d
  real(CHEM_KIND_R8), dimension(:, :, :), intent(in)  :: prl3d
  real(CHEM_KIND_R8), dimension(:, :, :), intent(in)  :: tk3d,u3d,v3d,ph3d,phl3d,ws3d
  real(CHEM_KIND_R8), dimension(:, :, :, :), intent(in)  :: tr3d_in8
  real(CHEM_KIND_R4), dimension(is:ie,js:je,nl,nchem) :: tr3d_inout4
  real(CHEM_KIND_R8), dimension(:,:), intent(in) :: rn2d,rc2d
  real(CHEM_KIND_R8), dimension(:, :, :, :), intent(out) :: tr3d_out
  real(CHEM_KIND_R8), dimension(:, :), intent(in) :: deg_lat
  real(CHEM_KIND_R8), dimension(:, :), intent(in) :: deg_lon
!  real(CHEM_KIND_R8), dimension(is:ie,js:je,nl) :: delp,diffozone
  real(CHEM_KIND_R8), dimension(is:ie,js:je,nl) :: diffozone
  integer jindx1(is:ie,js:je),jindx2(is:ie,js:je),i,j,k
  real(CHEM_KIND_R8), dimension(is:ie,js:je) :: ddy
  real(CHEM_KIND_R8) :: fhour
!  real(CHEM_KIND_R4),dimension(is:ie,js:je,nl) :: dens,difftrace
!  real(CHEM_KIND_R4),dimension(is:ie,js:je,nl) :: difftrace
!  real(CHEM_KIND_R4),dimension(is:ie,js:je) :: sorcpar
  real(CHEM_KIND_R4),dimension(is:ie,js:je) :: dz ! cm
  real(CHEM_KIND_R8),dimension(is:ie) :: draw2d
  real caldayfv
  logical  :: first=.true.
  logical  :: ldiag3d=.false.
  character *20 setall
  real*8 latat,lonat
  real dy_jdgmt
  real(CHEM_KIND_R8) :: chemutc
  type(chem_data_type) :: data
  integer ip,jp,mchem,nchem,juldayin
  type(chemlocaltype) :: chemlocal
  integer mbcurfv,mscurfv,mcdatefv,mcsecfv,localrc
  character *10 cwetdep
  logical dowetdep
  data dowetdep/.true./
  save dowetdep

  
!  if(first)then
!    if(mype.eq.0)then
!      write(6,*)'p_co',p_co,' p_no2 ',p_no2
!      call flush(6)
!      write(6,*)'raqmsustar',maxval(ustar),minval(ustar)
!     write(6,*)'raqms cmdrag',maxval(cmdrag),minval(cmdrag)
!     write(6,*)'raqms slmsk2d',maxval(slmsk2d),minval(slmsk2d)
!     endif
!  endif
  julday=juldayin
  nstepat=advancecount
  call chem_comm_get(localpe=mype)
  npts=(je-js+1)*(ie-is+1)
!  do k=1,nl
!    do j=js,je
!      jp=j-js+1
!      do i=is,ie
!        ip=i-is+1
!        delp(i,j,k)=pr3d(ip,jp,k)-pr3d(ip,jp,k+1)
!        dens(i,j,k)=7.2431122e+18*prl3d(ip,jp,k)/tk3d(ip,jp,k)
!      end do
!    end do
!  end do
!  numi=ie-is+1
!  numj=je-js+1
!  call initraqmschem_pmgrid(numi,numj,nl,js,je,12,mype)
!  call allocatechemcommperm
!  call allocatechemcommtemp
    idate(1)=h
    idate(2)=mm
    idate(3)=dd
    idate(4)=yy
!  write(6,*)'shape4',shape4
!  write(6,*)'ntra',ntra,'nchem',nchem
!  call flush(6)
  write(cdate,'(i4.4,3i2.2)')yy,mm,dd,h
  if(mype.eq.0)then
  write(6,*)'cdate',cdate
  endif
! first step might be after timestep
  fhour=float(advancecount+1)*dts/3600.
  if(mype.eq.0)then
  write(6,*)'fhour',fhour,'idate',idate
  endif
  nstep=advancecount+1
  dtime=dts
  if(first)then
    cwetdep=' '
    call getenv('WETDEP',cwetdep)
    if(cwetdep.ne.' ')then
      if(cwetdep.eq.'NO')then
        write(6,*)'turn off wetdep'
        call flush(6)
        dowetdep=.false.
      endif
    endif
    nnbdat=(mod(yy,100)*100+mm)*100+dd
    nnbsec=fhour*3600.

    mbdate=(yy*100+mm)*100+dd
!  write(6,*)'mbdate ',mbdate,' yy ',yy,mm,dd,'julday',julday
!  call flush(6)
    mbsec=0
  endif
! ajl try caldyi here to get date we need since ours is one dts greater than
! input
  call caldyi(nstep,dtime,nnbdat,nnbsec,mbdate,mbsec, &
    mbcurfv,mscurfv,mcdatefv,mcsecfv,caldayfv)
  if(mype.eq.0)then
      write(6,*)'mcdatefv',mcdatefv,'mcsecfv',mcsecfv,'caldayfv',caldayfv
      write(6,*)'hour',mcsecfv/3600.,'julday',julday
      julday=int(caldayfv)+1
      write(6,*)'julday chem',julday
  endif
  write(cdate,'(i8.8,i2.2)')mcdatefv,int(mcsecfv/3600.)
  if(mype.eq.0)then
     write(6,*)'cdatenew ',cdate
  endif
! we assumed that were dts later in raqms
!  chemutc=float(h)+float(m)/60.+float(s)/3600.
  chemutc=float(mcsecfv)/3600.
  gmt=chemutc
  dy_jdgmt=float(julday)+chemutc/24.0
  if(mype.eq.0)then
    write(6,*)'chemutc',chemutc,'dtime',dtime,'gmt',gmt,'dy_jdgmt',dy_jdgmt
  endif
  call chem_driver(julday)
!  write(6,*)'after chem_drive ',mbdate
!  call flush(6)
!  write(200+iam,*)'lbound prl3d',lbound(prl3d),ubound(prl3d),'js',js,je,'begj',begj,endj
!  call flush(200+iam)
!  write(6,*)'js',js,je
!  call flush(6)
!  write(6,*)'rn2d',maxval(rn2d),' rc2d ',maxval(rc2d)
!  call flush(6)
  tr3d_inout4=tr3d_in8
! can update in in4 in wetdep
! rn2d and rc2d are cumulative precip even though dictionary says  instantantaneous
! could do wetdep separate or in j loop and before chemistry or after chemistry
  if(dowetdep)then
    call wetdep_advance(advancecount,dts,rc2d,rn2d,ph3d,phl3d,pr3d,prl3d,tk3d,u3d,v3d,ws3d,tr3d_inout4, &
    nchem,num_moist,ntra,is, ie, js, je, 1, nl, &
    is, ie, js, je, 1, nl, &
    rc=localrc)
  endif
  if(iprn>0)then
    write(300+iam,*)'after wet',tr3d_inout4(iprnin,jprnin,1,p_bry)
    call flush(300+iam)
  endif
  do j=js,je
!    call maptochem(j,chemlocal,tr3d_in4,is,ie,nl,nchem)
!    write(6,*)'call maptochem ',j,'js',js,je
!    call flush(6)
!    write(200+iam,*)'call maptochem',j,'prl3d ',maxval(prl3d(:,j-js+1,:))
!    call flush(200+iam)
!    write(200+iam,*)'lbound',lbound(prl3d),ubound(prl3d)
!    call flush(200+iam)
    call maptochem(j,chemlocal,tr3d_inout4,nchem,pr3d,prl3d,tk3d,u3d,v3d,ph3d,phl3d)
!    call mapfromchem(j,chemlocal,tr3d_out,is,ie,nl,nchem)
!    write(6,*)j,'shape',shape(draw2d),'ie-is+1',ie-is+1,'is',is,'ie',ie
!    call flush(6)
!    write(6,*)'call chemdraw',j
!    !call flush(6)
!    call chemdraw(j,draw2d)
    call doloop(j,dy_jdgmt,chemlocal)
    call mapfromchem(j,chemlocal,tr3d_out,nchem)
    call deallocatechemlocal(chemlocal)
    call deallocatechemlocal2(chemlocal)
  end do
  call chem_diagnostics
  do mchem=1,ninput
    if(.not.ischemfull(mchem))then
      tr3d_out(:,:,:,icheminputpt(mchem))=tr3d_in8(:,:,:,icheminputpt(mchem))
    endif
  end do
  tr3d_out(:,:,:,1)=tr3d_in8(:,:,:,1)
!  do mchem=1,ninput
!    write(6,*)'trace ',mchem,icheminputpt(mchem),maxval(tr3d_out(:,:,:,icheminputpt(mchem)))
!  end do
!    difftrace=tr3d_out(:,:,:,icheminputpt(mchem))-tr3d_in(:,:,:,icheminputpt(mchem))
!    write(6,*)'diff trace ',mchem,icheminputpt(mchem),cheminput(mchem),maxval(difftrace),minval(difftrace)
  !end do
!  if(first)then
!    do mchem=1,ninput
!      tr3d_out(:,:,:,icheminputpt(mchem))=tr3d_in(:,:,:,icheminputpt(mchem))
!    end do
!    if(p_co>0)then
!       if(mype.eq.0)then
!         write(6,*)'lbound ',lbound(tr3d_out),lbound(data%srcco_totl),lbound(dens)
!         write(6,*)'co begin 11 20 ',tr3d_out(11,20,1,p_co),' inc ',data%srcco_totl(11,20)*dts/dens(11,20,1)/dz(11,20)
!       endif
!      tr3d_out(:,:,1,p_co)=tr3d_out(:,:,1,p_co)+data%srcco_totl(:,:)*dts/dens(:,:,1)/dz(:,:)
!      tr3d_out(:,:,1,p_co)=tr3d_out(:,:,1,p_co)+srcco_totl(:,:,mm)*dts/dens(:,:,1)/dz(:,:)
!      write(6,*)'lbound ',lbound(tr3d_out), 'src',lbound(srcco_totl)
!      call flush(6)
!    endif
!    if(mype.eq.0)then
!        write(6,*)'srcco',data%srcco_totl(11,20),'dts',dts,'dz',dz(11,20),'dens',dens(11,20,1)
!    endif
!    if(p_no2>0)then
!      tr3d_out(:,:,1,p_no2)=tr3d_out(:,:,1,p_no2)+data%srcnind(:,:)*dts/dens(:,:,1)/dz(:,:)
!      tr3d_out(:,:,1,p_no2)=tr3d_out(:,:,1,p_no2)+srcnind(:,:,mm)*dts/dens(:,:,1)/dz(:,:)
!    endif
!  else
    
!    do mchem=1,ninput
!      tr3d_out(:,:,:,icheminputpt(mchem))=tr3d_in(:,:,:,icheminputpt(mchem))
!    end do
!    if(p_co>0)then
!       if(mype.eq.0)then
!         write(6,*)'lbound ',lbound(tr3d_out),lbound(data%srcco_totl),lbound(dens)
!         write(6,*)'co begin ',tr3d_out(11,20,1,p_co),' inc ',data%srcco_totl(11,20)*dts/dens(11,20,1)/dz(11,20)
!       endif
!      tr3d_out(:,:,1,p_co)=tr3d_out(:,:,1,p_co)+data%srcco_totl(:,:)*dts/dens(:,:,1)/dz(:,:)
!      tr3d_out(:,:,1,p_co)=tr3d_out(:,:,1,p_co)+srcco_totl(:,:,mm)*dts/dens(:,:,1)/dz(:,:)
!    endif
!    if(p_no2>0)then

!      tr3d_out(:,:,1,p_no2)=tr3d_out(:,:,1,p_no2)+data%srcnind(:,:)*dts/dens(:,:,1)/dz(:,:)
      !tr3d_out(:,:,1,p_no2)=tr3d_out(:,:,1,p_no2)+srcnind(:,:,mm)*dts/dens(:,:,1)/dz(:,:)
!    endif
!  endif
!  if(p_par>0)then
!    sorcpar=(data%srcbutane(:,:)+5.*data%srcpentane(:,:)+6.*data%srchexane(:,:)+data%srcpropene(:,:))*dts/dens(:,:,1)
!    sorcpar=(srcbutane(:,:,mm)+5.*srcpentane(:,:,mm)+6.*srchexane(:,:,mm)+srcpropene(:,:,mm))*dts/dens(:,:,1)
!    tr3d_out(:,:,1,p_par)=tr3d_out(:,:,1,p_par)+sorcpar/dz(:,:)
!  endif
!  if(p_eth>0)then
!     tr3d_out(:,:,1,p_eth)=tr3d_out(:,:,1,p_eth)+srcethene(:,:,mm)*dts/dens(:,:,1)/dz(:,:)
!  endif
!  if(p_olet>0)then
!     tr3d_out(:,:,1,p_olet)=tr3d_out(:,:,1,p_olet)+srcpropene(:,:,mm)*dts/dens(:,:,1)/dz(:,:)
!     tr3d_out(:,:,1,p_olet)=tr3d_out(:,:,1,p_olet)+srcpropene(:,:,mm)*dts/dens(:,:,1)/dz(:,:)
!  endif 
!if(p_isop>0)then 
!  tr3d_out(:,:,1,p_isop)=tr3d_out(:,:,1,p_isop)+srcisop(:,:,mm)*dts/dens(:,:,1)/dz(:,:) endif
!  if(p_propar>0)then
!     tr3d_out(:,:,1,p_propar)=tr3d_out(:,:,1,p_propar)+srcpropane(:,:,mm)*dts/dens(:,:,1)/dz(:,:)
!  endif
    first=.false.
!  call deallocatechemcommperm
!  call deallocatechemcommtemp
  
  end subroutine raqms_advance
end module raqms_mod
