subroutine calcaodgsi(tr3d,dzaod,data,nchem,ids,ide,jds,jde,nl,pri3d,prl3d,tk3d,de)
  use computeaod
  use raqmschem_data_mod
  use raqmschem_model_mod
  use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
  use raqmschem_species_mod, only : idaod,nsolmax,idaodgsi
  use raqmschem_const_mod, only : epsqs,rd
  use chem_const_mod, only : grvity
  use raqmschem_pmgrid_mod, only : iam,tile
  use funcphys, only : fpvs,gfuncphys
  implicit none
  integer ids,ide,jds,jde,nl,ihs,ihe,jhs,jhe,nchem,n
!  real(CHEM_KIND_R4),intent(in) :: tr3d(ids:ide,jds:jde,nl,nchem)
  real(CHEM_KIND_R8),intent(in) :: tr3d(ids:ide,jds:jde,nl,nchem)
  real(CHEM_KIND_R8),intent(in) :: pri3d(ids:ide,jds:jde,nl+1)
  real(CHEM_KIND_R8),intent(in) :: prl3d(ids:ide,jds:jde,nl)
  real(CHEM_KIND_R4),intent(in) :: dzaod(ids:ide,jds:jde,nl)
  real(CHEM_KIND_R8),intent(in) :: tk3d(ids:ide,jds:jde,nl)
  real(CHEM_KIND_R4) :: aero(ids:ide,jds:jde,nl,nsolmax),sum
  real(CHEM_KIND_R4) :: es1,rh(nl,ids:ide,jds:jde),qs
  real(chem_kind_R4) :: ugkg_kgm2(ids:ide,jds:jde,nl)
  integer :: irh(ids:ide,jds:jde,nl),ii,ncid,ierr
  real    :: drh(ids:ide,jds:jde,nl),daod,dext,kecrtm,kecrtmmax(14)
  real(CHEM_KIND_R4),parameter :: threshold=1.e-6*18./28.97
  integer i,j,k,m,nc,nr,de,localrc,maxrh
  real(chem_kind_r4),allocatable,dimension(:,:) :: kecrtmrh
  type(chem_data_type),   pointer :: data
  character *8 caod(15)
  data caod/'sulf','bc1','bc2','oc1','oc2','dust1','dust2','dust3','dust4','dust5', &
  'seas1','seas2','seas3','seas4','seas5'/
  integer iprn,jprn,tileprn
  integer,parameter :: maxrhdim=540
  iprn=43
  jprn=51
  tileprn=3
  call raqmschem_model_get(de=de, data=data, rc=localrc)
  nc=ide-ids+1
  nr=jde-jds+1
!  do n=1,nsolmax
  do n=1,11 ! up to seas1
    if(idaodgsi(n)>0)then
      aero(:,:,:,n)=tr3d(:,:,:,idaodgsi(n))
    else
      aero(:,:,:,n)=0.0
    endif
  end do
  if(idaodgsi(12)>0)then ! conbine seas1 and seas2 both fine
    aero(:,:,:,11)=aero(:,:,:,11)+tr3d(:,:,:,idaodgsi(12))
  endif
  do n=12,14
    aero(:,:,:,n)=tr3d(:,:,:,idaodgsi(n+1))
  end do
! need rh so need temperature and density
  if(.not.allocated(kecrtmrh))then
    allocate(kecrtmrh(0:maxrhdim,14))
  endif
  call readkecrtmrh(maxrh,maxrhdim,kecrtmrh)
  do k=1,nl
    do j=jds,jde
      do i=ids,ide
        es1=fpvs(tk3d(i,j,k)) ! pa
        es1=min(es1,prl3d(i,j,k))
!       here qs is specific humidity
        qs=epsqs*es1/(prl3d(i,j,k)-(1.-epsqs)*es1)
        rh(k,i,j)=max(threshold,tr3d(i,j,k,1))/qs
        rh(k,i,j)=max(0.,min(.99,rh(k,i,j)))
        if(rh(k,i,j)<=.80)then
          irh(i,j,k)=int(rh(k,i,j)*200.)
          drh(i,j,k)=rh(k,i,j)*200.-float(irh(i,j,k))
        elseif(rh(k,i,j)>=.99)then
          irh(i,j,k)=maxrh
        else
          irh(i,j,k)=int((rh(k,i,j)-.80)*2000.)
          drh(i,j,k)=(rh(k,i,j)-.80)*2000.-float(irh(i,j,k))
        endif
        ugkg_kgm2(i,j,k)=(pri3d(i,j,k)-pri3d(i,j,k+1))*1.e-9/grvity
      end do
    end do
  end do
  if(.not.allocated(data%aod))then
    allocate (data%aod(nc,jds:jde))
  endif
  if(.not.allocated(data%aodg5))then
    allocate (data%aodg5(nc,jds:jde))
  endif
  if(.not.allocated(data%ext_3d))then
    allocate (data%ext_3d(nc,jds:jde,nl,nsolmax))
  endif
  if(.not.allocated(data%ext_3dg))then
    allocate (data%ext_3dg(nc,jds:jde,nl,nsolmax))
  endif
 data%aodg5=0.0
 data%ext_3dg=0.0
! kecrtmmax=0.0
!  write(6,*)'kecrtmrh',lbound(kecrtmrh),ubound(kecrtmrh),'maxrh',maxrh
!  flush(6)
  do k=1,nl
    do j=jds,jde
      do i=ids,ide
        ii=i-ids+1
        do n=1,14
          if(irh(i,j,k)>=maxrh)then
            kecrtm=kecrtmrh(maxrh,n)
          else
            kecrtm=kecrtmrh(irh(i,j,k),n)*(1.-drh(i,j,k))+kecrtmrh(irh(i,j,k)+1,n)*drh(i,j,k)
          endif
!          kecrtmmax(n)=max(kecrtmmax(n),kecrtm)
          daod=kecrtm*aero(i,j,k,n)*ugkg_kgm2(i,j,k)
          data%aodg5(ii,j)=data%aodg5(ii,j)+daod
          dext=daod/dzaod(i,j,k)
          data%ext_3dg(ii,j,k,n)=data%ext_3dg(ii,j,k,n)+dext
        end do
      end do
    end do
  end do
          
end subroutine calcaodgsi
subroutine readkecrtmrh(maxrh,maxrhdim,kecrtmrh)
use netcdf
 use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
 use raqmschem_pmgrid_mod, only : iam
implicit none
integer ncid,maxrh,dimid,idke,ierr,n,k,maxrhdim
real(chem_kind_r4):: kecrtmrh(0:maxrhdim,14)
character *256 file
file='/lfs4/BMC/rcm3/UFS_FV3_GSD_RAQMS_CHEM/GSI/ke.gocart.geos5.nc'
ierr=nf90_open(file,0,ncid)
if(ierr/=nf90_noerr)then
   write(6,*)'error open ',trim(file)
   write(6,*)trim(nf90_strerror(ierr))
endif
ierr=nf90_inq_dimid(ncid,'nrh',dimid)
if(ierr/=nf90_noerr)then
   write(6,*)'error inq dimid nrh '
   write(6,*)trim(nf90_strerror(ierr))
endif
ierr=nf90_inquire_dimension(ncid,dimid,len=maxrh)
if(ierr/=nf90_noerr)then
   write(6,*)'error getdim marrh '
   write(6,*)trim(nf90_strerror(ierr))
endif
maxrh=maxrh-1
ierr=nf90_inq_varid(ncid,'KECRTM',idke)
if(ierr/=nf90_noerr)then
   write(6,*)'error inq varid kecrtm'
   write(6,*)trim(nf90_strerror(ierr))
endif
ierr=nf90_get_var(ncid,idke,kecrtmrh)
if(ierr/=nf90_noerr)then
   write(6,*)'error get var kecrtmrh'
   write(6,*)trim(nf90_strerror(ierr))
endif
ierr=nf90_close(ncid)
return
end subroutine readkecrtmrh

