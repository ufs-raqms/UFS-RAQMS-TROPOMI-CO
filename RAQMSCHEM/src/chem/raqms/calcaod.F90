subroutine calcaod(tr3d,dzaod,data,nchem,ids,ide,jds,jde,nl,prl3d,tk3d,de)
  use computeaod
  use raqmschem_data_mod
  use raqmschem_model_mod
  use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
  use raqmschem_species_mod, only : idaod,nsolmax
  use raqmschem_const_mod, only : epsqs,rd
  use raqmschem_pmgrid_mod, only : iam,tile
  use funcphys, only : fpvs,gfuncphys
  implicit none
  integer ids,ide,jds,jde,nl,ihs,ihe,jhs,jhe,nchem,n
!  real(CHEM_KIND_R4),intent(in) :: tr3d(ids:ide,jds:jde,nl,nchem)
  real(CHEM_KIND_R8),intent(in) :: tr3d(ids:ide,jds:jde,nl,nchem)
  real(CHEM_KIND_R8),intent(in) :: prl3d(ids:ide,jds:jde,nl),tk3d(ids:ide,jds:jde,nl)
  real(CHEM_KIND_R4),intent(in) :: dzaod(ids:ide,jds:jde,nl)
  real(CHEM_KIND_R4) :: aero(ids:ide,jds:jde,nl,nsolmax),sum
  real(CHEM_KIND_R4) :: den(ids:ide,jds:jde,nl),es1,rh(nl,ids:ide,jds:jde),qs
  real(CHEM_KIND_R4),parameter :: threshold=1.e-6*18./28.97
  integer i,j,k,m,nc,nr,de,localrc
  type(chem_data_type),   pointer :: data
  integer iprn,jprn,tileprn
  iprn=43
  jprn=51
  tileprn=3
!  write(6,*)'top calcaod nchem ',nchem,'ids',ids,ide,jds,jde,'nl',nl,'de',de
!  call flush(6)
  call raqmschem_model_get(de=de, data=data, rc=localrc)
  nc=ide-ids+1
  nr=jde-jds+1
!  do n=1,nsolmax
  do n=1,11 ! include up to sea1
    if(idaod(n)>0)then
      aero(:,:,:,n)=tr3d(:,:,:,idaod(n))
!      write(200+iam,*)'aero',n,maxval(aero(:,:,:,n)),minval(aero(:,:,:,n))
      
!      write(6,*)'aero',n,maxval(aero(:,:,:,n)),minval(aero(:,:,:,n))
    else
      aero(:,:,:,n)=0.0
    endif
  end do
  if(idaod(12)>0)then
    aero(:,:,:,11)=aero(:,:,:,11)+tr3d(:,:,:,idaod(12))
  endif
  do n=12,14
    if(idaod(n+1)>0)then
      aero(:,:,:,n)=tr3d(:,:,:,idaod(n+1))
    else
      aero(:,:,:,n)=0.0
    endif
  end do
  call flush(200+iam)
! need rh so need temperature and density
  do k=1,nl
    do j=jds,jde
      do i=ids,ide
        den(i,j,k)=prl3d(i,j,k)/(rd*tk3d(i,j,k))
        es1=fpvs(tk3d(i,j,k)) ! pa
        es1=min(es1,prl3d(i,j,k))
        qs=epsqs*es1/(prl3d(i,j,k)-(1.-epsqs)*es1)
        rh(k,i,j)=100.*max(threshold,tr3d(i,j,k,1))/qs
!       if(i.eq.ids.and.j.eq.jds)then
!       write(400+iam,*)'i',i,j,'k',k,'rh',rh(k,i,j),'q',tr3d(i,j,k,1),'qs',qs
!       write(400+iam,*)'es1',es1,'prl3d',prl3d(i,j,k),'den',den(i,j,k)
!       write(400+iam,*)'tk3d',tk3d(i,j,k)
!       call flush(400+iam)
!        endif
        rh(k,i,j)=max(0.,min(99.,rh(k,i,j)))
      end do
    end do
  end do
!  write(200+iam,*)'den',maxval(den)
!  write(200+iam,*)'rh',maxval(rh),minval(rh)
!  write(200+iam,*)'dzaod',maxval(dzaod),minval(dzaod)
!  call flush(200+iam)
  if(.not.allocated(data%aod))then
    allocate (data%aod(nc,jds:jde))
  endif
  if(.not.allocated(data%ext_3d))then
!    write(600+iam,*)'allocate data%ext_3d de',de
!    flush(600+iam)
    allocate (data%ext_3d(nc,jds:jde,nl,nsolmax))
  endif
!  write(6,*)'call compute_aod',nl,nc,nsolmax,kind(data%aod),kind(data%ext_3d)
!  write(6,*)'size',size(data%aod),size(data%ext_3d)
!  call flush(6)
!  write(300+iam,*)'dzaod',kind(dzaod),'lb',lbound(dzaod),' ub ',ubound(dzaod)
!  write(300+iam,*)'nl',nl,'nc',nc,'nsolmax',nsolmax
!  call flush(300+iam)
  call compute_aod_raqms(nl,nc,nsolmax,den,rh,dzaod,aero,data%aod,nl,data%ext_3d)
  if(iam.eq.0)then
    write(6,*)'aod ',maxval(data%aod)
  endif
!  if(tile.eq.tileprn)then
!    if(ids<=iprn.and.ide>=iprn.and.jds<=jprn.and.jde>=jprn)then
!      do n=1,nsolmax
!        if(n.eq.6)cycle
        !sum=0.0
!        do k=1,nl
!          write(6,*)'ext',n,k,data%ext_3d(iprn,jprn,k,n)*dzaod(iprn,jprn,k)
!          sum=sum+data%ext_3d(iprn,jprn,k,n)*dzaod(iprn,jprn,k)
!        end do
!        write(6,*)'sum',n,sum
      !end do
!    endif
!  endif
!  write(6,*)'did aod'
!  call flush(6)
!  write(200+iam,*)'aod',maxval(data%aod),minval(data%aod)
!  call flush(6)
!  do n=1,nsolmax
!    do k=1,nl
!      write(200+iam,*)'ext_3d',n,k,maxval(data%ext_3d(:,:,k,n))
!      call flush(200+iam)
    !end do
!   end do
!  write(6,*)'den',maxval(den),minval(den)
!  write(6,*)'rh',maxval(rh),minval(rh)
end subroutine calcaod
