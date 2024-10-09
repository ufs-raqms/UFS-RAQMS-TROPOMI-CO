module vertmx_driver_mod
contains
  subroutine vertmx_driver(dt,tr3d,num_chem,exch_h,dryrho,z_at_w,zmid,ids,ide,jds,jde,kms,kme,kts,kte)
  use chem_types_mod
  use raqmschem_pmgrid_mod, only : iam,tile
  use raqmschem_species_mod,only : p_ox,p_co,lraqmschem
  use raqmschemcomm_mod, only : covermx,oxvermx
  use raqmschem_const_mod, only : epsilc2
  use vertmxv_mod, only : vertmxv
  implicit none
  
  integer ids,ide,jds,jde,i,j,k,numvermx,num_chem,m,kts,kte,kms,kme
  real(CHEM_KIND_R4) :: dt
  real(CHEM_KIND_R4) ,dimension(ids:ide,kms:kme,jds:jde) :: zmid,dryrho,exch_h
  real(CHEM_KIND_R4), dimension(ids:ide,jds:jde,kts:kte,num_chem) :: tr3d
  real(CHEM_KIND_R4), dimension(ids:ide,kms:kme,jds:jde) :: z_at_w ! is kte+1
  real, dimension(kts:kte+1) :: zzfull,ekmfull
  real, dimension(kts:kte) :: zz,dryrhopt
  real vd
  real,dimension(kts:kte,num_chem) :: pblst
  integer numvertmx,ipoint(num_chem),ii
  vd=0.
  numvermx=0
  do j=jds,jde
    do i=ids,ide
      ii=i-ids+1
      do k=kts,kte
        zz(k)=zmid(i,k,j)-z_at_w(i,kts,j)
        zzfull(k)=z_at_w(i,k,j)-z_at_w(i,kts,j)
        dryrhopt(k)=dryrho(i,k,j)
      end do
      zzfull(kte+1)=z_at_w(i,kte+1,j)-z_at_w(i,kts,j)
      ekmfull(kts)=0.
      do k=kts+1,kte
        ekmfull(k)=exch_h(i,k,j)
      end do
      ekmfull(kte+1)=0.0
      numvermx=0
      do m=1,num_chem
        if(lraqmschem(m))then
          numvermx=numvermx+1
          pblst(:,numvermx)=max(epsilc2,tr3d(i,j,:,m))
          ipoint(numvermx)=m
        endif
      end do
!      call vertmxv(dt,numvermx,pblst,ekmfull,dryrho,zzfull,zz,vd,1,kte-kts+1)
!     fix 2/15/2021 ajl
      call vertmxv(dt,numvermx,pblst,ekmfull,dryrhopt,zzfull,zz,vd,1,kte-kts+1)
      do m=1,numvermx
        if(ipoint(m).eq.p_co)then
          covermx(ii,j,:)=max(epsilc2,pblst(:,m))-tr3d(i,j,:,p_co)
        endif
        if(ipoint(m).eq.p_ox)then
          oxvermx(ii,j,:)=max(epsilc2,pblst(:,m))-tr3d(i,j,:,p_ox)
        endif
        tr3d(i,j,:,ipoint(m))=max(epsilc2,pblst(:,m))
      end do
    end do
  end do

  end subroutine vertmx_driver
end module vertmx_driver_mod
