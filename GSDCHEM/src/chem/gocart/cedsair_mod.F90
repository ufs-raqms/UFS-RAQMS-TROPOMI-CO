module cedsair_mod
use chem_rc_mod
use chem_types_mod
use_chem_tracers_mod,only : p_bc_1,p_bc_2,p_oc_1,p_oc_2,p_so2
implicit none
integer(chem_kind_i4),dimension(:,:,:),allocatable :: ktopcedsair
real(chem_kind_r4),dimension(:,:,:),allocatable :: intaltair
real(chem_kind_r4),dimension(:,:,:,:),allocatable :: intcedsair
contains
  subroutine cedsair_driver(dt,chem,dz8w,z_at_w,rri, &
  ims,ime,jms,jme,kms,kme,num_chem)
  real(chem_kind_r4),intent(in) :: dt
  real(chem_kind_r4),dimension(ims:ime,,kms:kme,jms:jme,num_chem),intent(inout) :: chem
  real(chem_kind_r4),dimension(ims:ime,kms:kme,jms:jme),intent(in) :: z_at_w,dz8w
  integer,intent(in) :: ims,ime,jms,jme,kms,kme,num_chem
  real(chem_kind_r4) :: factor
  real(chem_kind_r4),dimension(kms:kme) ::afup,aflo
  integer i,j,k
  integer(chem_kind_i4),dimension(kms:kme) ::  ktop,kbot
  factor=dt*rri(i,k,j)/dz8w(i,k,j)
  do j=jms,jme
    do i=ims,ime
      call cedsaircoeff(z_at_w(i,:,j),dz82(i,:,j),afup,aflo,kbot,ktop)
      do k=kms,kme
      end do
    end do
  end do
  end subroutine cedsair_driver
end module cedsair_mod
