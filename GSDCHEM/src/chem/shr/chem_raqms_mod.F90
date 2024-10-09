module chem_raqms_mod
use chem_state_mod, only : chem_state_type
type(chem_state_type),pointer :: chem_pass_state
!real*4,allocatable :: chem_pass(:,:,:,:)
end module chem_raqms_mod
