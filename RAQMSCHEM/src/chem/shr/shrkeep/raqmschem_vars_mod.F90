module raqmschem_vars_mod
  IMPLICIT NONE
  REAL, ALLOCATABLE :: chem( :, :, :, : )
  real,pointer     :: tr3d  (:,:,:) ! 1=pot.temp,2=water vapor,3=cloud water,4=ozone
end module raqmschem_vars_mod
