module chem_types_mod

  implicit none

#ifdef __LEGACY__
  integer, parameter :: CHEM_KIND_I4  = 4
  integer, parameter :: CHEM_KIND_I8  = 8
  integer, parameter :: CHEM_KIND_R4  = 4
  integer, parameter :: CHEM_KIND_R8  = 8
  integer, parameter :: CHEM_KIND_C8  = 8
  integer, parameter :: CHEM_KIND_C16 = 16
#else
  integer, parameter :: CHEM_KIND_I4  = selected_int_kind(9)
  integer, parameter :: CHEM_KIND_I8  = selected_int_kind(18)
  integer, parameter :: CHEM_KIND_R4  = selected_real_kind(3,25)
  integer, parameter :: CHEM_KIND_R8  = selected_real_kind(6,45)
  integer, parameter :: CHEM_KIND_C8  = selected_real_kind(3,25)
  integer, parameter :: CHEM_KIND_C16 = selected_real_kind(6,45)
#endif
  integer, parameter :: CHEM_MAXSTR  = 256

  public

end module chem_types_mod
