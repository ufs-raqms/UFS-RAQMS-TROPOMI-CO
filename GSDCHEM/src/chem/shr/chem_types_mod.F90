module chem_types_mod

  use ESMF, only : ESMF_KIND_R8

  implicit none

  ! -- local types
  integer, parameter :: CHEM_KIND_I4 = kind(1)
  integer, parameter :: CHEM_KIND_R4 = kind(1.0)
  integer, parameter :: CHEM_KIND_R8 = kind(1.d0)
  integer, parameter :: CHEM_KIND_C4 = kind(cmplx(1.,1.,kind=CHEM_KIND_R4))
  integer, parameter :: CHEM_KIND_C8 = kind(cmplx(1.,1.,kind=CHEM_KIND_R8))

  ! -- ESMF-based coupling types
  integer, parameter :: CHEM_KIND_F8 = ESMF_KIND_R8

  ! -- internal limits
  integer, parameter :: CHEM_MAXSTR  = 256

  private

  public :: CHEM_KIND_I4, &
            CHEM_KIND_R4, &
            CHEM_KIND_R8, &
            CHEM_KIND_C4, &
            CHEM_KIND_C8, &
            CHEM_KIND_F8, &
            CHEM_MAXSTR

end module chem_types_mod
