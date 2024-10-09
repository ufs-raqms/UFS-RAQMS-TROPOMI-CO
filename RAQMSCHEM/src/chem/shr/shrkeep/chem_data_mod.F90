module chem_data_mod
  use chem_rc_mod
  use chem_types_mod

  implicit none

  type chem_data_type
    ! -- input
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: tr3d
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcco_totl
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcnind
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcethane
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcpropane
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcbutane
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcpentane
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srchexane
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcethene
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcpropene
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcch2o
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcald
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcalkanone
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcisop
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcterp
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: ch4clim
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcbc
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcoc
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcnair
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcbcair
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcocair
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srccoair
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: srcso2air
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: griddx
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: griddy
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: griddxwithhalo
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: griddywithhalo
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: vort
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: absvort
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: potvort

  end type chem_data_type

  private

  public :: chem_data_type
  public :: chem_data_destroy
contains
  subroutine chem_data_destroy(data, rc) 
    type(chem_data_type)           :: data
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS
    if (allocated(data % tr3d)) then
      deallocate(data % tr3d, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % srcco_totl)) then
      deallocate(data % srcco_totl, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % srcnind)) then
      deallocate(data % srcnind, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % griddx)) then
      deallocate(data % griddx, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % griddy)) then
      deallocate(data % griddy, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % griddxwithhalo)) then
      deallocate(data % griddxwithhalo, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % griddywithhalo)) then
      deallocate(data % griddywithhalo, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % vort)) then
      deallocate(data % vort, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % absvort)) then
      deallocate(data % absvort, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % potvort)) then
      deallocate(data % potvort, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
  end subroutine chem_data_destroy

end module chem_data_mod
