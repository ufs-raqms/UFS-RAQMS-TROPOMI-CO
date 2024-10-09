module chem_data_mod

  use chem_rc_mod
  use chem_types_mod

  implicit none

  type chem_data_type
    ! -- input
    real(CHEM_KIND_R4), dimension(:),     allocatable :: p_gocart          ! GOCART pressure levels
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: clayfrac          ! clay fraction (AFWA & FENGSHA dust scheme)
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: dm0               ! dms reference emissions
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: emiss_ab          ! emissions for all available species
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: emiss_abu         ! emissions for all available species
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_ash_dt      ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_ash_height  ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_ash_mass    ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_tr_dt       ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_tr_height   ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_tr_mass     ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: ero1              ! dust erosion factor
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: ero2              ! dust erosion factor
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: ero3              ! dust erosion factor
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: rdrag             ! Drag Partition Map (FENGSHA)
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: uthr              ! Threshold Friction Velocity (FENGSHA)
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: ssm               ! Sediment Supply Map (FENGSHA)
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: h2o2_backgd       ! H2O2 background for GOCART
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: no3_backgd        ! NO3 background for GOCART
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: oh_backgd         ! OH background for GOCART
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: plume             ! fire info - MODIS & GBBEPx
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: plume_hr             ! fire info - MODIS & GBBEPx
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: emiss_abu_hr         ! emissions for all available species
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: sandfrac          ! sand fraction (AFWA & FENGSHA dust scheme)
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: th_pvsrf
    ! -- output
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: pm10
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: pm25
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: ebu_oc
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: oh_bg
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: h2o2_bg
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: no3_bg
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: wet_dep
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: ext_cof
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: sscal
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: asymp
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: tr3d
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: trdp
    ! -- internal buffers
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: aod2d
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: eburn
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

    if (allocated(data % p_gocart)) then
      deallocate(data % p_gocart, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % clayfrac)) then
      deallocate(data % clayfrac, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % dm0)) then
      deallocate(data % dm0, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % emiss_ab)) then
      deallocate(data % emiss_ab, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % emiss_abu)) then
      deallocate(data % emiss_abu, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % emiss_ash_dt)) then
      deallocate(data % emiss_ash_dt, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % emiss_ash_height)) then
      deallocate(data % emiss_ash_height, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % emiss_ash_mass)) then
      deallocate(data % emiss_ash_mass, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % emiss_tr_dt)) then
      deallocate(data % emiss_tr_dt, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % emiss_tr_height)) then
      deallocate(data % emiss_tr_height, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % emiss_tr_mass)) then
      deallocate(data % emiss_tr_mass, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % ero1)) then
      deallocate(data % ero1, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % ero2)) then
      deallocate(data % ero2, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % ero3)) then
      deallocate(data % ero3, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % rdrag)) then
      deallocate(data % rdrag, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % uthr)) then
      deallocate(data % uthr, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % ssm)) then
      deallocate(data % ssm, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % h2o2_backgd)) then
      deallocate(data % h2o2_backgd, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % no3_backgd)) then
      deallocate(data % no3_backgd, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % oh_backgd)) then
      deallocate(data % oh_backgd, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % plume)) then
      deallocate(data % plume, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % sandfrac)) then
      deallocate(data % sandfrac, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % th_pvsrf)) then
      deallocate(data % th_pvsrf, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % pm10)) then
      deallocate(data % pm10, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % pm25)) then
      deallocate(data % pm25, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % ebu_oc)) then
      deallocate(data % ebu_oc, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % oh_bg)) then
      deallocate(data % oh_bg, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % h2o2_bg)) then
      deallocate(data % h2o2_bg, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % no3_bg)) then
      deallocate(data % no3_bg, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % wet_dep)) then
      deallocate(data % wet_dep, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % ext_cof)) then
      deallocate(data % ext_cof, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % sscal)) then
      deallocate(data % sscal, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % asymp)) then
      deallocate(data % asymp, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % tr3d)) then
      deallocate(data % tr3d, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % trdp)) then
      deallocate(data % trdp, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % aod2d)) then
      deallocate(data % aod2d, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % eburn)) then
      deallocate(data % eburn, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if

  end subroutine chem_data_destroy

end module chem_data_mod
