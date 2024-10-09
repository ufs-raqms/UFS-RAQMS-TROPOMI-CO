module dust_data_mod

  use chem_types_mod,   only : CHEM_KIND_R8
  use chem_tracers_mod, only : p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5, &
                               p_edust1, p_edust2, p_edust3, p_edust4, p_edust5


  implicit none

  integer, parameter :: ndust = 5
  integer, parameter :: ndcls = 3
  integer, parameter :: ndsrc = 1
  integer, parameter :: maxstypes = 100
  integer, parameter :: nsalt = 9

  real,    parameter :: dyn_visc = 1.5E-5

  ! -- dust parameters
  integer,            dimension(ndust), parameter :: ipoint    = (/ 3, 2, 2, 2, 2 /)
  real(CHEM_KIND_R8), dimension(ndust), parameter :: den_dust  = (/   2500.,  2650.,  2650.,  2650.,  2650. /)
  real(CHEM_KIND_R8), dimension(ndust), parameter :: reff_dust = (/ 0.73D-6, 1.4D-6, 2.4D-6, 4.5D-6, 8.0D-6 /)
  real(CHEM_KIND_R8), dimension(ndust), parameter :: frac_s    = (/     0.1,   0.25,   0.25,   0.25,   0.25 /)
  real(CHEM_KIND_R8), dimension(ndust), parameter :: lo_dust   = (/  0.1D-6, 1.0D-6, 1.8D-6, 3.0D-6, 6.0D-6 /)
  real(CHEM_KIND_R8), dimension(ndust), parameter :: up_dust   = (/  1.0D-6, 1.8D-6, 3.0D-6, 6.0D-6,10.0D-6 /)
  real(CHEM_KIND_R8), dimension(ndust, 12)        :: ch_dust   = 0.8e-09_CHEM_KIND_R8

  ! -- default dust parameters
  ! -- AFWA & GOCART
  ! -----------+----------+-----------+
  ! Parameter  | FIM-Chem | HRRR-Chem |
  ! -----------+----------+-----------+
  ! alpha      |      1.0 |       0.5 |
  ! gamma      |      1.6 |       1.0 |
  ! -----------+----------+-----------+
  real, parameter :: afwa_alpha    = 0.2
  real, parameter :: afwa_gamma    = 1.3
  real, parameter :: gocart_alpha  = 0.3
  real, parameter :: gocart_gamma  = 1.3
  ! -- FENGSHA
  real, parameter :: fengsha_alpha = 0.3
  real, parameter :: fengsha_gamma = 1.3
  ! -- FENGSHA threshold velocities based on Dale A. Gillette's data
  integer, parameter :: fengsha_maxstypes = 13
  real, dimension(fengsha_maxstypes) :: dust_uthres = &
    (/ 0.065,   & ! Sand            - 1
       0.20,    & ! Loamy Sand      - 2
       0.52,    & ! Sandy Loam      - 3
       0.50,    & ! Silt Loam       - 4
       0.50,    & ! Silt            - 5
       0.60,    & ! Loam            - 6
       0.73,    & ! Sandy Clay Loam - 7
       0.73,    & ! Silty Clay Loam - 8
       0.80,    & ! Clay Loam       - 9
       0.95,    & ! Sandy Clay      - 10
       0.95,    & ! Silty Clay      - 11
       1.00,    & ! Clay            - 12
       9.999 /)   ! Other           - 13
  ! -- FENGSHA uses precalculated drag partition from ASCAT. See: Prigent et al. (2012,2015)
  integer :: dust_calcdrag = 1

  ! -- values set at initialization
  real :: dust_alpha = 0.
  real :: dust_gamma = 0.


  ! -- sea salt parameters
  integer,            dimension(nsalt), parameter :: spoint    = (/ 1, 2, 2, 2, 2, 2, 3, 3, 3 /)  ! 1 Clay, 2 Silt, 3 Sand
  real(CHEM_KIND_R8), dimension(nsalt), parameter :: reff_salt = &
    (/ 0.71D-6, 1.37D-6, 2.63D-6, 5.00D-6, 9.50D-6, 18.1D-6, 34.5D-6, 65.5D-6, 125.D-6 /)
  real(CHEM_KIND_R8), dimension(nsalt), parameter :: den_salt  = &
    (/   2500.,   2650.,   2650.,   2650.,   2650.,   2650.,   2650.,   2650.,   2650. /)
  real(CHEM_KIND_R8), dimension(nsalt), parameter :: frac_salt = &
    (/      1.,     0.2,     0.2,     0.2,     0.2,     0.2,   0.333,   0.333,   0.333 /)


  ! -- soil vagatation parameters
  integer, parameter :: max_soiltyp = 30
  real, dimension(max_soiltyp) :: &
    maxsmc = (/ 0.421, 0.464, 0.468, 0.434, 0.406, 0.465, &
                0.404, 0.439, 0.421, 0.000, 0.000, 0.000, &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)

  ! -- other soil parameters
  real, dimension(maxstypes) :: porosity

  public

end module dust_data_mod
