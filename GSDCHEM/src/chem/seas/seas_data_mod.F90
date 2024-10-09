module seas_data_mod

  use chem_types_mod, only : fp => CHEM_KIND_R8

  ! -- parameters from NGAC v2.4.0 (rev. d48932c)
  integer,                             parameter :: number_ss_bins  = 5
  ! -- lower/upper particle radii (um) for each bin
  real,     dimension(number_ss_bins), parameter :: ra =             (/  0.03,   0.1,   0.5,   1.5,   5.0 /)
  real,     dimension(number_ss_bins), parameter :: rb =             (/   0.1,   0.5,   1.5,   5.0,  10.0 /)
  ! -- sea salt density
  real(fp), dimension(number_ss_bins), parameter :: den_seas  = (/    2200._fp,    2200._fp,    2200._fp,    2200._fp,    2200._fp /)
  ! -- particle effective radius (m)
  real(fp), dimension(number_ss_bins), parameter :: reff_seas = (/ 0.079e-6_fp, 0.316e-6_fp, 1.119e-6_fp, 2.818e-6_fp, 7.772e-6_fp /)

  ! -- default values for tuning parameters
  ! -- NGAC sea salt mass emission method: 1 = Gong 2003, 2 = Gong 1997, 3 = GEOS5 2012 (default)
  integer                             :: emission_scheme = 3
  ! -- global scaling factors for sea salt emissions (originally 0.875 in NGAC namelist)
  real,     dimension(number_ss_bins) :: emission_scale = (/ 0.100, 0.100, 0.100, 0.100, 0.100 /)

  public

end module seas_data_mod
