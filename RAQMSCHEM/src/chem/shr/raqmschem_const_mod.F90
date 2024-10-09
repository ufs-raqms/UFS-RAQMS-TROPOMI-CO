module raqmschem_const_mod

  use chem_types_mod

  implicit none

  real,    parameter :: kappa=2./7.     
  real,    parameter :: p00=1.e5
  real,    parameter :: con_ttp=273.16
! try bigger may be too small
!  real,    parameter :: epsilc=1.e-30
!  real,    parameter :: epsilc2=1.e-34
  real,    parameter :: epsilc=1.e-22
  real,    parameter :: epsilc2=1.e-22
  REAL,    PARAMETER :: mw_so4_aer = 96.066
  real,    parameter :: mw_CO=28.01
  real,    parameter :: avgro=6.022E23
  real,    parameter :: cp=1004.6855     ! specific heat at const pres
  real,    parameter :: p1000=100000.       ! p at 1000mb (pascals)
  real*4,    parameter :: rd     = 2.8705e+2
  real*4,    parameter :: rdgsdchem     = 287.0586 
  real,    parameter :: rv     = 4.6150e+2
  real,    parameter :: eps=rd/rv
  real,    parameter :: epsqs=0.622
  real,    parameter :: grvity=9.80616
  real*4,    parameter :: grvityfms=9.80665
  public

end module raqmschem_const_mod
