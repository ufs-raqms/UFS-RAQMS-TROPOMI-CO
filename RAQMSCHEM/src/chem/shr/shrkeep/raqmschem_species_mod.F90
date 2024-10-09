module raqmschem_species_mod

  implicit none

  type raqmschem_species_type
    sequence
    ! -- atmospheric tracers
    integer :: p_atm_shum    = 0
    integer :: p_atm_cldq    = 0
    integer :: p_atm_o3mr    = 0
    ! -- chemical tracers
    integer :: p_o3vm2       = 0
    integer :: p_qc          = 0 
    integer :: p_qi          = 0 
    integer :: p_qv          = 0 
  end type raqmschem_species_type

  private

  public :: raqmschem_species_type

end module raqmschem_species_mod
