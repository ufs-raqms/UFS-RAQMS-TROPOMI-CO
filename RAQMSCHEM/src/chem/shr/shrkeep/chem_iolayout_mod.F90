module chem_iolayout_mod

  use mpi

  implicit none

  type chem_iolayout_type
    integer :: mpicomm     = MPI_COMM_NULL
    integer :: modelComm   = MPI_COMM_NULL
    integer :: tileComm    = MPI_COMM_NULL
    logical :: localIOflag = .false.
  end type chem_iolayout_type

  private

  public :: chem_iolayout_type

end module chem_iolayout_mod
