module raqmschem_config_mod

  use chem_rc_mod
  use chem_types_mod,   only : CHEM_MAXSTR, CHEM_KIND_R4
  use chem_comm_mod,    only : chem_comm_bcast, chem_comm_isroot
  use raqmschem_species_mod, only : raqmschem_species_type
!  use machine, only : kind_phys

  implicit none

  ! -- Currently available modules
   integer, parameter :: CHEM_OPT_NONE         = 0
   integer, parameter :: CHEM_OPT_MAX          = 500

  character(len=*), parameter :: chem_file_nml = 'input.nml'

  ! -- data structure for configuration options
  type raqmschem_config_type
    sequence
    character(len=CHEM_MAXSTR) :: emi_inname         = ''
    character(len=CHEM_MAXSTR) :: emi_outname        = ''
    character(len=CHEM_MAXSTR) :: input_chem_inname  = ''
    character(len=CHEM_MAXSTR) :: input_chem_outname = ''
    character(len=CHEM_MAXSTR) :: chem_hist_outname  = 'raqms_out_'
    logical :: readrestart        = .false.
    integer :: archive_step
    ! -- control variables
    integer :: ntra                = 3      ! # of tracers advected on small dt: FV3 only has 3 tracers
    integer :: ntrb                = 0      ! # of tracers advected on large dt: will include chemistry
    integer :: num_chem            = 0
    integer :: num_moist           = 0
    integer :: numgas              = 0
    integer :: nbegin              = 0
    integer :: raqmschem_opt       = 0
    real(CHEM_KIND_R4) :: raqmschemdt
    ! -- parameters
  ! -- configuration variables for output:
  !  . grid level defined in FIM Makefile
    integer :: glvl            = 0
  ! -- control variables
  type(raqmschem_species_type), pointer :: species => null()


  end type raqmschem_config_type

  private

  public :: raqmschem_config_type
  ! -- provide subtypes

  public :: raqmschem_config_read
  public :: chem_config_control_init
  public :: chem_config_species_init
  public :: CHEM_OPT_NONE,         &
            CHEM_OPT_MAX

contains

  subroutine raqmschem_config_read(config, rc)

    type(raqmschem_config_type), intent(inout) :: config
    integer, optional,      intent(out)   :: rc

    ! -- local variables
    integer, parameter :: unit = 200

    integer                :: localrc, iostat
    integer                :: raqmschem_opt
    integer                :: buffer(1)
    real(CHEM_KIND_R4)     :: raqmschemdt
    real(CHEM_KIND_R4)     :: rbuffer(1)
    character(CHEM_MAXSTR) :: sbuffer(3)

    ! -- variables in input namelist
    character(len=CHEM_MAXSTR) :: emi_inname
    character(len=CHEM_MAXSTR) :: fireemi_inname
    character(len=CHEM_MAXSTR) :: emi_outname
    character(len=CHEM_MAXSTR) :: fireemi_outname
    character(len=CHEM_MAXSTR) :: input_chem_inname
    character(len=CHEM_MAXSTR) :: input_chem_outname
    character(len=CHEM_MAXSTR) :: chem_hist_outname


    namelist /raqmschem_nml/          &
      emi_inname,                &
      fireemi_inname,            &
      emi_outname,               &
      fireemi_outname,           &
      input_chem_inname,         &
      input_chem_outname,        &
      chem_hist_outname         

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    ! -- set defaults
    emi_inname         = ""
    fireemi_inname     = ""
    emi_outname        = ""
    fireemi_outname    = ""
    input_chem_inname  = ""
    input_chem_outname = ""
    chem_hist_outname  = "chem_out_"

    ! -- read chem configuration namelist

    if (chem_comm_isroot()) then
      open(unit, file=chem_file_nml, form='formatted', status='old', iostat=iostat)
    end if
    call chem_comm_bcast(iostat, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    if (chem_rc_test((iostat /= 0), msg="Failed to open namelist file: "//chem_file_nml, &
        file=__FILE__, line=__LINE__, rc=rc)) return
    if (chem_comm_isroot()) then
      rewind(unit)
      read(unit, nml=raqmschem_nml, iostat=iostat )
      close(unit)
      write(6, nml=raqmschem_nml)
    end if
    call chem_comm_bcast(iostat, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    if (chem_rc_test((iostat /= 0), msg="Failed to read &chem namelist in: "//chem_file_nml, &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- pack integer values in buffer
    buffer = (/ &
      raqmschem_opt          &
      /)
    ! -- broadcast integer buffer
    call chem_comm_bcast(buffer, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    ! -- load integer values in local config
    config % raqmschem_opt          = buffer( 1 )

    ! -- pack real variables in buffer
    rbuffer = (/ raqmschemdt/)
    ! -- broadcast real buffer
    call chem_comm_bcast(rbuffer, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    ! -- local real values in local buffer
    config % raqmschemdt     = rbuffer(1)

    ! -- pack strings into buffer
    sbuffer = (/ chem_hist_outname, emi_inname, emi_outname /)
    ! -- broadcast string variable
    call chem_comm_bcast(sbuffer, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    ! -- set string values to config
    config % chem_hist_outname = sbuffer(1)
    config % emi_inname        = sbuffer(2)
    config % emi_outname       = sbuffer(3)
!    write(6,*)'ajl bottom raqmschem_config_read'
!    call flush(6)

  end subroutine raqmschem_config_read

  subroutine chem_config_control_init(config, rc)

    type(raqmschem_config_type), intent(inout) :: config
    integer, optional,      intent(out)   :: rc

    ! -- local variables

    ! -- begin
!    write(6,*)'ajl top chemconfigcontrolinit'
!    call flush(6)
    if (present(rc)) rc = CHEM_RC_SUCCESS

    config % ntrb = 0

    config % numgas     = 1

          config % num_chem            = 1
          config % num_moist           = 3
          ! compute total # of tracers - no ice variable transported
          config % ntrb = config % ntrb + config % num_moist + config % num_chem - 3

    ! -- nbegin is the start address (-1) of the first chem variable in tr3d
    if (config % num_moist > 3) then
      config % nbegin = config % ntra + config % num_moist - 2
    else
      config % nbegin = config % ntra + config % num_moist - 3
    end if

  end subroutine chem_config_control_init


  subroutine chem_config_species_init(config, rc)

    type(raqmschem_config_type), intent(inout) :: config
    integer, optional,      intent(out)   :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin
!    write(6,*)'ajl top  config_species_init'
!    call flush(6)
    if (.not.associated(config % species)) then
      allocate(config % species, stat=localrc)
      if (chem_rc_test((localrc /= 0), msg="Failed to allocate species data container", &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if

    ! -- set pointers to predefined atmospheric tracers
    ! -- NOTE: this is model-dependent
    config % species % p_atm_shum = 1
    config % species % p_atm_cldq = 2
    config % species % p_atm_o3mr = 3


        if (chem_rc_test((config % num_chem    /= 1), &
          msg="num_chem is not equal to 1", &
          file=__FILE__, line=__LINE__, rc=rc)) return
        config % species % p_qv=1
        config % species % p_qc=2
        config % species % p_qi=3
        config % species % p_o3vm2=1
        config % numgas=4
!    write(6,*)'ajl bot  config_species_init'
!    call flush(6)
  end subroutine chem_config_species_init






end module raqmschem_config_mod
