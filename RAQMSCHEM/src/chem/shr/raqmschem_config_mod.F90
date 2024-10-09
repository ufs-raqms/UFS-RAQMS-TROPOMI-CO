module raqmschem_config_mod

  use chem_rc_mod
  use chem_types_mod,   only : CHEM_MAXSTR, CHEM_KIND_R4
  use chem_comm_mod,    only : chem_comm_bcast, chem_comm_isroot

  implicit none

  ! -- Currently available modules
   integer, parameter :: CHEM_OPT_NONE         = 0
   integer, parameter :: CHEM_OPT_MAX          = 500
!  add for plume rise
  ! -- biomass burning emissions
  integer, parameter :: BURN_OPT_NONE   = 0  
  integer, parameter :: BURN_OPT_ENABLE = 1
  integer, parameter :: FIRE_OPT_NONE   = 0  
  integer, parameter :: FIRE_OPT_MODIS  = 1  
  integer, parameter :: FIRE_OPT_GBBEPx = 2
  ! -- subgrid convective transport
  integer, parameter :: CTRA_OPT_NONE  = 0  
  integer, parameter :: CTRA_OPT_GRELL = 2
  ! -- large scale wet deposition
  integer, parameter :: WDLS_OPT_NONE  = 0  
  integer, parameter :: WDLS_OPT_GSD   = 1  
  integer, parameter :: WDLS_OPT_NGAC  = 2  
  integer, parameter :: WDLS_OPT_NGAC_BOTH  = 3  


  character(len=*), parameter :: chem_file_nml = 'input.nml'
!  logical :: dofirehour=.false.

  ! -- data structure for configuration options
  type raqmschem_config_type
    sequence
    character(len=CHEM_MAXSTR) :: emi_inname         = ''
    character(len=CHEM_MAXSTR) :: fireemi_inname     = ''
    character(len=CHEM_MAXSTR) :: emi_outname        = ''
    character(len=CHEM_MAXSTR) :: gsi_path        = ''
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
    integer :: wetdep_ls_opt       =  WDLS_OPT_GSD
    integer :: plumerise_flag      = 0   ! default is none1g
    integer :: plumerisefire_frq  = 0   ! no plume rise
    INTEGER :: num_plume_data  = 0
    INTEGER :: chem_conv_tr     = CTRA_OPT_NONE
    real(CHEM_KIND_R4) :: raqmschemdt
    ! -- parameters
  ! -- configuration variables for output:
  !  . grid level defined in FIM Makefile
    integer :: glvl            = 0
  ! -- control variables


  end type raqmschem_config_type

  private

  public :: raqmschem_config_type
  ! -- provide subtypes
!  public :: dofirehour
  public :: raqmschem_config_read
  public :: raqmschem_config_control_init
  public :: FIRE_OPT_NONE,    &
            FIRE_OPT_MODIS,   &
            FIRE_OPT_GBBEPx
  public :: WDLS_OPT_NONE,    &
            WDLS_OPT_GSD,     &
            WDLS_OPT_NGAC,WDLS_OPT_NGAC_BOTH
  public :: CTRA_OPT_NONE,    &
            CTRA_OPT_GRELL
!  public :: CHEM_OPT_NONE,         &
!            CHEM_OPT_MAX

contains

  subroutine raqmschem_config_read(config, rc)

    type(raqmschem_config_type), intent(inout) :: config
    integer, optional,      intent(out)   :: rc

    ! -- local variables
    integer, parameter :: unit = 200

    integer                :: localrc, iostat
    integer                :: raqmschem_opt
    integer                :: wetdep_ls_opt
    integer                :: plumerise_flag
    integer                :: plumerisefire_frq 
    integer                :: buffer(5)
    integer                :: chem_conv_tr
    real(CHEM_KIND_R4)     :: raqmschemdt
    real(CHEM_KIND_R4)     :: rbuffer(1)
    character(CHEM_MAXSTR) :: sbuffer(5)

    ! -- variables in input namelist
    character(len=CHEM_MAXSTR) :: emi_inname
    character(len=CHEM_MAXSTR) :: fireemi_inname
    character(len=CHEM_MAXSTR) :: emi_outname
    character(len=CHEM_MAXSTR) :: gsi_path
    character(len=CHEM_MAXSTR) :: fireemi_outname
    character(len=CHEM_MAXSTR) :: input_chem_outname
    character(len=CHEM_MAXSTR) :: chem_hist_outname


    namelist /raqmschem_nml/          &
      emi_inname,                &
      fireemi_inname,            &
      emi_outname,               &
      gsi_path,                  &
      fireemi_outname,           &
      input_chem_outname,        &
      wetdep_ls_opt,             &
      chem_conv_tr,              &
      plumerise_flag,            &
      plumerisefire_frq,        &
      chem_hist_outname         

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    ! -- set defaults
    emi_inname         = ""
    fireemi_inname     = ""
    emi_outname        = ""
    gsi_path           = ''
    fireemi_outname    = ""
    input_chem_outname = ""
    chem_hist_outname  = "chem_out_"
    raqmschem_opt      = 0
    wetdep_ls_opt      = 0
    plumerise_flag     = 0
    raqmschemdt        = 3._CHEM_KIND_R4
    wetdep_ls_opt      = WDLS_OPT_GSD
    plumerisefire_frq  = 0
    chem_conv_tr       = CTRA_OPT_NONE

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
      raqmschem_opt,wetdep_ls_opt,plumerise_flag,plumerisefire_frq,chem_conv_tr  &
      /)
    ! -- broadcast integer buffer
    call chem_comm_bcast(buffer, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    ! -- load integer values in local config
    config % raqmschem_opt          = buffer( 1 )
    config % wetdep_ls_opt          = buffer( 2 )
    config % plumerise_flag         = buffer (3)
    config % plumerisefire_frq      = buffer (4)
    config % chem_conv_tr           = buffer (5)

    ! -- pack real variables in buffer
    rbuffer = (/ raqmschemdt/)
    ! -- broadcast real buffer
    call chem_comm_bcast(rbuffer, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    ! -- local real values in local buffer
    config % raqmschemdt     = rbuffer(1)

    ! -- pack strings into buffer
    sbuffer = (/ chem_hist_outname, emi_inname, fireemi_inname, emi_outname, gsi_path /)
    ! -- broadcast string variable
    call chem_comm_bcast(sbuffer, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    ! -- set string values to config
    config % chem_hist_outname = sbuffer(1)
    config % emi_inname        = sbuffer(2)
    config % fireemi_inname        = sbuffer(3)
    config % emi_outname       = sbuffer(4)
    config % gsi_path          = sbuffer(5)
!    write(6,*)'bottom raqmschem_config_read'
!    call flush(6)

  end subroutine raqmschem_config_read

  subroutine raqmschem_config_control_init(config, rc)

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
!   add fire
    ! -- fire options
!    write(6,*)'config%plumerise_flag',config%plumerise_flag
!    write(6,*)'opt',FIRE_OPT_NONE ,FIRE_OPT_MODIS,FIRE_OPT_GBBEPx
!    call flush(6)
    select case (config % plumerise_flag)
      case (FIRE_OPT_NONE)
        ! -- valid option
      case (FIRE_OPT_MODIS)
        config % num_plume_data = 8
      case (FIRE_OPT_GBBEPx)
        ! -- valid option
        config % num_plume_data = 1
      case (4) ! ajl till get code to make these
      case default
        call chem_rc_set(CHEM_RC_FAILURE, msg="plumerise_flag not implemented", &
          file=__FILE__, line=__LINE__, rc=rc)
        return
    end select
!    write(6,*)'num_plume_data',config%num_plume_data
!    call flush(6)

  end subroutine raqmschem_config_control_init

end module raqmschem_config_mod
