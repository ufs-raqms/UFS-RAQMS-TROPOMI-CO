module chem_iodata_mod

  use chem_rc_mod
  use chem_types_mod
  use chem_io_mod
  use chem_comm_mod
  use chem_model_mod
  use chem_config_mod

  implicit none

  private

  public :: chem_buffer_init
  public :: chem_backgd_init
  public :: chem_backgd_read
  public :: chem_backgd_write
  public :: chem_output_init
  public :: chem_output_write
  public :: chem_restart_read
  public :: chem_restart_write
  public :: dofirehour
  logical :: dofirehour=.false.

contains

  subroutine chem_buffer_init(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount
    integer :: ids, ide, jds, jde, ni
    type(chem_data_type),   pointer :: data   => null()
    type(chem_config_type), pointer :: config => null()

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, ni=ni, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      ! -- biomass burning emission buffers
      if (.not.allocated(data % eburn)) then
        allocate(data % eburn(ids:ide,1:ni,jds:jde,config % num_ebu), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % eburn = 0._CHEM_KIND_R4
      end if

    end do

  end subroutine chem_buffer_init


  subroutine chem_backgd_init(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount
    integer :: ids, ide, jds, jde
    type(chem_data_type),   pointer :: data   => null()
    type(chem_config_type), pointer :: config => null()

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      ! -- pressure levels for GOCART
      if (.not.allocated(data % p_gocart)) then
        allocate(data % p_gocart(config % nvl_gocart), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % p_gocart = (/ 1000., 992.5, 985., 977.5, 970., 955., 940., 925., 910.,    &
          895., 880., 865., 850., 825., 800., 775., 750., 712.5,  675., 637.5, 600.,     &
          562.5, 525., 487.5, 450., 412.5, 375., 337.5, 288.08, 244.88, 208.15, 176.93,  &
          150.39, 127.84, 108.66, 92.37, 78.51, 66.6, 56.39, 47.64, 40.18, 33.81, 28.37, &
          23.73, 19.79,  16.46, 13.64, 11.28, 9.29, 7.62, 6.22, 5.05, 4.08, 3.28, 2.62,  &
          2.08, 1.65, 1.3, 1.02, 0.8, 0.62, 0.48, 0.37, 0.28 /)
      end if

      ! -- dust
      if (.not.allocated(data % dm0)) then
        allocate(data % dm0(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % dm0 = 0._CHEM_KIND_R4
      end if

      ! -- dust erosion factors
      if (.not.allocated(data % ero1)) then
        allocate(data % ero1(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % ero1 = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % ero2)) then
        allocate(data % ero2(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % ero2 = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % ero3)) then
        allocate(data % ero3(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % ero3 = 0._CHEM_KIND_R4
      end if

      ! -- additional dust quantities for AFWA
      if (.not.allocated(data % clayfrac)) then
        allocate(data % clayfrac(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % clayfrac = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % sandfrac)) then
        allocate(data % sandfrac(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % sandfrac = 0._CHEM_KIND_R4
      end if

      ! -- drag partition map (FENGSHA)
      if (.not.allocated(data % rdrag)) then
        allocate(data % rdrag(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % rdrag = 0._CHEM_KIND_R4
      end if

      ! -- threshold friction velocity map (FENGSHA)
      if (.not.allocated(data % uthr)) then
        allocate(data % uthr(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % uthr = 0._CHEM_KIND_R4
      end if

      ! -- PJZ sediment supply map
      if (.not.allocated(data % ssm)) then
        allocate(data % ssm(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % ssm = 0._CHEM_KIND_R4
      end if

      ! -- chemical species background
      if (.not.allocated(data % h2o2_backgd)) then
        allocate(data % h2o2_backgd(ids:ide,jds:jde,config % nvl_gocart), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % h2o2_backgd = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % no3_backgd)) then
        allocate(data % no3_backgd(ids:ide,jds:jde,config % nvl_gocart), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % no3_backgd = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % oh_backgd)) then
        allocate(data % oh_backgd(ids:ide,jds:jde,config % nvl_gocart), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % oh_backgd = 0._CHEM_KIND_R4
      end if

      ! -- emissions
      if (.not.allocated(data % emiss_ab)) then
        allocate(data % emiss_ab(ids:ide,jds:jde,config % num_emis_ant), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % emiss_ab = 0._CHEM_KIND_R4
      end if

      ! -- volcanic ash
      if (.not.allocated(data % emiss_ash_height)) then
        allocate(data % emiss_ash_height(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % emiss_ash_height = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % emiss_ash_mass)) then
        allocate(data % emiss_ash_mass(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % emiss_ash_mass = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % emiss_ash_dt)) then
        allocate(data % emiss_ash_dt(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % emiss_ash_dt = 0._CHEM_KIND_R4
      end if

      ! -- emission tr
      if (.not.allocated(data % emiss_tr_height)) then
        allocate(data % emiss_tr_height(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % emiss_tr_height = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % emiss_tr_mass)) then
        allocate(data % emiss_tr_mass(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % emiss_tr_mass = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % emiss_tr_dt)) then
        allocate(data % emiss_tr_dt(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % emiss_tr_dt = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % th_pvsrf)) then
        allocate(data % th_pvsrf(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % th_pvsrf = 0._CHEM_KIND_R4
      end if

      ! -- emission from burning biomass
      if (.not.allocated(data % emiss_abu)) then
        allocate(data % emiss_abu(ids:ide,jds:jde,config % num_ebu_in), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % emiss_abu = 0._CHEM_KIND_R4
      end if
      if (.not.allocated(data % plume)) then
        allocate(data % plume(ids:ide,jds:jde,config % num_plume_data), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % plume = 0._CHEM_KIND_R4
      end if

    end do

  end subroutine chem_backgd_init

  subroutine chem_backgd_read(verbose, rc)
    logical, optional, intent(in)  :: verbose
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount, localpe, tile
    logical :: isVerbose
    type(chem_data_type),   pointer :: data => null()
    type(chem_config_type), pointer :: config => null()
    logical :: emishr=.false.
!   ajl for hourly
    integer ids,ide,jds,jde,its,ite,jts,jte
    real(chem_kind_r4),allocatable,dimension(:,:,:)  :: srcinhr

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    isVerbose = .false.
    if (present(verbose)) isVerbose = verbose

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_comm_get(localpe=localpe, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, tile=tile, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if ((config % chem_opt == CHEM_OPT_RACM_SOA_VBS) .or.  &
          (config % chem_opt >= CHEM_OPT_GOCART)       .and. &
          (config % chem_opt < 500)) then

        ! -- dust erosion factors
        call chem_io_read('dm0.dat', data % dm0, path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," dm0 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % dm0), maxval(data % dm0)

        ! -- dust erosion factors
        call chem_io_read('erod1.dat', data % ero1, path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ero1 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % ero1), maxval(data % ero1)
        call chem_io_read('erod2.dat', data % ero2, path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ero2 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % ero2), maxval(data % ero2)
        call chem_io_read('erod3.dat', data % ero3, path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ero3 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % ero3), maxval(data % ero3)

        ! -- background values for chemical species
        call chem_io_read('h2o2.dat', data % h2o2_backgd, path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," h2o2 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % h2o2_backgd), maxval(data % h2o2_backgd)
        call chem_io_read('no3.dat', data % no3_backgd, path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," no3 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % no3_backgd), maxval(data % no3_backgd)
        call chem_io_read('oh.dat', data % oh_backgd, path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," oh - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % oh_backgd), maxval(data % oh_backgd)

        ! -- emissions
        call chem_io_read('e_bc.dat', data % emiss_ab(:,:,config % species % p_e_bc), &
          path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_bc - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_bc)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_bc))

        call chem_io_read('e_oc.dat', data % emiss_ab(:,:,config % species % p_e_oc), &
          path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_oc - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_oc)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_oc))

!       call chem_io_read('e_pm_10.dat', data % emiss_ab(:,:,config % species % p_e_pm_10), &
!         path=trim(config % emi_inname), de=de, rc=localrc)
!       if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!       if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_pm_10 - min/max = "2g16.6)') &
!         localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_pm_10)), &
!         maxval(data % emiss_ab(:,:,config % species % p_e_pm_10))

        call chem_io_read('e_pm_25.dat', data % emiss_ab(:,:,config % species % p_e_pm_25), &
          path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_pm_25 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_pm_25)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_pm_25))

        call chem_io_read('e_so2.dat', data % emiss_ab(:,:,config % species % p_e_so2), &
          path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_so2 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_so2)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_so2))

        call chem_io_read('e_sulf.dat', data % emiss_ab(:,:,config % species % p_e_sulf), &
          path=trim(config % emi_inname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_sulf - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_sulf)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_sulf))

        select case (config % dust_opt)
          case (DUST_OPT_AFWA, DUST_OPT_FENGSHA)
            call chem_io_read('clay.dat', data % clayfrac, path=trim(config % dust_inname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," clayfrac - min/max = "2g16.6)') &
              localpe, de, tile, minval(data % clayfrac), maxval(data % clayfrac)
            call chem_io_read('sand.dat', data % sandfrac, path=trim(config % dust_inname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," sandfrac - min/max = "2g16.6)') &
              localpe, de, tile, minval(data % sandfrac), maxval(data % sandfrac)
        end select

        if (config % dust_opt == DUST_OPT_FENGSHA) then
          ! -- PJZ sediment supply map
          call chem_io_read('ssm.dat', data % ssm, path=trim(config % dust_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ssm - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % ssm), maxval(data % ssm)
          ! -- threshold friction velocity map
          call chem_io_read('uthr.dat', data % uthr, path=trim(config % dust_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," uthr - min/max = "2g16.6)') &
               localpe, de, tile, minval(data % uthr), maxval(data % uthr)
          if (config % dust_calcdrag == 1) then
            ! -- drag partition map
            call chem_io_read('rdrag.dat', data % rdrag, path=trim(config % dust_inname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," rdrag - min/max = "2g16.6)') &
              localpe, de, tile, minval(data % rdrag), maxval(data % rdrag)
          end if
        end if

        if ((config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
            (config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then

          call chem_io_read('e_ald.dat', data % emiss_ab(:,:,config % species % p_e_ald), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_ald - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_ald)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_ald))

          call chem_io_read('e_co.dat', data % emiss_ab(:,:,config % species % p_e_co), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_co - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_co)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_co))

          call chem_io_read('e_csl.dat', data % emiss_ab(:,:,config % species % p_e_csl), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_csl - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_csl)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_csl))

          call chem_io_read('e_dms.dat', data % emiss_ab(:,:,config % species % p_e_dms), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_dms - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_dms)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_dms))

          call chem_io_read('e_eth.dat', data % emiss_ab(:,:,config % species % p_e_eth), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_eth - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_eth)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_eth))

          call chem_io_read('e_hc3.dat', data % emiss_ab(:,:,config % species % p_e_hc3), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_hc3 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_hc3)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_hc3))

          call chem_io_read('e_hc5.dat', data % emiss_ab(:,:,config % species % p_e_hc5), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_hc5 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_hc5)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_hc5))

          call chem_io_read('e_hc8.dat', data % emiss_ab(:,:,config % species % p_e_hc8), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_hc8 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_hc8)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_hc8))

          call chem_io_read('e_hcho.dat', data % emiss_ab(:,:,config % species % p_e_hcho), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_hcho - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_hcho)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_hcho))

          call chem_io_read('e_iso.dat', data % emiss_ab(:,:,config % species % p_e_iso), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_iso - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_iso)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_iso))

          call chem_io_read('e_ket.dat', data % emiss_ab(:,:,config % species % p_e_ket), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_ket - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_ket)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_ket))

          call chem_io_read('e_nh3.dat', data % emiss_ab(:,:,config % species % p_e_nh3), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_nh3 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_nh3)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_nh3))

          call chem_io_read('e_no2.dat', data % emiss_ab(:,:,config % species % p_e_no2), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_no2 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_no2)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_no2))

          call chem_io_read('e_no.dat', data % emiss_ab(:,:,config % species % p_e_no), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_no - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_no)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_no))

          call chem_io_read('e_oli.dat', data % emiss_ab(:,:,config % species % p_e_oli), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_oli - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_oli)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_oli))

          call chem_io_read('e_olt.dat', data % emiss_ab(:,:,config % species % p_e_olt), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_olt - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_olt)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_olt))

          call chem_io_read('e_ora2.dat', data % emiss_ab(:,:,config % species % p_e_ora2), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_ora2 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_ora2)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_ora2))

          call chem_io_read('e_tol.dat', data % emiss_ab(:,:,config % species % p_e_tol), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_tol - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_tol)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_tol))

          call chem_io_read('e_xyl.dat', data % emiss_ab(:,:,config % species % p_e_xyl), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," e_xyl - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_xyl)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_xyl))

        end if
      end if

      ! -- emissions from burning biomass
      if (config % biomass_burn_opt == BURN_OPT_ENABLE) then
        ! -- emissions
            if(localpe==0)then
              write(6,*)'read ebu_bc.dat from fireemi_inname ',trim(config % fireemi_inname)
            endif
        emishr=.false.
        if(inquire_file_var('ebu_bc.tile1.nc','ebu_bc',trim(config%fireemi_inname)))then 
          if(inquire_file_dim('ebu_bc.tile1.nc','nhr',trim(config%fireemi_inname)))then
            emishr=.true.
          endif
        endif
        if(emishr)then
            dofirehour=.true. 
            call chem_model_domain_get(ids=ids, ide=ide, jds=jds, jde=jde, & 
            its=its,ite=ite,jts=jts,jte=jte)
            allocate (srcinhr(ide-ids+1,jds:jde,24))
            allocate (data%emiss_abu_hr(ide-ids+1,jds:jde,24,config%num_ebu_in))
!            write(200+localpe,*)'allocate srcinhr',ide-ids+1,jds,jde
!            call flush(200+localpe)
!            write(200+localpe,*)'p_e_pm_25',config % species % p_e_pm_25, &
!               config%num_ebu_in
!            call flush(200+localpe)
!          call chem_io_read('ebu_pm_25.dat', srcinhr, &
          call chem_io_read('ebu_pm_25.dat', data % emiss_abu_hr(:,:,:,config % species % p_e_pm_25), &
          path=trim(config % fireemi_inname), de=de, varname='ebu_pm25',rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!          data % emiss_abu_hr(:,:,:,config % species % p_e_pm_25)=srcinhr
!          write(200+localpe,*)'done read pm_25'
!          call flush(200+localpe)
!          write(6,*)'read bc',localpe
!          write(200+localpe,*)'read bc',localpe
!            call flush(200+localpe)
          !call flush(6)
!          call chem_io_read('ebu_bc.dat', srcinhr, &
          call chem_io_read('ebu_bc.dat', data % emiss_abu_hr(:,:,:,config % species % p_e_bc), &
          path=trim(config % fireemi_inname), de=de, varname='ebu_bc',rc=localrc)
!          write(200+localpe,*)'did read bc',localpe
!            call flush(200+localpe)
!          write(6,*)'did read bc',localpe,localrc
!          call flush(6)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!          data % emiss_abu_hr(:,:,:,config % species % p_e_bc)=srcinhr
!          call chem_io_read('ebu_oc.dat', srcinhr, &
          call chem_io_read('ebu_oc.dat', data % emiss_abu_hr(:,:,:,config % species % p_e_oc), &
          path=trim(config % fireemi_inname), de=de, varname='ebu_oc',rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!          data % emiss_abu_hr(:,:,:,config % species % p_e_oc)=srcinhr
!          write(200+localpe,*)'done read oc'
!          call flush(200+localpe)
!          call chem_io_read('ebu_co.dat', srcinhr, &
!          path=trim(config % fireemi_inname), de=de, varname='ebu_co',rc=localrc)
!          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!          data % emiss_abu_hr(:,:,:,config % species % p_e_co)=srcinhr
!          call chem_io_read('ebu_so2.dat', srcinhr, &
          call chem_io_read('ebu_so2.dat', data % emiss_abu_hr(:,:,:,config % species % p_e_so2), &
          path=trim(config % fireemi_inname), de=de, varname='ebu_so2',rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!          data % emiss_abu_hr(:,:,:,config % species % p_e_co)=srcinhr
!          write(200+localpe,*)'done read so2'
!          call flush(6)
        else
!          write(200+localpe,*)'no emishr'
!          call flush(200+localpe)

          call chem_io_read('ebu_bc.dat', data % emiss_abu(:,:,config % species % p_e_bc), &
          path=trim(config % fireemi_inname), de=de, varname='ebu_bc',rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_bc - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_bc)), &
          maxval(data % emiss_abu(:,:,config % species % p_e_bc))

          call chem_io_read('ebu_oc.dat', data % emiss_abu(:,:,config % species % p_e_oc), &
          path=trim(config % fireemi_inname), de=de, varname='ebu_oc',rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_oc - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_oc)), &
          maxval(data % emiss_abu(:,:,config % species % p_e_oc))

          call chem_io_read('ebu_pm_25.dat', data % emiss_abu(:,:,config % species % p_e_pm_25), &
          path=trim(config % fireemi_inname), de=de, varname='ebu_pm25',rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_pm_25 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_pm_25)), &
          maxval(data % emiss_abu(:,:,config % species % p_e_pm_25))

          call chem_io_read('ebu_so2.dat', data % emiss_abu(:,:,config % species % p_e_so2), &
          path=trim(config % fireemi_inname), de=de, varname='ebu_so2',rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_so2 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_so2)), &
          maxval(data % emiss_abu(:,:,config % species % p_e_so2))
!        write(6,*)'plumerise_flag',config % plumerise_flag,'MODIS',FIRE_OPT_MODIS,'GBBEPX',FIRE_OPT_GBBEPx
!        call flush(6)
          endif
!        write(200+localpe,*)'flag',config % plumerise_flag
!        call flush(200+localpe)
        select case (config % plumerise_flag)
          case (FIRE_OPT_MODIS)
            call chem_io_read('ebu_pm_10.dat', data % emiss_abu(:,:,config % species % p_e_pm_10), &
              path=trim(config % fireemi_inname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_pm_10 - min/max = "2g16.6)') &
              localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_pm_10)), &
              maxval(data % emiss_abu(:,:,config % species % p_e_pm_10))
            call chem_io_read('plumestuff.dat', data % plume, recrange=(/ 1, config % num_plume_data /), &
              path=trim(config % fireemi_inname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," plume - min/max = "2g16.6)') &
              localpe, de, tile, minval(data % plume), maxval(data % plume)
          case (FIRE_OPT_GBBEPx)
            if (config % plumerisefire_frq > 0) then
            if(localpe==0)then
              write(6,*)'read plumefrp.dat from',trim(config % fireemi_inname)
              call flush(6)
            endif
            if(emishr)then
!              write(6,*)'read plumefrp',localpe
!              call flush(6)
!              write(200+localpe,*)'read plumefrp',localpe,allocated(srcinhr)
!              call flush(200+localpe)
              if(.not.allocated(data%plume_hr))then
                allocate (data%plume_hr(ide-ids+1,jds:jde,24))
              endif
!              write(200+localpe,*)'plume_hr',shape(data%plume_hr)
!              call flush(200+localpe)
              call chem_io_read('plumefrp.dat',data%plume_hr , &
              path=trim(config % fireemi_inname), de=de, varname='plumefrp',rc=localrc)
!              write(6,*)'did read plumefrp',localpe
!              call flush(6)
!              write(200+localpe,*)'did read plumefrp',localpe,allocated(srcinhr)
!              call flush(200+localpe)
!              data%plume_hr(:,:,:)=srcinhr
!              write(6,*)'deallocated srcinhr'
!              write(200+localpe,*)'deallocated srcinhr',allocated(srcinhr)
!              call flush(200+localpe)
!              if(allocated(srcinhr))then
!                deallocate(srcinhr)
!              endif
!              write(200+localpe,*)'did deallocated srcinhr'
!              call flush(200+localpe)
            else
!              write(200+localpe,*)'should not be here'
!              call flush(200+localpe)
              call chem_io_read('plumefrp.dat', data % plume, &
              path=trim(config % fireemi_inname), de=de, varname='plumefrp',rc=localrc)
              if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
              if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," plume - min/max = "2g16.6)') &
                localpe, de, tile, minval(data % plume), maxval(data % plume)
              endif
            endif
          case default
            ! -- no further options available
        end select
!        write(200+localpe,*)'afer select'
!        call flush(200+localpe)

        if ((config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
            (config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then

          call chem_io_read('ebu_ald.dat', data % emiss_abu(:,:,config % species % p_e_ald), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_ald - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_ald)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_ald))

          call chem_io_read('ebu_co.dat', data % emiss_abu(:,:,config % species % p_e_co), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_co - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_co)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_co))

          call chem_io_read('ebu_csl.dat', data % emiss_abu(:,:,config % species % p_e_csl), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_csl - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_csl)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_csl))

          call chem_io_read('ebu_dms.dat', data % emiss_abu(:,:,config % species % p_e_dms), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_dms - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_dms)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_dms))

          call chem_io_read('ebu_eth.dat', data % emiss_abu(:,:,config % species % p_e_eth), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_eth - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_eth)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_eth))

          call chem_io_read('ebu_hc3.dat', data % emiss_abu(:,:,config % species % p_e_hc3), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_hc3 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_hc3)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_hc3))

          call chem_io_read('ebu_hc5.dat', data % emiss_abu(:,:,config % species % p_e_hc5), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_hc5 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_hc5)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_hc5))

          call chem_io_read('ebu_hc8.dat', data % emiss_abu(:,:,config % species % p_e_hc8), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_hc8 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_hc8)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_hc8))

          call chem_io_read('ebu_hcho.dat', data % emiss_abu(:,:,config % species % p_e_hcho), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_hcho - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_hcho)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_hcho))

          call chem_io_read('ebu_iso.dat', data % emiss_abu(:,:,config % species % p_e_iso), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_iso - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_iso)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_iso))

          call chem_io_read('ebu_ket.dat', data % emiss_abu(:,:,config % species % p_e_ket), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_ket - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_ket)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_ket))

          call chem_io_read('ebu_nh3.dat', data % emiss_abu(:,:,config % species % p_e_nh3), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_nh3 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_nh3)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_nh3))

          call chem_io_read('ebu_no2.dat', data % emiss_abu(:,:,config % species % p_e_no2), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_no2 - min/max = "2g16.6)') &
             localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_no2)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_no2))

          call chem_io_read('ebu_no.dat', data % emiss_abu(:,:,config % species % p_e_no), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_no - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_no)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_no))

          call chem_io_read('ebu_oli.dat', data % emiss_abu(:,:,config % species % p_e_oli), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_oli - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_oli)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_oli))

          call chem_io_read('ebu_olt.dat', data % emiss_abu(:,:,config % species % p_e_olt), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_olt - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_olt)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_olt))

          call chem_io_read('ebu_ora2.dat', data % emiss_abu(:,:,config % species % p_e_ora2), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_ora2 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_ora2)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_ora2))

          call chem_io_read('ebu_tol.dat', data % emiss_abu(:,:,config % species % p_e_tol), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_tol - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_tol)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_tol))

          call chem_io_read('ebu_xyl.dat', data % emiss_abu(:,:,config % species % p_e_xyl), &
            path=trim(config % fireemi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," ebu_xyl - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_xyl)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_xyl))

        end if
      end if

      ! -- volcanic stuff
      if (config % chem_opt == CHEM_OPT_GOCART) then
         ! -- also for chem_opt = 316, 317, 502

        if (config % ash_mass > -900._CHEM_KIND_R4) then
          ! -- TODO
          call chem_io_read('volcanic.dat', data % emiss_ash_mass,recrange=(/ 4, 4 /), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," emiss_ash_mass - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ash_mass), maxval(data % emiss_ash_mass)

          call chem_io_read('volcanic.dat', data % emiss_ash_height,recrange=(/ 5, 5 /), &
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," emiss_ash_height - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ash_height), maxval(data % emiss_ash_height)

          call chem_io_read('volcanic.dat', data % emiss_ash_dt, recrange=(/ 6, 6 /),&
            path=trim(config % emi_inname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," emiss_ash_dt - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ash_dt), maxval(data % emiss_ash_dt)

        end if
        ! -- overwrite ash_mass if namelist value exists
        if (config % ash_mass > -100._CHEM_KIND_R4)then
          write(0,*)'using namelist value for ash_mass'
          where(data % emiss_ash_mass > 0._CHEM_KIND_R4) data % emiss_ash_mass = config % ash_mass
        endif
        ! -- overwrite ash_height if namelist value exists
        if (config % ash_height > 0._CHEM_KIND_R4) then

          write(0,*)'using namelist value for ash_height'
          where(data % emiss_ash_height > 0._CHEM_KIND_R4) data % emiss_ash_height = config % ash_height
          where(data % emiss_ash_height < 1._CHEM_KIND_R4) data % emiss_ash_dt     = 0._CHEM_KIND_R4

        else if (config % ash_height < -990._CHEM_KIND_R4) then

          write(0,*)'resetting all ash variables to zero'
          data % emiss_ash_mass   = 0._CHEM_KIND_R4
          data % emiss_ash_height = 0._CHEM_KIND_R4
          data % emiss_ash_dt     = 0._CHEM_KIND_R4

        endif
      end if

    end do
!        write(200+localpe,*)'bottom chem_backgd_read'
!        call flush(200+localpe)

  end subroutine chem_backgd_read


  subroutine chem_backgd_write(verbose, rc)
    logical, optional, intent(in)  :: verbose
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount, localpe, tile
    logical :: isVerbose
    type(chem_data_type),   pointer :: data => null()
    type(chem_config_type), pointer :: config => null()

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    isVerbose = .false.
    if (present(verbose)) isVerbose = verbose

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_comm_get(localpe=localpe, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, tile=tile, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if ((config % chem_opt == CHEM_OPT_RACM_SOA_VBS) .or.  &
          (config % chem_opt >= CHEM_OPT_GOCART)       .and. &
          (config % chem_opt < 500)) then

        ! -- dust erosion factors
        call chem_io_write('dm0.dat', data % dm0, path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," dm0 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % dm0), maxval(data % dm0)

        ! -- dust erosion factors
        call chem_io_write('erod1.dat', data % ero1, path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ero1 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % ero1), maxval(data % ero1)
        call chem_io_write('erod2.dat', data % ero2, path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ero2 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % ero2), maxval(data % ero2)
        call chem_io_write('erod3.dat', data % ero3, path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ero3 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % ero3), maxval(data % ero3)

        ! -- background values for chemical species
        call chem_io_write('h2o2.dat', data % h2o2_backgd, path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," h2o2 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % h2o2_backgd), maxval(data % h2o2_backgd)
        call chem_io_write('no3.dat', data % no3_backgd, path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," no3 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % no3_backgd), maxval(data % no3_backgd)
        call chem_io_write('oh.dat', data % oh_backgd, path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," oh - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % oh_backgd), maxval(data % oh_backgd)

        ! -- emissions
        call chem_io_write('e_bc.dat', data % emiss_ab(:,:,config % species % p_e_bc), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_bc - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_bc)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_bc))

        call chem_io_write('e_oc.dat', data % emiss_ab(:,:,config % species % p_e_oc), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_oc - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_oc)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_oc))

        call chem_io_write('e_pm_10.dat', data % emiss_ab(:,:,config % species % p_e_pm_10), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_pm_10 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_pm_10)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_pm_10))

        call chem_io_write('e_pm_25.dat', data % emiss_ab(:,:,config % species % p_e_pm_25), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_pm_25 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_pm_25)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_pm_25))

        call chem_io_write('e_so2.dat', data % emiss_ab(:,:,config % species % p_e_so2), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_so2 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_so2)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_so2))

        call chem_io_write('e_sulf.dat', data % emiss_ab(:,:,config % species % p_e_sulf), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_sulf - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_sulf)), &
          maxval(data % emiss_ab(:,:,config % species % p_e_sulf))

        select case (config % dust_opt)
          case (DUST_OPT_AFWA, DUST_OPT_FENGSHA)
            call chem_io_write('clay.dat', data % clayfrac, path=trim(config % emi_outname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," clayfrac - min/max = "2g16.6)') &
              localpe, de, tile, minval(data % clayfrac), maxval(data % clayfrac)
            call chem_io_write('sand.dat', data % sandfrac, path=trim(config % emi_outname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," sandfrac - min/max = "2g16.6)') &
              localpe, de, tile, minval(data % sandfrac), maxval(data % sandfrac)
        end select

        if (config % dust_opt == DUST_OPT_FENGSHA) then
          ! -- PJZ sediment supply map
          call chem_io_write('ssm.dat', data % ssm, path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ssm - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % ssm), maxval(data % ssm)
          ! -- threshold friction velocity map
          call chem_io_write('uthr.dat', data % uthr, path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," uthr - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % uthr), maxval(data % uthr)
          if (config % dust_calcdrag == 1) then
            ! -- drag partition map
            call chem_io_write('rdrag.dat', data % rdrag, path=trim(config % emi_outname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," rdrag - min/max = "2g16.6)') &
              localpe, de, tile, minval(data % rdrag), maxval(data % rdrag)
          end if
        end if

        if ((config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
            (config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then

          call chem_io_write('e_ald.dat', data % emiss_ab(:,:,config % species % p_e_ald), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_ald - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_ald)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_ald))

          call chem_io_write('e_co.dat', data % emiss_ab(:,:,config % species % p_e_co), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_co - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_co)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_co))

          call chem_io_write('e_csl.dat', data % emiss_ab(:,:,config % species % p_e_csl), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_csl - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_csl)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_csl))

          call chem_io_write('e_dms.dat', data % emiss_ab(:,:,config % species % p_e_dms), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_dms - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_dms)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_dms))

          call chem_io_write('e_eth.dat', data % emiss_ab(:,:,config % species % p_e_eth), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_eth - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_eth)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_eth))

          call chem_io_write('e_hc3.dat', data % emiss_ab(:,:,config % species % p_e_hc3), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_hc3 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_hc3)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_hc3))

          call chem_io_write('e_hc5.dat', data % emiss_ab(:,:,config % species % p_e_hc5), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_hc5 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_hc5)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_hc5))

          call chem_io_write('e_hc8.dat', data % emiss_ab(:,:,config % species % p_e_hc8), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_hc8 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_hc8)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_hc8))

          call chem_io_write('e_hcho.dat', data % emiss_ab(:,:,config % species % p_e_hcho), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_hcho - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_hcho)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_hcho))

          call chem_io_write('e_iso.dat', data % emiss_ab(:,:,config % species % p_e_iso), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_iso - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_iso)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_iso))

          call chem_io_write('e_ket.dat', data % emiss_ab(:,:,config % species % p_e_ket), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_ket - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_ket)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_ket))

          call chem_io_write('e_nh3.dat', data % emiss_ab(:,:,config % species % p_e_nh3), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_nh3 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_nh3)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_nh3))

          call chem_io_write('e_no2.dat', data % emiss_ab(:,:,config % species % p_e_no2), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_no2 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_no2)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_no2))

          call chem_io_write('e_no.dat', data % emiss_ab(:,:,config % species % p_e_no), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_no - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_no)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_no))

          call chem_io_write('e_oli.dat', data % emiss_ab(:,:,config % species % p_e_oli), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_oli - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_oli)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_oli))

          call chem_io_write('e_olt.dat', data % emiss_ab(:,:,config % species % p_e_olt), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_olt - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_olt)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_olt))

          call chem_io_write('e_ora2.dat', data % emiss_ab(:,:,config % species % p_e_ora2), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_ora2 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_ora2)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_ora2))

          call chem_io_write('e_tol.dat', data % emiss_ab(:,:,config % species % p_e_tol), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_tol - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_tol)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_tol))

          call chem_io_write('e_xyl.dat', data % emiss_ab(:,:,config % species % p_e_xyl), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," e_xyl - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_ab(:,:,config % species % p_e_xyl)), &
            maxval(data % emiss_ab(:,:,config % species % p_e_xyl))

        end if
      end if

      ! -- emissions from burning biomass
      if (config % biomass_burn_opt == BURN_OPT_ENABLE) then
        ! -- emissions
        call chem_io_write('ebu_bc.dat', data % emiss_abu(:,:,config % species % p_e_bc), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_bc - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_bc)), &
          maxval(data % emiss_abu(:,:,config % species % p_e_bc))

        call chem_io_write('ebu_oc.dat', data % emiss_abu(:,:,config % species % p_e_oc), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_oc - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_oc)), &
          maxval(data % emiss_abu(:,:,config % species % p_e_oc))

        if (config % plumerise_flag == FIRE_OPT_MODIS) then
          call chem_io_write('ebu_pm_10.dat', data % emiss_abu(:,:,config % species % p_e_pm_10), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_pm_10 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_pm_10)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_pm_10))
        end if

        call chem_io_write('ebu_pm_25.dat', data % emiss_abu(:,:,config % species % p_e_pm_25), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_pm_25 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_pm_25)), &
          maxval(data % emiss_abu(:,:,config % species % p_e_pm_25))

        call chem_io_write('ebu_so2.dat', data % emiss_abu(:,:,config % species % p_e_so2), &
          path=trim(config % emi_outname), de=de, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_so2 - min/max = "2g16.6)') &
          localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_so2)), &
          maxval(data % emiss_abu(:,:,config % species % p_e_so2))

        if ((config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
            (config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then

          call chem_io_write('ebu_ald.dat', data % emiss_abu(:,:,config % species % p_e_ald), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_ald - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_ald)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_ald))

          call chem_io_write('ebu_co.dat', data % emiss_abu(:,:,config % species % p_e_co), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_co - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_co)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_co))

          call chem_io_write('ebu_csl.dat', data % emiss_abu(:,:,config % species % p_e_csl), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_csl - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_csl)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_csl))

          call chem_io_write('ebu_dms.dat', data % emiss_abu(:,:,config % species % p_e_dms), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_dms - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_dms)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_dms))

          call chem_io_write('ebu_eth.dat', data % emiss_abu(:,:,config % species % p_e_eth), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_eth - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_eth)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_eth))

          call chem_io_write('ebu_hc3.dat', data % emiss_abu(:,:,config % species % p_e_hc3), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_hc3 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_hc3)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_hc3))

          call chem_io_write('ebu_hc5.dat', data % emiss_abu(:,:,config % species % p_e_hc5), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_hc5 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_hc5)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_hc5))

          call chem_io_write('ebu_hc8.dat', data % emiss_abu(:,:,config % species % p_e_hc8), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_hc8 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_hc8)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_hc8))

          call chem_io_write('ebu_hcho.dat', data % emiss_abu(:,:,config % species % p_e_hcho), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_hcho - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_hcho)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_hcho))

          call chem_io_write('ebu_iso.dat', data % emiss_abu(:,:,config % species % p_e_iso), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_iso - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_iso)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_iso))

          call chem_io_write('ebu_ket.dat', data % emiss_abu(:,:,config % species % p_e_ket), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_ket - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_ket)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_ket))

          call chem_io_write('ebu_nh3.dat', data % emiss_abu(:,:,config % species % p_e_nh3), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_nh3 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_nh3)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_nh3))

          call chem_io_write('ebu_no2.dat', data % emiss_abu(:,:,config % species % p_e_no2), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_no2 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_no2)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_no2))

          call chem_io_write('ebu_no.dat', data % emiss_abu(:,:,config % species % p_e_no), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_no - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_no)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_no))

          call chem_io_write('ebu_oli.dat', data % emiss_abu(:,:,config % species % p_e_oli), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_oli - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_oli)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_oli))

          call chem_io_write('ebu_olt.dat', data % emiss_abu(:,:,config % species % p_e_olt), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_olt - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_olt)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_olt))

          call chem_io_write('ebu_ora2.dat', data % emiss_abu(:,:,config % species % p_e_ora2), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_ora2 - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_ora2)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_ora2))

          call chem_io_write('ebu_tol.dat', data % emiss_abu(:,:,config % species % p_e_tol), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_tol - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_tol)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_tol))

          call chem_io_write('ebu_xyl.dat', data % emiss_abu(:,:,config % species % p_e_xyl), &
            path=trim(config % emi_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          if (isVerbose) write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," ebu_xyl - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % emiss_abu(:,:,config % species % p_e_xyl)), &
            maxval(data % emiss_abu(:,:,config % species % p_e_xyl))

        end if
      end if

    end do

  end subroutine chem_backgd_write

  subroutine chem_output_init(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount
    integer :: ids, ide, jds, jde, nvl, nt
    type(chem_data_type),   pointer :: data   => null()
    type(chem_config_type), pointer :: config => null()

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, nl=nvl, nt=nt, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (.not.allocated(data % aod2d)) then
        allocate(data % aod2d(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % aod2d = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % ext_cof)) then
        allocate(data % ext_cof(ids:ide,jds:jde,nvl,config % nbands), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % ext_cof = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % sscal)) then
        allocate(data % sscal(ids:ide,jds:jde,nvl,config % nbands), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % sscal = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % asymp)) then
        allocate(data % asymp(ids:ide,jds:jde,nvl,config % nbands), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % asymp = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % pm10)) then
        allocate(data % pm10(ids:ide,jds:jde,nvl), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % pm10 = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % pm25)) then
        allocate(data % pm25(ids:ide,jds:jde,nvl), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % pm25 = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % ebu_oc)) then
        allocate(data % ebu_oc(ids:ide,jds:jde,nvl), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % ebu_oc = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % oh_bg)) then
        allocate(data % oh_bg(ids:ide,jds:jde,nvl), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % oh_bg = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % h2o2_bg)) then
        allocate(data % h2o2_bg(ids:ide,jds:jde,nvl), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % h2o2_bg = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % no3_bg)) then
        allocate(data % no3_bg(ids:ide,jds:jde,nvl), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % no3_bg = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % wet_dep)) then
        allocate(data % wet_dep(ids:ide,jds:jde,config % num_chem), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % wet_dep = 0._CHEM_KIND_R4
      end if

      if (.not.allocated(data % trdp)) then
        allocate(data % trdp(ids:ide,jds:jde,nvl,nt), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
        data % trdp = 0._CHEM_KIND_R4
      end if

    end do

  end subroutine chem_output_init

  subroutine chem_output_write(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount
    integer :: n, p
    integer :: advanceCount
    type(chem_config_type),  pointer :: config   => null()
    type(chem_species_type), pointer :: s        => null()
    type(chem_data_type),    pointer :: data     => null()
    type(chem_state_type),   pointer :: stateOut => null()

    character(len=*), parameter :: filepos = 'append'
    character(len=*), parameter :: emnames(10) = (/ &
      'emd1', 'emd2', 'emd3', 'emd4', 'emd5', &
      'ems1', 'ems2', 'ems3', 'ems4', 'ems5'  &
      /)
    character(len=*), parameter :: cbnames(6) = (/   &
      'cbae', 'cbbc', 'cboc', 'cbsf', 'cbdt', 'cbss' &
      /)

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount < 1) return

    if (config % archive_step < 0) return

    call chem_model_clock_get(advanceCount=advanceCount, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (mod(advanceCount, config % archive_step) == 0) then

      s => config % species

      do de = 0, deCount-1
        call chem_model_get(de=de, data=data, stateOut=stateOut, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

        select case (config % chem_opt)
          case (CHEM_OPT_GOCART, CHEM_OPT_GOCART_RACM)

            call chem_io_write('d1st', stateOut % tr3d(:,:,:,config % ntra + s % p_dust_1), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('d2st', stateOut % tr3d(:,:,:,config % ntra + s % p_dust_2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('d3st', stateOut % tr3d(:,:,:,config % ntra + s % p_dust_3), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('d4st', stateOut % tr3d(:,:,:,config % ntra + s % p_dust_4), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('d5st', stateOut % tr3d(:,:,:,config % ntra + s % p_dust_5), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

            call chem_io_write('dms1', stateOut % tr3d(:,:,:,config % ntra + s % p_dms), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('pmsa', stateOut % tr3d(:,:,:,config % ntra + s % p_msa), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('s1ea', stateOut % tr3d(:,:,:,config % ntra + s % p_seas_1), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('s2ea', stateOut % tr3d(:,:,:,config % ntra + s % p_seas_2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('s3ea', stateOut % tr3d(:,:,:,config % ntra + s % p_seas_3), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('s4ea', stateOut % tr3d(:,:,:,config % ntra + s % p_seas_4), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('s5ea', stateOut % tr3d(:,:,:,config % ntra + s % p_seas_5), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

            call chem_io_write('pbc1', stateOut % tr3d(:,:,:,config % ntra + s % p_bc1), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('pbc2', stateOut % tr3d(:,:,:,config % ntra + s % p_bc2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('poc1', stateOut % tr3d(:,:,:,config % ntra + s % p_oc1), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('poc2', stateOut % tr3d(:,:,:,config % ntra + s % p_oc2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('pso2', stateOut % tr3d(:,:,:,config % ntra + s % p_so2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('sulf', stateOut % tr3d(:,:,:,config % ntra + s % p_sulf), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('pp25', stateOut % tr3d(:,:,:,config % ntra + s % p_p25), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('pp10', stateOut % tr3d(:,:,:,config % ntra + s % p_p10), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

            call chem_io_write('pm10', data % pm10, &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('pm25', data % pm25, &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

            ! -- write emissions
            if (associated(stateOut % truf)) then
              do n = 1, size(StateOut % truf, dim=3)
                call chem_io_write(emnames(n), stateOut % truf(:,:,n), &
                  path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
                if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
              end do
            end if

            ! -- write column mass density
            if (associated(stateOut % trcm)) then
              do p = 1, size(StateOut % trcm, dim=3)
                call chem_io_write(cbnames(p), stateOut % trcm(:,:,p), &
                  path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
                if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
              end do
            end if

            call chem_io_write('ao2D', data % aod2d, &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

            call chem_io_write('wpre', data % wet_dep(:,:,s % p_so2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            call chem_io_write('wbc2', data % wet_dep(:,:,s % p_bc2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('woc2', data % wet_dep(:,:,s % p_oc2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wp10', data % wet_dep(:,:,s % p_p10), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wp25', data % wet_dep(:,:,s % p_p25), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wso4', data % wet_dep(:,:,s % p_sulf), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wdt1', data % wet_dep(:,:,s % p_dust_1), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wdt2', data % wet_dep(:,:,s % p_dust_2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wdt3', data % wet_dep(:,:,s % p_dust_3), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wdt4', data % wet_dep(:,:,s % p_dust_4), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wdt5', data % wet_dep(:,:,s % p_dust_5), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wse1', data % wet_dep(:,:,s % p_seas_1), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wse2', data % wet_dep(:,:,s % p_seas_2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wse3', data % wet_dep(:,:,s % p_seas_3), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wse4', data % wet_dep(:,:,s % p_seas_4), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('wse5', data % wet_dep(:,:,s % p_seas_5), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

            call chem_io_write('aso2', data % emiss_ab(:,:,s % p_e_so2), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('anbc', data % emiss_ab(:,:,s % p_e_bc), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('anoc', data % emiss_ab(:,:,s % p_e_oc), &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

            if (config % chem_opt == CHEM_OPT_GOCART) then

              call chem_io_write('ohbg', data % oh_bg, &
                path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
              if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
              call chem_io_write('hobg', data % h2o2_bg, &
                path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
              if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
              call chem_io_write('no3b', data % no3_bg, &
                path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
              if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
              call chem_io_write('ocbb', data % ebu_oc, &
                path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
              if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

            end if

          case (CHEM_OPT_RACM_SOA_VBS)

            call chem_io_write('ao2D', data % aod2d, &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('pm10', data % pm10, &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
            call chem_io_write('pm25', data % pm25, &
              path=trim(config % emi_outname), pos=filepos, de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

          case default

            ! -- no output

        end select

      end do

    end if

  end subroutine chem_output_write

  ! -- restart methods

  subroutine chem_restart_read(verbose, rc)
    logical, optional, intent(in)  :: verbose
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount, item, localpe, tile
    logical :: isVerbose
    character(len=CHEM_MAXSTR)       :: filename
    type(chem_config_type),  pointer :: config   => null()
    type(chem_data_type),    pointer :: data     => null()

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    isVerbose = .false.
    if (present(verbose)) isVerbose = verbose

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount < 1) return

    call chem_comm_get(localpe=localpe, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    filename = "chm_data"

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, tile=tile, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      select case (config % chem_opt)
        case (CHEM_OPT_GOCART)
          do item = 1, config % num_ebu
            call chem_io_read(filename, data % eburn, order="ikj", &
              path=trim(config % restart_inname), de=de, rc=localrc)
            if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
          end do
          if (isVerbose) write(6,'("chem_restart_read: PET:",i4," DE:",i2," tile=",i2," eburn - min/max = "2g16.6)') &
            localpe, de, tile, minval(data % eburn), maxval(data % eburn)
        case default
          ! -- not implemented yet
        end select
    end do

  end subroutine chem_restart_read


  subroutine chem_restart_write(timeStamp, rc)
    character(len=*), optional, intent(in)  :: timeStamp
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount, item
    character(len=CHEM_MAXSTR)       :: filename
    type(chem_config_type),  pointer :: config   => null()
    type(chem_data_type),    pointer :: data     => null()

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount < 1) return

    filename = "chm_data"
    if (present(timeStamp)) then
      filename = trim(timeStamp) // "." // filename
    end if

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      select case (config % chem_opt)
        case (CHEM_OPT_GOCART)
          call chem_io_write(filename, data % eburn, order="ikj", &
            path=trim(config % restart_outname), de=de, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        case default
          ! -- not implemented yet
        end select
    end do

  end subroutine chem_restart_write

end module chem_iodata_mod
