module chem_model_mod

  use chem_rc_mod
  use chem_types_mod
  use chem_config_mod
  use chem_clock_mod
  use chem_data_mod,     only : chem_data_type, chem_data_destroy
  use chem_domain_mod,   only : chem_domain_type
  use chem_state_mod,    only : chem_state_type
  use chem_iolayout_mod, only : chem_iolayout_type
  use chem_comm_mod

  implicit none

  type chem_model_type
    type(chem_clock_type),   pointer :: clock   => null()
    type(chem_config_type),  pointer :: config  => null()
    type(chem_domain_type)   :: domain
    type(chem_state_type)    :: stateIn, stateOut
    type(chem_data_type)     :: data
    type(chem_iolayout_type) :: iolayout
  end type chem_model_type

  type(chem_model_type), dimension(:), allocatable, target :: chem_model

  private

  public :: CHEM_OPT_GOCART
  public :: CHEM_OPT_GOCART_RACM
  public :: CHEM_OPT_RACM_SOA_VBS

  public :: chem_model_type
  ! -- also provide subtypes
  public :: chem_clock_type
  public :: chem_config_type
  public :: chem_domain_type
  public :: chem_state_type
  public :: chem_data_type
  public :: chem_iolayout_type

  public :: chem_model_create
  public :: chem_model_destroy
  public :: chem_model_init
  public :: chem_model_get
  public :: chem_model_set
  public :: chem_model_clock_create
  public :: chem_model_clock_get
  public :: chem_model_clock_set
  public :: chem_model_config_init
  public :: chem_model_domain_get
  public :: chem_model_domain_set
  public :: chem_model_domain_coord_set

contains

  subroutine chem_model_create(deCount, rc)

    integer, intent(in),  optional :: deCount
    integer, intent(out), optional :: rc

    !-- local variables
    integer :: localDeCount, localrc

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localDeCount = 1
    if (present(deCount)) localDeCount = deCount

    if (localDeCount > 0) then
      allocate(chem_model(0:localDeCount-1), stat=localrc)
      if (chem_rc_test((localrc /= 0), msg="Failure to create model", &
          file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    
  end subroutine chem_model_create


  subroutine chem_model_destroy(rc)

    integer, intent(out), optional :: rc

    !-- local variables
    integer :: localrc
    integer :: de

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    if (allocated(chem_model)) then
      do de = 1, size(chem_model)-1
        call chem_data_destroy(chem_model(de) % data, rc=localrc)
        if (chem_rc_test((localrc /= 0), &
          msg="Failure to allocate model data memory", &
          file=__FILE__, line=__LINE__, rc=rc)) return
        nullify(chem_model(de) % config)
        nullify(chem_model(de) % clock)
      end do
      de = 0
      call chem_data_destroy(chem_model(de) % data, rc=localrc)
      if (associated(chem_model(de) % config)) then
        if (associated(chem_model(de) % config % species)) then
          deallocate(chem_model(de) % config % species, stat=localrc)
          if (chem_rc_test((localrc /= 0), &
            msg="Failure to allocate model species memory", &
            file=__FILE__, line=__LINE__, rc=rc)) return
          nullify(chem_model(de) % config % species)
        end if
        deallocate(chem_model(de) % config, stat=localrc)
        if (chem_rc_test((localrc /= 0), &
          msg="Failure to allocate model config memory", &
          file=__FILE__, line=__LINE__, rc=rc)) return
        nullify(chem_model(de) % config)
      end if
      if (associated(chem_model(de) % clock)) then
        deallocate(chem_model(de) % clock, stat=localrc)
        if (chem_rc_test((localrc /= 0), &
          msg="Failure to allocate model clock memory", &
          file=__FILE__, line=__LINE__, rc=rc)) return
        nullify(chem_model(de) % clock)
      end if
      deallocate(chem_model, stat=localrc)
      if (chem_rc_test((localrc /= 0), msg="Failure to allocate model memory", &
          file=__FILE__, line=__LINE__, rc=rc)) return
    end if

  end subroutine chem_model_destroy


  subroutine chem_model_init(rc, comm, isolate)
    integer, optional, intent(in)  :: comm
    logical, optional, intent(in)  :: isolate
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: color, de, deCount, localrc
    integer :: modelComm

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_local_get(deCount=deCount, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- initialize communications, if needed
    call chem_comm_init(comm=comm, isolate=isolate, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to initialize communications", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- create new communicators to separate PETs with model from PETs with no model
    color = 0
    if (deCount > 0) color = 1

    call chem_comm_create(modelComm, color, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- set new communicator
    call chem_comm_set(comm=modelComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- store new communicator in model's data structure
    do de = 0, deCount-1
      call chem_model_set(de=de, modelComm=modelComm, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    end do

  end subroutine chem_model_init


  subroutine chem_model_clock_create(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer                        :: localrc
    integer                        :: de, deCount
    type(chem_model_type), pointer :: model

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    nullify(model)

    ! -- NOTE: config is allocated only on DE 0 on each PET
    ! --       config points to DE 0 on other DEs on the same PET
    deCount = 0
    call chem_model_local_get(de=0, deCount=deCount, model=model, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount > 0) then
      if (.not.associated(model % clock)) allocate(model % clock)
      call chem_clock_set(model % clock, &
        julday=0, yy=0, mm=0, dd=0, h=0, m=0, s=0, dts=0._CHEM_KIND_R8, &
        advanceCount=0, rc=localrc)
      if (chem_rc_check(localrc, msg="Failure to set model clock", &
        file=__FILE__, line=__LINE__, rc=rc)) return
      ! -- populate pointers on other local DEs
      do de = 1, deCount-1
        call chem_model_set(de=de, clock=model % clock, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      end do
    end if

  end subroutine chem_model_clock_create

  subroutine chem_model_clock_get(julday, yy, mm, dd, h, m, s, tz, dts, &
    advanceCount, tStamp, rc)
    integer,            optional, intent(out) :: julday
    integer,            optional, intent(out) :: yy, mm, dd, h, m, s, tz
    real(CHEM_KIND_R8), optional, intent(out) :: dts
    integer,            optional, intent(out) :: advanceCount
    character(len=12),  optional, intent(out) :: tStamp
    integer,            optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: deCount
    type(chem_clock_type), pointer :: clock

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=0, deCount=deCount, clock=clock, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount > 0) then
      call chem_clock_get(clock, julday=julday, yy=yy, mm=mm, dd=dd, &
        h=h, m=m, s=s, tz=tz, dts=dts, advanceCount=advanceCount, &
        tStamp=tStamp, rc=localrc)
      if (chem_rc_check(localrc, msg="Failure to get model clock info", &
        file=__FILE__, line=__LINE__, rc=rc)) return
    else
      if (present(julday)) julday = 0
      if (present(yy))     yy     = 0
      if (present(mm))     mm     = 0
      if (present(dd))     dd     = 0
      if (present(h))      h      = 0
      if (present(m))      m      = 0
      if (present(s))      s      = 0
      if (present(tz))     tz     = 0
      if (present(dts))    dts    = 0._CHEM_KIND_R8
      if (present(advanceCount)) advanceCount = 0
      if (present(tStamp)) tStamp = ""
    end if

  end subroutine chem_model_clock_get


  subroutine chem_model_clock_set(julday, yy, mm, dd, h, m, s, tz, dts, advanceCount, rc)
    integer,            optional, intent(in)  :: julday
    integer,            optional, intent(in)  :: yy, mm, dd, h, m, s, tz
    real(CHEM_KIND_R8), optional, intent(in)  :: dts
    integer,            optional, intent(in)  :: advanceCount
    integer,            optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: deCount
    type(chem_clock_type), pointer :: clock

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=0, deCount=deCount, clock=clock, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount > 0) then
      call chem_clock_set(clock, julday=julday, yy=yy, mm=mm, dd=dd, &
        h=h, m=m, s=s, tz=tz, dts=dts, advanceCount=advanceCount, rc=localrc)
      if (chem_rc_check(localrc, msg="Failure to set model clock", &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if

  end subroutine chem_model_clock_set


  subroutine chem_model_config_init(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount
    type(chem_model_type), pointer :: model

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    ! -- NOTE: config is allocated only on DE 0 on each PET
    ! --       config points to DE 0 on other DEs on the same PET
    deCount = 0
    call chem_model_local_get(de=0, deCount=deCount, model=model, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount > 0) then
      if (.not.associated(model % config)) allocate(model % config)
      ! -- read in namelist options
      call chem_config_read(model % config, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      ! -- setup internal parameters based on input options
      call chem_config_control_init(model % config, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      ! -- setup internal indexes for chemical species
      call chem_config_species_init(model % config, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      ! -- populate pointers on other local DEs
      do de = 1, deCount-1
        call chem_model_set(de=de, config=model % config, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      end do
    end if

  end subroutine chem_model_config_init


  subroutine chem_model_domain_set(minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
    minIndexLocal, maxIndexLocal, tile, tileCount, de, rc)

    integer, dimension(2), optional, intent(in)  :: minIndexPDe,   maxIndexPDe
    integer, dimension(2), optional, intent(in)  :: minIndexPTile, maxIndexPTile
    integer, dimension(2), optional, intent(in)  :: minIndexLocal, maxIndexLocal
    integer,               optional, intent(in)  :: tile, tileCount
    integer,               optional, intent(in)  :: de
    integer,               optional, intent(out) :: rc

    !-- local variables
    integer                        :: deCount, localrc
    type(chem_model_type), pointer :: model

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_local_get(de=de, deCount=deCount, model=model, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount < 1) return

    model % domain % tile = 0
    if (present(tile)) model % domain % tile = tile

    model % domain % tileCount = 0
    if (present(tileCount)) model % domain % tileCount = tileCount

    if (present(minIndexPDe)) then
      model % domain % ids = minIndexPDe(1)
      model % domain % jds = minIndexPDe(2)
    end if
    if (present(maxIndexPDe)) then
      model % domain % ide = maxIndexPDe(1)
      model % domain % jde = maxIndexPDe(2)
    end if
    if (present(minIndexPTile)) then
      model % domain % its = minIndexPTile(1)
      model % domain % jts = minIndexPTile(2)
    end if
    if (present(maxIndexPTile)) then
      model % domain % ite = maxIndexPTile(1)
      model % domain % jte = maxIndexPTile(2)
    end if
    if (present(minIndexLocal)) then
      model % domain % ims = minIndexLocal(1)
      model % domain % jms = minIndexLocal(2)
    end if
    if (present(maxIndexLocal)) then
      model % domain % ime = maxIndexLocal(1)
      model % domain % jme = maxIndexLocal(2)
    end if

  end subroutine chem_model_domain_set


  subroutine chem_model_domain_coord_set(coordDim, coord, de, rc)

    integer,            intent(in)  :: coordDim
    real(CHEM_KIND_R8), pointer     :: coord(:,:)
    integer, optional,  intent(in)  :: de
    integer, optional,  intent(out) :: rc

    !-- local variables
    integer                        :: deCount, localrc
    type(chem_model_type), pointer :: model

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_local_get(de=de, deCount=deCount, model=model, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount < 1) return

    select case (coordDim)
      case(1)
        model % domain % lon => coord
      case(2)
        model % domain % lat => coord
      case default
        call chem_rc_set(CHEM_RC_FAILURE, &
          msg="coordDim can only be 1 or 2.", &
          file=__FILE__, line=__LINE__, rc=rc)
        return
    end select

  end subroutine chem_model_domain_coord_set


  subroutine chem_model_domain_get(de, &
   ids, ide, jds, jde,  &
   its, ite, jts, jte,  &
   ims, ime, jms, jme,  &
   ni,  nl,  nt,  tile, &
   lon, lat, rc)
    integer, optional, intent(in)  :: de
    integer, optional, intent(out) :: ids, ide, jds, jde
    integer, optional, intent(out) :: its, ite, jts, jte
    integer, optional, intent(out) :: ims, ime, jms, jme
    integer, optional, intent(out) :: ni, nl, nt
    integer, optional, intent(out) :: tile
    real(CHEM_KIND_R8), optional, pointer :: lon(:,:)
    real(CHEM_KIND_R8), optional, pointer :: lat(:,:)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer                        :: localrc
    type(chem_model_type), pointer :: model

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_local_get(de=de, model=model, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (present(ids)) ids = model % domain % ids
    if (present(ide)) ide = model % domain % ide
    if (present(jds)) jds = model % domain % jds
    if (present(jde)) jde = model % domain % jde
    if (present(its)) its = model % domain % its
    if (present(ite)) ite = model % domain % ite
    if (present(jts)) jts = model % domain % jts
    if (present(jte)) jte = model % domain % jte
    if (present(ims)) ims = model % domain % ims
    if (present(ime)) ime = model % domain % ime
    if (present(jms)) jms = model % domain % jms
    if (present(jme)) jme = model % domain % jme
    if (present(ni))  ni  = model % domain % ni
    if (present(nl))  nl  = model % domain % nl
    if (present(nt))  nt  = model % domain % nt
    if (present(tile)) tile = model % domain % tile
    if (present(lon)) lon => model % domain % lon
    if (present(lat)) lat => model % domain % lat

  end subroutine chem_model_domain_get


  subroutine chem_model_set(de, numIntLayers, numModLayers, numSoilLayers, numTracers, modelComm, tileComm, &
    localIOflag, clock, config, rc)

    integer, optional, intent(in)  :: de
    integer, optional, intent(in)  :: numIntLayers
    integer, optional, intent(in)  :: numModLayers
    integer, optional, intent(in)  :: numSoilLayers
    integer, optional, intent(in)  :: numTracers
    integer, optional, intent(in)  :: modelComm
    integer, optional, intent(in)  :: tileComm
    logical, optional, intent(in)  :: localIOflag
    type(chem_clock_type),  optional, pointer :: clock
    type(chem_config_type), optional, pointer :: config
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer                        :: deCount, localrc
    type(chem_model_type), pointer :: model
    
    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_local_get(de=de, deCount=deCount, model=model, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount > 0) then
      if (present(numIntLayers)) model % domain % ni = numIntLayers
      if (present(numModLayers)) model % domain % nl = numModLayers
      if (present(numSoilLayers)) model % domain % ns = numSoilLayers
      if (present(numTracers))   model % domain % nt = numTracers
      if (present(modelComm))    model % iolayout % modelComm = modelComm
      if (present(tileComm))     model % iolayout % tileComm  = tileComm
      if (present(localIOflag))  model % iolayout % localIOflag = localIOflag
      if (present(config))       model % config => config
      if (present(clock))        model % clock  => clock
    end if

  end subroutine chem_model_set


  subroutine chem_model_get(de, deCount, stateIn, stateOut, clock, config, &
    data, domain, tile, tileCount, tileComm, modelComm, localIOflag, rc)

    integer,               optional,  intent(in)  :: de
    integer,               optional,  intent(out) :: deCount
    type(chem_state_type), optional,  pointer     :: stateIn, stateOut
    type(chem_clock_type), optional,  pointer     :: clock
    type(chem_config_type),optional,  pointer     :: config
    type(chem_data_type),  optional,  pointer     :: data
    type(chem_domain_type),optional,  pointer     :: domain
    integer,               optional,  intent(out) :: tile
    integer,               optional,  intent(out) :: tileCount, tileComm
    integer,               optional,  intent(out) :: modelComm
    logical,               optional,  intent(out) :: localIOflag
    integer,               optional,  intent(out) :: rc

    !-- local variables
    integer                        :: localrc, localDeCount
    type(chem_model_type), pointer :: model
    
    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localDeCount = 0
    call chem_model_local_get(de=de, deCount=localDeCount, model=model, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure to retrieve model on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (present(deCount)) deCount = localDeCount

    if (localDeCount > 0) then
      if (present(stateIn))     stateIn     => model % stateIn
      if (present(stateOut))    stateOut    => model % stateOut
      if (present(clock))       clock       => model % clock
      if (present(config))      config      => model % config
      if (present(data))        data        => model % data
      if (present(domain))      domain      => model % domain
      if (present(tile))        tile        =  model % domain % tile
      if (present(tileCount))   tileCount   =  model % domain % tileCount
      if (present(modelComm))   modelComm   =  model % iolayout % modelComm
      if (present(tileComm))    tileComm    =  model % iolayout % tileComm
      if (present(localIOflag)) localIOflag =  model % iolayout % localIOflag
    end if

  end subroutine chem_model_get


  subroutine chem_model_local_get(de, deCount, model, rc)

    integer,               optional,  intent(in)  :: de
    integer,               optional,  intent(out) :: deCount
    type(chem_model_type), optional,  pointer     :: model
    integer,               optional,  intent(out) :: rc

    !-- local variables
    integer :: localDe

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    if (present(model)) nullify(model)

    localDe = 0
    if (present(de)) localDe = de

    if (chem_rc_test((localDe < 0), &
      msg="DE must be >= 0", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (allocated(chem_model)) then
      if (present(deCount)) deCount = size(chem_model)
      if (present(model)) then
        if (localDe <= ubound(chem_model, dim=1)) then
          model => chem_model(localDe)
        else
          call chem_rc_set(CHEM_RC_FAILURE, msg="Model undefined on local DE", &
            file=__FILE__, line=__LINE__, rc=rc)
          return
        end if
      end if
    else
      if (present(deCount)) deCount = 0
    end if

  end subroutine chem_model_local_get


end module chem_model_mod
