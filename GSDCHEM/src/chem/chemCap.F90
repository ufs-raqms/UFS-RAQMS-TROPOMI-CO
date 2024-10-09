module CHM
  !-----------------------------------------------------------------------------
  ! CHM Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, inheritModel => SetServices

  use chem_methods
  
  implicit none

  ! -- import fields
  integer, parameter :: importFieldCount = 26
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "inst_pres_interface                  ", &
      "inst_pres_levels                     ", &
      "inst_geop_interface                  ", &
      "inst_geop_levels                     ", &
      "inst_temp_levels                     ", &
      "inst_zonal_wind_levels               ", &
      "inst_merid_wind_levels               ", &
      "inst_omega_levels                    ", &
      "inst_tracer_mass_frac                ", &
      "soil_type                            ", &
      "inst_pbl_height                      ", &
      "surface_cell_area                    ", &
      "inst_convective_rainfall_amount      ", &
      "inst_exchange_coefficient_heat_levels", &
      "inst_spec_humid_conv_tendency_levels ", &
      "inst_friction_velocity               ", &
      "inst_rainfall_amount                 ", &
      "inst_soil_moisture_content           ", &
      "inst_down_sw_flx                     ", &
      "inst_land_sea_mask                   ", &
      "inst_temp_height_surface             ", &
      "inst_up_sensi_heat_flx               ", &
      "inst_lwe_snow_thickness              ", &
      "vegetation_type                      ", &
      "inst_vegetation_area_frac            ", &
      "inst_surface_roughness               "  &
    /)
  ! -- export fields
  integer, parameter :: exportFieldCount = 5
  character(len=*), dimension(exportFieldCount), parameter :: &
    exportFieldNames = (/ &
      "inst_tracer_mass_frac                ", &
      "inst_tracer_up_surface_flx           ", &
      "inst_tracer_down_surface_flx         ", &
      "inst_tracer_clmn_mass_dens           ", &
      "inst_tracer_anth_biom_flx            "  &
    /)

  ! -- verbosity
  integer :: verbosity

  private
  
  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! begin
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, inheritModel, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Provide InitializeP0 to switch from default IPDv00 to IPDv03
    call ESMF_GridCompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(model, specLabel=label_DataInitialize, &
      specRoutine=DataInitialize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_Finalize, &
      specRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    character(len=5)     :: value

    ! begin
    rc = ESMF_SUCCESS

    ! get component verbosity
    call ESMF_AttributeGet(model, name="Verbosity", value=value, &
      defaultValue="min", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! convert value to verbosity
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min ","high", "max "/), &
      specialValueList=(/0,255,255/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! write to log files
    call ESMF_LogWrite("CHM: Verbosity = " // trim(value), &
      ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(model, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! begin
    rc = ESMF_SUCCESS

    ! -- advertise imported fields
    if (importFieldCount > 0) then
      call NUOPC_Advertise(importState, importFieldNames, &
        TransferOfferGeomObject="cannot provide", &
        SharePolicyField="share", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    ! -- advertise exported fields
    if (exportFieldCount > 0) then
      call NUOPC_Advertise(exportState, exportFieldNames, &
        TransferOfferGeomObject="cannot provide", &
        SharePolicyField="share", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine DataInitialize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)              :: importState, exportState
    type(ESMF_Field)              :: field
    type(ESMF_Clock)              :: clock
    type(ESMF_Alarm)              :: alarm
    type(ESMF_Grid)               :: grid
    type(ESMF_VM)                 :: vm
    type(ESMF_GeomType_flag)      :: geomtype
    type(ESMF_DistGrid)           :: distgrid
    type(ESMF_Array)              :: array
    type(ESMF_Config)             :: cf
    integer                       :: de, item, localrc, localDe, tile
    integer                       :: comm
    real(ESMF_KIND_R8), dimension(:,:), pointer :: coord

    integer :: dimCount, tileCount, deCount, localDeCount
    integer, dimension(:),   allocatable :: deToTileMap, localDeToDeMap
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile
    integer, dimension(:,:), allocatable :: computationalLBound, computationalUBound

    logical                    :: restart
    integer                    :: restart_interval
    integer                    :: yy, mm, dd, h, m
    integer(ESMF_KIND_I8)      :: advanceCount
    real(ESMF_KIND_R8)         :: dts
    type(ESMF_Time)            :: startTime, currTime
    type(ESMF_TimeInterval)    :: TimeStep
    character(len=255)         :: msgString

    character(len=*), parameter :: config_file = "model_configure"

    ! begin
    rc = ESMF_SUCCESS

    ! -- initialize chemistry model
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- check if import fields are defined
    if (importFieldCount < 1) then 
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="This component requires import fields to be defined.", &
        line=__LINE__, file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- check if export fields are defined
    if (exportFieldCount < 1) then 
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="This component requires export fields to be defined.", &
        line=__LINE__, file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get coordinates from Grid object
    ! assume all fields on same grid
    ! use first field 
    call ESMF_StateGet(importState, field=field, &
      itemName=trim(importFieldNames(1)), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_FieldGet(field, geomtype=geomtype, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (geomtype == ESMF_GEOMTYPE_GRID) then
      call ESMF_FieldGet(field, grid=grid, array=array, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_ArrayGet(array, deCount=deCount, dimCount=dimCount, &
        tileCount=tileCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(minIndexPDe(dimCount, deCount), maxIndexPDe(dimCount, deCount),  &
        minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount), &
        computationalLBound(dimCount, localDeCount), computationalUBound(dimCount, localDeCount), &
        deToTileMap(deCount), localDeToDeMap(localDeCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_ArrayGet(array, distgrid=distgrid, &
        deToTileMap=deToTileMap, localDeToDeMap=localDeToDeMap, &
        computationalLBound=computationalLBound, &
        computationalUBound=computationalUBound, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_DistGridGet(distgrid, &
        minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe, &
        minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! -- init chemistry model on local DEs
      call chem_model_create(deCount=localDeCount, rc=rc)
      if (chem_rc_check(rc)) then
        call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
          msg="Failed to initialize chemistry model for localDeCount", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return  ! bail out
      end if

      do localDe = 0, localDeCount-1
        de   = localDeToDeMap(localDe+1) + 1
        tile = deToTileMap(de)

        if (btest(verbosity,0)) then
          write(msgString,'("CHM: localDe: ",i4," DE: ",i4, " tile=",i2," minIndexPDe=",2i4,2x," maxIndexPDe=",2i4," minIndexPTile=",2i4," maxIndexPTile=",2i4,4i4)') &
            localDe, de-1, tile, minIndexPDe(:,de), maxIndexPDe(:,de), minIndexPTile(:,tile), maxIndexPTile(:,tile), &
            computationalLBound(:,localDe+1), computationalUBound(:,localDe+1)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        end if

        ! -- set model domains for local DEs
        call chem_model_domain_set(minIndexPDe=minIndexPDe(:,de), maxIndexPDe=maxIndexPDe(:,de), &
          minIndexPTile=minIndexPTile(:,tile), maxIndexPTile=maxIndexPTile(:,tile), &
          minIndexLocal=computationalLBound(:,localDe+1), maxIndexLocal=computationalUBound(:,localDe+1), &
          tile=deToTileMap(de), tileCount=tileCount, de=localDe, rc=rc)
        if (chem_rc_check(rc)) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="Failed to set domain for chemistry model", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if

        ! -- get local coordinate arrays
        do item = 1, 2
          call ESMF_GridGetCoord(grid, coordDim=item, staggerloc=ESMF_STAGGERLOC_CENTER, &
            localDE=localDe, farrayPtr=coord, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call chem_model_domain_coord_set(item, coord, de=localDe, rc=rc)
          if (chem_rc_check(rc)) then
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
              msg="Failed to set coordinates for chemistry model", &
              line=__LINE__, &
              file=__FILE__, &
              rcToReturn=rc)
            return  ! bail out
          end if
        end do

      end do
      deallocate(minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
        computationalLBound, computationalUBound, &
        deToTileMap, localDeToDeMap, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="Imported fields can only be defined on Grid objects.", &
        line=__LINE__, file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- init model communication subsystem over ESMF communicator
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_VMGet(vm, mpiCommunicator=comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- initialize model after creation to setup correct communication
    call chem_model_init(comm=comm, isolate=.true., rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize chemistry model", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- initialize model I/O
    call chem_io_init(verbose=btest(verbosity,0), rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize I/O model subsystem", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- initialize internal clock
    ! -- get clock information
    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, &
      timeStep=timeStep, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    ! -- set restart flag
    restart = (currTime > startTime)
    ! -- reset advanceCount if model is restarting
    if (restart) then
      advanceCount = (currTime - startTime) / timeStep
      call ESMF_ClockSet(clock, advanceCount=advanceCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if
    ! -- get forecast starting time
    call ESMF_TimeGet(currTime, yy=yy, mm=mm, dd=dd, h=h, m=m, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    ! -- get time step
    call ESMF_TimeIntervalGet(timeStep, s_r8=dts, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    ! -- set internal clock
    call chem_model_clock_create(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to create model clock", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    call chem_model_clock_set(yy=yy, mm=mm, dd=dd, h=h, m=m, &
      dts=real(dts,kind=CHEM_KIND_R8), tz=0, rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize model clock", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- read-in model configuration options
    cf = ESMF_ConfigCreate(rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    call ESMF_ConfigLoadFile(cf, config_file, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    call ESMF_ConfigGetAttribute(cf, restart_interval, &
      label="restart_interval:", default=0, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! -- if restart_interval > 0, create alarm to write restart files
    if (restart_interval > 0) then
      call ESMF_TimeIntervalSet(timeStep, h=restart_interval, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      alarm = ESMF_AlarmCreate(clock, ringInterval=timeStep, &
        ringTime=(startTime + timeStep), name="chemistry_restart", rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

    ! -- read-in namelist options and setup internal parameters
    call chem_model_config_init(restart=restart, rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize model configuration", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- allocate memory for background fields
    call chem_backgd_init(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize I/O model storage", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- read-in background fields
    call chem_backgd_read(verbose=btest(verbosity,0), rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to read background/emission values", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- connect import fields to model
    ! -- this can be done only once since remote fields are accessed by reference
    call chem_comp_connect('import', importState, importFieldNames, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- connect export fields to model
    ! -- this can be done only once since remote fields are accessed by reference
    call chem_comp_connect('export', exportState, exportFieldNames, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- allocate memory for internal workspace
    call chem_buffer_init(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize I/O model buffers", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- restore internal buffers after restart
    if (restart) then
      call chem_restart_read(verbose=btest(verbosity,0), rc=rc)
      if (chem_rc_check(rc)) then
        call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
          msg="Failed to read restart data", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return  ! bail out
      end if
    end if

    ! -- initialize internal component (GOCART)
    call chem_comp_init(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- allocate arrays for tracer output
    call chem_output_init(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to allocate output arrays", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! indicate that data initialization is complete (breaking out of init-loop)
    call NUOPC_CompAttributeSet(model, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine DataInitialize

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState
    type(ESMF_Time)               :: currTime
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Field)              :: field
    type(ESMF_VM)                 :: vm
    type(ESMF_Alarm), pointer     :: alarmList(:)
    integer                       :: alarmCount, item, localrc
    integer                       :: yy, mm, dd, h, m, s
    character(len=ESMF_MAXSTR)    :: alarmName
    character(len=15)             :: timeStamp
    logical :: chemoutput=.false. ! ajl
    integer mype                  ! ajl

    ! begin
    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    call chem_comm_get(localpe=mype) ! ajl
    if(mype==0)then    
      call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing CHM from: ", rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    endif
    
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    if(mype==0)then    
     call ESMF_TimePrint(currTime + timeStep, &
      preString="---------------------> to: ", rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    endif

    ! print field diagnostics
    if (btest(verbosity,0)) then
      call ESMF_GridCompGet(model, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      do item = 1, importFieldCount
        call ESMF_StateGet(importState, field=field, &
          itemName=trim(importFieldNames(item)), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call fieldPrintMinMax(field, vm=vm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end do
    end if

    ! advance model
    call chem_comp_advance(clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! write output if it is time
    if(chemoutput)then
    call chem_output_write(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to write chemistry output", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    endif

    ! write restart file if requested
    call ESMF_ClockGetAlarmList(clock, ESMF_ALARMLIST_NEXTRINGING, &
      alarmCount=alarmCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (alarmCount > 0) then
      allocate(alarmList(alarmCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_ClockGetAlarmList(clock, ESMF_ALARMLIST_NEXTRINGING, &
        alarmList=alarmList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      do item = 1, alarmCount
        call ESMF_AlarmGet(alarmList(item), name=alarmName, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (trim(alarmName) == "chemistry_restart") then
          call ESMF_TimeGet(currTime + timeStep, yy=yy, mm=mm, dd=dd, &
            h=h, m=m, s=s, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          write(timeStamp, '(i4.4,2i2.2,".",3i2.2)') yy, mm, dd, h, m, s
          call chem_restart_write(timeStamp=timeStamp, rc=rc)
          if (chem_rc_check(rc)) then
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
              msg="Failed to write chemistry restart file", &
              line=__LINE__, &
              file=__FILE__, &
              rcToReturn=rc)
            return  ! bail out
          end if
        end if
      end do
      deallocate(alarmList, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

  end subroutine ModelAdvance

  !-----------------------------------------------------------------------------

  subroutine ModelFinalize(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc

    ! begin
    rc = ESMF_SUCCESS

    ! finalize model
    call chem_comp_finalize(rc)

  end subroutine ModelFinalize

  !-----------------------------------------------------------------------------

end module CHM
