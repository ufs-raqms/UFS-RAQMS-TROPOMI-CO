module raqmschem_methods

  use ESMF
  use NUOPC
  use chem_rc_mod
  use chem_comm_mod
  use chem_types_mod, only : CHEM_MAXSTR,CHEM_KIND_R8
  use raqmschem_model_mod
!  use chem_io_mod
  use raqmschem_iodata_mod
  use raqmschem_config_mod, only : raqmschem_config_type

  implicit none

  public 

contains

  subroutine raqmschem_comp_init(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: deCount
    type(raqmschem_config_type), pointer :: config

    ! -- begin
    if (present(rc)) rc = ESMF_FAILURE
!    write(6,*)'top raqmschem_comp_init',ESMF_FAILURE
    !call flush(6)

    call raqmschem_model_get(deCount=deCount, config=config, rc=localrc)
!     write(6,*)'ajl raqmschemmodelget localrc',localrc,decount
!     call flush(6)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__)) return

    if (deCount > 0) then
!      select case (config % raqmschem_opt)
!        case(CHEM_OPT_GOCART, CHEM_OPT_GOCART_RACM, CHEM_OPT_RACM_SOA_VBS)
!          call gocart_model_init(rc=localrc)
!          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__)) then
!            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Failed to initialize model", &
!              line=__LINE__, file=__FILE__, rcToReturn=rc)
!            return  ! bail out
!          end if
!        case default
!          return
!      end select
    end if

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine raqmschem_comp_init


  subroutine raqmschem_comp_advance(clock, rc)
    use chem_types_mod
    use raqms_model_mod

    type(ESMF_Clock), intent(in) :: clock
    integer,         intent(out) :: rc

    ! -- local variables
    integer                 :: deCount
    integer                 :: julday, yy, mm, dd, h, m, s
    integer(ESMF_KIND_I8)   :: advanceCount
    real(ESMF_KIND_R8)      :: dts
    character(len=CHEM_MAXSTR) :: tStamp
    type(ESMF_Time)         :: currTime
    type(ESMF_TimeInterval) :: timeStep
    type(raqmschem_config_type), pointer :: config
    real(CHEM_KIND_R8), dimension(:,:), pointer :: lat, lon
    integer :: is, ie, js, je, ni, nl
    integer :: localrc,de
    integer :: mype
    integer :: entry
    save entry
    data entry/0/

    ! -- begin
    rc = ESMF_SUCCESS
    entry=entry+1
!    write(6,*)'find ajl top raqmschem_comp_advance',entry
!    call flush(6)

    ! -- check if model is active on this PET, bail out if not
    call raqmschem_model_get(deCount=deCount, config=config, rc=rc)
    if (chem_rc_check(rc, file=__FILE__, line=__LINE__)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Failed to get model info", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    end if

    if (deCount < 1) return
    call chem_comm_get(localpe=mype)
!    write(6,*)'find ajl chem methods ',mype
!    call flush(6)
!    write(200+mype,*)'find ajl top raqms advance'
!    call flush(200+mype)

    ! -- get current time and set model's internal clock
    if(mype.eq.0)then
      call ESMF_ClockPrint(clock, &
      preString="RCHM chem: run(): time step : ", rc=rc)
!    else
!      call ESMF_ClockPrint(clock) 
    endif
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!    if(mype.eq.0)then
!      call ESMF_ClockPrint(clock, options="currTime", &
!      preString="RCHM chem: run(): time stamp: ", rc=rc)
!    else
!      call ESMF_ClockPrint(clock, options="currTime") 
!    endif
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, &
      advanceCount=advanceCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeIntervalGet(timeStep, s_r8=dts, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(currTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
      dayOfYear=julday, timeString=tStamp, rc=rc)
!     write(6,*)'yy',yy,mm,dd,h,m,s,'julday',julday
!     call flush(6)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    if(mype.eq.0)then
     !write(6,*)'araqms clock set',yy,mm,dd,h,m,s,'dts',dts,'adv',advancecount
!    endif
    call chem_model_clock_set(julday=julday, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, dts=dts, &
      advanceCount=int(advanceCount), rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Failed to set model's internal clock", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    end if
!     write(6,*)'find ajl call raqms_model_advance'
!     call flush(6)
        call raqms_model_advance(rc=rc)
!      write(6,*)'find ajl bottom did call raqms_model_advance',rc
!     call flush(6)
        if (chem_rc_check(rc, file=__FILE__, line=__LINE__)) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Failed to advance model", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
        end if
!      write(200+mype,*)'bottom raqms advance',entry
!      call flush(200+mype)

  end subroutine raqmschem_comp_advance


  subroutine raqmschem_comp_finalize(rc)
    integer, intent(out) :: rc

    ! -- begin
    rc = ESMF_SUCCESS
!    write(6,*)'chem_model_destroy'
!    call flush(6)

    call chem_model_destroy()
!    write(6,*)'did chem_model_destroy'
!    call flush(6)

  end subroutine raqmschem_comp_finalize

  !-----------------------------------------------------------------------------

  subroutine raqmschem_comp_connect(stateType, state, fieldNames, rc)
    character(len=*),  intent(in)  :: stateType
    type(ESMF_State),  intent(in)  :: state
    character(len=*),  intent(in)  :: fieldNames(:)
    integer,           intent(out) :: rc

    ! -- begin
    rc = ESMF_RC_NOT_IMPL

    select case (trim(stateType))
      case('import','i')
!        write(6,*)'call raqmschem_comp_import',size(fieldnames),fieldnames
!        call flush(6)
        call raqmschem_comp_import(state, fieldNames, rc)
      case('export','e')
!        write(6,*)'call raqmschem_comp_export'
!        call flush(6)
        call raqmschem_comp_export(state, fieldNames, rc)
!        write(6,*)'call raqmschem_comp_export',rc
!        call flush(6)
      case default
        ! not implemented
    end select

  end subroutine raqmschem_comp_connect


  subroutine raqmschem_comp_export(state, fieldNames, rc)
    type(ESMF_State),               intent(in) :: state
    character(len=*), dimension(:), intent(in) :: fieldNames
    integer, intent(out) :: rc

    ! -- local variables
    type(chem_state_type), pointer :: stateOut
    type(ESMF_Field)               :: field
    integer                        :: item, localDe, localDeCount

    ! -- begin
    rc = ESMF_SUCCESS

    ! -- check if model is active on this PET, bail out if not
    call raqmschem_model_get(deCount=localDeCount, rc=rc)
    if (chem_rc_check(rc, file=__FILE__, line=__LINE__)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Failed to get model info", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    end if

    if (localDeCount < 1) return

    do item = 1, size(fieldNames)

      call ESMF_StateGet(state, field=field, &
        itemName=trim(fieldNames(item)), rc=rc)
!      write(6,*)'esmf_stateget ',trim(fieldnames(item)),rc
!      call flush(6)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail

      do localDe = 0, localDeCount-1

        call raqmschem_model_get(stateOut=stateOut, de=localDe, rc=rc)
        if (chem_rc_check(rc)) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="Failed to retrieve model's export state", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if

        select case (trim(fieldNames(item)))
          case ("raqms_inst_tracer_mass_frac")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateOut % tr3d, rc=rc)
!            write(6,*)'esmf_fieldget ',kind(stateout%tr3d),rc
!            call flush(6)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case default
            ! -- unused field
        end select

      end do
      call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail
    end do

  end subroutine raqmschem_comp_export

  subroutine raqmschem_comp_import(state, fieldNames, rc)
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: fieldNames(:)
    integer,          intent(out) :: rc

    ! -- local variables
    type(chem_state_type), pointer :: stateIn
    type(ESMF_Field)               :: field
    integer                        :: item, localDe, localDeCount

    ! -- begin
    rc = ESMF_SUCCESS

    ! -- check if model is active on this PET, bail out if not
    call raqmschem_model_get(deCount=localDeCount, rc=rc)
    if (chem_rc_check(rc, file=__FILE__, line=__LINE__)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Failed to get model info", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    end if
    if (localDeCount < 1) return

    do item = 1, size(fieldNames)

      call ESMF_StateGet(state, field=field, &
        itemName=trim(fieldNames(item)), rc=rc)
!      write(6,*)'item ',item,'rc',rc,fieldnames(item)
!      call flush(6)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail

      do localDe = 0, localDeCount-1

        call raqmschem_model_get(stateIn=stateIn, de=localDe, rc=rc)
        if (chem_rc_check(rc)) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="Failed to retrieve model's import state", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
!        write(6,*)'item',item,trim(fieldNames(item))
!        call flush(6)

        select case (trim(fieldNames(item)))
          case ("inst_pres_interface")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % pr3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
            call raqmschem_model_set(numIntLayers=size(stateIn % pr3d,dim=3), de=localDe)
          case ("inst_pres_levels")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % prl3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
            call raqmschem_model_set(numModLayers=size(stateIn % prl3d,dim=3), de=localDe)
          case ("inst_temp_levels")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % tk3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("soil_type")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % stype2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_pbl_height")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % pb2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("surface_cell_area")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % area, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_convective_rainfall_amount")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % rc2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_exchange_coefficient_heat_levels")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % exch, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_friction_velocity")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % us2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_geop_interface")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % ph3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_geop_levels")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % phl3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_cldtau")

!            call flush(6)
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % cldtau, rc=rc)
!            write(6,*)'raqmschem cldtau shape',shape(statein%cldtau),'kind', &
!            kind(statein%cldtau)
!            write(6,*)'cldtau rc ',rc
!            call flush(6)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_vort")

!            call flush(6)
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % vort, rc=rc)
!             write(6,*)'raqms vort in',maxval(statein%vort),minval(statein%vort)
!            write(6,*)'raqmschem cldtau shape',shape(statein%cldtau),'kind', &
!            kind(statein%cldtau)
!            write(6,*)'cldtau rc ',rc
!            call flush(6)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_tracer_mass_frac")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % tr3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
            call raqmschem_model_set(numTracers=size(stateIn % tr3d, dim=4), de=localDe)
          case ("inst_omega_levels")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % ws3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
            
              return  ! bail
          case ("inst_rainfall_amount")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % rn2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_soil_moisture_content")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % sm3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
            call raqmschem_model_set(numSoilLayers=size(stateIn % sm3d, dim=3), de=localDe)
          case ("inst_down_sw_flx")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % rsds, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_land_sea_mask")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % slmsk2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_temp_height_surface")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % ts2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_up_sensi_heat_flx")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % hf2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_lwe_snow_thickness")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % snwdph2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("vegetation_type")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % vtype2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_vegetation_area_frac")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % vfrac2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_zonal_wind_levels")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % us3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_merid_wind_levels")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % vs3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("inst_surface_roughness")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % zorl2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("cm_drag")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % cmdrag, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
!            write(6,*)'cmdrag',shape(stateIn % cmdrag),kind(stateIn % cmdrag)
!            call flush(6)
          case ("ktop_kbot_cnv")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % cktop_kbot_cnv, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
!            write(6,*)'cktop_kbot_cnv',shape(statein%cktop_kbot_cnv),kind(stateIn% cktop_kbot_cnv)
!            call flush(6)
          case ("inst_spec_humid_conv_tendency_levels")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stateIn % dqdt, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case default
            ! -- unused field
        end select
!        write(6,*)'did ',item,trim(fieldNames(item))
!        call flush(6)
      end do
    end do

  end subroutine raqmschem_comp_import

  !-----------------------------------------------------------------------------

  subroutine fieldPrintMinMax(field, vm, global, rc)
    type(ESMF_Field),        intent(in)  :: field
    type(ESMF_VM), optional, intent(in)  :: vm
    logical,       optional, intent(in)  :: global
    integer,       optional, intent(out) :: rc

    ! local variables
    type(ESMF_VM)               :: localVM
    real(ESMF_KIND_R8), pointer :: fp1d(:), fp2d(:,:), fp3d(:,:,:), fp4d(:,:,:,:)
    real(ESMF_KIND_R8)          :: fieldMaxValue, fieldMinValue, maxValue, minValue
    real(ESMF_KIND_R8)          :: globalMaxValue(1), globalMinValue(1)
    integer                     :: localDe, localDeCount, localPet, localrc, rank
    logical                     :: addGlobal
    character(len=ESMF_MAXSTR)  :: fieldName

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    addGlobal = .false.
    if (present(global)) addGlobal = global

    if (present(vm)) then
      localVM = vm
    else
      call ESMF_VMGetCurrent(localVM, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return  ! bail out
    end if

    call ESMF_VMGet(localVM, localPet=localPet, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return  ! bail out

    call ESMF_FieldGet(field, rank=rank, localDeCount=localDeCount, &
      name=fieldName, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return  ! bail out

    fieldMinValue = huge(1.0_ESMF_KIND_R8)
    fieldMaxValue = -fieldMinValue

    do localDe = 0, localDeCount - 1
      select case(rank)
        case(1)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp1d, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
          minValue = minval(fp1d)
          maxValue = maxval(fp1d)
        case(2)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp2d, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
          minValue = minval(fp2d)
          maxValue = maxval(fp2d)
        case(3)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp3d, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
          minValue = minval(fp3d)
          maxValue = maxval(fp3d)
        case(4)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp4d, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
          minValue = minval(fp4d)
          maxValue = maxval(fp4d)
        case default
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Field rank not implemented.", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=localrc)
          return ! bail out
      end select
      fieldMinValue = min(fieldMinValue, minValue)
      fieldMaxValue = max(fieldMaxValue, maxValue)
!      write(6,'(a,":",i0,2x,"DE: ",i0,2x,a," - raq checking  - min/max = ",2g16.6)') 'PET', &
!         localPet, localDe, trim(fieldName), minValue, maxValue
!      call flush(6)
    end do

    if (addGlobal) then

      globalMinValue(1) = 0._ESMF_KIND_R8
      globalMaxValue(1) = 0._ESMF_KIND_R8

      call ESMF_VMReduce(localVM, (/ fieldMinValue /), globalMinValue, 1, &
        reduceflag=ESMF_REDUCE_MIN, rootPet=0, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return  ! bail out
        return  ! bail out
      call ESMF_VMReduce(localVM, (/ fieldMaxValue /), globalMaxValue, 1, &
        reduceflag=ESMF_REDUCE_MAX, rootPet=0, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return  ! bail out

!      if (localPet == 0) then
!         write(6,'(a,":",a," - checking  - min/max = ",2g16.6)') 'Field', &
!           trim(fieldName), globalMinValue, globalMaxValue
!      end if

    end if

  end subroutine fieldPrintMinMax
! subroutine raqms_model_advance2(rc)
! use machine,only : kind_phys
! use ozne_def,only : oz_lat,kozpl,kozc
! use chem_types_mod
! implicit none
! type(chem_state_type),  pointer :: stateIn, stateOut
! integer rc,localrc,de,is,ie,js,je,ni,nl,i,decount
! real(CHEM_KIND_R8), dimension(:,:), pointer :: lat, lon 
! type(raqmschem_config_type), pointer :: config
! type(chem_data_type),   pointer :: data
! write(6,*)'oz_lat',shape(oz_lat),kozpl,kozc
! write(6,*)'max',maxval(oz_lat)
! call flush(6)
! call raqmschem_model_get(deCount=deCount, rc=localrc)
! if (deCount < 1) return
! do de = 0, deCount-1
!     call raqmschem_model_get(de=de, config=config, data=data, &
!       stateIn=stateIn, stateOut=stateOut, rc=localrc)
!  call chem_model_domain_get(de=de, ids=is, ide=ie, jds=js, jde=je, ni=ni, nl=nl, &
!       lon=lon, lat=lat, rc=localrc) 
!!  write(6,*)'ajl ni',ni,nl,'tr3d',shape(statein%tr3d)
!!  write(6,*)'ajl is',is,ie,js,je
!!  do i=1,4
!    !write(6,*)'ajl i',i,maxval(statein%tr3d(:,:,:,i)),minval(statein%tr3d(:,:,:,i))
!!  end do
! end do
! return
! end subroutine raqms_model_advance2


  !-----------------------------------------------------------------------------

end module raqmschem_methods
