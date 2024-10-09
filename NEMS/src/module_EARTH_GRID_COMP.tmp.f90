




!-----------------------------------------------------------------------
!
      MODULE module_EARTH_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the EARTH component.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  2010-03-24  Black - Created Earth component module.
!  2010-04     Yang  - Added Ensemble capability.
!  2011-05-11  Theurich & Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!  2011-10-04  Yang - Modified for using the ESMF 5.2.0r library.
!  2012-02     Tripp - Added ESMF superstructure to support an OCN model
!  2013-06     Theurich - Reworked OCN dependency to be NUOPC based
!  2013-07     Theurich - Macro based ESMF error handling
!-----------------------------------------------------------------------
!
!***  The EARTH component lies in the hierarchy seen here:
!
!          Main program
!               |
!               |
!          NEMS component
!               |     |________________________.
!               |                              |
!          EARTH component        Ensemble Coupler component
!              /|!             / | !          ATM/OCN/ICE/WAV/LND/IPM/HYD .. components
!          |    |   |
!          |    |   (CICE, etc.)
!          |    |
!          |    (MOM5, MOM6, HYCOM, POM, etc.)
!          |
!          CORE component (GSM, NMM, FV3, etc.)
!
!-----------------------------------------------------------------------
!
      USE ESMF

      use NUOPC
      use NUOPC_Driver, &
        Driver_routine_SS             => SetServices, &
        Driver_label_SetModelServices => label_SetModelServices, &
        Driver_label_SetRunSequence   => label_SetRunSequence, &
        Driver_label_SetRunClock      => label_SetRunClock, &
        Driver_label_Finalize         => label_Finalize
      use NUOPC_Connector, only: conSS => SetServices
  ! - Handle build time ATM options:
      use RAQMSCHM,    only: RAQMSCHEM_SS  => SetServices
      use fv3gfs_cap_mod,        only: FV3_SS   => SetServices
  ! - Handle build time OCN options:
  ! - Handle build time ICE options:
  ! - Handle build time WAV options:
  ! - Handle build time LND options:
  ! - Handle build time IPM options:
  ! - Handle build time HYD options:
      use CHM,    only: GSDCHEM_SS  => SetServices
  ! - Mediator
      use module_MEDIATOR,        only: MED_SS     => SetServices
      use module_MEDSpaceWeather, only: MEDSW_SS   => SetServices

      USE module_EARTH_INTERNAL_STATE,ONLY: EARTH_INTERNAL_STATE        &
                                           ,WRAP_EARTH_INTERNAL_STATE
!
!      USE module_ATM_GRID_COMP
!
      USE module_NEMS_UTILS,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: EARTH_REGISTER
      PUBLIC :: VERBOSE_DIAGNOSTICS
!
!-----------------------------------------------------------------------
!

      LOGICAL, PRIVATE :: flag_verbose_diagnostics = .false.


      CONTAINS

      logical function verbose_diagnostics(set)
        !! Mutator for the verbose diagnostics flag; returns true if
        !! verbose diagnostics should be used, and false otherwise.
        !! If the "set" argument is present, then the flag is set to
        !! the given value.
        logical, optional :: set
        if(present(set)) then
           flag_verbose_diagnostics=set
        endif
        verbose_diagnostics=flag_verbose_diagnostics
      end function verbose_diagnostics

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE EARTH_REGISTER(EARTH_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: EARTH_GRID_COMP                               !<-- The EARTH component
!
      INTEGER,INTENT(OUT) :: RC_REG                                        !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
      type(ESMF_Config)             :: config
      
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_REG = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      ! Derive from NUOPC_Driver
      call NUOPC_CompDerive(EARTH_GRID_COMP, Driver_routine_SS, rc=RC)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=233, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return

      ! specializations:

      call NUOPC_CompSpecialize(EARTH_GRID_COMP, &
        specLabel=Driver_label_SetModelServices, specRoutine=SetModelServices, &
        rc=RC)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=240, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return
      
      call NUOPC_CompSpecialize(EARTH_GRID_COMP, &
        specLabel=Driver_label_SetRunSequence, specRoutine=SetRunSequence, &
        rc=RC)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=245, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return

      ! The NEMS Earth component is currently the top-level driver and
      ! does not need to coordinate Clocks with its parent.
      call ESMF_MethodRemove(EARTH_GRID_COMP, Driver_label_SetRunClock, rc=RC_REG)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=250, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return
      call NUOPC_CompSpecialize(EARTH_GRID_COMP, &
        specLabel=Driver_label_SetRunClock, specRoutine=NUOPC_NoOp, rc=RC_REG)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=253, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return
      
      call NUOPC_CompSpecialize(EARTH_GRID_COMP, &
        specLabel=Driver_label_Finalize, specRoutine=Finalize, &
        rc=RC)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=258, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return
      
      ! register an internal initialization method
      call NUOPC_CompSetInternalEntryPoint(EARTH_GRID_COMP, ESMF_METHOD_INITIALIZE, &
        phaseLabelList=(/"IPDv04p2"/), userRoutine=ModifyCplLists, rc=rc)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=263, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return

      ! create, open, and set the config
      config = ESMF_ConfigCreate(rc=RC)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=267, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return
      call ESMF_ConfigLoadFile(config, "nems.configure", rc=RC)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=269, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return
      call ESMF_GridCompSet(EARTH_GRID_COMP, config=config, rc=RC)
      if (ESMF_LogFoundError(RC, msg="Breaking out of subroutine", line=271, file="module_EARTH_GRID_COMP.F90", rcToReturn=RC_REG)) return
      
      ! Added the following Field Dictionary block to the EARTH component level
      ! in order to prevent different dictionary definitions in the lower
      ! components. Doing this here isn't without problems because it
      ! potentially makes the components (ATM & OCN) depend on this environment,
      ! which lowers their transferability to other coupled systems. However,
      ! extending the Field Dictionary is a temporary solution anyway (see the
      ! TODO: below), so this isn't going to stay for ever this way.
      
      ! Extend the NUOPC Field Dictionary to hold required entries.
      !TODO: In the long run this section will not be needed when we have
      !TODO: absorbed the needed standard names into the default dictionary.
      ! -> 20 fields identified as exports by the GSM component
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "air_density_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="air_density_height_lowest", &
          canonicalUnits="kg m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=292, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_zonal_moment_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_zonal_moment_flx", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=303, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_merid_moment_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_merid_moment_flx", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=314, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=325, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx_atm_into_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx_atm_into_ice", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=336, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx_atm_into_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx_atm_into_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=347, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=358, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_atm_into_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_atm_into_ice", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=369, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_atm_into_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_atm_into_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=380, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_down_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=391, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_flx")) then
!        write(6,*)'add dictionary meandownswflux'
!        call flush(6)
        !write(0,*)'add dictionary meandownswflux'
!        call flush(0)
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=406, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_fprec_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_fprec_rate", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=417, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_prec_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_prec_rate", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=428, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_evap_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_evap_rate", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=439, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_evap_rate_atm_into_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_evap_rate_atm_into_ice", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=450, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_evap_rate_atm_into_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_evap_rate_atm_into_ocn", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=461, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_zonal_moment_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_zonal_moment_flx", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=472, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_merid_moment_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_merid_moment_flx", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=483, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_sensi_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_sensi_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=494, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_laten_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_laten_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=505, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_down_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=516, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=527, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_temp_height2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_temp_height2m", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=538, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_spec_humid_height2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_spec_humid_height2m", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=549, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_u_wind_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_u_wind_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=560, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_v_wind_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_v_wind_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=571, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_zonal_wind_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_zonal_wind_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=582, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_merid_wind_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_merid_wind_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=593, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      !For MOM6 and WW3 variables to match:
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_partitioned_stokes_drift_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_partitioned_stokes_drift_1", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=605, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_partitioned_stokes_drift_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_partitioned_stokes_drift_1", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=617, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_partitioned_stokes_drift_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_partitioned_stokes_drift_2", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=628, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_partitioned_stokes_drift_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_partitioned_stokes_drift_2", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=639, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
     if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_partitioned_stokes_drift_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_partitioned_stokes_drift_3", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=650, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_partitioned_stokes_drift_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_partitioned_stokes_drift_3", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=662, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      ! end of MOM6 and WW3 variables to match
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_temp_height_surface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_temp_height_surface", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=674, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_pres_height_surface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pres_height_surface", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=685, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_surface_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_surface_height", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=696, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      ! -> Additional fields identified as needed by MOM5 and others...
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=708, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=719, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=730, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=741, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=752, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=763, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=774, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=785, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=796, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=807, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=818, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=829, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=840, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=851, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=862, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=873, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_salt_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_salt_rate", &
          canonicalUnits="kg psu m-2 s", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=884, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_runoff_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_runoff_rate", &
          canonicalUnits="kg m-2 s", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=895, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_calving_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_calving_rate", &
          canonicalUnits="kg m-2 s", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=906, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_runoff_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_runoff_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=917, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry(  &
        "mean_calving_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_calving_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=928, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ice_fraction")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ice_fraction", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=939, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_sw_pen_to_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sw_pen_to_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=950, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_up_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_up_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=961, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mass_of_overlying_sea_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mass_of_overlying_sea_ice", &
          canonicalUnits="kg", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=972, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "s_surf")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="s_surf", &
          canonicalUnits="psu", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=983, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "freezing_melting_potential")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="freezing_melting_potential", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=994, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      ! following two added for export from MOM6
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "accum_heat_frazil")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="accum_heat_frazil", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1006, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_melt_potential")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_melt_potential", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1017, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "u_surf")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="u_surf", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1028, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "v_surf")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="v_surf", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1039, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_lev")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_lev", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1050, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "wind_stress_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wind_stress_zonal", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1061, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "wind_stress_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wind_stress_merid", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1072, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocn_current_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocn_current_zonal", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1083, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocn_current_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocn_current_merid", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1094, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocn_current_idir")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocn_current_idir", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1105, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocn_current_jdir")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocn_current_jdir", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1116, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_surface_slope_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_surface_slope_zonal", &
          canonicalUnits="m m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1127, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_surface_slope_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_surface_slope_merid", &
          canonicalUnits="m m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1138, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_surface_slope_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_surface_slope_zonal", &
          canonicalUnits="m m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1149, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_surface_slope_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_surface_slope_merid", &
          canonicalUnits="m m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1160, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_air_ice_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_air_ice_zonal", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1171, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_air_ice_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_air_ice_merid", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1182, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_air_ocn_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_air_ocn_zonal", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1193, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_air_ocn_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_air_ocn_merid", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1204, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_ocn_ice_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_ocn_ice_zonal", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1215, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_ocn_ice_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_ocn_ice_merid", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1226, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_ocn_ice_idir")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_ocn_ice_idir", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1237, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_ocn_ice_jdir")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_ocn_ice_jdir", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1248, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mixed_layer_depth")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mixed_layer_depth", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1259, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1270, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1281, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_up_lw_flx_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_up_lw_flx_ice", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1292, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_up_lw_flx_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_up_lw_flx_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1303, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1314, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1325, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ir_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ir_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1336, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ir_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ir_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1347, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_vis_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_vis_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1358, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_vis_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_vis_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1369, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ocn_ir_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ocn_ir_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1380, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ocn_ir_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ocn_ir_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1391, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ocn_vis_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ocn_vis_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1402, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ocn_vis_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ocn_vis_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1413, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ice_ir_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ice_ir_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1424, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ice_ir_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ice_ir_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1435, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ice_vis_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ice_vis_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1446, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ice_vis_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ice_vis_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1457, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_land_sea_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_land_sea_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1468, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_temp_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_temp_height_lowest", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1479, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_spec_humid_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_spec_humid_height_lowest", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1490, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "humidity_2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="humidity_2m", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1501, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_zonal_wind_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_zonal_wind_height_lowest", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1512, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_merid_wind_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_merid_wind_height_lowest", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1523, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_pres_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pres_height_lowest", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1534, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_height_lowest", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1545, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocean_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocean_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1556, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ice_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ice_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1567, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "land_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="land_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1578, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      ! special HYCOM exports
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "surface_downward_eastward_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_downward_eastward_stress", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1590, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "surface_downward_northward_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_downward_northward_stress", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1601, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "wind_speed_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wind_speed_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1612, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "wind_speed_squared_10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wind_speed_squared_10m", &
          canonicalUnits="m2 s-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1623, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "friction_speed")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="friction_speed", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1634, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_lat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_lat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1645, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_sens_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sens_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1656, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "water_flux_into_sea_water")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="water_flux_into_sea_water", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1667, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "frozen_water_flux_into_sea_water")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="frozen_water_flux_into_sea_water", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1678, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "surface_temperature")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_temperature", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1689, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "air_surface_temperature")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="air_surface_temperature", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1700, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "temperature_2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_2m", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1711, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "upward_sea_ice_basal_available_heat_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="upward_sea_ice_basal_available_heat_flux", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1722, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      ! special HYCOM imports
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_area_fraction")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_area_fraction", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1734, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_x_stress_at_sea_ice_base")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_x_stress_at_sea_ice_base", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1745, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_y_stress_at_sea_ice_base")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_y_stress_at_sea_ice_base", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1756, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_sea_ice_basal_solar_heat_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_sea_ice_basal_solar_heat_flux", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1767, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "upward_sea_ice_basal_heat_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="upward_sea_ice_basal_heat_flux", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1778, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_sea_ice_basal_salt_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_sea_ice_basal_salt_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1789, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_sea_ice_basal_water_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_sea_ice_basal_water_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1800, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_surface_temperature")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_surface_temperature", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1811, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_temperature")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_temperature", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1822, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_thickness")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_thickness", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1833, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_x_velocity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_x_velocity", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1844, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_y_velocity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_y_velocity", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1855, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "net_heat_flx_to_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="net_heat_flx_to_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1867, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_fresh_water_to_ocean_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_fresh_water_to_ocean_rate", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1878, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_ice_volume")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_ice_volume", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1889, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_snow_volume")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_snow_volume", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1900, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
     
      !Mass flux of liquid runoff
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "Foxx_rofl")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="Foxx_rofl", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1913, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      !Mass flux of frozen runoff
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "Foxx_rofi")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="Foxx_rofi", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1925, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      !Ocean surface boundary layer depth
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "So_bldepth")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="So_bldepth", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1937, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
 
      ! Synonyms for HYCOM fields
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"surface_downward_eastward_stress",&
                          "mean_zonal_moment_flx           "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=1947, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"surface_downward_northward_stress",&
                          "mean_merid_moment_flx            "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=1954, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"mean_lat_flx       ",&
                          "mean_laten_heat_flx"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=1961, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"mean_sens_flx      ",&
                          "mean_sensi_heat_flx"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=1968, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out

      ! DCR - Fields added for Regional Application
      ! ATM-OCN-ICE-LND-HYD
      ! List of exisitng fields
      ! ice_mask, inst_down_lw_flx, inst_down_sw_flx, inst_height_lowest,
      ! inst_merid_wind_height_lowest, inst_pres_height_lowest,
      ! inst_pres_height_surface, inst_spec_humid_height_lowest,
      ! inst_temp_height_lowest, inst_temp_height_surface,
      ! inst_zonal_wind_height_lowest, mean_down_lw_flx, mean_down_sw_flx,
      ! mean_fprec_rate, mean_laten_heat_flx, mean_net_lw_flx, mean_net_sw_flx,
      ! mean_prec_rate, mean_sensi_heat_flx

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "aerodynamic_roughness_length")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="aerodynamic_roughness_length", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1990, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "canopy_moisture_storage")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="canopy_moisture_storage", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2001, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "carbon_dioxide")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="carbon_dioxide", &
          canonicalUnits="ppmv", & ! Units must be clarified
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2012, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "cosine_zenith_angle")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="cosine_zenith_angle", &
          canonicalUnits="degree", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2023, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "exchange_coefficient_heat")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="exchange_coefficient_heat", &
          canonicalUnits="W m-2 K-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2034, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "exchange_coefficient_heat_height2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="exchange_coefficient_heat_height2m", &
          canonicalUnits="W m-2 K-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2045, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "exchange_coefficient_moisture_height2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="exchange_coefficient_moisture_height2m", &
          canonicalUnits="kg m-2 s-1 Pa-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2056, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_wind_speed_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_wind_speed_height_lowest", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2067, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_cprec_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_cprec_rate", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2078, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_grnd_sensi_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_grnd_sensi_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2089, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_kinematic")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_kinematic", &
          canonicalUnits="Kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2100, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_surface_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_surface_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2111, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_surface_skin_temp")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_surface_skin_temp", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2122, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mixing_ratio_surface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mixing_ratio_surface", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2133, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "root_moisture")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="root_moisture", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2144, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "saturated_mixing_ratio")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="saturated_mixing_ratio", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2155, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_snow_area_fraction")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_snow_area_fraction", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2166, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_snow_thickness")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_snow_thickness", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2177, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_snow_melt_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_snow_melt_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2188, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_surface_snow")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_surface_snow", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2199, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_depth")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_depth", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2210, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_hydraulic_conductivity_at_saturation")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_hydraulic_conductivity_at_saturation", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2221, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2232, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer_1", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2243, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer_2", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2254, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer_3", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2265, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer_4")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer_4", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2276, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_porosity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_porosity", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2287, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2298, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer_1", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2309, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer_2", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2320, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer_3", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2331, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer_4")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer_4", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2342, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_temperature_bottom")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_temperature_bottom", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2353, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_type")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_type", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2364, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_moisture_content")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_moisture_content", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2375, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "subsurface_basin_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="subsurface_basin_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2386, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "subsurface_runoff_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="subsurface_runoff_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2397, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_microwave_emissivity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_microwave_emissivity", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2408, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_runoff_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_runoff_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2419, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "vegetation_type")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="vegetation_type", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2430, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_frozen_water_in_soil")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_frozen_water_in_soil", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2441, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2452, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer_1", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2463, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer_2", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2474, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer_3", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2485, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer_4")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer_4", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2496, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_total_water_in_soil")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_total_water_in_soil", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2507, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_total_water_in_soil_at_critical_point")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_total_water_in_soil_at_critical_point", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2518, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_total_water_in_soil_at_field_capacity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_total_water_in_soil_at_field_capacity", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2529, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_total_water_in_soil_at_wilting_point")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_total_water_in_soil_at_wilting_point", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2540, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "water_surface_height_above_reference_datum")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="water_surface_height_above_reference_datum", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2551, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx_atm_into_lnd")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx_atm_into_lnd", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2563, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_atm_into_lnd")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_atm_into_lnd", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2575, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      ! Fields from and to WW3

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wind_at_10m_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wind_at_10m_height", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2589, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"eastward_wind_at_10m_height",&
                          "inst_zonal_wind_height10m  "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=2597, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wind_at_10m_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wind_at_10m_height", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2608, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"northward_wind_at_10m_height",&
                          "inst_merid_wind_height10m   "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=2616, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out

      !For MOM6 and WW3 variables to match: 
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"surface_eastward_sea_water_velocity",&
                          "ocn_current_zonal                  "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=2625, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"surface_northward_sea_water_velocity",&
                          "ocn_current_merid                   "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=2632, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_stokes_drift_current")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_stokes_drift_current", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2643, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_stokes_drift_current")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_stokes_drift_current", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2655, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wave_bottom_current")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wave_bottom_current", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2667, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wave_bottom_current")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wave_bottom_current", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2679, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "wave_bottom_current_radian_frequency")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wave_bottom_current_radian_frequency", &
          canonicalUnits="rad s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2691, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wave_radiation_stress_gradient")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wave_radiation_stress_gradient", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2703, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wave_radiation_stress_gradient")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wave_radiation_stress_gradient", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2715, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2727, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_northward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_northward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2739, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2751, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "wave_induced_charnock_parameter")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wave_induced_charnock_parameter", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2763, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "wave_z0_roughness_length")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wave_z0_roughness_length", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2775, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "wave_bottom_current_period")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wave_bottom_current_period", &
          canonicalUnits="s", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2787, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif


      ! Fields from WAM to IPE

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wind_neutral")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wind_neutral", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2802, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wind_neutral")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wind_neutral", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2814, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "upward_wind_neutral")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="upward_wind_neutral", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2826, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temp_neutral")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temp_neutral", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2838, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "O_Density")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="O_Density", &
          canonicalUnits="m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2850, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "O2_Density")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="O2_Density", &
          canonicalUnits="m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2862, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "N2_Density")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="N2_Density", &
          canonicalUnits="m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2874, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="height", &
          canonicalUnits="km", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2886, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      ! Chemistry fields

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_pres_interface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pres_interface", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2900, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_pres_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pres_levels", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2912, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_geop_interface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_geop_interface", &
          canonicalUnits="m2 s-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2924, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_geop_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_geop_levels", &
          canonicalUnits="m2 s-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2936, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_temp_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_temp_levels", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2948, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_zonal_wind_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_zonal_wind_levels", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2960, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_merid_wind_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_merid_wind_levels", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2972, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_omega_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_omega_levels", &
          canonicalUnits="Pa s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2984, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_mass_frac")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_mass_frac", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=2996, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "raqms_inst_tracer_mass_frac")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="raqms_inst_tracer_mass_frac", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3008, &
          file="module_EARTH_GRID_COMP.F90")) &
          return ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_up_surface_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_up_surface_flx", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3020, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_down_surface_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_down_surface_flx", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3032, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_clmn_mass_dens")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_clmn_mass_dens", &
          canonicalUnits="g m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3044, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_anth_biom_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_anth_biom_flx", &
          canonicalUnits="ug m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3056, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_pbl_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pbl_height", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3068, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_cell_area")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_cell_area", &
          canonicalUnits="m2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3080, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_convective_rainfall_amount")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_convective_rainfall_amount", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3092, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_exchange_coefficient_heat_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_exchange_coefficient_heat_levels", &
          canonicalUnits="W m-2 K-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3104, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_spec_humid_conv_tendency_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_spec_humid_conv_tendency_levels", &
          canonicalUnits="kg kg-1 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3116, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_friction_velocity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_friction_velocity", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3128, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif


      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_rainfall_amount")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_rainfall_amount", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3141, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_soil_moisture_content")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_soil_moisture_content", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3153, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

! ajl add
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "cm_drag")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="cm_drag", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3166, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "ktop_bot_cnv")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ktop_kbot_cnv", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3177, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
!        "cldtau_frac_liq")) then
        "inst_cldtau")) then
        call NUOPC_FieldDictionaryAddEntry( &
!          standardName="cldtau_frac_liq", &
          standardName="inst_cldtau", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3190, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
!        "cldtau_frac_liq")) then
        "inst_vort")) then
        call NUOPC_FieldDictionaryAddEntry( &
!          standardName="cldtau_frac_liq", &
          standardName="inst_vort", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3203, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
! end ajl add

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_up_sensi_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_up_sensi_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3216, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_lwe_snow_thickness")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_lwe_snow_thickness", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3228, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_vegetation_area_frac")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_vegetation_area_frac", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3240, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_surface_roughness")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_surface_roughness", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3252, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

      ! Dummy fields

      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "dummyfield")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="dummyfield", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3266, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "dummyfield1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="dummyfield1", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3277, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "dummyfield2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="dummyfield2", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3288, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
      endif

!-----------------------------------------------------------------------
!
      END SUBROUTINE EARTH_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!

      subroutine SetModelServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        ! local variables
        integer                         :: localrc, stat, i, j, petCount
        character(ESMF_MAXSTR)          :: name
        type(WRAP_EARTH_INTERNAL_STATE) :: is
        type(ESMF_GridComp)             :: comp
        type(ESMF_Config)               :: config
        character(len=32), allocatable  :: compLabels(:)
        integer, allocatable            :: petList(:)
        character(len=10)               :: value
        character(len=20)               :: model, prefix
        character(len=160)              :: msg
        integer                         :: petListBounds(2)
        integer                         :: componentCount
        type(NUOPC_FreeFormat)          :: attrFF, fdFF

        rc = ESMF_SUCCESS

        ! query the Component for info
        call ESMF_GridCompGet(driver, name=name, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3326, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out

        ! allocate memory for the internal state and store in Component
        allocate(is%EARTH_INT_STATE, stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Allocation of internal state memory failed.", &
          line=3332, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90", rcToReturn=rc)) &
          return  ! bail out
        call ESMF_GridCompSetInternalState(driver, is, rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3336, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        
        ! get petCount and config
        call ESMF_GridCompGet(driver, petCount=petCount, config=config, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3341, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        
        ! read and ingest free format driver attributes
        attrFF = NUOPC_FreeFormatCreate(config, label="EARTH_attributes::", &
          relaxedflag=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3347, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        call NUOPC_CompAttributeIngest(driver, attrFF, addFlag=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3350, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3353, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        
        ! dump the current field dictionary into the Log file
        call ESMF_AttributeGet(driver, name="DumpFieldDictionary", &
          value=value, defaultValue="false", &
          convention="NUOPC", purpose="Instance", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3360, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        if (trim(value)=="true") then
          call ESMF_LogWrite( &
            "===>===>===>===> Begin Dumping Field Dictionary <===<===<===<===",&
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3366, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          call NUOPC_FieldDictionaryEgest(fdFF, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3369, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          call NUOPC_FreeFormatLog(fdFF, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3372, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          call ESMF_LogWrite( &
            "===>===>===>===> Done Dumping Field Dictionary <===<===<===<===", &
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3377, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        endif
        
        ! determine the generic component labels
        componentCount = ESMF_ConfigGetLen(config, &
          label="EARTH_component_list:", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3384, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        allocate(compLabels(componentCount), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Allocation of compLabels failed.", &
          line=3388, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90", rcToReturn=rc)) &
          return  ! bail out
        call ESMF_ConfigGetAttribute(config, valueList=compLabels, &
          label="EARTH_component_list:", count=componentCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=3393, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
        
        ! determine information for each component and add to the driver
        do i=1, componentCount
          ! construct component prefix
          prefix=trim(compLabels(i))
          ! read in petList bounds
          call ESMF_ConfigGetAttribute(config, petListBounds, &
            label=trim(prefix)//"_petlist_bounds:", default=-1, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3403, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          ! handle the default situation
          if (petListBounds(1)==-1 .or. petListBounds(2)==-1) then
            petListBounds(1) = 0
            petListBounds(2) = petCount - 1
          endif
          ! read in model instance name
          call ESMF_ConfigGetAttribute(config, model, &
            label=trim(prefix)//"_model:", default="none", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3413, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          ! check that there was a model instance specified
          if (trim(model) == "none") then
            ! Error condition: no model was specified
            write (msg, *) "No model was specified for component: ",trim(prefix)
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3418, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          endif
          ! set petList for this component
          allocate(petList(petListBounds(2)-petListBounds(1)+1))
          do j=petListBounds(1), petListBounds(2)
            petList(j-petListBounds(1)+1) = j ! PETs are 0 based
          enddo
          
          if (trim(model) == "satm") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3437, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "xatm") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3450, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "datawam") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3463, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "gsm") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3476, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "nmmb") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3489, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "fv3") then
            call NUOPC_DriverAddComp(driver, trim(prefix), FV3_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=3498, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          elseif (trim(model) == "datm") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3515, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "socn") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3528, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "xocn") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3541, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "hycom") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3554, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "mom5") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3567, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "mom6") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3580, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "pom") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3593, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "sice") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3606, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "xice") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3619, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "cice") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3632, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "swav") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3645, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "xwav") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3658, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "ww3") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3671, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "slnd") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3684, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "xlnd") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3697, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "noah") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3710, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "lis") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3723, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "sipm") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3736, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "xipm") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3749, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "ipe") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3762, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "dataipe") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3775, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "shyd") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3788, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "xhyd") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3801, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "wrfhydro") then
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3814, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          elseif (trim(model) == "gsdchem") then
            call NUOPC_DriverAddComp(driver, trim(prefix), GSDCHEM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=3823, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  !  bail out
         elseif (trim(model) == "raqmschem")then
            call NUOPC_DriverAddComp(driver, trim(prefix), RAQMSCHEM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=3836, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  !  bail out
          ! - Two mediator choices currently built into NEMS from internal
          elseif (trim(model) == "nems") then
            call NUOPC_DriverAddComp(driver, trim(prefix), MED_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=3851, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          elseif (trim(model) == "spaceweather") then
            call NUOPC_DriverAddComp(driver, trim(prefix), MEDSW_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=3856, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          else
            ! Error condition: unknown model requested
            write (msg, *) "The requested model '", trim(model), &
              "' is an invalid choice!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=3861, &
              file="module_EARTH_GRID_COMP.F90", rcToReturn=rc)
            return  ! bail out
          endif
          
          ! read and ingest free format component attributes
          attrFF = NUOPC_FreeFormatCreate(config, &
            label=trim(prefix)//"_attributes::", relaxedflag=.true., rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3870, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          call NUOPC_CompAttributeIngest(comp, attrFF, addFlag=.true., rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3873, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3876, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
          
          ! clean-up
          deallocate(petList)
          
        enddo
        

        ! clean-up
        deallocate(compLabels)
        
      end subroutine

  !-----------------------------------------------------------------------------
  
  subroutine SetRunSequence(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    character(ESMF_MAXSTR)          :: name
    type(ESMF_Config)               :: config
    type(NUOPC_FreeFormat)          :: runSeqFF

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=3914, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out

    ! read free format run sequence from config
    call ESMF_GridCompGet(driver, config=config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=3920, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
    runSeqFF = NUOPC_FreeFormatCreate(config, label="runSeq::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=3923, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out

    ! ingest FreeFormat run sequence
    call NUOPC_DriverIngestRunSequence(driver, runSeqFF, &
      autoAddConnectors=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=3929, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
    
    ! Diagnostic output
    if(verbose_diagnostics()) then
       call NUOPC_DriverPrint(driver, orderflag=.true., rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=3941, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
    endif
    
  end subroutine
    
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------

  subroutine Finalize(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    integer                         :: localrc, stat
    type(WRAP_EARTH_INTERNAL_STATE) :: is
    logical                         :: existflag
    character(ESMF_MAXSTR)          :: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=4240, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
    
    ! query Component for this internal State
    nullify(is%EARTH_INT_STATE)
    call ESMF_GridCompGetInternalState(driver, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=4246, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90")) return  ! bail out
      
    ! deallocate internal state memory
    deallocate(is%EARTH_INT_STATE, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of internal state memory failed.", &
      line=4252, file=trim(name)//":"//"module_EARTH_GRID_COMP.F90", rcToReturn=rc)) &
      return  ! bail out
      
  end subroutine
      
  !-----------------------------------------------------------------------------
  
  recursive subroutine ModifyCplLists(driver, importState, exportState, clock, &
    rc)
    type(ESMF_GridComp)  :: driver
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=160)              :: name, msg
    type(ESMF_CplComp), pointer     :: connectorList(:)
    integer                         :: i, j, cplListSize
    character(len=160), allocatable :: cplList(:)
    character(len=160)              :: value
    type(WRAP_EARTH_INTERNAL_STATE) :: is

    rc = ESMF_SUCCESS
    
    ! query Component for this internal State
    nullify(is%EARTH_INT_STATE)
    call ESMF_GridCompGetInternalState(driver, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=4279, &
      file="module_EARTH_GRID_COMP.F90")) &
      return  ! bail out

    call ESMF_LogWrite("Driver is in ModifyCplLists()", ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=4285, &
      file="module_EARTH_GRID_COMP.F90")) &
      return  ! bail out
    
    nullify(connectorList)
    call NUOPC_DriverGetComp(driver, compList=connectorList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=4292, &
      file="module_EARTH_GRID_COMP.F90")) &
      return  ! bail out
    
    write (msg,*) "Found ", size(connectorList), " Connectors."// &
      " Modifying CplList Attribute...."
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=4300, &
      file="module_EARTH_GRID_COMP.F90")) &
      return  ! bail out
      
    do i=1, size(connectorList)
      ! query Connector i for its name
      call ESMF_CplCompGet(connectorList(i), name=name, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=4308, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out
      ! access CplList for Connector i
      call NUOPC_CompAttributeGet(connectorList(i), name="CplList", &
        itemCount=cplListSize, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=4315, &
        file="module_EARTH_GRID_COMP.F90")) &
        return  ! bail out
      if (cplListSize>0) then
        allocate(cplList(cplListSize))
        call NUOPC_CompAttributeGet(connectorList(i), name="CplList", &
          valueList=cplList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=4323, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
        ! go through all of the entries in the cplList and add options
        do j=1, cplListSize
          cplList(j) = trim(cplList(j))//":DumpWeights=true"
          cplList(j) = trim(cplList(j))//":SrcTermProcessing=1:TermOrder=SrcSeq"
          ! add connection options read in from configuration file
          call ESMF_AttributeGet(connectorList(i), name="ConnectionOptions", &
            value=value, defaultValue="", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=4334, &
            file="module_EARTH_GRID_COMP.F90")) &
            return  ! bail out
          cplList(j) = trim(cplList(j))//trim(value)
        enddo
        ! store the modified cplList in CplList attribute of connector i
        call NUOPC_CompAttributeSet(connectorList(i), &
          name="CplList", valueList=cplList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=4343, &
          file="module_EARTH_GRID_COMP.F90")) &
          return  ! bail out
        deallocate(cplList)
      endif
    enddo
      
    deallocate(connectorList)
    
  end subroutine

  !-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END MODULE module_EARTH_GRID_COMP
!
!-----------------------------------------------------------------------
