# 1 "GFS_layer/GFS_driver.F90"
# 1 "/usr/include/stdc-predef.h" 1
# 1 "GFS_layer/GFS_driver.F90" 2
module GFS_driver

  use machine,                  only: kind_phys
  use GFS_typedefs,             only: GFS_init_type,                       &
                                      GFS_statein_type, GFS_stateout_type, &
                                      GFS_sfcprop_type, GFS_coupling_type, &
                                      GFS_control_type, GFS_grid_type,     &
                                      GFS_tbd_type,     GFS_cldprop_type,  &
                                      GFS_radtend_type, GFS_diag_type



  use module_radiation_driver,  only: GFS_radiation_driver, radupdate
  use module_physics_driver,    only: GFS_physics_driver
  use funcphys,                 only: gfuncphys
  use gfdl_cloud_microphys_mod, only: gfdl_cloud_microphys_init

  use physcons,                 only: gravit => con_g,    rair    => con_rd, &
                                      rh2o   => con_rv,                      &
                                      tmelt  => con_ttp,  cpair   => con_cp, &
                                      latvap => con_hvap, latice  => con_hfus

  implicit none

  private

!--------------------------------------------------------------------------------
! GFS_init_type
!--------------------------------------------------------------------------------
!   This container is the minimum set of data required from the dycore/atmosphere
!   component to allow proper initialization of the GFS physics
!  
!   Type is defined in GFS_typedefs.F90
!--------------------------------------------------------------------------------
! type GFS_init_type
!   public
!   integer :: me                                !< my MPI-rank
!   integer :: master                            !< master MPI-rank
!   integer :: isc                               !< starting i-index for this MPI-domain
!   integer :: jsc                               !< starting j-index for this MPI-domain
!   integer :: nx                                !< number of points in i-dir for this MPI rank
!   integer :: ny                                !< number of points in j-dir for this MPI rank
!   integer :: levs                              !< number of vertical levels
!   integer :: cnx                               !< number of points in i-dir for this cubed-sphere face
!                                                !< equal to gnx for lat-lon grids
!   integer :: cny                               !< number of points in j-dir for this cubed-sphere face
!                                                !< equal to gny for lat-lon grids
!   integer :: gnx                               !< number of global points in x-dir (i) along the equator
!   integer :: gny                               !< number of global points in y-dir (j) along any meridian
!   integer :: nlunit                            !< fortran unit number for file opens
!   integer :: logunit                           !< fortran unit number for writing logfile
!   integer :: dt_dycore                         !< dynamics time step in seconds
!   integer :: dt_phys                           !< physics  time step in seconds
!   integer :: bdat(8)                           !< model begin date in GFS format   (same as idat)
!   integer :: cdat(8)                           !< model current date in GFS format (same as jdat)
!   !--- blocking data
!   integer, pointer :: blksz(:)                 !< for explicit data blocking
!                                                !< default blksz(1)=[nx*ny]
!   !--- ak/bk for pressure level calculations
!   integer, pointer :: ak(:)                    !< from surface (k=1) to TOA (k=levs)
!   integer, pointer :: bk(:)                    !< from surface (k=1) to TOA (k=levs)
!   !--- grid metrics
!   real(kind=kind_phys), pointer :: xlon(:,:)   !< column longitude for MPI rank
!   real(kind=kind_phys), pointer :: xlat(:,:)   !< column latitude  for MPI rank
!   real(kind=kind_phys), pointer :: area(:,:)   !< column area for length scale calculations
!
!   character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
!                                                 !< based on name location in array
!   character(len=65) :: fn_nml                  !< namelist filename
!   character(len=*), pointer :: input_nml_file(:) !< character string containing full namelist
!                                                  !< for use with internal file reads
! end type GFS_init_type
!--------------------------------------------------------------------------------
    
!------------------
! Module parameters
!------------------

!----------------------------
! Module variable definitions
!----------------------------
  real(kind=kind_phys), parameter :: con_24  =   24.0_kind_phys
  real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
  real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
  real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys
  real(kind=kind_phys), parameter :: qmin    =    1.0e-10

  integer, allocatable :: blksz(:)

!----------------
! Public entities
!----------------
  public  GFS_initialize              !< GFS initialization routine

  public  GFS_time_vary_step          !< perform operations needed prior radiation or physics
  public  GFS_radiation_driver        !< radiation_driver (was grrad)
  public  GFS_physics_driver          !< physics_driver (was gbphys)
  public  GFS_stochastic_driver       !< stochastic physics


  CONTAINS
!*******************************************************************************************


!--------------
! GFS initialze
!--------------
  subroutine GFS_initialize (Model, Statein, Stateout, Sfcprop,     &
                             Coupling, Grid, Tbd, Cldprop, Radtend, & 




                             Diag, Init_parm)







!   use module_microphysics, only: gsmconst
    use cldwat2m_micro,      only: ini_micro
    use micro_mg2_0,         only: micro_mg_init2_0 => micro_mg_init
    use micro_mg3_0,         only: micro_mg_init3_0 => micro_mg_init
    use aer_cloud,           only: aer_cloud_init
    use module_ras,          only: ras_init
    use module_mp_thompson,  only: thompson_init
    use module_mp_wsm6,      only: wsm6init
    use cires_ugwp_module,   only: cires_ugwp_init


    !--- interface variables
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(inout) :: Statein(:)
    type(GFS_stateout_type),  intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),  intent(inout) :: Coupling(:)
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),   intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),   intent(inout) :: Radtend(:)
    type(GFS_diag_type),      intent(inout) :: Diag(:)





    type(GFS_init_type),      intent(in)    :: Init_parm


    !--- local variables
    integer :: nb
    integer :: nblks





    integer :: ntrac
    integer :: ix

    integer :: blocksize
    real(kind=kind_phys), allocatable :: si(:)
    real(kind=kind_phys), parameter   :: p_ref = 101325.0d0


    nblks = size(Init_parm%blksz)
    ntrac = size(Init_parm%tracer_names)
    allocate (blksz(nblks))
    blksz(:) = Init_parm%blksz(:)

# 180 "GFS_layer/GFS_driver.F90"

    !--- set control properties (including namelist read)
    call Model%init (Init_parm%nlunit, Init_parm%fn_nml,           &
                     Init_parm%me, Init_parm%master,               &
                     Init_parm%logunit, Init_parm%isc,             &
                     Init_parm%jsc, Init_parm%nx, Init_parm%ny,    &
                     Init_parm%levs, Init_parm%cnx, Init_parm%cny, &
                     Init_parm%gnx, Init_parm%gny,                 &
                     Init_parm%dt_dycore, Init_parm%dt_phys,       &
                     Init_parm%iau_offset, Init_parm%bdat,         &
                     Init_parm%cdat, Init_parm%tracer_names,       &
                     Init_parm%input_nml_file, Init_parm%tile_num  &





                     )

! For CCPP,  these are called automatically in GFS_phys_time_vary_init as part of CCPP physics init.
! The reason why these are in GFS_phys_time_vary_init and not in ozphys/h2ophys is that the ozone
! and h2o interpolation of the data read here is done in GFS_phys_time_vary_run, i.e. all work
! related to the ozone/h2o input data is in GFS_phys_time_vary, while ozphys/h2ophys are applying
! ozone/h2o forcing to the model state.

    call read_o3data  (Model%ntoz, Model%me, Model%master)
    call read_h2odata (Model%h2o_phys, Model%me, Model%master)
    if (Model%aero_in) then
      call read_aerdata (Model%me,Model%master,Model%iflip,Model%idate)
    endif
    if (Model%iccn) then
      call read_cidata  ( Model%me, Model%master)
    endif


! For CCPP,  stochastic_physics_init is called automatically as part of CCPP physics init

    !--- initializing stochastic physics
    call init_stochastic_physics(Model,Init_parm,nblks)
    if(Model%me == Model%master) print*,'do_skeb=',Model%do_skeb


    do nb = 1,nblks
      ix = Init_parm%blksz(nb)
!     write(0,*)' ix in gfs_driver=',ix,' nb=',nb
      call Statein  (nb)%create (ix, Model)
      call Stateout (nb)%create (ix, Model)
      call Sfcprop  (nb)%create (ix, Model)
      call Coupling (nb)%create (ix, Model)
      call Grid     (nb)%create (ix, Model)

      call Tbd      (nb)%create (ix, nb, Model)



      call Cldprop  (nb)%create (ix, Model)
      call Radtend  (nb)%create (ix, Model)
!--- internal representation of diagnostics
      call Diag     (nb)%create (ix, Model)
    enddo

# 272 "GFS_layer/GFS_driver.F90"

    !--- populate the grid components
    call GFS_grid_populate (Grid, Init_parm%xlon, Init_parm%xlat, Init_parm%area)

! For CCPP, stochastic_physics_sfc_init is called automatically as part of CCPP physics init

!   get land surface perturbations here (move to GFS_time_vary if wanting to
!   update each time-step
    call run_stochastic_physics_sfc(nblks,Model,Grid,Coupling)


! For CCPP, these are called automatically in GFS_phys_time_vary_init as part of CCPP physics init

    !--- read in and initialize ozone and water
    if (Model%ntoz > 0) then
      do nb = 1, nblks
          call setindxoz (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_o3, &
                        Grid(nb)%jindx2_o3, Grid(nb)%ddy_o3)
      enddo
    endif

    !--- read in and initialize IN and CCN
    if (Model%iccn) then
      do nb = 1, nblks
        call setindxci (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_ci, &
                        Grid(nb)%jindx2_ci, Grid(nb)%ddy_ci, Grid(nb)%xlon_d,     &
                        Grid(nb)%iindx1_ci,Grid(nb)%iindx2_ci,Grid(nb)%ddx_ci)
      enddo
    endif

    !--- read in and initialize aerosols
    if (Model%aero_in) then
      do nb = 1, nblks
        call setindxaer (Init_parm%blksz(nb),Grid(nb)%xlat_d,Grid(nb)%jindx1_aer, &
                        Grid(nb)%jindx2_aer, Grid(nb)%ddy_aer, Grid(nb)%xlon_d,   &
                        Grid(nb)%iindx1_aer,Grid(nb)%iindx2_aer,Grid(nb)%ddx_aer, &
                        Init_parm%me, Init_parm%master )
      enddo
    endif

    if (Model%h2o_phys) then
      do nb = 1, nblks
        call setindxh2o (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_h, &
                         Grid(nb)%jindx2_h, Grid(nb)%ddy_h)
      enddo
    endif


! For CCPP, this is called automatically in GFS_time_vary_pre_init as part of CCPP physics init

    !--- Call gfuncphys (funcphys.f) to compute all physics function tables.
    call gfuncphys ()


!   call gsmconst (Model%dtp, Model%me, .TRUE.) ! This is for Ferrier microphysics - notused - moorthi





    !--- define sigma level for radiation initialization 
    !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
    !--- ps is replaced with p0. The value of p0 uses that in http:
    !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
    allocate(si(Model%levr+1))
    si = (Init_parm%ak + Init_parm%bk * p_ref - Init_parm%ak(Model%levr+1)) &
             / (p_ref - Init_parm%ak(Model%levr+1))

    call rad_initialize (si,  Model%levr,         Model%ictm,    Model%isol,      &
           Model%ico2,        Model%iaer,         Model%ialb,    Model%iems,      &
           Model%ntcw,        Model%num_p2d,      Model%num_p3d, Model%npdf3d,    &
           Model%ntoz,        Model%iovr_sw,      Model%iovr_lw, Model%isubc_sw,  &
           Model%isubc_lw,    Model%icliq_sw,     Model%crick_proof, Model%ccnorm,&
           Model%imp_physics, Model%norad_precip, Model%idate,   Model%iflip,  Model%me)
    deallocate (si)


!   microphysics initialization calls
!   --------------------------------- 

    if (Model%imp_physics == Model%imp_physics_mg) then          !--- initialize Morrison-Gettelman microphysics

      if (Model%fprcp <= 0) then
        call ini_micro (Model%mg_dcs, Model%mg_qcvar, Model%mg_ts_auto_ice(1))
      elseif (Model%fprcp == 1) then
        call micro_mg_init2_0(kind_phys, gravit, rair, rh2o, cpair,              &
                              tmelt, latvap, latice, 1.01_kind_phys,             &
                              Model%mg_dcs, Model%mg_ts_auto_ice,                &
                              Model%mg_qcvar,                                    &
                              Model%microp_uniform, Model%do_cldice,             &
                              Model%hetfrz_classnuc,                             &
                              Model%mg_precip_frac_method,                       &
                              Model%mg_berg_eff_factor,                          &
                              Model%sed_supersat,    Model%do_sb_physics,        &
                              Model%mg_do_ice_gmao,  Model%mg_do_liq_liu,        &
                              Model%mg_nccons,       Model%mg_nicons,            &
                              Model%mg_ncnst,        Model%mg_ninst)
      elseif (Model%fprcp == 2) then
        call micro_mg_init3_0(kind_phys, gravit, rair, rh2o, cpair,              &
                              tmelt, latvap, latice, 1.01_kind_phys,             &
                              Model%mg_dcs, Model%mg_ts_auto_ice,                &
                              Model%mg_qcvar,                                    &
                              Model%mg_do_hail,       Model%mg_do_graupel,       &
                              Model%microp_uniform,   Model%do_cldice,           &
                              Model%hetfrz_classnuc,                             &
                              Model%mg_precip_frac_method,                       &
                              Model%mg_berg_eff_factor,                          &
                              Model%sed_supersat,    Model%do_sb_physics,        &
                              Model%mg_do_ice_gmao,  Model%mg_do_liq_liu,        &
                              Model%mg_nccons,       Model%mg_nicons,            &
                              Model%mg_ncnst,        Model%mg_ninst,             &
                              Model%mg_ngcons,       Model%mg_ngnst)
      else
        write(0,*)' Model%fprcp = ',Model%fprcp,' is not a valid option - aborting'
        stop
      
      endif
      call aer_cloud_init ()

!
    elseif (Model%imp_physics == Model%imp_physics_thompson) then       !--- initialize Thompson Cloud microphysics
      if(Model%do_shoc) then
        print *,'SHOC is not currently compatible with Thompson MP -- shutting down'
        stop 
      endif 
! For CCPP the Thompson MP init is called automatically as part of CCPP physics init

      call thompson_init()                     !--- add aerosol version later
      if(Model%ltaerosol) then 
        print *,'Aerosol awareness is not included in this version of Thompson MP -- shutting down'
        stop 
      endif 
!
    elseif(Model%imp_physics == Model%imp_physics_wsm6) then        !--- initialize WSM6 Cloud microphysics
      if(Model%do_shoc) then 
        print *,'SHOC is not currently compatible with WSM6 -- shutting down'
        stop 
      endif 
      call  wsm6init()

!
    else if(Model%imp_physics == Model%imp_physics_gfdl) then      !--- initialize GFDL Cloud microphysics
! For CCPP the GFDL MP init is called automatically as part of CCPP physics init

      if(Model%do_shoc) then 
         print *,'SHOC is not currently compatible with GFDL MP -- shutting down'
         stop 
      endif 
       call gfdl_cloud_microphys_init (Model%me, Model%master, Model%nlunit, Model%input_nml_file, &
                                       Init_parm%logunit, Model%fn_nml)

    endif 


    !--- initialize ras
    if (Model%ras) call ras_init (Model%levs, Model%me)


! DH* Even though this gets called through CCPP in lsm_noah_init, we also
! need to do this here as long as FV3GFS_io.F90 is calculating Sfcprop%sncovr
! when reading restart files (which it shouldn't, this should be moved to physics).
    !--- initialize soil vegetation
    call set_soilveg(Model%me, Model%isot, Model%ivegsrc, Model%nlunit)
! *DH

    !--- lsidea initialization
    if (Model%lsidea) then
      print *,' LSIDEA is active but needs to be reworked for FV3 - shutting down'
      stop
      !--- NEED TO get the logic from the old phys/gloopb.f initialization area
    endif


!----  initialization of cires_ugwp .
!     if ( Model%me == Model%master) print *,  ' VAY-nml ',  Model%fn_nml
!     if ( Model%me == Model%master) print *,  ' VAY-nml2 ', Model%input_nml_file
    if (Model%do_ugwp) then
!     if ( Model%me == Model%master) print *,  ' VAY-nml ',  Model%fn_nml,
!     Model%input_nml_file
      call cires_ugwp_init(Model%me,      Model%master, Model%nlunit,  Init_parm%logunit, &
                           Model%fn_nml,  Model%lonr,   Model%latr,    Model%levs,        &
                           Init_parm%ak,  Init_parm%bk, p_ref,         Model%dtp,         &
                           Model%cdmbgwd, Model%cgwf,   Model%prslrd0, Model%ral_ts)
    endif


! For CCPP, this is called automatically as part of CCPP_physics_init

    !--- Initialize cellular automata
    if(Model%do_ca)then
    blocksize=size(Grid(1)%xlon)
    call cellular_automata(Model%kdt, Statein, Coupling, Diag, nblks, Model%levs,     &
                           Model%nca, Model%ncells, Model%nlives, Model%nfracseed,    &
                           Model%nseed, Model%nthresh, Model%ca_global, Model%ca_sgs, &
                           Model%iseed_ca, Model%ca_smooth, Model%nspinup, blocksize)
    endif


    !--- sncovr may not exist in ICs from chgres.
    !--- FV3GFS handles this as part of the IC ingest
    !--- this note is placed here to alert users to study
    !--- the FV3GFS_io.F90 module

  end subroutine GFS_initialize



!-------------------------------------------------------------------------
! time_vary_step
!-------------------------------------------------------------------------
!    routine called prior to radiation and physics steps to handle:
!      1) sets up various time/date variables
!      2) sets up various triggers
!      3) defines random seed indices for radiation (in a reproducible way)
!      5) interpolates coefficients for prognostic ozone calculation
!      6) performs surface data cycling via the GFS gcycle routine
!-------------------------------------------------------------------------
  subroutine GFS_time_vary_step (Model, Statein, Stateout, Sfcprop, Coupling, & 
                                 Grid, Tbd, Cldprop, Radtend, Diag)

    implicit none

    !--- interface variables
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(inout) :: Statein(:)
    type(GFS_stateout_type),  intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),  intent(inout) :: Coupling(:)
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),   intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),   intent(inout) :: Radtend(:)
    type(GFS_diag_type),      intent(inout) :: Diag(:)


    !--- local variables
    integer :: nb, nblks, k, kdt_rad, kdt_iau, blocksize
    logical :: iauwindow_center
    real(kind=kind_phys) :: rinc(5)
    real(kind=kind_phys) :: sec, sec_zero
    real(kind=kind_phys), parameter :: cn_hr     = 3600._kind_phys

    nblks = size(blksz)
    !--- Model%jdat is being updated directly inside of FV3GFS_cap.F90
    !--- update calendars and triggers
    rinc(1:5)   = 0
    call w3difdat(Model%jdat,Model%idat,4,rinc)
    sec = rinc(4)
    Model%phour = sec/con_hr
    !--- set current bucket hour
    Model%zhour = Model%phour
    Model%fhour = (sec + Model%dtp)/con_hr
    Model%kdt   = nint((sec + Model%dtp)/Model%dtp)

    Model%ipt    = 1
    Model%lprnt  = .false.
    Model%lssav  = .true.

    !--- radiation triggers
    Model%lsswr  = (mod(Model%kdt, Model%nsswr) == 1)
    Model%lslwr  = (mod(Model%kdt, Model%nslwr) == 1)
    !--- allow for radiation to be called on every physics time step, if needed
    if (Model%nsswr == 1)  Model%lsswr = .true.
    if (Model%nslwr == 1)  Model%lslwr = .true.

    !--- set the solar hour based on a combination of phour and time initial hour
    Model%solhr  = mod(Model%phour+Model%idate(1),con_24)

    if ((Model%debug) .and. (Model%me == Model%master)) then
      print *,'   sec ', sec
      print *,'   kdt ', Model%kdt
      print *,' nsswr ', Model%nsswr
      print *,' nslwr ', Model%nslwr
      print *,' nscyc ', Model%nscyc
      print *,' lsswr ', Model%lsswr
      print *,' lslwr ', Model%lslwr
      print *,' fhour ', Model%fhour
      print *,' phour ', Model%phour
      print *,' solhr ', Model%solhr
    endif

    !--- radiation time varying routine
    if (Model%lsswr .or. Model%lslwr) then
      call GFS_rad_time_vary (Model, Statein, Tbd, sec)
    endif

    !--- physics time varying routine
    call GFS_phys_time_vary (Model, Grid, Tbd, Statein)

    !--- repopulate specific time-varying sfc properties for AMIP/forecast runs
    if (Model%nscyc >  0) then
      if (mod(Model%kdt,Model%nscyc) == 1) THEN
        call gcycle (nblks, Model, Grid(:), Sfcprop(:), Cldprop(:))
      endif
    endif

    !--- determine if diagnostics buckets need to be cleared
    sec_zero = nint(Model%fhzero*con_hr)
    if (sec_zero >= nint(max(Model%fhswr,Model%fhlwr))) then
      if (mod(Model%kdt,Model%nszero) == 1) then
        do nb = 1,nblks
          call Diag(nb)%rad_zero  (Model)
          call Diag(nb)%phys_zero (Model)
    !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
        enddo
      endif
    else
      if (mod(Model%kdt,Model%nszero) == 1) then
        do nb = 1,nblks
          call Diag(nb)%phys_zero (Model)
    !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
        enddo
      endif
      kdt_rad = nint(min(Model%fhswr,Model%fhlwr)/Model%dtp)
      if (mod(Model%kdt, kdt_rad) == 1) then
        do nb = 1,nblks
          call Diag(nb)%rad_zero  (Model)
    !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
        enddo
      endif
    endif
!
    if (Model%iau_offset > 0) then
      kdt_iau = nint(Model%iau_offset*con_hr/Model%dtp)
      if (Model%kdt == kdt_iau+1) then
        iauwindow_center = .true.
        do nb = 1,nblks
          call Diag(nb)%rad_zero  (Model)
          call Diag(nb)%phys_zero (Model,iauwindow_center=iauwindow_center)
        enddo
        if(Model%me == Model%master) print *,'in gfs_driver, at iau_center, zero out rad/phys accumulated diag fields, kdt=',Model%kdt,'kdt_iau=',kdt_iau,'iau_offset=',Model%iau_offset
      endif
    endif
    call run_stochastic_physics(nblks,Model,Grid(:),Coupling(:))

    if(Model%do_ca)then
      blocksize = size(Grid(1)%xlon)
      call cellular_automata(Model%kdt,Statein,Coupling,Diag,nblks,Model%levs,      &
              Model%nca,Model%ncells,Model%nlives,Model%nfracseed,                  &
              Model%nseed,Model%nthresh,Model%ca_global,Model%ca_sgs,Model%iseed_ca,&
              Model%ca_smooth,Model%nspinup,blocksize)
    endif


! kludge for output
    if (Model%do_skeb) then
      do nb = 1,nblks
        do k=1,Model%levs
          Diag(nb)%skebu_wts(:,k) = Coupling(nb)%skebu_wts(:,Model%levs-k+1)
          Diag(nb)%skebv_wts(:,k) = Coupling(nb)%skebv_wts(:,Model%levs-k+1)
        enddo
      enddo
    endif
    !if (Model%do_sppt) then
    !  do nb = 1,nblks
    !    do k=1,Model%levs
    !      Diag(nb)%sppt_wts(:,k) = Coupling(nb)%sppt_wts(:,Model%levs-k+1)
    !    enddo
    !  enddo
    !endif
    if (Model%do_shum) then
      do nb = 1,nblks
        do k=1,Model%levs
          Diag(nb)%shum_wts(:,k)=Coupling(nb)%shum_wts(:,Model%levs-k+1)
        enddo
      enddo
    endif

  end subroutine GFS_time_vary_step


!-------------------------------------------------------------------------
! GFS stochastic_driver
!-------------------------------------------------------------------------
!    routine called prior to radiation and physics steps to handle:
!      1) sets up various time/date variables
!      2) sets up various triggers
!      3) defines random seed indices for radiation (in a reproducible way)
!      5) interpolates coefficients for prognostic ozone calculation
!      6) performs surface data cycling via the GFS gcycle routine
!-------------------------------------------------------------------------
  subroutine GFS_stochastic_driver (Model, Statein, Stateout, Sfcprop, Coupling, &
                                    Grid, Tbd, Cldprop, Radtend, Diag)

    implicit none

    !--- interface variables
! DH* gfortran correctly throws an error if the intent() declarations
! for arguments differ between the actual routine (here) and the dummy
! interface routine (IPD_func0d_proc in IPD_typedefs.F90):
!
! Error: Interface mismatch in procedure pointer assignment at (1): INTENT mismatch in argument 'control'
!
! Since IPD_func0d_proc declares all arguments as intent(inout), we
! need to do the same here - however, this way we are loosing the
! valuable information on the actual intent to this routine. *DH
# 680 "GFS_layer/GFS_driver.F90"
    type(GFS_control_type),   intent(in   ) :: Model
    type(GFS_statein_type),   intent(in   ) :: Statein
    type(GFS_stateout_type),  intent(in   ) :: Stateout
    type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
    type(GFS_coupling_type),  intent(inout) :: Coupling
    type(GFS_grid_type),      intent(in   ) :: Grid
    type(GFS_tbd_type),       intent(in   ) :: Tbd
    type(GFS_cldprop_type),   intent(in   ) :: Cldprop
    type(GFS_radtend_type),   intent(in   ) :: Radtend
    type(GFS_diag_type),      intent(inout) :: Diag

    !--- local variables
    integer :: k, i
    real(kind=kind_phys) :: upert, vpert, tpert, qpert, qnew,sppt_vwt
    !real(kind=kind_phys),dimension(size(Statein%tgrs,1),size(Statein%tgrs,2)) :: tconvtend, &
    !                     qconvtend,uconvtend,vconvtend

     if (Model%do_sppt) then
       do k = 1,size(Statein%tgrs,2)
         do i = 1,size(Statein%tgrs,1)
           sppt_vwt=1.0
           if (Diag%zmtnblck(i).EQ.0.0) then
              sppt_vwt=1.0
           else 
              if (k.GT.Diag%zmtnblck(i)+2) then
                 sppt_vwt=1.0
              endif
              if (k.LE.Diag%zmtnblck(i)) then
                 sppt_vwt=0.0
              endif
              if (k.EQ.Diag%zmtnblck(i)+1) then
                 sppt_vwt=0.333333
              endif
              if (k.EQ.Diag%zmtnblck(i)+2) then
                 sppt_vwt=0.666667
              endif
           endif
           if (Model%use_zmtnblck)then
              Coupling%sppt_wts(i,k)=(Coupling%sppt_wts(i,k)-1)*sppt_vwt+1.0
           endif
           Diag%sppt_wts(i,Model%levs-k+1)=Coupling%sppt_wts(i,k)

           
          ! if(Model%isppt_deep)then

           ! tconvtend(i,k)=Coupling%tconvtend(i,k)
           ! qconvtend(i,k)=Coupling%qconvtend(i,k)
           ! uconvtend(i,k)=Coupling%uconvtend(i,k)
           ! vconvtend(i,k)=Coupling%vconvtend(i,k)           


           ! upert = (Stateout%gu0(i,k)   - Statein%ugrs(i,k) - uconvtend(i,k)) + uconvtend(i,k) * Coupling%sppt_wts(i,k)
           ! vpert = (Stateout%gv0(i,k)   - Statein%vgrs(i,k) - vconvtend(i,k)) + vconvtend(i,k) * Coupling%sppt_wts(i,k)
           ! tpert = (Stateout%gt0(i,k)   - Statein%tgrs(i,k) - Tbd%dtdtr(i,k) - tconvtend(i,k)) + tconvtend(i,k) * Coupling%sppt_wts(i,k)
           ! qpert = (Stateout%gq0(i,k,1) - Statein%qgrs(i,k,1) - qconvtend(i,k)) + qconvtend(i,k) * Coupling%sppt_wts(i,k)

           !else
           
            upert = (Stateout%gu0(i,k)   - Statein%ugrs(i,k))   * Coupling%sppt_wts(i,k)
            vpert = (Stateout%gv0(i,k)   - Statein%vgrs(i,k))   * Coupling%sppt_wts(i,k)
            tpert = (Stateout%gt0(i,k)   - Statein%tgrs(i,k) - Tbd%dtdtr(i,k)) * Coupling%sppt_wts(i,k)
            qpert = (Stateout%gq0(i,k,1) - Statein%qgrs(i,k,1)) * Coupling%sppt_wts(i,k)
 
           !endif

           Stateout%gu0(i,k)  = Statein%ugrs(i,k)+upert
           Stateout%gv0(i,k)  = Statein%vgrs(i,k)+vpert
 
           !negative humidity check
           qnew = Statein%qgrs(i,k,1)+qpert
           if (qnew >= 1.0e-10) then
              Stateout%gq0(i,k,1) = qnew
              Stateout%gt0(i,k)   = Statein%tgrs(i,k) + tpert + Tbd%dtdtr(i,k)
           endif
         enddo
       enddo

       !if(Model%isppt_deep == .true.)then
       !  Sfcprop%tprcp(:) = Sfcprop%tprcp(:) + (Coupling%sppt_wts(:,15) - 1 )*Diag%rainc(:)
       !  Diag%totprcp(:)  = Diag%totprcp(:)  + (Coupling%sppt_wts(:,15) - 1 )*Diag%rainc(:) 
       !  Diag%cnvprcp(:)  = Diag%cnvprcp(:)  + (Coupling%sppt_wts(:,15) - 1 )*Diag%rainc(:)
       !!  ! bucket precipitation adjustment due to sppt                                                                                                                                                                     
       !  Diag%totprcpb(:) = Diag%totprcpb(:)  + (Coupling%sppt_wts(:,15) - 1 )*Diag%rainc(:)
       !  Diag%cnvprcpb(:) = Diag%cnvprcpb(:)  + (Coupling%sppt_wts(:,15) - 1 )*Diag%rainc(:)


       !  if (Model%cplflx) then !Need to make proper adjustments for deep convection only perturbations
       !    Coupling%rain_cpl(:) = Coupling%rain_cpl(:) + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%drain_cpl(:)
       !    Coupling%snow_cpl(:) = Coupling%snow_cpl(:) + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%dsnow_cpl(:)
       !  endif

       !else

        ! instantaneous precip rate going into land model at the next time step
        Sfcprop%tprcp(:) = Coupling%sppt_wts(:,15)*Sfcprop%tprcp(:)
        Diag%totprcp(:)      = Diag%totprcp(:)      + (Coupling%sppt_wts(:,15) - 1 )*Diag%rain(:)
        ! acccumulated total and convective preciptiation
        Diag%cnvprcp(:)      = Diag%cnvprcp(:)      + (Coupling%sppt_wts(:,15) - 1 )*Diag%rainc(:)
        ! bucket precipitation adjustment due to sppt
        Diag%totprcpb(:)      = Diag%totprcpb(:)      + (Coupling%sppt_wts(:,15) - 1 )*Diag%rain(:)
        Diag%cnvprcpb(:)      = Diag%cnvprcpb(:)      + (Coupling%sppt_wts(:,15) - 1 )*Diag%rainc(:)
      

        if (Model%cplflx) then
           Coupling%rain_cpl(:) = Coupling%rain_cpl(:) + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%drain_cpl(:)
           Coupling%snow_cpl(:) = Coupling%snow_cpl(:) + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%dsnow_cpl(:)
        endif

       !endif 

     endif

     if (Model%do_shum) then
       Stateout%gq0(:,:,1) = Stateout%gq0(:,:,1)*(1.0 + Coupling%shum_wts(:,:))
     endif

     if (Model%do_skeb) then
       do k = 1,size(Statein%tgrs,2)
           Stateout%gu0(:,k) = Stateout%gu0(:,k)+Coupling%skebu_wts(:,k)*(Statein%diss_est(:,k))
           Stateout%gv0(:,k) = Stateout%gv0(:,k)+Coupling%skebv_wts(:,k)*(Statein%diss_est(:,k))
       !    print*,'in do skeb',Coupling%skebu_wts(1,k),Statein%diss_est(1,k)
       enddo
     endif

  end subroutine GFS_stochastic_driver



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------
! GFS_rad_time_vary
!-----------------------------------------------------------------------
!
!  Routine containing all of the setup logic originally in phys/gloopr.f
!
!-----------------------------------------------------------------------
  subroutine GFS_rad_time_vary (Model, Statein, Tbd, sec)

    use physparam,        only: ipsd0, ipsdlim, iaerflg
    use mersenne_twister, only: random_setseed, random_index, random_stat

    implicit none

    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(in)    :: Statein(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    real(kind=kind_phys),     intent(in)    :: sec
    !--- local variables
    type (random_stat) :: stat
    integer :: ix, nb, j, i, nblks, ipseed
    integer :: numrdm(Model%cnx*Model%cny*2)

    nblks = size(blksz,1)

    call radupdate (Model%idat, Model%jdat, Model%fhswr, Model%dtf,  Model%lsswr, &
                    Model%me,   Model%slag, Model%sdec,  Model%cdec, Model%solcon)

    !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
    if ((Model%isubc_lw==2) .or. (Model%isubc_sw==2)) then
      ipseed = mod(nint(con_100*sqrt(sec)), ipsdlim) + 1 + ipsd0
      call random_setseed (ipseed, stat)
      call random_index (ipsdlim, numrdm, stat)

      !--- set the random seeds for each column in a reproducible way
      ix = 0
      nb = 1
      do j = 1,Model%ny
        do i = 1,Model%nx
          ix = ix + 1
          if (ix > blksz(nb)) then
            ix = 1
            nb = nb + 1
          endif
          !--- for testing purposes, replace numrdm with '100'
          Tbd(nb)%icsdsw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx)
          Tbd(nb)%icsdlw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + Model%cnx*Model%cny)
        enddo
      enddo
    endif  ! isubc_lw and isubc_sw

    if (Model%imp_physics == 99) then
      if (Model%kdt == 1) then
        do nb = 1,nblks
          Tbd(nb)%phy_f3d(:,:,1) = Statein(nb)%tgrs
          Tbd(nb)%phy_f3d(:,:,2) = max(qmin,Statein(nb)%qgrs(:,:,1))
          Tbd(nb)%phy_f3d(:,:,3) = Statein(nb)%tgrs
          Tbd(nb)%phy_f3d(:,:,4) = max(qmin,Statein(nb)%qgrs(:,:,1))
          Tbd(nb)%phy_f2d(:,1)   = Statein(nb)%prsi(:,1)
          Tbd(nb)%phy_f2d(:,2)   = Statein(nb)%prsi(:,1)
        enddo
      endif
    endif

  end subroutine GFS_rad_time_vary


!-----------------------------------------------------------------------
! GFS_phys_time_vary
!-----------------------------------------------------------------------
!
!  Routine containing all of the setup logic originally in phys/gloopb.f
!
!-----------------------------------------------------------------------
  subroutine GFS_phys_time_vary (Model, Grid, Tbd, Statein)
    use mersenne_twister, only: random_setseed, random_number

    implicit none
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    type(GFS_statein_type),   intent(in)    :: Statein(:)
    !--- local variables
    integer :: nb, ix, k, j, i, nblks, iseed, iskip
    real(kind=kind_phys) :: wrk(1)
    real(kind=kind_phys) :: rannie(Model%cny)
    real(kind=kind_phys) :: rndval(Model%cnx*Model%cny*Model%nrcm)

    nblks = size(blksz,1)

    !--- switch for saving convective clouds - cnvc90.f 
    !--- aka Ken Campana/Yu-Tai Hou legacy
    if ((mod(Model%kdt,Model%nsswr) == 0) .and. (Model%lsswr)) then
      !--- initialize,accumulate,convert
      Model%clstp = 1100 + min(Model%fhswr/con_hr,Model%fhour,con_99)
    elseif (mod(Model%kdt,Model%nsswr) == 0) then
      !--- accumulate,convert
      Model%clstp = 0100 + min(Model%fhswr/con_hr,Model%fhour,con_99)
    elseif (Model%lsswr) then
      !--- initialize,accumulate
      Model%clstp = 1100
    else
      !--- accumulate
      Model%clstp = 0100
    endif

    !--- random number needed for RAS and old SAS and when cal_pre=.true.
    !    Model%imfdeepcnv < 0 when Model%ras = .true.
    if ( (Model%imfdeepcnv <= 0 .or. Model%cal_pre) .and. Model%random_clds ) then
      iseed = mod(con_100*sqrt(Model%fhour*con_hr),1.0d9) + Model%seed0
      call random_setseed(iseed)
      call random_number(wrk)
      do i = 1,Model%cnx*Model%nrcm
        iseed = iseed + nint(wrk(1)*1000.0) * i
        call random_setseed(iseed)
        call random_number(rannie)
        rndval(1+(i-1)*Model%cny:i*Model%cny) = rannie(1:Model%cny)
      enddo

      do k = 1,Model%nrcm
        iskip = (k-1)*Model%cnx*Model%cny
        ix = 0
        nb = 1
        do j = 1,Model%ny
          do i = 1,Model%nx
            ix = ix + 1
            if (ix > blksz(nb)) then
              ix = 1
              nb = nb + 1
            endif
            Tbd(nb)%rann(ix,k) = rndval(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + iskip)
          enddo
        enddo
      enddo
    endif  ! imfdeepcnv, cal_re, random_clds

    !--- o3 interpolation
    if (Model%ntoz > 0) then
      do nb = 1, nblks
        call ozinterpol (Model%me, blksz(nb), Model%idate, Model%fhour,     &
                         Grid(nb)%jindx1_o3,  Grid(nb)%jindx2_o3,           &
                         Tbd(nb)%ozpl,        Grid(nb)%ddy_o3)
      enddo
    endif

    !--- h2o interpolation
     if (Model%h2o_phys) then
       do nb = 1, nblks
         call h2ointerpol (Model%me, blksz(nb), Model%idate, Model%fhour,   &
                           Grid(nb)%jindx1_h,   Grid(nb)%jindx2_h,          &
                           Tbd(nb)%h2opl,       Grid(nb)%ddy_h)
       enddo
     endif

    !--- ICCN interpolation
    if (Model%ICCN ) then
      do nb = 1, nblks
        call ciinterpol (Model%me, blksz(nb), Model%idate, Model%fhour, &
                         Grid(nb)%jindx1_ci, Grid(nb)%jindx2_ci,        &
                         Grid(nb)%ddy_ci,Grid(nb)%iindx1_ci,            &
                         Grid(nb)%iindx2_ci,Grid(nb)%ddx_ci,            &
                         Model%levs,Statein(nb)%prsl,                   &
                         Tbd(nb)%in_nm, Tbd(nb)%ccn_nm)
      enddo
    endif

    !--- aerosol interpolation
     if (Model%aero_in ) then
      do nb = 1, nblks
        call aerinterpol (Model%me, Model%master, blksz(nb),            &
                         Model%idate, Model%fhour,                      &
                         Grid(nb)%jindx1_aer, Grid(nb)%jindx2_aer,      &
                         Grid(nb)%ddy_aer,Grid(nb)%iindx1_aer,          &
                         Grid(nb)%iindx2_aer,Grid(nb)%ddx_aer,          &
                         Model%levs,Statein(nb)%prsl,                   &
                         Tbd(nb)%aer_nm)
      enddo
     endif

  end subroutine GFS_phys_time_vary



!------------------
! GFS_grid_populate
!------------------
  subroutine GFS_grid_populate (Grid, xlon, xlat, area)
    use physcons,                 only: pi => con_pi

    implicit none

    type(GFS_grid_type)              :: Grid(:)
    real(kind=kind_phys), intent(in) :: xlon(:,:)
    real(kind=kind_phys), intent(in) :: xlat(:,:)
    real(kind=kind_phys), intent(in) :: area(:,:)
    real(kind=kind_phys), parameter  :: rad2deg = 180.0_kind_phys/pi

    !--- local variables
    integer :: nb, ix, blksz, i, j

    blksz = size(Grid(1)%xlon)

    nb = 1
    ix = 0
    do j = 1,size(xlon,2)
      do i = 1,size(xlon,1)
        ix=ix+1
        if (ix > blksz) then
          nb = nb + 1
          ix = 1
        endif
        Grid(nb)%xlon(ix)   = xlon(i,j)
        Grid(nb)%xlat(ix)   = xlat(i,j)
        Grid(nb)%xlat_d(ix) = xlat(i,j) * rad2deg
        Grid(nb)%xlon_d(ix) = xlon(i,j) * rad2deg
        Grid(nb)%sinlat(ix) = sin(Grid(nb)%xlat(ix))
        Grid(nb)%coslat(ix) = sqrt(1.0_kind_phys - Grid(nb)%sinlat(ix)*Grid(nb)%sinlat(ix))
        Grid(nb)%area(ix)   = area(i,j)
        Grid(nb)%dx(ix)     = sqrt(area(i,j))
      enddo
    enddo

  end subroutine GFS_grid_populate

end module GFS_driver

