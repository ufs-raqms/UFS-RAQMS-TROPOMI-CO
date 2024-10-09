!!!!!  ==========================================================  !!!!!
!!!!!            'module_radiation_aerosols' description           !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   this module contains climatological atmospheric aerosol schemes for!
!   radiation computations.                                            !
!                                                                      !
!   in the module, the externally callable subroutines are :           !
!                                                                      !
!      'aer_init'   -- initialization                                  !
!         inputs:                                                      !
!           ( NLAY, me )                                               !
!         outputs:                                                     !
!           ( none )                                                   !
!                                                                      !
!      'aer_update' -- updating aerosol data                           !
!         inputs:                                                      !
!           ( iyear, imon, me )                                        !
!         outputs:                                                     !
!           ( none )                                                   !
!                                                                      !
!      'setaer'     -- mapping aeros profile, compute aeros opticals   !
!         inputs:                                                      !
!           (prsi,prsl,prslk,tvly,rhlay,slmsk,tracer,xlon,xlat,        !
!            IMAX,NLAY,NLP1, lsswr,lslwr,                              !
!         outputs:                                                     !
!           (aerosw,aerolw)                                            !
!!          (aerosw,aerolw,aerodp)                                     !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!       'module physpara'                in 'physpara.f'               !
!       'module physcons'                in 'physcons.f'               !
!       'module module_radsw_parameters' in 'radsw_xxxx#_param.f'      !
!       'module module_radlw_parameters' in 'radlw_xxxx#_param.f'      !
!       'module module_radlw_cntr_para'  in 'radsw_xxxx#_param.f'      !
!                                                                      !
!   output variable definitions:                                       !
!       aerosw(IMAX,NLAY,NBDSW,1) - aerosols optical depth for sw      !
!       aerosw(IMAX,NLAY,NBDSW,2) - aerosols single scat albedo for sw !
!       aerosw(IMAX,NLAY,NBDSW,3) - aerosols asymmetry parameter for sw!
!                                                                      !
!       aerolw(IMAX,NLAY,NBDLW,1) - aerosols optical depth for lw      !
!       aerolw(IMAX,NLAY,NBDLW,2) - aerosols single scattering albedo  !
!       aerolw(IMAX,NLAY,NBDLW,3) - aerosols asymetry parameter        !
!                                                                      !
!                                                                      !
!   program history:                                                   !
!     apr     2003  ---  y.-t. hou     created                         !
!     nov 04, 2003  ---  y.-t. hou     modified version                !
!     apr 15, 2005  ---  y.-t. hou     modified module structure       !
!     jul     2006  ---  y.-t. hou     add volcanic forcing            !
!     feb     2007  ---  y.-t. hou     add generalized spectral band   !
!                   interpolation for sw aerosol optical properties    !
!     mar     2007  ---  y.-t. hou     add generalized spectral band   !
!                   interpolation for lw aerosol optical properties    !
!     aug     2007  ---  y.-t. hou     change clim-aer vert domain     !
!                   from pressure reference to sigma reference         !
!     sep     2007  ---  y.-t. hou     moving temporary allocatable    !
!                   module variable arrays to subroutine dynamically   !
!                   allocated arrays (eliminate deallocate calls)      !
!     jan     2010  ---  sarah lu      add gocart option               !
!     may     2010  ---  sarah lu      add geos4-gocart climo          !
!     jul     2010  --   s. moorthi - merged NEMS version with new GFS !
!                        version                                       !
!     oct 23, 2010  ---  Hsin-mu lin   modified subr setclimaer to     !
!        interpolate the 5 degree aerosol data to small domain based on!
!        the nearby 4 points instead of previous nearby assignment by  !
!        using the 5 degree data. this process will eliminate the dsw  !
!        jagged edges in the east conus where aerosol effect are lagre.!
!     dec     2010  ---  y.-t. hou     modified and optimized bi-linear!
!        horizontal interpolation in subr setclimaer. added safe guard !
!        measures in lat/lon indexing and added sea/land mask variable !
!        slmsk as input field to help aerosol profile selection.       !
!     jan     2011  ---  y.-t. hou     divided the program into two    !
!        separated interchangeable modules: a climatology aerosol      !
!        module, and a gocart aerosol scheme module. the stratospheric !
!        volcanic aerosol part is still within the two driver modules, !
!        and may also become a separate one in the further development.!
!        unified in/out argument list for both clim and gocart types of!
!        schemes and added vertically integrated aer-opt-dep, aerodp,  !
!        to replace tau_gocart as optional output for various species. !
!     aug     2012  ---  y.-t. hou     changed the initialization subr !
!        'aerinit' into two parts: 'aer_init' is called at the start   !
!        of run to set up module parameters; and 'aer_update' is       !
!        called within the time loop to check and update data sets.    !
!     nov     2012  ---  y.-t. hou     modified control parameters thru!
!        module 'physpara'.                                            !
!     jan     2013  ---  sarah lu and y.-t. hou   reintegrate both     !
!        opac-clim and gocart schemes into one module to make the      !
!        program best utilize common components. added aerosol model   !
!        scheme selection control variable iaer_mdl to the namelist.   !
!                                                                      !
!   references for opac climatological aerosols:                       !
!     hou et al. 2002  (ncep office note 441)                          !
!     hess et al. 1998 - bams v79 831-844                              !
!                                                                      !
!   references for gocart interactive aerosols:                        !
!     chin et al., 2000 - jgr, v105, 24671-24687                       !
!                                                                      !
!   references for stratosperic volcanical aerosols:                   !
!     sato et al. 1993 - jgr, v98, d12, 22987-22994                    !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!=============================================!
      module module_radiation_aerosols_nmmb   !
!.............................................!
!
      USE ESMF

      use physparam,only : iaermdl, iaerflg, lavoflg, lalwflg, laswflg, &
     &                     lalw1bd, aeros_file, ivflip, kind_phys
      use physcons, only : con_pi, con_rd, con_g, con_t0c, con_c,       &
     &                     con_boltz, con_plnk, con_fvirt

      use module_iounitdef,        only : NIAERCM
      use module_radsw_parameters, only : NBDSW,  wvnsw1=>wvnum1,       &
     &                                    NSWSTR, wvnsw2=>wvnum2
      use module_radlw_parameters, only : NBDLW,  wvnlw1, wvnlw2
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGAER='NCEP-Radiation_aerosols  v5.2  Jan 2013 '
!    &   VTAGAER='NCEP-Radiation_aerosols  v5.1  Nov 2012 '
!    &   VTAGAER='NCEP-Radiation_aerosols  v5.0  Aug 2012 '

!  ---  general use parameter constants:
      integer, parameter, public :: NF_AESW = 3     ! num of output fields for sw rad
      integer, parameter, public :: NF_AELW = 3     ! num of output fields for lw rad
      integer, parameter, public :: NLWSTR  = 1     ! starting band number in ir region
      integer, parameter, public :: NSPC    = 5     ! num of species for output aod (opnl)
      integer, parameter, public :: NSPC1   = NSPC + 1   ! total + species

      real (kind=kind_phys), parameter :: f_zero = 0.0
      real (kind=kind_phys), parameter :: f_one  = 1.0

!  ---  module control parameters set in subroutine "aer_init"
      integer, save :: NSWBND  = NBDSW       ! number of actual bands for sw aerosols
                                             ! calculated according to laswflg setting
      integer, save :: NLWBND  = NBDLW       ! number of actual bands for lw aerosols
                                             ! calculated according to lalwflg and lalw1bd settings
      integer, save :: NSWLWBD = NBDSW+NBDLW ! total number of bands for sw+lw aerosols

! --------------------------------------------------------------------- !
!   section-1 : module variables for spectral band interpolation        !
!               similar to gfdl-sw treatment (2000 version)             !
! --------------------------------------------------------------------- !

!  ---  parameter constants:
      integer, parameter, public :: NWVSOL  = 151   ! num of wvnum regions where solar
                                                    ! flux is constant
      integer, parameter, public :: NWVTOT  = 57600 ! total num of wvnum included
      integer, parameter, public :: NWVTIR  = 4000  ! total num of wvnum in ir range

!  ---  number of wavenumbers in each region where the solar flux is constant
      integer, dimension(NWVSOL), save :: nwvns0

      data nwvns0   / 100,  11,  14,  18,  24,  33,  50,  83,  12,  12, &
     &  13,  15,  15,  17,  18,  20,  21,  24,  26,  30,  32,  37,  42, &
     &  47,  55,  64,  76,  91, 111, 139, 179, 238, 333,  41,  42,  45, &
     &  46,  48,  51,  53,  55,  58,  61,  64,  68,  71,  75,  79,  84, &
     &  89,  95, 101, 107, 115, 123, 133, 142, 154, 167, 181, 197, 217, &
     & 238, 263, 293, 326, 368, 417, 476, 549, 641, 758, 909, 101, 103, &
     & 105, 108, 109, 112, 115, 117, 119, 122, 125, 128, 130, 134, 137, &
     & 140, 143, 147, 151, 154, 158, 163, 166, 171, 175, 181, 185, 190, &
     & 196, 201, 207, 213, 219, 227, 233, 240, 248, 256, 264, 274, 282, &
     & 292, 303, 313, 325, 337, 349, 363, 377, 392, 408, 425, 444, 462, &
     & 483, 505, 529, 554, 580, 610, 641, 675, 711, 751, 793, 841, 891, &
     & 947,1008,1075,1150,1231,1323,1425,1538,1667,1633,14300 /

!  ---  solar flux (w/m**2) in each wvnumb region where it is constant
      real (kind=kind_phys), dimension(NWVSOL), save :: s0intv

      data  s0intv(  1: 50)       /                                     &
     &     1.60000E-6, 2.88000E-5, 3.60000E-5, 4.59200E-5, 6.13200E-5,  &
     &     8.55000E-5, 1.28600E-4, 2.16000E-4, 2.90580E-4, 3.10184E-4,  &
     &     3.34152E-4, 3.58722E-4, 3.88050E-4, 4.20000E-4, 4.57056E-4,  &
     &     4.96892E-4, 5.45160E-4, 6.00600E-4, 6.53600E-4, 7.25040E-4,  &
     &     7.98660E-4, 9.11200E-4, 1.03680E-3, 1.18440E-3, 1.36682E-3,  &
     &     1.57560E-3, 1.87440E-3, 2.25500E-3, 2.74500E-3, 3.39840E-3,  &
     &     4.34000E-3, 5.75400E-3, 7.74000E-3, 9.53050E-3, 9.90192E-3,  &
     &     1.02874E-2, 1.06803E-2, 1.11366E-2, 1.15830E-2, 1.21088E-2,  &
     &     1.26420E-2, 1.32250E-2, 1.38088E-2, 1.44612E-2, 1.51164E-2,  &
     &     1.58878E-2, 1.66500E-2, 1.75140E-2, 1.84450E-2, 1.94106E-2 /
      data  s0intv( 51:100)       /                                     &
     &     2.04864E-2, 2.17248E-2, 2.30640E-2, 2.44470E-2, 2.59840E-2,  &
     &     2.75940E-2, 2.94138E-2, 3.13950E-2, 3.34800E-2, 3.57696E-2,  &
     &     3.84054E-2, 4.13490E-2, 4.46880E-2, 4.82220E-2, 5.22918E-2,  &
     &     5.70078E-2, 6.19888E-2, 6.54720E-2, 6.69060E-2, 6.81226E-2,  &
     &     6.97788E-2, 7.12668E-2, 7.27100E-2, 7.31610E-2, 7.33471E-2,  &
     &     7.34814E-2, 7.34717E-2, 7.35072E-2, 7.34939E-2, 7.35202E-2,  &
     &     7.33249E-2, 7.31713E-2, 7.35462E-2, 7.36920E-2, 7.23677E-2,  &
     &     7.25023E-2, 7.24258E-2, 7.20766E-2, 7.18284E-2, 7.32757E-2,  &
     &     7.31645E-2, 7.33277E-2, 7.36128E-2, 7.33752E-2, 7.28965E-2,  &
     &     7.24924E-2, 7.23307E-2, 7.21050E-2, 7.12620E-2, 7.10903E-2 /
      data  s0intv(101:151)       /                        7.12714E-2,  &
     &     7.08012E-2, 7.03752E-2, 7.00350E-2, 6.98639E-2, 6.90690E-2,  &
     &     6.87621E-2, 6.52080E-2, 6.65184E-2, 6.60038E-2, 6.47615E-2,  &
     &     6.44831E-2, 6.37206E-2, 6.24102E-2, 6.18698E-2, 6.06320E-2,  &
     &     5.83498E-2, 5.67028E-2, 5.51232E-2, 5.48645E-2, 5.12340E-2,  &
     &     4.85581E-2, 4.85010E-2, 4.79220E-2, 4.44058E-2, 4.48718E-2,  &
     &     4.29373E-2, 4.15242E-2, 3.81744E-2, 3.16342E-2, 2.99615E-2,  &
     &     2.92740E-2, 2.67484E-2, 1.76904E-2, 1.40049E-2, 1.46224E-2,  &
     &     1.39993E-2, 1.19574E-2, 1.06386E-2, 1.00980E-2, 8.63808E-3,  &
     &     6.52736E-3, 4.99410E-3, 4.39350E-3, 2.21676E-3, 1.33812E-3,  &
     &     1.12320E-3, 5.59000E-4, 3.60000E-4, 2.98080E-4, 7.46294E-5  /

! --------------------------------------------------------------------- !
!   section-2 : module variables for stratospheric volcanic aerosols    !
!               from historical data (sato et al. 1993)                 !
! --------------------------------------------------------------------- !

!  ---  parameter constants:
      integer, parameter :: MINVYR = 1850    ! lower lim (year) data available
      integer, parameter :: MAXVYR = 1999    ! upper lim (year) data available

!  ---  monthly, 45-deg lat-zone aerosols data set in subroutine 'aer_init'
      integer, allocatable, save :: ivolae(:,:,:)

!  ---  static control variables:
      integer :: kyrstr, kyrend, kyrsav, kmonsav

! --------------------------------------------------------------------- !
!   section-3 : module variables for opac climatological aerosols       !
!               optical properties (hess et al. 1989)                   !
! --------------------------------------------------------------------- !

!  ---  parameters and constants:
      integer, parameter :: NXC = 5    ! num of max componets in a profile
      integer, parameter :: NAE = 7    ! num of aerosols profile structures
      integer, parameter :: NDM = 5    ! num of atmos aerosols domains
      integer, parameter :: IMXAE = 72 ! num of lon-points in glb aeros data set
      integer, parameter :: JMXAE = 37 ! num of lat-points in glb aeros data set
      integer, parameter :: NAERBND=61 ! num of bands for clim aer data (opac)
      integer, parameter :: NRHLEV =8  ! num of rh levels for rh-dep components
      integer, parameter :: NCM1 = 6   ! num of rh independent aeros species
      integer, parameter :: NCM2 = 4   ! num of rh dependent aeros species
      integer, parameter :: NCM  = NCM1+NCM2

      real (kind=kind_phys), dimension(NRHLEV), save :: rhlev
      data  rhlev (:) / 0.0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99 /

!  ---  the following arrays are for climatological data that are
!           allocated and read in subroutine 'clim_aerinit'.
!   - global aerosol distribution:
!      haer  (NDM,NAE)  - scale height of aerosols (km)
!      prsref(NDM,NAE)  - ref pressure lev (sfc to toa) in mb (100Pa)
!      sigref(NDM,NAE)  - ref sigma lev (sfc to toa)

      real (kind=kind_phys), save, dimension(NDM,NAE) :: haer, prsref,  &
     &                                                   sigref

!  ---  the following arrays are allocate and setup in subr 'clim_aerinit'
!   - for relative humidity independent aerosol optical properties:
!      species : insoluble        (inso); soot             (soot);
!                mineral nuc mode (minm); mineral acc mode (miam);
!                mineral coa mode (micm); mineral transport(mitr).
!      extrhi(NCM1,NSWLWBD) - extinction coefficient for sw+lw spectral band
!      scarhi(NCM1,NSWLWBD) - scattering coefficient for sw+lw spectral band
!      ssarhi(NCM1,NSWLWBD) - single scattering albedo for sw+lw spectral band
!      asyrhi(NCM1,NSWLWBD) - asymmetry parameter for sw+lw spectral band
!   - for relative humidity dependent aerosol optical properties:
!      species : water soluble    (waso); sea salt acc mode(ssam);
!                sea salt coa mode(sscm); sulfate droplets (suso).
!      rh level: 00%, 50%, 70%, 80%, 90%, 95%, 98%, 99%
!      extrhd(NRHLEV,NCM2,NSWLWBD) - extinction coefficient for sw+lw band
!      scarhd(NRHLEV,NCM2,NSWLWBD) - scattering coefficient for sw+lw band
!      ssarhd(NRHLEV,NCM2,NSWLWBD) - single scattering albedo for sw+lw band
!      asyrhd(NRHLEV,NCM2,NSWLWBD) - asymmetry parameter for sw+lw band
!   - for stratospheric aerosols optical properties:
!      extstra(NSWLWBD)            - extinction coefficient for sw+lw band

      real (kind=kind_phys), allocatable, save, dimension(:,:)   ::     &
     &       extrhi, scarhi, ssarhi, asyrhi
      real (kind=kind_phys), allocatable, save, dimension(:,:,:) ::     &
     &       extrhd, scarhd, ssarhd, asyrhd
      real (kind=kind_phys), allocatable, save, dimension(:)     ::     &
     &       extstra

!  ---  the following arrays are calculated in subr 'clim_aerinit'
!   - for topospheric aerosol profile distibution:
!      kprfg (    IMXAE*JMXAE)   - aeros profile index
!      idxcg (NXC*IMXAE*JMXAE)   - aeros component index
!      cmixg (NXC*IMXAE*JMXAE)   - aeros component mixing ratio
!      denng ( 2 *IMXAE*JMXAE)   - aerosols number density

      real (kind=kind_phys), dimension(NXC,IMXAE,JMXAE), save :: cmixg
      real (kind=kind_phys), dimension( 2 ,IMXAE,JMXAE), save :: denng
      integer,               dimension(NXC,IMXAE,JMXAE), save :: idxcg
      integer,               dimension(    IMXAE,JMXAE), save :: kprfg

! --------------------------------------------------------------------- !
!   section-4 : module variables for gocart aerosol optical properties  !
! --------------------------------------------------------------------- !

!  ---  parameters and constants:
!   - KCM, KCM1, KCM2 are determined from subroutine 'set_aerspc'
!     integer, parameter :: KAERBND=61 ! num of bands for aer data (gocart)
!     integer, parameter :: KRHLEV =36 ! num of rh levels for rh-dep components
!*    integer, parameter :: KCM1 = 8   ! num of rh independent aer !species
!*    integer, parameter :: KCM2 = 5   ! num of rh dependent aer species
!*    integer, parameter :: KCM  = KCM1 + KCM2
!     integer, save      :: KCM1 = 0   ! num of rh indep aerosols (set in subr set_aerspc)
!     integer, save      :: KCM2 = 0   ! num of rh dep aerosols   (set in subr set_aerspc)
!     integer, save      :: KCM        ! =KCM1+KCM2               (set in subr set_aerspc)

!     real (kind=kind_phys), dimension(KRHLEV) :: rhlev_grt             &
!     data  rhlev_grt (:)/ .00, .05, .10, .15, .20, .25, .30, .35,      &
!    &      .40, .45, .50, .55, .60, .65, .70, .75, .80, .81, .82,      &
!    &      .83, .84, .85, .86, .87, .88, .89, .90, .91, .92, .93,      &
!    &      .94, .95, .96, .97, .98, .99 /

!  ---  the following arrays are allocate and setup in subr 'gocrt_aerinit'
!  ------  gocart aerosol specification    ------
!  =>  transported aerosol species:
!      DU (5-bins)
!      SS (4 bins for climo mode and 5 bins for fcst mode)
!      SU (dms, so2, so4, msa)
!      OC (phobic, philic) and BC (phobic, philic)
!  =>  species and lumped species for aerosol optical properties
!      DU (5-bins, with 4 sub-groups in the submicron bin )
!      SS (ssam for submicron, sscm for coarse mode)
!      SU (so4)
!      OC (phobic, philic) and BC (phobic, philic)
!  =>  specification used for aerosol optical properties luts
!      DU (8 bins)
!      SS (ssam, sscm)
!      SU (suso)
!      OC (waso) and BC (soot)
!
!   - spectral band structure:
!      iendwv_grt(KAERBND)      - ending wavenumber (cm**-1) for each band
!   - relative humidity independent aerosol optical properties:
!   ===> species : dust (8 bins)
!      rhidext0_grt(KAERBND,KCM1) - extinction coefficient
!      rhidssa0_grt(KAERBND,KCM1) - single scattering albedo
!      rhidasy0_grt(KAERBND,KCM1) - asymmetry parameter
!   - relative humidity dependent aerosol optical properties:
!   ===> species : soot, suso, waso, ssam, sscm
!      rhdpext0_grt(KAERBND,KRHLEV,KCM2) - extinction coefficient
!      rhdpssa0_grt(KAERBND,KRHLEV,KCM2) - single scattering albedo
!      rhdpasy0_grt(KAERBND,KRHLEV,KCM2) - asymmetry parameter

!     integer,               allocatable, dimension(:) :: iendwv_grt
!     real (kind=kind_phys), allocatable, dimension(:,:)  ::            &
!    &                       rhidext0_grt, rhidssa0_grt, rhidasy0_grt
!     real (kind=kind_phys), allocatable, dimension(:,:,:)::            &
!    &                       rhdpext0_grt, rhdpssa0_grt, rhdpasy0_grt

!   - relative humidity independent aerosol optical properties:
!      extrhi_grt(KCM1,NSWLWBD) - extinction coefficient for sw+lw spectral band
!      ssarhi_grt(KCM1,NSWLWBD) - single scattering albedo for sw+lw spectral band
!      asyrhi_grt(KCM1,NSWLWBD) - asymmetry parameter for sw+lw spectral band
!   - relative humidity dependent aerosol optical properties:
!      extrhd_grt(KRHLEV,KCM2,NSWLWBD) - extinction coefficient for sw+lw band
!      ssarhd_grt(KRHLEV,KCM2,NSWLWBD) - single scattering albedo for sw+lw band
!      asyrhd_grt(KRHLEV,KCM2,NSWLWBD) - asymmetry parameter for sw+lw band

!     real (kind=kind_phys), allocatable, save, dimension(:,:)   ::     &
!    &       extrhi_grt, ssarhi_grt, asyrhi_grt
!     real (kind=kind_phys), allocatable, save, dimension(:,:,:) ::     &
!    &       extrhd_grt, ssarhd_grt, asyrhd_grt

! --------------------------------------------------------------------- !
!   section-5 : module variables for gocart aerosol climo data set      !
! --------------------------------------------------------------------- !
!     This version only supports geos3-gocart data set (Jan 2010)
!     Modified to support geos4-gocart data set        (May 2010)
!
!  geos3-gocart vs geos4-gocart
!  (1) Use the same module variables
!      IMXG,JMXG,KMXG,NMXG,psclmg,dmclmg,geos_rlon,geos_rlat
!  (2) Similarity between geos3 and geos 4:
!      identical lat/lon grids and aerosol specification;
!      direction of vertical index is bottom-up (sfc to toa)
!  (3) Difference between geos3 and geos4
!      vertical coordinate (sigma for geos3/hybrid_sigma_pressure for geos4)
!      aerosol units (mass concentration for geos3/mixing ratio for geos4)

!     integer, parameter :: IMXG = 144 ! num of lon-points in geos dataset
!     integer, parameter :: JMXG = 91  ! num of lat-points in geos dataset
!     integer, parameter :: KMXG = 30  ! num of vertical layers in geos dataset
!*    integer, parameter :: NMXG = 12  ! num of gocart aer spec for opt calc
!     integer, save      :: NMXG       ! to be determined by set_aerspc

!     real (kind=kind_phys), parameter :: dltx = 360.0 / float(IMXG)
!     real (kind=kind_phys), parameter :: dlty = 180.0 / float(JMXG-1)

!  --- the following arrays are allocated and setup in 'rd_gocart_clim'
!   - geos-gocart climo data (input dataset)
!     psclmg  - pressure in cb                   IMXG*JMXG*KMXG
!     dmclmg  - aerosol dry mass in g/m3         IMXG*JMXG*KMXG*NMXG
!               or aerosol mixing ratio in mol/mol or Kg/Kg

!     real (kind=kind_phys),allocatable, save:: psclmg(:,:,:),          &
!    &                                          dmclmg(:,:,:,:)

!   - geos-gocart lat/lon arrays
!     real (kind=kind_phys), allocatable, save, dimension(:)   ::       &
!    &                       geos_rlon, geos_rlat

!   - control flag for gocart climo data set
!     xxxx as default; ver3 for geos3; ver4 for geos4; 0000 for unknown data
!     character*4, save  :: gocart_climo = 'xxxx'

!   - molecular wght of gocart aerosol species
!     real (kind=kind_io4), allocatable :: molwgt(:)

!! ---  the following are for diagnostic purpose to output aerosol optical depth
!       aod from 10 components are grouped into 5 major different species:
!      1:dust (inso,minm,miam,micm,mitr); 2:black carbon (soot)
!      3:water soluble (waso);            4:sulfate (suso);      5:sea salt (ssam,sscm)
!
!      idxspc (NCM)         - index conversion array
!      lspcaod              - logical flag for aod from individual species
!
!     integer, dimension(NCM) :: idxspc
!     data  idxspc / 1, 2, 1, 1, 1, 1, 3, 5, 5, 4 /
!     logical, save :: lspcaod = .false.
!
!   - wvn550 is the wavenumber (1/cm) of wavelenth 550nm for diagnostic aod output
!     nv_aod is the sw spectral band covering wvn550 (comp in aer_init)
!
!     real (kind=kind_phys), parameter :: wvn550 = 1.0e4/0.55
!     integer, save      :: nv_aod = 1

!  ---  public interfaces

      public aer_init, aer_update, setaer


! =================
      contains
! =================

!-----------------------------------
      subroutine aer_init                                               &
!...................................
!  ---  inputs:
     &     ( NLAY, me )
!  ---  outputs: ( to module variables )

!  ==================================================================  !
!                                                                      !
!  aer_init is the initialization program to set up necessary          !
!    parameters and working arrays.                                    !
!                                                                      !
!  inputs:                                                             !
!     NLAY    - number of model vertical layers  (not used)            !
!     me      - print message control flag                             !
!                                                                      !
!  outputs: (to module variables)                                      !
!                                                                      !
!  external module variables: (in physpara)                            !
!     iaermdl - tropospheric aerosol model scheme flag                 !
!               =0 opac-clim; =1 gocart-clim, =2 gocart-prognostic     !
!     lalwflg - logical lw aerosols effect control flag                !
!               =t compute lw aerosol optical prop                     !
!     laswflg - logical sw aerosols effect control flag                !
!               =t compute sw aerosol optical prop                     !
!     lavoflg - logical stratosphere volcanic aerosol control flag     !
!               =t include volcanic aerosol effect                     !
!     lalw1bd = logical lw aeros propty 1 band vs multi-band cntl flag !
!               =t use 1 broad band optical property                   !
!               =f use multi bands optical property                    !
!                                                                      !
!  module constants:                                                   !
!     NWVSOL  - num of wvnum regions where solar flux is constant      !
!     NWVTOT  - total num of wave numbers used in sw spectrum          !
!     NWVTIR  - total num of wave numbers used in the ir region        !
!     NSWBND  - total number of sw spectral bands                      !
!     NLWBND  - total number of lw spectral bands                      !
!                                                                      !
!  usage:    call aer_init                                             !
!                                                                      !
!  subprograms called:  clim_aerinit, gcrt_aerinit,                    !
!                       wrt_aerlog, set_volcaer, set_spectrum,         !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      integer,  intent(in) :: NLAY, me

!  ---  output: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(NWVTOT) :: solfwv        ! one wvn sol flux
      real (kind=kind_phys), dimension(NWVTIR) :: eirfwv        ! one wvn ir flux
!
!===>  ...  begin here
!
      kyrstr  = 1
      kyrend  = 1
      kyrsav  = 1
      kmonsav = 1

!  --- ...  write aerosol parameter configuration to output logs

      if ( me == 0 ) then

        call wrt_aerlog      ! write aerosol param info to log file
!  ---  inputs:   (in scope variables)
!  ---  outputs:  ( none )

      endif

      if ( iaerflg == 0 ) return      ! return without any aerosol calculations

!  --- ...  in sw, aerosols optical properties are computed for each radiation
!           spectral band; while in lw, optical properties can be calculated
!           for either only one broad band or for each of the lw radiation bands

      if ( laswflg ) then
        NSWBND = NBDSW
      else
        NSWBND = 0
      endif

      if ( lalwflg ) then
        if ( lalw1bd ) then
          NLWBND = 1
        else
          NLWBND = NBDLW
        endif
      else
        NLWBND = 0
      endif

      NSWLWBD = NSWBND + NLWBND

      if ( iaerflg /= 100 ) then

!  --- ...  set up spectral one wavenumber solar/ir fluxes

        call set_spectrum
!  ---  inputs:   (module constants)
!  ---  outputs:  (in-scope variables)

!  --- ...  invoke tropospheric aerosol initialization

        if ( iaermdl == 0 ) then                    ! opac-climatology scheme

          call clim_aerinit                                             &
!  ---  inputs:
     &     ( solfwv, eirfwv, me                                         &
!  ---  outputs:
     &     )

!       elseif ( iaermdl == 1 ) then                ! gocart-climatology scheme
!       elseif ( iaermdl==1 .or. iaermdl==2 ) then  ! gocart-clim/prog scheme

!         call gcrt_climinit

!       elseif ( iaermdl == 2 ) then                ! gocart-prognostic scheme

!         call gcrt_aerinit

        else
          if ( me == 0 ) then
            print *,'  !!! ERROR in aerosol model scheme selection',    &
     &              ' iaermdl =',iaermdl
            ! stop
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
        endif

      endif    ! end if_iaerflg_block

!  --- ...  invoke stratosperic volcanic aerosol initialization

      if ( lavoflg ) then

        call set_volcaer
!  ---  inputs:  (module variables)
!  ---  outputs: (module variables)

      endif    ! end if_lavoflg_block


! =================
      contains
! =================

!--------------------------------
      subroutine wrt_aerlog
!................................
!  ---  inputs:    (in scope variables)
!  ---  outputs:   ( none )

!  ==================================================================  !
!                                                                      !
!  subprogram : wrt_aerlog                                             !
!                                                                      !
!    write aerosol parameter configuration to run log file.            !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  external module variables:  (in physpara)                           !
!   iaermdl  - aerosol scheme flag: 0:opac-clm; 1:gocart-clim;         !
!                                   2:gocart-prog                      !
!   iaerflg  - aerosol effect control flag: 3-digits (volc,lw,sw)      !
!   lalwflg  - toposphere lw aerosol effect: =f:no; =t:yes             !
!   laswflg  - toposphere sw aerosol effect: =f:no; =t:yes             !
!   lavoflg  - stratospherer volcanic aeros effect: =f:no; =t:yes      !
!                                                                      !
!  outputs: ( none )                                                   !
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call wrt_aerlog                                           !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: ( none )
!  ---  output: ( none )
!  ---  locals:

!
!===>  ...  begin here
!
      print *, VTAGAER    ! print out version tag

      if ( iaermdl == 0 ) then
        print *,' - Using OPAC-seasonal climatology for tropospheric',  &
     &          ' aerosol effect'
      elseif ( iaermdl == 1 ) then
        print *,' - Using GOCART-climatology for tropospheric',         &
     &          ' aerosol effect'
      elseif ( iaermdl == 2 ) then
        print *,' - Using GOCART-prognostic aerosols for tropospheric', &
     &          ' aerosol effect'
      else
        print *,' !!! ERROR in selection of aerosol model scheme',      &
     &          ' IAER_MDL =',iaermdl
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif   ! end_if_iaermdl_block

      print *,'   IAER=',iaerflg,'  LW-trop-aer=',lalwflg,              &
     &        '  SW-trop-aer=',laswflg,'  Volc-aer=',lavoflg

      if ( iaerflg <= 0 ) then        ! turn off all aerosol effects
        print *,' - No tropospheric/volcanic aerosol effect included'
        print *,'      Input values of aerosol optical properties to'   &
     &         ,' both SW and LW radiations are set to zeros'
      else
        if ( iaerflg >= 100 ) then    ! incl stratospheric volcanic aerosols
          print *,' - Include stratospheric volcanic aerosol effect'
        else                       ! no stratospheric volcanic aerosols
          print *,' - No stratospheric volcanic aerosol effect'
        endif

        if ( laswflg ) then          ! chcek for sw effect
          print *,'   - Compute multi-band aerosol optical'             &
     &           ,' properties for SW input parameters'
        else
          print *,'   - No SW radiation aerosol effect, values of'      &
     &           ,' aerosol properties to SW input are set to zeros'
        endif

        if ( lalwflg ) then          ! check for lw effect
          if ( lalw1bd ) then
            print *,'   - Compute 1 broad-band aerosol optical'         &
     &           ,' properties for LW input parameters'
          else
            print *,'   - Compute multi-band aerosol optical'           &
     &           ,' properties for LW input parameters'
          endif
        else
          print *,'   - No LW radiation aerosol effect, values of'      &
     &           ,' aerosol properties to LW input are set to zeros'
        endif
      endif     ! end if_iaerflg_block
!
      return
!................................
      end subroutine wrt_aerlog
!--------------------------------


!--------------------------------
      subroutine set_spectrum
!................................
!  ---  inputs:   (module constants)
!  ---  outputs:  (in-scope variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : set_spectrum                                           !
!                                                                      !
!    define the one wavenumber solar fluxes based on toa solar spectral!
!    distrobution, and define the one wavenumber ir fluxes based on    !
!    black-body emission distribution at a predefined temperature.     !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  inputs:  (module constants)                                         !
!   NWVTOT           - total num of wave numbers used in sw spectrum   !
!   NWVTIR           - total num of wave numbers used in the ir region !
!                                                                      !
!  outputs: (in-scope variables)                                       !
!   solfwv(NWVTOT)   - solar flux for each individual wavenumber (w/m2)!
!   eirfwv(NWVTIR)   - ir flux(273k) for each individual wavenum (w/m2)!
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call set_spectrum                                         !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: (module constants)
!     integer :: NWVTOT, NWVTIR

!  ---  output: (in-scope variables)
!     real (kind=kind_phys), dimension(NWVTOT) :: solfwv        ! one wvn sol flux
!     real (kind=kind_phys), dimension(NWVTIR) :: eirfwv        ! one wvn ir flux

!  ---  locals:
      real (kind=kind_phys) :: soltot, tmp1, tmp2, tmp3

      integer :: nb, nw, nw1, nw2, nmax, nmin
!
!===>  ...  begin here
!
!     nmax = min( NWVTOT, nint( maxval(wvnsw2) ))
!     nmin = max( 1,      nint( minval(wvnsw1) ))

!  ---  check print
!     print *,' MINWVN, MAXWVN = ',nmin, nmax
!  --- ...  define the one wavenumber solar fluxes based on toa solar
!           spectral distribution

!     soltot1 = f_zero
!     soltot  = f_zero
      do nb = 1, NWVSOL
        if ( nb == 1 ) then
          nw1 = 1
        else
          nw1 = nw1 + nwvns0(nb-1)
        endif

        nw2 = nw1 + nwvns0(nb) - 1

        do nw = nw1, nw2
          solfwv(nw) = s0intv(nb)
!         soltot1 = soltot1 + s0intv(nb)
!         if ( nw >= nmin .and. nw <= nmax ) then
!           soltot = soltot + s0intv(nb)
!         endif
        enddo
      enddo

!  --- ...  define the one wavenumber ir fluxes based on black-body
!           emission distribution at a predefined temperature

      tmp1 = 2.0 * con_pi * con_plnk * (con_c**2)
      tmp2 = con_plnk * con_c / (con_boltz * con_t0c)

      do nw = 1, NWVTIR
        tmp3 = 100.0 * nw
        eirfwv(nw) = (tmp1 * tmp3**3) / (exp(tmp2*tmp3) - 1.0)
      enddo
!
      return
!................................
      end subroutine set_spectrum
!--------------------------------


!-----------------------------
      subroutine set_volcaer
!.............................
!  ---  inputs:   ( none )
!  ---  outputs:  (module variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : set_volcaer                                            !
!                                                                      !
!    this is the initialization progrmam for stratospheric volcanic    !
!    aerosols.                                                         !
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call set_volcaer                                          !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: (none)

!  ---  output: (module variables)
!     integer :: ivolae(:,:,:)

!  ---  locals:
!
!===>  ...  begin here
!
!  ---  allocate data space

      if ( .not. allocated(ivolae) ) then
        allocate ( ivolae(12,4,10) )   ! for 12-mon,4-lat_zone,10-year
      endif
!
      return
!................................
      end subroutine set_volcaer
!--------------------------------
!
!...................................
      end subroutine aer_init
!-----------------------------------


!-----------------------------------
      subroutine clim_aerinit                                           &
!...................................
!  ---  inputs:
     &     ( solfwv, eirfwv, me                                         &
!  ---  outputs:
     &     )

!  ==================================================================  !
!                                                                      !
!  clim_aerinit is the opac-climatology aerosol initialization program !
!  to set up necessary parameters and working arrays.                  !
!                                                                      !
!  inputs:                                                             !
!   solfwv(NWVTOT)   - solar flux for each individual wavenumber (w/m2)!
!   eirfwv(NWVTIR)   - ir flux(273k) for each individual wavenum (w/m2)!
!   me               - print message control flag                      !
!                                                                      !
!  outputs: (to module variables)                                      !
!                                                                      !
!  external module variables: (in physpara)                            !
!     iaerflg - abc 3-digit integer aerosol flag (abc:volc,lw,sw)      !
!               a: =0 use background stratospheric aerosol             !
!                  =1 incl stratospheric vocanic aeros (MINVYR-MAXVYR) !
!               b: =0 no topospheric aerosol in lw radiation           !
!                  =1 include tropspheric aerosols for lw radiation    !
!               c: =0 no topospheric aerosol in sw radiation           !
!                  =1 include tropspheric aerosols for sw radiation    !
!     lalwflg - logical lw aerosols effect control flag                !
!               =t compute lw aerosol optical prop                     !
!     laswflg - logical sw aerosols effect control flag                !
!               =t compute sw aerosol optical prop                     !
!     lalw1bd = logical lw aeros propty 1 band vs multi-band cntl flag !
!               =t use 1 broad band optical property                   !
!               =f use multi bands optical property                    !
!                                                                      !
!  module constants:                                                   !
!     NWVSOL  - num of wvnum regions where solar flux is constant      !
!     NWVTOT  - total num of wave numbers used in sw spectrum          !
!     NWVTIR  - total num of wave numbers used in the ir region        !
!     NSWBND  - total number of sw spectral bands                      !
!     NLWBND  - total number of lw spectral bands                      !
!     NAERBND - number of bands for climatology aerosol data           !
!     NCM1    - number of rh independent aeros species                 !
!     NCM2    - number of rh dependent aeros species                   !
!                                                                      !
!  usage:    call clim_aerinit                                         !
!                                                                      !
!  subprograms called:  set_aercoef, optavg                            !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      real (kind=kind_phys), dimension(:) :: solfwv        ! one wvn sol flux
      real (kind=kind_phys), dimension(:) :: eirfwv        ! one wvn ir flux

      integer,  intent(in) :: me

!  ---  output: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(NAERBND,NCM1)       ::           &
     &       rhidext0, rhidsca0, rhidssa0, rhidasy0
      real (kind=kind_phys), dimension(NAERBND,NRHLEV,NCM2)::           &
     &       rhdpext0, rhdpsca0, rhdpssa0, rhdpasy0
      real (kind=kind_phys), dimension(NAERBND)            :: straext0

      real (kind=kind_phys), dimension(NSWBND,NAERBND) :: solwaer
      real (kind=kind_phys), dimension(NSWBND)         :: solbnd
      real (kind=kind_phys), dimension(NLWBND,NAERBND) :: eirwaer
      real (kind=kind_phys), dimension(NLWBND)         :: eirbnd

      integer, dimension(NSWBND) :: nv1, nv2
      integer, dimension(NLWBND) :: nr1, nr2
!
!===>  ...  begin here
!
!  --- ...  invoke tropospheric aerosol initialization

      call set_aercoef
!  ---  inputs:   (in-scope variables, module constants)
!  ---  outputs:  (module variables)


! =================
      contains
! =================

!--------------------------------
      subroutine set_aercoef
!................................
!  ---  inputs:   (in-scope variables, module constants)
!  ---  outputs:  (module variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : set_aercoef                                            !
!                                                                      !
!    this is the initialization progrmam for climatological aerosols   !
!                                                                      !
!    the program reads and maps the pre-tabulated aerosol optical      !
!    spectral data onto corresponding sw radiation spectral bands.     !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  inputs:  (in-scope variables, module constants)                     !
!   solfwv(:)    - real, solar flux for individual wavenumber (w/m2)   !
!   eirfwv(:)    - real, lw flux(273k) for individual wavenum (w/m2)   !
!   me           - integer, select cpu number as print control flag    !
!                                                                      !
!  outputs: (to the module variables)                                  !
!                                                                      !
!  external module variables:  (in physpara)                           !
!   lalwflg   - module control flag for lw trop-aer: =f:no; =t:yes     !
!   laswflg   - module control flag for sw trop-aer: =f:no; =t:yes     !
!   aeros_file- external aerosol data file name                        !
!                                                                      !
!  internal module variables:                                          !
!     IMXAE   - number of longitude points in global aeros data set    !
!     JMXAE   - number of latitude points in global aeros data set     !
!     wvnsw1,wvnsw2 (NSWSTR:NSWEND)                                    !
!             - start/end wavenumbers for each of sw bands             !
!     wvnlw1,wvnlw2 (     1:NBDLW)                                     !
!             - start/end wavenumbers for each of lw bands             !
!     NSWLWBD - total num of bands (sw+lw) for aeros optical properties!
!     NSWBND  - number of sw spectral bands actually invloved          !
!     NLWBND  - number of lw spectral bands actually invloved          !
!     NIAERCM - unit number for reading input data set                 !
!     extrhi  - extinction coef for rh-indep aeros         NCM1*NSWLWBD!
!     scarhi  - scattering coef for rh-indep aeros         NCM1*NSWLWBD!
!     ssarhi  - single-scat-alb for rh-indep aeros         NCM1*NSWLWBD!
!     asyrhi  - asymmetry factor for rh-indep aeros        NCM1*NSWLWBD!
!     extrhd  - extinction coef for rh-dep aeros    NRHLEV*NCM2*NSWLWBD!
!     scarhd  - scattering coef for rh-dep aeros    NRHLEV*NCM2*NSWLWBD!
!     ssarhd  - single-scat-alb for rh-dep aeros    NRHLEV*NCM2*NSWLWBD!
!     asyrhd  - asymmetry factor for rh-dep aeros   NRHLEV*NCM2*NSWLWBD!
!                                                                      !
!  major local variables:                                              !
!   for handling spectral band structures                              !
!     iendwv   - ending wvnum (cm**-1) for each band  NAERBND          !
!   for handling optical properties of rh independent species (NCM1)   !
!         1. insoluble        (inso); 2. soot             (soot);      !
!         3. mineral nuc mode (minm); 4. mineral acc mode (miam);      !
!         5. mineral coa mode (micm); 6. mineral transport(mitr).      !
!     rhidext0 - extinction coefficient             NAERBND*NCM1       !
!     rhidsca0 - scattering coefficient             NAERBND*NCM1       !
!     rhidssa0 - single scattering albedo           NAERBND*NCM1       !
!     rhidasy0 - asymmetry parameter                NAERBND*NCM1       !
!   for handling optical properties of rh ndependent species (NCM2)    !
!         1. water soluble    (waso); 2. sea salt acc mode(ssam);      !
!         3. sea salt coa mode(sscm); 4. sulfate droplets (suso).      !
!         rh level (NRHLEV): 00%, 50%, 70%, 80%, 90%, 95%, 98%, 99%    !
!     rhdpext0 - extinction coefficient             NAERBND,NRHLEV,NCM2!
!     rhdpsca0 - scattering coefficient             NAERBND,NRHLEV,NCM2!
!     rhdpssa0 - single scattering albedo           NAERBND,NRHLEV,NCM2!
!     rhdpasy0 - asymmetry parameter                NAERBND,NRHLEV,NCM2!
!   for handling optical properties of stratospheric bkgrnd aerosols   !
!     straext0 - extingction coefficients             NAERBND          !
!                                                                      !
!  usage:    call set_aercoef                                          !
!                                                                      !
!  subprograms called:  optavg                                         !
!                                                                      !
!  ==================================================================  !
!
!  ---  inputs:  ( none )
!  ---  output: ( none )

!  ---  locals:
      integer, dimension(NAERBND) :: iendwv

      integer :: i, j, k, m, mb, ib, ii, id, iw, iw1, iw2

      real (kind=kind_phys) :: sumsol, sumir

      logical :: file_exist
      character :: cline*80, ctyp*3
!
!===>  ...  begin here
!
!  --- ...  reading climatological aerosols data

      inquire (file=aeros_file, exist=file_exist)

      if ( file_exist ) then
        close (NIAERCM)
        open  (unit=NIAERCM,file=aeros_file,status='OLD',               &
     &        form='FORMATTED')
        rewind (NIAERCM)
      else
        print *,'    Requested aerosol data file "',aeros_file,         &
     &          '" not found!'
        print *,'    *** Stopped in subroutine aero_init !!'
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif     ! end if_file_exist_block

!  --- ...  skip monthly global distribution

      do m = 1, 12
        read (NIAERCM,12) cline
  12    format(a80/)

        do j = 1, JMXAE
          do i = 1, IMXAE
            read(NIAERCM,*) id
          enddo
        enddo
      enddo   ! end do_m_block

!  --- ...  aloocate and input aerosol optical data

      if ( .not. allocated( extrhi ) ) then
        allocate ( extrhi (       NCM1,NSWLWBD) )
        allocate ( scarhi (       NCM1,NSWLWBD) )
        allocate ( ssarhi (       NCM1,NSWLWBD) )
        allocate ( asyrhi (       NCM1,NSWLWBD) )
        allocate ( extrhd (NRHLEV,NCM2,NSWLWBD) )
        allocate ( scarhd (NRHLEV,NCM2,NSWLWBD) )
        allocate ( ssarhd (NRHLEV,NCM2,NSWLWBD) )
        allocate ( asyrhd (NRHLEV,NCM2,NSWLWBD) )
        allocate ( extstra(            NSWLWBD) )
      endif

      read(NIAERCM,21) cline   ! ending wave num for 61 aeros spectral bands
  21  format(a80)
      read(NIAERCM,22) iendwv(:)
  22  format(13i6)

      read(NIAERCM,21) cline   ! atmos scale height for 5 domains, 7 profs
      read(NIAERCM,24) haer(:,:)
  24  format(20f4.1)

      read(NIAERCM,21) cline   ! reference pressure for 5 domains, 7 profs
      read(NIAERCM,26) prsref(:,:)
  26  format(10f7.2)

      read(NIAERCM,21) cline   ! rh indep ext coef for 61 bands, 6 species
      read(NIAERCM,28) rhidext0(:,:)
  28  format(8e10.3)

      read(NIAERCM,21) cline   ! rh indep sca coef for 61 bands, 6 species
      read(NIAERCM,28) rhidsca0(:,:)

      read(NIAERCM,21) cline   ! rh indep ssa coef for 61 bands, 6 species
      read(NIAERCM,28) rhidssa0(:,:)

      read(NIAERCM,21) cline   ! rh indep asy coef for 61 bands, 6 species
      read(NIAERCM,28) rhidasy0(:,:)

      read(NIAERCM,21) cline   ! rh dep ext coef for 61 bands, 8 rh lev, 4 species
      read(NIAERCM,28) rhdpext0(:,:,:)

      read(NIAERCM,21) cline   ! rh dep sca coef for 61 bands, 8 rh lev, 4 species
      read(NIAERCM,28) rhdpsca0(:,:,:)

      read(NIAERCM,21) cline   ! rh dep ssa coef for 61 bands, 8 rh lev, 4 species
      read(NIAERCM,28) rhdpssa0(:,:,:)

      read(NIAERCM,21) cline   ! rh dep asy coef for 61 bands, 8 rh lev, 4 species
      read(NIAERCM,28) rhdpasy0(:,:,:)

      read(NIAERCM,21) cline   ! stratospheric background aeros for 61 bands
      read(NIAERCM,28) straext0(:)

      close (NIAERCM)

!  --- ...  convert pressure reference level (in mb) to sigma reference level
!           assume an 1000mb reference surface pressure

      sigref(:,:) = 0.001 * prsref(:,:)

!  --- ...  compute solar flux weights and interval indices for mapping
!           spectral bands between sw radiation and aerosol data

      if ( laswflg ) then
        solbnd (:)   = f_zero
        solwaer(:,:) = f_zero

        do ib = 1, NSWBND
          mb = ib + NSWSTR - 1
          ii = 1
          iw1 = nint(wvnsw1(mb))
          iw2 = nint(wvnsw2(mb))

!         if ( wvnsw2(mb)>=wvn550 .and. wvn550>=wvnsw1(mb) ) then
!           nv_aod = ib                  ! sw band number covering 550nm wavelenth
!         endif

          Lab_swdowhile : do while ( iw1 > iendwv(ii) )
            if ( ii == NAERBND ) exit Lab_swdowhile
            ii = ii + 1
          enddo  Lab_swdowhile

          sumsol = f_zero
          nv1(ib) = ii

          do iw = iw1, iw2
            solbnd(ib) = solbnd(ib) + solfwv(iw)
            sumsol = sumsol + solfwv(iw)

            if ( iw == iendwv(ii) ) then
              solwaer(ib,ii) = sumsol

              if ( ii < NAERBND ) then
                sumsol = f_zero
                ii = ii + 1
              endif
            endif
          enddo

          if ( iw2 /= iendwv(ii) ) then
            solwaer(ib,ii) = sumsol
          endif

          nv2(ib) = ii
!         frcbnd(ib) = solbnd(ib) / soltot
        enddo     ! end do_ib_block for sw
      endif    ! end if_laswflg_block

!  --- ...  compute lw flux weights and interval indices for mapping
!           spectral bands between lw radiation and aerosol data

      if ( lalwflg ) then
        eirbnd (:)   = f_zero
        eirwaer(:,:) = f_zero

        do ib = 1, NLWBND
          ii = 1
          if ( NLWBND == 1 ) then
!           iw1 = 250                   ! corresponding 40 mu
            iw1 = 400                   ! corresponding 25 mu
            iw2 = 2500                  ! corresponding 4  mu
          else
            mb = ib + NLWSTR - 1
            iw1 = nint(wvnlw1(mb))
            iw2 = nint(wvnlw2(mb))
          endif

          Lab_lwdowhile : do while ( iw1 > iendwv(ii) )
            if ( ii == NAERBND ) exit Lab_lwdowhile
            ii = ii + 1
          enddo  Lab_lwdowhile

          sumir = f_zero
          nr1(ib) = ii

          do iw = iw1, iw2
            eirbnd(ib) = eirbnd(ib) + eirfwv(iw)
            sumir  = sumir  + eirfwv(iw)

            if ( iw == iendwv(ii) ) then
              eirwaer(ib,ii) = sumir

              if ( ii < NAERBND ) then
                sumir = f_zero
                ii = ii + 1
              endif
            endif
          enddo

          if ( iw2 /= iendwv(ii) ) then
            eirwaer(ib,ii) = sumir
          endif

          nr2(ib) = ii
        enddo     ! end do_ib_block for lw
      endif    ! end if_lalwflg_block

!  ---  compute spectral band mean properties for each species

      call optavg
!  ---  inputs:  (in-scope variables, module variables)
!  ---  outputs: (module variables)

!  ---  check print
!     do ib = 1, NSWBND
!       print *,' After optavg, for sw band:',ib
!       print *,'  extrhi:', extrhi(:,ib)
!       print *,'  scarhi:', scarhi(:,ib)
!       print *,'  ssarhi:', ssarhi(:,ib)
!       print *,'  asyrhi:', asyrhi(:,ib)
!       mb = ib + NSWSTR - 1
!       print *,'  wvnsw1,wvnsw2 :',wvnsw1(mb),wvnsw2(mb)
!       do i = 1, NRHLEV
!         print *,'  extrhd for rhlev:',i
!         print *,extrhd(i,:,ib)
!         print *,'  scarhd for rhlev:',i
!         print *,scarhd(i,:,ib)
!         print *,'  ssarhd for rhlev:',i
!         print *,ssarhd(i,:,ib)
!         print *,'  asyrhd for rhlev:',i
!         print *,asyrhd(i,:,ib)
!       enddo
!       print *,' extstra:', extstra(ib)
!     enddo
!     print *,'  wvnlw1 :',wvnlw1
!     print *,'  wvnlw2 :',wvnlw2
!     do ib = 1, NLWBND
!       ii = NSWBND + ib
!       print *,' After optavg, for lw band:',ib
!       print *,'  extrhi:', extrhi(:,ii)
!       print *,'  scarhi:', scarhi(:,ii)
!       print *,'  ssarhi:', ssarhi(:,ii)
!       print *,'  asyrhi:', asyrhi(:,ii)
!       do i = 1, NRHLEV
!         print *,'  extrhd for rhlev:',i
!         print *,extrhd(i,:,ii)
!         print *,'  scarhd for rhlev:',i
!         print *,scarhd(i,:,ii)
!         print *,'  ssarhd for rhlev:',i
!         print *,ssarhd(i,:,ii)
!         print *,'  asyrhd for rhlev:',i
!         print *,asyrhd(i,:,ii)
!       enddo
!       print *,' extstra:', extstra(ii)
!     enddo
!
      return
!................................
      end subroutine set_aercoef
!--------------------------------

!--------------------------------
      subroutine optavg
!................................
!  ---  inputs:  (in-scope variables, module variables
!  ---  outputs: (module variables)

! ==================================================================== !
!                                                                      !
! subprogram: optavg                                                   !
!                                                                      !
!   compute mean aerosols optical properties over each sw radiation    !
!   spectral band for each of the species components.  This program    !
!   follows gfdl's approach for thick cloud opertical property in      !
!   sw radiation scheme (2000).                                        !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
! major input variables:                                               !
!   nv1,nv2 (NSWBND) - start/end spectral band indices of aerosol data !
!                      for each sw radiation spectral band             !
!   nr1,nr2 (NLWBND) - start/end spectral band indices of aerosol data !
!                      for each ir radiation spectral band             !
!   solwaer (NSWBND,NAERBND)                                           !
!                    - solar flux weight over each sw radiation band   !
!                      vs each aerosol data spectral band              !
!   eirwaer (NLWBND,NAERBND)                                           !
!                    - ir flux weight over each lw radiation band      !
!                      vs each aerosol data spectral band              !
!   solbnd  (NSWBND) - solar flux weight over each sw radiation band   !
!   eirbnd  (NLWBND) - ir flux weight over each lw radiation band      !
!   NSWBND           - total number of sw spectral bands               !
!   NLWBND           - total number of lw spectral bands               !
!                                                                      !
! external module variables:  (in physpara)                            !
!   laswflg          - control flag for sw spectral region             !
!   lalwflg          - control flag for lw spectral region             !
!                                                                      !
! output variables: (to module variables)                              !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
!  ---  output:

!  ---  locals:
      real (kind=kind_phys) :: sumk, sums, sumok, sumokg, sumreft,      &
     &       sp, refb, reft, rsolbd, rirbd

      integer :: ib, nb, ni, nh, nc
!
!===> ...  begin here
!
!  --- ...  loop for each sw radiation spectral band

      if ( laswflg ) then

        do nb = 1, NSWBND
          rsolbd = f_one / solbnd(nb)

!  ---  for rh independent aerosol species

          do nc = 1, NCM1
            sumk    = f_zero
            sums    = f_zero
            sumok   = f_zero
            sumokg  = f_zero
            sumreft = f_zero

            do ni = nv1(nb), nv2(nb)
              sp   = sqrt( (f_one - rhidssa0(ni,nc))                    &
     &             / (f_one - rhidssa0(ni,nc)*rhidasy0(ni,nc)) )
              reft = (f_one - sp) / (f_one + sp)
              sumreft = sumreft + reft*solwaer(nb,ni)

              sumk    = sumk    + rhidext0(ni,nc)*solwaer(nb,ni)
              sums    = sums    + rhidsca0(ni,nc)*solwaer(nb,ni)
              sumok   = sumok   + rhidssa0(ni,nc)*solwaer(nb,ni)        &
     &                * rhidext0(ni,nc)
              sumokg  = sumokg  + rhidssa0(ni,nc)*solwaer(nb,ni)        &
     &                * rhidext0(ni,nc)*rhidasy0(ni,nc)
            enddo

            refb = sumreft * rsolbd

            extrhi(nc,nb) = sumk   * rsolbd
            scarhi(nc,nb) = sums   * rsolbd
            asyrhi(nc,nb) = sumokg / (sumok + 1.0e-10)
            ssarhi(nc,nb) = 4.0*refb                                    &
     &         / ( (f_one+refb)**2 - asyrhi(nc,nb)*(f_one-refb)**2 )
          enddo   ! end do_nc_block for rh-ind aeros

!  ---  for rh dependent aerosols species

          do nc = 1, NCM2
            do nh = 1, NRHLEV
              sumk    = f_zero
              sums    = f_zero
              sumok   = f_zero
              sumokg  = f_zero
              sumreft = f_zero

              do ni = nv1(nb), nv2(nb)
                sp   = sqrt( (f_one - rhdpssa0(ni,nh,nc))               &
     &               / (f_one - rhdpssa0(ni,nh,nc)*rhdpasy0(ni,nh,nc)) )
                reft = (f_one - sp) / (f_one + sp)
                sumreft = sumreft + reft*solwaer(nb,ni)

                sumk    = sumk    + rhdpext0(ni,nh,nc)*solwaer(nb,ni)
                sums    = sums    + rhdpsca0(ni,nh,nc)*solwaer(nb,ni)
                sumok   = sumok   + rhdpssa0(ni,nh,nc)*solwaer(nb,ni)   &
     &                  * rhdpext0(ni,nh,nc)
                sumokg  = sumokg  + rhdpssa0(ni,nh,nc)*solwaer(nb,ni)   &
     &                  * rhdpext0(ni,nh,nc)*rhdpasy0(ni,nh,nc)
              enddo

              refb = sumreft * rsolbd

              extrhd(nh,nc,nb) = sumk   * rsolbd
              scarhd(nh,nc,nb) = sums   * rsolbd
              asyrhd(nh,nc,nb) = sumokg / (sumok + 1.0e-10)
              ssarhd(nh,nc,nb) = 4.0*refb                               &
     &         / ( (f_one+refb)**2 - asyrhd(nh,nc,nb)*(f_one-refb)**2 )
            enddo   ! end do_nh_block
          enddo   ! end do_nc_block for rh-dep aeros

!  ---  for stratospheric background aerosols

          sumk = f_zero
          do ni = nv1(nb), nv2(nb)
            sumk = sumk + straext0(ni)*solwaer(nb,ni)
          enddo

          extstra(nb) = sumk * rsolbd

!  ---  check print
!         if ( nb > 6 .and. nb < 10) then
!           print *,' in optavg for sw band',nb
!           print *,'  nv1, nv2:',nv1(nb),nv2(nb)
!           print *,'  solwaer:',solwaer(nb,nv1(nb):nv2(nb))
!           print *,'  extrhi:', extrhi(:,nb)
!           do i = 1, NRHLEV
!             print *,'  extrhd for rhlev:',i
!             print *,extrhd(i,:,nb)
!           enddo
!           print *,'  sumk, rsolbd, extstra:',sumk,rsolbd,extstra(nb)
!         endif

        enddo   !  end do_nb_block for sw
      endif   !  end if_laswflg_block

!  --- ...  loop for each lw radiation spectral band

      if ( lalwflg ) then

        do nb = 1, NLWBND

          ib = NSWBND + nb
          rirbd = f_one / eirbnd(nb)

!  ---  for rh independent aerosol species

          do nc = 1, NCM1
            sumk    = f_zero
            sums    = f_zero
            sumok   = f_zero
            sumokg  = f_zero
            sumreft = f_zero

            do ni = nr1(nb), nr2(nb)
              sp   = sqrt( (f_one - rhidssa0(ni,nc))                    &
     &             / (f_one - rhidssa0(ni,nc)*rhidasy0(ni,nc)) )
              reft = (f_one - sp) / (f_one + sp)
              sumreft = sumreft + reft*eirwaer(nb,ni)

              sumk    = sumk    + rhidext0(ni,nc)*eirwaer(nb,ni)
              sums    = sums    + rhidsca0(ni,nc)*eirwaer(nb,ni)
              sumok   = sumok   + rhidssa0(ni,nc)*eirwaer(nb,ni)        &
     &                * rhidext0(ni,nc)
              sumokg  = sumokg  + rhidssa0(ni,nc)*eirwaer(nb,ni)        &
     &                * rhidext0(ni,nc)*rhidasy0(ni,nc)
            enddo

            refb = sumreft * rirbd

            extrhi(nc,ib) = sumk   * rirbd
            scarhi(nc,ib) = sums   * rirbd
            asyrhi(nc,ib) = sumokg / (sumok + 1.0e-10)
            ssarhi(nc,ib) = 4.0*refb                                       &
     &         / ( (f_one+refb)**2 - asyrhi(nc,ib)*(f_one-refb)**2 )
          enddo   ! end do_nc_block for rh-ind aeros

!  ---  for rh dependent aerosols species

          do nc = 1, NCM2
            do nh = 1, NRHLEV
              sumk    = f_zero
              sums    = f_zero
              sumok   = f_zero
              sumokg  = f_zero
              sumreft = f_zero

              do ni = nr1(nb), nr2(nb)
                sp   = sqrt( (f_one - rhdpssa0(ni,nh,nc))               &
     &             / (f_one - rhdpssa0(ni,nh,nc)*rhdpasy0(ni,nh,nc)) )
                reft = (f_one - sp) / (f_one + sp)
                sumreft = sumreft + reft*eirwaer(nb,ni)

                sumk    = sumk    + rhdpext0(ni,nh,nc)*eirwaer(nb,ni)
                sums    = sums    + rhdpsca0(ni,nh,nc)*eirwaer(nb,ni)
                sumok   = sumok   + rhdpssa0(ni,nh,nc)*eirwaer(nb,ni)   &
     &                  * rhdpext0(ni,nh,nc)
                sumokg  = sumokg  + rhdpssa0(ni,nh,nc)*eirwaer(nb,ni)   &
     &                  * rhdpext0(ni,nh,nc)*rhdpasy0(ni,nh,nc)
              enddo

              refb = sumreft * rirbd

              extrhd(nh,nc,ib) = sumk   * rirbd
              scarhd(nh,nc,ib) = sums   * rirbd
              asyrhd(nh,nc,ib) = sumokg / (sumok + 1.0e-10)
              ssarhd(nh,nc,ib) = 4.0*refb                               &
     &         / ( (f_one+refb)**2 - asyrhd(nh,nc,ib)*(f_one-refb)**2 )
            enddo   ! end do_nh_block
          enddo   ! end do_nc_block for rh-dep aeros

!  ---  for stratospheric background aerosols

          sumk = f_zero
          do ni = nr1(nb), nr2(nb)
            sumk = sumk + straext0(ni)*eirwaer(nb,ni)
          enddo

          extstra(ib) = sumk * rirbd

!  ---  check print
!         if ( nb >= 1 .and. nb < 5) then
!           print *,' in optavg for ir band:',nb
!           print *,'  nr1, nr2:',nr1(nb),nr2(nb)
!           print *,'  eirwaer:',eirwaer(nb,nr1(nb):nr2(nb))
!           print *,'  extrhi:', extrhi(:,ib)
!           do i = 1, NRHLEV
!             print *,'  extrhd for rhlev:',i
!             print *,extrhd(i,:,ib)
!           enddo
!           print *,'  sumk, rirbd, extstra:',sumk,rirbd,extstra(ib)
!         endif

        enddo   !  end do_nb_block for lw
      endif   !  end if_lalwflg_block
!
      return
!................................
      end subroutine optavg
!--------------------------------
!
!...................................
      end subroutine clim_aerinit
!-----------------------------------



!-----------------------------------
      subroutine aer_update                                             &
!...................................
!  ---  inputs:
     &     ( iyear, imon, me )
!  ---  outputs: ( to module variables )

!  ==================================================================  !
!                                                                      !
!  aer_update checks and update time varying climatology aerosol       !
!    data sets.                                                        !
!                                                                      !
!  inputs:                                          size               !
!     iyear   - 4-digit calender year                 1                !
!     imon    - month of the year                     1                !
!     me      - print message control flag            1                !
!                                                                      !
!  outputs: ( none )                                                   !
!                                                                      !
!  external module variables: (in physpara)                            !
!     lalwflg     - control flag for tropospheric lw aerosol           !
!     laswflg     - control flag for tropospheric sw aerosol           !
!     lavoflg     - control flag for stratospheric volcanic aerosol    !
!                                                                      !
!  usage:    call aero_update                                          !
!                                                                      !
!  subprograms called:  trop_update, volc_update                       !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      integer,  intent(in) :: iyear, imon, me

!  ---  output: ( none )

!  ---  locals: ( none )
!
!===> ...  begin here
!
      if ( imon < 1 .or. imon > 12 ) then
        print *,' ***** ERROR in specifying requested month !!! ',      &
     &          'imon=', imon
        print *,' ***** STOPPED in subroutinte aer_update !!!'
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif

      if ( lalwflg .or. laswflg ) then ! update monthly tropspheric aerosol data
        call trop_update
      endif

      if ( lavoflg ) then              ! update yearly stratospheric volcanic aerosol data
        call volc_update
      endif


! =================
      contains
! =================

!--------------------------------
      subroutine trop_update
!................................
!  ---  inputs:    (in scope variables, module variables)
!  ---  outputs:   (module variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : trop_update                                            !
!                                                                      !
!    updates the  monthly global distribution of aerosol profiles in   !
!    five degree horizontal resolution.                                !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  inputs:  (in-scope variables, module constants)                     !
!   imon     - integer, month of the year                              !
!   me       - integer, print message control flag                     !
!                                                                      !
!  outputs: (module variables)                                         !
!                                                                      !
!  external module variables: (in physpara)                            !
!    aeros_file   - external aerosol data file name                    !
!                                                                      !
!  internal module variables:                                          !
!    kprfg (    IMXAE*JMXAE)   - aeros profile index                   !
!    idxcg (NXC*IMXAE*JMXAE)   - aeros component index                 !
!    cmixg (NXC*IMXAE*JMXAE)   - aeros component mixing ratio          !
!    denng ( 2 *IMXAE*JMXAE)   - aerosols number density               !
!                                                                      !
!    NIAERCM      - unit number for input data set                     !
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call trop_update                                          !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: ( none )
!  ---  output: ( none )

!  ---  locals:
!     real (kind=kind_io8) :: cmix(NXC), denn, tem
      real (kind=kind_phys) :: cmix(NXC), denn, tem
      integer              :: idxc(NXC), kprf

      integer :: i, id, j, k, m, nc
      logical :: file_exist

      character :: cline*80, ctyp*3
!
!===>  ...  begin here
!
!  --- ...  reading climatological aerosols data

      inquire (file=aeros_file, exist=file_exist)

      if ( file_exist ) then
        close(NIAERCM)
        open (unit=NIAERCM,file=aeros_file,status='OLD',                &
     &        form='FORMATTED')
        rewind (NIAERCM)

        if ( me == 0 ) then
          print *,'   Opened aerosol data file: ',aeros_file
        endif
      else
        print *,'    Requested aerosol data file "',aeros_file,         &
     &          '" not found!'
        print *,'    *** Stopped in subroutine trop_update !!'
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif      ! end if_file_exist_block

      do j = 1, JMXAE
        do i = 1, IMXAE
          do m = 1, NXC
            idxcg(m,i,j) = 0
            cmixg(m,i,j) = f_zero
          enddo
        enddo
      enddo

      do j = 1, JMXAE
        do i = 1, IMXAE
          do m = 1, 2
            denng(m,i,j) = f_zero
          enddo
        enddo
      enddo

!  --- ...  loop over 12 month global distribution

      Lab_do_12mon : do m = 1, 12

        read(NIAERCM,12) cline
  12    format(a80/)

        if ( m /= imon ) then
!         if ( me == 0 ) print *,'  *** Skipped ',cline

          do j = 1, JMXAE
          do i = 1, IMXAE
            read(NIAERCM,*) id
          enddo
          enddo
        else
!         if ( me == 0 ) print *,'  --- Reading ',cline

          do j = 1, JMXAE
          do i = 1, IMXAE
            read(NIAERCM,14) (idxc(k),cmix(k),k=1,NXC),kprf,denn,nc,ctyp
  14        format(5(i2,e11.4),i2,f8.2,i3,1x,a3)

            kprfg(i,j)     = kprf
            denng(1,i,j)   = denn       ! num density of 1st layer
            if ( kprf >= 6 ) then
              denng(2,i,j) = cmix(NXC)  ! num density of 2dn layer
            else
              denng(2,i,j) = f_zero
            endif

            tem = f_one
            do k = 1, NXC-1
              idxcg(k,i,j) = idxc(k)    ! component index
              cmixg(k,i,j) = cmix(k)    ! component mixing ratio
              tem          = tem - cmix(k)
            enddo
            idxcg(NXC,i,j) = idxc(NXC)
            cmixg(NXC,i,j) = tem        ! to make sure all add to 1.
          enddo
          enddo

          close (NIAERCM)
          exit  Lab_do_12mon
        endif     ! end if_m_block

      enddo  Lab_do_12mon

!  --  check print

!     print *,'  IDXCG :'
!     print 16,idxcg
! 16  format(40i3)
!     print *,'  CMIXG :'
!     print 17,cmixg
!     print *,'  DENNG :'
!     print 17,denng
!     print *,'  KPRFG :'
!     print 17,kprfg
! 17  format(8e16.9)
!
      return
!................................
      end subroutine trop_update
!--------------------------------


!--------------------------------
      subroutine volc_update
!................................
!  ---  inputs:    (in scope variables, module variables)
!  ---  outputs:   (module variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : volc_update                                            !
!                                                                      !
!    searches historical volcanic data sets to find and read in        !
!    monthly 45-degree lat-zone band data of optical depth.            !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  inputs:  (in-scope variables, module constants)                     !
!   iyear    - integer, 4-digit calender year                 1        !
!   imon     - integer, month of the year                     1        !
!   me       - integer, print message control flag            1        !
!   NIAERCM  - integer, unit number for input data set        1        !
!                                                                      !
!  outputs: (module variables)                                         !
!   ivolae   - integer, monthly, 45-deg lat-zone volc odp      12*4*10 !
!   kyrstr   - integer, starting year of data in the input file        !
!   kyrend   - integer, ending   year of data in the input file        !
!   kyrsav   - integer, the year of data in use in the input file      !
!   kmonsav  - integer, the month of data in use in the input file     !
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call volc_aerinit                                         !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: (in-scope variables, module constants)
!     integer :: iyear, imon, me, NIAERCM

!  ---  output: (module variables)
!     integer :: ivolae(:,:,:), kyrstr, kyrend, kyrsav, kmonsav

!  ---  locals:
      integer :: i, j, k
      logical :: file_exist

      character :: cline*80, volcano_file*32
      data volcano_file / 'volcanic_aerosols_1850-1859.txt ' /
!
!===>  ...  begin here
!
      kmonsav = imon

      if ( kyrstr<=iyear .and. iyear<=kyrend ) then   ! use previously input data
        kyrsav = iyear
        return
      else                                            ! need to input new data
        kyrsav = iyear
        kyrstr = iyear - mod(iyear,10)
        kyrend = kyrstr + 9

!  ---  check print
!       print *,'  kyrstr, kyrend, kyrsav, kmonsav =',                  &
!    &          kyrstr,kyrend,kyrsav,kmonsav

        if ( iyear < MINVYR .or. iyear > MAXVYR ) then
          ivolae(:,:,:) = 1            ! set as lowest value
          if ( me == 0 ) then
            print *,'   Request volcanic date out of range,',           &
     &              ' optical depth set to lowest value'
          endif
        else
          write(volcano_file(19:27),60) kyrstr,kyrend
  60      format(i4.4,'-',i4.4)

          inquire (file=volcano_file, exist=file_exist)
          if ( file_exist ) then
            close(NIAERCM)
            open (unit=NIAERCM,file=volcano_file,status='OLD',          &
     &            form='FORMATTED')

            read(NIAERCM,62) cline
  62        format(a80)

!  ---  check print
            if ( me == 0 ) then
              print *,'   Opened volcanic data file: ',volcano_file
              print *, cline
            endif

            do k = 1, 10
              do j = 1, 4
                read(NIAERCM,64) (ivolae(i,j,k),i=1,12)
  64            format(12i5)
              enddo
            enddo

            close (NIAERCM)
          else
            print *,'   Requested volcanic data file "',                &
     &              volcano_file,'" not found!'
            print *,'   *** Stopped in subroutine VOLC_AERINIT !!'
            ! stop
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif   ! end if_file_exist_block

        endif   ! end if_iyear_block
      endif   ! end if_kyrstr_block

!  ---  check print
      if ( me == 0 ) then
        k = mod(kyrsav,10) + 1
        print *,' CHECK: Sample Volcanic data used for month, year:',   &
     &           imon, iyear
        print *,  ivolae(kmonsav,:,k)
      endif
!
      return
!................................
      end subroutine volc_update
!--------------------------------
!
!...................................
      end subroutine aer_update
!-----------------------------------



!-----------------------------------
      subroutine setaer                                                 &
!...................................

!  ---  inputs:
     &     ( prsi,prsl,prslk,tvly,rhlay,slmsk,tracer,xlon,xlat,         &
     &       IMAX,NLAY,NLP1, lsswr,lslwr,                               &
!  ---   extra 2 vars for old code that has been replaced by 'tvly'
!        can be eliminated after aggrements being reached
     &       tlay,qlay,                                                 & 
!  ---  outputs:
     &       aerosw,aerolw                                              &
!    &       aerosw,aerolw,aerodp                                       &
     &     )

!  ==================================================================  !
!                                                                      !
!  setaer computes aerosols optical properties                         !
!                                                                      !
!  inputs:                                                   size      !
!     prsi    - pressure at interface              mb      IMAX*NLP1   !
!     prsl    - layer mean pressure                mb      IMAX*NLAY   !
!     prslk   - exner function = (p/p0)**rocp              IMAX*NLAY   !
!     tvly    - layer virtual temperature          k       IMAX*NLAY   !
!     rhlay   - layer mean relative humidity               IMAX*NLAY   !
!     slmsk   - sea/land mask (sea:0,land:1,sea-ice:2)       IMAX      !
!     tracer  - aerosol tracer concentration           IMAX*NLAY*NTRAC !
!     xlon    - longitude of given points in radiance        IMAX      !
!               ok for both 0->2pi or -pi->+pi ranges                  !
!     xlat    - latitude of given points in radiance         IMAX      !
!               default to pi/2 -> -pi/2, otherwise see in-line comment!
!     IMAX    - horizontal dimension of arrays                  1      !
!     NLAY,NLP1-vertical dimensions of arrays                   1      !
!     lsswr,lslwr                                                      !
!             - logical flags for sw/lw radiation calls         1      !
!                                                                      !
!  outputs:                                                            !
!     aerosw - aeros opt properties for sw      IMAX*NLAY*NBDSW*NF_AESW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!     aerolw - aeros opt properties for lw      IMAX*NLAY*NBDLW*NF_AELW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!!    aerodp - vertically integrated optical depth         IMAX*NSPC1  !
!                                                                      !
!  external module variable: (in physpara)                             !
!     iaerflg - aerosol effect control flag (volc,lw,sw, 3-dig)        !
!     laswflg - tropospheric aerosol control flag for sw radiation     !
!               =f: no sw aeros calc.  =t: do sw aeros calc.           !
!     lalwflg - tropospheric aerosol control flag for lw radiation     !
!               =f: no lw aeros calc.  =t: do lw aeros calc.           !
!     lavoflg - control flag for stratospheric vocanic aerosols        !
!               =t: add volcanic aerosols to the background aerosols   !
!     ivflip  - control flag for direction of vertical index           !
!               =0: index from toa to surface                          !
!               =1: index from surface to toa                          !
!                                                                      !
!  internal module variable: (set by subroutine aer_init)              !
!     ivolae  - stratosphere volcanic aerosol optical depth (fac 1.e4) !
!                                                     12*4*10          !
!  usage:    call setaer                                               !
!                                                                      !
!  subprograms called:  aer_property                                   !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      integer, intent(in) :: IMAX, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: prsi, prsl,  &
     &       prslk, tvly, rhlay, tlay, qlay
      real (kind=kind_phys), dimension(:),   intent(in) :: xlon, xlat,  &
     &       slmsk
      real (kind=kind_phys), dimension(:,:,:),intent(in):: tracer

      logical, intent(in) :: lsswr, lslwr

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:,:), intent(out) ::         &
     &       aerosw, aerolw
!     real (kind=kind_phys), dimension(:,:)    , intent(out) :: aerodp

!  ---  locals:
      real (kind=kind_phys), parameter :: psrfh = 5.0    ! ref press (mb) for upper bound

      real (kind=kind_phys), dimension(IMAX) :: alon,alat,volcae,rdelp
!     real (kind=kind_phys), dimension(IMAX) :: sumodp
      real (kind=kind_phys) :: prsln(NLP1),hz(IMAX,NLP1),dz(IMAX,NLAY)
      real (kind=kind_phys) :: tmp1, tmp2, psrfl

      integer               :: kcutl(IMAX), kcuth(IMAX)
      integer               :: i, i1, j, k, m, mb, kh, kl

      logical               :: laddsw=.false.,  laersw=.false.
      logical               :: laddlw=.false.,  laerlw=.false.
!$omp threadprivate(laddsw,laersw,laddlw,laerlw)
!  ---  conversion constants
      real (kind=kind_phys), parameter :: rdg  = 180.0 / con_pi
      real (kind=kind_phys), parameter :: rovg = 0.001 * con_rd / con_g

!===>  ...  begin here

      do m = 1, NF_AESW
        do j = 1, NBDSW
          do k = 1, NLAY
            do i = 1, IMAX
              aerosw(i,k,j,m) = f_zero
            enddo
          enddo
        enddo
      enddo

      do m = 1, NF_AELW
        do j = 1, NBDLW
          do k = 1, NLAY
            do i = 1, IMAX
              aerolw(i,k,j,m) = f_zero
            enddo
          enddo
        enddo
      enddo

!     aerodp = f_zero
!     sumodp = f_zero

      if ( .not. (lsswr .or. lslwr) ) then
        return
      endif

      if ( iaerflg <= 0 ) then
        return
      endif

      laersw = lsswr .and. laswflg
      laerlw = lslwr .and. lalwflg

!  ---  ...  convert lat/lon from radiance to degree

      do i = 1, IMAX
        alon(i) = xlon(i) * rdg
        if (alon(i) < f_zero) alon(i) = alon(i) + 360.0
        alat(i) = xlat(i) * rdg          ! if xlat in pi/2 -> -pi/2 range
        alat(i)=min(max(alat(i),-90.0),90.0)
!       alat(i) = 90.0 - xlat(i)*rdg     ! if xlat in 0 -> pi range
      enddo

!  ---  ...  compute level height and layer thickness

      if ( laswflg .or. lalwflg ) then

        lab_do_IMAX : do i = 1, IMAX

          lab_if_flip : if (ivflip == 1) then       ! input from sfc to toa

            do k = 1, NLAY
              prsln(k) = log(prsi(i,k))
            enddo
            prsln(NLP1)= log(prsl(i,NLAY))

            do k = NLAY, 1, -1
              dz(i,k) = rovg * (prsln(k) - prsln(k+1))                  &
     &                * tlay(i,k) * (f_one + con_fvirt*qlay(i,k))

              !  dz(i,k) = rovg * (prsln(k) - prsln(k+1)) * tvly(i,k)
            enddo
            dz(i,NLAY)  = 2.0 * dz(i,NLAY)

            hz(i,1) = f_zero
            do k = 1, NLAY
              hz(i,k+1) = hz(i,k) + dz(i,k)
            enddo

          else  lab_if_flip                         ! input from toa to sfc

            prsln(1) = log(prsl(i,1))
            do k = 2, NLP1
              prsln(k) = log(prsi(i,k))
            enddo

            do k = 1, NLAY
              dz(i,k) = rovg * (prsln(k+1) - prsln(k))                  &
     &                * tlay(i,k) * (f_one + con_fvirt*qlay(i,k))

              !  dz(i,k) = rovg * (prsln(k+1) - prsln(k)) * tvly(i,k)
            enddo
            dz(i,1) = 2.0 * dz(i,1)

            hz(i,NLP1) = f_zero
            do k = NLAY, 1, -1
              hz(i,k) = hz(i,k+1) + dz(i,k)
            enddo

          endif  lab_if_flip

        enddo  lab_do_IMAX

!  ---  ...  calculate sw aerosol optical properties for the corresponding
!            frequency bands

        call aer_property                                               &
!  ---  inputs:
     &     ( prsi,prsl,prslk,tvly,rhlay,dz,hz,tracer,                   &
     &       alon,alat,slmsk, laersw,laerlw,                            &
     &       IMAX,NLAY,NLP1,                                            &
!    &       IMAX,NLAY,NLP1,NSPC1,                                      &
!  ---  outputs:
     &       aerosw,aerolw                                              &
!    &       aerosw,aerolw,aerodp                                       &
     &     )

!  ---  check print
!       do m = 1, NBDSW
!         print *,'  ***  CHECK AEROSOLS PROPERTIES FOR SW BAND =',m,   &
!    &            ' ***'
!         do k = 1, 10
!           print *,'  LEVEL :',k
!           print *,'  TAUAER:',aerosw(:,k,m,1)
!           print *,'  SSAAER:',aerosw(:,k,m,2)
!           print *,'  ASYAER:',aerosw(:,k,m,3)
!         enddo
!       enddo
!       print *,'  ***  CHECK AEROSOLS OPTICAL DEPTH FOR 550nm REGION'
!       print *, aerodp(:,1)
!       if ( laod_out ) then
!         do m = 1, NSPC1
!           print *,'  ***  CHECK AEROSOLS OPTICAL DEPTH FOR SPECIES:', &
!    &              m
!           print *, aerodp(:,m)
!           sumodp(:) = sumodp(:) + aerodp(:,m)
!         enddo
!
!         print *,'  ***  CHECK AEROSOLS OPTICAL DEPTH FOR ALL SPECIES:'
!         print *, sumodp(:)
!       endif
!       do m = 1, NBDLW
!         print *,'  ***  CHECK AEROSOLS PROPERTIES FOR LW BAND =',m,   &
!    &            ' ***'
!         do k = 1, 10
!           print *,'  LEVEL :',k
!           print *,'  TAUAER:',aerolw(:,k,m,1)
!           print *,'  SSAAER:',aerolw(:,k,m,2)
!           print *,'  ASYAER:',aerolw(:,k,m,3)
!         enddo
!       enddo

      endif   ! end if_laswflg_or_lalwflg_block

!  ---  ...  stratosphere volcanic forcing

      if ( lavoflg ) then

        if ( iaerflg == 100 ) then
          laddsw = lsswr
          laddlw = lslwr
        else
          laddsw = lsswr .and. laswflg
          laddlw = lslwr .and. lalwflg
        endif

        i1 = mod(kyrsav, 10) + 1

!  ---  select data in 4 lat bands, interpolation at the boundaires

        do i = 1, IMAX
          if      ( alat(i) > 46.0 ) then
            volcae(i) = 1.0e-4 * ivolae(kmonsav,1,i1)
          else if ( alat(i) > 44.0 ) then
            volcae(i) = 5.0e-5                                          &
     &                * (ivolae(kmonsav,1,i1) + ivolae(kmonsav,2,i1))
          else if ( alat(i) >  1.0 ) then
            volcae(i) = 1.0e-4 * ivolae(kmonsav,2,i1)
          else if ( alat(i) > -1.0 ) then
            volcae(i) = 5.0e-5                                          &
     &                * (ivolae(kmonsav,2,i1) + ivolae(kmonsav,3,i1))
          else if ( alat(i) >-44.0 ) then
            volcae(i) = 1.0e-4 * ivolae(kmonsav,3,i1)
          else if ( alat(i) >-46.0 ) then
            volcae(i) = 5.0e-5                                          &
     &                * (ivolae(kmonsav,3,i1) + ivolae(kmonsav,4,i1))
          else
            volcae(i) = 1.0e-4 * ivolae(kmonsav,4,i1)
          endif
        enddo

        if ( ivflip == 0 ) then         ! input data from toa to sfc

!  ---  find lower boundary of stratosphere

          do i = 1, IMAX

            tmp1 = abs( alat(i) )
            if ( tmp1 > 70.0 ) then          ! polar, fixed at 25000pa (250mb)
              psrfl = 250.0
            elseif ( tmp1 < 20.0 ) then      ! tropic, fixed at 15000pa (150mb)
              psrfl = 150.0
            else                             ! mid-lat, interpolation
              psrfl = 110.0 + 2.0*tmp1
            endif

            kcuth(i) = NLAY - 1
            kcutl(i) = 2
            rdelp(i) = f_one / prsi(i,2)

            lab_do_kcuth0 : do k = 2, NLAY-2
              if ( prsi(i,k) >= psrfh ) then
                kcuth(i) = k - 1
                exit lab_do_kcuth0
              endif
            enddo  lab_do_kcuth0

            lab_do_kcutl0 : do k = 2, NLAY-2
              if ( prsi(i,k) >= psrfl ) then
                kcutl(i) = k - 1
                rdelp(i) = f_one / (prsi(i,k) - prsi(i,kcuth(i)))
                exit lab_do_kcutl0
              endif
            enddo  lab_do_kcutl0
          enddo

!  ---  sw: add volcanic aerosol optical depth to the background value

          if ( laddsw ) then
            do m = 1, NBDSW
              mb = NSWSTR + m - 1

              if     ( wvnsw1(mb) > 20000 ) then   ! range of wvlth < 0.5mu
                tmp2 = 0.74
              elseif ( wvnsw2(mb) < 20000 ) then   ! range of wvlth > 0.5mu
                tmp2 = 1.14
              else                                 ! range of wvlth in btwn
                tmp2 = 0.94
              endif
              tmp1 = (0.275e-4 * (wvnsw2(mb)+wvnsw1(mb))) ** tmp2

              do i = 1, IMAX
                kh = kcuth(i)
                kl = kcutl(i)
                do k = kh, kl
                  tmp2 = tmp1 * ((prsi(i,k+1) - prsi(i,k)) * rdelp(i))
                  aerosw(i,k,m,1) = aerosw(i,k,m,1) + tmp2*volcae(i)
                enddo

!  ---  smoothing profile at boundary if needed

                if ( aerosw(i,kl,m,1) > 10.*aerosw(i,kl+1,m,1) ) then
                  tmp2 = aerosw(i,kl,m,1) + aerosw(i,kl+1,m,1)
                  aerosw(i,kl  ,m,1) = 0.8 * tmp2
                  aerosw(i,kl+1,m,1) = 0.2 * tmp2
                endif
              enddo    ! end do_i_block
            enddo      ! end do_m_block

!  ---  check print
!           do i = 1, IMAX
!             print *,' LEV  PRESS      TAU      FOR PROFILE:',i,       &
!    &                '  KCUTH, KCUTL =',kcuth(i),kcutl(i)
!             kh = kcuth(i) - 1
!             kl = kcutl(i) + 10
!             do k = kh, kl
!               write(6,71) k, prsl(i,k), aerosw(i,k,1,1)
! 71            format(i3,2e11.4)
!             enddo
!           enddo

          endif        ! end if_laddsw_block

!  ---  lw: add volcanic aerosol optical depth to the background value

          if ( laddlw ) then
            if ( NLWBND == 1 ) then

              tmp1 = (0.55 / 11.0) ** 1.2
              do i = 1, IMAX
                kh = kcuth(i)
                kl = kcutl(i)
                do k = kh, kl
                  tmp2 = tmp1 * ((prsi(i,k+1) - prsi(i,k)) * rdelp(i))  &
     &                 * volcae(i)
                  do m = 1, NBDLW
                    aerolw(i,k,m,1) = aerolw(i,k,m,1) + tmp2
                  enddo
                enddo
              enddo    ! end do_i_block

            else

              do m = 1, NBDLW
                tmp1 = (0.275e-4 * (wvnlw2(m) + wvnlw1(m))) ** 1.2

                do i = 1, IMAX
                  kh = kcuth(i)
                  kl = kcutl(i)
                  do k = kh, kl
                    tmp2 = tmp1 * ((prsi(i,k+1)-prsi(i,k)) * rdelp(i))
                    aerolw(i,k,m,1) = aerolw(i,k,m,1) + tmp2*volcae(i)
                  enddo
                enddo    ! end do_i_block
              enddo      ! end do_m_block

            endif      ! end if_NLWBND_block
          endif        ! end if_laddlw_block

        else                            ! input data from sfc to toa

!  ---  find lower boundary of stratosphere

          do i = 1, IMAX

            tmp1 = abs( alat(i) )
            if ( tmp1 > 70.0 ) then          ! polar, fixed at 25000pa (250mb)
              psrfl = 250.0
            elseif ( tmp1 < 20.0 ) then      ! tropic, fixed at 15000pa (150mb)
              psrfl = 150.0
            else                             ! mid-lat, interpolation
              psrfl = 110.0 + 2.0*tmp1
            endif

            kcuth(i) = 2
            kcutl(i) = NLAY - 1
            rdelp(i) = f_one / prsi(i,NLAY-1)

            lab_do_kcuth1 : do k = NLAY-1, 2, -1
              if ( prsi(i,k) >= psrfh ) then
                kcuth(i) = k
                exit lab_do_kcuth1
              endif
            enddo  lab_do_kcuth1

            lab_do_kcutl1 : do k = NLAY, 2, -1
              if ( prsi(i,k) >= psrfl ) then
                kcutl(i) = k
                rdelp(i) = f_one / (prsi(i,k) - prsi(i,kcuth(i)+1))
                exit lab_do_kcutl1
              endif
            enddo  lab_do_kcutl1
          enddo

!  ---  sw: add volcanic aerosol optical depth to the background value

          if ( laddsw ) then
            do m = 1, NBDSW
              mb = NSWSTR + m - 1

              if     ( wvnsw1(mb) > 20000 ) then   ! range of wvlth < 0.5mu
                tmp2 = 0.74
              elseif ( wvnsw2(mb) < 20000 ) then   ! range of wvlth > 0.5mu
                tmp2 = 1.14
              else                                 ! range of wvlth in btwn
                tmp2 = 0.94
              endif
              tmp1 = (0.275e-4 * (wvnsw2(mb)+wvnsw1(mb))) ** tmp2

              do i = 1, IMAX
                kh = kcuth(i)
                kl = kcutl(i)
                do k = kl, kh
                  tmp2 = tmp1 * ((prsi(i,k) - prsi(i,k+1)) * rdelp(i))
                  aerosw(i,k,m,1) = aerosw(i,k,m,1) + tmp2*volcae(i)
                enddo

!  ---  smoothing profile at boundary if needed

                if ( aerosw(i,kl,m,1) > 10.*aerosw(i,kl-1,m,1) ) then
                  tmp2 = aerosw(i,kl,m,1) + aerosw(i,kl-1,m,1)
                  aerosw(i,kl  ,m,1) = 0.8 * tmp2
                  aerosw(i,kl-1,m,1) = 0.2 * tmp2
                endif
              enddo    ! end do_i_block
            enddo      ! end do_m_block

!  ---  check print
!           do i = 1, IMAX
!             print *,' LEV  PRESS      TAU      FOR PROFILE:',i,       &
!    &                '  KCUTH, KCUTL =',kcuth(i),kcutl(i)
!             kh = kcuth(i) + 1
!             kl = kcutl(i) - 10
!             do k = kh, kl, -1
!               write(6,71) NLP1-k,prsl(i,k),aerosw(i,k,1,1)
!             enddo
!           enddo

          endif        ! end if_laddsw_block

!  ---  lw: add volcanic aerosol optical depth to the background value

          if ( laddlw ) then
            if ( NLWBND == 1 ) then

              tmp1 = (0.55 / 11.0) ** 1.2
              do i = 1, IMAX
                kh = kcuth(i)
                kl = kcutl(i)
                do k = kl, kh
                  tmp2 = tmp1 * ((prsi(i,k) - prsi(i,k+1)) * rdelp(i))  &
     &                 * volcae(i)
                  do m = 1, NBDLW
                    aerolw(i,k,m,1) = aerolw(i,k,m,1) + tmp2
                  enddo
                enddo
              enddo    ! end do_i_block

            else

              do m = 1, NBDLW
                tmp1 = (0.275e-4 * (wvnlw2(m) + wvnlw1(m))) ** 1.2

                do i = 1, IMAX
                  kh = kcuth(i)
                  kl = kcutl(i)
                  do k = kl, kh
                    tmp2 = tmp1 * ((prsi(i,k)-prsi(i,k+1)) * rdelp(i))
                    aerolw(i,k,m,1) = aerolw(i,k,m,1) + tmp2*volcae(i)
                  enddo
                enddo    ! end do_i_block
              enddo      ! end do_m_block

            endif      ! end if_NLWBND_block
          endif        ! end if_laddlw_block

        endif                           ! end if_ivflip_block

      endif   ! end if_lavoflg_block
!
      return
!...................................
      end subroutine setaer
!-----------------------------------



!-----------------------------------
      subroutine aer_property                                           &
!...................................

!  ---  inputs:
     &     ( prsi,prsl,prslk,tvly,rhlay,dz,hz,tracer,                   &
     &       alon,alat,slmsk, laersw,laerlw,                            &
     &       IMAX,NLAY,NLP1,                                            &
!    &       IMAX,NLAY,NLP1,NSPC,                                       &
!  ---  outputs:
     &       aerosw,aerolw                                              &
!    &       aerosw,aerolw,aerodp                                       &
     &     )

!  ==================================================================  !
!                                                                      !
!  aer_property maps the 5 degree global climatological aerosol data   !
!  set onto model grids, and compute aerosol optical properties for sw !
!  and lw radiations.                                                  !
!                                                                      !
!  inputs:                                                             !
!     prsi    - pressure at interface              mb      IMAX*NLP1   !
!     prsl    - layer mean pressure         (not used)     IMAX*NLAY   !
!     prslk   - exner function=(p/p0)**rocp (not used)     IMAX*NLAY   !
!     tvly    - layer virtual temperature   (not used)     IMAX*NLAY   !
!     rhlay   - layer mean relative humidity               IMAX*NLAY   !
!     dz      - layer thickness                    m       IMAX*NLAY   !
!     hz      - level high                         m       IMAX*NLP1   !
!     tracer  - aer tracer concentrations   (not used)  IMAX*NLAY*NTRAC!
!     alon, alat                                             IMAX      !
!             - longitude and latitude of given points in degree       !
!     slmsk   - sea/land mask (sea:0,land:1,sea-ice:2)       IMAX      !
!     laersw,laerlw                                             1      !
!             - logical flag for sw/lw aerosol calculations            !
!     IMAX    - horizontal dimension of arrays                  1      !
!     NLAY,NLP1-vertical dimensions of arrays                   1      !
!!    NSPC    - num of species for optional aod output fields   1      !
!                                                                      !
!  outputs:                                                            !
!     aerosw - aeros opt properties for sw      IMAX*NLAY*NBDSW*NF_AESW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!     aerolw - aeros opt properties for lw      IMAX*NLAY*NBDLW*NF_AELW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!!    aerodp - vertically integrated aer-opt-depth         IMAX*NSPC+1 !
!                                                                      !
!  module parameters and constants:                                    !
!     NSWBND  - total number of actual sw spectral bands computed      !
!     NLWBND  - total number of actual lw spectral bands computed      !
!     NSWLWBD - total number of sw+lw bands computed                   !
!                                                                      !
!  external module variables: (in physpara)                            !
!     ivflip  - control flag for direction of vertical index           !
!               =0: index from toa to surface                          !
!               =1: index from surface to toa                          !
!                                                                      !
!  module variable: (set by subroutine aer_init)                       !
!     kprfg   - aerosols profile index                IMXAE*JMXAE      !
!               1:ant  2:arc  3:cnt  4:mar  5:des  6:marme 7:cntme     !
!     idxcg   - aerosols component index              NXC*IMXAE*JMXAE  !
!               1:inso    2:soot    3:minm    4:miam    5:micm         !
!               6:mitr    7:waso    8:ssam    9:sscm   10:suso         !
!     cmixg   - aerosols component mixing ratio       NXC*IMXAE*JMXAE  !
!     denng   - aerosols number density                2 *IMXAE*JMXAE  !
!               1:for domain-1   2:domain-2 (prof marme/cntme only)    !
!                                                                      !
!  usage:    call aer_property                                         !
!                                                                      !
!  subprograms called:  radclimaer                                     !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      integer, intent(in) :: IMAX, NLAY, NLP1
!     integer, intent(in) :: IMAX, NLAY, NLP1, NSPC
      logical, intent(in) :: laersw, laerlw

      real (kind=kind_phys), dimension(:,:), intent(in) :: prsi, prsl,  &
     &       prslk, tvly, rhlay, dz, hz
      real (kind=kind_phys), dimension(:),   intent(in) :: alon, alat,  &
     &       slmsk
      real (kind=kind_phys), dimension(:,:,:),intent(in):: tracer

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:,:), intent(out) ::         &
     &       aerosw, aerolw
!     real (kind=kind_phys), dimension(:,:)    , intent(out) :: aerodp

!  ---  locals:
      real (kind=kind_phys), dimension(NCM) :: cmix
      real (kind=kind_phys), dimension(  2) :: denn
!     real (kind=kind_phys), dimension(NSPC) :: spcodp

      real (kind=kind_phys), dimension(NLAY) :: delz, rh1, dz1
      integer,               dimension(NLAY) :: idmaer

      real (kind=kind_phys), dimension(NLAY,NSWLWBD):: tauae,ssaae,asyae
!test real (kind=kind_phys), dimension(IMAX,NLAY) :: aersav

      real (kind=kind_phys) :: tmp1, tmp2, rps, dtmp, h1
      real (kind=kind_phys) :: wi, wj, w11, w12, w21, w22

      integer :: i, ii, i1, i2, i3,  j1, j2, j3,  k, m, m1,             &
     &           kp, kpa, kpi, kpj

!  ---  conversion constants
      real (kind=kind_phys), parameter :: dltg = 360.0 / float(IMXAE)
      real (kind=kind_phys), parameter :: hdlt = 0.5 * dltg
      real (kind=kind_phys), parameter :: rdlt = 1.0 / dltg

!
!===>  ...  begin here
!
!  ---  map aerosol data to model grids

      i1 = 1
      i2 = 2
      j1 = 1
      j2 = 2

      lab_do_IMAX : do i = 1, IMAX

!  ---  map grid in longitude direction, lon from 0 to 355 deg resolution

!       print *,' Seeking lon index for point i =',i
        i3 = i1
        lab_do_IMXAE : do while ( i3 <= IMXAE )
          tmp1 = dltg * (i3 - 1)
          dtmp = alon(i) - tmp1
!         print *,'   alon, i3, tlon, dlon =',alon(i),i3,tmp1,dtmp

          if ( dtmp > dltg ) then
            i3 = i3 + 1
            if ( i3 > IMXAE ) then
              print *,' ERROR! In setclimaer alon>360. ipt =',i,        &
     &           ',  dltg,alon,tlon,dlon =',dltg,alon(i),tmp1,dtmp
              ! stop
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            endif
          elseif ( dtmp >= f_zero ) then
            i1 = i3
            i2 = mod(i3,IMXAE) + 1
            wi = dtmp * rdlt
            if ( dtmp <= hdlt ) then
              kpi = i3
            else
              kpi = i2
            endif
!           print *,'   found i1, i2, wi =',i1,i2,wi
            exit lab_do_IMXAE
          else
            i3 = i3 - 1
            if ( i3 < 1 ) then
              print *,' ERROR! In setclimaer alon< 0. ipt =',i,         &
     &           ',  dltg,alon,tlon,dlon =',dltg,alon(i),tmp1,dtmp
              ! stop
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            endif
          endif
        enddo  lab_do_IMXAE

!  ---  map grid in latitude direction, lat from 90n to 90s in 5 deg resolution

!       print *,' Seeking lat index for point i =',i
        j3 = j1
        lab_do_JMXAE : do while ( j3 <= JMXAE )
          tmp2 = 90.0 - dltg * (j3 - 1)
          dtmp = tmp2 - alat(i)
!         print *,'   alat, j3, tlat, dlat =',alat(i),j3,tmp2,dtmp

          if ( dtmp > dltg ) then
            j3 = j3 + 1
            if ( j3 >= JMXAE ) then
              print *,' ERROR! In setclimaer alat<-90. ipt =',i,        &
     &           ',  dltg,alat,tlat,dlat =',dltg,alat(i),tmp2,dtmp
              ! stop
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            endif
          elseif ( dtmp >= f_zero ) then
            j1 = j3
            j2 = j3 + 1
            wj = dtmp * rdlt
            if ( dtmp <= hdlt ) then
              kpj = j3
            else
              kpj = j2
            endif
!           print *,'   found j1, j2, wj =',j1,j2,wj
            exit lab_do_JMXAE
          else
            j3 = j3 - 1
            if ( j3 < 1 ) then
              print *,' ERROR! In setclimaer alat>90. ipt =',i,         &
     &           ',  dltg,alat,tlat,dlat =',dltg,alat(i),tmp2,dtmp
              ! stop
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            endif
          endif
        enddo  lab_do_JMXAE

!  ---  determin the type of aerosol profile (kp) and scale hight for domain 1 (h1)
!       to be used at this grid point

        kp = kprfg(kpi,kpj)                     ! nearest typical aeros profile as default
        kpa = max( kprfg(i1,j1),kprfg(i1,j2),kprfg(i2,j1),kprfg(i2,j2) )
        h1 = haer(1,kp)
        denn(2) = f_zero
        ii = 1

        if ( kp /= kpa ) then
          if ( kpa == 6 ) then                  ! if ocean prof with mineral aeros overlay
            ii = 2                              ! need 2 types of densities
            if ( slmsk(i) > f_zero ) then       ! but actually a land/sea-ice point
              kp = 7                            ! reset prof index to land
              h1 = 0.5*(haer(1,6) + haer(1,7))  ! use a transition scale hight
            else
              kp = kpa
              h1 = haer(1,6)
            endif
          elseif ( kpa == 7 ) then              ! if land prof with mineral aeros overlay
            ii = 2                              ! need 2 types of densities
            if ( slmsk(i) <= f_zero ) then      ! but actually an ocean point
              kp = 6                            ! reset prof index to ocean
              h1 = 0.5*(haer(1,6) + haer(1,7))  ! use a transition scale hight
            else
              kp = kpa
              h1 = haer(1,7)
            endif
          else                                  ! lower atmos without mineral aeros overlay
!           h1 = 0.5*(haer(1,kp) + haer(1,kpa)) ! use a transition scale hight
            h1 = haer(1,kpa)
            kp = kpa
          endif
        endif

!  ---  compute horizontal bi-linear interpolation weights

        w11 = (f_one-wi) * (f_one-wj)
        w12 = (f_one-wi) *       wj
        w21 =        wi  * (f_one-wj)
        w22 =        wi  * wj

!  ---  check print
!       print *,'  Grid pt', i,',   alon, alat =',alon(i),alat(i),      &
!    &                       ',   tlon, tlat =',tmp1,tmp2
!       print *,'   lon grid index i1, i2 =',i1,i2,',  weight wi =',wi
!       print *,'   lat grid index j1, j2 =',j1,j2,',  weight wj =',wj
!       print *,'   bi-linear weights w11,w21,w12,w22 =',w11,w21,w12,w22
!       print *,'   kp,kpa,slmsk,h1 =',kp,m1,slmsk(i),h1

!  ---  do horizontal bi-linear interpolation on aerosol partical density (denn)

        do m = 1, ii                            ! ii=1 for domain 1; =2 for domain 2.
          denn(m) = w11*denng(m,i1,j1) + w12*denng(m,i1,j2)             &
     &            + w21*denng(m,i2,j1) + w22*denng(m,i2,j2)
        enddo  ! end_do_m_loop

!  ---  do horizontal bi-linear interpolation on mixing ratios

        cmix(:) = f_zero
        do m = 1, NXC
          ii = idxcg(m,i1,j1)
          if ( ii > 0 ) then
            cmix(ii) = cmix(ii) + w11*cmixg(m,i1,j1)
          endif
          ii = idxcg(m,i1,j2)
          if ( ii > 0 ) then
            cmix(ii) = cmix(ii) + w12*cmixg(m,i1,j2)
          endif
          ii = idxcg(m,i2,j1)
          if ( ii > 0 ) then
            cmix(ii) = cmix(ii) + w21*cmixg(m,i2,j1)
          endif
          ii = idxcg(m,i2,j2)
          if ( ii > 0 ) then
            cmix(ii) = cmix(ii) + w22*cmixg(m,i2,j2)
          endif
        enddo  ! end_do_m_loop

!  ---  check print
!       print *,'   denn =',denn(:)
!       print *,'   cmix =',cmix(:)

!  ---  prepare to setup domain index array and effective layer thickness
!       also convert pressure level to sigma level to follow the terrain

        do k = 1, NLAY
          rh1(k) = rhlay(i,k)
          dz1(k) = dz   (i,k)
        enddo

        lab_if_flip : if (ivflip == 1) then       ! input from sfc to toa

          if ( prsi(i,1) > 100.0 ) then
            rps = f_one / prsi(i,1)
          else
            print *,' !!! Error in subr radiation_aerosols:',           &
     &              ' unrealistic surface pressure =', prsi(i,1)
            ! stop
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif

          ii = 1
          do k = 1, NLAY
            if (prsi(i,k+1)*rps < sigref(ii,kp)) then
              ii = ii + 1
              if (ii == 2 .and. prsref(2,kp) == prsref(3,kp)) then
                ii = 3
              endif
            endif
            idmaer(k) = ii

            if ( ii > 1 ) then
              tmp1 = haer(ii,kp)
            else
              tmp1 = h1
            endif

            if (tmp1 > f_zero) then
              tmp2 = f_one / tmp1
              delz(k) = tmp1 * (exp(-hz(i,k)*tmp2)-exp(-hz(i,k+1)*tmp2))
            else
              delz(k) = dz1(k)
            endif
          enddo

        else  lab_if_flip                         ! input from toa to sfc

          if ( prsi(i,NLP1) > 100.0 ) then
            rps =  1.0 / prsi(i,NLP1)
          else
            print *,' !!! Error in subr radiation_aerosols:',           &
     &              ' unrealistic surface pressure =', prsi(i,NLP1)
            ! stop
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif

          ii = 1
          do k = NLAY, 1, -1
            if (prsi(i,k)*rps < sigref(ii,kp)) then
              ii = ii + 1
              if (ii == 2 .and. prsref(2,kp) == prsref(3,kp)) then
                ii = 3
              endif
            endif
            idmaer(k) = ii

            if ( ii > 1 ) then
              tmp1 = haer(ii,kp)
            else
              tmp1 = h1
            endif

            if (tmp1 > f_zero) then
              tmp2   = f_one / tmp1
              delz(k) = tmp1 * (exp(-hz(i,k+1)*tmp2)-exp(-hz(i,k)*tmp2))
            else
              delz(k) = dz1(k)
            endif
          enddo

        endif  lab_if_flip

!  ---  check print

!       print *,' in setclimaer, profile:',i
!       print *,'  rh   :',rh1
!       print *,'  dz   :',dz1
!       print *,'  delz :',delz
!       print *,'  idmaer:',idmaer

!  ---  calculate sw/lw aerosol optical properties for the
!       corresponding frequency bands

        call radclimaer
!  ---  inputs:  (in-scope variables)
!  ---  outputs: (in-scope variables)

        if ( laersw ) then

          do m = 1, NBDSW
            do k = 1, NLAY
              aerosw(i,k,m,1) = tauae(k,m)
              aerosw(i,k,m,2) = ssaae(k,m)
              aerosw(i,k,m,3) = asyae(k,m)
            enddo
          enddo

!  ---  total aod (optional)
!         do k = 1, NLAY
!           aerodp(i,1) = aerodp(i,1) + tauae(k,nv_aod)
!         enddo

!  ---  for diagnostic output (optional)
!         if ( lspcaod ) then
!           do m = 1, NSPC
!             aerodp(i,m+1) = spcodp(m)
!           enddo
!         endif

        endif     ! end if_larsw_block

        if ( laerlw ) then

          if ( NLWBND == 1 ) then
            m1 = NSWBND + 1
            do m = 1, NBDLW
              do k = 1, NLAY
                aerolw(i,k,m,1) = tauae(k,m1)
                aerolw(i,k,m,2) = ssaae(k,m1)
                aerolw(i,k,m,3) = asyae(k,m1)
              enddo
            enddo
          else
            do m = 1, NBDLW
              m1 = NSWBND + m
              do k = 1, NLAY
                aerolw(i,k,m,1) = tauae(k,m1)
                aerolw(i,k,m,2) = ssaae(k,m1)
                aerolw(i,k,m,3) = asyae(k,m1)
              enddo
            enddo
          endif

        endif     ! end if_laerlw_block

      enddo  lab_do_IMAX

! =================
      contains
! =================

!--------------------------------
      subroutine radclimaer
!................................

!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

!  ==================================================================  !
!                                                                      !
!  compute aerosols optical properties in NSWLWBD bands. there are     !
!  seven different vertical profile structures. in the troposphere,    !
!  aerosol distribution at each grid point is composed from up to      !
!  six components out of a total of ten different substances.          !
!                                                                      !
!  ref: wmo report wcp-112 (1986)                                      !
!                                                                      !
!  input variables:                                                    !
!     cmix   - mixing ratioes of aerosol components  -     NCM         !
!     denn   - aerosol number densities              -     2           !
!     rh1    - relative humidity                     -     NLAY        !
!     delz   - effective layer thickness             km    NLAY        !
!     idmaer - aerosol domain index                  -     NLAY        !
!     NXC    - number of different aerosol components-     1           !
!     NLAY   - vertical dimensions                   -     1           !
!                                                                      !
!  output variables:                                                   !
!     tauae  - optical depth                         -     NLAY*NSWLWBD!
!     ssaae  - single scattering albedo              -     NLAY*NSWLWBD!
!     asyae  - asymmetry parameter                   -     NLAY*NSWLWBD!
!!    aerodp - vertically integrated aer-opt-depth   -     IMAX*NSPC+1 !
!                                                                      !
!  ==================================================================  !
!
      real (kind=kind_phys) :: crt1, crt2
      parameter (crt1=30.0, crt2=0.03333)

!  ---  inputs:
!  ---  outputs:

!  ---  locals:
      real (kind=kind_phys) :: cm, hd, hdi, sig0u, sig0l, ratio, tt0,   &
     &      ex00, sc00, ss00, as00, ex01, sc01, ss01, as01,     tt1,    &
     &      ex02, sc02, ss02, as02, ex03, sc03, ss03, as03,     tt2,    &
     &      ext1, sca1, ssa1, asy1, drh0, drh1, rdrh

      integer :: ih1, ih2, kk, idom, icmp, ib, ii, ic, ic1

!===> ...  begin here

!     spcodp = f_zero

!===> ... loop over vertical layers from top to surface

      lab_do_layer : do kk = 1, NLAY

! --- linear interp coeffs for rh-dep species

        ih2 = 1
        do while ( rh1(kk) > rhlev(ih2) )
          ih2 = ih2 + 1
          if ( ih2 > NRHLEV ) exit
        enddo
        ih1 = max( 1, ih2-1 )
        ih2 = min( NRHLEV, ih2 )

        drh0 = rhlev(ih2) - rhlev(ih1)
        drh1 = rh1(kk) - rhlev(ih1)
        if ( ih1 == ih2 ) then
          rdrh = f_zero
        else
          rdrh = drh1 / drh0
        endif

! --- assign optical properties in each domain

        idom = idmaer(kk)

        lab_if_idom : if (idom == 5) then
! --- 5th domain - upper stratosphere assume no aerosol

          do ib = 1, NSWLWBD
            tauae(kk,ib) = f_zero
            if ( ib <= NSWBND ) then
              ssaae(kk,ib) = 0.99
              asyae(kk,ib) = 0.696
            else
              ssaae(kk,ib) = 0.5
              asyae(kk,ib) = 0.3
            endif
          enddo

        elseif (idom == 4) then    lab_if_idom
! --- 4th domain - stratospheric layers

          do ib = 1, NSWLWBD
            tauae(kk,ib) = extstra(ib) * delz(kk)
            if ( ib <= NSWBND ) then
              ssaae(kk,ib) = 0.99
              asyae(kk,ib) = 0.696
            else
              ssaae(kk,ib) = 0.5
              asyae(kk,ib) = 0.3
            endif
          enddo

! --- compute aod from individual species' contribution (optional)
!         idx = idxspc(10)             ! for sulfate
!         if ( lspcaod ) then
!           spcodp(idx) = spcodp(idx) + tauae(kk,nv_aod)
!         endif

        elseif (idom == 3) then    lab_if_idom
! --- 3rd domain - free tropospheric layers
!   1:inso 0.17e-3; 2:soot 0.4; 7:waso 0.59983; n:730

          do ib = 1, NSWLWBD
            ex01 = extrhi(1,ib)
            sc01 = scarhi(1,ib)
            ss01 = ssarhi(1,ib)
            as01 = asyrhi(1,ib)

            ex02 = extrhi(2,ib)
            sc02 = scarhi(2,ib)
            ss02 = ssarhi(2,ib)
            as02 = asyrhi(2,ib)

            ex03 = extrhd(ih1,1,ib)                                     &
     &           + rdrh * (extrhd(ih2,1,ib) - extrhd(ih1,1,ib))
            sc03 = scarhd(ih1,1,ib)                                     &
     &           + rdrh * (scarhd(ih2,1,ib) - scarhd(ih1,1,ib))
            ss03 = ssarhd(ih1,1,ib)                                     &
     &           + rdrh * (ssarhd(ih2,1,ib) - ssarhd(ih1,1,ib))
            as03 = asyrhd(ih1,1,ib)                                     &
     &           + rdrh * (asyrhd(ih2,1,ib) - asyrhd(ih1,1,ib))

            ext1 = 0.17e-3*ex01 + 0.4*ex02 + 0.59983*ex03
            sca1 = 0.17e-3*sc01 + 0.4*sc02 + 0.59983*sc03
            ssa1 = 0.17e-3*ss01*ex01 + 0.4*ss02*ex02 + 0.59983*ss03*ex03
            asy1 = 0.17e-3*as01*sc01 + 0.4*as02*sc02 + 0.59983*as03*sc03

            tauae(kk,ib) = ext1 * 730.0 * delz(kk)
            ssaae(kk,ib) = min(f_one, ssa1/ext1)
            asyae(kk,ib) = min(f_one, asy1/sca1)

! --- compute aod from individual species' contribution (optional)
!           if ( lspcaod .and. ib==nv_aod ) then
!             spcodp(1) = spcodp(1) + 0.17e-3*ex01*730.0*delz(kk)   ! dust (inso)   #1
!             spcodp(2) = spcodp(2) + 0.4    *ex02*730.0*delz(kk)   ! black carbon  #2
!             spcodp(3) = spcodp(3) + 0.59983*ex03*730.0*delz(kk)   ! water soluble #7
!           endif

          enddo

        elseif (idom == 1) then    lab_if_idom
! --- 1st domain - mixing layer

          lab_do_ib : do ib = 1, NSWLWBD
            ext1 = f_zero
            sca1 = f_zero
            ssa1 = f_zero
            asy1 = f_zero

            lab_do_icmp : do icmp = 1, NCM
              ic = icmp
!             idx = idxspc(icmp)

              cm = cmix(icmp)
              lab_if_cm : if ( cm > f_zero ) then

                lab_if_ic : if ( ic <= NCM1 ) then        ! component withour rh dep
                  tt0  = cm * extrhi(ic,ib)
                  ext1 = ext1 + tt0
                  sca1 = sca1 + cm * scarhi(ic,ib)
                  ssa1 = ssa1 + cm * ssarhi(ic,ib) * extrhi(ic,ib)
                  asy1 = asy1 + cm * asyrhi(ic,ib) * scarhi(ic,ib)
                else  lab_if_ic                           ! component with rh dep
                  ic1 = ic - NCM1

                  ex00 = extrhd(ih1,ic1,ib)                             &
     &               + rdrh * (extrhd(ih2,ic1,ib) - extrhd(ih1,ic1,ib))
                  sc00 = scarhd(ih1,ic1,ib)                             &
     &               + rdrh * (scarhd(ih2,ic1,ib) - scarhd(ih1,ic1,ib))
                  ss00 = ssarhd(ih1,ic1,ib)                             &
     &               + rdrh * (ssarhd(ih2,ic1,ib) - ssarhd(ih1,ic1,ib))
                  as00 = asyrhd(ih1,ic1,ib)                             &
     &               + rdrh * (asyrhd(ih2,ic1,ib) - asyrhd(ih1,ic1,ib))

                  tt0  = cm * ex00
                  ext1 = ext1 + tt0
                  sca1 = sca1 + cm * sc00
                  ssa1 = ssa1 + cm * ss00 * ex00
                  asy1 = asy1 + cm * as00 * sc00
                endif  lab_if_ic

! --- compute aod from individual species' contribution (optional)
!               if ( lspcaod .and. ib==nv_aod ) then
!                 spcodp(idx) = spcodp(idx) + tt0*denn(1)*delz(kk)   ! idx for dif species
!               endif

              endif  lab_if_cm
            enddo  lab_do_icmp

            tauae(kk,ib) = ext1 * denn(1) * delz(kk)
            ssaae(kk,ib) = min(f_one, ssa1/ext1)
            asyae(kk,ib) = min(f_one, asy1/sca1)
          enddo  lab_do_ib

        elseif (idom == 2) then    lab_if_idom
! --- 2nd domain - mineral transport layers

          do ib = 1, NSWLWBD
            tauae(kk,ib) = extrhi(6,ib) * denn(2) * delz(kk)
            ssaae(kk,ib) = ssarhi(6,ib)
            asyae(kk,ib) = asyrhi(6,ib)
          enddo

! --- compute aod from individual species' contribution (optional)
!         if ( laersw ) then
!            spcodp(1) = spcodp(1) + tauae(kk,nv_aod)            ! dust
!         endif

        else  lab_if_idom
! --- domain index out off range, assume no aerosol

          do ib = 1, NSWLWBD
            tauae(kk,ib) = f_zero
            ssaae(kk,ib) = f_one
            asyae(kk,ib) = f_zero
          enddo

!         write(6,19) kk,idom
! 19      format(/'  ***  ERROR in sub AEROS: domain index out'         &
!    &,            ' of range!  K, IDOM =',3i5,' ***')
!         stop 19

        endif  lab_if_idom

      enddo  lab_do_layer
!
!===> ... smooth profile at domain boundaries
!
      if ( ivflip == 0 ) then    ! input from toa to sfc

        do ib = 1, NSWLWBD
        do kk = 2, NLAY
          if ( tauae(kk,ib) > f_zero ) then
            ratio = tauae(kk-1,ib) / tauae(kk,ib)
          else
            ratio = f_one
          endif

          tt0 = tauae(kk,ib) + tauae(kk-1,ib)
          tt1 = 0.2 * tt0
          tt2 = tt0 - tt1

          if ( ratio > crt1 ) then
            tauae(kk,ib)   = tt1
            tauae(kk-1,ib) = tt2
          endif

          if ( ratio < crt2 ) then
            tauae(kk,ib)   = tt2
            tauae(kk-1,ib) = tt1
          endif
        enddo   ! do_kk_loop
        enddo   ! do_ib_loop

      else                      ! input from sfc to toa

        do ib = 1, NSWLWBD
        do kk = NLAY-1, 1, -1
          if ( tauae(kk,ib) > f_zero ) then
            ratio = tauae(kk+1,ib) / tauae(kk,ib)
          else
            ratio = f_one
          endif

          tt0 = tauae(kk,ib) + tauae(kk+1,ib)
          tt1 = 0.2 * tt0
          tt2 = tt0 - tt1

          if ( ratio > crt1 ) then
            tauae(kk,ib)   = tt1
            tauae(kk+1,ib) = tt2
          endif

          if ( ratio < crt2 ) then
            tauae(kk,ib)   = tt2
            tauae(kk+1,ib) = tt1
          endif
        enddo   ! do_kk_loop
        enddo   ! do_ib_loop

      endif

!
      return
!................................
      end subroutine radclimaer
!--------------------------------
!
!...................................
      end subroutine aer_property
!-----------------------------------


!...............................................!
      end module module_radiation_aerosols_nmmb !
!===============================================!
