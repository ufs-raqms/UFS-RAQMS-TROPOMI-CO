!!!!!  ==============================================================  !!!!!
!!!!!              sw-rrtm3 radiation package description              !!!!!
!!!!!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtm-sw radiation    !
!   code from aer inc.                                                     !
!                                                                          !
!   the sw-rrtm3 package includes these parts:                             !
!                                                                          !
!      'radsw_rrtm3_param.f'                                               !
!      'radsw_rrtm3_datatb.f'                                              !
!      'radsw_rrtm3_main.f'                                                !
!                                                                          !
!   the 'radsw_rrtm3_param.f' contains:                                    !
!                                                                          !
!      'module_radsw_parameters'  -- band parameters set up                !
!                                                                          !
!   the 'radsw_rrtm3_datatb.f' contains:                                   !
!                                                                          !
!      'module_radsw_ref'         -- reference temperature and pressure    !
!      'module_radsw_cldprtb'     -- cloud property coefficients table     !
!      'module_radsw_sflux'       -- spectral distribution of solar flux   !
!      'module_radsw_kgbnn'       -- absorption coeffients for 14          !
!                                    bands, where nn = 16-29               !
!                                                                          !
!   the 'radsw_rrtm3_main.f' contains:                                     !
!                                                                          !
!      'module_radsw_main'        -- main sw radiation transfer            !
!                                                                          !
!   in the main module 'module_radsw_main' there are only two              !
!   externally callable subroutines:                                       !
!                                                                          !
!      'swrad'      -- main sw radiation routine                           !
!         inputs:                                                          !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                         !
!            clouds,icseed,aerosols,sfcalb,                                !
!            cosz,solcon,NDAY,idxday,                                      !
!            npts, nlay, nlp1, lprnt,                                      !
!         outputs:                                                         !
!            hswc,topflx,sfcflx,                                           !
!!        optional outputs:                                                !
!            HSW0,HSWB,FLXPRF,FDNCMP)                                      !
!           )                                                              !
!                                                                          !
!      'rswinit'    -- initialization routine                              !
!         inputs:                                                          !
!           ( me )                                                         !
!         outputs:                                                         !
!           (none)                                                         !
!                                                                          !
!   all the sw radiation subprograms become contained subprograms          !
!   in module 'module_radsw_main' and many of them are not directly        !
!   accessable from places outside the module.                             !
!                                                                          !
!    derived data type constructs used:                                    !
!                                                                          !
!     1. radiation flux at toa: (from module 'module_radsw_parameters')    !
!          topfsw_type   -  derived data type for toa rad fluxes           !
!            upfxc              total sky upward flux at toa               !
!            dnfxc              total sky downward flux at toa             !
!            upfx0              clear sky upward flux at toa               !
!                                                                          !
!     2. radiation flux at sfc: (from module 'module_radsw_parameters')    !
!          sfcfsw_type   -  derived data type for sfc rad fluxes           !
!            upfxc              total sky upward flux at sfc               !
!            dnfxc              total sky downward flux at sfc             !
!            upfx0              clear sky upward flux at sfc               !
!            dnfx0              clear sky downward flux at sfc             !
!                                                                          !
!     3. radiation flux profiles(from module 'module_radsw_parameters')    !
!          profsw_type    -  derived data type for rad vertical prof       !
!            upfxc              level upward flux for total sky            !
!            dnfxc              level downward flux for total sky          !
!            upfx0              level upward flux for clear sky            !
!            dnfx0              level downward flux for clear sky          !
!                                                                          !
!     4. surface component fluxes(from module 'module_radsw_parameters'    !
!          cmpfsw_type    -  derived data type for component sfc flux      !
!            uvbfc              total sky downward uv-b flux at sfc        !
!            uvbf0              clear sky downward uv-b flux at sfc        !
!            nirbm              surface downward nir direct beam flux      !
!            nirdf              surface downward nir diffused flux         !
!            visbm              surface downward uv+vis direct beam flx    !
!            visdf              surface downward uv+vis diffused flux      !
!                                                                          !
!   external modules referenced:                                           !
!                                                                          !
!       'module physparam'                                                 !
!       'module physcons'                                                  !
!       'mersenne_twister'                                                 !
!                                                                          !
!   compilation sequence is:                                               !
!                                                                          !
!      'radsw_rrtm3_param.f'                                               !
!      'radsw_rrtm3_datatb.f'                                              !
!      'radsw_rrtm3_main.f'                                                !
!                                                                          !
!   and all should be put in front of routines that use sw modules         !
!                                                                          !
!==========================================================================!
!                                                                          !
!   the original program declarations:                                     !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
!  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  !
!  This software may be used, copied, or redistributed as long as it is    !
!  not sold and this copyright notice is reproduced on each copy made.     !
!  This model is provided as is without any express or implied warranties. !
!                       (http://www.rtweb.aer.com/)                        !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
! ************************************************************************ !
!                                                                          !
!                              rrtmg_sw                                    !
!                                                                          !
!                                                                          !
!                   a rapid radiative transfer model                       !
!                    for the solar spectral region                         !
!            atmospheric and environmental research, inc.                  !
!                        131 hartwell avenue                               !
!                        lexington, ma 02421                               !
!                                                                          !
!                           eli j. mlawer                                  !
!                        jennifer s. delamere                              !
!                         michael j. iacono                                !
!                         shepard a. clough                                !
!                                                                          !
!                                                                          !
!                       email:  miacono@aer.com                            !
!                       email:  emlawer@aer.com                            !
!                       email:  jdelamer@aer.com                           !
!                                                                          !
!        the authors wish to acknowledge the contributions of the          !
!        following people:  steven j. taubman, patrick d. brown,           !
!        ronald e. farren, luke chen, robert bergstrom.                    !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    references:                                                           !
!    (rrtm_sw/rrtmg_sw):                                                   !
!      clough, s.a., m.w. shephard, e.j. mlawer, j.s. delamere,            !
!      m.j. iacono, k. cady-pereira, s. boukabara, and p.d. brown:         !
!      atmospheric radiative transfer modeling: a summary of the aer       !
!      codes, j. quant. spectrosc. radiat. transfer, 91, 233-244, 2005.    !
!                                                                          !
!    (mcica):                                                              !
!      pincus, r., h. w. barker, and j.-j. morcrette: a fast, flexible,    !
!      approximation technique for computing radiative transfer in         !
!      inhomogeneous cloud fields, j. geophys. res., 108(d13), 4376,       !
!      doi:10.1029/2002jd003322, 2003.                                     !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    aer's revision history:                                               !
!     this version of rrtmg_sw has been modified from rrtm_sw to use a     !
!     reduced set of g-point intervals and a two-stream model for          !
!     application to gcms.                                                 !
!                                                                          !
! --  original version (derived from rrtm_sw)                              !
!        2002: aer. inc.                                                   !
! --  conversion to f90 formatting; addition of 2-stream radiative transfer!
!        feb 2003: j.-j. morcrette, ecmwf                                  !
! --  additional modifications for gcm application                         !
!        aug 2003: m. j. iacono, aer inc.                                  !
! --  total number of g-points reduced from 224 to 112.  original          !
!     set of 224 can be restored by exchanging code in module parrrsw.f90  !
!     and in file rrtmg_sw_init.f90.                                       !
!        apr 2004: m. j. iacono, aer, inc.                                 !
! --  modifications to include output for direct and diffuse               !
!     downward fluxes.  there are output as "true" fluxes without          !
!     any delta scaling applied.  code can be commented to exclude         !
!     this calculation in source file rrtmg_sw_spcvrt.f90.                 !
!        jan 2005: e. j. mlawer, m. j. iacono, aer, inc.                   !
! --  revised to add mcica capability.                                     !
!        nov 2005: m. j. iacono, aer, inc.                                 !
! --  reformatted for consistency with rrtmg_lw.                           !
!        feb 2007: m. j. iacono, aer, inc.                                 !
! --  modifications to formatting to use assumed-shape arrays.             !
!        aug 2007: m. j. iacono, aer, inc.                                 !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!   ncep modifications history log:                                        !
!                                                                          !
!       sep 2003,  yu-tai hou        -- received aer's rrtm-sw gcm version !
!                    code (v224)                                           !
!       nov 2003,  yu-tai hou        -- corrected errors in direct/diffuse !
!                    surface alabedo components.                           !
!       jan 2004,  yu-tai hou        -- modified code into standard modular!
!                    f9x code for ncep models. the original three cloud    !
!                    control flags are simplified into two: iflagliq and   !
!                    iflagice. combined the org subr sw_224 and setcoef    !
!                    into radsw (the main program); put all kgb##together  !
!                    and reformat into a separated data module; combine    !
!                    reftra and vrtqdr as swflux; optimized taumol and all !
!                    taubgs to form a contained subroutines.               !
!       jun 2004,  yu-tai hou        -- modified code based on aer's faster!
!                    version rrtmg_sw (v2.0) with 112 g-points.            !
!       mar 2005,  yu-tai hou        -- modified to aer v2.3, correct cloud!
!                    scaling error, total sky properties are delta scaled  !
!                    after combining clear and cloudy parts. the testing   !
!                    criterion of ssa is saved before scaling. added cloud !
!                    layer rain and snow contributions. all cloud water    !
!                    partical contents are treated the same way as other   !
!                    atmos particles.                                      !
!       apr 2005,  yu-tai hou        -- modified on module structures (this!
!                    version of code was given back to aer in jun 2006)    !
!       nov 2006,  yu-tai hou        -- modified code to include the       !
!                    generallized aerosol optical property scheme for gcms.!
!       apr 2007,  yu-tai hou        -- added spectral band heating as an  !
!                    optional output to support the 500km model's upper    !
!                    stratospheric radiation calculations. restructure     !
!                    optional outputs for easy access by different models. !
!       oct 2008,  yu-tai hou        -- modified to include new features   !
!                    from aer's newer release v3.5-v3.61, including mcica  !
!                    sub-grid cloud option and true direct/diffuse fluxes  !
!                    without delta scaling. added rain/snow opt properties !
!                    support to cloudy sky calculations. simplified and    !
!                    unified sw and lw sub-column cloud subroutines into   !
!                    one module by using optional parameters.              !
!       mar 2009,  yu-tai hou        -- replaced the original random number!
!                    generator coming with the original code with ncep w3  !
!                    library to simplify the program and moved sub-column  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       mar 2009,  yu-tai hou        -- replaced random number generator   !
!                    programs coming from the original code with the ncep  !
!                    w3 library to simplify the program and moved sub-col  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       nov 2009,  yu-tai hou        -- updated to aer v3.7-v3.8 version.  !
!                    notice the input cloud ice/liquid are assumed as      !
!                    in-cloud quantities, not grid average quantities.     !
!       aug 2010,  yu-tai hou        -- uptimized code to improve efficiency
!                    splited subroutine spcvrt into two subs, spcvrc and   !
!                    spcvrm, to handling non-mcica and mcica type of calls.!
!       apr 2012,  b. ferrier and y. hou -- added conversion factor to fu's!
!                    cloud-snow optical property scheme.                   !
!       jul 2012,  s. moorthi and Y. hou  -- eliminated the pointer array  !
!                     in subr 'spcvrt' for multi-threading issue running   !
!                     under intel's fortran compiler.                      !
!       nov 2012,  yu-tai hou        -- modified control parameters thru   !
!                     module 'physparam'.                                  !
!       jun 2013,  yu-tai hou        -- moving band 9 surface treatment    !
!                     back as in the rrtm2 version, spliting surface flux  !
!                     into two spectral regions (vis & nir), instead of    !
!                     designated it in nir region only.                    !
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!


!========================================!
      module module_radsw_main_nmmb      !
!........................................!
!
      USE ESMF

      use physparam,        only : iswrate, iswrgas, iswcliq, iswcice,  &
     &                             isubcsw, icldflg, iovrsw,  ivflip,   &
     &                             iswmode, kind_phys, kind_taum
      use physcons,         only : con_g, con_cp, con_avgd, con_amd,    &
     &                             con_amw, con_amo3

      use module_radsw_parameters
      use mersenne_twister, only : random_setseed, random_number,       &
     &                             random_stat
      use module_radsw_ref, only : preflog, tref
      use module_radsw_sflux
!
      implicit none
!
      private
!
!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGSW='NCEP SW v5.1  Nov 2012 -RRTMG-SW v3.8   '
!    &   VTAGSW='NCEP SW v5.0  Aug 2012 -RRTMG-SW v3.8   '
!    &   VTAGSW='RRTMG-SW v3.8   Nov 2009'
!    &   VTAGSW='RRTMG-SW v3.7   Nov 2009'
!    &   VTAGSW='RRTMG-SW v3.61  Oct 2008'
!    &   VTAGSW='RRTMG-SW v3.5   Oct 2008'
!    &   VTAGSW='RRTM-SW 112v2.3 Apr 2007'
!    &   VTAGSW='RRTM-SW 112v2.3 Mar 2005'
!    &   VTAGSW='RRTM-SW 112v2.0 Jul 2004'

!  ---  constant values
      real (kind=kind_phys), parameter :: eps     = 1.0e-6
      real (kind=kind_phys), parameter :: oneminus= 1.0 - eps
      real (kind=kind_phys), parameter :: bpade   = 1.0/0.278  ! pade approx constant
      real (kind=kind_phys), parameter :: stpfac  = 296.0/1013.0
      real (kind=kind_phys), parameter :: ftiny   = 1.0e-12
      real (kind=kind_phys), parameter :: s0      = 1368.22    ! internal solar const
                                                               ! adj through input
      real (kind=kind_phys), parameter :: f_zero  = 0.0
      real (kind=kind_phys), parameter :: f_one   = 1.0

!  ---  atomic weights for conversion from mass to volume mixing ratios
      real (kind=kind_phys), parameter :: amdw    = con_amd/con_amw
      real (kind=kind_phys), parameter :: amdo3   = con_amd/con_amo3

!  ---  band indices
      integer, dimension(nblow:nbhgh) :: nspa, nspb, idxebc, idxsfc

      data nspa(:) /  9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1 /
      data nspb(:) /  1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1 /

!     data idxsfc(:) / 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1 /  ! band index for sfc flux
      data idxsfc(:) / 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1 /  ! band index for sfc flux
      data idxebc(:) / 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1, 5 /  ! band index for cld prop

!  ---  band wavenumber intervals
!     real (kind=kind_phys), dimension(nblow:nbhgh):: wavenum1,wavenum2
!     data wavenum1(:)  /                                               &
!    &         2600.0, 3250.0, 4000.0, 4650.0, 5150.0, 6150.0, 7700.0,  &
!    &         8050.0,12850.0,16000.0,22650.0,29000.0,38000.0,  820.0 /
!     data wavenum2(:)  /                                               &
!              3250.0, 4000.0, 4650.0, 5150.0, 6150.0, 7700.0, 8050.0,  &
!    &        12850.0,16000.0,22650.0,29000.0,38000.0,50000.0, 2600.0 /
!     real (kind=kind_phys), dimension(nblow:nbhgh) :: delwave
!     data delwave(:)   /                                               &
!    &          650.0,  750.0,  650.0,  500.0, 1000.0, 1550.0,  350.0,  &
!    &         4800.0, 3150.0, 6650.0, 6350.0, 9000.0,12000.0, 1780.0 /

      integer, parameter :: nuvb = 27            !uv-b band index

!! ---  logical flags for optional output fields

      logical :: lhswb  = .false.
      logical :: lhsw0  = .false.
      logical :: lflxprf= .false.
      logical :: lfdncmp= .false.
!$omp threadprivate(lhswb,lhsw0,lflxprf,lfdncmp)

!  ---  those data will be set up only once by "rswinit"

      real (kind=kind_phys) :: exp_tbl(0:NTBMX)

!  ...  heatfac is the factor for heating rates
!       (in k/day, or k/sec set by subroutine 'rswinit')

      real (kind=kind_phys) :: heatfac

!  ---  the following variables are used for sub-column cloud scheme

      integer, parameter :: ipsdsw0 = 1          ! initial permutation seed
!
!rv   integer, parameter :: CHK= CHNK_RRTM      ! CHNK_RRTM is macro defined in configure.nems
      integer, parameter :: CHK= 1              ! CHNK_RRTM is macro defined in configure.nems


        logical buggal

!  ---  public accessible subprograms

      public swrad, rswinit
      public buggal, buggalon, buggaloff, CHK

! =================
      contains
! =================

      subroutine buggalon
          buggal=.true.
      end subroutine buggalon
      subroutine buggaloff
          buggal=.false.
      end subroutine buggaloff


!-----------------------------------
      subroutine swrad                                                  &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icseed,aerosols,sfcalb,                             &
     &       cosz,solcon,NDAY,idxday,                                   &
     &       npts, nlay, nlp1, lprnt,                                   &
!  ---  outputs:
     &       hswc,topflx,sfcflx                                         &
!! ---  optional:
     &,      HSW0,HSWB,FLXPRF,FDNCMP                                    &
     &     )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!   plyr (npts,nlay) : model layer mean pressure in mb                  !
!   plvl (npts,nlp1) : model level pressure in mb                       !
!   tlyr (npts,nlay) : model layer mean temperature in k                !
!   tlvl (npts,nlp1) : model level temperature in k    (not in use)     !
!   qlyr (npts,nlay) : layer specific humidity in gm/gm   *see inside   !
!   olyr (npts,nlay) : layer ozone concentration in gm/gm               !
!   gasvmr(npts,nlay,:): atmospheric constent gases:                    !
!                      (check module_radiation_gases for definition)    !
!      gasvmr(:,:,1)  - co2 volume mixing ratio                         !
!      gasvmr(:,:,2)  - n2o volume mixing ratio                         !
!      gasvmr(:,:,3)  - ch4 volume mixing ratio                         !
!      gasvmr(:,:,4)  - o2  volume mixing ratio                         !
!      gasvmr(:,:,5)  - co  volume mixing ratio        (not used)       !
!      gasvmr(:,:,6)  - cfc11 volume mixing ratio      (not used)       !
!      gasvmr(:,:,7)  - cfc12 volume mixing ratio      (not used)       !
!      gasvmr(:,:,8)  - cfc22 volume mixing ratio      (not used)       !
!      gasvmr(:,:,9)  - ccl4  volume mixing ratio      (not used)       !
!   clouds(npts,nlay,:): cloud profile                                  !
!                      (check module_radiation_clouds for definition)   !
!                ---  for  iswcliq > 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer in-cloud liq water path   (g/m**2)     !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)     !
!       clouds(:,:,4)  -   layer in-cloud ice water path   (g/m**2)     !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)     !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)     !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)     !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)     !
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)     !
!                ---  for  iswcliq = 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer cloud optical depth                    !
!       clouds(:,:,3)  -   layer cloud single scattering albedo         !
!       clouds(:,:,4)  -   layer cloud asymmetry factor                 !
!     icseed(npts)   : auxiliary special cloud related array            !
!                      when module variable isubcsw=2, it provides      !
!                      permutation seed for each column profile that    !
!                      are used for generating random numbers.          !
!                      when isubcsw /=2, it will not be used.           !
!   aerosols(npts,nlay,nbdsw,:) : aerosol optical properties            !
!                      (check module_radiation_aerosols for definition) !
!         (:,:,:,1)   - optical depth                                   !
!         (:,:,:,2)   - single scattering albedo                        !
!         (:,:,:,3)   - asymmetry parameter                             !
!   sfcalb(npts, : ) : surface albedo in fraction                       !
!                      (check module_radiation_surface for definition)  !
!         ( :, 1 )    - near ir direct beam albedo                      !
!         ( :, 2 )    - near ir diffused albedo                         !
!         ( :, 3 )    - uv+vis direct beam albedo                       !
!         ( :, 4 )    - uv+vis diffused albedo                          !
!   cosz  (npts)     : cosine of solar zenith angle                     !
!   solcon           : solar constant                      (w/m**2)     !
!   NDAY             : num of daytime points                            !
!   idxday(npts)     : index array for daytime points                   !
!   npts             : number of horizontal points                      !
!   nlay,nlp1        : vertical layer/lavel numbers                     !
!   lprnt            : logical check print flag                         !
!                                                                       !
!  output variables:                                                    !
!   hswc  (npts,nlay): total sky heating rates (k/sec or k/day)         !
!   topflx(npts)     : radiation fluxes at toa (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at toa                   !
!     dnflx            - total sky downward flux at toa                 !
!     upfx0            - clear sky upward flux at toa                   !
!   sfcflx(npts)     : radiation fluxes at sfc (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at sfc                   !
!     dnfxc            - total sky downward flux at sfc                 !
!     upfx0            - clear sky upward flux at sfc                   !
!     dnfx0            - clear sky downward flux at sfc                 !
!                                                                       !
!!optional outputs variables:                                           !
!   hswb(npts,nlay,nbdsw): spectral band total sky heating rates        !
!   hsw0  (npts,nlay): clear sky heating rates (k/sec or k/day)         !
!   flxprf(npts,nlp1): level radiation fluxes (w/m**2), components:     !
!                      (check module_radsw_parameters for definition)   !
!     dnfxc            - total sky downward flux at interface           !
!     upfxc            - total sky upward flux at interface             !
!     dnfx0            - clear sky downward flux at interface           !
!     upfx0            - clear sky upward flux at interface             !
!   fdncmp(npts)     : component surface downward fluxes (w/m**2):      !
!                      (check module_radsw_parameters for definition)   !
!     uvbfc            - total sky downward uv-b flux at sfc            !
!     uvbf0            - clear sky downward uv-b flux at sfc            !
!     nirbm            - downward surface nir direct beam flux          !
!     nirdf            - downward surface nir diffused flux             !
!     visbm            - downward surface uv+vis direct beam flux       !
!     visdf            - downward surface uv+vis diffused flux          !
!                                                                       !
!  external module variables:  (in physparam)                           !
!   iswrgas - control flag for rare gases (ch4,n2o,o2, etc.)            !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   iswcliq - control flag for liq-cloud optical properties             !
!           =0: input cloud optical depth, fixed ssa, asy               !
!           =1: use hu and stamnes(1993) method for liq cld             !
!           =2: not used                                                !
!   iswcice - control flag for ice-cloud optical properties             !
!           *** if iswcliq==0, iswcice is ignored                       !
!           =1: use ebert and curry (1992) scheme for ice clouds        !
!           =2: use streamer v3.0 (2001) method for ice clouds          !
!           =3: use fu's method (1996) for ice clouds                   !
!   iswmode - control flag for 2-stream transfer scheme                 !
!           =1; delta-eddington    (joseph et al., 1976)                !
!           =2: pifm               (zdunkowski et al., 1980)            !
!           =3: discrete ordinates (liou, 1973)                         !
!   isubcsw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   iovrsw  - cloud overlapping control flag                            !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud                                   !
!   ivflip  - control flg for direction of vertical index               !
!           =0: index from toa to surface                               !
!           =1: index from surface to toa                               !
!                                                                       !
!  module parameters, control variables:                                !
!     nblow,nbhgh      - lower and upper limits of spectral bands       !
!     maxgas           - maximum number of absorbing gaseous            !
!     ngptsw           - total number of g-point subintervals           !
!     ng##             - number of g-points in band (##=16-29)          !
!     ngb(ngptsw)      - band indices for each g-point                  !
!     bpade            - pade approximation constant (1/0.278)          !
!     nspa,nspb(nblow:nbhgh)                                            !
!                      - number of lower/upper ref atm's per band       !
!     ipsdsw0          - permutation seed for mcica sub-col clds        !
!                                                                       !
!  major local variables:                                               !
!     pavel  (nlay)         - layer pressures (mb)                      !
!     delp   (nlay)         - layer pressure thickness (mb)             !
!     tavel  (nlay)         - layer temperatures (k)                    !
!     coldry (nlay)         - dry air column amount                     !
!                                   (1.e-20*molecules/cm**2)            !
!     cldfrc (nlay)         - layer cloud fraction (norm by tot cld)    !
!     cldfmc (nlay,ngptsw)  - layer cloud fraction for g-point          !
!     taucw  (nlay,nbdsw)   - cloud optical depth                       !
!     ssacw  (nlay,nbdsw)   - cloud single scattering albedo (weighted) !
!     asycw  (nlay,nbdsw)   - cloud asymmetry factor         (weighted) !
!     tauaer (nlay,nbdsw)   - aerosol optical depths                    !
!     ssaaer (nlay,nbdsw)   - aerosol single scattering albedo          !
!     asyaer (nlay,nbdsw)   - aerosol asymmetry factor                  !
!     colamt (nlay,maxgas)  - column amounts of absorbing gases         !
!                             1 to maxgas are for h2o, co2, o3, n2o,    !
!                             ch4, o2, co, respectively (mol/cm**2)     !
!     facij  (nlay)         - indicator of interpolation factors        !
!                             =0/1: indicate lower/higher temp & height !
!     selffac(nlay)         - scale factor for self-continuum, equals   !
!                          (w.v. density)/(atm density at 296K,1013 mb) !
!     selffrac(nlay)        - factor for temp interpolation of ref      !
!                             self-continuum data                       !
!     indself(nlay)         - index of the lower two appropriate ref    !
!                             temp for the self-continuum interpolation !
!     forfac (nlay)         - scale factor for w.v. foreign-continuum   !
!     forfrac(nlay)         - factor for temp interpolation of ref      !
!                             w.v. foreign-continuum data               !
!     indfor (nlay)         - index of the lower two appropriate ref    !
!                             temp for the foreign-continuum interp     !
!     laytrop               - layer at which switch is made from one    !
!                             combination of key species to another     !
!     jp(nlay),jt(nlay),jt1(nlay)                                       !
!                           - lookup table indexes                      !
!     flxucb(nlp1,nbdsw)    - spectral bnd total-sky upward flx (w/m2)  !
!     flxdcb(nlp1,nbdsw)    - spectral bnd total-sky downward flx (w/m2)!
!     flxu0b(nlp1,nbdsw)    - spectral bnd clear-sky upward flx (w/m2)  !
!     flxd0b(nlp1,nbdsw)    - spectral b d clear-sky downward flx (w/m2)!
!                                                                       !
!                                                                       !
!  =====================    end of definitions    ====================  !

!  ---  inputs:
      integer, intent(in) :: npts, nlay, nlp1, NDAY

      integer, dimension(:), intent(in) :: idxday, icseed

      logical, intent(in) :: lprnt

      real (kind=kind_phys), dimension(CHK ,nlp1), intent(in) ::        &
     &       plvl, tlvl
      real (kind=kind_phys), dimension(CHK ,nlay), intent(in) ::        &
     &       plyr, tlyr, qlyr, olyr
      real (kind=kind_phys), dimension(CHK ,4),    intent(in) :: sfcalb

      real (kind=kind_phys), dimension(CHK ,nlay,9),intent(in):: gasvmr
      real (kind=kind_phys), dimension(CHK ,nlay,9),intent(in):: clouds
      real (kind=kind_phys), dimension(CHK ,nlay,nbdsw,3),intent(in)::  &
     &       aerosols

      real (kind=kind_phys), intent(in) :: cosz(CHK ), solcon


!  ---  outputs:
      real (kind=kind_phys), dimension(CHK ,nlay), intent(out) :: hswc

      type (topfsw_type),    dimension(CHK ), intent(out) :: topflx
      type (sfcfsw_type),    dimension(CHK ), intent(out) :: sfcflx

!! ---  optional outputs:
      real (kind=kind_phys), dimension(CHK ,nlay,nbdsw), optional,      &
     &       intent(out) :: hswb

      real (kind=kind_phys), dimension(CHK ,nlay),       optional,      &
     &       intent(out) :: hsw0
      type (profsw_type),    dimension(CHK ,nlp1),       optional,      &
     &       intent(out) :: flxprf
      type (cmpfsw_type),    dimension(CHK ),            optional,      &
     &       intent(out) :: fdncmp

!  ---  locals:
      real (kind=kind_phys), dimension( CHK     ,nlay,ngptsw) ::        &
     &       cldfmc,                                                    &
     &       taug, taur
      real (kind=kind_phys), dimension( CHK     ,nlp1,nbdsw)::          &
     &       fxupc, fxdnc,                                              &
     &       fxup0, fxdn0

      real (kind=kind_phys), dimension( CHK     ,nlay,nbdsw)  ::        &
     &       tauae, ssaae, asyae, taucw, ssacw, asycw

      real (kind=kind_phys), dimension( CHK     ,ngptsw) :: sfluxzen

      real (kind=kind_phys), dimension( CHK     ,nlay) :: cldfrc,delp,  &
     &       pavel, tavel, coldry, colmol, h2ovmr, o3vmr, temcol,       &
     &       cliqp, reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4,    &
     &       cfrac, fac00, fac01, fac10, fac11, forfac, forfrac,        &
     &       selffac, selffrac, rfdelp

      real (kind=kind_phys), dimension(CHK,nlp1) ::                     &
     &       fnet, flxdc, flxuc,                                        &
     &       flxd0, flxu0

      real (kind=kind_phys), dimension(CHK,2) ::                        &
     &       albbm, albdf, sfbmc,                                       &
     &       sfbm0, sfdfc, sfdf0

      real (kind=kind_phys),DIMENSION(CHK) ::                           &
     &       cosz1, sntz1,                                              &
     &       ssolar, zcf0, zcf1, ftoau0, ftoauc, ftoadc,                &
     &       fsfcu0, fsfcuc, fsfcd0, fsfcdc, suvbfc, suvbf0
      real (kind=kind_phys) :: tem0, tem1, tem2, s0fac,                 &
     &       fp, fp1, ft, ft1, zdpgcp


!  ---  column amount of absorbing gases:
!       (:,m) m = 1-h2o, 2-co2, 3-o3, 4-n2o, 5-ch4, 6-o2, 7-co
      real (kind=kind_phys) ::  colamt(CHK       ,nlay,maxgas)

      integer, dimension(CHK) :: ipseed
      integer, dimension(CHK       ,nlay) :: indfor,indself,jp,jt,jt1

      integer, dimension(CHK       ) :: laytrop
      integer :: i, ib, ipt, j1, k, kk, mb
!
!===> ... begin here

!
      IF ( NDAY .GT.  CHK      ) THEN
        write(0,*)'radsw_main: swrad aborting NDAY ',NDAY,              &
     &            '  > CHK      ', CHK     
        CALL abort
      ENDIF

      lhswb  = present ( hswb )
      lhsw0  = present ( hsw0 )
      lflxprf= present ( flxprf )
      lfdncmp= present ( fdncmp )

!  --- ...  compute solar constant adjustment factor according to solcon.
!      ***  s0, the solar constant at toa in w/m**2, is hard-coded with
!           each spectra band, the total flux is about 1368.22 w/m**2.

      s0fac = solcon / s0

!  --- ...  initial output arrays

      hswc(:,:) = f_zero
      topflx = topfsw_type ( f_zero, f_zero, f_zero )
      sfcflx = sfcfsw_type ( f_zero, f_zero, f_zero, f_zero )

!! --- ...  initial optional outputs
      if ( lflxprf ) then
        flxprf = profsw_type ( f_zero, f_zero, f_zero, f_zero )
      endif

      if ( lfdncmp ) then
        fdncmp = cmpfsw_type (f_zero,f_zero,f_zero,f_zero,f_zero,f_zero)
      endif

      if ( lhsw0 ) then
        hsw0(:,:) = f_zero
      endif

      if ( lhswb ) then
        hswb(:,:,:) = f_zero
      endif

!  --- ...  change random number seed value for each radiation invocation

      if     ( isubcsw == 1 ) then     ! advance prescribed permutation seed
        do i = 1, npts
          ipseed(i) = ipsdsw0 + i
        enddo
      elseif ( isubcsw == 2 ) then     ! use input array of permutaion seeds
        do i = 1, npts
          ipseed(i) = icseed(i)
        enddo
      endif

      if ( lprnt ) then
        write(0,*)'  In radsw, isubcsw, ipsdsw0,ipseed =',              &
     &           isubcsw, ipsdsw0, ipseed
      endif

!  --- ...  loop over each daytime grid point  (this loop has been demoted to inner)

!jm      lab_do_ipt : do ipt = 1, NDAY
!jm        j1 = idxday(ipt)
!#define RUNTIME_CHECKING
#ifndef RUNTIME_CHECKING
# define NDAY CHK
#endif

      DO j1 = 1,NDAY
        cosz1(j1)  = cosz(j1)
        sntz1(j1)  = f_one / cosz(j1)
        ssolar(j1) = s0fac * cosz(j1)

!  --- ...  surface albedo: bm,df - dir,dif;  1,2 - nir,uvv
        albbm(j1,1) = sfcalb(j1,1)
        albdf(j1,1) = sfcalb(j1,2)
        albbm(j1,2) = sfcalb(j1,3)
        albdf(j1,2) = sfcalb(j1,4)
      ENDDO

!  --- ...  prepare atmospheric profile for use in rrtm
!           the vertical index of internal array is from surface to top

        if (ivflip == 0) then       ! input from toa to sfc

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd

          do k = 1, nlay
            kk = nlp1 - k
      DO j1 = 1,NDAY
            pavel(j1,k) = plyr(j1,kk)
            tavel(j1,k) = tlyr(j1,kk)
            delp (j1,k) = plvl(j1,kk+1) - plvl(j1,kk)

!  --- ...  set absorber amount
!test use
!           h2ovmr(k)= max(f_zero,qlyr(j1,kk)*amdw)                     ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qlyr(j1,kk))                          ! input vol mixing ratio
!           o3vmr (k)= max(f_zero,olyr(j1,kk))                          ! input vol mixing ratio
!ncep model use
            h2ovmr(j1,k)=                                               &
     &                max(f_zero,qlyr(j1,kk)*amdw/(f_one-qlyr(j1,kk))) ! input specific humidity
            o3vmr (j1,k)=max(f_zero,olyr(j1,kk)*amdo3)                    ! input mass mixing ratio

            tem0 = (f_one - h2ovmr(j1,k))*con_amd + h2ovmr(j1,k)*con_amw
            coldry(j1,k) = tem2*delp(j1,k)/(tem1*tem0*(f_one            &
     &                      + h2ovmr(j1,k)))
            temcol(j1,k) = 1.0e-12 * coldry(j1,k)

            colamt(j1,k,1) = max(f_zero,    coldry(j1,k)*h2ovmr(j1,k))         ! h2o
            colamt(j1,k,2) =                                            &
     &                   max(temcol(j1,k), coldry(j1,k)*gasvmr(j1,kk,1))    ! co2
            colamt(j1,k,3) = max(f_zero,    coldry(j1,k)*o3vmr(j1,k))          ! o3
            colmol(j1,k)   = coldry(j1,k) + colamt(j1,k,1)
      ENDDO
          enddo

!  --- ...  set up gas column amount, convert from volume mixing ratio
!           to molec/cm2 based on coldry (scaled to 1.0e-20)

          if (iswrgas > 0) then
            do k = 1, nlay
              kk = nlp1 - k
      DO j1 = 1,NDAY
              colamt(j1,k,4)=                                           &
     &           max(temcol(j1,k),coldry(j1,k)*gasvmr(j1,kk,2))  ! n2o
              colamt(j1,k,5)=                                           &
     &           max(temcol(j1,k),coldry(j1,k)*gasvmr(j1,kk,3))  ! ch4
              colamt(j1,k,6)=                                           &
     &           max(temcol(j1,k),coldry(j1,k)*gasvmr(j1,kk,4))  ! o2
!             colamt(j1,k,7)=                                           &
!    &           max(temcol(j1,k),coldry(j1,k)*gasvmr(j1,kk,5))  ! co - notused
      ENDDO
            enddo
          else
            do k = 1, nlay
      DO j1 = 1,NDAY
              colamt(j1,k,4) = temcol(j1,k)                                  ! n2o
              colamt(j1,k,5) = temcol(j1,k)                                  ! ch4
              colamt(j1,k,6) = temcol(j1,k)                                  ! o2
!             colamt(j1,k,7) = temcol(j1,k)                                  ! co - notused
      ENDDO
            enddo
          endif

!  --- ...  set aerosol optical properties

          do k = 1, nlay
            kk = nlp1 - k
            do ib = 1, nbdsw
      DO j1 = 1,NDAY
              tauae(j1,k,ib) = aerosols(j1,kk,ib,1)
              ssaae(j1,k,ib) = aerosols(j1,kk,ib,2)
              asyae(j1,k,ib) = aerosols(j1,kk,ib,3)
      ENDDO
            enddo
          enddo

          if (iswcliq > 0) then    ! use prognostic cloud method
            do k = 1, nlay
              kk = nlp1 - k
      DO j1 = 1,NDAY
              cfrac(j1,k) = clouds(j1,kk,1)      ! cloud fraction
              cliqp(j1,k) = clouds(j1,kk,2)      ! cloud liq path
              reliq(j1,k) = clouds(j1,kk,3)      ! liq partical effctive radius
              cicep(j1,k) = clouds(j1,kk,4)      ! cloud ice path
              reice(j1,k) = clouds(j1,kk,5)      ! ice partical effctive radius
              cdat1(j1,k) = clouds(j1,kk,6)      ! cloud rain drop path
              cdat2(j1,k) = clouds(j1,kk,7)      ! rain partical effctive radius
              cdat3(j1,k) = clouds(j1,kk,8)      ! cloud snow path
              cdat4(j1,k) = clouds(j1,kk,9)      ! snow partical effctive radius
      ENDDO
            enddo
          else                     ! use diagnostic cloud method
            do k = 1, nlay
              kk = nlp1 - k
      DO j1 = 1,NDAY
              cfrac(j1,k) = clouds(j1,kk,1)      ! cloud fraction
              cdat1(j1,k) = clouds(j1,kk,2)      ! cloud optical depth
              cdat2(j1,k) = clouds(j1,kk,3)      ! cloud single scattering albedo
              cdat3(j1,k) = clouds(j1,kk,4)      ! cloud asymmetry factor
      ENDDO
            enddo
          endif                    ! end if_iswcliq

        else                        ! input from sfc to toa

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd

          do k = 1, nlay
      DO j1 = 1,NDAY
            pavel(j1,k) = plyr(j1,k)
            tavel(j1,k) = tlyr(j1,k)
            delp (j1,k) = plvl(j1,k) - plvl(j1,k+1)


!  --- ...  set absorber amount
!test use
!           h2ovmr(k)= max(f_zero,qlyr(j1,k)*amdw)                    ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qlyr(j1,k))                         ! input vol mixing ratio
!           o3vmr (k)= max(f_zero,olyr(j1,k))                         ! input vol mixing ratio
!ncep model use
            h2ovmr(j1,k)= max(f_zero,qlyr(j1,k)*amdw/(f_one-qlyr(j1,k))) ! input specific humidity
            o3vmr (j1,k)= max(f_zero,olyr(j1,k)*amdo3)                   ! input mass mixing ratio

            tem0=(f_one - h2ovmr(j1,k))*con_amd + h2ovmr(j1,k)*con_amw
            coldry(j1,k) = tem2 * delp(j1,k) /                          &
     &                              (tem1*tem0*(f_one + h2ovmr(j1,k)))
            temcol(j1,k) = 1.0e-12 * coldry(j1,k)

            colamt(j1,k,1) =                                            &
     &                  max(f_zero,    coldry(j1,k)*h2ovmr(j1,k))         ! h2o
            colamt(j1,k,2) =                                            &
     &                  max(temcol(j1,k), coldry(j1,k)*gasvmr(j1,k,1))    ! co2
            colamt(j1,k,3) =                                            &
     &                  max(f_zero,    coldry(j1,k)*o3vmr(j1,k))             ! o3
            colmol(j1,k)   = coldry(j1,k) + colamt(j1,k,1)
      ENDDO
          enddo


       if (lprnt) then
         write(0,*)' pavel=',pavel(1,:)
         write(0,*)' tavel=',tavel(1,:)
         write(0,*)' delp=',delp(1,:)
         write(0,*)' h2ovmr=',h2ovmr(1,:)*1000
         write(0,*)' o3vmr=',o3vmr(1,:)*1000000
       endif

!  --- ...  set up gas column amount, convert from volume mixing ratio
!           to molec/cm2 based on coldry (scaled to 1.0e-20)

          if (iswrgas > 0) then
            do k = 1, nlay
      DO j1 = 1,NDAY
              colamt(j1,k,4) =                                          &
     &            max(temcol(j1,k), coldry(j1,k)*gasvmr(j1,k,2))  ! n2o
              colamt(j1,k,5) =                                          &
     &            max(temcol(j1,k), coldry(j1,k)*gasvmr(j1,k,3))  ! ch4
              colamt(j1,k,6) =                                          &
     &            max(temcol(j1,k), coldry(j1,k)*gasvmr(j1,k,4))  ! o2
!             colamt(j1,k,7) =                                          &
!    &            max(temcol(j1,k), coldry(j1,k)*gasvmr(j1,k,5))  ! co - notused
            enddo
      ENDDO
          else
            do k = 1, nlay
      DO j1 = 1,NDAY
              colamt(j1,k,4) = temcol(j1,k)                                 ! n2o
              colamt(j1,k,5) = temcol(j1,k)                                 ! ch4
              colamt(j1,k,6) = temcol(j1,k)                                 ! o2
!             colamt(j1,k,7) = temcol(j1,k)                                 ! co - notused
      ENDDO
            enddo
          endif

!  --- ...  set aerosol optical properties

          do ib = 1, nbdsw
            do k = 1, nlay
      DO j1 = 1,NDAY
              tauae(j1,k,ib) = aerosols(j1,k,ib,1)
              ssaae(j1,k,ib) = aerosols(j1,k,ib,2)
              asyae(j1,k,ib) = aerosols(j1,k,ib,3)
      ENDDO
            enddo
          enddo

          if (iswcliq > 0) then    ! use prognostic cloud method
            do k = 1, nlay
      DO j1 = 1,NDAY
              cfrac(j1,k) = clouds(j1,k,1)       ! cloud fraction
              cliqp(j1,k) = clouds(j1,k,2)       ! cloud liq path
              reliq(j1,k) = clouds(j1,k,3)       ! liq partical effctive radius
              cicep(j1,k) = clouds(j1,k,4)       ! cloud ice path
              reice(j1,k) = clouds(j1,k,5)       ! ice partical effctive radius
              cdat1(j1,k) = clouds(j1,k,6)       ! cloud rain drop path
              cdat2(j1,k) = clouds(j1,k,7)       ! rain partical effctive radius
              cdat3(j1,k) = clouds(j1,k,8)       ! cloud snow path
              cdat4(j1,k) = clouds(j1,k,9)       ! snow partical effctive radius
      ENDDO
            enddo
          else                     ! use diagnostic cloud method
            do k = 1, nlay
      DO j1 = 1,NDAY
              cfrac(j1,k) = clouds(j1,k,1)       ! cloud fraction
              cdat1(j1,k) = clouds(j1,k,2)       ! cloud optical depth
              cdat2(j1,k) = clouds(j1,k,3)       ! cloud single scattering albedo
              cdat3(j1,k) = clouds(j1,k,4)       ! cloud asymmetry factor
      ENDDO
            enddo
          endif                    ! end if_iswcliq

        endif                       ! if_ivflip

!  --- ...  compute fractions of clear sky view

        zcf0   = f_one
        zcf1   = f_one

        if (iovrsw == 0) then                    ! random overlapping
          do k = 1, nlay
      DO j1 = 1,NDAY
            zcf0(j1) = zcf0(j1) * (f_one - cfrac(j1,k))
      ENDDO
          enddo
        else if (iovrsw == 1) then               ! max/ran overlapping
          do k = 1, nlay
       if (buggal) write(0,*)'cfrac ',k,cfrac(1,k)
      DO j1 = 1,NDAY
            if (cfrac(j1,k) > ftiny) then                ! cloudy layer
              zcf1(j1) = min ( zcf1(j1), f_one-cfrac(j1,k) )
            elseif (zcf1(j1) < f_one) then                ! clear layer
              zcf0(j1) = zcf0(j1) * zcf1(j1)
              zcf1(j1) = f_one
            endif
      ENDDO
          enddo
          zcf0 = zcf0 * zcf1
        else if (iovrsw == 2) then               ! maximum overlapping
          do k = 1, nlay
      DO j1 = 1,NDAY
            zcf0(j1) = min ( zcf0(j1), f_one-cfrac(j1,k) )
      ENDDO
          enddo
        endif

        WHERE (zcf0 <= ftiny) zcf0 = f_zero
        WHERE (zcf0 > oneminus) zcf0 = f_one
        zcf1 = f_one - zcf0

!  --- ...  compute cloud optical properties

      DO j1 = 1,NDAY
        if (zcf1(j1) > f_zero) then     ! cloudy sky column
! the copy and transpose implied here is not good but will fix later.... JM
          call cldprop                                                  &
!  ---  inputs:
     &     ( cfrac(j1,:),cliqp(j1,:),reliq(j1,:),cicep(j1,:),           &
     &       reice(j1,:),cdat1(j1,:),cdat2(j1,:),cdat3(j1,:),           &
     &       cdat4(j1,:),                                               &
     &       zcf1(j1), nlay, ipseed(j1),                                &
!  ---  outputs:
     &       taucw(j1,:,:), ssacw(j1,:,:), asycw(j1,:,:),               &
     &       cldfrc(j1,:), cldfmc(j1,:,:)                               &
     &     )
        else                        ! clear sky column
          cldfrc(j1,:)  = f_zero
          cldfmc(j1,:,:)= f_zero
          taucw(j1,:,:) = f_zero
          ssacw(j1,:,:) = f_zero
          asycw(j1,:,:) = f_zero
        endif   ! end if_zcf1_block

        call setcoef                                                    &
!  ---  inputs:
     &     ( pavel(j1,:),tavel(j1,:),h2ovmr(j1,:), nlay,nlp1,           &
!  ---  outputs:
     &       laytrop(j1),jp(j1,:),jt(j1,:),jt1(j1,:),                   &
     &       fac00(j1,:),fac01(j1,:),fac10(j1,:),fac11(j1,:),           &
     &       selffac(j1,:),selffrac(j1,:),indself(j1,:),                &
     &       forfac(j1,:),forfrac(j1,:),indfor(j1,:)                    &
     &     )

      ENDDO

!  --- ...  calculate optical depths for gaseous absorption and Rayleigh
!           scattering

        call taumol                                                     &
!  ---  inputs:
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself,NLAY,nday,  &
!  ---  outputs:
     &       sfluxzen, taug, taur                                       &
     &     )

!  --- ...  call the 2-stream radiation transfer model

        if ( isubcsw <= 0 ) then     ! use standard cloud scheme
          call spcvrtc                                                  &
!  ---  inputs:
     &     ( ssolar,cosz1,sntz1,albbm,albdf,sfluxzen,cldfrc,            &
     &       zcf1,zcf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,   &
     &       nlay, nlp1,nday,                                           &
!  ---  outputs:
     &       fxupc,fxdnc,fxup0,fxdn0,                                   &
     &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
     &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
     &     )

#ifdef TODO
        else                         ! use mcica cloud scheme

!rv       call spcvrtm                                                  &
!  ---  inputs:
!rv  &     ( ssolar,cosz1,sntz1,albbm,albdf,sfluxzen,cldfmc,            &
!rv  &       zcf1,zcf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,   &
!rv  &       nlay, nlp1,                                                &
!  ---  outputs:
!rv  &       fxupc,fxdnc,fxup0,fxdn0,                                   &
!rv  &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
!rv  &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
!rv  &     )
#endif

        endif

!  --- ...  sum up total spectral fluxes for total-sky

        do k = 1, nlp1
          flxuc(1:NDAY,k) = f_zero
          flxdc(1:NDAY,k) = f_zero

          do ib = 1, nbdsw
            flxuc(1:NDAY,k) = flxuc(1:NDAY,k) + fxupc(1:NDAY,k,ib)
            flxdc(1:NDAY,k) = flxdc(1:NDAY,k) + fxdnc(1:NDAY,k,ib)
          enddo
        enddo

!! --- ...  optional clear sky fluxes

        if ( lhsw0 .or. lflxprf ) then
          do k = 1, nlp1
            flxu0(1:NDAY,k) = f_zero
            flxd0(1:NDAY,k) = f_zero

            do ib = 1, nbdsw
              flxu0(1:NDAY,k) = flxu0(1:NDAY,k) + fxup0(1:NDAY,k,ib)
              flxd0(1:NDAY,k) = flxd0(1:NDAY,k) + fxdn0(1:NDAY,k,ib)
            enddo
          enddo
        endif

!  --- ...  prepare for final outputs

        do k = 1, nlay
          rfdelp(1:NDAY,k) = heatfac / delp(1:NDAY,k)
        enddo

        if ( lfdncmp ) then
!! --- ...  optional uv-b surface downward flux
      DO j1 = 1,NDAY
          fdncmp(j1)%uvbf0 = suvbf0(j1)
          fdncmp(j1)%uvbfc = suvbfc(j1)

!! --- ...  optional beam and diffuse sfc fluxes
          fdncmp(j1)%nirbm = sfbmc(j1,1)
          fdncmp(j1)%nirdf = sfdfc(j1,1)
          fdncmp(j1)%visbm = sfbmc(j1,2)
          fdncmp(j1)%visdf = sfdfc(j1,2)
       ENDDO
        endif    ! end if_lfdncmp

!  --- ...  toa and sfc fluxes

      DO j1 = 1,NDAY
        topflx(j1)%upfxc = ftoauc(j1)
        topflx(j1)%dnfxc = ftoadc(j1)
        topflx(j1)%upfx0 = ftoau0(j1)

        sfcflx(j1)%upfxc = fsfcuc(j1)
        sfcflx(j1)%dnfxc = fsfcdc(j1)
        sfcflx(j1)%upfx0 = fsfcu0(j1)
        sfcflx(j1)%dnfx0 = fsfcd0(j1)
      ENDDO

        if (ivflip == 0) then       ! output from toa to sfc

!  --- ...  compute heating rates

          fnet(1:NDAY,1) = flxdc(1:NDAY,1) - flxuc(1:NDAY,1)

          do k = 2, nlp1
            kk = nlp1 - k + 1
!DEC$ VECTOR ALWAYS
      DO j1 = 1,NDAY
            fnet(j1,k) = flxdc(j1,k) - flxuc(j1,k)
            hswc(j1,kk) = (fnet(j1,k)-fnet(j1,k-1)) * rfdelp(j1,k-1)
      ENDDO
          enddo

!! --- ...  optional flux profiles

          if ( lflxprf ) then
            do k = 1, nlp1
              kk = nlp1 - k + 1
!DEC$ VECTOR ALWAYS
      DO j1 = 1,NDAY
              flxprf(j1,kk)%upfxc = flxuc(j1,k)
              flxprf(j1,kk)%dnfxc = flxdc(j1,k)
              flxprf(j1,kk)%upfx0 = flxu0(j1,k)
              flxprf(j1,kk)%dnfx0 = flxd0(j1,k)
      ENDDO
            enddo
          endif

!! --- ...  optional clear sky heating rates

          if ( lhsw0 ) then
            fnet(1:NDAY,1) = flxd0(1:NDAY,1) - flxu0(1:NDAY,1)

            do k = 2, nlp1
              kk = nlp1 - k + 1
!DEC$ VECTOR ALWAYS
      DO j1 = 1,NDAY
              fnet(j1,k) = flxd0(j1,k) - flxu0(j1,k)
              hsw0(j1,kk) = (fnet(j1,k)-fnet(j1,k-1)) * rfdelp(j1,k-1)
      ENDDO
            enddo
          endif

!! --- ...  optional spectral band heating rates

          if ( lhswb ) then
            do mb = 1, nbdsw
              fnet(1:NDAY,1) = fxdnc(1:NDAY,1,mb) - fxupc(1:NDAY,1,mb)

              do k = 2, nlp1
                kk = nlp1 - k + 1
!DEC$ VECTOR ALWAYS
      DO j1 = 1,NDAY
                fnet(j1,k) = fxdnc(j1,k,mb) - fxupc(j1,k,mb)
                hswb(j1,kk,mb)=(fnet(j1,k)-fnet(j1,k-1))*rfdelp(j1,k-1)
      ENDDO
              enddo
            enddo
          endif

        else                        ! output from sfc to toa

!  --- ...  compute heating rates

          fnet(1:NDAY,1) = flxdc(1:NDAY,1) - flxuc(1:NDAY,1)

          do k = 2, nlp1
            fnet(1:NDAY,k) = flxdc(1:NDAY,k) - flxuc(1:NDAY,k)
            hswc(1:NDAY,k-1) =                                          &
     &             (fnet(1:NDAY,k)-fnet(1:NDAY,k-1))*rfdelp(1:NDAY,k-1)
          enddo

!! --- ...  optional flux profiles

          if ( lflxprf ) then
            do k = 1, nlp1
!DEC$ VECTOR ALWAYS
      DO j1 = 1,NDAY
              flxprf(j1,k)%upfxc = flxuc(j1,k)
              flxprf(j1,k)%dnfxc = flxdc(j1,k)
              flxprf(j1,k)%upfx0 = flxu0(j1,k)
              flxprf(j1,k)%dnfx0 = flxd0(j1,k)
      ENDDO
            enddo
          endif

!! --- ...  optional clear sky heating rates

          if ( lhsw0 ) then
            fnet(1:NDAY,1) = flxd0(1:NDAY,1) - flxu0(1:NDAY,1)

            do k = 2, nlp1
              fnet(1:NDAY,k) = flxd0(1:NDAY,k) - flxu0(1:NDAY,k)
              hsw0(1:NDAY,k-1) =                                        &
     &            (fnet(1:NDAY,k)-fnet(1:NDAY,k-1)) * rfdelp(1:NDAY,k-1)
            enddo
          endif

!! --- ...  optional spectral band heating rates

          if ( lhswb ) then
            do mb = 1, nbdsw
              fnet(1:NDAY,1) = fxdnc(1:NDAY,1,mb) - fxupc(1:NDAY,1,mb)

              do k = 1, nlay
                fnet(1:NDAY,k+1) =                                      &
     &              fxdnc(1:NDAY,k+1,mb) - fxupc(1:NDAY,k+1,mb)
                hswb(1:NDAY,k,mb) =                                     &
     &              (fnet(1:NDAY,k+1)-fnet(1:NDAY,k)) * rfdelp(1:NDAY,k)
              enddo
            enddo
          endif

        endif                       ! if_ivflip

!jm      enddo   lab_do_ipt

      return
!...................................
      end subroutine swrad
!-----------------------------------


!-----------------------------------
      subroutine rswinit                                                &
!...................................

!  ---  inputs:
     &     ( me )
!  ---  outputs: (none)

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  initialize non-varying module variables, conversion factors,!
! and look-up tables.                                                   !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!    me       - print control for parallel process                      !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  external module variables:  (in physparam)                           !
!   iswrate - heating rate unit selections                              !
!           =1: output in k/day                                         !
!           =2: output in k/second                                      !
!   iswrgas - control flag for rare gases (ch4,n2o,o2, etc.)            !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   iswcliq - liquid cloud optical properties contrl flag               !
!           =0: input cloud opt depth from diagnostic scheme            !
!           >0: input cwp,rew, and other cloud content parameters       !
!   isubcsw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   icldflg - cloud scheme control flag                                 !
!           =0: diagnostic scheme gives cloud tau, omiga, and g.        !
!           =1: prognostic scheme gives cloud liq/ice path, etc.        !
!   iovrsw  - clouds vertical overlapping control flag                  !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud                                   !
!   iswmode - control flag for 2-stream transfer scheme                 !
!           =1; delta-eddington    (joseph et al., 1976)                !
!           =2: pifm               (zdunkowski et al., 1980)            !
!           =3: discrete ordinates (liou, 1973)                         !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
! definitions:                                                          !
!     arrays for 10000-point look-up tables:                            !
!     tau_tbl  clear-sky optical depth                                  !
!     exp_tbl  exponential lookup table for transmittance               !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: me

!  ---  outputs: none

!  ---  locals:
      real (kind=kind_phys), parameter :: expeps = 1.e-20

      integer :: i

      real (kind=kind_phys) :: tfn, tau

!
!===> ... begin here
!
      if ( iovrsw<0 .or. iovrsw>2 ) then
        print *,'  *** Error in specification of cloud overlap flag',   &
     &          ' IOVRSW=',iovrsw,' in RSWINIT !!'
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif

      if (me == 0) then
        print *,' - Using AER Shortwave Radiation, Version: ',VTAGSW

        if (iswmode == 1) then
          print *,'   --- Delta-eddington 2-stream transfer scheme'
        else if (iswmode == 2) then
          print *,'   --- PIFM 2-stream transfer scheme'
        else if (iswmode == 3) then
          print *,'   --- Discrete ordinates 2-stream transfer scheme'
        endif

        if (iswrgas <= 0) then
          print *,'   --- Rare gases absorption is NOT included in SW'
        else
          print *,'   --- Include rare gases N2O, CH4, O2, absorptions',&
     &            ' in SW'
        endif

        if ( isubcsw == 0 ) then
          print *,'   --- Using standard grid average clouds, no ',     &
     &            'sub-column clouds approximation applied'
        elseif ( isubcsw == 1 ) then
          print *,'   --- Using MCICA sub-colum clouds approximation ', &
     &            'with a prescribed sequence of permutation seeds'
        elseif ( isubcsw == 2 ) then
          print *,'   --- Using MCICA sub-colum clouds approximation ', &
     &            'with provided input array of permutation seeds'
        else
          print *,'  *** Error in specification of sub-column cloud ',  &
     &            ' control flag isubcsw =',isubcsw,' !!'
          ! stop
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif
      endif

!  --- ...  check cloud flags for consistency

      if ((icldflg == 0 .and. iswcliq /= 0) .or.                        &
     &    (icldflg == 1 .and. iswcliq == 0)) then
        print *,'  *** Model cloud scheme inconsistent with SW',        &
     &          ' radiation cloud radiative property setup !!'
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif

!  --- ...  setup constant factors for heating rate
!           the 1.0e-2 is to convert pressure from mb to N/m**2

      if (iswrate == 1) then
!       heatfac = 8.4391
!       heatfac = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
        heatfac = con_g * 864.0 / con_cp            !   (in k/day)
      else
        heatfac = con_g * 1.0e-2 / con_cp           !   (in k/second)
      endif

!  --- ...  define exponential lookup tables for transmittance. tau is
!           computed as a function of the tau transition function, and
!           transmittance is calculated as a function of tau.  all tables
!           are computed at intervals of 0.0001.  the inverse of the
!           constant used in the Pade approximation to the tau transition
!           function is set to bpade.

      exp_tbl(0) = 1.0
      exp_tbl(NTBMX) = expeps

      do i = 1, NTBMX-1
        tfn = float(i) / float(NTBMX-i)
        tau = bpade * tfn
        exp_tbl(i) = exp( -tau )
      enddo

      return
!...................................
      end subroutine rswinit
!-----------------------------------


!-----------------------------------
      subroutine cldprop                                                &
!...................................
!  ---  inputs:
     &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     &
     &       cf1, nlay, ipseed,                                         &
!  ---  output:
     &       taucw, ssacw, asycw, cldfrc, cldfmc                        &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! Purpose: Compute the cloud optical properties for each cloudy layer   !
! and g-point interval.                                                 !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!    cfrac - real, layer cloud fraction                            nlay !
!        .....  for  iswcliq > 0 (prognostic cloud sckeme)  - - -       !
!    cliqp - real, layer in-cloud liq water path (g/m**2)          nlay !
!    reliq - real, mean eff radius for liq cloud (micron)          nlay !
!    cicep - real, layer in-cloud ice water path (g/m**2)          nlay !
!    reice - real, mean eff radius for ice cloud (micron)          nlay !
!    cdat1 - real, layer rain drop water path (g/m**2)             nlay !
!    cdat2 - real, effective radius for rain drop (micron)         nlay !
!    cdat3 - real, layer snow flake water path(g/m**2)             nlay !
!    cdat4 - real, mean eff radius for snow flake(micron)          nlay !
!        .....  for iswcliq = 0  (diagnostic cloud sckeme)  - - -       !
!    cdat1 - real, layer cloud optical depth                       nlay !
!    cdat2 - real, layer cloud single scattering albedo            nlay !
!    cdat3 - real, layer cloud asymmetry factor                    nlay !
!    cdat4 - real, optional use                                    nlay !
!    cliqp - real, not used                                        nlay !
!    cicep - real, not used                                        nlay !
!    reliq - real, not used                                        nlay !
!    reice - real, not used                                        nlay !
!                                                                       !
!    cf1   - real, effective total cloud cover at surface           1   !
!    nlay  - integer, vertical layer number                         1   !
!    ipseed- permutation seed for generating random numbers (isubcsw>0) !
!                                                                       !
!  outputs:                                                             !
!    taucw  - real, cloud optical depth, w/o delta scaled    nlay*nbdsw !
!    ssacw  - real, weighted cloud single scattering albedo  nlay*nbdsw !
!                             (ssa = ssacw / taucw)                     !
!    asycw  - real, weighted cloud asymmetry factor          nlay*nbdsw !
!                             (asy = asycw / ssacw)                     !
!    cldfrc - real, cloud fraction of grid mean value              nlay !
!    cldfmc - real, cloud fraction for each sub-column       nlay*ngptsw!
!                                                                       !
!                                                                       !
!  explanation of the method for each value of iswcliq, and iswcice.    !
!  set up in module "physparam"                                         !
!                                                                       !
!     iswcliq=0  : input cloud optical property (tau, ssa, asy).        !
!                  (used for diagnostic cloud method)                   !
!     iswcliq>0  : input cloud liq/ice path and effective radius, also  !
!                  require the user of 'iswcice' to specify the method  !
!                  used to compute aborption due to water/ice parts.    !
!  ...................................................................  !
!                                                                       !
!     iswcliq=1  : liquid water cloud optical properties are computed   !
!                  as in hu and stamnes (1993), j. clim., 6, 728-742.   !
!                                                                       !
!     iswcice used only when iswcliq > 0                                !
!                  the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are inputs.                                !
!     iswcice=1  : ice cloud optical properties are computed as in      !
!                  ebert and curry (1992), jgr, 97, 3831-3836.          !
!     iswcice=2  : ice cloud optical properties are computed as in      !
!                  streamer v3.0 (2001), key, streamer user's guide,    !
!                  cooperative institude for meteorological studies,95pp!
!     iswcice=3  : ice cloud optical properties are computed as in      !
!                  fu (1996), j. clim., 9.                              !
!                                                                       !
!  other cloud control module variables:                                !
!     isubcsw =0: standard cloud scheme, no sub-col cloud approximation !
!             >0: mcica sub-col cloud scheme using ipseed as permutation!
!                 seed for generating rundom numbers                    !
!                                                                       !
!  ======================  end of description block  =================  !
!
      use module_radsw_cldprtb

!  ---  inputs:
      integer, intent(in) :: nlay, ipseed
      real (kind=kind_phys), intent(in) :: cf1

      real (kind=kind_phys), dimension(nlay), intent(in) :: cliqp,      &
     &       reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4, cfrac

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay,ngptsw), intent(out) ::     &
     &       cldfmc
      real (kind=kind_phys), dimension(nlay,nbdsw),  intent(out) ::     &
     &       taucw, ssacw, asycw
      real (kind=kind_phys), dimension(nlay), intent(out) :: cldfrc

!  ---  locals:
      real (kind=kind_phys), dimension(nblow:nbhgh) :: tauliq, tauice,  &
     &       ssaliq, ssaice, ssaran, ssasnw, asyliq, asyice,            &
     &       asyran, asysnw, tausnow
      real (kind=kind_phys), dimension(nlay)       :: cldf

      real (kind=kind_phys) :: dgeice, factor, fint, tauran, tausnw,    &
     &       cldliq, refliq, cldice, refice, cldran, cldsnw, refsnw,    &
     &       extcoliq, ssacoliq, asycoliq, extcoice, ssacoice, asycoice,&
     &       extcosno, ssacosno, asycosno,                              &
     &       dgesnw, dummy_tau

      logical :: lcloudy(nlay,ngptsw)
      integer :: ia, ib, ig, jb, k, index

!
!===> ...  begin here
!
      do ib = 1, nbdsw
        do k = 1, nlay
          taucw (k,ib) = f_zero
          ssacw (k,ib) = f_one
          asycw (k,ib) = f_zero
        enddo
      enddo

!  --- ...  compute cloud radiative properties for a cloudy column

      lab_if_iswcliq : if (iswcliq > 0) then

        lab_do_k : do k = 1, nlay
          lab_if_cld : if (cfrac(k) > ftiny) then

!  --- ...  optical properties for rain and snow
            cldran = cdat1(k)
!           refran = cdat2(k)
            cldsnw = cdat3(k)
            refsnw = cdat4(k)
            dgesnw = 1.0315 * refsnw        ! for fu's snow formula

            if (iswcice < 4) then

              tauran = cldran * a0r
!  ---  if use fu's formula it needs to be normalized by snow/ice density
!       !not use snow density = 0.1 g/cm**3 = 0.1 g/(mu * m**2)
!       use ice density = 0.9167 g/cm**3 = 0.9167 g/(mu * m**2)
!       1/0.9167 = 1.09087
!       factor 1.5396=8/(3*sqrt(3)) converts reff to generalized ice particle size
!       use newer factor value 1.0315
              if (cldsnw>f_zero .and. refsnw>10.0_kind_phys) then
!               tausnw = cldsnw * (a0s + a1s/refsnw)
                tausnw = cldsnw*1.09087*(a0s + a1s/dgesnw)     ! fu's formula
              else
                tausnw = f_zero
              endif

              do ib = nblow, nbhgh
                ssaran(ib) = tauran * (f_one - b0r(ib))
                ssasnw(ib) = tausnw * (f_one - (b0s(ib)+b1s(ib)*dgesnw))
                asyran(ib) = ssaran(ib) * c0r(ib)
                asysnw(ib) = ssasnw(ib) * c0s(ib)
              enddo

            endif

            cldliq = cliqp(k)
            cldice = cicep(k)
            refliq = reliq(k)
            refice = reice(k)

!  --- ...  calculation of absorption coefficients due to water clouds.

            if ( cldliq <= f_zero ) then
              do ib = nblow, nbhgh
                tauliq(ib) = f_zero
                ssaliq(ib) = f_zero
                asyliq(ib) = f_zero
              enddo
            else
              if ( iswcliq == 1 ) then
                factor = refliq - 1.5
                index  = max( 1, min( 57, int( factor ) ))
                fint   = factor - float(index)

                do ib = nblow, nbhgh
                  extcoliq = max(f_zero,            extliq1(index,ib)   &
     &              + fint*(extliq1(index+1,ib)-extliq1(index,ib)) )
                  ssacoliq = max(f_zero, min(f_one, ssaliq1(index,ib)   &
     &              + fint*(ssaliq1(index+1,ib)-ssaliq1(index,ib)) ))

                  asycoliq = max(f_zero, min(f_one, asyliq1(index,ib)   &
     &              + fint*(asyliq1(index+1,ib)-asyliq1(index,ib)) ))
!                 forcoliq = asycoliq * asycoliq

                  tauliq(ib) = cldliq     * extcoliq
                  ssaliq(ib) = tauliq(ib) * ssacoliq
                  asyliq(ib) = ssaliq(ib) * asycoliq
                enddo
              endif   ! end if_iswcliq_block
            endif   ! end if_cldliq_block

!  --- ...  calculation of absorption coefficients due to ice clouds.

            if ( cldice <= f_zero ) then
              do ib = nblow, nbhgh
                tauice(ib) = f_zero
                ssaice(ib) = f_zero
                asyice(ib) = f_zero
              enddo
            else

!  --- ...  ebert and curry approach for all particle sizes though somewhat
!           unjustified for large ice particles

              if ( iswcice == 1 ) then
                refice = min(130.0_kind_phys,max(13.0_kind_phys,refice))

                do ib = nblow, nbhgh
                  ia = idxebc(ib)           ! eb_&_c band index for ice cloud coeff

                  extcoice = max(f_zero, abari(ia)+bbari(ia)/refice )
                  ssacoice = max(f_zero, min(f_one,                     &
     &                             f_one-cbari(ia)-dbari(ia)*refice ))
                  asycoice = max(f_zero, min(f_one,                     &
     &                                   ebari(ia)+fbari(ia)*refice ))
!                 forcoice = asycoice * asycoice

                  tauice(ib) = cldice     * extcoice
                  ssaice(ib) = tauice(ib) * ssacoice
                  asyice(ib) = ssaice(ib) * asycoice
                enddo

!  --- ...  streamer approach for ice effective radius between 5.0 and 131.0 microns

              elseif ( iswcice == 2 ) then
                refice = min(131.0_kind_phys,max(5.0_kind_phys,refice))

                factor = (refice - 2.0) / 3.0
                index  = max( 1, min( 42, int( factor ) ))
                fint   = factor - float(index)

                do ib = nblow, nbhgh
                  extcoice = max(f_zero,            extice2(index,ib)   &
     &                + fint*(extice2(index+1,ib)-extice2(index,ib)) )
                  ssacoice = max(f_zero, min(f_one, ssaice2(index,ib)   &
     &                + fint*(ssaice2(index+1,ib)-ssaice2(index,ib)) ))
                  asycoice = max(f_zero, min(f_one, asyice2(index,ib)   &
     &                + fint*(asyice2(index+1,ib)-asyice2(index,ib)) ))
!                 forcoice = asycoice * asycoice

                  tauice(ib) = cldice     * extcoice
                  ssaice(ib) = tauice(ib) * ssacoice
                  asyice(ib) = ssaice(ib) * asycoice
                enddo

!  --- ...  fu's approach for ice effective radius between 4.8 and 135 microns
!           (generalized effective size from 5 to 140 microns)

              elseif ( iswcice == 3 ) then
                dgeice = max( 5.0, min( 140.0, 1.0315*refice ))

                factor = (dgeice - 2.0) / 3.0
                index  = max( 1, min( 45, int( factor ) ))
                fint   = factor - float(index)

                do ib = nblow, nbhgh
                  extcoice = max(f_zero,            extice3(index,ib)   &
     &                + fint*(extice3(index+1,ib)-extice3(index,ib)) )
                  ssacoice = max(f_zero, min(f_one, ssaice3(index,ib)   &
     &                + fint*(ssaice3(index+1,ib)-ssaice3(index,ib)) ))
                  asycoice = max(f_zero, min(f_one, asyice3(index,ib)   &
     &                + fint*(asyice3(index+1,ib)-asyice3(index,ib)) ))
!                 fdelta   = max(f_zero, min(f_one, fdlice3(index,ib)   &
!    &                + fint*(fdlice3(index+1,ib)-fdlice3(index,ib)) ))
!                 forcoice = min( asycoice, fdelta+0.5/ssacoice )           ! see fu 1996 p. 2067

                  tauice(ib) = cldice     * extcoice
                  ssaice(ib) = tauice(ib) * ssacoice
                  asyice(ib) = ssaice(ib) * asycoice
                enddo

              elseif ( iswcice == 4 ) then
                dgeice = max( 10.1, min( 129.9, refice ))

                factor = (dgeice - 2.0) / 3.0
                index  = max( 1, min( 45, int( factor ) ))
                fint   = factor - float(index)

                do ib = nblow, nbhgh
                  extcoice = max(f_zero,            extice3(index,ib)   &
     &                + fint*(extice3(index+1,ib)-extice3(index,ib)) )
                  ssacoice = max(f_zero, min(f_one, ssaice3(index,ib)   &
     &                + fint*(ssaice3(index+1,ib)-ssaice3(index,ib)) ))
                  asycoice = max(f_zero, min(f_one, asyice3(index,ib)   &
     &                + fint*(asyice3(index+1,ib)-asyice3(index,ib)) ))

                  tauice(ib) = cldice     * extcoice
                  ssaice(ib) = tauice(ib) * ssacoice
                  asyice(ib) = ssaice(ib) * asycoice
                enddo

              endif   ! end if_iswcice_block
            endif   ! end if_cldice_block

!  ===================================================================
!  --- Extra from Greg

!+---+-----------------------------------------------------------------+
!..Calculation of coefficients due to snow, new if iflagice=4.
!.. Treated entirely same as cloud ice in effective radius and band
!.. mode, and combined later as a linear sum.  G. Thompson
!+---+-----------------------------------------------------------------+

            if (iswcice == 4) then
              if (cldsnw <= f_zero) then
                do ib = nblow, nbhgh
                  tausnow(ib) = f_zero
                  ssasnw (ib) = f_zero
                  asysnw (ib) = f_zero
                enddo
              else
                dgesnw = max( 25.1, min( 129.9, refsnw ))

                factor = (dgesnw - 2.0) / 3.0
                index  = max( 1, min( 45, int( factor ) ))
                fint   = factor - float(index)

                do ib = nblow, nbhgh
                  extcosno = max(f_zero,            extice3(index,ib)   &
     &                + fint*(extice3(index+1,ib)-extice3(index,ib)) )
                  ssacosno = max(f_zero, min(f_one, ssaice3(index,ib)   &
     &                + fint*(ssaice3(index+1,ib)-ssaice3(index,ib)) ))
                  asycosno = max(f_zero, min(f_one, asyice3(index,ib)   &
     &                + fint*(asyice3(index+1,ib)-asyice3(index,ib)) ))

                  tausnow(ib) = MIN(35., cldsnw * extcosno)
                  ssasnw (ib) = tausnow(ib) * ssacosno
                  asysnw (ib) = ssasnw (ib) * asycosno
                enddo

              endif
            endif

            if (iswcice==4 ) then

              do ib = 1, nbdsw
                jb = nblow + ib - 1
                tauliq(jb) = MIN(35., tauliq(jb))
                tauice(jb) = MIN(35., tauice(jb))
                dummy_tau = MIN(tauliq(jb)+tauice(jb)+tausnow(jb),35.)
                taucw(k,ib)=MAX(1.0E-3,dummy_tau)
                ssacw(k,ib) = ssaliq(jb) + ssaice(jb) + ssasnw(jb)
                asycw(k,ib) = asyliq(jb) + asyice(jb) + asysnw(jb)
              enddo

!  ===================================================================

            else

              do ib = 1, nbdsw
                jb = nblow + ib - 1
                taucw(k,ib) = tauliq(jb)+tauice(jb)+tauran+tausnw
                ssacw(k,ib) = ssaliq(jb)+ssaice(jb)+ssaran(jb)+ssasnw(jb)
                asycw(k,ib) = asyliq(jb)+asyice(jb)+asyran(jb)+asysnw(jb)
              enddo

            endif

          endif  lab_if_cld
        enddo  lab_do_k

      else  lab_if_iswcliq

        do k = 1, nlay
          if (cfrac(k) > ftiny) then
            do ib = 1, nbdsw
              taucw(k,ib) = cdat1(k)
              ssacw(k,ib) = cdat1(k)    * cdat2(k)
              asycw(k,ib) = ssacw(k,ib) * cdat3(k)
            enddo
          endif
        enddo

      endif  lab_if_iswcliq

!  ---  distribute cloud properties to each g-point

      if ( isubcsw > 0 ) then      ! mcica sub-col clouds approx

        cldf(:) = cfrac(:)
        where (cldf(:) < ftiny)
          cldf(:) = f_zero
        end where

!  --- ...  call sub-column cloud generator

        call mcica_subcol                                               &
!  ---  inputs:
     &     ( cldf, nlay, ipseed,                                        &
!  ---  outputs:
     &       lcloudy                                                    &
     &     )

        do ig = 1, ngptsw
          do k = 1, nlay
            if ( lcloudy(k,ig) ) then
              cldfmc(k,ig) = f_one
            else
              cldfmc(k,ig) = f_zero
            endif
          enddo
        enddo

      else                         ! non-mcica, normalize cloud

        do k = 1, nlay
          cldfrc(k) = cfrac(k) / cf1
        enddo
      endif   ! end if_isubcsw_block

      return
!...................................
      end subroutine cldprop
!-----------------------------------


! ----------------------------------
      subroutine mcica_subcol                                           &
! ..................................
!  ---  inputs:
     &    ( cldf, nlay, ipseed,                                         &
!  ---  outputs:
     &      lcloudy                                                     &
     &    )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                size !
!   cldf    - real, layer cloud fraction                           nlay !
!   nlay    - integer, number of model vertical layers               1  !
!   ipseed  - integer, permute seed for random num generator         1  !
!    ** note : if the cloud generator is called multiple times, need    !
!              to permute the seed between each call; if between calls  !
!              for lw and sw, use values differ by the number of g-pts. !
!                                                                       !
!  output variables:                                                    !
!   lcloudy - logical, sub-colum cloud profile flag array    nlay*ngptsw!
!                                                                       !
!  other control flags from module variables:                           !
!     iovrsw    : control flag for cloud overlapping method             !
!                 =0:random; =1:maximum/random; =2:maximum              !
!                                                                       !
!                                                                       !
!  =====================    end of definitions    ====================  !

      implicit none

!  ---  inputs:
      integer, intent(in) :: nlay, ipseed

      real (kind=kind_phys), dimension(nlay), intent(in) :: cldf

!  ---  outputs:
      logical, dimension(nlay,ngptsw), intent(out):: lcloudy

!  ---  locals:
      real (kind=kind_phys) :: cdfunc(nlay,ngptsw), tem1,               &
     &       rand2d(nlay*ngptsw), rand1d(ngptsw)

      type (random_stat) :: stat          ! for thread safe random generator

      integer :: k, n, k1
!
!===> ...  begin here
!
!  --- ...  advance randum number generator by ipseed values

      call random_setseed                                               &
!  ---  inputs:
     &    ( ipseed,                                                     &
!  ---  outputs:
     &      stat                                                        &
     &    )

!  --- ...  sub-column set up according to overlapping assumption

      select case ( iovrsw )

        case( 0 )        ! random overlap, pick a random value at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptsw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(k,n) = rand2d(k1)
            enddo
          enddo

        case( 1 )        ! max-ran overlap

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptsw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(k,n) = rand2d(k1)
            enddo
          enddo

!  ---  first pick a random number for bottom/top layer.
!       then walk up the column: (aer's code)
!       if layer below is cloudy, use the same rand num in the layer below
!       if layer below is clear,  use a new random number

!  ---  from bottom up
          do k = 2, nlay
            k1 = k - 1
            tem1 = f_one - cldf(k1)

            do n = 1, ngptsw
              if ( cdfunc(k1,n) > tem1 ) then
                cdfunc(k,n) = cdfunc(k1,n)
              else
                cdfunc(k,n) = cdfunc(k,n) * tem1
              endif
            enddo
          enddo

!  ---  then walk down the column: (if use original author's method)
!       if layer above is cloudy, use the same rand num in the layer above
!       if layer above is clear,  use a new random number

!  ---  from top down
!         do k = nlay-1, 1, -1
!           k1 = k + 1
!           tem1 = f_one - cldf(k1)

!           do n = 1, ngptsw
!             if ( cdfunc(k1,n) > tem1 ) then
!               cdfunc(k,n) = cdfunc(k1,n)
!             else
!               cdfunc(k,n) = cdfunc(k,n) * tem1
!             endif
!           enddo
!         enddo

        case( 2 )        ! maximum overlap, pick same random numebr at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand1d, stat )

          do n = 1, ngptsw
            tem1 = rand1d(n)

            do k = 1, nlay
              cdfunc(k,n) = tem1
            enddo
          enddo

      end select

!  --- ...  generate subcolumns for homogeneous clouds

      do k = 1, nlay
        tem1 = f_one - cldf(k)

        do n = 1, ngptsw
          lcloudy(k,n) = cdfunc(k,n) >= tem1
        enddo
      enddo

      return
! ..................................
      end subroutine mcica_subcol
! ----------------------------------


! ----------------------------------
      subroutine setcoef                                                &
! ..................................
!  ---  inputs:
     &     ( pavel,tavel,h2ovmr, nlay,nlp1,                             &
!  ---  outputs:
     &       laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,                 &
     &       selffac,selffrac,indself,forfac,forfrac,indfor             &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute various coefficients needed in radiative transfer   !
!    calculations.                                                      !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       -size- !
!   pavel     - real, layer pressures (mb)                         nlay !
!   tavel     - real, layer temperatures (k)                       nlay !
!   h2ovmr    - real, layer w.v. volum mixing ratio (kg/kg)        nlay !
!   nlay/nlp1 - integer, total number of vertical layers, levels    1   !
!                                                                       !
!  outputs:                                                             !
!   laytrop   - integer, tropopause layer index (unitless)          1   !
!   jp        - real, indices of lower reference pressure          nlay !
!   jt, jt1   - real, indices of lower reference temperatures      nlay !
!                 at levels of jp and jp+1                              !
!   facij     - real, factors multiply the reference ks,           nlay !
!                 i,j=0/1 for lower/higher of the 2 appropriate         !
!                 temperatures and altitudes.                           !
!   selffac   - real, scale factor for w. v. self-continuum        nlay !
!                 equals (w. v. density)/(atmospheric density           !
!                 at 296k and 1013 mb)                                  !
!   selffrac  - real, factor for temperature interpolation of      nlay !
!                 reference w. v. self-continuum data                   !
!   indself   - integer, index of lower ref temp for selffac       nlay !
!   forfac    - real, scale factor for w. v. foreign-continuum     nlay !
!   forfrac   - real, factor for temperature interpolation of      nlay !
!                 reference w.v. foreign-continuum data                 !
!   indfor    - integer, index of lower ref temp for forfac        nlay !
!                                                                       !
!  ======================    end of definitions    ===================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(:), intent(in) :: pavel, tavel,  &
     &       h2ovmr

!  ---  outputs:
      integer, dimension(nlay), intent(out) :: indself, indfor,         &
     &       jp, jt, jt1
      integer, intent(out) :: laytrop

      real (kind=kind_phys), dimension(nlay), intent(out) :: fac00,     &
     &       fac01, fac10, fac11, selffac, selffrac, forfac, forfrac

!  ---  locals:
      real (kind=kind_phys) :: plog, fp, fp1, ft, ft1, tem1, tem2

      integer :: i, k, jp1
!
!===> ... begin here
!
      laytrop= nlay

      do k = 1, nlay

        forfac(k) = pavel(k)*stpfac / (tavel(k)*(f_one + h2ovmr(k)))

!  --- ...  find the two reference pressures on either side of the
!           layer pressure.  store them in jp and jp1.  store in fp the
!           fraction of the difference (in ln(pressure)) between these
!           two values that the layer pressure lies.

        plog  = log(pavel(k))
        jp(k) = max(1, min(58, int(36.0 - 5.0*(plog+0.04)) ))
        jp1   = jp(k) + 1
        fp    = 5.0 * (preflog(jp(k)) - plog)

!  --- ...  determine, for each reference pressure (jp and jp1), which
!          reference temperature (these are different for each reference
!          pressure) is nearest the layer temperature but does not exceed it.
!          store these indices in jt and jt1, resp. store in ft (resp. ft1)
!          the fraction of the way between jt (jt1) and the next highest
!          reference temperature that the layer temperature falls.

        tem1 = (tavel(k) - tref(jp(k))) / 15.0
        tem2 = (tavel(k) - tref(jp1  )) / 15.0
        jt (k) = max(1, min(4, int(3.0 + tem1) ))
        jt1(k) = max(1, min(4, int(3.0 + tem2) ))
        ft  = tem1 - float(jt (k) - 3)
        ft1 = tem2 - float(jt1(k) - 3)

!  --- ...  we have now isolated the layer ln pressure and temperature,
!           between two reference pressures and two reference temperatures
!           (for each reference pressure).  we multiply the pressure
!           fraction fp with the appropriate temperature fractions to get
!           the factors that will be needed for the interpolation that yields
!           the optical depths (performed in routines taugbn for band n).

        fp1 = f_one - fp
        fac10(k) = fp1 * ft
        fac00(k) = fp1 * (f_one - ft)
        fac11(k) = fp  * ft1
        fac01(k) = fp  * (f_one - ft1)

!  --- ...  if the pressure is less than ~100mb, perform a different
!           set of species interpolations.

        if ( plog > 4.56 ) then

          laytrop =  k

!  --- ...  set up factors needed to separately include the water vapor
!           foreign-continuum in the calculation of absorption coefficient.

          tem1 = (332.0 - tavel(k)) / 36.0
          indfor (k) = min(2, max(1, int(tem1)))
          forfrac(k) = tem1 - float(indfor(k))

!  --- ...  set up factors needed to separately include the water vapor
!           self-continuum in the calculation of absorption coefficient.

          tem2 = (tavel(k) - 188.0) / 7.2
          indself (k) = min(9, max(1, int(tem2)-7))
          selffrac(k) = tem2 - float(indself(k) + 7)
          selffac (k) = h2ovmr(k) * forfac(k)

        else

!  --- ...  set up factors needed to separately include the water vapor
!           foreign-continuum in the calculation of absorption coefficient.

          tem1 = (tavel(k) - 188.0) / 36.0
          indfor (k) = 3
          forfrac(k) = tem1 - f_one

          indself (k) = 0
          selffrac(k) = f_zero
          selffac (k) = f_zero

        endif

      enddo    ! end_do_k_loop

      return
! ..................................
      end subroutine setcoef
! ----------------------------------


!-----------------------------------
      subroutine spcvrtc                                                &
!...................................
!  ---  inputs:
     &     ( ssolar,cosz,sntz,albbm,albdf,sfluxzen,cldfrc,              &
     &       cf1,cf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,     &
     &       nlay, nlp1, nday,                                          &
!  ---  outputs:
     &       fxupc,fxdnc,fxup0,fxdn0,                                   &
     &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
     &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
!   purpose:  computes the shortwave radiative fluxes using two-stream  !
!             method                                                    !
!                                                                       !
!   subprograms called:  swflux                                         !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!    ssolar  - real, incoming solar flux at top                    1    !
!    cosz    - real, cosine solar zenith angle                     1    !
!    sntz    - real, secant solar zenith angle                     1    !
!    albbm   - real, surface albedo for direct beam radiation      2    !
!    albdf   - real, surface albedo for diffused radiation         2    !
!    sfluxzen- real, spectral distribution of incoming solar flux ngptsw!
!    cldfrc  - real, layer cloud fraction                         nlay  !
!    cf1     - real, >0: cloudy sky, otherwise: clear sky          1    !
!    cf0     - real, =1-cf1                                        1    !
!    taug    - real, spectral optical depth for gases        nlay*ngptsw!
!    taur    - real, optical depth for rayleigh scattering   nlay*ngptsw!
!    tauae   - real, aerosols optical depth                  nlay*nbdsw !
!    ssaae   - real, aerosols single scattering albedo       nlay*nbdsw !
!    asyae   - real, aerosols asymmetry factor               nlay*nbdsw !
!    taucw   - real, weighted cloud optical depth            nlay*nbdsw !
!    ssacw   - real, weighted cloud single scat albedo       nlay*nbdsw !
!    asycw   - real, weighted cloud asymmetry factor         nlay*nbdsw !
!    nlay,nlp1 - integer,  number of layers/levels                 1    !
!                                                                       !
!  output variables:                                                    !
!    fxupc   - real, tot sky upward flux                     nlp1*nbdsw !
!    fxdnc   - real, tot sky downward flux                   nlp1*nbdsw !
!    fxup0   - real, clr sky upward flux                     nlp1*nbdsw !
!    fxdn0   - real, clr sky downward flux                   nlp1*nbdsw !
!    ftoauc  - real, tot sky toa upwd flux                         1    !
!    ftoau0  - real, clr sky toa upwd flux                         1    !
!    ftoadc  - real, toa downward (incoming) solar flux            1    !
!    fsfcuc  - real, tot sky sfc upwd flux                         1    !
!    fsfcu0  - real, clr sky sfc upwd flux                         1    !
!    fsfcdc  - real, tot sky sfc dnwd flux                         1    !
!    fsfcd0  - real, clr sky sfc dnwd flux                         1    !
!    sfbmc   - real, tot sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdfc   - real, tot sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    sfbm0   - real, clr sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdf0   - real, clr sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    suvbfc  - real, tot sky sfc dnwd uv-b flux                    1    !
!    suvbf0  - real, clr sky sfc dnwd uv-b flux                    1    !
!                                                                       !
!  internal variables:                                                  !
!    zrefb   - real, direct beam reflectivity for clear/cloudy    nlp1  !
!    zrefd   - real, diffuse reflectivity for clear/cloudy        nlp1  !
!    ztrab   - real, direct beam transmissivity for clear/cloudy  nlp1  !
!    ztrad   - real, diffuse transmissivity for clear/cloudy      nlp1  !
!    zldbt   - real, layer beam transmittance for clear/cloudy    nlp1  !
!    ztdbt   - real, lev total beam transmittance for clr/cld     nlp1  !
!                                                                       !
!  control parameters in module "physparam"                             !
!    iswmode - control flag for 2-stream transfer schemes               !
!              = 1 delta-eddington    (joseph et al., 1976)             !
!              = 2 pifm               (zdunkowski et al., 1980)         !
!              = 3 discrete ordinates (liou, 1973)                      !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  method:                                                              !
!  -------                                                              !
!     standard delta-eddington, p.i.f.m., or d.o.m. layer calculations. !
!     kmodts  = 1 eddington (joseph et al., 1976)                       !
!             = 2 pifm (zdunkowski et al., 1980)                        !
!             = 3 discrete ordinates (liou, 1973)                       !
!                                                                       !
!  modifications:                                                       !
!  --------------                                                       !
!   original: h. barker                                                 !
!   revision: merge with rrtmg_sw: j.-j.morcrette, ecmwf, feb 2003      !
!   revision: add adjustment for earth/sun distance:mjiacono,aer,oct2003!
!   revision: bug fix for use of palbp and palbd: mjiacono, aer, nov2003!
!   revision: bug fix to apply delta scaling to clear sky: aer, dec2004 !
!   revision: code modified so that delta scaling is not done in cloudy !
!             profiles if routine cldprop is used; delta scaling can be !
!             applied by swithcing code below if cldprop is not used to !
!             get cloud properties. aer, jan 2005                       !
!   revision: uniform formatting for rrtmg: mjiacono, aer, jul 2006     !
!   revision: use exponential lookup table for transmittance: mjiacono, !
!             aer, aug 2007                                             !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: zcrit = 0.9999995 ! thresold for conservative scattering
      real (kind=kind_phys), parameter :: zsr3  = sqrt(3.0)
      real (kind=kind_phys), parameter :: od_lo = 0.06
      real (kind=kind_phys), parameter :: eps1  = 1.0e-8


!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, nday

      real (kind=kind_phys),dimension(CHK,nlay,ngptsw),intent(inout) ::      &
     &       taug, taur
      real (kind=kind_phys),dimension(CHK,nlay,nbdsw), intent(inout) ::      &
     &       taucw, ssacw, asycw, tauae, ssaae, asyae

      real (kind=kind_phys),dimension(CHK,ngptsw),intent(inout)::       &
     &       sfluxzen
      real (kind=kind_phys),dimension(CHK,nlay),  intent(inout)::cldfrc

      real (kind=kind_phys),dimension(CHK,2), intent(inout)::albbm,albdf

      real (kind=kind_phys), dimension(CHK), intent(inout) ::           &
     &       cosz, sntz, cf1, cf0, ssolar


!  ---  outputs:
      real (kind=kind_phys), dimension(CHK,nlp1,nbdsw), intent(out) ::      &
     &       fxupc, fxdnc, fxup0, fxdn0

      real (kind=kind_phys),dimension(CHK,2),intent(out):: sfbmc, sfdfc, &
     &       sfbm0, sfdf0

      real (kind=kind_phys), dimension(CHK), intent(out) ::             &
     &       suvbfc, suvbf0, ftoadc,                                    &
     &       ftoauc, ftoau0, fsfcuc, fsfcu0, fsfcdc, fsfcd0

!  ---  locals:
      real (kind=kind_phys), dimension(CHK,nlay) :: ztaus, zssas, zasys,    &
     &       zldbt0

      real (kind=kind_phys), dimension(CHK,nlp1) :: zrefb, zrefd, ztrab,    &
     &       ztrad, ztdbt, zldbt, zfu, zfd

      real (kind=kind_phys),dimension(CHK) ::                           &
     &       zsolar, ztdbt0
      real (kind=kind_phys) :: ztau1, zssa1, zasy1, ztau0, zssa0,       &
     &       zasy0, zasy3, zssaw, zasyw, zgam1, zgam2, zgam3, zgam4,    &
     &       zc0, zc1, za1, za2, zb1, zb2, zrk, zrk2, zrp, zrp1, zrm1,  &
     &       zrpp, zrkg1, zrkg3, zrkg4, zexp1, zexm1, zexp2, zexm2,     &
     &       zexp3, zexp4, zden1, ze1r45, ftind,         zrefb1,        &
     &       zrefd1, ztrab1, ztrad1,         zr1, zr2, zr3, zr4, zr5,   &
     &       zt1, zt2, zt3, zf1, zf2

      integer :: ib, ibd, jb, jg, k, kp, itind, j1, i, j
      logical, dimension(CHK) :: lmask
!
!===> ...  begin here
!
!  --- ... initialization of output fluxes

      do ib = 1, nbdsw
        do k = 1, nlp1
          fxdnc(:,k,ib) = f_zero
          fxupc(:,k,ib) = f_zero
          fxdn0(:,k,ib) = f_zero
          fxup0(:,k,ib) = f_zero
        enddo
      enddo

      ftoadc = f_zero
      ftoauc = f_zero
      ftoau0 = f_zero
      fsfcuc = f_zero
      fsfcu0 = f_zero
      fsfcdc = f_zero
      fsfcd0 = f_zero

!! --- ...  uv-b surface downward fluxes
      suvbfc  = f_zero
      suvbf0  = f_zero

!! --- ...  output surface flux components
      sfbmc = f_zero
      sfdfc = f_zero
      sfbm0 = f_zero
      sfdf0 = f_zero

!  --- ...  loop over all g-points in each band

      lab_do_jg : do jg = 1, ngptsw

        jb = NGB(jg)
        ib = jb + 1 - nblow
        ibd = idxsfc(jb)

        zsolar = ssolar * sfluxzen(:,jg)

!  --- ...  set up toa direct beam and surface values (beam and diff)

        ztdbt(:,nlp1) = f_one
        ztdbt0   = f_one

        zldbt(:,1) = f_zero
        if (ibd /= 0) then
          zrefb(:,1) = albbm(:,ibd)
          zrefd(:,1) = albdf(:,ibd)
        else
          zrefb(:,1) = 0.5 * (albbm(:,1) + albbm(:,2))
          zrefd(:,1) = 0.5 * (albdf(:,1) + albdf(:,2))
        endif
        ztrab(:,1) = f_zero
        ztrad(:,1) = f_zero

!  --- ...  compute clear-sky optical parameters, layer reflectance and transmittance

        do k = nlay, 1, -1
          kp = k + 1
         DO j1=1,NDAY
          ztau0 = max(ftiny, taur(j1,k,jg)+taug(j1,k,jg)+tauae(j1,k,ib))
          zssa0 = taur(j1,k,jg) + tauae(j1,k,ib)*ssaae(j1,k,ib)
          zasy0 = asyae(j1,k,ib)*ssaae(j1,k,ib)*tauae(j1,k,ib)
          zssaw = min( oneminus, zssa0 / ztau0 )
          zasyw = zasy0 / max( ftiny, zssa0 )

!  --- ...  saving clear-sky quantities for later total-sky usage
          ztaus(j1,k) = ztau0
          zssas(j1,k) = zssa0
          zasys(j1,k) = zasy0

!  --- ...  delta scaling for clear-sky condition
          za1 = zasyw * zasyw
          za2 = zssaw * za1

          ztau1 = (f_one - za2) * ztau0
          zssa1 = (zssaw - za2) / (f_one - za2)
!org      zasy1 = (zasyw - za1) / (f_one - za1)   ! this line is replaced by the next
          zasy1 = zasyw / (f_one + zasyw)         ! to reduce truncation error
          zasy3 = 0.75 * zasy1

!  --- ...  general two-stream expressions
          if ( iswmode == 1 ) then
            zgam1 = 1.75 - zssa1 * (f_one + zasy3)
            zgam2 =-0.25 - zssa1 * (f_one - zasy3)
            zgam3 = 0.5  - zasy3 * cosz(j1)
          elseif ( iswmode == 2 ) then               ! pifm
            zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
            zgam2 = 0.75* zssa1 * (f_one- zasy1)
            zgam3 = 0.5 - zasy3 * cosz(j1)
          elseif ( iswmode == 3 ) then               ! discrete ordinates
            zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
            zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
            zgam3 = (1.0 - zsr3 * zasy1 * cosz(j1)) * 0.5
          endif
          zgam4 = f_one - zgam3

!  --- ...  compute homogeneous reflectance and transmittance

          if ( zssaw >= zcrit ) then    ! for conservative scattering
            za1 = zgam1 * cosz(j1) - zgam3
            za2 = zgam1 * ztau1

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

            zb1 = min ( ztau1*sntz(j1) , 500.0 )
            if ( zb1 <= od_lo ) then
              zb2 = f_one - zb1 + 0.5*zb1*zb1
            else
              ftind = zb1 / (bpade + zb1)
              itind = ftind*NTBMX + 0.5
              zb2 = exp_tbl(itind)
            endif

!      ...  collimated beam
            zrefb(j1,kp) = max(f_zero, min(f_one,                          &
     &                  (za2 - za1*(f_one - zb2))/(f_one + za2) ))
            ztrab(j1,kp) = max(f_zero, min(f_one, f_one-zrefb(j1,kp) ))

!      ...  isotropic incidence
            zrefd(j1,kp) = max(f_zero, min(f_one, za2/(f_one + za2) ))
            ztrad(j1,kp) = max(f_zero, min(f_one, f_one-zrefd(j1,kp) ))

          else                          ! for non-conservative scattering
            za1 = zgam1*zgam4 + zgam2*zgam3
            za2 = zgam1*zgam3 + zgam2*zgam4
            zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
            zrk2= 2.0 * zrk

            zrp  = zrk * cosz(j1)
            zrp1 = f_one + zrp
            zrm1 = f_one - zrp
            zrpp = f_one - zrp*zrp
            if ( abs(zrpp) < eps1 )then
              zrpp = sign ( eps1 , zrpp )
            endif
            zrkg1= zrk + zgam1
            zrkg3= zrk * zgam3
            zrkg4= zrk * zgam4

            zr1  = zrm1 * (za2 + zrkg3)
            zr2  = zrp1 * (za2 - zrkg3)
            zr3  = zrk2 * (zgam3 - za2*cosz(j1))
            zr4  = zrpp * zrkg1
            zr5  = zrpp * (zrk - zgam1)

            zt1  = zrp1 * (za1 + zrkg4)
            zt2  = zrm1 * (za1 - zrkg4)
            zt3  = zrk2 * (zgam4 + za1*cosz(j1))

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

            zb1 = min ( zrk*ztau1, 500.0 )
            if ( zb1 <= od_lo ) then
              zexm1 = f_one - zb1 + 0.5*zb1*zb1
            else
              ftind = zb1 / (bpade + zb1)
              itind = ftind*NTBMX + 0.5
              zexm1 = exp_tbl(itind)
            endif
            zexp1 = f_one / zexm1

            zb2 = min ( sntz(j1)*ztau1, 500.0 )
            if ( zb2 <= od_lo ) then
              zexm2 = f_one - zb2 + 0.5*zb2*zb2
            else
              ftind = zb2 / (bpade + zb2)
              itind = ftind*NTBMX + 0.5
              zexm2 = exp_tbl(itind)
            endif
            zexp2 = f_one / zexm2
            ze1r45 = zr4*zexp1 + zr5*zexm1

!      ...  collimated beam
            if (ze1r45>=-eps1 .and. ze1r45<=eps1) then
              zrefb(j1,kp) = eps1
              ztrab(j1,kp) = zexm2
            else
              zden1 = zssa1 / ze1r45
              zrefb(j1,kp) = max(f_zero, min(f_one,                        &
     &                    (zr1*zexp1 - zr2*zexm1 - zr3*zexm2)*zden1 ))
              ztrab(j1,kp) = max(f_zero, min(f_one, zexm2*(f_one           &
     &                  - (zt1*zexp1 - zt2*zexm1 - zt3*zexp2)*zden1) ))
            endif

!      ...  diffuse beam
            zden1 = zr4 / (ze1r45 * zrkg1)
            zrefd(j1,kp) = max(f_zero, min(f_one,                          &
     &                  zgam2*(zexp1 - zexm1)*zden1 ))
            ztrad(j1,kp) = max(f_zero, min(f_one, zrk2*zden1 ))
          endif    ! end if_zssaw_block

!  --- ...  direct beam transmittance. use exponential lookup table
!           for transmittance, or expansion of exponential for low
!           optical depth

          zr1 = ztau1 * sntz(j1)
!        if(j1.eq.1)write(0,*)"zap zr1 ",kp,zr1
!        if(j1.eq.1)write(0,*)"zap ztau1 ",kp,ztau1
!        if(j1.eq.1)write(0,*)"zap sntz ",kp,sntz(j1)
          if ( zr1 <= od_lo ) then
            zexp3 = f_one - zr1 + 0.5*zr1*zr1
!        if(j1.eq.1)write(0,*)"zap 1 zexp3",kp,zexp3
          else
            ftind = zr1 / (bpade + zr1)
            itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
            zexp3 = exp_tbl(itind)
!        if(j1.eq.1)write(0,*)"zap 2 zexp3",kp,zexp3
          endif

!        if(j1.eq.1)write(0,*)"zap ztdbt",kp,ztdbt(j1,kp)
          ztdbt(j1,k)  = zexp3 * ztdbt(j1,kp)
          zldbt(j1,kp) = zexp3

!  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
!           (must use 'orig', unscaled cloud optical depth)

          zr1 = ztau0 * sntz(j1)
          if ( zr1 <= od_lo ) then
            zexp4 = f_one - zr1 + 0.5*zr1*zr1
          else
            ftind = zr1 / (bpade + zr1)
            itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
            zexp4 = exp_tbl(itind)
          endif

          zldbt0(j1,k) = zexp4
          ztdbt0(j1) = zexp4 * ztdbt0(j1)
       ENDDO  ! end of j1 loop
        enddo    ! end do_k_loop

        lmask = .true.
        call swflux                                                     &
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1, nday, lmask,                                         &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!  --- ...  compute upward and downward fluxes at levels
        do k = 1, nlp1
          fxup0(:,k,ib) = fxup0(:,k,ib) + zsolar*zfu(:,k)
          fxdn0(:,k,ib) = fxdn0(:,k,ib) + zsolar*zfd(:,k)
        enddo

!! --- ...  surface downward beam/diffused flux components
         DO j1=1,NDAY
        zb1 = zsolar(j1)*ztdbt0(j1)
        zb2 = zsolar(j1)*(zfd(j1,1) - ztdbt0(j1))

        if (ibd /= 0) then
          sfbm0(j1,ibd) = sfbm0(j1,ibd) + zb1
          sfdf0(j1,ibd) = sfdf0(j1,ibd) + zb2
        else
          zf1 = 0.5 * zb1
          zf2 = 0.5 * zb2
          sfbm0(j1,1) = sfbm0(j1,1) + zf1
          sfdf0(j1,1) = sfdf0(j1,1) + zf2
          sfbm0(j1,2) = sfbm0(j1,2) + zf1
          sfdf0(j1,2) = sfdf0(j1,2) + zf2
        endif
!       sfbm0(ibd) = sfbm0(ibd) + zsolar*ztdbt0
!       sfdf0(ibd) = sfdf0(ibd) + zsolar*(zfd(1) - ztdbt0)
         ENDDO

!  --- ...  compute total sky optical parameters, layer reflectance and transmittance

        lmask = ( cf1 > eps )
        if ( ANY(lmask) ) then

!  --- ...  set up toa direct beam and surface values (beam and diff)
!          WHERE(lmask)
          ztdbt0 = f_one
          zldbt(:,1) = f_zero
!          ENDWHERE

          do k = nlay, 1, -1
            kp = k + 1
!DEC$ VECTOR ALWAYS
          DO j1=1,NDAY
!            IF ( lmask(j1) ) THEN
            zc0 = f_one - cldfrc(j1,k)
            zc1 = cldfrc(j1,k)
            if ( zc1 > ftiny ) then          ! it is a cloudy-layer

              ztau0 = ztaus(j1,k) + taucw(j1,k,ib)
              zssa0 = zssas(j1,k) + ssacw(j1,k,ib)
              zasy0 = zasys(j1,k) + asycw(j1,k,ib)
              zssaw = min(oneminus, zssa0 / ztau0)
              zasyw = zasy0 / max(ftiny, zssa0)

!  --- ...  delta scaling for total-sky condition
              za1 = zasyw * zasyw
              za2 = zssaw * za1

              ztau1 = (f_one - za2) * ztau0
              zssa1 = (zssaw - za2) / (f_one - za2)
!org          zasy1 = (zasyw - za1) / (f_one - za1)
              zasy1 = zasyw / (f_one + zasyw)
              zasy3 = 0.75 * zasy1

!  --- ...  general two-stream expressions
              if ( iswmode == 1 ) then
                zgam1 = 1.75 - zssa1 * (f_one + zasy3)
                zgam2 =-0.25 - zssa1 * (f_one - zasy3)
                zgam3 = 0.5  - zasy3 * cosz(j1)
              elseif ( iswmode == 2 ) then               ! pifm
                zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
                zgam2 = 0.75* zssa1 * (f_one- zasy1)
                zgam3 = 0.5 - zasy3 * cosz(j1)
              elseif ( iswmode == 3 ) then               ! discrete ordinates
                zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
                zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
                zgam3 = (1.0 - zsr3 * zasy1 * cosz(j1)) * 0.5
              endif
              zgam4 = f_one - zgam3

              zrefb1 = zrefb(j1,kp)
              zrefd1 = zrefd(j1,kp)
              ztrab1 = ztrab(j1,kp)
              ztrad1 = ztrad(j1,kp)

!  --- ...  compute homogeneous reflectance and transmittance

              if ( zssaw >= zcrit ) then    ! for conservative scattering
                za1 = zgam1 * cosz(j1) - zgam3
                za2 = zgam1 * ztau1

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

                zb1 = min ( ztau1*sntz(j1) , 500.0 )
                if ( zb1 <= od_lo ) then
                  zb2 = f_one - zb1 + 0.5*zb1*zb1
                else
                  ftind = zb1 / (bpade + zb1)
                  itind = ftind*NTBMX + 0.5
                  zb2 = exp_tbl(itind)
                endif

!      ...  collimated beam
                zrefb(j1,kp)=max(f_zero, min(f_one,                      &
     &                      (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                ztrab(j1,kp)=max(f_zero, min(f_one, f_one-zrefb(j1,kp)))

!      ...  isotropic incidence
                zrefd(j1,kp)=max(f_zero,min(f_one,za2 / (f_one+za2) ))
                ztrad(j1,kp)=max(f_zero,min(f_one,f_one - zrefd(j1,kp)))

              else                          ! for non-conservative scattering
                za1 = zgam1*zgam4 + zgam2*zgam3
                za2 = zgam1*zgam3 + zgam2*zgam4
                zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
                zrk2= 2.0 * zrk

                zrp  = zrk * cosz(j1)
                zrp1 = f_one + zrp
                zrm1 = f_one - zrp
                zrpp = f_one - zrp*zrp
                if ( abs(zrpp) < eps1 )then
                  zrpp = sign ( eps1 , zrpp )
                endif
                zrkg1= zrk + zgam1
                zrkg3= zrk * zgam3
                zrkg4= zrk * zgam4

                zr1  = zrm1 * (za2 + zrkg3)
                zr2  = zrp1 * (za2 - zrkg3)
                zr3  = zrk2 * (zgam3 - za2*cosz(j1))
                zr4  = zrpp * zrkg1
                zr5  = zrpp * (zrk - zgam1)

                zt1  = zrp1 * (za1 + zrkg4)
                zt2  = zrm1 * (za1 - zrkg4)
                zt3  = zrk2 * (zgam4 + za1*cosz(j1))

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

                zb1 = min ( zrk*ztau1, 500.0 )
                if ( zb1 <= od_lo ) then
                  zexm1 = f_one - zb1 + 0.5*zb1*zb1
                else
                  ftind = zb1 / (bpade + zb1)
                  itind = ftind*NTBMX + 0.5
                  zexm1 = exp_tbl(itind)
                endif
                zexp1 = f_one / zexm1

                zb2 = min ( ztau1*sntz(j1), 500.0 )
                if ( zb2 <= od_lo ) then
                  zexm2 = f_one - zb2 + 0.5*zb2*zb2
                else
                  ftind = zb2 / (bpade + zb2)
                  itind = ftind*NTBMX + 0.5
                  zexm2 = exp_tbl(itind)
                endif
                zexp2 = f_one / zexm2
                ze1r45 = zr4*zexp1 + zr5*zexm1

!      ...  collimated beam
                if ( ze1r45>=-eps1 .and. ze1r45<=eps1 ) then
                  zrefb(j1,kp) = eps1
                  ztrab(j1,kp) = zexm2
                else
                  zden1 = zssa1 / ze1r45
                  zrefb(j1,kp) = max(f_zero, min(f_one,                    &
     &                        (zr1*zexp1-zr2*zexm1-zr3*zexm2)*zden1 ))
                  ztrab(j1,kp) = max(f_zero, min(f_one, zexm2*(f_one -     &
     &                        (zt1*zexp1-zt2*zexm1-zt3*zexp2)*zden1) ))
                endif

!      ...  diffuse beam
                zden1 = zr4 / (ze1r45 * zrkg1)
                zrefd(j1,kp) = max(f_zero, min(f_one,                      &
     &                      zgam2*(zexp1 - zexm1)*zden1 ))
                ztrad(j1,kp) = max(f_zero, min(f_one, zrk2*zden1 ))
              endif    ! end if_zssaw_block

!  --- ...  combine clear and cloudy contributions for total sky
!           and calculate direct beam transmittances

              zrefb(j1,kp) = zc0*zrefb1 + zc1*zrefb(j1,kp)
              zrefd(j1,kp) = zc0*zrefd1 + zc1*zrefd(j1,kp)
              ztrab(j1,kp) = zc0*ztrab1 + zc1*ztrab(j1,kp)
              ztrad(j1,kp) = zc0*ztrad1 + zc1*ztrad(j1,kp)

!  --- ...  direct beam transmittance. use exponential lookup table
!           for transmittance, or expansion of exponential for low
!           optical depth

              zr1 = ztau1 * sntz(j1)
              if ( zr1 <= od_lo ) then
                zexp3 = f_one - zr1 + 0.5*zr1*zr1
              else
                ftind = zr1 / (bpade + zr1)
                itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
                zexp3 = exp_tbl(itind)
              endif

              zldbt(j1,kp) = zc0*zldbt(j1,kp) + zc1*zexp3
              ztdbt(j1,k) = zldbt(j1,kp) * ztdbt(j1,kp)

!  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
!           (must use 'orig', unscaled cloud optical depth)

              zr1 = ztau0 * sntz(j1)
              if ( zr1 <= od_lo ) then
                zexp4 = f_one - zr1 + 0.5*zr1*zr1
              else
                ftind = zr1 / (bpade + zr1)
                itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
                zexp4 = exp_tbl(itind)
              endif

              ztdbt0(j1) = (zc0*zldbt0(j1,k) + zc1*zexp4) * ztdbt0(j1)

            else     ! if_zc1_block  ---  it is a clear layer

!  --- ...  direct beam transmittance
              ztdbt(j1,k) = zldbt(j1,kp) * ztdbt(j1,kp)

!  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
              ztdbt0(j1) = zldbt0(j1,k) * ztdbt0(j1)

            endif    ! end if_zc1_block
!            endif    ! mask
           enddo  ! j1=1,nday
          enddo   ! end do_k_loop

!  --- ...  perform vertical quadrature

          call swflux                                                   &
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1,nday,lmask,                                     &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!  --- ...  compute upward and downward fluxes at levels
          do k = 1, nlp1
           WHERE (lmask)
            fxupc(:,k,ib) = fxupc(:,k,ib) + zsolar*zfu(:,k)
            fxdnc(:,k,ib) = fxdnc(:,k,ib) + zsolar*zfd(:,k)
           ENDWHERE
          enddo

!! --- ...  surface downward beam/diffused flux components

           DO j1 = 1,NDAY
            if ( lmask(j1) ) then
          zb1 = zsolar(j1)*ztdbt0(j1)
          zb2 = zsolar(j1)*(zfd(j1,1) - ztdbt0(j1))

          if (ibd /= 0) then
            sfbmc(j1,ibd) = sfbmc(j1,ibd) + zb1
            sfdfc(j1,ibd) = sfdfc(j1,ibd) + zb2
          else
            zf1 = 0.5 * zb1
            zf2 = 0.5 * zb2
            sfbmc(j1,1) = sfbmc(j1,1) + zf1
            sfdfc(j1,1) = sfdfc(j1,1) + zf2
            sfbmc(j1,2) = sfbmc(j1,2) + zf1
            sfdfc(j1,2) = sfdfc(j1,2) + zf2
          endif
!         sfbmc(ibd) = sfbmc(ibd) + zsolar*ztdbt0
!         sfdfc(ibd) = sfdfc(ibd) + zsolar*(zfd(1) - ztdbt0)
            endif
           ENDDO

        endif      ! end if_cf1_block

      enddo  lab_do_jg

!  --- ...  end of g-point loop

      do ib = 1, nbdsw
!DEC$ VECTOR ALWAYS
      DO i=1,CHK
        ftoadc(i) = ftoadc(i) + fxdn0(i,nlp1,ib)
        ftoau0(i) = ftoau0(i) + fxup0(i,nlp1,ib)
        fsfcu0(i) = fsfcu0(i) + fxup0(i,1,ib)
        fsfcd0(i) = fsfcd0(i) + fxdn0(i,1,ib)
      ENDDO
      enddo

!! --- ...  uv-b surface downward flux
      ibd = nuvb - nblow + 1
      suvbf0 = fxdn0(:,1,ibd)

!      if ( cf1 <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
        do ib = 1, nbdsw
          do k = 1, nlp1
!DEC$ VECTOR ALWAYS
      DO i=1,CHK
        IF ( cf1(i) <= eps ) THEN
             fxupc(i,k,ib) = fxup0(i,k,ib)
             fxdnc(i,k,ib) = fxdn0(i,k,ib)
        ENDIF
      ENDDO
          enddo
        enddo

!DEC$ VECTOR ALWAYS
      DO i=1,CHK
        IF ( cf1(i) <= eps ) THEN
        ftoauc(i) = ftoau0(i)
        fsfcuc(i) = fsfcu0(i)
        fsfcdc(i) = fsfcd0(i)

!! --- ...  surface downward beam/diffused flux components
        sfbmc(i,1) = sfbm0(i,1)
        sfdfc(i,1) = sfdf0(i,1)
        sfbmc(i,2) = sfbm0(i,2)
        sfdfc(i,2) = sfdf0(i,2)
        ENDIF
      ENDDO

!! --- ...  uv-b surface downward flux
        suvbfc = suvbf0
!      else                        ! cloudy column, compute total-sky fluxes
        do ib = 1, nbdsw
          do k = 1, nlp1
!DEC$ VECTOR ALWAYS
      DO i=1,CHK
        IF ( cf1(i) > eps ) THEN
            fxupc(i,k,ib) = cf1(i)*fxupc(i,k,ib) + cf0(i)*fxup0(i,k,ib)
            fxdnc(i,k,ib) = cf1(i)*fxdnc(i,k,ib) + cf0(i)*fxdn0(i,k,ib)
        ENDIF
      ENDDO
          enddo
        enddo

        do ib = 1, nbdsw
!DEC$ VECTOR ALWAYS
      DO i=1,CHK
        IF ( cf1(i) > eps ) THEN
          ftoauc(i) = ftoauc(i) + fxupc(i,nlp1,ib)
          fsfcuc(i) = fsfcuc(i) + fxupc(i,1,ib)
          fsfcdc(i) = fsfcdc(i) + fxdnc(i,1,ib)
        ENDIF
      ENDDO
        enddo

!! --- ...  uv-b surface downward flux
!DEC$ VECTOR ALWAYS
      DO i=1,CHK
        IF ( cf1(i) > eps ) THEN
        suvbfc(i) = fxdnc(i,1,ibd)

!! --- ...  surface downward beam/diffused flux components
        sfbmc(i,1) = cf1(i)*sfbmc(i,1) + cf0(i)*sfbm0(i,1)
        sfbmc(i,2) = cf1(i)*sfbmc(i,2) + cf0(i)*sfbm0(i,2)
        sfdfc(i,1) = cf1(i)*sfdfc(i,1) + cf0(i)*sfdf0(i,1)
        sfdfc(i,2) = cf1(i)*sfdfc(i,2) + cf0(i)*sfdf0(i,2)
        ENDIF
      ENDDO
!      endif    ! end if_cf1_block

      return
!...................................
      end subroutine spcvrtc
!-----------------------------------


#ifdef TODO
!-----------------------------------
!rv   subroutine spcvrtm
!rv   end subroutine spcvrtm
!-----------------------------------
#endif


!-----------------------------------
      subroutine swflux                                                 &
!...................................
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       NLAY, NLP1, nday,lmask,                                    &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
!   purpose:  computes the upward and downward radiation fluxes         !
!                                                                       !
!   interface:  "swflux" is called by "spcvrc" and "spcvrm"             !
!                                                                       !
!   subroutines called : none                                           !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    zrefb(NLP1)     - layer direct beam reflectivity                   !
!    zrefd(NLP1)     - layer diffuse reflectivity                       !
!    ztrab(NLP1)     - layer direct beam transmissivity                 !
!    ztrad(NLP1)     - layer diffuse transmissivity                     !
!    zldbt(NLP1)     - layer mean beam transmittance                    !
!    ztdbt(NLP1)     - total beam transmittance at levels               !
!    NLAY, NLP1      - number of layers/levels                          !
!                                                                       !
!  output variables:                                                    !
!    zfu  (NLP1)     - upward flux at layer interface                   !
!    zfd  (NLP1)     - downward flux at layer interface                 !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1,nday
      logical, dimension(CHK), intent(in) :: lmask

      real (kind=kind_phys), dimension(CHK,nlp1), intent(in) :: zrefb,      &
     &       zrefd, ztrab, ztrad, ztdbt, zldbt

!  ---  outputs:
      real (kind=kind_phys), dimension(CHK,nlp1),intent(out) :: zfu, zfd

!  ---  locals:
      real (kind=kind_phys), dimension(CHK,nlp1)::zrupb,zrupd,zrdnd,ztdn

      real (kind=kind_phys)                 :: zden1

      integer :: k, kp, i
!
!===> ... begin here
!

!  --- ...  link lowest layer with surface

#ifdef RUNTIME_CHECKING
# define CHK nday
#endif


!DEC$ VECTOR ALWAYS
      DO i=1,CHK
        zrupb(i,1) = zrefb(i,1)        ! direct beam
        zrupd(i,1) = zrefd(i,1)        ! diffused
      ENDDO

!  --- ...  pass from bottom to top

        do k = 1, nlay
          kp = k + 1

!DEC$ VECTOR ALWAYS
      DO i=1,CHK
          zden1 = f_one / ( f_one - zrupd(i,k)*zrefd(i,kp) )
          zrupb(i,kp) = zrefb(i,kp) + ( ztrad(i,kp) *                         &
     &                ( (ztrab(i,kp) - zldbt(i,kp))*zrupd(i,k) +              &
     &                zldbt(i,kp)*zrupb(i,k)) ) * zden1
          zrupd(i,kp)=zrefd(i,kp)+ztrad(i,kp)*ztrad(i,kp)*zrupd(i,k)    &
     &                *zden1
      ENDDO
        enddo

!  --- ...  upper boundary conditions

!DEC$ VECTOR ALWAYS
      DO i=1,CHK
        ztdn (i,nlp1) = f_one
        zrdnd(i,nlp1) = f_zero
        ztdn (i,nlay) = ztrab(i,nlp1)
        zrdnd(i,nlay) = zrefd(i,nlp1)
      ENDDO

!  --- ...  pass from top to bottom

        do k = nlay, 2, -1
!DEC$ VECTOR ALWAYS
      DO i=1,CHK
          zden1 = f_one / (f_one - zrefd(i,k)*zrdnd(i,k))
          ztdn (i,k-1) = ztdbt(i,k)*ztrab(i,k) + ( ztrad(i,k) *                 &
     &                 ( (ztdn(i,k) - ztdbt(i,k)) + ztdbt(i,k) *              &
     &                 zrefb(i,k)*zrdnd(i,k) )) * zden1
          zrdnd(i,k-1)=zrefd(i,k)+ztrad(i,k)*ztrad(i,k)*zrdnd(i,k)*zden1
      ENDDO
        enddo

!  --- ...  up and down-welling fluxes at levels

        do k = 1, nlp1
!DEC$ VECTOR ALWAYS
      DO i=1,CHK
          zden1 = f_one / (f_one - zrdnd(i,k)*zrupd(i,k))
          zfu(i,k) = ( ztdbt(i,k)*zrupb(i,k) +                                &
     &             (ztdn(i,k) - ztdbt(i,k))*zrupd(i,k) ) * zden1
          zfd(i,k) = ztdbt(i,k) + ( ztdn(i,k) - ztdbt(i,k) +                    &
     &             ztdbt(i,k)*zrupb(i,k)*zrdnd(i,k) ) * zden1
      ENDDO
        enddo

#ifdef RUNTIME_CHECKING
# undef CHK
#endif
      return
!...................................
      end subroutine swflux
!-----------------------------------


!-----------------------------------
      subroutine taumol                                                 &
!...................................
!  ---  inputs:
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay,nday, &
!  ---  outputs:
     &       sfluxzen, taug, taur                                       &
     &     )

!  ==================   program usage description   ==================  !
!                                                                       !
!  description:                                                         !
!    calculate optical depths for gaseous absorption and rayleigh       !
!    scattering.                                                        !
!                                                                       !
!  subroutines called: taugb## (## = 16 - 29)                           !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                         size !
!    colamt  - real, column amounts of absorbing gases the index        !
!                    are for h2o, co2, o3, n2o, ch4, and o2,            !
!                    respectively (molecules/cm**2)          nlay*maxgas!
!    colmol  - real, total column amount (dry air+water vapor)     nlay !
!    facij   - real, for each layer, these are factors that are         !
!                    needed to compute the interpolation factors        !
!                    that multiply the appropriate reference k-         !
!                    values.  a value of 0/1 for i,j indicates          !
!                    that the corresponding factor multiplies           !
!                    reference k-value for the lower/higher of the      !
!                    two appropriate temperatures, and altitudes,       !
!                    respectively.                                 naly !
!    jp      - real, the index of the lower (in altitude) of the        !
!                    two appropriate ref pressure levels needed         !
!                    for interpolation.                            nlay !
!    jt, jt1 - integer, the indices of the lower of the two approp      !
!                    ref temperatures needed for interpolation (for     !
!                    pressure levels jp and jp+1, respectively)    nlay !
!    laytrop - integer, tropopause layer index                       1  !
!    forfac  - real, scale factor needed to foreign-continuum.     nlay !
!    forfrac - real, factor needed for temperature interpolation   nlay !
!    indfor  - integer, index of the lower of the two appropriate       !
!                    reference temperatures needed for foreign-         !
!                    continuum interpolation                       nlay !
!    selffac - real, scale factor needed to h2o self-continuum.    nlay !
!    selffrac- real, factor needed for temperature interpolation        !
!                    of reference h2o self-continuum data          nlay !
!    indself - integer, index of the lower of the two appropriate       !
!                    reference temperatures needed for the self-        !
!                    continuum interpolation                       nlay !
!    nlay    - integer, number of vertical layers                    1  !
!                                                                       !
!  output:                                                              !
!    sfluxzen- real, spectral distribution of incoming solar flux ngptsw!
!    taug    - real, spectral optical depth for gases        nlay*ngptsw!
!    taur    - real, opt depth for rayleigh scattering       nlay*ngptsw!
!                                                                       !
!  ===================================================================  !
!  ************     original subprogram description    ***************  !
!                                                                       !
!                  optical depths developed for the                     !
!                                                                       !
!                rapid radiative transfer model (rrtm)                  !
!                                                                       !
!            atmospheric and environmental research, inc.               !
!                        131 hartwell avenue                            !
!                        lexington, ma 02421                            !
!                                                                       !
!                                                                       !
!                           eli j. mlawer                               !
!                         jennifer delamere                             !
!                         steven j. taubman                             !
!                         shepard a. clough                             !
!                                                                       !
!                                                                       !
!                                                                       !
!                       email:  mlawer@aer.com                          !
!                       email:  jdelamer@aer.com                        !
!                                                                       !
!        the authors wish to acknowledge the contributions of the       !
!        following people:  patrick d. brown, michael j. iacono,        !
!        ronald e. farren, luke chen, robert bergstrom.                 !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  taumol                                                               !
!                                                                       !
!    this file contains the subroutines taugbn (where n goes from       !
!    16 to 29).  taugbn calculates the optical depths and Planck        !
!    fractions per g-value and layer for band n.                        !
!                                                                       !
!  output:  optical depths (unitless)                                   !
!           fractions needed to compute planck functions at every layer !
!           and g-value                                                 !
!                                                                       !
!  modifications:                                                       !
!                                                                       !
! revised: adapted to f90 coding, j.-j.morcrette, ecmwf, feb 2003       !
! revised: modified for g-point reduction, mjiacono, aer, dec 2003      !
! revised: reformatted for consistency with rrtmg_lw, mjiacono, aer,    !
!          jul 2006                                                     !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nday
      integer, intent(inout) :: laytrop( CHK     )

      integer, dimension(CHK    ,nlay), intent(inout):: indfor, indself,&
     &       jp, jt, jt1

      real (kind=kind_phys),dimension( CHK     ,nlay),intent(inout) ::  &
     &       colmol, fac00, fac01, fac10, fac11, forfac, forfrac,       &
     &       selffac, selffrac

      real (kind=kind_phys), dimension( CHK     ,nlay,maxgas),          &
     &      intent(inout) :: colamt

!  ---  outputs:
      real (kind=kind_phys), dimension( CHK     ,ngptsw),intent(out) :: &
     &         sfluxzen

      real (kind=kind_phys), dimension( CHK     ,nlay,ngptsw),          &
     &      intent(out) :: taug, taur

!  ---  locals:
      real (kind=kind_phys) :: fs, speccomb, specmult, colm1, colm2

      integer, dimension(CHK,nlay,nblow:nbhgh) :: id0, id1

      integer :: ibd, j, jb, js, k, klow, khgh, klim, ks, njb, ns, j1
!
!===> ... begin here

!
!  --- ...  loop over each spectral band

      do jb = nblow, nbhgh

!  --- ...  indices for layer optical depth

      do k = 1, nlay
       WHERE( k .LE. laytrop(1: nday     ) )
          id0(:,k,jb) = ((jp(:,k)-1)*5 + (jt (:,k)-1)) * nspa(jb)
          id1(:,k,jb) = ( jp(:,k)   *5 + (jt1(:,k)-1)) * nspa(jb)
       ENDWHERE

       WHERE( k .GT. laytrop(1: nday     ) )
          id0(:,k,jb) = ((jp(:,k)-13)*5 + (jt (:,k)-1)) * nspb(jb)
          id1(:,k,jb) = ((jp(:,k)-12)*5 + (jt1(:,k)-1)) * nspb(jb)
       ENDWHERE
      ENDDO

!  --- ...  calculate spectral flux at toa

        ibd = ibx(jb)
        njb = ng (jb)
        ns  = ngs(jb)

        select case (jb)

          case (16, 20, 23, 25, 26, 29)

            do j = 1, njb
              sfluxzen(:,ns+j) = sfluxref01(j,1,ibd)
            enddo

          case (27)

            do j = 1, njb
              sfluxzen(:,ns+j) = scalekur * sfluxref01(j,1,ibd)
            enddo

          case default

      DO j1 = 1,nday
            if (jb==17 .or. jb==28) then

              ks = nlay
              lab_do_k1 : do k = laytrop(j1), nlay-1
                if (jp(j1,k)<layreffr(jb)                               &
     &                          .and.jp(j1,k+1)>=layreffr(jb)) then
                  ks = k + 1
                  exit lab_do_k1
                endif
              enddo  lab_do_k1

              colm1 = colamt(j1,ks,ix1(jb))
              colm2 = colamt(j1,ks,ix2(jb))
              speccomb = colm1 + strrat(jb)*colm2
              specmult = specwt(jb) * min( oneminus, colm1/speccomb )
              js = 1 + int( specmult )
              fs = mod(specmult, f_one)

              do j = 1, njb
                sfluxzen(j1,ns+j) = sfluxref02(j,js,ibd)                   &
     &           + fs * (sfluxref02(j,js+1,ibd) - sfluxref02(j,js,ibd))
              enddo

            else

              ks = laytrop(j1)
              lab_do_k2 : do k = 1, laytrop(j1)-1
                if (jp(j1,k)<layreffr(jb)                               &
     &                          .and. jp(j1,k+1)>=layreffr(jb)) then
                  ks = k + 1
                  exit lab_do_k2
                endif
              enddo  lab_do_k2

              colm1 = colamt(j1,ks,ix1(jb))
              colm2 = colamt(j1,ks,ix2(jb))
              speccomb = colm1 + strrat(jb)*colm2
              specmult = specwt(jb) * min( oneminus, colm1/speccomb )
              js = 1 + int( specmult )
              fs = mod(specmult, f_one)

              do j = 1, njb
                sfluxzen(j1,ns+j) = sfluxref03(j,js,ibd)                   &
     &           + fs * (sfluxref03(j,js+1,ibd) - sfluxref03(j,js,ibd))
              enddo

            endif
      ENDDO

        end select

      enddo

!  --- ...  call taumol## to calculate layer optical depth

      call taumol16
      call taumol17
      call taumol18
      call taumol19
      call taumol20
      call taumol21
      call taumol22
      call taumol23
      call taumol24
      call taumol25
      call taumol26
      call taumol27
      call taumol28
      call taumol29

! =================
      contains
! =================

!-----------------------------------
      subroutine taumol16
!...................................

!  ------------------------------------------------------------------  !
!     band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb16

!  ---  locals:

      real (kind=kind_taum),dimension( CHK     ) ::                     &
     &       speccomb, specmult, tauray, fs, fs1,                       &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer,dimension( CHK     ) ::                                   &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1

!
!===> ... begin here
!

!  --- ... compute the optical depth by interpolating in ln(pressure),
!          temperature, and appropriate species.  below laytrop, the water
!          vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
!NOP$ SIMD
        tauray(1: CHK     ) = colmol(1: CHK     ,k) * rayl

        do j = 1, NG16
          taur(1: CHK     ,k,NS16+j) = tauray(1: CHK     )
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1: CHK     ))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1: CHK     ) )
        speccomb = colamt(1:CHK,k,1) + strrat(16)*colamt(1:CHK,k,5)
        specmult = 8.0 * min( oneminus, colamt(1:CHK,k,1)/speccomb )

        js = 1 + int( specmult )
        fs = mod( specmult, f_one )
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSA16,id0(1:CHK,k,16) + js))
        ind02 =       MIN(MSA16,ind01 + 1)
        ind03 =       MIN(MSA16,ind01 + 9)
        ind04 =       MIN(MSA16,ind01 + 10)
        ind11 = MAX(1,MIN(MSA16,id1(1:CHK,k,16) + js))
        ind12 =       MIN(MSA16,ind11 + 1)
        ind13 =       MIN(MSA16,ind11 + 9)
        ind14 =       MIN(MSA16,ind11 + 10)
        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG16
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK))
          taug(1:CHK,k,NS16+j) = speccomb                               &
     &        *( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)        &
     &        +  fac010 * absa(ind03,j) + fac110 * absa(ind04,j)        &
     &        +  fac001 * absa(ind11,j) + fac101 * absa(ind12,j)        &
     &        +  fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )      &
     &  + colamt(1:CHK,k,1) * (selffac(1:CHK,k) * (selfref(inds,j)      &
     &  + selffrac(1:CHK,k) * (selfref(indsp,j)-selfref(inds,j)))       &
     &  + forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)          &
     &        * (forref(indfp,j) - forref(indf,j))))
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1,nlay
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSB16,id0(1:CHK,k,16) + 1))
        ind02 =       MIN(MSB16,ind01 + 1)
        ind11 = MAX(1,MIN(MSB16,id1(1:CHK,k,16) + 1))
        ind12 =       MIN(MSB16,ind11 + 1)
       ENDWHERE

        do j = 1, NG16
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS16+j) = colamt(1:CHK,k,5)                      &
     &      *(fac00(1:CHK,k)*absb(ind01,j)+fac10(1:CHK,k)*absb(ind02,j) &
     &      + fac01(1:CHK,k)*absb(ind11,j)+fac11(1:CHK,k)*absb(ind12,j))
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol16
!-----------------------------------


!-----------------------------------
      subroutine taumol17
!...................................

!  ------------------------------------------------------------------  !
!     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)         !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb17

!  ---  locals:

      real (kind=kind_taum),dimension(CHK) ::                           &
     &       speccomb, specmult, tauray, fs, fs1,                       &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer,dimension(CHK) ::                                         &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1

!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
!NOP$ SIMD
        tauray(1:CHK) = colmol(1:CHK,k) * rayl

        do j = 1, NG17
          taur(1:CHK,k,NS17+j) = tauray(1:CHK)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        speccomb = colamt(1:CHK,k,1) + strrat(17)*colamt(1:CHK,k,2)
        specmult = 8.0 * min(oneminus, colamt(1:CHK,k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSA17,id0(1:CHK,k,17) + js))
        ind02 =       MIN(MSA17,ind01 + 1)
        ind03 =       MIN(MSA17,ind01 + 9)
        ind04 =       MIN(MSA17,ind01 + 10)
        ind11 = MAX(1,MIN(MSA17,id1(1:CHK,k,17) + js))
        ind12 =       MIN(MSA17,ind11 + 1)
        ind13 =       MIN(MSA17,ind11 + 9)
        ind14 =       MIN(MSA17,ind11 + 10)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG17
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS17+j) = speccomb                                     &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &  + colamt(1:CHK,k,1) * (selffac(1:CHK,k) * (selfref(inds,j)      &
     &  + selffrac(1:CHK,k) * (selfref(indsp,j)-selfref(inds,j)))       &
     &  + forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)          &
     &        * (forref(indfp,j) - forref(indf,j))))
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        speccomb = colamt(1:CHK,k,1) + strrat(17)*colamt(1:CHK,k,2)
        specmult = 4.0 * min(oneminus, colamt(1:CHK,k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSB17,id0(1:CHK,k,17) + js))
        ind02 =       MIN(MSB17,ind01 + 1)
        ind03 =       MIN(MSB17,ind01 + 5)
        ind04 =       MIN(MSB17,ind01 + 6)
        ind11 = MAX(1,MIN(MSB17,id1(1:CHK,k,17) + js))
        ind12 =       MIN(MSB17,ind11 + 1)
        ind13 =       MIN(MSB17,ind11 + 5)
        ind14 =       MIN(MSB17,ind11 + 6)

        indf = indfor(1:CHK,k)
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG17
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS17+j) = speccomb                                     &
     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
     &  + colamt(1:CHK,k,1) * forfac(1:CHK,k) * (forref(indf,j)         &
     &  + forfrac(1:CHK,k) * (forref(indfp,j) - forref(indf,j)))
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol17
!-----------------------------------


!-----------------------------------
      subroutine taumol18
!...................................

!  ------------------------------------------------------------------  !
!     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb18

!  ---  locals:
      real (kind=kind_taum),dimension(CHK) ::                           &
     &       speccomb, specmult, tauray, fs, fs1,                       &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer,dimension(CHK) ::                                         &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1
!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
!NOP$ SIMD
        tauray(1:CHK) = colmol(1:CHK,k) * rayl

        do j = 1, NG18
          taur(1:CHK,k,NS18+j) = tauray(1:CHK)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        speccomb = colamt(1:CHK,k,1) + strrat(18)*colamt(1:CHK,k,5)
        specmult = 8.0 * min(oneminus, colamt(1:CHK,k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSA18,id0(1:CHK,k,18) + js))
        ind02 =       MIN(MSA18,ind01 + 1)
        ind03 =       MIN(MSA18,ind01 + 9)
        ind04 =       MIN(MSA18,ind01 + 10)
        ind11 = MAX(1,MIN(MSA18,id1(1:CHK,k,18) + js))
        ind12 =       MIN(MSA18,ind11 + 1)
        ind13 =       MIN(MSA18,ind11 + 9)
        ind14 =       MIN(MSA18,ind11 + 10)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG18
!!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
!!NOP$ IVDEP
          taug(1:CHK,k,NS18+j) = speccomb                               &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &  + colamt(1:CHK,k,1) * (selffac(1:CHK,k) * (selfref(inds,j)      &
     &  + selffrac(1:CHK,k) * (selfref(indsp,j)-selfref(inds,j)))       &
     &  + forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)          &
     &        * (forref(indfp,j) - forref(indf,j))))
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSB18,id0(:,k,18) + 1))
        ind02 =       MIN(MSB18,ind01 + 1)
        ind11 = MAX(1,MIN(MSB18,id1(:,k,18) + 1))
        ind12 =       MIN(MSB18,ind11 + 1)
       ENDWHERE

        do j = 1, NG18
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS18+j) = colamt(1:CHK,k,5)                      &
     &   * (fac00(1:CHK,k)*absb(ind01,j)+fac10(1:CHK,k)*absb(ind02,j)   &
     &   +  fac01(1:CHK,k)*absb(ind11,j)+fac11(1:CHK,k)*absb(ind12,j) )
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol18
!-----------------------------------


!-----------------------------------
      subroutine taumol19
!...................................

!  ------------------------------------------------------------------  !
!     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb19

!  ---  locals:
      real (kind=kind_taum),dimension(CHK) ::                           &
     &       speccomb, specmult, tauray, fs, fs1,                       &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer,dimension(CHK) ::                                         &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      do k = 1, nlay
!NOP$ SIMD
        tauray(1:CHK) = colmol(1:CHK,k) * rayl

        do j = 1, NG19
          taur(1:CHK,k,NS19+j) = tauray(1:CHK)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        speccomb = colamt(1:CHK,k,1) + strrat(19)*colamt(1:CHK,k,2)
        specmult = 8.0 * min(oneminus, colamt(1:CHK,k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSA19,id0(1:CHK,k,19) + js))
        ind02 =       MIN(MSA19,ind01 + 1)
        ind03 =       MIN(MSA19,ind01 + 9)
        ind04 =       MIN(MSA19,ind01 + 10)
        ind11 = MAX(1,MIN(MSA19,id1(1:CHK,k,19) + js))
        ind12 =       MIN(MSA19,ind11 + 1)
        ind13 =       MIN(MSA19,ind11 + 9)
        ind14 =       MIN(MSA19,ind11 + 10)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG19
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS19+j) = speccomb                               &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     & + colamt(1:CHK,k,1) * (selffac(1:CHK,k) * (selfref(inds,j)       &
     & + selffrac(1:CHK,k) * (selfref(indsp,j)-selfref(inds,j)))        &
     & + forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)           &
     &        * (forref(indfp,j) - forref(indf,j))))
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSB19,id0(1:CHK,k,19) + 1))
        ind02 =       MIN(MSB19,ind01 + 1)
        ind11 = MAX(1,MIN(MSB19,id1(1:CHK,k,19) + 1))
        ind12 =       MIN(MSB19,ind11 + 1)
       ENDWHERE

        do j = 1, NG19
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS19+j) = colamt(1:CHK,k,2)                      &
     & * ( fac00(1:CHK,k)*absb(ind01,j) + fac10(1:CHK,k)*absb(ind02,j)  &
     & +   fac01(1:CHK,k)*absb(ind11,j) + fac11(1:CHK,k)*absb(ind12,j) ) 
       ENDWHERE
        enddo
      enddo

!...................................
      end subroutine taumol19
!-----------------------------------


!-----------------------------------
      subroutine taumol20
!...................................

!  ------------------------------------------------------------------  !
!     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)                 !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb20

!  ---  locals:

      real (kind=kind_taum),dimension(CHK) :: tauray

      integer,dimension(CHK) ::                                         &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
!NOP$ SIMD
        tauray(1:CHK) = colmol(1:CHK,k) * rayl

        do j = 1, NG20
          taur(1:CHK,k,NS20+j) = tauray(1:CHK)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSA20,id0(1:CHK,k,20) + 1))
        ind02 =       MIN(MSA20,ind01 + 1)
        ind11 = MAX(1,MIN(MSA20,id1(1:CHK,k,20) + 1))
        ind12 =       MIN(MSA20,ind11 + 1)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG20
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS20+j) = colamt(1:CHK,k,1)                      &
     & * ( (fac00(1:CHK,k)*absa(ind01,j) + fac10(1:CHK,k)*absa(ind02,j) &
     & +    fac01(1:CHK,k)*absa(ind11,j) + fac11(1:CHK,k)*absa(ind12,j))&
     & +   selffac(1:CHK,k) * (selfref(inds,j) + selffrac(1:CHK,k)      &
     & *   (selfref(indsp,j) - selfref(inds,j)))                        &
     & +   forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)         &
     & *   (forref(indfp,j) - forref(indf,j))) )                        &
     & + colamt(1:CHK,k,5) * absch4(j)
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSB20,id0(1:CHK,k,20) + 1))
        ind02 =       MIN(MSB20,ind01 + 1)
        ind11 = MAX(1,MIN(MSB20,id1(1:CHK,k,20) + 1))
        ind12 =       MIN(MSB20,ind11 + 1)

        indf = indfor(1:CHK,k)
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG20
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS20+j) = colamt(1:CHK,k,1)                      &
     & * ( fac00(1:CHK,k)*absb(ind01,j) + fac10(1:CHK,k)*absb(ind02,j)  &
     & +   fac01(1:CHK,k)*absb(ind11,j) + fac11(1:CHK,k)*absb(ind12,j)  &
     & +   forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)         &
     & *   (forref(indfp,j) - forref(indf,j))) )                        &
     & + colamt(1:CHK,k,5) * absch4(j)
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol20
!-----------------------------------


!-----------------------------------
      subroutine taumol21
!...................................

!  ------------------------------------------------------------------  !
!     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)         !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb21

!  ---  locals:

      real (kind=kind_taum),dimension(CHK) ::                           &
     &       speccomb, specmult, tauray, fs, fs1,                       &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer,dimension(CHK) ::                                         &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
!NOP$ SIMD
        tauray(1:CHK) = colmol(1:CHK,k) * rayl

        do j = 1, NG21
          taur(1:CHK,k,NS21+j) = tauray(1:CHK)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        speccomb = colamt(1:CHK,k,1) + strrat(21)*colamt(1:CHK,k,2)
        specmult = 8.0 * min(oneminus, colamt(1:CHK,k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSA21,id0(1:CHK,k,21) + js))
        ind02 =       MIN(MSA21,ind01 + 1)
        ind03 =       MIN(MSA21,ind01 + 9)
        ind04 =       MIN(MSA21,ind01 + 10)
        ind11 = MAX(1,MIN(MSA21,id1(1:CHK,k,21) + js))
        ind12 =       MIN(MSA21,ind11 + 1)
        ind13 =       MIN(MSA21,ind11 + 9)
        ind14 =       MIN(MSA21,ind11 + 10)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG21
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS21+j) = speccomb                               &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     & + colamt(1:CHK,k,1) * (selffac(1:CHK,k) * (selfref(inds,j)       &
     & + selffrac(1:CHK,k) * (selfref(indsp,j) - selfref(inds,j)))      &
     & + forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)           &
     & * (forref(indfp,j) - forref(indf,j))))
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
	speccomb = colamt(1:CHK,k,1) + strrat(21)*colamt(1:CHK,k,2)
        specmult = 4.0 * min(oneminus, colamt(1:CHK,k,1) / speccomb)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSB21,id0(:,k,21) + js))
        ind02 =       MIN(MSB21,ind01 + 1)
        ind03 =       MIN(MSB21,ind01 + 5)
        ind04 =       MIN(MSB21,ind01 + 6)
        ind11 = MAX(1,MIN(MSB21,id1(:,k,21) + js))
        ind12 =       MIN(MSB21,ind11 + 1)
        ind13 =       MIN(MSB21,ind11 + 5)
        ind14 =       MIN(MSB21,ind11 + 6)

        indf = indfor(1:CHK,k)
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG21
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS21+j) = speccomb                                     &
     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
     &        + colamt(1:CHK,k,1) * forfac(1:CHK,k) * (forref(indf,j)               &
     &        + forfrac(1:CHK,k) * (forref(indfp,j) - forref(indf,j)))
       ENDWHERE
        enddo
      enddo


!...................................
      end subroutine taumol21
!-----------------------------------


!-----------------------------------
      subroutine taumol22
!...................................

!  ------------------------------------------------------------------  !
!     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)               !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb22

!  ---  locals:
      real (kind=kind_taum),dimension(CHK) ::                           &
     &       speccomb, specmult, tauray, fs, fs1,                       &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer,dimension(CHK) ::                                         &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1
      real (kind=kind_taum) :: o2adj, o2tem,o2cont(CHK)

!
!===> ... begin here
!
!  --- ...  the following factor is the ratio of total o2 band intensity (lines
!           and mate continuum) to o2 band intensity (line only). it is needed
!           to adjust the optical depths since the k's include only lines.

      o2adj = 1.6
      o2tem = 4.35e-4 / (350.0*2.0)
      
!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
!NOP$ SIMD
        tauray(1:CHK) = colmol(1:CHK,k) * rayl

        do j = 1, NG22
          taur(1:CHK,k,NS22+j) = tauray(1:CHK)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        o2cont   = o2tem * colamt(:,k,6)
        speccomb = colamt(1:CHK,k,1) + strrat(22)*colamt(1:CHK,k,6)
        specmult = 8.0 * min(oneminus, colamt(1:CHK,k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSA22,id0(1:CHK,k,22) + js))
        ind02 =       MIN(MSA22,ind01 + 1)
        ind03 =       MIN(MSA22,ind01 + 9)
        ind04 =       MIN(MSA22,ind01 + 10)
        ind11 = MAX(1,MIN(MSA22,id1(1:CHK,k,22) + js))
        ind12 =       MIN(MSA22,ind11 + 1)
        ind13 =       MIN(MSA22,ind11 + 9)
        ind14 =       MIN(MSA22,ind11 + 10)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG22
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS22+j) = speccomb                               &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &        + colamt(1:CHK,k,1) * (selffac(1:CHK,k) * (selfref(inds,j)&
     &        + selffrac(1:CHK,k) * (selfref(indsp,j)-selfref(inds,j))) &
     &        + forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)    &
     &        * (forref(indfp,j) - forref(indf,j)))) + o2cont
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSB22,id0(1:CHK,k,22) + 1))
        ind02 =       MIN(MSB22,ind01 + 1)
        ind11 = MAX(1,MIN(MSB22,id1(1:CHK,k,22) + 1))
        ind12 =       MIN(MSB22,ind11 + 1)
       ENDWHERE

        do j = 1, NG22
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          o2cont = o2tem * colamt(:,k,6)
          taug(:,k,NS22+j) = colamt(:,k,6) * o2adj                          &
     &        * ( fac00(:,k)*absb(ind01,j) + fac10(:,k)*absb(ind02,j)       &
     &        +   fac01(:,k)*absb(ind11,j) + fac11(:,k)*absb(ind12,j) )     &
     &        + o2cont
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol22
!-----------------------------------


!-----------------------------------
      subroutine taumol23
!...................................

!  ------------------------------------------------------------------  !
!     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)            !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb23

!  ---  locals:
      integer, dimension(1:CHK) :: ind01, ind02, ind11, ind12
      integer, dimension(1:CHK) :: inds, indf, indsp, indfp
      integer ::  j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        do j = 1, NG23
          taur(1:CHK,k,NS23+j) = colmol(1:CHK,k) * rayl(j)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSA23,id0(1:CHK,k,23) + 1))
        ind02 =       MIN(MSA23,ind01 + 1)
        ind11 = MAX(1,MIN(MSA23,id1(1:CHK,k,23) + 1))
        ind12 =       MIN(MSA23,ind11 + 1)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG23
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS23+j) = colamt(1:CHK,k,1) * (givfac            &
     & * ( fac00(1:CHK,k)*absa(ind01,j) + fac10(1:CHK,k)*absa(ind02,j)  &
     & +   fac01(1:CHK,k)*absa(ind11,j) + fac11(1:CHK,k)*absa(ind12,j) )&
     & + selffac(1:CHK,k) * (selfref(inds,j) + selffrac(1:CHK,k)        &
     & * (selfref(indsp,j) - selfref(inds,j)))                          &
     & + forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)           &
     & * (forref(indfp,j) - forref(indf,j))))
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
        do j = 1, NG23
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(:,k,NS23+j) = f_zero
       ENDWHERE
        enddo
      enddo

!...................................
      end subroutine taumol23
!-----------------------------------


!-----------------------------------
      subroutine taumol24
!...................................

!  ------------------------------------------------------------------  !
!     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb24

!  ---  locals:
      real (kind=kind_taum),dimension(CHK) ::                           &
     &       speccomb, specmult, tauray, fs, fs1,                       &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer,dimension(CHK) ::                                         &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        speccomb = colamt(1:CHK,k,1) + strrat(24)*colamt(1:CHK,k,6)
        specmult = 8.0 * min(oneminus, colamt(1:CHK,k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSA24,id0(1:CHK,k,24) + js))
        ind02 =       MIN(MSA24,ind01 + 1)
        ind03 =       MIN(MSA24,ind01 + 9)
        ind04 =       MIN(MSA24,ind01 + 10)
        ind11 = MAX(1,MIN(MSA24,id1(1:CHK,k,24) + js))
        ind12 =       MIN(MSA24,ind11 + 1)
        ind13 =       MIN(MSA24,ind11 + 9)
        ind14 =       MIN(MSA24,ind11 + 10)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG24
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS24+j) = speccomb                               &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &        + colamt(1:CHK,k,3) * abso3a(j) + colamt(1:CHK,k,1)       &
     &        * (selffac(1:CHK,k) * (selfref(inds,j) + selffrac(1:CHK,k)&
     &        * (selfref(indsp,j) - selfref(inds,j)))                   &
     &        + forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)    &
     &        * (forref(indfp,j) - forref(indf,j))))

          taur(1:CHK,k,NS24+j) = colmol(1:CHK,k)                        &
     &           * (rayla(j,js) + fs*(rayla(j,js+1) - rayla(j,js)))
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSB24,id0(1:CHK,k,24) + 1))
        ind02 =       MIN(MSB24,ind01 + 1)
        ind11 = MAX(1,MIN(MSB24,id1(1:CHK,k,24) + 1))
        ind12 =       MIN(MSB24,ind11 + 1)
       ENDWHERE

        do j = 1, NG24
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS24+j) = colamt(1:CHK,k,6)                      &
     & * ( fac00(1:CHK,k)*absb(ind01,j) + fac10(1:CHK,k)*absb(ind02,j)  &
     & +   fac01(1:CHK,k)*absb(ind11,j) + fac11(1:CHK,k)*absb(ind12,j) )&
     &        + colamt(1:CHK,k,3) * abso3b(j)

          taur(1:CHK,k,NS24+j) = colmol(1:CHK,k) * raylb(j)
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol24
!-----------------------------------


!-----------------------------------
      subroutine taumol25
!...................................

!  ------------------------------------------------------------------  !
!     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)           !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb25

!  ---  locals:
      integer, dimension(1:CHK) :: ind01, ind02, ind11, ind12
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        do j = 1, NG25
!NOP$ SIMD
          taur(1:CHK,k,NS25+j) = colmol(1:CHK,k) * rayl(j)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSA25,id0(1:CHK,k,25) + 1))
        ind02 =       MIN(MSA25,ind01 + 1)
        ind11 = MAX(1,MIN(MSA25,id1(1:CHK,k,25) + 1))
        ind12 =       MIN(MSA25,ind11 + 1)
       ENDWHERE

        do j = 1, NG25
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS25+j) = colamt(1:CHK,k,1)                      &
     & * ( fac00(1:CHK,k)*absa(ind01,j) + fac10(1:CHK,k)*absa(ind02,j)  &
     & +   fac01(1:CHK,k)*absa(ind11,j) + fac11(1:CHK,k)*absa(ind12,j) )&
     &        + colamt(1:CHK,k,3) * abso3a(j) 
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
        do j = 1, NG25
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS25+j) = colamt(1:CHK,k,3) * abso3b(j) 
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol25
!-----------------------------------


!-----------------------------------
      subroutine taumol26
!...................................

!  ------------------------------------------------------------------  !
!     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)       !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb26

!  ---  locals:
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        do j = 1, NG26
!NOP$ SIMD
          taug(1:CHK,k,NS26+j) = f_zero
!NOP$ SIMD
          taur(1:CHK,k,NS26+j) = colmol(1:CHK,k) * rayl(j) 
        enddo
      enddo

      return
!...................................
      end subroutine taumol26
!-----------------------------------


!-----------------------------------
      subroutine taumol27
!...................................

!  ------------------------------------------------------------------  !
!     band 27:  29000-38000 cm-1 (low - o3; high - o3)                 !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb27
!
!  ---  locals:
      integer, dimension(1:CHK) :: ind01, ind02, ind11, ind12
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        do j = 1, NG27
!NOP$ SIMD
          taur(1:CHK,k,NS27+j) = colmol(1:CHK,k) * rayl(j)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSA27,id0(1:CHK,k,27) + 1))
        ind02 =       MIN(MSA27,ind01 + 1)
        ind11 = MAX(1,MIN(MSA27,id1(1:CHK,k,27) + 1))
        ind12 =       MIN(MSA27,ind11 + 1)
       ENDWHERE

        do j = 1, NG27
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS27+j) = colamt(1:CHK,k,3)                                  &
     & * ( fac00(1:CHK,k)*absa(ind01,j) + fac10(1:CHK,k)*absa(ind02,j)       &
     & +   fac01(1:CHK,k)*absa(ind11,j) + fac11(1:CHK,k)*absa(ind12,j) )
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSB27,id0(1:CHK,k,27) + 1))
        ind02 =       MIN(MSB27,ind01 + 1)
        ind11 = MAX(1,MIN(MSB27,id1(1:CHK,k,27) + 1))
        ind12 =       MIN(MSB27,ind11 + 1)
       ENDWHERE

        do j = 1, NG27
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS27+j) = colamt(1:CHK,k,3)                                  &
     & * ( fac00(1:CHK,k)*absb(ind01,j) + fac10(1:CHK,k)*absb(ind02,j)       &
     & +   fac01(1:CHK,k)*absb(ind11,j) + fac11(1:CHK,k)*absb(ind12,j) )
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol27
!-----------------------------------


!-----------------------------------
      subroutine taumol28
!...................................

!  ------------------------------------------------------------------  !
!     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)           !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb28

!  ---  locals:

      real (kind=kind_taum),dimension(CHK) ::                           &
     &       speccomb, specmult, tauray, fs, fs1,                       &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer,dimension(CHK) ::                                         &
     & ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14,          &
     & inds, indf, indsp, indfp, js
      integer :: j, k, j1

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
!NOP$ SIMD
        tauray(1:CHK) = colmol(1:CHK,k) * rayl

        do j = 1, NG28
          taur(1:CHK,k,NS28+j) = tauray(1:CHK)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        speccomb = colamt(1:CHK,k,3) + strrat(28)*colamt(1:CHK,k,6)
        specmult = 8.0 * min(oneminus, colamt(1:CHK,k,3) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSA28,id0(1:CHK,k,28) + js))
        ind02 =       MIN(MSA28,ind01 + 1)
        ind03 =       MIN(MSA28,ind01 + 9)
        ind04 =       MIN(MSA28,ind01 + 10)
        ind11 = MAX(1,MIN(MSA28,id1(1:CHK,k,28) + js))
        ind12 =       MIN(MSA28,ind11 + 1)
        ind13 =       MIN(MSA28,ind11 + 9)
        ind14 =       MIN(MSA28,ind11 + 10)
       ENDWHERE

        do j = 1, NG28
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS28+j) = speccomb                               &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        speccomb = colamt(1:CHK,k,3) + strrat(28)*colamt(1:CHK,k,6)
        specmult = 4.0 * min(oneminus, colamt(1:CHK,k,3) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(1:CHK,k)
        fac010 = fs1 * fac10(1:CHK,k)
        fac100 = fs  * fac00(1:CHK,k)
        fac110 = fs  * fac10(1:CHK,k)
        fac001 = fs1 * fac01(1:CHK,k)
        fac011 = fs1 * fac11(1:CHK,k)
        fac101 = fs  * fac01(1:CHK,k)
        fac111 = fs  * fac11(1:CHK,k)

        ind01 = MAX(1,MIN(MSB28,id0(1:CHK,k,28) + js))
        ind02 =       MIN(MSB28,ind01 + 1)
        ind03 =       MIN(MSB28,ind01 + 5)
        ind04 =       MIN(MSB28,ind01 + 6)
        ind11 = MAX(1,MIN(MSB28,id1(1:CHK,k,28) + js))
        ind12 =       MIN(MSB28,ind11 + 1)
        ind13 =       MIN(MSB28,ind11 + 5)
        ind14 =       MIN(MSB28,ind11 + 6)
       ENDWHERE

        do j = 1, NG28
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS28+j) = speccomb                                     &
     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol28
!-----------------------------------


!-----------------------------------
      subroutine taumol29
!...................................

!  ------------------------------------------------------------------  !
!     band 29:  820-2600 cm-1 (low - h2o; high - co2)                  !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb29

!  ---  locals:
      real (kind=kind_taum), dimension(1:CHK) :: tauray

      integer, dimension(1:CHK) :: ind01, ind02, ind11, ind12
      integer, dimension(1:CHK) :: inds, indf, indsp, indfp
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
!NOP$ SIMD
        tauray(1:CHK) = colmol(1:CHK,k) * rayl

        do j = 1, NG29
          taur(1:CHK,k,NS29+j) = tauray(1:CHK)
        enddo
      enddo

      do k = 1, MAXVAL(laytrop(1:CHK))
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSA29,id0(1:CHK,k,29) + 1))
        ind02 =       MIN(MSA29,ind01 + 1)
        ind11 = MAX(1,MIN(MSA29,id1(1:CHK,k,29) + 1))
        ind12 =       MIN(MSA29,ind11 + 1)

        inds = indself(1:CHK,k)
        indf = indfor (1:CHK,k)
        indsp= inds + 1
        indfp= indf + 1
       ENDWHERE

        do j = 1, NG29
!NOP$ SIMD
       WHERE( k .LE. laytrop(1:CHK) )
          taug(1:CHK,k,NS29+j) = colamt(1:CHK,k,1)                      &
     & * ( (fac00(1:CHK,k)*absa(ind01,j) + fac10(1:CHK,k)*absa(ind02,j) &
     & +    fac01(1:CHK,k)*absa(ind11,j) + fac11(1:CHK,k)*absa(ind12,j))&
     & +  selffac(1:CHK,k) * (selfref(inds,j) + selffrac(1:CHK,k)       &
     & *  (selfref(indsp,j) - selfref(inds,j)))                         &
     & +  forfac(1:CHK,k) * (forref(indf,j) + forfrac(1:CHK,k)          &
     & *  (forref(indfp,j) - forref(indf,j))))                          &
     & +  colamt(1:CHK,k,2) * absco2(j)
       ENDWHERE
        enddo
      enddo

      do k = MINVAL(laytrop(1:CHK))+1, NLAY
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
        ind01 = MAX(1,MIN(MSB29,id0(1:CHK,k,29) + 1))
        ind02 =       MIN(MSB29,ind01 + 1)
        ind11 = MAX(1,MIN(MSB29,id1(1:CHK,k,29) + 1))
        ind12 =       MIN(MSB29,ind11 + 1)
       ENDWHERE

        do j = 1, NG29
!NOP$ SIMD
       WHERE( k .GE. laytrop(1:CHK) )
          taug(1:CHK,k,NS29+j) = colamt(1:CHK,k,2)                      &
     & * ( fac00(1:CHK,k)*absb(ind01,j) + fac10(1:CHK,k)*absb(ind02,j)  &
     & +   fac01(1:CHK,k)*absb(ind11,j) + fac11(1:CHK,k)*absb(ind12,j) )&
     & + colamt(1:CHK,k,1) * absh2o(j) 
       ENDWHERE
        enddo
      enddo

      return
!...................................
      end subroutine taumol29
!-----------------------------------

!...................................
      end subroutine taumol
!-----------------------------------

!
!........................................!
      end module module_radsw_main_nmmb  !
!========================================!

