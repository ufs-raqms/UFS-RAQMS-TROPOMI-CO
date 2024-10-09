module gocart_mod

  use chem_rc_mod
  use chem_types_mod
  use chem_const_mod,  only : cp, grvity,rv,xlv, mwdry, mw_so2_aer, mw_so4_aer, &
                              p1000, rd, epsilc
  use chem_tracers_mod
  use chem_config_mod, only : BURN_OPT_ENABLE,       &
                              CHEM_OPT_NONE,         &
                              CHEM_OPT_GOCART,       &
                              CHEM_OPT_GOCART_RACM,  &
                              CHEM_OPT_RACM_SOA_VBS, &
                              CHEM_OPT_MAX,          &
                              CTRA_OPT_GRELL,        &
                              DMSE_OPT_ENABLE,       &
                              DUST_OPT_AFWA,         &
                              DUST_OPT_FENGSHA,      &
                              DUST_OPT_GOCART,       &
                              DUST_OPT_NONE,         &
                              FIRE_OPT_GBBEPx,       &
                              FIRE_OPT_MODIS,        &
                              SEAS_OPT_NONE,         &
                              WDLS_OPT_GSD,          &
                              WDLS_OPT_NGAC
  use chem_iodata_mod,only : dofirehour
  use gocart_prep_mod
  use gocart_settling_mod
  use gocart_aerosols_mod
  use gocart_dmsemis_mod
  use gocart_chem_mod
  use gocart_diag_mod
  use plume_rise_mod
  use vash_settling_mod
  use dep_mod
  use dust_mod
  use seas_mod
  use opt_mod
  use chem_comm_mod

  implicit none

  public :: gocart_advance

contains

  subroutine gocart_init(config, rc)

   type(chem_config_type), intent(in) :: config
   integer,     optional, intent(out) :: rc

   ! -- local variables
   integer :: localrc, n

   ! -- begin
   if (present(rc)) rc = CHEM_RC_SUCCESS

   ! -- initialize dust module
   ! -- set default values
   select case (config % dust_opt)
     case (DUST_OPT_AFWA   )
       dust_alpha = afwa_alpha
       dust_gamma = afwa_gamma
     case (DUST_OPT_FENGSHA)
       dust_alpha    = fengsha_alpha
       dust_gamma    = fengsha_gamma
       dust_calcdrag = config % dust_calcdrag
     case (DUST_OPT_GOCART )
       dust_alpha = gocart_alpha
       dust_gamma = gocart_gamma
     case default
       call chem_rc_set(CHEM_RC_FAILURE, msg="Dust option not implemented", &
         file=__FILE__, line=__LINE__, rc=rc)
       return
   end select
   ! -- replace with input values if available
   if (config % dust_alpha > 0._CHEM_KIND_R4) dust_alpha = config % dust_alpha
   if (config % dust_gamma > 0._CHEM_KIND_R4) dust_gamma = config % dust_gamma

   ! -- initialize sea salt module
   ! -- replace first default nbins parameters with input values if available
   if (any(config % seas_emis_scale > 0._CHEM_KIND_R4)) then
     n = min(number_ss_bins, size(config % seas_emis_scale))
     emission_scale(1:n) = config % seas_emis_scale(1:n)
   end if
   if (config % seas_emis_scheme > 0) emission_scheme = config % seas_emis_scheme

   ! -- initialize wet deposition module
   call dep_wet_ls_init(config, localrc)
   if (chem_rc_check(localrc, msg="Failure initializing wet deposition", &
     file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine gocart_init


  subroutine gocart_advance(readrestart, chem_opt, chem_in_opt, chem_conv_tr, &
    biomass_burn_opt, seas_opt, dust_opt, dmsemis_opt, wetdep_ls_opt, call_chemistry, &
    aer_ra_feedback, aer_ra_frq, plumerise_flag, plumerisefire_frq, &
    kemit, ktau, dts, current_month, tz, julday,      &
    p_gocart, clayfrac, dm0, emiss_ab, emiss_abu,                         &
    emiss_ash_dt, emiss_ash_height, emiss_ash_mass, &
    emiss_tr_dt, emiss_tr_height, emiss_tr_mass, ero1, ero2, ero3, rdrag, uthr, ssm, &
    h2o2_backgd, no3_backgd, oh_backgd, plume, sandfrac, th_pvsrf,  &
    area, hf2d, pb2d, rc2d, rn2d, rsds, slmsk2d, snwdph2d, stype2d,       &
    ts2d, us2d, vtype2d, vfrac2d, zorl2d, dqdt, exch, ph3d, phl3d, pr3d, prl3d, &
    sm3d, tk3d, us3d, vs3d, ws3d, tr3d_in, tr3d_out, trcm, trab, truf, trdf, trdp, &
    ext_cof, sscal, asymp, &
    p10, pm25, ebu_oc, oh_bg, h2o2_bg, no3_bg, wet_dep, aod2d, ebu, &
    nvl, nvi, ntra, ntrb, nvl_gocart, nbands, numgas, num_ebu, num_ebu_in, &
    num_plume_data, num_soil_layers, num_chem, num_moist, num_emis_vol, &
    num_emis_ant, num_emis_dust, num_emis_seas, &
    num_asym_par, num_bscat_coef, num_ext_coef, deg_lon, deg_lat, &
    its, ite, jts, jte, kts, kte, &
    ims, ime, jms, jme, kms, kme, tile, verbose, rc)

    logical,            intent(in) :: readrestart
    integer,            intent(in) :: chem_opt
    integer,            intent(in) :: chem_in_opt
    integer,            intent(in) :: chem_conv_tr
    integer,            intent(in) :: biomass_burn_opt
    integer,            intent(in) :: seas_opt
    integer,            intent(in) :: dust_opt
    integer,            intent(in) :: dmsemis_opt
    integer,            intent(in) :: wetdep_ls_opt
    integer,            intent(in) :: call_chemistry
    integer,            intent(in) :: aer_ra_feedback
    integer,            intent(in) :: aer_ra_frq
    integer,            intent(in) :: plumerise_flag
    integer,            intent(in) :: plumerisefire_frq
    integer,            intent(in) :: kemit
    integer,            intent(in) :: ktau
    integer,            intent(in) :: current_month
    integer,            intent(in) :: julday
    integer,            intent(in) :: tz
    integer,            intent(in) :: its, ite, jts, jte, kts, kte
    integer,            intent(in) :: ims, ime, jms, jme, kms, kme
    integer,            intent(in) :: nvl, nvi, ntra, ntrb, nvl_gocart, &
                                      nbands, num_ebu, num_ebu_in, num_plume_data, &
                                      num_soil_layers, num_chem, num_moist, &
                                      num_emis_vol, num_emis_ant, &
                                      num_emis_dust, num_emis_seas
    integer,            intent(in) :: num_asym_par, num_bscat_coef, num_ext_coef
    integer,            intent(in) :: numgas
    integer,            intent(in) :: tile
    logical, optional,  intent(in) :: verbose
    integer, optional, intent(out) :: rc

    real(CHEM_KIND_R8), intent(in) :: dts

    real(CHEM_KIND_R4), dimension(nvl_gocart),       intent(in) :: p_gocart
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: clayfrac
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: dm0
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero1
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero2
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero3
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: rdrag
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: uthr
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ssm
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_dt
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_height
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_mass
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_ash_dt
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(inout) :: emiss_ash_height
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(inout) :: emiss_ash_mass
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: sandfrac
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: th_pvsrf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart),     intent(in) :: h2o2_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart),     intent(in) :: no3_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart),     intent(in) :: oh_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_emis_ant),   intent(in) :: emiss_ab
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_ebu_in),     intent(in) :: emiss_abu
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_plume_data), intent(in) :: plume

    ! -- input pointers: indexing must always start from 1
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: area
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: hf2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: pb2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: rc2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: rn2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: rsds
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: slmsk2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: snwdph2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: stype2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: ts2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: us2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: vtype2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: vfrac2d
    real(CHEM_KIND_F8), dimension(:, :), intent(in) :: zorl2d
    real(CHEM_KIND_R8), dimension(:, :), intent(in) :: deg_lat
    real(CHEM_KIND_R8), dimension(:, :), intent(in) :: deg_lon

    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: dqdt
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: exch
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: ph3d
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: phl3d
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: pr3d
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: prl3d
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: sm3d
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: tk3d
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: us3d
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: vs3d
    real(CHEM_KIND_F8), dimension(:, :, :), intent(in) :: ws3d

    real(CHEM_KIND_F8), dimension(:, :, :, :), intent(in)  :: tr3d_in
    real(CHEM_KIND_F8), dimension(:, :, :, :), intent(out) :: tr3d_out

    ! -- output diagnostics
    real(CHEM_KIND_F8), dimension(:, :, :),    intent(out) :: trcm
    real(CHEM_KIND_F8), dimension(:, :, :),    intent(out) :: trab
    real(CHEM_KIND_F8), dimension(:, :, :),    intent(out) :: truf
    real(CHEM_KIND_F8), dimension(:, :, :, :), intent(out) :: trdf

    ! -- output tracers
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem), intent(out) :: wet_dep
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: p10
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: pm25
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: ebu_oc
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: oh_bg
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: h2o2_bg
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: no3_bg
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl, 1:nbands),  intent(out) :: ext_cof
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl, 1:nbands),  intent(out) :: sscal
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl, 1:nbands),  intent(out) :: asymp
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl, ntra+ntrb), intent(out) :: trdp

    ! -- buffers
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),                     intent(inout) :: aod2d
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_ebu), intent(inout) :: ebu

    ! -- local variables

    ! -- work arrays
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: aod
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: ash_fall
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: clayf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: dep_vel_o3
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: e_co
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: dms_0
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: dusthelp
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: dxy
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: firesize_agef
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: firesize_aggr
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: firesize_agsv
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: firesize_agtf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: gsw
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: hfx
    integer,            dimension(ims:ime, jms:jme) :: isltyp
    integer,            dimension(ims:ime, jms:jme) :: ivgtyp
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: mean_fct_agef
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: mean_fct_aggr
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: mean_fct_agsv
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: mean_fct_agtf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: pbl
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: raincv_b
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: precc
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rcav
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rnav
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rmol
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: sandf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: seashelp
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: snowh
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: tcosz
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: tsk
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: ttday
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: u10
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: ust
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: v10
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: vegfra
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: xland
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: xlat
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: xlong
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: znt

    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_ebu_in)    :: ebu_in
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:3)           :: erod
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_frp_plume) :: plume_frp

    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: ac3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: ahno3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: anh3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: asulf
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: backg_h2o2
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: backg_no3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: backg_oh
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: cor3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: dz8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: dqdti
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: exch_h
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: h2o2_t
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: h2oai
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: h2oaj
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: no3_t
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: nu3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: oh_t
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: p8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: pm10
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: pm2_5_dry
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: pm2_5_dry_ec
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: relhum
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: rho_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: rri
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: t8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: t_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: u_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: v_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: vvel
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: z_at_w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: zmid
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: convfac
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: cu_co_ten
    real(CHEM_KIND_R4), dimension(ims:ime, 1:num_soil_layers, jms:jme) :: smois

    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_asym_par)   :: asym_par
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:nbands)         :: asympar
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_bscat_coef) :: bscat_coeff
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: bscoefsw
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem)       :: chem
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kemit, jms:jme, 1:num_emis_ant) :: emis_ant
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_emis_vol)   :: emis_vol
    real(CHEM_KIND_R4), dimension(ims:ime,   1:1  , jms:jme, 1:num_emis_dust)  :: emis_dust
    real(CHEM_KIND_R4), dimension(ims:ime,   1:1  , jms:jme, 1:num_emis_seas)  :: emis_seas
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_ext_coef)   :: ext_coeff
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:nbands)         :: extt
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: gaersw
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l2aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l3aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l4aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l5aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l6aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l7aer
    real(CHEM_KIND_R4), dimension(ims:ime,   1:1  , jms:jme, 1:5)              :: srce_dust
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:nbands)         :: ssca

    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem)       :: var_rmv
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem)       :: tr_fall
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem)       :: dry_fall
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem)       :: sedim
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_moist)      :: moist
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:16)             :: tauaerlw
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: tauaersw
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: waersw

    real(CHEM_KIND_R4), dimension(1:num_chem) :: ppm2ugkg

    ! -- local variables
    logical :: call_gocart, call_plume, call_radiation, scale_fire_emiss
    logical :: store_arrays

    integer :: localrc
    integer :: nbegin, nv, nvv
    integer :: i, ip, j, jp, k, kp
    integer :: ids, ide, jds, jde, kds, kde
    real(CHEM_KIND_R4) :: dt
    real(CHEM_KIND_R8) :: curr_secs

    real(CHEM_KIND_R4) :: factor, factor2, factor3
    real(CHEM_KIND_R4) :: dtstep, gmt

    real(CHEM_KIND_R4), parameter :: m2mm = 1.e+03_CHEM_KIND_R4
    real(CHEM_KIND_R4), parameter :: frpc = 1.e+06_CHEM_KIND_R4
    real(CHEM_KIND_R4), parameter :: ugkg = 1.e-09_CHEM_KIND_R4 !lzhang
    integer localpe

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    ! -- initialize output arrays
    wet_dep  = 0._CHEM_KIND_R4
    p10      = 0._CHEM_KIND_R4
    pm25     = 0._CHEM_KIND_R4
    ebu_oc   = 0._CHEM_KIND_R4
    oh_bg    = 0._CHEM_KIND_R4
    h2o2_bg  = 0._CHEM_KIND_R4
    no3_bg   = 0._CHEM_KIND_R4
    ext_cof  = 0._CHEM_KIND_R4
    sscal    = 0._CHEM_KIND_R4
    asymp    = 0._CHEM_KIND_R4
    trdp     = 0._CHEM_KIND_R4
    tr3d_out = 0._CHEM_KIND_F8

    ! -- initialize output diagnostics
    trcm = 0._CHEM_KIND_F8
    trab = 0._CHEM_KIND_F8
    truf = 0._CHEM_KIND_F8
    trdf = 0._CHEM_KIND_F8

    if (chem_opt == CHEM_OPT_NONE) return

    ! -- initialize local arrays
    aod           = 0._CHEM_KIND_R4
    ash_fall      = 0._CHEM_KIND_R4
    clayf         = 0._CHEM_KIND_R4
    dep_vel_o3    = 0._CHEM_KIND_R4
    e_co          = 0._CHEM_KIND_R4
    dms_0         = 0._CHEM_KIND_R4
    dusthelp      = 0._CHEM_KIND_R4
    dxy           = 0._CHEM_KIND_R4
    firesize_agef = 0._CHEM_KIND_R4
    firesize_aggr = 0._CHEM_KIND_R4
    firesize_agsv = 0._CHEM_KIND_R4
    firesize_agtf = 0._CHEM_KIND_R4
    gsw           = 0._CHEM_KIND_R4
    hfx           = 0._CHEM_KIND_R4
    isltyp        = 0
    ivgtyp        = 0
    mean_fct_agef = 0._CHEM_KIND_R4
    mean_fct_aggr = 0._CHEM_KIND_R4
    mean_fct_agsv = 0._CHEM_KIND_R4
    mean_fct_agtf = 0._CHEM_KIND_R4
    pbl           = 0._CHEM_KIND_R4
    raincv_b      = 0._CHEM_KIND_R4
    precc         = 0._CHEM_KIND_R4
    rcav          = 0._CHEM_KIND_R4
    rnav          = 0._CHEM_KIND_R4
    rmol          = 0._CHEM_KIND_R4
    sandf         = 0._CHEM_KIND_R4
    seashelp      = 0._CHEM_KIND_R4
    snowh         = 0._CHEM_KIND_R4
    tcosz         = 0._CHEM_KIND_R4
    tsk           = 0._CHEM_KIND_R4
    ttday         = 0._CHEM_KIND_R4
    u10           = 0._CHEM_KIND_R4
    ust           = 0._CHEM_KIND_R4
    v10           = 0._CHEM_KIND_R4
    vegfra        = 0._CHEM_KIND_R4
    xland         = 0._CHEM_KIND_R4
    xlat          = 0._CHEM_KIND_R4
    xlong         = 0._CHEM_KIND_R4
    znt           = 0._CHEM_KIND_R4
    ebu_in        = 0._CHEM_KIND_R4
    erod          = 0._CHEM_KIND_R4
    ac3           = 0._CHEM_KIND_R4
    ahno3         = 0._CHEM_KIND_R4
    anh3          = 0._CHEM_KIND_R4
    asulf         = 0._CHEM_KIND_R4
    backg_h2o2    = 0._CHEM_KIND_R4
    backg_no3     = 0._CHEM_KIND_R4
    backg_oh      = 0._CHEM_KIND_R4
    cor3          = 0._CHEM_KIND_R4
    dz8w          = 0._CHEM_KIND_R4
    dqdti         = 0._CHEM_KIND_R4
    exch_h        = 0._CHEM_KIND_R4
    h2o2_t        = 0._CHEM_KIND_R4
    h2oai         = 0._CHEM_KIND_R4
    h2oaj         = 0._CHEM_KIND_R4
    no3_t         = 0._CHEM_KIND_R4
    nu3           = 0._CHEM_KIND_R4
    oh_t          = 0._CHEM_KIND_R4
    p8w           = 0._CHEM_KIND_R4
    p_phy         = 0._CHEM_KIND_R4
    pm10          = 0._CHEM_KIND_R4
    pm2_5_dry     = 0._CHEM_KIND_R4
    pm2_5_dry_ec  = 0._CHEM_KIND_R4
    relhum        = 0._CHEM_KIND_R4
    rho_phy       = 0._CHEM_KIND_R4
    rri           = 0._CHEM_KIND_R4
    t8w           = 0._CHEM_KIND_R4
    t_phy         = 0._CHEM_KIND_R4
    u_phy         = 0._CHEM_KIND_R4
    v_phy         = 0._CHEM_KIND_R4
    vvel          = 0._CHEM_KIND_R4
    z_at_w        = 0._CHEM_KIND_R4
    zmid          = 0._CHEM_KIND_R4
    convfac       = 0._CHEM_KIND_R4
    cu_co_ten     = 0._CHEM_KIND_R4
    smois         = 0._CHEM_KIND_R4
    asym_par      = 0._CHEM_KIND_R4
    asympar       = 0._CHEM_KIND_R4
    bscat_coeff   = 0._CHEM_KIND_R4
    bscoefsw      = 0._CHEM_KIND_R4
    chem          = 0._CHEM_KIND_R4
    emis_ant      = 0._CHEM_KIND_R4
    emis_vol      = 0._CHEM_KIND_R4
    emis_dust     = 0._CHEM_KIND_R4
    emis_seas     = 0._CHEM_KIND_R4
    ext_coeff     = 0._CHEM_KIND_R4
    extt          = 0._CHEM_KIND_R4
    gaersw        = 0._CHEM_KIND_R4
    l2aer         = 0._CHEM_KIND_R4
    l3aer         = 0._CHEM_KIND_R4
    l4aer         = 0._CHEM_KIND_R4
    l5aer         = 0._CHEM_KIND_R4
    l6aer         = 0._CHEM_KIND_R4
    l7aer         = 0._CHEM_KIND_R4
    srce_dust     = 0._CHEM_KIND_R4
    ssca          = 0._CHEM_KIND_R4
    var_rmv       = 0._CHEM_KIND_R4
    tr_fall       = 0._CHEM_KIND_R4
    dry_fall      = 0._CHEM_KIND_R4
    sedim         = 0._CHEM_KIND_R4
    moist         = 0._CHEM_KIND_R4
    tauaerlw      = 0._CHEM_KIND_R4
    tauaersw      = 0._CHEM_KIND_R4
    waersw        = 0._CHEM_KIND_R4

    ! -- set domain
    ids = ims
    ide = ime
    jds = jms
    jde = jme
    kds = kms
    kde = kme

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._CHEM_KIND_R4
  !  ppm2ugkg(p_so2 ) = 1.e+03_CHEM_KIND_R4 * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_CHEM_KIND_R4 * mw_so4_aer / mwdry
!  

    ! -- set run parameters
    ! -- nbegin is the start address (-1) of the first chem variable in tr3d
    if (num_moist > 3) then
      nbegin = ntra + num_moist - 2
    else
      nbegin = ntra + num_moist - 3
    end if

    ! -- get time & time step
    dt = real(dts, kind=CHEM_KIND_R4)
    curr_secs = ktau * dts
    gmt = real(tz)
!    write(6,*)'gmt',gmt,'tz',tz,'ktau',ktau,ktau*dts/3600.
!    call flush(6)

    ! -- set control flags
    call_plume       = (biomass_burn_opt == BURN_OPT_ENABLE) .and. (plumerisefire_frq > 0)
    call chem_comm_get(localpe=localpe)
!    if(localpe==0)then
!      write(800,*)'dofirehour',dofirehour,'call_plume',call_plume
!    endif
    if(dofirehour)then
      if (call_plume) &
       call_plume    = (mod(int(curr_secs-dts), max(1, 60*plumerisefire_frq)) == 0) &
                        .or. (ktau == 1)
    else
      if (call_plume) &
       call_plume    = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0) &
                        .or. (ktau == 1)
    endif
    if(localpe==0)then
      if(call_plume)then
        write(800,*)'gocart call_plume',call_plume,'ktau',ktau,'curr_secs',curr_secs,curr_secs/3600.
      else
        write(800,*)'gocart not call pliume ',ktau,mod(int(curr_secs-dts), max(1, 60*plumerisefire_frq))
      endif
    endif
    call_gocart      = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)
    call_radiation   = (mod(int(curr_secs), max(1, 60*aer_ra_frq)) == 0) .or. (ktau == 1)
    scale_fire_emiss = .false.

    if (present(verbose)) then
      if (verbose) &
       call gocart_diag_output(ktau, plumerise_flag, plumerisefire_frq, &
         call_gocart, call_plume, call_radiation, readrestart)
    end if

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if
    if(call_plume)then
      if(localpe==0)then
        write(6,*)'gocart call_plume ktau',ktau,'curr_secs',curr_secs,'dts',dts
        call flush(6)
        write(800,*)'gocart ktau',ktau,'curr_secs',curr_secs,'dts',dts
        call flush(800)
      endif
    endif

    do j = jts, jte
      jp = j - jts + 1
      do i = its, ite
        ip = i - its + 1
        ! -- compute incremental convective and large-scale rainfall
        rcav(i,j)  = max(m2mm * rc2d(ip,jp)                , 0._CHEM_KIND_R4)
        rnav(i,j)  = max(m2mm * (rn2d(ip,jp) - rc2d(ip,jp)), 0._CHEM_KIND_R4)
      end do
    end do

    ! -- get ready for chemistry run
!    do i=1,30
!      write(6,*)'gocart',i,maxval(tr3d_in(:,:,:,i)),minval(tr3d_in(:,:,:,i))
!    end do
    call gocart_prep(readrestart,chem_opt,chem_in_opt,ktau,dt,tr3d_in,tk3d,sm3d,   &
                   ts2d,us2d,rsds,pr3d,prl3d,ph3d,phl3d,emiss_ash_mass,emiss_ash_height, &
                   emiss_ash_dt,dm0,emiss_tr_mass,emiss_tr_height,      &
                   emiss_tr_dt,snwdph2d,vfrac2d,vtype2d,stype2d,us3d,vs3d,ws3d,  &
                   slmsk2d,zorl2d,dqdt,exch,pb2d,hf2d,clayfrac,clayf,sandfrac,sandf,th_pvsrf, &
                   oh_backgd,h2o2_backgd,no3_backgd,backg_oh,backg_h2o2,backg_no3,p_gocart,   &
                   nvl_gocart,ttday,tcosz,gmt,julday,area,ero1,   &
                   ero2,ero3,rcav,raincv_b,deg_lat,deg_lon,ntra, &
                   relhum,rri,t_phy,moist,u_phy,v_phy,p_phy,chem,tsk, &
                   grvity,rd,erod,emis_ant,emis_vol,e_co,dms_0,        &
                   u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,xland,dxy, &
                   t8w,p8w,dqdti,exch_h,pbl,hfx,snowh,xlat,xlong,convfac,z_at_w,zmid,dz8w,vvel,&
                   rho_phy,smois,num_soil_layers,num_chem,num_moist,        &
                   emiss_abu,ebu_in,emiss_ab,num_ebu_in,num_plume_data,num_emis_ant,     &
                   num_emis_vol,kemit,call_gocart,plumerise_flag, &
                   plume,plume_frp, &
                   mean_fct_agtf,mean_fct_agef,mean_fct_agsv, &
                   mean_fct_aggr,firesize_agtf,firesize_agef, &
                   firesize_agsv,firesize_aggr, &
                   ppm2ugkg, &
                   ids,ide, jds,jde, kds,kde, &
                   ims,ime, jms,jme, kms,kme, &
                   its,ite, jts,jte, kts,kte, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure in gocart_prep", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- compute sea salt
    if (seas_opt /= SEAS_OPT_NONE) then
      call gocart_seasalt_driver(ktau,dt,rri,t_phy,moist, &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,ust,p8w,tsk,&
        xland,xlat,xlong,dxy,grvity,emis_seas,           &
        seashelp,num_emis_seas,num_moist,num_chem,seas_opt,&
        ids,ide, jds,jde, kds,kde,                        &
        ims,ime, jms,jme, kms,kme,                        &
        its,ite, jts,jte, kts,kte)
      truf(:,:,num_emis_dust+1:num_emis_dust+num_emis_seas) = &
        emis_seas(its:ite,1,jts:jte,1:num_emis_seas) !kg/m2/s
    endif
    store_arrays = .false.
    select case (dust_opt)
      case (DUST_OPT_AFWA)
        call gocart_dust_afwa_driver(ktau,dt,rri,t_phy,moist,u_phy,    &
          v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,ivgtyp,isltyp,&
          vegfra,xland,xlat,xlong,gsw,dxy,grvity,emis_dust,srce_dust,  &
          dusthelp,ust,znt,clayf,sandf,                                &
          num_emis_dust,num_moist,num_chem,num_soil_layers,            &
          ids,ide, jds,jde, kds,kde,                                   &
          ims,ime, jms,jme, kms,kme,                                   &
          its,ite, jts,jte, kts,kte)
        store_arrays = .true.
      case (DUST_OPT_FENGSHA)
       call gocart_dust_fengsha_driver(dt,chem,rho_phy,smois,p8w,ssm,  &
            isltyp,vegfra,snowh,xland,dxy,grvity,emis_dust,ust,znt,    &
            clayf,sandf,rdrag,uthr,                                    &
            num_emis_dust,num_moist,num_chem,num_soil_layers,          &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte)
        store_arrays = .true.
      case (DUST_OPT_GOCART)
        call gocart_dust_driver(chem_opt,dt,u_phy,v_phy,chem,rho_phy,  &
          dz8w,smois,u10,v10,p8w,erod,isltyp,xland,dxy,grvity,         &
          emis_dust,srce_dust,dusthelp,                                &
          num_emis_dust,num_chem,num_soil_layers,                      &
          current_month,                                               &
          ids,ide, jds,jde, kds,kde,                                   &
          ims,ime, jms,jme, kms,kme,                                   &
          its,ite, jts,jte, kts,kte)
        store_arrays = .true.
    end select
    store_arrays = store_arrays .and. (chem_opt >= CHEM_OPT_GOCART)

    ! -- set output arrays
    if (store_arrays) then
      truf(:,:,1:num_emis_dust) = ugkg*emis_dust(its:ite,1,jts:jte,1:num_emis_dust) !kg/m2/s
    end if
    if (call_plume) then
      call plumerise_driver (ktau,dtstep,num_chem,num_ebu,num_ebu_in, &
        ebu,ebu_in,mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr, &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr, &
        'GOCART','BIOMASSB', t_phy,moist(:,:,:,p_qv),            &
        rho_phy,vvel,u_phy,v_phy,p_phy,                          &
        z_at_w,scale_fire_emiss,plume_frp,plumerise_flag,        &
        ids,ide, jds,jde, kds,kde,                               &
        ims,ime, jms,jme, kms,kme,                               &
        its,ite, jts,jte, kts,kte                                )
    end if

    if (dmsemis_opt == DMSE_OPT_ENABLE) then
      call gocart_dmsemis(dt,rri,t_phy,u_phy,v_phy,     &
         chem,rho_phy,dz8w,u10,v10,p8w,dms_0,tsk,       &
         ivgtyp,isltyp,xland,dxy,grvity,mwdry,          &
         num_chem,p_dms,                                &
         ids,ide, jds,jde, kds,kde,                     &
         ims,ime, jms,jme, kms,kme,                     &
         its,ite, jts,jte, kts,kte)
    endif

    if ((dust_opt /= DUST_OPT_NONE) .or. &
        (seas_opt /= SEAS_OPT_NONE)) then
      call gocart_settling_driver(dt,t_phy,moist,  &
        chem,rho_phy,dz8w,p8w,p_phy,sedim, &
        dusthelp,seashelp,dxy,grvity,  &
        num_moist,num_chem,            &
        ids,ide, jds,jde, kds,kde,     &
        ims,ime, jms,jme, kms,kme,     &
        its,ite, jts,jte, kts,kte)
    end if

#if 0
    if (chem_opt == 316) then
      ! -- 10 volcanic size bins
      call vash_settling_driver(dt,t_phy,moist, &
         chem,rho_phy,dz8w,p8w,p_phy,dxy,       &
         ash_fall,grvity,num_moist,num_chem,    &
         ids,ide, jds,jde, kds,kde,             &
         ims,ime, jms,jme, kms,kme,             &
         its,ite, jts,jte, kts,kte)
      ashfall = ash_fall
    else
#endif
      ! -- 4 volcanic size bins
      call vashshort_settling_driver(dt,t_phy,moist, &
           chem,rho_phy,dz8w,p8w,p_phy,dxy,          &
           ash_fall,grvity,num_moist,num_chem,       &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte) 
!     ashfall = ash_fall !!!!!!!!!!!!!!!!!! WARNING: rewrites ashfall if chem_opt = 316

#if 0
    end if
#endif
    ! -- add biomass burning emissions at every timestep
    if (biomass_burn_opt == BURN_OPT_ENABLE) then
      jp = jte
      factor3 = 0._CHEM_KIND_R4
      select case (plumerise_flag)
        case (FIRE_OPT_MODIS)
          factor3 = 4.828e-04_CHEM_KIND_R4/60
          kp = kte    ! full column
        case (FIRE_OPT_GBBEPx)
          factor3 = 1.e-03_CHEM_KIND_R4 * mwdry / mw_so2_aer
          if (plumerisefire_frq > 0) then
            kp = kte  ! full column
          else
            kp = kts  ! surface only
          end if
        case default
          ! -- no further options available, skip this step
          jp = jts - 1
      end select

      if (kp == kts) then
        ! -- only include surface emissions
        k = kts
        do j = jts, jp
          do i = its, ite
            ! -- factor for pm emissions, factor2 for burn emissions
            factor  = dt*rri(i,k,j)/dz8w(i,k,j)
            factor2 = factor * factor3
            chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu_in(i,j,p_ebu_in_oc  )
            chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu_in(i,j,p_ebu_in_bc  )
            chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu_in(i,j,p_ebu_in_pm25)
            chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu_in(i,j,p_ebu_in_pm10)
            chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu_in(i,j,p_ebu_in_so2 )
          end do
        end do

      else
        ! -- use full-column emissions
        do j = jts, jp
          do k = kts, kp
            do i = its, ite
              ! -- factor for pm emissions, factor2 for burn emissions
              factor  = dt*rri(i,k,j)/dz8w(i,k,j)
              factor2 = factor * factor3
!              if(ebu(i,k,j,p_ebu_oc)/=0.0)then
!                write(6,*)'add ebuoc k',k,'i',i,j,ebu(i,k,j,p_ebu_oc  ),ebu(i,k,j,p_ebu_bc  )
!              endif
              chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu(i,k,j,p_ebu_oc  )
              chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu(i,k,j,p_ebu_bc  )
              chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu(i,k,j,p_ebu_pm25)
              chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu(i,k,j,p_ebu_pm10)
              chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu(i,k,j,p_ebu_so2 )
            end do
          end do
        end do
      end if

    end if
    ! -- subgrid convective transport
    if (chem_conv_tr == CTRA_OPT_GRELL) then
      call grelldrvct(dt,ktau,                  &
        rho_phy,rcav,chem,tr_fall,              &
        u_phy,v_phy,t_phy,moist,dz8w,p_phy,p8w, &
        pbl,xlv,cp,grvity,rv,z_at_w,cu_co_ten,  &
        numgas,chem_opt,                        &
        num_chem,num_moist,tile,                &
        ids,ide, jds,jde, kds,kde,              &
        ims,ime, jms,jme, kms,kme,              &
        its,ite, jts,jte, kts,kte)
     endif

     call dry_dep_driver(ktau,dt,julday,current_month,t_phy,p_phy,&
       moist,p8w,rmol,rri,gmt,t8w,rcav,                           &
       chem,rho_phy,dz8w,exch_h,hfx,                              &
       ivgtyp,tsk,gsw,vegfra,pbl,ust,znt,zmid,z_at_w,             &
       xland,xlat,xlong,h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,     &
       anh3,dry_fall,dep_vel_o3,grvity,                           &
       e_co,kemit,snowh,numgas,                                   &
       num_chem,num_moist,                                        &
       ids,ide, jds,jde, kds,kde,                                 &
       ims,ime, jms,jme, kms,kme,                                 &
       its,ite, jts,jte, kts,kte)

     ! -- ls wet deposition
     select case (wetdep_ls_opt)
       case (WDLS_OPT_GSD)
         call wetdep_ls(dt,chem,rnav,moist,rho_phy,var_rmv,xlat, &
                        num_moist,num_chem,p_qc,p_qi,dz8w,vvel, &
                        ids,ide, jds,jde, kds,kde,              &
                        ims,ime, jms,jme, kms,kme,              &
                        its,ite, jts,jte, kts,kte)
       case (WDLS_OPT_NGAC)
         call WetRemovalGOCART(its,ite, jts,jte, kts,kte, 1,1, dt, &
                               num_chem,var_rmv,chem,p_phy,t_phy,  &
                               rho_phy,dqdti,rcav,rnav,            &
                               ims,ime, jms,jme, kms,kme, localrc)
         if (chem_rc_check(localrc, msg="Failure in NGAC wet removal scheme", &
           file=__FILE__, line=__LINE__, rc=rc)) return
       case default
         ! -- no further option implemented
    end select

    if (call_gocart) then
      call gocart_chem_driver(ktau,dt,dtstep,gmt,julday,    &
           t_phy,moist,chem,rho_phy,dz8w,p8w,backg_oh,oh_t, &
           backg_h2o2,h2o2_t,backg_no3,no3_t,               &
           dxy,grvity,xlat,xlong,ttday,tcosz,               &
           chem_opt,num_chem,num_moist,                     &
           ids,ide, jds,jde, kds,kde,                       &
           ims,ime, jms,jme, kms,kme,                       &
           its,ite, jts,jte, kts,kte                        )
      call gocart_aerosols_driver(ktau,dtstep,t_phy,moist,  &
           chem,rho_phy,dz8w,p8w,dxy,grvity,                &
           chem_opt,num_chem,num_moist,                     &
           ids,ide, jds,jde, kds,kde,                       &
           ims,ime, jms,jme, kms,kme,                       &
           its,ite, jts,jte, kts,kte                        )
    endif

    if (call_radiation) then
      store_arrays = .false.
      select case (aer_ra_feedback)
        case (1) 
          call optical_driver(curr_secs,dtstep,             &
               chem,dz8w,rri,relhum,                        &
               h2oai,h2oaj,                                 &
               tauaersw,gaersw,waersw,bscoefsw,tauaerlw,    &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,         &
               num_chem,chem_opt,ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme,                   &
               its,ite, jts,jte, kts,kte)
          call aer_opt_out(aod,dz8w,                        &
               ext_coeff,bscat_coeff,asym_par,              &
               tauaersw,gaersw,waersw,tauaerlw,             &
               num_ext_coef,num_bscat_coef,num_asym_par,    &
               ids,ide, jds,jde, kds,kde,                   &
               ims,ime, jms,jme, kms,kme,                   &
               its,ite, jts,jte, kts,kte)
          call aer_ra(dz8w                                  &
               ,extt,ssca,asympar,nbands                    &
               ,tauaersw,gaersw,waersw,tauaerlw             &
               ,ids,ide, jds,jde, kds,kde                   &
               ,ims,ime, jms,jme, kms,kme                   &
               ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case (2) 
          call aero_opt('sw',dz8w,chem         &
                   ,rri,relhum,aod             &
                   ,extt,ssca,asympar,num_chem &
                   ,ids,ide, jds,jde, kds,kde  &
                   ,ims,ime, jms,jme, kms,kme  &
                   ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case default
          ! -- no feedback
      end select
      
      ! -- store aerosol optical variables for feedback in radiation
      if (store_arrays) then
        do nv = 1, nbands
          do j = jts, jte
            do k = kts, kte
              do i = its, ite
                ext_cof(i,j,k,nv) = extt   (i,k,j,nv)
                sscal  (i,j,k,nv) = ssca   (i,k,j,nv)
                asymp  (i,j,k,nv) = asympar(i,k,j,nv)
              end do
            end do
          end do
        end do
        aod2d(its:ite,jts:jte) = aod(its:ite,jts:jte)
      end if
    endif

    call sum_pm_gocart (                              &
         rri, chem,pm2_5_dry, pm2_5_dry_ec, pm10,     &
         num_chem,chem_opt,                           &
         ids,ide, jds,jde, kds,kde,                   &
         ims,ime, jms,jme, kms,kme,                   &
         its,ite, jts,jte, kts,kte)

    ! -- pm25 and pm10 for output , not for tracer options
    do j = jts, jte
      do k = kts, kte
        do i = its, ite
          pm25  (i,j,k) = pm2_5_dry(i,k,j)
          p10   (i,j,k) = pm10     (i,k,j)
          ebu_oc(i,j,k) = ebu      (i,k,j,p_ebu_oc)
        end do
      end do
    end do

    if (call_gocart) then
      do j = jts, jte
        do k = kts, kte
          do i = its, ite
            oh_bg  (i,j,k) = max(0., oh_t  (i,k,j))
            h2o2_bg(i,j,k) = max(0., h2o2_t(i,k,j))
            no3_bg (i,j,k) = max(0., no3_t (i,k,j))
          end do
        end do
      end do
    end if

    ! -- put chem stuff back into tracer array
    do nv = 1, num_chem
      nvv = nbegin + nv
      jp = 0
      do j = jts, jte
        jp = jp + 1
        kp = 0
        do k = kts, kte
          kp = kp + 1
          ip = 0
          do i = its, ite
            ip = ip + 1
            ! -- export updated tracers
            tr3d_out(ip,jp,kp,nvv) = real(ppm2ugkg(nv) * max(epsilc,chem(i,k,j,nv)), kind=CHEM_KIND_F8)
            ! -- compute auxiliary array trdp
            trdp(i,j,k,nvv) = tr3d_out(ip,jp,kp,nvv)*(pr3d(ip,jp,kp)-pr3d(ip,jp,kp+1))
          end do
        end do
      end do
    end do

    ! -- calculate column mass density
    call gocart_diag_cmass(chem_opt, nbegin, grvity, pr3d, tr3d_out, trcm)

    ! -- output anthropogenic emissions, kg/m2/s
    trab(:,:,1) = ugkg*emis_ant(its:ite, kts, jts:jte, p_e_bc )
    trab(:,:,2) = ugkg*emis_ant(its:ite, kts, jts:jte, p_e_oc )
    trab(:,:,3) = ugkg*emis_ant(its:ite, kts, jts:jte, p_e_so2)

    ! -- output biomass burning emissions,kg/m2/s
    trab(:,:,4) = ugkg*ebu_in(its:ite, jts:jte, p_ebu_in_bc )
    trab(:,:,5) = ugkg*ebu_in(its:ite, jts:jte, p_ebu_in_oc )
    trab(:,:,6) = ugkg*ebu_in(its:ite, jts:jte, p_ebu_in_so2)

    ! -- output aod !lzhang
    trab(:,:,7) = aod2d(its:ite, jts:jte)

    ! -- output sedimentation and dry/wet deposition
    call gocart_diag_store(1, sedim, trdf)
    ! -- output dry deposition
    call gocart_diag_store(2, dry_fall, trdf)
    ! -- output large-scale wet deposition
    call gocart_diag_store(3, var_rmv, trdf)
    ! -- output convective-scale wet deposition
    if (chem_conv_tr == CTRA_OPT_GRELL) then
      where (tr_fall > 0._CHEM_KIND_R4) wet_dep = tr_fall
      call gocart_diag_store(4, tr_fall, trdf)
    end if

  end subroutine gocart_advance

  subroutine gocart_diag_output(ktau, plumerise_flag, plumerisefire_frq, &
    call_gocart, call_plume, call_radiation, readrestart)

    ! -- arguments
    integer, intent(in) :: ktau
    integer, intent(in) :: plumerise_flag
    integer, intent(in) :: plumerisefire_frq
    logical, intent(in) :: call_gocart
    logical, intent(in) :: call_plume
    logical, intent(in) :: call_radiation
    logical, intent(in) :: readrestart

    ! -- local parameters
    character(len=3), dimension(0:1), parameter :: switch = (/"OFF", "ON "/)

    ! -- local variables
    integer       :: state
    logical, save :: firstfire = .true.

    ! -- begin
    write(6,'(38("-"))')
    if (readrestart) then
      write(6,'(1x,"Restarting GSDCHEM @ step: ",i0)') ktau
    else
      write(6,'(1x,"Running GSDCHEM @ step: ",i0)') ktau
    end if
    write(6,'(14("-")," Modules: ",14("-"))')
    state = 0
    if (call_gocart) state = 1
    write(6,'(1x,"Chemistry and Aerosols",11x,a)') switch(state)
    state = 0
    if (call_radiation) state = 1
    write(6,'(1x,"Radiation",24x,a)') switch(state)
    state = 0
    if (call_plume) state = 1
    write(6,'(1x,"Fires    ",24x,a)') switch(state)
    if (call_plume) then
      if (firstfire) then
        write(6,'(1x,"* Initializing...")')
        firstfire = .false.
      end if
      select case (plumerise_flag)
        case (FIRE_OPT_MODIS)
          if (plumerisefire_frq > 0) write(6,'(1x,"* Using MODIS with plume-rise")')
        case (FIRE_OPT_GBBEPx)
          if (plumerisefire_frq > 0) then
            write(6,'(1x,"* Using GBBEPx + FRP with plume-rise")')
          else
            write(6,'(1x,"* Using surface GBBEPx only")')
          end if
        case default
          ! -- no further option
      end select
    end if
    write(6,'(38("-"))')

  end subroutine gocart_diag_output

end module gocart_mod
