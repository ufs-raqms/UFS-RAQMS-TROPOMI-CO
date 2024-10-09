module seas_mod

  use chem_const_mod,   only : pi => pi_r8
  use chem_rc_mod,      only : chem_rc_test
  use chem_types_mod,   only : CHEM_KIND_R8
  use chem_tracers_mod, only : p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_seas_5, &
                               p_eseas1, p_eseas2, p_eseas3, p_eseas4, p_eseas5, &
                               config => chem_config
  use chem_config_mod,  only : SEAS_OPT_GOCART, SEAS_OPT_NGAC

  use seas_data_mod
  use seas_ngac_mod

  implicit none

  private

  public :: number_ss_bins
  public :: emission_scale, emission_scheme
  public :: gocart_seasalt_driver

CONTAINS

  subroutine gocart_seasalt_driver(ktau,dt,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,u10,v10,ustar,p8w,tsk,            &
         xland,xlat,xlong,area,g,emis_seas, &
         seashelp,num_emis_seas,num_moist,num_chem,seas_opt,  &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )

     INTEGER,      INTENT(IN   ) :: ktau,num_emis_seas,num_moist,num_chem,   &
                                    ids,ide, jds,jde, kds,kde,               &
                                    ims,ime, jms,jme, kms,kme,               &
                                    its,ite, jts,jte, kts,kte,seas_opt
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
           INTENT(IN ) ::                                   moist
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
           INTENT(INOUT ) ::                                   chem
     REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_seas),                    &
           INTENT(OUT   ) ::                                                 &
           emis_seas
     REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
            INTENT(IN   ) ::                                                 &
                                                       u10,                  &
                                                       v10,                  &
                                                       ustar,tsk,            &
                                                       xland,                &
                                                       xlat,                 &
                                                       xlong,area
     REAL,  DIMENSION( ims:ime , jms:jme ),                        &
            INTENT(OUT   ) :: seashelp
     REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
            INTENT(IN   ) ::                                                 &
                                                          alt,               &
                                                        t_phy,               &
                                                       dz8w,p8w,             &
                                                u_phy,v_phy,rho_phy

    REAL, INTENT(IN   ) :: dt,g
!
! local variables
!
    integer :: ipr,i,j,imx,jmx,lmx,n,rc
    integer,dimension (1,1) :: ilwi
    real               :: fsstemis, memissions, nemissions, tskin_c, ws10m
    real(CHEM_KIND_R8) :: delp
    real(CHEM_KIND_R8), DIMENSION (number_ss_bins) :: tc,bems
    real(CHEM_KIND_R8), dimension (1,1) ::w10m,airmas,tskin
    real(CHEM_KIND_R8), dimension (1) :: dxy

    real(CHEM_KIND_R8), dimension(1,1,1) :: airmas1
    real(CHEM_KIND_R8), dimension(1,1,1,number_ss_bins) :: tc1
    real(CHEM_KIND_R8), dimension(1,1,number_ss_bins) :: bems1
    
!
! local parameters
!
    real(CHEM_KIND_R8), parameter :: conver  = 1.e-9_CHEM_KIND_R8
    real(CHEM_KIND_R8), parameter :: converi = 1.e+9_CHEM_KIND_R8
!
! number of dust bins
!
    imx=1
    jmx=1
    lmx=1

    emis_seas = 0.

    select case (config % chem_opt)

      case (304, 316, 317)

        seashelp(:,:)=0.
        do j=jts,jte
          do i=its,ite
!
! donṫ do dust over water!!!
!
            if(xland(i,j).lt.0.5)then
              ilwi(1,1)=0
              tc(1)=chem(i,kts,j,p_seas_1)*conver
              tc(2)=1.e-30_CHEM_KIND_R8
              tc(3)=chem(i,kts,j,p_seas_2)*conver
              tc(4)=1.e-30_CHEM_KIND_R8
              w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
              tskin(1,1)=tsk(i,j)
              delp = p8w(i,kts,j)-p8w(i,kts+1,j)
              airmas(1,1)=area(i,j) * delp / g
!
! we donṫ trust the u10,v10 values, is model layers are very thin near surface
!
              if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
!
              dxy(1)=area(i,j)
              ipr=0

              airmas1(1,1,1) = airmas(1,1)
              tc1(1,1,1,:) = tc
              bems1(1,1,:) = bems
              call source_ss( imx, jmx, lmx, number_ss_bins, dt, tc1,ilwi, dxy, w10m, airmas1, bems1,ipr)
              tc = tc1(1,1,1,:)
              chem(i,kts,j,p_seas_1)=(tc(1)+.75*tc(2))*converi
              chem(i,kts,j,p_seas_2)=(tc(3)+.25*tc(2))*converi
              seashelp(i,j)=tc(2)*converi
            endif
          enddo
        enddo

      case default

        select case (seas_opt)
          case (SEAS_OPT_GOCART)
            ! -- original GOCART sea salt scheme
            do j = jts, jte
              do i = its, ite

                ! -- only use sea salt scheme over water
                if (xland(i,j) < 0.5) then

                  ! -- compute auxiliary variables
                  delp = p8w(i,kts,j)-p8w(i,kts+1,j)
                  if (dz8w(i,kts,j) < 12.) then
                    w10m = sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
                  else
                    w10m = sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
                  end if

                  ilwi(1,1)=0
                  tc = 0.
                  tskin(1,1)=tsk(i,j)
                  airmas(1,1)=area(i,j) * delp / g
                  dxy(1)=area(i,j)
                  ipr=0

                  airmas1(1,1,1) = airmas(1,1)
                  tc1(1,1,1,:) = tc
                  bems1(1,1,:) = bems
                  call source_ss( imx,jmx,lmx,number_ss_bins, dt, tc1, ilwi, dxy, w10m, airmas1, bems1,ipr)
                  tc   = tc1(1,1,1,:)
                  bems = bems1(1,1,:)

                  ! -- add sea salt emission increments to existing airborne concentrations
                  chem(i,kts,j,p_seas_1) = chem(i,kts,j,p_seas_1) + tc(1)*converi
                  chem(i,kts,j,p_seas_2) = chem(i,kts,j,p_seas_2) + tc(2)*converi
                  chem(i,kts,j,p_seas_3) = chem(i,kts,j,p_seas_3) + tc(3)*converi
                  chem(i,kts,j,p_seas_4) = chem(i,kts,j,p_seas_4) + tc(4)*converi
                  chem(i,kts,j,p_seas_5) = chem(i,kts,j,p_seas_5) + tc(5)*converi

                  ! for output diagnostics
                  emis_seas(i,1,j,p_eseas1) = bems(1)
                  emis_seas(i,1,j,p_eseas2) = bems(2)
                  emis_seas(i,1,j,p_eseas3) = bems(3)
                  emis_seas(i,1,j,p_eseas4) = bems(4)
                  emis_seas(i,1,j,p_eseas5) = bems(5)

                end if

              end do
            end do

          case (SEAS_OPT_NGAC)
            ! -- NGAC sea salt scheme
            do j = jts, jte
              do i = its, ite

                ! -- only use sea salt scheme over water
                if (xland(i,j) < 0.5) then

                  ! -- compute auxiliary variables
                  delp = p8w(i,kts,j)-p8w(i,kts+1,j)
                  if (dz8w(i,kts,j) < 12.) then
                    ws10m = sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
                  else
                    ws10m = sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
                  end if

                  ! -- compute NGAC SST correction
                  tskin_c  = tsk(i,j) - 273.15
                  tskin_c  = min(max(tskin_c, -0.1), 36.0)    ! temperature range (0, 36) C

                  fsstemis = -1.107211 &
                             - tskin_c*(0.010681+0.002276*tskin_c) &
                             + 60.288927/(40.0 - tskin_c)
                  fsstemis = min(max(fsstemis, 0.0), 7.0)

                  do n = 1, number_ss_bins
                    memissions = 0.
                    nemissions = 0.
                    call SeasaltEmission( ra(n), rb(n), emission_scheme, &
                                          ws10m, ustar(i,j), memissions, nemissions, rc )
                    if (chem_rc_test((rc /= 0), msg="Error in NGAC sea salt scheme", &
                      file=__FILE__, line=__LINE__)) return

                    bems(n) = emission_scale(n) * fsstemis * memissions
                    tc(n) = bems(n) * dt * g / delp
                  end do

                  ! -- add sea salt emission increments to existing airborne concentrations
                  chem(i,kts,j,p_seas_1) = chem(i,kts,j,p_seas_1) + tc(1)*converi
                  chem(i,kts,j,p_seas_2) = chem(i,kts,j,p_seas_2) + tc(2)*converi
                  chem(i,kts,j,p_seas_3) = chem(i,kts,j,p_seas_3) + tc(3)*converi
                  chem(i,kts,j,p_seas_4) = chem(i,kts,j,p_seas_4) + tc(4)*converi
                  chem(i,kts,j,p_seas_5) = chem(i,kts,j,p_seas_5) + tc(5)*converi

                  ! for output diagnostics kg/m2/s
                  emis_seas(i,1,j,p_eseas1) = bems(1)
                  emis_seas(i,1,j,p_eseas2) = bems(2)
                  emis_seas(i,1,j,p_eseas3) = bems(3)
                  emis_seas(i,1,j,p_eseas4) = bems(4)
                  emis_seas(i,1,j,p_eseas5) = bems(5)
                end if

              end do
            end do

          case default
          ! -- no sea salt scheme

        end select

    end select

  end subroutine gocart_seasalt_driver

  SUBROUTINE source_ss(imx,jmx,lmx,nmx, dt1, tc, &
                       ilwi, dxy, w10m, airmas, &
                       bems,ipr)

! ****************************************************************************
! *  Evaluate the source of each seasalt particles size classes  (kg/m3) 
! *  by soil emission.
! *  Input:
! *         SSALTDEN  Sea salt density                               (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *      
! *  Output:
! *         DSRC      Source of each sea salt bins       (kg/timestep/cell) 
! *
! *
! * Number flux density: Original formula by Monahan et al. (1986) adapted
! * by Sunling Gong (JGR 1997 (old) and GBC 2003 (new)).  The new version is
! * to better represent emission of sub-micron sea salt particles.
!
! * dFn/dr = c1*u10**c2/(r**A) * (1+c3*r**c4)*10**(c5*exp(-B**2))
! * where B = (b1 -log(r))/b2
! * see c_old, c_new, b_old, b_new below for the constants.
! * number fluxes are at 80% RH.
! *
! * To calculate the flux:
! * 1) Calculate dFn based on Monahan et al. (1986) and Gong (2003)
! * 2) Assume that wet radius r at 80% RH = dry radius r_d *frh
! * 3) Convert particles flux to mass flux :
! *    dFM/dr_d = 4/3*pi*rho_d*r_d^3 *(dr/dr_d) * dFn/dr
! *             = 4/3*pi*rho_d*r_d^3 * frh * dFn/dr
! *               where rho_p is particle density [kg/m3]
! *    The factor 1.e-18 is to convert in micro-meter r_d^3
! ****************************************************************************
   

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: nmx,imx,jmx,lmx,ipr
    INTEGER, INTENT(IN)    :: ilwi(imx,jmx)
    REAL(CHEM_KIND_R8),    INTENT(IN)    :: dxy(jmx), w10m(imx,jmx)
    REAL(CHEM_KIND_R8),    INTENT(IN)    :: airmas(imx,jmx,lmx)
    REAL(CHEM_KIND_R8),    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
    REAL(CHEM_KIND_R8),    INTENT(OUT)   :: bems(imx,jmx,nmx)

    REAL(CHEM_KIND_R8) :: c0(5), b0(2)
!  REAL(CHEM_KIND_R8), PARAMETER :: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/) 
!  REAL(CHEM_KIND_R8), PARAMETER :: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
    ! Change suggested by MC
    REAL(CHEM_KIND_R8), PARAMETER :: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/) 
    REAL(CHEM_KIND_R8), PARAMETER :: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
    REAL(CHEM_KIND_R8), PARAMETER :: b_old(2)=(/0.380, 0.650/)
    REAL(CHEM_KIND_R8), PARAMETER :: b_new(2)=(/0.433, 0.433/)
    REAL(CHEM_KIND_R8), PARAMETER :: dr=5.0D-2 ! um   
    REAL(CHEM_KIND_R8), PARAMETER :: theta=30.0
    ! Swelling coefficient frh (d rwet / d rd)
!!!  REAL(CHEM_KIND_R8),    PARAMETER :: frh = 1.65
    REAL(CHEM_KIND_R8),    PARAMETER :: frh = 2.d0
    LOGICAL, PARAMETER :: old=.TRUE., new=.FALSE.
    REAL(CHEM_KIND_R8) :: rho_d, r0, r1, r, r_w, a, b, dfn, r_d, dfm, src
    INTEGER :: i, j, n, nr, ir
    REAL :: dt1,fudge_fac


    REAL(CHEM_KIND_R8)    :: tcmw(nmx), ar(nmx), tcvv(nmx)
    REAL(CHEM_KIND_R8)    :: ar_wetdep(nmx), kc(nmx)
    CHARACTER(LEN=20)     :: tcname(nmx), tcunits(nmx)
    LOGICAL               :: aerosol(nmx)


    REAL(CHEM_KIND_R8) :: tc1(imx,jmx,lmx,nmx)
    REAL(CHEM_KIND_R8), TARGET :: tcms(imx,jmx,lmx,nmx) ! tracer mass (kg; kgS for sulfur case)
    REAL(CHEM_KIND_R8), TARGET :: tcgm(imx,jmx,lmx,nmx) ! g/m3

    !-----------------------------------------------------------------------  
    ! sea salt specific
    !-----------------------------------------------------------------------  
! REAL(CHEM_KIND_R8), DIMENSION(nmx) :: ra, rb
! REAL(CHEM_KIND_R8) :: ch_ss(nmx,12)

    !-----------------------------------------------------------------------  
    ! emissions (input)
    !-----------------------------------------------------------------------  
    REAL(CHEM_KIND_R8) :: e_an(imx,jmx,2,nmx), e_bb(imx,jmx,nmx), &
            e_ac(imx,jmx,lmx,nmx)

    !-----------------------------------------------------------------------  
    ! diagnostics (budget)
    !-----------------------------------------------------------------------
!  ! tendencies per time step and process
!  REAL(CHEM_KIND_R8), TARGET :: bems(imx,jmx,nmx), bdry(imx,jmx,nmx), bstl(imx,jmx,nmx)
!  REAL(CHEM_KIND_R8), TARGET :: bwet(imx,jmx,nmx), bcnv(imx,jmx,nmx)!

!  ! integrated tendencies per process
!  REAL(CHEM_KIND_R8), TARGET :: tems(imx,jmx,nmx), tstl(imx,jmx,nmx)
!  REAL(CHEM_KIND_R8), TARGET :: tdry(imx,jmx,nmx), twet(imx,jmx,nmx), tcnv(imx,jmx,nmx)

    ! global mass balance per time step 
    REAL(CHEM_KIND_R8) :: tmas0(nmx), tmas1(nmx)
    REAL(CHEM_KIND_R8) :: dtems(nmx), dttrp(nmx), dtdif(nmx), dtcnv(nmx)
    REAL(CHEM_KIND_R8) :: dtwet(nmx), dtdry(nmx), dtstl(nmx)
    REAL(CHEM_KIND_R8) :: dtems2(nmx), dttrp2(nmx), dtdif2(nmx), dtcnv2(nmx)
    REAL(CHEM_KIND_R8) :: dtwet2(nmx), dtdry2(nmx), dtstl2(nmx)

    ! detailed integrated budgets for individual emissions
    REAL(CHEM_KIND_R8), TARGET :: ems_an(imx,jmx,nmx),    ems_bb(imx,jmx,nmx), ems_tp(imx,jmx)
    REAL(CHEM_KIND_R8), TARGET :: ems_ac(imx,jmx,lmx,nmx)
    REAL(CHEM_KIND_R8), TARGET :: ems_co(imx,jmx,nmx)

    ! executable statements
! decrease seasalt emissions (Colarco et al. 2010)
!
    !fudge_fac= 1. !.5
    !fudge_fac= .5 !lzhang
    fudge_fac= .25 !lzhang
!
    DO n = 1,nmx
       bems(:,:,n) = 0.0
       rho_d = den_seas(n)
       r0 = ra(n)*frh
       r1 = rb(n)*frh
       r = r0
       nr = INT((r1-r0)/dr+.001)
       DO ir = 1,nr
          r_w = r + dr*0.5
          r = r + dr
          IF (new) THEN
             a = 4.7*(1.0 + theta*r_w)**(-0.017*r_w**(-1.44))
             c0 = c_new
             b0 = b_new
          ELSE
             a = 3.0
             c0 = c_old
             b0 = b_old
          END IF
          !
          b = (b0(1) - LOG10(r_w))/b0(2)
          dfn = (c0(1)/r_w**a)*(1.0 + c0(3)*r_w**c0(4))* &
               10**(c0(5)*EXP(-(b**2)))
          
          r_d = r_w/frh*1.0D-6  ! um -> m
          dfm = 4.0/3.0*pi*r_d**3*rho_d*frh*dfn*dr*dt1 ! 3600 !dt1
          DO i = 1,imx
             DO j = 1,jmx
!              IF (water(i,j) > 0.0) THEN
                IF (ilwi(i,j) == 0) THEN
!                 src = dfm*dxy(j)*water(i,j)*w10m(i,j)**c0(2)
                   src = dfm*dxy(j)*w10m(i,j)**c0(2)
!                 src = ch_ss(n,dt(1)%mn)*dfm*dxy(j)*w10m(i,j)**c0(2)
                   tc(i,j,1,n) = tc(i,j,1,n) + fudge_fac*src/airmas(i,j,1)
                ELSE
                   src = 0.0
                END IF
                bems(i,j,n) = bems(i,j,n) + src*fudge_fac/(dxy(j)*dt1) !kg/m2/s
             END DO  ! i
          END DO ! j
       END DO ! ir
    END DO ! n

  END SUBROUTINE source_ss

end module seas_mod
