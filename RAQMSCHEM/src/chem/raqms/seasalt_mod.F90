module seasalt_mod
  IMPLICIT NONE
  integer, parameter :: nmx=5
  PUBLIC
  integer, parameter :: idss1=1,idss2=2,idss3=3,idss4=4
 
  INTEGER, PARAMETER :: msub = 4
  REAL, DIMENSION(nmx) :: ssaltden, ssaltreff, ra, rb
  REAL :: ch_ss(nmx,12)
  REAL,    PARAMETER :: pi = 3.1415926
  REAL,    PARAMETER :: avogad   = 6.023E23 ! new value: 6.0221415E23
  REAL,    PARAMETER :: airmw    = 28.97
  real sorcaero(nmx)
  contains
  SUBROUTINE get_aerosource_ss(i,j,ilwi,us3d,vs3d)
!  use raqmschemcomm_mod, only : w10
  use raqmschem_pmgrid_mod, only : tile,iam,iamprn
  IMPLICIT NONE
  integer ilwi
  integer i,j
  real tc(nmx),bems(nmx),tcmw(nmx)
  integer,parameter :: idum=1,jdum=1,ldum=1
  integer,parameter :: ndtdum=1
  real,dimension(1,1),parameter :: airmas=1
  real,dimension(1),parameter :: dxydum=1.
  real,dimension(1,1) :: w10m
  real*4 us3d,vs3d
  w10m=sqrt(us3d**2.+vs3d**2.)
!  w10m=w10(i,j)
!  write(6,*)'get aerosource ss',i,j,ilwi,'w10m',w10m
  call set_bcond_ss
  tc=0.0
  call source_ss(idum,jdum,ldum,nmx,ndtdum,tc, &
         ilwi,dxydum,w10m,airmas,bems)
!  Convert from kg/m2/s to molec/cm2/sec for raqms: 
!    (1 m2/1e4 cm2)*(1e3 g/1 kg)*(1 mol/MW g)*(avogad molec/mol)
!   = avogad/(10*MW)
!
  tcmw=airmw 
  sorcaero(idss1) = tc(1)*avogad/(10.*tcmw(idss1))
  sorcaero(idss2) = tc(2)*avogad/(10.*tcmw(idss2))
  sorcaero(idss3) = tc(3)*avogad/(10.*tcmw(idss3))
  sorcaero(idss4) = tc(4)*avogad/(10.*tcmw(idss4)) 
  if(tc(1)<0.0)then
    write(6,*)'tc neg',tc(1),'idss1',idss1,'w10m',w10m,'ilwi',ilwi
    call flush(6)
  endif

!  if(iam.eq.iamprn.and.ilwi.eq.2)then
!    write(6,*)'i',i,j,'ilwi',ilwi,'tc1',tc(1),'sorcaero',idss1,sorcaero(idss1)
!  endif
  end SUBROUTINE get_aerosource_ss

  SUBROUTINE set_bcond_ss


  IMPLICIT NONE



     ! Density of sea salt (kg/m3)

     ssaltden(1:nmx) = 2200.0
     
     ! Main effective radius (m)

!     ssaltreff(1) = 0.26E-6
!     ssaltreff(2) = 1.19E-6
!!     ssaltreff(3) = 2.43E-6
!     ! probably a typo ...
!     ssaltreff(3) = 3.43E-6
!     ssaltreff(4) = 7.57E-6

     ssaltreff(1) = 0.30E-6
     ssaltreff(2) = 1.00E-6
     ssaltreff(3) = 3.25E-6
     ssaltreff(4) = 7.50E-6
!    ajl
     ssaltreff(5) = 15.00E-6

     ! Bin edges in um:

     ra(1) = 0.1
     ra(2) = 0.5
     ra(3) = 1.5
     ra(4) = 5.0
!    ajl
     ra(5)=10.0

     rb(1) = 0.5
     rb(2) = 1.5
     rb(3) = 5.0
     rb(4) = 10.0
!    ajl
     rb(5)=20.0
     
     ch_ss(:,:) = &
          1.0
 !         RESHAPE(SOURCE=ch_ss_aerocom_1(:),  SHAPE=(/4,12/)) / &
 !         RESHAPE(SOURCE=ch_ss_geos4_mod_1(:),SHAPE=(/4,12/)) 


  END SUBROUTINE set_bcond_ss
  SUBROUTINE source_ss(imx, jmx, lmx, nmx, ndt1, tc, &
                     ilwi, dxy, w10m, airmas, &
                     bems)

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

!raqms Following variable is not even used
!raqms  USE mo_time_control, ONLY: dt

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: imx, jmx, lmx, nmx, ndt1
!  INTEGER, INTENT(IN)    :: ilwi(imx,jmx)
  INTEGER, INTENT(IN)    :: ilwi
  REAL,    INTENT(IN)    :: dxy(jmx), w10m(imx,jmx)
  REAL,    INTENT(IN)    :: airmas(imx,jmx,lmx)
  REAL,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL,    INTENT(OUT)   :: bems(imx,jmx,nmx)

  REAL :: c0(5), b0(2)
!  REAL, PARAMETER :: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/) 
!  REAL, PARAMETER :: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
  ! Change suggested by MC
  REAL, PARAMETER :: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/) 
  REAL, PARAMETER :: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
  REAL, PARAMETER :: b_old(2)=(/0.380, 0.650/)
  REAL, PARAMETER :: b_new(2)=(/0.433, 0.433/)
  REAL, PARAMETER :: dr=5.0E-2 ! um   
  REAL, PARAMETER :: theta=30.0
  ! Swelling coefficient frh (d rwet / d rd)
!!!  REAL,    PARAMETER :: frh = 1.65
  REAL,    PARAMETER :: frh = 2.0
  LOGICAL, PARAMETER :: old=.FALSE., new=.TRUE.
  REAL :: rho_d, r0, r1, r, r_w, a, b, dfn, r_d, dfm, src
  INTEGER :: i, j, n, nr, ir

  ! executable statements

  DO n = 1,nmx
     bems(:,:,n) = 0.0
     rho_d = ssaltden(n)
     r0 = ra(n)*frh
     r1 = rb(n)*frh
     r = r0
     nr = INT((r1-r0)/dr)
!     write(6,*)'nr',nr,'rho',rho_d
!     call flush(6)
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
        
        r_d = r_w/frh*1.0E-6  ! um -> m
        dfm = 4.0/3.0*pi*r_d**3*rho_d*frh*dfn*dr*ndt1
        DO i = 1,imx
           DO j = 1,jmx
!              IF (ilwi(i,j) == 0) THEN
              IF (ilwi == 0) THEN
                 src = dfm*dxy(j)*w10m(i,j)**c0(2)
!                 write(6,*)'dfm',dfm,'dxy',dxy(j),'w10',w10m(i,j),'c0',c0(2),'airmas',airmas(i,j,1)
!                 src = ch_ss(n,dt(1)%mn)*dfm*dxy(j)*w10m(i,j)**c0(2)
                 tc(i,j,1,n) = tc(i,j,1,n) + src/airmas(i,j,1)
              ELSE
                 src = 0.0
              END IF
              bems(i,j,n) = bems(i,j,n) + src
           END DO  ! i
        END DO ! j
     END DO ! ir
  END DO ! n

  END SUBROUTINE source_ss


end  module seasalt_mod
