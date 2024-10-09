module wetdep_alpha_mod
  use chem_types_mod
  implicit none
  real(CHEM_KIND_R8), dimension(:,:,:),allocatable :: cldice,cldliq,c_h2o
  private
  public initwetdep,solalpha

contains
   subroutine initwetdep(its,ite,jts,jte,kts,kte,kms,kme,tkin,pmid)
   use chem_types_mod
   implicit none
   integer its,ite,jts,jte,kts,kte,i,j,k,kms,kme
   real(CHEM_KIND_R4), dimension(its:ite,kms:kme,jts:jte), intent(in) :: tkin,pmid
   real(CHEM_KIND_R8) :: tc,pl,tk
!   write(6,*)'initwetdep kts',kts,kte,'kms',kms,kme,'ub tk',ubound(tkin)
!   call flush(6)
   if(.not.allocated(cldice))then
   allocate (cldice(its:ite,kts:kte,jts:jte),cldliq(its:ite,kts:kte,jts:jte))
   allocate (c_h2o(its:ite,kts:kte,jts:jte))
   endif
   do k=kts,kte
     do j=jts,jte
       do i=its,ite
        !==============================================================
         ! CLDLIQ, the cloud liquid water content [cm3 H2O/cm3 air], 
         ! is a function of the local Kelvin temperature:
         !    
         !    CLDLIQ = 2e-6                    [     T >= 268 K    ]    
         !    CLDLIQ = 2e-6 * ((T - 248) / 20) [ 248 K < T < 268 K ]
         !    CLDLIQ = 0                       [     T <= 248 K    ]    
         !==============================================================
         tk=tkin(i,k,j)
         TC=tk-273.15d0
         pl=pmid(i,k,j)
         IF ( TK >= 268d0 ) THEN 
            CLDLIQ(I,k,j) = 2d-6 

         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN 
            CLDLIQ(I,k,j) = 2d-6 * ( ( TK - 248d0 ) / 20d0 )

         ELSE 
            CLDLIQ(I,k,j) = 0d0
     
         ENDIF
     
         !=============================================================
         ! CLDICE, the cloud ice content [cm3 ice/cm3 air] is given by:
         !    
         !    CLDICE = 2e-6 - CLDLIQ
         !=============================================================
         CLDICE(I,k,j) = 2d-6 - CLDLIQ(I,k,j)
         !=============================================================
         ! C_H2O is the mixing ratio of water computed from Eice(Tc)
         ! (the sat. vapor pressure of ice at the local Celsius temp.)
         !
         ! C_H2O is given by Dalton's Law as:
         !
         !       C_H2O = Eice( Tc(I,J,L) ) / P(I,J,L)
         !
         ! where P(L) = pressure in grid box (I,J,L)
         !
         ! and   Tc(I,J,L) is the Celsius temp. of grid box (I,J,L).
         !
         ! and   Eice( Tc(I,J,L) ) is the saturation vapor pressure
         !       of ice at temperature Tc(I,J,L) -- computed in
         !       routine E_ICE above.
         !
         ! NOTE: Do not call E_ICE for T > 0C, since C_H2O only needs
         !       to be defined for T <= 0 anyway.
         !==============================================================
         IF ( TC > 0d0 .or. TC < -120d0 ) THEN
            C_H2O(I,k,j) = 0d0
         ELSE
            C_H2O(I,k,j) = E_ICE( TC ) / PL
         ENDIF
       end do
     end do
   end do
   return
   end subroutine initwetdep
!   subroutine solalpha(its,ite,kts,kte,mask,kbot,ktop,j,t,qcld,alphasoluable)
   subroutine solalpha(its,ite,kts,kte,mask,kbot,ktop,j,t,alphasoluable)
   use raqmschem_species_mod
   use raqmschem_const_mod, only : con_ttp
   use raqmschem_pmgrid_mod, only : iam
   implicit none
   integer :: its,ite,kts,kte,i,k,l,nn,n,j
   integer, dimension(its:ite), intent(in) :: kbot,ktop
   real,intent(in) :: t(its:ite,kts:kte)
!   real,intent(in) :: qcld(its:ite,kts:kte)
   real,intent(inout) :: alphasoluable(its:ite,kts:kte,nsol)
!   real :: fice,c_h2o(its:ite,kts:kte),I2G
   real :: fice,I2G
!   real*8 :: cldice(its:ite,kts:kte),cldliq(its:ite,kts:kte)
   real*8 :: r81,r82
   real*8 :: C_TOT,F_L,F_I,L2G
   REAL*8, PARAMETER    :: CONV = 8.27042925126d-1
   logical mask(its:ite)
!   write(6,*)'top sol ',shape(alphasoluable),kind(alphasoluable)
!   call flush(6)
#ifdef OLDWAY
   do i=its,ite
     do k=kbot(i),ktop(i)
       cldliq(i,k)=qcld(i,k)
       if(t(i,k)>=con_ttp)then
         c_h2o(i,k)=qcld(i,k)
         cldice(i,k)=0.
       elseif(t(i,k)<=con_ttp-20.)then
         c_h2o(i,k)=0.0
         cldice(i,k)=qcld(i,k)
       else
         fice=max(0.,min(1.,(con_ttp-t(i,k))*.05))
         cldice(i,k)=fice*qcld(i,k)
         c_h2o(i,k)=qcld(i,k)-cldice(i,k)
       endif
     end do
   end do
#endif
   alphasoluable=0.0
   do n=1,nsol
     nn=idwetd(n)
!     if(iam.eq.0)then
!       write(6,*)'top set alpha n',n,'nn',nn,'p_h2o2',p_h2o2,'p-vrp',p_vrp
!       call flush(6)
!     endif
     if(nn.eq.p_h2o2)then
       do i=its,ite
         if(mask(i))then
           do l=max(2,kbot(i)),ktop(i)

            ! Compute ice to gas ratio for H2O2 by co-condensation
            ! (Eq. 9, Jacob et al, 2000)
            IF ( C_H2O(I,L,j) > 0d0 ) THEN 
               I2G = ( CLDICE(I,L,j) / C_H2O(I,L,j) ) * CONV
            ELSE
               I2G = 0d0
            ENDIF

            ! Compute liquid to gas ratio for H2O2, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            r81=8.3d4
            r82=-7.4d3
            call compute_l2g(r81,r82,T(I,L), CLDLIQ(I,L,j), L2G )

            ! Fraction of H2O2 in liquid & ice phases
            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1d0 + L2G + I2G
            F_L   = L2G / C_TOT
            F_I   = I2G / C_TOT

            ! Compute the rate constant K.  The retention factor for 
            ! liquid H2O2 is 0.05 for 248 K < T < 268 K and 1.0 for 
            ! T >= 268 K. (Eq. 1, Jacob et al, 2000)
            IF ( T(I,L) >= 268d0 ) THEN
               alphasoluable(i,l,n)= F_L +F_I
 
            ELSE IF ( T(I,L) > 248d0  .and. T(I,L) < 268d0 ) THEN
               alphasoluable(i,l,n)=  ( 5d-2 * F_L ) + F_I  

            ELSE
               alphasoluable(i,l,n)=F_I
                  
            ENDIF

           ENDDO
         endif
         ENDDO
            


      !----------------------------
      ! CH3OOH (liquid phase only)
      !----------------------------
!      ELSE IF ( N == IDTMP ) THEN
      elseif(nn.eq.p_ch3ooh)then

         ! No scavenging at the surface
        do i=its,ite
         if(mask(i))then
           do l=max(2,kbot(i)),ktop(i)


            ! Compute liquid to gas ratio for CH3OOH, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            r81=3.1d2
            r82=-5.2d3
            CALL COMPUTE_L2G( r81,r82,T(I,L), CLDLIQ(I,L,j), L2G )

            ! Fraction of CH3OOH in liquid phase
            ! NOTE: CH3OOH does not exist in the ice phase!
            ! (Eq. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1d0 + L2G
            F_L   = L2G / C_TOT

            ! Compute the rate constant K.  The retention factor  
            ! for liquid CH3OOH is 0.0 for T <= 248 K and 0.02 for 
            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
            IF ( T(I,L) >= 268d0 ) THEN
               alphasoluable(i,l,n)=F_L

            ELSE IF ( T(I,L) > 248d0 .and. T(I,L) < 268d0 ) THEN
               alphasoluable(i,l,n)= 2d-2 * F_L  
                  
            ELSE
!               write(200+iam,*)'cold t zero ',i,'k',l,'n',n,'nn',nn,'temp',t(i,l)
!               call flush(200+iam)
               alphasoluable(i,l,n)=0.0

            ENDIF
               

           ENDDO
         endif
         ENDDO



      !----------------------------
!      ! HCl, HF, HBr  (fully soluble, treat like aerosol)
      !----------------------------
!      ELSE IF ( N == idthcl .or. N == idthf .or. N == idthbr ) THEN
      elseif(nn.eq.p_hcl.or.nn.eq.p_hbr.or.nn.eq.p_hno3t)then
        do i=its,ite
         if(mask(i))then
          do l=max(2,kbot(i)),ktop(i)
!           alphasoluable(i,:,n)=1.
           alphasoluable(i,l,n)=1.
         end do
        endif
       end do


      !----------------------------
      ! HNO4 (assuming liquid phase only)
      !----------------------------
!      ELSE IF ( N == idthno4 ) THEN
      elseif(nn.eq.p_hno4)then


       do i=its,ite
         if(mask(i))then
         ! Start scavenging at level 2
         do l=max(2,kbot(i)),ktop(i)


            ! Compute liquid to gas ratio using
            ! the appropriate parameters for Henry's law
            r81=4.0d3
            r82= 0.d0
            CALL COMPUTE_L2G( r81,r82,T(I,L), CLDLIQ(I,L,j), L2G )

            ! Fraction in liquid phase
            ! NOTE: Assuming it does not exist in the ice phase!
            C_TOT = 1d0 + L2G
            F_L   = L2G / C_TOT

            ! Compute the rate constant K using retention factor  
            ! assumed for other species that don't exist in ice phase
            IF ( T(I,L) >= 268d0 ) THEN
               alphasoluable(i,l,n)=F_L
            ELSE IF ( T(I,L) > 248d0 .and. T(I,L) < 268d0 ) THEN
               alphasoluable(i,l,n)= 2d-2 * F_L  
            ELSE
               alphasoluable(i,l,n)= 0d0
            ENDIF
               

         ENDDO
         endif
         ENDDO


      !----------------------------
      ! ethooh, rooh (assuming liquid phase only)
      !----------------------------
!      ELSE IF ( N == idtethooh .OR. N == idtrooh .OR.
!     &          N == idtrip .OR. N == idtprdp .OR.
!     &          N == idtxmrp .OR. N == idtvrp ) THEN
      elseif(nn.eq.p_ethooh.or.nn.eq.p_rooh.or.nn.eq.p_rip.or.nn.eq.p_prdp &
         .or.nn.eq.p_xmrp.or.nn.eq.p_vrp)then

         ! No scavenging at the surface
!       if(iam.eq.0)then
!         write(6,*)'ethooh rooh prdp vrp',nn,'mask',mask,maxval(t),maxval(cldliq)
!         call flush(6)
!       endif
       do i=its,ite
         if(mask(i))then
         do l=max(2,kbot(i)),ktop(i)


            ! Compute liquid to gas ratio using
            ! the appropriate parameters for Henry's law
            r81=3.4d2
            r82=-6.0d3
            CALL COMPUTE_L2G( r81,r82,T(I,L), CLDLIQ(I,L,j), L2G )

            ! Fraction in liquid phase
            ! NOTE: Assuming it does not exist in the ice phase!
            C_TOT = 1d0 + L2G
            F_L   = L2G / C_TOT

            ! Compute the rate constant K using retention factor  
            ! assumed for other species that don't exist in ice phase
            IF ( T(I,L) >= 268d0 ) THEN
               alphasoluable(i,l,n)=F_L
            ELSE IF ( T(I,L) > 248d0 .and. T(I,L) < 268d0 ) THEN
               alphasoluable(i,l,n)= 2d-2 * F_L  
            ELSE
               alphasoluable(i,l,n)= 0d0
            ENDIF
               

         ENDDO
         endif
         ENDDO


      !----------------------------
      ! mgly (assuming liquid phase only)
      !----------------------------
!      ELSE IF ( N == idtmgly ) THEN
      elseif(nn.eq.p_mgly)then

         ! No scavenging at the surface
!       if(iam.eq.0)then
!          write(6,*)'its  mgly',its,ite
!          call flush(6)
!       endif
       do i=its,ite
         if(mask(i))then
!           if(iam.eq.0)then
!              write(6,*)'i',i,j,'kbot',kbot(i),ktop(i)
!              call flush(6)
!           endif
         do l=max(2,kbot(i)),ktop(i)


            ! Compute liquid to gas ratio using
            ! the appropriate parameters for Henry's law
            r81=3.7d3
            r82=-7.5d3
            CALL COMPUTE_L2G( r81,r82,T(I,L), CLDLIQ(I,L,j), L2G )

            ! Fraction in liquid phase
            ! NOTE: Assuming it does not exist in the ice phase!
            C_TOT = 1d0 + L2G
            F_L   = L2G / C_TOT

            ! Compute the rate constant K using retention factor  
            ! assumed for other species that don't exist in ice phase
!           if(iam.eq.0)then
!             write(6,*)'t',l,t(i,l),'cldliq',cldliq(i,l,j)
!             call flush(6)
!           endif
            IF ( T(I,L) >= 268d0 ) THEN
               alphasoluable(i,l,n)=F_L
            ELSE IF ( T(I,L) > 248d0 .and. T(I,L) < 268d0 ) THEN
               alphasoluable(i,l,n)= 2d-2 * F_L  
            ELSE
               alphasoluable(i,l,n)= 0d0
            ENDIF
!           if(iam.eq.0)then
!             write(6,*)'i,l,n',i,l,n,'alpha',alphasoluable(i,l,n)
!             call flush(6)
!           endif
               

         ENDDO
         endif
         ENDDO


      ! CH2O (liquid phase only)
      !----------------------------
!      ELSE IF ( N == IDTCH2O ) THEN 
      elseif(nn.eq.p_ch2o)then

         ! No scavenging at the surface

       do i=its,ite
         if(mask(i))then
         do l=max(2,kbot(i)),ktop(i)


            ! Compute liquid to gas ratio for CH2O, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            r81=3.d3
            r82=-7.2d3
            CALL COMPUTE_L2G( r81,    r82,T(I,L), CLDLIQ(I,L,j), L2G )

            ! Fraction of CH2O in liquid phase 
            ! NOTE: CH2O does not exist in the ice phase!
            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1d0 + L2G
            F_L   = L2G / C_TOT

            ! Compute the rate constant K.  The retention factor 
            ! for liquid CH2O is 0.0 for T <= 248K and 0.02 for 
            ! 248 K < T < 268 K. (Eq. 1, Jacob et al, 2000)
            IF ( T(I,L) >= 268d0 ) THEN
               alphasoluable(i,l,n)=F_L

            ELSE IF ( T(I,L) > 248d0 .and. T(I,L) < 268d0 ) THEN
               alphasoluable(i,l,n)= 2d-2 * F_L  

            ELSE
               alphasoluable(i,l,n)= 0d0

            ENDIF


         ENDDO
         endif
       ENDDO
     endif
   end do
   return
   end subroutine solalpha

      SUBROUTINE COMPUTE_L2G( Kstar298, H298_R, TK, H2OLIQ, L2G )
!
!******************************************************************************
!  Subroutine COMPUTE_L2G computes the ratio L2G = Cliq / Cgas, which is 
!  the mixing ratio of tracer in the liquid phase, divided by the mixing 
!  ratio of tracer in the gas phase.  (bmy, 2/23/00, 11/8/02)
!
!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate 
!  values of Kstar298 and H298_R must be supplied for each tracer.  
!  (cf Jacob et al 2000, p. 3)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) Kstar298 (REAL*8) : Eff. Henry's law constant @ 298 K [moles/atm]
!  (2 ) H298_R   (REAL*8) : Henry's law coefficient           [K]
!  (3 ) TK       (REAL*8) : Temperature at grid box (I,J,L)   [K]
!  (4 ) H2OLIQ   (REAL*8) : Liquid water content at (I,J,L)   [cm3 H2O/cm3 air]
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) L2G      (REAL*8) : Cliq/Cgas ratio for given tracer  [unitless]
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Jacob et al, 2000
!
!  NOTES:
!  (1 ) Bundled into "wetscav_mod.f" (bmy, 11/8/02)
!******************************************************************************
!
      ! Arguments
      use chem_types_mod
      implicit none
      REAL*8, INTENT(IN)  :: KStar298, H298_R
!      REAL*4, INTENT(IN)  :: TK,H2OLIQ
      REAL*4, INTENT(IN)  :: TK
      real(CHEM_KIND_R8),intent(in) :: H2OLIQ
      REAL*8, INTENT(OUT) :: L2G
      
      ! Local variables
      REAL*8              :: Kstar

      ! R = universal gas constant [atm/moles/K]
      REAL*8, PARAMETER   :: R = 8.32d-2

      ! INV_T0 = 1/298 K
      REAL*8, PARAMETER   :: INV_T0 = 1d0 / 298d0

      !=================================================================
      ! COMPUTE_L2G begins here!
      !=================================================================

      ! Get Kstar, the effective Henry's law constant for temperature TK
      Kstar = Kstar298 * EXP( -H298_R * ( ( 1d0 / TK ) - INV_T0 ) )

      ! Use Henry's Law to get the ratio:
      ! [ mixing ratio in liquid phase / mixing ratio in gas phase ]
      L2G   = Kstar * H2OLIQ * R * TK

      ! Return to calling program
      END SUBROUTINE COMPUTE_L2G

      FUNCTION E_ICE( TC ) 
!
!******************************************************************************
!  Subroutine E_ICE computes Eice(T), the saturation vapor pressure of ice
!  at a given Celsius temperature.  (bmy, djj, 2/10/00, 11/8/02)
!  
!  Arguments as Input:
!  ===========================================================================
!  (1 ) TC (REAL*8) : Ambient temperature [C]
!
!  NOTES:
!  (1 ) The saturation vapor pressure of ice is generated by a polynomial
!        formulation from -50 C < T < 0 C.  From -120 C to -50 C, log10
!        interpolations are done in 10 C bins.  This should be good enough
!        for most purposes. (bmy, djj, 2/10/00) 
!  (2 ) Use the "D" exponent (e.g. 1.0d+00) to explicitly force double
!        precision for all floating point numbers (bmy, 2/10/00)
!  (3 ) Eliminate double declaration of E_ICE as REAL*8 (bmy, 7/16/01)
!  (4 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run.  Also updated comments
!        and made cosmetic changes.  Now bundled into "wetscav_mod.f". 
!        (bmy, 11/8/02)
!******************************************************************************
!
      ! References to F90 modules
!jaa      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
     
      ! Arguments 
      REAL*8, INTENT(IN) :: TC

      ! Local variables
      REAL*8             :: E_ICE

      ! Parameters for the polynomial (Lowe & Ficke, 1974)
      REAL*8, PARAMETER  :: A0 = 6.109177956d+00
      REAL*8, PARAMETER  :: A1 = 5.034698970d-01
      REAL*8, PARAMETER  :: A2 = 1.886013408d-02
      REAL*8, PARAMETER  :: A3 = 4.176223716d-04
      REAL*8, PARAMETER  :: A4 = 5.824720280d-06
      REAL*8, PARAMETER  :: A5 = 4.838803174d-08
      REAL*8, PARAMETER  :: A6 = 1.838826904d-10

      !=================================================================
      ! E_ICE begins here!
      !
      ! The polynomial formulation for Eice(T), -50 C <= T <= 0 C, 
      ! is as follows:
      !
      !   Eice(T) = a0 + T(a1 + T(a2 + T(a3 + T(a4 + T(a5 + a6T)))))
      !
      ! with coefficients a0...a6 given above by Lowe & Ficke, 1974.
      !
      ! The log10 interpolations for temperatures below -50 C are 
      ! done in 10 degree C bins.  The equation is as follows:
      !
      !  T0  = Temperature in C at lower bound of temperature bin
      !  T1  = Temperature in C at upper bound of temperature bin
      !  T   = Temperature for which Eice(T) is desired
      !  E0  = Eice(T0) = value of Eice at T0 (obtained from literature)
      !  E1  = Eice(T1) = value of Eice at T1 (obtained from literature)
      !  E   = Eice(T)  = value of Eice at desired temperature T
      !
      !                         {  T - T0                              }
      ! log10(E) =  log10(E0) + { -------- * [ log10(E1) - log10(E0) ] }
      !                         {    10                                }
      !
      !  The relevant values of T0, T1, E0, E1 for each bin are: 
      !
      !  Bin      T0       T1      E0 = Eice(T0)    E1 = Eice(T1)     
      ! ----------------------------------------------------------------
      !   1     -60 C    -50 C     1.0765e-2 mb     3.947e-2  mb
      !   2     -70 C    -60 C     2.577e-3  mb     1.0765e-2 mb
      !   3     -80 C    -70 C     5.307e-4  mb     2.577e-3  mb
      !   4     -90 C    -80 C     9.17e-5   mb     5.307e-4  mb
      !   5    -100 C    -90 C     1.29e-5   mb     9.17e-5   mb
      !   6    -110 C   -100 C     1.4e-6    mb     1.29e-5   mb
      !   7    -120 C   -110 C     1.0e-7    mb     1.4e-6    mb
      !
      ! For computational expediency, the following terms were also 
      ! computed for use in the code below:
      !  
      !  Bin           log10(E0)         [ log10(E1) - log10(E0) ]  
      ! ----------------------------------------------------------------
      !   1        -1.96798596584            5.6425309224e-01            
      !   2        -2.58888558145            6.2089961561e-01            
      !   3        -3.27515091237            6.8626533092e-01 
      !   4        -4.03763066433            7.6247975196e-01  
      !   5        -4.88941028970            8.5177962537e-01 
      !   6        -5.85387196432            9.6446167462e-01
      !   7        -7.00000000000            1.14612803568e+00
      !=================================================================

      ! Use a polynomial formulation for the saturation vapor 
      ! pressure over ice from -50 C to 0 C (Lowe & Ficke, 1974)
      IF ( TC >= -50d0 .and. TC <= 0d0 ) THEN
         E_ICE = A0 + TC * ( A1 + TC * ( A2 + TC * ( A3 +  &
                     TC * ( A4 + TC * ( A5 + A6 * TC ) ) ) ) )

      ! Bin 1 -- Log10 interpolation 
      ELSE IF ( TC >= -60d0 .and. TC < -50d0 ) THEN
         E_ICE = -1.96798596584d0 +  &
                ( ( ( TC + 60d0 ) / 10d0 ) * 5.6425309224d-1 )

         E_ICE = 10d0**E_ICE

      ! Bin 2 -- Log10 interpolation
      ELSE IF ( TC >= -70d0 .and. TC < -60d0 ) THEN
         E_ICE = -2.58888558145d0 + &
                 ( ( ( TC + 70d0 ) / 10d0 ) * 6.2089961561d-1 ) 
                 
         E_ICE = 10d0**E_ICE

      ! Bin 3 -- Log10 interpolation
      ELSE IF ( TC >= -80d0 .and. TC < -70d0 ) THEN
         E_ICE = -3.27515091237d0 + &
                ( ( ( TC + 80d0 ) / 10d0 ) * 6.8626533092d-1 ) 

         E_ICE = 10d0**E_ICE

      ! Bin 4 -- Log10 interpolation
      ELSE IF ( TC >= -90d0 .and. TC < -80d0 ) THEN
         E_ICE = -4.03763066433d0 + &
                ( ( ( TC + 90d0 ) / 10d0 ) * 7.6247975196d-1 ) 

         E_ICE = 10d0**E_ICE

      ! Bin 5 -- Log10 interpolation
      ELSE IF ( TC >= -100d0 .and. TC < -90d0 ) THEN 
         E_ICE =  -4.88941028970d0  + &
                ( ( ( TC + 100d0 ) / 10d0 ) * 8.5177962537d-1 ) 

         E_ICE = 10d0**E_ICE

      ! Bin 6 -- Log10 interpolation
      ELSE IF ( TC >= -110d0 .and. TC < -100d0 ) THEN
         E_ICE =  -5.85387196432d0 +  &
                 ( ( ( TC + 110d0 ) / 10d0 ) * 9.6446167462d-1 ) 

         E_ICE = 10d0**E_ICE

      ! Bin 7 -- Log10 interpolation
      ELSE IF ( TC >= -120d0 .and. TC < -110d0 ) THEN
         E_ICE = -7d0 + ( ( ( TC + 120d0 ) / 10d0 ) * 1.14612803568d0 ) 

         E_ICE = 10d0**E_ICE

      ! Below -120 C stop with an error message
      ELSE
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '( ''Temperature: '', f10.5, '' C'' )' ) TC
         WRITE( 6, '( ''T must be in the range -120 C <= T <= 0C'' )' )
         WRITE( 6, '( ''STOP in e_ice.f'' )' )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!jaa         CALL GEOS_CHEM_STOP
         stop 'E_ICE in wetscav code'
!jaa
      ENDIF

      ! Return to calling program
      END FUNCTION E_ICE
end module wetdep_alpha_mod
