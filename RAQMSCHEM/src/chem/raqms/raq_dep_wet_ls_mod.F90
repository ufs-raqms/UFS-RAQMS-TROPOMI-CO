module raq_dep_wet_ls_mod
#include <options.h>


  implicit none

  private

  public :: wetdep_ls,WetRemovalGOCART

contains

subroutine wetdep_ls(dt,var,rain,moist,t,rho,var_rmv,num_moist, &
         num_chem,p_qc,p_qi,dz8w,vvel,             &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  use chem_types_mod
  use raqmschem_species_mod, only : chemsoluablefull,ichemfullpt,nsol, &
     chempointsoluablefull,lraqmschem,chemname,chemsoluable,cheminputlist

  use wetdep_alpha_mod, only : solalpha
  use wetdep_fam_mod, only :wetdeplsfam
  use raqmschem_pmgrid_mod, only : iam,iprnin,jprnin,tile,iamprn,kprnin
  use raqmschem_species_mod, only : p_co
#ifdef WETDEPDIAG
  use raqmschemcomm_mod, only : lpaveij
#endif
  IMPLICIT NONE
!  ids,ide,jds,jde not used by ls
   INTEGER,      INTENT(IN   ) :: num_chem,num_moist,p_qc,p_qi,            &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   real(CHEM_KIND_R4), INTENT(IN ) :: dt
    REAL(CHEM_KIND_R4), DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL(CHEM_KIND_R4),  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) :: rho,dz8w,vvel,t        
!   our var uses nl not ni like rest
   REAL(CHEM_KIND_R4),  DIMENSION( ims:ime , jms:jme , kts:kte ,1:num_chem),                        &
          INTENT(INOUT) :: var        
!   REAL,  DIMENSION( jms:jme ),                                  &
   REAL(CHEM_KIND_R4),  DIMENSION( ims:ime,jms:jme ),                                  &
          INTENT(IN   ) :: rain
   REAL(CHEM_KIND_R4),  DIMENSION( ims:ime ,  jms:jme,num_chem ),                                  &
          INTENT(INOUT   ) :: var_rmv
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: var_sum
   REAL(CHEM_KIND_R4),  DIMENSION( its:ite ,  kts:kte, jts:jte ) :: var_rmvl
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: frc,var_sum_clw,rain_clw     
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: var_sum_clw_wrf
!    real :: dvar,factor,clsum,alpha,rho_water
    real :: dvar,factor,clsum,rho_water
   real alpha(kts:kte),alphasoluable(its:ite,kts:kte,nsol),newfrc
   real t2d(its:ite,kts:kte)
!   real qcld(its:ite,kts:kte)
   integer :: nv,i,j,k,km,kb,kbeg,nmask,kbot(its:ite),ktop(its:ite)
   logical maskalpha(its:ite)
   character *10 cenv
   logical first,dowetdeplswrf
   save first
   data first/.true./
!  in ls all three indexes are the same and full limits dont have to start at 1
!   write(400+iam,*)'ls its',its,ite,' jts ',jts,jte
!   write(400+iam,*)'ls ims',ims,ime,' jms ',jms,jme
!   write(6,*)'kind ls ',kind(qcld) real*4
!   write(6,*)'alphasol ',shape(alphasoluable),'kind',kind(alphasoluable)
!   call flush(6)
!   write(6,*)'top wetdep'
!   call flush(6)
    alphasoluable=.5 ! ajl for now
    if(first)then
      cenv=' '
      call getenv('WETDEPLSWRF',cenv)
      if(cenv.eq.'YES')then
        dowetdeplswrf=.true.
        if(iam.eq.0)then
          write(6,*)'raqms dowetdeplswrfchem'
        endif
      else
        dowetdeplswrf=.false.
      endif
    endif
    if(dowetdeplswrf)then
      do j=jts,jte
        do i=its,ite
          var_sum_clw(i,j)=0.
          var_sum_clw_wrf(i,j)=0.
          var_sum(i,j)=0.
          var_rmvl(i,:,j)=0.
          rain_clw(i,j)=0.
          frc(i,j)=0.0
          if(rain(i,j).gt.1.e-3)then
! convert rain back to rate
!
             rain_clw(i,j)=rain(i,j)/dt
             do k=1,kte-1
               dvar=max(0.,(moist(i,k,j,p_qc)+moist(i,k,j,p_qi))*rho(i,k,j)*vvel(i,k,j)*dz8w(i,k,j))
               var_sum_clw_wrf(i,j)=var_sum_clw_wrf(i,j)+dvar
             end do
             if(var_sum_clw_wrf(i,j)>1.e-10)then
               frc(i,j)=rain_clw(i,j)/var_sum_clw_wrf(i,j)
!               write(6,*)'raqms frc ', frc(i,j),var_sum_clw_wrf(i,j),'rain_clw',rain_clw(i,j)
               frc(i,j)=max(1.e-6,min(frc(i,j),.005))
             endif
! total cloud water
!
             do k=1,kte
                dvar=max(0.,(moist(i,k,j,p_qc)+moist(i,k,j,p_qi)))
                var_sum_clw(i,j)=var_sum_clw(i,j)+dvar
             enddo
          endif
        enddo
      enddo
    else
      cenv=' '
      call getenv('WETDEPLSFRC',cenv)
      if(cenv.eq.' ')then
        frc(:,:)=.1 ! new for 8.8
      else
        read(cenv,*)newfrc
        frc(:,:)=newfrc
        if(first)then
           write(6,*)'raqms newwetdepls frc',newfrc
           first=.false.
        endif
      endif
      do j=jts,jte
        do i=its,ite
          var_sum_clw(i,j)=0.
          var_sum(i,j)=0.
          var_rmvl(i,:,j)=0.
!          frc(i,j)=.1 ! new for 8.8
          rain_clw(i,j)=0.
          if(rain(i,j).gt.1.e-3)then
! convert rain back to rate
!
!             rain_clw(i,j)=rain(j)/dt
             rain_clw(i,j)=rain(i,j)/dt
! total cloud water
!
!             do k=1,kte-1
             do k=1,kte
                dvar=max(0.,(moist(i,k,j,p_qc)+moist(i,k,j,p_qi)))
                var_sum_clw(i,j)=var_sum_clw(i,j)+dvar
             enddo
          endif
        enddo
      enddo
    endif
    rho_water = 1000.
    var_rmv (:,:,:)=0.
#define NODIAGWETDEP
    do nv=1,num_chem
!      if(lraqmschem(nv))then
!     only do for chemicals that are soluable
#ifdef DIAGWETDEP
      if(iam.eq.iamprn.and.chemsoluablefull(nv))then
        write(6,*)' ls nv ',nv,' cehmsoluablefull ',chemsoluablefull(nv)
        call flush(6)
      endif
#endif
      if(.not.chemsoluablefull(nv))cycle
        
!
! get rid of it
!
      kbot=kts
      ktop=kte
!     write(6,*)'var_sum',maxval(var_sum),'var_sum_clw',maxval(var_sum_clw),'qc', &
!     maxval(moist(:,:,:,p_qc))
!     call flush(6)
      do j=jts,jte
        maskalpha=0
        nmask=0
        do i=its,ite
          if(rain(i,j).gt.1.e-6)then
           maskalpha(i)=.true.
           nmask=nmask+1
           endif
        end do
        if(nmask>0)then
!          if(iam.eq.iamprn.and.j.eq.jprnin)then
!          write(6,*)'nmask large ',nmask
!          call flush(6)
!          endif
          do k=kts,kte
            do i=its,ite
!              qcld(i,k)=moist(i,k,j,p_qc)
              t2d(i,k)=t(i,k,j)
            end do
          end do
!          write(6,*)'wetls ',shape(alphasoluable),kind(alphasoluable)
!          call flush(6)
!          call solalpha(its,ite,kts,kte,maskalpha,kbot,ktop,j,t2d,qcld,alphasoluable)
          call solalpha(its,ite,kts,kte,maskalpha,kbot,ktop,j,t2d,alphasoluable)
        endif
        do i=its,ite
#ifdef DIAGWETDEP
         if(iam.eq.iamprn.and.i.eq.iprnin.and.j.eq.jprnin)then
           if(rain(i,j)>1.e-3)then
           write(6,*)'rain',rain(i,j),'var_sum  nv ',nv,var_sum(i,j),'var_sum_clw ',var_sum_clw(i,j)
           call flush(6)
           endif
         endif
#endif
!          if(iam.eq.iamprn.and.i.eq.iprnin.and.j.eq.jprnin)then
!             write(6,*)'rain 2 ',rain(i,j),'var_sum',var_sum(i,j),var_sum_clw(i,j),' nv ',nv,' p_co',p_co
!          endif
         if(rain(i,j).gt.1.e-6)then
!           if(iam.eq.iamprn.and.j.eq.jprnin)then
!             write(6,*)'rain at ',i,j,rain(i,j)
!           endif
           if(var_sum_clw(i,j).gt.1.e-10)then
             do k=kts,kte
               alpha(k)=alphasoluable(i,k,chempointsoluablefull(nv))
             end do
             do k=kts,kte-2
!               if(i.eq.iprnin.and.j.eq.jprnin.and.nv.eq.p_co)then
!                 write(6,*)'co ls bef',nv,k,var(i,j,k,nv),'rain',rain(i,j)
!               endif
               if(var(i,j,k,nv).gt.1.e-19 .and.  (moist(i,k,j,p_qc)+moist(i,k,j,p_qi)).gt.1.e-8)then
                 factor = max(0.,frc(i,j)*rho(i,k,j)*dz8w(i,k,j)*vvel(i,k,j))
                 dvar=max(0.,alpha(k)*factor/(1+factor)*var(i,j,k,nv))
                 dvar=min(dvar,var(i,j,k,nv))
                 if((var(i,j,k,nv)-dvar).lt.1.e-24)then
                   dvar=var(i,j,k,nv)-1.e-24
                   var(i,j,k,nv)=var(i,j,k,nv)-dvar
                 else
                   var(i,j,k,nv)=var(i,j,k,nv)-dvar
                 endif
                 var_rmvl(i,k,j)=dvar
!               if(i.eq.iprnin.and.j.eq.jprnin.and.nv.eq.p_co)then
!                 write(6,*)'co ls aft',nv,k,var(i,j,k,nv),'rain',rain(i,j),'dvar',dvar
!               endif
#ifdef DIAGWETDEP
         if(iam.eq.iamprn.and.i.eq.iprnin.and.j.eq.jprnin.and.k.eq.kprnin)then 
           write(6,*)'var_rmvl ',nv,dvar
           call flush(6)
         endif
#endif
            
              
                 var_rmv(i,j,nv)=var_rmv(i,j,nv)+var_rmvl(i,k,j)
               endif
             enddo
!           var_rmv(i,j)=var_rmv(i,j)+var_rmvl(i,j)
           endif
          endif
        enddo
      enddo
!     handle wetdep ls  family update and diagnostics
!      write(400+iam,*)'call wetdeplsfam',its,ite,jts,jte
!      call flush(400+iam)
!     its etc are full indices dont have to start at 1
!     if(iam.eq.iamprn)then
!       write(6,*)'var_rmv',nv,var_rmvl(iprnin,1,jprnin)
!     endif
!      write(300+iam,*)'call wetdeplsfam nv ',nv,num_chem
!      call flush(300+iam)
#ifdef DIAGWETDEP
      if(iam.eq.iamprn.and.var_rmvl(iprnin,kprnin,jprnin).ne.0.0)then
        write(6,*)'call wetdepslam nv',nv,'num_chem',num_chem,'var_rmvl ',chemsoluablefull(nv),var_rmvl(iprnin,kprnin,jprnin)
        call flush(6)
      endif
#endif
!      write(6,*)'num_chem',num_chem,'nv',nv,'shape var',shape(var)
!      call flush(6)
      call wetdeplsfam(nv,num_chem,its,ite,jts,jte,kts,kte,var,var_rmvl)
    enddo
#ifdef WETDEPDIAG
    do j=jts,jte
      do i=its,ite
        lpaveij(i,j)=lpaveij(i,j)+rain_clw(i,j)
      end do
    end do
!      write(300+iam,*)'bottom rain_clw',maxval(rain_clw),'dt',dt
!      call flush(300+iam)
#endif
#ifdef DIAGWETDEP
    if(iam.eq.iamprn)then
    write(6,*)'num_chem ls',num_chem,'co out 1',var(iprnin,jprnin,1,p_co)
    write(6,*)'num_chem ls',num_chem,'co out 64',var(iprnin,jprnin,1,p_co)
    call flush(6)
    endif
#endif
    first=.false.
END SUBROUTINE WETDEP_LS


!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: WetRemovalGOCART - Calculate aerosol wet removal due
!                               to large scale processes.
!
! !INTERFACE:
!

  subroutine WetRemovalGOCART ( i1, i2, j1, j2, k1, k2, n1, n2, cdt, &
                                num_chem, var_rmv, chem, ple, tmpu,  &
                                rhoa, dqcond, precl,         &
                                ims, ime, jms, jme, kms, kme, rc,precc )

! !USES:
   use chem_rc_mod

   use chem_types_mod
   use chem_const_mod, only : grav => grvity
   use raqmschem_species_mod, only : chemsoluablefull,ichemfullpt,nsol, &
     chempointsoluablefull,lraqmschem,chemname,chemsoluable,cheminputlist

   use wetdep_alpha_mod, only : solalpha
   use raqmschem_pmgrid_mod, only : iam,tile
   use wetdep_fam_mod, only :wetdeplsfam
#ifdef WETDEPDIAG
   use raqmschemcomm_mod, only : lpaveij
#endif
   IMPLICIT NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, k1, k2, n1, n2, num_chem, &
                          ims, ime, jms, jme, kms, kme
   real, intent(in)    :: cdt
   REAL(CHEM_KIND_R4),  DIMENSION( i1:i2 , j1:j2, k1:k2 , 1:num_chem),&
          INTENT(INOUT) :: chem
   REAL(CHEM_KIND_R4),  DIMENSION( ims:ime ,  jms:jme,num_chem ), &
          INTENT(INOUT   ) :: var_rmv !! tracer loss flux [kg m-2 s-1]
   real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme),&
          INTENT(IN)     :: ple, tmpu, rhoa, dqcond
   real(CHEM_KIND_R4), dimension(ims:ime ,  jms:jme) , &
          INTENT(IN)      ::  precl    ! cv, ls precip [mm day-1]
   real(CHEM_KIND_R4), dimension(ims:ime ,  jms:jme) , &
          INTENT(IN),optional    :: precc    ! cv, ls precip [mm day-1]
   real(CHEM_KIND_R4) :: alpha(k1:k2),alphasoluable(i1:i2,k1:k2,nsol)

! !OUTPUT PARAMETERS:
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 -

! !DESCRIPTION: Calculates the updated species concentration due to wet
!               removal.  As written, intended to function for large
!               scale (not convective) wet removal processes

!
! !REVISION HISTORY:
!
!  08Jan2010 - Colarco, based on GOCART implementation, does not
!                       include any size dependent term
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   character(len=*), parameter :: myname = 'WetRemovalGOCART'
   integer  ::  i, j, k, n, nbins, LH, kk, ios,nv
   real(CHEM_KIND_R4) :: pdog(i1:i2,k1:k2,j1:j2)      ! air mass factor dp/g [kg m-2]
   real :: pls, pcv, pac             ! ls, cv, tot precip [mm day-1]
   real :: qls(k1:k2), qcv(k1:k2)          ! ls, cv portion dqcond [kg m-3 s-1]
   real :: qmx, qd, A                ! temporary variables on moisture
   real :: F, B, BT                  ! temporary variables on cloud, freq.
   real, allocatable :: fd(:,:)      ! flux across layers [kg m-2]
   real, allocatable :: DC(:)        ! scavenge change in mass mixing ratio
!  Rain parameters from Liu et al.
   real, parameter :: B0_ls = 1.0e-4
   real, parameter :: F0_ls = 1.0
   real, parameter :: XL_ls = 5.0e-4
   real, parameter :: B0_cv = 1.5e-3
   real, parameter :: F0_cv = 0.3
   real, parameter :: XL_cv = 2.0e-3
!  Duration of rain: ls = model timestep, cv = 1800 s (<= cdt)
   real            :: Td_ls
   real            :: Td_cv
   logical bottop
   logical maskalpha(i1:i2)
   integer nmask,kbota(i1:i2),ktopa(i1:i2),kbot,ktop,kinc
   real t2d(i1:i2,k1:k2)
   REAL,  DIMENSION( i1:i2 ,  k1:k2, j1:j2 ) :: var_rmvl
   real dval


!  Efficiency of dust wet removal (since dust is really not too hygroscopic)
!  Applied only to in-cloud scavenging
   real :: effRemoval
!  real,dimension(20) ::fwet
!  tracer: p_so2=1 p_sulf=2 p_dms=3 p_msa=4 p_p25=5 p_bc1=6 p_bc2=7 p_oc1=8
!  p_oc2=9 p_dust_1=10 p_dust_2=11 p_dust_3=12 p_dust_4=13 p_dust_5=14
!  p_seas_1=15 p_seas_2=16 p_seas_3=17 p_seas_4=18 p_seas_5=19 p_p10  =20
!   data fwet /0.,1.5,0.,0.,1.,0.7,0.7,0.4,0.4,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./
!  rc=0.

!  Initialize local variables
!  --------------------------
   rc = CHEM_RC_SUCCESS

   Td_ls = cdt
   Td_cv = cdt
   nbins = n2-n1+1
   var_rmv = 0.0

!  Allocate the dynamic arrays
   allocate(fd(k1:k2,nbins),stat=ios)
   if (chem_rc_test((ios .ne. 0), msg="Failed to allocate memory", &
     file=__FILE__, line=__LINE__, rc=rc)) return
   allocate(dc(nbins),stat=ios)
   if (chem_rc_test((ios .ne. 0), msg="Failed to allocate memory", &
     file=__FILE__, line=__LINE__, rc=rc)) return

!  Accumulate the 3-dimensional arrays of rhoa and pdog
   if(ple(i1,k1,j1)>ple(i1,k2,j2))then
!    k1 is bottom
     kbot=k1
     ktop=k2
     kinc=-1
     bottop=.true.
   else 
     kbot=k2
     ktop=k1
     kinc=1
     bottop=.false.
   end if
   do j = j1, j2
    do i = i1, i2
      if(bottop)then
        do k=k1,k2
         pdog(i,k,j)=(ple(i,k,j)-ple(i,k+1,j))/grav
!         if(i.eq.57.and.j.eq.15)then
!           write(6,*)'pdog',k,pdog(i,k,j)
!         endif
        end do
      else 
        do k=k1,k2
         pdog(i,k,j)=(ple(i,k+1,j)-ple(i,k,j))/grav
        end do
      endif
!      pdog(i,k1:k2,j) = (ple(i,k1+1:k2+1,j)-ple(i,k1:k2,j)) / grav 
    enddo
   enddo
#ifdef WETDEPDIAG
    do j=j1,j2
      do i=i1,i2
        lpaveij(i,j)=lpaveij(i,j)+precl(i,j)
      end do
    end do
!      write(300+iam,*)'bottom rain_clw',maxval(rain_clw),'dt',dt
!      call flush(300+iam)
#endif

   do nv=1, num_chem
!  Loop over spatial indices
    if(.not.chemsoluablefull(nv))cycle
   var_rmvl=0.
   do j = j1, j2
     maskalpha=0
     nmask=0
     if(present(precc))then
       do i=i1,i2
         if(precl(i,j)+precc(i,j).gt.0.0)then
           maskalpha(i)=.true.
           nmask=nmask+1
         endif
       end do
     else
       do i=i1,i2
         if(precl(i,j).gt.0.0)then
           maskalpha(i)=.true.
           nmask=nmask+1
         endif
       end do
     endif
     if(nmask==0)cycle
     do k=k1,k2
       do i=i1,i2
         t2d(i,k)=tmpu(i,k,j)
       end do
     end do
     kbota=k1
     ktopa=k2
     call solalpha(i1,i2,k1,k2,maskalpha,kbota,ktopa,j,t2d,alphasoluable)
    do i = i1, i2

!    Check for total precipitation amount
!    Assume no precip in column if precl+precc = 0
     if(present(precc))then
       pac = precl(i,j) + precc(i,j)
       pcv = precc(i,j)
     else
       pac=precl(i,j)
       pcv=0.0
     endif
     if(pac .le. 0.) cycle
     pls = precl(i,j)

!    Initialize the precipitation fields
     qls(:)  = 0.
     qcv(:)  = 0.
     fd(:,:) = 0.

!    Find the highest model layer experiencing rainout.  Assumes no
!    scavenging if T < 258 K
     LH = 0
     do k = ktop,kbot,kinc
      if(dqcond(i,k,j) .lt. 0. .and. tmpu(i,k,j) .gt. 258.) then
       LH = k
       goto 15
      endif
     end do
 15  continue
!     !if(tile.eq.3.and.i.eq.94.and.j.eq.65)then
!        write(6,*)'lh',lh
!     endif
     if(LH .lt. 1) cycle
!     if(tile.eq.3)then
!       write(6,*)'i',i,j,pls,'lh',lh
     !endif

!    convert dqcond from kg water/kg air/s to kg water/m3/s and reverse
!    sign so that dqcond < 0. (positive precip) means qls and qcv > 0.
     do k = LH, kbot,kinc
      qls(k) = -dqcond(i,k,j)*pls/pac*rhoa(i,k,j)
      qcv(k) = -dqcond(i,k,j)*pcv/pac*rhoa(i,k,j)
     end do
     do k=LH,kbot,kinc
       alpha(k)=alphasoluable(i,k,chempointsoluablefull(nv))
     end do

!    Loop over vertical to do the scavenging!
!     do k = LH, k2
     do k = LH, kbot,kinc

!-----------------------------------------------------------------------------
!   (1) LARGE-SCALE RAINOUT:
!       Tracer loss by rainout = TC0 * F * exp(-B*dt)
!         where B = precipitation frequency,
!               F = fraction of grid box covered by precipitating clouds.
!       We assume that tracer scavenged by rain is falling down to the
!       next level, where a fraction could be re-evaporated to gas phase
!       if Qls is less then 0 in that level.
!-----------------------------------------------------------------------------
!     if(tile.eq.3.and.i.eq.94.and.j.eq.65)then
!         write(6,*)'qls',k,qls(k)
!     endif
      if (qls(k) .gt. 0.) then
!        if(iam.eq.41.and.nv.eq.47)then
!          write(200+iam,*)'qls',k,i,j,'nv',nv,qls(k),chem(i,j,k,nv)
!          call flush(200+iam)
!        endif
       F  = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qls(k)*cdt/Td_ls))
       B  = B0_ls/F0_ls +1./(F0_ls*XL_ls/qls(k))
       BT = B * Td_ls
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >
!      Adjust du level:
       do n = 1, nbins
         effRemoval = alpha(k)
         DC(n) = chem(i,j,k,nv) * F * effRemoval *(1.-exp(-BT))
         if (DC(n).lt.0.) DC(n) = 0.
         dval=DC(n)
         if (chem(i,j,k,nv)-DC(n) .lt. 1.0E-32) then
!          fix subscript below had i,k,j 2/25/2021
           dval=chem(i,j,k,nv)-1.0e-32
           chem(i,j,k,nv) = 1.0E-32
         else
           chem(i,j,k,nv) = chem(i,j,k,nv)-DC(n)
         endif
         var_rmvl(i,k,j)=var_rmvl(i,k,j)+dval
       end do
!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = DC(n)*pdog(i,k,j)
       end do

      end if                                    ! if Qls > 0  >>>

!-----------------------------------------------------------------------------
! * (2) LARGE-SCALE WASHOUT:
! *     Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------
      if(k /= LH .and. qls(k) .ge. 0.) then
       if(qls(k) .lt. qls(k-kinc)) then
!       Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-kinc,LH,-kinc
         if (Qls(kk).gt.0.) then
          Qmx = max(Qmx,Qls(kk))
         else
          goto 333
         end if
        end do

 333    continue
        F = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(Qmx*cdt/Td_ls))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx /rhoa(i,k,j)*pdog(i,k,j)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Adjust du level:
        do n = 1, nbins
         DC(n) = chem(i,j,k,nv) * F * (1.-exp(-BT))
         if (DC(n).lt.0.) DC(n) = 0.
         chem(i,j,k,nv) = chem(i,j,k,nv)-DC(n)
         dval=DC(n)
         if (chem(i,j,k,nv) .lt. 1.0E-32) then 
           dval=dval-(1.0e-32-chem(i,j,k,nv))
           chem(i,j,k,nv) = 1.0E-32
         endif
          var_rmv(i,j,nv) = var_rmv(i,j,nv)+DC(n)*pdog(i,k,j)/cdt
          var_rmvl(i,k,j)=var_rmvl(i,k,j)+dval
        end do

       end if
      end if                                    ! if ls washout  >>>
!-----------------------------------------------------------------------------
!  (3) CONVECTIVE RAINOUT:
!      Tracer loss by rainout = dd0 * F * exp(-B*dt)
!        where B = precipitation frequency,
!              F = fraction of grid box covered by precipitating clouds.
!-----------------------------------------------------------------------------
      if(present(precc))then
       if (qcv(k) .gt. 0.) then
        F  = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qcv(k)*cdt/Td_cv))
        B  = B0_cv
        BT = B * Td_cv
        if (BT.gt.10.) BT = 10.               !< Avoid overflow >

!      Adjust du level:
       do n = 1, nbins
        effRemoval = alpha(k)
        DC(n) = chem(i,j,k,nv) * F * effRemoval * (1.-exp(-BT))
        if (DC(n).lt.0.) DC(n) = 0.
        dval=dc(n)
        chem(i,j,k,nv) = chem(i,j,k,nv)-DC(n)
        if (chem(i,j,k,nv)<1.0e-32)then
          dval=dval-(1.0e-32-chem(i,j,k,nv))
          chem(i,j,k,nv) = 1.0E-32
        endif
         var_rmvl(i,k,j)=var_rmvl(i,k,j)+dval
       end do

!------  Flux down:  unit is kg. Including both ls and cv.
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + DC(n)*pdog(i,k,j)
       end do

      end if                                  ! if Qcv > 0   >>>

!-----------------------------------------------------------------------------
!  (4) CONVECTIVE WASHOUT:
!      Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------

      if (k/=LH .and. Qcv(k).ge.0.) then
       if (Qcv(k).lt.Qcv(k-kinc)) then
!-----  Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-kinc, LH, -kinc
         if (Qcv(kk).gt.0.) then
          Qmx = max(Qmx,Qcv(kk))
         else
          goto 444
         end if
        end do

 444    continue
        F = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qmx*cdt/Td_cv))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx / rhoa(i,k,j)*pdog(i,k,j)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Adjust du level:
        do n = 1, nbins
         DC(n) = chem(i,j,k,nv) * F * (1.-exp(-BT))
         if (DC(n).lt.0.) DC(n) = 0.
         chem(i,j,k,nv) = chem(i,j,k,nv)-DC(n)
         dval=dc(n)
         if (chem(i,j,k,nv) .lt. 1.0E-32)then
           dval=dval-(1.e-32-chem(i,j,k,nv))
           chem(i,j,k,nv) = 1.0E-32
         endif
          var_rmv(i,j,nv) = var_rmv(i,j,nv)+DC(n)*pdog(i,k,j)/cdt
          var_rmvl(i,k,j)=var_rmvl(i,k,j)+dval
        end do

       end if
      end if                                    ! if cv washout  >>>
      end if ! present precc
!-----------------------------------------------------------------------------
!  (5) RE-EVAPORATION.  Assume that SO2 is re-evaporated as SO4 since it
!      has been oxidized by H2O2 at the level above.
!-----------------------------------------------------------------------------
!     Add in the flux from above, which will be subtracted if reevaporation occurs
      if(k /= LH) then
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + Fd(k-kinc,n)
       end do

!      Is there evaporation in the currect layer?
       if (-dqcond(i,k,j) .lt. 0.) then
!       Fraction evaporated = H2O(k)evap / H2O(next condensation level).
        if (-dqcond(i,k-kinc,j) .gt. 0.) then

          A =  abs(  dqcond(i,k,j) * pdog(i,k,j)    &
            /      ( dqcond(i,k-kinc,j) * pdog(i,k-kinc,j))  )
          if (A .gt. 1.) A = 1.

!         Adjust tracer in the level
          do n = 1, nbins
           DC(n) =  Fd(k-kinc,n) / pdog(i,k,j) * A
           chem(i,j,k,nv) = chem(i,j,k,nv) + DC(n)
           dval=-dc(n)
           if(chem(i,j,k,nv)<1.e-32)then
             dval=dval+1.e-32-chem(i,j,k,nv)
           endif
           chem(i,j,k,nv) = max(chem(i,j,k,nv),1.e-32)
!          Adjust the flux out of the bottom of the layer
           Fd(k,n)  = Fd(k,n) - DC(n)*pdog(i,k,j)
           var_rmvl(i,k,j)=var_rmvl(i,k,j)+dval

          end do

        endif
       endif                                   ! if -moistq < 0
      endif
     end do  ! k

     do n = 1, nbins
       var_rmv(i,j,nv) = var_rmv(i,j,nv)+Fd(kbot,n)/cdt
     end do
!     if(iam.eq.41.and.nv.eq.47)then
!       write(200+iam,*)i,j,'var_rmv',var_rmv(i,j,nv)
!     endif
!     if(tile.eq.3.and.i.eq.80.and.j.eq.93)then
!       write(6,*)'k1',k1,k2,kinc,'kbot',kbot,ktop,'nv',nv
!       do k=ktop,kbot,kinc
!         if(var_rmvl(i,k,j)/=0.0)then
!           write(6,*)'var_rmvl',nv,'k',k,var_rmvl(i,k,j)
!         endif
!       end do
     !endif

    end do   ! i
   end do    ! j
!   if(tile.eq.3)then
!       write(200+iam,*)'var_rmvl',nv,maxval(var_rmv(:,:,nv))
!   endif
   call wetdeplsfam(nv,num_chem,i1,i2,j1,j2,k1,k2,chem,var_rmvl)
   end do    !nv for num_chem

   deallocate(fd,DC,stat=ios)
   if (chem_rc_test((ios .ne. 0), msg="Failed to deallocate memory", &
     file=__FILE__, line=__LINE__, rc=rc)) return
!   write(6,*)'bot wetdeplsnew'
!   call flush(6)
!   write(200+iam,*)'bot setdep harv'
!   call flush(200+iam)

   end subroutine WetRemovalGOCART

end module raq_dep_wet_ls_mod
#ifdef WETDEPDIAG
subroutine savewetdep(i1,i2,j1,j2,k1,k2,ivar,num_chem,var_rm,var_rml)
use chem_types_mod
use raqmschem_pmgrid_mod,only :iam
use raqmschemcomm_mod,only : gsdcol,gsd3d
integer i1,i2,j1,j2,k1,k2,ivar,num_chem,ii
real(CHEM_KIND_R4),dimension(i1:i2,j1:j2,num_chem) :: var_rm
real(CHEM_KIND_R4),dimension(i1:i2,k1:k2,j1:j2)          :: var_rml
select case(ivar)
case (6)
  mvar=1
case (7)
  mvar=2
case (8)
  mvar=3
case (9)
  mvar=4
case default
  mvar=0
end select
if (mvar/=0)then
 do j=j1,j2
  do k=k1,k2
    do i=i1,i2
      ii=i-i1+1
      gsd3d(ii,j,k,mvar)=gsd3d(ii,j,k,mvar)+var_rml(i,k,j)
    enddo
  enddo
  do i=i1,i2
    ii=i-i1+1
    gsdcol(ii,j,mvar)=gsdcol(ii,j,mvar)+var_rm(i,j,ivar)
  end do
 enddo
! write(6,*)'mvar',mvar,maxval(gsdcol(:,:,mvar)),minval(gsdcol(:,:,mvar))
endif
return
end subroutine savewetdep
subroutine savewetdepr(i1,i2,j1,j2,k1,k2,ivar,num_chem,var_rm,var_rml_r)
use chem_types_mod
use raqmschem_pmgrid_mod,only :iam
use raqmschemcomm_mod,only : gsdcol,gsd3d
integer i1,i2,j1,j2,k1,k2,ivar,num_chem,ii
real(CHEM_KIND_R4),dimension(i1:i2,j1:j2,num_chem) :: var_rm
real(CHEM_KIND_R4),dimension(i1:i2,k1:k2,j1:j2)          :: var_rml
real(CHEM_KIND_R4),dimension(i1:i2,k1:k2,j1:j2)          :: var_rml_r
select case(ivar)
case (6)
  mvar=1
case (7)
  mvar=2
case (8)
  mvar=3
case (9)
  mvar=4
case default
  mvar=0
end select
if (mvar/=0)then
 do k=k1,k2
   var_rml(:,k,:)=var_rml_r(:,k2-k+1,:)
 end do
 do j=j1,j2
  do k=k1,k2
    do i=i1,i2
      ii=i-i1+1
      gsd3d(ii,j,k,mvar)=gsd3d(ii,j,k,mvar)+var_rml(i,k,j)
    enddo
  enddo
  do i=i1,i2
    ii=i-i1+1
    gsdcol(ii,j,mvar)=gsdcol(ii,j,mvar)+var_rm(i,j,ivar)
  end do
 enddo
! write(6,*)'mvar',mvar,maxval(gsdcol(:,:,mvar)),minval(gsdcol(:,:,mvar))
endif
return
end subroutine savewetdepr
#endif
