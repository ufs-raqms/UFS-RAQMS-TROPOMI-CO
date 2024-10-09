module dep_ctrans_grell_mod

  use chem_const_mod,  only : epsilc
  use chem_tracers_mod

  IMPLICIT NONE

  private
  public :: grelldrvct

CONTAINS

!-------------------------------------------------------------
   SUBROUTINE grelldrvct(DT,itimestep,                          &
              rho_phy,RAINCV,chem,trfall,                       &
              U,V,t_phy,moist,dz8w,p_phy,p8w,                   &
              pbl,XLV,CP,G,r_v,z,cu_co_ten,                     &
              numgas,chem_opt,                                  &
              num_chem,num_moist,tile,                          &
              ids,ide, jds,jde, kds,kde,                        &
              ims,ime, jms,jme, kms,kme,                        &
              its,ite, jts,jte, kts,kte                         )
!-------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------
   INTEGER,      INTENT(IN   ) ::                               &
                                  numgas,chem_opt,              &
                                  num_chem,num_moist, tile,     &
                                  ids,ide, jds,jde, kds,kde,    & 
                                  ims,ime, jms,jme, kms,kme,    & 
                                  its,ite, jts,jte, kts,kte

  
   INTEGER,      INTENT(IN   ) :: ITIMESTEP

   REAL,         INTENT(IN   ) :: XLV, R_v
   REAL,         INTENT(IN   ) :: CP,G

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme,num_moist ),   &
          INTENT(IN   ) ::                            moist
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                          U,    &
                                                          V,    &
                                                      t_phy,    &
                                                          z,    &
                                                      p_phy,    &
                                                        p8w,    &
                                                       dz8w,    &
                                                    rho_phy
!
! on output for control only, purely diagnostic
  REAL,  DIMENSION( ims:ime , jms:jme )                    ,    &
          INTENT(INOUT) ::                                      &
                                                        pbl
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(INOUT   ) ::                                   &
                                                  cu_co_ten


!
   REAL, INTENT(IN   ) :: DT
!
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme, num_chem ),    &
         INTENT(INOUT) ::                                       &
                                                       chem
   REAL, DIMENSION( ims:ime ,           jms:jme, num_chem ),    &
         INTENT(INOUT) ::                                       &
                                                     trfall

   REAL, DIMENSION( ims:ime , jms:jme ),                        &
         INTENT(IN) ::                               RAINCV

! LOCAL VARS
     real,    dimension (its:ite,kts:kte) ::                    &
        OUTT,OUTQ,OUTQC,rho,zz,dz
     real,    dimension (its:ite)         ::                    &
        pret, ter11

!
! basic environmental input includes moisture convergence (mconv)
! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
! convection for this call only and at that particular gridpoint
!
     real,    dimension (its:ite,kts:kte) ::                    &
        T,TN,q,qo,PO,P,US,VS,hstary
     real,    dimension (its:ite,kts:kte,num_chem) ::           &
           tracer,tracert
     real,    dimension (its:ite,num_chem) ::                   &
           trdep
     real, dimension (its:ite)            ::                    &
        xmb,Z1,PSUR,AAEQ
     integer, dimension (its:ite)            ::                 &
        kpbli,ktop,csum,ipr

   INTEGER :: nv,i,j,k,ICLDCK,jpr,npr
   REAL    :: tcrit,dp,dq,trmax
   INTEGER :: itf,jtf,ktf,iopt

   jpr=0
   npr=10
   tcrit=258.
   iopt=1
   itf=ite !lzhang
   ktf=kte !lzhang
   jtf=jte !lzhang
!                                                                      
!
     trfall(:,:,:)=0.
     DO 100 J = jts,jtf  
     DO I=ITS,ITF
         xmb(i)=0.
         csum(i)=0
         ktop(i)=0
         kpbli(i)=0
         !PSUR(I)=p_phy(I,kts,J)*.01 !lzhang
         PSUR(I)=p8w(I,kts,J)*.01
         TER11(I)=z(i,kts,j)
         aaeq(i)=0.
!
!   rainrate is input for this transport/wet-deposition routine
!
         pret(i)=raincv(i,j)/dt
         if(pret(i).le.0.)aaeq(i)=20.
     ENDDO
     DO K=kts,ktf
     DO I=ITS,ITF
         po(i,k)=p_phy(i,k,j)*.01
         P(I,K)=PO(i,k)
         US(I,K) =u(i,k,j)
         VS(I,K) =v(i,k,j)
         zz(I,K) =z(i,k,j)
         T(I,K)=t_phy(i,k,j)
         rho(I,K)=rho_phy(i,k,j)
         dz(I,K)=dz8w(i,k,j)
         q(I,K)=moist(i,k,j,p_qv)
         IF(Q(I,K).LT.1.E-08)Q(I,K)=1.E-08
         if (pbl(i,j).ge.z(i,k,j).and.pbl(i,j).lt.z(i,k+1,j)) then
         kpbli(i)=k
         endif
     ENDDO
     ENDDO
     do nv=1,num_chem
     DO I=ITS,ITF
      trdep(i,nv)=0.
     ENDDO
     DO K=kts,ktf
     DO I=ITS,ITF
         tracer(i,k,nv)=max(epsilc,chem(i,k,j,nv))
         tracert(i,k,nv)=0.
     ENDDO
     ENDDO
     ENDDO
     DO K=kts,ktf
     DO I=ITS,ITF
         cu_co_ten(i,k,j)=0.
!        hstary(i,k)=hstar4(nv)*exp(dhr(nv)*(1./t(i,k)-1./298.))
     ENDDO
     ENDDO
!    ENDDO
!
!---- CALL NON_RESOLVED CONVECTIVE TRANSPORT
      ipr(:)=0
#if 0
      DO I=ITS,ITF
         !if(pret(i)*3600. .gt.4) then
          !write(20,*)'zl0',tile,j,i,itimestep,pret(i)
           !write(20,*)'zl0'
         if(pret(i)*3600. .gt.4.and.j.eq.76 .and. (i.eq.87 ) .and. (itimestep .eq.8).and.(tile.eq.2)) then
           write(20,*)'zltile,j,i,pret(i),ter11(i) = ',tile,j,i,itimestep,pret(i)
           DO K=kts,ktf
              write(20,123)k,p(i,k),t(i,k),q(i,k),tracer(i,k,16)
           ENDDO
           ipr(i)=1
         endif
         !endif

      ENDDO
#endif

123  format(1x,i3,f7.1,1x,f6.1,2(1x,e13.4))
124  format(1x,i3,f7.1,1x,2(1x,e13.4))
!
      CALL CUP_gf(itimestep,tile,kpbli,ktop,tracer,trdep,aaeq,t,q,ter11,zz,pret,p,tracert,    &
           hstary,dt,psur,us,vs,tcrit,num_chem,chem_opt,numgas,      &
           rho,xlv,r_v,cp,g,ipr,its,ite,itf,kts,kte,ktf,csum,xmb )
#if 0
      DO I=ITS,ITF
      if(pret(i)*3600. .gt. 4.and.j.eq.76.and.i.eq.87.and.tile.eq.2.and.(itimestep .eq.8))then
           write(20,*)'zl9,j,i,pret(i),xmb(i) =,ktop(i)',j,i,itimestep,tile,xmb(i),ktop(i),pret(i)
           DO K=kts,ktf
              write(20,124)k,p(i,k),tracert(i,k,p_bc2),trdep(i,16)
           ENDDO
       endif
       enddo
#endif
       do nv=1,num_chem
         DO I=its,itf
           if(pret(i).le.0.)then
             DO K=kts,ktf
               tracert(i,k,nv)=0.
             ENDDO
           endif
         enddo
       enddo
!     CALL neg_check_ct(pret,ktop,epsilc,dt,tracer,tracert,iopt,num_chem,   &
!                       its,ite,kts,kte,itf,ktf,jpr,jpr,npr,j)
!lzhang diag trdep:
!      do i=its,itf
!        trdep(i,1)=pret(i)
!        trdep(i,2)=maxval(tracert(i,:,p_seas_2))
!      enddo

       do nv=1,num_chem
         do I=its,itf
           if(pret(i).gt.0.)then
             !wet deposition into ug/m2/s
             !trfall(i,j,nv)=trdep(i,nv) !lzhang
             do k=kts,kte
             trfall(i,j,nv)=trfall(i,j,nv)+(tracert(i,k,nv)*rho(i,k)*dz(i,k)) !lzhang 
             enddo
             trfall(i,j,nv)=max(0.,trfall(i,j,nv))  !lzhang
             trmax = maxval(chem(i,kts:ktf,j,nv))
             do k=kts,ktop(i)
               chem(i,k,j,nv)=min(trmax,chem(i,k,j,nv)+tracert(i,k,nv)*dt)
               chem(i,k,j,nv)=max(epsilc,chem(i,k,j,nv))
               if(nv.eq.npr)then
                 cu_co_ten(i,k,j)=tracert(i,k,npr)*dt
               endif
              
             enddo
           else
             do K=kts,ktf-1
               tracert(i,k,nv)=0.
               if(nv.eq.npr)cu_co_ten(i,k,j)=0.
             enddo
           endif
         enddo
       enddo

 100   continue

   END SUBROUTINE GRELLDRVCT

!this routine derived from the GF scheme, consistent when using GF, but can also
!be used with any other convective parameterization
!
   SUBROUTINE CUP_gf(ktau,tile,kpbli,ktop,chem,chem_psum,aaeq,t,q,z1,z,pre,p,outc,    &
              hstary,DTIME,PSUR,US,VS,TCRIT,numc,chem_opt,numgas, &
              rho,xlv,r_v,cp,g,ipr,its,ite,itf,kts,kte,ktf,csum,xmb )
   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        numc,chem_opt,numgas,its,ite,itf,kts,kte,ktf,tile,ktau
  !
  ! 
  !
  ! outc  = output c tendency (per s)
  ! pre    = input precip
     real,    dimension (its:ite,kts:kte,numc)                          &
        ,intent (inout  )                   ::                           &
         chem,outc
     real,    dimension (its:ite,numc)                              &
        ,intent (inout  )                   ::                           &
        chem_psum
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                           &
        pre,aaeq
     integer,    dimension (its:ite)                                   &
        ,intent (inout  )                   ::                           &
        ktop
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,Q,T,P,US,VS,hstary
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout   )                   ::                           &
        z
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        Z1,PSUR
       real                                                            &
        ,intent (in   )                   ::                           &
        xlv,r_v,cp,g,dtime,tcrit
     integer,    dimension (its:ite), intent(in) :: ipr,kpbli,csum

!
!  local ensemble dependent variables in this routine
!
     real,    dimension (its:ite,kts:kte,numc)   ::                           &
         chem_cup,chem_up,chem_down,dellac,chem_c,chem_pw,chem_pwd,chem_pwout
     real,    dimension (its:ite,numc)   ::                           &
         chem_pwav
     real,    dimension (its:ite,kts:kte,numc)   :: chem_pwdout
     real,    dimension (its:ite,1)  ::                         &
        xaa0_ens
     real,    dimension (1)  ::                                 &
        mbdt_ens
     real,    dimension (1) ::                                 &
        edt_ens
     real,    dimension (its:ite,1) ::                         &
        edtc
     real,    dimension (its:ite,kts:kte,1) ::                 &
        dellat_ens,dellaqc_ens,dellaq_ens,pwo_ens
!
!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level

     real,    dimension (its:ite,kts:kte) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,                      &
        c0t,pwdper,cd,cdd,                                                 & 
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        hcot,dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       edt,HKB,QKB,XMB,PWAV,PWEV,PWEVO,BU,BUD,cap_max,                 &
       cap_max_increment,closure_n,psum,psumh,sig,sigd,zuhe,           &
       ztexec,zqexec
     real,    dimension (its:ite) ::                                   &
        edtmax,edtmin,entr_rate
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,K22,KB,JMIN,kstabi,kstabm,        &   !-lxz
       kbcon,ktopdby,KBMAX,ierr

     integer                              ::                           &
       iloop,nall,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,dzo,mbdt,radius, &
      zcutdown,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,frh,radiusd,frhd
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: zusafe,xxx,xx1,xx2,zutop,power_entr,dzm1,dzp1


     integer :: k1,k2,kbegzu,kdefi,kfinalzu,kstart,jmini
     logical :: keep_going,flg(its:ite)
     real tot_time_hr,xff_shal(9),blqe,xkshal
     character*50 :: ierrc(its:ite) 
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr,c1d 
     real,    dimension (kts:kte) :: smth
     real :: xff_mid(its:ite,2)
! rainevap from sas
      real, dimension (its:ite) :: ccn,deltv,rntot,delqev,delq2,delt,delq,qevap,rn,qcond
      real :: alpha,trcc,rain,t1,q1,elocp,evef,el2orc,evfact,evfactl,g_rain,e_dn,c_up
!srf- begin - DICYCLE
      real :: beta,denom,h_entr,c1,ccnclean
      integer :: autoconv,aeroevap,use_excess,nv,masscon
      autoconv=1
      aeroevap=1
      c1=.001
      ccnclean=150.
      ccn(:)=150.
      use_excess=0
      ztexec(:)=0.
      zqexec(:)=0.
!srf- end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!proportionality constant to estimate pressure gradient of updraft (Zhang and Wu, 2003, JAS
      day=86400.
! ! no restriction if pre incoming is non zero
      cap_maxs=350.
      do i=its,itf
        edt(i)=0.
        xmb(i)=0.
        ierr(i)=0
        cap_max(i)=cap_maxs
        cap_max_increment(i) = 20.
        ierrc(i)=" "
      enddo
!
!--- initial entrainment rate (these may be changed later on in the
!--- program
!
      do i=its,itf
         c1d(i,:)= c1 ! 0. ! c1 ! max(.003,c1+float(csum(i))*.0001)
         entr_rate(i)=7.e-5
      enddo

      
!
!--- initial detrainmentrates
!
      do k=kts,ktf
      do i=its,itf
        zu(i,k)=0.
        zd(i,k)=0.
        hcot(i,k)=0.
        he(i,k)=0.
        cd(i,k)=1.e-9 ! 1.*entr_rate
        cdd(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!
      edtmax(:)=1.
      edtmin(:)=.1
!
!--- minimum depth (m), clouds must have
!
      depth_min=1000.
!
!--- maximum depth (mb) of capping 
!--- inversion (larger cap = no convection)
!
      DO i=its,itf
        kbmax(i)=1
        edt(i)=0.
        kstabm(i)=ktf-1
      enddo

!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1000.
!
      edt_ens(1)=.94
!
!--- environmental conditions, FIRST HEIGHTS
!
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,-1,   &
           xlv,r_v,cp,g,itf,ktf, &
           its,ite,  kts,kte)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup, &
           he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,  &
           xlv,r_v,cp,g,ierr,z1,itf,ktf,its,ite,  kts,kte)
      do i=its,itf
        if(ierr(i).eq.0)then
!          if (pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)then
!            write(20,*)'zl10,he,hes ',i
!            do k=kts,ktf
!              write(20,*)'zl11',k,z(i,k),he_cup(i,k),hes_cup(i,k),ierr(i)
!            enddo
!          endif

        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!
      do k=kts,ktf
        if(z_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
!--- level where detrainment for downdraft starts
!
      do k=kts,ktf
        if(z_cup(i,k).gt.z_detr+z1(i))then
          kdet(i)=k
          go to 26
        endif
      enddo
 26   continue
!
      endif
      enddo
!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      k22(:)=1
      ktopdby(:)=0
      CALL cup_MAXIMI(HE_CUP,1,KBMAX,K22,ierr, &
           itf,ktf, its,ite,  kts,kte)
       DO 36 i=its,itf
          !if (pre(i)*3600. .gt.4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)then
          !  write(20,*)'zl12,i,k22,kbmax = ',i,k22(i),kbmax(i),ierr(i)
          !endif
         IF(ierr(I).eq.0.)THEN
         IF(K22(I).GE.KBMAX(i))then
           ierr(i)=2
           ierrc(i)="could not find k22"
           ktop(i)=0
           k22(i)=0
           kbcon(i)=0
         endif
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      do i=its,itf
       IF(ierr(I).eq.0.)THEN
         hkb(i)=he_cup(i,k22(i))
       endif ! ierr
      enddo
      iloop=1
      call cup_kbcon(ierrc,cap_max_increment,iloop,k22,kbcon,he_cup,hes_cup, &
           hkb,ierr,kbmax,p_cup,cap_max, &
           ztexec,zqexec,use_excess,       &
           itf,ktf, &
           its,ite,  kts,kte, &
           z_cup,entr_rate,he,0)
      !do i=its,itf
      !    if (pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)then
      !      write(20,*)'zl13,i,k22,kbcon = ',i,k22(i),kbcon(i),ierr(i)
      !    endif
      ! enddo
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEs_cup,Kbcon,kstabm,kstabi,ierr,  &
           itf,ktf, its,ite,  kts,kte)
      DO i=its,itf
         if(kstabi(i).lt.kbcon(i))then
           kbcon(i)=1
           ierr(i)=42
         endif
         do k=kts,ktf
            entr_rate_2d(i,k)=entr_rate(i)
         enddo
         IF(ierr(I).eq.0.)THEN
            kdefi=0
            do k=kts,ktf
               frh = min(q_cup(i,k)/qes_cup(i,k),1.)
               entr_rate_2d(i,k)=entr_rate(i) *(1.3-frh)
            enddo
          endif
       enddo
      DO i=its,itf
         IF(ierr(I).eq.0.)THEN
           hkb(i)=he_cup(i,k22(i))
           if(kbcon(i).lt.2)then
!              write(0,*)"kbcon",kbcon(i),k22(i),ierr(i)
              kbcon(i)=2
              ierr(i)=43
          endif
         endif
       enddo
!
! get entrainment and detrainmentrates for updraft
!
      call rates_up_pdf('deep',ktop,ierr,p_cup,entr_rate_2d,hkb,he,hes_cup,z_cup, &
           kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zu,kpbli,ktopdby,csum)
      do i=its,itf
          !if (ipr(i).eq.1)then
          if (ipr(i).eq.1.and.pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)then
            write(20,*)'ktop,ierr,kpbli,csum,kstabi = ',ktop(i),ierr(i),kpbli(i),csum(i),kstabi(i)
            do k=kts,ktf
              write(20,*)k,z_cup(i,k),zu(i,k)
            enddo
          endif
       enddo

!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         hc(i,k)=0.
         DBY(I,K)=0.
         up_massentr(i,k)=0.
         up_massdetr(i,k)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0.)THEN
         do k=1,max(1,k22(i) -1)
            hc(i,k)=he_cup(i,k)
         enddo
       endif ! ierr
      enddo
!
!
      do i=its,itf
         if(ierr(i).eq.0)then
         hc(i,k22(i))=hkb(i)
!-- formulation 2
         if(k22(i).gt.1)then
            do k=1,k22(i) -1
              zu(i,k)=0.
            enddo
         endif
         do k=ktop(i)+1,ktf
           zu(i,k)=0.
          enddo

         !- mass entrainment and detrinament is defined on model levels
        
         do k=2,maxloc(zu(i,:),1)
         !=> below maximum value zu -> change entrainment
           dz=z_cup(i,k)-z_cup(i,k-1)
        
           up_massdetr(i,k-1)=cd(i,k-1)*dz*zu(i,k-1)
           up_massentr(i,k-1)=zu(i,k)-zu(i,k-1)+up_massdetr(i,k-1)
           if(up_massentr(i,k-1).lt.0.)then
              up_massentr(i,k-1)=0.
              up_massdetr(i,k-1)=zu(i,k-1)-zu(i,k)
              if(zu(i,k-1).gt.0.)cd(i,k-1)=up_massdetr(i,k-1)/(dz*zu(i,k-1))
           endif
           if(zu(i,k-1).gt.0.)entr_rate_2d(i,k-1)=(up_massentr(i,k-1))/(dz*zu(i,k-1))
         enddo
         do k=maxloc(zu(i,:),1)+1,ktop(i)
         !=> above maximum value zu -> change detrainment
           dz=z_cup(i,k)-z_cup(i,k-1)
           up_massentr(i,k-1)=entr_rate_2d(i,k-1)*dz*zu(i,k-1)
           up_massdetr(i,k-1)=zu(i,k-1)+up_massentr(i,k-1)-zu(i,k)
           if(up_massdetr(i,k-1).lt.0.)then
              up_massdetr(i,k-1)=0.
              up_massentr(i,k-1)=zu(i,k)-zu(i,k-1)
              if(zu(i,k-1).gt.0.)entr_rate_2d(i,k-1)=(up_massentr(i,k-1))/(dz*zu(i,k-1))
           endif
        
           if(zu(i,k-1).gt.0.)cd(i,k-1)=up_massdetr(i,k-1)/(dz*zu(i,k-1))
         enddo
        
!
!   NOTE: Ktop here already includes overshooting, ktopdby is without
!   overshooting
!
        do k=k22(i)  +1,ktop(i)  !mass cons option
          denom=zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1)
          if(denom.lt.1.e-8)then
           ierr(i)=52
           exit
          endif

          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)
         enddo
41       continue
         if(ktop(i).lt.kbcon(i)+2)then
            ierr(i)=5
            ierrc(i)='ktop too small'
            ktop(i)=0
         endif
         do k=ktop(i)+1,ktf
           HC(i,K)=hes_cup(i,k)
           DBY(I,K)=0.
           zu(i,k)=0.
           cd(i,k)=0.
           entr_rate_2d(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
         enddo
      endif
      enddo
!
      DO 37 i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(z_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            ! if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)'zl14',zktop,z_cup(i,ktop(i))
            do k=kts,ktf
              if(z_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 kzdown(i)=min(kzdown(i),kstabi(i)-1)  !
                 go to 37
              endif
              enddo
         endif
 37   CONTINUE
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
      call cup_minimi(HEs_cup,K22,kzdown,JMIN,ierr, &
           itf,ktf, its,ite,  kts,kte)
      DO 100 i=its,itf
        ! if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)'zl15,minim',kzdown(i),jmin(i),ierr(i)
         IF(ierr(I).eq.0.)THEN
!
!--- check whether it would have buoyancy, if there where
!--- no entrainment/detrainment
!
         jmini = jmin(i)
         keep_going = .TRUE.
         do while ( keep_going )
           keep_going = .FALSE.
           if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
           if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
           ki = jmini
           hcd(i,ki)=hes_cup(i,ki)
           DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
           dh=0.
           do k=ki-1,1,-1
             hcd(i,k)=hes_cup(i,jmini)
             DZ=Z_cup(i,K+1)-Z_cup(i,K)
             dh=dh+dz*(HCD(i,K)-hes_cup(i,k))
             if(dh.gt.0.)then
               jmini=jmini-1
               if ( jmini .gt. 5 ) then
                 keep_going = .TRUE.
               else
                 ierr(i) = 9
                 ierrc(i) = "could not find jmini9"
                 exit
               endif
             endif
           enddo
         enddo
         jmin(i) = jmini 
         if ( jmini .le. 5 ) then
           ierr(i)=4
           ierrc(i) = "could not find jmini4"
         endif
       ENDIF
100   continue
!
! - Must have at least depth_min m between cloud convective base
!     and cloud top.
!
      do i=its,itf
        ! if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)'zl16,jmin = ',jmin(i),i
         IF(ierr(I).eq.0.)THEN
            if ( jmin(i) - 1 .lt. kdet(i)   ) kdet(i) = jmin(i)-1
            IF(-z_cup(I,KBCON(I))+z_cup(I,KTOP(I)).LT.depth_min)then
               ierr(i)=6
               ierrc(i)="cloud depth very shallow"
            endif
         endif
      enddo
!
!--- normalized downdraft mass flux profile,also work on bottom detrainment
!--- in this routine
!
      do i=its,itf
      do k=kts,ktf
       zd(i,k)=0.
       cdd(i,k)=0.
       dd_massentr(i,k)=0.
       dd_massdetr(i,k)=0.
       hcd(i,k)=hes_cup(i,k)
       dbyd(i,k)=0.
       mentrd_rate_2d(i,k)=entr_rate(i)
        !if(pre(i)*3600. .gt. 0.1.and.i.eq.64)write(0,*)'zl,go into pdf routine1',kdet(i),jmin(i),kpbli(i),ierr(i)
      enddo
        !if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)'zlp00,go into pdf routine1',i,pre(i),kdet(i),jmin(i),kpbli(i),ierr(i)
      enddo 

      do i=its,itf
        beta=max(.035,.055-float(csum(i))*.0015)  !.02
        bud(i)=0.
        !if(ipr(i).eq.1)write(0,*)'zl,go into pdf routine2',kdet(i),jmin(i),beta,kpbli(i),ierr(i)
        !if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)'zl0,go into pdf routine2',kdet(i),jmin(i),beta,kpbli(i),ierr(i)
        IF(ierr(I).eq.0)then
        cdd(i,1:jmin(i))=1.e-9
        cdd(i,jmin(i))=0.
        !dd_massdetr(i,:)=0.
        !dd_massentr(i,:)=0.
           call get_zu_zd_pdf_fim("DOWN",ierr(i),kdet(i),jmin(i),zd(i,:),kts,kte,ktf,beta,kpbli(i),csum(i))
         !if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)then
         !  do k=kts,ktop(i)
         !  write(20,*)'zl1',kdet(i),jmin(i),kpbli(i),zd(i,jmin(i)),zd(i,jmin(i)-1)
         !  enddo
         !endif

        if(zd(i,jmin(i)) .lt.1.e-8)then
          zd(i,jmin(i))=0.
          jmin(i)=jmin(i)-1
          if(zd(i,jmin(i)) .lt.1.e-8)then
             ierr(i)=876
             cycle
          endif
        endif
        
        do ki=jmin(i)  ,maxloc(zd(i,:),1),-1
!-srf mass cons
          !=> from jmin to maximum value zd -> change entrainment
          dzo=z_cup(i,ki+1)-z_cup(i,ki)
          dd_massdetr(i,ki)=cdd(i,ki)*dzo*zd(i,ki+1)
          dd_massentr(i,ki)=zd(i,ki)-zd(i,ki+1)+dd_massdetr(i,ki)
          if(dd_massentr(i,ki).lt.0.)then
             dd_massentr(i,ki)=0.
             dd_massdetr(i,ki)=zd(i,ki+1)-zd(i,ki)
             if(zd(i,ki+1).gt.0.)cdd(i,ki)=dd_massdetr(i,ki)/(dzo*zd(i,ki+1))
          endif
          if(zd(i,ki+1).gt.0.)mentrd_rate_2d(i,ki)=dd_massentr(i,ki)/(dzo*zd(i,ki+1))
        enddo
        mentrd_rate_2d(i,1)=0.
        do ki=maxloc(zd(i,:),1)-1,1,-1
          !=> from maximum value zd to surface -> change detrainment
          dzo=z_cup(i,ki+1)-z_cup(i,ki)
          dd_massentr(i,ki)=mentrd_rate_2d(i,ki)*dzo*zd(i,ki+1)
          dd_massdetr(i,ki) = zd(i,ki+1)+dd_massentr(i,ki)-zd(i,ki)
          if(dd_massdetr(i,ki).lt.0.)then
            dd_massdetr(i,ki)=0.
            dd_massentr(i,ki)=zd(i,ki)-zd(i,ki+1)
            if(zd(i,ki+1).gt.0.)mentrd_rate_2d(i,ki)=dd_massentr(i,ki)/(dzo*zd(i,ki+1))
          endif
          if(zd(i,ki+1).gt.0.)cdd(i,ki)= dd_massdetr(i,ki)/(dzo*zd(i,ki+1))
        enddo

        !do k=kts,jmin(i)
        ! if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)'zl2',k,jmin(i),zd(i,k),hcd(i,jmin(i)),hes_cup(i,jmin(i))
        !enddo
        do k=kts,ktop(i)+1
          if(p_cup(i,k).gt.600.)then
            c1d(i,k)=0.
          elseif(p_cup(i,k).gt.400. .and. p_cup(i,k).le.600.)then
            c1d(i,k)=c1*(600.-p_cup(i,k))/200.
            c1d(i,k)=max(0.,c1d(i,k))
          elseif(p_cup(i,k).le.400.)then
            c1d(i,k)=c1
          endif
        enddo


! downdraft moist static energy + moisture budget
            dbyd(i,jmin(i))=hcd(i,jmin(i))-hes_cup(i,jmin(i))
            bud(i)=dbyd(i,jmin(i))*(z_cup(i,jmin(i)+1)-z_cup(i,jmin(i)))

       ! if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)'zlts1', dbyd(i,jmin(i)),bud(i),z_cup(i,jmin(i)+1),z_cup(i,jmin(i))

            do ki=jmin(i)  ,1,-1
!-srf mass cons
             dzo=z_cup(i,ki+1)-z_cup(i,ki)
!             h_entr=heo(i,ki)
             h_entr=.5*(he(i,ki)+.5*(hc(i,ki)+hc(i,ki+1)))
             hcd(i,ki)=(hcd(i,ki+1)*zd(i,ki+1)                       &
                         -.5*dd_massdetr(i,ki)*hcd(i,ki+1)+ &
                        dd_massentr(i,ki)*h_entr)   /            &
                        (zd(i,ki+1)-.5*dd_massdetr(i,ki)+dd_massentr(i,ki))
             dbyd(i,ki)=hcd(i,ki)-hes_cup(i,ki)
             !if(ipr(i).eq.1)write(0,*)'zl,ki,bud = ',ki,bud(i),hcd(i,ki)
        !     if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)'zl3,ki,bud = ',ki,bud(i),hcd(i,ki)
             bud(i)=bud(i)+dbyd(i,ki)*dzo
            enddo
          endif

        if(bud(i).gt.0)then
          ierr(i)=7
          ierrc(i)='downdraft is not negatively buoyant '
        endif
      enddo
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture_new(ierrc,zd,hcd,hes_cup,qcd,qes_cup, &
           pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,gamma_cup, &
           pwev,bu,qrcd,q,he,t_cup,1,xlv,itf,ktf, its,ite,  kts,kte)

!
!--- calculate moisture properties of updraft
!
      call cup_up_moisture('deep',ierr,z_cup,qc,qrc,pw,pwav, &
           ccnclean,p_cup,kbcon,ktop,cd,dby,clw_all, &
           t_cup,q,GAMMA_cup,zu,qes_cup,k22,q_cup,        &
           ZQEXEC,use_excess,ccn,rho,c1d,up_massentr,up_massdetr,psum,psumh,&
           xlv,c0t,autoconv,aeroevap,1,itf,ktf,ipr, &
           its,ite,  kts,kte)
!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
           pw,ccn,pwev,edtmax,edtmin,1,edtc,psum,psumh, &
           ccnclean,rho,aeroevap,itf,ktf,ipr, &
           its,ite,  kts,kte)
! initialize tracers if they exist
!
           chem_pwav(:,:)=0.
           chem_psum(:,:)=0.
           chem_pw(:,:,:)=0.
           chem_pwd(:,:,:)=0.
           chem_pwout(:,:,:)=0.
           chem_pwdout(:,:,:)=0.
           pwdper(:,:)=0.
           chem_down(:,:,:)=0.
           chem_up(:,:,:)=0.
           chem_c(:,:,:)=0. !lzhang
           if(numc.gt.0)then
             do i=its,itf
               if(ierr(i).eq.0)then
                 !if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)"zl4,pwdper",pwev(i),pwav(i),edtc(i,1)
                do k=kts,jmin(i)
                 pwdper(i,k)=-edtc(i,1)*pwd(i,k)/pwav(i)
                 !if(pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)write(20,*)"zl5,pwdper",k,pwd(i,k)
                enddo
                edt(i)=edtc(i,1)
              do nv=1,numc
!                alpha=0.5
!                if(nv.eq.p_bc1 .or. nv.eq.p_oc1 .or. nv.eq.p_dms)alpha=0.
!                !if(nv.eq.p_bc2 .or. nv.eq.p_oc2)alpha=0.6
!                if(nv.eq.p_bc2 .or. nv.eq.p_oc2)alpha=0.5  !lzhang
!                if(nv.eq.p_sulf .or. nv.eq.p_seas_1 .or. nv.eq.p_seas_2 .or. &
!                nv.eq.p_seas_3 .or. nv.eq.p_seas_4)then
!                   alpha=1.
!                endif

                    alpha = 0.
   if(chem_opt >= 300 .and. chem_opt < 500)then
# if 0
    if(nv.gt. numgas .or. nv.eq.p_sulf) then
     alpha = .5    ! scavenging factor
       !if(nv.eq.p_bc1 .or. nv.eq.p_oc1 .or. nv.eq.p_dms) alpha=0.
       if(nv.eq.p_bc1 .or. nv.eq.p_oc1) alpha=0.4 !lzhang
       if(nv.eq.p_dms) alpha=0.
       if(nv.eq.p_sulf .or. nv.eq.p_seas_1 .or. nv.eq.p_seas_2 .or. &
          nv.eq.p_seas_3 .or. nv.eq.p_seas_4 .or. nv.eq.p_seas_5)alpha=1.
       !if(nv.eq.p_bc2 .or. nv.eq.p_oc2)alpha=0.8
       if(nv.eq.p_bc2 .or. nv.eq.p_oc2)alpha=0.5  !lzhang
     endif
# endif
    alpha = 1. 
   endif

!    if(chem_opt == 301 .or. chem_opt==108 ) then
!       if(nv.eq.p_hno3) alpha=1.
!       if(nv .eq.p_h2o2) alpha=0.5
!    endif

   if (chem_opt==108) then
    if(nv.gt. numgas .or. nv.eq.p_sulf) then
        alpha = .5
     if (nv.eq.p_naai.or.nv.eq.p_naaj.or.nv.eq.p_claj.or.nv.eq.p_clai.or. &
           nv.eq.p_nh4aj.or.nv.eq.p_nh4ai.or.nv.eq.p_no3aj.or.nv.eq.p_no3ai.or.& 
      nv.eq.p_so4ai.or.nv.eq.p_so4aj.or.nv.eq.p_seas) alpha=1.

     endif
   endif

    if(chem_opt == 301 .or. chem_opt==108 ) then
       if(nv.eq.p_hno3) alpha=1.
       if(nv .eq.p_h2o2) alpha=0.5
    endif



                do k=kts+1,ktf
                   chem_cup(i,k,nv)=.5*(chem(i,k-1,nv)+chem(i,k,nv))
                enddo
                chem_cup(i,kts,nv)=chem(i,kts,nv)
  !
! in updraft  
!
                do k=1,k22(i)
                   chem_up(i,k,nv)=chem_cup(i,k,nv)
                enddo
                do k=k22(i)+1,ktop(i)
                   chem_up(i,k,nv)=(chem_up(i,k-1,nv)*zu(i,k-1) &
                             -.5*up_massdetr(i,k-1)*chem_up(i,k-1,nv)+ &
                         up_massentr(i,k-1)*chem(i,k-1,nv))   /            & !lzhang chem
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
                   chem_c(i,k,nv)=alpha*chem_up(i,k,nv)
               !    if(nv.eq.p_sulf .or. nv.eq.p_seas_1 .or. nv.eq.p_seas_2 .or. &
               !       nv.eq.p_seas_3 .or. nv.eq.p_seas_4)then
               !        alpha=1.
               !    endif
                   dz=z_cup(i,K)-z_cup(i,K-1)
                   trcc=chem_c(i,k,nv)/(1.+c0t(i,k)*dz)
                   chem_pw(i,k,nv)=c0t(i,k)*dz*trcc*zu(i,k)
                   chem_pwav(i,nv)=chem_pwav(i,nv)+chem_pw(i,k,nv)
                   chem_up(i,k,nv)=trcc+chem_up(i,k,nv)-chem_c(i,k,nv)

                  chem_pwout(i,k,nv)=chem_pw(i,k,nv)
                enddo
                do k=ktop(i)+1,ktf
                 chem_up(i,k,nv)=chem_cup(i,k,nv)
                enddo
!
! in downdraft
!
                chem_down(i,jmin(i)+1,nv)=chem_cup(i,jmin(i)+1,nv)
                do ki=jmin(i)  ,1,-1
                  chem_down(i,ki,nv)=(chem_down(i,ki+1,nv)*zd(i,ki+1)  &
                         -.5*dd_massdetr(i,ki)*chem_down(i,ki+1,nv)+ &
                        dd_massentr(i,ki)*chem(i,ki,nv))   /            &
                        (zd(i,ki+1)-.5*dd_massdetr(i,ki)+dd_massentr(i,ki))
                  chem_down(i,ki,nv)=chem_down(i,ki,nv)+pwdper(i,ki)*chem_pwav(i,nv) 
                  chem_pwd(i,ki,nv)=max(0.,pwdper(i,ki)*chem_pwav(i,nv))
                  chem_psum(i,nv)=chem_psum(i,nv)+chem_pw(i,ki,nv)-chem_pwd(i,ki,nv)
                  chem_pwdout(i,ki,nv)=chem_pwd(i,ki,nv)
                enddo
!
              enddo ! numc
#if 0
                 if(pre(i)*3600. .gt. 4.and.i.eq.87 .and.tile.eq.2.and.ktau.eq.8)then
                   nv=p_bc2
                   do k=kts,ktop(i)
                   write(20,*)"zl6,chem_up",k,chem(i,k,16),chem_up(i,k,16),chem_pw(i,k,16)
                   enddo
                   do k=kts,jmin(i)
                   write(20,*)"zl7,chem_down",k,chem_down(i,k,16),chem_pwd(i,k,16)
                   enddo
                 endif
#endif
            endif ! ierr=0
          enddo ! i
        endif ! numc=0
!
!
!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level 3
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1

     if (numc.gt.0)then
      do nv=1,numc
        do k=kts,ktf
        do i=its,itf
          dellac(i,k,nv)=0.
        enddo
        enddo
        do i=its,itf
          if(ierr(i).eq.0)then
            dp=100.*(p_cup(i,1)-p_cup(i,2))
            G_rain=  0.5*(chem_pw (i,1,nv)+chem_pw (i,2,nv))*g/dp
            E_dn  =  0.5*(chem_pwd(i,1,nv)+chem_pwd(i,2,nv))*g/dp
            dellac(i,1,nv)=(edt(i)*zd(i,2)*chem_down(i,2,nv)   &
                     -edt(i)*zd(i,2)*chem_cup(i,2,nv)&
                     -zu(i,2)*(chem_up(i,2,nv)-chem_cup(i,2,nv)) )*g/dp &
                         - g_rain + E_dn
            do k=kts+1,ktop(i)
               dp=100.*(p_cup(i,k)-p_cup(i,k+1))
            ! these three are only used at or near mass detrainment and/or
            ! entrainment levels
	    G_rain=  0.5*(chem_pw (i,k,nv)+chem_pw (i,k+1,nv))*g/dp
	    E_dn  =  0.5*(chem_pwd(i,k,nv)+chem_pwd(i,k+1,nv))*g/dp  ! pwdo < 0 and E_dn must > 0
            !-- condensation source term = detrained + flux divergence of
            !-- cloud liquid water (qrco) + converted to rain
	
!            C_up = dellaqc(i,k)+(zu(i,k+1)* qrco(i,k+1) -       &
!                                 zu(i,k  )* qrco(i,k  )  )*g/dp + G_rain

               dellac(i,k,nv) =-(zu(i,k+1)*(chem_up(i,k+1,nv)-chem_cup(i,k+1,nv) ) - &
                    zu(i,k  )*(chem_up (i,k,nv  )-chem_cup(i,k,nv  ) ) )*g/dp &
                    +(zd(i,k+1)*(chem_down(i,k+1,nv)-chem_cup(i,k+1,nv) ) - &
                    zd(i,k  )*(chem_down(i,k,nv  )-chem_cup(i,k,nv)))*g/dp*edt(i) & 
                         - g_rain + E_dn
#if 0
      if(pre(i)*3600. .gt. 4.and.i.eq.87 .and.tile.eq.2.and.ktau.eq.8.and.nv.eq.16)then
      write (20,*),'zln8',zu(i,k+1),chem_up(i,k+1,nv),chem_cup(i,k+1,nv)
      write (20,*),'zln9',zu(i,k  ),chem_up(i,k,nv),chem_cup(i,k,nv),zd(i,k+1)
      write (20,*),'zln10',chem_down(i,k+1,nv),chem_cup(i,k+1,nv),zd(i,k  )
      write (20,*),'zln11',chem_down(i,k,nv  ),chem_cup(i,k,nv),g/dp,edt(i)
      write (20,*),'zln12',dellac(i,k,nv),g_rain,E_dn
      write (20,*),'zln13',pre(i),pwav(i)+edt(i)*pwev(i),xmb(i)
      endif
#endif

            enddo

         endif ! ierr
       enddo ! i
       enddo ! numc loop
     endif ! numc, chemstry check

444   format(1x,i2,1x,6e12.4) !,1x,f7.2,2x,e13.5)
!
!--- using dellas, calculate changed environmental profiles
!
             if(numc.gt.0)then
              do i=its,itf
                 if(ierr(i).eq.0)then
                    denom=max(pre(i),pwav(i)+edt(i)*pwev(i))
                    xmb(i)=pre(i)/max(epsilc,denom)
                    do nv=1,numc
                      do k=kts,ktop(i)
                          outc(i,k,nv)=dellac(i,k,nv)*xmb(i)
                      enddo
                      chem_psum(i,nv)=chem_psum(i,nv)*xmb(i) !lzhang
                    enddo
                 endif ! ierr
              enddo
             endif ! numc
      do i=its,itf
         ! if (pre(i)*3600. .gt. 4.and.i.eq.87.and.tile.eq.2.and.ktau.eq.8)then
         !   write(20,*)'zl8,ierr(i) = ',ierr(i),xmb(i)
         ! endif
       enddo


!
!---------------------------done------------------------------
!

   END SUBROUTINE CUP_gf


   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pw,ccn,pwev,edtmax,edtmin,maxens2,edtc,psum2,psumh, &
              ccnclean,rho,aeroevap,itf,ktf,ipr,         &
              its,ite,  kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        aeroevap,itf,ktf,           &
        its,ite,  kts,kte
     integer, intent (in   )              ::                           &
        maxens2
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,us,vs,z,p,pw
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev,ccn,psum2,psumh,edtmax,edtmin
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon,ipr
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer i,k,kk
     real    einc,pef,pefb,prezk,zkbc
     real,    dimension (its:ite)         ::                           &
      vshear,sdp,vws
     real :: prop_c,pefc,aeroadd,alpha3,beta3,rhoc
     prop_c=8. !10.386
     alpha3 = 1.9
     beta3  = -1.13
     pefc=0.

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
! */ calculate an average wind shear over the depth of the cloud
!
       do i=its,itf
        edt(i)=0.
        vws(i)=0.
        sdp(i)=0.
        vshear(i)=0.
       enddo
       do k=1,maxens2
       do i=its,itf
        edtc(i,k)=0.
       enddo
       enddo
       do kk = kts,ktf-1
         do 62 i=its,itf
          IF(ierr(i).ne.0)GO TO 62
          if (kk .le. min0(ktop(i),ktf) .and. kk .ge. kbcon(i)) then
             vws(i) = vws(i)+ &
              (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
          +   abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
              (p(i,kk) - p(i,kk+1))
            sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
          endif
          if (kk .eq. ktf-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
   62   continue
       end do
      do i=its,itf
         IF(ierr(i).eq.0)then
            pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
               -.00496*(VSHEAR(I)**3))
            if(pef.gt.0.9)pef=0.9
            if(pef.lt.0.1)pef=0.1
!
!--- cloud base precip efficiency
!
            zkbc=z(i,kbcon(i))*3.281e-3
            prezk=.02
            if(zkbc.gt.3.)then
               prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
               *(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
            endif
            if(zkbc.gt.25)then
               prezk=2.4
            endif
            pefb=1./(1.+prezk)
!           write(11,*)pefb,prezk,zkbc
            if(pefb.gt.0.9)pefb=0.9
            if(pefb.lt.0.1)pefb=0.1
            EDT(I)=1.-.5*(pefb+pef)
            if(aeroevap.gt.1)then
               aeroadd=(ccnclean**beta3)*((psumh(i))**(alpha3-1)) !*1.e6
!              prop_c=.9/aeroadd
               prop_c=.5*(pefb+pef)/aeroadd
               aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1)) !*1.e6
               aeroadd=prop_c*aeroadd
               pefc=aeroadd
               if(pefc.gt.0.9)pefc=0.9
               if(pefc.lt.0.1)pefc=0.1
               EDT(I)=1.-pefc
               if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
            endif


!--- edt here is 1-precipeff!
            einc=.2*edt(i)
            do k=1,maxens2
                edtc(i,k)=edt(i)+float(k-2)*einc
            enddo
         endif
      enddo
      do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               EDTC(I,K)=-EDTC(I,K)*pwav(I)/PWEV(I)
               IF(EDTC(I,K).GT.edtmax(i))EDTC(I,K)=edtmax(i)
               IF(EDTC(I,K).LT.edtmin(i))EDTC(I,K)=edtmin(i)
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_moisture_new(ierrc,zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,           &
              xlv,itf,ktf,                     &
              its,ite,  kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,ktf,           &
                                  its,ite,  kts,kte
  ! cdd= detrainment function 
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
        dd_massentr,dd_massdetr,gamma_cup,q,he 
     real,    intent (in   ) :: xlv
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qcd,qrcd,pwd
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwev,bu
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dh,dz,dqeva

      do i=its,itf
         bu(i)=0.
         pwev(i)=0.
      enddo
      do k=kts,ktf
      do i=its,itf
         qcd(i,k)=0.
         qrcd(i,k)=0.
         pwd(i,k)=0.
      enddo
      enddo
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      DZ=Z_cup(i,K+1)-Z_cup(i,K)
      qcd(i,k)=q_cup(i,k)
      DH=HCD(I,k)-HES_cup(I,K)
      if(dh.lt.0)then
        QRCD(I,K)=(qes_cup(i,k)+(1./XLV)*(GAMMA_cup(i,k) &
                  /(1.+GAMMA_cup(i,k)))*DH)
        else
          qrcd(i,k)=qes_cup(i,k)
        endif
      pwd(i,jmin(i))=zd(i,jmin(i))*min(0.,qcd(i,k)-qrcd(i,k))
      qcd(i,k)=qrcd(i,k)
      pwev(i)=pwev(i)+pwd(i,jmin(i)) ! *dz
!
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
!        QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki+1)*DZ) &
!                 +entr*DZ*q(i,Ki) &
!                )/(1.+entr*DZ-.5*CDD(i,Ki+1)*DZ)
!        dz=qcd(i,ki)
         qcd(i,ki)=(qcd(i,ki+1)*zd(i,ki+1)                          &
                  -.5*dd_massdetr(i,ki)*qcd(i,ki+1)+ &
                  dd_massentr(i,ki)*q(i,ki))   /            &
                  max(epsilc,zd(i,ki+1)-.5*dd_massdetr(i,ki)+dd_massentr(i,ki))
!        write(0,*)'qcd in dd_moi = ',qcd(i,ki)

!
!--- to be negatively buoyant, hcd should be smaller than hes!
!--- ideally, dh should be negative till dd hits ground, but that is not always
!--- the case
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XLV)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)then
          dqeva=0.
          qrcd(i,ki)=qcd(i,ki)
         endif
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki) ! *dz
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem with buoy in cup_dd_moisture"
       endif
       if(BU(I).GE.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem2 with buoy in cup_dd_moisture"
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture_new

   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,                   &
              xlv,r_v,cp,g,itf,ktf,                     &
              its,ite,  kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,ktf,           &
        its,ite,  kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t,q
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        he,z
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real, intent (in   )                   ::                           &
        xlv,r_v,cp,g
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k,iph
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar
!      real, external :: satvap
!      real :: satvap


      HT(1)=XLV/CP
      HT(2)=2.834E6/CP
      BE(1)=.622*HT(1)/.286
      AE(1)=BE(1)/273.+ALOG(610.71)
      BE(2)=.622*HT(2)/.286
      AE(2)=BE(2)/273.+ALOG(610.71)
!      print *, 'TCRIT = ', tcrit,its,ite
      DO k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
        IPH=1
        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
!       E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
!       QES(I,K)=.622*E/(100.*P(I,K)-E)
        e=satvap(t(i,k))
        qes(i,k)=0.622*e/max(1.e-8,(p(i,k)-e))
        IF(QES(I,K).LE.1.E-16)QES(I,K)=1.E-16
        IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
!       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.eq.1 .or. itest.eq.0)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                 ALOG(PSUR(I)))*287.*TV(I,1)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
      else if(itest.eq.2)then
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      else if(itest.eq.-1)then
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.le.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo
 
   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              xlv,r_v,cp,g,ierr,z1,                                &
              itf,ktf,                       &
              its,ite,  kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,ktf,           &
        its,ite,  kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real,    intent (in   )                   ::                           &
        xlv,r_v,cp,g
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k


      do k=kts,ktf
      do i=its,itf
        qes_cup(i,k)=0.
        q_cup(i,k)=0.
        hes_cup(i,k)=0.
        he_cup(i,k)=0.
        z_cup(i,k)=0.
        p_cup(i,k)=0.
        t_cup(i,k)=0.
        gamma_cup(i,k)=0.
      enddo
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xlv/cp)*(xlv/(r_v*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
!       hes_cup(i,1)=hes(i,1)
!       he_cup(i,1)=he(i,1)
        hes_cup(i,1)=9.81*z1(i)+1004.*t(i,1)+2.5e6*qes(i,1)
        he_cup(i,1)=9.81*z1(i)+1004.*t(i,1)+2.5e6*q(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        z_cup(i,1)=z1(i)
        p_cup(i,1)=psur(i)
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xlv/cp*(xlv/(r_v*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup_env_clev


   SUBROUTINE cup_kbcon(ierrc,cap_inc,iloop_in,k22,kbcon,he_cup,hes_cup, &
              hkb,ierr,kbmax,p_cup,cap_max,                         &
              ztexec,zqexec,use_excess,       &
              itf,ktf,                        &
              its,ite, kts,kte, &
              z_cup,entr_rate,heo,imid                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        use_excess,itf,ktf,imid,           &
        its,ite, kts,kte
  ! 
  ! 
  ! 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        entr_rate,ztexec,zqexec,cap_inc,cap_max
     real,    dimension (its:ite)                                      &
        ,intent (inout   )                   ::                           &
        hkb !,cap_max
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop_in
     character*50 :: ierrc(its:ite)
     real, dimension (its:ite,kts:kte),intent (in) :: z_cup,heo
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,k1,k2,iloop
     real                                 ::                           &
        pbcdif,plus,hetest,dz
     real, dimension (its:ite,kts:kte) ::hcot
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      iloop=iloop_in
       DO 27 i=its,itf
      kbcon(i)=1
!
! reset iloop for mid level convection
      if(cap_max(i).gt.200 .and. imid.eq.1)iloop=5
!
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)+1
      if(iloop.eq.5)KBCON(I)=K22(I)
       !== including entrainment for hetest
        hcot(i,1:k22(i)) = HKB(I)
        do k=k22(i)+1,KBMAX(i)+3
           dz=z_cup(i,k)-z_cup(i,k-1)
           hcot(i,k)= ( (1.-0.5*entr_rate(i)*dz)*hcot(i,k-1)   &
                         + entr_rate(i)*dz*heo(i,k-1)       )/ &
                      (1.+0.5*entr_rate(i)*dz)
        enddo
       !==

      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
         endif
        GO TO 27
      ENDIF
 32   CONTINUE
      hetest=hcot(i,kbcon(i)) !hkb(i) ! HE_cup(I,K22(I))
!      if(iloop.eq.5)then
!       hetest=HKB(I)
!      endif
      IF(HETEST.LT.HES_cup(I,KBCON(I)))then
        GO TO 31
      endif

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.1)go to 27
      if(iloop.eq.5 .and. (KBCON(I)-K22(I)).le.2)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
!
! for shallow convection, if cap_max is greater than 25, it is the pressure at pbltop
      if(iloop.eq.5)plus=150.
      if(iloop.eq.5.and.cap_max(i).gt.200)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
      IF(PBCDIF.GT.plus)THEN
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)+1
       !== including entrainment for hetest
        hkb(i)=HE_cup(I,K22(I)) !lzhang
        hcot(i,1:k22(i)) = HKB(I)
        do k=k22(i)+1,KBMAX(i)+3
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)= ( (1.-0.5*entr_rate(i)*dz)*hcot(i,k-1)   &
                         + entr_rate(i)*dz*heo(i,k-1)       )/ &
                      (1.+0.5*entr_rate(i)*dz)
        enddo
       !==

        if(iloop.eq.5)KBCON(I)=K22(I)
        IF(KBCON(I).GT.KBMAX(i)+2)THEN
            if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
            endif
            GO TO 27
        ENDIF
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon


   SUBROUTINE cup_ktop(ierrc,ilo,dby,kbcon,ktop,ierr,              &
              itf,ktf,                     &
              its,ite,  kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,ktf,           &
        its,ite,  kts,kte
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! ilo  = flag
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        ilo
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
!
        DO 42 i=its,itf
        ktop(i)=1
         IF(ierr(I).EQ.0)then
          DO 40 K=KBCON(I)+1,ktf-1
            IF(DBY(I,K).LE.0.)THEN
                KTOP(I)=K-1
                GO TO 41
             ENDIF
  40      CONTINUE
          if(ilo.eq.1)ierr(i)=5
          if(ilo.eq.1)ierrc(i)="problem with defining ktop"
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         if(kbcon(i).eq.ktop(i))then
            ierr(i)=55
            ierrc(i)="kbcon == ktop "
         endif
         endif
  42     CONTINUE

   END SUBROUTINE cup_ktop


   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              itf,ktf,                     &
              its,ite,  kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,ktf,                                    &
         its,ite,  kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         maxx
     real,    dimension (its:ite)         ::                           &
         x
     real                                 ::                           &
         xar
     integer                              ::                           &
         i,k

       DO 200 i=its,itf
       MAXX(I)=KS
       if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS)
!
       DO 100 K=KS,KE(i)
         XAR=ARRAY(I,K)
         IF(XAR.GE.X(I)) THEN
            X(I)=XAR
            MAXX(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MAXIMI


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
              itf,ktf,                     &
              its,ite,  kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,ktf,                                    &
         its,ite,  kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         kt
     real,    dimension (its:ite)         ::                           &
         x
     integer                              ::                           &
         i,k,kstop

       DO 200 i=its,itf
      KT(I)=KS(I)
      if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS(I))
       KSTOP=MAX(KS(I)+1,KEND(I))
!
       DO 100 K=KS(I)+1,KSTOP
         IF(ARRAY(I,K).LT.X(I)) THEN
              X(I)=ARRAY(I,K)
              KT(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MINIMI

!====================================================================
   SUBROUTINE neg_check(ktop,name,j,dt,q,outq,outt,outu,outv,   &
                        outqc,outc,pret,its,ite,kts,kte,itf,ktf,numc)

   INTEGER,      INTENT(IN   ) ::  numc,j,its,ite,kts,kte,itf,ktf

     real, dimension (its:ite,kts:kte  )                    ,                 &
      intent(inout   ) ::                                                     &
       outq,outt,outqc,outu,outv
     real, dimension (its:ite,kts:kte,numc  )                    ,            &
      intent(inout   ) ::                                                     &
        outc
     real, dimension (its:ite,kts:kte  )                    ,                 &
      intent(inout   ) ::                                                     &
       q
     real, dimension (its:ite  )                            ,                 &
      intent(inout   ) ::                                                     &
       pret
      character *(*), intent (in)         ::                           &
       name
     real                                                                     &
        ,intent (in  )                   ::                                   &
        dt
     real :: names,scalef,thresh,qmem,qmemf,qmem2,qtest,qmem1
     integer :: i,k,nv,icheck
     integer, dimension (its:ite  )                            , &
        intent (in  )                   ::                                   &
         ktop
!
! first do check on vertical heating rate
!
      icheck=0
      thresh=300.01
      thresh=200.01	!ss
      thresh=250.01
      names=1.
      if(name == 'shallow')then
        thresh=48.01
        names=2.
      endif
      scalef=86400.
      do i=its,itf
      qmemf=1.
      qmem=0.
      do k=kts,ktop(i)
         qmem=(outt(i,k))*86400.
         if(qmem.gt.thresh)then
           qmem2=thresh/qmem
           qmemf=min(qmemf,qmem2)
      icheck=1
!
!
!          print *,'1',' adjusted massflux by factor ',i,j,k,qmem,qmem2,qmemf,dt
         endif
         if(qmem.lt.-.5*thresh*names)then
           qmem2=-.5*names*thresh/qmem
           qmemf=min(qmemf,qmem2)
      icheck=2
!
!
!!sms$ignore begin
!          write(0,*)name, icheck,' adjusted massflux by factor ',qmem2, 'at j,k = ',j,k
!!sms$ignore end
         endif
      enddo
!    if(name == 'shallow')then
 !    if(j.eq.586981)then
!!sms$ignore begin
!      if(qmemf.lt.0.5)then
!           write(12,*)'1',' adjusted massflux by factor ',j,qmemf,icheck,name
!        do k=kts,ktf-5
!         write(12,*)k,scalef*(outq(i,k)),scalef*qmemf*(outt(i,k)+)
!        enddo
!      endif
!!sms$ignore end
!!     endif
!      endif
      do k=kts,ktop(i)
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outu(i,k)=outu(i,k)*qmemf
         outv(i,k)=outv(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      if(numc.gt.0)then
         do nv=1,numc
            do k=kts,ktop(i)
              outc(i,k,nv)=outc(i,k,nv)*qmemf
            enddo
         enddo
      endif
      pret(i)=pret(i)*qmemf 
      enddo
!
! check whether routine produces negative q's. This can happen, since 
! tendencies are calculated based on forced q's. This should have no
! influence on conservation properties, it scales linear through all
! tendencies
!
      thresh=1.e-10
      do i=its,itf
      qmemf=1.
      do k=kts,ktop(i)
         qmem=outq(i,k)
         if(abs(qmem).gt.0.)then
         qtest=q(i,k)+(outq(i,k))*dt
         if(qtest.lt.thresh)then
!
! qmem2 would be the maximum allowable tendency
!
           qmem1=outq(i,k)
           qmem2=(thresh-q(i,k))/dt
           qmemf=min(qmemf,qmem2/qmem1)
           qmemf=max(0.,qmemf)
!          write(0,*)'4 adjusted tendencies ',i,k,qmem,qmem2,qmemf
!          write(0,*)'4 adjusted tendencies ',i,j,k,q(i,k),qmem1,qmemf
         endif
         endif
      enddo
!     if(qmemf.lt.1.)write(0,*)'4 adjusted tendencies ',i,j,qmemf
      do k=kts,ktop(i)
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outu(i,k)=outu(i,k)*qmemf
         outv(i,k)=outv(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      if(numc.gt.0)then
         do nv=1,numc
            do k=kts,ktop(i)
              outc(i,k,nv)=outc(i,k,nv)*qmemf
            enddo
         enddo
      endif
      pret(i)=pret(i)*qmemf 
      enddo

   END SUBROUTINE neg_check

!-------------------------------------------------------
   SUBROUTINE cup_up_moisture(name,ierr,z_cup,qc,qrc,pw,pwav,     &
              ccnclean,p_cup,kbcon,ktop,cd,dby,clw_all,&
              t_cup,q,GAMMA_cup,zu,qes_cup,k22,qe_cup,         &
              ZQEXEC,use_excess,ccn,rho,c1d, &
              up_massentr,up_massdetr,psum,psumh,                 &
              xlv,c0t,autoconv,aeroevap,itest,itf,ktf,ipr,                &
              its,ite,  kts,kte                     )

   IMPLICIT NONE
  real, parameter :: BDISPM = 0.366       !Berry--size dispersion (martime)
  REAL, PARAMETER :: BDISPC = 0.146       !Berry--size dispersion (continental)
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  use_excess,itest,autoconv,aeroevap,itf,ktf,           &
                                  its,ite, kts,kte
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        t_cup,p_cup,rho,q,zu,gamma_cup,qe_cup,                         &
        up_massentr,up_massdetr,dby,qes_cup,z_cup,cd
     real,    dimension (its:ite)                              &
        ,intent (in   )                   ::                           &
        zqexec
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean,xlv
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22,ipr
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all,c0t
     real,    dimension (its:ite,kts:kte) ::                           &
        qch,qrcb,pwh,clw_allh,c1d
     real,    dimension (its:ite)         ::                           &
        pwavh
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav,psum,psumh
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                           &
        ccn
!
!  local variables in this routine
!

     integer                              ::                           &
        iounit,iprop,iall,i,k,k1,k2
     real                                 ::                           &
        prop_ave,qrcb_h,bdsp,dp,g,rhoc,dh,qrch,dz,radius,berryc0,q1,berryc
     real                                 ::                           &
        denom,c0
     real,    dimension (kts:kte)         ::                           &
        prop_b
!
        prop_b(kts:kte)=0
        c0=.002
        !c0t(:,:)=c0
        iall=0
        g=9.81
        bdsp=BDISPM
!
!--- no precip for small clouds
!
        if(name.eq.'shallow')then
!           iall=1
        endif
        do i=its,itf
          pwav(i)=0.
          pwavh(i)=0.
          psum(i)=0.
          psumh(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          pwh(i,k)=0.
          qc(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qe_cup(i,k)
          if(ierr(i).eq.0)qch(i,k)=qe_cup(i,k)
          clw_all(i,k)=0.
          clw_allh(i,k)=0.
          qrc(i,k)=0.
          qrcb(i,k)=0.
!from Haiqing Li
         if(t_cup(i,k) > 273.16) then
           c0t(i,k) = c0
         else
           c0t(i,k) = c0 * exp(0.07 * (t_cup(i,k) - 273.16))
         endif

        enddo
        enddo
!     if(use_excess < 2 ) then
      do i=its,itf
      if(ierr(i).eq.0.)then
!
!  initialize below originating air
!
      do k=2,k22(i)
        DZ=Z_cup(i,K)-Z_cup(i,K-1)
        qc(i,k)=qe_cup(i,k)!+float(use_excess)*zqexec(i)
        qch(i,k)=qe_cup(i,k)!+float(use_excess)*zqexec(i)
      enddo
      endif
      enddo

       DO 100 i=its,itf
         IF(ierr(i).eq.0)then

! below LFC, but maybe above LCL
!
            DO k=k22(i)+1,kbcon(i)-1
              qc(i,k)=   (qc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
!              QRCH=QES_cup(I,K)
               QRCH=QES_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k) &
                 /(1.+GAMMA_cup(i,k)))*DBY(I,K)
              if(qc(i,k).gt.qrch)then
                DZ=Z_cup(i,K)-Z_cup(i,K-1)
                QRC(I,K)=(QC(I,K)-QRCH)/(1.+(c0t(i,k)+c1d(i,k))*DZ)
                PW(i,k)=c0t(i,k)*dz*QRC(I,K)*zu(i,k)
                qc(i,k)=qrch+qrc(i,k)
              endif
            enddo
!
!now do the rest
!
            DO k=kbcon(i),ktop(i)
               denom=zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1)
               if(denom.lt.1.e-12)then
           write(14,*)k,denom,zu(i,k-1),up_massdetr(i,k-1),up_massentr(i,k-1)
                     ierr(i)=51
                exit
               endif

   
               rhoc=.5*(rho(i,k)+rho(i,k-1))
               DZ=Z_cup(i,K)-Z_cup(i,K-1)
               DP=p_cup(i,K)-p_cup(i,K-1)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
               QRCH=QES_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k) &
                 /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
               qc(i,k)=   (qc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
               qch(i,k)= (qch(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*qch(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

               if(qc(i,k).le.qrch)then
                 qc(i,k)=qrch
               endif
               if(qch(i,k).le.qrch)then
                 qch(i,k)=qrch
               endif
!
!------- Total condensed water before rainout
!
               clw_all(i,k)=max(0.,QC(I,K)-QRCH)
               QRC(I,K)=max(0.,(QC(I,K)-QRCH)) ! /(1.+C0*DZ*zu(i,k))
               clw_allh(i,k)=max(0.,QCH(I,K)-QRCH)
               QRCB(I,K)=max(0.,(QCH(I,K)-QRCH)) ! /(1.+C0*DZ*zu(i,k))
               IF(autoconv.eq.2) then


! 
! normalized berry
!
! first calculate for average conditions, used in cup_dd_edt!
! this will also determine proportionality constant prop_b, which, if applied,
! would give the same results as c0 under these conditions
!
                 q1=1.e3*rhoc*qrcb(i,k)  ! g/m^3 ! g[h2o]/cm^3
                 berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCNclean/ &
                    ( q1 * BDSP)  ) ) !/(
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',k,rhoc,rho(i,k)
!        qrcb_h=qrcb(i,k)/(1.+c0*dz)
                 qrcb_h=((QCH(I,K)-QRCH)*zu(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1)+c0t(i,k)*dz*zu(i,k))
                 prop_b(k)=c0t(i,k)*qrcb_h*zu(i,k)/(1.e-3*berryc0)
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',berryc0,prop_b(k),qrcb_h
                 pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k) ! 2.
                 berryc=qrcb(i,k)
                 qrcb(i,k)=((QCh(I,K)-QRCH)*zu(i,k)-pwh(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                       (zu(i,k)+.5*up_massdetr(i,k-1))
                 if(qrcb(i,k).lt.0.)then
                   berryc0=(qrcb(i,k-1)*(.5*up_massdetr(i,k-1))-(QCh(I,K)-QRCH)*zu(i,k))/zu(i,k)*1.e-3*dz*prop_b(k)
                   pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k)
                   qrcb(i,k)=0.
                 endif
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',zu(i,k),pwh(i,k),dz,qrch,qrcb(i,k),clw_allh(i,k)
                 QCh(I,K)=QRCb(I,K)+qrch
                 PWAVH(I)=PWAVH(I)+pwh(I,K)
                 Psumh(I)=Psumh(I)+clw_allh(I,K)*zu(i,k) *dz
        !
! then the real berry
!
                 q1=1.e3*rhoc*qrc(i,k)  ! g/m^3 ! g[h2o]/cm^3
                 berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCN(i)/ &
                    ( q1 * BDSP)  ) ) !/(
                 berryc0=1.e-3*berryc0*dz*prop_b(k) ! 2.
                 berryc=qrc(i,k)
                 qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-zu(i,k)*berryc0-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                       (zu(i,k)+.5*up_massdetr(i,k-1))
                 if(qrc(i,k).lt.0.)then
                    berryc0=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/zu(i,k)
                    qrc(i,k)=0.
                 endif
                 pw(i,k)=berryc0*zu(i,k)
                 QC(I,K)=QRC(I,K)+qrch
!
!  if not running with berry at all, do the following
!
               ELSE       !c0=.002
                 if(iall.eq.1)then
                   qrc(i,k)=0.
                   pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
                   if(pw(i,k).lt.0.)pw(i,k)=0.
                 else
                   QRC(I,K)=(QC(I,K)-QRCH)/(1.+(c1d(i,k)+c0t(i,k))*DZ)
                   PW(i,k)=c0t(i,k)*dz*QRC(I,K)*zu(i,k)
          !        qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
!                  (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
!        PW(i,k)=c0*dz*qrc(i,k)*zu(i,k)
                   if(qrc(i,k).lt.0)then
                     qrc(i,k)=0.
                     pw(i,k)=0.
                   endif
                 endif
                 QC(I,K)=QRC(I,K)+qrch
               endif !autoconv
               PWAV(I)=PWAV(I)+PW(I,K)
               Psum(I)=Psum(I)+clw_all(I,K)*zu(i,k) *dz
            enddo ! k=kbcon,ktop
      endif ! ierr
!
!--- integrated normalized ondensate
!
 100     CONTINUE
       prop_ave=0.
       iprop=0
       do k=kts,kte
        prop_ave=prop_ave+prop_b(k)
        if(prop_b(k).gt.0)iprop=iprop+1
       enddo
       iprop=max(iprop,1)
!      write(11,*)'prop_ave = ',prop_ave/float(iprop)
!      print *,'pwav = ',pwav(1)

   END SUBROUTINE cup_up_moisture
!====================================================================

!--------------------------------------------------------------------

      real function satvap(temp2)
      implicit none
      real :: temp2, temp, toot, toto, eilog, tsot,  &
     &        ewlog, ewlog2, ewlog3, ewlog4
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
        toot = 273.16 / temp2
        toto = 1 / toot
        eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
     &    log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
        satvap = 10 ** eilog
      else
        tsot = 373.16 / temp2
        ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
     &             (log(tsot) / log(10.))
        ewlog2 = ewlog - 1.3816e-07 * &
     &             (10 ** (11.344 * (1 - (1 / tsot))) - 1)
        ewlog3 = ewlog2 + .0081328 * &
     &             (10 ** (-3.49149 * (tsot - 1)) - 1)
        ewlog4 = ewlog3 + (log(1013.246) / log(10.))
        satvap = 10 ** ewlog4
      end if
      return
      end function
 SUBROUTINE rates_up_pdf(name,ktop,ierr,p_cup,entr_rate_2d,hkbo,heo,heso_cup,z_cup, &
               kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zu,kpbli,ktopdby,csum)
     implicit none
     integer, intent(in) :: its,ite,itf,kts,kte,ktf
     real, dimension (its:ite,kts:kte),intent (inout) :: entr_rate_2d,zu
     real, dimension (its:ite,kts:kte),intent (in) ::p_cup, heo,heso_cup,z_cup
     real, dimension (its:ite),intent (in) :: hkbo
     integer, dimension (its:ite),intent (in) :: kstabi,k22,kbcon,kpbli,csum
     integer, dimension (its:ite),intent (inout) :: ierr,ktop,ktopdby
     real, dimension (its:ite,kts:kte) :: hcot
      character *(*), intent (in)         ::                           &
       name
     real :: dz,dh, &
             dbythresh,dzh2
     real :: dby(kts:kte)
      real zh2(40)
     integer :: kk,kbegin,i,k,ipr,kdefi,kstart,kbegzu,kfinalzu
     !integer :: lev_adj,power_entr
     dbythresh=0.6
     dbythresh=0.85
     dbythresh=0.95
     dby(:)=0.
!     dbythresh=1.
!
!
     DO i=its,itf
      zu(i,:)=0.
      dby(:)=0.
      ktop(i)=0
      if(ierr(i).eq.0)then
        hcot(i,kbcon(i))=hkbo(i)
        dz=z_cup(i,kbcon(i))-z_cup(i,kbcon(i)-1)
        dby(kbcon(i))=(hcot(i,kbcon(i))-heso_cup(i,kbcon(i)))*dz
        do k=kbcon(i)+1,ktf-2
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)=( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1) &
                      + entr_rate_2d(i,k-1)*dz*heo(i,k-1))/ &
                      (1.+0.5*entr_rate_2d(i,k-1)*dz)
           dby(k)=dby(k-1)+(hcot(i,k)-heso_cup(i,k))*dz
       enddo
       ktopdby(i)=maxloc(dby(:),1)
       do k=maxloc(dby(:),1)+1,ktf-2
          if(dby(k).lt.dbythresh*maxval(dby))then
              kfinalzu=k  - 1
              ktop(i)=kfinalzu
              go to 412
          endif
       enddo
       kfinalzu=ktf-2
       ktop(i)=kfinalzu
412    continue
!
! at least overshoot by one level
!
       kfinalzu=min(max(kfinalzu,ktopdby(i)+1),ktopdby(i)+2)
       ktop(i)=kfinalzu
       if(kfinalzu.le.kbcon(i)+2)then
           ierr(i)=41
           ktop(i)= 0
       else
           call get_zu_zd_pdf_fim("UP",ierr(i),k22(i),kfinalzu,zu(i,kts:kte),kts,kte,ktf,0.,kpbli(i),csum(i))
       endif
      ENDIF
     ENDDO

  END SUBROUTINE rates_up_pdf
!-------------------------------------------------------------------------
subroutine get_zu_zd_pdf_fim(draft,ierr,kb,kt,zu,kts,kte,ktf,max_mass,kpbli,csum)

implicit none
integer, intent(in) ::kb,kt,kts,kte,ktf,kpbli,csum
real, intent(in) ::max_mass
real, intent(inout) :: zu(kts:kte)
real  :: zuh(kts:kte)
integer, intent(inout) :: ierr
character*(*), intent(in) ::draft

!- local var
integer :: kk,add,i,nrec=0,k,kb_adj,kpbli_adj
real ::krmax,zumax,ztop_adj
real ::beta, alpha,kratio,tunning,FZU
!- kb cannot be at 1st level

   !-- fill zu with zeros
   zu=0.0
   zuh=0.0

IF(draft == "UP") then
   kb_adj=max(kb,2)
  !beta=4.  !=> must larger than 1
            !=> higher makes the profile sharper
            !=> around the maximum zu 
  !- 2nd approach for beta and alpha parameters
  !- the tunning parameter must be between 0.5 (low  level max zu)
  !-                                   and 1.5 (high level max zu)
  tunning = 1.2-min(0.6,float(csum)*.05)
  beta    = 2. ! 2.0/tunning
  alpha   = tunning*beta
  fzu=1.
  do k=kts,min(kte,kt+1)
      kratio= float(k)/float(kt+1)
      zu(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
   enddo
   if(maxval(zu(kts:min(ktf,kt+1))).gt.0.)  &
      zu(kts:min(ktf,kt+1))= zu(kts:min(ktf,kt+1))/maxval(zu(kts:min(ktf,kt+1)))
     do k=maxloc(zu(:),1),1,-1
       if(zu(k).lt.1.e-6)then
         kb_adj=k+1
         exit
       endif
     enddo
     kb_adj=max(2,kb_adj)
     do k=kts,kb_adj-1
       zu(k)=0.
     enddo

ELSEIF(draft == "DOWN" .or. draft == "DOWNM") then

  tunning = 0.8
  beta    =2.0/tunning
  alpha   = tunning*beta
  fzu=1.
  zuh(:)=0.
  zu(:)=0.
  do k=kts,min(kt+1,ktf)
      kratio= float(k)/float(kt+1)
      zuh(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
!      write(0,*) k,zuh(k)
   enddo
   if(maxloc(zuh(:),1).ge.kpbli)then
      do k=maxloc(zuh(:),1),1,-1
         kk=kpbli+k-maxloc(zuh(:),1)
         if(kk.gt.1)zu(kk)=zuh(k)
!         if(kk.gt.1)write(0,*)kk,zu(kk)
      enddo
      do k=maxloc(zuh(:),1)+1,kt
         kk=kpbli+k-maxloc(zuh(:),1)
         if(kk.le.kt)zu(kk)=zuh(k)
!         if(kk.ge.1)write(0,*)kk,zu(kk)
      enddo
   else
      do k=2,kt ! maxloc(zuh(:),1)
        zu(k)=zuh(k-1)
!         write(0,*)k,zu(k)
      enddo
   endif
   fzu=maxval(zu(kts:min(ktf,kt+1)))
   if(fzu.gt.0.)  &
      zu(kts:min(ktf,kt+1))= zu(kts:min(ktf,kt+1))/fzu
    if(zu(2).gt.max_mass)fzu=max_mass/zu(2) ! max(0.,zu(2)-max_mass)
     do k=2,kt+1
       zu(k)=fzu*zu(k)
     enddo
     zu(1)=0.

ENDIF
   !- normalize ZU
!  if(maxval(zu(kts:min(ktf,kt+1))).gt.0.)  &
!     zu(kts:min(ktf,kt+1))= zu(kts:min(ktf,kt+1))/ maxval(zu(kts:min(ktf,kt+1)))



return


end subroutine get_zu_zd_pdf_fim



!---------------------------------------------------------------------- 



   SUBROUTINE neg_check_ct(pret,ktop,epsilc,dt,q,outq,iopt,num_chem,    &
                           its,ite,kts,kte,itf,ktf,ipr,jpr,npr,j)

   INTEGER,      INTENT(IN   ) ::   iopt,num_chem,its,ite,kts,kte,itf,ktf,ipr,jpr,npr,j

     real, dimension (its:ite,kts:kte,num_chem  )          ,                  &
      intent(inout   ) ::                                                     &
       q,outq
     real, dimension (its:ite  )          ,                                   &
      intent(in      ) ::                                                     &
       pret
     integer, dimension (its:ite  )          ,                                &
      intent(in   ) ::                                                        &
      ktop
     real                                                                     &
        ,intent (in  )                   ::                                   &
        dt,epsilc
     real :: tracermin,tracermax,thresh,qmem,qmemf,qmem2,qtest,qmem1

     integer :: i,k,nv
! check whether routine produces negative q's. This can happen, since 
! tendencies are calculated based on forced q's. This should have no
! influence on conservation properties, it scales linear through all
! tendencies. Use iopt=0 to test for each tracer seperately, iopt=1
! for a more severe limitation...
!
      thresh=epsilc
!     thresh=1.e-30
      if(iopt.eq.0)then
      do nv=1,num_chem
      do 100 i=its,itf
         if(pret(i).le.0.)go to 100
         tracermin=q(i,kts,nv)
         tracermax=q(i,kts,nv)
         do k=kts+1,ktop(i)
           tracermin=min(tracermin,q(i,k,nv))
           tracermax=max(tracermax,q(i,k,nv))
         enddo
         tracermin=max(tracermin,thresh)
         qmemf=1.
!
! first check for minimum restriction
!
         do k=kts,ktop(i)
!
! tracer tendency
!
            qmem=outq(i,k,nv)
!
! only necessary if there is a tendency
!
            if(qmem.lt.0.)then
               qtest=q(i,k,nv)+outq(i,k,nv)*dt
               if(qtest.lt.tracermin)then
!
! qmem2 would be the maximum allowable tendency
!
                    qmem1=outq(i,k,nv)
                    qmem2=(tracermin-q(i,k,nv))/dt
                    qmemf=min(qmemf,qmem2/qmem1)
                    if(qmemf.gt.1.)print *,'something wrong in negct_1',qmem2,qmem1
                    if(i.eq.ipr.and.j.eq.jpr.and.nv.eq.npr)then
                      print *,k,qtest,qmem2,qmem1,qmemf
                    endif
                    qmemf=max(qmemf,0.)
               endif
            endif
         enddo
         do k=kts,ktop(i)
            outq(i,k,nv)=outq(i,k,nv)*qmemf
         enddo
!
! now check max
!
         qmemf=1.
         do k=kts,ktop(i)
!
! tracer tendency
!
            qmem=outq(i,k,nv)
!
! only necessary if there is a tendency
!
            if(qmem.gt.0.)then
               qtest=q(i,k,nv)+outq(i,k,nv)*dt
               if(qtest.gt.tracermax)then
!
! qmem2 would be the maximum allowable tendency
!
                    qmem1=outq(i,k,nv)
                    qmem2=(tracermax-q(i,k,nv))/dt
                    qmemf=min(qmemf,qmem2/qmem1)
                    if(qmemf.gt.1.)print *,'something wrong in negct_2',qmem2,qmem1
                    if(i.eq.ipr.and.j.eq.jpr.and.nv.eq.npr)then
                      print *,'2',k,qtest,qmem2,qmem1,qmemf
                    endif
                    qmemf=max(qmemf,0.)
               endif
            endif
         enddo
         do k=kts,ktop(i)
            outq(i,k,nv)=outq(i,k,nv)*qmemf
         enddo
 100  continue
      enddo
!
! ELSE
!
      elseif(iopt.eq.1)then
      do i=its,itf
      qmemf=1.
      do k=kts,ktop(i)
      do nv=1,num_chem
!
! tracer tendency
!
         qmem=outq(i,k,nv)
!
! only necessary if tendency is larger than zero
!
         if(qmem.lt.0.)then
         qtest=q(i,k,nv)+outq(i,k,nv)*dt
         if(qtest.lt.thresh)then
!
! qmem2 would be the maximum allowable tendency
!
           qmem1=outq(i,k,nv)
           qmem2=(thresh-q(i,k,nv))/dt
           qmemf=min(qmemf,qmem2/qmem1)
           qmemf=max(0.,qmemf)
         endif
         endif
      enddo
      enddo
      do nv=1,num_chem
      do k=kts,ktop(i)
         outq(i,k,nv)=outq(i,k,nv)*qmemf
      enddo
      enddo
      enddo
      endif

   END SUBROUTINE neg_check_ct


!-------------------------------------------------------
END MODULE dep_ctrans_grell_mod
