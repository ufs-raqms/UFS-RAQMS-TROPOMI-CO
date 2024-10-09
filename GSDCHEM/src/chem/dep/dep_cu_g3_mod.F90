MODULE dep_cu_g3_mod


CONTAINS
!-------------------------------------------------------------


   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pwev,edtmax,edtmin,maxens2,edtc,               &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        maxens2
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        us,vs,z,p
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev
     real                                                              &
        ,intent (in   )                   ::                           &
        edtmax,edtmin
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
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
          if (kk .eq. ktf)vshear(i) = 1.e3 * vws(i) / sdp(i)
   62   continue
       end do
      do i=its,itf
         IF(ierr(i).eq.0)then
            pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
               -.00496*(VSHEAR(I)**3))
            if(pef.gt.1.)pef=1.
            if(pef.lt.0.)pef=0.
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
            if(pefb.gt.1.)pefb=1.
            if(pefb.lt.0.)pefb=0.
            EDT(I)=1.-.5*(pefb+pef)
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
               EDTC(I,K)=-EDTC(I,K)*PWAV(I)/PWEV(I)
               IF(EDTC(I,K).GT.edtmax)EDTC(I,K)=edtmax
               IF(EDTC(I,K).LT.edtmin)EDTC(I,K)=edtmin
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,entr,       &
              jmin,ierr,he,dby,he_cup,                       &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! hcd = downdraft moist static energy
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cdd= detrainment function 
  ! z_cup = heights of model cloud levels 
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cdd,zd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hcd,dby
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dz


      do k=kts+1,ktf
      do i=its,itf
      dby(i,k)=0.
      IF(ierr(I).eq.0)then
         hcd(i,k)=hes_cup(i,k)
      endif
      enddo
      enddo
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      hcd(i,k)=hes_cup(i,k)
      dby(i,k)=hcd(i,jmin(i))-hes_cup(i,k)
!
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         HCD(i,Ki)=(HCD(i,Ki+1)*(1.-.5*CDD(i,Ki)*DZ) &
                  +entr*DZ*HE(i,Ki) &
                  )/(1.+entr*DZ-.5*CDD(i,Ki)*DZ)
         dby(i,ki)=HCD(i,Ki)-hes_cup(i,ki)
      enddo
!
      endif
!--- end loop over i
100    continue


   END SUBROUTINE cup_dd_he


   SUBROUTINE cup_dd_moisture_3d(zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,cdd,entr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,xl,high_resolution,           &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte,high_resolution
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
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,cdd,gamma_cup,q,he 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr,xl
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
      if(high_resolution.eq.1)qcd(i,k)=.5*(qes_cup(i,k)+q_cup(i,k))
      qrcd(i,k)=qes_cup(i,k)
      pwd(i,jmin(i))=min(0.,qcd(i,k)-qrcd(i,k))
      pwev(i)=pwev(i)+pwd(i,jmin(i))
      qcd(i,k)=qes_cup(i,k)
!
      DH=HCD(I,k)-HES_cup(I,K)
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki)*DZ) &
                  +entr*DZ*q(i,Ki) &
                  )/(1.+entr*DZ-.5*CDD(i,Ki)*DZ)
!
!--- to be negatively buoyant, hcd should be smaller than hes!
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XL)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)dqeva=0.
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki)
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
       endif
       if(BU(I).GE.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture_3d


   SUBROUTINE cup_dd_nms(zd,z_cup,cdd,entr,jmin,ierr,        &
              itest,kdet,z1,                                 &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! z_cup = height of cloud model level
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! itest = flag to whether to calculate cdd
  
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        z1 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr 
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin,kdet
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
                                                                 ierr
   ! zd is the normalized downdraft mass flux
   ! cdd is the downdraft detrainmen function

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
                                                             zd
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
                                                             cdd
!
!  local variables in this routine
!

     integer                              ::                           &
                                                  i,k,ki
     real                                 ::                           &
                                            a,perc,dz

!
!--- perc is the percentage of mass left when hitting the ground
!
      perc=.03

      do k=kts,ktf
      do i=its,itf
         zd(i,k)=0.
         if(itest.eq.0)cdd(i,k)=0.
      enddo
      enddo
      a=1.-perc
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      zd(i,jmin(i))=1.
!
!--- integrate downward, specify detrainment(cdd)!
!
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         if(ki.le.kdet(i).and.itest.eq.0)then
           cdd(i,ki)=entr+(1.- (a*(z_cup(i,ki)-z1(i)) &
                     +perc*(z_cup(i,kdet(i))-z1(i)) ) &
                         /(a*(z_cup(i,ki+1)-z1(i)) &
                      +perc*(z_cup(i,kdet(i))-z1(i))))/dz
         endif
         zd(i,ki)=zd(i,ki+1)*(1.+(entr-cdd(i,ki))*dz)
      enddo
!
      endif
!--- end loop over i
100    continue

   END SUBROUTINE cup_dd_nms



   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,xl,cp,                   &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
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
        p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z,q
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,cp
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


      HT(1)=XL/CP
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
        E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
        QES(I,K)=.622*E/(100.*P(I,K)-E)
        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.ne.2)then
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
      else
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.eq.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              itf,jtf,ktf,                       &
              its,ite, jts,jte, kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
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
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,rv,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k


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
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
        hes_cup(i,1)=hes(i,1)
        he_cup(i,1)=he(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup_env_clev


   SUBROUTINE cup_forcing_ens_3d(closure_n,xland,aa0,aa1,xaa0,mbdt,dtime,ierr,ierr2,ierr3,&
              xf_ens,j,name,axx,maxens,iens,iedt,maxens2,maxens3,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,massflx,      &
              iact_old_gr,dir,ensdim,massfln,icoic,edt_out,            &
              high_resolution,itf,jtf,ktf,               &
              its,ite, jts,jte, kts,kte,ens4,ktau                )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte,ens4,high_resolution,ktau
     integer, intent (in   )              ::                           &
        j,ensdim,maxens,iens,iedt,maxens2,maxens3
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens,massfln
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (inout   )                   ::                           &
        edt_out
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        massflx
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,zu,p_cup
     real,    dimension (its:ite,kts:kte,1:ens4)                              &
        ,intent (in   )                   ::                           &
        omeg
     real,    dimension (its:ite,1:maxens)                             &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt,dir,xland
     real,    dimension (its:ite,1:ens4)                                      &
        ,intent (in   )                   ::                           &
        mconv,axx
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        aa0,closure_n
     real,    dimension (1:maxens)                                     &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        iact_old_gr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        icoic
      character *(*), intent (in)         ::                           &
       name
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3,iresult,iresultd,iresulte,mkxcrt,kclim
     parameter (mkxcrt=15)
     real                                 ::                           &
       fens4,a1,massfld,a_ave,xff0,xff00,xxx,xomg,aclim1,aclim2,aclim3,aclim4
     real,    dimension(1:mkxcrt)         ::                           &
       pcrit,acrit,acritt

     integer :: nall2,ixxx,irandom
     integer,  dimension (12) :: seed


      DATA PCRIT/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,    &
                 350.,300.,250.,200.,150./
      DATA ACRIT/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,       &
                 .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!  GDAS DERIVED ACRIT
      DATA ACRITT/.203,.515,.521,.566,.625,.665,.659,.688,             &
                  .743,.813,.886,.947,1.138,1.377,1.896/
!
       seed=0
       seed(2)=j
       seed(3)=ktau
       nens=0
       irandom=1
       if(high_resolution.eq.1)irandom=0
       irandom=0
       fens4=float(ens4)

!--- LARGE SCALE FORCING
!
       DO 100 i=its,itf
          if(name.eq.'deeps'.and.ierr(i).gt.995)then
           aa0(i)=0.
           ierr(i)=0
          endif
          IF(ierr(i).eq.0)then
!
!---
!
             if(name.eq.'deeps')then
!
                a_ave=0.
                do ne=1,ens4
                  a_ave=a_ave+axx(i,ne)
                enddo
                a_ave=max(0.,a_ave/fens4)
                a_ave=min(a_ave,aa1(i))
                a_ave=max(0.,a_ave)
                do ne=1,16
                  xff_ens3(ne)=0.
                enddo
                xff0= (AA1(I)-AA0(I))/DTIME
                if(high_resolution.eq.1)xff0= (a_ave-AA0(I))/DTIME
                xff_ens3(1)=(AA1(I)-AA0(I))/dtime
                xff_ens3(2)=(a_ave-AA0(I))/dtime
                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(3)=(axx(i,ixxx)-AA0(I))/dtime
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(13)=(axx(i,ixxx)-AA0(I))/dtime
                else
                   xff_ens3(3)=(AA1(I)-AA0(I))/dtime
                   xff_ens3(13)=(AA1(I)-AA0(I))/dtime
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(1)=(a_ave-AA0(I))/dtime
                   xff_ens3(2)=(a_ave-AA0(I))/dtime
                   xff_ens3(3)=(a_ave-AA0(I))/dtime
                   xff_ens3(13)=(a_ave-AA0(I))/dtime
                endif
!   
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
!
!--- omeg is in bar/s, mconv done with omeg in Pa/s
!     more like Brown (1979), or Frank-Cohen (199?)
!
                xff_ens3(14)=0.
                do ne=1,ens4
                  xff_ens3(14)=xff_ens3(14)-omeg(i,k22(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(14).lt.0.)xff_ens3(14)=0.
                xff_ens3(5)=0.
                do ne=1,ens4
                  xff_ens3(5)=xff_ens3(5)-omeg(i,kbcon(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(5).lt.0.)xff_ens3(5)=0.
!  
! minimum below kbcon
!
                if(high_resolution.eq.0)then
                   xff_ens3(4)=-omeg(i,2,1)/9.81
                   do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.lt.xff_ens3(4))xff_ens3(4)=xomg
                   enddo
                   enddo
                   if(xff_ens3(4).lt.0.)xff_ens3(4)=0.
!
! max below kbcon
                   xff_ens3(6)=-omeg(i,2,1)/9.81
                   do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.gt.xff_ens3(6))xff_ens3(6)=xomg
                   enddo
                   enddo
                   if(xff_ens3(6).lt.0.)xff_ens3(6)=0.
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(5)=min(xff_ens3(5),xff_ens3(14))
                   xff_ens3(4)=xff_ens3(5)
                   xff_ens3(6)=xff_ens3(5)
                endif
!
!--- more like Krishnamurti et al.; pick max and average values
!
                xff_ens3(7)=mconv(i,1)
                xff_ens3(8)=mconv(i,1)
                xff_ens3(9)=mconv(i,1)
                if(ens4.gt.1)then
                   do ne=2,ens4
                      if (mconv(i,ne).gt.xff_ens3(7))xff_ens3(7)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      if (mconv(i,ne).lt.xff_ens3(8))xff_ens3(8)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      xff_ens3(9)=xff_ens3(9)+mconv(i,ne)
                   enddo
                   xff_ens3(9)=xff_ens3(9)/fens4
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(7)=xff_ens3(9)
                   xff_ens3(8)=xff_ens3(9)
                   xff_ens3(15)=xff_ens3(9)
                endif
!
                if(high_resolution.eq.0)then
                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(15)=mconv(i,ixxx)
                else
                   xff_ens3(15)=mconv(i,1)
                endif
                endif
!
!--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
!
                xff_ens3(10)=A_AVE/(60.*40.)
                xff_ens3(11)=AA1(I)/(60.*40.)
                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(12)=AXX(I,ixxx)/(60.*40.)
                else
                   xff_ens3(12)=AA1(I)/(60.*40.)
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(11)=xff_ens3(10)
                   xff_ens3(12)=xff_ens3(10)
                endif
!  
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
!               edt_out(i,j)=xff0
                if(icoic.eq.0)then
                if(xff0.lt.0.)then
                     xff_ens3(1)=0.
                     xff_ens3(2)=0.
                     xff_ens3(3)=0.
                     xff_ens3(13)=0.
                     xff_ens3(10)=0.
                     xff_ens3(11)=0.
                     xff_ens3(12)=0.
                endif
                endif



                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(2)
                   if(xk(nens).le.0.and.xk(nens).gt.-1.e-6) &
                           xk(nens)=-1.e-6
                   if(xk(nens).gt.0.and.xk(nens).lt.1.e-6) &
                           xk(nens)=1.e-6
                enddo
!
!--- add up all ensembles
!
                do 350 ne=1,maxens
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
                   iresult=0
                   iresultd=0
                   iresulte=0
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(ne-1)*maxens3
!
! over water, enfor!e small cap for some of the closures
!
                if(xland(i).lt.0.1)then
                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
                      xff_ens3(1) =0.
                      massfln(i,j,nall+1)=0.
                      xff_ens3(2) =0.
                      massfln(i,j,nall+2)=0.
                      xff_ens3(3) =0.
                      massfln(i,j,nall+3)=0.
                      xff_ens3(10) =0.
                      massfln(i,j,nall+10)=0.
                      xff_ens3(11) =0.
                      massfln(i,j,nall+11)=0.
                      xff_ens3(12) =0.
                      massfln(i,j,nall+12)=0.
                      xff_ens3(7) =0.
                      massfln(i,j,nall+7)=0.
                      xff_ens3(8) =0.
                      massfln(i,j,nall+8)=0.
                      xff_ens3(9) =0.
                      massfln(i,j,nall+9)=0.
                      closure_n(i)=closure_n(i)-1.
                      xff_ens3(13) =0.
                      massfln(i,j,nall+13)=0.
                      xff_ens3(15) =0.
                      massfln(i,j,nall+15)=0.
                endif
                endif
!
! end water treatment
!
!
!--- check for upwind convection
!                  iresult=0
                   massfld=0.

!                  call cup_direction2(i,j,dir,iact_old_gr, &
!                       massflx,iresult,1,                  &
!                       massfld,                            &
!                       itf,jtf,ktf,          &
!                       ims,ime, jms,jme, kms,kme,          &
!                       its,ite, jts,jte, kts,kte          )
!                  if(i.eq.ipr.and.j.eq.jpr.and.iedt.eq.1.and.ne.eq.1)then
!                  if(iedt.eq.1.and.ne.eq.1)then
!                   print *,massfld,ne,iedt,iens
!                   print *,xk(ne),xff_ens3(1),xff_ens3(2),xff_ens3(3)
!                  endif
!                  print *,i,j,massfld,aa0(i),aa1(i)
                   IF(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                   iresulte=max(iresult,iresultd)
                   iresulte=1
                   if(iresulte.eq.1)then
!
!--- special treatment for stability closures
!

                      if(xff0.ge.0.)then
                         xf_ens(i,j,nall+1)=massfld
                         xf_ens(i,j,nall+2)=massfld
                         xf_ens(i,j,nall+3)=massfld
                         xf_ens(i,j,nall+13)=massfld
                         if(xff_ens3(1).gt.0)xf_ens(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne)) &
                                        +massfld
                         if(xff_ens3(2).gt.0)xf_ens(i,j,nall+2)=max(0.,-xff_ens3(2)/xk(ne)) &
                                        +massfld
                         if(xff_ens3(3).gt.0)xf_ens(i,j,nall+3)=max(0.,-xff_ens3(3)/xk(ne)) &
                                        +massfld
                         if(xff_ens3(13).gt.0)xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne)) &
                                        +massfld
!                       endif
                      else
                         xf_ens(i,j,nall+1)=massfld
                         xf_ens(i,j,nall+2)=massfld
                         xf_ens(i,j,nall+3)=massfld
                         xf_ens(i,j,nall+13)=massfld
                      endif
!
!--- if iresult.eq.1, following independent of xff0
!
                         xf_ens(i,j,nall+4)=max(0.,xff_ens3(4) &
                            +massfld)
                         xf_ens(i,j,nall+5)=max(0.,xff_ens3(5) &
                                        +massfld)
                         xf_ens(i,j,nall+6)=max(0.,xff_ens3(6) &
                                        +massfld)
                         xf_ens(i,j,nall+14)=max(0.,xff_ens3(14) &
                                        +massfld)
                         a1=max(1.e-3,pr_ens(i,j,nall+7))
                         xf_ens(i,j,nall+7)=max(0.,xff_ens3(7) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+8)=max(0.,xff_ens3(8) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+9))
                         xf_ens(i,j,nall+9)=max(0.,xff_ens3(9) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+15))
                         xf_ens(i,j,nall+15)=max(0.,xff_ens3(15) &
                                     /a1)
                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+10)=max(0., &
                                        -xff_ens3(10)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+11)=max(0., &
                                        -xff_ens3(11)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+12)=max(0., &
                                        -xff_ens3(12)/xk(ne)) &
                                        +massfld
                         else
                            xf_ens(i,j,nall+10)=massfld
                            xf_ens(i,j,nall+11)=massfld
                            xf_ens(i,j,nall+12)=massfld
                         endif
                      if(icoic.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,j,nall+1)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+2)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+3)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+4)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+5)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+6)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+7)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+8)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+9)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+13)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+14)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+15)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall+icoic)
                      endif
!
! 16 is a randon pick from the oher 15
!
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(15,max(1,int(15.*xxx+1.e-8)))
                   xf_ens(i,j,nall+16)=xf_ens(i,j,nall+ixxx)
                else
                   xf_ens(i,j,nall+16)=xf_ens(i,j,nall+1)
                endif
!
!
!--- store new for next time step
!
                      do nens3=1,maxens3
                        massfln(i,j,nall+nens3)=edt(i) &
                                                *xf_ens(i,j,nall+nens3)
                        massfln(i,j,nall+nens3)=max(0., &
                                              massfln(i,j,nall+nens3))
                      enddo
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!
!
                if(ne.eq.2.and.ierr2(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                      massfln(i,j,nall+1)=0.
                      massfln(i,j,nall+2)=0.
                      massfln(i,j,nall+3)=0.
                      massfln(i,j,nall+4)=0.
                      massfln(i,j,nall+5)=0.
                      massfln(i,j,nall+6)=0.
                      massfln(i,j,nall+7)=0.
                      massfln(i,j,nall+8)=0.
                      massfln(i,j,nall+9)=0.
                      massfln(i,j,nall+10)=0.
                      massfln(i,j,nall+11)=0.
                      massfln(i,j,nall+12)=0.
                      massfln(i,j,nall+13)=0.
                      massfln(i,j,nall+14)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                endif
                if(ne.eq.3.and.ierr3(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                      massfln(i,j,nall+1)=0.
                      massfln(i,j,nall+2)=0.
                      massfln(i,j,nall+3)=0.
                      massfln(i,j,nall+4)=0.
                      massfln(i,j,nall+5)=0.
                      massfln(i,j,nall+6)=0.
                      massfln(i,j,nall+7)=0.
                      massfln(i,j,nall+8)=0.
                      massfln(i,j,nall+9)=0.
                      massfln(i,j,nall+10)=0.
                      massfln(i,j,nall+11)=0.
                      massfln(i,j,nall+12)=0.
                      massfln(i,j,nall+13)=0.
                      massfln(i,j,nall+14)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                endif

                   endif
 350            continue
! ne=1, cap=175
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3
! ne=2, cap=100
!
                   nall2=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(2-1)*maxens3
                      xf_ens(i,j,nall+4) = xf_ens(i,j,nall2+4)
                      xf_ens(i,j,nall+5) =xf_ens(i,j,nall2+5)
                      xf_ens(i,j,nall+6) =xf_ens(i,j,nall2+6)
                      xf_ens(i,j,nall+14) =xf_ens(i,j,nall2+14)
                      xf_ens(i,j,nall+7) =xf_ens(i,j,nall2+7)
                      xf_ens(i,j,nall+8) =xf_ens(i,j,nall2+8)
                      xf_ens(i,j,nall+9) =xf_ens(i,j,nall2+9)
                      xf_ens(i,j,nall+15) =xf_ens(i,j,nall2+15)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall2+10)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall2+11)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall2+12)
                go to 100
             endif
          elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
             do n=1,ensdim
               xf_ens(i,j,n)=0.
               massfln(i,j,n)=0.
             enddo
          endif
 100   continue

   END SUBROUTINE cup_forcing_ens_3d


   SUBROUTINE cup_kbcon(cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              ierr,kbmax,p_cup,cap_max,                         &
              itf,jtf,ktf,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
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
        cap_max,cap_inc
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
!
!  local variables in this routine
!

     integer                              ::                           &
        i
     real                                 ::                           &
        pbcdif,plus
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
       DO 27 i=its,itf
      kbcon(i)=1
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)
      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.lt.4)ierr(i)=3
!        if(iloop.lt.4)ierr(i)=997
        GO TO 27
      ENDIF
 32   CONTINUE
      IF(HE_cup(I,K22(I)).LT.HES_cup(I,KBCON(I)))GO TO 31

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.1)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
      IF(PBCDIF.GT.plus)THEN
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon


   SUBROUTINE cup_ktop(ilo,dby,kbcon,ktop,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
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
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         endif
  42     CONTINUE

   END SUBROUTINE cup_ktop


   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
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
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
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


   SUBROUTINE cup_output_ens_3d(xf_ens,ierr,dellat,dellaq,dellaqc,  &
              subt_ens,subq_ens,subt,subq,outtem,outq,outqc,     &
              zu,sub_mas,pre,pw,xmb,ktop,                 &
              j,name,nx,nx2,iens,ierr2,ierr3,pr_ens,             &
              maxens3,ensdim,massfln,                            &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                 &
              APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,    &
              itf,jtf,ktf, &
              its,ite, jts,jte, kts,kte)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ensdim,nx,nx2,iens,maxens3
  ! xf_ens = ensemble mass fluxes
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
       xf_ens,pr_ens,massfln
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (inout)                   ::                           &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,            &
               APR_CAPME,APR_CAPMI 

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc,subt,subq,sub_mas
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        zu
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        closure_n,xland1
     real,    dimension (its:ite,kts:kte,1:nx)                     &
        ,intent (in   )                   ::                           &
       subt_ens,subq_ens,dellat,dellaqc,dellaq,pw
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,n,ncount
     real                                 ::                           &
        outtes,ddtes,dtt,dtq,dtqc,dtpw,tuning,prerate,clos_wei,xmbhelp
     real                                 ::                           &
        dtts,dtqs
     real,    dimension (its:ite)         ::                           &
       xfac1,xfac2
     real,    dimension (its:ite)::                           &
       xmb_ske,xmb_ave,xmb_std,xmb_cur,xmbweight
     real,    dimension (its:ite)::                           &
       pr_ske,pr_ave,pr_std,pr_cur
     real,    dimension (its:ite,jts:jte)::                           &
               pr_gr,pr_w,pr_mc,pr_st,pr_as,pr_capma,     &
               pr_capme,pr_capmi
     real, dimension (5) :: weight,wm,wm1,wm2,wm3
     real, dimension (its:ite,5) :: xmb_w

!
      character *(*), intent (in)        ::                           &
       name
!
     weight(1) = -999.  !this will turn off weights
     wm(1)=-999.

     tuning=0.
!
!
      DO k=kts,ktf
      do i=its,itf
        outtem(i,k)=0.
        outq(i,k)=0.
        outqc(i,k)=0.
        subt(i,k)=0.
        subq(i,k)=0.
        sub_mas(i,k)=0.
      enddo
      enddo
      do i=its,itf
        pre(i)=0.
        xmb(i)=0.
         xfac1(i)=0.
         xfac2(i)=0.
        xmbweight(i)=1.
      enddo
      do i=its,itf
        IF(ierr(i).eq.0)then
        do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
           if(pr_ens(i,j,n).le.0.)then
             xf_ens(i,j,n)=0.
           endif
        enddo
        endif
      enddo
!
!--- calculate ensemble average mass fluxes
!
       call massflx_stats(xf_ens,ensdim,nx2,nx,maxens3,      &
            xmb_ave,xmb_std,xmb_cur,xmb_ske,j,ierr,1,    &
            APR_GR,APR_W,APR_MC,APR_ST,APR_AS,           &
            APR_CAPMA,APR_CAPME,APR_CAPMI,               &
            pr_gr,pr_w,pr_mc,pr_st,pr_as,                &
            pr_capma,pr_capme,pr_capmi,                  &
            itf,jtf,ktf,                   &
            its,ite, jts,jte, kts,kte                   )
       xmb_w=0.
       call massflx_stats(pr_ens,ensdim,nx2,nx,maxens3,  &
            pr_ave,pr_std,pr_cur,pr_ske,j,ierr,2,        &
            APR_GR,APR_W,APR_MC,APR_ST,APR_AS,           &
            APR_CAPMA,APR_CAPME,APR_CAPMI,               &
            pr_gr,pr_w,pr_mc,pr_st,pr_as,                &
            pr_capma,pr_capme,pr_capmi,                  &
            itf,jtf,ktf,                   &
            its,ite, jts,jte, kts,kte                   )
!
!-- now do feedback
!
      ddtes=100.
      do i=its,itf
        if(ierr(i).eq.0)then
         if(xmb_ave(i).le.0.)then
              ierr(i)=13
              xmb_ave(i)=0.
         endif
         xmb(i)=max(.1*xmb_ave(i),xmb_ave(i)-tuning*xmb_std(i))
! --- Now use proper count of how many closures were actually
!       used in cup_forcing_ens (including screening of some
!       closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))
           if (xland1(i).lt.0.5)xmb(i)=xmb(i)*clos_wei
           if(xmb(i).eq.0.)then
              ierr(i)=19
           endif
           if(xmb(i).gt.100.)then
              ierr(i)=19
           endif
           xfac1(i)=xmb(i)
           xfac2(i)=xmb(i)

        endif
!       if(weight(1).lt.-100.)xfac1(i)=xmb_ave(i)
!       if(weight(1).lt.-100.)xfac2(i)=xmb_ave(i)
      ENDDO
      DO k=kts,ktf
      do i=its,itf
            dtt=0.
            dtts=0.
            dtq=0.
            dtqs=0.
            dtqc=0.
            dtpw=0.
        IF(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,nx
              dtt=dtt+dellat(i,k,n)
              dtts=dtts+subt_ens(i,k,n)
              dtq=dtq+dellaq(i,k,n)
              dtqs=dtqs+subq_ens(i,k,n)
              dtqc=dtqc+dellaqc(i,k,n)
              dtpw=dtpw+pw(i,k,n)
           enddo
           OUTTEM(I,K)=XMB(I)*dtt/float(nx)
           SUBT(I,K)=XMB(I)*dtts/float(nx)
           OUTQ(I,K)=XMB(I)*dtq/float(nx)
           SUBQ(I,K)=XMB(I)*dtqs/float(nx)
           OUTQC(I,K)=XMB(I)*dtqc/float(nx)
           PRE(I)=PRE(I)+XMB(I)*dtpw/float(nx)
           sub_mas(i,k)=zu(i,k)*xmb(i)
        endif
      enddo
      enddo

      do i=its,itf
        if(ierr(i).eq.0)then
        do k=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
          massfln(i,j,k)=massfln(i,j,k)*xfac1(i)
          xf_ens(i,j,k)=xf_ens(i,j,k)*xfac1(i)
        enddo
        endif
      ENDDO

   END SUBROUTINE cup_output_ens_3d


   SUBROUTINE cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,       &
              kbcon,ktop,ierr,                               &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,da
!
        do i=its,itf
         aa0(i)=0.
        enddo
        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.LE.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z(I,K)-Z(I,K-1)
         da=zu(i,k)*DZ*(9.81/(1004.*( &
                (T_cup(I,K)))))*DBY(I,K-1)/ &
             (1.+GAMMA_CUP(I,K))
         IF(K.eq.KTOP(I).and.da.le.0.)go to 100
         AA0(I)=AA0(I)+da
         if(aa0(i).lt.0.)aa0(i)=0.
100     continue

   END SUBROUTINE cup_up_aa0


   SUBROUTINE cup_up_he(k22,hkb,z_cup,cd,entr,he_cup,hc,     &
              kbcon,ierr,dby,he,hes_cup,                     &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! z_cup = heights of model cloud levels 
  ! entr = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hc,dby
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        hkb
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz
!
!--- moist static energy inside cloud
!
      do k=kts,ktf
      do i=its,itf
       hc(i,k)=0.
       DBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
       hkb(i)=0.
      enddo
      do i=its,itf
      if(ierr(i).eq.0.)then
      hkb(i)=he_cup(i,k22(i))
      do k=1,k22(i)
        hc(i,k)=he_cup(i,k)
!       DBY(I,K)=0.
      enddo
      do k=k22(i),kbcon(i)-1
        hc(i,k)=hkb(i)
!       DBY(I,K)=0.
      enddo
        k=kbcon(i)
        hc(i,k)=hkb(i)
        DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
      endif
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(k.gt.kbcon(i).and.ierr(i).eq.0.)then
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
           HC(i,K)=(HC(i,K-1)*(1.-.5*CD(i,K)*DZ)+entr* &
                DZ*HE(i,K-1))/(1.+entr*DZ-.5*cd(i,k)*dz)
           DBY(I,K)=HC(I,K)-HES_cup(I,K)
        endif
      enddo

      enddo

   END SUBROUTINE cup_up_he


   SUBROUTINE cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav,     &
              kbcon,ktop,cd,dby,mentr_rate,clw_all,                  &
              q,GAMMA_cup,zu,qes_cup,k22,qe_cup,xl,          &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        q,zu,gamma_cup,qe_cup,dby,qes_cup,z_cup,cd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        mentr_rate,xl
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav
!
!  local variables in this routine
!

     integer                              ::                           &
        iall,i,k
     real                                 ::                           &
        dh,qrch,c0,dz,radius
!
        iall=0
        c0=.002
!
!--- no precip for small clouds
!
        if(mentr_rate.gt.0.)then
          radius=.2/mentr_rate
          if(radius.lt.900.)c0=0.
!         if(radius.lt.900.)iall=0
        endif
        do i=its,itf
          pwav(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          qc(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qes_cup(i,k)
          clw_all(i,k)=0.
          qrc(i,k)=0.
        enddo
        enddo
      do i=its,itf
      if(ierr(i).eq.0.)then
      do k=k22(i),kbcon(i)-1
        qc(i,k)=qe_cup(i,k22(i))
      enddo
      endif
      enddo

        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Lt.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z_cup(i,K)-Z_cup(i,K-1)
!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
        QC(i,K)=(QC(i,K-1)*(1.-.5*CD(i,K)*DZ)+mentr_rate* &
                DZ*Q(i,K-1))/(1.+mentr_rate*DZ-.5*cd(i,k)*dz)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
         QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k) &
              /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------- LIQUID WATER CONTENT IN cloud after rainout
!
        clw_all(i,k)=QC(I,K)-QRCH
        QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ)
        if(qrc(i,k).lt.0.)then
          qrc(i,k)=0.
        endif
!
!-------   3.Condensation
!
         PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
        if(iall.eq.1)then
          qrc(i,k)=0.
          pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
          if(pw(i,k).lt.0.)pw(i,k)=0.
        endif
!
!----- set next level
!
         QC(I,K)=QRC(I,K)+qrch
!
!--- integrated normalized ondensate
!
         PWAV(I)=PWAV(I)+PW(I,K)
 100     CONTINUE

   END SUBROUTINE cup_up_moisture


   SUBROUTINE cup_up_nms(zu,z_cup,entr,cd,kbcon,ktop,ierr,k22,  &
              itf,jtf,ktf,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE

!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! cd= detrainment function 
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         z_cup,cd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
         entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
         ierr
   ! zu is the normalized mass flux

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
         zu
!
!  local variables in this routine
!

     integer                              ::                           &
         i,k
     real                                 ::                           &
         dz
!
!   initialize for this go around
!
       do k=kts,ktf
       do i=its,itf
         zu(i,k)=0.
       enddo
       enddo
!
! do normalized mass budget
!
       do i=its,itf
          IF(ierr(I).eq.0)then
             do k=k22(i),kbcon(i)
               zu(i,k)=1.
             enddo
             DO K=KBcon(i)+1,KTOP(i)
               DZ=Z_cup(i,K)-Z_cup(i,K-1)
               ZU(i,K)=ZU(i,K-1)*(1.+(entr-cd(i,k))*DZ)
             enddo
          endif
       enddo

   END SUBROUTINE cup_up_nms

!====================================================================
   SUBROUTINE g3init(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,           &
                        MASS_FLUX,cp,restart,                       &
                        P_QC,P_QI,P_FIRST_SCALAR,                   &
                        RTHFTEN, RQVFTEN,                           &
                        APR_GR,APR_W,APR_MC,APR_ST,APR_AS,          &
                        APR_CAPMA,APR_CAPME,APR_CAPMI,              &
                        cugd_tten,cugd_ttens,cugd_qvten,            &
                        cugd_qvtens,cugd_qcten,                     &
                        allowed_to_read,                            &
                        ids, ide, jds, jde, kds, kde,               &
                        ims, ime, jms, jme, kms, kme,               &
                        its, ite, jts, jte, kts, kte               )
!--------------------------------------------------------------------   
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           ::  restart,allowed_to_read
   INTEGER , INTENT(IN)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   INTEGER , INTENT(IN)           ::  P_FIRST_SCALAR, P_QI, P_QC
   REAL,     INTENT(IN)           ::  cp

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          CUGD_TTEN,         &
                                                          CUGD_TTENS,        &
                                                          CUGD_QVTEN,        &
                                                          CUGD_QVTENS,       &
                                                          CUGD_QCTEN
   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHCUTEN, &
                                                          RQVCUTEN, &
                                                          RQCCUTEN, &
                                                          RQICUTEN   

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHFTEN,  &
                                                          RQVFTEN

   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(OUT) ::        &
                                APR_GR,APR_W,APR_MC,APR_ST,APR_AS,  &
                                APR_CAPMA,APR_CAPME,APR_CAPMI,      &
                                MASS_FLUX

   INTEGER :: i, j, k, itf, jtf, ktf
 
   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)
 
   IF(.not.restart)THEN
     DO j=jts,jte
     DO k=kts,kte
     DO i=its,ite
        RTHCUTEN(i,k,j)=0.
        RQVCUTEN(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO
     DO j=jts,jte
     DO k=kts,kte
     DO i=its,ite
       cugd_tten(i,k,j)=0.
       cugd_ttens(i,k,j)=0.
       cugd_qvten(i,k,j)=0.
       cugd_qvtens(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RTHFTEN(i,k,j)=0.
        RQVFTEN(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

     IF (P_QC .ge. P_FIRST_SCALAR) THEN
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           RQCCUTEN(i,k,j)=0.
           cugd_qcten(i,k,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

     IF (P_QI .ge. P_FIRST_SCALAR) THEN
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           RQICUTEN(i,k,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

     DO j=jts,jtf
     DO i=its,itf
        mass_flux(i,j)=0.
     ENDDO
     ENDDO

   ENDIF
     DO j=jts,jtf
     DO i=its,itf
        APR_GR(i,j)=0.
        APR_ST(i,j)=0.
        APR_W(i,j)=0.
        APR_MC(i,j)=0.
        APR_AS(i,j)=0.
        APR_CAPMA(i,j)=0.
        APR_CAPME(i,j)=0.
        APR_CAPMI(i,j)=0.
     ENDDO
     ENDDO

   END SUBROUTINE g3init


   SUBROUTINE massflx_stats(xf_ens,ensdim,maxens,maxens2,maxens3, &
              xt_ave,xt_std,xt_cur,xt_ske,j,ierr,itest,           &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                  &
              APR_CAPMA,APR_CAPME,APR_CAPMI,                      &
              pr_gr,pr_w,pr_mc,pr_st,pr_as,                       &
              pr_capma,pr_capme,pr_capmi,                         &
              itf,jtf,ktf,  &
              its,ite, jts,jte, kts,kte)

   IMPLICIT NONE

   integer, intent (in   )              ::                                    &
                     j,ensdim,maxens3,maxens,maxens2,itest
   INTEGER,      INTENT(IN   ) ::                                             &
                                  itf,jtf,ktf,                  &
                                  its,ite, jts,jte, kts,kte


     real, dimension (its:ite)                                                &
         , intent(inout) ::                                                   &
           xt_ave,xt_cur,xt_std,xt_ske
     integer, dimension (its:ite), intent (in) ::                             &
           ierr
     real, dimension (its:ite,jts:jte,1:ensdim)                               &
         , intent(in   ) ::                                                   &
           xf_ens
     real, dimension (its:ite,jts:jte)                                        &
         , intent(inout) ::                                                   &
           APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                                 &
           APR_CAPMA,APR_CAPME,APR_CAPMI
     real, dimension (its:ite,jts:jte)                                        &
         , intent(inout) ::                                                   &
           pr_gr,pr_w,pr_mc,pr_st,pr_as,                                      &
           pr_capma,pr_capme,pr_capmi

!
! local stuff
!
     real, dimension (its:ite , 1:maxens3 )       ::                          &
           x_ave,x_cur,x_std,x_ske
     real, dimension (its:ite , 1:maxens  )       ::                          &
           x_ave_cap


      integer, dimension (1:maxens3) :: nc1
      integer :: i,k
      integer :: num,kk,num2,iedt
      real :: a3,a4

      num=ensdim/maxens3
      num2=ensdim/maxens
      if(itest.eq.1)then
      do i=its,ite
       pr_gr(i,j) =  0.
       pr_w(i,j) =  0.
       pr_mc(i,j) = 0.
       pr_st(i,j) = 0.
       pr_as(i,j) = 0.
       pr_capma(i,j) =  0.
       pr_capme(i,j) = 0.
       pr_capmi(i,j) = 0.
      enddo
      endif

      do k=1,maxens
      do i=its,ite
        x_ave_cap(i,k)=0.
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        x_ave(i,k)=0.
        x_std(i,k)=0.
        x_ske(i,k)=0.
        x_cur(i,k)=0.
      enddo
      enddo
      do i=its,ite
        xt_ave(i)=0.
        xt_std(i)=0.
        xt_ske(i)=0.
        xt_cur(i)=0.
      enddo
      do kk=1,num
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave(i,k)=x_ave(i,k)+xf_ens(i,j,maxens3*(kk-1)+k)
        endif
      enddo
      enddo
      enddo
      do iedt=1,maxens2
      do k=1,maxens
      do kk=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave_cap(i,k)=x_ave_cap(i,k)                               &
            +xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+kk)
        endif
      enddo
      enddo
      enddo
      enddo
      do k=1,maxens
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave_cap(i,k)=x_ave_cap(i,k)/float(num2)
        endif
      enddo
      enddo

      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave(i,k)=x_ave(i,k)/float(num)
        endif
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        xt_ave(i)=xt_ave(i)+x_ave(i,k)
        endif
      enddo
      enddo
      do i=its,ite
        if(ierr(i).eq.0)then
        xt_ave(i)=xt_ave(i)/float(maxens3)
        endif
      enddo
!
!--- now do std, skewness,curtosis
!
      do kk=1,num
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.x_ave(i,k).gt.0.)then
!       print *,i,j,k,kk,x_std(i,k),xf_ens(i,j,maxens3*(kk-1)+k),x_ave(i,k)
        x_std(i,k)=x_std(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**2
        x_ske(i,k)=x_ske(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**3
        x_cur(i,k)=x_cur(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**4
        endif
      enddo
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.xt_ave(i).gt.0.)then
        xt_std(i)=xt_std(i)+(x_ave(i,k)-xt_ave(i))**2
        xt_ske(i)=xt_ske(i)+(x_ave(i,k)-xt_ave(i))**3
        xt_cur(i)=xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
        endif
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.x_std(i,k).gt.0.)then
           x_std(i,k)=x_std(i,k)/float(num)
           a3=max(1.e-6,x_std(i,k))
           x_std(i,k)=sqrt(a3)
           a3=max(1.e-6,x_std(i,k)**3)
           a4=max(1.e-6,x_std(i,k)**4)
           x_ske(i,k)=x_ske(i,k)/float(num)/a3
           x_cur(i,k)=x_cur(i,k)/float(num)/a4
        endif
!       print*,'                               '
!       print*,'Some statistics at gridpoint i,j, ierr',i,j,ierr(i)
!       print*,'statistics for closure number ',k
!       print*,'Average= ',x_ave(i,k),'  Std= ',x_std(i,k)
!       print*,'Skewness= ',x_ske(i,k),' Curtosis= ',x_cur(i,k)
!       print*,'                               '

      enddo
      enddo
      do i=its,ite
        if(ierr(i).eq.0.and.xt_std(i).gt.0.)then
           xt_std(i)=xt_std(i)/float(maxens3)
           a3=max(1.e-6,xt_std(i))
           xt_std(i)=sqrt(a3)
           a3=max(1.e-6,xt_std(i)**3)
           a4=max(1.e-6,xt_std(i)**4)
           xt_ske(i)=xt_ske(i)/float(maxens3)/a3
           xt_cur(i)=xt_cur(i)/float(maxens3)/a4
!       print*,'                               '
!       print*,'Total ensemble independent statistics at i =',i
!       print*,'Average= ',xt_ave(i),'  Std= ',xt_std(i)
!       print*,'Skewness= ',xt_ske(i),' Curtosis= ',xt_cur(i)
!       print*,'                               '
!
!  first go around: store massflx for different closures/caps
!
      if(itest.eq.1)then
       pr_gr(i,j) = .25*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)+x_ave(i,13))
       pr_w(i,j) = .25*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,14))
       pr_mc(i,j) = .25*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)+x_ave(i,15))
       pr_st(i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
       pr_as(i,j) = x_ave(i,16)
       pr_capma(i,j) = x_ave_cap(i,1)
       pr_capme(i,j) = x_ave_cap(i,2)
       pr_capmi(i,j) = x_ave_cap(i,3)
!
!  second go around: store preciprates (mm/hour) for different closures/caps
!
        else if (itest.eq.2)then
       APR_GR(i,j)=.25*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)+x_ave(i,13))*      &
                  3600.*pr_gr(i,j) +APR_GR(i,j)
       APR_W(i,j)=.25*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,14))*       &
                  3600.*pr_w(i,j) +APR_W(i,j)
       APR_MC(i,j)=.25*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)+x_ave(i,15))*      &
                  3600.*pr_mc(i,j) +APR_MC(i,j)
       APR_ST(i,j)=.333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))*   &
                  3600.*pr_st(i,j) +APR_ST(i,j)
       APR_AS(i,j)=x_ave(i,16)*                       &
                  3600.*pr_as(i,j) +APR_AS(i,j)
       APR_CAPMA(i,j) = x_ave_cap(i,1)*                          &
                  3600.*pr_capma(i,j) +APR_CAPMA(i,j)
       APR_CAPME(i,j) = x_ave_cap(i,2)*                          &
                  3600.*pr_capme(i,j) +APR_CAPME(i,j)
       APR_CAPMI(i,j) = x_ave_cap(i,3)*                          &
                  3600.*pr_capmi(i,j) +APR_CAPMI(i,j)
        endif
        endif
      enddo

   END SUBROUTINE massflx_stats

   SUBROUTINE cup_axx(tcrit,kbmax,z1,p,psur,xl,rv,cp,tx,qx,axx,ierr,    &
           cap_max,cap_max_increment,entr_rate,mentr_rate,&
           j,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte,ens4)
   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) ::                                             &
                                  j,itf,jtf,ktf,                &
                                  its,ite, jts,jte, kts,kte,ens4
     real, dimension (its:ite,kts:kte,1:ens4)                                 &
         , intent(inout) ::                                                   &
           tx,qx
     real, dimension (its:ite,kts:kte)                                 &
         , intent(in) ::                                                   &
           p
     real, dimension (its:ite)                                 &
         , intent(in) ::                                                   &
           z1,psur,cap_max,cap_max_increment
     real, intent(in) ::                                                   &
           tcrit,xl,rv,cp,mentr_rate,entr_rate
     real, dimension (its:ite,1:ens4)                                 &
         , intent(out) ::                                                   &
           axx
     integer, dimension (its:ite), intent (in) ::                             &
           ierr,kbmax
     integer, dimension (its:ite) ::                             &
           ierrxx,k22xx,kbconxx,ktopxx,kstabm,kstabi
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: e,tvbar
     integer n,i,k,iph
     real,    dimension (its:ite,kts:kte) ::                           &
        he,hes,qes,z,                                                  &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        tn_cup,                                                        &
        dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,cd

     real,    dimension (its:ite) ::                                   &
       AA0,HKB,QKB,          &
       PWAV,BU
      do n=1,ens4
      do i=its,ite
       axx(i,n)=0.
      enddo
      enddo
     HT(1)=XL/CP
     HT(2)=2.834E6/CP
     BE(1)=.622*HT(1)/.286
     AE(1)=BE(1)/273.+ALOG(610.71)
     BE(2)=.622*HT(2)/.286
     AE(2)=BE(2)/273.+ALOG(610.71)
!
!
     do 100 n=1,ens4

      do k=kts,ktf
      do i=its,itf
        cd(i,k)=0.1*entr_rate
      enddo
      enddo


      do i=its,itf
        ierrxx(i)=ierr(i)
        k22xx(i)=1
        kbconxx(i)=1
        ktopxx(i)=1
        kstabm(i)=ktf-1
      enddo
      DO k=kts,ktf
      do i=its,itf
        if(ierrxx(i).eq.0)then
        IPH=1
        IF(Tx(I,K,n).LE.TCRIT)IPH=2
        E=EXP(AE(IPH)-BE(IPH)/TX(I,K,N))
        QES(I,K)=.622*E/(100.*P(I,K)-E)
        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(Qx(I,K,N).GT.QES(I,K))Qx(I,K,N)=QES(I,K)
        TV(I,K)=Tx(I,K,N)+.608*Qx(I,K,N)*Tx(I,K,N)
        endif
      enddo
      enddo
!
         do i=its,itf
           if(ierrxx(i).eq.0)then
             Z(I,KTS)=max(0.,Z1(I))-(ALOG(P(I,KTS))- &
                 ALOG(PSUR(I)))*287.*TV(I,KTS)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierrxx(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierrxx(i).eq.0)then
         HE(I,K)=9.81*Z(I,K)+1004.*Tx(I,K,n)+2.5E06*Qx(I,K,n)
         HES(I,K)=9.81*Z(I,K)+1004.*Tx(I,K,n)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

! cup levels
!
      do k=kts+1,ktf
      do i=its,itf
        if(ierrxx(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(qx(i,k-1,n)+qx(i,k,n))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(tx(i,k-1,n)+tx(i,k,n))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierrxx(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=qx(i,1,n)
        hes_cup(i,1)=hes(i,1)
        he_cup(i,1)=he(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        t_cup(i,1)=tx(i,1,n)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      CALL cup_MAXIMI(HE_CUP,3,KBMAX,K22XX,ierrxx, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
       DO 36 i=its,itf
         IF(ierrxx(I).eq.0.)THEN
         IF(K22xx(I).GE.KBMAX(i))ierrxx(i)=2
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      call cup_kbcon(cap_max_increment,1,k22xx,kbconxx,he_cup,hes_cup, &
           ierrxx,kbmax,p_cup,cap_max, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEs_cup,Kbconxx,kstabm,kstabi,ierrxx,  &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
      IF(ierrxx(I).eq.0.)THEN
        if(kstabm(i)-1.gt.kstabi(i))then
           do k=kstabi(i),kstabm(i)-1
             cd(i,k)=cd(i,k-1)+1.5*entr_rate
             if(cd(i,k).gt.10.0*entr_rate)cd(i,k)=10.0*entr_rate
           enddo
        ENDIF
      ENDIF
      ENDDO
!
!--- calculate incloud moist static energy
!
      call cup_up_he(k22xx,hkb,z_cup,cd,mentr_rate,he_cup,hc, &
           kbconxx,ierrxx,dby,he,hes_cup, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)

!--- DETERMINE CLOUD TOP - KTOP
!
      call cup_ktop(1,dby,kbconxx,ktopxx,ierrxx, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!c--- normalized updraft mass flux profile
!
      call cup_up_nms(zu,z_cup,mentr_rate,cd,kbconxx,ktopxx,ierrxx,k22xx, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbconxx,ktopxx,ierrxx,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
       if(ierrxx(i).eq.0)axx(i,n)=aa0(i)
      enddo
100   continue
     END SUBROUTINE cup_axx

      SUBROUTINE conv_grell_spread3d(rthcuten,rqvcuten,raincv,            &
     &         cugd_avedx,cugd_tten,cugd_qvten,cugd_ttens,                &
     &         cugd_qvtens,cugd_qcten,pi_phy,moist_qv,pratec,dt,num_tiles,&
     &         imomentum,                                                 &
     &         ids, ide, jds, jde, kds, kde,                              &
     &         ips, ipe, jps, jpe, kps, kpe,                              &
     &         ims, ime, jms, jme, kms, kme,                              &
     &         its, ite, jts, jte, kts, kte,                              &
     &         fqc, fqi, rqccuten, rqicuten                               )

!

   INTEGER,      INTENT(IN   )    ::   num_tiles,imomentum
   INTEGER,  INTENT(IN) ::                       &
     &           its,ite,jts,jte
   INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde,&
                                       ims,ime, jms,jme, kms,kme, &
                                       ips,ipe, jps,jpe, kps,kpe, &
                                       kts,kte,cugd_avedx
   REAL, DIMENSION (ims:ime,kms:kme,jms:jme), optional,INTENT (INOUT) ::  &
     &  rthcuten,rqvcuten,cugd_tten,                              &
     &  cugd_qvten,cugd_ttens,cugd_qvtens,cugd_qcten
   REAL, DIMENSION (ims:ime,kms:kme,jms:jme),INTENT (IN) ::       &
          moist_qv
   REAL, DIMENSION (ims:ime,kms:kme,jms:jme), INTENT (IN) ::      &
          PI_PHY
   REAL, DIMENSION (ims:ime,jms:jme), INTENT (INOUT) ::           &
          RAINCV,PRATEC
   REAL,                              INTENT(IN) ::   dt
   REAL, DIMENSION (ims:ime,kms:kme,jms:jme), optional,INTENT (INOUT) ::  &
     &  rqccuten, rqicuten

   LOGICAL, OPTIONAL,                 INTENT(IN) :: fqc, fqi

   INTEGER                        :: ikk1,ikk2,ikk11,i,j,k,kk,nn,smoothh,smoothv
   INTEGER                        :: ifs,ife,jfs,jfe,ido,jdo,cugd_spread
   LOGICAL                        :: new, f_qc, f_qi
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   REAL, DIMENSION (ips-2:ipe+2,kps:kpe,jps-2:jpe+2) ::     &
          rthcutent,rqvcutent
   real, dimension (ips-2:ipe+2,jps-2:jpe+2) :: qmem
   real, dimension (ips-1:ipe+1,jps-1:jpe+1) :: smtt,smtq
   real, dimension (kps:kpe) :: conv_trasht,conv_trashq
   REAL                           :: qmem1,qmem2,qmemf,thresh

      f_qc = .false.
      if (present(fqc)) f_qc = fqc
      f_qi = .false.
      if (present(fqi)) f_qi = fqi

      smoothh=1
      smoothv=1
      cugd_spread=cugd_avedx/2
! SET START AND END POINTS FOR TILES
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij ,ifs,ife,jfs,jfe,its,ite,jts,jte, i,j,k,kk,nn,ikk1,ikk2,ikk11) &
      !$OMP PRIVATE ( ido,jdo,qmemf,qmem1,qmem2,qmem,thresh,conv_trasht,conv_trashq,smtt,smtq)


        do j=jts-2,jte+2
        do i=its-2,ite+2
          qmem(i,j)=1.
        enddo
        enddo
        do j=jts-1,jte+1
        do i=its-1,ite+1
          smtt(i,j)=0.
          smtq(i,j)=0.
        enddo
        enddo
        do j=jts,jte              
        do k=kts,kte              
        do i=its,ite
          rthcuten(i,k,j)=0. 
          rqvcuten(i,k,j)=0.      
        enddo
        enddo
        enddo
        do j=jts-2,jte+2              
        do k=kts,kte              
        do i=its-2,ite+2
          rthcutent(i,k,j)=0. 
          rqvcutent(i,k,j)=0.      
        enddo
        enddo
        enddo
!     
        ifs=max(its,ids)
        jfs=max(jts,jds)
        ife=min(ite,ide-1)
        jfe=min(jte,jde-1)
!
!       
!  
! prelims finished, now go real for every grid point
!  
   ifs=max(its,ids)
   ife=min(ite,ide-1)
   jfs=max(jts,jds)
   jfe=min(jte,jde-1)
   if(cugd_spread.gt.0.or.smoothh.eq.1)then
      if(its.eq.ips)ifs=max(its-1,ids)
      if(ite.eq.ipe)ife=min(ite+1,ide-1)
      if(jts.eq.jps)jfs=max(jts-1,jds)
      if(jte.eq.jpe)jfe=min(jte+1,jde-1)
   endif
        do j=jfs,jfe
        do i=ifs,ife
!
        do k=kts,kte
        rthcutent(i,k,j)=cugd_tten(i,k,j)
        rqvcutent(i,k,j)=cugd_qvten(i,k,j)
        enddo
!
! for high res run, spread the subsidence
! this is tricky......only consider grid points where there was no rain,
! so cugd_tten and such are zero!
!
        if(cugd_spread.gt.0)then
        do k=kts,kte
          do nn=-1,1,1
            jdo=max(j+nn,jds)
            jdo=min(jdo,jde-1)
            do kk=-1,1,1
             ido=max(i+kk,ids)
             ido=min(ido,ide-1)
             rthcutent(i,k,j)=rthcutent(i,k,j)     &
                                  +qmem(ido,jdo)*cugd_ttens(ido,k,jdo)
             rqvcutent(i,k,j)=rqvcutent(i,k,j)     &
                                  +qmem(ido,jdo)*cugd_qvtens(ido,k,jdo)
          enddo
          enddo
        enddo
        endif
!       
! end spreading
       
        if(cugd_spread.eq.0)then
        do k=kts,kte
           rthcutent(i,k,j)=rthcutent(i,k,j)+cugd_ttens(i,k,j)
           rqvcutent(i,k,j)=rqvcutent(i,k,j)+cugd_qvtens(i,k,j)
        enddo
        endif
        enddo  ! end j
        enddo  ! end i
! smooth
        do k=kts,kte
          if(smoothh.eq.0)then
             ifs=max(its,ids+4)
             ife=min(ite,ide-5)
             jfs=max(jts,jds+4)
             jfe=min(jte,jde-5)
             do i=ifs,ife
             do j=jfs,jfe
               rthcuten(i,k,j)=rthcutent(i,k,j)
               rqvcuten(i,k,j)=rqvcutent(i,k,j)
             enddo  ! end j
             enddo  ! end j
          else if(smoothh.eq.1)then   ! smooth
             ifs=max(its,ids)
             ife=min(ite,ide-1)
             jfs=max(jts,jds)
             jfe=min(jte,jde-1)
! we need an extra row for j (halo comp)
             if(jts.eq.jps)jfs=max(jts-1,jds)
             if(jte.eq.jpe)jfe=min(jte+1,jde-1)
             do i=ifs,ife
             do j=jfs,jfe
                smtt(i,j)=.25*(rthcutent(i-1,k,j)+2.*rthcutent(i,k,j)+rthcutent(i+1,k,j))
                smtq(i,j)=.25*(rqvcutent(i-1,k,j)+2.*rqvcutent(i,k,j)+rqvcutent(i+1,k,j))
             enddo  ! end j
             enddo  ! end j
             ifs=max(its,ids+4)
             ife=min(ite,ide-5)
             jfs=max(jts,jds+4)
             jfe=min(jte,jde-5)
             do i=ifs,ife
               do j=jfs,jfe
                 rthcuten(i,k,j)=.25*(smtt(i,j-1)+2.*smtt(i,j)+smtt(i,j+1))
                 rqvcuten(i,k,j)=.25*(smtq(i,j-1)+2.*smtq(i,j)+smtq(i,j+1))
              enddo  ! end j
              enddo  ! end i
           endif  ! smoothh

        enddo  ! end k
!       
! check moistening rates
!         
          ifs=max(its,ids+4)
          ife=min(ite,ide-5)
          jfs=max(jts,jds+4)
          jfe=min(jte,jde-5)
        do j=jfs,jfe
        do i=ifs,ife
         qmemf=1.
        thresh=1.e-20
        do k=kts,kte              
         if(rqvcuten(i,k,j).lt.0.)then

           qmem1=moist_qv(i,k,j)+rqvcuten(i,k,j)*dt
           if(qmem1.lt.thresh)then
              qmem1=rqvcuten(i,k,j)
              qmem2=(thresh-moist_qv(i,k,j))/dt
              qmemf=min(qmemf,qmem2/qmem1)
              qmemf=max(0.,qmemf)
              qmemf=min(1.,qmemf)
           endif

         endif
         enddo
        do k=kts,kte
         rqvcuten(i,k,j)=rqvcuten(i,k,j)*qmemf
         rthcuten(i,k,j)=rthcuten(i,k,j)*qmemf
        enddo
        if(present(rqccuten))then
          if(f_qc) then
           do k=kts,kte
               rqccuten(i,k,j)=rqccuten(i,k,j)*qmemf
           enddo
          endif
        endif
        if(present(rqicuten))then
          if(f_qi) then
           do k=kts,kte
               rqicuten(i,k,j)=rqicuten(i,k,j)*qmemf
           enddo
          endif
        endif
        RAINCV(I,J)=RAINCV(I,J)*qmemf
        PRATEC(I,J)=PRATEC(I,J)*qmemf
!
! check heating rates

!
        thresh=200.
        qmemf=1.
        qmem1=0.
        do k=kts,kte
        qmem1=abs(rthcuten(i,k,j))*86400. 

        if(qmem1.gt.thresh)then
          qmem2=thresh/qmem1
          qmemf=min(qmemf,qmem2)
          qmemf=max(0.,qmemf) 
        endif

        enddo
        RAINCV(I,J)=RAINCV(I,J)*qmemf
        PRATEC(I,J)=PRATEC(I,J)*qmemf
        do k=kts,kte
           rqvcuten(i,k,j)=rqvcuten(i,k,j)*qmemf
           rthcuten(i,k,j)=rthcuten(i,k,j)*qmemf
        enddo
        if(present(rqccuten))then
          if(f_qc) then
             do k=kts,kte
               rqccuten(i,k,j)=rqccuten(i,k,j)*qmemf
             enddo
          endif
        endif
        if(present(rqicuten))then
          if(f_qi) then
             do k=kts,kte
               rqicuten(i,k,j)=rqicuten(i,k,j)*qmemf
             enddo
          endif
        endif
        if(smoothv.eq.1)then
! 
! smooth for now
!
           do k=kts+2,kte-2
               conv_trasht(k)= .25*(rthcuten(i,k-1,j)+2.*rthcuten(i,k,j)+rthcuten(i,k+1,j))
               conv_trashq(k)= .25*(rqvcuten(i,k-1,j)+2.*rqvcuten(i,k,j)+rqvcuten(i,k+1,j))
           enddo
           do k=kts+2,kte-2
                rthcuten(i,k,j)=conv_trasht(k)
                rqvcuten(i,k,j)=conv_trashq(k)
           enddo
        endif
        do k=kts,kte
            rthcuten(i,k,j)=rthcuten(i,k,j)/pi_phy(i,k,j)
        enddo
        enddo  ! end j
        enddo  ! end i


   END SUBROUTINE CONV_GRELL_SPREAD3D
!-------------------------------------------------------
END MODULE dep_cu_g3_mod
