module wetdep_prep_mod

  use chem_rc_mod
  use chem_types_mod

contains

   subroutine wetdep_prep(tr3d,tk3d,pr3d,prl3d,ph3d,phl3d,          &
                       us3d,vs3d,ws3d,exch,dqdt,           &
                       pb2d,&
                       rcav,raincv_b,     &
!                       t_phy,moist,u_phy,v_phy,p_phy,chem,       &
                       t_phy,moist,u_phy,v_phy,p_phy,       &
                       g,rd,p1000,cp,              &
!                       t8w,p8w,pbl,z_at_w,zmid,dz8w,vvel,kt_turb,&
                       p8w,pbl,z_at_w,zmid,dz8w,vvel,kt_turb,dqdti, &
                       rho_phy,nsol,num_moist,ntra,            &
                       ids,ide, jds,jde, kds,kde,                                    &
                       ims,ime, jms,jme, kms,kme,                                    &
                       its,ite, jts,jte, kts,kte, rc)
    use chem_const_mod, only : rgasuniv
    use raqmschem_const_mod, only : epsilc
    use raqmschem_species_mod, only : p_atm_shum,p_atm_cldq,idwetd
    use raqmschem_pmgrid_mod,only : iam,iamprn,iprnin,jprnin,tile
    IMPLICIT NONE

    ! -- input variables

    INTEGER,      INTENT(IN) :: nsol,num_moist,ntra
    INTEGER,      INTENT(IN) ::   ids,ide, jds,jde, kds,kde,      &
                                   ims,ime, jms,jme, kms,kme,                          &
                                   its,ite, jts,jte, kts,kte
    REAL(CHEM_KIND_R4), INTENT(IN) :: g,rd,p1000,cp

    ! -- input pointers: indexing must always start from 1
    real(CHEM_KIND_R4), dimension(:, :), intent(in) :: pb2d

    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: ph3d
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: phl3d
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: pr3d
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: prl3d
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: tk3d
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: us3d
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: vs3d
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: ws3d
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: exch
    real(CHEM_KIND_R8), dimension(:, :, :), intent(in) :: dqdt

    real(CHEM_KIND_R4), dimension(:, :, :, :), intent(inout)  :: tr3d

    ! -- I/O arrays
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: rcav
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: raincv_b
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: t_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: u_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: v_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: p_phy
!    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, nsol), intent(out) :: chem
!    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: p8w,t8w ! dont do kme+1
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: p8w
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: pbl
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w ! need kte+1 which is kme
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: zmid
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: dz8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: vvel
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: kt_turb 
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: dqdti
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: rho_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme)  :: rri,convfac

    INTEGER, OPTIONAL, INTENT(OUT) :: rc

    ! -- local variables
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n
    real(CHEM_KIND_R4) ::  maxv,factor,factor2,pu,pl,aln,pwant,rlat

    ! .. Intrinsic Functions ..
    INTRINSIC max, min, float

    ! -- begin
!    print *,'chem_prep: begin...'
    if (present(rc)) rc = CHEM_RC_SUCCESS
!    write(6,*)'raqms prep kds',kds,kde,'kms',kms,kme,'kts',kts,kte
!    call flush(6)


    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,jp,1)/g)
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=(ph3d(ip,jp,kp+1)-ph3d(ip,jp,kp))/g
          if (dz8w(i,k,j) < 0.) dz8w(i,k,j)=-dz8w(i,k,j)
!          if(k.ne.kte)then
            z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
!          endif
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
!      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          p8w(i,k,j)=pr3d(ip,jp,kp)
!         if(iam.eq.iamprn)then
!         if(i.eq.iprnin.and.j.eq.jprnin)then
!           write(250+iam,*)'p8w',k,p8w(i,k,j),'ip',ip,'jp',jp,'kp',kp
!           write(250+iam,*)'kts',kts,'kte',kte,'shape',shape(p8w),'pr3d',shape(pr3d)
!           call flush(250+iam)
!         endif
!         endif

        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
        ip = i - its + 1
        raincv_b(i,j)=rcav(i,j)
        pbl(i,j)=pb2d(ip,jp)
      enddo
    enddo


    factor=0.
    jmax=0
    jmaxi=0
    k=1

    do j=jts,jte
      jp = j - jts + 1
!      do k=kts,kte+1
      do k=kts,kte
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          ip = i - its + 1
          zmid(i,k,j)=phl3d(ip,jp,kkp)/g
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(ip,jp,kkp)
          p_phy(i,k,j)=prl3d(ip,jp,kkp)
          u_phy(i,k,j)=us3d(ip,jp,kkp)
          v_phy(i,k,j)=vs3d(ip,jp,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(RD*T_phy(i,k,j)*(1.+.608*tr3d(ip,jp,kkp,p_atm_shum)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-ws3d(ip,jp,kkp)*rri(i,k,j)/g
          kt_turb(i,k,j)=exch(ip,jp,kkp)
          dqdti(i,k,j)=dqdt(ip,jp,kkp)
!          if(tile.eq.3.and.i.eq.16.and.j.eq.11)then
!            write(200+iam,*)'dqdt',ip,jp,kkp,dqdt(ip,jp,kkp),'t',t_phy(i,k,j)
!          endif
          convfac(i,k,j)=p_phy(i,k,j)/rgasuniv/t_phy(i,k,j)
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=tr3d(ip,jp,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=tr3d(ip,jp,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=tr3d(ip,jp,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
        enddo
      enddo
    enddo

!    do j=jts,jte
!      do k=2,kte
!        do i=its,ite
!          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
!        enddo
!      enddo
    !enddo

    ! -- only used in phtolysis....
!    do j=jts,jte
!      do i=its,ite
!        t8w(i,1,j)=t_phy(i,1,j)
!        t8w(i,kte+1,j)=t_phy(i,kte,j)
!      enddo
!    enddo
!    write(6,*)'nsol',nsol,'num_moist',num_moist
!    call flush(6)
!#ifdef CHEMLOCAL
    do nv=1,nsol
      do j=jts,jte
        jp = j - jts + 1
!        do k=kts,kte+1
        do k=kts,kte ! ajl
          kk=min(k,kte)
          kkp = kk - kts + 1
          do i=its,ite
            ip = i - its + 1
!            chem(i,k,j,nv)=max(epsilc,tr3d(ip,jp,kkp,idwetd(nv)))
            tr3d(ip,jp,kkp,ntra+nv)=max(epsilc,tr3d(ip,jp,kkp,ntra+nv))
          enddo
        enddo
      enddo
    enddo
!#endif
   end subroutine wetdep_prep
end module wetdep_prep_mod
