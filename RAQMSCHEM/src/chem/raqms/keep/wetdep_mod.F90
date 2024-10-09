module wetdep_mod
contains
  subroutine wetdep_advance(advancecount,dts,rc2d,rn2d,ph3d,phl3d,pr3d,prl3d,tk3d,us3d,vs3d,ws3d,tr3d_inout4, &
    num_chem,num_moist,ntra,its, ite, jts, jte, kts, kte, &
    ims, ime, jms, jme, kms, kme, &
    rc)
    use chem_types_mod
    use chem_rc_mod
    use raqmschemcomm_mod, only : rcav,rnav,pblht
    use wetdep_prep_mod, only : wetdep_prep
    use chem_const_mod,  only : cp, grvity,rv,xlv, mwdry, p1000, rd, epsilc
    use dep_ctrans_grell_mod, only : grelldrvct
    use raqmschem_species_mod, only : nsol,nchemfull,chemname,ichemfullpt,p_qc,chemfull
    use raqmschem_species_mod, only : p_brcl,p_brno3
    use dep_wet_ls_mod
    use wetdep_alpha_mod, only : initwetdep
    use wetdep_fam_mod, only : initsubpointers
    use raqmschem_pmgrid_mod,only : iam,ibeg,tile,iamprn,iprn,jprn,iprnin,jprnin
    use raqmschem_pmgrid_mod,only : jbeg
    implicit none
    integer, intent(in) :: advancecount,its,ite,jts,jte,kts,kte
    integer, intent(in) :: ims,ime,jms,jme,kms,kme,num_moist,num_chem,ntra
    integer, intent(out), optional :: rc
    real(CHEM_KIND_R8), intent(in) :: dts
    real(CHEM_KIND_R8), dimension(:,:), intent(in) :: rn2d,rc2d
    real(CHEM_KIND_R8), dimension(:,:,:), intent(in) :: ph3d,phl3d,pr3d,prl3d,tk3d 
    real(CHEM_KIND_R8), dimension(:,:,:), intent(in) :: us3d,vs3d,ws3d
    real(CHEM_KIND_R4), dimension(:,:,:,:), intent(inout) :: tr3d_inout4
!   tr3d_in4 will be inout so wetdep can update tracers at end to pass to chemistry
!   these need to be in a module so are saved
!    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rcav
!    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rnav
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: dz8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme+1, jms:jme) :: p8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: rho_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: rri
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme+1, jms:jme) :: t8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: t_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: u_phy,raincv_v
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: v_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: vvel
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme+1, jms:jme) :: z_at_w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: zmid
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_moist) :: moist
    real(CHEM_KIND_R4) :: dt
    integer :: ids, ide, jds, jde, kds, kde
    integer :: jp,j,ip,i,k
!   temporary variables
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: pbl 
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: raincv_b
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: relhum
!    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:nsol) :: chem
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem) :: tr_fall !  convective wetdep
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem) :: var_rmv !  large scale wetdep
    real(CHEM_KIND_R4), dimension(kms:kme) :: brclold,brno3old
    logical first
    save first
    data first/.true./
    if(first)then
      first=.false.
      call initsubpointers
    endif
    if (present(rc)) rc = 0

    ! -- set domain
    ids = ims
    ide = ime
    jds = jms
    jde = jme
    kds = kms
    kde = kme
    ! -- initialize local arrays
    raincv_b   = 0._CHEM_KIND_R4

    dt = real(dts, kind=CHEM_KIND_R4)
!    write(6,*)'advnacecount',advancecount
    if(advancecount <= 0 )then
        jp=0
        do j=jts, jte
        jp = jp + 1
          ip = 0
         do i= its, ite
            ip = ip + 1
      rcav(ip,j) = rc2d(ip,jp)*1000.
      rnav(ip,j) = (rn2d(ip,jp)-rc2d(ip,jp))*1000.
         enddo
        enddo
    else
        jp=0
        do j=jts, jte
        jp = jp + 1
          ip = 0
         do i= its, ite
            ip = ip + 1
      rcav(ip,j) = max(0.,rc2d(ip,jp)*1000.-rcav(ip,j))
      rnav(ip,j) = max(0.,(rn2d(ip,jp)-rc2d(ip,jp))*1000.-rnav(ip,j))
         enddo
        enddo
    end if
!    write(6,*)'rc2d',maxval(rc2d),maxval(rn2d)
!    write(6,*)'rcav',maxval(rcav),' rnav ',maxval(rnav)
!   now call wetdep_prep
!    write(6,*)'ntra',ntra
!    call flush(6)
    if(iam.eq.iamprn)then
      do k=1,63
        write(250+iam,*)'pr3d',k,pr3d(iprn,jprn,k),prl3d(iprn,jprn,k),'tk3d',tk3d(iprn,jprn,k)
        call flush(250+iam)
      end do
      write(250+iam,*)'lbound',lbound(pr3d),'ubound',ubound(pr3d)
      call flush(250+iam)
      write(250+iam,*)'ims',ims,ime,jms,jme,'ibeg',ibeg,'jbeg',jbeg
        
    endif
    call wetdep_prep(tr3d_inout4,tk3d,pr3d,prl3d,ph3d,phl3d, &
      us3d,vs3d,ws3d, &
      pblht, &
      rcav,raincv_v, &
!      t_phy,moist,u_phy,v_phy,p_phy,chem, &
      t_phy,moist,u_phy,v_phy,p_phy, &
      grvity,rd,p1000,cp, &
      t8w,p8w,pbl,z_at_w,zmid,dz8w,vvel, &
      rho_phy, nsol, num_moist,ntra, &
      ids,ide, jds,jde, kds,kde, &
      ims,ime, jms,jme, kms,kme, &
      its,ite, jts,jte, kts,kte,rc=rc) 
!    write(6,*)'call grell nsol',nsol,'num_chem',num_chem
!    write(6,*)'shape tr3d_inout4',shape(tr3d_inout4)
!    call flush(6)
!    write(6,*)'nchemfull',nchemfull,'num_moist',num_moist
!    call flush(6)
    call initwetdep(its,ite,jts,jte,kts,kte,t_phy,p_phy)
    if(iam.eq.iamprn)then
      do k=1,63
        write(250+iam,*)'brcl in',k,tr3d_inout4(iprn,jprn,k,p_brcl),tr3d_inout4(iprn,jprn,k,p_brno3)
        call flush(250+iam)
      end do
      brclold=tr3d_inout4(iprn,jprn,:,p_brcl)
      brno3old=tr3d_inout4(iprn,jprn,:,p_brno3)
      do k=1,63
        write(250+iam,*)k,'rho',rho_phy(iprnin,k,jprnin),'u',u_phy(iprnin,k,jprnin),v_phy(iprnin,k,jprnin)
        call flush(250+iam)
      end do
      do k=1,63
        write(250+iam,*)k,'t',t_phy(iprnin,k,jprnin),moist(iprnin,k,jprnin,:)
        call flush(250+iam)
      end do
      do k=1,63
        write(250+iam,*)k,'dz',dz8w(iprnin,k,jprnin),'p',p_phy(iprnin,k,jprnin)
        call flush(250+iam)
      end do
      do k=1,64
        write(250+iam,*)k,'zatw',z_at_w(iprnin,k,jprnin)
      end do
      write(6,*)'pbl',pbl(iprnin,jprnin)
      call flush(6)
      write(250+iam,*)'call grellshape p8w',shape(p8w),'lb',lbound(p8w),'ub',ubound(p8w)
      call flush(250+iam)
      write(250+iam,*)'kms',kms,kme+1
      do k=1,64
        write(250+iam,*)'p8w before grell',k,p8w(iprnin,k,jprnin)
      end do
      call flush(250+iam)

    endif
    call grelldrvct(dt,advancecount, &
     rho_phy,raincv_v,tr3d_inout4,tr_fall, &
     u_phy,v_phy,t_phy,moist,dz8w,p_phy,p8w, &
     pbl,xlv,cp,grvity,rv,z_at_w, &
     num_chem, num_moist,nchemfull, &
     ids,ide, jds,jde, kds,kde, &
     ims,ime, jms,jme, kms,kme, &
     its,ite, jts,jte, kts,kte)
    if(iam.eq.iamprn)then
!      do k=1,63
!        write(6,*)'brcl grell',k,tr3d_inout4(iprn,jprn,k,p_brcl),tr3d_inout4(iprn,jprn,k,p_brno3)
!      end do
      do k=1,63
        if(brclold(k).ne.tr3d_inout4(iprn,jprn,k,p_brcl))then
          write(250+iam,*)'diff grell',k,tr3d_inout4(iprn,jprn,k,p_brcl)-brclold(k), &
          tr3d_inout4(iprn,jprn,k,p_brno3)-brno3old(k)
        endif
      end do
      brclold=tr3d_inout4(iprn,jprn,:,p_brcl)
      brno3old=tr3d_inout4(iprn,jprn,:,p_brno3)
    endif
     call wetdep_ls(dt,tr3d_inout4,rnav,moist,t_phy,rho_phy,var_rmv,num_moist, &
         num_chem,p_qc,dz8w,vvel,        &   
         ids,ide, jds,jde, kds,kde,                               &   
         ims,ime, jms,jme, kms,kme,                               &   
         its,ite, jts,jte, kts,kte)
    if(iam.eq.iamprn)then
!      do k=1,63
!        write(6,*)'brcl after ls',k,tr3d_inout4(iprn,jprn,k,p_brcl),tr3d_inout4(iprn,jprn,k,p_brno3)
!      end do
      do k=1,63
        if(brclold(k).ne.tr3d_inout4(iprn,jprn,k,p_brcl))then
          write(250+iam,*)'diff ls',k,tr3d_inout4(iprn,jprn,k,p_brcl)-brclold(k), &
          tr3d_inout4(iprn,jprn,k,p_brno3)-brno3old(k)
        endif
      end do
    endif
!     do i=1,nchemfull
!       if(maxval(tr_fall(:,:,ichemfullpt(i)))>0)then
!         write(6,*)'trfall',i,maxval(tr_fall(:,:,ichemfullpt(i))),chemfull(i)
!       endif
!     end do
!     do i=1,nchemfull
!        if(maxval(var_rmv(:,:,ichemfullpt(i)))>0.0)then
!          write(6,*)'var_rmv',i,maxval(var_rmv(:,:,ichemfullpt(i))),chemfull(i)
!        endif
!     end do
!     do i=1,nsol
!       write(6,*)'wetdep ',i,maxval(tr_fall(:,:,i))
!     end do
!    write(6,*)'size tr3d',size(tr3d_inout4),' chem ',size(chem)
!    call updateqsol(tr3d_inout4,chem,nsol)
    return


  end subroutine wetdep_advance
  subroutine updateqsol(tr3d_inout4,chem,nsol)
  use chem_types_mod
  use raqmschem_species_mod, only : idwetd
  
  integer i,j,k,m
  real(CHEM_KIND_R4), dimension(:,:,:,:), intent(inout) :: tr3d_inout4
  real(CHEM_KIND_R4), dimension(:,:,:,:), intent(in) :: chem
  integer isize(4),isize2(4)
  isize=shape(tr3d_inout4)
  isize2=shape(chem)
!  write(6,*)'isize',isize,' isize2 ',isize2
!  call flush(6)
  do m=1,nsol
    do k=1,isize(3)
      do j=1,isize(2)
        do i=1,isize(1)
          tr3d_inout4(i,j,k,idwetd(m))=chem(i,k,j,m)
        end do
      end do
    end do
  end do
  end subroutine updateqsol
end module wetdep_mod

