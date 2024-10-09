module raqms_plume_mod
!  use raqmschem_config_mod, only : FIRE_OPT_GBBEPx, FIRE_OPT_MODIS
!  use raqmschem_const_mod,  only : g => grvity, cp, &
!                              r_d => rd, r_v => rv, p1000mb => p1000
!  use plume_data_mod,  only : nveg_agreg
  use raq_plume_zero_mod
  use raq_plume_scalar_mod

  implicit none

  private

  public :: num_frp_plume
  public :: raqms_plumerise_driver
contains
  subroutine raqms_plumerise_driver(tk3d,rv3d,ws3d,u3d,v3d,p3d,ph3d,vgtyp, &
    num_plume_data,plume,ims,ime,jms,jme,nl )
    use chem_types_mod 
    use raqmschem_pmgrid_mod, only : iam,tile
    use raqmschem_const_mod, only : rd,rdgsdchem,grvityfms
    use raqmschem_const_mod, only : avgro,mw_co
!    use plume_data_mod, only : catb,flaming,msize,frp_mean
    use raq_plume_data_mod
    use raqmschemcomm_mod, only : bbco_d,emisco3d,ktopco,kbotco,ltb
    use raqmschemcomm_mod, only : colemisco,colemisco_chem,colemisoc
    use raqmschem_pmgrid_mod, only : iam,tile
    use raqmschemcomm_mod,only : plumerisefire_frq
    implicit none
    integer ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte,nl
    integer i,j,k,ii
    integer, intent(in) :: num_plume_data
    real(CHEM_KIND_R8),intent(in),dimension(ims:ime,jms:jme) :: vgtyp
    integer,dimension(ims:ime,jms:jme) :: ivgtyp
    real(CHEM_KIND_R8), dimension(ims:ime,jms:jme,nl),intent(in) :: tk3d,ws3d,u3d,v3d,p3d
    real(CHEM_KIND_R8), dimension(ims:ime,jms:jme,nl+1),intent(in) :: ph3d
    real(CHEM_KIND_R4), dimension(ims:ime,jms:jme,nl),intent(in) :: rv3d
    real(CHEM_KIND_R4), dimension(ims:ime,jms:jme,nl)  :: rho,vvel
    real(CHEM_KIND_R4), dimension(nl) :: u_in,v_in,w_in,qv_in,pi_in,zmid,z_lev,rho_phyin,theta_in
    real(CHEM_KIND_R4), dimension(ims:ime,jms:jme,nl+1) :: z_at_w
    real(chem_KIND_R4) :: dz,zk,rcp
    real(CHEM_KIND_R4) :: plume_frp(num_frp_plume),eburnout(nl,1)
    real(CHEM_KIND_R4), dimension(ims:ime,jms:jme,num_plume_data) :: plume
    real(CHEM_KIND_R4), parameter :: frp2plume          = 1.e+06_CHEM_KIND_R4  !  FRP-to-plume conversion factor
    real, dimension(nveg_agreg)        :: firesize,mean_fct
    real(chem_kind_r4) :: colemiscoadd(ime-ims+1,jms:jme),colemiscodiff(ime-ims+1,jms:jme)
    integer,parameter :: plumerise_flag=2
    integer iat,jat,imax(2)
    common/ij/iat,jat
    rcp=rdgsdchem/cp
!    write(6,*)'rcp',rcp,'rd',rd,'cp',cp,'p100',p1000mb
!    write(6,*)'top raqms_plumerise',ims,ime,jms,jme,nl
!    call flush(6)
    rho(:,:,:)=p3d(:,:,:)/(RDgsdchem*tk3d(:,:,:)*(1.+.608*rv3d(:,:,:)))
!   added g 3/31/2022
    vvel(:,:,:)=-ws3d(:,:,:)/rho(:,:,:)/grvityfms
    if(.not.allocated(emisco3d))then
      allocate(emisco3d(ime-ims+1,jms:jme,nl))
    endif
    colemiscoadd=0.0
    if(.not.allocated(colemisco))then

      allocate(colemisco(ime-ims+1,jms:jme),colemisco_chem(ime-ims+1,jms:jme))
      allocate(colemisoc(ime-ims+1,jms:jme))
      colemisoc=0.0
      colemisco=0.0
      colemisco_chem=0.0
    else
      colemiscodiff=colemisco_chem-colemisco
!      write(300+iam,*)'colemisco_chem',maxval(colemisco_chem),minval(colemisco_chem)
!      write(300+iam,*)'diff colemisco',maxval(colemiscodiff),minval(colemiscodiff)
!      write(300+iam,*)'colemisoc',maxval(colemisoc),minval(colemisoc)
      imax= maxloc(colemisoc)
!      write(300+iam,*)'imax oc',imax
      colemisoc=0.0
      colemisco=0.0
      colemisco_chem=0.0
    endif
    if(allocated(ktopco))then
      deallocate(ktopco)
    endif
    if(.not.allocated(ktopco))then
      allocate(ktopco(ime-ims+1,jms:jme))
    endif
    if(.not.allocated(kbotco))then
      allocate(kbotco(ime-ims+1,jms:jme))
    endif
    emisco3d=0.0
    ktopco=0
    kbotco=0

   
    do j=jms,jme
      do i=ims,ime
        ii=i-ims+1
!        if(i==145.and.j==73)then
!          write(500+iam,*)'p_phy',p3d(i,j,1),'t',tk3d(i,j,1)
!          write(500+iam,*)'ws3d',ws3d(i,j,1),'rho',rho(i,j,1),'vvel',vvel(i,j,1)
!        endif
        if(bbco_d(ii,j)==0.0)cycle
        u_in=u3d(i,j,:)
        v_in=v3d(i,j,:)
        w_in=vvel(i,j,:)
        qv_in=rv3d(i,j,:)
        pi_in=cp*(p3d(i,j,:)/p1000mb)**rcp
        z_at_w(i,j,1)=max(0.,ph3d(i,j,1)/g)
        do k=1,nl
          dz=abs(ph3d(i,j,k+1)-ph3d(i,j,k))/g
          z_at_w(i,j,k+1)=z_at_w(i,j,k)+dz
!        if(i==145.and.j==73)then
!           write(500+iam,*)'ph13d',k,ph3d(i,j,k:k+1),'g',g
!           write(500+iam,*)'dz',dz,k,'z_at_w',z_at_w(i,j,k:k+1)
!        endif
        end do
        do k=1,nl
          zk=.5*(z_at_w(i,j,k+1)+z_at_w(i,j,k))
          zmid(k)=zk-z_at_w(i,j,1)
          z_lev(k)=z_at_w(i,j,k)-z_at_w(i,j,1)
!          write(300+iam,*)k,'p3d',p3d(i,j,k),'z_lev',z_lev(k)
!          call flush(300+iam)
        end do
        rho_phyin=rho(i,j,:)
        theta_in=tk3d(i,j,:)/pi_in*cp
        ivgtyp(i,j)=nint(vgtyp(i,j))
              !if(i==ims.and.j==jms)then
!       if(i==145.and.j==73)then
!         write(500+iam,*)'i',i,j,kind(rdgsdchem)
!              do k=1,nl
!              write(500+iam,*)k,'P',real(p3d(i,j,k),4),'cp',cp,'rcp',rcp
!              write(500+iam,*)k,'T',real(tk3d(i,j,k),4),'rd',rdgsdchem
!              write(500+iam,*)k,'u',u_in(k),v_in(k),'w',w_in(k),'qv',qv_in(k)
!              write(500+iam,*)k,'z_at_w',z_at_w(i,j,k),'rho',rho_phyin(k)
!              write(500+iam,*)k,'pi',pi_in(k),'z_lev',z_lev(k)
!              write(500+iam,*)k,'theta',theta_in(k)
!              write(500+iam,*)k,'zmid',zmid(k)
!              end do
!               write(500+iam,*)'ivg',ivgtyp(i,j)
!             endif
        plume_frp(p_frp_flam_frac) = flaming(catb(ivgtyp(i,j)))
        plume_frp(p_frp_mean     ) = frp2plume * plume(i,j,1)
        plume_frp(p_frp_std      ) = 0.3_CHEM_KIND_R4   * frp2plume * plume(i,j,1)
        plume_frp(p_frp_mean_size) = msize(ivgtyp(i,j)) * frp2plume * plume(i,j,1)
        plume_frp(p_frp_std_size ) = 0.5_CHEM_KIND_R4 * plume_frp(p_frp_mean_size )
!        if(tile.eq.2.and.i.eq.83.and.j.eq.134)then
!          write(6,*)'plume1 ',plume(i,j,1),'ebrunout',eburnout(1,1),'ivg',ivgtyp(i,j)
!          write(6,*)'bbco_d',ii,j,bbco_d(ii,j)
!        endif
        if(plume(i,j,1)/=0.0)then
          iat=i
          jat=j
!          write(500+iam,*)'firesize',firesize,'mean_fct',mean_fct
!          write(500+iam,*)'plume_frp',plume_frp,plumerise_flag
          call plumerise(nl,firesize,mean_fct,1,bbco_d(ii,j),eburnout,u_in,v_in,w_in,theta_in,pi_in, &
          rho_phyin,qv_in,zmid,z_lev,plume_frp,plumerise_flag)
          if(eburnout(1,1)/=0.0)then
            emisco3d(ii,j,ltb(1))=eburnout(1,1)
            colemiscoadd(ii,j)=colemiscoadd(ii,j)+emisco3d(ii,j,ltb(1))
!            write(6,*)i,j,'bbco-d',bbco_d(ii,j),'eburnout(1,1)',eburnout(1,1),'diff',bbco_d(ii,j)-eburnout(1,1)
            do k=2,nl
              if(eburnout(k,1)/=0.0)then
!               ajl need to convert back to kg/m2 to pass to chem_driver like
!               gsdchem does 
!                emisco3d(ii,j,ltb(k))=eburnout(k,1)
!                if(k==2)then
                  !write(6,*)'error burn emis at k=2'
!                  call flush(6)
!                endif
                emisco3d(ii,j,ltb(k))=eburnout(k,1)*(z_at_w(i,j,k+1)-z_at_w(i,j,k))
                colemiscoadd(ii,j)=colemiscoadd(ii,j)+emisco3d(ii,j,ltb(k))
!                if(emisco3d(ii,j,ltb(k))<0.0)then
!                  write(6,*)'raq neg dz',(z_at_w(i,j,k+1)-z_at_w(i,j,k)),'eburnout',eburnout(k,1)
!                  write(6,*)i,j,'z_at_w',k+1,z_at_w(i,j,k+1),k,z_at_w(i,j,k)
!                endif

!                write(6,*)i,j,k,'emisco3d',emisco3d(ii,j,ltb(k)),'dz',z_at_w(i,k+1,j)-z_at_w(i,k,j)
!                write(500+iam,'("ratio",3i3,f10.6)')i,j,k,emisco3d(ii,j,ltb(k))/emisco3d(ii,j,ltb(1))
                ktopco(ii,j)=ltb(k) ! make for raqms which has 1 at top for chem
                if(kbotco(ii,j)==0)then
                  kbotco(ii,j)=ltb(k)
                endif
              endif
!              if(tile.eq.2.and.i.eq.83.and.j.eq.134)then
!                write(6,*)'emisco3d',ii,j,ltb(k),eburnout(k,1),'ktopco',ktopco(ii,j)
!              endif
            enddo
          endif
        endif
      end do
    end do
!   write(6,*)'z_at_w',kts,minval(z_at_w(:,kts,:))
!    write(6,*)'z_at_w kts+1',minval(z_at_w(:,kts+1,:)),maxval(z_at_w(:,kts+1,:))
!    write(6,*)'dz8w',minval(dz8w(:,kts,:)),maxval(dz8w(:,kts,:))
!   here convert Kg/m3/sec to molec/cm3/sec which chemdriver needs multiplied by
!   dz
!   here convert Kg/m2/sec to molec/cm3/sec which chemdriver needs
!   change per Brad 3/28/2022
    emisco3d=emisco3d*.1*avgro/mw_co ! make molec/ccm2 since per moles/km2/sec layer
    colemiscoadd=colemiscoadd*.1*avgro/mw_co
    colemisco=colemisco+colemiscoadd*60*plumerisefire_frq
!    colemisco=colemisco*.1*avgro/mw_co
    if(maxval(abs(colemisco))>1.e40)then
      write(300+iam,*)'inf colemsico'
        !ii=i-ims+1
      do j=jms,jme
      do ii=1,i-ims+1
        if(abs(colemisco(ii,j))>1.e40)then
          write(300+iam,*)'colemisco',ii,colemisco(ii,j)
        endif
      end do
      end do
    endif
!    if(maxval(abs(colemisco))>0.0)then
!    write(300+iam,*)'colemisco hr',maxval(colemisco),minval(colemisco)
!    call flush(300+iam)
!  endif
!    write(6,*)'emisco3d',maxval(emisco3d),minval(emisco3d)
!    write(6,*)'avgro',avgro,'mw_co',mw_co
    call flush(6)
    return
  end subroutine raqms_plumerise_driver
end module raqms_plume_mod
