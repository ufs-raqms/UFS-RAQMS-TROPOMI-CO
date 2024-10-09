module dep_wet_ls_mod


  implicit none

  private

  public :: wetdep_ls

contains

subroutine wetdep_ls(dt,var,rain,moist,t,rho,var_rmv,num_moist, &
         num_chem,p_qc,dz8w,vvel,             &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  use raqmschem_species_mod, only : chemsoluablefull,ichemfullpt,nsol, &
     chempointsoluablefull,lraqmschem,chemname,chemsoluable,cheminputlist
  use wetdep_alpha_mod, only : solalpha
  use wetdep_fam_mod, only :wetdeplsfam
  use raqmschem_pmgrid_mod, only : iam,iprnin,jprnin,tile,iamprn
  IMPLICIT NONE
!  ids,ide,jds,jde not used by ls
   INTEGER,      INTENT(IN   ) :: num_chem,num_moist,p_qc,                          &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   real, INTENT(IN ) :: dt
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) :: rho,dz8w,vvel,t        
   REAL,  DIMENSION( ims:ime , jms:jme , kms:kme ,1:num_chem),                        &
          INTENT(INOUT) :: var        
!   REAL,  DIMENSION( jms:jme ),                                  &
   REAL,  DIMENSION( ims:ime,jms:jme ),                                  &
          INTENT(IN   ) :: rain
   REAL,  DIMENSION( ims:ime ,  jms:jme,num_chem ),                                  &
          INTENT(INOUT   ) :: var_rmv
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: var_sum
   REAL,  DIMENSION( its:ite ,  kts:kte, jts:jte ) :: var_rmvl
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: frc,var_sum_clw,rain_clw     
!    real :: dvar,factor,clsum,alpha,rho_water
    real :: dvar,factor,clsum,rho_water
   real alpha(kts:kte),alphasoluable(its:ite,kts:kte,nsol)
   real qcld(its:ite,kts:kte),t2d(its:ite,kts:kte)
   integer :: nv,i,j,k,km,kb,kbeg,nmask,kbot(its:ite),ktop(its:ite)
   logical maskalpha(its:ite)
!  in ls all three indexes are the same and full limits dont have to start at 1
!   write(400+iam,*)'ls its',its,ite,' jts ',jts,jte
!   write(400+iam,*)'ls ims',ims,ime,' jms ',jms,jme
   alphasoluable=.5 ! ajl for now
    rho_water = 1000.
    var_rmv (:,:,:)=0.
!   write(6,*) 'in wetdepls, p_qc = ',p_qc
!   nv=p_bc1
!    write(6,*)'num_chem ls',num_chem
    do nv=1,num_chem
!      if(lraqmschem(nv))then
!     only do for chemicals that are soluable
      if(.not.chemsoluablefull(nv))cycle
!      if(iam.eq.0)then
!         write(6,*)'nv',nv,'soluable'
!         call flush(6)
!      endif
!            if(chemsoluable(nv))then
!               write(6,*)nv,chemname(nv),maxval(var(:,:,:,nv)),minval(var(:,:,:,nv))
!               call flush(6)
!            endif
        
!        if(tile.eq.2)then
!          write(6,*)'ls its',its,ite,jts,jte,'nv',nv,cheminputlist(nv)
          !write(6,*)'iprnin',iprnin,jprnin
!          call flush(6)
!        endif
        do j=jts,jte
        do i=its,ite
          var_sum_clw(i,j)=0.
          var_sum(i,j)=0.
          var_rmvl(i,:,j)=0.
          frc(i,j)=0.
          rain_clw(i,j)=0.
!          if(rain(j).gt.1.e-3)then
          if(rain(i,j).gt.1.e-3)then
! convert rain back to rate
!
!             rain_clw(i,j)=rain(j)/dt
             rain_clw(i,j)=rain(i,j)/dt
! total cloud water
!
             do k=1,kte-1
                dvar=max(0.,moist(i,k,j,p_qc)*rho(i,k,j)*vvel(i,k,j)*dz8w(i,k,j))
                var_sum_clw(i,j)=var_sum_clw(i,j)+dvar
                var_sum(i,j)=var_sum(i,j)+var(i,j,k,nv)*rho(i,k,j)
             enddo
!            ajl way too big for gases and not aerosols
             if(var_sum(i,j).gt.1.e-20 .and. var_sum_clw(i,j).gt.1.e-5 ) then
!             assuming that frc is onstant, it is my conversion factor 
!            (just like in convec. parameterization
                frc(i,j)=rain_clw(i,j)/var_sum_clw(i,j)
!         print *,'frc ', frc(i,j),var_sum_clw(i,j),var_sum(i,j)
                frc(i,j)=max(1.e-6,min(frc(i,j),.004))
             endif
          endif
        enddo
        enddo
!
! get rid of it
!
        kbot=kts
        ktop=kte
!        write(6,*)'var_sum',maxval(var_sum),'var_sum_clw',maxval(var_sum_clw),'qc', &
!        maxval(moist(:,:,:,p_qc))
        call flush(6)
        do j=jts,jte
          maskalpha=0
          nmask=0
          do i=its,ite
            if(rain(i,j).gt.1.e-3)then
             maskalpha(i)=.true.
             nmask=nmask+1
             endif
          end do
          if(nmask>0)then
!            write(6,*)'nmask large ',nmask
!            call flush(6)
            do k=kts,kte
              do i=its,ite
                qcld(i,k)=moist(i,k,j,p_qc)
                t2d(i,k)=t(i,k,j)
              end do
            end do
            call solalpha(its,ite,kts,kte,maskalpha,kbot,ktop,j,t2d,qcld,alphasoluable)
          endif
          do i=its,ite
#ifdef DIAGWET
            if(iam.eq.iamprn)then
              if(i.eq.iprnin.and.j.eq.jprnin)then
                write(6,*)'rain',i,j,rain(i,j)
                write(300+iam,*)'rain',i,j,rain(i,j)
                call flush(6)
              endif
           endif
#endif
           if(rain(i,j).gt.1.e-3)then
!           if(var_sum(i,j).gt.1.e-6 .and. var_sum_clw(i,j).gt.1.e-5)then
#ifdef DIAGWET
            if(iam.eq.iamprn)then
              if(i.eq.iprnin.and.j.eq.jprnin)then
                write(6,*)'var_sum',var_sum(i,j),var_sum_clw(i,j)
                write(300+iam,*)'var_sum',var_sum(i,j),var_sum_clw(i,j)
                call flush(6)
              endif
           endif
#endif
           if(var_sum(i,j).gt.1.e-20 .and. var_sum_clw(i,j).gt.1.e-5)then
!            if(chemsoluablefull(nv))then
              do k=kts,kte-2
                alpha(k)=alphasoluable(i,k,chempointsoluablefull(nv))
#ifdef DIAGWET
            if(iam.eq.iamprn)then
              if(i.eq.iprnin.and.j.eq.jprnin)then
                  write(6,*)'alpha',k,alpha(k),nv,chempointsoluablefull(nv)
                  write(300+iam,*)'alpha',k,alpha(k),nv,chempointsoluablefull(nv)
                  call flush(6)
              endif
            endif
#endif
!                if(tile.eq.2)then
#ifdef DIAGWET
                if(i.eq.iprnin.and.j.eq.jprnin.and.tile.eq.2)then
                  write(300+iam,*)'chempointsol',i,j,'nv',nv,chempointsoluablefull(nv),'k',k,alpha(k)
                  call flush(300+iam)
                endif
#endif
!                if(iam.eq.0)then
!                write(6,*)'alpha nv',nv,'k',k,alpha(k),'point',chempointsoluablefull(nv)
!                endif

!                if(isnan(alpha(k)))then
!                   write(6,*)'nan alpha ls ',i,j,k,nv
!                   call flush(6)
!                endif

              end do
!            else
!              alpha=0.
!            endif
            do k=kts,kte-2
!              if(var(i,j,k,nv).gt.1.e-08 .and. moist(i,k,j,p_qc).gt.1.e-8)then
              if(var(i,j,k,nv).gt.1.e-19 .and. moist(i,k,j,p_qc).gt.1.e-8)then
                factor = max(0.,frc(i,j)*rho(i,k,j)*dz8w(i,k,j)*vvel(i,k,j))
                dvar=max(0.,alpha(k)*factor/(1+factor)*var(i,j,k,nv))
#ifdef DIAGWET
                if(i.eq.iprnin.and.j.eq.jprnin.and.tile.eq.2)then
!                if(tile.eq.2)then
                  write(300+iam,*)i,j,'dvar',k,dvar,'alpha',alpha(k),'factor ',factor,'var',nv,var(i,j,k,nv)
                  call flush(300+iam)
                endif
            if(iam.eq.iamprn)then
              if(i.eq.iprnin.and.j.eq.jprnin)then
                  write(6,*)'dvar',dvar,'factor',factor
                  write(300+iam,*)'dvar',dvar,'factor',factor
                  call flush(6)
              endif
            endif
#endif
!                if(isnan(dvar))then
!                   write(6,*)'dvar nan',i,j,k,nv
!                   call flush(6)
!                endif
                dvar=min(dvar,var(i,j,k,nv))
                if((var(i,j,k,nv)-dvar).lt.1.e-24)then
                  dvar=var(i,j,k,nv)-1.e-24
                  var(i,j,k,nv)=var(i,j,k,nv)-dvar
                else
                  var(i,j,k,nv)=var(i,j,k,nv)-dvar
                endif
                var_rmvl(i,k,j)=dvar
!                if(isnan(var(i,j,k,nv)))then
!                  write(6,*)'nan ls var',i,j,k,nv
!                  call flush(6)
!                endif
              
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
#ifdef DIAGWET
      if(iam.eq.iamprn)then
        write(6,*)'var_rmv',nv,var_rmvl(iprnin,1,jprnin)
      endif
#endif
      call wetdeplsfam(nv,num_chem,its,ite,jts,jte,kts,kte,var,var_rmvl)
    enddo
END SUBROUTINE WETDEP_LS

end module dep_wet_ls_mod
