module wetdep_fam_mod
#include <options.h>
  use raqmschem_species_mod
  implicit none
! variables for wetdep families
  integer,parameter :: nfamily=4
  integer,parameter :: nsubfam=2
  integer :: pointsubls(nfamily,nsubfam),nsub(nfamily),pointfam(nfamily)
  integer :: pointsubcum(nfamily,nsubfam),pointfamcum(nfamily)
  integer :: ifamls(30+nchemfulldim,2)
  integer,parameter :: nvar=4
  integer :: pointfamls(nvar,2),mfamls(nvar,2)
  integer ::ipointvar(nvar)
  logical :: isfamvar(30+nchemfulldim)
  integer :: varnum(30+nchemfulldim)
  data nsub/2,2,1,1/
  character*8 famname(nfamily)
  data famname/'ox','noy','cly','bry'/
contains
  subroutine initsubpointers
  use raqmschem_species_mod
  use raqmschem_pmgrid_mod, only : iam,iprnin,jprnin,iamprn
  implicit none
  integer nf,ns,ivar,n,nn
  ifamls=-1
  isfamvar=.false.
  varnum=0
  ipointvar(1)=p_hno3t
  ipointvar(2)=p_hno4
  ipointvar(3)=p_hcl
  ipointvar(4)=p_hbr
  do ivar=1,nvar
    isfamvar(ipointvar(ivar))=.true.
    varnum(ipointvar(ivar))=ivar
  end do
  pointfam(1)=p_ox
  pointfam(2)=p_noy
  pointfam(3)=p_cly
  pointfam(4)=p_bry 
  pointfamls(1,1)=p_ox
  pointfamls(1,2)=p_noy
  pointfamls(2,1)=p_ox
  pointfamls(2,2)=p_noy
  pointfamls(3,1)=p_cly
  pointfamls(4,1)=p_bry
! ajl added like 9/21/2018 to put into correct family
  mfamls(1,1)=1
  mfamls(1,2)=2
  mfamls(2,1)=1
  mfamls(2,2)=2
  mfamls(3,1)=3
  mfamls(4,1)=4
  pointsubls(1,1)=p_hno3t
  pointsubls(1,2)=p_hno4
  pointsubls(2,1)=p_hno3t
  pointsubls(2,2)=p_hno4
  pointsubls(3,1)=p_hcl
  pointsubls(4,1)=p_hbr
  do nf=1,nfamily
    if(pointfam(nf)>0)then
!      pointfamcum(nf)=chempointcum(pointfam(nf))
      pointfamcum(nf)=pointfam(nf)
    endif
    do ns=1,nsub(nf)
      if(pointsubls(nf,ns)>0)then
!      pointsubcum(nf,ns)=chempointcum(pointsubls(nf,ns))
        pointsubcum(nf,ns)=chempointsoluablefull(pointsubls(nf,ns))
!       ifamls new like 9/21/2018
        ifamls(pointsubcum(nf,ns),ns)=nf
#ifdef DIAGWETDEP
      if(iam.eq.iamprn)then
        write(6,*)'nf ',nf,' ns ',ns,' pointsubls ',pointsubls(nf,ns)
        write(6,*)'ifamls pointsubcum ',pointsubcum(nf,ns),ns,' )= ',nf
        write(6,*)'mfamls ',mfamls(nf,ns)
      endif
      if(iam.eq.0)then
        write(6,*)'nf',nf,'ns',ns,'ifamls( ',pointsubcum(nf,ns),' )=',nf
        call flush(6)
      endif
#endif
!      pointsubcum(nf,ns)=pointsubls(nf,ns)
#ifdef DIAGWETDEP
       if(iam.eq.0)then
       do n=1,nsol
         nn=idwetd(n)
         if(pointsubls(nf,ns).eq.nn)then
           write(6,*)'write nn',nn,'pointsubcum','n',n,pointsubcum(nf,ns),'nf',nf,ns
         endif
        end do
        endif
#endif
      endif
    end do
  end do
#ifdef DIAGWETDEP
  if(iam.eq.0)then
    do nf=1,nfamily
      write(6,*)'pointfamcum',nf,pointfamcum(nf)
      do ns=1,nsub(nf)
        write(6,*)'family ',nf,' pointfam ',pointfamls(nf,ns)
        write(6,*)'pointlssub nf',nf,' ns ',ns,pointsubls(nf,ns),' pointcumsub ',pointsubcum(nf,ns)
      end do
    end do
    call flush(6)
  endif
#endif
#ifdef DIAGWETDEP
  if(iam.eq.iamprn)then
    do nf=1,nfamily
      write(6,*)'pointfamcum',nf,pointfamcum(nf)
      do ns=1,nsub(nf)
        write(6,*)'family ',nf,' pointfam ',pointfamls(nf,ns)
        write(6,*)'pointlssub nf',nf,' ns ',ns,pointsubls(nf,ns),' pointcumsub ',pointsubcum(nf,ns)
      end do
    end do
    call flush(6)
  endif
#endif
  end subroutine initsubpointers
  subroutine wetdeplsfam(nv,num_chem,its,ite,jts,jte,kts,kte,var,var_rmvl)
  use raqmschem_pmgrid_mod, only : iam,tile,ibeg,jbeg,iprnin,jprnin,iprn,jprn,kprnin
  use raqmschem_species_mod, only : p_brcl,p_bry
  use raqmschem_pmgrid_mod, only : iamprn
#ifdef WETDEPDIAG
  use raqmschemcomm_mod, only : wetc,wetl,wetcf,wetlf,nwetdepdiag,wetccol,wetlcol,dpmgrd
#endif
  implicit none
  integer nv,num_chem,its,ite,jts,jte,kts,kte,ivar
  real,intent(inout) :: var(its:ite,jts:jte,kts:kte,num_chem)
  real, intent(in) :: var_rmvl(its:ite,kts:kte,jts:jte) ! is amount removed so need to subtract
  integer nf,ns,j,k,jj,jat,iat,m
!  write(300+iam,*)'its lsfam',its,ite,jts,jte,kts,kte,'nv',nv,'tile',tile
!  call flush(300+iam)
!  iat=iprnin-ibeg+1
   iat=iprnin
!  write(400+iam,*)'iat',iat,'iprnin',iprnin,'ibeg',ibeg,'jts',jts,jte,'jprnin',jprnin
!  write(400+iam,*)'nv',nv,'isfamvar',isfamvar(nv)
!  call flush(400+iam)
!  all indices are full
!  write(6,*)'top wetdeplsfam ',nv
!  call flush(6)
!  write(300+iam,*)'wetdeplsfam',nv,jts,jte
!  call flush(300+iam)
  if(isfamvar(nv))then
    ivar=varnum(nv) ! point to soluable species list which map to families 1-4 hno3t,hno4,hcl,hbr
!   handle nsub families that this variable maps into
    do nf=1,nsub(ivar)

!      write(6,*)'ls var_rmvl',nv,nf,maxval(var_rmvl),minval(var_rmvl)
!      call flush(6)
      do k=kts,kte
        do j=jts,jte
!         if(pointfamls(ifam,nf).eq.p_bry)then
!         if(iat>=its.and.iat<=ite.and.tile.eq.2.and.j.eq.jprnin.and.k.eq.1)then
!           write(300+iam,*)'nv',nv,'ls ifam',ifam,'nf',nf,var(iat,j,k,pointfamls(ifam,nf))
!           write(300+iam,*)'var_rmvl ',var_rmvl(iat,k,j)
!           call flush(300+iam)
!         endif
!         endif
!         fill in family values in var list (inout4)
          var(:,j,k,pointfamls(ivar,nf))=var(:,j,k,pointfamls(ivar,nf))-var_rmvl(:,k,j)
!          var(:,j,k,nv)=var(:,j,k,nv)-var_rmvl(:,k,j)
#ifdef WETDEPDIAG
#ifdef DIAGWETDEP
          if(j.eq.jprnin.and.k.eq.kprnin.and.iam.eq.iamprn.and.var_rmvl(iprnin,k,j).ne.0.0)then
            write(6,*)'nv ',nv,' ifamls ',ifamls(ivar,nf),' nf ',nf,'ivar',ivar
            write(6,*)' pointfamls ',pointfamls(ivar,nf),'var_rmvl ',var_rmvl(iprnin,k,j)
            write(6,*)'mfamls',mfamls(ivar,nf)
            call flush(6)
          endif
#endif
!         there are 4 families add into correct one
!         new like 9/21/2018
          wetlf(:,j,k,mfamls(ivar,nf))=wetlf(:,j,k,mfamls(ivar,nf))-var_rmvl(:,k,j)
          !wetlf(:,j,k,ifam)=wetlf(:,j,k,ifam)-var_rmvl(:,k,j)
#endif
        end do
      end do
    end do
  endif
#ifdef WETDEPDIAG
#ifdef DIAGWETDEP
  if(iam.eq.iamprn)then
    write(6,*)'nv',nv,'chemsoluable ',chemsoluable(nv),' chemsoluablefull ',chemsoluablefull(nv)
  endif
#endif
!  if(chemsoluable(nv))then
  if(chemsoluablefull(nv))then
!    m=chempointsoluable(nv)
!   point to soluable species list 1-14
    m=chempointsoluablefull(nv)
#ifdef DIAGWETDEP
    if(iam.eq.iamprn)then
    write(6,*)j,'m at',m,'nv',nv,'chemsoluablefull',chemsoluablefull(nv),'chemsoluable',chemsoluable(nv)
    write(6,*)'chempointsoluable',chempointsoluable(nv)
    call flush(6)
    endif
#endif
    if(m>14)then
      write(6,*)'error m',m
      call flush(6)
    endif
    do k=kts,kte
      do j=jts,jte
#ifdef DIAGWETDEP
          if(j.eq.jprnin.and.k.eq.kprnin.and.iam.eq.iamprn)then
            write(6,*)'m ',m,' nv ',nv,' wetl ',wetl(iprnin,j,k,m),' var_rmvl ',var_rmvl(iprnin,k,j)
            call flush(6)
          endif
#endif
        wetl(:,j,k,m)=wetl(:,j,k,m)-var_rmvl(:,k,j)
!       make it tendency due rainout
        wetlcol(:,j,m)=wetlcol(:,j,m)-var_rmvl(:,k,j)*dpmgrd(:,j,k)/9.80616 !  want mb
#ifdef DIAGWETDEP
          if(j.eq.jprnin.and.k.eq.kprnin.and.iam.eq.iamprn)then
              write(6,*)'ajl wetl ',m,' nv ',nv,wetl(iprnin,j,k,m)
              call flush(6)
          endif
#endif
      end do
    end do
!    if(minval(var_rmvl)<0.0)then
!      write(200+iam,*)'neg var_rmvl ',m,minval(var_rmvl)
!    endif
  endif
#endif
!  write(300+iam,*)'btottom wetdeplsfam',nv
!  call flush(300+iam)
!  write(6,*)'bottom wetdeplsfam'
!  call flush(6)
  return
  end subroutine wetdeplsfam
!  subroutine wetdepcumfam(num_chem,its,ite,jts,jte,kts,kte,dt,chem,oracert,j)
  subroutine wetdepcumfam(num_chem,its,ite,jts,jte,kts,kte,dt,chem,wetdepsol,j)
  use chem_const_mod, only : epsilc
  use raqmschem_pmgrid_mod,only : iam,tile,ibeg,jbeg,iprnin,jprnin,tileprn,kprnin,iamprn,tile
  use raqmschem_species_mod
#ifdef WETDEPDIAG
  use raqmschemcomm_mod, only : wetc,wetl,wetcf,wetlf,nwetdepdiag,wetccol,dpmgrd
#endif
  implicit none
  integer nv,num_chem,its,ite,jts,jte,kts,kte,j,k
  integer nf,ns,iat,jat,m
  real,intent(in) :: dt
  real, intent(inout) :: chem(its:ite,jts:jte,kts:kte,num_chem)
!  real, intent(in) :: tracert(its:ite,kts:kte,num_chem) ! is tendency of tracer
  real, intent(in) :: wetdepsol(its:ite,kts:kte,nsol) ! is tendency of tracer
  integer i,ii
!  write(300+iam,*)'cum',its,ite,jts,jte,kts,kte,'j',j,'tile',tile
!  call flush(300+iam)
  iat=iprnin
!  jat=j-jts+jbeg
!  write(400+iam,*)'cumfam its',its,ite,'jts',jts,jte,'ibeg',ibeg,jbeg
!  write(400+iam,*)'ls iat',iat,'jat',jat,'iprnin',iprnin,'jprnin',jprnin ,'nfamily',nfamily
!  call flush(400+iam)
#ifdef WETDEPDIAG
!  write(400+iam,*)'cumfam wetdepsol ',shape(wetdepsol)
!  write(400+iam,*)'wetc',shape(wetc),'wetcf',shape(wetcf)
!  call flush(400+iam)
  if(j.eq.jts)then
    nwetdepdiag=nwetdepdiag+1
  endif 
!  write(300+iam,*)'lb dmpgrd ',lbound(dpmgrd),' ub ',ubound(dpmgrd)
!  write(300+iam,*)'lb wetccol ',lbound(wetccol),' ub ',ubound(wetccol)
!  write(300+iam,*)'lb wetdepsol',lbound(wetdepsol),' ub ',ubound(wetdepsol)
!  call flush(300+iam)
  do m=1,nsol
    do k=1,kts,kte
      do ii=its,ite
        i=ii-its+1
        if(m.eq.13.and.ii.eq.iprnin.and.j.eq.jprnin.and.iam.eq.iamprn)then
        if(wetdepsol(ii,k,m).ne.0.0)then
          write(300+iam,*)'wetccol ',m,wetccol(i,j,m),' wetdepsol ii ',ii,wetdepsol(ii,k,m),'dpm',i,j,k,dpmgrd(i,j,k)
          write(300+iam,*)'nwetdepdiag',nwetdepdiag,soluablechemname(m)
          call flush(300+iam)
        endif
        endif
        wetccol(i,j,m)=wetccol(i,j,m)+wetdepsol(ii,k,m)*dpmgrd(i,j,k)/9.80616 !  want mb
      end do
    end do
    wetc(:,j,:,m)=wetc(:,j,:,m)+wetdepsol(:,:,m) ! is a tendency alreaDY
  end do
!  write(300+iam,*)'did loop'
!  call flush(300+iam)
#endif
  do nf=1,nfamily
    do ns=1,nsub(nf)
!     if(pointfamcum(nf).eq.p_bry)then
!     if(iat>=its.and.iat<=ite.and.j.eq.jprnin.and.tile.eq.2)then
!       write(300+iam,*)'nf',nf,'ns',ns,'cum add to ',pointfamcum(nf),chem(iat,j,1,pointfamcum(nf))
!       write(300+iam,*)'add ',pointsubcum(nf,ns),tracert(iat,1,pointsubcum(nf,ns))
!       call flush(300+iam)
!     endif
!     endif
!      write(6,*)'tracert',nf,ns,maxval(tracert(:,:,pointsubcum(nf,ns))),minval(tracert(:,:,pointsubcum(nf,ns)))
!      call flush(6)
!      if(j.eq.jprnin)then
!        write(300+iam,*)'pointfamcum',nf,ns,pointfamcum(nf),chem(iprnin,j,1,pointfamcum(nf))
!        call flush(300+iam)
!      endif
!      chem(:,j,:,pointfamcum(nf))=max(epsilc,chem(:,j,:,pointfamcum(nf))+tracert(:,:,pointsubcum(nf,ns))*dt)
#ifdef DIAGWETDEP
     if(iat>=its.and.iat<=ite.and.j.eq.jprnin.and.tile.eq.tileprn)then
!       write(500+iam,*)'chem nf ',nf,' ns ',ns,pointfamcum(nf),' msol ',pointsubcum(nf,ns),'dt',dt
!       call flush(500+iam)
!       write(500+iam,*)'iat',iat,'j',j
!       do k=1,10
!       write(500+iam,*)'chem fam to add to',k,chem(iat,j,k,pointfamcum(nf)),wetdepsol(iat,k,pointsubcum(nf,ns))
!       end do
!       call flush(500+iam)
        write(6,*)'nf ',nf,' ns ',ns,' pointfamcum ',pointfamcum(nf)
        write(6,*)'pointsubcum ',pointsubcum(nf,ns)
        write(6,*)'chem before ',chem(iprnin,j,kprnin,pointfamcum(nf)),' add ', &
          wetdepsol(iprnin,kprnin, pointsubcum(nf,ns) )*dt
        call flush(6)
      endif
#endif
      if(pointfamcum(nf)<1.or.pointfamcum(nf)>num_chem)then
         write(6,*)'ajl error pointjamcum',nf,pointfamcum(nf)
         call flush(6)
      endif
      if(pointsubcum(nf,ns)<1.or.pointsubcum(nf,ns)>nsol)then
        write(6,*)'ajl error pointsubcum',nf,ns,pointsubcum(nf,ns)
        call flush(6)
      endif
      chem(:,j,:,pointfamcum(nf))=max(epsilc,chem(:,j,:,pointfamcum(nf))+wetdepsol(:,:,pointsubcum(nf,ns))*dt)
#ifdef WETDEPDIAG
       wetcf(:,j,:,nf)=wetcf(:,j,:,nf)+wetdepsol(:,:,pointsubcum(nf,ns))
#endif
!     if(iat>=its.and.iat<=ite.and.j.eq.jprnin.and.tile.eq.tileprn)then
!       do k=1,10
!       write(500+iam,*)'chem fam after',k,chem(iat,j,k,pointfamcum(nf))
!       end do
       !call flush(500+iam)
!      endif
      
    end do
  end do
!  write(500+iam,*)'did wetdepcumfam'
!  call flush(500+iam)
!  write(300+iam,*)'did wetdepcumfam'
!  call flush(300+iam)
  return
  end subroutine wetdepcumfam

end module wetdep_fam_mod
