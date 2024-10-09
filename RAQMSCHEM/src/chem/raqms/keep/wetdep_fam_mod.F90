module wetdep_fam_mod
  use raqmschem_species_mod
  implicit none
! variables for wetdep families
  integer,parameter :: nfamily=4
  integer,parameter :: nsubfam=2
  integer :: pointsubls(nfamily,nsubfam),nsub(nfamily),pointfam(nfamily)
  integer :: pointsubcum(nfamily,nsubfam),pointfamcum(nfamily)
  integer,parameter :: nvar=4
  integer :: pointfamls(nvar,2)
  integer ::ipointvar(nvar)
  logical :: isfamvar(30+nchemfull)
  integer :: varnum(30+nchemfull)
  data nsub/2,2,1,1/
contains
  subroutine initsubpointers
  use raqmschem_species_mod
  use raqmschem_pmgrid_mod, only : iam,iprnin,jprnin
  implicit none
  integer nf,ns,ivar
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
  pointsubls(1,1)=p_hno3t
  pointsubls(1,2)=p_hno4
  pointsubls(2,1)=p_hno3t
  pointsubls(2,2)=p_hno4
  pointsubls(3,1)=p_hcl
  pointsubls(4,1)=p_hbr
  do nf=1,nfamily
    if(pointfam(nf)>0)then
      pointfamcum(nf)=chempointcum(pointfam(nf))
    endif
    do ns=1,nsub(nf)
      if(pointsubls(nf,ns)>0)then
      pointsubcum(nf,ns)=chempointcum(pointsubls(nf,ns))
      endif
    end do
  end do
  if(iam.eq.0)then
    do nf=1,nfamily
      do ns=1,nsub(nf)
        write(6,*)'family ',nf,' pointfam ',pointfamls(nf,ns)
        write(6,*)'pointlssub',ns,pointsubls(nf,ns),' pointcumsub ',pointsubcum(nf,ns)
      end do
    end do
    call flush(6)
  endif
  end subroutine initsubpointers
  subroutine wetdeplsfam(nv,num_chem,its,ite,jts,jte,kts,kte,var,var_rmvl)
  use raqmschem_pmgrid_mod, only : iam,tile,ibeg,jbeg,iprnin,jprnin
  use raqmschem_species_mod, only : p_brcl,p_bry
  implicit none
  integer nv,num_chem,its,ite,jts,jte,kts,kte,ivar
  real,intent(inout) :: var(its:ite,jts:jte,kts:kte,num_chem)
  real, intent(in) :: var_rmvl(its:ite,kts:kte,jts:jte)
  integer nf,ns,j,k,jj,jat,iat
!  write(300+iam,*)'its lsfam',its,ite,jts,jte,kts,kte,'nv',nv,'tile',tile
!  call flush(300+iam)
!  iat=iprnin-ibeg+1
   iat=iprnin
!  write(400+iam,*)'iat',iat,'iprnin',iprnin,'ibeg',ibeg,'jts',jts,jte,'jprnin',jprnin
!  write(400+iam,*)'nv',nv,'isfamvar',isfamvar(nv)
!  call flush(400+iam)
!  all indices are full
  if(isfamvar(nv))then
    ivar=varnum(nv)
    do nf=1,nsub(ivar)
      do k=kts,kte
        do j=jts,jte
          if(pointfamls(ivar,nf).eq.p_bry)then
          if(iat>=its.and.iat<=ite.and.tile.eq.2.and.j.eq.jprnin.and.k.eq.1)then
            write(300+iam,*)'nv',nv,'ls ivar',ivar,'nf',nf,var(iat,j,k,pointfamls(ivar,nf))
            write(300+iam,*)'var_rmvl ',var_rmvl(iat,k,j)
            call flush(300+iam)
          endif
          endif
          var(:,j,k,pointfamls(ivar,nf))=var(:,j,k,pointfamls(ivar,nf))+var_rmvl(:,k,j)
        end do
      end do
    end do
  endif
  return
  end subroutine wetdeplsfam
  subroutine wetdepcumfam(num_chem,its,ite,jts,jte,kts,kte,dt,chem,tracert,j)
  use chem_const_mod, only : epsilc
  use raqmschem_pmgrid_mod,only : iam,tile,ibeg,jbeg,iprnin,jprnin
  use raqmschem_species_mod
  implicit none
  integer nv,num_chem,its,ite,jts,jte,kts,kte,j
  integer nf,ns,iat,jat
  real,intent(in) :: dt
  real, intent(inout) :: chem(its:ite,jts:jte,kts:kte,num_chem)
  real, intent(in) :: tracert(its:ite,kts:kte,num_chem)
!  write(300+iam,*)'cum',its,ite,jts,jte,kts,kte,'j',j,'tile',tile
!  call flush(300+iam)
  iat=iprnin
!  jat=j-jts+jbeg
!  write(400+iam,*)'cumfam its',its,ite,'jts',jts,jte,'ibeg',ibeg,jbeg
!  write(400+iam,*)'ls iat',iat,'jat',jat,'iprnin',iprnin,'jprnin',jprnin ,'nfamily',nfamily
!  call flush(400+iam)
  do nf=1,nfamily
    do ns=1,nsub(nf)
      if(pointfamcum(nf).eq.p_bry)then
      if(iat>=its.and.iat<=ite.and.j.eq.jprnin.and.tile.eq.2)then
        write(300+iam,*)'nf',nf,'ns',ns,'cum add to ',pointfamcum(nf),chem(iat,j,1,pointfamcum(nf))
        write(300+iam,*)'add ',pointsubcum(nf,ns),tracert(iat,1,pointsubcum(nf,ns))
        call flush(300+iam)
      endif
      endif
      chem(:,j,:,pointfamcum(nf))=max(epsilc,chem(:,j,:,pointfamcum(nf))+tracert(:,:,pointsubcum(nf,ns))*dt)
    end do
  end do
  return
  end subroutine wetdepcumfam

end module wetdep_fam_mod
