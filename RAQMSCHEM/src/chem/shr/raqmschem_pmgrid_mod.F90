module raqmschem_pmgrid_mod
public
integer nc,plat,plev,begj,endj,nmonHTAP,iam,nr,nlev,plon,beglat,endlat,jbeg,plond
integer ibeg,iend,tile
logical masterproc
character *240 bbedirf,bbedira,prefi
character *19 cdatestr
character *10 cdateat,cdatestart
integer nhtuw
integer iprn,jprn,kprn,iamprn,iprnin,jprnin,kprnin,pecount,tileprn
logical :: arcpac=.false.
logical :: calnex=.false.
logical :: twolayer_emissions=.false.
integer iatchem,jatchem,katchem,nstepat,iamprint
logical :: skipphotoj=.false.
! move these to pmgrid out of omp loop
logical dopr94,doareafix,usescale_factor
real*8 fracday,juldayat
real*4 forecast_hr,timestep
integer numgsicem
character *8,dimension(:),allocatable :: gsicem,gsivarr
integer,dimension(:),allocatable,target :: loccem
logicaL raqms_localIOflag
logical :: aerosol_ugpkg=.true.
logical :: aerosol_in_old=.false.
logical :: ompsgsi=.false.
  type gsidatatype
    integer ngsivar,nsecgsi
    real*4,allocatable,dimension(:,:,:,:) :: gsiinc
    integer,allocatable,dimension(:) :: loccem
    character *8,allocatable,dimension(:) :: gsicem
    character *8 :: gsivarr
  end type
  type(gsidatatype),allocatable,target :: gsiarray(:)
  character *8 cemaod(15)
  data cemaod/'sulf','bc1','bc2','oc1','oc2','dust1','dust2','dust3','dust4','dust5', &
     'seas1','seas2','seas3','seas4','seas5'/

contains
  subroutine initraqmschem_pmgrid(ncin,platin,plevin,begjin,endjin,ibegin,iendin,nmonHTAPin,iamin)
!  use mpi
  implicit none
  integer ncin,platin,plevin,begjin,endjin,nmonHTAPin,iamin,ierr
  integer ibegin,iendin
  character *20 ctest
  juldayat=0
  fracday=0.0
  ctest=' '
  call getenv('AEROSOL_UNITS',ctest)
  if(ctest=='raqms')then
    aerosol_ugpkg=.false.
  endif
  if(tile.eq.1.and.iamin.eq.0)then
    write(6,*)'aerosol_ugpkg',aerosol_ugpkg
  endif
 
  ctest=' '
  call getenv('AEROSOL_IN_OLD',ctest)
  if(ctest=='YES')then
    aerosol_in_old=.true.
  endif
  ctest=' '
  call getenv('USEPR94',ctest)
  if(ctest.eq.'YES'.or.ctest.eq.'yes')then
    dopr94=.true.
    if(tile.eq.1.and.iamin.eq.0)then
      write(6,*)'do pr94 for lightning'
      call flush(6)
    endif
  else
    dopr94=.false.
    if(tile.eq.1.and.iamin.eq.0)then
      write(6,*)'dont do pr94 for lightning'
      call flush(6)
    endif
  endif
  ctest=' ' 
  call getenv("AREAFIX",ctest)
  if(ctest.eq.'YES'.or.ctest.eq.'yes')then
    doareafix=.true.
    if(tile.eq.1.and.iamin.eq.0)then
      write(6,*)'doa reafix for lightning'
    endif
  else
    doareafix=.false.
    if(tile.eq.1.and.iamin.eq.0)then
      write(6,*)'dont do areafix for lightning'
    endif
  endif
  ctest=' '
  call getenv('USESCALEFACTOR',ctest)
  if(ctest.eq.'NO')then
    usescale_factor=.false.
    if(tile.eq.1.and.iamin.eq.0)then
      write(6,*)'dont use new fire scale factor'
    endif
  else
    usescale_factor=.true.
    if(tile.eq.1.and.iamin.eq.0)then
      write(6,*)'do use new fire scale factor'
    endif
  endif
  ctest=' '
  iatchem=0
  jatchem=0
  katchem=0
  call getenv('SKIPPHOTO',ctest)
  if(ctest.ne.' ' )then
    if(ctest.eq.'YES')then
      write(6,*)'tile',tile,'iamin',iamin
      call flush(6)
      if(tile.eq.1.and.iamin.eq.0)then
         write(6,*)'SKIP PHOTOJ'
         call flush(6)
      else
         write(6,*)'problem SKIP'
         call flush(6)
      endif
      skipphotoj=.true.
    endif
  else
    if(tile.eq.1.and.iamin.eq.0)then
      write(6,*)'NO SKIP PHOTOJ'
      call flush(6)
    endif
  endif
  if(skipphotoj)then
      write(6,*)'skipphotoj ',iamin
      call flush(6)
  endif
  ctest=' '
  call getenv('TWOLAYER_EMISSIONS',ctest)
  if(ctest.ne.' ')then
    if(ctest.eq.'YES'.or.ctest.eq.'yes')then
      twolayer_emissions=.true.
      if(tile.eq.1.and.iamin.eq.0)then
!      if(masterproc)then
        write(6,*)'use twolayer_emsissions'
      endif
    endif
  else
   if(tile.eq.1.and.iamin.eq.0)then
      write(6,*)'dont use twolayer emissions'
   endif
  endif
  call getenv('NHTUW',ctest)
  if(ctest.ne.' ')then
    read(ctest,*)nhtuw
  else
    nhtuw=6
  endif
  nc=ncin
  plon=nc
  plond=nc
  nr=platin
!  write(6,*)'set nr and plat ',platin,'nc',nc,'begj',begjin,endjin
!  call flush(6)
  nlev=plevin
  plat=platin
  plev=plevin
  begj=begjin
  endj=endjin
  jbeg=begj
  beglat=begj
  endlat=endj
  ibeg=ibegin
  iend=iendin
  nmonHTAP=nmonHTAPin
  iam=iamin
  masterproc=.false.
  if(iam.eq.0)then
    write(6,*)'initraqmchem_pmgrid',nc,plat,plev,begj,endj,nmonhtap,'ibeg',ibeg,iend
    masterproc=.true.
  endif
  return
  end subroutine initraqmschem_pmgrid
!  subroutine gsivar(numvar,cemgsi,loccem)
  subroutine gsivar
  use field_manager_mod, only : find_field_index
  use raqmschem_species_mod,only : p_ox,idaodgsi
  implicit none
!  character *(*) cemgsi(*)
!  integer numvar,is(30),ie(30),igot,num,i,loccem(*)
  integer numvar,is(30),ie(30),igot,num,i,n,nsecgsi(30),is2(30),ie2(30),igot2,num2
  character *256 cinput,cinput2
  logical :: first=.true.
  save first
  cinput=' '
  call getenv('GSIVAR',cinput)
  cinput2=' '
  call getenv('GSISEC',cinput2)
  igot2=0
  if(cinput2/=' ')then
    num2=30
    call strip(cinput2,is2,ie2,num2,igot2)
    do i=1,igot2
      read(cinput2(is2(i):ie2(i)),*)nsecgsi(i)
      if(first)then
      if(iam==0)then
      write(6,*)'nsecgsi',i,nsecgsi(i)
      endif
      first=.false.
      endif
    end do
  endif

  if(cinput/=' ')then
    num=30
    call strip(cinput,is,ie,num,igot)
    numgsicem=igot
!    if(iam==0)then
!      write(6,*)'numgsicem',numgsicem,'allocated',allocated(loccem)
!      call flush(6)
!    endif
    if(.not.allocated(gsiarray))then
      allocate(gsiarray(numgsicem))
    endif
    if(.not.allocated(loccem))then
      allocate (loccem(numgsicem),gsicem(numgsicem),gsivarr(numgsicem))
      if(iam==0)then
        write(6,*)'GSIVAR',cinput
      endif
    else
       return
    endif
    do i=1,numgsicem
      gsicem(i)=cinput(is(i):ie(i))
      gsivarr(i)=cinput(is(i):ie(i))
      if(iam==0)then
        write(6,*)'gsicem',i,gsicem(i),'p_ox',p_ox
        call flush(6)
      endif
      gsiarray(i)%gsivarr=gsivarr(i)
      if(gsicem(i).eq.'o3vmr')then
        allocate (gsiarray(i)%loccem(1),gsiarray(i)%gsicem(1))
        loccem(i)=p_ox
        gsiarray(i)%ngsivar=1
        gsiarray(i)%loccem(1)=loccem(i)
        gsiarray(i)%gsicem(1)=gsicem(i)
        if(igot2>=i)then
          gsiarray(i)%nsecgsi=nsecgsi(i)
        elseif(igot2==1)then
          gsiarray(i)%nsecgsi=nsecgsi(1)
        else
          gsiarray(i)%nsecgsi=-1
        endif
      elseif(gsicem(i).eq.'AOD')then
        gsiarray(i)%ngsivar=15
        allocate (gsiarray(i)%loccem(15),gsiarray(i)%gsicem(15))
        do n=1,15
          gsiarray(i)%loccem(n)=idaodgsi(n)
          gsiarray(i)%gsicem(n)=cemaod(n)
!          if(iam==0)then
!            write(6,*)'gsicem',i,gsicem(i)
!            write(6,*)'cem',n,cemaod(n),'loccem',idaodgsi(n)
!          endif
        end do
        if(igot2>=i)then
          gsiarray(i)%nsecgsi=nsecgsi(i)
        elseif(igot2==1)then
          gsiarray(i)%nsecgsi=nsecgsi(1)
        else
          gsiarray(i)%nsecgsi=-1
        endif
      else
        allocate (gsiarray(i)%loccem(1),gsiarray(i)%gsicem(1))
        loccem(i)=find_field_index(1,gsicem(i))
        gsiarray(i)%loccem(1)=loccem(i)
        gsiarray(i)%ngsivar=1
        gsiarray(i)%gsicem(1)=gsicem(i)
        if(igot2>=i)then
          gsiarray(i)%nsecgsi=nsecgsi(i)
        elseif(igot2==1)then
          gsiarray(i)%nsecgsi=nsecgsi(1)
        else
!          if(igot2==0)then
!           write(6,*)'igot2 zero ',i,gsicem(i)
!           call flush(6)
!          endif
          gsiarray(i)%nsecgsi=-1
        endif
      endif
!      write(6,*)'loccem',i,loccem(i),gsicem(i)
!      call flush(6)
!      if(iam==0)then
!        write(6,*)'i',i,'gsicem',gsicem(i),'ngsivar',gsiarray(i)%ngsivar
!        do n=1,gsiarray(i)%ngsivar
!          write(6,*)'loccem',gsiarray(i)%loccem(n),'gsichem',gsiarray(i)%gsicem(n)
!        end do
!      endif
    end do
  else
    numgsicem=0
  endif
  return
  end subroutine gsivar

      subroutine strip(input,is,ie,num,igot)
      character *(*) input
      integer is(num),ie(num)
      logical quote,skipblank
      quote=.false.
      skipblank=.true.
      ilen=len(input)
      id=1
      igot=0
      ifound=0
      do i=1,ilen
          if(input(i:i).eq."'")then
            if(quote)then
!             quote was on so this is end of character string
              ie(id)=i-1
              if(id.eq.num)go to 999
              id=id+1
              skipblank=.true.
              quote=.false.
            else
!             new quote so start string
              is(id)=i+1
              quote=.true.
              igot=igot+1
              ifound=1
            endif
          else
            if(.not.quote)then
!           if quote not on then continue parse else keep moving till
!           hit another quote
              if(input(i:i).eq.' '.or.input(i:i).eq.',')then
                if(.not.skipblank)then
!                 hit a blank after a number so set end
!                 print *,'hit blank previous was ',input(i-1:i-1),i
                  ie(id)=i-1
                  igot=igot+1
                  if(id.eq.num)go to 999
                  id=id+1
                  skipblank=.true.
                endif
              else
!               must be a number so continue
                if(skipblank)then
                  is(id)=i
!              print *,'hit a number input=',input(i:i)
                  skipblank=.false.
                endif
!               print *,'hit a number set ie input=',input(i:i)
                ie(id)=i
              endif
            endif
          endif
        end do
999     continue
        return
        endsubroutine strip
      subroutine cloweruw(text,cout)
      character *(*) text,cout*(*)
      lenc=len(text)
      lenc=min(lenc,len(cout))
      icapa=ichar('A')
      ismalla=ichar('a')
      cout=' '
      do i=1,lenc
        if(text(i:i).ge.'A'.and.text(i:i).le.'Z')then
          ic=ichar(text(i:i))
          cout(i:i)=char(ic-icapa+ismalla)
        else
          cout(i:i)=text(i:i)
        endif
      end do
      return
      end subroutine cloweruw
! SGI version of routine
!
!-----------------------------------------------------------------------
subroutine rsearch(im,km1,ixz1,kxz1,z1,km2,ixz2,kxz2,z2,ixl2,kxl2,l2)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    rsearch     search for a surrounding real interval
!   prgmmr: iredell          org: np23                date: 1998-05-01
!
! abstract: This subprogram searches monotonic sequences of real numbers
!   for intervals that surround another given search set of real numbers.
!   The sequences may be monotonic in either direction; the real numbers
!   may be single or double precision; the input sequences and sets
!   and the output locations may be arbitrarily dimensioned.
!
! program history log:
!   1999-01-05  mark iredell
!   2008-04-16  safford - rm unused vars
!
!   input argument list:
!     im           integer number of sequences to search
!     km1          integer number of points in each sequence
!     ixz1         integer sequence skip number for z1
!     kxz1         integer point skip number for z1
!     z1           real (1+(im-1)*ixz1+(km1-1)*kxz1)
!                  sequence values to search
!                  (z1 must be monotonic in either direction)
!     km2          integer number of points to search for
!                  in each respective sequence
!     ixz2         integer sequence skip number for z2
!     kxz2         integer point skip number for z2
!     z2           real (1+(im-1)*ixz2+(km2-1)*kxz2)
!                  set of values to search for
!                  (z2 need not be monotonic)
!     ixl2         integer sequence skip number for l2
!     kxl2         integer point skip number for l2
!     
!   output argument list:
!     l2           real (1+(im-1)*ixl2+(km2-1)*kxl2)
!                  interval locations of the set of values
!                  (z2 will be between z1(l2) and z1(l2+1).)
!
! remarks:
!   If the array z1 is dimensioned (im,km1), then the skip numbers are
!   ixz1=1 and kxz1=im; if it is dimensioned (km1,im), then the skip
!   numbers are ixz1=km1 and kxz1=1; if it is dimensioned (im,jm,km1),
!   then the skip numbers are ixz1=1 and kxz1=im*jm; etcetera.
!   similar examples apply to the skip numbers for z2 and l2.
!
!   Returned values of 0 or km1 indicate that the given search value
!   is outside the range of the sequence.
!
!   If a search value is identical to one of the sequence values
!   then the location returned points to the identical value.
!   If the sequence is not strictly monotonic and a search value is
!   identical to more than one of the sequence values, then the
!   location returned may point to any of the identical values.
!
!   To be exact, for each i from 1 to im and for each k from 1 to km2,
!   z=z2(1+(i-1)*ixz2+(k-1)*kxz2) is the search value and
!   l=l2(1+(i-1)*ixl2+(k-1)*kxl2) is the location returned.
!   If l=0, then z is less than the start point z1(1+(i-1)*ixz1)
!   for ascending sequences (or greater than for descending sequences).
!   If l=km1, then z is greater than or equal to the end point
!   z1(1+(i-1)*ixz1+(km1-1)*kxz1) for ascending sequences
!   (or less than or equal to for descending sequences).
!   otherwise z is between the values z1(1+(i-1)*ixz1+(l-1)*kxz1) and
!   z1(1+(i-1)*ixz1+(l-0)*kxz1) and may equal the former.
!
! attributes:
!   language:  f90
!   machine:   RS/6000 SP
!
!$$$ end documentation block

  implicit none
 
  integer,intent(in   ) :: im,km1,ixz1,kxz1,km2,ixz2,kxz2,ixl2,kxl2
  real*4   ,intent(in   ) :: z1(1+(im-1)*ixz1+(km1-1)*kxz1)
  real*4   ,intent(in   ) :: z2(1+(im-1)*ixz2+(km2-1)*kxz2)
  integer,intent(  out) :: l2(1+(im-1)*ixl2+(km2-1)*kxl2)

  integer i,k2,l
  real*4 z
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Find the surrounding input interval for each output point.
  do i=1,im
     if(z1(1+(i-1)*ixz1)<=z1(1+(i-1)*ixz1+(km1-1)*kxz1)) then

!       Input coordinate is monotonically ascending.
        do k2=1,km2
           z=z2(1+(i-1)*ixz2+(k2-1)*kxz2)
           l=0
           do
              
              if(z<z1(1+(i-1)*ixz1+l*kxz1)) exit
              l=l+1
              if(l==km1) exit
           enddo
           l2(1+(i-1)*ixl2+(k2-1)*kxl2)=l
        enddo

     else

!        write(200+iam,*)'descend'
!    Input coordinate is monotonically descending.
        do k2=1,km2
           z=z2(1+(i-1)*ixz2+(k2-1)*kxz2)
           l=0
           do
!              write(200+iam,*)'z',z,'k2',k2,'z1','l=',l,z1(1+(i-1)*ixz1+l*kxz1)
              if(z>z1(1+(i-1)*ixz1+l*kxz1)) exit
              l=l+1
              if(l==km1) exit
           enddo
           l2(1+(i-1)*ixl2+(k2-1)*kxl2)=l
        enddo
     endif

  enddo
end subroutine rsearch
end module raqmschem_pmgrid_mod
