module cedsair_data_mod
use chem_types_mod
use chem_const_mod, only : mw_so2_aer,airmw
implicit none
integer(chem_kind_i4),allocatable,dimension(:,:) :: ktopcedsair,kbotcedsair
integer(chem_kind_i4),parameter :: p_cedsair_bc=1,p_cedsair_oc=2,p_cedsair_so2=3
integer(chem_kind_i4),parameter :: nl_cedsair=25
real(chem_kind_r4),allocatable,dimension(:,:,:,:) :: intcedsair
real(chem_kind_r4),allocatable,dimension(:) :: intaltair
integer(chem_kind_i4),parameter :: num_cedsair=3
character*8 cem_cedsair(num_cedsair)
data cem_cedsair/'BC','OC','SO2'/
contains
  subroutine cedsair_init_data(cedsemi_inname,is,ie,js,je,ks,ke,tile,localpe,imon)
  implicit none
  logical :: first=.true.
  logical exist
  character *(*) cedsemi_inname
  character*8 file*100,cmonth2(2)*2,ctile*1
  real*4 cedsair_in(192,192,nl_cedsair),altair(nl_cedsair)
  integer(chem_kind_i4) :: kktopcedsair(is:ie,js:je,3),kkbotcedsair(is:ie,js:je,3)
  integer is,ie,js,je,ks,ke,i,j,k,n,imon,nlevin
  integer(chem_kind_i4) :: tile,localpe
  real*4 :: mw_bc=12.01
  real*4 :: mw_oc=16.80
  real*4 :: mw_air=28.97
  if(.not.first)return
  write(cmonth2(1),'(i2.2)')imon
  first=.false.
  allocate(ktopcedsair(is:ie,js:je),kbotcedsair(is:ie,js:je))
  allocate(intcedsair(is:ie,js:je,0:nl_cedsair,3))
  allocate(intaltair(0:nl_cedsair))
  kktopcedsair=0
  kkbotcedsair=0
  intaltair(0)=0.0
  intcedsair=0.0
  write(ctile,'(i1.1)')tile
  do n=1,3
      file='tile'//ctile//'/emi_'//trim(cem_cedsair(n))//'-AIR_ugperm2persec.dat'
    inquire(file=trim(cedsemi_inname)//trim(cmonth2(1))//'/'//trim(file),exist=exist)
    if(.not.exist)then
      write(6,*)'file not exist',trim(cedsemi_inname)//trim(cmonth2(1))//'/'//trim(file)
      flush(6)
!      caLL KILLIT('FILE')
    endif
    open(10,file=trim(cedsemi_inname)//trim(cmonth2(1))//'/'//trim(file), &
    convert='big_endian',form='unformatted')
    read(10)nlevin
    read(10)altair
!   make meters for gsdchem code
    altair=altair*1.e3
    if(n==1)then
      intaltair(1)=2.*altair(1)
      do k=2,nl_cedsair
        intaltair(k)=altair(k)+(altair(k)-intaltair(k-1))
      end do
    endif
    read(10)cedsair_in
    if(cem_cedsair(n)=='BC')then
      cedsair_in=cedsair_in ! since ug/m2 not ppv
    elseif(cem_cedsair(n)=='OC')then
      cedsair_in=cedsair_in  ! since ug/m2 not ppv
    elseif(cem_cedsair(n)=='SO2')then
!      cedsair_in=cedsair_in*mw_so2_aer/airmw*.001 ! make umoleso2/m2 /airmw since want ppmv
      cedsair_in=cedsair_in*airmw/mw_so2_aer*.001 ! make umoleso2/m2 *airmw since want ppmv
    endif
    do j=js,je
      do i=is,ie
        if(maxval(cedsair_in(i,j,:))/=0.0)then
        do k=1,nl_cedsair
!          if(i==190.and.j==47.and.tile==1.and.n==3)then
!            write(6,*)'k',k,'cedsair_in',cedsair_in(i,j,k),'kbot',kkbotcedsair(i,j,n)
!          endif
          intcedsair(i,j,k,n)=intcedsair(i,j,k-1,n)+cedsair_in(i,j,k)
          if(cedsair_in(i,j,k)/=0.0)then
            if(kkbotcedsair(i,j,n)==0)then
              kkbotcedsair(i,j,n)=k
            endif
              kktopcedsair(i,j,n)=k
          endif
!          if(i==141.and.j==192.and.tile==1.and.cem_cedsair(n)=='SO2')then
!            write(6,*)'cedsair_in',k,cedsair_in(i,j,k),'intcedsair',intcedsair(i,j,k,n)
!            write(6,*)'kkbot',kkbotcedsair(i,j,n),kktopcedsair(i,j,n)
!          endif
        end do
        endif
      end do
    end do
    close(10)
    
  end do ! n
  do j=js,je
    do i=is,ie
        ktopcedsair(i,j)=maxval(kktopcedsair(i,j,:))
        kbotcedsair(i,j)=minval(kkbotcedsair(i,j,:))
    end do
  enddo
  end subroutine cedsair_init_data
end module cedsair_data_mod
