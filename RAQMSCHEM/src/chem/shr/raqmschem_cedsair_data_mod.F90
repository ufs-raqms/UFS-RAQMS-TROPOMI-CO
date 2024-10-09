module raqmschem_cedsair_data_mod
use raqmschem_cedsair_mod
contains
  subroutine cedsair_init_data(m,ids,ide,jds,jde,tile,pathin,cmonth)
  use raqmschem_config_mod, only : cedsairemis_opt
  use raqmschem_cedsair_mod, only : lcedsair,cedsairstep,numcedsair
  use raqmschem_cedsair_mod, only : p_ceds_co,p_ceds_nox
  use raqmschem_cedsair_mod, only : p_ceds_n2o
  use raqmschem_pmgrid_mod, only : iam,strip
  implicit none 
  character *256 file,ctile*1
  character *(*) pathin,cmonth
  integer m,nlevin,k,i,j,ids,ide,jds,jde,igot
  integer is(10),ie(10),tile,n,ii
  real*4,allocatable :: air3d(:,:,:)
  integer,parameter :: nlnairceds=25
  real*4 :: altair(nlnairceds)

  if(cedsairemis_opt==1)then
    numcedsair=3
    lcedsair=.true.
    ccedsair=(/'CO','NOx','N2O'/)
    do k=1,numcedsair
      select caSe (ccedsair(k))
      case ('CO')
        p_ceds_co=k
      case ('NOx')
        p_ceds_nox=k
      case ('N2O')
        p_ceds_n2o=k
      case default
        write(6,*)'bad ceds cehm',ccedsair(k)
        call killit('stop cedsin')
      end select
    end do
  endif
  if(iam==0)then
    write(6,*)'numcedsair',numcedsair
    write(6,*)'P_ceds_co',p_ceds_co,'nox',p_ceds_nox,'n2o',p_ceds_n2o
    do k=1,numcedsair
      write(6,*)'cedsair',k,ccedsair(k)
    end do
  endif
  write(ctile,'(i1.1)')tile
  if(.not.allocated(intcedsair).and.numcedsair/=0)then
    allocate(intcedsair(ide-ids+1,jds:jde,0:nlnairceds,numcedsair))
    intcedsair=0.0
    allocate(ktopcedsair(ide-ids+1,jds:jde,numcedsair))
    allocate(kbotcedsair(ide-ids+1,jds:jde,numcedsair))
    ktopcedsair=0
    kbotcedsair=0
    if(.not.allocated(air3d))then
      allocate(air3d(192,192,nlnairceds))
    endif
    if(.not.allocated(intaltair))then
      allocate (intaltair(0:nlnairceds))
    endif
    do n=1,numcedsair
      file='tile'//ctile//'/emi_'//trim(ccedsair(n))//'-AIR_molpercm2persec.dat'
      open(10,file=trim(pathin)//trim(cmonth)//'/'//trim(file), &
      convert='big_endian',form='unformatted')
      read(10)nlevin
      if(iam==75.and.n==1)then
        write(6,*)'nlevin',nlevin
        call flush(6)
      endif
      read(10)altair
      if(iam==75.and.n==1)then
        write(6,*)'altair',altair
        call flush(6)
      endif
      if(n==1)then
        intaltair(0)=0
        intaltair(1)=2.*altair(1)
        do k=2,nlnairceds
          intaltair(k)=altair(k)+(altair(k)-intaltair(k-1))
        end do
        if(iam==75)then
          do k=1,nlnair
            write(6,*)k,altair(k),intaltair(k)
          end do
        endif
      endif
      read(10)air3d
      do j=jds,jde
        do i=ids,ide
          ii=i-ids+1
          do k=1,nlnairceds
            intcedsair(i-ids+1,j,k,n)=intcedsair(i-ids+1,j,k-1,n)+air3d(i,j,k)
            if(kbotcedsair(ii,j,n)==0.and.air3d(i,j,k)/=0.)then
              kbotcedsair(ii,j,n)=k
            endif
            if(air3d(i,j,k)/=0.0)then
              ktopcedsair(ii,j,n)=k
            endif
          end do
        end do ! i
      end do ! j
      close(10)
    end do ! n
  endif
  end subroutine cedsair_init_data
end module raqmschem_cedsair_data_mod
