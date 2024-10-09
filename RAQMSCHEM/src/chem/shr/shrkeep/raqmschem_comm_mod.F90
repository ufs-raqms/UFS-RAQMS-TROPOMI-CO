module raqmschem_comm_mod
  use mpi 
  use chem_types_mod
  use chem_rc_mod
!  use chem_comm_mod, only : chem_comm_get,CHEM_COMM_SUM
  use chem_comm_mod
  use chem_model_mod, only : chem_model_get => raqmschem_model_get,chem_model_domain_get
  implicit none
  interface chem_reducetile_pushwithhalo

    module procedure chem_reducetile_pushwithhalo_2dr4
    module procedure chem_reducetile_pushwithhalo_2dr8
    module procedure chem_reducetile_pushwithhalo_3dr8

  end interface chem_reducetile_pushwithhalo

  private
  public :: chem_reducetile_pushwithhalo

contains

    subroutine chem_reducetile_pushwithhalo_2dr4(datain,ids,ide,jds,jde,datawithhalo, &
    ihs,ihe,jhs,jhe,de,rc)
    use chem_types_mod,only :CHEM_KIND_R4
    real(CHEM_KIND_R4),dimension(ids:ide,jds:jde),intent(in) :: datain
    real(CHEM_KIND_R4),dimension(ihs:ihe,jhs:jhe),intent(out) :: datawithhalo
    ! -- local variables
    integer,optional :: rc,de
    integer :: localrc,tileCOMM
    integer :: localpe, tile 
    integer ids,ide,jds,jde,ihs,ihe,jhs,jhe,its,ite,jts,jte,mype,i
!    type(chem_data_type),   pointer :: data => null()
!    type(chem_config_type), pointer :: config => null()
    real(CHEM_KIND_R4),allocatable,dimension(:,:) :: recvbuf,buf2d
    logical localIOflag
!    write(6,*)'top reducetile'
!    call flush(6)
    call chem_comm_get(localpe=mype)
     

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS


    call chem_comm_get(localpe=localpe, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_get(de=de, tile=tile, tileCOMM=tileCOMM, localIOflag=localIOflag,rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!   need to get domain limits and define limits with haloe
    call chem_model_domain_get(de=de,  &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
!    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!    ihs=max(ids-1,its)
! ihe=min(ide+1,ite)
!    jhs=max(jds-1,jts)
!    jhe=min(jde+1,jte)
!   allocate sendbufr
    ! -- check size consistency
    if (chem_rc_test((size(datain) /= (ide-ids+1)*(jde-jds+1)), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'allocate buf2d',its,ite,jts,jte

!    call flush(6)

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buf2d = 0._CHEM_KIND_R4

    buf2d(ids:ide, jds:jde) = datain
!    if(tile.eq.1)then
!      write(6,*)mype,'buf2d tile',ids,ide,jds,jde,'maxval',maxval(datain),minval(datain)
!      call flush(6)
!    endif
!   allocate reduce bufr
    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4
!   do reduction
!    write(70+mype,*)'call chem_comm_reduce',shape(buf2d),shape(recvbuf)
!    call flush(70+mype)
    call chem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!    if(tile.eq.1)then
!      write(6,*)'chem_comm_reduce tile 1 ',mype,maxval(recvbuf),minval(recvbuf),'localIOflag',localIOflag
!      call flush(6)
!    endif
!    write(6,*)'bcast ajl',shape(buf2d),lbound(buf2d),ubound(buf2d)
!    call flush(6)
    !write(70+mype,*)'bcast ajl',shape(buf2d),lbound(buf2d),ubound(buf2d)
!    call flush(70+mype)
!   now need to bcast and then put out data with tile inside halo included
!    if(mype==0)then
!      do i=its,ite
!        write(6,*)'ajl recvbuf push',i,recvbuf(i,48:49)
!      end do
!    endif
    call chem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
!    if(tile.eq.1)then
!      write(6,*)mype,'ihs',ihs,ihe,jhs,jhe
!      write(6,*)mype,'recvbuf',maxval(recvbuf),minval(recvbuf)
!    !endif
    datawithhalo=recvbuf(ihs:ihe,jhs:jhe)
    deallocate (buf2d,recvbuf)
!    write(6,*)'bottom reducetile',tile,maxval(datawithhalo),minval(datawithhalo)
    !call flush(6)
    return
    end subroutine chem_reducetile_pushwithhalo_2dr4

    subroutine chem_reducetile_pushwithhalo_2dr8(datain,ids,ide,jds,jde,datawithhalo, &
    ihs,ihe,jhs,jhe,de,rc)
    use chem_types_mod,only :CHEM_KIND_R4
    real(CHEM_KIND_R8),dimension(ids:ide,jds:jde),intent(in) :: datain
    real(CHEM_KIND_R8),dimension(ihs:ihe,jhs:jhe),intent(out) :: datawithhalo
    ! -- local variables
    integer,optional :: rc,de
    integer :: localrc,tileCOMM
    integer :: localpe, tile 
    integer ids,ide,jds,jde,ihs,ihe,jhs,jhe,its,ite,jts,jte,mype,i
!    type(chem_data_type),   pointer :: data => null()
!    type(chem_config_type), pointer :: config => null()
    real(CHEM_KIND_R8),allocatable,dimension(:,:) :: recvbuf,buf2d
    logical localIOflag
!    write(6,*)'top reducetile'
!    call flush(6)
    call chem_comm_get(localpe=mype)
     

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS


    call chem_comm_get(localpe=localpe, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_get(de=de, tile=tile, tileCOMM=tileCOMM, localIOflag=localIOflag,rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!   need to get domain limits and define limits with haloe
    call chem_model_domain_get(de=de,  &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
!    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!    ihs=max(ids-1,its)
! ihe=min(ide+1,ite)
!    jhs=max(jds-1,jts)
!    jhe=min(jde+1,jte)
!   allocate sendbufr
    ! -- check size consistency
    if (chem_rc_test((size(datain) /= (ide-ids+1)*(jde-jds+1)), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'allocate buf2d',its,ite,jts,jte

!    call flush(6)

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buf2d = 0._CHEM_KIND_R4

    buf2d(ids:ide, jds:jde) = datain
!    if(tile.eq.1)then
!      write(6,*)mype,'buf2d tile',ids,ide,jds,jde,'maxval',maxval(datain),minval(datain)
!      call flush(6)
!    endif
!   allocate reduce bufr
    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4
!   do reduction
!    write(70+mype,*)'call chem_comm_reduce',shape(buf2d),shape(recvbuf)
!    call flush(70+mype)
    call chem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!    if(tile.eq.1)then
!      write(6,*)'chem_comm_reduce tile 1 ',mype,maxval(recvbuf),minval(recvbuf),'localIOflag',localIOflag
!      call flush(6)
!    endif
!    write(6,*)'bcast ajl',shape(buf2d),lbound(buf2d),ubound(buf2d)
!    call flush(6)
    !write(70+mype,*)'bcast ajl',shape(buf2d),lbound(buf2d),ubound(buf2d)
!    call flush(70+mype)
!   now need to bcast and then put out data with tile inside halo included
!    if(mype==0)then
!      do i=its,ite
!        write(6,*)'ajl recvbuf push',i,recvbuf(i,48:49)
!      end do
!    endif
    call chem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
!    if(tile.eq.1)then
!      write(6,*)mype,'ihs',ihs,ihe,jhs,jhe
!      write(6,*)mype,'recvbuf',maxval(recvbuf),minval(recvbuf)
!    !endif
    datawithhalo=recvbuf(ihs:ihe,jhs:jhe)
    deallocate (buf2d,recvbuf)
!    write(6,*)'bottom reducetile',tile,maxval(datawithhalo),minval(datawithhalo)
    !call flush(6)
    return
    end subroutine chem_reducetile_pushwithhalo_2dr8
    subroutine chem_reducetile_pushwithhalo_3dr8(datain,ids,ide,jds,jde,nl,datawithhalo, &
    ihs,ihe,jhs,jhe,de,rc)
    use chem_types_mod,only :CHEM_KIND_R4
    real(CHEM_KIND_R8),dimension(ids:ide,jds:jde,nl),intent(in) :: datain
    real(CHEM_KIND_R8),dimension(ihs:ihe,jhs:jhe,nl),intent(out) :: datawithhalo
    ! -- local variables
    integer,optional :: rc,de
    integer :: localrc,tileCOMM
    integer :: localpe, tile 
    integer ids,ide,jds,jde,ihs,ihe,jhs,jhe,its,ite,jts,jte,mype,nl,k,i
!    type(chem_data_type),   pointer :: data => null()
!    type(chem_config_type), pointer :: config => null()
    real(CHEM_KIND_R8),allocatable,dimension(:,:) :: recvbuf,buf2d
    logical localIOflag
!    write(6,*)'top reducetile 3d'
!    call flush(6)
    call chem_comm_get(localpe=mype)
     

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS


    call chem_comm_get(localpe=localpe, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_get(de=de, tile=tile, tileCOMM=tileCOMM, localIOflag=localIOflag,rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!   need to get domain limits and define limits with haloe
    call chem_model_domain_get(de=de,  &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
!    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!    ihs=max(ids-1,its)
! ihe=min(ide+1,ite)
!    jhs=max(jds-1,jts)
!    jhe=min(jde+1,jte)
!   allocate sendbufr
    ! -- check size consistency
    if (chem_rc_test((size(datain) /= (ide-ids+1)*(jde-jds+1)*nl), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'allocate buf2d',its,ite,jts,jte
!    call flush(6)

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    do k=1,nl
      buf2d = 0._CHEM_KIND_R4

      buf2d(ids:ide, jds:jde) = datain(:,:,k)
#if 0
      if(tile.eq.1)then
        write(6,*)mype,'buf2d tile',ids,ide,jds,jde,'maxval',maxval(datain),minval(datain)
        call flush(6)
      endif
#endif
!     allocate reduce bufr

      recvbuf = 0._CHEM_KIND_R4
!     do reduction
      call chem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
#if 0
      if(tile.eq.1)then
        if(mype.eq.0)then
          if(k.eq.63)then
            do i=its,ite
              write(6,*)'ajl reduce 3d i',i,recvbuf(i,48)
            end do
          endif
        endif
        write(6,*)'chem_comm_reduce tile 1 ',mype,maxval(recvbuf),minval(recvbuf),'localIOflag',localIOflag
        call flush(6)
      endif
#endif
!      write(6,*)'bcast ajl',shape(buf2d),lbound(buf2d),ubound(buf2d)
!      call flush(6)
!      write(70+mype,*)'bcast ajl',shape(buf2d),lbound(buf2d),ubound(buf2d)
!      call flush(70+mype)
!     now need to bcast and then put out data with tile inside halo included
      call chem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
#if 0
      if(tile.eq.1)then
        if(mype<=1)then
          if(k.eq.63)then
            do i=its,ite
              write(6,*)'ajl reduce 3d afger bcast i',i,recvbuf(i,48)
            end do
          endif
        endif
        write(6,*)mype,'ihs',ihs,ihe,jhs,jhe
        write(6,*)mype,'recvbuf',maxval(recvbuf),minval(recvbuf)
      endif
#endif
      datawithhalo(:,:,k)=recvbuf(ihs:ihe,jhs:jhe)
    end do
    deallocate (buf2d,recvbuf)
!    write(6,*)'bottom reducetile',tile,maxval(datawithhalo),minval(datawithhalo)
!    call flush(6)
    return
    end subroutine chem_reducetile_pushwithhalo_3dr8
end module raqmschem_comm_mod
