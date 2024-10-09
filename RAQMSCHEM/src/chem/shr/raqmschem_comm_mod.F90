module raqmschem_comm_mod
  use mpi 
  use chem_types_mod
  use chem_rc_mod
!  use chem_comm_mod, only : chem_comm_get,CHEM_COMM_SUM
  use chem_comm_mod, only : chem_comm_rootpe
  use raqmschem_model_mod, only : chem_model_get => raqmschem_model_get,chem_model_domain_get
  implicit none
  integer :: mpi_comm_chem    = MPI_COMM_NULL
  interface chem_reducetile_pushwithhalo

    module procedure chem_reducetile_pushwithhalo_2dr4
    module procedure chem_reducetile_pushwithhalo_2dr8
    module procedure chem_reducetile_pushwithhalo_3dr8

  end interface chem_reducetile_pushwithhalo
  interface raqmschem_comm_all_bcast
    module procedure chem_comm_all_bcast_0dr4
    module procedure chem_comm_all_bcast_1dr4
    module procedure chem_comm_all_bcast_2dr4
    module procedure chem_comm_all_bcast_3dr4
    module procedure chem_comm_all_bcast_4dr4
    module procedure chem_comm_all_bcast_0di4
    module procedure chem_comm_all_bcast_1di4
    module procedure chem_comm_all_bcast_char
    module procedure chem_comm_all_bcast_char1
    module procedure chem_comm_all_bcast_char2
    module procedure chem_comm_all_bcast_0dr8
    module procedure chem_comm_all_bcast_1dr8
    module procedure chem_comm_all_bcast_2dr8
    module procedure chem_comm_all_bcast_3dr8
    module procedure chem_comm_all_bcast_log
  end interface raqmschem_comm_all_bcast
  interface raqmschem_comm_bcast
    module procedure chem_comm_bcast_r2
    module procedure chem_comm_bcast_d2
  end interface raqmschem_comm_bcast
  interface raqmschem_comm_reduce
    module procedure chem_comm_reduce_r2
    module procedure chem_comm_reduce_d2
  end interface raqmschem_comm_reduce

  private
  public :: chem_reducetile_pushwithhalo
  public :: raqmschem_comm_abort
  public :: raqmschem_comm_all_bcast
  public :: raqmschem_comm_bcaSt
  public :: raqmschem_comm_reduce

contains


  subroutine chem_comm_all_bcast_char(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb
  Character(*)  buffer
  lenb=len(buffer)
  call mpi_bcast(buffer,lenb,MPI_CHARACTER,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char
  subroutine chem_comm_all_bcast_char1(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb,lena
  Character(*)  buffer(:)
  lenb=len(buffer)
  lena=size(buffer)
  call mpi_bcast(buffer,lenb*lena,MPI_CHARACTER,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char1
  subroutine chem_comm_all_bcast_char2(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb,lena
  Character(*)  buffer(:,:)
  lenb=len(buffer)
  lena=size(buffer)
  call mpi_bcast(buffer,lenb*lena,MPI_CHARACTER,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char2
  subroutine chem_comm_all_bcast_log(buffer,rc)
  integer, optional :: rc
  logical :: buffer
  integer localrc
  call mpi_bcast(buffer,1,MPI_LOGICAL,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_log
  subroutine chem_comm_all_bcast_0dr4(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4) :: buffer
  call mpi_bcast(buffer,1,MPI_REAL,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0dr4
  subroutine chem_comm_all_bcast_1dr4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4),dimension(:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1dr4
  subroutine chem_comm_all_bcast_2dr4(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4),dimension(:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_2dr4
  subroutine chem_comm_all_bcast_3dr4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4),dimension(:,:,:) :: buffer
!  write(6,*)'3dr4 ',size(buffer)
!  call flush(6)
!  write(300+iam,*)'all_bcast_3dr4 ',size(buffer)
!  call flush(300+iam)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_world,localrc)
!  write(300+iam,*)'did all_b uffer 3dr4 ',localrc
!  call flush(300+iam)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_3dr4
  subroutine chem_comm_all_bcast_4dr4(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4),dimension(:,:,:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_4dr4
  subroutine chem_comm_all_bcast_0di4(buffer,rc)
  integer, optional :: rc
  integer localrc
  integer :: buffer
  call mpi_bcast(buffer,1,MPI_INTEGER,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0di4
  subroutine chem_comm_all_bcast_1di4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  integer :: buffer(:)
!  write(300+iam,*)'1di4 ',size(buffer)
!  call flush(300+iam)
  call mpi_bcast(buffer,size(buffer),MPI_INTEGER,0,mpi_comm_world,localrc)
!  write(300+iam,*)'1di4 local',localrc
!  call flush(300+iam)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1di4
  subroutine chem_comm_all_bcast_0dr8(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R8) :: buffer
  call mpi_bcast(buffer,1,MPI_DOUBLE,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0dr8
  subroutine chem_comm_all_bcast_1dr8(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R8),dimension(:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1dr8
  subroutine chem_comm_all_bcast_2dr8(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R8),dimension(:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_2dr8
  subroutine chem_comm_all_bcast_3dr8(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R8),dimension(:,:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,0,mpi_comm_world,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_3dr8
  subroutine raqmschem_comm_abort(errcode, msg)
    use raqmschem_model_mod, only : raqmschem_model_get
    integer,          optional, intent(in) :: errcode
    character(len=*), optional, intent(in) :: msg 

    ! -- local variables
    integer :: ierr, localerrcode
    integer :: mpi_comm_chem

    ! -- begin
    localerrcode = CHEM_RC_FAILURE
    call chem_model_get(modelcomm=mpi_comm_chem,rc=ierr)
    if (present(msg)) then
       write(0,'("chem_comm_abort:",a)') trim(msg)
       call flush(0)
    endif
    if (present(errcode)) localerrcode = errcode
    call mpi_abort(mpi_comm_chem, localerrcode, ierr)

  end subroutine raqmschem_comm_abort

    subroutine chem_reducetile_pushwithhalo_2dr4(datain,ids,ide,jds,jde,datawithhalo, &
    ihs,ihe,jhs,jhe,de,rc)
    use chem_types_mod,only :CHEM_KIND_R4
    use chem_comm_mod
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
    call raqmschem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
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
    call raqmschem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
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
    use chem_comm_mod
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
    call raqmschem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
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
    call raqmschem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
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
    use chem_comm_mod
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
      call raqmschem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
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
      call raqmschem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
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


  subroutine chem_comm_bcast_r2(buffer, count, rootpe, comm, rc)
    real(CHEM_KIND_R4), intent(inout) :: buffer(:,:)
    integer, optional,  intent(in)    :: count
    integer, optional,  intent(in)    :: rootpe
    integer, optional,  intent(in)    :: comm
    integer, optional,  intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: localcomm, localcount, root

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS
    if (present(count)) then
      if (count > size(buffer)) return
      localcount = count
    else
      localcount = size(buffer)
    end if
    root = chem_comm_rootpe
    if (present(rootpe)) root = rootpe
    localcomm = mpi_comm_chem
    if (present(comm)) localcomm = comm
!    write(6,*)'call mpi_bcast ',shape(buffer),'localcount',localcount,'root',root
!    call flush(6)
!    call mpi_barrier(localcomm,localrc) ! ajl
    call mpi_bcast(buffer, localcount, MPI_REAL, root, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_bcast_r2
  subroutine chem_comm_bcast_d2(buffer, count, rootpe, comm, rc)
    real(CHEM_KIND_R8), intent(inout) :: buffer(:,:)
    integer, optional,  intent(in)    :: count
    integer, optional,  intent(in)    :: rootpe
    integer, optional,  intent(in)    :: comm
    integer, optional,  intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: localcomm, localcount, root

    ! -- begin
!    write(6,*)'top bcast_d2',shape(buffer)
!    call flush(6)
    if (present(rc)) rc = CHEM_RC_SUCCESS
    if (present(count)) then
      if (count > size(buffer)) return
      localcount = count
    else
      localcount = size(buffer)
    end if
!    write(6,*)'localcount',localcount
!    call flush(6)
    root = chem_comm_rootpe
    if (present(rootpe)) root = rootpe
    localcomm = mpi_comm_chem
    if (present(comm)) localcomm = comm
!    write(6,*)'root for bcast_d2',root
!    call flush(6)
!    write(6,*)'call mpi_bcast ',shape(buffer),'localcount',localcount,'root',root
!    call flush(6)
!    call mpi_barrier(localcomm,localrc) ! ajl
    call mpi_bcast(buffer, localcount, MPI_DOUBLE, root, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'bottom bcast_d2',maxval(buffer),minval(buffer),shape(buffer)
!    call flush(6)
    
  end subroutine chem_comm_bcast_d2


  subroutine chem_comm_reduce_r2(sendbuf, recvbuf, op, rootpe, comm, rc)
    real(CHEM_KIND_R4), intent(in)    :: sendbuf(:,:)
    real(CHEM_KIND_R4), intent(inout) :: recvbuf(:,:)
    integer,            intent(in)    :: op
    integer, optional,  intent(in)    :: rootpe
    integer, optional,  intent(in)    :: comm
    integer, optional,  intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: localcomm, localroot

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localcomm = mpi_comm_chem
    if (present(comm)) localcomm = comm

    localroot = chem_comm_rootpe
    if (present(rootpe)) localroot = rootpe

    call mpi_reduce(sendbuf, recvbuf, size(sendbuf), MPI_REAL, op, localroot, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_comm_reduce_r2
  subroutine chem_comm_reduce_d2(sendbuf, recvbuf, op, rootpe, comm, rc)
    real(CHEM_KIND_R8), intent(in)    :: sendbuf(:,:)
    real(CHEM_KIND_R8), intent(inout) :: recvbuf(:,:)
    integer,            intent(in)    :: op
    integer, optional,  intent(in)    :: rootpe
    integer, optional,  intent(in)    :: comm
    integer, optional,  intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: localcomm, localroot

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localcomm = mpi_comm_chem
    if (present(comm)) localcomm = comm

    localroot = chem_comm_rootpe
    if (present(rootpe)) localroot = rootpe

    call mpi_reduce(sendbuf, recvbuf, size(sendbuf), MPI_DOUBLE, op, localroot, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
   end subroutine chem_comm_reduce_d2
end module raqmschem_comm_mod
