module chem_comm_mod

  use mpi
  use chem_types_mod
  use chem_rc_mod

  implicit none

  integer :: mpi_comm_chem    = MPI_COMM_NULL
  integer :: mpi_group_chem   = MPI_GROUP_NULL
  integer :: chem_comm_rootpe = 0

  integer, parameter :: CHEM_COMM_MAX     = MPI_MAX     ! return the maximum
  integer, parameter :: CHEM_COMM_MIN     = MPI_MIN     ! return the minumum
  integer, parameter :: CHEM_COMM_SUM     = MPI_SUM     ! return the sum
  integer, parameter :: CHEM_COMM_PROD    = MPI_PROD    ! return the product
  integer, parameter :: CHEM_COMM_LAND    = MPI_LAND    ! return the logical and
  integer, parameter :: CHEM_COMM_BAND    = MPI_BAND    ! return the bitwise and
  integer, parameter :: CHEM_COMM_LOR     = MPI_LOR     ! return the logical or
  integer, parameter :: CHEM_COMM_BOR     = MPI_BOR     ! return the bitwise of
  integer, parameter :: CHEM_COMM_LXOR    = MPI_LXOR    ! return the logical exclusive or
  integer, parameter :: CHEM_COMM_BXOR    = MPI_BXOR    ! return the bitwise exclusive or
  integer, parameter :: CHEM_COMM_MINLOC  = MPI_MINLOC  ! return the minimum and the location (actually, the value of
                                                        ! the second element of the structure where the minimum of
                                                        ! the first is found)
  integer, parameter :: CHEM_COMM_MAXLOC  = MPI_MAXLOC  ! return the maximum and the location
  integer, parameter :: CHEM_COMM_REPLACE = MPI_REPLACE ! replace b with a
  integer, parameter :: CHEM_COMM_NO_OP   = MPI_NO_OP   ! perform no operation


  interface chem_comm_bcast
    module procedure chem_comm_bcast_i0
    module procedure chem_comm_bcast_i1
    module procedure chem_comm_bcast_r0
    module procedure chem_comm_bcast_r1
    module procedure chem_comm_bcast_r2
    module procedure chem_comm_bcast_d2
    module procedure chem_comm_bcast_s0
    module procedure chem_comm_bcast_s1
  end interface chem_comm_bcast

  interface chem_comm_allgather
    module procedure chem_comm_allgather_i1
  end interface chem_comm_allgather

  interface chem_comm_reduce
    module procedure chem_comm_reduce_r1
    module procedure chem_comm_reduce_r2
    module procedure chem_comm_reduce_d2
    module procedure chem_comm_reduce_r3
  end interface chem_comm_reduce

  interface chem_comm_create
    module procedure chem_comm_create_group
    module procedure chem_comm_create_comm
  end interface chem_comm_create
#if 0
  interface chem_comm_reducetile_pushwithhalo
    module procedure chem_reducetile_pushwithhalo_2dr4
    module procedure chem_reducetile_pushwithhalo_2dr8
    module procedure chem_reducetile_pushwithhalo_3dr8
  end interface chem_comm_reducetile_pushwithhalo
#endif

  private

  public :: chem_comm_rootpe
  public :: CHEM_RC_SUCCESS
  public :: CHEM_RC_FAILURE
  public :: CHEM_COMM_MAX
  public :: CHEM_COMM_MIN
  public :: CHEM_COMM_SUM
  public :: CHEM_COMM_PROD
  public :: CHEM_COMM_LAND
  public :: CHEM_COMM_BAND
  public :: CHEM_COMM_LOR
  public :: CHEM_COMM_BOR
  public :: CHEM_COMM_LXOR
  public :: CHEM_COMM_BXOR
  public :: CHEM_COMM_MINLOC
  public :: CHEM_COMM_MAXLOC
  public :: CHEM_COMM_REPLACE
  public :: CHEM_COMM_NO_OP

  public :: chem_comm_abort
  public :: chem_comm_allgather
  public :: chem_comm_bcast
  public :: chem_comm_create
  public :: chem_comm_get
  public :: chem_comm_init
  public :: chem_comm_inquire
  public :: chem_comm_isroot
  public :: chem_comm_log
  public :: chem_comm_reduce
  public :: chem_comm_set
!  public :: chem_comm_reducetile_pushwithhalo

contains

  subroutine chem_comm_init(rc, comm, isolate)
    integer,           intent(out) :: rc
    integer, optional, intent(in)  :: comm
    logical, optional, intent(in)  :: isolate

    ! -- local variables
    integer :: ierr
    logical :: flag

    ! -- begin
    rc = CHEM_RC_FAILURE

    call mpi_initialized(flag, ierr)
    if (.not.flag) then
      call mpi_init(ierr)
      if (ierr /= MPI_SUCCESS) return
    end if

    flag = .false.
    if (present(isolate)) flag = isolate
      
    if (present(comm)) then
      if (flag) then
        call mpi_comm_dup(comm, mpi_comm_chem, ierr)
        if (ierr /= MPI_SUCCESS) return
      else
        mpi_comm_chem = comm
      end if
    end if

    ! -- change MPI default error handler to return
    call mpi_errhandler_set(mpi_comm_chem, MPI_ERRORS_RETURN, ierr)
    if (ierr /= MPI_SUCCESS) return

    ! -- get group handle
    call mpi_comm_group(mpi_comm_chem, mpi_group_chem, ierr)
    if (ierr /= MPI_SUCCESS) return

    rc = CHEM_RC_SUCCESS
   
  end subroutine chem_comm_init

  logical function chem_comm_isroot()
    ! -- local variables
    integer :: ierr, rank
    ! -- begin
    chem_comm_isroot = .false.
    call mpi_comm_rank(mpi_comm_chem, rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    
    chem_comm_isroot = (rank == chem_comm_rootpe)

  end function chem_comm_isroot

  subroutine chem_comm_get(localpe, pecount, comm, group, rc)
    integer, optional, intent(out) :: localpe
    integer, optional, intent(out) :: pecount
    integer, optional, intent(out) :: comm
    integer, optional, intent(out) :: group
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    if (present(localpe)) then
      localpe = -1
      call mpi_comm_rank(mpi_comm_chem, localpe, ierr) 
      if (chem_rc_test((ierr /= MPI_SUCCESS), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (present(pecount)) then
      pecount = -1
      call mpi_comm_size(mpi_comm_chem, pecount, ierr) 
      if (chem_rc_test((ierr /= MPI_SUCCESS), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (present(comm)) comm = mpi_comm_chem

    if (present(group)) then
      call mpi_comm_group(mpi_comm_chem, group, ierr)
      if (chem_rc_test((ierr /= MPI_SUCCESS), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if

  end subroutine chem_comm_get

  subroutine chem_comm_inquire(comm, localpe, pecount, rc)
    integer,           intent(in)  :: comm
    integer, optional, intent(out) :: localpe
    integer, optional, intent(out) :: pecount
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    if (present(localpe)) then
      localpe = -1
      call mpi_comm_rank(comm, localpe, localrc)
      if (chem_rc_test((localrc /= MPI_SUCCESS), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (present(pecount)) then
      pecount = -1
      call mpi_comm_size(comm, pecount, localrc)
      if (chem_rc_test((localrc /= MPI_SUCCESS), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if

  end subroutine chem_comm_inquire

  subroutine chem_comm_set(comm, rc)
    integer, optional, intent(in)  :: comm
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    if (present(comm)) then
      mpi_comm_chem = comm
      ! -- update group handle
      call mpi_comm_group(mpi_comm_chem, mpi_group_chem, localrc)
      if (chem_rc_test((localrc /= MPI_SUCCESS), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if

  end subroutine chem_comm_set

  subroutine chem_comm_log(msg, sync)
    character(len=*),  intent(in) :: msg
    logical, optional, intent(in) :: sync

    write(0,'("chem_comm_log:",a)') trim(msg)
    if (present(sync)) then
      if (sync) call flush(0)
    end if

  end subroutine chem_comm_log

  subroutine chem_comm_abort(errcode, msg)
    integer,          optional, intent(in) :: errcode
    character(len=*), optional, intent(in) :: msg

    ! -- local variables
    integer :: ierr, localerrcode

    ! -- begin
    localerrcode = CHEM_RC_FAILURE
    if (present(msg)) write(0,'("chem_comm_abort:",a)') trim(msg)
    if (present(errcode)) localerrcode = errcode
    call mpi_abort(mpi_comm_chem, localerrcode, ierr)

  end subroutine chem_comm_abort

  ! -- allgather

  subroutine chem_comm_allgather_i1(sendbuf, recvbuf, count, comm, rc)
    integer,           intent(in)    :: sendbuf(:)
    integer,           intent(inout) :: recvbuf(:)
    integer, optional, intent(in)    :: count
    integer, optional, intent(in)    :: comm
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: localcomm, localcount

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localcomm = mpi_comm_chem
    if (present(comm)) localcomm = comm

    if (present(count)) then
      if (count > size(sendbuf)) return
      localcount = count
    else
      localcount = size(sendbuf)
    end if

    call mpi_allgather(sendbuf, localcount, MPI_INTEGER, recvbuf, localcount, &
      MPI_INTEGER, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_allgather_i1

  ! -- reduce

  subroutine chem_comm_reduce_r1(sendbuf, recvbuf, op, rootpe, comm, rc)
    real(CHEM_KIND_R4), intent(in)    :: sendbuf(:)
    real(CHEM_KIND_R4), intent(inout) :: recvbuf(:)
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

  end subroutine chem_comm_reduce_r1

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

  subroutine chem_comm_reduce_r3(sendbuf, recvbuf, op, rootpe, comm, rc)
    real(CHEM_KIND_R4), intent(in)    :: sendbuf(:,:,:)
    real(CHEM_KIND_R4), intent(inout) :: recvbuf(:,:,:)
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

  end subroutine chem_comm_reduce_r3

  ! -- communicators and groups

  subroutine chem_comm_create_comm(newcomm, color, comm, rc)
    integer,           intent(out) :: newcomm
    integer,           intent(in)  :: color
    integer, optional, intent(in)  :: comm
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: localcomm

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localcomm  = mpi_comm_chem
    if (present(comm)) localcomm = comm

    call mpi_comm_split(localcomm, color, 0, newcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_comm_create_comm

  subroutine chem_comm_create_group(newcomm, peList, comm, rc)
    integer,           intent(out) :: newcomm
    integer,           intent(in)  :: peList(:)
    integer, optional, intent(in)  :: comm
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: localcomm, localgroup, newgroup

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localcomm  = mpi_comm_chem
    localgroup = mpi_group_chem
    if (present(comm)) then
      localcomm = comm
      call mpi_comm_group(localcomm, localgroup, localrc)
      if (chem_rc_test((localrc /= MPI_SUCCESS), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    end if

    call mpi_group_incl(localgroup, size(peList), peList, newgroup, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return

    call mpi_comm_create_group(localcomm, newgroup, 0, newcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_create_group

  ! -- broadcast

  subroutine chem_comm_bcast_i0(data, rootpe, comm, rc)
    integer,           intent(inout) :: data
    integer, optional, intent(in)    :: rootpe
    integer, optional, intent(in)    :: comm
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: buffer(1)

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS
    buffer(1) = data
    call chem_comm_bcast_i1(buffer, rootpe=rootpe, comm=comm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    data = buffer(1)
    
  end subroutine chem_comm_bcast_i0

  subroutine chem_comm_bcast_i1(buffer, count, rootpe, comm, rc)
    integer,           intent(inout) :: buffer(:)
    integer, optional, intent(in)    :: count
    integer, optional, intent(in)    :: rootpe
    integer, optional, intent(in)    :: comm
    integer, optional, intent(out)   :: rc

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
    call mpi_bcast(buffer, localcount, MPI_INTEGER, root, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_bcast_i1

  subroutine chem_comm_bcast_r0(data, rootpe, comm, rc)
    real(CHEM_KIND_R4), intent(inout) :: data
    integer, optional,  intent(in)    :: rootpe
    integer, optional,  intent(in)    :: comm
    integer, optional,  intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    real    :: buffer(1)

    ! -- begin
    buffer(1) = data
    call chem_comm_bcast_r1(buffer, rootpe=rootpe, comm=comm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_bcast_r0

  subroutine chem_comm_bcast_r1(buffer, count, rootpe, comm, rc)
    real(CHEM_KIND_R4), intent(inout) :: buffer(:)
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
    call mpi_barrier(localcomm,localrc) ! ajl
    call mpi_bcast(buffer, localcount, MPI_REAL, root, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_bcast_r1
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
    call mpi_barrier(localcomm,localrc) ! ajl
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
    call mpi_barrier(localcomm,localrc) ! ajl
    call mpi_bcast(buffer, localcount, MPI_DOUBLE, root, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'bottom bcast_d2',maxval(buffer),minval(buffer),shape(buffer)
!    call flush(6)
    
  end subroutine chem_comm_bcast_d2

  subroutine chem_comm_bcast_s0(data, rootpe, comm, rc)
    character(len=*),  intent(inout) :: data
    integer, optional, intent(in)    :: rootpe
    integer, optional, intent(in)    :: comm
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: localcomm, root

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS
    root = chem_comm_rootpe
    if (present(rootpe)) root = rootpe
    localcomm = mpi_comm_chem
    if (present(comm)) localcomm = comm
    call mpi_bcast(data, len(data), MPI_CHARACTER, root, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_bcast_s0

  subroutine chem_comm_bcast_s1(buffer, count, rootpe, comm, rc)
    character(len=*),  intent(inout) :: buffer(:)
    integer, optional, intent(in)    :: count
    integer, optional, intent(in)    :: rootpe
    integer, optional, intent(in)    :: comm
    integer, optional, intent(out)   :: rc

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
    call mpi_bcast(buffer, localcount*len(buffer(1)), MPI_CHARACTER, root, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_bcast_s1

#if 0
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
#endif

end module chem_comm_mod
