module raqmschem_comm_mod
  use mpi 
  use chem_types_mod
  use chem_rc_mod
!  use chem_comm_mod, only : chem_comm_get,CHEM_COMM_SUM
  use chem_comm_mod, only : chem_comm_rootpe
  use raqmschem_model_mod, only : chem_model_get => raqmschem_model_get,chem_model_domain_get
  implicit none
  integer :: mpi_comm_chem    = MPI_COMM_NULL
  integer :: bcast_comm = MPI_COMM_NULL
  integer :: tilecomm = MPI_COMM_NULL
  logical :: usechemcomm
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
  interface raqmschem_comm_all_bcast_tile
    module procedure chem_comm_all_bcast_char_tile
    module procedure chem_comm_all_bcast_char1_tile
    module procedure chem_comm_all_bcast_char2_tile
    module procedure chem_comm_all_bcast_0dr4_tile
    module procedure chem_comm_all_bcast_1dr4_tile
    module procedure chem_comm_all_bcast_2dr4_tile
    module procedure chem_comm_all_bcast_3dr4_tile
    module procedure chem_comm_all_bcast_4dr4_tile
    module procedure chem_comm_all_bcast_0dr8_tile
    module procedure chem_comm_all_bcast_1dr8_tile
    module procedure chem_comm_all_bcast_2dr8_tile
    module procedure chem_comm_all_bcast_3dr8_tile
    module procedure chem_comm_all_bcast_4dr8_tile
    module procedure chem_comm_all_bcast_0di4_tile
    module procedure chem_comm_all_bcast_1di4_tile
    module procedure chem_comm_all_bcast_log_tile
  end interface raqmschem_comm_all_bcast_tile
  interface raqmschem_comm_bcast
    module procedure chem_comm_bcast_r2
    module procedure chem_comm_bcast_d2
  end interface raqmschem_comm_bcast
  interface raqmschem_comm_reduce
    module procedure chem_comm_reduce_r0
    module procedure chem_comm_reduce_r1
    module procedure chem_comm_reduce_r2
    module procedure chem_comm_reduce_d2
  end interface raqmschem_comm_reduce

  private
  public :: chem_reducetile_pushwithhalo
  public :: raqmschem_comm_abort
  public :: raqmschem_comm_all_bcast,setbcastcomm
  public :: raqmschem_comm_bcaSt
  public :: raqmschem_comm_reduce
  public :: raqmschem_comm_all_bcast_tile

contains

  subroutine setbcastcomm
  use raqmschem_pmgrid_mod, only : iam,tile
  integer ierr
  call chem_model_get(modelcomm=bcast_comm,tilecomm=tilecomm,rc=ierr)
  if(usechemcomm)then
    bcast_comm=bcast_comm
  else
    bcast_comm=mpi_comm_world
  endif
  return
  end subroutine setbcastcomm
  subroutine chem_comm_all_bcast_char(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb
  Character(*)  buffer
  lenb=len(buffer)
  call mpi_bcast(buffer,lenb,MPI_CHARACTER,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char
  subroutine chem_comm_all_bcast_char1(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb,lena
  Character(*)  buffer(:)
  lenb=len(buffer)
  lena=size(buffer)
  call mpi_bcast(buffer,lenb*lena,MPI_CHARACTER,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char1
  subroutine chem_comm_all_bcast_char2(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb,lena
  Character(*)  buffer(:,:)
  lenb=len(buffer)
  lena=size(buffer)
  call mpi_bcast(buffer,lenb*lena,MPI_CHARACTER,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char2
  subroutine chem_comm_all_bcast_log(buffer,rc)
  integer, optional :: rc
  logical :: buffer
  integer localrc
  call mpi_bcast(buffer,1,MPI_LOGICAL,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_log
  subroutine chem_comm_all_bcast_0dr4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4) :: buffer
!  write(6,*)iam,'bcast_1d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,1,MPI_REAL,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0dr4
  subroutine chem_comm_all_bcast_1dr4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4),dimension(:) :: buffer
!  write(6,*)iam,'bcast_1d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1dr4
  subroutine chem_comm_all_bcast_2dr4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4),dimension(:,:) :: buffer
!  write(6,*)iam,'bcast_2d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_2dr4
  subroutine chem_comm_all_bcast_3dr4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R4),dimension(:,:,:) :: buffer
!  mpi_comm_chem=0
!  tilecomm=0
!  call chem_model_get(modelcomm=mpi_comm_chem,tilecomm=tilecomm,rc=ierr)
!  write(6,*)'tile',tile,'mpi_comm_chem',mpi_comm_chem,'bcast_comm',mpi_comm_world,'tilec',tilecomm
!  flush(6)
!  write(6,*)iam,'bcast_3d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,bcast_comm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_3dr4
  subroutine chem_comm_all_bcast_4dr4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R4),dimension(:,:,:,:) :: buffer
  write(6,*)iam,'bcast_4d ',shape(buffer)
  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_4dr4
  subroutine chem_comm_all_bcast_0di4(buffer,rc)
  integer, optional :: rc
  integer localrc
  integer :: buffer
  call mpi_bcast(buffer,1,MPI_INTEGER,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0di4
  subroutine chem_comm_all_bcast_1di4(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  integer :: buffer(:)
  call mpi_bcast(buffer,size(buffer),MPI_INTEGER,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1di4
  subroutine chem_comm_all_bcast_0dr8(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R8) :: buffer
  call mpi_bcast(buffer,1,MPI_DOUBLE,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0dr8
  subroutine chem_comm_all_bcast_1dr8(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R8),dimension(:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1dr8
  subroutine chem_comm_all_bcast_2dr8(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R8),dimension(:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_2dr8
  subroutine chem_comm_all_bcast_3dr8(buffer,rc)
  integer, optional :: rc
  integer localrc
  real(CHEM_KIND_R8),dimension(:,:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,0,bcast_comm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_3dr8
  subroutine chem_comm_all_bcast_log_tile(buffer,rc)
  integer, optional :: rc
  logical :: buffer
  integer localrc
  call mpi_bcast(buffer,1,MPI_LOGICAL,0,tilecomm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_log_tile
  subroutine chem_comm_all_bcast_char_tile(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb
  Character(*)  buffer
  lenb=len(buffer)
  call mpi_bcast(buffer,lenb,MPI_CHARACTER,0,tilecomm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char_tile
  subroutine chem_comm_all_bcast_char1_tile(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb,lena
  Character(*)  buffer(:)
  lenb=len(buffer)
  lena=size(buffer)
  call mpi_bcast(buffer,lenb*lena,MPI_CHARACTER,0,tilecomm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char1_tile
  subroutine chem_comm_all_bcast_char2_tile(buffer,rc)
  integer, optional :: rc
  integer localrc,lenb,lena
  Character(*)  buffer(:,:)
  lenb=len(buffer)
  lena=size(buffer)
  call mpi_bcast(buffer,lenb*lena,MPI_CHARACTER,0,tilecomm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_char2_tile
  subroutine chem_comm_all_bcast_0di4_tile(buffer,rc)
  integer, optional :: rc
  integer localrc
  integer :: buffer
  call mpi_bcast(buffer,1,MPI_INTEGER,0,tilecomm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0di4_tile
  subroutine chem_comm_all_bcast_1di4_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam
  integer, optional :: rc
  integer localrc
  integer :: buffer(:)
  call mpi_bcast(buffer,size(buffer),MPI_INTEGER,0,tilecomm,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1di4_tile
  subroutine chem_comm_all_bcast_0dr4_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R4) :: buffer
!  mpi_comm_chem=0
!  tilecomm=0
!  call chem_model_get(modelcomm=mpi_comm_chem,tilecomm=tilecomm,rc=ierr)
!  write(6,*)'tile',tile,'mpi_comm_chem',mpi_comm_chem,'bcast_comm',mpi_comm_world,'tilec',tilecomm
!  flush(6)
!  write(6,*)iam,'bcast_3d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,1,MPI_REAL,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0dr4_tile
  subroutine chem_comm_all_bcast_1dr4_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R4),dimension(:) :: buffer
!  mpi_comm_chem=0
!  tilecomm=0
!  call chem_model_get(modelcomm=mpi_comm_chem,tilecomm=tilecomm,rc=ierr)
!  write(6,*)'tile',tile,'mpi_comm_chem',mpi_comm_chem,'bcast_comm',mpi_comm_world,'tilec',tilecomm
!  flush(6)
!  write(6,*)iam,'bcast_3d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1dr4_tile
  subroutine chem_comm_all_bcast_2dr4_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R4),dimension(:,:) :: buffer
!  mpi_comm_chem=0
!  tilecomm=0
!  call chem_model_get(modelcomm=mpi_comm_chem,tilecomm=tilecomm,rc=ierr)
!  write(6,*)'tile',tile,'mpi_comm_chem',mpi_comm_chem,'bcast_comm',mpi_comm_world,'tilec',tilecomm
!  flush(6)
!  write(6,*)iam,'bcast_3d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_2dr4_tile
  subroutine chem_comm_all_bcast_3dr4_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R4),dimension(:,:,:) :: buffer
!  mpi_comm_chem=0
!  tilecomm=0
!  call chem_model_get(modelcomm=mpi_comm_chem,tilecomm=tilecomm,rc=ierr)
!  write(6,*)'tile',tile,'mpi_comm_chem',mpi_comm_chem,'bcast_comm',mpi_comm_world,'tilec',tilecomm
!  flush(6)
!  write(6,*)iam,'bcast_3d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_3dr4_tile
  subroutine chem_comm_all_bcast_4dr4_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R4),dimension(:,:,:,:) :: buffer
!  mpi_comm_chem=0
!  tilecomm=0
!  call chem_model_get(modelcomm=mpi_comm_chem,tilecomm=tilecomm,rc=ierr)
!  write(6,*)'tile',tile,'mpi_comm_chem',mpi_comm_chem,'bcast_comm',mpi_comm_world,'tilec',tilecomm
!  flush(6)
!  write(6,*)iam,'bcast_3d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_REAL,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_4dr4_tile
  subroutine chem_comm_all_bcast_0dr8_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R8) :: buffer
!  mpi_comm_chem=0
!  tilecomm=0
!  call chem_model_get(modelcomm=mpi_comm_chem,tilecomm=tilecomm,rc=ierr)
!  write(6,*)'tile',tile,'mpi_comm_chem',mpi_comm_chem,'bcast_comm',mpi_comm_world,'tilec',tilecomm
!  flush(6)
!  write(6,*)iam,'bcast_3d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,1,MPI_DOUBLE,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_0dr8_tile
  subroutine chem_comm_all_bcast_1dr8_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R8),dimension(:) :: buffer
!  mpi_comm_chem=0
!  tilecomm=0
!  call chem_model_get(modelcomm=mpi_comm_chem,tilecomm=tilecomm,rc=ierr)
!  write(6,*)'tile',tile,'mpi_comm_chem',mpi_comm_chem,'bcast_comm',mpi_comm_world,'tilec',tilecomm
!  flush(6)
!  write(6,*)iam,'bcast_3d ',shape(buffer)
!  flush(6)
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_1dr8_tile
  subroutine chem_comm_all_bcast_2dr8_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R8),dimension(:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_2dr8_tile
  subroutine chem_comm_all_bcast_3dr8_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R8),dimension(:,:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_3dr8_tile
  subroutine chem_comm_all_bcast_4dr8_tile(buffer,rc)
  use raqmschem_pmgrid_mod, only : iam,tile,masterproct
  integer, optional :: rc
  integer localrc,ierr,mpi_comm_chem
  real(CHEM_KIND_R8),dimension(:,:,:,:) :: buffer
  call mpi_bcast(buffer,size(buffer),MPI_DOUBLE,chem_comm_rootpe,tilecomm,localrc)
!  call mpi_bcast(buffer,size(buffer),MPI_REAL,0,mpi_comm_chem,localrc)
  if(present(rc))rc=localrc
  return
  end subroutine chem_comm_all_bcast_4dr8_tile
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

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buf2d = 0._CHEM_KIND_R4

    buf2d(ids:ide, jds:jde) = datain
    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4
!   do reduction
    call raqmschem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    call raqmschem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
    datawithhalo=recvbuf(ihs:ihe,jhs:jhe)
    deallocate (buf2d,recvbuf)
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

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buf2d = 0._CHEM_KIND_R4

    buf2d(ids:ide, jds:jde) = datain
!   allocate reduce bufr
    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4
!   do reduction
    call raqmschem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!   now need to bcast and then put out data with tile inside halo included
    call raqmschem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
    datawithhalo=recvbuf(ihs:ihe,jhs:jhe)
    deallocate (buf2d,recvbuf)
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
!     allocate reduce bufr

      recvbuf = 0._CHEM_KIND_R4
!     do reduction
      call raqmschem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!     now need to bcast and then put out data with tile inside halo included
      call raqmschem_comm_bcast(recvbuf, comm=tileComm, rc=localrc)
      datawithhalo(:,:,k)=recvbuf(ihs:ihe,jhs:jhe)
    end do
    deallocate (buf2d,recvbuf)
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
    integer :: localrc,ierr
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
    call mpi_bcast(buffer, localcount, MPI_DOUBLE, root, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
    
  end subroutine chem_comm_bcast_d2
  subroutine chem_comm_reduce_r0(sendbuf, recvbuf, op, rootpe, comm, rc)
    real(CHEM_KIND_R4), intent(in)    :: sendbuf
    real(CHEM_KIND_R4), intent(inout) :: recvbuf
    integer,            intent(in)    :: op
    integer, optional,  intent(in)    :: rootpe
    integer, optional,  intent(in)    :: comm
    integer, optional,  intent(out)   :: rc

    ! -- local variables
    integer :: localrc,ierr
    integer :: localcomm, localroot

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(modelcomm=mpi_comm_chem,rc=ierr)
    localcomm = mpi_comm_chem
!    localcomm=mpi_comm_world
    if (present(comm)) localcomm = comm

    localroot = chem_comm_rootpe
    if (present(rootpe)) localroot = rootpe

    call mpi_reduce(sendbuf, recvbuf, 1 , MPI_REAL, op, localroot, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_comm_reduce_r0
  subroutine chem_comm_reduce_r1(sendbuf, recvbuf, op, rootpe, comm, rc)
    real(CHEM_KIND_R4), intent(in)    :: sendbuf(:)
    real(CHEM_KIND_R4), intent(inout) :: recvbuf(:)
    integer,            intent(in)    :: op
    integer, optional,  intent(in)    :: rootpe
    integer, optional,  intent(in)    :: comm
    integer, optional,  intent(out)   :: rc

    ! -- local variables
    integer :: localrc,ierr
    integer :: localcomm, localroot

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(modelcomm=mpi_comm_chem,rc=ierr)
    localcomm = mpi_comm_chem
!    localcomm=mpi_comm_world
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
    integer :: localrc,ierr
    integer :: localcomm, localroot

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(modelcomm=mpi_comm_chem,rc=ierr)
    localcomm = mpi_comm_chem
!    localcomm=mpi_comm_world
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
    integer :: localrc,ierr
    integer :: localcomm, localroot

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(modelcomm=mpi_comm_chem,rc=ierr)
    localcomm = mpi_comm_chem
    if (present(comm)) localcomm = comm

    localroot = chem_comm_rootpe
    if (present(rootpe)) localroot = rootpe

    call mpi_reduce(sendbuf, recvbuf, size(sendbuf), MPI_DOUBLE, op, localroot, localcomm, localrc)
    if (chem_rc_test((localrc /= MPI_SUCCESS), &
      file=__FILE__, line=__LINE__, rc=rc)) return
   end subroutine chem_comm_reduce_d2
end module raqmschem_comm_mod
