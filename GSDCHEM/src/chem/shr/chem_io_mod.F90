module chem_io_mod

  use mpi
  use chem_rc_mod
  use chem_types_mod
  use chem_comm_mod
  use chem_model_mod

  implicit none

  integer, parameter :: ioUnit = 100

  logical :: chem_io_verbose = .false.

  interface chem_io_read
    module procedure chem_io_read_2DR4
    module procedure chem_io_read_3DR4
    module procedure chem_io_read_4DR4
  end interface chem_io_read

  interface chem_io_write
    module procedure chem_io_write_2DR4
    module procedure chem_io_write_2DR8
    module procedure chem_io_write_3DR4
    module procedure chem_io_write_3DR8
    module procedure chem_io_write_4DR4
  end interface chem_io_write

  private

  public :: chem_io_init
  public :: chem_io_read
  public :: chem_io_write
  public :: inquire_file_var
  public :: inquire_file_dim

contains

  subroutine chem_io_init(verbose, rc)
    logical, optional, intent(in)  :: verbose
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount, tile, tileCount
    integer :: pe, peCount
    integer :: i, localpe, npe
    integer :: comm, tileComm
    integer, dimension(:), allocatable :: localTile, tileToPet, pes

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    chem_io_verbose = .false.
    if (present(verbose)) chem_io_verbose = verbose

    call chem_model_get(deCount=deCount, tileCount=tileCount, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- if no model on this PET, bail out
    if (deCount < 1) return

    call chem_comm_get(comm=comm, localpe=localpe, pecount=peCount, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(localTile(tileCount), tileToPet(tileCount*peCount), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Unable to allocate memory", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- store which tiles are assigned to this PET
    localTile = -1
    do de = 0, deCount-1
      call chem_model_get(de=de, tile=tile, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      localTile(tile) = localpe
    end do

    ! -- build a global tile-to-PET map
    call chem_comm_allgather(localTile, tileToPet, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    deallocate(localTile, stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Unable to free memory", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- extract the list of PETs assigned to each tile and create MPI groups
    allocate(pes(peCount), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Unable to allocate memory", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- gather PET list for each tile and create tile-specific communicator
    pes = -1
    do tile = 1, tileCount
      npe = 0
      do i = tile, tileCount*peCount, tileCount
        if (tileToPet(i) > -1) then
          npe = npe + 1
          pes(npe) = tileToPet(i)
        end if
      end do

      ! -- create new communicator
      call chem_comm_create(tileComm, pes(1:npe), rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      do de = 0, deCount-1
        call chem_model_get(de=de, tile=i, rc=localrc)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        if (tile == i) then
          call chem_model_set(de=de, tileComm=tileComm, rc=localrc)
          if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        end if
      end do
    end do

    deallocate(pes, tileToPet, stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Unable to free memory", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- store I/O layout into model
    do de = 0, deCount-1
      call chem_model_get(de=de, tile=tile, tileComm=tileComm, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      ! -- get local PET for tile 
      call chem_comm_inquire(tileComm, localpe=pe, pecount=npe, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      ! -- mark local root PET as I/O PET
      call chem_model_set(de=de, localIOflag=(pe == 0), rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      if (chem_io_verbose) &
        write(6,'("chem_io_init: PET:",i2," DE:",i02," tile=",i0," - comm=",i0," PE:",i0,"/",i0)') &
          localpe, de, tile, tileComm, pe, npe
    end do

  end subroutine chem_io_init


  

  subroutine chem_io_file_name(fullname, filename, tile, pathname,doread)
    character(len=*),           intent(out) :: fullname
    character(len=*),           intent(in)  :: filename
    integer,                    intent(in)  :: tile
    character(len=*), optional, intent(in)  :: pathname
    logical, optional,intent(in)            :: doread

    ! -- local variables
    integer :: lstr,lstr2
    character(len=CHEM_MAXSTR) :: fname
    logical exist


    ! -- begin
    fname = ""
    fullname = ""

    lstr = len_trim(filename)
    if (lstr > 4) then
      if (filename(lstr-3:lstr) == ".dat") then
        write(fname, '("tile",i0,"/",a)') tile, trim(filename)
      elseif(filename(lstr-2:lstr) == '.nc')then
         write(fname,'(a,".tile",i0,".nc")')filename(1:lstr-3),tile
      else
        write(fname, '(a,".tile",i0,".dat")') trim(filename), tile
      end if
    else
      write(fname, '(a,".tile",i0,".dat")') trim(filename), tile
    end if

    if (present(pathname)) then
      lstr = len_trim(pathname)
      if (lstr > 0) then
!        write(6,*)'pathname',trim(pathname)
!        write(6,*)'lstr',lstr
!        write(6,*)'len',len(fullname)
!        call flush(6)
!        write(6,*)'fname',trim(fname)
!        call flush(6)
!        write(6,*)'lstr',pathname(lstr:lstr)
!        call flush(6)
        if (pathname(lstr:lstr) == "/") then
!          write(6,*)'176'
!          call flush(6)
!          write(6,*)'lenpathname',len_trim(pathname)
!          write(6,*)'lenfname',len_trim(fname)
!          call flush(6)
          fullname = trim(pathname) // trim(fname)
!          write(6,*)'fullname',len_trim(fullname),trim(fullname)
!          call flush(6)
        else
!          write(6,*)'180'
!          call flush(6)
          fullname = trim(pathname) // "/" // trim(fname)
        end if
      else
        fullname = trim(fname)
      end if
    else
      fullname = trim(fname)
    end if
    
    if(present(doread))then
      inquire(file=fullname,exist=exist)
      if(exist)return
!     if not there see if netcdf file is    
      lstr2=len_trim(fullname)
      if (fullname(lstr2-3:lstr2) == '.dat')then
        fullname(lstr2-3:lstr2)='.nc '
        inquire(file=fullname,exist=exist)
!       write(6,*)'inquire 2 ',exist,'fullname',trim(fullname)
        if(exist)return
        write(6,*)'gsdchemneither .dat nor .nc are there',fullname,'fname',fname
        call killit('neither')
      elseif(fullname(lstr2-2:lstr2)=='.nc')then
        write(6,*)'.nc not there',fullname
        call killit('.nc')
      endif
    endif
    return

  end subroutine chem_io_file_name

  subroutine chem_io_file_read(datafile, buffer, recrange, recsize, recstride, rc)
    character(len=*),   intent(in)  :: datafile
    real(CHEM_KIND_R4), intent(out) :: buffer(:)
    integer, optional,  intent(in)  :: recrange(2)
    integer, optional,  intent(in)  :: recsize
    integer, optional,  intent(in)  :: recstride
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: irec, is, ie
    integer :: rcount, rsize, rstride
    integer :: rrange(2)

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    buffer = 0._CHEM_KIND_R4

    if (present(recrange)) then
      rrange = recrange
    else
      rrange = 1
    end if
    rcount = rrange(2) - rrange(1) + 1

    if (present(recsize)) then
      rsize = recsize
    else
      rsize = size(buffer) / rcount
    end if

    if (present(recstride)) then
      rstride = max(rsize, recstride)
    else
      rstride = rsize
    end if

    if (chem_rc_test((size(buffer) < rcount * rstride), &
        msg="insufficient buffer size", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- open file
    open(unit=ioUnit, file=trim(datafile), form='unformatted', action='read', position='rewind', iostat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Failure opening file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- advance to first record
    do irec = 2, rrange(1)
      read(unit=ioUnit, iostat=localrc)
      if (chem_rc_test((localrc /= 0), msg="Unable to locate record in file: "//trim(datafile), &
          file=__FILE__, line=__LINE__, rc=rc)) then
        close(unit=iounit)
        return
      end if
    end do

    ! -- read records
    is = 1
    ie = rsize
    do irec = 1, rcount
      read(unit=iounit, iostat=localrc) buffer(is:ie)
      if (chem_rc_test((localrc /= 0), msg="Failure reading data from file: "//trim(datafile), &
          file=__FILE__, line=__LINE__, rc=rc)) then
        close(unit=iounit)
        return
      end if
      is = is + rstride
      ie = ie + rstride
    end do

    ! -- close file
    close(unit=ioUnit)

  end subroutine chem_io_file_read


  subroutine chem_io_file_write(datafile, buffer, pos, reccount, recsize, recstride, rc)
    character(len=*),           intent(in)  :: datafile
    real(CHEM_KIND_R4),         intent(in)  :: buffer(:)
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(in)  :: reccount
    integer,          optional, intent(in)  :: recsize
    integer,          optional, intent(in)  :: recstride
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: irec, is, ie
    integer :: rcount, rsize, rstride
    character(len=6) :: filepos

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    filepos = 'rewind'
    if (present(pos)) then
      select case (trim(pos))
        case ('a', 'append')
          filepos = 'append'
        case default
          filepos = 'rewind'
      end select
    end if

    if (present(reccount)) then
      rcount = reccount
    else
      rcount = 1
    end if

    if (present(recsize)) then
      rsize = recsize
    else
      rsize = size(buffer) / rcount
    end if

    if (present(recstride)) then
      rstride = max(rsize, recstride)
    else
      rstride = rsize
    end if

    ! -- adjust record count
    rcount = size(buffer) / rstride

    if (chem_rc_test((size(buffer) < rcount * rstride), &
        msg="insufficient buffer size", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- open file
    open(unit=ioUnit, file=trim(datafile), form='unformatted', action='write', &
      position=trim(filepos), iostat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Failure opening file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- write records
    is = 1
    ie = rsize
    do irec = 1, rcount
      write(unit=ioUnit, iostat=localrc) buffer(is:ie)
      if (chem_rc_test((localrc /= 0), msg="Failure writing data to file: "//trim(datafile), &
          file=__FILE__, line=__LINE__, rc=rc)) then
        close(unit=iounit)
        return
      end if
      is = is + rstride
      ie = ie + rstride
    end do

    ! -- close file
    close(unit=ioUnit)

  end subroutine chem_io_file_write




  subroutine chem_io_read_2DR4(filename, farray, path, recrange, recsize, recstride, de, varname, rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(out) :: farray(:,:)
    character(len=*), optional, intent(in)  :: path
    integer,          optional, intent(in)  :: recrange(2)
    integer,          optional, intent(in)  :: recsize
    integer,          optional, intent(in)  :: recstride
    integer,          optional, intent(in)  :: de
    character(len=*), optional, intent(in)  :: varname
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: ids, ide, jds, jde, its, ite, jts, jte
    integer :: bsize(2)
    logical :: localIOflag
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:), allocatable :: buffer
    real(CHEM_KIND_R4), dimension(:,:), allocatable, target :: buf2d
    logical :: doread=.true.
    integer lenstr,bsiz1d(1)

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- check size consistency
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    bsize = (/ ite-its+1, jte-jts+1 /)
    bsiz1d(1)=bsize(1)*bsize(2)
    allocate(buffer(bsize(1)*bsize(2)), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buffer = 0._CHEM_KIND_R4
    allocate(buf2d(its:ite,jts:jte), stat=localrc)

    if (localIOflag) then

      call chem_io_file_name(datafile, filename, tile, pathname=path,doread=doread)
      lenstr=len_trim(datafile)
!      write(6,*)'lenstr',lenstr,'datafile',trim(datafile)
!      call flush(6)
      if(datafile(lenstr-2:lenstr)=='.nc')then
        if(present(varname))then
          call chem_io_file_read_2d_nc(datafile, buf2d, recrange=recrange, &
          recsize=recsize, recstride=recstride, varname=varname,rc=localrc)
          buffer=reshape(buf2d,bsiz1d)
        else 
          write(6,*)'netcdf file without varname',trim(filename)
        endif
      else 
        call chem_io_file_read(datafile, buffer, recrange=recrange, &
        recsize=recsize, recstride=recstride, rc=localrc)
      endif
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (chem_io_verbose) &
        write(6,'("chem_io_read: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
          trim(datafile), minval(buffer), maxval(buffer)
    end if

    call chem_comm_bcast(buffer, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    buf2d = reshape(buffer, bsize)

    farray = buf2d(ids:ide, jds:jde) 
     
    deallocate(buffer, buf2d, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_read_2DR4


  subroutine chem_io_read_3DR4(filename, farray, path, recrange, recsize, recstride, de, varname, rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(out) :: farray(:,:,:)
    character(len=*), optional, intent(in)  :: path
    integer,          optional, intent(in)  :: recrange(2)
    integer,          optional, intent(in)  :: recsize
    integer,          optional, intent(in)  :: recstride
    integer,          optional, intent(in)  :: de
    character(len=*), optional, intent(in)  :: varname
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: ids, ide, jds, jde, its, ite, jts, jte
    integer :: bsize(3)
    logical :: localIOflag
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:),     allocatable :: buffer
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable, target :: buf3d
    integer bsiz1d(1),lenstr
    logical :: doread=.true.
    integer localpe

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS
    call chem_comm_get(localpe=localpe)
!    if(present(varname))then
!    write(200+localpe,*)'top read3dr4',trim(varname),'size',size(farray)
!    write(200+localpe,*)'shape',shape(farray)
    !call flush(200+localpe)
!    endif

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    bsize = (/ ite-its+1, jte-jts+1, size(farray, dim=3) /)
!    write(200+localpe,*)'bsize',bsize
!    call flush(200+localpe)
    bsiz1d=bsize(1)*bsize(2)*bsize(3)

    ! -- check size consistency
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)*bsize(3)), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buffer(bsize(1)*bsize(2)*bsize(3)), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buffer = 0._CHEM_KIND_R4
    allocate(buf3d(its:ite,jts:jte,bsize(3)), stat=localrc)

    if (localIOflag) then
!      write(6,*)'len datafile',len(datafile)
!      call flush(6)
!      write(6,*)'len filename',len_trim(filename)
!      write(6,*)'len path',len_trim(path)
!      call flush(6)

      call chem_io_file_name(datafile, filename, tile, pathname=path,doread=doread)
      lenstr=len_trim(datafile)
!      write(6,*)'lenstr3d',lenstr,'datafile',trim(datafile)
!      call flush(6)
!      write(200+localpe,*)'datafile',trim(datafile)
!      call flush(200+localpe)
      if(datafile(lenstr-2:lenstr)=='.nc')then
        if(present(varname))then
          call chem_io_file_read_3d_nc(datafile, buf3d, recrange=recrange, &
          recsize=recsize, recstride=recstride, varname=varname,rc=localrc)
          buffer=reshape(buf3d,bsiz1d)
        else
          write(6,*)'netcdf file without varname',trim(filename)
        endif
      else

        call chem_io_file_read(datafile, buffer, recrange=recrange, &
        recsize=recsize, recstride=recstride, rc=localrc)
      endif
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (chem_io_verbose) &
        write(6,'("chem_io_read: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
          trim(datafile), minval(buffer), maxval(buffer)
    end if
!    write(200+localpe,*)'buffer3dr4',shape(buffer),'var',trim(varname)
!    call flush(200+localpe)

    call chem_comm_bcast(buffer, comm=tileComm, rc=localrc)
    !write(200+localpe,*)'localrc',localrc
!    call flush(200+localpe)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    buf3d = reshape(buffer, bsize)

    farray = buf3d(ids:ide, jds:jde, :)
     
    deallocate(buffer, buf3d, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
!    if(present(varname))then
!    write(200+localpe,*)'bottom read3dr4',trim(varname)
!    call flush(200+localpe)
    !endif

  end subroutine chem_io_read_3DR4

  subroutine chem_io_read_4DR4(filename, farray, order, path, recsize, recstride, de, rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(out) :: farray(:,:,:,:)
    character(len=*), optional, intent(in)  :: order
    character(len=*), optional, intent(in)  :: path
    integer,          optional, intent(in)  :: recsize
    integer,          optional, intent(in)  :: recstride
    integer,          optional, intent(in)  :: de
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: i, j, k, l, id, jd, lbuf
    integer :: ids, ide, jds, jde, its, ite, jts, jte, nk, nl
    logical :: localIOflag
    character(len=3) :: localOrder
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:),       allocatable :: buffer
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable, target :: buf4d

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    localOrder = "ijk"
    if (present(order)) localOrder = order

    select case (trim(localOrder))
      case("ikj")
        ! -- (i,k,j)
        nk = size(farray, dim=2)
      case default
        ! -- default to (i,j,k)
        nk = size(farray, dim=3)
    end select

    nl = size(farray, dim=4)

    ! -- check size consistency
    lbuf = (ide-ids+1) * (jde-jds+1) * nk * nl
    if (chem_rc_test((size(farray) /= lbuf), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    lbuf = (ite-its+1) * (jte-jts+1) * nk * nl
    allocate(buffer(lbuf), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    buffer = 0._CHEM_KIND_R4

    if (localIOflag) then

      call chem_io_file_name(datafile, filename, tile, pathname=path)

      print *,'READ4D: '//trim(filename)//" -> "//trim(datafile)

      call chem_io_file_read(datafile, buffer, recrange=(/ 1, nl /), &
        recsize=recsize, recstride=recstride, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (chem_io_verbose) &
        write(6,'("chem_io_read: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
          trim(datafile), minval(buffer), maxval(buffer)
    end if

    call chem_comm_bcast(buffer, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf4d(its:ite,jts:jte,nk,nl), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    buf4d = reshape(buffer, shape(buf4d))

    select case (trim(localOrder))
      case("ikj")
        do l = 1, nl
          do k = 1, nk
            j = 0
            do jd = jds, jde
              j = j + 1
              i = 0
              do id = ids, ide
                i = i + 1
                farray(i, k, j, l) = buf4d(id, jd, k, l)
              end do
            end do
          end do
        end do
      case default
        farray = buf4d(ids:ide, jds:jde, :, :)
    end select

    deallocate(buffer, buf4d, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_read_4DR4



  subroutine chem_io_write_2DR4(filename, farray, path, pos, de, rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(in)  :: farray(:,:)
    character(len=*), optional, intent(in)  :: path
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(in)  :: de
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: ids, ide, jds, jde, its, ite, jts, jte
    logical :: localIOflag
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:,:), allocatable, target :: buf2d, recvbuf

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- check size consistency
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buf2d = 0._CHEM_KIND_R4

    buf2d(ids:ide, jds:jde) = farray

    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4

    call chem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (localIOflag) then

      call chem_io_file_name(datafile, filename, tile, pathname=path)

      call chem_io_file_write(datafile, reshape(recvbuf, (/size(buf2d)/)), &
        pos=pos, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (chem_io_verbose) &
        write(6,'("chem_io_write: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
          trim(datafile), minval(recvbuf), maxval(recvbuf)
    end if

    deallocate(buf2d, recvbuf, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_write_2DR4


  subroutine chem_io_write_2DR8(filename, farray, path, pos, de, rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R8),         intent(in)  :: farray(:,:)
    character(len=*), optional, intent(in)  :: path
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(in)  :: de
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: ids, ide, jds, jde, its, ite, jts, jte
    logical :: localIOflag
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:,:), allocatable, target :: buf2d, recvbuf

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- check size consistency
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buf2d = 0._CHEM_KIND_R4

    buf2d(ids:ide, jds:jde) = real(farray, kind=CHEM_KIND_R4)

    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4

    call chem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (localIOflag) then

      call chem_io_file_name(datafile, filename, tile, pathname=path)

      call chem_io_file_write(datafile, reshape(recvbuf, (/size(buf2d)/)), &
        pos=pos, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (chem_io_verbose) &
        write(6,'("chem_io_write: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
          trim(datafile), minval(recvbuf), maxval(recvbuf)
    end if

    deallocate(buf2d, recvbuf, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_write_2DR8


  subroutine chem_io_write_3DR4(filename, farray, order, path, pos, de, rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(in)  :: farray(:,:,:)
    character(len=*), optional, intent(in)  :: order
    character(len=*), optional, intent(in)  :: path
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(in)  :: de
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: i, j, k, id, jd, lbuf
    integer :: ids, ide, jds, jde, its, ite, jts, jte, nk
    logical :: localIOflag
    character(len=3) :: localOrder
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:),     allocatable :: recvbuf
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: buf3d

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    localOrder = "ijk"
    if (present(order)) localOrder = order

    select case (trim(localOrder))
      case("ikj")
        ! -- (i,k,j)
        nk = size(farray,dim=2)
      case default
        ! -- default to (i,j,k)
        nk = size(farray,dim=3)
    end select

    ! -- check consistency in decomposition
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)*nk), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf3d(its:ite,jts:jte,nk), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot allocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    buf3d = 0._CHEM_KIND_R4

    select case (trim(localOrder))
      case("ikj")
        do k = 1, nk
          j = 0
          do jd = jds, jde
            j = j + 1
            i = 0
            do id = ids, ide
              i = i + 1
              buf3d(id, jd, k) = farray(i, k, j)
            end do
          end do
        end do
      case default
        buf3d(ids:ide, jds:jde, 1:nk) = farray
    end select

    lbuf = (ite-its+1)*(jte-jts+1)*nk
    allocate(recvbuf(lbuf), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot allocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4

    call chem_comm_reduce(reshape(buf3d, (/lbuf/)), recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (localIOflag) then

      call chem_io_file_name(datafile, filename, tile, pathname=path)

      call chem_io_file_write(datafile, recvbuf, pos=pos, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (chem_io_verbose) &
        write(6,'("chem_io_write: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
          trim(datafile), minval(recvbuf), maxval(recvbuf)
    end if

    deallocate(buf3d, recvbuf, stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot deallocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_write_3DR4

  subroutine chem_io_write_3DR8(filename, farray, order, path, pos, de, rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R8),         intent(in)  :: farray(:,:,:)
    character(len=*), optional, intent(in)  :: order
    character(len=*), optional, intent(in)  :: path
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(in)  :: de
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: i, j, k, id, jd, lbuf
    integer :: ids, ide, jds, jde, its, ite, jts, jte, nk
    logical :: localIOflag
    character(len=3) :: localOrder
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:),     allocatable :: recvbuf
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: buf3d

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    localOrder = "ijk"
    if (present(order)) localOrder = order

    select case (trim(localOrder))
      case("ikj")
        ! -- (i,k,j)
        nk = size(farray,dim=2)
      case default
        ! -- default to (i,j,k)
        nk = size(farray,dim=3)
    end select

    ! -- check consistency in decomposition
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)*nk), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf3d(its:ite,jts:jte,nk), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot allocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    buf3d = 0._CHEM_KIND_R4

    select case (trim(localOrder))
      case("ikj")
        do k = 1, nk
          j = 0
          do jd = jds, jde
            j = j + 1
            i = 0
            do id = ids, ide
              i = i + 1
              buf3d(id, jd, k) = real(farray(i, k, j), kind=CHEM_KIND_R4)
            end do
          end do
        end do
      case default
        buf3d(ids:ide, jds:jde, 1:nk) = real(farray, kind=CHEM_KIND_R4)
    end select

    lbuf = (ite-its+1)*(jte-jts+1)*nk
    allocate(recvbuf(lbuf), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot allocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4

    call chem_comm_reduce(reshape(buf3d, (/lbuf/)), recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (localIOflag) then

      call chem_io_file_name(datafile, filename, tile, pathname=path)

      call chem_io_file_write(datafile, recvbuf, pos=pos, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (chem_io_verbose) &
        write(6,'("chem_io_write: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
          trim(datafile), minval(recvbuf), maxval(recvbuf)
    end if

    deallocate(buf3d, recvbuf, stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot deallocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_write_3DR8


  subroutine chem_io_write_4DR4(filename, farray, order, path, pos, recsize, recstride, de, rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(in)  :: farray(:,:,:,:)
    character(len=*), optional, intent(in)  :: order
    character(len=*), optional, intent(in)  :: path
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(in)  :: recsize
    integer,          optional, intent(in)  :: recstride
    integer,          optional, intent(in)  :: de
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: i, j, k, l, id, jd, lbuf
    integer :: ids, ide, jds, jde, its, ite, jts, jte, nk, nl
    logical :: localIOflag
    character(len=3) :: localOrder
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:),       allocatable :: recvbuf
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: buf4d

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    localOrder = "ijk"
    if (present(order)) localOrder = order

    select case (trim(localOrder))
      case("ikj")
        ! -- (i,k,j)
        nk = size(farray, dim=2)
      case default
        ! -- default to (i,j,k)
        nk = size(farray, dim=3)
    end select

    nl = size(farray, dim=4)

    ! -- check consistency in decomposition
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)*nk*nl), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf4d(its:ite,jts:jte,nk,nl), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot allocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    buf4d = 0._CHEM_KIND_R4

    select case (trim(localOrder))
      case("ikj")
        do l = 1, nl
          do k = 1, nk
            j = 0
            do jd = jds, jde
              j = j + 1
              i = 0
              do id = ids, ide
                i = i + 1
                buf4d(id, jd, k, l) = farray(i, k, j, l)
              end do
            end do
          end do
        end do
      case default
        buf4d(ids:ide, jds:jde, 1:nk, 1:nl) = farray
    end select

    lbuf = (ite-its+1)*(jte-jts+1)*nk*nl

    allocate(recvbuf(lbuf), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot allocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4

    call chem_comm_reduce(reshape(buf4d, (/lbuf/)), recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (localIOflag) then

      call chem_io_file_name(datafile, filename, tile, pathname=path)
      print *,'WRITE4D: '//trim(filename)//" -> "//trim(datafile)

      call chem_io_file_write(datafile, recvbuf, pos=pos, reccount=nl, &
        recsize=recsize, recstride=recstride, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      if (chem_io_verbose) &
        write(6,'("chem_io_write: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
          trim(datafile), minval(recvbuf), maxval(recvbuf)
    end if

    deallocate(buf4d, recvbuf, stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot deallocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_write_4DR4


! add for netcdf files ajl

  subroutine chem_io_file_read_2d_nc(datafile, buffer, recrange, recsize, recstride, varname,rc)
    use chem_comm_mod, only : chem_comm_get
    use netcdf
    character(len=*),   intent(in)  :: datafile
    real(CHEM_KIND_R4), intent(out) :: buffer(:,:)
    integer, optional,  intent(in)  :: recrange(2)
    integer, optional,  intent(in)  :: recsize
    integer, optional,  intent(in)  :: recstride
    character(len=*), optional,  intent(in)  :: varname
    character*120 namevar
    integer, optional,  intent(out) :: rc
     

    ! -- local variables
    integer :: localrc
    integer :: irec, is, ie,i
    integer :: rcount, rsize, rstride
    integer :: rrange(2)
    integer :: nc,nr
    integer :: mype
    integer :: ncid ! ncid for netcdf file
    integer :: ierr ! netcdf error code
    integer :: varid ! netcdf varid
    if(.not.present(varname))then
      rc=1
      write(6,*)'varname must be present for reading netcdf file input'
      call flush(6)
      return
    endif
!    write(6,*)'read_2d_nc var',trim(varname)
!    call flush(6)
    call chem_comm_get(localpe=mype)
!    write(6,*)'chem_io_file_read_2d_nc mype',mype,shape(buffer)
!    call flush(6)

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS
!    return ! for now

    buffer = 0._CHEM_KIND_R4

    if (present(recrange)) then
      rrange = recrange
    else
      rrange = 1
    end if
    rcount = rrange(2) - rrange(1) + 1

    if (present(recsize)) then
      rsize = recsize
    else
      rsize = size(buffer) / rcount
    end if

    if (present(recstride)) then
      rstride = recstride
    else
      rstride = rsize
    end if
!    write(6,*)'sizebuffer',size(buffer),'rcount',rcount,'rsize',rsize,'rstride',rstride
!    write(70+mype,*)'sizebuffer',size(buffer),'rcount',rcount,'rsize',rsize,'rstride',rstride
!    call flush(70+mype)
!    write(6,*)'prod ',rcount * max(rsize, rstride)
!    call flush(6)

    if (chem_rc_test((size(buffer) < rcount * max(rsize, rstride)), &
        msg="insufficient buffer size", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- open file
!    write(70+mype,*)'open',trim(datafile)
!    call flush(70+mype)
    ierr=nf90_open(trim(datafile),0,ncid)
    if (ierr /= nf90_noerr)then
       write(6,*)'error open read2dnc',trim(datafile)
       call flush(6)
       rc=1
       return
    else
!      if(mype.eq.0)then
!        write(6,*)'open 2d nc ',trim(datafile)
!        call flush(6)
!      endif
    endif
    ierr=nf90_inq_varid(ncid,trim(varname),varid)
    if (ierr /= nf90_noerr)then
      write(6,*)'error getting varid for ',trim(varname),' file ',trim(datafile)
      write(6,*)trim(nf90_strerror(ierr))
      call flush(6)
      do i=1,3
        ierr=nf90_inquire_variable(ncid,i,name=namevar)
        write(6,*)'namevar',i,trim(namevar)
        call flush(6)
      end do
      rc=1
      return
!    else
!      write(6,*)'read id ',trim(varname),varid
!      call flush(6)
    endif
    ierr=nf90_get_var(ncid,varid,buffer)
    if (ierr /= nf90_noerr)then
      write(6,*)'error getting var size buffer',shape(buffer),'var',trim(varname),'file',trim(datafile)
      write(6,*)trim(nf90_strerror(ierr))
      call flush(6)
      rc=1
      return
    else
!      if(mype.eq.0)then
!        write(6,*)'read var',maxval(buffer),minval(buffer)
!        call flush(6)
!      endif
    endif


    ! -- close file
     ierr=nf90_close(ncid)

  end subroutine chem_io_file_read_2d_nc
  subroutine chem_io_file_read_3d_nc(datafile, buffer, recrange, recsize, recstride, varname,rc)
    use chem_comm_mod, only : chem_comm_get
    use netcdf
    character(len=*),   intent(in)  :: datafile
    real(CHEM_KIND_R4), intent(out) :: buffer(:,:,:)
    integer, optional,  intent(in)  :: recrange(2)
    integer, optional,  intent(in)  :: recsize
    integer, optional,  intent(in)  :: recstride
    character(len=*), optional,  intent(in)  :: varname
    character *256 datafileopen
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: irec, is, ie
    integer :: rcount, rsize, rstride
    integer :: rrange(2)
    integer :: nc,nr
    integer :: mype
    integer :: ncid ! ncid for netcdf file
    integer :: ierr ! netcdf error code
    integer :: varid ! netcdf varid
    if(.not.present(varname))then
      rc=1
      write(6,*)'varname must be present for reading netcdf file input'
      call flush(6)
      return
    endif
    call chem_comm_get(localpe=mype)
    ! -- open file
    datafileopen=datafile
!    write(200+mype,*)'open',trim(datafileopen),' varname',trim(varname)
!    call flush(200+mype)
    ierr=nf90_open(trim(datafileopen),0,ncid)
!    write(200+mype,*)'ierr',ierr
!    call flush(200+mype)
    if (ierr /= nf90_noerr)then
       write(6,*)'error open read3dnc',trim(datafile)
       write(200+mype,*)'error open read3dnc',trim(datafile)
       call flush(200+mype)
       call flush(6)
       rc=1
       return
!    else
!      write(200+mype,*)'success open',trim(datafile)
!      call flush(200+mype)
    endif

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    buffer = 0._CHEM_KIND_R4

    if (present(recrange)) then
      rrange = recrange
    else
      rrange = 1
    end if
    rcount = rrange(2) - rrange(1) + 1

    if (present(recsize)) then
      rsize = recsize
    else
      rsize = size(buffer) / rcount
    end if

    if (present(recstride)) then
      rstride = recstride
    else
      rstride = rsize
    end if

    if (chem_rc_test((size(buffer) < rcount * max(rsize, rstride)), &
        msg="insufficient buffer size", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ierr=nf90_inq_varid(ncid,trim(varname),varid)
    if (ierr /= nf90_noerr)then
      write(6,*)'error getting varid for ',trim(varname),' file ',trim(datafile)
      write(200+mype,*)'error getting varid for ',trim(varname),' file ',trim(datafile)
      call flush(200+mype)
      call flush(6)
      write(6,*)trim(nf90_strerror(ierr))
      rc=1
      return
!    else
!      write(200+mype,*)'got varid',trim(varname)
!      call flush(200+mype)
!      write(200+mype,*)'buffer',shape(buffer)
!      call flush(200+mype)
    endif
    ierr=nf90_get_var(ncid,varid,buffer)
    if (ierr /= nf90_noerr)then
      write(6,*)'error getting var size buffer',shape(buffer),'var',trim(varname),'file',trim(datafile)
      write(6,*)trim(nf90_strerror(ierr))
      write(200+mype,*)'error get var',trim(nf90_strerror(ierr))
      call flush(200+mype)
      call flush(6)
      rc=1
      return
!    else
!      write(200+mype,*)'got buffer',trim(varname)
!      call flush(200+mype)

    endif
!     write(200+mype,*)'nf90_close bottom filre_read_3d_nc',trim(varname)
!     call flush(200+mype)
     ierr=nf90_close(ncid)
!     write(200+mype,*)'bottom filre_read_3d_nc',trim(varname)
!     call flush(200+mype)

  end subroutine chem_io_file_read_3d_nc


  subroutine chem_io_writenc_2DR4(filename, farray, path, pos, de, time,varname,units,rc)
    use netcdf
    use chem_comm_mod, only : chem_comm_get
!    use raqmschem_pmgrid_mod, only : iam
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(in)  :: farray(:,:)
    character(len=*), optional, intent(in)  :: path
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(in)  :: de
    integer,          optional, intent(in)  :: time
    character(len=*),optional,  intent(in)  :: varname
    character(len=*),optional,  intent(in)  :: units
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: ids, ide, jds, jde, its, ite, jts, jte
    logical :: localIOflag
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:,:), allocatable, target :: buf2d, recvbuf
    real(chem_kind_R8), pointer, dimension(:,: ) :: lat2dr8,lon2dr8
    real(chem_kind_R4), pointer, dimension(:,: ) :: lat2d,lon2d,latbuf,lonbuf
    integer ncid,ierr,lenf,mype
    lenf=len_trim(filename)

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS
    call chem_comm_get(localpe=mype)

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, lon=lon2dr8,lat=lat2dr8,rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- check size consistency
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buf2d = 0._CHEM_KIND_R4

    buf2d(ids:ide, jds:jde) = farray

    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    allocate(lonbuf(its:ite,jts:jte), stat=localrc)
    allocate(latbuf(its:ite,jts:jte), stat=localrc)
    allocate(lat2d(its:ite,jts:jte))
    allocate(lon2d(its:ite,jts:jte))
    lat2d=0.
    lon2d=0.
    lat2d(ids:ide,jds:jde)=lat2dr8
    lon2d(ids:ide,jds:jde)=lon2dr8
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4
    latbuf=0._CHEM_KIND_R4
    lonbuf=0._CHEM_KIND_R4

    call chem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    call chem_comm_reduce(lat2d, latbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    call chem_comm_reduce(lon2d, lonbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (localIOflag) then

      call chem_io_file_name(datafile, filename, tile, pathname=path)
        if(present(varname))then
          if(present(units))then
            call chem_io_file_writenc2d(datafile, filename,recvbuf,latbuf,lonbuf, &
            time=time, varname=varname,units=units,rc=localrc)
          else
            call chem_io_file_writenc2d(datafile, filename,recvbuf,latbuf,lonbuf, &
            time=time, varname=varname,rc=localrc)
          endif
        else
          if(present(units))then
            call chem_io_file_writenc2d(datafile, filename,recvbuf,latbuf,lonbuf, &
            varname=varname,units=units,rc=localrc)
          else
            call chem_io_file_writenc2d(datafile, filename,recvbuf,latbuf,lonbuf, &
            varname=varname,rc=localrc)
          endif
        endif
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    end if

    deallocate(buf2d, recvbuf, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_writenc_2DR4
  subroutine chem_io_file_writenc2d(datafile, filename,buffer, lat,lon,time,varname, units,rc)
    use netcdf
    character(len=*),           intent(in)  :: datafile
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(in)  :: buffer(:,:)
    real(CHEM_KIND_R4),         intent(in)  :: lat(:,:),lon(:,:)
    integer,optional,intent(in)             :: time
    integer,          optional, intent(out) :: rc
    character (len=*), optional, intent(in) :: varname
    character (len=*), optional, intent(in) :: units
    integer ncid,ierr,varid,dims(2),ilat,ilon,idimzt,lenf,i
    real(chem_kind_r4),allocatable :: out(:,:,:)

    ! -- local variables
    integer :: localrc
    real(CHEM_KIND_R4) :: atest
    integer idimxt,idimyt,idims(3),itime,idims2(2),count(3),start(3)
    logical dodefine,exist
    integer idgridxt,idgridyt,idtime
    real(CHEM_KIND_R8),allocatable :: gridxy(:)
    real(CHEM_KIND_R8) :: rtime(1)
    

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    inquire(file=trim(datafile),exist=exist)
    if(.not.present(time).and. .not.exist)then
      ierr=nf90_create(trim(datafile),nf90_clobber+nf90_netcdf4,ncid)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error create no time',trim(datafile)
           write(6,*)trim(nf90_strerror(ierr))
      endif
      dodefine=.true.
    else
      if(exist)then
        dodefine=.false.
        ierr=nf90_open(trim(datafile),nf90_write,ncid)
        ierr=nf90_inq_varid(ncid,varname,varid)
        if (ierr /= nf90_noerr)then
          ierr=nf90_inq_dimid(ncid,'grid_xt',idimxt)
          if (ierr /= nf90_noerr)then
              write(6,*)'error grid_xt inq'
              call flush(6)
          endif
          ierr=nf90_inq_dimid(ncid,'grid_yt',idimyt)
          if (ierr /= nf90_noerr)then
              write(6,*)'error grid_yt inq'
              call flush(6)
          endif
          ierr=nf90_inq_dimid(ncid,'time',itime)
          if (ierr /= nf90_noerr)then
              write(6,*)'error time inq'
              call flush(6)
          endif
          idims(1)=idimxt
          idims(2)=idimyt
          idims(3)=itime
          ierr=nf90_redef(ncid)
          if (ierr /= nf90_noerr)then
             write(6,*)'error redef'
              write(6,*)trim(nf90_strerror(ierr))
          endif
          ierr=nf90_def_var(ncid,trim(varname),NF90_FLOAT,idims,varid,chunksizes=[dims(1),28,1], &
              shuffle=.true.,deflate_level=1)
          if (ierr /= nf90_noerr)then
              write(6,*)'erroro define ',varname
              write(6,*)trim(nf90_strerror(ierr))
          endif
          if(present(units))then
            ierr=nf90_put_att(ncid,varid,'units',trim(units))
            if (ierr /= nf90_noerr)then
              write(6,*)'error put_att 1,trim(varname),trim(units)'
              write(6,*)trim(nf90_strerror(ierr))
            endif
          endif
          ierr=nf90_enddef(ncid)
          if (ierr /= nf90_noerr)then
                write(6,*)'error enddef ',varname
                write(6,*)trim(nf90_strerror(ierr))
          endif
        endif
      elseif(time<=1)then
        inquire(file=trim(datafile),exist=exist)
        if(exist)then
          dodefine=.false.
          ierr=nf90_open(trim(datafile),nf90_write,ncid)
          ierr=nf90_inq_varid(ncid,varname,varid)
          if (ierr /= nf90_noerr)then
            ierr=nf90_inq_dimid(ncid,'grid_xt',idimxt)
            if (ierr /= nf90_noerr)then
                write(6,*)'error grid_xt inq'
                call flush(6)
            endif
            ierr=nf90_inq_dimid(ncid,'grid_yt',idimyt)
            if (ierr /= nf90_noerr)then
                write(6,*)'error grid_yt inq'
                call flush(6)
            endif
            ierr=nf90_inq_dimid(ncid,'time',itime)
            if (ierr /= nf90_noerr)then
                write(6,*)'error time inq'
                call flush(6)
            endif
            idims(1)=idimxt
            idims(2)=idimyt
            idims(3)=itime
            ierr=nf90_redef(ncid)
            if (ierr /= nf90_noerr)then
               write(6,*)'error redef'
                write(6,*)trim(nf90_strerror(ierr))
            endif
            ierr=nf90_def_var(ncid,varname,nf90_float,idims,varid)
            if (ierr /= nf90_noerr)then
                write(6,*)'erroro define ',varname
                write(6,*)trim(nf90_strerror(ierr))
            endif
            if(present(units))then
               ierr=nf90_put_att(ncid,varid,'units',trim(units))
              if (ierr /= nf90_noerr)then
                write(6,*)'error put_att 2',trim(varname),trim(units)
                write(6,*)trim(nf90_strerror(ierr))
              endif
            endif
            ierr=nf90_enddef(ncid)
            if (ierr /= nf90_noerr)then
                write(6,*)'error enddef ',varname
                write(6,*)trim(nf90_strerror(ierr))
            endif
          endif
        else
          write(6,*)'nf90_create',trim(datafile)
          ierr=nf90_create(trim(datafile),NF90_clobber+nf90_netcdf4,ncid)
          dodefine=.true.
        endif
        if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error create time',trim(datafile)
           write(6,*)trim(nf90_strerror(ierr))
        endif
      else
        ierr=nf90_open(trim(datafile),NF90_WRITE,ncid)
        if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error open time',trim(datafile)
           write(6,*)trim(nf90_strerror(ierr))
        endif
        dodefine=.false.
      endif
    endif
    if (ierr /= nf90_noerr)then
      write(6,*)'zzzz ajl error open ',trim(datafile)
      call flush(6)
    endif
    dims=shape(buffer)
    allocate (out(dims(1),dims(2),1))
    out(:,:,1)=buffer(:,:)
    if(dodefine)then
      write(6,*)'dims',dims(1:2)
      ierr=nf90_def_dim(ncid,'grid_xt',dims(1),idimxt)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_dim grid_xt'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_def_dim(ncid,'grid_yt',dims(2),idimyt)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_dim grid_yt'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_def_dim(ncid,'time',NF90_UNLIMITED,itime)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_dim time'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_def_var(ncid,'time',nf90_double,itime,idtime)
      ierr=nf90_def_var(ncid,'grid_xt',nf90_double,idimxt,idgridxt)
      ierr=nf90_def_var(ncid,'grid_yt',nf90_double,idimyt,idgridyt)
      ierr=nf90_put_att(ncid,idtime,'long_name','time')
!      ierr=nf90_put_att(ncid,idtime,'units','days since '//trim(cdatestr))
      ierr=nf90_put_att(ncid,idtime,'cartesian_axis','T')
      ierr=nf90_put_att(ncid,idtime,'calendar_type','JULIAN')
      ierr=nf90_put_att(ncid,idtime,'calendar','JULIAN')
      ierr=nf90_put_att(ncid,idgridxt,'long_name','T-cell longitude')
      ierr=nf90_put_att(ncid,idgridxt,'units','degrees_E')
      ierr=nf90_put_att(ncid,idgridxt,'cartesian_axis','X')
      ierr=nf90_put_att(ncid,idgridyt,'long_name','T-cell latitude')
      ierr=nf90_put_att(ncid,idgridyt,'units','degrees_N')
      ierr=nf90_put_att(ncid,idgridyt,'cartesian_axis','Y')
      idims(1)=idimxt
      idims(2)=idimyt
      idims(3)=itime
      idims2(1)=idimxt
      idims2(2)=idimyt

      ierr=nf90_def_var(ncid,'lat',nf90_float,idims2,ilat) 
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_var write3dnc lat'
           write(6,*)trim(nf90_strerror(ierr))
           call flush(6)
      endif
      ierr=nf90_def_var(ncid,'lon',nf90_float,idims2,ilon) 
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_var lon write3dnc'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      if(present(varname))then
        ierr=nf90_def_var(ncid,trim(varname),NF90_FLOAT,idims,varid,chunksizes=[dims(1),28,1],shuffle=.true.,deflate_level=1)
      else
        ierr=nf90_def_var(ncid,trim(filename),NF90_FLOAT,idims,varid)
        write(6,*)'no varname'
      endif
      if (ierr /= nf90_noerr)then
           if(present(varname))then
           write(6,*)'zzzz ajl error def_var varname',varname
           endif
           write(6,*)trim(nf90_strerror(ierr))
      endif
      if(present(units))then
        ierr=nf90_put_att(ncid,varid,'units',trim(units))
        if (ierr /= nf90_noerr)then
          write(6,*)'error put_att 3',trim(units)
          write(6,*)trim(nf90_strerror(ierr))
        endif
      endif
      ierr=nf90_enddef(ncid)
      if (ierr /= nf90_noerr)then
       write(6,*)'zzzz ajl enddef'
       write(6,*)trim(nf90_strerror(ierr))
      endif
!      rtime=fracday
!      ierr=nf90_put_var(ncid,idtime,rtime)
      allocate (gridxy(dims(1)))
      do i=1,dims(1)
        gridxy(i)=i
      end do
      ierr=nf90_put_var(ncid,idgridxt,gridxy)
      ierr=nf90_put_var(ncid,idgridyt,gridxy)
      deallocate(gridxy)
    endif
    if(present(time))then
      start=1
      start(3)=time
      count(1)=dims(1)
      count(2)=dims(2)
      count(3)=1
      if(present(varname))then
        ierr=nf90_inq_varid(ncid,trim(varname),varid)
      else
        ierr=nf90_inq_varid(ncid,trim(filename),varid)
      endif
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error inq_varid bb'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_put_var(ncid,varid,out,start=start,count=count)
      if (ierr /= nf90_noerr)then
        write(6,*)'put var outtime ',time
        write(6,*)trim(nf90_strerror(ierr))
      endif
    else
      ierr=nf90_put_var(ncid,varid,out)
      if (ierr /= nf90_noerr)then
         write(6,*)'zzzz ajl put var out'
         write(6,*)trim(nf90_strerror(ierr))
      endif
    endif
!    if(dodefine)then
!      ierr=nf90_put_var(ncid,ilon,lon)
!      if (ierr /= nf90_noerr)then
!        write(6,*)'put_var lon'
!        write(6,*)trim(nf90_strerror(ierr))
      !endif
!      ierr=nf90_put_var(ncid,ilat,lat)
!      if (ierr /= nf90_noerr)then
!        write(6,*)'put_var lat'
        !write(6,*)trim(nf90_strerror(ierr))
!      endif
!    endif
    deallocate (out)
    ierr=nf90_close(ncid)
    if (ierr /= nf90_noerr)then
      write(6,*)'zzzz ajl error close'
      write(6,*)trim(nf90_strerror(ierr))
      call flush(6)
    endif

  end subroutine chem_io_file_writenc2d

  logical function inquire_file_var(datafile,varname,path)
    use netcdf
!    use raqmschem_pmgrid_mod, only : iam
    integer ncid,idvar,ierr
    character(len=*),   intent(in)  :: datafile,path
    character(len=*),  intent(in)  :: varname
    logical filevar,localIOflag
    integer tilecomm,localrc
    call chem_model_get(tileComm=tileComm,localIOflag=localIOflag)
!    write(200+iam,*)'inquire_file_var'
!    call flush(200+iam)
!    write(200+iam,*)'path',trim(path)
!    call flush(200+iam)
!    write(200+iam,*)'datafile',trim(datafile)
!    call flush(200+iam)
  
    if(localIOflag)then      
      ierr=nf90_open(trim(path)//trim(datafile),0,ncid)
      if (ierr /= nf90_noerr)then
        write(6,*)'path',trim(path)
        write(6,*)'datafile',trim(datafile)
        write(6,*)'error open inquire_file_var ',trim(path)//trim(datafile)
        call flush(6)
        filevar=.false.
      else
        ierr=nf90_inq_varid(ncid,varname,idvar)
        if(ierr==nf90_noerr)then
          filevar=.true.
        else
          filevar=.false.
        endif
        ierr=nf90_close(ncid)
      endif
    endif
    call chem_comm_bcast(filevar, comm=tileComm,rc=localrc)
    inquire_file_var=filevar
  end function inquire_file_var
  logical function inquire_file_dim(datafile,dimname,path)
    use netcdf
!    use raqmschem_pmgrid_mod, only : iam
    integer ncid,idvar,ierr,tilecomm,localrc
    character(len=*),   intent(in)  :: datafile,path
    character(len=*),  intent(in)  :: dimname
    logical :: localIOflag,inquire_file_dimin
    call chem_model_get(tileComm=tileComm,localIOflag=localIOflag)
!    write(200+iam,*)'inquire_file_var'
!    call flush(200+iam)
!    write(200+iam,*)'path',trim(path)
!    call flush(200+iam)
!    write(200+iam,*)'datafile',trim(datafile)
!    call flush(200+iam)
  
    if(localIOflag)then      
      ierr=nf90_open(trim(path)//trim(datafile),0,ncid)
      if (ierr /= nf90_noerr)then
        write(6,*)'path',trim(path)
        write(6,*)'datafile',trim(datafile)
        write(6,*)'error open file_dim ',trim(path)//trim(datafile)
        call flush(6)
        inquire_file_dimin=.false.
!       return
!       write(6,*)'did open',iam
      else
        ierr=nf90_inq_dimid(ncid,dimname,idvar)
        if(ierr==nf90_noerr)then
          inquire_file_dimin=.true.
        else
          inquire_file_dimin=.false.
        endif
      endif
      ierr=nf90_close(ncid)
    endif
    call chem_comm_bcast(inquire_file_dimin, comm=tileComm,rc=localrc)
    inquire_file_dim=inquire_file_dimin
  end function inquire_file_dim

end module chem_io_mod
