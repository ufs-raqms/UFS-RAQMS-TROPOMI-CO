module chem_io_mod

  use mpi
  use chem_rc_mod
  use chem_types_mod
  use chem_comm_mod
  use chem_model_mod,only : chem_model_get => raqmschem_model_get
  use chem_model_mod,only : chem_model_domain_get 
  use chem_model_mod,only : chem_model_set => raqmschem_model_set

  implicit none

  integer, parameter :: ioUnit = 100

  interface chem_io_read
    module procedure chem_io_read_2DR4
    module procedure chem_io_read_3DR4
  end interface chem_io_read

  interface chem_io_write
    module procedure chem_io_write_2DR4
    module procedure chem_io_write_3DR4
    module procedure chem_io_write_3DR8
  end interface chem_io_write

  private

  public :: raqmschem_io_init
  public :: chem_io_read
  public :: chem_io_write

contains

  subroutine raqmschem_io_init(rc)
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
      write(6,'("chem_io_init: PET:",i2," DE:",i02," tile=",i0," - comm=",i0," PE:",i0,"/",i0)') &
        localpe, de, tile, tileComm, pe, npe
      flush(6)
    end do

  end subroutine raqmschem_io_init


  subroutine raqmschem_io_file_name(fullname, filename, tile, pathname)
    character(len=*),           intent(out) :: fullname
    character(len=*),           intent(in)  :: filename
    integer,                    intent(in)  :: tile
    character(len=*), optional, intent(in)  :: pathname

    ! -- local variables
    integer :: lstr
    character(len=CHEM_MAXSTR) :: fname

    ! -- begin
    fname = ""
    fullname = ""

    lstr = len_trim(filename)
    if (lstr > 4) then
!      print *, 'filename = ', filename(lstr-3:lstr)
      if (filename(lstr-3:lstr) == ".dat") then
!        we have a different path for now
         write(fname,'(a,"_tile",i0,".dat")')filename(1:lstr-4),tile
!        write(fname, '("tile",i0,"/",a)') tile, trim(filename)
      elseif(filename(lstr-2:lstr) == '.nc')then
         write(fname,'(a,"_tile",i0,".nc")')filename(1:lstr-3),tile

      else
        write(fname, '(a,".tile",i0,".dat")') trim(filename), tile
      end if
    else
      write(fname, '(a,".tile",i0,".dat")') trim(filename), tile
    end if
!    write(6,*)'fname full',trim(fname)

    if (present(pathname)) then
      lstr = len_trim(pathname)
      if (pathname(lstr:lstr) == "/") then
        fullname = trim(pathname) // trim(fname)
      else
        fullname = trim(pathname) // "/" // trim(fname)
      end if
    else
      fullname = trim(fname)
    end if
!    write(6,*)'bottom raqmschem_io_file_name fullname',trim(fullname)
!    call flush(6)

  end subroutine raqmschem_io_file_name
  subroutine bigchem_io_file_read(datafile, buffer, recrange, recsize, recstride, rc)
    use chem_comm_mod, only : chem_comm_get
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
    integer :: nc,nr
    integer :: mype
    call chem_comm_get(localpe=mype)
!    write(6,*)'chem_io_filie_read mype',mype,shape(buffer)
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
    open(unit=ioUnit, file=trim(datafile), form='unformatted', action='read', position='rewind', convert='big_endian',iostat=localrc)
!    write(6,*)'open file ',trim(datafile),'localrc',localrc
!    call flush(6)
!    write(6,*)'rrange',rrange,'rsize',rsize
!    call flush(6)
!    read(iounit)nc,nr
!    write(6,*)'nc',nc,nr
!    call flush(6)
!    return
    if (chem_rc_test((localrc /= 0), msg="Failure opening file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- advance to first record
!    write(6,*)'rrange',rrange
!    call flush(6)
!    write(72+mype,*)'rrange',rrange
!    call flush(72+mype)
!   we need to skip 2 records so rrange(1) should be 3
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
!    write(6,*)'rcount',rcount,'rstride',rstride
!    call flush(6)
    do irec = 1, rcount
!      write(6,*)'read buffer',shape(buffer),'rsize',rsize
!      call flush(6)
!      write(72+mype,*)'read buffer',shape(buffer),ie
!      call flush(72+mype)
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

  end subroutine bigchem_io_file_read
  

  subroutine chem_io_file_read(datafile, buffer, recrange, recsize, recstride, rc)
    use chem_comm_mod, only : chem_comm_get
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
    integer :: nc,nr
    integer :: mype
    call chem_comm_get(localpe=mype)
!    write(6,*)'chem_io_filie_read mype',mype,shape(buffer)
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
    !call flush(70+mype)
!    write(6,*)'prod ',rcount * max(rsize, rstride)
!    call flush(6)

    if (chem_rc_test((size(buffer) < rcount * max(rsize, rstride)), &
        msg="insufficient buffer size", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- open file
!    write(70+mype,*)'open',trim(datafile)
!    call flush(70+mype)
      open(unit=ioUnit, file=trim(datafile), form='unformatted', action='read', position='rewind', iostat=localrc)
!    write(6,*)'open file ',trim(datafile),'localrc',localrc
!    call flush(6)
!    write(6,*)'rrange',rrange,'rsize',rsize
!    call flush(6)
!    read(iounit)nc,nr
!    write(6,*)'nc',nc,nr
!    call flush(6)
!    return
    if (chem_rc_test((localrc /= 0), msg="Failure opening file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- advance to first record
!    write(6,*)'rrange',rrange
!    call flush(6)
!    write(72+mype,*)'rrange',rrange
!    call flush(72+mype)
!   we need to skip 2 records so rrange(1) should be 3
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
!    write(6,*)'rcount',rcount,'rstride',rstride
!    call flush(6)
    do irec = 1, rcount
!      write(6,*)'read buffer',shape(buffer),'rsize',rsize
!      call flush(6)
!      write(72+mype,*)'read buffer',shape(buffer),ie
!      call flush(72+mype)
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


  subroutine chem_io_file_write(datafile, buffer, pos, rc)
    character(len=*),           intent(in)  :: datafile
    real(CHEM_KIND_R4),         intent(in)  :: buffer(:)
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
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

    open(unit=ioUnit, file=trim(datafile), form='unformatted', action='write', &
      position=trim(filepos), iostat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Failure opening file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'chem_io_file_write buffer',shape(buffer),trim(datafile)
!    call flush(6)
    write(unit=ioUnit, iostat=localrc) buffer
    if (chem_rc_test((localrc /= 0), msg="Failure reading data from file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    close(unit=ioUnit)

  end subroutine chem_io_file_write

  subroutine chem_io_file_writenc(datafile, filename,buffer, lat,lon,pos, rc)
    use netcdf
    character(len=*),           intent(in)  :: datafile
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(in)  :: buffer(:,:)
    real(CHEM_KIND_R4),         intent(in)  :: lat(:,:),lon(:,:)
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(out) :: rc
    integer ncid,ierr,itest,dims(2),ilat,ilon
    real(chem_kind_r4),allocatable :: out(:,:,:)

    ! -- local variables
    integer :: localrc
    character(len=6) :: filepos
    real(CHEM_KIND_R4) :: atest
    integer idimxt,idimyt,idims(3),itime,idims2(2)

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    filepos = 'rewind'
!    write(6,*)'datafile writenc',trim(datafile),'pos',pos
    if (present(pos)) then
      select case (trim(pos))
        case ('a', 'append')
          filepos = 'append'
        case default
          filepos = 'rewind'
      end select
    end if
!    write(6,*)'ajl nf90_create ',trim(datafile)
!    call flush(6)
    ierr=nf90_create(trim(datafile),0,ncid)
    if (ierr /= nf90_noerr)then
      write(6,*)'error open ',trim(datafile)
      call flush(6)
    endif
    dims=shape(buffer)
    allocate (out(dims(1),dims(2),1))
    out(:,:,1)=buffer(:,:)
!    write(6,*)'dims',dims
!    call flush(6)
    ierr=nf90_def_dim(ncid,'grid_xt',dims(1),idimxt)
    ierr=nf90_def_dim(ncid,'grid_yt',dims(2),idimyt)
    ierr=nf90_def_dim(ncid,'time',NF90_UNLIMITED,itime)
    idims(1)=idimxt
    idims2(1)=idimxt
    idims(2)=idimyt
    idims2(2)=idimyt
    idims(3)=itime

    ierr=nf90_def_var(ncid,'lat',nf90_float,idims2,ilat) 
    ierr=nf90_def_var(ncid,'lon',nf90_float,idims2,ilon) 
    ierr=nf90_def_var(ncid,trim(filename),NF90_FLOAT,idims,itest)
    if (ierr /= nf90_noerr)then
         write(6,*)'error def_var'
    endif
    ierr=nf90_enddef(ncid)
    if (ierr /= nf90_noerr)then
     write(6,*)'enddef'
    endif
    ierr=nf90_put_var(ncid,itest,out)
    ierr=nf90_put_var(ncid,ilon,lon)
    ierr=nf90_put_var(ncid,ilat,lat)
    deallocate (out)
    if (ierr /= nf90_noerr)then
      write(6,*)'put_var'
    endif
    ierr=nf90_close(ncid)
    if (ierr /= nf90_noerr)then
      write(6,*)'close'
    endif
  

#if 0
    open(unit=ioUnit, file=trim(datafile), form='unformatted', action='write', &
      position=trim(filepos), iostat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Failure opening file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'chem_io_file_write buffer',shape(buffer),trim(datafile)
!    call flush(6)
    write(unit=ioUnit, iostat=localrc) buffer
    if (chem_rc_test((localrc /= 0), msg="Failure reading data from file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    close(unit=ioUnit)
#endif

  end subroutine chem_io_file_writenc
  subroutine chem_io_file_writenc3d(datafile, filename,buffer, lat,lon,pos,time,varname, rc)
    use netcdf
    character(len=*),           intent(in)  :: datafile
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(in)  :: buffer(:,:,:)
    real(CHEM_KIND_R4),         intent(in)  :: lat(:,:),lon(:,:)
    character(len=*), optional, intent(in)  :: pos
    integer,optional,intent(in)             :: time
    integer,          optional, intent(out) :: rc
    character (len=*), optional, intent(in) :: varname
    integer ncid,ierr,itest,dims(3),ilat,ilon,idimzt,lenf
    real(chem_kind_r4),allocatable :: out(:,:,:,:)

    ! -- local variables
    integer :: localrc
    character(len=6) :: filepos
    real(CHEM_KIND_R4) :: atest
    integer idimxt,idimyt,idims(4),itime,idims2(2),count(4),start(4)
    logical dodefine

    ! -- begin
!    write(6,*)'writenc3d',maxval(buffer),minval(buffer)
    if (present(rc)) rc = CHEM_RC_SUCCESS

    filepos = 'rewind'
!    write(6,*)'datafile writenc',trim(datafile),'pos',pos
    if (present(pos)) then
      select case (trim(pos))
        case ('a', 'append')
          filepos = 'append'
        case default
          filepos = 'rewind'
      end select
    end if
!    write(6,*)'ajl nf90_create ',trim(datafile)
!    call flush(6)
    if(.not.present(time))then
      ierr=nf90_create(trim(datafile),0,ncid)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error create no time',trim(datafile)
           write(6,*)trim(nf90_strerror(ierr))
      endif
      dodefine=.true.
    else
!      write(6,*)'ajl time present',time
!      call flush(6)
      if(time<=1)then
        ierr=nf90_create(trim(datafile),0,ncid)
        if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error create time',trim(datafile)
           write(6,*)trim(nf90_strerror(ierr))
        endif
        dodefine=.true.
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
    allocate (out(dims(1),dims(2),dims(3),1))
    out(:,:,:,1)=buffer(:,:,:)
!    write(6,*)'dims',dims
!    call flush(6)
    if(dodefine)then
      ierr=nf90_def_dim(ncid,'grid_xt',dims(1),idimxt)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_dim'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_def_dim(ncid,'grid_yt',dims(2),idimyt)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_dim'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_def_dim(ncid,'grid_zt',dims(3),idimzt)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_dim'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_def_dim(ncid,'time',NF90_UNLIMITED,itime)
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_dim'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      idims(1)=idimxt
      idims2(1)=idimxt
      idims(2)=idimyt
      idims(3)=idimzt
      idims2(2)=idimyt
      idims(4)=itime

      ierr=nf90_def_var(ncid,'lat',nf90_float,idims2,ilat) 
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_var'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_def_var(ncid,'lon',nf90_float,idims2,ilon) 
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_var'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      if(present(varname))then
        ierr=nf90_def_var(ncid,trim(varname),NF90_FLOAT,idims,itest)
      else
        ierr=nf90_def_var(ncid,trim(filename),NF90_FLOAT,idims,itest)
      endif
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error def_var'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_enddef(ncid)
      if (ierr /= nf90_noerr)then
       write(6,*)'zzzz ajl enddef'
       write(6,*)trim(nf90_strerror(ierr))
      endif
    endif
    if(present(time))then
      start=1
      start(4)=time
      count(1)=dims(1)
      count(2)=dims(2)
      count(3)=dims(3)
      count(4)=1
      if(present(varname))then
        ierr=nf90_inq_varid(ncid,trim(varname),itest)
      else
        ierr=nf90_inq_varid(ncid,trim(filename),itest)
      endif
      if (ierr /= nf90_noerr)then
           write(6,*)'zzzz ajl error inq_varid'
           write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_put_var(ncid,itest,out,start=start,count=count)
      if (ierr /= nf90_noerr)then
        write(6,*)'put var outtime ',time
        write(6,*)trim(nf90_strerror(ierr))
      endif
    else
      ierr=nf90_put_var(ncid,itest,out)
      if (ierr /= nf90_noerr)then
         write(6,*)'zzzz ajl put var out'
         write(6,*)trim(nf90_strerror(ierr))
      endif
    endif
    if(dodefine)then
      ierr=nf90_put_var(ncid,ilon,lon)
      if (ierr /= nf90_noerr)then
        write(6,*)'put_var'
        write(6,*)trim(nf90_strerror(ierr))
      endif
      ierr=nf90_put_var(ncid,ilat,lat)
      if (ierr /= nf90_noerr)then
        write(6,*)'put_var'
        write(6,*)trim(nf90_strerror(ierr))
      endif
    endif
    deallocate (out)
    ierr=nf90_close(ncid)
    if (ierr /= nf90_noerr)then
      write(6,*)'zzzz ajl error close'
      write(6,*)trim(nf90_strerror(ierr))
    endif
  

#if 0
    open(unit=ioUnit, file=trim(datafile), form='unformatted', action='write', &
      position=trim(filepos), iostat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Failure opening file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'chem_io_file_write buffer',shape(buffer),trim(datafile)
!    call flush(6)
    write(unit=ioUnit, iostat=localrc) buffer
    if (chem_rc_test((localrc /= 0), msg="Failure reading data from file: "//trim(datafile), &
        file=__FILE__, line=__LINE__, rc=rc)) return
    close(unit=ioUnit)
#endif

  end subroutine chem_io_file_writenc3d


  subroutine chem_io_read_2DR4(filename, farray, path, recrange, recsize, recstride, de, bigendian,rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(out) :: farray(:,:)
    character(len=*), optional, intent(in)  :: path
    integer,          optional, intent(in)  :: recrange(2)
    integer,          optional, intent(in)  :: recsize
    integer,          optional, intent(in)  :: recstride
    integer,          optional, intent(in)  :: de
    logical,          optional, intent(in)  :: bigendian
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
!    write(6,*)'bsize',bsize
!    call flush(6)
    allocate(buffer(bsize(1)*bsize(2)), stat=localrc)
!    write(6,*)'localrc',localrc
!    call flush(6)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buffer = 0._CHEM_KIND_R4
!    write(6,*)'localIOflag',localIOflag
!    call flush(6)

    if (localIOflag) then
!      write(6,*)'call raqmschem_io_file_name',trim(filename),'path',trim(path)
!      call flush(6)
      call raqmschem_io_file_name(datafile, filename, tile, pathname=path)
!     write(6,*)'datafile',trim(datafile)
!     call flush(6)
!     write(6,*)'present recrange',present(recrange)
!     call flush(6)
!     write(6,*)'present recsize',present(recsize)
!     call flush(6)
!     write(6,*)'present recstride',present(recstride)
!     call flush(6)
!!     write(6,*)'recrange',recrange,'recsize',recsize,'recstride',recstride
!!     call flush(6)
!     write(6,*)'call chem_io_file_read',trim(datafile)
!     call flush(6)
      if(present(bigendian))then
        call bigchem_io_file_read(datafile, buffer, recrange=recrange, recsize=recsize, recstride=recstride, rc=localrc)
      else
        call chem_io_file_read(datafile, buffer, recrange=recrange, recsize=recsize, recstride=recstride, rc=localrc)
      endif
!      write(6,*)'chem_io_file_read ',localrc
!      call flush(6)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      write(6,'("chem_data_read: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
        trim(datafile), minval(buffer), maxval(buffer)
    end if
!    write(6,*)'call chem_comm_bcast',shape(buffer),'localIOflag',localIOflag,'bsize',bsize
!    call flush(6)
!    call mpi_barrier(tilecomm,localrc)

    call chem_comm_bcast(buffer, comm=tileComm, rc=localrc)
!     write(6,*)'after bcast',localrc
!     call flush(6)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf2d(its:ite,jts:jte), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    buf2d = reshape(buffer, bsize)

    farray = buf2d(ids:ide, jds:jde) 
!    write(6,*)'farray',shape(farray)
!    call flush(6)
     
    deallocate(buffer, buf2d, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_read_2DR4


  subroutine chem_io_read_3DR4(filename, farray, path, recrange, recsize, recstride, de, bigendian,rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(out) :: farray(:,:,:)
    character(len=*), optional, intent(in)  :: path
    integer,          optional, intent(in)  :: recrange(2)
    integer,          optional, intent(in)  :: recsize
    integer,          optional, intent(in)  :: recstride
    integer,          optional, intent(in)  :: de
    logical,          optional, intent(in)  :: bigendian
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

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    bsize = (/ ite-its+1, jte-jts+1, size(farray, dim=3) /)

    ! -- check size consistency
    if (chem_rc_test((size(farray) /= (ide-ids+1)*(jde-jds+1)*bsize(3)), &
      msg="size of input array inconsistent with domain decomposition", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buffer(bsize(1)*bsize(2)*bsize(3)), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return
    buffer = 0._CHEM_KIND_R4

    if (localIOflag) then

      call raqmschem_io_file_name(datafile, filename, tile, pathname=path)
      if(present(bigendian))then
        call bigchem_io_file_read(datafile, buffer, recrange=recrange, recsize=recsize, recstride=recstride, rc=localrc)
      else
        call chem_io_file_read(datafile, buffer, recrange=recrange, recsize=recsize, recstride=recstride, rc=localrc)
      endif
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      write(6,'("chem_data_read: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
        trim(datafile), minval(buffer), maxval(buffer)
    end if

    call chem_comm_bcast(buffer, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    allocate(buf3d(its:ite,jts:jte,bsize(3)), stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    buf3d = reshape(buffer, bsize)

    farray = buf3d(ids:ide, jds:jde, :)
     
    deallocate(buffer, buf3d, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_read_3DR4


  subroutine chem_io_write_2DR4(filename, farray, path, pos, de, rc)
    use netcdf
    use chem_comm_mod, only : chem_comm_get
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
!    write(6,*)'get lat2d,lon2d de',de
!    call flush(6)
    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, lon=lon2dr8,lat=lat2dr8,rc=localrc)
!    write(6,*)'subde',de,'pos',pos,'ids',ids,ide,'jds',jds,jde,'its',its,ite,'jts',jts,jte
!    call flush(6)
!     write(6,*)'got latlon',localrc
!     call flush(6)
!    write(6,*)'lbound ',lbound(lat2dr8)
!    write(6,*)'ubound ',ubound(lat2dr8)
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
!    write(6,*)'buf2d',ids,ide,'jds',jds,jde,'shapefarray',shape(farray)
!    write(6,*)'its',its,ite,'jts',jts,jte
!    call flush(6)

    allocate(recvbuf(its:ite,jts:jte), stat=localrc)
    allocate(lonbuf(its:ite,jts:jte), stat=localrc)
    allocate(latbuf(its:ite,jts:jte), stat=localrc)
    allocate(lat2d(its:ite,jts:jte))
    allocate(lon2d(its:ite,jts:jte))
    lat2d=0.
    lon2d=0.
    lat2d(ids:ide,jds:jde)=lat2dr8
    lon2d(ids:ide,jds:jde)=lon2dr8
!    write(6,*)'lat2d',maxval(lat2d),minval(lat2d),'shape',shape(lat2d)
!    write(6,*)'lon2d',maxval(lon2d),minval(lon2d),'shape',shape(lon2d)
!    call flush(6)
!    write(6,*)'shape recvbuf',shape(recvbuf)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot allocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4
    latbuf=0._CHEM_KIND_R4
    lonbuf=0._CHEM_KIND_R4
!     write(6,*)'shape latbuf',shape(latbuf),'recvbuf',shape(recvbuf)
!    write(70+mype,*)'de',de,ids,ide,jds,jde,'local',localIoflag
!    call flush(70+mype)

    call chem_comm_reduce(buf2d, recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    call chem_comm_reduce(lat2d, latbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    call chem_comm_reduce(lon2d, lonbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (localIOflag) then
!      write(6,*)'write_io filename',trim(filename),' path ',trim(path)
!      call flush(6)
!      write(6,*)'lonbuf',maxval(lonbuf),minval(lonbuf),'lat',maxval(latbuf),minval(latbuf)
!      write(6,*)'shape',shape(lonbuf)

      call raqmschem_io_file_name(datafile, filename, tile, pathname=path)
!      write(6,*)'datafile write ',trim(datafile)
!      call flush(6)
!      write(6,*)'shape recvbuf',shape(recvbuf),'size bud2d',size(buf2d)
      if(filename(lenf-2:lenf).eq.'.nc')then
!        write(6,*)'do netcdf',trim(datafile)
!        call flush(6)
!        call chem_io_file_writenc(datafile, reshape(recvbuf, (/size(buf2d)/)), &
        call chem_io_file_writenc(datafile, filename,recvbuf,latbuf,lonbuf, &
        pos=pos, rc=localrc)
      else
!        write(6,*)'do binary'
!        call flush(6)
        call chem_io_file_write(datafile, reshape(recvbuf, (/size(buf2d)/)), &
        pos=pos, rc=localrc)
      endif
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      write(6,'("chem_data_write: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
        trim(datafile), minval(recvbuf), maxval(recvbuf)
    end if

    deallocate(buf2d, recvbuf, stat=localrc)
    if (chem_rc_test((localrc /= 0), &
      msg="Cannot deallocate read buffer", &
      file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_write_2DR4


  subroutine chem_io_write_3DR4(filename, farray, order, path, pos, de,time, varname,rc)
    character(len=*),           intent(in)  :: filename
    real(CHEM_KIND_R4),         intent(in)  :: farray(:,:,:)
    character(len=*), optional, intent(in)  :: order
    character(len=*), optional, intent(in)  :: path
    character(len=*), optional, intent(in)  :: pos
    integer,          optional, intent(in)  :: de
    integer,          optional, intent(in)  :: time
    integer,          optional, intent(out) :: rc
    character(len=*),optional,  intent(in)  :: varname

    ! -- local variables
    integer :: localrc
    integer :: tile, tileComm
    integer :: i, j, k, id, jd, lbuf
    integer :: ids, ide, jds, jde, its, ite, jts, jte, nk,lenf
    logical :: localIOflag
    character(len=3) :: localOrder
    character(len=CHEM_MAXSTR) :: datafile
    real(CHEM_KIND_R4), dimension(:),     allocatable :: recvbuf
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: buf3d
    real(chem_kind_R8), pointer, dimension(:,: ) :: lat2dr8,lon2dr8
    real(chem_kind_R4), pointer, dimension(:,: ) :: lat2d,lon2d,latbuf,lonbuf
    integer bsize(3)
!    write(6,*)'write3dr4',maxval(farray),minval(farray)

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(de=de, tile=tile, tileComm=tileComm, &
      localIOflag=localIOflag, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, &
      its=its, ite=ite, jts=jts, jte=jte, lon=lon2dr8,lat=lat2dr8,rc=localrc)
!      its=its, ite=ite, jts=jts, jte=jte, rc=localrc)
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
    bsize = (/ ite-its+1, jte-jts+1,nk /)
    allocate(recvbuf(lbuf), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot allocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    recvbuf = 0._CHEM_KIND_R4

    call chem_comm_reduce(reshape(buf3d, (/lbuf/)), recvbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    allocate(lonbuf(its:ite,jts:jte), stat=localrc)
    allocate(latbuf(its:ite,jts:jte), stat=localrc)
    allocate(lat2d(its:ite,jts:jte))
    allocate(lon2d(its:ite,jts:jte))
    lat2d=0.
    lon2d=0.
    lat2d(ids:ide,jds:jde)=lat2dr8
    lon2d(ids:ide,jds:jde)=lon2dr8
    latbuf=0.
    lonbuf=0.
    call chem_comm_reduce(lat2d, latbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
    call chem_comm_reduce(lon2d, lonbuf, CHEM_COMM_SUM, comm=tileComm, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    if (localIOflag) then

      call raqmschem_io_file_name(datafile, filename, tile, pathname=path)
      lenf=len_trim(filename)
      if(filename(lenf-2:lenf).eq.'.nc')then
        buf3d=reshape(recvbuf,bsize)
!        write(6,*)'buf3d write',maxval(buf3d),minval(buf3d),'shape',shape(buf3d)
!        call flush(6)
        call chem_io_file_writenc3d(datafile, filename,buf3d,latbuf,lonbuf, pos=pos, time=time, varname=varname, rc=localrc)
      else
        call chem_io_file_write(datafile, recvbuf, pos=pos, rc=localrc)
      endif
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

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

      call raqmschem_io_file_name(datafile, filename, tile, pathname=path)

      call chem_io_file_write(datafile, recvbuf, pos=pos, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      write(6,'("chem_io_write: tile=",i2,2x,a," - min/max = "2g16.6)') tile, &
        trim(datafile), minval(recvbuf), maxval(recvbuf)
    end if

    deallocate(buf3d, recvbuf, stat=localrc)
    if (chem_rc_test((localrc /= 0), msg="Cannot deallocate read buffer", &
        file=__FILE__, line=__LINE__, rc=rc)) return

  end subroutine chem_io_write_3DR8



end module chem_io_mod
