module raqmschem_iodata_mod
  use chem_rc_mod
  use chem_types_mod
  use chem_io_mod
  use chem_comm_mod
  use chem_model_mod, only : chem_model_get => raqmschem_model_get
  use chem_model_mod, only : chem_config_type => raqmschem_config_type
  use chem_model_mod, only : chem_model_domain_get
  use raqmschem_config_mod
  use chem_data_mod, only : chem_data_type

  implicit none 

  private

  public :: raqmschem_backgd_init
  public :: raqmschem_backgd_read
  public :: raqmschem_backgd_write
contains
  subroutine raqmschem_backgd_init(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount
    integer :: ids, ide, jds, jde
    type(chem_data_type),   pointer :: data   => null()
    type(chem_config_type), pointer :: config => null()

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

      ! -- srcco_totl
      if (.not.allocated(data % srcco_totl)) then 
!        write(6,*)'allocate srcco_totl',ids,ide,'jds',jds,jde
!        call flush(6)
        allocate(data % srcco_totl(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
      endif
      ! -- srcnind
      if (.not.allocated(data % srcnind)) then 
        allocate(data % srcnind(ids:ide,jds:jde), stat=localrc)
        if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
      endif
    end do
    return
  end subroutine raqmschem_backgd_init
  subroutine raqmschem_backgd_read(bigendian,rc)
    integer, optional, intent(out) :: rc
    logical,optional,intent(in) :: bigendian

    ! -- local variables
    integer :: localrc
    integer :: de, deCount, localpe, tile
    type(chem_data_type),   pointer :: data => null()
    type(raqmschem_config_type), pointer :: config => null()
!    real(CHEM_KIND_R4),allocatable :: srcco_totl(:,:),srcnind(:,:)
    integer ids,ide,jds,jde,recrange(2)
    recrange=3

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_comm_get(localpe=localpe, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
!    write(6,*)'decount',decount

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, tile=tile, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      call chem_model_domain_get(de=de, ids=ids, ide=ide, jds=jds, jde=jde)
!      write(6,*)'srcco_totl',shape(data%srcco_totl),'ids',ids,ide,'jds',jds,jde
!      call flush(6)


        ! -- dust erosion factors
        call chem_io_read('C96_edgar_HTAP_CO_emi_molpercm2persec_2010_9.dat', data%srcco_totl, path=trim(config % emi_inname), de=de, recrange=recrange,bigendian=bigendian,rc=localrc)
!        write(6,*)'chem_io_read',localrc
!        call flush(6)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," co - min/max = "2g16.6)') localpe, de, &
          tile, minval(data%srcco_totl), maxval(data%srcco_totl)
        call chem_io_read('C96_edgar_HTAP_NOX_emi_molpercm2persec_2010_9.dat', data%srcnind, path=trim(config % emi_inname), de=de, recrange=recrange,bigendian=bigendian,rc=localrc)
!        write(6,*)'chem_io_read',localrc
!        call flush(6)
        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        write(6,'("chem_backgd_read: PET:",i4," DE:",i2," tile=",i2," nox - min/max = "2g16.6)') localpe, de, &
          tile, minval(data%srcnind), maxval(data%srcnind)
    end do 
    return
    end subroutine raqmschem_backgd_read
    subroutine raqmschem_backgd_write(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: de, deCount, localpe, tile 
    type(chem_data_type),   pointer :: data => null()
    type(chem_config_type), pointer :: config => null()

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    call chem_comm_get(localpe=localpe, rc=localrc)
    if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return

    do de = 0, deCount-1
      call chem_model_get(de=de, data=data, tile=tile, rc=localrc)
      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return


!      call chem_io_writenc('srcco_totl.nc', data % srcco_totl, path=trim(config % emi_outname), de=de, rc=localrc)
      call chem_io_write('srcco_totl.nc', data % srcco_totl, path=trim(config % emi_outname), de=de, rc=localrc)

      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," co - min/max = "2g16.6)') localpe, de, &
        tile, minval(data % srcco_totl), maxval(data % srcco_totl) 
      call chem_io_write('srcnind.nc', data % srcnind, path=trim(config % emi_outname), de=de, rc=localrc)

      if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
      write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," nox - min/max = "2g16.6)') localpe, de, &
        tile, minval(data % srcnind), maxval(data % srcnind) 
#if 0
!      write(6,*)'allocated griddx ',allocated(data%griddx)
!      call flush(6)
      if(allocated(data%griddx))then
        call chem_io_write('griddx.nc', data % griddx, path=trim(config % emi_outname), de=de, rc=localrc)

        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," griddx - min/max = "2g16.6)') localpe, de, &
          tile, minval(data % griddx), maxval(data % griddx) 
        call chem_io_write('griddy.nc', data % griddy, path=trim(config % emi_outname), de=de, rc=localrc)

        if (chem_rc_check(localrc, file=__FILE__, line=__LINE__, rc=rc)) return
        write(6,'("chem_backgd_write: PET:",i4," DE:",i2," tile=",i2," griddy - min/max = "2g16.6)') localpe, de, &
         tile, minval(data % griddy), maxval(data % griddy) 
       endif
#endif
     end do
     return
     end subroutine raqmschem_backgd_write

end module raqmschem_iodata_mod
