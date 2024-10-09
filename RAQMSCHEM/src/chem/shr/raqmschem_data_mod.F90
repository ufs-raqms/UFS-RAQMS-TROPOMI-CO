module raqmschem_data_mod
  use chem_rc_mod
  use chem_types_mod

  implicit none

  type chem_data_type
    ! -- input
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: griddx
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: griddy
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: griddxwithhalo
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: griddywithhalo
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: area
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: fcor
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: vort
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: absvort
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: potvort
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: aod
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: aodg5
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: aodgsi
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: aodgsibcoc
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: aodgsidust
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: aodgsissalt
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: aodgsisulf
    real(CHEM_KIND_R4), dimension(:,:), allocatable :: aodincgsi
!    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: ugrid,vgrid
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: ext_bcoc,ext_dust,ext_seas,ext_tot,ext_sulf
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: ext_3d,ext_3dg
!   for GBBEPX
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: plume             ! fire info - MODIS & GBBEPx
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: plume_hr             ! fire info - MODIS & GBBEPx
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: bbco_hr ! fire info GbbEPX
    integer num_plume_data

  end type chem_data_type

  private

  public :: chem_data_type
  public :: chem_data_destroy
contains
  subroutine chem_data_destroy(data, rc) 
    use raqmschem_pmgrid_mod, only : iam 
    type(chem_data_type)           :: data
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS
!    if(iam.eq.0)then
!      write(6,*)'at a chem destroy'
!      call flush(6)
!    endif
    if (allocated(data%fcor))then
      deallocate(data%fcor)
    endif
    if (allocated(data % area)) then
      deallocate(data % area, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
#if 0
    if (allocated(data % griddx)) then
      deallocate(data % griddx, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % griddy)) then
      deallocate(data % griddy, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % griddxwithhalo)) then
      deallocate(data % griddxwithhalo, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % griddywithhalo)) then
      deallocate(data % griddywithhalo, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
#endif
    if (allocated(data % vort)) then
      deallocate(data % vort, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % absvort)) then
      deallocate(data % absvort, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
    if (allocated(data % potvort)) then
      deallocate(data % potvort, stat=localrc)
      if (chem_rc_test((localrc /= 0), file=__FILE__, line=__LINE__, rc=rc)) return
    end if
!    if (allocated(data % ugrid))then
!      deallocate(data%ugrid)
!    endif
!    if (allocated(data % vgrid))then
!      deallocate(data%vgrid)
!    endif
!    if(iam.eq.0)then
!      write(6,*)'at b chem destroy'
!      call flush(6)
!    endif
     if (allocated(data%ext_bcoc))then
       deallocate(data%ext_bcoc)
     endif
     if (allocated(data%ext_dust))then
       deallocate(data%ext_dust)
     endif
     if(allocated(data%ext_seas))then
       deallocate(data%ext_seas)
     endif
     if(allocated(data%ext_tot))then
       deallocate(data%ext_tot)
     endif
     if(allocated(data%ext_sulf))then
       deallocate(data%ext_sulf)
     endif
!    if(iam.eq.0)then
!      write(6,*)'at c chem destroy'
!      call flush(6)
!    endif
     if(allocated(data%aod))then
       deallocate(data%aod)
     endif
     if(allocated(data%aodg5))then
       deallocate(data%aodg5)
     endif
     if(allocated(data%aodgsi))then
       deallocate(data%aodgsi)
     endif
     if(allocated(data%aodgsibcoc))then
       deallocate(data%aodgsibcoc)
     endif
     if(allocated(data%aodgsidust))then
       deallocate(data%aodgsidust)
     endif
     if(allocated(data%aodgsissalt))then
       deallocate(data%aodgsissalt)
     endif
     if(allocated(data%aodgsisulf))then
       deallocate(data%aodgsisulf)
     endif
     if(allocated(data%aodincgsi))then
       deallocate(data%aodincgsi)
     endif
    if(iam.eq.0)then
      write(300+iam,*)'at d chem destroy data%ext_3d',allocated(data%ext_3d)
      call flush(300+iam)
    endif
     if(allocated(data%ext_3d))then
       deallocate(data%ext_3d)
     endif
     if(allocated(data%ext_3dg))then
       deallocate(data%ext_3dg)
     endif
    if(iam.eq.0)then
      write(6,*)'at e chem destroy plume'
      call flush(6)
    endif
     if (allocated(data%plume))then
       deallocate(data%plume)
     endif
    if(iam.eq.0)then
      write(6,*)'at f chem destroy bottom'
      call flush(6)
    endif
  end subroutine chem_data_destroy

end module raqmschem_data_mod
