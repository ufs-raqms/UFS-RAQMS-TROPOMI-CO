module chem_domain_mod

  use chem_types_mod

!------ legacy
! use mpp_mod,         only : mpp_npes
! use mpp_domains_mod, only : domain2d, mpp_define_domains, mpp_get_compute_domain
! use mpp_domains_mod, only : domain1d
!------ legacy end

  implicit none

  type chem_domain_type
    integer :: tile = 0
    integer :: tileCount = 0
    ! -- DE bounds
    integer :: ids   = 0
    integer :: ide   = 0
    integer :: jds   = 0
    integer :: jde   = 0
    ! -- memory bounds
    integer :: ims  = 0
    integer :: ime  = 0
    integer :: jms  = 0
    integer :: jme  = 0
    ! -- tile boundaries
    integer :: its  = 0
    integer :: ite  = 0
    integer :: jts  = 0
    integer :: jte  = 0
    ! -- number of vertical levels (i: interface, l: model)
    integer :: ni   = 0
    integer :: nl   = 0
    integer :: ns   = 0
    integer :: nt   = 0
    real(CHEM_KIND_R8), pointer :: lon(:,:)  => null()
    real(CHEM_KIND_R8), pointer :: lat(:,:)  => null()
    real(CHEM_KIND_R8), pointer :: area(:,:) => null()
  end type chem_domain_type

! private

! public :: chem_domain_type


  ! -- atmospheric grid parameters
  integer :: nip                 =  0      ! # of icosaedral cells
  integer :: nvl                 =  0      ! Number of vertical native levels
  integer :: nvlp                =  0      ! # of isobaric vertical levels - ex.  1000-25 hPa
  integer :: nvlp1               =  0      ! # of vertical levels ( = layers+1)

  integer, parameter :: nvl_chem =  55     ! Number of vertical native levels

  ! -- domain decomposition
! type(domain2d) :: chem_domain

  integer :: ids, ide, jds, jde, kds, kde
  integer :: ims, ime, jms, jme, kms, kme
  integer :: its, ite, jts, jte, kts, kte

  ! -- PET list
  integer, dimension(:), allocatable :: chem_pelist

  integer, allocatable :: indx(:), map(:)

  public

contains

  subroutine chem_domain_get(domain, &
   ids, ide, jds, jde,  &
   its, ite, jts, jte,  &
   ims, ime, jms, jme,  &
   ni,  nl,  ns, nt,  tile)
    type(chem_domain_type), intent(in)  :: domain
    integer, optional,      intent(out) :: ids, ide, jds, jde
    integer, optional,      intent(out) :: its, ite, jts, jte
    integer, optional,      intent(out) :: ims, ime, jms, jme
    integer, optional,      intent(out) :: ni, nl, ns, nt
    integer, optional,      intent(out) :: tile

    ! -- local variables

    ! -- begin
    if (present(ids))  ids  = domain % ids
    if (present(ide))  ide  = domain % ide
    if (present(jds))  jds  = domain % jds
    if (present(jde))  jde  = domain % jde
    if (present(its))  its  = domain % its
    if (present(ite))  ite  = domain % ite
    if (present(jts))  jts  = domain % jts
    if (present(jte))  jte  = domain % jte
    if (present(ims))  ims  = domain % ims
    if (present(ime))  ime  = domain % ime
    if (present(jms))  jms  = domain % jms
    if (present(jme))  jme  = domain % jme
    if (present(ni))   ni   = domain % ni
    if (present(nl))   nl   = domain % nl
    if (present(ns))   ns   = domain % ns
    if (present(nt))   nt   = domain % nt
    if (present(tile)) tile = domain % tile

  end subroutine chem_domain_get


  subroutine chem_domain_setup()
  end subroutine chem_domain_setup

#if 0
  subroutine chem_domain_setup(domain)
    type(domain2d), optional, intent(in) :: domain

    type(domain1d) :: domain_1d
    ! -- begin

    if (present(domain)) then
      chem_domain = domain
    else
      call mpp_define_domains( (/ 1, nvl, 1, nip /), (/ 1, mpp_npes() /), chem_domain)
    end if

    ! set up trivial indices
    call mpp_get_compute_domain(chem_domain, kts, kte, jts, jte)

    ids =   1 ; ide =   1 ; jds =   1 ; jde = nip ; kds =   1 ; kde = nvl + 1
    ims =   1 ; ime =   1 ; jms = jts ; jme = jte ; kms =   1 ; kme = nvl + 1
    its =   1 ; ite =   1 

  end subroutine chem_domain_setup

#endif

end module chem_domain_mod
