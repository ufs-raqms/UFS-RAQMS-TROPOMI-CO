      module exchange_o3vmr_mod

      implicit none
      real*4, allocatable  :: o3vmr_dyn(:,:,:)

      contains

      subroutine init_o3vmr_dyn
      use raqmschem_pmgrid_mod, only : plon, plat, plev
      allocate (o3vmr_dyn(plon,plat,plev))
      return
      end subroutine init_o3vmr_dyn

      end module exchange_o3vmr_mod
