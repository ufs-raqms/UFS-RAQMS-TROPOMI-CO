!#**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! miscellaneous debuging routines for CBMZ and MOSAIC
!**********************************************************************************  
! Sep 5, 2017 - Adapted from module_peg_util in FIM, R. Montuoro

module opt_peg_util_mod

  use chem_comm_mod, only : chem_comm_abort, chem_comm_log

  implicit none

  public


contains


!-----------------------------------------------------------------------
        subroutine peg_debugmsg( lun, level, str )
!
! when lun >  0, writes "str" to unit "lun"
! when lun <= 0, passes "str" on to wrf_debug
!
        implicit none
! subr arguments
        integer, intent(in) :: lun, level
        character(len=*), intent(in) :: str
! local variables
        integer n

        n = max( 1, len_trim(str) )
        if (lun .gt. 0) then
            write(lun,'(a)') str(1:n)
        else
          call chem_comm_log( "WARNING: " // str(1:n) )
        end if
        return
        end subroutine peg_debugmsg


!-----------------------------------------------------------------------
        subroutine peg_message( lun, str )
!
! when lun >  0, writes "str" to unit "lun"
! when lun <= 0, passes "str" on to wrf_message
!
        implicit none
! subr arguments
        integer, intent(in) :: lun
        character(len=*), intent(in) :: str
! local variables
        integer n

        n = max( 1, len_trim(str) )
        if (lun .gt. 0) then
            write(lun,'(a)') str(1:n)
        else
          call chem_comm_log( "WARNING: " // str(1:n) )
        end if
        return
        end subroutine peg_message


!-----------------------------------------------------------------------
        subroutine peg_error_fatal( lun, str )
!
! when lun >  0, writes "str" to unit "lun"
! then (always)  passes "str" on to wrf_error_fatal
!
        implicit none
! subr arguments
        integer, intent(in) :: lun
        character(len=*), intent(in) :: str
! local variables
        integer n

        n = max( 1, len_trim(str) )
        if (lun .gt. 0) write(lun,'(a)') str(1:n)
        call chem_comm_abort( msg=str(1:n) )
        return
        end subroutine peg_error_fatal


!-----------------------------------------------------------------------
end module opt_peg_util_mod
