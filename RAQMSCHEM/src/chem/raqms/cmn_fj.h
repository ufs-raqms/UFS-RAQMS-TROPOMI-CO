!
!******************************************************************************
!  cmn_fj.h -- Header file containing parameters and common
!  blocks used to interface between RAQMS chemistry and UC-Irvine 
!  Fast-J photolysis programs.
!  
!  Based on code from Oliver Wild (9 July 1999)
!
!  NOTES:
!  (1) Uses Fortran 77 declarations for parameters and variables
!  
!  hyl, 11/13/02
!
!******************************************************************************
!
! Array size in altitude
      integer  lpar

! Variables for number of layers (jpnl), and number of photolysis rxns (jppj),
! and max # of photolytic reactions (jpmax)
      integer  jpnl, jppj, jpmax
      
!set lpar=jpnl=nl (see comm_3d_chem for nl)
!use number here, don't use nl (nl is used in both jv_mie.h and comm_3d_chem!!)
!(hyl,02/28/03)
      parameter(lpar=63,jpnl=63)   !Number of levels in CTM and 
!      parameter(lpar=64,jpnl=64)   !Number of levels in CTM and 
                                   !max number of levels requiring chemistry
                                   !set jpnl=lpar for RAQMS regional (hyl)
      parameter(jppj=53)           !Number of photolytic reactions supplied
      parameter(jpmax=55)          !max # of photolytic reactions

