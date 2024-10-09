c
c $Id: commss.h,v 1.3 1995/11/03 21:03:17 rosinski Exp $
c $Author: rosinski $
c
C
C Character variables associated with Mass Store pathnames
C
      common/commsc/nhpath  ,nrpath  ,ncdata  ,bndti   ,bndtvs  ,
     $              bndtvo  ,nrfil   ,nsmvn   ,nrmvn   ,msscom  ,
     $              nswrps  ,lcroot
C
c     character*72  nhpath   ! MSS pathname for history tapes
c     character*72  nrpath   ! MSS pathname for restart files
c      character*120  nhpath   ! MSS pathname for history tapes
      character*240  nhpath   ! MSS pathname for history tapes
c      character*120  nrpath   ! MSS pathname for restart files
      character*240  nrpath   ! MSS pathname for restart files
c      integer lnhpath
c      parameter (lnhpath=120)
c      parameter (lnhpath=240)
      character*80  ncdata   ! MSS pathname for initial dataset
      character*80  bndti    ! MSS path for time-inv boundary dataset
C
      character*120  bndtvs   ! MSS path for time-variant sst dataset
      character*120  bndtvo   ! MSS path for time-variant ozone dataset
c     character*80  bndtvs   ! MSS path for time-variant sst dataset
c     character*80  bndtvo   ! MSS path for time-variant ozone dataset
      character*22   nrfil    ! Current file name for regen dataset
      character*8   nsmvn    ! Virtual volume name for history tapes
C
      character*8   nrmvn    ! Virtual volume name for restart data
      character*80  msscom   ! MSS comment field
      character*8   nswrps   ! MSS write password
      character*40  lcroot   ! Prepend to MSS paths for local disk name
C
C Non-character variables associated with Mass Store pathnames
C
      common/commss/nhpthl  ,nrpthl  , cnvrt2nc,
     *lsminithoj
C
      integer  nhpthl        ! Length of nhpath,history tape pathname
      integer  nrpthl        ! Length of nrpath,restart data pathname
      integer  irt           ! Mass Store retention period, history tapes
      integer  rirt          ! Mass Store retention time, restart data
C
      logical cnvrt2nc       ! Flag to convert history tapes to netCDF
      logical lsminithoj      ! flag if true lsminidata is hoj if false
c                              ascii
C
#ifdef LOGENTRY 
$Log: commss.h,v $
#endif








