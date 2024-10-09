      subroutine caldyi(nstep   ,dtime   ,mdbase  ,msbase  ,mbdate  ,
     $                  mbsec   ,mdcur   ,mscur   ,mcdate  ,mcsec   ,
     $                  calday  )
C-----------------------------------------------------------------------
C
C Compute current date and day information for history tape header.
C Compute current julian day (including fraction) for sun angle 
C computations and time interpolation of boundary dataset information.
C One year is defined as 365 days.
C
C Computational notes: 
C
C 86400 is the number of seconds in 1 day.
C
C Dividing an integer by 10**n has the effect of right-shifting the           
C decimal digits n positions (ex: 861231/100 = 008612).
C
C mod(integer,10**n) has the effect of extracting the rightmost     
C n decimal digits of the integer (ex: mod(861231,10000) = 1231).
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Boville, August 1992
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: caldyi.F,v 1.2 1995/10/14 21:10:08 rosinski Exp $
c $Author: rosinski $
c
      use raqmschem_pmgrid_mod, only : iam,masterproc
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer nstep        ! current time step
      real*8 dtime           ! length of time step (seconds)
      integer mdbase       ! base day of run
      integer msbase       ! base seconds of base day
      integer mbdate       ! base date (yymmdd format) of run
      integer mbsec        ! base seconds of base date
C
C Output arguments
C
      integer mdcur        ! current day
      integer mscur        ! current seconds of current day
      integer mcdate       ! current date (yymmdd format)
      integer mcsec        ! current seconds of current date
      real calday          ! current julian day including fraction
C
C---------------------------Local workspace-----------------------------
C
      integer nsecs        ! run time in seconds
      integer ndays        ! ndays + nyears = number of day changes since start
      integer nyears       ! ndays + nyears = number of day changes since start
      integer mcyear       ! current year of current date (0-99)
      integer mcmnth       ! current month of current date (1-12)
      integer mcday        ! current day of current date (1-31)
      integer jday         ! Julian day (1-365)
      integer mbmnth       ! base month (1-12) validity check
      integer mbday        ! base day of base date (1-31) validity check
      integer ndm(12,2)      ! number of days per month
      integer jdcon(12)    ! convert month index to julian day
      integer il  ! 1 if not leap 2 if leap year
c
      save ndm,jdcon
      data ndm/31,28,31,30,31,30,31,31,30,31,30,31
     *,31,29,31,30,31,30,31,31,30,31,30,31/
      data jdcon/0,31,59,90,120,151,181,212,243,273,304,334/
C
C-----------------------------------------------------------------------
C
C Check validity of input data
C
!      if(masterproc)print *,'mbdate',mbdate,' /10000 ',mbdate/10000
      if(mod(mbdate/10000,4).eq.0)then
        il=2
!         if(masterproc)print *,'leap year'
      else
        il=1
!         if(masterproc)print *,'not leap year'
      endif
      mbmnth = mod(mbdate,10000)/100
      mbday = mod(mbdate,100)
!      if(masterproc)then
!      write(6,*)'mbmnth',mbmnth,'mbday',mbday
      !endif
      if (mbmnth.lt.1 .or. mbmnth.gt.12) then
        write(6,*)' CALDYI: Invalid base month input:',mbmnth
        write(6,*)'mbdate=',mbdate,'nstep=',nstep,'dtime=',dtime
        call endrun
      end if
      if (mbday.lt.1 .or. mbday.gt.ndm(mbmnth,il)) then
        write(6,*)' CALDYI: Invalid base day of base date input:',
     *mbday,' il ',il,ndm(mbmnth,il)
        call endrun
      end if
      if (msbase.lt.0 .or. msbase.ge.86400) then
        write(6,*)' CALDYI: Invalid base seconds(msbase):',msbase
        call endrun
      end if
      if (mbsec.lt.0 .or. mbsec.ge.86400) then
        write(6,*)' CALDYI: Invalid base seconds(mbsec):',mbsec
        call endrun
      end if
C
      nsecs = nstep*nint(dtime)
C
C First: current day, seconds
C
      mdcur = mdbase + (nsecs+msbase)/86400
      mscur = mod((nsecs+msbase),86400)
C
C Next: current date, seconds.
C Uncommenting the next line will have the effect of modifying  
C nsecs to include base day and base seconds.  This is the way             
C the current date was computed in CCM1.
C
C      nsecs = nsecs + mdbase*86400 + msbase

      ndays = (nsecs+mbsec)/86400
      nyears = ndays/365
      ndays = mod(ndays,365)
C
C Current seconds of current date
C
      mcsec = mod(nsecs+mbsec,86400)
C
C Initialize current year, month, day.
C
      mcyear = mbdate/10000 + nyears
      mcmnth = mod(mbdate,10000)/100
      mcday = mod(mbdate,100) + ndays
C
C Now loop through months, converting yy, mm, and ddd to yymmdd.
C ex: 791235 becomes 800104. 190001370 becomes 19010105.
C
   10 if (mcday.gt.ndm(mcmnth,il)) then
        mcday = mcday - ndm(mcmnth,il)
        mcmnth = mcmnth + 1
        if (mcmnth.eq.13) then             ! add a year
          mcyear = mcyear + 1
          mcmnth = 1
        end if
        go to 10
      end if
      mcdate = mcyear*10000 + mcmnth*100 + mcday
c      if(masterproc)print *,'caldyi mcdate out ',mcdate
C
C Convert current month, day, seconds to Julian day + fraction
C
      jday = jdcon(mcmnth) + mcday
      calday = float(jday) + float(mcsec)/86400.
C
      return
      end

      subroutine caldyi2(nstep   ,dtime   ,mdbase  ,msbase  ,mbdate  ,
c    $                  mbsec   ,mdcur   ,mscur   ,mcdate  ,mcsec   ,
     $                  mbsec   ,mcdate  ,
     $                  mcyear,mcmnth,mcday)
C-----------------------------------------------------------------------
C
C Compute current date and day information for history tape header.
C Compute current julian day (including fraction) for sun angle 
C computations and time interpolation of boundary dataset information.
C One year is defined as 365 days.
C
C Computational notes: 
C
C 86400 is the number of seconds in 1 day.
C
C Dividing an integer by 10**n has the effect of right-shifting the           
C decimal digits n positions (ex: 861231/100 = 008612).
C
C mod(integer,10**n) has the effect of extracting the rightmost     
C n decimal digits of the integer (ex: mod(861231,10000) = 1231).
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Boville, August 1992
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: caldyi.F,v 1.2 1995/10/14 21:10:08 rosinski Exp $
c $Author: rosinski $
c
      use raqmschem_pmgrid_mod, only : masterproc
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer nstep        ! current time step
      real dtime           ! length of time step (seconds)
      integer mdbase       ! base day of run
      integer msbase       ! base seconds of base day
      integer mbdate       ! base date (yymmdd format) of run
      integer mbsec        ! base seconds of base date
C
C Output arguments
C
      integer mdcur        ! current day
      integer mscur        ! current seconds of current day
      integer mcdate       ! current date (yymmdd format)
      integer mcsec        ! current seconds of current date
      real calday          ! current julian day including fraction
C
C---------------------------Local workspace-----------------------------
C
      integer nsecs        ! run time in seconds
      integer ndays        ! ndays + nyears = number of day changes since start
      integer nyears       ! ndays + nyears = number of day changes since start
      integer mcyear       ! current year of current date (0-99)
      integer mcmnth       ! current month of current date (1-12)
      integer mcday        ! current day of current date (1-31)
      integer jday         ! Julian day (1-365)
      integer mbmnth       ! base month (1-12) validity check
      integer mbday        ! base day of base date (1-31) validity check
      integer ndm(12,2)      ! number of days per month
      integer jdcon(12)    ! convert month index to julian day
      integer il ! 2 leap year 1 not leap year
c
      save ndm,jdcon
      data ndm/31,28,31,30,31,30,31,31,30,31,30,31,
     *31,29,31,30,31,30,31,31,30,31,30,31/
      data jdcon/0,31,59,90,120,151,181,212,243,273,304,334/
C
C-----------------------------------------------------------------------
C
C Check validity of input data
C
      if(mod(mbdate/10000,4).eq.0)then
        il=2
!         if(masterproc)print *,'leap year'
      else
        il=1
      endif
      mbmnth = mod(mbdate,10000)/100
      mbday = mod(mbdate,100)
      if (mbmnth.lt.1 .or. mbmnth.gt.12) then
        write(6,*)' CALDYI2: Invalid base month input:',mbmnth
        write(6,*)'mbdate=',mbdate,'nstep=',nstep,'dtime=',dtime
        call endrun
      end if
      if (mbday.lt.1 .or. mbday.gt.ndm(mbmnth,il)) then
        write(6,*)' CALDYI: Invalid base day of base date input:',
     *mbday,' il ',il,ndm(mbmnth,il)
        call endrun
      end if
      if (msbase.lt.0 .or. msbase.ge.86400) then
        write(6,*)' CALDYI: Invalid base seconds(msbase):',msbase
        call endrun
      end if
      if (mbsec.lt.0 .or. mbsec.ge.86400) then
        write(6,*)' CALDYI: Invalid base seconds(mbsec):',mbsec
        call endrun
      end if
C
      nsecs = nstep*nint(dtime)
C
C First: current day, seconds
C
      mdcur = mdbase + (nsecs+msbase)/86400
      mscur = mod((nsecs+msbase),86400)
C
C Next: current date, seconds.
C Uncommenting the next line will have the effect of modifying  
C nsecs to include base day and base seconds.  This is the way             
C the current date was computed in CCM1.
C
C      nsecs = nsecs + mdbase*86400 + msbase

      ndays = (nsecs+mbsec)/86400
      nyears = ndays/365
      ndays = mod(ndays,365)
C
C Current seconds of current date
C
      mcsec = mod(nsecs+mbsec,86400)
C
C Initialize current year, month, day.
C
      mcyear = mbdate/10000 + nyears
      mcmnth = mod(mbdate,10000)/100
      mcday = mod(mbdate,100) + ndays
C
C Now loop through months, converting yy, mm, and ddd to yymmdd.
C ex: 791235 becomes 800104. 190001370 becomes 19010105.
C
   10 if (mcday.gt.ndm(mcmnth,il)) then
        mcday = mcday - ndm(mcmnth,il)
        mcmnth = mcmnth + 1
        if (mcmnth.eq.13) then             ! add a year
          mcyear = mcyear + 1
          mcmnth = 1
        end if
        go to 10
      end if
      mcdate = mcyear*10000 + mcmnth*100 + mcday
!      if(masterproc)print *,'caldyi2 mcdate out ',mcdate
C
C Convert current month, day, seconds to Julian day + fraction
C
      jday = jdcon(mcmnth) + mcday
      calday = float(jday) + float(mcsec)/86400.
C
      return
      end
      subroutine endrun
      use raqmschem_comm_mod
      call raqmschem_comm_abort(1,' date error')
      return
      end subroutine endrun
      subroutine killit(text)
      use raqmschem_comm_mod
      character*(*) text
      call raqmschem_comm_abort(2,text)
      return
      end subroutine killit
