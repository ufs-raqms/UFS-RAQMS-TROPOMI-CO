
c tks,  I made some changes on 4-20-09 following code from Murali's verson
c of chem_driver that included updates for bromine. All of my changes
c are marked with date and 'tks'
#include <options.h>
      subroutine chem_driver(jday)
      use raqmschem_data_mod

#include <choosechem.h>
#include <choosevan.h>

! This module used to pass o3vmr to dynamics
      use exchange_o3vmr_mod

#ifdef USEHARVARDWETDEP
      use raqms_wetdep_interface_mod, only:
#ifdef WETD_USE
     *   nsol_fam, wetd_cd, wetd_rw, aphys_sum, nowet_sum,
#else
     *   nsol_fam,  aphys_sum, nowet_sum,
#endif
     *   cmflx_sum, pcpdiag2d, navgwet
      use tracerid_mod, only:
     *   idthno3, idth2o2, idtch2o, idtmp, idthcl, idthf,
     *   idthno4, idthbr, idtethooh, idtrooh, idtmgly,
     *   idtrip, idtprdp, idtxmrp, idtvrp,
     *   idtox, idtnoy, idtcly, idtbry,
     *   idtbcpo, idtbcpi, idtocpo, idtocpi,
     *   idtsal1, idtsal2, idtsal3, idtsal4,
     *   idtdms, idtso2, idtso4, idtmsa,
     *   idtdst1, idtdst2, idtdst3, idtdst4, idtdst5
      USE WETSCAV_MOD,  ONLY : NSOL
#endif

#ifdef DO_AEROSOL
      use mo_raqms_gocart_interface,  only : startup_gocart
#endif

c
c Addition to write out pbl height JAA 2005/03/14
      use raqmschemcomm_mod, only : pblht
      use raqmschemcomm_mod, only : nl_clim,nl_clim_dyn,nlnair,nlnlgt
c
       use raqmschem_pmgrid_mod, only : nhtuw
!      use mpimcomm, only :rad,dxdeg,dydeg,aarea,areasum,nhtuw
#ifdef MD5SUM
      use mpimcomm, only : checkmd5sum
#endif
!      use mpimcomm, only : rsynchist
#ifdef USEPKAPPA
      use raqmschem_pmgrid_mod, only : pdashk
#endif
      use raqmschem_pmgrid_mod, only : beglat,endlat,iam,masterproc,jbeg,nr=>plat,nc,nl=>plev
      use raqmschem_pmgrid_mod, only : bbedirf,bbedira,nstepat
      
#ifdef RES11FROM22
      use modnter2to1
#endif
c 1-30-08
#ifdef DOVDCHEM
      use vdchemmod, only : outputvdchem
#endif
      use raqmschemcomm_mod, only : o3ften_save,o3dten_save,o3dep_save,ho2_inst
      use raqmschemcomm_mod, only : phytendco,phytendn2o,phytendch4
      use raqmschemcomm_mod, only : noydep_save,codep_save,o3vmr_inst,oh_inst
      use raqmschemcomm_mod, only : no_inst,bro_inst,jo1d_inst,jno_inst
#ifdef CLD6HR
      use raqmschemcomm_mod, only : cld6hr
#endif
      use raqmschemcomm_mod, only : noyften_save,coften_save,codten_save
#ifdef CLDTAU
      use raqmschemcomm_mod, only : tau6hr,lnox_save,noydten_save,colnox_save
#else
      use raqmschemcomm_mod, only : lnox_save,noydten_save,colnox_save
#endif
!      use raqmschemcomm_mod, only : cloud3d,tauxcl3d,tauxci3d,navg,pcolnox
      use raqmschemcomm_mod, only : navg,pcolnox
#ifdef WETDEPDIAG
      use raqmschemcomm_mod, only : ch4ften_save,ch4dten_save,n2often_save,n2odten_save
#endif
#ifdef NEWTERM
!      use raqmschemcomm_mod, only : o3ften_save_new,o3dten_save_new,
      use raqmschemcomm_mod, only : ch4ften_save,ch4dten_save,n2often_save,n2odten_save
      use raqmschemcomm_mod, only : ch4sfcten_save,n2osfcten_save
      use raqmschemcomm_mod, only : o3ten_part
#endif
      use raqmschemcomm_mod, only : no2ten
      use raqmschemcomm_mod, only : xgrid,ygrid,bbco_d,bbco_n,firehtkm_d,firehtkm_n
      use raqmschemcomm_mod, only : bbco_hr
#ifdef DOMPI
      use mpishorthand
#endif
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
      use raqmschemcomm_mod, only : nl_clim,nl_clim_dyn,nlnair,nlnlgt,julday,gmt
      use raqmschemcomm_mod, only : bb_nox_adjust
c
c include following com files to get output directory and file label
#include <commss.h>
#include <comcdate.h>

c
c include following com file to get date for selecting daily fire emissions
#include <comtim.h>

c
#include <comm_3d_chem>
      include 'combined_constants'
      include 'combined_model_logic'

c save Ox photochemical tendency
#ifndef CLD6HR
      real *4 cld6hr(nc,beglat:endlat,nl)
#endif
      character *2 clvl
      character *32 checksum
      real lnoxout(nc,beglat:endlat,nl),noxsource,noxsourcej,areatot
      real areaearth
      type(chem_data_type),   pointer :: data

#ifdef DAILYFIRES
c
c Additions for reading daily fire emissions
      character cmon*2,cday*2,cyr*4,cyyyymmdd*8,cyyyymmddold*8,
     *newdate*8
!      common /todayfires/ bbco_d(nc,nr), firehtkm_d(nc,nr),
!     &                    bbco_n(nc,nr), firehtkm_n(nc,nr)
!      real *4 bbco_d, firehtkm_d, bbco_n, firehtkm_n
      integer nlondum,nlatdum
      real *4 dumlon(nc,nr),dumlat(nc,nr)
c      common/bbinfo/ bbedir
c      character bbedir*132
      save cyyyymmddold
      data cyyyymmddold/' '/
#endif

#ifdef DO_AEROSOL
      real xlon1d(nc), xlat1d(nr)
      save xlon1d, xlat1d
#endif
      common /chemnfo_24/ imon
      character*256 filename


c
c
c flag to read auxiliary files at model startup
      logical dochemread,printchem,lopenrsync
      save dochemread,printchem
!      data dochemread/.true./,printchem/.true./
      data dochemread/.true./,printchem/.false./
      integer nt6hr,imonold,itry
      save imonold
      data imonold/-99999/
!      character firehoj*10,fireoutdir*256
      character fireoutdir*256
      integer lenpath
      

#ifdef USEHARVARDWETDEP
c To save space, write wetdep diagnostics for some selected species only
      parameter(noutdiag=2)
      integer ioutdiag(noutdiag)
      save ioutdiag
      real *4 r4out(nc,beglat:endlat,nl,noutdiag)
#endif
      integer rclocal
c 4-20-09 tks

c define time step
c
#ifndef CLD6HR
      cld6hr=0.0
#endif
!      firehoj=' '
!      write(6,*)'top chem_deiver',mbdate
!      call flush(6)
!      write(300+iam,*)'top chem_driver ',nstepat
!      call flush(300+iam)
      fireoutdir=' '
!      call getenv('FIREHOJ',firehoj)
      call getenv('FIREOUTDIR',fireoutdir)
      dt=dtime
      ntspd=86400/nint(dt)
      gmtinc=dt/3600.
      nt6hr=nint(6.*3600./dtime)
!      write(6,*)'top chemdriver ',nstepat
!      call flush(6)
!      !if(masterproc)then
!        write(6,*)'nhtuw',nhtuw,'nt6hr',nt6hr,'nstep',nstep
!        call flush(6)
!      endif
!      if(mod(nstep,nt6hr).eq.0)then
!        if(masterproc)print *,'nt6hr',nt6hr,'nstep',nstep,'6hr','nhtuw',nhtuw
!      endif
c
c define date
c
      julday=jday
!      if(masterproc)then
!        if(mod(nstep,12).eq.0)then
!          print *,'julday,gmt=',julday,gmt
!        endif
!      !endif
c 4-20-09 tks
      dy_jdgmt=float(julday)+gmt/24.0
!      write(6,*)'dy_jdgmt',dy_jdgmt,'julday',julday,'gmt',gmt
!      call flush(6)
!      !if(masterproc)then
!        write(6,*)'dy_jdgmt',dy_jdgmt,'julday',julday,'gmt',gmt,'nstep',nstep
!      endif
c
c read auxiliary files at model startup
c
!      write(100+iam,*)'dochemread',dochemread
!      call flush(100+iam)
!      write(6,*)'dochemread',dochemread
!      call flush(6)
      if(dochemread)then
       dochemread=.false.

c Run/model parameters
       hydrostatic=.true.
       dohet=.true.
       fixtend=.false.
       usefastj=.true.
       cloudphot=.true.         !photolysis cloud flag: false for clear-sky
       if(masterproc)then
!       print *,'nhpath(nhpthl): ',trim(nhpath)
!       print *,'cdate(lencdate): ',trim(cdate)
c
c  read albedo data
crbp currently call at start up should call each day
c
c       print *,'reading albedo data',julday
c       call rd_albedo(julday)
c
c  read larc datasets
c
       endif
!      if(masterproc)then
!         write(6,*)'caldyi one ',nstep,'mbdate',mbdate,mbsec,'julday',julday
!         call flush(6)
!      endif
!      write(6,*)'nstep',nstep,' dtime fire ',dtime
      call caldyi(nstep,dtime,nnbdat,nnbsec,mbdate,mbsec,
     *mbcurfv,mscurfv,mcdatefv,mcsecfv,caldayfv)
      if(masterproc)then
         write(6,*)'at fires nstep',nstep,'mcsecfv',mcsecfv
         call flush(6)
      endif
!      if(masterproc)then
!           write(6,*)'for fires caldayfv',caldayfv,'nstep',nstep,'mcdateafv',mcdatefv
!      endif
        iyr=mcdatefv/10000
        if(iyr.lt.40)then
          iyr=2000+iyr
        else 
          iyr=1900+iyr
        endif
       write(cyr,'(i4.4)')iyr
!       write(100+iam,*),'init_larc'
!       call flush(100+iam)
       call init_larc
!       write(100+iam,*),'did init_larc'
!       call flush(100+iam)

       !Initial call to Fast-J to set things up
       if (usefastj) then
         nlevdum=nl
         nphojdum=iphot
         if(masterproc)then
         write(6,*)'Initializing fastj2 ',nlevdum,nphojdum
         call flush(6)
         endif
         CALL INPHOT (nlevdum,nphojdum)
       endif
c
c save Ox photochemical tendency
       navg=0
!       if(mype.eq.0)then
!         write(6,*)'navg=0 aa'
!       endif
       do i=1,nc
        do j=beglat,endlat
         o3dep_save(i,j) = 0.
         if(allocated(bb_nox_adjust))then
           bb_nox_adjust(i,j)=0.0
         endif
         noydep_save(i,j) = 0.
         codep_save(i,j) = 0.
         colnox_save(i,j)=0.0
#ifdef NEWTERM
         ch4sfcten_save(i,j)=0.0
         n2osfcten_save(i,j)=0.0
#endif
         do l=1,nl
          no2ten(i,j,l)=0.0
          o3ften_save(i,j,l) = 0.
          o3dten_save(i,j,l) = 0.
          phytendco(i,l,j)=0.
          phytendn2o(i,l,j)=0.
          phytendch4(i,l,j)=0.
#ifdef WETDEPDIAG
          ch4ften_save(i,j,l)=0.0
          ch4dten_save(i,j,l)=0.0
          n2often_save(i,j,l)=0.0
          n2odten_save(i,j,l)=0.0
#endif
#ifdef NEWTERM
!          o3ften_save_new(i,j,l) = 0.
!          o3dten_save_new(i,j,l) = 0.
          ch4ften_save(i,j,l)=0.0
          ch4dten_save(i,j,l)=0.0
          n2often_save(i,j,l)=0.0
          n2odten_save(i,j,l)=0.0
          o3ten_part(i,j,l,:)=0.0
#endif
!        ajl this is filled does not need to be zero
!          o3vmr_inst(i,j,l) = 0.
          oh_inst(i,j,l) = 0.
          ho2_inst(i,j,l) = 0.
          no_inst(i,j,l) = 0.
          bro_inst(i,j,l) = 0.
          jo1d_inst(i,j,l) = 0.
          jno_inst(i,j,l) = 0.
          noyften_save(i,j,l) = 0.
          noydten_save(i,j,l) = 0.
          coften_save(i,j,l) = 0.
          codten_save(i,j,l) = 0.
#ifdef CLD6HR
          cld6hr(i,j,l) = 0.
#endif
#ifdef CLDTAU
          tau6hr(i,j,l) = 0.
#endif
          lnox_save(i,j,l) = 0.
         enddo
        enddo
       enddo

       call init_o3vmr_dyn

c       open(unit=666,file='oxtend.dat',
c     &  form='unformatted',status='unknown')
c       write(666) nc,nr,nl

c       open(unit=555,file='colo3.dat',
c     &  form='unformatted',status='unknown')
c       write(555) nc,nr,nl

c Define longitudes (xgrid) and latitudes (ygrid)
c Caution: xgrid and ygrid are both 2D (nc,nr) for modularity
!      define xgrid and ygrid outside
!       call getgrid_uwhyb
!       write(100+iam,*)'did getgrid_uwhyb'
!       call flush(100+iam)

#ifdef USEHARVARDWETDEP
c To save space, write wetdep diagnostics for some selected species only
       if(masterproc)then
       print*,'In chem driver, wetdep species IDs:'
       print*,'idthno3, idth2o2, idtch2o, idtmp, idthcl, idthf'
       print*, idthno3, idth2o2, idtch2o, idtmp, idthcl, idthf
       print*,'idthno4, idthbr, idtethooh, idtrooh, idtmgly'
       print*, idthno4, idthbr, idtethooh, idtrooh, idtmgly
       print*,'idtrip, idtprdp, idtxmrp, idtvrp'
       print*, idtrip, idtprdp, idtxmrp, idtvrp
       print*,'idtox, idtnoy, idtcly, idtbry'
       print*, idtox, idtnoy, idtcly, idtbry
       print*,'idtbcpo, idtbcpi, idtocpo, idtocpi'
       print*, idtbcpo, idtbcpi, idtocpo, idtocpi
       print*,'idtsal1, idtsal2, idtsal3, idtsal4'
       print*, idtsal1, idtsal2, idtsal3, idtsal4
       print*, 'idtdms, idtso2, idtso4, idtmsa'
       print*, idtdms, idtso2, idtso4, idtmsa
       print*, 'idtdst1, idtdst2, idtdst3, idtdst4, idtdst5'
       print*, idtdst1, idtdst2, idtdst3, idtdst4, idtdst5
       call flush(6)
       endif
c Ensure that number of species defined here in ioutdiag = noutdiag
       ioutdiag(1) = idtox
       ioutdiag(2) = idtnoy
#endif

!       write(100+iam,*)'use xgrid',shape(xgrid),shape(ygrid)
!       write(100+iam,*)'lb',lbound(xgrid),ubound(xgrid)
!       call flush(100+iam)
!       call mpibarrier(mpicom)
#ifdef DO_AEROSOL
       do i=1,nc
         xlon1d(i) = dble(xgrid(i,1))
       enddo
       do j=1,nr
         xlat1d(j) = dble(ygrid(1,j))
       enddo
!       write(6,*)'startup_gocart'
!       call flush(6)
!       call mpibarrier(mpicom)
!       write(100+iam,*)'call startup_gocart'
!       call flush(100+iam)
       call startup_gocart(emdeplabel,xlon1d,xlat1d,cyr)
#endif

      endif        !End of model startup section
!      write(6,*)'end of model startup'
!      call flush(6)
!      write(100+iam,*)'end of model startup'
!      call flush(100+iam)
!       call mpibarrier(mpicom)
cc
c JAA 20080225 Modification to build LWI flag based on daily sea ice.
c zlwigrd now built in getlwi using model physics data; prior to this it was
c built in getdiag from climatology.  Also store tskin and sea ice thickness.
c Also, move rd_albedo outside of "model startup" loop so that surface
c albedo is temporally interpolated daily instead of only at startup.
c
!      write(6,*)'getlwi'
!      call flush(6)
!       call mpibarrier(mpicom)
!      call getlwi
!      if(masterproc)then
!        if(printchem)then
          !write(6,*)'reading albedo data'
!          call flush(6)
!        endif
!#ifdef PRINTJULDAY
!        print *,'rd_albedo read ',julday
!#endif
!      endif
!      if(masterproc)then
!          write(6,*) 'for albedo ajl julday',julday
!      endif
!      write(6,*)'call rd_ALBEDO',julday
!      call flush(6)
!      write(200+iam,*)'call rd_albedo',julday
!      call flush(200+iam)
      call rd_albedo(julday)
!      write(6,*)'getdiag'
      !call flush(6)
      call getdiag
!      write(100+iam,*)'did getdiag',nstepat
!      write(300+iam,*)'did getdiag',nstepat
!      call flush(100+iam)
!      call flush(300+iam)
!      write(6,*)'did getdiag'
!      call flush(6)
!       call mpibarrier(mpicom)

c
c get date for daily fire emissions and read biomass burning CO source
c
!      if(masterproc)then
!         write(6,*)'caldyi two ',nstep,'mbdate',mbdate,mbsec,'julday',julday
!         call flush(6)
!      endif,nstepat
!      write(6,*)'caldyi two dtime',dtime
      call caldyi(nstep,dtime,nnbdat,nnbsec,mbdate,mbsec,
     *mbcurfv,mscurfv,mcdatefv,mcsecfv,caldayfv)
!      if(masterproc)then
!           write(6,*)'for fires caldayfv',caldayfv,'nstep',nstep,'mcdateafv',mcdatefv
!      endif
        iyr=mcdatefv/10000
        if(iyr.lt.40)then
          iyr=2000+iyr
        else
          iyr=1900+iyr
        endif
        imon=mod(mcdatefv,10000)/100
        iday=mod(mcdatefv,100)
        if(masterproc)then
          if(mod(nstep,12).eq.0)then
            print*,'Actual date/secs: ',iyr,imon,iday,mcsecfv
          endif
        endif
#ifdef DO_AEROSOL
#ifdef AEROSOL_SULFATE
c This subroutine fills emission arrays with values for this month
      if(imon.ne.imonold)then
        imonold=imon
        call set_bcond_su(.TRUE.,imon,xlon1d,xlat1d)
      endif
#endif
#endif
c     
      if(masterproc)then
        if(printchem)then
          write(6,*)'calling chemistry ',julday,gmt
          call flush(6)
        endif
      endif
!      write(6,*)'call chem_driver_3d',dy_jdgmt
!      call flush(6)
!      write(6,*)'bottom  chemdriver ',nstepat
!      call flush(6)
!      write(300+iam,*)'bottom  chemdriver ',nstepat
!      call flush(300+iam)
      return
      end subroutine chem_driver
!     call chem_driver_3d(dy_jdgmt)
      subroutine chem_diagnostics
      use exchange_o3vmr_mod
      use raqmschem_data_mod
      use raqmschem_pmgrid_mod, only : beglat,endlat,iam,masterproc,jbeg,nr=>plat,nc,nl=>plev
      use raqmschem_pmgrid_mod, only : nhtuw,nstepat
#ifdef DOVDCHEM
      use vdchemmod, only : outputvdchem
#endif
      use raqmschemcomm_mod, only : o3ften_save,o3dten_save,o3dep_save,ho2_inst
      use raqmschemcomm_mod, only : phytendco,phytendn2o,phytendch4
      use raqmschemcomm_mod, only : noydep_save,codep_save,o3vmr_inst,oh_inst
      use raqmschemcomm_mod, only : no_inst,bro_inst,jo1d_inst,jno_inst
#ifdef CLD6HR
      use raqmschemcomm_mod, only : cld6hr
#endif
      use raqmschemcomm_mod, only : noyften_save,coften_save,codten_save
#ifdef CLDTAU
      use raqmschemcomm_mod, only : tau6hr,lnox_save,noydten_save,colnox_save
#else
      use raqmschemcomm_mod, only : lnox_save,noydten_save,colnox_save
#endif
!      use raqmschemcomm_mod, only : cloud3d,tauxcl3d,tauxci3d,navg,pcolnox
      use raqmschemcomm_mod, only : navg,pcolnox
#ifdef WETDEPDIAG
      use raqmschemcomm_mod, only : ch4ften_save,ch4dten_save,n2often_save,n2odten_save
#endif
#ifdef NEWTERM
!      use raqmschemcomm_mod, only : o3ften_save_new,o3dten_save_new,
      use raqmschemcomm_mod, only : ch4ften_save,ch4dten_save,n2often_save,n2odten_save
      use raqmschemcomm_mod, only : ch4sfcten_save,n2osfcten_save
      use raqmschemcomm_mod, only : o3ten_part
#endif
      use raqmschemcomm_mod, only : xgrid,ygrid,bbco_d,bbco_n,firehtkm_d,firehtkm_n
#ifdef DOMPI
      use mpishorthand
#endif
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
      use raqmschemcomm_mod, only : nl_clim,nl_clim_dyn,nlnair,nlnlgt
!      use raqmschemcomm_mod, only : pblht,taucldfrc_liq,taucldfrc_ice
      use raqmschemcomm_mod, only : pblht,taucldfrc
      use raqmschemcomm_mod, only : no2ten
      integer i,j,l
      real *4 pbldiag(nc,beglat:endlat)
c
c include following com files to get output directory and file label
#include <commss.h>
#include <comcdate.h>

c
c include following com file to get date for selecting daily fire emissions
#include <comtim.h>

c
#include <comm_3d_chem>
      include 'combined_constants'
      include 'combined_model_logic'

#include <choosechem.h>
#include <choosevan.h>


cc accumulate mean cloud and optical depth diagnostics
cc also, map o3vmr into dynamics array
      do l=1,nl
!        if(maxval(taucldfrc(:,:,l))/=0.0)then
!          write(200+iam,*)' ave ',k,maxval(taucldfrc(:,:,l))
!          call flush(200+iam)
!        endif
        do j=beglat,endlat

          do i=1,nc
#ifdef CLD6HR
          cld6hr(i,j,l) = cld6hr(i,j,l)+cloud3d(i,j,l)
#endif
c          tau6hr(i,j,l) = tau6hr(i,j,l)+tauxcl3d(i,j,l)+tauxci3d(i,j,l)
c         fix ajl 7/31/2013 tauxcl3d and tauxci3d are defined everywhere even if no cloud 
c         so need to multiply by cloud fraction to get correct cloud tau
#ifdef CLDTAU
!          tau6hr(i,j,l) = tau6hr(i,j,l)+(tauxcl3d(i,j,l)+tauxci3d(i,j,l))*cloud3d(i,j,l)
!          tau6hr(i,j,l) = tau6hr(i,j,l)+taucldfrc_liq(i,j,l)+taucldfrfc_ice(i,j,l)
          tau6hr(i,j,l) = tau6hr(i,j,l)+taucldfrc(i,j,l)
#endif
c          if(i.eq.15.and.j.eq.109)then
c            if(tau6hr(i,j,l).ne.0.0)then
c              write(6,*)'tau6hr=',l,tau6hr(i,j,l),'tauxc3d',tauxcl3d(i,j,l),tauxci3d(i,j,l),'cloud3d',cloud3d(i,j,l)
c            endif
c          endif
          o3vmr_dyn(i,j-beglat+jbeg,l) = o3vmr_inst(i,j,l)
        enddo
       enddo
      enddo


cc save Ox photochemical tendency and other diagnostics
      navg = navg+1
#ifdef USEHARVARDWETDEP
      if (navg .ne. navgwet) then
        write(6,*) 'Error in diagnostics: navg ne navgwet',navg,navgwet
        call flush(6)
        call killit('useharv')
!        stop
      endif
#endif
c      testoxtend = mod(chemutc,6.)           !chemutc is hours
c      if (abs(testoxtend) .lt. 0.01) then
c       if(nstep<=1)then
c         call writewrffileo3vmr(o3vmr_inst)
c       endif
#ifdef SAVEALL
      if(.true.)then
#else
!      if(mod(nstep,nt6hr).eq.0)then
!     try using nhtuw instead so that each uwh has it
!      if(masterproc)then
!        write(6,*)'nstep',nstep,'nhtuw',nhtuw
!        call flush(6)
!      endif
      if(masterproc)then
      write(6,*)'nstep',nstep,'nhtuw',nhtuw
      call flush(6)
      endif
!      write(6,*)' nstep ',nstep,' nhtuw ',nhtuw
!      call flush(6)
      if(mod(nstep,nhtuw).eq.0.and.nstep.ne.0)then
#endif
#ifdef SAVELNOX
        noxsource=0.0
        areatot=0.0
        do j=beglat,endlat
          noxsourcej=0.0
          do i=1,nc
            noxsourcej=noxsourcej+colnox_save(i,j)*data@area(i,j)
          end do
          noxsource=noxsource+noxsourcej*aarea(j-beglat+jbeg)
          areatot=areatot+float(nc)*aarea(j-beglat+jbeg)
        end do
#endif
#if 0
#ifdef DOMPI
       call reduce0(noxsource)
       call reduce0(areatot)
#endif
       if(masterproc)then
         write(6,*)'global noxsource',noxsource
         radian=2.*pi/360.
         avenoxsource=noxsource/dt/float(navg)/areatot
         write(6,*)'avenoxsource',avenoxsource,'areatot',areatot
         areaearth=rad**2*4.*pi
         write(6,*)'areaearth',areaearth,'prod',
     *   avenoxsource*areaearth
         write(6,*)'prod tg/year',
     *   avenoxsource*areaearth*86400.*365./1.e9
         noxsource=noxsource/dt/float(navg)
!     *   (rad*radian)**2.*dydeg*dxdeg

!        now per second
!        make teragrams/year colnox in kg/m**2 ?
         noxsource=noxsource*86400.*365./1.e9
         write(6,*)'noxsource TGNOX/yr',noxsource
       endif
#endif
!       write(6,*)iam,'phytendco',shape(phytendco),shape(phytendn2o),shape(phytendch4)
!       call flush(6)
!       call mpibarrier(mpicom)
        if(iam.eq.0)then
        write(6,*)'chem_diag nc',nc,beglat,endlat,'navg',navg,
     * 'nhtuw',nhtuw,'nstep',nstep
        endif
!       call flush(6)
!       write(6,*)'navg',navg,nstepat
       if(navg<1)navg=1
       do i=1,nc
        do j=beglat,endlat
         o3dep_save(i,j) = o3dep_save(i,j)/float(navg)
         noydep_save(i,j) = noydep_save(i,j)/float(navg)
         codep_save(i,j) = codep_save(i,j)/float(navg)
         colnox_save(i,j)=colnox_save(i,j)/float(navg)
c         if(i.eq.271.and.j>130.and.j<140)then
c           write(6,*)'codep_save ',j,codep_save(i,j),
c     *   ' o3dep ',o3dep_save(i,j),'noydep ',noydep_save(i,j)
c         endif
         pbldiag(i,j) = pblht(i,j)
         do l=1,nl
          o3ften_save(i,j,l) = o3ften_save(i,j,l)/float(navg)
          o3dten_save(i,j,l) = o3dten_save(i,j,l)/float(navg)
          phytendco(i,l,j)=phytendco(i,l,j)/float(navg)
          phytendn2o(i,l,j)=phytendn2o(i,l,j)/float(navg)
          phytendch4(i,l,j)=phytendch4(i,l,j)/float(navg)
#ifdef WETDEPDIAG
          ch4ften_save(i,j,l)=ch4ften_save(i,j,l)/float(navg)
          ch4dten_save(i,j,l)=ch4dten_save(i,j,l)/float(navg)
          n2often_save(i,j,l)=n2often_save(i,j,l)/float(navg)
          n2odten_save(i,j,l)=n2odten_save(i,j,l)/float(navg)
#endif
#ifdef NEWTERM
!          o3ften_save_new(i,j,l) = o3ften_save_new(i,j,l)/float(navg)
!          o3dten_save_new(i,j,l) = o3dten_save_new(i,j,l)/float(navg)
          ch4ften_save(i,j,l)=ch4ften_save(i,j,l)/float(navg)
          ch4dten_save(i,j,l)=ch4dten_save(i,j,l)/float(navg)
          n2often_save(i,j,l)=n2often_save(i,j,l)/float(navg)
          n2odten_save(i,j,l)=n2odten_save(i,j,l)/float(navg)
          o3ten_part(i,j,l,:)=o3ten_part(i,j,l,:)/float(navg)
          
#endif
          noyften_save(i,j,l) = noyften_save(i,j,l)/float(navg)
          noydten_save(i,j,l) = noydten_save(i,j,l)/float(navg)
          coften_save(i,j,l) = coften_save(i,j,l)/float(navg)
          codten_save(i,j,l) = codten_save(i,j,l)/float(navg)
#ifdef CLD6HR
          cld6hr(i,j,l) = cld6hr(i,j,l)/float(navg)
#endif
#ifdef CLDTAU
          tau6hr(i,j,l) = tau6hr(i,j,l)/float(navg)
#endif
          lnox_save(i,j,l) = lnox_save(i,j,l)/float(navg)
          no2ten(i,j,l)=no2ten(i,j,l)/float(navg)
         enddo
        enddo
       enddo
!       do l=1,nl
!        write(6,*)l,'chemdriver n2often',maxval(n2often_save(:,:,l)),minval(n2often_save(:,:,l))
!        call flush(6)
!       end do
       return
!       write(6,*)'did tend ',iam
!       call flush(6)
!       call mpibarrier(mpicom)
c       if(masterproc)then
c         write(6,*)'before oxtendnc'
c       endif
#ifdef DONETCDFHISTOUT
c        if(masterproc)then
c          write(6,*)'calle oxtendnc in chemdriver'
c        endif
       call oxtendncsave
c       call mpibarrier(mpicom)
c       write(6,*)'call wediagncsave'
c       call flush(6)
c       call mpibarrier(mpicom)
c       call wetdiagncsave
c       call aodncwrite
c       write(6,*)'did call wediagncsave'
c       call flush(6)
c       call mpibarrier(mpicom)
#endif
!       write(6,*)'did oxtendncsave ',iam
!       call flush(6)
!       call mpibarrier(mpicom)
       navg=0
!       if(mype.eq.0)then
!          write(6,*)'navg = 0'
!          call flush(6)
!       endif
!       write(6,*)iam,'save chem ',beglat,endlat
!       call flush(6)
       do i=1,nc
        do j=beglat,endlat
         o3dep_save(i,j) = 0.
         noydep_save(i,j) = 0.
         codep_save(i,j) = 0.
#ifdef NEWTERM
         ch4sfcten_save(i,j)=0.0
         n2osfcten_save(i,j)=0.0
#endif
         do l=1,nl
          o3ften_save(i,j,l) = 0.
          o3dten_save(i,j,l) = 0.
#ifdef WETDEPDIAG
          ch4ften_save(i,j,l)=0.0
          ch4dten_save(i,j,l)=0.0
          n2often_save(i,j,l)=0.0
          n2odten_save(i,j,l)=0.0
#endif
#ifdef NEWTERM
!          o3ften_save_new(i,j,l) = 0.
!          o3dten_save_new(i,j,l) = 0.
          ch4ften_save(i,j,l)=0.0
          ch4dten_save(i,j,l)=0.0
          n2often_save(i,j,l)=0.0
          n2odten_save(i,j,l)=0.0
          o3ten_part(i,j,l,:)=0.0
#endif
!         ajl this is filled does not need to be zero
!          o3vmr_inst(i,j,l) = 0.
          oh_inst(i,j,l) = 0.
          ho2_inst(i,j,l) = 0.
          no_inst(i,j,l) = 0.
          bro_inst(i,j,l) = 0.
          jo1d_inst(i,j,l) = 0.
          jno_inst(i,j,l) = 0.
          noyften_save(i,j,l) = 0.
          noydten_save(i,j,l) = 0.
          coften_save(i,j,l) = 0.
          codten_save(i,j,l) = 0.
#ifdef CLD6HR
          cld6hr(i,j,l) = 0.
#endif
#ifdef CLDTAU
          tau6hr(i,j,l) = 0.
#endif
          lnox_save(i,j,l) = 0.
          no2ten(i,j,l)=0.0
         enddo
         colnox_save(i,j)=0.
        enddo
       enddo
!       write(6,*)iam,'did save chem ',beglat,endlat
!       call flush(6)
!       call mpibarrier(mpicom)

#ifdef USEHARVARDWETDEP
C Write 6-hr averaged wet-removal, convective, precip diagnostics
       do i=1,nc

        do j=beglat,endlat
         do l=1,nl
          do n=1,nsol+nsol_fam
#ifdef WETD_USE
           wetd_cd(i,j,l,n) = wetd_cd(i,j,l,n)/float(navgwet)
           wetd_rw(i,j,l,n) = wetd_rw(i,j,l,n)/float(navgwet)
#endif
           aphys_sum(i,j,l,n) = aphys_sum(i,j,l,n)/float(navgwet)
           nowet_sum(i,j,l,n) = nowet_sum(i,j,l,n)/float(navgwet)
          enddo
          do n=1,2
           cmflx_sum(i,j,l,n) = cmflx_sum(i,j,l,n)/float(navgwet)
          enddo
         enddo
         do n=1,2
           pcpdiag2d(i,j,n) = pcpdiag2d(i,j,n)/float(navgwet)
         enddo
        enddo
       enddo
       !write(6,*)iam,'didd save 2'
!       call flush(6)
!       call mpibarrier(mpicom)
#ifdef DONETCDFHISTOUT
c       if(masterproc)then
c         write(6,*)'calle oxtendnc'
c       endif
c       call oxtendncsave
c       call mpibarrier(mpicom)
c       write(6,*)'call wediagncsave'
c       call flush(6)
c       call mpibarrier(mpicom)
c       if(masterproc)then
c         write(6,*)'call wetdiagncsave in chemdriver'
c       endif
       call wetdiagncsave
c       call aodncwrite
c       write(6,*)'did call wediagncsave'
c       call flush(6)
c       call mpibarrier(mpicom)
#endif
       navgwet=0
       do i=1,nc
        do j=beglat,endlat
         do l=1,nl
          do n=1,nsol+nsol_fam
#ifdef WETD_USE
           wetd_cd(i,j,l,n) = 0.
           wetd_rw(i,j,l,n) = 0.
#endif
           aphys_sum(i,j,l,n) = 0.
           nowet_sum(i,j,l,n) = 0.
          enddo
          do n=1,2
           cmflx_sum(i,j,l,n) = 0.
          enddo
         enddo
         do n=1,2
           pcpdiag2d(i,j,n) = 0.
         enddo
        enddo
       enddo
#endif


      endif  !End of test for time to write out diagnostics
      printchem=.false.
!       write(6,*)iam,'bottom of chem_driver'
!       call flush(6)
!       call mpibarrier(mpicom)

      return 
      end subroutine chem_diagnostics
c-------------------------------------------------------------------------
      subroutine chem_driver_3d
      use raqmschem_pmgrid_mod, only : beglat,endlat,iam
      use raqmschemcomm_mod
      use raqmschemlocaltype_mod
#ifdef DOCHEMKK
      use raqmschem_pmgrid_mod,only : maxx1n
      use mpishorthand,only : mpimax
#endif
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
#include <choosechem.h>
#include <comm_3d_chem>
      include 'combined_constants'
c 4-20-09 tks
      logical lopenrysnc
      type(chemlocaltype) :: chemlocal
#ifdef DOCHEMKK
      real maxx1ntot
      maxx1n=0.0
#endif
c
c start chemistry loop now rewritten in preparation for parallelizing
c
ccccC$OMP parallel do
ccccC$OMP1 private(jjzap)
      do 160 jjzap=beglat,endlat                 ! latitude loop
c          print*,'160 loop jjzap=',jjzap
c          call flush(6)
c
c  map gridded dynamics to chemical arrays
c
c        print*,'calling map_to, ',jjzap
c        call flush(6)
!        write(100+iam,*)'call map_to_uwhyb',jjzap
!        call flush(100+iam)
!        call map_to_uwhyb(jjzap,chemlocal)
!        write(100+iam,*)'did call map_to_uwhyb',jjzap
!        call flush(100+iam)

!        write(6,*)'calling doloop, ',jjzap,'dy_jdgmt',dy_jdgmt,'kind',kind(dy_jdgmt)
!        call flush(6)
        call doloop(jjzap,chemlocal)
        !write(6,*)'did dollop',jjzap
!        call flush(6)
c
c map gridded chemistry back to dynamical arrays
c
!        write(100+iam,*)'call map_from_uwhyb',jjzap
!        call flush(100+iam)
!        write(6,*)'calling map_from, ',jjzap
!        call flush(6)
!        call map_from_uwhyb(jjzap,chemlocal)
        call deallocatechemlocal(chemlocal)
        call deallocatechemlocal2(chemlocal)

160   continue                          ! end of latitude loop
#ifdef DOCHEMKK
      call mpireduce0(maxx1n,maxx1ntot,mpimax)
      write(6,'(i3," maxx1n ",e12.4)')iam,maxx1n
      if(masterproc)then
          write(6,'("maxx1ntot ",e12.4)')maxx1ntot
      endif
#endif
      return
      end
c-------------------------------------------------------------------------
      subroutine doloop(jjzap,chemlocal)

#include <choosechem.h>
c
c add doloop subroutine for multiprocessor
c
      use seasalt_mod, only : get_aerosource_ss,idss1,idss2,idss3,idss4,sorcaero
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
#ifdef CALCLIGHTNINGNOX
!      use raqms_wetdep_interface_mod, only: ccthk, lconvcld
       use raqmschemcomm_mod, only : lconvcld,ccthk
#endif

c Additions to extract pbl height JAA 2004/09/02
      use raqmschem_pmgrid_mod, only : plond,plat,masterproc,iam,nc,nl=>plev
      use raqmschem_pmgrid_mod, only : nmonHTAP,nr=>plat,ibeg,begj,tile
      use raqmschem_pmgrid_mod, only : twolayer_emissions,dopr94,doareafix,usescale_factor

!      use buffer, only : pblht
      use raqmschemcomm_mod, only : pblht
      use raqmschem_pmgrid_mod,only : iam,iprnin,jprnin,skipphotoj,iamprn,nstepat,aerosol_ugpkg
!#ifdef DOCHEMJJ
      use raqmschem_pmgrid_mod,only : iatchem,jatchem,katchem,iprn,jprn,kprn,kprnin,tileprn
!#endif

#ifdef DO_AEROSOL
      use mo_raqms_gocart_interface, only : calc_depvel, get_aerosource
      use mo_parameters, only: nmx_raq
c
c module mo_raqms_gocart_vardef provides n_raq and idNNN
!      USE mo_raqms_gocart_vardef
!      use mo_tracer, only: tcmw
#endif
      use raqmschemcomm_mod, only: tcmw,idbc1,idbc2,idoc1,idoc2,idso4,ndust,naer
      use raqmschemcomm_mod, only: ns,nw,np,mx,ltb,no2ten
!      use ijprn
#ifdef COTRACERS
      use raqmschemcomm_mod, only : xcod25,sorccod5025
#ifdef XBBCOD25
     use raqmschemcomm_mod, only : bbcod25grd,xbbcod25,sorcbbcod5025
#endif
#endif
#ifdef COTRACERS50
      use raqmschemcomm_mod, only : bbcod50grd,xbbcod50
#endif
      use raqmschemcomm_mod, only : sorcbbcod25,sorccoanth25
      use raqmschemcomm_mod, only : xcoanth25,xbbcod25
      use raqmschemcomm_mod, only : o3ften_save,o3dten_save,o3dep_save,ho2_inst
      use raqmschemcomm_mod, only : phytendco,phytendn2o,phytendch4

      use raqmschemcomm_mod, only : noydep_save,codep_save,o3vmr_inst,oh_inst
      use raqmschemcomm_mod, only : no_inst,bro_inst,jo1d_inst,jno_inst
      use raqmschemcomm_mod, only : noyften_save,coften_save,codten_save
      use raqmschemcomm_mod, only : lnox_save,noydten_save,zlwigrd,colnox_save
      use raqmschemcomm_mod, only : spgrd,pgrd,slp,zsurf,tgrd,tskin,no2ten
      use raqmschemcomm_mod, only : zgrd,thgrd,zeupgrd,szagrd,albgrd,ugrd,vgrd
      use raqmschemcomm_mod, only : srcnind,srcnsoil,srcco_totl,srcethane,srcalkanone
      use raqmschemcomm_mod, only : srcpentane,srchexane,srcpropene,srcethene
      use raqmschemcomm_mod, only : srcisop,diurngrd,srcch2o,srcald,srcpropane
      use raqmschemcomm_mod, only : srcterp,bbefac,srcnair,rcraingrd
#ifdef AEROSOL_CARBON
      use raqmschemcomm_mod, only : srcbc,srcoc
#endif
      use raqmschemcomm_mod, only : scale_factor
!     ajl new aircraft sources
#ifdef AEROSOL_CARBON
      use raqmschemcomm_mod, only : srcbcair,srcocair
      use raqmschemcomm_mod, only : srcso2air
#endif
      use raqmschemcomm_mod, only : srccoair
      use raqmschemcomm_mod, only : dpo3,dphno3,dph2o2,dpco,ch4clim,srcbutane
      use raqmschemcomm_mod, only : qgrd,slghtng,dpno2
!      use raqmschemcomm_mod, only : cloud3d,tauxcl3d,tauxci3d
      use raqmschemcomm_mod, only : dpmgrd,pcolnox,navg,latgrd,longrd
      use raqmschemcomm_mod, only : xgrid,ygrid,bbco_d,bbco_n,firehtkm_d,firehtkm_n
      use raqmschemcomm_mod, only : bbco_hr
#ifdef WETDEPDIAG
      use raqmschemcomm_mod, only : ch4ften_save,ch4dten_save,n2often_save,n2odten_save
#endif
#ifdef NEWTERM
!      use raqmschemcomm_mod, only : o3ften_save_new,o3dten_save_new,prodlossnsave
      use raqmschemcomm_mod, only : ch4ften_save,ch4dten_save,n2often_save,n2odten_save
      use raqmschemcomm_mod, only : ch4sfcten_save,n2osfcten_save
!      use mpimcomm, only : nstepa
      use raqmschemcomm_mod, only : o3ten_part
#endif
      use raqmschemcomm_mod, only : julday,gmt
#ifdef DOEMISSION
      use emisnc, only : doemis,emisout,emis3dout
#endif
      use raqmschemlocaltype_mod
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
!      use raqmschemcomm_mod, only : nl_clim,nl_clim_dyn,nlnair,nlnlgt,taucldfrc_liq,taucldfrc_ice,areagrd
      use raqmschemcomm_mod, only : nl_clim,nl_clim_dyn,nlnair,nlnlgt,taucldfrc,areagrd
      use raqmschemcomm_mod, only : i1grd,j1grd,dist1grd,jppj,jpnl,dy_jdgmt
      use raqmschemcomm_mod, only : sourceno,bb_nox_scale_factor,bb_nox_adjust
      use raqmschemcomm_mod, only : emisco3d,ktopco,kbotco,doplumerise
      use raqmschemcomm_mod, only : colemisco_chem,colemisoc
      use raqmschemcomm_mod, only : docofire85percent,donofire85percent,use_fv3_aero

#include <comm_3d_chem>
      include 'combined_constants'
!      include 'cmn_fj.h'
      include 'ccwvl.h'
! ajl try 3/24/2019
!!!#include "jv_cmn.h"
c 4-20-09 tks

#ifdef CALCLIGHTNINGNOX
c Lightning flashrate param from Allen and Pickering, JGR 2002
      parameter(ap_a = 0.021, ap_b = -0.648, ap_c = 7.493,
     &          ap_d = -36.540, ap_e = 63.090)
!     actuall constanst fort Price and Rind 1993
!     fcg=0.02 dz > 14 km
!     fcg=fraction of cloud to ground
!     dz from cloud top to hightest layer T > 273K
!     min dz is 5.5km
      parameter(cgnox = 6.7e26)   ! molecules NO/flash
      parameter(ccnox = 6.7e25)   ! molecules NO/flash
#endif

      parameter(recm = 6.371e8)   ! radius of earth in cm
      common /chemnfo_24/ imon
      real(CHEM_KIND_R8) :: ox,xnoy,hno3,cly,xn2o5,h2o2,hcl,clno3,oclo,xn2o,
     1              f11,f12,ccl4,ch3cl,xmtcfm,bry,ch3br,h1301,h1211,
#ifdef ISOPRENE_PEROX
     2              rip,prdp,xmrp,ch4,hno4,hocl,vrp,
#else
     2              hf,cfclo,cf2o,ch4,hno4,hocl,h2otest,
#endif
     2              xno3,xno2,ch2o,ch3ooh,xco,hbr,brno3,hobr,brcl,cl2,
     3              xnat,xice,chem_tend,sednat_in,sedice_in
      common /blk2_24/ ox,xnoy,hno3,cly,xn2o5,h2o2,hcl,clno3,oclo,xn2o,
     1              f11,f12,ccl4,ch3cl,xmtcfm,bry,ch3br,h1301,h1211,
#ifdef ISOPRENE_PEROX
     2              rip,prdp,xmrp,ch4,hno4,hocl,vrp,
#else
     2              hf,cfclo,cf2o,ch4,hno4,hocl,h2otest,
#endif
     2              xno3,xno2,ch2o,ch3ooh,xco,hbr,brno3,hobr,brcl,cl2,
     3              xnat,xice,chem_tend(ikmult),sednat_in,sedice_in
c (JRO) add new common block for nmhc transported species
      common/blk_nmhc/c2h6,ald2,ethooh,pan,par,xonit,
     5              aone,rooh,xmgly,eth,xolet,xolei,xisop,
     6              xisoprd,prop_par,ch3oh,xmvk,xmacr,xmpan
      common/blk_nmhc2/hcn,ch3cn
#if (defined DO_AEROSOL || defined DO_AEROSOL_DEFINE)
      common/blk_aero0d/so2,so4aer,dms,xmsa,xno3aer,xnh3,xnh4aer,
     1      bc1,bc2,oc1,oc2,du1,du2,du3,du4,du5,ss1,ss2,ss3,ss4,ss5
      common/srfsnkaero/depso2,depso4aer,depdms,depxmsa,
     &      depxno3aer,depxnh3,depxnh4aer,depbc1,depbc2,
     &      depoc1,depoc2,depdu1,depdu2,depdu3,depdu4,
     &      depdu5,depss1,depss2,depss3,depss4,depss5,
     &      sorcso2,sorcso4aer,sorcdms,sorcxmsa,
     &      sorcxno3aer,sorcxnh3,sorcxnh4aer,sorcbc1,sorcbc2,
     &      sorcoc1,sorcoc2,sorcdu1,sorcdu2,sorcdu3,sorcdu4,
     &      sorcdu5,sorcss1,sorcss2,sorcss3,sorcss4,sorcss5
CCCMNN0508    include dv_hobr_ice
c 4-20-09 tks     &      ,dv_hobr_ice
!     &      ,dv_hobr_ice,br_ice_fctr
      common/aerodat/r_su,r_ssa,r_ssc,v_su,v_ssa,v_ssc,
     &      sa_su,sa_ssa,sa_ssc
#ifdef DO_AEROSOL
      real sorcaero(nmx_raq),ohgrd(nl),delz(nl),cldf(nl)
#endif
#endif
#ifdef DO_AEROSOL_BR2
      common/blk_aerobr2/br2
#endif
      common/srfbr/dv_hobr_ice,br_ice_fctr
      common /blk4_24/ r(ireac),hr(ihreac),phoj(iphot),dn(inum),
     1              zeta,pur,sad,prodloss(13)
c      common /prodlossnew/prodlossn(3),prodlossch4(3),prodlossn2o(3)
      real n2oprodloss(3),noyprodloss(3),coprodloss,ch4prodloss
      common /prodlossnew/oxprodloss(9),coprodloss(3),noyprodloss,ch4prodloss(3),n2oprodloss
      real oxtend,cotend,noytend,ch4tend,n2otend
      common /prodlosstend/oxtend,cotend,noytend,ch4tend,n2otend
      real oxsource,oxsink
      real(CHEM_KIND_R8) :: troplat,troptheta
      common /tropopause/ troplat(64),troptheta(64)
      real(CHEM_KIND_R8) :: xlsig,fixer,climacetone
      common /fixer_24/ xlsig(inl),fixer(3,ijg2,inl)
     1                 ,climacetone(19,17,2)
      real(CHEM_KIND_R8) :: depo3,dephno3,deph2o2,depch2o,
     1       depch3ooh,depno2,depco,sorcno,sorcco,drip,trht
c (JRO) additional species in srfsnk common block for nmhc species
     2       ,sorcc2h6,sorcpar,sorcaone,sorceth,
     3       sorcolet,sorcolei,sorcisop,sorcch2o,sorcald2,
     4       depald2,depethooh,deprooh,depmgly,drip_ch3ooh,
     5       sorcproppar,sorcch3oh
c Additional variables for biomass burning tracers
     6       ,dephcn,depch3cn,sorchcn,sorcch3cn

      common /srfsnk/ depo3,dephno3,deph2o2,depch2o,
     1       depch3ooh,depno2,depco,sorcno,sorcco,drip,trht
c (JRO) additional species in srfsnk common block for nmhc species
     2       ,sorcc2h6,sorcpar,sorcaone,sorceth,
     3       sorcolet,sorcolei,sorcisop,sorcch2o,sorcald2,
     4       depald2,depethooh,deprooh,depmgly,drip_ch3ooh,
     5       sorcproppar,sorcch3oh
c Additional variables for biomass burning tracers
     6       ,dephcn,depch3cn,sorchcn,sorcch3cn
c
c Add block chembc_var to pass surface concentrations of species 
c  that are not held constant, eg CH4 is now f(lat,time) JAA 200601
      real ch4_s
      common /chembc_var/ch4_s
      real,parameter :: airmw=28.97 ! added 7/22/2022

c cloud optical depth
cjaa      real*4 tau3d(nc,nr,nl)
c      tauxcl3d - water cloud extinction optical depth
c      tauxci3d - ice cloud extinction optical depth

c store 3d cloud info JAA 09/22/2004

      real sorclnox(nl)

      real rjs(nxsec)

      real zgeo(nl),zgeoh(nl),plog(nl),tc(nl),svo3(nl),
     1     svcolo3(nl),svcolo3v(nl)
!      integer,allocatable :: ltb(:)
      real dz2

      integer nmonth(12),ncount
      data nmonth /31,28,31,30,31,30,31,31,30,31,30,31/
      save nmonth
#ifdef DOINT
      real*8 gmt8,zsurf8,spgrd8
#ifndef SKIPISNAN
      logical isnan
#endif
#endif
#if (defined DOINT || defined DOINTCHEMB || defined DOINTCHEMM)
      real*8 t8,a8,p8,u8,v8,tskin8,tll8
#endif

c hyl, 03/11/03
      real*8 dbp(nl),dbz(nl),dbt(nl),dbox(nl),dbcolo3v(nl)
cjaa      real*8 doptd(nl)
CCCMNN0308      real*8 doptdl(nl),doptdi(nl),doptdbc(nl)
      real*8 doptdl(nl),doptdi(nl),aer_opd(naer,nl),dust_opd(ndust,nl)
      real odsclf(5),wwcon(5),dwcon(2),wcdu(5),rwet(5),qwet(5)
      real wwconr8(5),dwconr8(2)
      real wwconr4(5),dwconr4(2)
      real rh(5),aer_den(5),dust_den(5)
      real v_arsl(5,nl),s_arsl(5,nl),eff_r(5,nl)
      real zero8
      parameter (zero8=0.0)
      data aer_den/1.7e3,1.0e3,1.8e3,2.2e3,2.2e3/
      save aer_den
      data dust_den/2.5e3,2.65e3,2.65e3,2.65e3,2.65e3/
      save dust_den
      data rh/0.0,0.5,0.7,0.8,0.9/
      save rh
      integer indx(5)
      data indx/22,29,36,43,50/
      save indx
      real*8 zpj(jpnl,jppj)
      integer ii,jj
      real noxsource,noxsourcej,radian,avenoxsource
      real dob_to_den,dtor,zer_min,XMW_DRYAIR,XMW_H2O

      parameter(dob_to_den=2.637e19/1000.)     ! Dobson Units to column density
      parameter(dtor=3.1415926535/180.)        ! degrees to radians
      parameter(zer_min=1.)
      parameter(XMW_DRYAIR=28.964)
      parameter(XMW_H2O=18.0153)
      logical haveemis
CCCMNN 0318
      

c save Ox photochemical tendency

      common/dtchem/ timstep,tfac
      real(CHEM_KIND_R8) :: draw2d(nc)

#ifdef DAILYFIRES
c Additions for reading daily fire emissions
!      common /todayfires/ bbco_d(nc,nr), firehtkm_d(nc,nr),
!     &                    bbco_n(nc,nr), firehtkm_n(nc,nr)
!      real *4 bbco_d, firehtkm_d, bbco_n, firehtkm_n
#endif
C$OMP THREADPRIVATE (/blk2_24/)
C$OMP THREADPRIVATE (/blk_nmhc/)
C$OMP THREADPRIVATE (/blk4_24/)
C$OMP THREADPRIVATE (/srfsnk/)
C$OMP THREADPRIVATE (/chembc_var/)
C$OMP THREADPRIVATE (/blk_nmhc2/)
#ifdef DO_AEROSOL
C$OMP THREADPRIVATE (/blk_aero0d/)
C$OMP THREADPRIVATE (/srfsnkaero/)
C$OMP THREADPRIVATE (/aerodat/)
#endif
C$OMP THREADPRIVATE (/srfbr/)
      real rr,coefpr94
      integer errcount,idum
!      integer iprn,jprn
      type(chemlocaltype),target :: chemlocal
#include <chemlocaldefinepointerF.h>
#include <chemlocaldefinepointer2F.h>
!      write(6,*)'top doloop'
!      call flush(6)
!      write(200+iam,*)'top doloop ',jjzap,' nstepa ',nstepat
!      call flush(200+iam)
#include <chemlocalsetpointer.h>
#include <chemlocalsetpointer2.h>
      j=jjzap
      jatchem=j
      errcount=0

c
c
!       write(6,*)'top doloop ',jjzap
!       call flush(6)
!      write(6,*)'jjzap',jjzap,'tile',tile
!      call flush(6)
#ifdef FINDTOP
      if(j.eq.jprn+1)then
         write(6,*)'findtop',j
         call flush(6)
      endif
#endif
c
      timstep=dt
      tfac=1.
!      if(iam.eq.0)then
!         write(6,*)'timestep',timestep,kind(timestep)
!         write(6,*)'tfac',tfac,kind(tfac)
!      endif
#ifdef DO_AEROSOL
      ndtaero=nint(timstep)
      jaero=jjzap
#endif
c   imonhtap = imon if nmonHTAP = 12
      if(nmonHTAP.eq.12)then
        imonHTAP=imon
      else
        imonHTAP=1
      endif

c  ensure that altitude loop is top to bottom
c   define ltb indexing array - useful for accessing grid points
c   above/below current point, given arbitrary altitude ordering of grid data
c
       
c  get draw (entire latitude band)
!      write(6,*)'call doloop chemdraw',j
!      call flush(6)
      call chemdraw(j,draw2d)
!      write(100+iam,*)'did doloop chemdraw'
!      call flush(100+iam)

c
c check pbl height
c
c        write(6,*) 'PBL heights'
c        write(6,*) j,ygrid(1,j)
c        write(6,*) (pblht(i,j),i=1,nc)


!      write(6,*)'do 160 nc ',nc
!      call flush(6)
!      write(200+iam,*)'top nc ',jjzap,nc
!      call flush(200+iam)
      do 160 i=1,nc                      ! longitude loop
        iatchem=ibeg+i-1
!        if(tile.eq.1.and.iatchem.eq.83.and.jatchem.eq.25)then
!           write(6,*)'printpoint',iatchem,jatchem,'tile',tile
!        endif
!        if(nstepat==84.and.jjzap.eq.99)then
!          write(200+iam,*)'i',i
!          call flush(200+iam)
!        endif
#ifdef DIAGHTAP
        if(i.eq.120.and.j.eq.130)then
          write(6,*)'top of i j 120,130'
        endif
#endif
        rlon=xgrid(i,j)
!        write(6,*)'rlon',i,j,rlon,img
!        call flush(6)
c  find longitude index for nearest point in larc grid
        tol=999.
        do ii=1,img
          diff=abs(xlarc(ii)-rlon)
          if (diff .lt. tol) then
            klon=ii
            tol=diff
          endif
        enddo
c
        rlat = ygrid(i,j)
c  find latitude index for nearest point in larc grid
        tol=999.
        do jj=1,ijg2
          diff=abs(ylarc(jj)-rlat)
          if (diff .lt. tol) then
            jlat=jj
            tol=diff
          endif
        enddo

        sednat_in = 0.          !initialize NAT sedimentation into layer
        sedice_in = 0.          !initialize ice sedimentation into layer
        ktf = -1                !initialize tropopause flag
        ksf = -1                !initialize surface flag

        lbot=nl
        ltrop=1
        ptop=pgrd(i,j,ltb(1))

        do ip=1,nl
!          if(nstepat==84.and.jjzap.eq.99)then
!          if(nstepat==84.and.jjzap.eq.99)then
!            write(200+iam,*)'i',i,'ip',ip,'nl',nl
!            write(200+iam,*)'iatchem',iatchem,jatchem,katchem,' tile ',tile
!            call flush(200+iam)
!          endif
!          write(6,*)'ip',ip,nl
!          call flush(6)
c  fill T->B column arrays used in column ozone calculation
          l=ltb(ip)
!         here ip=1 is the top and l=nl at top
#ifdef DIAGHTAP
          if(i.eq.120.and.j.eq.130)then
             write(6,*)'ip',ip,'l',l,'p',pgrd(i,j,l)
          endif
#endif
          plog(ip)=alog(1000.*pgrd(i,j,l))
          tc(ip) = tgrd(i,j,l)
          dens=7.2431122e+18*pgrd(i,j,l)/tgrd(i,j,l)
          svo3(ip) = max(oxgrd(i,l)*dens,zer_min)
cJAA commented out because assim is putting negatives in some locations
c          if(svo3(ip).eq.zer_min)print *,'dens=',dens,'oxgrd=',
c     *     oxgrd(i,l),i,j,l,
c     *    'pgrd=',pgrd(i,j,l),'t=',tgrd(i,j,l)
c  find lowest grid point that is above ground
cjjj          if (pgrd(i,j,l) .le. spgrd(i,j)) lbot=ip
          if (zgrd(i,j,l) .ge. zsurf(i,j)) lbot=ip
!          if (iam.eq.7.and.i.eq.33.and.j.eq.49)then
!             write(6,*)'zgrd',l,zgrd(i,j,l),'zsurf',zsurf(i,j),'ip',ip,'lbot',lbot
!          Qendif
#ifdef DIAGFV3
        if(ibeg+i-1.eq.1.and.j.eq.1.and.tile.eq.1)then
          if (zgrd(i,j,l) .ge. zsurf(i,j))then
            write(6,*)'zgrd ',i,j,l,zgrd(i,j,l),'zsurf',zsurf(i,j),'lbot',lbot
          endif
          endif
#endif
c  find highest grid point that is in the troposphere (for lightning NOx)
!          if(i.eq.1.and.(j.eq.45.or.j.eq.1))then
!            write(6,*)'thgrd',thgrd(i,j,l),j,l,'troptheta',troptheta(jlat),jlat
!            call flush(6)
!          endif
          if (thgrd(i,j,l) .ge. troptheta(jlat)) ltrop=ip
!          if (thgrd(i,j,l) .ge. troptheta(jlat)) then
!            write(200+iam,*)'ip',ip,' l ',l,'thgrd',thgrd(i,j,l),'tropth',jlat,troptheta(jlat)
            !write(200+iam,*)'ltrop',ltrop
!          endif
          
!          if(ltrop.eq.nl)then
!            write(6,*)'thgrd',i,j,l,thgrd(i,j,l),'troptheta',troptheta(jlat),'jlat',jlat,'ip',ip,'ii',i+ibeg+1
!            write(6,*)'xgrid',i,j,xgrid(i,j),ygrid(i,j)
!            write(6,*)'lat',latgrd(i,j),'lon',longrd(i,j),'ltrop',ltrop
!            call flush(6)
!          endif
        enddo
        ltrop=ltrop+1
!          if(i.eq.1.and.j.eq.45)then
!              write(6,*)'ltrop',ltrop
!           endif

c        print *,'*******LOOK!!!*************'
c        print *,i,j,lbot,ptop,spgrd(i,j),zsurf(i,j)
c        print *,ltb(lbot)

c  zgeo is cell-center, zgeoh is cell-boundary (omitting top-of-column)
c   Define zgeoh and zgeo as T->B with units of km
c   Approximate zgeoh as midway between cell-centers
c   Update 200801: now have zgrd (cell-center) and zeupgrd (cell-edge, 
c    upper boundary, so does not include surface)
c    Note that zgeoh includes surface but not top of atmosphere!
c
        zgeo(1) = zgrd(i,j,ltb(1))/1000.
          
        do ip=2,lbot
          zgeo(ip) = zgrd(i,j,ltb(ip))/1000.
c          zgeoh(ip-1) = 0.5*(zgeo(ip)+zgeo(ip-1))
          zgeoh(ip-1) = zeupgrd(i,j,ltb(ip))/1000.
        enddo
        if(lbot<1)then
          write(6,*)'lbot error',lbot
          call flush(6)
        endif
        zgeoh(lbot) = zsurf(i,j)/1000.

cjaa Note: using cloud optical depth calculated in ccm3 instead of this
cjaa  calculation; just for diagnostics now, remove at cleanup
c  calculate cloud optical depth
#ifdef OLDTAUWATER
        call get_cloud_opdepth(i,j,ltb,zgeoh)
#endif

c  calculate overhead O3 column at each model level at this lon/lat
        sza_deg = szagrd(i,j)           ! solar zenith angle (deg)
        sza = sza_deg*dtor
        zeta = cos(sza)
!        if(i.eq.1.and.j.eq.1)then
!           write(6,*)'szagrd',szagrd(i,j),'zeta',zeta
!           call flush(6)
!        endif
c  for efficiency, only calculate in daytime (beware, not re-init at night!)
!        if(nstepat==84.and.jjzap.eq.99)then
!          write(200+iam,*)'0zeta',zeta
!          call flush(200+iam)
!        endif
        if (zeta .ge. -0.068) then
          ierr = 0

          call calc_colo3_col(sza,plog,tc,svo3,zgeo,ptop,lbot,
     +             klon,jlat,svcolo3,svcolo3v,chapno,ierr)
cJAA commented out because assim is putting negatives in some locations
c          if (ierr .gt. 0) then
c            print *,'h_o3',julday,sza*180./3.14159,rlat,rlon
c            print *,'svo3 ='
c            print *,svo3
c          endif
CCC MNN (9/8/2000)  COLUMN NO REQUIRED FOR OVERHEAD NO OPACITY.
CCC                 CONSTANT VERTICAL COLUMN AT THE TOP LEVEL OF THE MODEL
CCC                 IS ASSUMED  (3.0E14 MOLECULES/CM**2).  STRATOSPHERIC NO
CCC                 IS NOT INCLUDED AT PRESENT SINCE THIS DATA IS NOT AVAILABLE
CCC                 FOR RESTART. O3 SCALE HEIGHT IS USED FOR NO ALSO.
CCC
          slntno = 3.0e14 * chapno
CCC
        endif

! ------------------------------ Fast-J2 implementation ------------------------
! For the photolytic reactions included in the fastj2 scheme, use the rates
! from fastj2 instead of those from lookup tables (hyl, 11/12/2002)
        if(skipphotoj)then
          do ii=1,jpnl
            do jj=1,jppj
              zpj(ii,jj)=1.e-20
            enddo
          enddo
        else
        if (usefastj) then                                 !Begin fastj2 j-calc

! Initialize zpj (jpnl,jppj) array (hyl,02/19/03)
          do ii=1,jpnl
            do jj=1,jppj
              zpj(ii,jj)=0.D0
            enddo
          enddo

! Note fastj arrays are B->T (dbp,dbz,dbt,dbox,dbcolo3v) (hyl, 11/12/2002)
! Note svcolo3v is T->B and model "grd" arrays are B->T in raqms_global.
!  Because model "grd" arrays may be T->B in raqms_regional, use
!  indexing vector ltb (jaa 08/09/2007)

          do 12 ll = 1,nl
            lbt=nl-ltb(ll)+1
            dbp(ll) = dble(pgrd(i,j,lbt))
            dbz(ll) = dble(zgrd(i,j,lbt))
            dbt(ll) = dble(tgrd(i,j,lbt))
            dbox(ll)= dble(oxgrd(i,lbt))
            dbcolo3v(ll)= dble(svcolo3v(nl-ll+1))   !column ozone (molec cm-2)
!            if(tile.eq.1.and.iatchem.eq.139.and.jatchem.eq.101.and.doptdl(ll)/=0.0)then
!              write(6,*)'cloudphot',cloudphot
!            endif

            if (cloudphot) then
              !cloudy sky
cjaa          !Note that tau3d is B-->T (hyl, 04/12/03)
cjaa          doptd(ll) = dble ( tau3d(i,j,ll) )
cjaa          if ( doptd(ll) .lt. 0 .or. doptd(ll) .gt. 500 ) then
cjaa           write(6,*) 'doptd(ll),i,j,ll= ', doptd(ll),i,j,ll
cjaa           write(6,*) 'doptd: negative OR too large !!!'
cjaa           call flush(6)
cjaa           stop
cjaa          endif

cjaa Using modified ccm3 extinction optical depth calc
cjaa  1) liquid and ice are available separately, so use them
cjaa  2) ccm3 uses a factor of cloud^(3/2) to account for "random overlap"
cjaa     but this gives cot values MUCH too low (instead of just too low)
cjaa     in our current 6-hourly reinit cycle.  Instead use cloud^1
c
              lbtccm = nl-ll+1  ! ccm3 is 1->nl for T->B
!              doptdl(ll) = dble(tauxcl3d(i,j,lbtccm)*cloud3d(i,j,lbtccm))
!              doptdi(ll) = dble(tauxci3d(i,j,lbtccm)*cloud3d(i,j,lbtccm))
!              doptdl(ll) = dble(taucldfrc_liq(i,j,lbtccm))
!              doptdl(ll) = dble(taucldfrc(i,j,lbtccm))
!             lbt accounts for t-b or b-t
              doptdl(ll) = dble(taucldfrc(i,j,lbt))
!            if(tile.eq.1.and.iatchem.eq.139.and.jatchem.eq.101.and.doptdl(ll)/=0.0)then
!ji            if(tile.eq.1.and.jatchem.eq.101.and.doptdl(ll)/=0.0)then
!              write(6,*)iatchem,'ll',ll,'doptdl',doptdl(ll),'pres',dbp(ll),'lbt',lbt
!              call flush(6)
!            endif
!              if(doptdl(ll).ne.0.0)then
!                write(200+iam,*)i,j,'ll',ll,' lbt ',lbt,' lbtccm ',lbtccm,' pgrd ',pgrd(i,j,lbt)
!                write(200+iam,*)'zgrd',dbz(ll),'taucld',doptdl(ll)
!              endif
!              doptdi(ll) = dble(taucldfrc_ice(i,j,lbtccm))
              doptdi(ll) = 0.0
            
              alb = albgrd(i,j)     !RAQMS 2D
!              if(iam.eq.2)then
!              write(6,*)'i,j',i,j,'alb',alb,'doptdl',doptdl(ll),ll,doptdi(ll)
!              endif

cjaa Crude kludge, scale optical depths to partially correct for low
cjaa  bias resulting from 6-hourly reinitialization (Dec 08 2004)
c          doptdl(ll) = doptdl(ll)*4.0
c          doptdi(ll) = doptdi(ll)*4.0
cjaa 2005/06/14 Now using factor of 3.0 based on March comparison with
cjaa  Hongyu's MODIS monthly means
              doptdl(ll) = doptdl(ll)*3.0
              doptdi(ll) = doptdi(ll)*3.0

CCCMNN 2008/03/18  GET OPTICAL DEPTH DUE TO DIFFERENT AEROSOL TYPES.
CCCMNN             USE MODEL CALCULATED SULFATE,DUST,BC,OC, AND SEASALT
CCCMNN             ABSORPTION AND SCATTERING PARAMETERS FROM HARVARD
CCCMNN             (DUNCAN FAIRLIE).  EACH AEROSOL (SULFATE,BC,OC, AND
CCCMNN             SEASALT) AT DIFFERENT RELATIVE HUMIDITIES IS 
CCCMNN             CONSIDERED AS A SEPARATE CATEGORY OF AEROSOL.
CCCMNN             BASED ON THE LOCAL RELATIVE HUMIDITY, ONLY ONE
CCCMNN             CATEGORY IN EACH TYPE WILL HAVE NONZERO OPTICAL
CCCMNN             DEPTH).
CCCC     Mass concentration of aerosols
              dens = 7.2431122e+18*pgrd(i,j,lbt)/tgrd(i,j,lbt)
              dens_m=dens*28.97/6.023e20
!              write(100+iam,*)' at bc i',i,'lbt',lbt
!              call flush(100+iam)
#ifdef DOINT
      
              so4aergrdr8=so4aergrd(i,lbt)
              bc1grdr8=bc1grd(i,lbt)
              bc2grdr8=bc2grd(i,lbt)
              oc1grdr8=oc1grd(i,lbt)
              oc2grdr8=oc2grd(i,lbt)
#ifdef AEROSOL_SEASALT
              ss1grdr8=ss1grd(i,lbt)
              ss2grdr8=ss2grd(i,lbt)
              ss3grdr8=ss3grd(i,lbt)
              ss4grdr8=ss4grd(i,lbt)
#endif
#ifdef AEROSOL_SEASALT_AOD
              ss1grdr8=ss1grd(i,lbt)
              ss2grdr8=ss2grd(i,lbt)
              ss3grdr8=ss3grd(i,lbt)
              ss4grdr8=ss4grd(i,lbt)
                ss5grdr8=ss5grd(i,lbt)
#endif
              so4aergrdr8=max(so4aergrdr8,zero8)
              bc1grdr8=max(bc1grdr8,zero8)
              bc2grdr8=max(bc2grdr8,zero8)
              oc1grdr8=max(oc1grdr8,zero8)
              oc2grdr8=max(oc2grdr8,zero8)
#ifdef AEROSOL_SEASALT
              ss1grdr8=max(ss1grdr8,zero8)
              ss2grdr8=max(ss2grdr8,zero8)
              ss3grdr8=max(ss3grdr8,zero8)
              ss4grdr8=max(ss4grdr8,zero8)
#endif
#ifdef AEROSOL_SEASALT_AOD
              ss1grdr8=max(ss1grdr8,zero8)
              ss2grdr8=max(ss2grdr8,zero8)
              ss3grdr8=max(ss3grdr8,zero8)
              ss4grdr8=max(ss4grdr8,zero8)
              ss5grdr8=max(ss5grdr8,zero8)
#endif
! ajl 10/15/2021 here ppv for so4,bc,oc1 kg/kg for ss and dust findajlz
!   tr3d_inout4(:,:,:,p_bc1)=gsdchem_tr3d_inout8(:,:,:,p_bc1)*1.e-9*airmw/tcmw(idbc1)
!  tr3d_inout4(:,:,:,p_seas1)=gsdchem_tr3d_inout8(:,:,:,p_seas1)*1.e-9
              if(aerosol_ugpkg)then
!                wwcon(1)=so4aergrdr8*dens*tcmw(idso4)/6.023e20*1.e-9
                wwcon(1)=so4aergrdr8*dens*airmw/6.023e20*1.e-9
                dwcon(1)=bc1grdr8*dens*airmw/6.023e20*1.e-9
!                if(isnan(dwcon(1)))then
!                   write(6,*)'bc1grd',i,lbt,bc1grd(i,lbt),'dens',dens,'airmw',airmw
!                endif
                wwcon(2)=bc2grdr8*dens*airmw/6.023e20*1.e-9
                dwcon(2)=oc1grdr8*dens*airmw/6.023e20*1.e-9
!                if(isnan(dwcon(2)))then
!                   write(6,*)'oc1grd',i,lbt,oc1grd(i,lbt),'dens',dens,'airmw',airmw
!                endif
                wwcon(3)=oc2grdr8*dens*airmw/6.023e20*1.e-9
                if(use_fv3_aero)then
                  wwcon(4)=(ss1grdr8+ss2grdr8)*dens_m*1.e-9
                  wwcon(5)=(ss3grdr8+ss4grdr8+ss5grdr8)*dens_m*1.e-9
!                  if(isnan(wwcon(5)))then
!                     write(6,*)'ss3',ss3grdr8,ss4grdr8,ss5grdr8
!                   endif

                else
                  wwcon(4)=ss1grdr8*dens_m*1.e-9
                  wwcon(5)=(ss2grdr8+ss3grdr8+ss4grdr8+ss5grdr8)*dens_m*1.e-9
                endif
              else
                wwcon(1)=so4aergrdr8*dens*tcmw(idso4)/6.023e20
                dwcon(1)=bc1grdr8*dens*tcmw(idbc1)/6.023e20
                wwcon(2)=bc2grdr8*dens*tcmw(idbc2)/6.023e20
                dwcon(2)=oc1grdr8*dens*tcmw(idoc1)/6.023e20
                wwcon(3)=oc2grdr8*dens*tcmw(idoc2)/6.023e20
                if(use_fv3_aero)then
                  wwcon(4)=(ss1grdr8+ss2grdr8)*dens_m
                  wwcon(5)=(ss3grdr8+ss4grdr8+ss5grdr8)*dens_m
                else
                  wwcon(4)=ss1grdr8*dens_m
                  wwcon(5)=(ss2grdr8+ss3grdr8+ss4grdr8+ss5grdr8)*dens_m
                endif
              endif
#else
              if(aerosol_ugpkg)then
!                wwcon(1)=max(so4aergrd(i,lbt),0.)*dens*tcmw(idso4)/6.023e20*1.e-9
                wwcon(1)=max(so4aergrd(i,lbt),0.)*dens*airmw/6.023e20*1.e-9
                dwcon(1)=max(bc1grd(i,lbt),0.)*dens*airmw/6.023e20*1.e-9
!                if(isnan(dwcon(1)))then
!                   write(6,*)'bc1grd',i,lbt,bc1grd(i,lbt),'dens',dens,'airmw',airmw
!                endif
                wwcon(2)=max(bc2grd(i,lbt),0.)*dens*airmw/6.023e20*1.e-9
                dwcon(2)=max(oc1grd(i,lbt),0.)*dens*airmw/6.023e20*1.e-9
!                if(isnan(dwcon(2)))then
!                   write(6,*)'oc1grd',i,lbt,oc1grd(i,lbt),'dens',dens,'airmw',airmw
!                endif
                wwcon(3)=max(oc2grd(i,lbt),0.)*dens*airmw/6.023e20*1.e-9
                if(use_fv3_aero)then
                  wwcon(4)=(max(ss1grd(i,lbt),0.)+max(ss2grd(i,lbt,0.))*dens_m*1.e-9
                  wwcon(5)=(max(ss5grd(i,lbt),0.)+max(ss3grd(i,lbt),0.)+
     &            max(ss4grd(i,lbt),0.)+max(ss5grd(i,lbt),0.))*dens_m*1.e-9
                else
                  wwcon(4)=max(ss1grd(i,lbt),0.)*dens_m*1.e-9
                  wwcon(5)=(max(ss2grd(i,lbt),0.)+max(ss3grd(i,lbt),0.)+
     &            max(ss4grd(i,lbt),0.)+max(ss5grd(i,lbt),0.))*dens_m*1.e-9
                endif
              else
                wwcon(1)=max(so4aergrd(i,lbt),0.)*dens*tcmw(idso4)/6.023e20
                dwcon(1)=max(bc1grd(i,lbt),0.)*dens*tcmw(idbc1)/6.023e20
                wwcon(2)=max(bc2grd(i,lbt),0.)*dens*tcmw(idbc2)/6.023e20
                dwcon(2)=max(oc1grd(i,lbt),0.)*dens*tcmw(idoc1)/6.023e20
                wwcon(3)=max(oc2grd(i,lbt),0.)*dens*tcmw(idoc2)/6.023e20
                if(use_fv3_aero)then
                  wwcon(4)=(max(ss1grd(i,lbt),0.)+max(ss2grd(i,lbt),0.)*dens_m
                  wwcon(5)=(max(ss5grd(i,lbt),0.)+max(ss3grd(i,lbt),0.)+
     &            max(ss4grd(i,lbt),0.)+max(ss5grd(i,lbt),0.))*dens_m
                else
                  wwcon(4)=max(ss1grd(i,lbt),0.)*dens_m
                  wwcon(5)=(max(ss2grd(i,lbt),0.)+max(ss3grd(i,lbt),0.)+
     &            max(ss4grd(i,lbt),0.)+max(ss5grd(i,lbt),0.))*dens_m
                endif
              endif
#endif
!              write(6,*)'wwcon',wwcon,'dwcon',dwcon
!              call flush(6)
!            if(tile.eq.1.and.iatchem.eq.139.and.jatchem.eq.101)then
!              write(6,*)'ll',ll,'wwcon',wwcon,'lbt',lbt,'dens',dens,'p',pgrd(i,j,lbt)
!            endif
!        if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!          write(200+iam,*)'wwcon ll','i',i,'lbt',lbt,' tile ',tile
!          call flush(200+iam)
!        endif

              if (ll.eq.nl) then
                dz = 0.
              else
                if(nl-ll<1)then
                  write(6,*)'dz error nl-ll',nl-ll,nl,ll
                  call flush(6)
                endif
                dz = (zgeoh(nl-ll)-zgeoh(nl-ll+1)) * 1.e5   !cm
!                if(tile.eq.1.and.iatchem.eq.139.and.jatchem.eq.101)then
!                  write(6,*)'dz',dz,'zgeoh',nl-ll,zgeoh(nl-ll:nl-ll+1)
!                endif
              endif
C             GET RELATIVE HUMIDITY
              abshum = max(qgrd(i,j,lbt)*XMW_DRYAIR/XMW_H2O,3.0e-6)*dens
              trm1 = 17.2693882d0*(dbt(ll)-273.16d0)/(dbt(ll)-35.86d0)
              rel_hum=abshum * dbt(ll)/(6.1078e3*exp(trm1)/1.38054e-16)
c             RH bins are 0.0 to 0.9 so restrict to this range
              rel_hum = max(rel_hum,0.)
              rel_hum = min(rel_hum,0.9)
!              if(iam.eq.2)then
!                  write(6,*)'rel_hum',rel_hum,ll
!              endif
      
CCCCCCCCCCCC LOOP OVER AEROSOL TYPES
CC
              do ian = 1,5
CCCCCCCCCCCC INITIALIZE ODSCLF (OPTICAL DEPTH SCALE FACTOR)
                do jrh = 1,5
                  odsclf(jrh) = 0.0
                enddo
CCC             Humidity bin loop
                do jrh = 1,5
                  rwet(jrh) = raa(4,indx(ian)+jrh-1)
                  frac = (1.0 - (rwet(1)/rwet(jrh))**3)
                  qwet(jrh) = qaa(4,indx(ian)+jrh-1)*frac + qaa(4,indx(ian))*
     1               (1.0 - frac)
                end do
ccc
                if (rel_hum .le. rh(2)) then
                  irh = 1
                else if (rel_hum .le. rh(3)) then
                  irh = 2
                else if (rel_hum .le. rh(4)) then
                  irh = 3
                else if (rel_hum .le. rh(5)) then
                  irh = 4
                else
                  irh = 5
                endif

                if (irh .eq. 5) then
                  qsclf = qwet(irh)/qwet(1)
                  reff = rwet(irh)
                else
                  frac = (rel_hum - rh(irh))/(rh(irh+1)-rh(irh))
                  qsclf = (qwet(irh+1)*frac+qwet(irh)*(1.0-frac))/qwet(1)
                  reff = rwet(irh+1)*frac + rwet(irh)*(1.0-frac)
                endif
CCCC
                rsclf = reff/rwet(1)
                odsclf(irh) = qsclf*rsclf*rsclf
         
CCCC            calculate optical depth 
CCCC            OD  = (0.75*dz*M*Q)/(aerosol_density*r)

                do jrh = 1,5
                  irhn = (ian-1)*5 + jrh
                  aer_opd(irhn,ll) = odsclf(jrh)*0.75D0*wwcon(ian)*dz*
     1           qaa(4,indx(ian))/(aer_den(ian)*raa(4,indx(ian))*1.0D-4)
!                  if(isnan(wwcon(ian)).or.isnan(aer_opd(irhn,ll)))then
!                   write(6,*)'wwcon',ian,wwcon(ian),irhn,ll
!                  endif
                end do
CCCC            Aerosol volume and surface area density
CCCC
                v_arsl(ian,ll)=wwcon(ian)/aer_den(ian)
                s_arsl(ian,ll) = 3.0*v_arsl(ian,ll)*1.0e4*rsclf**3/reff
                eff_r(ian,ll) = reff*1.0e-4

              end do            !  end of loop over aerosol types
CCCCCCC*****  insert code for surface area density
CCCCCCCC  INCLUDE HYDROPHOBIC AEROSOLS  (BC AND OC; IAN=2,3)
CC
              do ian = 2,3
                irhn = (ian-1)*5 + 1
                aer_opd(irhn,ll) = aer_opd(irhn,ll) + 0.75D0 * dz *
     1          dwcon(ian-1) * qaa(4,indx(ian))/
     1          (aer_den(ian)*raa(4,indx(ian))*1.0D-4)
!                if(isnan(aer_opd(irhn,ll)))then
!                  write(6,*)'dwcon',ian-1,dwcon(ian-1),'qaa',qaa(4,indx(ian)),'den',aer_den(ian),'raa',raa(4,indx(ian))
!                 endif
              end do

CCCCCCCCCCC   DUST OPTICAL DEPTH

#ifdef DOINT
              du1r8=du1grd(i,lbt)
              du2r8=du2grd(i,lbt)
              du3r8=du3grd(i,lbt)
              du4r8=du4grd(i,lbt)
              du5r8=du5grd(i,lbt)
              wcdu(1)=max(du1r8,zero8)*dens_m
              wcdu(2)=max(du2r8,zero8)*dens_m
              wcdu(3)=max(du3r8,zero8)*dens_m
              wcdu(4)=max(du4r8,zero8)*dens_m
              wcdu(5)=max(du5r8,zero8)*dens_m
#else
              wcdu(1)=max(du1grd(i,lbt),0.)*dens_m
              wcdu(2)=max(du2grd(i,lbt),0.)*dens_m
              wcdu(3)=max(du3grd(i,lbt),0.)*dens_m
              wcdu(4)=max(du4grd(i,lbt),0.)*dens_m
              wcdu(5)=max(du5grd(i,lbt),0.)*dens_m
#endif
              if(aerosol_ugpkg)then
                wcdu(1:5)=wcdu(1:5)*1.e-9
              endif

c
              do ian = 1,ndust             !  loop over size bins
                indst = 17+ian
                if (ian .eq. 5) indst=21    ! use 4 micron data for the 5th bin
                dust_opd(ian,ll) = 0.75D0*dz*wcdu(ian)*qaa(4,indst)/
     1                      (dust_den(ian)*raa(4,indst)*1.0D-4)
              end do


            else           !end of "if(cloudphot)"
              !clear sky
              doptdl(ll) = 0.D0
              doptdi(ll) = 0.D0
              alb = 0.1     !use mean (tropospheric) value of 10% as in lookup table
              do idum=1,naer
                aer_opd(idum,ll) = 0.D0
              enddo
              do idum=1,ndust
                dust_opd(idum,ll) = 0.D0
              enddo
            endif

   12     continue
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!            write(200+iam,*)'after 12 loop'
!            call flush(200+iam)
!          endif
c         real*8 doptdl(nl),doptdi(nl),aer_opd(25,nl),dust_opd(5,nl)
 1010     format(7e11.4)

!         For efficiency, only calculate in daytime. See above for the
!         definition of daytime. (hyl,02/19/03)
#ifdef DIAGFV3
          if(tile.eq.1.and.ibeg+i-1.eq.1.and.j.eq.1)then
            write(6,*)'zeta',zeta
            call flush(6)
          endif
#endif
!          write(6,*)'zeta',zeta
!          call flush(6)
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!            write(200+iam,*)'azeta',zeta
!            call flush(200+iam)
!          endif
          if (zeta .ge. -0.068) then

!           Call subroutine PHOTOJ to compute J-values for the entire column.
!           note now lbot is passed into Fast-J2
#ifdef DOINT
            gmt8=gmt
!            if(iam.eq.0)then
!              write(6,*)'gmt',gmt,'gmt8',gmt8
!            endif
            spgrd8=spgrd(i,j)
            zsurf8=zsurf(i,j)
#ifdef DIAGFV3
            if(tile.eq.1.and.ibeg+i-1.eq.1.and.j.eq.1)then
!              write(200+iam,*)'lbot',lbot,'rlat',rlat,'gmt8',gmt8,'julday',julday,'imon',imon
              write(6,*)'gmt',gmt8,'julday',julday,'imon',imon,'call phojoj'
              write(6,*)'sza_deg',sza_deg,kind(sza_deg),'lbot',lbot,'rlat',rlat
              write(6,*)'sp',spgrd8,'zsurf',zsurf8,'alb',alb
              call flush(6)
            endif
#endif
#ifdef DIAGCHEM
        if(tile.eq.tileprn.and.ibeg+i-1.eq.iprnin.and.j.eq.jprnin)then
              write(300+iam,*)'tile',tile
              write(300+iam,*)'iprn',iprn,jprn,'iatchem',iatchem,jatchem,'iamprn',iamprn
              write(300+iam,*)'gmt8',gmt8,'julday',julday,'imon',imon
              write(300+iam,*)'spgrd8',spgrd8,'zsurf8',zsurf8
              write(300+iam,*)'dbp',dbp
              write(300+iam,*)'dbz',dbz
              write(300+iam,*)'dbt',dbt
              write(300+iam,*)'dbox',dbox
              write(300+iam,*)'colo',dbcolo3v
              write(300+iam,*)'doptdl',doptdl
              write(300+iam,*)'doptdi',doptdi
              write(300+iam,*)'dust',dust_opd
              call flush(300+iam)
            endif
#endif
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!              write(200+iam,*)'call photoj ',i,j,lbot
!              !call flush(200+iam)
!            !endif
!            if(iatchem.eq.139.and.jatchem.eq.101.and.tile.eq.1)then
!              write(6,*)'call photoj lbot',lbot,'rlat',rlat,'gmt',gmt8,'doptdl',nl-lbot+1,doptdl(nl-lbot+1)
!              call flush(6)
!            endif
            call photoj ( i,  j, lbot, rlat, gmt8, julday, imon,
     &             dble(sza_deg), spgrd8, zsurf8,
     &             dbp(nl-lbot+1), dbz(nl-lbot+1), dbt(nl-lbot+1),
     &             dbox(nl-lbot+1), dbcolo3v(nl-lbot+1),
cjaa     &             doptd(nl-lbot+1),dble(alb), zpj )
     &             doptdl(nl-lbot+1),doptdi(nl-lbot+1),
CCCMNN     &             doptdbc(nl-lbot+1),dble(alb), zpj )
     &             aer_opd,dust_opd,dble(alb),zpj)
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!              write(200+iam,*)'did call photoj ',i,j,lbot
!              call flush(200+iam)
!            !endif
#ifdef DIAGCHEM
        if(tile.eq.tileprn.and.ibeg+i-1.eq.iprnin.and.j.eq.jprnin)then
               write(300+iam,*)'zpj',zpj
               write(300+iam,*)'alb',alb
               write(300+iam,*)'sza_deg',sza_deg
               write(300+iam,*)'aer_opd',aer_opd
               write(300+iam,*)'lbot',lbot
               write(300+iam,*)'rlat',rlat
             endif
!            if(i.eq.1.and.j.eq.begj)then
!             write(200+iam,*)'zpj',zpj
!             call flush(200+iam)
!            endif
#endif
#else
            call photoj ( i,  j, lbot, dble(rlat), dble(gmt), julday, imon,
     &             dble(sza_deg), dble(spgrd(i,j)), dble(zsurf(i,j)),
     &             dbp(nl-lbot+1), dbz(nl-lbot+1), dbt(nl-lbot+1),
     &             dbox(nl-lbot+1), dbcolo3v(nl-lbot+1),
cjaa     &             doptd(nl-lbot+1),dble(alb), zpj )
     &             doptdl(nl-lbot+1),doptdi(nl-lbot+1),
CCCMNN     &             doptdbc(nl-lbot+1),dble(alb), zpj )
     &             aer_opd,dust_opd,dble(alb),zpj)
#endif

c      if (i.eq.1 .and. j.eq.46) then
c       print*,'After call to fastj2'
c       print*,'usefastj,cloudphot',usefastj,cloudphot
c       print*,i,j,lbot,rlat,gmt,julday,
c     &             sza_deg,spgrd(i,j),zsurf(i,j),alb
c       print*,'p, z, T, O3, OD_l, OD_i, jo1d'
c       do ll=1,nl
c        write(6,1010) dbp(ll),dbz(ll),dbt(ll),dbox(ll),
c     &  doptdl(ll),doptdi(ll),zpj(ll,3)
c       end do
c       call flush(6)
c      endif
c 1010 format(7e11.4)

          endif   !daytime

        endif                                              !End fastj2 j-calc
        endif ! skipphotoj
! ------------------------------------------------------------------------

c  construct lightning NOx distribution, B->T for consistency w/dataset
ccccmn0299   reduce the lightning source (the total value from the
cccc         datafile (t-32) is 12.39 tg/year) by a factor (0.242)
cccc         to make the total rate = 3 tg/year
!        write(100+iam,*)'do light i',i,imon
!        call flush(100+iam)
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!            write(200+iam,*)' lightning'
!            call flush(200+iam)
!          endif
        column=slghtng(i,j,imon)*0.242
        lnoxtop = ltrop
c
!        if(i.eq.nc) then
!          dlon = abs(xgrid(i,j)-xgrid(i-1,j))*dtor
!        else
!          dlon = abs(xgrid(i+1,j)-xgrid(i,j))*dtor
!        endif
        !if(j.eq.nr) then
!          dlat = abs(ygrid(i,j)-ygrid(i,j-1))*dtor
!        else
!          dlat = abs(ygrid(i,j+1)-ygrid(i,j))*dtor
!        endif
!        dlon=1.
!        dlat=1. ! fix later
        !areacm2 = recm*recm*dlon*dlat*cos(rlat*dtor)
        areacm2=areagrd(i,j)*1.e4
   
        area225=recm*recm*2.*2.5*dtor*dtor*cos(30.*dtor)
!        if(iam.eq.0)then
!          write(6,*)'areacm2',areacm2,area225
!        endif
c
#ifdef CALCLIGHTNINGNOX
c
c  Replace climatological column with instantaneous flashrate-based calc.
c  Top of profile is convective cloudtop height instead of tropopause
c
        column = 0.
!        !lnoxtop = lconvcld(i,j)           ! lconvcld index is T->B
!       fix for fv3
        dz = ccthk(i,j)
        pcolnox(i,j)=1.e32
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!            write(200+iam,*)'a dz',dz
!            call flush(200+iam)
!          endif
        sorclnox=0.0
        if (dz .gt. 0.) then
!         in raqms ccm3 1 was at top of conv so need to flip 
!         in fv3 phys 1 is at bottom so don't need to flip 
!          lnoxtop = ltb(lconvcld(i,j))           ! lconvcld index is T->B is b-t in fv3
          lnoxtop = lconvcld(i,j)           ! lconvcld index is T->B is b-t in fv3
          if(dopr94)then
            rr=dlon*dlat*cos(rlat*dtor)
            coefpr94=.97241*exp(.048203*rr)

          endif
          if (dz .ge. 14.) then
            fcg = 0.02
          else
            fcg = 1./(ap_a*dz**4+ap_b*dz**3+ap_c*dz**2+ap_d*dz+ap_e+1)
          endif
!         in raqms ccm3 1 was at top of conv so need to flip 
!         in fv3 phys 1 is at bottom so don't need to flip 
!          l = ltb(lconvcld(i,j))          ! lconvcld index is T->B
          l = lconvcld(i,j)          ! lconvcld index is T->B
          if (pgrd(i,j,l) .le. 440.) then
            pcolnox(i,j)=pgrd(i,j,lnoxtop)
          endif
          if (pgrd(i,j,l) .le. 440.) then
            ccldht = zgrd(i,j,l)/1000.     !km
            if (zlwigrd(i,j) .ge. 1.5 .and. zlwigrd(i,j) .le. 2.5) then
              ! Land
              flash = 3.44e-5*ccldht**4.9
            else
              ! Ocean (include ice)
              flash = 6.40e-4*ccldht**1.73
            endif
            !
            ! Scale flashrate according to Pickering [pers.comm. 2005]
            !
            flash = flash*0.648                   !global
            if (zlwigrd(i,j) .lt. 1.5) then
              flash = flash*5.145                 !oceanic
            else if (zlwigrd(i,j) .lt. 2.5) then
              if (rlat .ge. 0.) then
                flash = flash*1.044               !NH continental
              else
                flash = flash*0.485               !SH continental
              endif
            endif
            !
            ! Have flashrate/min/gridbox, calculate NOx molecules/cm^2/sec
            !
            column = flash*(cgnox*fcg+ccnox*(1.-fcg))
           
#ifdef DIAGSNO2
      if(i.eq.99.and.j.eq.95)then
        write(6,*)'column',column,'flash',flash,'cgnox',cgnox,'fcg',fcg
        write(6,*)'ccnox',ccnox,'zlwigrd',zlwigrd(i,j),'p',pgrd
     *  (i,j,l),'l',l,'zgrd',zgrd(i,j,l),'dz',dz
      endif
#endif
            if(doareafix)then
              column=column*areacm2/area225
            endif
            if(dopr94)then
!            area5kmcm2=(5.*1.e5)**2
            columnold = column/areacm2/60.
            column = column/areacm2/60.*coefpr94
!            print *,'columnold',columnold,'new',column,
!     *      'ratio',column/columnold
            else
            column = column/areacm2/60.
            endif
          endif
        endif
#endif
c
c       Determine whether midlatitude, tropical marine, or tropical land
c       land/water/ice flag: 1=ocean, 2=land, 3=ice
        iprof=1
        if (rlat.ge. -23.5 .and. rlat.le. 23.5) then
          if (zlwigrd(i,j) .lt. 1.5) then
            iprof=2
          else
            iprof=3
          endif
        endif
cc       test=0
        zhu=16.
        li=16
!        if(j.eq.46)then
!          write(6,*)'lbot',lbot,'zgeohlnoxttop',lnoxtop-1,
!     *    zgeoh(lnoxtop-1),'zhlbot',zgeoh(lbot)
!        endif
!        write(200+iam,*)'lbot',lbot,' ltb(1) ',ltb(1),'lnoxtop',lnoxtop
!        call flush(200+iam)
!        write(200+iam,*)'zgeohlbot',zgeoh(lbot),'lnoxtop-1',zgeoh(lnoxtop-1)
!            call flush(200+iam)
        do ip=1,lbot
          l=ltb(ip)
          colno2a=0.0
          colno2b=0.0
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!            write(200+iam,*)'ip',ip,'l',l
!            call flush(200+iam)
!          endif
          sorclnox(ip)=0.
!          write(200+iam,*)'ip',ip,'zgeoh',zgeoh(ip)
          if(column.gt.0. .and. ip.ge.lnoxtop .and. ip.le.lbot) then
c
c           Get normalized Z at cell edges corresponding to Pickering 0-16km grid
            zhl=16.*(zgeoh(ip)-zgeoh(lbot))/(zgeoh(lnoxtop-1)-zgeoh(lbot))
!            write(200+iam,*)'zhl',zhl,'ip',ip
!            call flush(200+iam)
            if(zhl<0.0)then
               write(6,*)'ip',ip,'zgeoh',zgeoh(ip),'lbot',zgeoh(lbot),
     *           'lnoxtop-1',lnoxtop-1,zgeoh(lnoxtop-1)
               call flush(6)
            endif
c
c           Sum/interpolate Pickering intervals to get distribution on our grid
            lip=int(zhl)
            wtav=0.
c           print*,ip,zhu,zhl,li,lip
            z1=zhu
            ll=li
            do while (ll .gt. lip+1)
              z2=float(ll)-1.
!              if(ll>6)then
!                 write(200+iam,*)'ll',ll,'lip',lip,'zhl',zhl,'ip',ip,'l',l,'zgeoh',zgeoh(ip)
!                 !write(200+iam,*)'lbot',lbot,'zgeoh',zgeoh(lbot),'lnoxtop',lnoxtop,zgeoh(lnoxtop-1)
!                 call flush(200+iam)
!              endif
              wtav=wtav+slnox_profile(ll,iprof)*(z1-z2)
#ifdef DIAGSNO2
              if(i.eq.99.and.j.eq.95.and.(ip.eq.18.or.ip.eq.19))then
                write(6,*)'wtava',wtav,'slnox_profile',slnox_profile(ll,iprof)
                 write(6,*)'ll',ll,'iprof',iprof,'z1',z1,z2
              endif
#endif
              z1=z2
              ll=ll-1
            enddo
            if(ll<1)then
              write(6,*)'i',i,j,'zhl',zhl,'ip',ip,zgeoh(ip),'lbot',lbot,
     *'lnoxtop',lnoxtop,zgeoh(lnoxtop-1),'l',l
              write(6,*)'zgeoh',zgeoh
              call flush(6)
            endif
            z2=zhl
            wtav=wtav+slnox_profile(ll,iprof)*(z1-z2)
c           Divide by layer thickness to get rate, convert km-1 to cm-1,
c           and divide by 100 to convert from percent
c           zgeoh is in km
            sorclnox(ip)=wtav*column*1.e-7/(zgeoh(ip-1)-zgeoh(ip))
!            if(iam.eq.iamprn.and.abs(i-iprn)<3 .and. abs(j-jprn)<3)then
!              write(6,*)'ip',ip,' sorclnox ',sorclnox(ip),' wtav ',wtav
!              write(6,*)'column',column,' zgeoh ',ip-1,zgeoh(ip-1:ip),'dz ',zgeoh(ip-1)-zgeoh(ip)
!              write(6,*)'ll ',ll,' iprof ',iprof,' z1 ',z1,z2,' z1-z2 ',z1-z2
            !endif
#ifdef DIAGLIGHT
            if(tile.eq.2)then
              write(300+iam,*)'i',i,j,'lat',latgrd(i,j),longrd(i,j)
              write(300+iam,*)'ip',ip,' sorclnox ',sorclnox(ip),' wtav ',wtav
              write(300+iam,*)'column',column,' zgeoh ',ip-1,zgeoh(ip-1:ip),'dz ',zgeoh(ip-1)-zgeoh(ip)
              write(300+iam,*)'ll ',ll,' iprof ',iprof,' z1 ',z1,z2,' z1-z2 ',z1-z2
              write(300+iam,*)'wtav/dz ',wtav/(zgeoh(ip-1)-zgeoh(ip))
              call flush(300+iam)
            endif
#endif
!              write(300+iam,*)i,j,'wtav',ip,wtav,'col',column,'zgeoh',zgeo(ip-1:ip),
!     *   'z1',z1,z2,'ll',ll,'l',l,'p',pgrd(i,j,l)
!!            endif
            zhu=zhl
            li=ll
cc            test=test+wtav
          endif
        enddo

c       land/water/ice flag: RAQMS, 1=ocean, 2=land, 3=ice; gocart (ilwi) 0=W, 1=L, 2=I
        ilwi = 0
        if (zlwigrd(i,j) .ge. 1.5) ilwi = ilwi+1
        if (zlwigrd(i,j) .ge. 2.5) ilwi = ilwi+1

#ifdef DO_AEROSOL
c
c Implementation of GOCART:  (JAA 2007)
c  Deposition velocities and emissions are needed prior to do_chem so
c   that additional gas phase source terms can be calculated (eg, Br from SS)
c  Emissions are applied in do_chem to make use of PBL logic
c  All other tendencies (settling, drydep, chemical transformation) applied
c   after do_chem (subroutine raqms_gocart_driver)
c  Wetdep is applied in physics as for gas phase species
c
c Here: interface with gocart to:
c  1) calculate deposition velocities for each aerosol related species
c  2) calculate source/emission as needed for particular species
c
c delz_sfc in meters, zgeoh is T->B (km)
c pblht in meters
c ts: T in lowest layer (K)
c land/water/ice flag: RAQMS, 1=ocean, 2=land, 3=ice; gocart, 0=W, 1=L, 2=I
c ps: p in lowest layer (hPa) (mb)
c zeta = cos(sza) is calculated above (exists in common block)
        if(lbot-1<1)then 
          write(6,*)'lbot error 2 lbot',lbot
          call flush(6)
        endif
        dz = 1.e3*(zgeoh(lbot-1)-zgeoh(lbot))
        l = ltb(lbot)
        zeta1 = zeta
#if (defined DOINT || defined DOINTCHEMB || defined DOINTCHEMM)
        t8=tgrd(i,j,l)
        a8=albgrd(i,j)
        p8=pgrd(i,j,l)
        call calc_depvel(i,j,dz,pblht(i,j),t8,
     &  a8,ilwi,p8,zeta1)
#else
#ifdef COMMR8
        call calc_depvel(i,j,dz,pblht(i,j),tgrd(i,j,l),
     &  dble(albgrd(i,j)),ilwi,pgrd(i,j,l),zeta1)
#else
        call calc_depvel(i,j,dz,pblht(i,j),dble(tgrd(i,j,l)),
     &  dble(albgrd(i,j)),ilwi,dble(pgrd(i,j,l)),zeta1)
#endif
#endif
c
c       dxy is cell area in square meters, areacm2 is calculated above
        dxy = areacm2*1.e-4
        sorcaero(:) = 0.
#if (defined DOINT || defined DOINTCHEMC || defined DOINTCHEMM)
        u8=ugrd(i,j,l)
        v8=vgrd(i,j,l)
        tskin8=tskin(i,j)
        tll8=tgrd(i,j,l)
        p8=pgrd(i,j,l)
c       find need to change ajl
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!        write(200+iam,*)'get aersource'
!         call flush(200+iam)
!         endif
        call get_aerosource(i,j,ndtaero,ilwi,dxy,u8,
     &  v8,tskin8,tll8,p8,imon,sorcaero)
#else
#ifdef COMMR8
        call get_aerosource(i,j,ndtaero,ilwi,dxy,dble(ugrd(i,j,l)),
     &  dble(vgrd(i,j,l)),dble(tskin(i,j)),tgrd(i,j,l),
     &  pgrd(i,j,l),imon,sorcaero)
#else
        call get_aerosource(i,j,ndtaero,ilwi,dxy,dble(ugrd(i,j,l)),
     &  dble(vgrd(i,j,l)),dble(tskin(i,j)),dble(tgrd(i,j,l)),
     &  dble(pgrd(i,j,l)),imon,sorcaero)
#endif
#endif
#endif
!        write(200+iam,*)'1917',lbot
!        call flush(200+iam)
        call get_aerosource_ss(i,j,ilwi,ugrd(i,j,1),vgrd(i,j,1))

        zratio = 1.

        sza_deg = szagrd(i,j)           ! solar zenith angle (deg)
#ifdef DIAGFV3
        if(ibeg+i-1.eq.1.and.j.eq.1.and.tile.eq.1)then
           write(6,*)'szagrd',i,j,szagrd(i,j)
        endif
#endif
        alb = albgrd(i,j)               ! surface albedo
c  start of altitude loop
        do 170 ip=1,lbot
          l=ltb(ip)
!         in fv3 ip=1 is top of atm and l=nl for ip=1
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!        write(200+iam,*)'ip',ip,' 170','lbot',lbot,' l ',l
!        call flush(200+iam)
!          endif

          if (ip .ge. 3) then
            dzp=zgeoh(ip-2)-zgeoh(ip-1)
            dz =zgeoh(ip-1)-zgeoh(ip)
            zratio=dzp/dz
          endif
          if (sednat_in .gt. 0.) sednat_in = sednat_in*zratio
          if (sedice_in .gt. 0.) sedice_in = sedice_in*zratio

          press = pgrd(i,j,l)             ! pressure (mb)
!          if(ccthk(i,j)>0.0)then
!          write(200+iam,*)'ip',ip,' l ',l,' pres ',press,' zgeoh ',zgeoh(ip),' t ',tgrd(i,j,l)
!          write(200+iam,*)'lconvcld',lconvcld(i,j)
!          call flush(200+iam)
!          endif
!          write(200+iam,*)'pgrs',i,j,l,pgrd(i,j,l),'lbound',lbound(pgrd),ubound(pgrd),maxval(pgrd(:,j,:))
!          call flush(200+iam)
          temp = tgrd(i,j,l)              ! temperature (K)
          theta = thgrd(i,j,l)            ! potential temperature (K)
cjjjcolo3        colo3_vert = colo3grd(i,j,l)    ! overhead O3 column (DU)
cjjj        xsad = surfgrd(i,j,l)           ! SAD (micron^2/cm^3)
c convert overhead ozone from dobsons to (cm)-2
cjjjcolo3        colo3_vert = colo3_vert*dob_to_den
          sza = sza_deg*dtor
          zeta = cos(sza)
c         convert from volume mixing ratio to number density
          dens = 7.2431122e+18*press/temp
          pur = press
          ox = max(oxgrd(i,l)*dens,zer_min)
!          if(i.eq.1.and.j.eq.1)then
!             write(6,*)'ox l',l,ox,'oxgrd',oxgrd(i,l)
!             call flush(6)
!          endif
          xnoy = max(xnoygrd(i,l)*dens,zer_min)
          hno3 = max(hno3tgrd(i,l)*dens,zer_min)
          cly = max(clygrd(i,l)*dens,zer_min)
          xn2o5 = max(xn2o5grd(i,l)*dens,zer_min)
          h2o2 = max(h2o2grd(i,l)*dens,zer_min)
          hcl = max(hclgrd(i,l)*dens,zer_min)
          if(isnan(hcl))then
            write(6,*)i,j,l,'hclgrd',hclgrd(i,l),'dens',dens
            call flush(6)
          endif
          clno3 = max(clno3grd(i,l)*dens,zer_min)
          oclo = max(oclogrd(i,l)*dens,zer_min)
          xn2o = max(xn2ogrd(i,l)*dens,zer_min)
!          if(i.eq.1.and.j.eq.1)then
!             write(6,*)'xn2o l',l,xn2o,'xn2ogrd',xn2ogrd(i,l)
!             call flush(6)
!          endif
          f11 = max(f11grd(i,l)*dens,zer_min)
          f12 = max(f12grd(i,l)*dens,zer_min)
          ccl4 = max(ccl4grd(i,l)*dens,zer_min)
          ch3cl = max(ch3clgrd(i,l)*dens,zer_min)
          xmtcfm = max(xmtcfmgrd(i,l)*dens,zer_min)
          bry = max(brygrd(i,l)*dens,zer_min)
!         if(iprn>0.and.(l.eq.1.or.l.eq.63.or.l.eq.kprnin.or.l.eq.64-kprnin))then
!           write(300+iam,*)'latgrd',latgrd(i,j),i,j,longrd(i,j),'iamprn',iamprn
!           write(300+iam,*)'bry',bry,l,'brygrd ','i',i,brygrd(i,l)
!           call flush(300+iam)
!         endif
          ch3br = max(ch3brgrd(i,l)*dens,zer_min)
          h1301 = max(f1301grd(i,l)*dens,zer_min)
          h1211 = max(f1211grd(i,l)*dens,zer_min)
#ifdef   ISOPRENE_PEROX
          rip = max(ripgrd(i,l)*dens,zer_min)
          prdp = max(prdpgrd(i,l)*dens,zer_min)
          xmrp = max(xmrpgrd(i,l)*dens,zer_min)
#else
          hf = max(hfgrd(i,l)*dens,zer_min)
          cfclo = max(cfclogrd(i,l)*dens,zer_min)
          cf2o = max(cf2ogrd(i,l)*dens,zer_min)
#endif
          ch4 = max(ch4grd(i,l)*dens,zer_min)
          hno4 = max(hno4grd(i,l)*dens,zer_min)
          hocl = max(hoclgrd(i,l)*dens,zer_min)
#ifdef   ISOPRENE_PEROX
          vrp = max(vrpgrd(i,l)*dens,zer_min)
#else
          h2otest = max(h2ogrd(i,l)*dens,zer_min)
#endif
          xno3 = max(xno3grd(i,l)*dens,zer_min)
          xno2 = max(xno2grd(i,l)*dens,zer_min)
!          if(i.eq.1.and.j.eq.1)then
!             write(6,*)'xno2 l',l,xno2,'xno2grd',xno2grd(i,l)
!             call flush(6)
!          endif
!          if (xno2grd(i,l).ne.chemlocal%xno2grd(i,l))then
!            if(errorcount<10)then
!            write(6,*)iam,'error local xno2grd',i,l,xno2grd(i,l),chemlocal%xno2grd(i,l)
!            call flush(6)
!            errorcount=errcount+1
!            endif
!          endif
          ch2o = max(ch2ogrd(i,l)*dens,zer_min)
          ch3ooh = max(ch3oohgrd(i,l)*dens,zer_min)
          xco = max(cogrd(i,l)*dens,zer_min)
!          write(200+iam,*)'2012'
!          call flush(200+iam)

#ifdef COTRACERS
          xcod25 = max(cod25grd(i,l)*dens,zer_min)
#ifdef XBBCOD25
          xbbcod25 = max(bbcod25grd(i,l)*dens,zer_min)
#endif
#endif
#ifdef COTRACERS50
          xcod50 = max(cod50grd(i,l)*dens,zer_min)
          xbbcod50 = max(bbcod50grd(i,l)*dens,zer_min)
#endif
!          xcoanth25=max(coanth25grd(i,l)*dens,zer_min)
          xcoanth25=max(coanth25grd(i,l)*dens,0.0)
!          xbbcod25=max(bbcod25grd(i,l)*dens,zer_min)
          xbbcod25=max(bbcod25grd(i,l)*dens,0.0)
          hbr = max(hbrgrd(i,l)*dens,zer_min)
          brno3 = max(brno3grd(i,l)*dens,zer_min)
          hobr = max(hobrgrd(i,l)*dens,zer_min)
          brcl = max(brclgrd(i,l)*dens,zer_min)
!         if(iprn>0.and.(l.eq.1.or.l.eq.63.or.l.eq.kprnin.or.l.eq.64-kprnin))then
!           write(300+iam,*)'brcl',brcl,l,'brclgrd ','i',i,brclgrd(i,l)
!           call flush(300+iam)
!         endif
          cl2 = max(cl2grd(i,l)*dens,zer_min)
c (JRO)   add nmhc species definitions here:
          c2h6 = max(c2h6grd(i,l)*dens,zer_min)
          ald2 = max(ald2grd(i,l)*dens,zer_min)
          ethooh = max(ethoohgrd(i,l)*dens,zer_min)
          pan = max(pangrd(i,l)*dens,zer_min)
          par = max(pargrd(i,l)*dens,zer_min)
          xonit = max(xonitgrd(i,l)*dens,zer_min)
          aone = max(aonegrd(i,l)*dens,zer_min)
          rooh = max(roohgrd(i,l)*dens,zer_min)
          xmgly = max(xmglygrd(i,l)*dens,zer_min)
          eth = max(ethgrd(i,l)*dens,zer_min)
          xolet = max(xoletgrd(i,l)*dens,zer_min)
          xolei = max(xoleigrd(i,l)*dens,zer_min)
          xisop = max(xisopgrd(i,l)*dens,zer_min)
          xisoprd = max(xisoprdgrd(i,l)*dens,zer_min)
          prop_par = max(prop_pargrd(i,l)*dens,zer_min)
          ch3oh = max(ch3ohgrd(i,l)*dens,zer_min)
          xmvk = max(xmvkgrd(i,l)*dens,zer_min)
          xmacr = max(xmacrgrd(i,l)*dens,zer_min)
          xmpan = max(xmpangrd(i,l)*dens,zer_min)
          xnat=0.
          xice=0.
cJAA      add hcn, ch3cn, aerosol species
          hcn = max(hcngrd(i,l)*dens,zer_min)
          ch3cn = max(ch3cngrd(i,l)*dens,zer_min)
#ifdef   DO_AEROSOL
          so2 = max(so2grd(i,l)*dens,zer_min)
          so4aer = max(so4aergrd(i,l)*dens,zer_min)
          dms = max(dmsgrd(i,l)*dens,zer_min)
          xmsa = max(msagrd(i,l)*dens,zer_min)
          xno3aer = max(no3aergrd(i,l)*dens,zer_min)
          xnh3 = max(nh3grd(i,l)*dens,zer_min)
          xnh4aer = max(nh4aergrd(i,l)*dens,zer_min)
          bc1 = max(bc1grd(i,l)*dens,zer_min)
          bc2 = max(bc2grd(i,l)*dens,zer_min)
          oc1 = max(oc1grd(i,l)*dens,zer_min)
          oc2 = max(oc2grd(i,l)*dens,zer_min)
          du1 = max(du1grd(i,l)*dens,zer_min)
          du2 = max(du2grd(i,l)*dens,zer_min)
          du3 = max(du3grd(i,l)*dens,zer_min)
          du4 = max(du4grd(i,l)*dens,zer_min)
          du5 = max(du5grd(i,l)*dens,zer_min)
#ifdef AEROSOL_SEASALT
          ss1 = max(ss1grd(i,l)*dens,zer_min)
          ss2 = max(ss2grd(i,l)*dens,zer_min)
          ss3 = max(ss3grd(i,l)*dens,zer_min)
          ss4 = max(ss4grd(i,l)*dens,zer_min)
          ss5 = max(ss5grd(i,l)*dens,zer_min)
#endif
#endif
#ifdef DO_AEROSOL_BR2
!          write(200+iam,*)'2079'
!          call flush(200+iam)
          br2 = max(br2grd(i,l)*dens,zer_min)
!          write(200+iam,*)'set br2',br2,'i',i,l
!          call flush(200+iam)
#endif

c         flag for tropopause: ktf +1 in troposphere, -1 otherwise
cJAA      20070208 begin
cJAA      Use 380K as switch between trop chem and strat chem (incl met h2o)
cJAA      NOTE: still using troptheta climo to limit BB inj height and lightning nox
cJAA
c         if (theta .lt. troptheta(jlat)) ktf = 1
          if (theta .lt. 380.) ktf = 1
cJAA      20070208 end
c
c         flag for surface: ksf +1 for bottom layer, -1 otherwise
          if (ip .eq. lbot) ksf = 1
c         flag for day or night conditions: kdn -1 for night, +1 for day
c         daytime at all levels if zeta > 0                       (SZA < 90)
c         daytime above level zshad if -0.068 < zeta < 0   (93.9 > SZA > 90)
          kdn = -1
!          write(200+iam,*)'zeta',zeta
!          call flush(200+iam)
#ifdef DIAGFV3
          if(ibeg+i-1.eq.1.and.j.eq.1.and.tile.eq.1)then
            write(6,*)'zeta',zeta,'theta',theta,'ktf',ktf,'ksf',ksf
            write(6,*)'radea',radea
            call flush(6)
          endif
#endif
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!           write(200+iam,*)'yzeta',zeta
!            call flush(200+iam)
!          endif
          if (zeta .ge. -0.068) then
            zshad = 0.0
            if (zeta .le. 0.)
     1        zshad = radea/1000.*(1.0/sqrt(1.0-zeta**2) - 1.0)
CCCMNN 0802   KDN SET TO -1 (NIGHT) FOR ZETA < 0. IN THE TROPOSPHERE
            if (ktf .eq. -1) then
              if (zeta .ge. 0.0 .or. (zgeo(ip) .ge. zshad)) kdn=1
            else
              if (zeta .ge. 0.0) kdn = 1
            endif
#ifdef DIAGFV3
            if(ibeg+i-1.eq.1.and.j.eq.1.and.tile.eq.1)then
              write(6,*)'zshad',zshad,'kdn',kdn,'zgeo',ip,zgeo(ip),'kdn',kdn
              call flush(6)
            endif
#endif
          endif

c         Get altitude B->T index (IMPACT grid) for clim. data sets
c         xlsig is ground to top
          tol=999.
          pressl=alog(press)
          do ll = 1,inl
            plarcl = -xlsig(ll)
            diff=abs(plarcl-pressl)
            if (diff .lt. tol) then
              mpr=ll
              tol=diff
            endif
          enddo

c         get climatologically-specified values (h2o, co, h2, aerosol SAD)
          rlwi=zlwigrd(i,j)
!          write(200+iam,*)'call clim_spec ',jlat,'ireac',ireac,ihreac,iphot,inum,zeta,pur
!          call flush(200+iam)
          call clim_spec(julday,jlat,dens,mpr,ktf,rlwi)
          if (ktf .eq. 1) sad=0.0

c         replace climatological values with met h2o values in troposphere
#ifdef ISOPRENE_PEROX
c         To convert water mass mixing ratio to vmr:
c         vmr = mmr * MW_DRYAIR / MW_H2O
c
cJAA 20070208 Now using met h2o up to 380K, use min threshhold of 3ppmv
          if (ktf .eq. 1)
cJAA           &    dn(13) = max(qgrd(i,j,l)*XMW_DRYAIR/XMW_H2O*dens,zer_min)
     &    dn(13) = max(qgrd(i,j,l)*XMW_DRYAIR/XMW_H2O,3.0e-6)*dens
#else
          if (ktf .eq. 1) dn(13) = h2otest
#endif

cjjj      dn(36) = max(h2grd(i,l)*dens,zer_min)
c         use climatological background SAD
c         if (xsad .lt. 1.e24) sad = xsad
c         if (xsad .lt. 1000.) sad = xsad

cjjj      colo3grd(i,j,l)=0. !cjjj check of colo3

c         Initialize photo array
          do ii=1,iphot
            phoj(ii)=0.
          enddo
          

          if (kdn .gt. 0) then
c
c  ***      daytime ***    Compute photolysis frequencies (J values)
c
            colo3 = svcolo3(ip)
            colo3_vert = svcolo3v(ip)
cjjj          colo3grd(i,j,l)=colo3_vert/dob_to_den !cjjj check of colo3

c           fill photolysis array with either fastj2 or kawa values
            if (usefastj) then
              ! Get Fast-J2 values
              do ii=1,iphot
                phoj(ii)=zpj(lbot-ip+1,ii)
                if(isnan(phoj(ii)))then
                    write(6,*)'NAN phoj ',lbot-ip+1,'lbot',lbot,'ip',ip,ii,phoj(ii)
                    write(6,*)'iatchem',iatchem,katchem,'j',j,'kdn',kdn,'ktf',ktf,'l',l
                    call flush(6)
                    call killit('nan phoj')
                endif
              end do
            else
              ! Get Kawa lookup-table values
c
c             calculate o2 slant column
              scalht=287.*temp*100./9.81
              z1=amax1(zgeo(ip),0.)
!             if sza le 60 degrees then use zeta, otherwise calculate chapman function
              if (sza .gt. 1.0472) then
                chapo2=chap(z1,sza,scalht)
              else
                chapo2=1./zeta
              end if
              colo2=0.209476*dens*scalht*chapo2
c
c             find pressure indices in IMPACT standard grid (sigma*1000) for
c             interpolation of radiative flux function to the local pressure
              do jlrc = 1,inl
                plarc = exp(-xlsig(jlrc))
                if (pur .ge. plarc) goto 222
              enddo
  222         continue
c
c             If actual pressure is greater than max tabulated pressure (ground,
c             jlrc=1) or less than min tabulated pressure (top of atm, jlrc=inl+1),
c             use nearest tabulated pressure without extrapolation (pfac=1).
c             Otherwise, calculate pressure interpolating factor pfac.
              mpru = jlrc
              pfac = 1.
              if (jlrc .ge. 2 .and. jlrc .le. inl) then
                pfac=(alog(pur)+xlsig(jlrc-1))/(xlsig(jlrc-1)-xlsig(jlrc))
              else if (jlrc .eq. inl+1) then
                mpru=inl
              endif
c
c             calculate photolysis rates using kawa lookup tables
CCCMOD        3/14/2002 TROPOSPHERE INDEX KTF INCLUDED IN THE CALL TO JCALC
CCC           ALBEDO = 10% FOR THE TROPOSPHERE, 30% FOR THE STRATOSPHERE.
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!             write(200+iam,*)'call jcalc '
!             call flush(200+iam)
!          endif
              call jcalc(mpru,pfac,sza_deg,colo3_vert,colo3,
     1               colo2,slntno,temp,dens,rjs,ktf)
c
              do ii=1,nxsec           !nxsec is #cross-secs defined in lookup table
                phoj(ii)=rjs(ii)
              end do

              phoj(33)=0.             !hbr photolysis ignored: no data in JPL97

c             equivalence some nmhc rates to those existing in lookup table
              phoj(47)=0. !for consistency with phoj(48),phoj(49),phoj(50) being 0.
              phoj(51)=phoj(28)
              phoj(52)=phoj(28)
              phoj(53)=9.64*phoj(24)

c             use nicolet (PSS, 33, 69-80, 1985) parameterization for J(CH4)
c             at lyman-alpha
              expfac=2.115e-18*(colo2**0.8855)
              phoj(41)=5.55e-6*exp(-expfac)
c             use nicolet (PSS, 1985) for H2O Lyman-alpha (quiet sun)
c             yield to H+OH channel is about 70% at lyman-alpha
              xjh2o_lya=5.7e-6*exp(-expfac)
              xjh2o_srb=1.2e-6*exp(-1.0e-7*(colo2**0.35))
              phoj(8)=0.7*xjh2o_lya+xjh2o_srb
c
            end if        !end of fill photolysis array, fastj2 or kawa
c
          end if
c
cc        get surface deposition rates, and tropospheric no and co emission rates.
c
c
          depo3 = 0.0
          dephno3 = 0.0
          deph2o2 = 0.0
          depch2o = 0.0
          depch3ooh = 0.0
          depno2 = 0.0
          depco = 0.0
          sorcno = 0.0
          sorcco = 0.0
#ifdef COTRACERS
          sorccod5025 = 0.0
          sorcbbcod5025 = 0.0
#endif
          sorccoanth25=0.0
          sorcbbcod25=0.0
          trht = 0.0
          drip = 0.0
c    (JRO) different rainout for ch3ooh in lower trop
          drip_ch3ooh=0.0
c    (JRO) new source and dep rates for nmhc chemistry
          sorcc2h6=0.0
          sorcpar=0.0
          sorcaone=0.0
          sorceth=0.0
          sorcolet=0.0
          sorcolei=0.0
          sorcisop=0.0
          sorcch2o=0.0
          sorcald2=0.0
          sorcproppar = 0.0
          sorcch3oh = 0.0
          depald2=0.0
          depethooh=0.0
          deprooh=0.0
          depmgly=0.0
c
          dephcn=0.0
          depch3cn=0.0
          sorchcn=0.0
          sorcch3cn=0.0
c
#if (defined DO_AEROSOL || defined DO_AEROSOL_DEFINE)
          depso2=0.0
          depso4aer=0.0
          depdms=0.0
          depxmsa=0.0
          depxno3aer=0.0
          depxnh3=0.0
          depxnh4aer=0.0
          depbc1=0.0
          depbc2=0.0
          depoc1=0.0
          depoc2=0.0
          depdu1=0.0
          depdu2=0.0
          depdu3=0.0
          depdu4=0.0
          depdu5=0.0
          depss1=0.0
          depss2=0.0
          depss3=0.0
          depss4=0.0
          depss5=0.0
c
          sorcso2=0.0
          sorcso4aer=0.0
          sorcdms=0.0
          sorcxmsa=0.0
          sorcxno3aer=0.0
          sorcxnh3=0.0
          sorcxnh4aer=0.0
          sorcbc1=0.0
          sorcbc2=0.0
          sorcoc1=0.0
          sorcoc2=0.0
          sorcdu1=0.0
          sorcdu2=0.0
          sorcdu3=0.0
          sorcdu4=0.0
          sorcdu5=0.0
          sorcss1=0.0
          sorcss2=0.0
          sorcss3=0.0
          sorcss4=0.0
          sorcss5=0.0
#endif

CCCMNN0508        ***  DEP. VELOCITY OF HOBr ON ICE.  USED FOR Br2 SOURCE
!          write(200+iam,*)'zero dv_hobr_ice'
!          call flush(200+iam)
          dv_hobr_ice = 0.0
c    4-20-09 tks
          br_ice_fctr = 0.0
c      goto 998  !jjj skip emissions, deposition
c
          if (ksf .eq. 1) then
C
cjjj        fcdr = 1.0e-2/draw(klon,jlat)
            if(draw2d(i).lt.1.e-20.or.draw2d(i).gt.10.)then
              write(6,*)'draw oob ',i,draw2d(i)
              call flush(6)
            endif
            fcdr = 1.0e-2/draw2d(i)

            dz = (zgeoh(lbot-1)-zgeoh(lbot))*1.0e5
CCCC
            ddum = dpo3(i,j,imon)
            depo3 = ddum/((1.0+ddum*fcdr)*dz)
!            write(200+iam,*)'depo3',depo3,'ddum',ddum,'fcdr',fcdr,'dz',dz,' cm'
!            write(200+iam,*)'ddum*fcdr',ddum*fcdr
!            dz2 = (zgeoh(lbot-2)-zgeoh(lbot))*1.0e5
!            write(200+iam,*)'dz2',dz2,'depo32 ',ddum/((1.0+ddum*fcdr)*dz2)
!            write(200+iam,*)'ratio',dz2/dz
!            call flush(200+iam)
!            write(6,*)'ddum',ddum,i,j,imon,'dz',dz,'fcdr',fcdr,
!     &      'zgeoh',lbot-1,lbot,zgeoh(lbot-1:lbot),'depo3',depo3
!            call flush(6)
CCC
            ddum = dphno3(i,j,imon)
            dephno3 = ddum/((1.0+ddum*fcdr)*dz)
CCC
            ddum = dph2o2(i,j,imon)
            deph2o2 = ddum/((1.0+ddum*fcdr)*dz)
CCC
            ddum = dpno2(i,j,imon)
            depno2 = ddum/((1.0+ddum*fcdr)*dz)
CCC
            ddum = dpco(i,j,imon)
            depco = ddum/((1.0+ddum*fcdr)*dz)
CCC
            depch2o = depo3
            depch3ooh = deph2o2*0.5
c     (JRO)  add deposition assumptions for nmhc aldehyde and peroxide species
            depald2=depch2o
            depmgly=depch2o
            depethooh=depch3ooh
            deprooh=depch3ooh
c
CCCMNN 20080605  DEPOSITION VELOCITY FOR HCN AND CH3CN OVER WATER
CCCMNN      DEPOSITION AT THE SURFACE IS ASSUMED TO BE 0.13 CM/SEC
c
c      
c 4-20-09 tks        if (ilwi .eq. 0) then
            if (zlwigrd(i,j) .lt. 1.5) then
              dephcn =  0.13/((1.0+0.13*fcdr)*dz)
              depch3cn = dephcn
            endif
c      
c 4-20-09 tks
CCCMNN0508   set dv_hobr_ice if ilwi = 2  (ice)
            if (zlwigrd(i,j) .gt. 2.5) then
              if(dz.eq.0)then
                write(6,*)'dz zero',i,j
                call flush(6)
              endif
!              write(6,*)'draq2d',i,draw2d(i),'j',j
!              call flush(6)
              dv_hobr_ice=draw2d(i)*1.0e2/dz
!              write(200+iam,*)'dv_hobr_ice ice ',dv_hobr_ice
!              call flush(200+iam)
CCCMNN0209    First year seaice contact correlates better with BrO measurements
CCC           at Barrow.  During 2008, BrO values are higher in March than in
CCC           April, and they seem to drop off in May  (Bill Simpson).
CCC           We approximate this effect by a factor which drops off
CCC           exponentially after April 15  (jday=106) with a 15-day time
CCC           constant.
CCC
              br_ice_fctr=min( exp(-(dy_jdgmt-106.)/15.),1.0)
CCCMNN0209    seaice at lat>78N is taken to be aged and therefore not
CCCC          included as a source bromine
              if (rlat .gt. 78.0) br_ice_fctr = 0.0
            endif
CCMNN0708   seasalt source only at the surface
#ifdef AEROSOL_SEASALT
            sorcss1 = sorcaero(idss1)/dz
            sorcss2 = sorcaero(idss2)/dz
            sorcss3 = sorcaero(idss3)/dz
            sorcss4 = sorcaero(idss4)/dz
#endif
c 4-20-09 tks end
           ch4_s = ch4clim(i,j,imon)
c
          endif ! ksf
c
c         Emissions through depth of PBL rather than just lowest layer
c         Handle biomass burning separately to use different injection ht
c
          pblkm = pblht(i,j)*1.e-3                 !pbl thickness (km)
          emdepth = pblkm
          zlo = zgeoh(ip)-zgeoh(lbot)              !height above surface (km)
c 
c JAA 20080829 Test, do emissions only in the lowest model layer
c JAA          So emdepth is dz (km) of lowest layer
          if(twolayer_emissions)then
            emdepth = zgeoh(lbot-2)-zgeoh(lbot)
          else
            emdepth = zgeoh(lbot-1)-zgeoh(lbot)
          endif
c JAA 20080829 End
          if (zlo.lt.emdepth .and. ip .gt. 1) then
            zup = zgeoh(ip-1)-zgeoh(lbot)           !height above surface (km)
            if (zup.le.emdepth) then
              zscale = 1.0
            else
              zscale = (emdepth-zlo)/(zup-zlo)
            endif
c           dz is now emission depth (cm) normalized by gridbox fraction within emdepth
            dz = emdepth*1.e5/zscale

            sorcno = (srcnind(i,j,imonHTAP)+srcnsoil(i,j,imon))/dz
!            if(i.eq.iprn.and.j.eq.jprn)then
!              write(6,*)'sorcno',sorcno,'ind',srcnind(i,j,imonHTAP),'srcsoil',srcnsoil(i,j,imon)
!              write(6,*)'imon',imon,'dz',dz
!            endif
#ifdef DOEMISSION
        if(doemis)then
          emisout(i,j,1)=srcald(i,j,imonHTAP)
          emisout(i,j,2)=srcalkanone(i,j,imonHTAP)
          emisout(i,j,4)=srcbutane(i,j,imonHTAP)
          emisout(i,j,5)=srcch2o(i,j,imonHTAP)
          emisout(i,j,6)=srcco_totl(i,j,imon)
          emisout(i,j,7)=srcethane(i,j,imonHTAP)
          emisout(i,j,8)=srcethene(i,j,imonHTAP)
          emisout(i,j,9)=srchexane(i,j,imonHTAP)
          emisout(i,j,10)=srcnind(i,j,imonHTAP)
          emisout(i,j,12)=srcpentane(i,j,imonHTAP)
          emisout(i,j,13)=srcpropane(i,j,imonHTAP)
          emisout(i,j,14)=srcpropene(i,j,imonHTAP)
        endif
#endif
        sorcco = (srcco_totl(i,j,imon))/dz

c
        sorcc2h6=srcethane(i,j,imonHTAP)/dz
c industrial par (4butane+5pentane+6hexane+1propene)
        sorcpar=(4.*srcbutane(i,j,imonHTAP)+
     1   5.*srcpentane(i,j,imonHTAP)+6.*srchexane(i,j,imonHTAP)+
     2   1.*srcpropene(i,j,imonHTAP))/dz
        sorceth=srcethene(i,j,imonHTAP)/dz
        sorcolet=srcpropene(i,j,imonHTAP)/dz
        sorcolei=0.0
C (MNN) *** include isoprene source. Source strength depends on
C       *** solar zenith angle (or time of day) and is set
C       *** equal to zero during night time
        sorcisop = srcisop(i,j,imon) * diurngrd(i,j)/dz
        sorcch2o=srcch2o(i,j,imonHTAP)/dz
        sorcald2=srcald(i,j,imonHTAP)/dz
c (JRO!!!) NEED TO PUT IN A SOURCE OF ACETONE!!!!!???
        sorcaone=0.0
        sorcproppar=3.0*srcpropane(i,j,imonHTAP)/dz
c mnn    set methanol source to zero
        sorcch3oh=0.0
c
c JAA 20050907 Add additional biogenic CO sources from methanol and
c     monoterpene oxidation.  Use diurnal scaling as for isoprene
c
c     Following GMI/Harvard recommendation, want to emit 100 tg/yr of CO
c     from methanol oxidation.  Scale from isoprene flux; assume that
c     our current isoprene emission is approx 500 tg/yr.
c
c     CO_from_ch3oh = (100tg CO)/(500tg isop) * srcisop(molec isop/cm2/s)
c      * (5 molec C/molec isop) * (12g C/Av molec C) * (Av molec CO/28g CO)
c      = (100/500) * sorcisop * 5 * (12/28)
c      = (3/7) * sorcisop with units of molec CO/cm2/s
c
c     Following GMI/Harvard recommendation, yield of CO from monoterpene
c     is 20% on carbon number basis.
c
c     CO_from_terp = srcterp(molec terp/cm2/s) * (10 molec C/molec terp) * 0.2
c      = 2 * srcterp with units of molec CO/cm2/s
c
        sorcco = sorcco + 3./7.*sorcisop
     &           + 2.*srcterp(i,j,imon)*diurngrd(i,j)/dz
#ifdef COTRACERS
c      ajl new source for passive tracers cod50 and cod25 do not include fires which were added to sorcco below
        sorccod5025=sorcco
!      if(sorccod5025>0.0)then
!        write(6,*)'sorccod5025',sorccod5025
!      endif
#endif
       sorccoanth25=sorcco
!       if(sorccoanth25>0)then
!         write(6,*)'set sorccoanth25',i,j,sorccoanth25
!       endif

c JAA 20050907 end

#ifdef AEROSOL_CARBON
#ifdef DOEMISSION
        if(doemis)then
          emisout(i,j,3)=srcbc(i,j,imonHTAP)
          emisout(i,j,11)=srcoc(i,j,imonHTAP)
        endif
#endif
          sorcbc1 = 0.8*srcbc(i,j,imonHTAP)/dz  !phobic
          sorcbc2 = 0.2*srcbc(i,j,imonHTAP)/dz  !philic
          sorcoc1 = 0.5*srcoc(i,j,imonHTAP)/dz  !phobic
          sorcoc2 = 0.5*srcoc(i,j,imonHTAP)/dz  !philic
#ifdef DIAGOC4
          if(i.eq.274.and.j.eq.128)then
            write(6,*)'sorcoc1',sorcoc1,'sorcoc2',sorcoc2
          endif
#endif
c
c   Following gocart, yield of aerosol oc2 from terpenes is 10% on mole/mole carbon basis
c     oc2_from_terp = srcterp(molec terp/cm2/s) * (10 molec C/molec terp) * 0.1
c      = 1 * srcterp with units of molec oc2/cm2/s
c
          sorcoc2 = sorcoc2+srcterp(i,j,imon)*diurngrd(i,j)/dz
#ifdef DIAGOC4
          if(i.eq.274.and.j.eq.128)then
           write(6,*)'srococ2 two',sorcoc2,'terp',srcterp(i,j,imon),'imon',imon
          endif
#endif
#endif
c 4-20-09 tks
c#ifdef AEROSOL_SEASALT
!  new code ajl 12/2/2018 to fix blowup over seaice
!          if(i.eq.iprn.and.j.eq.jprnin)then
!            write(6,*)'twolayer',twolayer_emissions,'ip',ip,'dz',dz
!          endif
          if(twolayer_emissions)then
#ifdef AEROSOL_SEASALT_AOD
!#ifdef AEROSOL_SEASALT
!            write(6,*)'i',i,j,'sorcaero',sorcaero(idss1),'dz',dz,'ip',ip
            sorcss1 = sorcaero(idss1)/dz
            sorcss2 = sorcaero(idss2)/dz
            sorcss3 = sorcaero(idss3)/dz
            sorcss4 = sorcaero(idss4)/dz
#endif
CCCMNN0508   set dv_hobr_ice if ilwi = 2  (ice)
c         if (ilwi .eq. 2) dv_hobr_ice =  draw(i,j)*1.0e2/dz
c fix error ajl 7/29/2008
            if (ilwi .eq. 2) dv_hobr_ice =  draw2d(i)*1.0e2/dz
!            write(6,*)i,j,'dy_jdgmt',dy_jdgmt,'kind',kind(dy_jdgmt)
!            call flush(6)
             br_ice_fctr=min( exp(-(dy_jdgmt-106.)/15.),1.0)
CCCMNN0209    seaice at lat>78N is taken to be aged and therefore not
CCCC          included as a source bromine
             if (rlat .gt. 78.0) br_ice_fctr = 0.0
!          if(i.eq.iprn.and.j.eq.jprnin)then
!            write(6,*)'ip',ip,'dv_hobr_ice',dv_hobr_ice,'br_ice_fctr',br_ice_fctr
!           endif        
          endif
c end 4-20-09 tks
#ifdef AEROSOL_SULFATE
          sorcso2 = sorcaero(idso2)/dz
          sorcso4aer = sorcaero(idso4)/dz
#ifdef DIAGOC4
          if(i.eq.274.and.j.eq.128)then
               write(6,*)'sorcso2',sorcso2,'sorcso4aer',sorcso4aer
          endif
#endif
          sorcdms = sorcaero(iddms)/dz
          sorcxmsa = sorcaero(idmsa)/dz
#endif
#ifdef AEROSOL_DUST
          sorcdu1 = sorcaero(iddu1)/dz
          sorcdu2 = sorcaero(iddu2)/dz
          sorcdu3 = sorcaero(iddu3)/dz
          sorcdu4 = sorcaero(iddu4)/dz
          sorcdu5 = sorcaero(iddu5)/dz
#endif

c      if (j.eq.70) then
c        write(6,*) 'sulfate sources',i,j,ip
c        write(6,*) sorcso2,sorcso4aer,sorcdms,sorcxmsa
c        call flush(6)
c      endif

       endif      !End of "surface" emissions through pbl depth

c Biomass burning emissions now separate to allow different injection ht

!       if(iatchem.eq.83.and.jatchem.eq.134.and.tile.eq.2)then
!          write(6,*)'chemdr doplumereise',doplumerise
!       endif
       if(doplumerise)then
!         if(iam==0)then
!            write(6,*)'doplumerise no scale_factor'
!         endif
!         xbbco=emisco3d(i,j,lbot)
c         don't inject above tropopause       JAA Dec 16, 2005
         if(ltrop>nl.or.lbot>nl)then
           ltrop=min(ltrop,nl)
           lbot=min(lbot,nl)
         endif
         if (emisco3d(i,j,lbot).gt.0.) then
           if(ip.eq.lbot.or.(ip>=ktopco(i,j).and.ip<=kbotco(i,j)))then
             dz=(zgeoh(ip-1)-zgeoh(ip))*1.e5 ! make cm
!             write(6,*)'chem emis',i,j,ip,emisco3d(i,j,ip),'dz',dz
             if(emisco3d(i,j,ip)<0.0)then
               write(6,*)'ktopco',ktopco(i,j),'emislbot',lbot,emisco3d(i,j,lbot)
             endif
             colemisco_chem(i,j)=colemisco_chem(i,j)+emisco3d(i,j,ip)*dt
             xbbco=emisco3d(i,j,ip)/dz ! put back in now to be like gdschem and old raqms
!             if(emisco3d(i,j,ip)/=0.0.and.ip==lbot-1)then
!               write(6,*)'emissions at level 2',ip,'lbot',lbot
!             endif
!            plume rise now has emisco3d in molec/cm3/sec so don't need
!            to divide by dz anymore per Brad 3/29/2022
!            with fix now molec/cm2/sec
             bbc = xbbco*bbefac(i,j,2)
             bbno = xbbco*bbefac(i,j,1)
             if(docofire85percent)then
                sorcco = sorcco+.85*xbbco
#ifdef COTRACERS
                 sorcbbcod5025=.85*xbbco
#endif
                 sorcbbcod25=.85*xbbco
             else 
               sorcco = sorcco+xbbco

#ifdef COTRACERS
               sorcbbcod5025=xbbco
#endif
               sorcbbcod25=xbbco
             endif
!            check if bb_nox_adjustment
             if(allocated(bb_nox_scale_factor))then
               if(bb_nox_scale_factor(i,j)>-1.01)then
!                 write(6,*)'bb_nox_scale_factor',i,j,bb_nox_scale_factor(i,j)
                 bb_nox_adjust(i,j)=bb_nox_adjust(i,j)+bbno*bb_nox_scale_factor(i,j)
                 bbno=bbno*(1.+bb_nox_scale_factor(i,j))
               endif
             endif
             if(donofire85percent)then
               sorcno = sorcno+.85*bbno
             else
               sorcno = sorcno+bbno
             endif
c            MNN,JAA 20060328 fix nmhc biomass burning emissions
c            Emission ratios are carbon mass, so ER's for species treated as
c            carbon-bond do NOT need multiplying by number of carbons while
c            ER's for species treated explicitly NEED to be divided by # carbons
c            sorcc2h6 = sorcc2h6+bbc*.00061
             sorcc2h6 = sorcc2h6+bbc*.00061/2.
c            biomass burning par(4butane+5pentane+propene)
c            sorcpar = sorcpar+ (4.*bbc*.00005+
c     1      5.*.00007*bbc+1.*.00066*bbc)
             sorcpar = sorcpar+(bbc*.00005+bbc*.00007+bbc*.00066)
c            sorceth = sorceth+.00123*bbc
             sorceth = sorceth+.00123*bbc/2.
             sorcolet = sorcolet+.00066*bbc
c            sorcproppar = sorcproppar+3.0*bbc*.00019
             sorcproppar = sorcproppar+bbc*.00019
CCCMNN 20080605  HCN AND CH3CN  FROM BIOMASS BURNING
CCCMNN       EMISSION RATIO WRT  CO FROM BIOMASS BURNING
CCCMNN       USE XBBCO IF DAILYFIRES AND BBC OTHERWISE
CCC
             sorchcn = xbbco * 0.0027
             sorcch3cn = xbbco * 0.002
c
             colemisoc(i,j)=colemisoc(i,j)+bbc*.018*8./9.*dz
#ifdef AEROSOL_CARBON
c
c Simple assumptions for now: have assumed Eaero/Eco2=0.02 (Soja, Cofer)
c  Assume nominal Eco2/Ecarbon=0.90, so Eaero/Ecarbon=.02*.90=0.018
c  Following gocart, assume OC/BC = 8 (mole ratio)
c  Following gocart, assume 80% of BC and 50% of OC emitted as hydrophobic
c
             bbc_aer = bbc*0.018
             bbbc = bbc_aer*1./9.
             bboc = bbc_aer*8./9.
#ifdef AEROSOL_19CARBON
             sorcbc1 = sorcbc1+0.1*bbbc !phobic
             sorcbc2 = sorcbc2+0.9*bbbc  !philic
             sorcoc1 = sorcoc1+0.1*bboc  !phobic
             sorcoc2 = sorcoc2+0.9*bboc  !philic
#else
             sorcbc1 = sorcbc1+0.8*bbbc  !phobic
             sorcbc2 = sorcbc2+0.2*bbbc  !philic
             sorcoc1 = sorcoc1+0.5*bboc  !phobic
             sorcoc2 = sorcoc2+0.5*bboc  !philic
#endif
#endif
           endif
         endif
       else ! do plumerise false
#ifdef DAILYFIRES
c Burning emissions obtained from daily CO files
!         write(6,*)'no plumeris usescale',usescale_factor

         if(usescale_factor)then
           if (kdn .gt. 0) then
!            new code to use Brad's new scale_factor
!            xbbco = bbco_d(i,j)
             xbbco = bbco_d(i,j)*scale_factor(i,j)
             xfirehtkm = firehtkm_d(i,j)
           else
!            xbbco = bbco_n(i,j)
             xbbco = bbco_n(i,j)*scale_factor(i,j)
             xfirehtkm = firehtkm_n(i,j)
           endif
         else
           if (kdn .gt. 0) then
!            new code to use Brad's new scale_factor
             xbbco = bbco_d(i,j)
             xfirehtkm = firehtkm_d(i,j)
           else
             xbbco = bbco_n(i,j)
             xfirehtkm = firehtkm_n(i,j)
           endif
         endif
c AJL new did not have in old code
c JAA 20080829 Test, reduce bb emission height
         if (xfirehtkm.gt.0.) xfirehtkm = xfirehtkm-1.
c JAA 20080829 End
       if (xbbco.gt.0.) then
         emdepth = pblkm+xfirehtkm
c don't inject above tropopause       JAA Dec 16, 2005
         if(ltrop>nl.or.lbot>nl)then
           ltrop=min(ltrop,nl)
           lbot=min(lbot,nl)
         endif
         emdepth = min(emdepth,zgeo(ltrop)-zgeoh(lbot))
         bbno = xbbco*bbefac(i,j,1)
!       ajl this could be done before l, ip loop since won't change
         bbc = xbbco*bbefac(i,j,2)
#else
c Burning emissions obtained from monthly NO climatology
c CO emissions already include climatological burning emissions
c
c (JRO)  Retrieve sources for nmhc species
c  Sources are read in as molec/cm(2)/s of the specified hydrocarbon
c  Also, scale NMHC biomass burning source from NOx
c  Convert NOx to total C burned: NOx*39=CO; CO/.0573=C)
c
       bbno = srcnbio(i,j,imon)
       if (bbno.gt.0.) then
         emdepth = pblkm
         bbc=srcnbio(i,j,imon)*39./.0573
         xbbco = 0.
#endif
         if (zlo.lt.emdepth .and. ip .gt. 1) then
             zup = zgeoh(ip-1)-zgeoh(lbot)           !height above surface (km)
             if (zup.le.emdepth) then
               zscale = 1.0
              else
               zscale = (emdepth-zlo)/(zup-zlo)
              endif
c dz is now emission depth (cm) normalized by gridbox fraction within emdepth
              dz = emdepth*1.e5/zscale
c
            if(docofire85percent)then
              sorcco = sorcco+.85*xbbco/dz

#ifdef COTRACERS
              sorcbbcod5025=.85*xbbco/dz
#endif
              sorcbbcod25=.85*xbbco/dz
            else
              sorcco = sorcco+xbbco/dz

#ifdef COTRACERS
              sorcbbcod5025=xbbco/dz
#endif
              sorcbbcod25=xbbco/dz
            endif
!           check if bb_nox_adjustment
            if(allocated(bb_nox_scale_factor))then
              if(bb_nox_scale_factor(i,j)>-1.01)then
                bbno=bbno*(1.+bb_nox_scale_factor(i,j))
              endif
            endif
            if(donofire85percent)then
              sorcno = sorcno+.85*bbno/dz
            else
              sorcno = sorcno+bbno/dz
            endif
c MNN,JAA 20060328 fix nmhc biomass burning emissions
c  Emission ratios are carbon mass, so ER's for species treated as
c  carbon-bond do NOT need multiplying by number of carbons while
c  ER's for species treated explicitly NEED to be divided by # carbons
c          sorcc2h6 = sorcc2h6+bbc*.00061/dz
            sorcc2h6 = sorcc2h6+bbc*.00061/2./dz
c biomass burning par(4butane+5pentane+propene)
c          sorcpar = sorcpar+ (4.*bbc*.00005+
c     1     5.*.00007*bbc+1.*.00066*bbc)/dz
          sorcpar = sorcpar+(bbc*.00005+bbc*.00007+bbc*.00066)/dz
c          sorceth = sorceth+.00123*bbc/dz
          sorceth = sorceth+.00123*bbc/2./dz
          sorcolet = sorcolet+.00066*bbc/dz
c          sorcproppar = sorcproppar+3.0*bbc*.00019/dz
          sorcproppar = sorcproppar+bbc*.00019/dz
CCCMNN 20080605  HCN AND CH3CN  FROM BIOMASS BURNING
CCCMNN    EMISSION RATIO WRT  CO FROM BIOMASS BURNING
CCCMNN    USE XBBCO IF DAILYFIRES AND BBC OTHERWISE
CCC
#ifdef DAILYFIRES
c
      sorchcn = xbbco * 0.0027/dz
      sorcch3cn = xbbco * 0.002/dz
c
#else
c
      sorchcn = bbc * 0.0573 * 0.0027/dz
      sorcch3cn = bbc * 0.0573 * 0.002/dz
c 
#endif

#ifdef AEROSOL_CARBON
c
c Simple assumptions for now: have assumed Eaero/Eco2=0.02 (Soja, Cofer)
c  Assume nominal Eco2/Ecarbon=0.90, so Eaero/Ecarbon=.02*.90=0.018
c  Following gocart, assume OC/BC = 8 (mole ratio)
c  Following gocart, assume 80% of BC and 50% of OC emitted as hydrophobic
c
          bbc_aer = bbc*0.018
          bbbc = bbc_aer*1./9.
          bboc = bbc_aer*8./9.
#ifdef AEROSOL_19CARBON
          sorcbc1 = sorcbc1+0.1*bbbc/dz  !phobic
          sorcbc2 = sorcbc2+0.9*bbbc/dz  !philic
          sorcoc1 = sorcoc1+0.1*bboc/dz  !phobic
          sorcoc2 = sorcoc2+0.9*bboc/dz  !philic
#else
          sorcbc1 = sorcbc1+0.8*bbbc/dz  !phobic
          sorcbc2 = sorcbc2+0.2*bbbc/dz  !philic
          sorcoc1 = sorcoc1+0.5*bboc/dz  !phobic
          sorcoc2 = sorcoc2+0.5*bboc/dz  !philic
#endif
#endif
           endif
         endif ! xbbco
       endif ! doplumerise

c Get aircraft nox source; sairnox_p is ground to top
c       ajl note sairnox_p are international standard atmoshpere pressures corresponding to 1 km height flight levels
c       0,1,2 .. 18 km  (1013 to 75mb)
       sacnox=0.
C ajl add in co,bc,oc,so2 for htap
       sacco=0.
       sacbc=0.
       sacoc=0.
       sacso2=0.

#ifdef DOEMISSION
       if(doemis)then
       if(i.eq.268.and.j.eq.130)then
          write(6,*)'ip',ip,'press',press,'l',l,'nlnair',nlnair,'sairnoxp',sairnox_p(nlnair)
          call flush(6)
          write(6,*)'shape',shape(emis3dout),'emis3dout',emis3dout(i,:,j,1)
          call flush(6)
       endif
       endif
#endif
C  ajl Note layers are meant to be 1 km thick and in units molecules/ccm2/s/ for 1 km
#ifdef DIAGCHEM
       if(i.eq.iprn.and.j.eq.jprnin)then
         write(300+iam,*)nstepat,katchem,'press ',press,' nlnair ',nlnair,' sairnoxp ',sairnox_p(nlnair) 
         call flush(300+iam)
       endif
#endif
    

       if (press .ge. sairnox_p(nlnair)) then
         do ll = 1,nlnair
!       if(i.eq.iprn.and.j.eq.jprnin)then
!         write(300+iam,*)'ll',l,' sairnox_p ',sairnox_p(ll)
!         call flush(300+iam)
!        endif
           if (press .ge. sairnox_p(ll)) then
#ifdef DOEMSSION
       if(doemis)then
       if(i.eq.268.and.j.eq.130)then
          write(6,*)'ll',ll,'sairnox_p',sairnox_p(ll)
          call flush(6)
       endif
       endif
#endif
 
             if (ll .eq. 1) then
               sacnox = srcnair(i,j,1)
c ajl add in htap
#ifndef OLDNAIR
               sacco = srccoair(i,j,1)
#ifdef AEROSOL_CARBON
               sacbc = srcbcair(i,j,1)
               sacoc = srcocair(i,j,1)
               sacso2 = srcso2air(i,j,1)
#endif
#endif
c ajl end add in htap
             else
#if (defined DOINT || defined DOINTCHEMD || defined ALOG10R8)
               sairnox_p8=sairnox_p(ll)
               plu = alog10(sairnox_p8)
               sairnox_p8=sairnox_p(ll-1)
               pll = alog10(sairnox_p8)
#else
c              ajl 2/28/2008 sairnox_p is r*4
c              serenity gives a r*4 alog10 back
               plu = alog10(sairnox_p(ll))
               pll = alog10(sairnox_p(ll-1))
#endif
               pfac = (alog10(press)-plu)/(pll-plu)
               sacnox = srcnair(i,j,ll)*(1.-pfac)+srcnair(i,j,ll-1)*pfac
#ifdef DOEMISSION
       if(doemis)then
       if(i.eq.268.and.j.eq.130)then
          write(6,*)'sacnox',sacnox,'pfac',pfac,'srcnair',srcnair(i,j,ll-1:ll)
          call flush(6)
       endif
       endif
#endif
#ifndef OLDNAIR
c ajl add in htap
               sacco = srccoair(i,j,ll)*(1.-pfac)+srccoair(i,j,ll-1)*pfac
!       if(i.eq.iprn.and.j.eq.jprnin)then
!          write(300+iam,*)'sacco3',sacco
!          call flush(300+iam)
!       endif
#ifdef AEROSOL_CARBON
               sacbc = srcbcair(i,j,ll)*(1.-pfac)+srcbcair(i,j,ll-1)*pfac
               sacoc = srcocair(i,j,ll)*(1.-pfac)+srcocair(i,j,ll-1)*pfac
               sacso2 = srcso2air(i,j,ll)*(1.-pfac)+srcso2air(i,j,ll-1)*pfac
#endif
c ajl end add in htap
#endif
             endif
c            Divide by layer thickness to get rate, convert km-1 to cm-1
!             if(ip.eq.1)then
!                write(6,*)'ll',ll,'press',press,'sairnox',sairnox_p(ll),'nlnair',nlnair,
!     * 'sairnox_p(nlnair)',sairnox_p(nlnair),'j=',j,'i',i
!                call flush(6)
!             endif
#ifdef DOEMISSION
        if(doemis)then
          if(i.eq.268.and.j.eq.130)then
            write(6,*)'sacnox out',ip,sacnox,'dz',zgeoh(ip-1)-zgeoh(ip)
          call flush(6)
          endif
          emis3dout(i,ip,j,1)=sacnox
#ifndef OLDNAIR
          emis3dout(i,ip,j,2)=sacbc
          emis3dout(i,ip,j,3)=sacco
          emis3dout(i,ip,j,4)=sacoc
          emis3dout(i,ip,j,5)=sacso2
#endif
        endif
#endif
#ifndef OLDNAIR
             sacnox = sacnox*1.e-5/(zgeoh(ip-1)-zgeoh(ip))
             sacco = sacco*1.e-5/(zgeoh(ip-1)-zgeoh(ip))
!       if(i.eq.iprn.and.j.eq.jprnin)then
!          write(300+iam,*)'sacco4',sacco
!          call flush(300+iam)
!       endif
#ifdef AEROSOL_CARBON
             sacbc = sacbc*1.e-5/(zgeoh(ip-1)-zgeoh(ip))
             sacoc = sacoc*1.e-5/(zgeoh(ip-1)-zgeoh(ip))
             sacso2 = sacso2*1.e-5/(zgeoh(ip-1)-zgeoh(ip))
#endif
#endif
! ajl I think delta z should not be there
! amount in a delta should be amount/km * dz
! thus amount in /ccm3 should be (amount/km *dz)/dz and thus dz drops out
#ifdef DOEMISSION
        if(doemis)then
       if(i.eq.268.and.j.eq.130)then
           write(6,*)'sacnox out/dz ',ip,sacnox
          call flush(6)
       endif
          emis3dout(i,ip,j,6)=sacnox
#ifndef OLDNAIR
#ifdef AEROSOL_CARBON
          emis3dout(i,ip,j,7)=sacbc
          emis3dout(i,ip,j,9)=sacoc
          emis3dout(i,ip,j,10)=sacso2
#endif
          emis3dout(i,ip,j,8)=sacco
#endif
        endif
#endif
             goto 223
           endif
         enddo
       endif
  223  continue

c Add lightning and aircraft nox sources
!       if(i.eq.iprn.and.j.eq.jprn)then
!          write(6,*)'add lnox and air',ip,sorclnox(ip),sacnox
!       endif
       sorcno=sorcno+sorclnox(ip)+sacnox
       if(ip.eq.lbot)then
         sourceno(i,j)=sorcno
       endif
#ifndef OLDNAIR
       sorcco=sorcco+sacco
!       if(i.eq.iprn.and.j.eq.jprnin)then
!          write(300+iam,*)'sacco',sacco
!          call flush(300+iam)
!       endif
#ifdef AEROSOL_CARBON
       sorcso2=sorcso2+sacso2
c  Brad said to use same ratios as surface emissions
       sorcbc1=sorcbc1+0.8*sacbc
       sorcbc2=sorcbc2+0.2*sacbc
       sorcoc1=sorcoc1+0.5*sacoc
       sorcoc2=sorcoc2+0.5*sacoc
#endif
#endif
c
c  tropospheric heterogeneous reaction rate
c
       if (ktf .eq. 1 .and.pur .ge. 100.) then
        ipr = int(pur/100.)
        if (ipr .lt. 10) then
         trht=trhet(ipr,jlat)+(trhet(ipr+1,jlat)-trhet(ipr,jlat))*
     1    (alog(pur*1.0e-2/ipr)/alog(float(ipr+1)/ipr))
        else
         trht = trhet(ipr,jlat)
        endif
       endif
c
cc  get rain out rate
       if (ktf .eq. 1) drip = rcraingrd(i,j,l)
       drip_ch3ooh = drip
! ajl change since don't set flag but are using wet removal 10/1/2018
!  #ifdef USEHARVARDWETDEP
      ! Wet removal will be handled in the moist physics routines
      drip = 0.0
      drip_ch3ooh = 0.0
!  #endif

  998 continue  !jjj skip emissions, deposition

#if (defined DO_AEROSOL || defined DO_AEROSOL_DEFINE)
CCCMNN0308    AEROSOL VOLUME, SURFACE AREA, AND EFECTIVE RADIUS
        lbt = lbot-ip+1
        r_su=eff_r(1,lbt)
        r_ssa=eff_r(4,lbt)
        r_ssc=eff_r(5,lbt)
        v_su = v_arsl(1,lbt)
        v_ssa = v_arsl(4,lbt)
        v_ssc = v_arsl(5,lbt)
        sa_su = s_arsl(1,lbt)
        sa_ssa = s_arsl(4,lbt)
        sa_ssc = s_arsl(5,lbt)
#endif
c
        ierr=0
        call chem_rate(temp,dens,kdn)
        do idum=1,ikmult
          chem_tend(idum)=0.
        enddo
!#ifdef DOCHEMJJ
        iatchem=i+ibeg-1
        jatchem=j
        katchem=l
!        if(iam.eq.6.and.iatchem.eq.38.and.jatchem.eq.79)then
!          write(6,*)iatchem,jatchem,'latgrd',latgrd(i,j),longrd(i,j),
!     *'zsurf',zsurf(i,j),'szagrd',szagrd(i,j),'zlwigrd',zlwigrd(i,j),
!     * 'ip',ip,'pres',press,'l',l,'ip',ip
!        !endif
!        if(iam.eq.4.and.iatchem.eq.31.and.jatchem.eq.32)then
!          write(6,*)iatchem,jatchem,'latgrd',latgrd(i,j),longrd(i,j),
!     *'zsurf',zsurf(i,j),'szagrd',szagrd(i,j),'zlwigrd',zlwigrd(i,j),
!     * 'ip',ip,'pres',press,'l',l,'ip',ip
!        endif
!#endif
        kprn=l
#ifdef DOCHEMJ
!       l=1 is bottom and nl is top of model
!        if(i.eq.1.and.j.eq.121)then
!         write(6,*)'call do_chem i',i,j,'l',l,'ksf',ksf
!        endif
#ifdef DIAGFV3
        if(ibeg+i-1.eq.1.and.j.eq.1.and.tile.eq.1)then
          write(6,*)'call dochem kdn',kdn,'ktf',ktf,'ksf',ksf,'l',l
          write(6,*)'lon lat',xgrid(i,j),ygrid(i,j)
          call flush(6)
        endif
#endif
!          if(nstepat==84.and.iatchem==168.and.jatchem==38.and.katchem==1.and.tile.eq.2)then
!           write(200+iam,*)'call dochem ',kdn,ktf,ksf,l
           !call flush(200+iam)
!          endif
!        if(iatchem.eq.139.and.jatchem.eq.101.and.l.eq.1.and.tile.eq.1)then
!          write(6,*)iam,'ibeg',ibeg,ibeg+i-1,'i',i,'j',j,'l',l
!          write(6,'("l ",i2," dens ",e12.5," rlat ",f10.2)')l,dens,rlat
!          call flush(6)
!        endif
!        write(6,*)'call do_chem',ibeg+i-1,j,l
!        call flush(6)
        call do_chem(temp,dens,rlat,kdn,ktf,ksf,ierr,ibeg+i-1,j,l)
!          if(nstepat==84.and.iatchem==168.and.jatchem==38.and.katchem==1.and.tile.eq.2)then
!           write(200+iam,*)'did call dochem ',kdn,ktf,ksf,l
!           call flush(200+iam)
!          endif
#else
!        write(200+iam,*)'call do_chem',i,j,l
!        call flush(200+iam)
!        write(200+iam,*)'cv_hobr_icd',dv_hobr_ice,br_ice_fctr
!        call flush(200+iam)
!          if(nstepat==84.and.iatchem==168.and.jatchem==38.and.katchem==1.and.tile.eq.2)then
!           write(200+iam,*)'call dochem ',kdn,ktf,ksf,l
!           call flush(200+iam)
!          endif
!        if(iatchem.eq.139.and.jatchem.eq.101.and.l.eq.1.and.tile.eq.1)then
!          write(6,*)iam,'ibeg',ibeg,ibeg+i-1,'i',i,'j',j,'l',l
!          write(6,'("l ",i2," dens ",e12.5," rlat ",f10.2)')l,dens,rlat
!          call flush(6)
!        endif
        call do_chem(temp,dens,rlat,kdn,ktf,ksf,ierr)
!          if(nstepat==84.and.iatchem==168.and.jatchem==38.and.katchem==1.and.tile.eq.2)then
!           write(200+iam,*)'didcall dochem ',kdn,ktf,ksf,l
!           call flush(200+iam)
!          endif
!        write(200+iam,*)'did do_chem'
!        call flush(200+iam)
#endif
        idid=0
#ifndef SKIPISNAN
        ncount=0
        do idum=1,ikmult
          if (isnan(chem_tend(idum))) then
!            if(idid .eq. 0) then
!            write(6,*) 'aerosol volume, etc'
!            write(6,*) r_su,r_ssa,r_ssc,v_su,v_ssa,v_ssc,sa_su,
!     &                 sa_ssa,sa_ssc
!            idid=1
!            endif
            ncount=ncount+1
            write(6,*) 'NaN val',idum,'i=',i,j,ip,ilwi,kdn,ktf,ksf
            write(6,*)'iprn',iprn,jprn,kprn,'tile',tile
            write(6,*)'iatchem',iatchem,jatchem,katchem
            call flush(6)
!            stop
          endif
        enddo
        if(ncount>0)then
          call killit('nan aerosol')
!          stop
        endif
#endif
       if (ierr .ne. 0) then
          print*,'jlat,klon,mp ',jlat,klon,mp
          print*,'pur,temp,sza ',pur,temp,acos(zeta)*180./3.14159
        endif
c      if (usefastj .and. zeta .ge. -0.068) then
c        print*,'phoj'
c        print*,phoj
c        print*,'dn'
c        print*,dn
c        call flush(6)
c      endif

cJAA 20050325 Store time-avg quantities in units of vmr
cJAA 20050325 Store these diagnostic quantities top->bottom (*,*,ip)
cJAA  instead of according to local grid (*,*,l)

c save Ox photochemical tendency
!        o3ften_save(i,j,ip) = o3ften_save(i,j,ip)+prodloss(1)/dens
!        o3dten_save(i,j,ip) = o3dten_save(i,j,ip)+(prodloss(2)
!     +       +prodloss(3)+prodloss(4)+prodloss(5)+prodloss(6))/dens
#ifdef WETDEPDIAG
        o3ften_save(i,j,ip) = o3ften_save(i,j,ip)+(oxprodloss(1)+oxprodloss(7))/dens
        o3dten_save(i,j,ip) = o3dten_save(i,j,ip)+(oxprodloss(2)
     +       +oxprodloss(3)+oxprodloss(4)+oxprodloss(5)+oxprodloss(6)+oxprodloss(8))/dens
#endif
#ifdef NEWTERM
        o3ften_save(i,j,ip) = o3ften_save(i,j,ip)+(oxprodloss(1)+oxprodloss(7))/dens
        o3dten_save(i,j,ip) = o3dten_save(i,j,ip)+(oxprodloss(2)
     +       +oxprodloss(3)+oxprodloss(4)+oxprodloss(5)+oxprodloss(6)+oxprodloss(8))/dens
!       nmhc o3 gain term
        o3ten_part(i,j,ip,1)=o3ten_part(i,j,ip,1)+oxprodloss(7)/dens
!       nmhc o3 loss term
        o3ten_part(i,j,ip,2)=o3ten_part(i,j,ip,2)+oxprodloss(8)/dens
!       cl o3 loss term
        o3ten_part(i,j,ip,3)=o3ten_part(i,j,ip,3)+oxprodloss(5)/dens
!       br o3 loss term
        o3ten_part(i,j,ip,4)=o3ten_part(i,j,ip,4)+oxprodloss(6)/dens
#ifdef DIAGRESIDUAL
        oxsource=oxprodloss(1)+oxprodloss(7)
        oxsink=oxprodloss(2)+oxprodloss(3)+oxprodloss(4)+oxprodloss(5)+oxprodloss(6)+oxprodloss(8)
        oxdiff=oxsource-oxsink
        o3residual=oxtend-(oxsource-oxsink)
        if(ksf.eq.1)then
!         minus a minus , here deposition was separater from the loss terms
          o3residual=o3residual+oxprodloss(9)
          oxdiff=oxdiff-oxprodloss(9)
        endif
        o3residualdens=o3residual/dens
        if(abs(o3residualdens)>1.e-20)then
!           if(i.eq.1.and.j.eq.121)then
           write(6,*)'ip',ip,'l',l,'ksf',ksf
           write(6,*)'resisual big',i,j,ip,o3residualdens,'o3ften_save',o3ften_save(i,j,ip),o3dten_save(i,j,ip)
           write(6,*)i,j,'oxsourcedens',oxsource/dens,'oxsinkdens',oxsink/dens
           if(oxprodloss(9).ne.00)write(6,*)'depdens',oxprodloss(9)/dens
           write(6,*)i,j,'oxtenddens',oxtend/dens,'oxdiffdens',oxdiff/dens
           write(6,*)i,j,'oxtend',oxtend,'oxdiff',oxdiff
           write(6,*)i,j,'oxresidual',o3residual
           write(6,*)i,j,'oxprodloss',oxprodloss
           call flush(6)
!           endif
        endif
#endif
#endif
#ifdef WETDEPDIAG
        coften_save(i,j,ip) = coften_save(i,j,ip)+coprodloss(1)/dens
        codten_save(i,j,ip) = codten_save(i,j,ip)+coprodloss(2)/dens
        ch4ften_save(i,j,ip)=ch4ften_save(i,j,ip)+ch4prodloss(1)/dens
        ch4dten_save(i,j,ip)=ch4dten_save(i,j,ip)+ch4prodloss(2)/dens
        n2often_save(i,j,ip)=n2often_save(i,j,ip)+n2oprodloss(1)/dens
        n2odten_save(i,j,ip)=n2odten_save(i,j,ip)+n2oprodloss(2)/dens
!        if(n2oprodloss(1)/=0.0.or.n2oprodloss(2)/=0.0)then
!          write(6,*)'n2oprodloss',n2oprodloss(1:2)
!          call flush(6)
!        endif
#endif
#ifdef NEWTERM
        noyften_save(i,j,ip) = noyften_save(i,j,ip)+noyprodloss(1)/dens
        noydten_save(i,j,ip) = noydten_save(i,j,ip)+noyprodloss(2)/dens
        coften_save(i,j,ip) = coften_save(i,j,ip)+coprodloss(1)/dens
        codten_save(i,j,ip) = codten_save(i,j,ip)+coprodloss(2)/dens
        ch4ften_save(i,j,ip)=ch4ften_save(i,j,ip)+ch4prodloss(1)/dens
        ch4dten_save(i,j,ip)=ch4dten_save(i,j,ip)+ch4prodloss(2)/dens
        n2often_save(i,j,ip)=n2often_save(i,j,ip)+n2oprodloss(1)/dens
        n2odten_save(i,j,ip)=n2odten_save(i,j,ip)+n2oprodloss(2)/dens
#endif
#ifdef DIAGCO
        coften_save(i,j,ip) = coften_save(i,j,ip)+coprodloss(1)/dens
        codten_save(i,j,ip) = codten_save(i,j,ip)+coprodloss(2)/dens
        coresidual=(cotend-(coprodloss(1)-coprodloss(2)-coprodloss(3)))/dens
        if(abs(coresidual)>1.e20)then
          write(6,*)'coresidual',i,j,ip,coresidual
        endif
#endif
#ifdef DIAGRESIDUAL
        if(ksf.eq.1)then
!         note noy,co,ch4,n2o deposition is part of sink prodloss(2) for each species
!          noyresidual=(noytend-(noyprodloss(1)-noyprodloss(2)))/dens
!          n2oresidual=(n2otend-(n2oprodloss(1)-n2oprodloss(2)))/dens
!         now term 3 is separate
          noyresidual=(noytend-(noyprodloss(1)-noyprodloss(2)-noyprodloss(3)))/dens
          n2oresidual=(n2otend-(n2oprodloss(1)-n2oprodloss(2)-n2oprodloss(3)))/dens
!         for ch4 and n2o have sfc newtonian term to balance tendendy
          ch4residual=(ch4tend-ch4prodloss(3))/dens
!          coresidual=(cotend-coprodloss(3))/dens
!         all in 1 and 2 3 is part of 2
!          coresidual=(cotend-(coprodloss(1)-coprodloss(2)))/dens
!         now 3 is not part of 2
          coresidual=(cotend-(coprodloss(1)-coprodloss(2)-coprodloss(3)))/dens
        else
          noyresidual=(noytend-(noyprodloss(1)-noyprodloss(2)))/dens
          n2oresidual=(n2otend-(n2oprodloss(1)-n2oprodloss(2)))/dens
          ch4residual=(ch4tend-(ch4prodloss(1)-ch4prodloss(2)))/dens
          coresidual=(cotend-(coprodloss(1)-coprodloss(2)))/dens
        endif
        if(abs(noyresidual)>1.e-20)then
          write(6,*)'noyresidual',i,j,ip,noyresidual
        endif
        if(abs(n2oresidual)>1.e-20)then
          write(6,*)'n2oresidual',i,j,ip,n2oresidual
        endif
        if(abs(coresidual)>1.e-20)then
          
          write(6,*)'coresidual',i,j,ip,'ksf',ksf,coresidual,'cotend',cotend,'coprodloss',coprodloss,'dens',dens
     * ,'cotend/dens',cotend/dens
        endif
        if(abs(ch4residual)>1.e-20)then
c          if(i.eq.1.and.j.eq.97)then
          write(6,*)'ch4residual',i,j,ip,ch4residual,'chtend',chtend,'ch4prodloss',ch4prodloss
c          endif
        endif
#endif
!        write(200+iam,*)'ksf 3049 ',ksf
!        call flush(200+iam)
        if (ksf .eq. 1) then
          o3dep_save(i,j) = o3dep_save(i,j)+oxprodloss(9)/dens
!          if(oxprodloss(9)/=0.0)then
!          write(6,*)'o3dep_save',i,j,o3dep_save(i,j),'add',oxprodloss(9)/dens,'dens',dens
!          call flush(6)
!          endif
          noydep_save(i,j) = noydep_save(i,j)+noyprodloss(3)/dens
          codep_save(i,j) = codep_save(i,j)+coprodloss(3)/dens
#ifdef NEWTERM
          ch4sfcten_save(i,j)=ch4sfcten_save(i,j)+ch4prodloss(3)/dens
          n2osfcten_save(i,j)=n2osfcten_save(i,j)+n2oprodloss(3)/dens
#endif
!          if(i.eq.271)then
!           write(6,*)'j',j,'prodloss12',prodloss(12)
!           call flush(6)
!          endif
        endif
!        write(200+iam,*)'i',i,j,ip,'save',lnox_save(i,j,ip),'sorc',sorclnox(ip),'dens',dens
!        call flush(200+iam)
        lnox_save(i,j,ip) = lnox_save(i,j,ip)+sorclnox(ip)/dens
!        if(sorclnox(ip)>0.0)then
!          write(300+iam,*)'sorc',i,j,ip,sorclnox(ip)
!          call flush(300+iam)
!        endif
        colnox_save(i,j) = 
     *    colnox_save(i,j)+sorclnox(ip)/dens*dpmgrd(i,j,ltb(ip))

c save instantaneous ozone (NOT Ox) volume mixing ratio
        o3vmr_inst(i,j,ip) = dn(1)/dens
!        if(tile.eq.1.and.i.eq.20.and.j.eq.20)then
!           write(300+iam,*)'o3vmr',ip,o3vmr_inst(i,j,ip),' p ',pgrd(i,j,ip)
!           call flush(300+iam)
!        endif
        oh_inst(i,j,ip) = dn(9)/dens
        ho2_inst(i,j,ip) = dn(10)/dens
        no_inst(i,j,ip) = dn(3)/dens
        bro_inst(i,j,ip) = dn(44)/dens
        jo1d_inst(i,j,ip) = phoj(3)
        jno_inst(i,j,ip) = phoj(4)

c  perform forward euler time integration (in units of mixing ratio)
!        if(i.eq.1.and.l.eq.31.and.j.eq.1)then
!         if(i.eq.30.and.j.eq.15)then
!          write(6,*)'chem_tend(1)',l,'ip',ip,chem_tend(1),dens,dt
!          call flush(6)
!        endif
!        write(200+iam,*)'chem_tend',chem_tend,'dens',dens,'dt',dt,shape(chem_tend)
!        call flush(200+iam)
!        write(200+iam,*)'chem_tend(78)',chem_tend(78)
#ifdef DIAGFV3
        if(i.eq.1.and.j.eq.1)then
          write(6,*)'chemtend1',ip,chem_tend(1),'dens',dens,'dt',dt,'xgrid',xgrid(i,j),ygrid(i,j)
          call flush(6)
        endif
#endif
#ifdef DIAGCHEM
        if(iprn>0)then
          do idum=1,ikmult
             if(l.eq.48.or.64-48.eq.l)then
             write(300+iam,*)'chemtend',idum,'l',l,chem_tend(idum)
             endif
          end do 
          call flush(300+iam)
        endif
#endif
!        if(tile.eq.5)then
!          if(iatchem.eq.39.and.jatchem.eq.74)then
!            do idum=1,ikmult
!              write(200+l,*)nstepat,'chemtend',idum,chem_tend(idum)
!            enddo
!          !endif
!        endif
!        if(iam.eq.iamprn.and.ibeg+i-1.eq.iprnin.and.j.eq.jprnin.and.l.eq.62)then
!          write(300+iam,*)'nstepat ',nstepat
!          do idum=1,ikmult
!            if(chem_tend(idum).ne.0.0)then
!              write(300+iam,'(i3," chemtend ",e10.3)')idum,chem_tend(idum)/dens*dt
!            endif
!          end do
!          call flush(300+iam)
!        endif
        o3ten = chem_tend(1)/dens*dt
        xnoyten = chem_tend(2)/dens*dt
        hno3tten = chem_tend(3)/dens*dt
        clyten = chem_tend(4)/dens*dt
        xn2o5ten = chem_tend(5)/dens*dt
        h2o2ten = chem_tend(6)/dens*dt
        hclten = chem_tend(7)/dens*dt
        clno3ten = chem_tend(8)/dens*dt
        ocloten = chem_tend(9)/dens*dt
        xn2oten = chem_tend(10)/dens*dt
        f11ten = chem_tend(11)/dens*dt
        f12ten = chem_tend(12)/dens*dt
        ccl4ten = chem_tend(13)/dens*dt
        ch3clten = chem_tend(14)/dens*dt
        xmtcfmten = chem_tend(15)/dens*dt
        bryten = chem_tend(16)/dens*dt
        ch3brten = chem_tend(17)/dens*dt
        f1301ten = chem_tend(18)/dens*dt
        f1211ten = chem_tend(19)/dens*dt
#ifdef ISOPRENE_PEROX
        ripten = chem_tend(20)/dens*dt
        prdpten = chem_tend(21)/dens*dt
        xmrpten = chem_tend(22)/dens*dt
#else
        hften = chem_tend(20)/dens*dt
        cfcloten = chem_tend(21)/dens*dt
        cf2oten = chem_tend(22)/dens*dt
#endif
        ch4ten = chem_tend(23)/dens*dt
        hno4ten = chem_tend(24)/dens*dt
        hoclten = chem_tend(25)/dens*dt
#ifdef ISOPRENE_PEROX
        vrpten = chem_tend(26)/dens*dt
#else
        h2oten = chem_tend(26)/dens*dt
#endif
        xno3ten = chem_tend(27)/dens*dt
        xno2ten = chem_tend(28)/dens*dt
        ch2oten = chem_tend(29)/dens*dt
        ch3oohten = chem_tend(30)/dens*dt
        coten = chem_tend(31)/dens*dt
#ifdef CO25D

#ifdef COTRACERS
        cod25ten = chem_tend(79)/dens*dt
#ifdef XBBCOD25
        bbcod25ten=chem_tend(80)/dens*dt
#endif
#ifdef DIAGCOTRACERS
        if(i.eq.86.and.j.eq.155)then
          write(6,*)'bbcod50ten',bbcod50ten
        endif
#endif
#endif
#ifdef COTRACERS50
        cod50ten = chem_tend(81)/dens*dt
        bbcod50ten=chem_tend(82)/dens*dt
#endif

#else
        coanth25ten=chem_tend(79)/dens*dt
        bbcod25ten=chem_tend(80)/dens*dt
!        if(coanth25ten/=0)then
!          write(6,*)'coanth25ten',i,j,coanth25ten
!        endif 
!        if(bbcod25ten/=0)then
!          write(6,*)'bbco25ten',i,j,bbcod25ten
!         endif
#endif
!        write(200+iam,*)'chem_tend(78)',chem_tend(78)
        hbrten = chem_tend(32)/dens*dt
        brno3ten = chem_tend(33)/dens*dt
        hobrten = chem_tend(34)/dens*dt
        brclten = chem_tend(35)/dens*dt
#ifdef DIAGCHEM
        if(iprn>0.and.(ip.eq.kprnin.or.ip.eq.64-kprnin))then
         write(300+iam,*)'brclten',brclten,'chemtend35',ip,chem_tend(35)
         call flush(300+iam)
        endif
#endif
        cl2ten = chem_tend(36)/dens*dt
        c2h6ten = chem_tend(37)/dens*dt
        ald2ten = chem_tend(38)/dens*dt
        ethoohten = chem_tend(39)/dens*dt
        panten = chem_tend(40)/dens*dt
        parten = chem_tend(41)/dens*dt
        xonitten = chem_tend(42)/dens*dt
        aoneten = chem_tend(43)/dens*dt
        roohten = chem_tend(44)/dens*dt
        xmglyten = chem_tend(45)/dens*dt
        ethten = chem_tend(46)/dens*dt
        xoletten = chem_tend(47)/dens*dt
        xoleiten = chem_tend(48)/dens*dt
        xisopten = chem_tend(49)/dens*dt
        xisoprdten = chem_tend(50)/dens*dt
        prop_parten = chem_tend(51)/dens*dt
        ch3ohten = chem_tend(52)/dens*dt
        xmvkten = chem_tend(53)/dens*dt
        xmacrten = chem_tend(54)/dens*dt
        xmpanten = chem_tend(55)/dens*dt
c
        hcnten = chem_tend(56)/dens*dt
        ch3cnten = chem_tend(57)/dens*dt
!        write(200+iam,*)'chem_tend(78)',chem_tend(78)
!        write(200+iam,*)'ch3cnten',ch3cnten
!        call flush(200+iam)
#ifdef DO_AEROSOL
        so2ten = chem_tend(58)/dens*dt
        so4aerten = chem_tend(59)/dens*dt
        dmsten = chem_tend(60)/dens*dt
        xmsaten = chem_tend(61)/dens*dt
        xno3aerten = chem_tend(62)/dens*dt
        xnh3ten = chem_tend(63)/dens*dt
        xnh4aerten = chem_tend(64)/dens*dt
        bc1ten = chem_tend(65)/dens*dt
        bc2ten = chem_tend(66)/dens*dt
#ifdef DIAGBC2TEN
        if(i.eq.185.and.abs(j-121)<2.and.ip.eq.20)then
          write(6,*)'bc2ten',bc2ten,'chemtend66',chem_tend(66),'dens',dens,'dt',dt,'ip',ip,'l',l
        endif
#endif
        oc1ten = chem_tend(67)/dens*dt
        oc2ten = chem_tend(68)/dens*dt
        du1ten = chem_tend(69)/dens*dt
        du2ten = chem_tend(70)/dens*dt
        du3ten = chem_tend(71)/dens*dt
        du4ten = chem_tend(72)/dens*dt
        du5ten = chem_tend(73)/dens*dt
#ifdef AEROSOL_SEASALT
        ss1ten = chem_tend(74)/dens*dt
        ss2ten = chem_tend(75)/dens*dt
        ss3ten = chem_tend(76)/dens*dt
        ss4ten = chem_tend(77)/dens*dt
        ss5ten = chem_tend(78)/dens*dt
#endif
#endif
#ifdef DO_AEROSOL_BR2
!        write(200+iam,*)'zz chem_tend(78)',chem_tend(78),dens,dt
!        call flush(200+iam)
        br2ten = chem_tend(78)/dens*dt
!        write(200+iam,*)'yy br2ten',br2ten
!        call flush(200+iam)
#endif

        if(fixtend) then
          xnoyten = - sednat_in/dens*dt
          write(6,*)'in fixtend **(*(*(*(*'
          clyten = 0.
          xn2oten = 0.
          f11ten = 0.
          f12ten = 0.
          ccl4ten = 0.
          ch3clten = 0.
          xmtcfmten = 0.
          bryten = 0.
          ch3brten = 0.
          f1301ten = 0.
          f1211ten = 0.
        endif
c
c overwrite haloe no with equilibrium value
cjjj        xnogrd(i,j,l) = dn(3)/dens
c Store diagnostic values
cjjj        ohgrd(i,j,l) = dn(9)/dens
cjjj        ho2grd(i,j,l) = dn(10)/dens
cjjj        pscgrd(i,j,l) = xnat/dens
cjjj        xicegrd(i,j,l)= xice/dens
cjjj        poxgrd(i,j,l) = prodloss(1)
cjjj        doxogrd(i,j,l) = prodloss(2)
cjjj        doxhgrd(i,j,l) = prodloss(3)
cjjj        doxngrd(i,j,l) = prodloss(4)
cjjj        doxcgrd(i,j,l) = prodloss(5)
cjjj        doxbgrd(i,j,l) = prodloss(6)
c
c  perform forward euler time integration (in units of mixing ratio)

c          if(isnan(o3ten))then
c            write(6,*) 'o3NaN val',idum,i,j,ip,dens,dt
c            call flush(6)
c          endif
!        if(ibeg+i-1.eq.80.and.j.eq.91.and.tile.eq.5)then
!          write(6,'(3i3," oxgrd ",i2,e12.4," o3ten ",e12.4)')tile,i,j,l,oxgrd(i,l),o3ten
!          call flush(6)
!        endif
#ifdef DIAGFV3
        if(ibeg+i-1.eq.1.and.j.eq.1.and.tile.eq.1)then
          write(6,'(3i3," oxgrd ",i2,e12.4," o3ten ",e12.4)')tile,i,j,l,oxgrd(i,l),o3ten
          call flush(6)
        endif
#endif
        oxgrd(i,l) = oxgrd(i,l) + o3ten
!        if(tile.eq.5)then
!          if(iatchem.eq.45.and.jatchem.eq.79.and.(l.eq.64.or.l.eq.1))then
!            write(6,*)nstepat,'oxgrd',l,oxgrd(i,l)
!          endif
!        endif
!        if(iatchem.eq.122.and.jatchem.eq.75.and.l.eq.1.and.tile.eq.5)then
!          write(6,*)'oxgrd',i,l,oxgrd(i,l),'o3ten',o3ten
!        endif
!        if(ibeg+i-1.eq.80.and.j.eq.91.and.tile.eq.5)then
!          write(6,'(2i2," oxgrd new ",e12.4," p ",e12.4," z ",e12.4)')i,j,oxgrd(i,l),pgrd(i,j,l),zgrd(i,j,l)
!          call flush(6)
!        endif
#ifdef DIAGFV3
        if(ibeg+i-1.eq.1.and.j.eq.1.and.tile.eq.1)then
          write(6,'(2i2," oxgrd new ",e12.4," p ",e12.4," z ",e12.4)')i,j,oxgrd(i,l),pgrd(i,j,l),zgrd(i,j,l)
          call flush(6)
        endif
#endif
        
!        if(i.eq.1.and.l.eq.31.and.j.eq.1)then
!          write(6,*)'o3ten',o3ten
!          call flush(6)
!        endif
        xnoygrd(i,l) = xnoygrd(i,l) + xnoyten
        hno3tgrd(i,l) = hno3tgrd(i,l) + hno3tten
        clygrd(i,l) = clygrd(i,l) + clyten
        xn2o5grd(i,l) = xn2o5grd(i,l) + xn2o5ten
        h2o2grd(i,l) = h2o2grd(i,l) + h2o2ten
        hclgrd(i,l) = hclgrd(i,l) + hclten
        clno3grd(i,l) = clno3grd(i,l) + clno3ten
        oclogrd(i,l) = oclogrd(i,l) + ocloten
        xn2ogrd(i,l) = xn2ogrd(i,l) + xn2oten 
        f11grd(i,l) = f11grd(i,l) + f11ten
        f12grd(i,l) = f12grd(i,l) + f12ten
        ccl4grd(i,l) = ccl4grd(i,l) + ccl4ten 
        ch3clgrd(i,l) = ch3clgrd(i,l) + ch3clten
        xmtcfmgrd(i,l) = xmtcfmgrd(i,l) + xmtcfmten 
#ifdef DIAGCHEM
        if(iprn>0.and.(l.eq.kprnin.or.l.eq.64-kprnin))then
!          write(250+iam,*)l,'brygrd',brygrd(i,l),' ten ',bryten,'per',bryten/brygrd(i,l)
!          call flush(250+iam)
          write(300+iam,*)l,'brygrd',brygrd(i,l),' ten ',bryten,'per',bryten/brygrd(i,l)
          call flush(300+iam)
        endif
#endif
        brygrd(i,l) = brygrd(i,l) + bryten
        ch3brgrd(i,l) = ch3brgrd(i,l) + ch3brten
        f1301grd(i,l) = f1301grd(i,l) + f1301ten
        f1211grd(i,l) = f1211grd(i,l) + f1211ten
#ifdef ISOPRENE_PEROX
!        write(200+iam,*)'ripten',ripten      
!        call flush(200+iam)
        ripgrd(i,l) = ripgrd(i,l) + ripten
        prdpgrd(i,l) = prdpgrd(i,l) + prdpten
        xmrpgrd(i,l) = xmrpgrd(i,l) + xmrpten
#else
        hfgrd(i,l) = hfgrd(i,l) + hften
        cfclogrd(i,l) = cfclogrd(i,l) + cfcloten
        cf2ogrd(i,l) = cf2ogrd(i,l) + cf2oten
#endif
        ch4grd(i,l) = ch4grd(i,l) + ch4ten
        hno4grd(i,l) = hno4grd(i,l) + hno4ten
        hoclgrd(i,l) = hoclgrd(i,l) + hoclten
#ifdef ISOPRENE_PEROX
        vrpgrd(i,l) = vrpgrd(i,l) + vrpten
#else
        h2ogrd(i,l) = h2ogrd(i,l) + h2oten 
#endif
        xno3grd(i,l) = xno3grd(i,l) + xno3ten
        colno2b=colno2b+dpmgrd(i,jatchem,l)*xno2grd(i,l)
        xno2grd(i,l) = xno2grd(i,l) + xno2ten
        colno2a=colno2a+dpmgrd(i,jatchem,l)*xno2grd(i,l)
        no2ten(i,j,l)=no2ten(i,j,l)+xno2ten
!        if(iatchem.eq.139.and.jatchem.eq.101.and.tile.eq.1)then
!          write(6,*)'dpm',l,dpmgrd(i,jatchem,l)
!        endif
!        if(iatchem.eq.139.and.jatchem.eq.101.and.l.eq.1.and.tile.eq.1)then
!          write(6,*)'xno2grd',l,xno2grd(i,l),'katchem',katchem,'xno2ten',xno2ten
!        endif
        ch2ogrd(i,l) = ch2ogrd(i,l) + ch2oten
        ch3oohgrd(i,l) = ch3oohgrd(i,l) + ch3oohten
        cogrd(i,l) = cogrd(i,l) + coten
!        if(iatchem.eq.122.and.jatchem.eq.75.and.(l.eq.1.or.l.eq.64).and.tile.eq.5)then
!          write(6,*)'cogrd',i,l,cogrd(i,l),'coten',coten
!          call flush(6)
!        endif
#ifdef CO25D

#ifdef COTRACERS
        cod25grd(i,l) = cod25grd(i,l) + cod25ten
#ifdef XBBCOD25
        bbcod25grd(i,l) = bbcod25grd(i,l) + bbcod25ten
#endif
#endif

#else
        coanth25grd(i,l)=coanth25grd(i,l)+coanth25ten
        bbcod25grd(i,l)=bbcod25grd(i,l)+bbcod25ten
#endif
#ifdef COTRACERS50
        cod50grd(i,l) = cod50grd(i,l) + cod50ten
        bbcod50grd(i,l) = bbcod50grd(i,l) + bbcod50ten
#endif
        hbrgrd(i,l) = hbrgrd(i,l) + hbrten 
        brno3grd(i,l) = brno3grd(i,l) + brno3ten
        hobrgrd(i,l) = hobrgrd(i,l) + hobrten
#ifdef DIAGCHEM
        if(iprn>0.and.(l.eq.kprnin.or.l.eq.64-kprnin))then
!         write(250+iam,*)l,'brclgrd',brclgrd(i,l),' ten ',brclten,'per',brclten/brclgrd(i,l)
!         call flush(250+iam)
         write(300+iam,*)l,'brclgrd',brclgrd(i,l),' ten ',brclten,'per',brclten/brclgrd(i,l)
         call flush(300+iam)
        endif
#endif
        brclgrd(i,l) = brclgrd(i,l) + brclten
        cl2grd(i,l) = cl2grd(i,l) + cl2ten
        c2h6grd(i,l) = c2h6grd(i,l) + c2h6ten
        ald2grd(i,l) = ald2grd(i,l) + ald2ten
        ethoohgrd(i,l) = ethoohgrd(i,l) + ethoohten
        pangrd(i,l) = pangrd(i,l) + panten
        pargrd(i,l) = pargrd(i,l) + parten
        xonitgrd(i,l) = xonitgrd(i,l) + xonitten
        aonegrd(i,l) = aonegrd(i,l) + aoneten
        roohgrd(i,l) = roohgrd(i,l) + roohten
        xmglygrd(i,l) = xmglygrd(i,l) + xmglyten
        ethgrd(i,l) = ethgrd(i,l) + ethten
        xoletgrd(i,l) = xoletgrd(i,l) + xoletten
        xoleigrd(i,l) = xoleigrd(i,l) + xoleiten
        if(xoleiten/=0.0)then
          write(6,*)'xoleigrd',i,l,xoleigrd(i,l),'ten',xoleiten
        endif
        xisopgrd(i,l) = xisopgrd(i,l) + xisopten
        xisoprdgrd(i,l) = xisoprdgrd(i,l) + xisoprdten
        prop_pargrd(i,l) = prop_pargrd(i,l) + prop_parten
        ch3ohgrd(i,l) = ch3ohgrd(i,l) + ch3ohten
        xmvkgrd(i,l) = xmvkgrd(i,l) + xmvkten
        xmacrgrd(i,l) = xmacrgrd(i,l) + xmacrten
        xmpangrd(i,l) = xmpangrd(i,l) + xmpanten
c
        hcngrd(i,l) = hcngrd(i,l) + hcnten
        ch3cngrd(i,l) = ch3cngrd(i,l) + ch3cnten
!        write(200+iam,*)'ch3cnten',ch3cnten
!        call flush(200+iam)
#ifdef DO_AEROSOL
        so2grd(i,l) = so2grd(i,l) + so2ten
        so4aergrd(i,l) = so4aergrd(i,l) + so4aerten
        dmsgrd(i,l) = dmsgrd(i,l) + dmsten
        msagrd(i,l) = msagrd(i,l) + xmsaten
        no3aergrd(i,l) = no3aergrd(i,l) + xno3aerten
        nh3grd(i,l) = nh3grd(i,l) + xnh3ten
        if(nh3grd(i,l).ne.0.0)then
         write(6,*)'n3grd',i,j,l,nh3grd(i,l),'tend',xnh3ten
         call flush(6)
        endif
!        if(i.eq.1.and.j.eq.30.and.l.eq.1)then
!          write(6,*)'nh3grd',nh3grd(i,l),'ten',xnh3ten
!        endif
        nh4aergrd(i,l) = nh4aergrd(i,l) + xnh4aerten
        bc1grd(i,l) = bc1grd(i,l) + bc1ten
#ifdef DIAGBC2TEN
        if(i.eq.185.and.abs(j-121)<2.and.ip.eq.20)then
          write(6,*)'b2grd',bc2grd(i,l),'bc2ten',bc2ten
        endif
#endif
        bc2grd(i,l) = bc2grd(i,l) + bc2ten
        oc1grd(i,l) = oc1grd(i,l) + oc1ten
        oc2grd(i,l) = oc2grd(i,l) + oc2ten
        write(6,*)'update du1grd',i,l,du1ten
        du1grd(i,l) = du1grd(i,l) + du1ten
        du2grd(i,l) = du2grd(i,l) + du2ten
        du3grd(i,l) = du3grd(i,l) + du3ten
        du4grd(i,l) = du4grd(i,l) + du4ten
        du5grd(i,l) = du5grd(i,l) + du5ten
#ifdef AEROSOL_SEASALT
        ss1grd(i,l) = ss1grd(i,l) + ss1ten
        ss2grd(i,l) = ss2grd(i,l) + ss2ten
        ss3grd(i,l) = ss3grd(i,l) + ss3ten
        ss4grd(i,l) = ss4grd(i,l) + ss4ten
        ss5grd(i,l) = ss5grd(i,l) + ss5ten
#endif
#endif
#ifdef DO_AEROSOL_BR2
!        write(200+iam,*)'add in ten ss5ten',ss5ten,ss5grd(i,l),i,l
!        call flush(200+iam)
        br2grd(i,l) = br2grd(i,l) + br2ten
#endif

c Apply NAT sedimentation as loss of NOy
c  (ONLY when xnoy tendency not applied above)
c       xnoygrd(i,l) = xnoygrd(i,l) + xnoyten  ! rbp 31 oct 01
c  .. trajectory sedimentation tendency arrays
c       xnoysed(i,j,l) = xnoysed(i,j,l)+xnoyten    ! rbp 31 oct 01

c end altitude loop
  170 continue

!          write(6,*)'colno2a',colno2a,colno2b,'diff',colno2a-colno2b
!        endif
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!        write(200+iam,*)'after 170'
!        call flush(200+iam)
!          endif

#ifdef DO_AEROSOL
c
c Note that "grd" arrays are defined B->T, gocart tracers are B->T
c oh_inst is T->B
c cloud3d is T->B
c gocart wants delz in meters, T->B
c zeta1=cos(sza) is stored above, prior to altitude loop
c
      do ip=1,nl
        l=nl-ip+1
        ohgrd(l) = dble(oh_inst(i,j,ip))
        cldf(l) = dble(cloud3d(i,j,ip))
        if (ip.eq.1) then
          delz(ip) = 2.*(zgeo(ip)-zgeoh(ip))*1000.
        else
          delz(ip) = (zgeoh(ip-1)-zgeoh(ip))*1000.
        endif
      enddo
      iaero = i
!      write(6,*) 'calling raqms_gocart_driver ',iaero,jaero
!      call flush(6)
      call raqms_gocart_driver(iaero,jaero,ndtaero,ltb,delz,ohgrd,
     &      cldf,zeta1,chemlocal)
c      write(6,*) 'back from raqms_gocart_driver ',iaero,jaero
c      call flush(6)
#endif

c Set below-ground values (if any) equal to surface values JAA 6/14/02
      ip=lbot+1
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!        write(200+iam,*)'after 170 ip',ip,'nl',nl,'lbot',lbot
!        call flush(200+iam)
!          endif
      do while (ip .le. nl)
        l=ltb(ip)
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!          write(200+iam,*)' while ',ip,nl,l
!          call flush(200+iam)
!          endif
        oxgrd(i,l) = oxgrd(i,lbot)
        xnoygrd(i,l) = xnoygrd(i,lbot)
        hno3tgrd(i,l) = hno3tgrd(i,lbot)
        clygrd(i,l) = clygrd(i,lbot)
        xn2o5grd(i,l) = xn2o5grd(i,lbot)
        h2o2grd(i,l) = h2o2grd(i,lbot)
        hclgrd(i,l) = hclgrd(i,lbot)
        clno3grd(i,l) = clno3grd(i,lbot)
        oclogrd(i,l) = oclogrd(i,lbot)
        xn2ogrd(i,l) = xn2ogrd(i,lbot)
        f11grd(i,l) = f11grd(i,lbot)
        f12grd(i,l) = f12grd(i,lbot)
        ccl4grd(i,l) = ccl4grd(i,lbot)
        ch3clgrd(i,l) = ch3clgrd(i,lbot)
        xmtcfmgrd(i,l) = xmtcfmgrd(i,lbot)
        brygrd(i,l) = brygrd(i,lbot)
        ch3brgrd(i,l) = ch3brgrd(i,lbot)
        f1301grd(i,l) = f1301grd(i,lbot)
        f1211grd(i,l) = f1211grd(i,lbot)
#ifdef ISOPRENE_PEROX
        ripgrd(i,l) = ripgrd(i,lbot)
        prdpgrd(i,l) = prdpgrd(i,lbot)
        xmrpgrd(i,l) = xmrpgrd(i,lbot)
#else
        hfgrd(i,l) = hfgrd(i,lbot)
        cfclogrd(i,l) = cfclogrd(i,lbot)
        cf2ogrd(i,l) = cf2ogrd(i,lbot)
#endif
        ch4grd(i,l) = ch4grd(i,lbot)
        hno4grd(i,l) = hno4grd(i,lbot)
        hoclgrd(i,l) = hoclgrd(i,lbot)
#ifdef ISOPRENE_PEROX
        vrpgrd(i,l) = vrpgrd(i,lbot)
#else
        h2ogrd(i,l) = h2ogrd(i,lbot)
#endif
        xno3grd(i,l) = xno3grd(i,lbot)
        xno2grd(i,l) = xno2grd(i,lbot)
        ch2ogrd(i,l) = ch2ogrd(i,lbot)
        ch3oohgrd(i,l) = ch3oohgrd(i,lbot)
!        if(iatchem.eq.122.and.jatchem.eq.75.and.(l.eq.1.or.l.eq.64).and.tile.eq.5)then
!          write(6,*)'store lbot ',cogrd(i,lbot),' into l',l,cogrd(i,lbot),'i',i
!          call flush(6)
!        endif
        cogrd(i,l) = cogrd(i,lbot)
#ifdef COTRACERS
        cod25grd(i,l) = cod25grd(i,lbot)
#ifdef XBBCOD25
        bbcod25grd(i,l) = bbcod25grd(i,lbot)
#endif
#ifdef DIAGCOTRACERS
        if(i.eq.86.and.j.eq.155)then
          write(6,*)'bbcod25grd bot',l,bbcod25grd(i,l)
        endif
#endif
#endif
#ifdef COTRACERS50
        cod50grd(i,l) = cod50grd(i,lbot)
        bbcod50grd(i,l) = bbcod50grd(i,lbot)
#endif
        hbrgrd(i,l) = hbrgrd(i,lbot)
        brno3grd(i,l) = brno3grd(i,lbot)
        hobrgrd(i,l) = hobrgrd(i,lbot)
        brclgrd(i,l) = brclgrd(i,lbot)
        cl2grd(i,l) = cl2grd(i,lbot)
        c2h6grd(i,l) = c2h6grd(i,lbot)
        ald2grd(i,l) = ald2grd(i,lbot)
        ethoohgrd(i,l) = ethoohgrd(i,lbot)
        pangrd(i,l) = pangrd(i,lbot)
        pargrd(i,l) = pargrd(i,lbot)
        xonitgrd(i,l) = xonitgrd(i,lbot)
        aonegrd(i,l) = aonegrd(i,lbot)
        roohgrd(i,l) = roohgrd(i,lbot)
        xmglygrd(i,l) = xmglygrd(i,lbot)
        ethgrd(i,l) = ethgrd(i,lbot)
        xoletgrd(i,l) = xoletgrd(i,lbot)
        xoleigrd(i,l) = xoleigrd(i,lbot)
        xisopgrd(i,l) = xisopgrd(i,lbot)
        xisoprdgrd(i,l) = xisoprdgrd(i,lbot)
        prop_pargrd(i,l) = prop_pargrd(i,lbot)
        ch3ohgrd(i,l) = ch3ohgrd(i,lbot)
        xmvkgrd(i,l) = xmvkgrd(i,lbot)
        xmacrgrd(i,l) = xmacrgrd(i,lbot)
        xmpangrd(i,l) = xmpangrd(i,lbot)
c
        hcngrd(i,l) = hcngrd(i,lbot)
        ch3cngrd(i,l) = ch3cngrd(i,lbot)
#ifdef DO_AEROSOL
        so2grd(i,l) = so2grd(i,lbot)
        so4aergrd(i,l) = so4aergrd(i,lbot)
        dmsgrd(i,l) = dmsgrd(i,lbot)
        msagrd(i,l) = msagrd(i,lbot)
        no3aergrd(i,l) = no3aergrd(i,lbot)
        nh3grd(i,l) = nh3grd(i,lbot)
        nh4aergrd(i,l) = nh4aergrd(i,lbot)
        bc1grd(i,l) = bc1grd(i,lbot)
        bc2grd(i,l) = bc2grd(i,lbot)
        oc1grd(i,l) = oc1grd(i,lbot)
        oc2grd(i,l) = oc2grd(i,lbot)
        du1grd(i,l) = du1grd(i,lbot)
        du2grd(i,l) = du2grd(i,lbot)
        du3grd(i,l) = du3grd(i,lbot)
        du4grd(i,l) = du4grd(i,lbot)
        du5grd(i,l) = du5grd(i,lbot)
#ifdef AEROSOL_SEASALT
        ss1grd(i,l) = ss1grd(i,lbot)
        ss2grd(i,l) = ss2grd(i,lbot)
        ss3grd(i,l) = ss3grd(i,lbot)
        ss4grd(i,l) = ss4grd(i,lbot)
        ss5grd(i,l) = ss5grd(i,lbot)
#endif
#endif
#ifdef DO_AEROSOL_BR2
        br2grd(i,l) = br2grd(i,lbot)
#endif
        ip=ip+1
      enddo

c end longitude loop
  160 continue
!          if(nstepat==84.and.iatchem.eq.168.and.jatchem.eq.38.and.tile.eq.2)then
!      write(200+iam,*)'bottom doloop ',j,' nstepat ',nstepat call flush(200+iam)
!      endif

      return
      end subroutine doloop


c-------------------------------------------------------------------------
      subroutine init_larc
      use raqmschem_pmgrid_mod, only : masterproc,beglat,endlat,iam
      use raqmschem_pmgrid_mod, only : bbedirf,bbedira,calnex,arcpac
      use raqmschem_pmgrid_mod, only : nmonHTAP,nr=>plat,nc,nl=>plev
      use raqmschem_pmgrid_mod, only : prefi
!      use mpimcomm, only : areasum,dydeg
      use raqmschemcomm_mod
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
#ifdef DOMPI
      use mpishorthand, only : mpicom,mpireal,mpir4
#endif
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
      implicit none

#include <choosechem.h>
#include <comm_3d_chem>
      include 'combined_constants'
      include 'combined_model_logic'

      character*40 xjunk
      real(CHEM_KIND_R8) :: xlat,surf_new,sa1,sa2
      common /aero_24/ xlat(ijg2),surf_new(31,37,2),
     1                 sa1(31,37,2),sa2(31,37,2)
      real(CHEM_KIND_R8) :: troplat,troptheta
      common /tropopause/ troplat(64),troptheta(64)
      real(CHEM_KIND_R8) :: xlsig,fixer,climacetone
      common /fixer_24/ xlsig(inl),fixer(3,ijg2,inl)
     1                 ,climacetone(19,17,2)
      real(CHEM_KIND_R8) :: xn2o_s,f11_s,f12_s,ccl4_s,ch3cl_s,xmtcfm_s,
     1                ch3br_s,h1301_s,h1211_s
      common /chembc/ xn2o_s,f11_s,f12_s,ccl4_s,ch3cl_s,xmtcfm_s,
c     1                ch3br_s,h1301_s,h1211_s,ch4_s
     1                ch3br_s,h1301_s,h1211_s
c
c
      integer nc1deg,nr1deg
      parameter(nc1deg=360,nr1deg=181)
      integer nlondum,nlatdum,nlndum
      real *4 dumlon(nc1deg,nr1deg),dumlat(nc1deg,nr1deg)
c
      integer*4 d1,d2,d3,kmain
      real *4 xdum(img),ydum(ijg2),plev_dum(nl_clim_dyn),s1,s2
      integer*2 parameter_id
      character *200 fileemis
      integer rclocal,ii,jj,ll,ndum,idum,nlm,l,j,i
      logical exist
      real(CHEM_KIND_R8) :: sina,weight,csj,rcsj,alat,rmg

444   continue
!      write(200+iam,*)'top initlarc '
!      call flush(200+iam)
      if(masterproc)then
        inquire(file=trim(prefi)//'trop_het_zmean_T32.d',exist=exist)
        if(exist)then
!        write(6,*)'file exists ',trim(prefi)//'trop_het_zmean_T32.d'
!        call flush(6)
        open(610,file=trim(prefi)//'trop_het_zmean_T32.d',
     1  status='old',form = 'unformatted')
!        write(6,*)'kind trhet ',kind(trhet),' shape ',shape(trhet)
!        write(6,*)'bounds ',lbound(trhet),ubound(trhet)
!        call flush(6)
        read (610) trhet
        close (610)
        else
          write(6,*)'does not exist'
          call flush(6)
        endif
      endif


      if (.not.usefastj) then
        if(masterproc)then
c  read in kawa photolysis routines lookup tables
          write(6,*)'reading kawa lookup tables'
          call flush(6)
          open(unit=305,file=trim(prefi)//
     1    'Kawa/xsect_47_79.dat',
     1    form='unformatted',status='old')
          read(305) xtab
          close(305)
        endif       
#ifdef DOMPI
!        call mpibcast(xtab,nlam*nts*nxsec,mpir4,0,mpicom)
        call raqmschem_comm_all_bcast(xtab,rc=rclocal)
#endif
        call rdstab
      endif

c  read file containing climatologies for h2o, co, h2
      if(masterproc)then
        write(6,*)'reading climatologies'
        call flush(6)
          write(6,*)'prefi',trim(prefi)
          call flush(6)
        open(243,file=trim(prefi)//'fixer_trop_24_64_3',status='old')
        read(243,34) fixer
   34   format(8e10.3)
        close(243)
     

c       read acetone climatology
        open(37,file=trim(prefi)//'acetone_dist.dat',status='old')
        read(37,950) xjunk
        read(37,950) xjunk
c        read Jan mean
        do ll=1,17                                   ! alt index 0-16km (1-17)
          read(37,*) (climacetone(jj,ll,1),jj=1,19)  ! lat index S->N (1-19)
        enddo
        read(37,950) xjunk
        read(37,950) xjunk
        read(37,950) xjunk
c       read Jul mean
        do ll=1,17                                   ! alt index 0-16km (1-17)
          read(37,*) (climacetone(jj,ll,2),jj=1,19)  ! lat index S->N (1-19)
        enddo
        close(37)

c  read in wmo91 sulfate aerosol surface area distribution and enhanced
c  surface areas (SA1, SA2) from Debra Weisenstein (AER)
        write(6,*)'reading surface area climotologies'
        call flush(6)
        open(37,file=trim(prefi)//'surface_area.dat',status='old')
        do ii=1,21
          read(37,950) xjunk
        end do
        read(37,*) ((surf_new(ii,jj,1),ii=1,31),jj=1,37)
        read(37,950) xjunk
        read(37,950) xjunk
        read(37,*) ((surf_new(ii,jj,2),ii=1,31),jj=1,37)
        read(37,950) xjunk
        read(37,950) xjunk
        read(37,*) ((sa1(ii,jj,1),ii=1,31),jj=1,37)
        read(37,950) xjunk
        read(37,950) xjunk
        read(37,*) ((sa1(ii,jj,2),ii=1,31),jj=1,37)
        read(37,950) xjunk
        read(37,950) xjunk
        read(37,*) ((sa2(ii,jj,1),ii=1,31),jj=1,37)
        read(37,950) xjunk
        read(37,950) xjunk
        read(37,*) ((sa2(ii,jj,2),ii=1,31),jj=1,37)
        close(37)
  950   format(a40)

c  read mixing ratios of long-lived species at the surface
       print *,'reading surface mixing ratios'
c  read mixing ratios of long-lived species at the surface
       open(38,file=trim(prefi)//'tracer_surf_conc_1990.dat',
     +     status='old')
        read(38,950) xjunk
        read(38,*) xn2o_s
        read(38,*) f11_s
        read(38,*) f12_s
        read(38,*) ccl4_s
        read(38,*) ch3cl_s
        read(38,*) xmtcfm_s
        read(38,*) ch3br_s
        read(38,*) h1301_s
        read(38,*) h1211_s
c JAA 20050907 CH4 is now specified by latitude and month and is read
c  in with emission/deposition files, see ch4clim below
c      read(38,*) ch4_s
        close(38)
      endif
!      write(6,*)'do raqmschem_comm_all_bcast'
!      call flush(6)
#ifdef DOMPI
!      call mpibcast(fixer,3*ijg2*inl,mpireal,0,mpicom)
!      write(6,*)'raqmschem_comm_all'
!      call flush(6)
      call raqmschem_comm_all_bcast(fixer,rc=rclocal)
!      write(6,*)'did it rclocal'
!      call flush(6)
!      call mpibcast(climacetone,19*17*2,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(climacetone,rc=rclocal)
!      call mpibcast(surf_new,31*37*2,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(surf_new,rc=rclocal)
!      call mpibcast(xn2o_s,1,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(xn2o_s,rc=rclocal)
      call raqmschem_comm_all_bcast(f11_s,rc=rclocal)
      call raqmschem_comm_all_bcast(f12_s,rc=rclocal)
      call raqmschem_comm_all_bcast(ccl4_s,rc=rclocal)
      call raqmschem_comm_all_bcast(ch3cl_s,rc=rclocal)
      call raqmschem_comm_all_bcast(xmtcfm_s,rc=rclocal)
      call raqmschem_comm_all_bcast(ch3br_s,rc=rclocal)
      call raqmschem_comm_all_bcast(h1301_s,rc=rclocal)
      call raqmschem_comm_all_bcast(h1211_s,rc=rclocal)
!      call mpibcast(f11_s,1,mpireal,0,mpicom)
!      call mpibcast(f12_s,1,mpireal,0,mpicom)
!      call mpibcast(ccl4_s,1,mpireal,0,mpicom)
!      call mpibcast(ch3cl_s,1,mpireal,0,mpicom)
      !call mpibcast(xmtcfm_s,1,mpireal,0,mpicom)
!      call mpibcast(ch3br_s,1,mpireal,0,mpicom)
!      call mpibcast(h1301_s,1,mpireal,0,mpicom)
!      call mpibcast(h1211_s,1,mpireal,0,mpicom)
#endif
!      write(6,*)'masterproc',masterproc
!      call flush(6)
      if(masterproc)then
        write(6,*)'xn2o_s',xn2o_s,f11_s,f12_s,ccl4_s,ch3cl_s,xmtcfm_s,ch3br_s,
     +  h1301_s,h1211_s
        call flush(6)
      

c  read tropopause potential temperature definition
        write(6,*)'reading tropopause potential temperature',kind(troplat)
        call flush(6)
        open(39,file=trim(prefi)//'tropopause_theta.dat',status='old')
        read(39,*) ndum
        read(39,*) troplat
        read(39,*) troptheta
!        do j=1,64
!          write(6,*)'troplat',j,troplat(j),'troptheta',troptheta(j)
        !end do
        close(39)

c  read heterogeneous reaction rate constant for N2O5 hydrolysis in the
c  troposphere (zonal averaged data from Dentner and Crutzen)
        write(6,*)'reading rate constant for N2O5'
        call flush(6)
        inquire(file=trim(prefi)//'trop_het_zmean_T32.d',exist=exist)
        if(exist)then
!        write(6,*)'file exists ',trim(prefi)//'trop_het_zmean_T32.d'
!        call flush(6)
        open(610,file=trim(prefi)//'trop_het_zmean_T32.d',
     1  status='old',form = 'unformatted')
!        write(6,*)'kind trhet ',kind(trhet),' shape ',shape(trhet)
!        write(6,*)'bounds ',lbound(trhet),ubound(trhet)
!        call flush(6)
        read (610) trhet
        close (610)
        else
          write(6,*)'does not exist'
          call flush(6)
        endif
      endif
!      write(6,*)'now bcast'
!      call flush(6)
#ifdef DOMPI
!      call mpibcast(troplat,64,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(troplat,rc=rclocal)
!      call mpibcast(troptheta,64,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(troptheta,rc=rclocal)
!      call mpibcast(trhet,10*ijg2,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(trhet,rc=rclocal)
#endif
!      write(300+iam,*)'masterproc',masterproc
!      call flush(300+iam)


10005 if(masterproc)then
        close (500)
!        write(300+iam,*)'open ',trim(prefi)//'upd_nox_air_hsrp_'//emdeplabel//'.ALL'
!        call flush(300+iam)

        open (500,file=trim(prefi)//
     1 'upd_nox_air_hsrp_'//emdeplabel//'.ALL',
     1  status='old',form='unformatted')
        read(500)nlondum,nlatdum,nlndum
!        write(300+iam,*)'nlon',nlondum,nlatdum,nlndum
!        call flush(300+iam)
c       ajl note sairnox_p are international standard atmoshpere
c       pressures corresponding to 1 km height flight levels
c       0,1,2 .. 18 km  (1013 to 75mb)
!        write(6,*)'sairnox_p',shape(sairnox_p)
!        call flush(6)
        sairnox_p=0.0
        read(500)dumlon,dumlat,sairnox_p
!        write(300+iam,*)'dumlon',maxval(dumlon),maxval(dumlat)
!        call flush(300+iam)
!        write(300+iam,*)kind(sairnox_p),'read sairnox_p',sairnox_p
!        call flush(300+iam)
!!        write(6,*)'sairnox_p',sairnox_p
        close(500)
!        write(6,*)'read lnox_frac'
!        call flush(6)
        open(500,file=trim(prefi)//'lnox_frac_profile.dat',
     1  status='old',form='formatted')
        do idum=1,4
          read(500,*)
        enddo
        read(500,*) (slnox_profile(idum,1),idum=1,nlnlgt)
        print*,'Midlat Contin:'
        print*,(slnox_profile(idum,1),idum=1,nlnlgt)
        read(500,*)
        read(500,*) (slnox_profile(idum,2),idum=1,nlnlgt)
        print*,'Tropic Marine:'
        print*,(slnox_profile(idum,2),idum=1,nlnlgt)
        read(500,*)
        read(500,*) (slnox_profile(idum,3),idum=1,nlnlgt)
        print*,'Tropic Contin:'
        print*,(slnox_profile(idum,3),idum=1,nlnlgt)
        close(500)
      endif
#ifdef DOMPI
!      call mpibcast(slnox_profile,nlnlgt*3,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(slnox_profile,rc=rclocal)
      call raqmschem_comm_all_bcast(sairnox_p,rc=rclocal)
#endif
!      write(300+iam,*)'sairnox_p',sairnox_p,'nlnair',nlnair
!      call flush(300+iam)


      nlm=inl-1
      if(masterproc)then
c       read and process sigma levels
        open(201,file=trim(prefi)//'sig_34_short',status='old')
        read(201,*) (sigmah(l),l=1,nlm)
        close(201)
      endif
#ifdef DOMPI
!      call mpibcast(sigmah,nlm,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(sigmah,rc=rclocal)
#endif
      sigmah(inl)=1.
      s1=0.
      do l=1,nlm
        s2=sigmah(l)
        sigma(l)=(s1+s2)*.5
        s1=s2
      end do
      sigma(inl)=.5*(1.+sigmah(nlm))
c  define log(sigma) on inverted altitude scale (xlsig is ground to top)
      do l=1,inl
        xlsig(l)=-log(sigma(inl-l+1)*1000.)
      end do
      if(masterproc)then
        print *,xlsig

c       read and process larc grid gaussian latitudes
        open(208,file=trim(prefi)//'gaulat_64',status='old')
        do j=1,ijg
          read(208,102) sina,weight
102       format(2e30.23)
          csj=1.-sina*sina
          rcsj=1./csj
          alat=atan(sina*sqrt(rcsj))*57.29578
c         ylarc goes from NP to SP
          ylarc(j)=alat
          ylarc(ijg2-j+1)=-alat
          xlat(j)=ylarc(j)
          xlat(ijg2-j+1)=ylarc(ijg2-j+1)
        end do
        close(208)
      endif
#ifdef DOMPI
!      call mpibcast(ylarc,ijg2,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(ylarc,rc=rclocal)
!      call mpibcast(xlat,ijg2,mpireal,0,mpicom)
      call raqmschem_comm_all_bcast(xlat,rc=rclocal)
#endif

c  define larc grid longitude array
      rmg=1./float(img)
      do i=1,img
        xlarc(i)=360.*float(i-1)*rmg
      enddo
c
c read LaRC-gridded ground info: draw still used in deposition param.
c
c Not used - now using online drag (draw) calculation JAA 07/08/03
c       print *,'reading LaRC air-sea interaction file'
c      open(unit=15,
c     & file='IMPACT_data/Dsigma_567361',
c     & form='unformatted',status='old')
c       read(15)cd,wmag,draw,gbphis
c       close(15)
c
c read LaRC-gridded monthly-mean pressure data for upper-strat column O3
c
!      call mpibarrier(mpicom)
      if(masterproc)then
        print*,'reading LaRC monthly-mean pressure files'
        open(unit=15,file=trim(prefi)//'IMPACT_data/Zpress_yr12_Mar',
     & form='unformatted',status='old')
        read(15)kmain
        read(15)parameter_id
        read(15)d1,d2,d3
        read(15)xdum,ydum,plev_dum
        read(15)clim_zpress
        close(15)
        open(unit=15,file=trim(prefi)//'IMPACT_data/Tpress_yr12_Mar',
     & form='unformatted',status='old')
        read(15)kmain
        read(15)parameter_id
        read(15)d1,d2,d3
        read(15)xdum,ydum,plev_dum
        read(15)clim_tpress
        close(15)
        open(unit=15,file=trim(prefi)//'IMPACT_data/oxpress_yr12_Mar',
     & form='unformatted',status='old')
        read(15)kmain
        read(15)parameter_id
        read(15)d1,d2,d3
        read(15)xdum,ydum,clim_plev
        read(15)clim_oxpress
        close(15)
      endif
#ifdef DOMPI
!      call mpibcast(clim_zpress,img*ijg2*nl_clim_dyn,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(clim_zpress,rc=rclocal)
!      call mpibcast(clim_tpress,img*ijg2*nl_clim_dyn,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(clim_tpress,rc=rclocal)
!      call mpibcast(clim_oxpress,img*ijg2*nl_clim,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(clim_oxpress,rc=rclocal)
!      call mpibcast(clim_plev,nl_clim,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(clim_plev,rc=rclocal)
#endif

!      write(200+iam,*)'ajl bottom initlarc '
!      call flush(200+iam)
      return
      end subroutine init_larc
c-------------------------------------------------------------------------
      subroutine clim_spec(julday,jlat,dens,mpr,ktf,rlwi)
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
      use raqmschem_pmgrid_mod, only : masterproc,iam
      include 'comm_defs_24t.chem'
      common /chemnfo_24/ imon
      common /blk4_24/ r(ireac),hr(ihreac),phoj(iphot),dn(inum),
     1              zeta,pur,sad,prodloss(13)
      real(CHEM_KIND_R8) :: xlat,surf_new,sa1,sa2
      common /aero_24/ xlat(ijg2),surf_new(31,37,2),
     1                 sa1(31,37,2),sa2(31,37,2)
      real*8 xlsig,fixer,climacetone
      common /fixer_24/ xlsig(inl),fixer(3,ijg2,inl)
     1                 ,climacetone(19,17,2)
c
C$OMP THREADPRIVATE (/blk4_24/)
c
      dimension jip(3)
      data jip/13,20,36/
      save jip
c
c The fixer array is stored (species,latitude,altitude)
c  species = H2O, H2 (CO  is now treated as a transported species)
c  latitude = Summer pole to Winter pole
c  altitude = bottom (ground) to top
c
!      write(200+iam,*)'clim_spec ',jlat,'ireac',ireac,ihreac,iphot,inum,zeta,pur
!      call flush(200+iam)
      dday=julday
      daytst=mod(dday,365.)
      frag=abs(182.5-daytst)/182.5
      jrev=(ijg2+1)-jlat
c  interpolate specified quantities to actual pressure (account for
c  topography)                                          Rich 3/26/96
c  first find sigma level just below current pressure
      xlogpur=-log(pur)
      do i=inl,1,-1
        if (xlsig(i) .le. xlogpur) goto 100
      end do
c  bottom level is less than xlsig(1) if we fall out of do loop
100   nlev=i
      if ((nlev .ne. inl) .and. (nlev .ne. 0)) then
        slope=(xlogpur-xlsig(nlev))/(xlsig(nlev+1)-xlsig(nlev))
      end if
c
c  pick up specified quantities  (only h2o and h2)
      do i=1,3,2
        ink=jip(i)
c  interpolate fixer to proper level                    Rich 3/26/96
        if ((nlev .ne. inl) .and. (nlev .ne. 0)) then
          fix_jrev=(fixer(i,jrev,nlev+1)-fixer(i,jrev,nlev))*slope +
     1             fixer(i,jrev,nlev)
          fix_jlat=(fixer(i,jlat,nlev+1)-fixer(i,jlat,nlev))*slope +
     1             fixer(i,jlat,nlev)
          apex=frag*(fix_jrev-fix_jlat)+fix_jlat
        else
          apex=frag*(fixer(i,jrev,mpr)-fixer(i,jlat,mpr))+
     1         fixer(i,jlat,mpr)
        end if	
        dn(ink)=dens*apex
      end do
c  --------------------------------------------------------
c  sulfate aerosol chemistry parameterization from WMO 91
c  using Debra Weisenstein's files.  ** Use SA0 profile for HSCTs. **
c  *** Background aerosols used in this model ***
c  first find index points in sulfate surface area array
      zstar=16.*alog10(1000./pur)
c  altitude (pressure) index (0-60 km: 1-31)
      i_aero=nint(zstar/2.)+1
      i_aero=max(min(i_aero,31),1)
c  latitude index (SP-NP: 1-37)
      j_aero=int(xlat(jlat)/5.)+19
c  month index
      k_aero=int(imon/7.)+1
c
c Avoid calculations in troposphere (H2SO4 fit is invalid, but no
c  aerosol there anyway)
cjjj Test to ensure we have aerosol down to the tropopause
cjjj      if (i_aero .le. 6) then
      if (i_aero .le. 5) then
        sad=0.
      else if (i_aero .eq. 6) then
        sad=surf_new(i_aero+1,j_aero,k_aero)*0.5
      else
        sad=surf_new(i_aero,j_aero,k_aero)
      endif
c
c  specify acetone concentration (troposphere only)
      dn(64) = 0.
c
c  20060203 JAA from MNN
c  Use all values in climatology, even above climatological tropopause,
c   to prevent "zero" values near tropopause (climo decays to 0 at 16km)
c 20060203      if (ktf .eq. 1) then
        ll = nint(zstar)+1              ! alt index 0-16km (1-17)
        if(ll<1)then
           write(6,*)'ll',ll,'zstar',zstar
           call flush(6)
        endif
        if (ll .le. 17) then
          jj = nint(xlat(jlat)/10.)+10  ! lat index S->N (1-19)
c   linear interp in time between mid-Jan and mid-Jul tabled values
c    Fixed dayoff bug 20060203 JAA
c          dayoff = daytst+16.
          dayoff = daytst-16.
          tf = abs(1.-(dayoff/181.))
          if (tf .gt. 1.) tf=2.-tf
          dn(64) = tf*climacetone(jj,ll,1)+(1.-tf)*climacetone(jj,ll,2)
          dn(64) = dn(64)*1.e-12*dens   ! climacetone units pptv
        endif
c 20060203      endif
c
cccmnn 0703  specify tropospheric methanol concentration
ccc         free trop 600 ppt, marine bl 900 pt, and cont. bl 2000 ppt
      dn(84) = 0.0
      if (ktf .eq. 1) then
        if (pur .gt. 800.) then
          if (rlwi .lt. 1.5) then
            dn(84) = 9.0e-10*dens
          else if (rlwi .lt. 2.5) then
            dn(84) = 2.0e-9*dens
          else
            dn(84) = 6.0e-10*dens
          endif
        else
          dn(84) = 6.0e-10*dens
        endif
      endif
c
      return
      end
c-------------------------------------------------------------------------
      subroutine jul_gmt (istime,julday,gmt)
      use raqmschem_pmgrid_mod, only : masterproc
c
c*****************************************************************c
c  From the start time, istime, compute the Julian day and GMT    c
c  istime is in "yymmddhh" format.                                c
c                                             Duncan 8/13/99      c
c*****************************************************************c
c
      DIMENSION MNO(12)        
      integer*4 istime,ictime,inch,incd,incm,incy
      integer*4 yy1,mm1,dd1,hh1,yy2,mm2,dd2,hh2
      data (MNO(I),I=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/    
      save mno
C                                                                       
      yy1=istime/1000000
      mm1=istime/10000-yy1*100
      dd1=istime/100-yy1*10000-mm1*100
      hh1=istime-yy1*1000000-mm1*10000-dd1*100
c
      if(masterproc)then
      print *,'Jul_gmt yy1,mm1,dd1,hh1 ',yy1,mm1,dd1,hh1
      endif
      julday = dd1
      if (mm1 .gt. 1) then
        do m = 1, mm1-1
          julday = julday + mno(m)
        end do
      endif
      gmt = 1.0 * hh1
      if(masterproc)then
      print *,'Jul_gmt julday,gmt ',julday,gmt
      endif
      return
      end
c-------------------------------------------------------------------------
#ifndef RDALBDOSLOW
      subroutine rd_albedo(jday)
      use raqmschem_pmgrid_mod, only : masterproc
      use raqmschem_pmgrid_mod, only : prefi,iam
      use raqmschemcomm_mod
#ifdef DOMPI
      use mpishorthand
#endif
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
c
c  read in monthly albedo data and interpolate to Julian day 
c 
#include <choosechem.h>
#include <comm_3d_chem>
      integer jmon(13),jday,jdayold
      character*3 amon(13)
c     ajl
c     real zalbmm1(ialb,jalb),zalbm(ialb,jalb),zalb0(ialb,jalb)
      real *4 zalbmm1(ialb,jalb),zalbm(ialb,jalb,12),zalb0(ialb,jalb)
      data jmon/1,32,60,91,121,152,182,213,244,274,305,335,366/       
      data amon/'jan','feb','mar','apr','may','jun',
     &          'jul','aug','sep','oct','nov','dec','jan'/
      save jmon,amon
      logical first
      integer rclocal
      data first/.true./,jdayold/0/
      save first,zalbm,jdayold
#ifdef PRINTJULDAY
      if(masterproc)print *,'top fast rd_albedo',jday,'first',first
#endif
!      write(6,*)'fast rdalb'
!      call flush(6)
      if(first)then
#ifdef NOTEMPARRAY
        if(masterproc)print *,'NOTEMPARRAY'
#else
        if(masterproc)print *,'TEMPARRAY'
#endif
        if(masterproc)then
          OPEN(10,FILE=trim(prefi)//'/GLA_albedo.f77',
     &    FORM='UNFORMATTED',STATUS='OLD')
c
c         read in GLA grid
c
          read(10)xalb,yalb
c         ajl 6/24/2008 yalb in data file off by 2 at north pole
          do j=1,jalb
            yalb(j)=-90.+180./float(jalb-1)*float(j-1)
c           write(6,*)'yalb',j,yalb(j)
          end do
c         print *,xalb,yalb
c
c         read in GLA land/water/ice flag
c
          read(10)zlwi
c         print *,zlwi(1,1)
c
c         read in albedo
c
c         read(10)zalbmm1
          do  m=1,12 
            read(10)zalbm(:,:,m)
          end do
          close(10)
        endif
#ifdef DOMPI
!        call mpibcast(xalb,ialb,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(xalb,rc=rclocal)
      if(rclocal/=0)then
        write(6,*)'error xalb ',rclocal
        call flush(6)
        call killit('bcast error')
!        stop
      endif
!        call mpibcast(yalb,jalb,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(yalb,rc=rclocal)
      if(rclocal/=0)then
        write(6,*)'error yalb ',rclocal
        call flush(6)
        call killit('bcast2 2 errror')
!        stop
      endif
!        call mpibcast(zlwi,ialb*jalb,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(zlwi,rc=rclocal)
      if(rclocal/=0)then
        write(6,*)'error zlwi',rclocal
        call flush(6)
        call killit('bcst3 error')
!        stop
      endif
!        call mpibcast(zalbm,ialb*jalb*12,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(zalbm,rc=rclocal)
      if(rclocal/=0)then
         write(6,*)'error zalbm',rclocal
        call flush(6)
        call killit('bcast4 error')
!         stop
      endif
#endif
        first=.false.
      endif
!      write(200+iam,*)'jday',jday,'jdayold',jdayold
!      call flush(200+iam)
      if(jday.eq.jdayold)return
#ifdef DOMPI
!      call mpibcasti0(jday)
      call raqmschem_comm_all_bcast(jday,rc=rclocal)
      if(rclocal/=0)then
        write(6,*)'error jday',rclocal
        call flush(6)
        call killit('bcsat 5 error')
!        stop
       endif
!       write(200+iam,*)' bcast jday'
!       call flush(200+iam)
#endif
c     print *,zalbmm1(1,1)
c
c map current month to january for save
c
c     do 5 i=1,ialb
c       do 5 j=1,jalb         
c        zalb0(i,j)=zalbmm1(i,j)
c5    continue
      do 30 m=2,12 
c          read(10)zalbm(:,:,m)
c         print *,zalbm(1,1)
c
c find proper month and interpolate
c
        if(jday.ge.jmon(m-1).and.jday.lt.jmon(m))then
          tscale=1.*(jday-jmon(m))/(jmon(m-1)-jmon(m))
          do i=1,ialb
            do j=1,jalb         
             zalb(i,j)=zalbm(i,j,m)+tscale*(zalbm(i,j,m-1)-zalbm(i,j,m))
            end do
          end do
          if(masterproc)then
            if(jday.ne.jdayold)then
              write(6,*)'Interpolated albedo between ',
     &        amon(m-1),' & ',amon(m),' for Julian day ',
     &        jday,' scale=',tscale
              call flush(6)
            endif
          endif
          goto 50
        endif
c
c map current month to previous month
c
c        do 20 i=1,ialb
c          do 20 j=1,jalb         
c           zalbmm1(i,j)=zalbm(i,j)
c20      continue
30      continue
c
c handle dec to january interpolation
c
        if(jday.ge.jmon(12)) then
          tscale=1.*(jday-jmon(13))/(jmon(12)-jmon(13))
          do i=1,ialb
            do j=1,jalb         
c            zalb(i,j)=zalb0(i,j)+tscale*(zalbm(i,j)-zalb0(i,j))
              zalb(i,j)=zalbm(i,j,1)+tscale*(zalbm(i,j,12)-zalbm(i,j,1))
            end do
          end do
          if(masterproc)then
            if(jday.ne.jdayold)then
              print *,'Interpolated albedo between ',
     &        amon(12),' & ',amon(1),' for Julian day ',
     &        jday,' scale=',tscale
            endif
          endif
        endif
50    continue
c      endif
#ifdef DOMPI
c      call mpibcast(xalb,ialb,mpir4,0,mpicom)
c      call mpibcast(yalb,jalb,mpir4,0,mpicom)
c      call mpibcast(zlwi,ialb*jalb,mpir4,0,mpicom)
c      call mpibcast(zalb,ialb*jalb,mpir4,0,mpicom)
#endif
      jdayold=jday
!      write(200+iam,*)'bottom new rd_aLBEDO',jdayold
!      CALL FLUSH(200+IAM)
!      write(6,*)'bottom new rd_albedo ',jday
!      call flush(6)
      return
      end
#else
      subroutine rd_albedo(jday)
      use raqmschem_pmgrid_mod, only : masterproc
      use raqmschem_pmgrid_mod, only : prefi
      use raqmschemcomm_mod
#ifdef DOMPI
      use mpishorthand
#endif
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
c
c  read in monthly albedo data and interpolate to Julian day 
c 
#include <choosechem.h>
#include <comm_3d_chem>
      integer jmon(13),rclocal
      character*3 amon(13)
c     ajl
c     real zalbmm1(ialb,jalb),zalbm(ialb,jalb),zalb0(ialb,jalb)
      real *4 zalbmm1(ialb,jalb),zalbm(ialb,jalb),zalb0(ialb,jalb)
      data jmon/1,32,60,91,121,152,182,213,244,274,305,335,366/       
      data amon/'jan','feb','mar','apr','may','jun',
     &          'jul','aug','sep','oct','nov','dec','jan'/
      save jmon,amon
      write(6,*)'new alb read'
      call flush(6)
      if(masterproc)then
#ifdef PRINTJULDAY
        print *,'top rd_alb'
#endif
        OPEN(10,FILE=trim(prefi)//'/GLA_albedo.f77',
     &  FORM='UNFORMATTED',STATUS='OLD')
c
c read in GLA grid
c
        read(10)xalb,yalb
c       ajl 6/24/2008 yalb in data file off by 2 at north pole
        do j=1,jalb
          yalb(j)=-90.+180./float(jalb-1)*float(j-1)
c         write(6,*)'yalb',j,yalb(j)
        end do
c       print *,xalb,yalb
c
c read in GLA land/water/ice flag
c
        read(10)zlwi
c       print *,zlwi(1,1)
c
c read in albedo
c
        read(10)zalbmm1
c       print *,zalbmm1(1,1)
c
c map current month to january for save
c
        do 5 i=1,ialb
          do 5 j=1,jalb         
           zalb0(i,j)=zalbmm1(i,j)
5       continue
        do 30 m=2,12 
          read(10)zalbm
c         print *,zalbm(1,1)
c
c find proper month and interpolate
c
          if(jday.ge.jmon(m-1).and.jday.lt.jmon(m))then
            tscale=1.*(jday-jmon(m))/(jmon(m-1)-jmon(m))
            do 10 i=1,ialb
            do 10 j=1,jalb         
              zalb(i,j)=zalbm(i,j)+tscale*(zalbmm1(i,j)-zalbm(i,j))
10          continue
            print *,'Interpolated albedo between ',
     &      amon(m-1),' & ',amon(m),' for Julian day ',
     &      jday,' scale=',tscale
            goto 50
          endif
c
c map current month to previous month
c
          do 20 i=1,ialb
            do 20 j=1,jalb         
             zalbmm1(i,j)=zalbm(i,j)
20        continue
30      continue
c
c handle dec to january interpolation
c
        if(jday.ge.jmon(12)) then
          tscale=1.*(jday-jmon(13))/(jmon(12)-jmon(13))
          do 40 i=1,ialb
          do 40 j=1,jalb         
            zalb(i,j)=zalb0(i,j)+tscale*(zalbm(i,j)-zalb0(i,j))
40        continue
          print *,'Interpolated albedo between ',
     &    amon(12),' & ',amon(1),' for Julian day ',
     &    jday,' scale=',tscale
        endif
50    continue
      close(10)
      endif
#ifdef DOMPI
!      call mpibcast(xalb,ialb,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(xalb,rc=rclocal)
!      call mpibcast(yalb,jalb,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(yalb,rc=rclocal)
!      call mpibcast(zlwi,ialb*jalb,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(zlwi,rc=rclocal)
!      call mpibcast(zalb,ialb*jalb,mpir4,0,mpicom)
      call raqmschem_comm_all_bcast(zalb,rc=rclocal)
#endif
      return
      end
#endif
c-------------------------------------------------------------------------
      subroutine getdiag
c 
c    get physical and dynamical diagnostics at the trajectory locations 
c    including solar zenith angle, 2D albedo, land-water-ice flag,
c    overhead column ozone   
      use raqmschem_pmgrid_mod, only : beglat,endlat,nc,endj,begj,ibeg,tile
      use raqmschemcomm_mod
      use raqmschem_pmgrid_mod, only : masterproc,iam
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
c
#include <choosechem.h>
#include <comm_3d_chem>
c
c albedo over water zenith angle correction arrays
c
      real alb_scale(10),zen_bin(10)
      data alb_scale/1.,1.079,1.197,1.329,1.513,1.750,2.118,
     &               2.671,3.526,4.395/
      data zen_bin/25.84,36.87,45.57,53.13,60.00,66.42,72.54,
     &             78.46,84.26,90.00/
      save alb_scale,zen_bin
!      write(6,*)iam,'getdiag nc',nc,'beglat',beglat,endlat,'lb',lbound(xgrid)
!      call flush(6)
!      write(200+iam,*)tile,'i1grd',maxval(i1grd),minval(i1grd)
!      write(200+iam,*)tile,'j1grd',maxval(j1grd),minval(j1grd)
!      call flush(200+iam)
#ifdef DIAGFV3
      if(tile.eq.1)then
        write(6,*)'nc',nc,'beglat',beglat,endlat,ibeg,'i1grd',i1grd(1,beglat),j1grd(1,beglat)
        call flush(6)
      endif
#endif
 
      do 80 ii=1,nc
       do 80 jj=beglat,endlat
c
        xp=xgrid(ii,jj)
        yp=ygrid(ii,jj)
c
c compute solar zenith angle for each grid point
c
        doy = 1. * julday
        if(masterproc.and.jj.eq.beglat.and.ii.eq.1)then
            write(6,*)iam,'ajl doy',doy,'solar gmt',gmt
        endif
        call zenith_angle(doy,gmt,yp,xp,chi,diurn_scal)
        szagrd(ii,jj) = chi
!        if(tile.eq.3.and.i1grd(ii,jj).eq.272.and.j1grd(ii,jj).eq.54)then
#ifdef DIAGFV3
        if(tile.eq.1.and.ibeg+ii-1.eq.1.and.jj.eq.1)then
!        if(iam>=9.and.iam<=12)then
           write(6,*)'chi',chi,ii,jj,'xp',xp,yp,'doy',doy,'solar gmt',gmt,'diurn',diurn_scal
!           write(200+iam,*)'chi',chi,ii,jj,'xp',xp,yp,'doy',doy,'solar gmt',gmt,'diurn',diurn_scal
           !call flush(200+iam)
           call flush(6)
         endif
#endif

!        if(ii.eq.30.and.jj.eq.15)then
!           write(6,*)iam,'szagrd',ii,jj,szagrd(ii,jj),'xgrid',xgrid(ii,jj),ygrid(ii,jj),'beglat',beglat,endlat
!           !write(6,*)iam,'doy',doy,'gmt',gmt,'nc',nc
!           write(6,*)iam,'xp',xp,yp,'chi',chi,diurn_scal
!           call flush(6)
!        !endif
!        if(ii.eq.1.and.jj.eq.1)then
!          write(6,*)iam,'szagrd',szagrd(ii,jj),'xgrid',xgrid(ii,jj),ygrid(ii,jj),'beglat',beglat,endlat
!          write(6,*)iam,'doy',doy,'gmt',gmt,'nc',nc
!          call flush(6)
!        endif
        diurngrd(ii,jj) = diurn_scal
c
c interpolate 2D albedo and land/water/ice flag separately
c due to fixed GLA grid note: reversal of i,j order
c
        do 70 i=1,ialb
        ip1=i+1
        if(i.eq.ialb)ip1=1
        xlon=xalb(i) 
        xlonp1=xalb(ip1)
        if(i.eq.ialb)xlonp1=360.+xalb(ip1)
        if((xp.ge.xlon).and.(xp.lt.xlonp1))then
        xscale=(xp-xlon)/(xlonp1-xlon)
c
        do 60 j=1,jalb-1
        jp1=j+1
        xlat0=yalb(j)
        xlatp1=yalb(jp1)
c
        if((yp.ge.xlat0).and.(yp.le.xlatp1))then
        yscale=(yp-xlat0)/(xlatp1-xlat0)
c
        zalbj=zalb(i,j)+
     &          xscale*(zalb(ip1,j)-zalb(i,j))
        zalbjp=zalb(i,jp1)
     +        +xscale*(zalb(ip1,jp1)-zalb(i,jp1))
        albgrd(ii,jj)=zalbj+yscale*(zalbjp-zalbj)
c
c JAA 20080225 zlwigrd has already been built in subroutine getlwi
c  and is passed through comm_3d_chem
c
c        zlwij=zlwi(i,j)+
c     &          xscale*(zlwi(ip1,j)-zlwi(i,j))
c        zlwijp=zlwi(i,jp1)
c     +        +xscale*(zlwi(ip1,jp1)-zlwi(i,jp1))
c        zlwigrd(ii,jj)=zlwij+yscale*(zlwijp-zlwij)
c
c determine whether over water (not sea ice) and adjust albedo to 
c account for zenith angle dependence
c
c Note: scale by normalized directional clear sky Albedo over ocean from
c       NASA Ref. Pub. 1184 Angular Radiation models for Earth-Atmosphere
c       System Table 7  
c
        if(zlwigrd(ii,jj).eq. 1.0 .and. albgrd(ii,jj).lt. 0.1) then 
        do 55 nbin = 1, 9
          if(szagrd(ii,jj).ge.zen_bin(nbin).and.
     &       szagrd(ii,jj).lt.zen_bin(nbin+1))then     
            zscale=(szagrd(ii,jj)-zen_bin(nbin))/
     &             (zen_bin(nbin+1)-zen_bin(nbin))
            ascale=alb_scale(nbin)+
     &             zscale*(alb_scale(nbin+1)-alb_scale(nbin))
            albgrd(ii,jj)=albgrd(ii,jj)*ascale
          endif
55      continue
        if(szagrd(ii,jj).ge.zen_bin(10))
     &  albgrd(ii,jj)=albgrd(ii,jj)*alb_scale(10)
        endif
        goto 75
        endif
60      continue
        endif
70      continue
75      continue
80    continue
      return
      end  
c-------------------------------------------------------------------------
      subroutine zenith_angle(doy,gmt,rlat,rlon,chi,diurn_scal)
c*****************************************************************c
c  Compute solar zenith angle for each parcel                     c
c                                                                 c
c INPUT:                                                          c
c  doy   : julian day (INPUT)                                     c
c  gmt   : greenwich mean time (INPUT)                            c
c  tlat  : true latitude  (INPUT)                                 c
c  tlon  : true longitude (INPUT)                                 c
c OUTPUT:                                                         c
c  chi   : zenith angle (degrees) at position (tlat,tlon)         c
c*****************************************************************c
c
c  zenith angle computation
      real *8 gmt
      pi=3.14159265
      dtor=pi/180.
      earinc=23.5
c
      sinlat=sin(rlat*dtor)
      coslat=sqrt(1.-sinlat**2)
      sinlon=sin(rlon*dtor)
      coslon=cos(rlon*dtor)
      soya=(doy-81.25)*pi/182.5           ! day angle
      soha=2.*pi*(gmt-12.)/24.            ! hour angle
      soha=-soha
      sininc=sin(earinc*dtor)
      sindec=sininc*sin(soya)
      cosdec= sqrt(1.-sindec**2)
      coszen=cos(soha)*coslon+sin(soha)*sinlon
      coszen=coszen*cosdec*coslat
      coszen=sindec*sinlat+coszen
!      if(abs(rlat-53.8186)<.2.and.abs(rlon-272.2149)<.2)then
!        write(6,*)'coszen',coszen,'doy',doy,'soha',soha,'cosdec',cosdec
!        write(6,*)'coslat',coslat,'gmt',gmt
!      endif
      coszen=min(max(coszen,-1.),1.)
      chi = acos(coszen)
      chi = chi/dtor

CCCMNN 8/2003
CCCMNN   SAVE PARAMETERS NEEDED FOR IMPOSING DIURNAL VARIATION
CCC      IN A SOURCE (E.G. ISOPRENE). SOURCE DURING DAYLIGHT IS GIVEN BY
CCC       S = SB *(PI/(SLSD*THDY+CLCD*SIN(THDY)) * COSZEN
CCC       WHERE, S=SOURCE AT ZENITH ANGLE (coszen=cos(sza))
CCC       SB=DIURNALLY AVERAGED SOURCE FROM GEIA DATABASE,
CCC       PI = 3.141592654, SLSD=SIN(LAT)*SIN(DEC), CLCD=COS(LAT)*COS(DEC),
CCC       THDY = HOUR ANGLE AT SUNRISE OR SUNSET (LENGTH OF HALF DAY),
CCC       AND COSZEN=COS(ZENITH ANGLE)
CCC
      diurn_scal = 0.0
      if (coszen .gt. 0.) then
        slsd = sinlat*sindec
        clcd = coslat*cosdec
        if(clcd.eq.0.0)then
#ifdef PRINTCLCD
          write(6,*)'clcd',clcd,'coslat',coslat,'sinlat',sinlat,
     * 'sindec',sindec,'cosdec',cosdec,'rlat',rlat
#endif
          if(slsd.gt.0.)then
            thdy=pi
            diurn_scal=coszen/slsd
          elseif(slsd.lt.0.0)then
            thdy=0.0
            diurn_scal=0.0
          else
            thdy=0.0
            diurn_scal=0.0
          endif
        else
          hdy = -slsd/clcd
          if (hdy .le. -1.0) then
            thdy = pi
            diurn_scal = coszen/slsd
          elseif (hdy .ge. 1.0) then
            thdy = 0.0
            diurn_scal = 0.0
          else
            thdy = acos(hdy)
            diurn_scal = pi * coszen/(slsd*thdy + clcd*sin(thdy))
          endif
        endif
      endif

      return
      end
c-------------------------------------------------------------------------

      subroutine getgrid_uwhyb
#include <choosechem.h>
      use raqmschemcomm_mod
      use raqmschem_pmgrid_mod, only : beglat,endlat,nr=>plat,nc
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
#include <comm_3d_chem>
      dx0=360./nc
      dy0=180./(nr-1)
      do i=1,nc
       do j=1,nr
        xgrid(i,j)=(i-1)*dx0
        ygrid(i,j)=-90+(j-1)*dy0
       enddo
      enddo
      return
      end
c-------------------------------------------------------------------------
      subroutine calc_colo3_col(sza,plog,tc,svo3,zgeo,ptop,lbot,
     +           klon,jlat,svcolo3,svcolo3v,chapno,ierr)

c  calculate and store the column of o3 overhead and slant column densities
c  (e.g., overhead column at each model level at this lon/lat)
      use raqmschem_pmgrid_mod, only : nl=>plev
      use raqmschemcomm_mod
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8

#include <choosechem.h>
#include <comm_3d_chem>
      real zgeo(nl),plog(nl),tc(nl),svo3(nl),svcolo3(nl),svcolo3v(nl)

      zeta = cos(sza)

c  Initialize column ozone at top level
c  If upper boundary is 1 mb or lower, use IMPACT climatology to initialize 
      if (ptop .ge. 1.) then
        call larco3clim(klon,jlat,zeta,sza,tc(1),ptop,svo3(1),zgeo(1),
     +                  colo3,colo3_vert,chapno,ierr)
      else
        call colo3top(zeta,sza,tc(1),svo3(1),svo3(2),zgeo(1),zgeo(2),
     +                colo3,colo3_vert,chapno,ierr)
      endif
      svcolo3(1) = colo3
      svcolo3v(1) = colo3_vert

c  Sum and store column O3 at all other levels
      do ip=2,lbot
        ia=ip-1
        call getchapo3(zeta,sza,tc(ip),zgeo(ip),svo3(ia),svo3(ip),
     +                 zgeo(ia),zgeo(ip),chapo3,h_o3,ierr)
        if(hydrostatic) then
          deltao3 = 1463.5891*(plog(ip)-plog(ia))*
     1     (svo3(ia)*tc(ia)+svo3(ip)*tc(ip))
        else
crbp      account for pressure increasing with height (non-hydrostatic)
          deltao3 = 1.e5*(svo3(ia)+svo3(ip))*(zgeo(ia)-zgeo(ip))/2.
        endif
        colo3 = colo3+deltao3*chapo3
        colo3_vert = colo3_vert+deltao3
        svcolo3(ip) = colo3
        svcolo3v(ip) = colo3_vert
      enddo

c  Set below-ground values (if any) equal to surface values
      ip=lbot+1
      do while (ip .le. nl)
        svcolo3(ip) = colo3
        svcolo3v(ip) = colo3_vert
        ip=ip+1
      enddo

      return
      end
c-------------------------------------------------------------------------
      subroutine colo3top(zeta,sza,temp,o3up,o3lo,zgup,zglo,
     +                    colo3,colo3_vert,chapno,ierr)
c Initialize column ozone at top of domain (valid in upper stratosphere)

cjjj crude accounting for diffs between ox and o3 at stratopause
      o3up=o3up*0.6

      call getchapo3(zeta,sza,temp,zgup,o3up,o3lo,
     +               zgup,zglo,chapo3,o3scalht,ierr)
      colo3 = o3up * o3scalht * chapo3
      colo3_vert = o3up * o3scalht
      chapno = chapo3

      return
      end
c-------------------------------------------------------------------------
      subroutine larco3clim(klon,jlat,zeta,sza,temp,ptop,o3top,ztop,
     +                      colo3,colo3_vert,chapno,ierr)
      use raqmschemcomm_mod
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
#include <choosechem.h>
#include <comm_3d_chem>

      real tc(nl_clim),zgeo(nl_clim),svo3(nl_clim),plog(nl_clim)
      parameter(zer_min=1.)

      loff=nl_clim_dyn-nl_clim
      do l=1,nl_clim
       if (clim_plev(l) .lt. ptop) then
        tc(l)=clim_tpress(klon,jlat,l+loff)
        zgeo(l)=clim_zpress(klon,jlat,l+loff)/1000.
        dens=7.2431122e+18*clim_plev(l)/tc(l)
        svo3(l)=max(clim_oxpress(klon,jlat,l)*dens,zer_min)
        plog(l)=alog(1000.*clim_plev(l))
       else
        tc(l)=temp
        zgeo(l)=ztop
        svo3(l)=o3top
        plog(l)=alog(1000.*ptop)
        lbot=l
        goto 99
       endif
      enddo
      print*,klon,jlat,temp,sza,ptop,o3top,ztop
      print*,nl_clim,clim_plev
      stop 'problem in larco3clim'

   99 continue

      call colo3top(zeta,sza,tc(1),svo3(1),svo3(2),zgeo(1),zgeo(2),
     +              colo3,colo3_vert,chapno,ierr)

c Sum rest of O3 column between top of LaRC grid and top of current grid
      do ip=2,lbot
        ia=ip-1
        call getchapo3(zeta,sza,temp,zgeo(ip),svo3(ia),svo3(ip),
     +                 zgeo(ia),zgeo(ip),chapo3,h_o3,ierr)
        deltao3 = 1463.5891*(plog(ip)-plog(ia))*
     1    (svo3(ia)*tc(ia)+svo3(ip)*tc(ip))
        colo3 = colo3+deltao3*chapo3
        colo3_vert = colo3_vert+deltao3
      enddo

      return
      end
c-------------------------------------------------------------------------
      subroutine getchapo3(zeta,sza,temp,zloc,o3up,o3lo,
     +                     zgup,zglo,chapo3,h_o3,ierr)

      scalht=287.*temp*100./9.81
      z1=amax1(zloc,0.)

      if (o3lo .ne. o3up) then
        h_o3=-(zglo-zgup)/log(o3lo/o3up)*1.0e5
        h_o3 = max(h_o3,0.1*scalht)                     !JAA 04/24/97
      else
        h_o3=1.0e9
        ierr = 1
      end if

c  if sza le 60 degrees then use zeta, otherwise calculate chapman function
      if (sza .gt. 1.0472) then
        chapo3=chap(z1,sza,abs(h_o3))
      else
        chapo3=1./zeta
      end if

      return
      end
c------------------------------------------------------------------------
      subroutine chemdraw(lat,draw)
      use raqmschem_pmgrid_mod
      use raqmschemcomm_mod
      use chem_types_mod, only : CHEM_KIND_R4,CHEM_KIND_R8
      implicit none
      integer lat,i
      real(CHEM_KIND_R8) :: draw(nc)
      real(CHEM_KIND_R8) ::  wmags,wmag,dts,cdrag,wmagc
!      write(200+iam,*)'top chemdraw lat',lat,'nc',nc,'lbound',lbound(cmdrag),ubound(cmdrag),'kind',kind(cmdrag)
!      call flush(200+iam)
      do i=1,nc
!        write(200+iam,*)'cmdrag ',i,lat,cmdrag(i,lat)
!        call flush(200+iam)
        cdrag=cmdrag(i,lat)
        wmags=sfcwind(i,lat)**2
        wmag=max(sfcwind(i,lat),0.001)
        dts=tskin(i,lat)-tgrd(i,lat,1)
        wmagc=wmags*wmag 
        if(dts.ge.0.)draw(i)=cdrag*(wmag+sqrt(dts))
        if(dts.lt.0.0)draw(i)=cdrag*wmagc/(wmags-7.*dts)
!        if(iam.eq.0.and.i.eq.1)then
!          write(6,*)'cdrag',cdrag,'wmags',wmags,'wmag',wmag
!          write(6,*)'tskin',tskin(i,lat),'tgrd',tgrd(i,lat,i)
!          write(6,*)'dts',dts,'draw',draw(i)
!        endif
      end do
        
      return
      end
c------------------------------------------------------------------------

      subroutine caldm1(iyr,imon,iday)
c  Subtract 1 day from iday, adjust month and year if necessary,
c   and return updated values

      dimension ndays(12)
c days in month
      data ndays/31,28,31,30,31,30,31,31,30,31,30,31/

      ist = iday-1
      if (ist .eq. 0) then
        imnth = imon-1
        if(imnth.eq.0)then
           iyr = iyr-1
           imon = 12
           iday = ndays(12)
        else
           imon = imnth
           iday = ndays(imnth)
c handle leap year
           if (imon.eq.2) then
             if (mod(iyr,4).ne.0) then
               iday = 28
             else if (mod(iyr,400).eq.0) then
               iday = 29
             else if (mod(iyr,100).eq.0) then
               iday = 28
             else
               iday = 29
             endif
           endif
c
        endif
      else
        iday = ist
      endif

      return
      end

c------------------------------------------------------------------------
#if 0
!     this routine is no longer used zlwigrd defined elsewhere
      subroutine getlwi
#if 0
      use raqmschem_pmgrid_mod
#ifdef LANDMASKONCE
      use raqmschemcomm_mod, only : zlwigrd,tskin,seaicethk,landlake
#else
      use raqmschemcomm_mod, only : zlwigrd,tskin,seaicethk
#endif

#ifdef DOMPI
      use mpishorthand, only : mpicom,mpiint
#endif
      use comsrf, only : ts,oro,sicthk
      
!     #include <comm_link_chem>
#ifndef LANDMASKONCE
      integer landlake(nc,plat)
      save landlake
#endif
      logical first
      save first
      data first/.true./

      integer i,j
!     need to tell where lakes are and set to land for chemisty
#ifdef LANDMASKONCE
      if(first)then
        allocate (landlake(nc,plat))
        first=.false.
        if(masterproc)then
          print *,'landmaskonce'
#ifdef RES11
          call hojofn('landmask.gfsmask.1x1.hoj',40,'OLD')
#else
          call hojofn('landmask.gfsmask.2x2.hoj',40,'OLD')
#endif
          call hojori(nc*plat,landlake,'LANDLAKE')
          call hojocl
        endif
#ifdef DOMPI
        call mpibcast(landlake,nc*plat,mpiint,0,mpicom)
#endif
      endif
#else
      if(masterproc)then
#ifdef RES11
        call hojofn('landmask.gfsmask.1x1.hoj',40,'OLD')
#else
        call hojofn('landmask.gfsmask.2x2.hoj',40,'OLD')
#endif
        call hojori(nc*plat,landlake,'LANDLAKE')
        call hojocl
      endif
#ifdef DOMPI
      call mpibcast(landlake,nc*plat,mpiint,0,mpicom)
#endif
#endif

c
c land/water/ice flag:
c RAQMS: 1=ocean, 2=land, 3=ice;    oro: 0=W, 1=L, 2=I
c
      do j=beglat,endlat
      do i=1,nc
        zlwigrd(i,j) = 1.
        if (nint(oro(i,j)).eq.1) then
          zlwigrd(i,j)=2.
        endif
        if (nint(oro(i,j)).eq.2) zlwigrd(i,j)=3.
        if (nint(oro(i,j)).eq.0)then
!         marked ocean make sure not a lake
!#ifdef NEWLANDLAKE
          if(landlake(i,j).ne.0)then
            zlwigrd(i,j)=2.
!            print *,'make ocean lake',i,j,'lon=',
!     *      float(i-1)*2.,-90.+float(j-1)*2.,'landlake',
!     *landlake(i,j),'oro',oro(i,j)
          endif
!#endif
        endif
        tskin(i,j) = ts(i,j)
        seaicethk(i,j) = sicthk(i,j)
      end do
      end do
#endif
      return
      end
#endif
      subroutine openrsync(bbdir,filename)
      character *(*) bbdir,filename
      character *256 cmd,exdir*240
      print *,'opernrsync bbdffir',bbdir

      if(bbdir(1:5).eq.'rsync')then
        exdir=' '
        call getenv('EXDIR',exdir)
        if(exdir.ne.' ')then
          cmd=trim(bbdir)//trim(filename)//' '//trim(exdir)//
     *    '/'//trim(filename)
        else
          cmd=trim(bbdir)//trim(filename)//' '//trim(filename)
        endif
        write(6,*)'cmd openrsync'
        write(6,*)cmd
        write(6,*)'filename',filename
        call flush(6)
        call system(cmd)
        open(500,file=filename,status='old',form='unformatted')
      else
        write(6,*)'open file',trim(bbdir)//trim(filename)
        call flush(6)
        open(500,file=trim(bbdir)//trim(filename),status='old',form=
     *  'unformatted')
!        write(6,*)'did open file',trim(bbdir)//trim(filename)
!        call flush(6)
      endif
  
      return
      end subroutine openrsync
      logical function lopenrsync(bbdir,filename)
      character *(*) bbdir,filename
      character *256 cmd,exdir*240
      logical exist
      print *,'lopernrsync bbdffir',bbdir

      if(bbdir(1:5).eq.'rsync')then
        exdir=' '
        call getenv('EXDIR',exdir)
        if(exdir.ne.' ')then
          cmd=trim(bbdir)//trim(filename)//' '//trim(exdir)//
     *    '/'//trim(filename)
        else
          cmd=trim(bbdir)//trim(filename)//' '//trim(filename)
        endif
        write(6,*)'cmd openrsync'
        write(6,*)cmd
        write(6,*)'filename',filename
        call flush(6)
        call system(cmd)
        inquire(file=filename,exist=exist)
        lopenrsync=exist
        if(exist)then
          print *,'fire file exists a ',filename
          open(500,file=filename,status='old',form='unformatted')
        endif
      else
        inquire(file=filename,exist=exist)
!       ajl add for iris
!       if given a path which is not scratch move to scratch
!       1/29/2015
        if(.not.exist)then
          exdir=' '
          call getenv('EXDIR',exdir)
          if(exdir.ne.' ')then
            cmd='cp -p '//trim(bbdir)//trim(filename)//' '//trim(exdir)//
     *      '/'//trim(filename)
          else
            cmd='cp -p '//trim(bbdir)//trim(filename)//' '//trim(filename)
          endif
          write(6,*)'cmd cp fire'
          write(6,*)cmd
          write(6,*)'filename',filename
          call flush(6)
          call system(cmd)
          inquire(file=filename,exist=exist)
          lopenrsync=exist
          if(exist)then
            print *,'fire file exists a ',filename
          else
            write(6,*)'file ',trim(filename),' not there try sleep 5'
            call system('sleep 5')
            inquire(file=filename,exist=exist)
            lopenrsync=exist
            write(6,*)'exist',exist
          endif
        endif
        lopenrsync=exist
        if(exist)then
          print *,'fire file exists b ',filename
          open(500,file=trim(bbdir)//trim(filename),status='old',form=
     *    'unformatted')
        endif
      endif
      return
      end function lopenrsync
!      subroutine firenetcdf(bbco_n,'bb_N_co_fuel.//cyyyymmdd//'.nc')
      subroutine firenetcdf(bbco,filefire)
      character *(*)filefire,fireoutdir*200
      real*4 bbco(360,181)
#if 0
      include 'netcdf.inc'
      integer ncid,dimlat,dimlon,dim(2),ifire,ierr
      write(6,*)'top firenetcdf ',filefire
      call flush(6)
      fireoutdir=' '
      call getenv('FIREOUTDIR',fireoutdir)
      write(6,*)'fireoutdif',trim(fireoutdir)
      call flush(6)

      ierr=nf_create(trim(fireoutdir)//filefire,0,ncid)
      ierr=nf_def_dim(ncid,'lat',181,dimlat)
      ierr=nf_def_dim(ncid,'lon',360,dimlon)
      dim(1)=dimlon
      dim(2)=dimlat
      write(6,*)'dim',dim
      call flush(6)
      ierr=nf_def_var(ncid,'bbco',nf_float,2,dim,ifire)
      ierr=nf_enddef(ncid)
      ierr=nf_put_var_real(ncid,ifire,bbco)
      write(6,*)trim(filefire),' bbco ',maxval(bbco)
      call flush(6)
      ierr=nf_close(ncid)
#endif
      return
      end subroutine firenetcdf
      subroutine readr4unit
      return
      end subroutine readr4unit

