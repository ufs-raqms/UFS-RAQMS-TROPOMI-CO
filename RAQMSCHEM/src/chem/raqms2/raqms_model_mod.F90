module raqms_model_mod
#include <options.h>
contains
  subroutine raqms_model_advance(rc)
  use mpi
  use machine,only : kind_phys
  use funcphys
  use chem_types_mod
  use raqmschem_state_mod
  use raqmschem_config_mod
  use raqmschem_data_mod
  use raqmschem_model_mod
  use chem_comm_mod, only : chem_comm_get
!  use raqmschem_model_mod, only : chem_model_get => raqmschem_model_get
  use raqmschem_io_mod, only : chem_io_write,chem_io_read
  use raqms_mod
  use calcdxdy, only : initdxdy
!  use calcvort_mod, only :calcvort
  use calcvort_dgrid_mod, only : calcabsvort
  use field_manager_mod, only : find_field_index
  use raqmschem_iodata_mod, only : raqmschem_chem_write,raqmschem_chem_write_var,raqmschem_chem_write_o3vmr
  use raqmschem_iodata_mod, only : raqmschem_chem_write_aodfcst
  use raqmschem_iodata_mod, only : raqmschem_gsi_write
  use raqmschem_iodata_mod, only : raqmschem_chem_write_var_2d,raqmschem_chem_write_debug
  use raqmschem_species_mod
  use chem_const_mod, only : grvity
  use raqmschemcomm_mod, only : allocatechemcommperm,pblht,slmsk2d,ustar,cmdrag
  !use raqmschemcomm_mod, only : szgrd,tskin,sfcwind,tgrd,lconvcld,ccthk,zlwigrd,thgrd
  use raqmschemcomm_mod, only : tskin,sfcwind,tgrd,lconvcld,ccthk,zlwigrd,thgrd,ltb
  use raqmschemcomm_mod, only : pgrd,ugrd,vgrd,spgrd,zsurf,zgrd,dpmgrd
  use raqmschemcomm_mod, only : allocatechemcommtemp,latgrd,longrd,xgrid,ygrid,areagrd
  use raqmschemcomm_mod, only : taucldfrc,i1grd,j1grd,dist1grd
!  use raqmschemcomm_mod, only : !  ztropraq,ptropraq,ktrop,ptropgsi,ztropgsi,ktropgsi
!  use raqmschemcomm_mod, only : taucldfrc_liq,taucldfrc_ice,i1grd,j1grd,dist1grd
!  use raqmschemcomm_mod, only : gsiinc
!  use raqmschemcomm_mod, only : incMLS3d,incOMI3d
!  use raqmschemcomm_mod, only : incNUCAPSco3d,incTROPOMIco3d
!  use raqmschemcomm_mod, only : percov,intco
!  use raqmschemcomm_mod, only : pergsiv,colgsiinc
  use raqmschemcomm_mod, only : deallocatechemcommtemp,o3dep_save,navg,noydep_save,codep_save,o3vmr_inst
  use raqmschem_pmgrid_mod, only : tile,ibeg,iam,iprn,jprn,iamprn,iprnin,jprnin,pecount,kprnin,tileprn,cdateat
  use raqmschem_pmgrid_mod, only : nstepat,fracday,timestep,forecaSt_hr,nhtuw,cdateat,gsivar
  use raqmschem_pmgrid_mod, only : gsivar,numgsicem,gsicem,loccem
  use raqmschem_pmgrid_mod, only : raqms_localIOflag,ompsgsi,cdatestart
  use chem_const_mod, only : omegx2,degrad
  use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
  use chem_raqms_mod, only : chem_pass_state
  use raqmschem_io_mod,only : inquire_file_var
  use chemgsimod
  implicit none
  type(chem_state_type),  pointer :: stateIn, stateOut
  integer rc,localrc,de,is,ie,js,je,ni,nl,i,decount,advancecount
  integer ihs,ihe,jhs,jhe
  integer its,ite,jts,jte,ims,ime,jms,jme,mype
  integer advnaceCount, julday, mm, tz,yy,dd,h,m,s
  real(CHEM_KIND_R8) :: dts,f
  real(CHEM_KIND_R8), dimension(:,:), pointer :: lat, lon 
  type(raqmschem_config_type), pointer :: config
  type(chem_data_type),   pointer :: data
  integer icall,isave
  integer j,k,ii,jj,dims(3),ktop,kbot,dims2(2)
  real(CHEM_KIND_R4),allocatable,dimension(:,:,:) :: varin,traceout
  real(CHEM_KIND_R4) :: af,fhmax
  REAL,    PARAMETER :: airmw    = 28.97
  real,    parameter :: mwo3=47.98
  real(CHEM_KIND_R4),allocatable,dimension(:,:) :: dz
!  real(CHEM_KIND_R4),allocatable,dimension(:,:,:) :: dzaod,coinc,no2inc
  real(CHEM_KIND_R4),allocatable,dimension(:,:,:) :: dzaod
!  real(CHEM_KIND_R4),allocatable,dimension(:,:,:,:) :: tempgsiinc
!  real(CHEM_KIND_R4),allocatable,dimension(:,:,:) :: tempincMLS,tempincOMI
!  real(CHEM_KIND_R4),allocatable,dimension(:,:,:) :: tempincNUCAPS,tempincTROPOMI
!  real(CHEM_KIND_R8) :: conew,no2new
!  real(CHEM_KIND_R8) :: gsinew,gsioxnew,hrout
  character *10 ctest,cstep,filename*100,cdateout,gsiremapsh*256,cmd*256
  character *256 remap05,crealtimerun*10
  character *256 gsiscratch,base_data_disk,pathdata,scratch_out
  character *10 ciat,cjat,ckat,ctileat,ctile*1
  real(CHEM_KIND_R8) :: wallwrite,wallwritecum,rtime1,rtimed0,rtimed1,rtimecum
  real(CHEM_KIND_R8) :: wallstart
  character *256 GSIDATAPATH,GSINAME,BASE_PATH,GSIINC_PATH,GSI_INQ_SCRIPT
  integer idate
  save wallwritecum,rtimecum,rtimed1,wallstart,gsiremapsh,gsiscratch
  save GSIDATAPATH,GSINAME,GSI_INQ_SCRIPT
  save fhmax
#define NOCLDTAUZERO
#ifdef CLDTAUZERO
  logical lcldtauzero
  data lcldtauzero/.false./
  save lcldtauzero
#endif
!  logical :: doall=.true.
  logical :: doall=.false.
  logical :: doaod

  real(chem_kind_r4) :: afsfc
  integer shape4(4),nchem
#include <comcdate.h>
  save icall
  data icall/0/
  character *256 nodelist
  logical first,useraqmso3vmr,exist,existt,existi,dogsiname,doaodfcst
  save first,useraqmso3vmr,dogsiname,doaodfcst
  data first/.true./
  data useraqmso3vmr/.false./,doaodfcst/.false./
  integer nsecsave,ihrsave,modelcomm,nsecgsi,ihrgsi
  data nsecsave/10800/
  save nsecsave,nsecgsi,base_data_disk
  character*10 cenv,expname*50,datascratch*256
!  integer numgsicem,loccem(100),countsleep,n,numneg
  integer countsleep,n,numneg
  integer*8 fsize
  integer*4 ibegdate,ienddate
  save ibegdate,ienddate
  real*4 latnegno2,lonnegno2,no2bef,no2inc
  real*4 negno2
!  character *20 gsicem(100)
  character *10 slurm_job_id
  character *10 date,time,zone
  integer ivalues(8),ihr,imod,idtsec

#define ALLOWSKIPCHEM
#ifdef ALLOWSKIPCHEM
  logical lskipchem
  data lskipchem/.false./
  save lskipchem
#endif
!  call raqmschem_model_get(deCount=deCount, tile=tile,rc=localrc)
  call raqmschem_model_get(deCount=deCount, tile=tile,localIoflag=raqms_localIOflag,rc=localrc) ! new
!  if (decount<1)then
!    write(200+iam,*)'decount zero return'
!  else
!    write(200+iam,*)'decount ',decount
!  endif
!  flush(200+iam)
  if (deCount < 1) return
  call chem_comm_get(localpe=mype,pecount=pecount)
!  write(200+iam,*)'localpe ',mype,'pecount',pecount,'tile',tile,'localIO',raqms_localIOflag
!  flush(200+iam)
  call raqmschem_model_get(modelcomm=modelcomm,rc=rc) ! ajl add new
  call chem_model_clock_get(advanceCount=advanceCount, dts=dts, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s,tz=tz, julday=julday, rc=localrc)
  advancecount=advancecount+1 ! ajl 4/6/2022 to sync plumerise calls
  timestep=dts
    icall=icall+1
  cenv=' '
  call getenv('AODFCST',cenv)
  if(cenv=='YES')then
    doaodfcst=.true.
    if(iam==0)then
      write(6,*)'doaodfcst'
    endif
  endif

  base_path=' '
  call getenv('BASE_PATH',base_path)
  if(first)then
  if(iam.eq.0.and.tile.eq.1)then
    write(6,*)'base_path',base_path
    write(6,*)'CHEMGSI DIR'
  endif
  endif
  slurm_job_id=' '
  call getenv('SLURM_JOB_ID',slurm_job_id)
  if(.not.first)then
    rtimecum=rtimecum+mpi_wtime()-rtimed1
      if(iam.eq.0.and.mod(nstepat,nhtuw).eq.0)then
      write(6,*)'walldyncum',rtimecum,' ',rtimecum/60.,' min '
    endif
  endif
  if(first)then
    if(iam.eq.0.and.tile.eq.1)then
      write(6,*)'SLURM_JOB_ID',trim(slurm_job_id)
      call flush(6)
    endif
    wallstart=mpi_wtime()
    base_data_disk=' '
    call getenv('base_data_disk',base_data_disk) 
    cdatestart=' '
    call getenv('CDATE',cdatestart)
      if(tile.eq.1.and.mype.eq.0)then
        write(6,*)'CDATESTART',cdatestart
        call flush(6)
      endif
    scratch_out=' '
    call getenv('SCRATCH_OUT',scratch_out)
    remap05=' '
    call getenv('REMAP05',remap05)
    crealtimerun=' '
    call getenv('REALTIMERUN',crealtimerun)

    gsiremapsh=' '
    call getenv('GSIREMAPSH',gsiremapsh)
    if(iam.eq.0.and.tile.eq.1)then
      write(6,*)'base_data_disk',trim(base_data_disk)
      write(6,*)'gsiremapsh',gsiremapsh
      call flush(6)
    endif

    nodelist=' '
  
    call getenv('SLURM_JOB_NODELIST',nodelist)
    if(nodelist/=' ')then
      if(tile.eq.1.and.mype.eq.0)then
         write(6,*)'nodelist',trim(nodelist)
         call flush(6)
      endif
    endif
    expname=' '
    call getenv('EXPNAME',expname)
    cenv=' '
    call getenv('FHMAX',cenv)
    if(cenv/=' ')then
      read(cenv,*)fhmax
      if(mype.eq.0)then
        write(6,*)'fhmax',fhmax
      endif
    endif
    cenv=' '
    call getenv('RAQMSO3MR',cenv)
    if(cenv.eq.'YES')then
      useraqmso3vmr=.true.
      if(mype.eq.0)then
        write(6,*)'use raqmso3vmr in o3phys'
        call flush(6)
      endif
    endif
    wallwritecum=0.0
    call setchempointers(mype)
    if(mype==0)then
    write(6,*)'GSIREMAPSH',trim(GSIREMAPSH)
    call flush(6)
    endif
    if(GSIREMAPSH/=' ' )then
!      if(mod(nint(dts*float(icall)),nsecgsi).eq.0.or.icall<=1)then
!      call gsivar(numgsicem,gsicem,loccem)
      call gsivar
     endif
#ifdef ALLOWSKIPCHEM
    ctest=' '
    call getenv('SKIPCHEM',ctest)
    if(ctest.ne.' ')then
      if(ctest.eq.'YES')then
        lskipchem=.true.
      endif
    endif
#endif
    ctest=' '
    call getenv('HRSAVE',ctest)
    if(ctest.ne.' ')then
      read(ctest,*)ihrsave
!      write(6,*)'ihrsave',ihrsave
      nsecsave=float(ihrsave)*3600.
    endif
    ctest=' '
    call getenv('DHRGSI',ctest)
    if(ctest.ne.' ')then
      read(ctest,*)ihrgsi
      nsecgsi=float(ihrgsi)*3600.
      if(iam==0)then
        write(6,*)'nsecgsi',nsecgsi,'DHRGSI',ctest
      endif
    endif
#ifdef CLDTAUZERO
    ctest=' '
    call getenv('CLDTAUZERO',ctest)
    if(ctest.eq.'YES')then
      lcldtauzero=.true.
      if(mype.eq.0)then
        write(6,*)'test CLDTAUZERO'
      endif
    endif
#endif
  endif
#ifdef ALLOWSKIPCHEM
  if(lskipchem)then
    first=.false.
    rtimed1=mpi_wtime()
    return
  endif
#endif
  do de = 0, deCount-1
      call raqmschem_model_get(de=de, config=config, data=data, &
        stateIn=stateIn, stateOut=stateOut, rc=localrc)
    shape4=shape(statein%tr3d)
    nchem=shape4(4)
    if(nchem<5)then
      if(mype.eq.0)then
        write(6,*)'no chem'
        call flush(6)
      endif
      first=.false.
      return
    endif
     call chem_model_domain_get(de=de, ids=is, ide=ie, jds=js, jde=je, ni=ni, nl=nl, &
     lon=lon, lat=lat, rc=localrc,its=its,ite=ite,jts=jts,jte=jte,ims=ims,ime=ime,jms=jms,jme=jme) 
!    ids,ide full limits of this subregion may not start at 1
!    its,ite full limits for this tile 1 to 96,192,384 etc
!    ims,ime limits start with index 1
     ihs=max(is-1,its)
     jhs=max(js-1,jts)
     ihe=min(ie+1,ite)
     jhe=min(je+1,jte)
!     write(6,*)'lat',maxval(lat),minval(lat)
!     write(6,*)'lon',maxval(lon),minval(lon)
     if(first)then
!      can calculate griddx,and griddy first time
       call initdxdy(is,ie,js,je,lon,lat,ihs,ihe,jhs,jhe,data%griddx,data%griddy,de,mype)
     endif
    if(.not.allocated(data%fcor))then
      allocate(data%fcor(is:ie,js:je))
      data%fcor=omegx2*sin(lat*degrad)
    endif
    if(.not.allocated(data%vort))then
      allocate(data%vort(is:ie,js:je,nl))
    endif
    call calcabsvort(statein%vort,statein%prl3d,statein%tk3d,lat,is,ie,js,je,nl,de,mype)
          
        
    
    if(first)then
      latgrd=lat
      longrd=lon
      xgrid=lon
      ygrid=lat
      i1grd=nint(lon)
      j1grd=nint(lat)
      dist1grd=sqrt((i1grd-lon)**2.+(j1grd-lat)**2.)
      data%area=statein%area
      areagrd=statein%area
!     will have to read in since ph3d 1 is zero
!      szgrd=statein%ph3d(:,:,1)/grvity
      slmsk2d=statein%slmsk2d
    endif
    call allocatechemcommtemp
    pblht=statein%pb2d
    spgrd=statein%pr3d(:,:,1)*.01 ! mb
    
    
!    zsurf=statein%phl3d(:,:,1)/9.80616 ! meters
    ustar=statein%us2d
    cmdrag=statein%cmdrag
    tskin=statein%ts2d
    tgrd=statein%tk3d
    ugrd=statein%us3d
    vgrd=statein%vs3d
    pgrd=statein%prl3d*.01 ! mb
    if(pgrd(1,js,1) > pgrd(1,js,nl))then
      do k=1,nl
        ltb(k)=nl-k+1
      end do
    else
      do k=1,nl
        ltb(k)=k
      end do
    endif
    thgrd(:,:,:)=tgrd(:,:,:)*(1000./pgrd(:,:,:))**.286
    do k=1,nl
!     zgrd is dashes
!      zgrd(:,:,k)=statein%phl3d(:,:,k)/9.80616+zsurf(:,:)
      zgrd(:,:,k)=statein%phl3d(:,:,k)/grvity+zsurf(:,:)
      dpmgrd(:,:,k)=(statein%pr3d(:,:,k)-statein%pr3d(:,:,k+1))*.01 ! mb
    end do
    call trop(is,ie,js,je,nl,spgrd,pgrd,zgrd,tgrd)
#ifdef DOPROF
    if(first)then
      if(mype.eq.9.and.tile.eq.3)then
        open(101,file='/odyssey/proxy/lenzen/prof.fv3.110.147',form='formatted')
        open(102,file='/odyssey/proxy/lenzen/prof.fv3.107.150',form='formatted')
        write(101,*)spgrd(18,18),zsurf(18,18)
        write(102,*)spgrd(15,20),zsurf(15,20)
!        do k=1,nl
!          write(101,*)k,pgrd(18,18,k),zgrd(18,18,k)
!          write(102,*)k,pgrd(15,20,k),zgrd(15,20,k)
!        end do
        call flush(101)
        call flush(102)
      endif
      if(tile.eq.3)then
        do j=js,je
          do i=1,48
            if(nint(latgrd(i,j)).eq.56.and.nint(longrd(i,j)).eq.109)then
              write(6,*)'latgrd',i,j,latgrd(i,j),longrd(i,j)
              call flush(6)
              write(6,*)'sp',i,j,spgrd(i,j),'zsurf',zsurf(i,j)
              do k=1,nl
                write(6,*)i,j,'p',pgrd(i,j,k),'zgrd',zgrd(i,j,k)
              end do
            endif
            if(nint(latgrd(i,j)).eq.59 .and.nint(longrd(i,j)).eq.106)then
              write(6,*)'latgrd',i,j,latgrd(i,j),longrd(i,j)
              call flush(6)
              write(6,*)'sp',spgrd(i,j),'zsurf',zsurf(i,j)
              do k=1,nl
                write(6,*)i,j,'p',pgrd(i,j,k),'zgrd',zgrd(i,j,k)
              end do
            endif
          end do
        end do
      endif
    endif
#endif
!    szgrd(:,:)=statein%phl3d(:,:,1)/9.80616
   
    do j=js,je
      jj=j-js+1
      do i=1,ie-is+1
        ii=is+i-1          
        afsfc=(statein%phl3d(i,jj,1)-statein%ph3d(i,jj,1))/(statein%ph3d(i,jj,2)-statein%ph3d(i,jj,1))
        sfcwind(i,j)=sqrt( (statein%us3d(i,jj,1)*(1.+afsfc)+statein%us3d(i,jj,2)*afsfc)**2+ &
        (statein%vs3d(i,jj,1)*(1.+afsfc)+statein%vs3d(i,jj,2)*afsfc)**2 )
        ktop=statein%cktop_kbot_cnv(i,jj)/1000
        kbot=statein%cktop_kbot_cnv(i,jj)-ktop*1000
        if(ktop>0)then
!          ccthk(i,j)=(statein%ph3d(i,jj,ktop)-statein%ph3d(i,jj,kbot))/9.80616/1000.
          ccthk(i,j)=(statein%ph3d(i,jj,ktop)-statein%ph3d(i,jj,kbot))/grvity/1000.
!          write(6,*)'ccthk=',i,j,ccthk(i,j),'ktop',ktop,kbot
!          flush(6)
          lconvcld(i,j)=ktop
        else
          ccthk(i,j)=0.0
          lconvcld(i,j)=0
        endif
        if(statein%slmsk2d(i,jj)<.001)then
!         ocean
          if(statein%ph3d(i,jj,1)<.0001)then
            zlwigrd(i,j)=1. ! ocean
          else
            zlwigrd(i,j)=2. ! make lake land
          endif
        else
          if(nint(statein%slmsk2d(i,jj)).eq.2)then
            zlwigrd(i,j)=3.
          else
!         land
            zlwigrd(i,j)=2.
          endif
        endif
      end do
    end do
    
!   put at hour
    if(mod(icall,2).eq.0)then
      isave=icall/2
    endif
!    dz=(statein%ph3d(:,:,2)-statein%ph3d(:,:,1))/grvity*100. ! make cm
    allocate(dzaod(is:ie,js:je,nl))
    do k=1,nl
      dzaod(:,:,k)=(statein%ph3d(:,:,k+1)-statein%ph3d(:,:,k))/grvity ! make m
    end do
#ifndef CLDTAUZERO
    taucldfrc=statein%cldtau
#else
    if(lcldtauzero)then
      taucldfrc=0.
    else
      taucldfrc=statein%cldtau
    endif
#endif
!    call wetdep(de, dts, Statein%pb2d, Statein%rc2d, Statein%rn2d, Statein%ph3d, Statein%phl3d, &
!      Statein%pr3d, Statein%prl3d, Statein%tk3d, Statein%us3d, Statein%vs3d, Statein%ws3d, &
!      Statein%tr3d, Stateout%tr3d, wet_dep, ntra, numgas, num_chem, num_moist, &
!      its, ite, jts, jte, kte, kts, &
!      ims, ime, jms, jme, kme, kms, rc)
    iprn=0
    jprn=0
     iprnin=63
     jprnin=57
     iprnin=52
     jprnin=49
     iprnin=80
     jprnin=33
     iprnin=30
     jprnin=64
     kprnin=31
    iprnin=93
    jprnin=89
    iprnin=73
    jprnin=9
    kprnin=48
     iprnin=51
     jprnin=53
     kprnin=28
     iprnin=31
     jprnin=38
     kprnin=12
     iprnin=88
     jprnin=80
     kprnin=20
!    iprnin=42
!    jprnin=33
    
!    endif
    iamprn=-1
    iprn=-1
    jprn=-1
    iprnin=47
    jprnin=31
    kprnin=62
    iprnin=36
    jprnin=31
    iprnin=69
    jprnin=56
    kprnin=61
    iprnin=66
    jprnin=54
    kprnin=1
    iprnin=122
    jprnin=75
    tileprn=5
    ciat=' '
    call getenv('IAT',ciat)
    if(ciat.ne.' ')then
      read(ciat,*)iprnin
      call getenv('JAT',cjat)
      read(cjat,*)jprnin
      call getenv('KAT',ckat)
      read(ckat,*)kprnin
      call getenv('TILEAT',ctileat)
      read(ctileat,*)tileprn
      if(iam.eq.0)then
        write(6,*)'iprnin',iprnin,jprnin,kprnin,tileprn
      endif
      if(tile.eq.tileprn)then
        if(is<=iprnin.and.ie>=iprnin.and.js<=jprnin.and.je>=jprnin)then
          iprn=iprnin-is+1
          jprn=jprnin-js+1
          iamprn=iam
        endif
      endif
    endif
    if(mod(nint(dts*float(icall)),nsecsave).eq.0.or.doall.or.icall.eq.1)then
      doaod=.true.
      if(iam.eq.0)then
        write(6,*)'set doaod true'
      endif
    else
      doaod=.false.
    endif
    call raqms_advance(de,dts,advanceCount, &
    StateIn%pr3d, &
    StateIn%prl3d, &
    Statein%tk3d, &
    Statein%us3d, &
    statein%vs3d, &
    Statein%ws3d, &
    statein%ph3d, &
    statein%phl3d, &
    Statein%rn2d, &
    Statein%rc2d, &
    StateIn%tr3d, &
    chem_pass_state%tr3d, &
    StateIn%vtype2d, &
    Statein%exch, &
    Statein%dqdt, &
    StateOut%tr3d, &
    config%ntra, &
    config%nbegin, &
    config%numgas, &
    config%num_moist, &
    nchem, &
    lon,lat, &
    yy,mm,dd,h,m,s,tz,julday, &
!   these are full limits for this subregion
    is, ie, js, je, 1, nl, &
    is, ie, js, je, 1, ni, &
    data, doaod,&
    rc=localrc)
! here aerosols from gsdchem are in chem_pass_state%tr3d
! stateout will be used by gsi increment and gsi write
  Stateout%tr3d(:,:,:,p_bc1)=chem_pass_state%tr3d(:,:,:,p_bc1)
  Stateout%tr3d(:,:,:,p_bc2)=chem_pass_state%tr3d(:,:,:,p_bc2)
  Stateout%tr3d(:,:,:,p_oc1)=chem_pass_state%tr3d(:,:,:,p_oc1)
  Stateout%tr3d(:,:,:,p_oc2)=chem_pass_state%tr3d(:,:,:,p_oc2)
  Stateout%tr3d(:,:,:,p_sulf)=chem_pass_state%tr3d(:,:,:,p_sulf)
  Stateout%tr3d(:,:,:,p_dust1)=chem_pass_state%tr3d(:,:,:,p_dust1)
  Stateout%tr3d(:,:,:,p_dust2)=chem_pass_state%tr3d(:,:,:,p_dust2)
  Stateout%tr3d(:,:,:,p_dust3)=chem_pass_state%tr3d(:,:,:,p_dust3)
  Stateout%tr3d(:,:,:,p_dust4)=chem_pass_state%tr3d(:,:,:,p_dust4)
  Stateout%tr3d(:,:,:,p_dust5)=chem_pass_state%tr3d(:,:,:,p_dust5)
  Stateout%tr3d(:,:,:,p_seas1)=chem_pass_state%tr3d(:,:,:,p_seas1)
  Stateout%tr3d(:,:,:,p_seas2)=chem_pass_state%tr3d(:,:,:,p_seas2)
  Stateout%tr3d(:,:,:,p_seas3)=chem_pass_state%tr3d(:,:,:,p_seas3)
  Stateout%tr3d(:,:,:,p_seas4)=chem_pass_state%tr3d(:,:,:,p_seas4)
  Stateout%tr3d(:,:,:,p_seas5)=chem_pass_state%tr3d(:,:,:,p_seas5)
  Stateout%tr3d(:,:,:,p_so2)=chem_pass_state%tr3d(:,:,:,p_so2)*1.e-6 ! new ajl 6/30/2023
    if(gsiremapsh/=' ')then
      if(iam==0)then
        write(6,*)'call chemgsi',nsecgsi
      endif
      call chemgsi(stateout,statein,is,ie,js,je,nl,de,icall,dts,nsecgsi,modelcomm,data)
    endif
    if(mod(nint(dts*float(icall)),nsecsave).eq.0.or.doall)then
      fracday=dts*float(icall)/86400.
      forecast_hr=dts*float(icall)/3600.
      call calcaod(stateout%tr3d,dzaod,data,nchem,is,ie,js,je,nl,statein%prl3d,statein%tk3d,de)
      call calcaodgsi(stateout%tr3d,dzaod,data,nchem,is,ie,js,je,nl,statein%pr3d,statein%prl3d,statein%tk3d,de)
!    if(mod(icall,6).eq.0)then
!      isave=icall/6
      isave=1
      dims(1)=ie-is+1
      dims(2)=je-js+1
      dims(3)=nl
      dims2=dims(1:2)
!      call calcaod(statein%tr3d,dzaod,data,nchem,is,ie,js,je,nl,statein%prl3d,statein%tk3d,de)
      if(mype.eq.0)then
        write(6,*)'call raqmschem_chem_write ',cdate
        call flush(6)
      endif
 
      rtime1=mpi_wtime()
      cdateout=cdate
      cdateat=cdate
!      write(cdateout(1:2),'(i2.2)')icall
      call raqmschem_model_get(modelcomm=modelcomm,rc=rc)
      if(iam==0)then
         write(6,*)'doaodfcst',doaodfcst
      endif
      if(doaodfcst)then
        if(iam==0)then
          write(6,*)'write_aodfcst ok ',doaodfcst,icall
          call flush(6)
        endif
        call raqmschem_chem_write_aodfcst(stateout%tr3d,data,statein,dims,itime=isave,cdate=cdateout,rc=rc)
      else
        if(iam==0)then
        write(6,*)'write_full bad  ',doaodfcst,icall
        call flush(6)
        endif
        if(iam==0)then
          write(6,*)'cdatestart',trim(cdatestart),'cdateout',trim(cdateout)
          write(6,*)'BASE_PATH',trim(BASE_PATH)
          write(6,*)'remap05',remap05
          write(6,*)'realtimerun',trim(crealtimerun)
        endif 
        call raqmschem_chem_write(stateout%tr3d,data,statein,dims,itime=isave,cdate=cdateout,rc=rc)
      if(iam==0)then
        if(crealtimerun=='YES')then
          call system('export CDATE='//trim(cdate)//' ; export CDATEBEG='//trim(cdatestart)// &
          '; export BASE_PATH='//trim(base_path)//'; export SCRATCH_OUT='//trim(scratch_out)// &
          ' ;  sbatch  -t 10 --export=ALL '//trim(remap05))
        endif
      endif
      call raqmschem_chem_write_o3vmr(dims,itime=isave,cdate=cdateout,rc=rc)
      endif
      wallwrite=mpi_wtime()-rtime1
      wallwritecum=wallwritecum+wallwrite
      if(iam.eq.0)then
        open(30,file='tracer.done'//trim(cdateout),form='formatted')
        write(30,'(a4)')'done'
        close(30)
      endif
      if(iam.eq.0.and.mod(nstepat,nhtuw).eq.0)then
        write(6,*)'wallwrite',wallwrite,' cum ',wallwritecum/60.,' min '
        call flush(6)
      endif
    endif
    if(useraqmso3vmr)then
!     convert to kg/kg from ppv 
      do k=1,nl
        stateout%tr3d(:,:,k,p_atm_o3mr)=mwo3/airmw*o3vmr_inst(:,:,nl+1-k)
      end do
    else
       stateout%tr3d(:,:,:,p_atm_o3mr)=statein%tr3d(:,:,:,p_atm_o3mr)
    endif
!    write(6,*)'o3mr from raqms',maxval(stateout%tr3d(:,:,:,3))
!      call raqmschem_chem_write_var(o3vmr_inst,dims,isave,'o3vmr',cdate=cdate,rc=rc)
!      call raqmschem_chem_write_var(thgrd,dims,isave,'theta',cdate=cdate,rc=rc)
!      o3dep_save=o3dep_save*86400.
!      call raqmschem_chem_write_var_2d(o3dep_save,dims2,isave,'o3dep',cdate=cdate,rc=rc)
!      noydep_save=noydep_save*86400.
!      call raqmschem_chem_write_var_2d(noydep_save,dims2,isave,'noydep',cdate=cdate,rc=rc)
!      !codep_save=codep_save*86400.
!      call raqmschem_chem_write_var_2d(codep_save,dims2,isave,'codep',cdate=cdate,rc=rc)
!      o3dep_save=0.0
!      noydep_save=0.0
!      codep_save=0.0
!      navg=0
!    if(first)then
!      dz=dz*.01 ! makd meters
!      call chem_io_write('slmsk2d.nc',slmsk2d,path=trim(config%emi_outname),de=de,rc=rc)
!      call chem_io_write('zlwi.nc',zlwigrd,path=trim(config%emi_outname),de=de,rc=rc)
!      call chem_io_write('dz.nc',dz,path=trim(config%emi_outname),de=de,rc=rc)
       
!    endif
!    if(first)then
!    call chem_io_write('o3dep.nc',o3dep_save,path=trim(config%emi_outname),de=de,rc=rc)
!    endif
    if(allocated(dzaod))then
      deallocate(dzaod)
    endif
    call deallocatechemcommtemp
  end do ! de
  if(first)rtimecum=0.0
  first=.false.
  rtimed1=mpi_wtime()
  if(iam.eq.0.and.tile.eq.1)then
  write(6,*)'wall cum ',(rtimed1-wallstart)/60.
  call flush(6)
  endif
!  write(200+iam,*)'bottom of raqms_model_advance'
!  flush(200+iam)
  return
  end subroutine raqms_model_advance
    subroutine getfsize(filename,size)
    use ifport
    use ifport_types
    character*(*) filename
    integer*8 size
    integer iresult
    integer(int_ptr_kind()) handle
    type(file$infoi8) info
    handle=file$first
    iresult=getfileinfoqq(filename,info,handle)
    if(iresult.ne.0)then
       size=info%LENGTH
    else
      size=0
    endif
    return
    end subroutine getfsize
      subroutine trop(is,ie,js,je,nl,psfc,prl,zm,tm)
      use raqmschem_pmgrid_mod, only : tile,ibeg,iam,rsearch
      use raqmschemcomm_mod, only : ztropraq,ptropraq,ktrop,ptropgsi,ztropgsi,ktropgsi
      implicit none
      integer i,j,k,ktopbot,is,ie,js,je,nl,kbot,ktop,ktropbot
      real*4,dimension(is:ie,js:je,nl) :: prl,zm,tm
      real*4,dimension(is:ie,js:je) :: psfc
      real*4 :: dz,gam,gamctrl,dzlay,p(nl),hphd(nl)
      integer :: iflag
      real*4,parameter :: ptplim(2)=(/500.,50./),hd=2000.,gamtp=2.e-3,one=1.0,half=0.5
      real*4 gamu,gamd,gami,td,wtp,ttrop,htp,ptp,rd_over_g,grav,rd,gamk(nl)
      integer kmm2,klim(2),ktp,kd(1)
    if(.not.allocated(ztropraq))then
      allocate(ztropraq(is:ie,js:je),ptropraq(is:ie,js:je),ktrop(is:ie,js:je))
      allocate(ztropgsi(is:ie,js:je),ptropgsi(is:ie,js:je),ktropgsi(is:ie,js:je))
    endif
      grav   = 9.80665
      rd     = 2.8705e+2
      rd_over_g  = rd/grav
      kmm2=nl-2
      gamctrl=-2.e-3
      do j=js,je
        do i=is,ie
          p=prl(i,j,:)
          call rsearch(1,kmm2,1,1,p(2:),2,1,1,ptplim,1,1,klim)
          klim(1)=klim(1)+1
          klim(2)=min(kmm2,klim(2))
          hphd=zm(i,j,:)+hd
          ptropraq(i,j)=-9999.
          ztropraq(i,j)=-9999.
          do k=1,nl
            if(prl(i,j,k)<psfc(i,j)-300.)then
              kbot=k
              exit
            endif
          end do
          do k=nl,1,-1
            if(prl(i,j,k)>50.)then
              ktop=k
              exit
            endif
          end do
          iflag=0
!          write(200+iam,*)i,j,'kbot',kbot,ktop,psfc(i,j),'prl1',prl(i,j,1),tm(i,j,1),zm(i,j,1)
          
          do k=kbot,ktop
            dz=zm(i,j,k+1)-zm(i,j,k)
            gam=(tm(i,j,k+1)-tm(i,j,k))/dz
            gamk(k)=gam
            if(gam>=gamctrl)then
              if(iflag==0)then
                ktropbot=k
                dzlay=dz
                iflag=1
                ptropraq(i,j)=prl(i,j,k)
                ztropraq(i,j)=zm(i,j,k)
                ktrop(i,j)=k
              elseif(iflag==1)then
                dzlay=dzlay+dz
                if(dzlay>2000.)then
                  exit
                endif
              endif
            else
              if(iflag==1)then
                if(dzlay<2000.)then
                  iflag=0
                  ptropraq(i,j)=-9999.
                  ztropraq(i,j)=-9999.
                  ktrop(i,j)=-100
                  iflag=0
                  dzlay=0.0
                endif
              endif
            endif
          end do
          if(ptropraq(i,j)<50.0)then
            ptropraq(i,j)=50.
          endif
          gamd=1.e9
          ktp=klim(2)
          wtp=0.0
          do k=klim(1),klim(2)
            dz=zm(i,j,k+1)-zm(i,j,k)
            gamu=-(tm(i,j,k+1)-tm(i,j,k))/dz
            if(gamu<=gamtp)then
              call rsearch(1,nl-k-1,1,1,zm(i,j,k+1:),1,1,1,hphd(k:),1,1,kd)
              if(k+kd(1)>nl.or.k+kd(1)+1>nl)then
               write(200+iam,*)'error k',k,'kd',kd(1),'sum',k+kd(1),'nl',nl
               call flush(200+iam)
               td=tm(i,j,k)
              else
                td=tm(i,j,k+kd(1))+(hphd(k)-zm(i,j,k+kd(1)))/(zm(i,j,k+kd(1)+1)-zm(i,j,k+kd(1)))* &
                (tm(i,j,k+kd(1)+1)-tm(i,j,k+kd(1)))
              endif
              gami=(tm(i,j,k)-td)/hd
              if(gami<=gamtp)then
                ktp=k
                wtp=(gamtp-gamu)/(max(gamd,gamtp+0.1e-3)-gamu)
                exit
              endif
            endif
            gamd=gamu
          end do
          if(ktp>0)then
            ttrop=tm(i,j,ktp)-wtp*(tm(i,j,ktp)-tm(i,j,ktp-1))
            htp=zm(i,j,ktp)-wtp*(zm(i,j,ktp)-zm(i,j,ktp-1))
            ptp=p(ktp)*exp((zm(i,j,ktp)-htp)*(one-half*(ttrop/tm(i,j,ktp)-one))/(rd_over_g*tm(i,j,ktp)))
            ptropgsi(i,j)=ptp*100. ! Pa
            ztropgsi(i,j)=htp
            ktropgsi(i,j)=ktp
         endif
        end do
      end do
      return
      end subroutine trop 
end module raqms_model_mod
