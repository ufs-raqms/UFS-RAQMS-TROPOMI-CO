module chemgsimod 
  contains 
  subroutine chemgsi(stateout,statein,is,ie,js,je,nl,de,icall,dts,nsecgsi,modelcomm,data)
  use mpi
  use chem_types_mod
  use raqmschem_state_mod
  use raqmschem_config_mod
  use raqmschem_data_mod
  use raqmschem_model_mod
  use chem_comm_mod, only : chem_comm_get
!  use raqmschem_model_mod, only : chem_model_get => raqmschem_model_get
  use raqmschem_io_mod, only : chem_io_write,chem_io_read
  use raqmschem_iodata_mod, only : raqmschem_gsi_write,flipz
  use raqms_mod
  use chem_raqms_mod, only : chem_pass_state
  use raqmschem_pmgrid_mod, only : tile,ibeg,iam,numgsicem,gsivar,gsicem,nstepat,cdateat,gsivarr
!  use raqmschem_pmgrid_mod, only : loccem,forecast_hr,ompsgsi,aerosol_ugpkg,raqms_localIOflag
  use raqmschem_pmgrid_mod, only : forecast_hr,ompsgsi,aerosol_ugpkg,raqms_localIOflag
  use raqmschem_pmgrid_mod, only : gsiarray
!  use raqmschemcomm_mod, only : gsiinc
  use raqmschemcomm_mod, only : incMLS3d,incOMI3d,o3vmr_inst
  use raqmschemcomm_mod, only : incNUCAPSco3d,incTROPOMIco3d
  use raqmschem_io_mod,only : inquire_file_var
  use raqmschem_species_mod
  use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
  implicit none
  type(chem_state_type),  pointer ::  stateOut,statein
  type(chem_data_type),   pointer :: data
  integer dims(3),is,ie,js,je,nl,icall,ibegdate,ienddate,nsecgsi,idate,mype
  integer imod,ihr,n,i,j,k,l,m,isave,countsleep,localrc,rc,jj,ii,nn,mm
  integer idtsec,julday,modelcomm,de,nsecgsi2
  character *256 GSIREMAPSH,cstep*4,base_data_disk,cdatestart*11,cdateout*10
  character *256 gsidatapath,gsiname,pathdata,datascratch,ctile*1,base_path
  character *256 gsiinc_path,gsi_inq_script
  character *20 filename,cnsecgsi2*10
  logical exist,exist2,first,dogsiname,existt,haveNUCAPS,checkthere
  real(chem_kind_r4),allocatable,dimension(:,:,:,:) :: tempgsiinc
  real(CHEM_KIND_R4),allocatable,dimension(:,:,:) :: tempincMLS,tempincOMI
  real(CHEM_KIND_R4),allocatable,dimension(:,:,:) :: tempincNUCAPS,tempincTROPOMI,o3vmr_inst_inc
  real(CHEM_KIND_R4),allocatable,dimension(:,:) :: tempincaodgsi
  real(CHEM_KIND_R8) :: gsinew,gsioxnew,hrout
  
  real(CHEM_KIND_R8) :: dts 
  real(CHEM_KIND_R4) :: af,fhmax
  integer,pointer,dimension(:) :: loccem
  integer,pointer :: ngsivar
  character*8,pointer,dimension(:) :: gsichemp
  real*4,pointer,dimension(:,:,:,:) :: gsiinc
  save first,gsiremapsh,cdatestart,GSIDATAPATH,GSINAME,GSI_INQ_SCRIPT,fhmax
  data first/.true./
#include <comcdate.h>
  cnsecgsi2=' '
  call getenv('NSECGSI2',cnsecgsi2)
  if(cnsecgsi2/=' ')then
    read(cnsecgsi2,*)nsecgsi2
  else
   nsecgsi2=-1
  endif
  

  base_data_disk=' '
  call getenv('base_data_disk',base_data_disk) 
  mype=iam
!  write(6,*)'top chemgsi ',tile,mype,'iam',iam
!  call flush(6)
  cdatestart=' '
  call getenv('CDATE',cdatestart)
  if(first)then
    if(tile.eq.1.and.mype.eq.0)then
      write(6,*)'CDATESTART',cdatestart
      call flush(6)
    endif
  endif
  gsiremapsh=' '
  call getenv('GSIREMAPSH',gsiremapsh)
  if(first)then
    if(iam.eq.0.and.tile.eq.1)then
      write(6,*)'base_data_disk',trim(base_data_disk)
      write(6,*)'gsiremapsh',gsiremapsh
      call flush(6)
    endif
  endif
  GSIDATAPATH=' '
!  write(6,*)mype,'getenv GSIDATAPATH'
!  call flush(6)
  call getenv('GSIDATAPATH',GSIDATAPATH)
  if(first)then
    if(iam.eq.0.and.tile.eq.1)then
      write(6,*)mype,'did  getenv GSIDATAPATH'
      call flush(6)
    endif
  endif
  if(GSIDATAPATH/=' ')then
    if(first)then
      if(mype.eq.0)then
        write(6,*)mype,'GSIDATAPATH',trim(GSIDATAPATH)
        call flush(6)
       endif
      endif
    endif
    GSINAME=' '
    call getenv('GSINAME',gsiname)
    if(GSINAME/=' ')then
      dogsiname=.true.
      if(first)then
      if(mype==0)then
        write(6,*)'GSINAME',trim(GSINAME)
        flush(6)
        write(6,*)'GSIDATAPATH',trim(GSIDATAPATH)
        flush(6)
      endif
      endif
!      gsiinqscript=' '
!      call getenv('GSIINQSCRIPT',gsiinqscript)
    else
      dogsiname=.false.
    endif
    if(dogsiname)then
      do  i=1,len_trim(gsiname)-4
        if(gsiname(i:i+4).eq.'CDATE')then
          ibegdate=i
          ienddate=i+4
      if(first)then
          if(mype==0)then
          write(6,*)'ibegdate',ibegdate,ienddate
          flush(6)
          endif
          endif
          exit
        endif
      end do
    endif
    base_path=' '
    call getenv('BASE_PATH',base_path)
    if(first)then
      if(iam.eq.0.and.tile.eq.1)then
        write(6,*)'chemgsi base_path',iam,tile,trim(base_path)
      write(6,*)'base_path',base_path
      write(6,*)'dogsiname',dogsiname
    endif
  endif ! datapath
#if 0
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
#endif

! ni should be second
! put at hour
  write(cstep,'(i4.4)')icall
  filename='each'//trim(cstep)//'.nc'
  dims(1)=ie-is+1
  dims(2)=je-js+1
  dims(3)=nl
  cdateat=cdate
  call gsivar
  if(first)then
    if(iam.eq.0)then
      write(6,*)'ajl dts',dts,'icall',icall,'nsecgsi',nsecgsi
      call flush(6)
    endif
  endif
! if not time return ?
  if(mod(nint(dts*float(icall)),nsecgsi)/=0)then
    if(nsecgsi2>=1)then
      if(mod(nint(dts*float(icall)),nsecgsi2)/=0)then
        return
      endif
    else
      return
    endif
  endif
  
!  if(mod(nint(dts*float(icall)),nsecgsi).eq.0)then
   cdateout=cdate
  checkthere=.false.
  if(checkthere)then
  if(iam.eq.0.and.tile.eq.1)then
    write(6,*)'findajl cdate',cdate
    write(6,*)'findajl icall',icall
    call flush(6)
    if(base_data_disk/=' ')then
!     check if OBSdir is there
      exist=.false.
      write(6,*)'numgsicem',numgsicem
      flush(6)
!     AOD
      if(numgsicem>=14)then
        write(6,*)'gsicem',gsicem(1:numgsicem)
        flush(6)
        if(gsicem(1).eq.'sulf')then
          if(dogsiname)then
            write(6,*)'set up data path aerosol',ibegdate,ienddate,len_trim(gsiname)
            flush(6)
            if(gsiname(1:ibegdate-1)=='viirs.aod.')then
!             make sure multiple of 3 hrs
              read(cdateout(1:10),'(i10)')idate
              write(6,*)'idate',idate
              ihr=mod(idate,100)
              imod=mod(ihr,3)
              write(6,*)'ihr',ihr,'imod',imod
              if(imod==1)then
                idtsec=-3600
                call advcdatehr(idate,idtsec,idate,julday,hrout)
                write(cdateout(1:10),'(i10)')idate
                write(6,*)'idate',idate
              elseif(imod==2)then
                  idtsec=3600
                  call advcdatehr(idate,idtsec,idate,julday,hrout)
                  write(cdateout(1:10),'(i10)')idate
                  write(6,*)'idate',idate
                endif 
              
                pathdata=trim(gsidatapath)//gsiname(1:ibegdate-1)//cdateout(1:10)//gsiname(ienddate+1:)
              else
                pathdata=trim(gsidatapath)//gsiname(1:ibegdate-1)//cdateout(1:10)//gsiname(ienddate+1:)
              endif
              write(6,*)'aerosols path',trim(pathdata)
              flush(6)
             inquire(file=pathdata,exist=existt)
             if(existt)then
               exist=.true.
             else
               exist=.false.
             endif
            endif
          endif
          write(6,*)'ajl gsi name',exist
          flush(6)
        else ! not AOD
          do n=1,numgsicem
            if(dogsiname)then
!              write(6,*)'set up data path ',ibegdate,ienddate,len_trim(gsiname)
!              flush(6)
              pathdata=trim(gsidatapath)//'/'//gsiname(1:ibegdate-1)//cdateout(1:10)// &
              gsiname(ienddate+1:)
!              write(6,*)'path data',trim(pathdata)
!              flush(6)
             inquire(file=pathdata,exist=existt)
!             write(6,*)'pathdata',pathdata,'exist',existt
             if(existt)then
               exist=.true.
             else
               exist=.false.
             endif
            else
              if(gsicem(n).eq.'o3vmr')then
                pathdata=trim(base_data_disk)//'/'//cdateout(1:4)//'/'//cdateout(1:6)//'/' &
                //cdateout(1:10)//'/'//'mlsbufr.'//cdateout(1:10)
              else
                pathdata=trim(base_data_disk)//'/'//cdateout(1:4)//'/'//cdateout(1:6)//'/' &
                //cdateout(1:10)//'/'//'tropomi.'//trim(gsicem(n))//'.bufr.'//cdateout(1:10)
              endif
            endif
!            write(6,*)'pathdata define ',trim(pathdata),'exist',exist
!            call flush(6)
            inquire(file=pathdata,exist=existt)
            if(existt)then
              write(6,*)'exist',existt,' ',trim(pathdata)
              exist=.true.
!            elseif(gsicem(m).eq.'o3vmr')then
!            fix 2/25/2021
            elseif(gsicem(n).eq.'o3vmr')then
              write(6,*)'check omi'
            call flush(6)
              pathdata=trim(base_data_disk)//'/'//cdateout(1:4)//'/'//cdateout(1:6)//'/' &
              //cdateout(1:10)//'/'//'omibufr.'//cdateout(1:10)
              write(6,*)'pathdata omi ',trim(pathdata)
            call flush(6)
              inquire(file=pathdata,exist=existt)
              if(existt)then
                write(6,*)'exist omi'
                exist=.true.
              endif
            endif
            if(gsicem(n).eq.'o3vmr')then
              write(6,*)'mlsomi exist',exist
              write(6,*)'pathdata',trim(pathdata)
            endif
          end do ! n
        endif ! end AOD NOT AOD
      else ! end base_data_disk
        write(6,*)'base_data_disk blank'
        call flush(6)
        exist=.true.
      endif ! end base data_disk
    endif
  endif ! checkthere
    call raqmschem_comm_all_bcast(exist,rc=localrc)
!  endif ! end check if do gsi and should have determined exist
! return if not exist ?
!  if(mod(nint(dts*float(icall)),nsecgsi).eq.0.and.exist)then
  exist=.true. ! for testing
  if(.not.exist)return
#ifdef PERGSIV
    if(.not.allocated(pergsiv))then
!      allocate(percov(is:ie,js:je,nl),intco(is:ie,js:je,nl))
      allocate(pergsiv(is:ie,js:je,nl,numgsicem))
      allocate(colgsiinc(is:ie,js:je,nl,numgsicem))
    endif
#endif
  cdateout=cdate
  isave=1
  dims(1)=ie-is+1
  dims(2)=je-js+1
  dims(3)=nl
! write(cdateout(1:2),'(i2.2)')icall
  call raqmschem_model_get(modelcomm=modelcomm,rc=rc)
  if(iam.eq.0.and.tile.eq.1)then
     write(6,*)'gsi_write',cdateout
     call flush(6)
  endif
  call raqmschem_gsi_write(stateout%tr3d,data,statein,dims,itime=isave,cdate=cdateout,rc=rc) 
!  call mpi_barrier(modelcomm,rc)
  if(iam.eq.0.and.tile.eq.1)then
    write(6,*)'done with raqmschem_gsi_write ',cdateout
    call flush(6)
  endif
!  datascratch='/scratch/users/lenzen/GSISCRATCH/'//trim(expname)//'/'//cdatestart//'/'
  datascratch='../../GSISCRATCH/'
 
  if(iam.eq.0.and.tile.eq.1)then
    write(6,*)'datascratch',trim(datascratch)
    write(6,*)'cdatestart',cdatestart
    call flush(6)
!    call system('mkdir -p ../../GSISCRATCH/'//trim(cdateout))
    call system('mkdir -p '//trim(datascratch))
    call system('/bin/rm -f '//trim(datascratch)//'gsi.done')
    call system('/bin/rm -f '//trim(datascratch)//'input.ready')
    call system('/bin/rm -f '//trim(datascratch)//'cdate.inc')
    write(6,*)'open cdate.inc',trim(datascratch)//'cdate.inc'
    call flush(6)
    write(6,*)'findajl cdatestart',cdatestart,'cdateout',cdateout
    call flush(6)
    open(40,file=trim(datascratch)//'/cdate.inc',form='formatted')
    write(40,'("export CDATE0=",a10)')trim(cdatestart)
    write(40,'("export CDATE=",a10)')trim(cdateout)
    forecast_hr=dts*float(icall)/3600.
    write(6,*)'forecast_hr',forecast_hr,'fhmax',fhmax
    call flush(6)
!    if(forecast_hr+.00001>fhmax)then
    if(cdateout(9:10).eq.'12')then
      write(6,*)'can end GSI part'
      call flush(6)
      write(40,'("export DONE=true")')
    else
      write(6,*)'keep running GSI'
      call flush(6)
      write(40,'("export DONE=false")')
    endif
    close(40)
    write(6,*)'open input.ready ',trim(datascratch)//'input.ready'
    call flush(6)
    call system('rm -f '//trim(datascratch)//'/gsi.done.'//trim(cdateout))
    call system('sleep 10')
    open(40,file=trim(datascratch)//'/input.ready',form='formatted')
    write(40,'(a10)')trim(cdateout)
    close(40)
    write(6,*)'input.ready'
    call flush(6)
    call system('date')
  endif
  countsleep=0
  exist=.false.
  do while (countsleep<100.and..not.exist)
    if(iam.eq.0.and.tile.eq.1)then
      write(6,*)'countsleep',countsleep
      write(6,*)'check gsi.done ',trim(datascratch)//'/gsi.done.'//trim(cdateout)
      call flush(6)
      inquire(file=trim(datascratch)//'/gsi.done.'//trim(cdateout),exist=exist)
    endif
    call raqmschem_comm_all_bcast(exist,rc=localrc)
    if(exist)then
      if(iam.eq.0.and.tile.eq.1)then
        write(6,*)'gsidone'  ,countsleep
      endif
    else
      !write(6,*)'gsi not done'
      call system('sleep 30')
      countsleep=countsleep+1
    endif
  end do
  call mpi_barrier(modelcomm,rc)
  if(iam.eq.0.and.tile.eq.1)then
    if(exist)then
      call system('rm -f '//trim(datascratch)//'/gsi.done.'//trim(cdateout))
    endif
    write(6,*)'gsi done exist',exist,'countsleep',countsleep
    call flush(6)
    call system('date')
  endif
  call system('sleep 10')
  if(iam.eq.0.and.tile.eq.1)then
!   make sure don't goof up next time period
    call system('rm -f '//trim(datascratch)//'/gsi.done.'//trim(cdateout))
  endif
  write(ctile,'(i1.1)')tile
  if(iam.eq.0.and.tile.eq.1)then
    call system('date')
 
    write(6,*)'read from ../../GSIINC/'//trim(cdatestart)//'/inc.' &
    //trim(cdateout)//'.tile'//ctile//'.nc'
    call flush(6)
  endif
  if(base_path/=' ')then
    inquire(file=trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/'//'inc.'//trim(cdateout)//'.tile'//ctile//'.nc',exist=exist)
    if(iam.eq.0)then
       write(6,*)'GSIINC',trim(base_path)//'/GSIINC/'//trim(cdatestart)// &
       '/'//'inc.'//trim(cdateout)//'.tile'//ctile//'.nc'
    endif
    if(exist)then
      GSIINC_PATH=trim(base_path)//'/GSIINC/'
    endif
  else
    inquire(file='../../GSIINC/'//trim(cdatestart)//'/'//'inc.'// &
     trim(cdateout)//'.tile'//ctile//'.nc',exist=exist)
    if(exist)then
      GSIINC_PATH='../../GSIINC/'
    endif
  endif
  if(.not.exist)then
    if(base_path/=' ')then
      write(6,*)'missing ',mype,trim(base_path)//'/GSIINC/'//trim(cdatestart)//  &
       '/inc.'//trim(cdateout)//'.tile'//ctile//'.nc'
    else
      write(6,*)'missing ',mype,'../../GSIINC/'//trim(cdatestart)//  &
           '/inc.'//trim(cdateout)//'.tile'//ctile//'.nc'
      call flush(6)
    endif
  endif
  do n=1,numgsicem
    loccem=>gsiarray(n)%loccem
    gsichemp=>gsiarray(n)%gsicem
    ngsivar=>gsiarray(n)%ngsivar
!    write(200+iam,*)'gsicem',shape(gsicem),'numgsicem',numgsicem
!    call flush(200+iam)
!    write(200+iam,*)'gsicem',gsicem,'n',n
!    call flush(200+iam)
!    write(200+iam,*)'gsicem',gsicem(n)
!    call flush(200+iam)
  
 
    allocate (tempgsiinc(is:ie,js:je,nl,ngsivar))
!    write(200+iam,*)'loccem',lbound(loccem),'ub',ubound(loccem),'n',n
!    write(200+iam,*)'ngsivar',ngsivar
!    call flush(200+iam)
!    write(6,*)'ngsivar',ngsivar
!    write(6,*)'allocated gsiarray%gsiinc ',allocated(gsiarray(n)%gsiinc)
!    call flush(6)
    if(.not.allocated(gsiarray(n)%gsiinc))then
!      if(mype==0)then
        !write(6,*)'allocated gsiinc ngsivar',ngsivar
!        write(6,*)'numgsicem',numgsicem,'is',is,ie,js,je,nl,'ngsivar',ngsivar
!        call flush(6)
!      endif
      allocate(gsiarray(n)%gsiinc(is:ie,js:je,nl,ngsivar))
      gsiarray(n)%gsiinc=0.0
      gsiinc=>gsiarray(n)%gsiinc
    else
      gsiinc=>gsiarray(n)%gsiinc
!      write(6,*)'gsiinc',shape(gsiinc),'n',n,'gsicem',gsicem(n)
!      call flush(6)
    endif
!    write(6,*)'GSIINC_PATH',trim(GSIINC_PATH)
!    call flush(6)
    if(gsicem(n).eq.'o3vmr')then
      if(.not.allocated(incMLS3d))then
         allocate(incMLS3d(is:ie,js:je,nl),incOMI3d(is:ie,js:je,nl))
         incMLS3d=0.0
         incOMI3d=0.0
      endif
      if(inquire_file_var('inc.'//trim(cdateout)//'.tile1.nc','incMLS', &
        trim(GSIINC_PATH)//trim(cdatestart)//'/' ))then
        if(.not.allocated(tempincMLS))then
          allocate (tempincMLS(is:ie,js:je,nl))
        endif
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempincMLS, &
        trim(GSIINC_PATH)//trim(cdatestart)//'/', &
        de=de,varname='incMLS',rc=localrc)
      elseif(inquire_file_var('inc.'//trim(cdateout)//'.tile1.nc','inco3lp', &
        trim(GSIINC_PATH)//trim(cdatestart)//'/' ))then
        if(.not.allocated(tempincMLS))then
          allocate (tempincMLS(is:ie,js:je,nl))
        endif
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempincMLS, &
        trim(GSIINC_PATH)//trim(cdatestart)//'/', &
        de=de,varname='inco3lp',rc=localrc)
        ompsgsi=.true.
      endif
      if(inquire_file_var('inc.'//trim(cdateout)//'.tile1.nc','incOMI', &
        trim(GSIINC_PATH)//trim(cdatestart)//'/' ))then
        if(.not.allocated(tempincOMI))then
          allocate (tempincOMI(is:ie,js:je,nl))
        endif
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempincOMI, &
        trim(GSIINC_PATH)//trim(cdatestart)//'/', &
        de=de,varname='incOMI',rc=localrc)
      elseif(inquire_file_var('inc.'//trim(cdateout)//'.tile1.nc','inco3nm', &
        trim(GSIINC_PATH)//trim(cdatestart)//'/' ))then
        if(.not.allocated(tempincOMI))then
          allocate (tempincOMI(is:ie,js:je,nl))
        endif
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempincOMI, &
        trim(GSIINC_PATH)//trim(cdatestart)//'/', &
        de=de,varname='inco3nm',rc=localrc)
        ompsgsi=.true.
      endif
!      write(6,*)'temp LP',maxval(tempincMLS),' NM ',maxval(tempincOMI)
!      write(6,*)'gsicem1',gsicem(n)
!      call flush(6)
!    endif ! gsicem==o3vmr
    elseif(gsicem(n).eq.'co')then
!     allocate (tempincTROPOMI(is:ie,js:je,nl),tempincNUCAPS(is:ie,js:je,nl))
!     if(.not.allocated(incNUCAPSco3d))then
!       !allocate(incNUCAPSco3d(is:ie,js:je,nl),incTROPOMIco3d(is:ie,js:je,nl))
!       incNUCAPSco3d=0.0
!       incTROPOMIco3d=0.0
!     endif
!     tempincNUCAPS=0.0
!     tempincTROPOMI=0.0
      haveNUCAPS=.false.
!      write(6,*)'path',trim(GSIINC_PATH)//trim(cdatestart)//"/"
      if(inquire_file_var('inc.'//trim(cdateout)//'.tile1.nc','incco_NU', &
        trim(GSIINC_PATH)//trim(cdatestart)//'/' ))then
        haveNUCAPS=.true.
        if(.not.allocated(tempincNUCAPS))then
          allocate (tempincNUCAPS(is:ie,js:je,nl))
        endif
        if(.not.allocated(incNUCAPSco3d))then
          allocate(incNUCAPSco3d(is:ie,js:je,nl))
          incNUCAPSco3d=0.0
        endif
        tempincNUCAPS=0.0
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempincNUCAPS, &
        trim(GSIINC_PATH)//trim(cdatestart)//'/', &
        de=de,varname='incco_NU',rc=localrc)
      endif
!    endif ! gsicem==co
!    if(gsicem(n)=='CO')then
      if(inquire_file_var('inc.'//trim(cdateout)//'.tile1.nc','incco_TR', &
        trim(GSIINC_PATH)//trim(cdatestart)//'/' ))then
        if(.not.allocated(tempincTROPOMI))then
          allocate (tempincTROPOMI(is:ie,js:je,nl))
        endif
        if(.not.allocated(incTROPOMIco3d))then
          allocate(incTROPOMIco3d(is:ie,js:je,nl))
         incTROPOMIco3d=0.0
        endif
        tempincTROPOMI=0.0
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempincTROPOMI, &
           trim(GSIINC_PATH)//trim(cdatestart)//'/', &
           de=de,varname='incco_TR',rc=localrc)
      endif ! inquire
!      if(havenucaps)then
        !write(6,*)'temp NUCAPS',maxval(tempincNUCAPS),' TROMI ',maxval(tempincTROPOMI)
!        call flush(6)
!      else
        !write(6,*)'temp TROMI ',maxval(tempincTROPOMI)
!        call flush(6)
!      endif
!    endif ! gsicem==co
    else if(gsicem(n)=='AOD')then
!     doing aod assimilation
      if(base_path/=' ')then
        if(.not.allocated(data%aodincgsi))then
          allocate(data%aodincgsi(is:ie,js:je))
!          allocate(data%aodgsi(is:ie,js:je),data%aodincgsi(is:ie,js:je))
          data%aodincgsi=0.0
        endif
        if(.not.allocated(tempincaodgsi))then
          allocate(tempincaodgsi(is:ie,js:je))
        endif
#ifdef DOALLGSIAOD
        if(.not.allocated(data%aodgsibcoc))then
          allocate(data%aodgsibcoc(is:ie,js:je))
        endif
        if(.not.allocated(data%aodgsidust))then
          allocate(data%aodgsidust(is:ie,js:je))
        endif
        if(.not.allocated(data%aodgsissalt))then
          allocate(data%aodgsissalt(is:ie,js:je))
        endif
        if(.not.allocated(data%aodgsisulf))then
          allocate(data%aodgsisulf(is:ie,js:je))
        endif
#endif
!        write(200+iam,*)'chem_io_read aod'
!        call flush(200+iam)
#ifdef DOALLGSIAOD
        call chem_io_read('inc.'//trim(cdateout)//'.nc',data%aodgsi, &
          trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/', &
          de=de,varname='aod',rc=localrc)
        if(localrc/=0)then
          write(6,*)'error reading aodgsi'
          deallocate(data%aodgsi)
        endif
        call chem_io_read('inc.'//trim(cdateout)//'.nc',data%aodgsibcoc, &
          trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/', &
          de=de,varname='aodbcoc',rc=localrc)
        if(localrc/=0)then
          write(6,*)'error reading aodgsibcoc'
          deallocate(data%aodgsibcoc)
        endif
        call chem_io_read('inc.'//trim(cdateout)//'.nc',data%aodgsidust, &
          trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/', &
        de=de,varname='aoddust',rc=localrc)
        if(localrc/=0)then
          write(6,*)'error reading aodgsidust'
          deallocate(data%aodgsidust)
        endif
        call chem_io_read('inc.'//trim(cdateout)//'.nc',data%aodgsissalt, &
           trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/', &
           de=de,varname='aodssalt',rc=localrc)
        if(localrc/=0)then
          write(6,*)'error reading aodgsissalt'
          deallocate(data%aodgsissalt)
        endif
        call chem_io_read('inc.'//trim(cdateout)//'.nc',data%aodgsisulf, &
           trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/', &
           de=de,varname='aodsulf',rc=localrc)
        if(localrc/=0)then
          write(6,*)'error reading aodgsisulf'
          deallocate(data%aodgsisulf)
        endif
#endif
!        write(200+iam,*)'chem_io_read inc aod'
!        call flush(200+iam)
!        call chem_io_read('inc.'//trim(cdateout)//'.nc',data%aodincgsi, &
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempincaodgsi, &
           trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/', &
           de=de,varname='incaod',rc=localrc)
!        write(200+iam,*)'read aodincgsi',maxval(tempincaodgsi),minval(tempincaodgsi)
!        call flush(200+iam)
        data%aodincgsi=data%aodincgsi+tempincaodgsi
        if(tile==5)then
           write(6,*)'increment aodincgsi',maxval(data%aodincgsi),minval(data%aodincgsi)
        endif
        deallocate(tempincaodgsi)
        if(localrc/=0)then
          write(6,*)'error reading incaodgsi'
          deallocate(data%aodincgsi)
        endif
      endif ! base_path
    endif ! aod

!    if(iam.eq.0.and.tile.eq.1)then
!      write(6,*)'ngsivar',ngsivar
!      write(6,*)'gsiinc',shape(gsiinc),kind(gsiinc)
      !write(6,*)'read ',trim(gsicem(n))//'inc',n
!      call flush(6)
!      write(6,*)'nsecgsi',gsiarray(n)%nsecgsi
!      write(6,*)'read inc.'//trim(cdateout)//'.nc'
!      if(base_path/=' ')then
!        write(6,*)'path '//trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/'
        !call flush(6)
!      else
!        write(6,*)'path ../../GSIINC/'//trim(cdatestart)//'/'
!        call flush(6)
!      endif
!    endif
!    write(200+iam,*)'n',n,'tempgisinc',shape(tempgsiinc)
!    write(200+iam,*)'lb',lbound(tempgsiinc),'ub',ubound(tempgsiinc)
!    write(200+iam,*)'ngsivar',ngsivar
!    call flush(200+iam)
!    write(200+iam,*)'shape gsiinc',shape(gsiinc),kind(gsiinc)
!    call flush(200+iam)
    do nn=1,ngsivar
      if(base_path/=' ')then
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempgsiinc(:,:,:,nn), &
          trim(base_path)//'/GSIINC/'//trim(cdatestart)//'/', &
          de=de,varname='inc'//trim(gsichemp(nn)),rc=localrc)
      else
        call chem_io_read('inc.'//trim(cdateout)//'.nc',tempgsiinc(:,:,:,nn),'../../GSIINC/'//trim(cdatestart)//'/', &
            de=de,varname='inc'//trim(gsichemp(nn)),rc=localrc)
      endif
      if(localrc/=0.0)then
        write(6,*)'localrc error',localrc
       flush(6)
      endif
#ifdef   PERGSIV
      pergsiv(:,:,:,nn)=0.0
#endif
!      write(6,*)'n ',n,gsivarr(n),'nn',nn,'tempgsiinc',maxval(tempgsiinc),minval(tempgsiinc)
!      call flush(6)
      do k=1,nl
        do j=js,je
          jj=j-js+1
          do i=is,ie
            ii=i-is+1
            if(tempgsiinc(i,j,k,nn)==0.0)then
            elseif(tempgsiinc(i,j,k,nn)<0.0)then ! neg
              if(gsicem(n).eq.'co')then
                gsinew=stateout%tr3d(ii,jj,k,loccem(nn))+tempgsiinc(i,j,k,nn)
                if(gsinew<1.e-10)then
                  tempgsiinc(i,j,k,nn)=-(stateout%tr3d(ii,jj,k,loccem(nn))-1.e-10)
                  gsinew=1.e-10
                endif
              elseif(gsicem(n).eq.'no2')then
                gsinew=stateout%tr3d(ii,jj,k,loccem(nn))+tempgsiinc(i,j,k,nn)
                if(gsinew<1.e-15)then
                  tempgsiinc(i,j,k,nn)=-(stateout%tr3d(ii,jj,k,loccem(nn))-1.e-15)
                  gsinew=1.e-15
!                  numneg=numneg+1
                endif
              elseif(gsicem(n).eq.'o3vmr')then
                gsinew=stateout%tr3d(ii,jj,k,p_ox)+tempgsiinc(i,j,k,nn)
                if(gsinew<1.e-11)then
!                 don't want to be lower than 1.e-11 unless input ozone is
!                 less than 1.e-11 but greater than zero
!                 gsinew2=max (1.e-15,min(1.e-11,input ozone))
!                  tempgsiinc(i,j,k,n)=tempgsiinc(i,j,k,n)+(1.e-11-gsinew)
                  tempgsiinc(i,j,k,nn)=-(stateout%tr3d(ii,jj,k,p_ox)-1.e-11)
                  gsinew=1.e-11
                endif
!               need to update ox also
!               gsioxnew=stateout%tr3d(ii,jj,k,p_ox)+tempgsiinc(i,j,k,n)
!               stateout%tr3d(ii,jj,k,p_ox)=gsioxnew
              elseif(gsicem(n).eq.'ch2o')then
                gsinew=stateout%tr3d(ii,jj,k,loccem(nn))+tempgsiinc(i,j,k,nn)
                if(gsinew<1.e-15)then
!                  tempgsiinc(i,j,k,n)=tempgsiinc(i,j,k,n)+(1.e-15-gsinew)
                  tempgsiinc(i,j,k,nn)=-(stateout%tr3d(ii,jj,k,loccem(nn))-1.e-15)
                  gsinew=1.e-15
                endif
              elseif(gsicem(n)=='AOD')then
!               need sub loop
!                write(200+iam,*)'ii',ii,jj,k,'nn',nn,'n',n
!                call flush(200+iam)
!                write(200+iam,*)'loccem',loccem(nn)
!                call flush(200+iam)
!                write(200+iam,*)'lbtemp',lbound(tempgsiinc),'ub',ubound(tempgsiinc)
!                call flush(200+iam)
!                write(200+iam,*)'lctr3d',lbound(stateout%tr3d),'ub',ubound(stateout%tr3d)
!                call flush(200+iam)
    
                gsinew=stateout%tr3d(ii,jj,k,loccem(nn))+tempgsiinc(i,j,k,nn)
!                write(200+iam,*)'gsinew',i,j,k,nn,gsinew
!                call flush(200+iam)
           
                if(aerosol_ugpkg)then
                  if(gsinew<1.e-6)then
                    tempgsiinc(i,j,k,nn)=-(stateout%tr3d(ii,jj,k,loccem(nn))-1.e-6)
                    gsinew=1.e-6
                  endif
                else
                  if(gsinew<1.e-15)then
                    tempgsiinc(i,j,k,nn)=-(stateout%tr3d(ii,jj,k,loccem(nn))-1.e-15)
                    gsinew=1.e-15
                  endif
                endif
              else
                write(6,*)'var not defined yet',n,gsicem(n)
                call flush(6)
                stop
              endif
#ifdef PERGSIV
              pergsiv(i,j,k,n)=tempgsiinc(i,j,k,nn)/max(1.e-15,stateout%tr3d(ii,jj,k,loccem(nn)))*100.
#endif
!              if(abs(pergsiv(i,j,k,n))>200.)then
!                write(6,*)'n',n,'loccem',loccem(n),'gsicem',gsicem(n)
!                write(6,*)'p_ox',p_ox,'gsinew',gsinew
!                write(6,*)'pergsiv',pergsiv(i,j,k,n),'temp',tempgsiinc(i,j,k,n),'orig', &
 
!                stateout%tr3d(ii,jj,k,loccem(n)),stateout%tr3d(ii,jj,k,p_ox)
!              endif
              stateout%tr3d(ii,jj,k,loccem(nn))=gsinew
              gsiinc(i,j,k,nn)=gsiinc(i,j,k,nn)+tempgsiinc(i,j,k,nn)
!              if(gsivarr(n)=='AOD')then
!              write(200+iam,*)'neg gsiincnew',nn,gsiinc(i,j,k,nn)
!              endif
            else ! positive
!              write(200+iam,*)'n',n,'gsiinc',i,j,k,nn
!              call flush(200+iam)
!              write(200+iam,*)'tempgsiinc',tempgsiinc(i,j,k,nn)
!              call flush(200+iam)
!              write(200+iam,*)'gsiinc',gsiinc(i,j,k,nn),'nn',nn,'loccem',loccem(nn)
!              call flush(200+iam)
              if(loccem(nn)>87)then
                write(6,*)iam,'error loccem',nn,loccem(nn),gsivarr(n)
                call flush(6)
                call killit('error 111')
              endif

              gsinew=stateout%tr3d(ii,jj,k,loccem(nn))+tempgsiinc(i,j,k,nn)
              gsiinc(i,j,k,nn)=gsiinc(i,j,k,nn)+tempgsiinc(i,j,k,nn)
!              if(gsivarr(n)=='AOD')then
              !write(200+iam,*)gsivarr(n),'pos gsiincnew',nn,gsiinc(i,j,k,nn)
!              endif
!              write(200+iam,*)'ii',ii,jj,k,'loccem',nn,loccem(nn)
              !call flush(200+iam)
              stateout%tr3d(ii,jj,k,loccem(nn))=gsinew
#ifdef PERGSIV
              pergsiv(i,j,k,nn)=tempgsiinc(i,j,k,nn)/max(1.e-15,stateout%tr3d(ii,jj,k,loccem(nn)))*100.
#endif
            endif ! neg positive
!           need ajl end do nn somewhere need to fix block
!           refix 2/26/2021 ajl
            if(gsicem(n).eq.'o3vmr')then
              incMLS3d(i,j,k)=incMLS3d(i,j,k)+tempincMLS(i,j,k)
              incOMI3d(i,j,k)=incOMI3d(i,j,k)+tempincOMI(i,j,k)
              allocate(o3vmr_inst_inc(is:ie,js:je,nl))
!              write(200+iam,*)'o3vmr_inst in',maxval(o3vmr_inst),minval(o3vmr_inst)
!              write(200+iam,*)'o3vmr_inst_inc',maxval(tempgsiinc(:,:,:,1)), &
!              minval(tempgsiinc(:,:,:,1))

              call flipz(tempgsiinc(:,:,:,1),o3vmr_inst_inc)
              o3vmr_inst=o3vmr_inst+o3vmr_inst_inc
!              write(200+iam,*)'new o3vmr_inst',maxval(o3vmr_inst),minval(o3vmr_inst)
!              call flush(200+iam)
              deallocate(o3vmr_inst_inc)
            else if(gsicem(n).eq.'co')then
              if(allocated(tempincNUCAPS))then
                incNUCAPSco3d(i,j,k)=incNUCAPSco3d(i,j,k)+tempincNUCAPS(i,j,k)
              endif
              if(allocated(tempincTROPOMI))then
                incTROPOMIco3d(i,j,k)=incTROPOMIco3d(i,j,k)+tempincTROPOMI(i,j,k)
              endif
            endif
         
#ifdef PERGSIV
            if(k.eq.1)then
 !             intco(i,j,k)=stateout%tr3d(ii,jj,k,p_co)*dpmgrd(ii,j,k)
              colgsiinc(i,j,k,n)=stateout%tr3d(ii,jj,k,loccem(nn))*dpmgrd(ii,j,k)
            else
!              intco(i,j,k)=intco(i,j,k-1)+stateout%tr3d(ii,jj,k,p_co)*dpmgrd(ii,j,k)
              colgsiinc(i,j,k,n)=colgsiinc(i,j,k-1,n)+stateout%tr3d(ii,jj,k,loccem(nn))*dpmgrd(ii,j,k)
            endif
#endif
           end do ! i
         end do ! j
       end do ! k
!       write(6,*)'gsiinc',nn,gsivarr(n),maxval(gsiinc(:,:,:,nn)),minval(gsiinc(:,:,:,nn))
    end do ! nn
!    endif
!   if aerosols
    if(gsicem(n)=='AOD')then
!     return aerosol back to gdshcem tr3dout
      chem_pass_state%tr3d(:,:,:,p_bc1)=stateout%tr3d(:,:,:,p_bc1)
      chem_pass_state%tr3d(:,:,:,p_bc2)=stateout%tr3d(:,:,:,p_bc2)
      chem_pass_state%tr3d(:,:,:,p_oc1)=stateout%tr3d(:,:,:,p_oc1)
      chem_pass_state%tr3d(:,:,:,p_oc2)=stateout%tr3d(:,:,:,p_oc2)
      chem_pass_state%tr3d(:,:,:,p_sulf)=stateout%tr3d(:,:,:,p_sulf)
      chem_pass_state%tr3d(:,:,:,p_dust1)=stateout%tr3d(:,:,:,p_dust1)
      chem_pass_state%tr3d(:,:,:,p_dust2)=stateout%tr3d(:,:,:,p_dust2)
      chem_pass_state%tr3d(:,:,:,p_dust3)=stateout%tr3d(:,:,:,p_dust3)
      chem_pass_state%tr3d(:,:,:,p_dust4)=stateout%tr3d(:,:,:,p_dust4)
      chem_pass_state%tr3d(:,:,:,p_dust5)=stateout%tr3d(:,:,:,p_dust5)
      chem_pass_state%tr3d(:,:,:,p_seas1)=stateout%tr3d(:,:,:,p_seas1)
      chem_pass_state%tr3d(:,:,:,p_seas2)=stateout%tr3d(:,:,:,p_seas2)
      chem_pass_state%tr3d(:,:,:,p_seas3)=stateout%tr3d(:,:,:,p_seas3)
      chem_pass_state%tr3d(:,:,:,p_seas4)=stateout%tr3d(:,:,:,p_seas4)
      chem_pass_state%tr3d(:,:,:,p_seas5)=stateout%tr3d(:,:,:,p_seas5)
    endif
    deallocate(tempgsiinc)
  end do ! n
!  endif
  first=.false.
  return
  end subroutine chemgsi
end module chemgsimod
