module raqms_model_mod
contains
  subroutine raqms_model_advance(rc)
  use machine,only : kind_phys
  use chem_types_mod
  use raqmschem_state_mod
  use raqmschem_config_mod
  use raqmschem_data_mod
  use raqmschem_model_mod
  use chem_comm_mod, only : chem_comm_get
  use raqmschem_io_mod, only : chem_io_write,chem_io_read
  use raqms_mod
  use calcdxdy, only : initdxdy
  use calcvort_mod, only :calcvort
  use calcvort_dgrid_mod, only :calcvortdgrid
  use field_manager_mod, only : find_field_index
  use raqmschem_iodata_mod, only : raqmschem_chem_write,raqmschem_chem_write_var
  use raqmschem_iodata_mod, only : raqmschem_chem_write_var_2d,raqmschem_chem_write_debug
  use raqmschem_species_mod
  use chem_const_mod, only : grvity
  use raqmschemcomm_mod, only : allocatechemcommperm,pblht,slmsk2d,ustar,cmdrag
  !use raqmschemcomm_mod, only : szgrd,tskin,sfcwind,tgrd,lconvcld,ccthk,zlwigrd,thgrd
  use raqmschemcomm_mod, only : tskin,sfcwind,tgrd,lconvcld,ccthk,zlwigrd,thgrd
  use raqmschemcomm_mod, only : pgrd,ugrd,vgrd,spgrd,zsurf,zgrd
  use raqmschemcomm_mod, only : allocatechemcommtemp,latgrd,longrd,xgrid,ygrid,areagrd
  use raqmschemcomm_mod, only : taucldfrc_liq,taucldfrc_ice,i1grd,j1grd,dist1grd
  use raqmschemcomm_mod, only : deallocatechemcommtemp,o3dep_save,navg,noydep_save,codep_save,o3vmr_inst
  use raqmschem_pmgrid_mod, only : tile,ibeg,iam,iprn,jprn,iamprn,iprnin,jprnin,pecount,kprnin
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
  real(CHEM_KIND_R4) :: af
  real(CHEM_KIND_R4),allocatable,dimension(:,:) :: dz
  character *10 ctest,cstep,filename*100
#define CLDTAUZERO
#ifdef CLDTAUZERO
  logical lcldtauzero
  data lcldtauzero/.false./
  save lcldtauzero
#endif

  real(chem_kind_r4) :: afsfc
  integer shape4(4),nchem
#include <comcdate.h>
  save icall
  data icall/0/
  logical first
  save first
  data first/.true./
  integer nsecsave,ihrsave
  data nsecsave/10800/
  save nsecsave
  call raqmschem_model_get(deCount=deCount, tile=tile,rc=localrc)
  if (deCount < 1) return
  call chem_comm_get(localpe=mype,pecount=pecount)
  call chem_model_clock_get(advanceCount=advanceCount, dts=dts, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s,tz=tz, julday=julday, rc=localrc)
!  if(mype.eq.0)then
    icall=icall+1
!  endif
  if(first)then
    call setchempointers(mype)
!    if(mype.eq.0)then
!    write(6,*)'p_co',p_co,'p_no2',p_no2,'decount',decount
!    call flush(6)
!    endif
    ctest=' '
    call getenv('HRSAVE',ctest)
    if(ctest.ne.' ')then
      read(ctest,*)ihrsave
!      write(6,*)'ihrsave',ihrsave
      nsecsave=float(ihrsave)*3600.
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
  do de = 0, deCount-1
      call raqmschem_model_get(de=de, config=config, data=data, &
        stateIn=stateIn, stateOut=stateOut, rc=localrc)
     call chem_model_domain_get(de=de, ids=is, ide=ie, jds=js, jde=je, ni=ni, nl=nl, &
     lon=lon, lat=lat, rc=localrc,its=its,ite=ite,jts=jts,jte=jte,ims=ims,ime=ime,jms=jms,jme=jme) 
!    ids,ide full limits of this subregion
!    its,ite full limits for this tile
!    ims,ime limits start with index 1
!     write(300+iam,*)'ids',is,ie,'jds',js,je
!     write(300+iam,*)'its',its,ite,'jts',jts,jte
!     write(300+iam,*)'ims',ims,ime,'jms',jms,jme
!     call flush(300+iam)
     ihs=max(is-1,its)
     jhs=max(js-1,jts)
     ihe=min(ie+1,ite)
     jhe=min(je+1,jte)
     if(first)then
!      can calculate griddx,and griddy first time
       call initdxdy(is,ie,js,je,lon,lat,ihs,ihe,jhs,jhe,data%griddx,data%griddy,de,mype)
     endif
    if(.not.allocated(data%vort))then
      allocate(data%vort(is:ie,js:je,nl))
    endif
    call calcvortdgrid(stateIn%us3d,stateIn%vs3d,statein%prl3d,statein%tk3d,statein%area, &
     lat,is,ie,js,je,nl,ihs,ihe,jhs,jhe,de,mype)
!   prl3d is on dashes
    if(.not.allocated(data%absvort))then
       allocate(data%absvort(is:ie,js:je,nl))
    endif
          
        
    
    if(first)then
      latgrd=lat
!      write(300+iam,*)'latgrd ',lbound(latgrd),'u',ubound(latgrd)
!      write(300+iam,*)'lat',lbound(lat),'u',ubound(lat)
      longrd=lon
      xgrid=lon
      ygrid=lat
      i1grd=nint(lon)
      j1grd=nint(lat)
      dist1grd=sqrt((i1grd-lon)**2.+(j1grd-lat)**2.)
#ifdef DIAGFV3
      write(6,*)'i1grd',maxval(i1grd),minval(i1grd)
      write(6,*)'j1grd',maxval(j1grd),minval(j1grd)
#endif
!      write(6,*)'dist1grd',maxval(dist1grd),minval(dist1grd)
      data%area=statein%area
      areagrd=statein%area
!     will have to read in since ph3d 1 is zero
!      szgrd=statein%ph3d(:,:,1)/grvity
      slmsk2d=statein%slmsk2d
    endif
    call allocatechemcommtemp
    pblht=statein%pb2d
    spgrd=statein%pr3d(:,:,1)*.01 ! mb
!    write(6,*)'ph3d 1',maxval(statein%ph3d(:,:,1))
!    zsurf=statein%phl3d(:,:,1)/9.80616 ! meters
    ustar=statein%us2d
    cmdrag=statein%cmdrag
    tskin=statein%ts2d
    tgrd=statein%tk3d
    ugrd=statein%us3d
    vgrd=statein%vs3d
    pgrd=statein%prl3d*.01 ! mb
    thgrd(:,:,:)=tgrd(:,:,:)*(1000./pgrd(:,:,:))**.286
    do k=1,nl
!     zgrd is dashes
      zgrd(:,:,k)=statein%phl3d(:,:,k)/9.80616+zsurf(:,:)
    end do
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
   
!    if(mype.eq.0.and.first)then
!      write(6,*)'pbhlht',pblht(1,js),'sz',zsurf(1,js),'sp',spgrd(1,js)
!      do k=1,nl
!        write(6,*)'u',k,ugrd(1,js,k),vgrd(1,js,k),'t',tgrd(1,js,k)
!        write(6,*)'thgrd',thgrd(1,js,k),'p',pgrd(1,js,k)
!        write(6,*)'zgrd',zgrd(1,js,k)
!      !end do
!    endif
!    write(6,*)'tk3d',lbound(statein%tk3d),' ub ',ubound(statein%tk3d)
!    call flush(6)
!    write(6,*)'tgrd',lbound(tgrd),' ub ',ubound(tgrd)
!    call flush(6)
!    write(6,*)'cldtau_frac_liq',maxval(statein%cldtau_frac_liq),minval(statein%cldtau_frac_liq),shape(statein%cldtau_frac_liq)
!    write(6,*)'cldtau_frac_ice',maxval(statein%cldtau_frac_ice),minval(statein%cldtau_frac_ice),shape(statein%cldtau_frac_ice)
#if 0
!    if(mype.eq.7)then
      do j=js,je
        jj=j-js+1
        do i=1,ie-is+1
          ii=is+i-1
          if(j.eq.80.and.ii.eq.54)then
             write(6,*)'mype',mype
            call flush(6)
            do k=1,nl
              write(6,*)'tin',ii,j,statein%tk3d(i,jj,k),'pin',statein%prl3d(i,jj,k),'zin',statein%phl3d(i,jj,k)
              call flush(6)
            end do
          endif
        end do
      end do
!    endif
#endif
!    write(6,*)'slmsk2d',maxval(statein%slmsk2d),minval(statein%slmsk2d)
!    call flush(6)
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
          ccthk(i,j)=(statein%ph3d(i,jj,ktop)-statein%ph3d(i,jj,kbot))/9.80616/1000.
          lconvcld(i,j)=ktop
        else
          ccthk(i,j)=0.0
        endif
        if(statein%slmsk2d(i,jj)<.001)then
!         ocean
          if(statein%ph3d(i,jj,1)<.0001)then
!            if(statein%ts2d(i,jj)<273.0)then
!              zlwigrd(i,j)=3. ! make sea ice for now
!              write(6,*)'sea ice ',i,j,'ts',statein%ts2d(i,jj),'ph',statein%ph3d(i,jj,1)
!              call flush(6)
!            else
            zlwigrd(i,j)=1. ! ocean
!            endif
          else
            zlwigrd(i,j)=2. ! make lake land
          endif
        else
          if(nint(statein%slmsk2d(i,jj)).eq.2)then
!            write(6,*)'sea ice at ',ii,j,'sz',statein%ph3d(i,jj,1),'ts',statein%ts2d(i,jj)
!            call flush(6)
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
    if(mype.eq.0)then
!    write(6,*)'shape ph3d',shape(statein%ph3d)
!    write(6,*)'shape phl3d',shape(statein%phl3d)
     write(6,*)'grvity',grvity
     do k=1,nl+1
!       write(6,*)'ph3d',statein%ph3d(1,1,k),'z',statein%ph3d(1,1,k)/grvity
!       if(k/=nl+1)then
!          write(6,*)'dz',k,(statein%ph3d(1,1,k+1)-statein%ph3d(1,1,k))/grvity,'p',statein%pr3d(1,1,k)
!       endif
     end do
    endif
    dz=(statein%ph3d(:,:,2)-statein%ph3d(:,:,1))/grvity*100. ! make cm
!    write(6,*)'dz',maxval(dz)*.01,minval(dz)*.01
#ifndef CLDTAUZERO
    taucldfrc_liq=statein%cldtau_frac_liq
    taucldfrc_ice=statein%cldtau_frac_ice
#else
    if(lcldtauzero)then
      taucldfrc_liq=0.
      taucldfrc_ice=0.
    else
      taucldfrc_liq=statein%cldtau_frac_liq
      taucldfrc_ice=statein%cldtau_frac_ice
    endif
#endif
!    write(6,*)'cmdrag',maxval(statein%cmdrag),minval(statein%cmdrag)
!    write(6,*)'phl3d bot',maxval(statein%phl3d(:,:,1)),minval(statein%phl3d(:,:,1))
!    call flush(6)
!    write(6,*)'ph3d bot',maxval(statein%ph3d(:,:,1)),minval(statein%ph3d(:,:,1))
!    call flush(6)
    if(mype.eq.0)then
      write(6,*)'call advance cdate',cdate,'h',h,'m',m,'s',s,'icall',icall,'advancecount',advancecount
      call flush(6)
    endif
!    call wetdep(de, dts, Statein%pb2d, Statein%rc2d, Statein%rn2d, Statein%ph3d, Statein%phl3d, &
!      Statein%pr3d, Statein%prl3d, Statein%tk3d, Statein%us3d, Statein%vs3d, Statein%ws3d, &
!      Statein%tr3d, Stateout%tr3d, wet_dep, ntra, numgas, num_chem, num_moist, &
!      its, ite, jts, jte, kte, kts, &
!      ims, ime, jms, jme, kme, kms, rc)
!    write(6,*)'rn2d',maxval(statein%rn2d),'rc2d',maxval(statein%rc2d)
!    call flush(6)
    shape4=shape(statein%tr3d)
    nchem=shape4(4)
!    if(mype.eq.0)then
!       write(6,*)'ntra',config%ntra,' nbegin ',config%nbegin, 'numgas ',config%numgas,' num_moist',config%num_moist
!    endif
    iprn=0
    jprn=0
!    write(6,*)'pecount',pecount
!    call flush(6)
!    if(pecount.eq.24)then
!    iprnin=74
!    jprnin=32
!    else
!    iprnin=45
!    jprnin=62
     iprnin=63
     jprnin=57
     iprnin=52
     jprnin=49
     iprnin=80
     jprnin=33
     kprnin=1
!    endif
    iamprn=-1
    if(tile.eq.2)then
      write(6,*)'is',is,js,js,je
      if(is<=iprnin.and.ie>=iprnin.and.js<=jprnin.and.je>=jprnin)then
        write(250+iam,*)'is',is,ie,'js',js,je
        iprn=iprnin-is+1
        jprn=jprnin-js+1
        iamprn=iam
        write(300+iam,*)'iamprn',iamprn,'iprn',iprn,jprn,'is',is,'js',js
        call flush(300+iam)
        write(250+iam,*)'iprn',iprn,jprn
        call flush(6)

      endif
    endif
    if(iprn>0)then
!      do k=1,63
      k=1
        write(250+iam,*)k,'brcladvbefore',statein%tr3d(iprn,jprn,k,p_brcl)
        write(250+iam,*)k,'bry',statein%tr3d(iprn,jprn,k,p_bry)
        write(300+iam,*)'lat iprn,jprn ',iprn,jprn,lat(iprn,jprn),lon(iprn,jprn)
        write(300+iam,*)'latgrd',iprn,'jprnin',jprnin,latgrd(iprn,jprnin)
        write(300+iam,*)'iamprn',iamprn
        write(300+iam,*)k,'brcladvbefore',statein%tr3d(iprn,jprn,k,p_brcl)
        write(300+iam,*)k,'bry',statein%tr3d(iprn,jprn,k,p_bry)
!      end do
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
    is, ie, js, je, 1, nl, &
    data, dz,&
    rc=localrc)
    if(iprn>0)then
!      do k=1,63
      k=1
        write(300+iam,*)k,'brclaft',stateout%tr3d(iprn,jprn,k,p_brcl)
        write(300+iam,*)k,'bry',stateout%tr3d(iprn,jprn,k,p_bry)
        write(6,*)k,'brclaft',stateout%tr3d(iprn,jprn,k,p_brcl)
        write(300+iam,*)'maxval',maxval(stateout%tr3d(:,:,k,p_brcl))
        call flush(300+iam)
!      end do
    endif
!   put at hour
    write(cstep,'(i4.4)')icall
    filename='each'//trim(cstep)//'.nc'
!    call raqmschem_chem_write_debug(stateout%tr3d,dims,filename=filename,rc=rc)
    dims(1)=ie-is+1
    dims(2)=je-js+1
    dims(3)=nl
    dims2=dims(1:2)
    call raqmschem_chem_write(stateout%tr3d,dims,filename=filename,itime=1,cdate=cdate,rc=rc)

    if(mod(nint(dts*float(icall)),nsecsave).eq.0)then
!    if(mod(icall,6).eq.0)then
!      isave=icall/6
      isave=1
      dims(1)=ie-is+1
      dims(2)=je-js+1
      dims(3)=nl
      dims2=dims(1:2)
      if(mype.eq.0)then
        write(6,*)'cdate',cdate,'h',h,'m',m,'s',s,'icall',icall
        write(6,*)'call raqmschem_chem_write ',cdate
        call flush(6)
        write(6,*)'output tracers'
        call flush(6)
      endif
!      call raqmschem_chem_write_var_2d(o3dep_save,dims2,isave,'o3dep',cdate=cdate,rc=rc)
 
      call raqmschem_chem_write(stateout%tr3d,dims,itime=isave,cdate=cdate,rc=rc)
!      call raqmschem_chem_write_var(thgrd,dims,isave,'theta',cdate=cdate,rc=rc)
!      call raqmschem_chem_write_var(o3vmr_inst,dims,isave,'o3vmr',cdate=cdate,rc=rc)
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
    endif
!    if(first)then
!      dz=dz*.01 ! makd meters
!      call chem_io_write('slmsk2d.nc',slmsk2d,path=trim(config%emi_outname),de=de,rc=rc)
!      call chem_io_write('zlwi.nc',zlwigrd,path=trim(config%emi_outname),de=de,rc=rc)
!      call chem_io_write('dz.nc',dz,path=trim(config%emi_outname),de=de,rc=rc)
       
!    endif
!    if(first)then
!    call chem_io_write('o3dep.nc',o3dep_save,path=trim(config%emi_outname),de=de,rc=rc)
!    endif
  end do
  call deallocatechemcommtemp
  first=.false.
  return
  end subroutine raqms_model_advance
end module raqms_model_mod
