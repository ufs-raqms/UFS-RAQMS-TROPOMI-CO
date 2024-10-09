      real*8 drip,dens,xbr,xbro,xhbr,xbrno3,xhobr,xbrcl,xbr2,sice_br2,ss_br2
      integer kdn,ktf,ihet_flag
      call brintg(kdn,ktf,ihet_flag,drip,dens,xbr,xbro,
     1  xhbr,xbrno3,xhobr,xbrcl,xbr2,sice_br2,ss_br2)
      stop 
      endif
      subroutine brintg (kdn,ktf,ihet_flag,drip,dens,xbr,xbro,
     1  xhbr,xbrno3,xhobr,xbrcl,xbr2,sice_br2,ss_br2)
      use raqmschem_pmgrid_mod, only : iam,iatchem,jatchem,katchem,
     *nstepat,iprn,jprn,kprn
c
      include 'comm_defs_24t.chem'
      common /blk4_24/ r(ireac),hr(ihreac),phoj(iphot),dn(inum),
     1              zeta,pur,sad,prodloss(13)
      common/pceinfo/timinc,tdenom,save_clo,save_bro,
     +       ipce_clno3,ipce_hocl,ipce_brno3,ipce_hobr,rbrobr
     +       ,equil_fam
      common/brx_com/pjbr2,rbr_ss,fr_hbrss,tau_ss,dvh,fctr_br
      real xhbrold,xbrno3old,xhobrold,xbrclold,xbr2old,xbroold
      real dtxhbr,dtxbrno3,dtxhobr,dtxbrcl,dtxbr2,dtxbro
      real xbrold,dtxbr


c
C$OMP THREADPRIVATE (/blk4_24/)
C$OMP THREADPRIVATE (/pceinfo/)
C$OMP THREADPRIVATE (/brx_com/)
c
c
      t = 0.
      n = 0
      dt = 0.2
c
cccccc  DAYTIME
cc
      if (kdn .eq. 1) then
!       if(iprn>0.and.kprn.eq.1)then
!         write(300+iam,*)'drip',drip,'xbr',xbr,'xbro',xbro
!         write(300+iam,*)'xhbr',xhbr,'xbrno3',xbrno3,'xhobr',xhobr
!         write(300+iam,*)'xbrcl',xbrcl,'xbr2',xbr2
!         write(300+iam,*)'sice_br2',sice_br2,'ss_br2',ss_br2
!         write(300+iam,*)'pjbr2',pjbr2,'rbr_ss',rbr_ss,'fr_hbrss',fr_hbrss
!         write(300+iam,*)'tau_ss',tau_ss,'dvh',dvh,'fctr_br',fctr_br
!         call flush(300+iam)
!       endif
10      continue
        t = t+dt
        n = n+1
        if (t .gt. timinc) then
          t = t - dt
          dt = timinc+0.1-t
          t = timinc + 0.1
        endif
        xhbrold=xhbr
        xbrno3old=xbrno3
        xhobrold=xhobr
        xbrclold=xbrcl
        xbr2old=xbr2
        xbroold=xbro
        xbrold=xbr

!      if(n>10)then
!        write(200+iam,*)'n',n,'t',t,'dt',dt,'timinc',timinc,'kind',kind(timinc),
!     *  'i',iatchem,jatchem,katchem,'xhobr',xhobr,'xbro',xbro,'xbr2',xbr2,
!     *   'xbrno3',xbrno3,'tau_ss',tau_ss,'fr_hbrss',fr_hbrss,'drip',
!     *   drip,'xbr',xbr,'xhbr',xhbr,'dvn',dvn,'sice_br2',sice_br2
!        call flush(200+iam)
!      endif
!     if(isnan(xhobr))then
!       write(6,*)'top xhobr ',xhobr,'n',n,'t',t,'dt',dt,iatchem,jatchem,katchem
!       call flush(6)
!     endif
c
      pbr = (phoj(32)+r(97)*dn(2)+r(99)*dn(3)+
     1       (r(101)+r(102))*dn(25)
     1       +2.*0.85*r(104)*xbro)*xbro + r(105)*dn(9)*xhbr
     1       +phoj(31)*xbrcl + phoj(34)*xhobr
     1       +0.29*phoj(35)*xbrno3+ 2.0*pjbr2*xbr2
      dbr = r(94)*dn(1)+r(95)*dn(10)+r(96)*dn(19)
c
      pbro = r(94)*xbr*dn(1)+(r(63)*xbrcl+r(89)*xhobr)*dn(2)
     1       +0.71*phoj(35)*xbrno3
#ifdef DIAGXBR
      if(pbro>2.e8)then
        write(200+iam,*)n,'xbr',xbr,'xbrcl',xbrcl,'xhobr',xhobr,'xbrno3',
     *  xbrno3,iatchem,jatchem,katchem
        write(200+iam,*)n,'r94',r(94),'dn1',dn(1),'r63',r(63),'r89',r(89),'dn2',dn(2)
        write(200+iam,*)'phoj35',phoj(35),'term1 ',r(94)*xbr*dn(1),' term2 ',
     *  r(63)*xbrcl*dn(2),' term3 ',r(89)*xhobr*dn(2)
        call flush(200+iam)
      endif
#endif
      dbro = phoj(32)+r(97)*dn(2)+r(98)*dn(10)+r(99)*dn(3)+
     1       r(100)*dn(4)*dens+(r(101)+r(102)+r(103))*dn(25)
     1       +2.*r(104)*xbro
      phbr = (r(95)*dn(10)+r(96)*dn(19))*xbr 
      dhbr = phoj(33)+r(105)*dn(9)+r(88)*dn(2)+hr(6)*xhobr+hr(8)*dn(38)
     1       +drip+fr_hbrss/tau_ss
c
      phobr = r(98)*xbro*dn(10)+hr(9)*xbrno3*dn(13)
      dhobr = phoj(34)+r(89)*dn(2)+hr(6)*xhbr+hr(7)*dn(14)
CCCMNN0508    loss of hobr on seaice
      dhobr = dhobr + dvh*fctr_br
c
      pbrno3 = r(100)*dn(4)*dens*xbro
      dbrno3 = phoj(35)+hr(9)*dn(13)
c
      pbrcl = r(103)*dn(25)*xbro+hr(7)*dn(14)*xhobr+
     1        hr(8)*dn(38)*xhbr
      dbrcl = phoj(31)+r(63)*dn(2)
c
      pbr2 = r(104)*0.15*xbro*xbro+hr(6)*xhobr*xhbr
      if(abs(pbr2)>1.e40)then
        write(200+iam,*)'pbr2',pbr2,iatchem,jatchem,katchem
        call flush(200+iam)
        write(200+iam,*)'xbro',xbro,'hr6',hr(6),'xhobr',xhobr,'xhbr',xhbr
        call flush(200+iam)
        write(200+iam,*)'t',t,'n',n,'timinc',timinc,'dt',dt
        call flush(200+iam)
      endif
!     if(isnan(pbr2))then
!        write(200+iam,*)'pbr2 aa ',pbr2,'t',t,n
!        call flush(200+iam)
!     endif
!CCMNN0508    Br2 source at seaice
      pbr2 = pbr2 + dvh*fctr_br*xhobr
!     if(abs(pbr2)>1.e40)then
!         write(200+iam,*)'iatchem',iatchem,jatchem,katchem
!         write(200+iam,*)'pbr2 b big ',pbr2,'dvh',dvh,fctr_br,xhobr
!         call flush(200+iam)
!     endif

!     if(isnan(pbr2))then
!        write(200+iam,*)'pbr2 bb ',pbr2,'t',t,n
!        call flush(200+iam)
!     endif
!CCMNN0508    Br2 from seasalt within the PBL
      if (xhobr .gt. dens*1.0e-15) then
       pbr2 = pbr2 + rbr_ss
!     if(abs(pbr2)>1.e40)then
!         write(200+iam,*)'iatchem',iatchem,jatchem,katchem
!       write(200+iam,*)'pbr2 c big ',pbr2,'rbr_ss',rbr_ss
!       call flush(200+iam)
!     endif
       ss_br2 = ss_br2 + rbr_ss*dt/timinc
       endif
      dbr2 = pjbr2
c
      sice_br2 = sice_br2+dvh*fctr_br*xhobr*dt/timinc
!     if(isnan(sice_br2).or.abs(sice_br2)>1.e40)then
!        write(6,*)'at b dvh',dvh,'fctr_br',fctr_br,iatchem,jatchem,
!    *   katchem,'xhobr',xhobr,'dt',dt,'timinc',timinc,'t',t,n
!     endif
c
CCCMNN 0608
      dtinv = 1.0/dt
cc      xbr = xbr+petend(xbr,pbr,dbr,dt,dtinv)*dt
cc      xbro = xbro+petend(xbro,pbro,dbro,dt,dtinv)*dt
cc      xbrcl = xbrcl+petend(xbrcl,pbrcl,dbrcl,dt,dtinv)*dt
cc      xbr2 = xbr2+petend(xbr2,pbr2,dbr2,dt,dtinv)*dt
cc      xhbr = xhbr+petend(xhbr,phbr,dhbr,dt,dtinv)*dt
cc      xhobr = xhobr+petend(xhobr,phobr,dhobr,dt,dtinv)*dt
CCCMNN  Limit brno3 production by available NO2
cc      xbrno3 = xbrno3+min(petend(xbrno3,pbrno3,dbrno3,dt,dtinv),
cc     1         dn(4)/timinc)*dt
cc      xbrno3 = xbrno3+petend(xbrno3,pbrno3,dbrno3,dt,dtinv)*dt


      xbr=pbr/dbr + (xbr - pbr/dbr)*exp(-dbr*dt)
      xbro=pbro/dbro + (xbro - pbro/dbro)*exp(-dbro*dt)
      xbrcl=pbrcl/dbrcl + (xbrcl-pbrcl/dbrcl)*exp(-dbrcl*dt)
!     if(isnan(xbr2))then
!         write(200+iam,*)'iatchem',iatchem,jatchem,katchem
!      write(200+iam,*)'abover br2'
!      call flush(200+iam)
!     endif
      xbr2=pbr2/dbr2 + (xbr2-pbr2/dbr2)*exp(-dbr2*dt)
!     if(isnan(xbr2))then
!         write(200+iam,*)'iatchem',iatchem,jatchem,katchem
!       write(200+iam,*)' bb xbr2',xbr2,'pbr2',pbr2,'dbr2',dbr2,'dt',dt,'xbr',xbr,'dbr',dbr
!       call flush(200+iam)
!     endif
      xhbr=phbr/dhbr + (xhbr - phbr/dhbr)*exp(-dhbr*dt)
      xhobr=phobr/dhobr + (xhobr - phobr/dhobr)*exp(-dhobr*dt)
      xtmp=pbrno3/dbrno3 + (xbrno3 - pbrno3/dbrno3)*exp(-dbrno3*dt)
      if ((xtmp-xbrno3)/dt .gt. dn(4)/timinc) then
       xbrno3=xbrno3+dn(4)*dt/timinc
       else
       xbrno3 = xtmp
       endif

c
      if (t .lt. timinc) then
       dt = dt*1.2
       dt = min(dt,10.)
       if(n>25)then
!      real dtxhbr,dtxbrno3,dtxhobr,dtxbrcl,dtxbr2,dtxbro
       dtxbr=abs(xbr-xbrold)/max(1.e-20,xbrold)
       dtxhbr=abs(xhbr-xhbrold)/max(1.e-20,xhbrold)
       dtxbro=abs(xbro-xbroold)/max(1.e-20,xbroold)
       dtxbrcl=abs(xbrcl-xbrclold)/max(1.e-20,xbrclold)
       dtxbr2=abs(xbr2-xbr2old)/max(1.e-20,xbr2old)
!       write(500+iam,*)'xhobr',xhobr,'old',xhobrold
!       write(500+iam,*)'xbrno3',xbrno3,xbrno3old
!       call flush(500+iam)
       dtxhobr=abs(xhobr-xhobrold)/max(1.e-20,xhobrold)
       dtxbrno3=abs(xbrno3-xbrno3old)/max(1.e-20,xbrno3old)
!       write(500+iam,*)'dtxhobr',dtxhobr,dtxbrno3,iatchem,jatchem,katchem
!       call flush(500+iam)
!      real xhbrold,xbrno3old,xhobrold,xbrclold,xbr2old,xbroold
       if(dtxhbr>.90.or.dtxbro>.90.or.dtxbrcl>.90.or.dtxbr2>.90.or.
     *    dtxhobr>.90.or.dtxbrno3>.90.or.dtxbr>.90)then
          dt=dt*.5
#ifdef DIAGXBR
         write(200+iam,*)'phobr/dhobr',phobr/dhobr,'xhobr',xhobr,'dhobr',dhobr
         write(200+iam,*)'nstepat',nstepat,' dt ',dt,'dtxbr',dtxbr
         write(200+iam,*)n,dt,'iatchem',iatchem,jatchem,katchem,
     *   ' dtxhbr ',dtxhbr,' dtxbro ',dtxbro,
    
     *   ' dtxbrcl ',dtxbrcl,' dtxbr2 ',dtxbr2,' xhobr ',
     *   dtxhobr,' dtbrno3 ',dtxbrno3
         write(200+iam,*)'xhobrold ',xhobr,xhobrold,
     *  ' xbrno3old ',xbrno3,xbrno3old
         call flush(200+iam)
!      xbro=pbro/dbro + (xbro - pbro/dbro)*exp(-dbro*dt)
        write(200+iam,*)'n',n,'xbro',xbro,'pbro',pbro,'dbro',dbro
        write(200+iam,*)'ratio ',pbro/dbro,'exp',-dbro*dt
        write(200+iam,*)'xbrno3',xbrno3
        write(200+iam,*)'xbr',xbr,'xhobr',xhobr
        call flush(200+iam)
#endif
       endif
       endif
       go to 10
       endif
c
c
cccccc     NIGHTTIME
      else
c
20    continue
      t = t+dt
      n = n+1
      if (t .gt. timinc) then
       t = t - dt
       dt = timinc+0.1-t
       t = timinc + 0.1
       endif
c
      pbr = (r(99)*dn(3)+
     1       (r(101)+r(102))*dn(25)
     1       +2.*0.85*r(104)*xbro)*xbro + r(105)*dn(9)*xhbr
      dbr = r(94)*dn(1)+r(95)*dn(10)+r(96)*dn(19)
c
      pbro = r(94)*xbr*dn(1)
      dbro = r(98)*dn(10)+r(99)*dn(3)+
     1       r(100)*dn(4)*dens+(r(101)+r(102)+r(103))*dn(25)
     1       +2.*r(104)*xbro
      phbr = (r(95)*dn(10)+r(96)*dn(19))*xbr
      dhbr = r(105)*dn(9)+hr(6)*xhobr+hr(8)*dn(38)
     1     +drip + fr_hbrss/tau_ss
c
c
      phobr = r(98)*xbro*dn(10)+hr(9)*xbrno3*dn(13)
      dhobr = hr(6)*xhbr+hr(7)*dn(14)
CCCMNN0508   loss on seaice
      dhobr = dhobr + dvh*fctr_br
c
      pbrno3 = r(100)*dn(4)*dens*xbro
      dbrno3 = hr(9)*dn(13)
c
      pbr2 = r(104)*0.15*xbro*xbro+hr(6)*xhbr*xhobr
CCCMNN0508   source on seaice
      pbr2 = pbr2 + dvh*fctr_br*xhobr
CCCMNN0508   Br2 from seasalt
      if (xhobr .gt. dens*1.0e-15) then
       pbr2  = pbr2 + rbr_ss
       ss_br2 = ss_br2 + rbr_ss*dt/timinc
       endif
      dbr2=0.0
c
      pbrcl = r(103)*dn(25)*xbro+hr(7)*dn(14)*xhobr+
     1        hr(8)*dn(38)*xhbr
      dbrcl = 0.0
c
      sice_br2 = sice_br2+dvh*fctr_br*xhobr*dt/timinc
!     if(isnan(sice_br2).or.abs(sice_br2)>1.e40)then
!        write(6,*)'at a dvh',iatchem,jatchem,katchem,dvh,'fctr_br',
!    *fctr_br,'xhobr',xhobr,'dt',dt,'timinc',timinc,'t',t,'n',n
!        call flush(6)
!     endif
c
CCCMNN 0608
      dtinv = 1.0/dt
cc      xbr = xbr+petend(xbr,pbr,dbr,dt,dtinv)*dt
cc      xbro = xbro+petend(xbro,pbro,dbro,dt,dtinv)*dt
cc      xbrcl = xbrcl+petend(xbrcl,pbrcl,dbrcl,dt,dtinv)*dt
cc      xbr2 = xbr2+petend(xbr2,pbr2,dbr2,dt,dtinv)*dt
cc      xhbr = xhbr+petend(xhbr,phbr,dhbr,dt,dtinv)*dt
cc      xhobr = xhobr+petend(xhobr,phobr,dhobr,dt,dtinv)*dt
CCCMNN  Limit brno3 production by available NO2
cc      xbrno3 = xbrno3+min(petend(xbrno3,pbrno3,dbrno3,dt,dtinv),
cc     1         dn(4)/timinc)*dt
cc      xbrno3 = xbrno3+petend(xbrno3,pbrno3,dbrno3,dt,dtinv)*dt


      xbr=pbr/dbr + (xbr - pbr/dbr)*exp(-dbr*dt)
      xbro=pbro/dbro + (xbro - pbro/dbro)*exp(-dbro*dt)
      if (dhbr .gt. 0.0) then
       xhbr=phbr/dhbr + (xhbr - phbr/dhbr)*exp(-dhbr*dt)
       else
       xhbr = xhbr + phbr*dt
       endif
      if (dhobr .gt. 0.0) then
       xhobr=phobr/dhobr + (xhobr - phobr/dhobr)*exp(-dhobr*dt)
       else
       xhobr = xhobr + phobr*dt
       endif
      if (dbrno3 .gt. 0.0) then
cc       xbrno3=pbrno3/dbrno3 + (xbrno3 - pbrno3/dbrno3)*exp(-dbrno3*dt)
cc       else
cc       xbrno3 = xbrno3 + pbrno3*dt
cc       endif
       xtmp=pbrno3/dbrno3 + (xbrno3 - pbrno3/dbrno3)*exp(-dbrno3*dt)
       else
       xtmp = xbrno3 + pbrno3*dt
       endif
      if ((xtmp-xbrno3)/dt .gt. dn(4)/timinc) then
       xbrno3=xbrno3+dn(4)*dt/timinc
       else
       xbrno3=xtmp
       endif

      xbrcl = xbrcl + pbrcl * dt
      xbr2 = xbr2 + pbr2 * dt

c
      if (t .lt. timinc) then
       dt = dt*1.2
       dt = min(dt,10.)
       go to 20
       endif
c
c
      endif
      return
      end
