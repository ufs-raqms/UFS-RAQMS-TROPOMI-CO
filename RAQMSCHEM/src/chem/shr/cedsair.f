
    subroutine definecedsair
    use raqmschem_pmgrid_mod, only : iam,strip
    character *256 :: ccedsairin
    integer is(20),ie(20),k,igot
    call getenv('CCEDSAIR',ccedsairin)
    call strip(ccedsairin,is,ie,9,igot)
    numcedsair=igot
    do k=1,igot
      ccedsair(k)=ccedsairin(is(k):ie(k))
      select caSe (ccedsair(k))
      case ('CO')
        p_ceds_co=k
      case ('NOx')
        p_ceds_nox=k
      case ('SO2')
        p_ceds_so2=k
      case ('BC')
        p_ceds_bc=k
      case ('OC')
        p_ceds_oc=k
      case ('N2O')
        p_ceds_n2o=k
      caSe ('NH3')
        p_ceds_nh3=k
      caSe ('CH4')
        p_ceds_ch4=k
      case default
        write(6,*)'bad ceds cehm',ccedsair(k)
        call killit('stop cedsin')
      end select
    end do
    if(iam==0)then
      write(6,*)'numcedsair',numcedsair
      write(6,*)'P_ceds_co',p_ceds_co,'nox',p_ceds_nox,'so2',p_ceds_so2
      do k=1,numcedsair
        write(6,*)'cedsair',k,ccedsair(k)
      end do
    endif

    end subroutine definecedsair
      subroutine docedsair(i,j,ip,lbot,ltb,dens,zgeoh,sacco,sacnox,dt,chemlocal)
      use raqmschemlocaltype_mod
      use raqmschem_pmgrid_mod, only : iam,tile,ibeg,beglat
!      use raqmschemcomm_mod, only : intcedsair,ccedsair,p_ceds_co,p_ceds_nox,p_ceds_so2
!      use raqmschemcomm_mod, only : p_ceds_bc,p_ceds_oc,p_ceds_n2o,p_ceds_ch4
!      use raqmschemcomm_mod, only : numcedsair,ktopcedsair
!      use raqmschemcomm_mod, only : srccocedsair,srcso2cedsair
!      use raqmschemcomm_mod, only : lcedsair,intaltair
      implicit none
      type(chemlocaltype),target :: chemlocal


      integer lbot,klo,kup,k,i,j,ip,m,l,ltb(:)
      real dens,sacco,sacnox,zgeoh(:),aflo,afup
      real*4 dt
      real intvallo,intvalup,zlo,zup,sacso2,sacbc,sacoc,sacnh3,sacn2o
      logical :: zloset,zupset
#include <chemlocaldefinepointer.h>
#include <chemlocaldefinepointer2.h>
#include <chemlocalsetpointer.h>
#include <chemlocalsetpointer2.h>
!      write(6,*)'cogrd in commmon cogrd',shape(cogrd)
!      write(6,*)'lb cogrd',lbound(cogrd),'ub',ubound(cogrd)
      zloset=.false.
      zupset=.false.
!     need to find depth of layer 1 in altitude
      zlo = zgeoh(ip)-zgeoh(lbot)              !height above surface (km)
      klo=-1
      kup=-1
!     if(tile==5.and.ibeg+i-1==168.and.j==192)then
!       write(6,*)'ip',ip,'lbot',lbot,'ktopcessair',ktopcedsair
!       write(6,*)'zlo',zlo,'intaltair',intaltair(ktopcedsair)
!     endif
!      write(6,*)'ktopcedsair',ktopcedsair
!      flush(6)
!      write(6,*)'allocated',allocated(intaltair)
!      flush(6)
      !if(iam==75)then
!      write(6,*)'newktopcedsair chemmoc ',maxval(ktopcedsair(:,:,1)),'shape intaltair',shape(intaltair)
!      call flush(6)
      !endif
      if(zlo<intaltair(maxval(ktopcedsair(i,j,:))).and.ip/=1)then
        zup = zgeoh(ip-1)-zgeoh(lbot)           !height above surface (km)
!       if(tile==5.and.ibeg+i-1==168.and.j==192)then
!             write(6,*)ip,'zlo',zlo,zup
!             call flush(6)
!         do k=1,ktopcedsair
!           write(6,*)'intaltair',k,intaltair(k),'intcoair',intcoair(i,j,k)
!        end do
!           endif
           if(zup>=intaltair(maxval(ktopcedsair(i,j,:))))then
             kup=ktopcedsair(i,j,1)
             afup=1.
             zupset=.true.
           endif
           do k=1,ktopcedsair(i,j,1)
             if(zlo<intaltair(k).and.zlo>=intaltair(k-1))then
               klo=k-1
               aflo=(zlo-intaltair(k-1))/(intaltair(k)-intaltair(k-1))
!             if(iam==75.and.i==1.and.j==49)then
               !write(6,*)'newsetklo',klo,'k',k,'zlo',zlo
!               write(6,*)'intaltair',k-1,k,intaltair(k-1:k)
!               write(6,*)'aflo',aflo
!             endif
!         if(tile==5.and.i==168.and.j==192)then
!       if(tile==5.and.ibeg+i-1==168.and.j==192)then
!               if(coiplo/=0.0)then
!               write(6,*)'iam',iam,i,j,'ip',ip,'k1',k1
               !write(6,*)'zlo',zlo,'af',af,'coiplo',coiplo
!               write(6,*)'intaltair',k-1,intaltair,k,intaltair(k)
               !call flush(6)
!               endif
               zloset=.true.
             endif
             if(zup<intaltair(k).and.zup>=intaltair(k-1))then
               kup=k-1
               afup=(zup-intaltair(k-1))/(intaltair(k)-intaltair(k-1))
!         if(tile==5.and.i==168.and.j==192)then
!       if(tile==5.and.ibeg+i-1==168.and.j==192)then
!               if(coiphi/=0.0)then
!               write(6,*)'iam',iam,i,j,'k1',k1
!               write(6,*)'zup',zup,'af',af,'coiphi',coiphi
!               write(6,*)'intaltair',k-1,intaltair,k,intaltair(k)
!               call flush(6)
!               endif
               zupset=.true.
             endif
             if(zupset.and.zloset)exit
           end do ! k=1,ktop
          if(klo<0)then
!             if(klo<0.and.iam==75.and.i==1.and.j==49)then
!             !write(6,*)'new klo',klo,'i',i,j,ip,'zlo',zlo,'zup',zup,'kup',kup
!             write(6,*)'new kbotcedsair',kbotcedsair(i,j,1),ktopcedsair(i,j,1)
!             flush(6)
             klo=0
!             do k=1,25
!               write(6,*)'new intaltair',k,intaltair(k),'intcedsair',intcedsair(i,j,k,1)
!
!             end do
             
           endif
           if(kup<0)then
!             if(kup<0.and.iam==75.and.i==1.and.j==49)then
!             write(6,*)'new kup',kup,'i',i,j,ip,'zlo',zlo,'zup',zup,'klo',klo
!             write(6,*)'new kbotcedsair',kbotcedsair(i,j,1),ktopcedsair(i,j,1)
!             flush(6)
             kup=0
           endif
!           if(coiplo/=0.0)then
!         if(tile==5.and.i==168.and.j==192)then
!       if(tile==5.and.ibeg+i-1==168.and.j==192)then
!             write(6,*)'iam',iam,i,j,'zgeoh',ip-1,zgeoh(ip-1),ip,zgeoh(ip)
!             write(6,*)'coiphi-coiplo',coiphi-coiplo
!             write(6,*)'ip',ip,'zup-zlo',zup-zlo
!                     !write(6,*)'sacco',ip,(coiphi-coiplo)/(zup-zlo)*1.e-5
        !           endif
                 l=ltb(ip)
                 do m=1,numcedsair
                   select case (ccedsair(m))
                   case ('CO')
!                     if(iam==75.and.i==1.and.j==49)then
!                       write(6,*)i,j,ip,'new codeaflo',aflo,'up',afup,'klo',klo,'kup',kup,'P-ceds_co',p_ceds_co
!                       flush(6)
!                       do k=kbotcedsair(i,j,1),ktopcedsair(i,j,1)
!                         write(6,*)'newcode intcedsair',k,intcedsair(i,j,k,p_ceds_co),intaltair(k)
!                       end do
!                     endif
                     intvallo=intcedsair(i,j,klo+1,p_ceds_co)*aflo+intcedsair(i,j,klo,p_ceds_co)*(1.-aflo)
                     intvalup=intcedsair(i,j,kup+1,p_ceds_co)*afup+intcedsair(i,j,kup,p_ceds_co)*(1.-afup)
!                    if(iam==75.and.i==1.and.j==49)then
!                      write(6,*)'new bottom klo',klo,'kup',kup
!                      write(6,*)'new bottom aflo',aflo,'afup',afup
!                      write(6,*)'intvallo',intvallo,'up',intvalup
!                      write(6,*)'intcedsair',klo,intcedsair(i,j,klo:klo+1,p_ceds_co)
!                      write(6,*)'intcedsair',kup,intcedsair(i,j,kup:kup+1,p_ceds_co)
!                      write(6,*)i,j,'new kbotcedsair',kbotcedsair(i,j,1),'ktopceddair',ktopcedsair(i,j,1)
!                      write(6,*)i,j,ip,'newintvallo',intvallo
!                      write(6,*)'intvalup',intvalup,'zup',zup,'zlo',zlo,'diff',zup-zlo
!                    flush(6)
!                    endif
                     sacco=(intvalup-intvallo)/(zup-zlo)*1.e-5
                     sacco=max(sacco,0.0)
                     ratiococeds(i,j,l)=sacco/dens*dt/max(cogrd(i,l),1.e-20)
        !             if(iam==75.and.i==1.and.j==49)then
        !             if(iam==75)then
        !             write(6,*)'newsacco',sacco
        !             flush(6)
        !             write(6,*)'lb',lbound(srccocedsair),'ub',ubound(srccocedsair)
        !             flush(6)
        !             endif
        !             srccocedsair(i,j,ip)=sacco
                   case ('NOx')
                     intvallo=intcedsair(i,j,klo+1,p_ceds_nox)*aflo+intcedsair(i,j,klo,p_ceds_nox)*(1.-aflo)
                     intvalup=intcedsair(i,j,kup+1,p_ceds_nox)*afup+intcedsair(i,j,kup,p_ceds_nox)*(1.-afup)
                     sacnox=(intvalup-intvallo)/(zup-zlo)*1.e-5
                     sacnox=max(sacnox,0.0)
                   case ('SO2')
                     intvallo=intcedsair(i,j,klo+1,p_ceds_so2)*aflo+intcedsair(i,j,klo,p_ceds_so2)*(1.-aflo)
                     intvalup=intcedsair(i,j,kup+1,p_ceds_so2)*afup+intcedsair(i,j,kup,p_ceds_so2)*(1.-afup)
                     sacso2=(intvalup-intvallo)/(zup-zlo)*1.e-5
                     sacso2=max(sacso2,0.0)
#ifdef DIAGCEDS
                     if(tile==1)then
                       if(i+ibeg-1==178.and.j==31)then
                         write(6,*)'do so2grd'
                         flush(6)
                         write(6,*)'docedsair so2grd in',l,so2grd(i,l),'sacso2',sacso2,'tend',sacso2/dens*dt
                         flush(6)
                       endif
                     endif
#endif
                     so2grd(i,l)=so2grd(i,l)+sacso2/dens*dt
#ifdef DIAGCEDS
                     if(tile==1)then
                       if(i+ibeg-1==178.and.j==31)then
                         write(6,*)'docedsair so2grd out',l,so2grd(i,l)
                       endif
                     endif
#endif
                     ratioso2ceds(i,j,l)=sacso2/dens*dt/max(so2grd(i,l),1.e-20)
                     so2grd(i,l)=so2grd(i,l)+sacso2/dens*dt
                     srcso2cedsair(i,j,l)=sacso2
                   case ('BC')
                     intvallo=intcedsair(i,j,klo+1,p_ceds_bc)*aflo+intcedsair(i,j,klo,p_ceds_bc)*(1.-aflo)
                     intvalup=intcedsair(i,j,kup+1,p_ceds_bc)*afup+intcedsair(i,j,kup,p_ceds_bc)*(1.-afup)
                     sacbc=(intvalup-intvallo)/(zup-zlo)*1.e-5
                     sacbc=max(sacbc,0.0)*1.e9 ! now ug/kg
#ifdef DIAGCEDS
                     if(tile==1)then
                       if(i+ibeg-1==178.and.j==31)then
                         write(6,*)'do bcgrd'
                         flush(6)
                         write(6,*)'bc1grd in',l,bc1grd(i,l),'sacbc',sacbc,'tend',sacbc/dens*dt
                         flush(6)
                       endif
                     endif
#endif
                     bc1grd(i,l)=bc1grd(i,l)+sacbc/dens*dt*.8
                     bc2grd(i,l)=bc2grd(i,l)+sacbc/dens*dt*.2
#ifdef DIAGCEDS
                     if(tile==1)then
                     if(tile==1)then
                       if(i+ibeg-1==178.and.j==31)then
                         write(6,*)'bc1grd out',l,bc1grd(i,l)
                       endif
                     endif
#endif
                   case ('OC')
                     intvallo=intcedsair(i,j,klo+1,p_ceds_oc)*aflo+intcedsair(i,j,klo,p_ceds_oc)*(1.-aflo)
                     intvalup=intcedsair(i,j,kup+1,p_ceds_oc)*afup+intcedsair(i,j,kup,p_ceds_oc)*(1.-afup)
                     sacoc=(intvalup-intvallo)/(zup-zlo)*1.e-5
                     sacoc=max(sacoc,0.0)*1.e9 ! noew ug/kg
#ifdef DIAGCEDS
                     if(tile==1)then
                       if(i+ibeg-1==178.and.j==31)then
                         write(6,*)'do ocgrd'
                         flush(6)
                         write(6,*)'oc1grd in',l,oc1grd(i,l),'sacoc',sacoc,'tend',sacoc/dens*dt
                         flush(6)
                       endif
                     endif
#endif
                     oc1grd(i,l)=oc1grd(i,l)+sacoc/dens*dt*.5
                     oc2grd(i,l)=oc2grd(i,l)+sacoc/dens*dt*.5
#ifdef DIAGCEDS
                     if(tile==1)then
                     if(tile==1)then
                       if(i+ibeg-1==178.and.j==31)then
                         write(6,*)'oc1grd out',l,oc1grd(i,l)
                       endif
                     endif
#endif
                   caSe ('N2O')
                     intvallo=intcedsair(i,j,klo+1,p_ceds_n2o)*aflo+intcedsair(i,j,klo,p_ceds_n2o)*(1.-aflo)
                     intvalup=intcedsair(i,j,kup+1,p_ceds_n2o)*afup+intcedsair(i,j,kup,p_ceds_n2o)*(1.-afup)
                     sacn2o=(intvalup-intvallo)/(zup-zlo)*1.e-5
                     sacn2o=max(sacn2o,0.0)
                     ration2oceds(i,j,l)=sacn2o/dens*dt/max(xn2ogrd(i,l),1.e-20)
                     xn2ogrd(i,l)=xn2ogrd(i,l)+sacn2o/dens*dt
                     srcn2ocedsair(I,J,l)=sacn2o
#ifdef DIAGCEDSAIRN2O
                     if(ration2oceds(i,j,l)>.02)then
!                     if(tile==5)then
!                       if(i+ibeg-1==18.and.j==114)then
                         write(6,*)'tile=',tile,i,j,'ii',i+ibeg-1,'docedsair n2ogrd out',l,xno2grd(i,l),'tend',sacn2o/dens*dt
                         write(6,*)i,j,'ii',i+ibeg-1,'ratio',l,ration2oceds(i,j,l)
                       endif
                     if(tile==5)then
                       if(i+ibeg-1==18.and.j==114)then
                         write(6,*)'tile=',tile,i,j,'ii',i+ibeg-1,'docedsair n2ogrd out',l,xno2grd(i,l),'tend',sacn2o/dens*dt
                         write(6,*)i,j,'ii',i+ibeg-1,'ratio',l,ration2oceds(i,j,l)
                       endif
                     endif
#endif
!                   caSe ('NH3')
!                     intvallo=intcedsair(i,j,klo+1,p_ceds_nh3)*aflo+intcedsair(i,j,klo,p_ceds_nh3)*(1.-aflo)
!             intvalup=intcedsair(i,j,kup+1,p_ceds_nh3)*afup+intcedsair(i,j,kup,p_ceds_nh3)*(1.-afup)
             !sacnh3=(intvalup-intvallo)/(zup-zlo)*1.e-5
!             sacnh3=max(sacnh3,0.0)
!             nh3grd(i,l)=nh3grd(i,l)+sacnh3/dens*dt
           case default
           end select
         end do ! M  
           

         endif ! zlo< ip/=1
      end subroutine docedsair
