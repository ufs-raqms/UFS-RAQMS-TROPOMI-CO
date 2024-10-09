      subroutine ztrop(is,ie,js,je,nl,pri,prl,zm,tm,ztrop,ptrop,ktrop)
      integer i,j,k,ktopbot,is,ie,js,je,nl
      real*4,dimension(is:ie,js:je,nl) :: pri,prl,zm,tm
      real*4,dimension(is:ie,js:je) :: psfc,ztrop,ptrop
      integer,dimension(is:ie,js:je) ktrop
      real :: dz,gam,gamctrl
      integer :: iflag
      gamctrl=-2.e-3
      do j=1,nrl
        do i=1,ncl
          ptrop(i,j)=-9999.
          ztrop(i,j)=-9999.
          do k=1,nl
            if(prl(i,j,k)<pri(i,j,1)-300.)then
              kbot=k
              exit
            end do
          end do
          do k=nl,1,-1
            if(prl(i,j,k)>50.)then
              ktop=k
              exit
            endif
          end do
          iflag=0
          do k=kbot,ktop
            dz=zm(i,j,k+1)-zm(i,j,k)
            gam=(tm(i,j,k+1)-tm(i,j,k))/dz
            if(gam>=gamctrl)then
              if(iflag==0)then
                ktropbot=k
                dzlay=dz
                iflag=1
                ptrop(i,j)=prl(i,j,k)
                ztrop(i,j)=am(i,j,k)
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
                  ptrop(i,j)=-9999.
                  ztrop(i,j)=-9999.
                  iflag=0
                  dzlay=0.0
                endif
              endif
            endif
          end do
          if(ptrop(i,j)<0.0)then
            write(6,*)'trop missing at',i,j
          endif
        end do
      end do
      return
      end subroutine ztrop 
