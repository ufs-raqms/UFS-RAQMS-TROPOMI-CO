module calcdxdy
use chem_types_mod
real(CHEM_KIND_R8), public, parameter :: RADIUS = 6.3712e+6 ! for GFS physics
real(CHEM_KIND_R8), public, parameter :: PI     = 3.1415926535897931 ! for GFS physics
real(CHEM_KIND_R8), public, parameter :: DEGTORAD=pi/180.
contains
 subroutine initdxdy(ids,ide,jds,jde,lonin,latin,ihs,ihe,jhs,jhe,griddx,griddy,de,mype)
 use raqmschem_comm_mod, only : chem_reducetile_pushwithhalo
 implicit none
 integer,intent(in) :: ids,ide,jds,jde,mype
 integer :: is,ie,js,je
 integer,intent(in) :: ihs,ihe,jhs,jhe
 real(CHEM_KIND_R8) :: q1(2),q2(2)
 real(CHEM_KIND_R8),intent(in) :: lonin(ids:ide,jds:jde),latin(ids:ide,jds:jde)
 real(CHEM_KIND_R8) :: londeg(ihs:ihe,jhs:jhe),latdeg(ihs:ihe,jhs:jhe)
 real(CHEM_KIND_R4),allocatable,dimension(:,:),intent(out) :: griddx,griddy
 real(CHEM_KIND_R8) :: lon(ihs:ihe,jhs:jhe),lat(ihs:ihe,jhs:jhe)
 real(CHEM_KIND_R8) :: lonwithhalo(ihs:ihe,jhs:jhe),latwithhalo(ihs:ihe,jhs:jhe)
 integer localrc,de,j,i

 call chem_reducetile_pushwithhalo(lonin,ids,ide,jds,jde,lonwithhalo,ihs,ihe,jhs,jhe,de=de,rc=localrc)
 call chem_reducetile_pushwithhalo(latin,ids,ide,jds,jde,latwithhalo,ihs,ihe,jhs,jhe,de=de,rc=localrc)

 lon=lonwithhalo*degtorad
 lat=latwithhalo*degtorad
   
 allocate (griddx(ids:ide,jds:jde),griddy(ids:ide,jds:jde))
 griddx=0.0
 griddy=0.0
 is=min(ihs,ids)
 js=min(jhs,jds)
 ie=max(ihe,ide)
 je=max(jhe,jde)
 
 do j=jds,jde
   do i=is+1,ie-1
     q1(1)=lon(i-1,j)
     q1(2)=lat(i-1,j)
     q2(1)=lon(i+1,j)
     q2(2)=lat(i+1,j)
!    will fix later to calculate exact centers
     griddx(i,j)=great_circle_distr(q1,q2)*.5
   end do
   if(is.eq.ids)then
     q1(1)=lon(is,j)
     q1(2)=lat(is,j)
     q2(1)=lon(is+1,j)
     q2(2)=lat(is+1,j)
     griddx(ids,j)=great_circle_distr(q1,q2)
   endif
   if(ie.eq.ide)then
     q1(1)=lon(ie-1,j)
     q1(2)=lat(ie-1,j)
     q2(1)=lon(ie,j)
     q2(2)=lat(ie,j)
     griddx(ide,j)=great_circle_distr(q1,q2)
   endif
 end do
 do j=js+1,je-1
   do i=ids,ide
     q1(1)=lon(i,j-1)
     q1(2)=lat(i,j-1)
     q2(1)=lon(i,j+1)
     q2(2)=lat(i,j+1)
     griddy(i,j)=great_circle_distr(q1,q2)*.5
   end do
 end do
 do i=ids,ide
   if(js.eq.jds)then
     q1(1)=lon(i,js)
     q1(2)=lat(i,js)
     q2(1)=lon(i,js+1)
     q2(2)=lat(i,js+1)
     griddy(i,jds)=great_circle_distr(q1,q2)
   endif
   if(je.eq.jde)then
     q1(1)=lon(i,je-1)
     q1(2)=lat(i,je-1)
     q2(1)=lon(i,je)
     q2(2)=lat(i,je)
     griddy(i,jde)=great_circle_distr(q1,q2)
   endif
 end do
 end subroutine initdxdy
 real function great_circle_distr( q1, q2 )
! ajl 1 is lon 2 is lat
! ajl q1 is end q2 is beginning
   use chem_types_mod, only : CHEM_KIND_R8,CHEM_KIND_R4
      real(CHEM_KIND_R8), intent(IN)           :: q1(2), q2(2)
 
      real (CHEM_KIND_R8):: p1(2), p2(2)
      real (CHEM_KIND_R8):: beta 
      integer n

      do n=1,2
         p1(n) = q1(n)
         p2(n) = q2(n)
      enddo

      beta = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
                         sin((p1(1)-p2(1))/2.)**2 ) ) * 2. 

      great_circle_distr = radius * beta 

  end function great_circle_distr
end module calcdxdy

