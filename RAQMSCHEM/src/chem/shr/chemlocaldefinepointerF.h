      real*4,pointer,dimension(:,:) :: oxgrd,xno2grd,ch4grd,hclgrd, 
     &  xnoygrd,hno3tgrd,clygrd,xn2o5grd,h2o2grd,clno3grd, 
     &  xn2ogrd,f11grd,f12grd,ccl4grd,ch3clgrd,xmtcfmgrd, 
     &  brygrd,ch3brgrd,f1301grd,f1211grd,hno4grd, 
     &  hoclgrd,cogrd,oclogrd,xno3grd, 
     &  ch2ogrd,ch3oohgrd,hbrgrd,brno3grd,hobrgrd,brclgrd, 
     &  cl2grd,xnoysed, 
#ifdef ISOPRENE_PEROX
     &  vrpgrd,ripgrd,xmrpgrd,prdpgrd, 
#else
     &  h2ogrd,hfgrd,cf2ogrd,cfclogrd, 
#endif
     &  c2h6grd,ald2grd,ethoohgrd,pangrd,pargrd,xonitgrd, 
     &  aonegrd,roohgrd,xmglygrd,ethgrd,xoletgrd,xoleigrd, 
     &  xisopgrd,xisoprdgrd,prop_pargrd,ch3ohgrd,xmvkgrd, 
     &  xmacrgrd,xmpangrd
      real*4,pointer,dimension(:,:) :: cod50grd,cod25grd
      real*4,pointer,dimension(:,:) :: coanth25grd,bbcod25grd
