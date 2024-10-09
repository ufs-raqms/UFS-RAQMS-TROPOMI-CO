module raqmschem_species_mod
#include <options.h>
#include <choosechem.h>
  implicit none
    public 
  integer nchemtable
!  parameter (nchemtable=80)
  parameter (nchemtable=81)
  logical :: didsetchempointers=.false.


    ! -- atmospheric tracers
    integer :: p_atm_shum    = 0
    integer :: p_atm_cldq    = 0
    integer :: p_atm_o3mr    = 0
    ! -- chemical tracers
    integer :: p_qc          = 0 
    integer :: p_qi          = 0 
    integer :: p_qv          = 0 
    integer :: p_brcl        = 0
    integer :: p_brno3       = 0
    integer :: p_bry         = 0
    integer :: p_ccl4        = 0
    integer :: p_xmrp        = 0
    integer :: p_prdp        = 0
    integer :: p_ch2h6       = 0
    integer :: p_ch2o        = 0
    integer :: p_ch3br       = 0
    integer :: p_ch3cl       = 0
    integer :: p_ch3ooh      = 0
    integer :: p_ch4         = 0
    integer :: p_cl2         = 0
    integer :: p_clno3       = 0
    integer :: p_cly         = 0
    integer :: p_co          = 0
    integer :: p_f11         = 0
    integer :: p_f1211       = 0
    integer :: p_f12         = 0
    integer :: p_f1301       = 0
    integer :: p_h2o2        = 0
    integer :: p_vrp         = 0
    integer :: p_hbr         = 0
    integer :: p_hcl         = 0
    integer :: p_rip         = 0
    integer :: p_hno3t       = 0
    integer :: p_hno4        = 0
    integer :: p_hobr        = 0
    integer :: p_hocl        = 0
    integer :: p_mtcfm       = 0
    integer :: p_n2o5        = 0
    integer :: p_n2o         = 0
    integer :: p_no2         = 0
    integer :: p_no3         = 0
    integer :: p_noy         = 0
    integer :: p_oclo        = 0
    integer :: p_ox          = 0
    integer :: p_c2h6        = 0
    integer :: p_ald2        = 0
    integer :: p_ethooh      = 0
    integer :: p_pan         = 0
    integer :: p_par         = 0
    integer :: p_onit        = 0
    integer :: p_aone        = 0
    integer :: p_rooh        = 0
    integer :: p_mgly        = 0
    integer :: p_eth         = 0
    integer :: p_olet        = 0
    integer :: p_olei        = 0
    integer :: p_isop        = 0
    integer :: p_isoprd      = 0
    integer :: p_propar      = 0
    integer :: p_ch3oh       = 0
    integer :: p_mvk         = 0
    integer :: p_macr        = 0
    integer :: p_mpan        = 0
    integer :: p_hcn         = 0
    integer :: p_ch3cn       = 0
    integer :: p_so2         = 0
!    integer :: p_so4aer      = 0
    integer :: p_sulf      = 0
    integer :: p_dms         = 0
    integer :: p_msa         = 0
    integer :: p_no3aer      = 0
    integer :: p_nh3         = 0
    integer :: p_nh4aer      = 0
    integer :: p_bc1         = 0
    integer :: p_bc2         = 0
    integer :: p_oc1         = 0
    integer :: p_oc2         = 0
    integer :: p_dust1         = 0
    integer :: p_dust2         = 0
    integer :: p_dust3         = 0
    integer :: p_dust4         = 0
    integer :: p_dust5         = 0
    integer :: p_seas1         = 0
    integer :: p_seas2         = 0
    integer :: p_seas3         = 0
    integer :: p_seas4         = 0
    integer :: p_seas5         = 0
    integer :: p_br2           = 0
#ifdef CO2550
    integer :: p_co50d       = 0
    integer :: p_co25d       = 0
#endif
    integer :: p_coanth25    = 0
    integer :: p_bbcod25    = 0
!    integer,parameter :: nchemfulldim     =58
!    integer,parameter :: nchemfulldim     =57
    integer,parameter :: nchemfulldim     =59
    integer :: nchemfull
  logical ischemfull(nchemtable)
  integer, parameter :: nsolmax=15 ! for fv3
  integer :: idwetd(nsolmax),nsol
  integer idaod(nsolmax),idaodgsi(nsolmax)
  data idaod/15*-1/,idaodgsi/15*-1/
  character *8 chemname(nchemtable),cheminput(nchemtable),chemfull(nchemfulldim)
  character *8 soluablechemname(nchemtable),cheminputlist(30+nchemtable)
  character *8 familyname(4)
  data familyname/'ox','noy','cly','bry'/
  integer chempoint(30+nchemtable),chempointsoluable(30+nchemtable)
  integer chempointcum(30+nchemtable)
  integer chempointsoluablefull(30+nchemtable)
  logical chemsoluable(30+nchemtable),lraqmschem(30+nchemtable)
  logical chemsoluablefull(30+nchemtable)
  data chemname/'brcl','brno3','bry','ccl4','xmrp', &
                'prdp','ch2o','ch3br','ch3cl','ch3ooh', &
                'ch4','cl2','clno3','cly','co', &
                'f11','f1211','f12','f1301','h2o2', &
                'vrp','hbr','hcl','rip','hno3t', &
                'hno4','hobr','hocl','mtcfm','n2o5', &
                'n2o','no2','no3','noy','oclo', & 
                'ox','c2h6','ald2','ethooh','pan', &
                'par','onit','aone','rooh','mgly', &
                'eth','olet','olei','isop','isoprd', &
                'propar','ch3oh','mvk','macr','mpan', &
                'hcn','ch3cn','so2','sulf','DMS','msa', &
                'no3aer','nh3','nh4aer','bc1','bc2', &
                'oc1','oc2','dust1','dust2','dust3', &
                'dust4','dust5','seas1','seas2','seas3','seas4', &
#ifdef CO2550
                 'seas5','co50d', 'co25d','br2'/
#else
                 'seas5','coanth25', 'bbcod25','br2'/
#endif
  integer ichempt(nchemtable),icheminputpt(nchemtable),ninput
  integer ichemfullpt(nchemfulldim)
  data chemfull/'brcl','brno3','bry','ccl4','xmrp', &
                'prdp','ch2o','ch3br','ch3cl','ch3ooh', &
                'ch4','cl2','clno3','cly','co', &
                'f11','f1211','f12','f1301','h2o2', &
                'vrp','hbr','hcl','rip','hno3t', &
                'hno4','hobr','hocl','mtcfm','n2o5', &
                'n2o','no2','no3','noy','oclo', & 
                'ox','c2h6','ald2','ethooh','pan', &
                'par','onit','aone','rooh','mgly', &
                'eth','olet','isop','isoprd', &
                'propar','ch3oh','mvk','macr','mpan', &
                'hcn','ch3cn','br2','coanth25','bbcod25'/
!                'hcn','ch3cn','br2','olei'/

  contains
    subroutine setchempointers(mype)
    use field_manager_mod, only : find_field_index
    implicit none
    integer ipoint,m,msol,nfull
    integer, optional :: mype
    logical :: debug=.false.
    didsetchempointers=.true.
    ischemfull=.false.
    ninput=0
    ischemfull(1:nchemfulldim-1)=.true.
    chempoint=0
    p_atm_shum=find_field_index(1,'sphum')
    p_qv=p_atm_shum
    p_atm_cldq=find_field_index(1,'liq_wat')
!    p_qc=p_atm_cldq
    p_qc=2 ! ajl new
    p_qi=3 ! ajl new
!    write(6,*)'p_atm_shum',p_atm_shum,'p_atm_cldq',p_atm_cldq,'P-qc',p_qc
!    call flush(6)
    p_atm_o3mr=find_field_index(1,'o3mr')
!    write(6,*)'p_atm_shum',p_atm_shum,p_atm_cldq
    lraqmschem=.false.
    nchemfull=nchemfulldim
    do m=1,nchemfulldim
      ichemfullpt(m)=find_field_index(1,chemfull(m))
      if(ichemfullpt(m)>0)then
        lraqmschem(ichemfullpt(m))=.true.
        chempointcum(ichemfullpt(m))=m
        if(present(mype))then
          if(mype.eq.0)then
            write(6,*)m,'ichemfullpt',ichemfullpt(m),' chemfull ',chemfull(m)
            call flush(6)
          endif
        endif
      endif
    end do
    cheminputlist=' '
    do m=1,nchemtable
      ichempt(m)=find_field_index(1,chemname(m))
      if(ichempt(m)>0)then
        chempoint(ichempt(m))=m
!        if(chemname(m).eq.'seas5')ischemfull(m)=.true.
        if(chemname(m).eq.'br2')ischemfull(m)=.true.
        ninput=ninput+1
!        if(debug)then
        if(present(mype))then
          if(mype.eq.0)then
            write(6,*)m,'ichempt',ichempt(m),'chemname',trim(chemname(m))
            call flush(6)
          endif
!        endif
        endif
        cheminput(ninput)=chemname(m)
        cheminputlist(ichempt(m))=chemname(m)
        icheminputpt(ninput)=ichempt(m)
        select case (chemname(m))
          case ('brcl')
            p_brcl=ichempt(m)
          case ('brno3')
            p_brno3=ichempt(m)
          case ('bry')
            p_bry=ichempt(m)
          case ('ccl4')
            p_ccl4=ichempt(m)
          case ('xmrp')
            p_xmrp=ichempt(m)
          case ('prdp')
            p_prdp=ichempt(m)
          case ('ch2o')
            p_ch2o=ichempt(m)
          case ('ch3br')
            p_ch3br=ichempt(m)
          case ('ch3cl')
           p_ch3cl=ichempt(m)
          case ('ch3ooh')
           p_ch3ooh=ichempt(m)
          case ('ch4')
           p_ch4=ichempt(m)
          case ('cl2')
           p_cl2=ichempt(m)
          case ('clno3')
           p_clno3=ichempt(m)
          case ('cly')
            p_cly=ichempt(m)
          case ('co')
            p_co=ichempt(m)
!                'f11','f1211','f12','f1301','h2o2', &
          case ('f11')
            p_f11=ichempt(m)
          case ('f1211')
            p_f1211=ichempt(m)
          case ('f12')
            p_f12=ichempt(m)
          case ('f1301')
            p_f1301=ichempt(m)
          case ('h2o2')
            p_h2o2=ichempt(m)
!                'vrp','hbr','hcl','rip','hno3t', &
          case ('vrp')
            p_vrp=ichempt(m)
          case ('hbr')
            p_hbr=ichempt(m)
          case ('hcl')
            p_hcl=ichempt(m)
          case ('rip')
            p_rip=ichempt(m)
          case ('hno3t')
            p_hno3t=ichempt(m)
!                'hno4','hobr','hocl','mtcfm','n2o5', &
          case ('hno4')
            p_hno4=ichempt(m)
          case ('hobr')
            p_hobr=ichempt(m)
          case ('hocl')
            p_hocl=ichempt(m)
          case ('mtcfm')
            p_mtcfm=ichempt(m)
          case ('n2o5')
            p_n2o5=ichempt(m)
!                'n2o','no2','no3','noy','oclo', & 
          case ('n2o')
            p_n2o=ichempt(m)
          case ('no2')
            p_no2=ichempt(m)
          case ('no3')
            p_no3=ichempt(m)
          case ('noy')
            p_noy=ichempt(m)
          case ('oclo')
            p_oclo=ichempt(m)
!                'ox','c2h5','ald2','ethooh','pan', &
          case ('ox')
            p_ox=ichempt(m)
          case ('c2h6')
            p_c2h6=ichempt(m)
          case ('ald2')
            p_ald2=ichempt(m)
          case ('ethooh')
            p_ethooh=ichempt(m)
          case ('pan')
            p_pan=ichempt(m)
!                'par','onit','aone','rooh','mgly', &
          case ('par')
            p_par=ichempt(m)
          case ('onit')
            p_onit=ichempt(m)
          case ('aone')
            p_aone=ichempt(m)
          case ('rooh')
            p_rooh=ichempt(m)
          case ('mgly')
            p_mgly=ichempt(m)
                !'eth','olet','olei','isop','isoprd', &
          case ('eth')
            p_eth=ichempt(m)
          case ('olet')
            p_olet=ichempt(m)
          case ('olei')
            p_olei=ichempt(m)
          case ('isop')
            p_isop=ichempt(m)
          case ('isoprd')
            p_isoprd=ichempt(m)
!                'propar','ch3oh','mvk','macr','mpan', &
          case ('propar')
            p_propar=ichempt(m)
          case ('ch3oh')
            p_ch3oh=ichempt(m)
          case ('mvk')
            p_mvk=ichempt(m)
          case ('macr')
            p_macr=ichempt(m)
          case ('mpan')
            p_mpan=ichempt(m)
!                'hcn','ch3cn','so2','sulf','DMS', &
          case ('hcn')
            p_hcn=ichempt(m)
          case ('ch3cn')
            p_ch3cn=ichempt(m)
          case ('so2')
            p_so2=ichempt(m)
          case ('sulf')
            p_sulf=ichempt(m)
            idaod(1)=p_sulf
            idaodgsi(1)=p_sulf
          case ('DMS')
            p_dms=ichempt(m)
!                'no3aer','nh3','nh4aer','bc1','bc2', &
          case ('msa')
            p_msa=ichempt(m)
          case ('no3aer')
            p_no3aer=ichempt(m)
          case ('nh3')
            p_nh3=ichempt(m)
          case ('nh4aer')
            p_nh4aer=ichempt(m)
          case ('bc1')
            p_bc1=ichempt(m)
            idaod(2)=p_bc1
            idaodgsi(2)=p_bc1
          case ('bc2')
            p_bc2=ichempt(m)
            idaod(3)=p_bc2
            idaodgsi(3)=p_bc2
          case ('oc1')
            p_oc1=ichempt(m)
            idaod(4)=p_oc1
            idaodgsi(4)=p_oc1
          case ('oc2')
            p_oc2=ichempt(m)
            idaod(5)=p_oc2
            idaodgsi(5)=p_oc2
          case ('dust1')
            p_dust1=ichempt(m)
            idaod(7)=p_dust1
            idaodgsi(6)=p_dust1
!            write(6,*)'dust1',m,p_dust1
          case ('dust2')
            p_dust2=ichempt(m)
            idaod(8)=p_dust2
            idaodgsi(7)=p_dust2
          case ('dust3')
            p_dust3=ichempt(m)
            idaod(9)=p_dust3
            idaodgsi(8)=p_dust3
          case ('dust4')
            p_dust4=ichempt(m)
            idaod(10)=p_dust4
            idaodgsi(9)=p_dust4
!                'dust4', 'seas1','seas2','seas3','seas4', &
          case ('dust5')
            p_dust5=ichempt(m)
!            write(6,*)'dsut5',m,p_dust5
!            idaod(10)=p_dust5
            idaodgsi(10)=p_dust5
          case ('seas1')
            p_seas1=ichempt(m)
            idaod(11)=p_seas1
            idaodgsi(11)=p_seas1
          case ('seas2')
            p_seas2=ichempt(m)
            idaod(12)=p_seas2
            idaodgsi(12)=p_seas2
          case ('seas3')
            p_seas3=ichempt(m)
            idaod(13)=p_seas3
            idaodgsi(13)=p_seas3
          case ('seas4')
            p_seas4=ichempt(m)
            idaod(14)=p_seas4
            idaodgsi(14)=p_seas4
          case ('seas5')
            p_seas5=ichempt(m)
            idaod(15)=p_seas5
            idaodgsi(15)=p_seas5
!                 'seas5','co50d', 'co25d'/
#ifdef CO2550
          case ('co50d')
            p_co50d=ichempt(m)
          case ('co25d')
            p_co25d=ichempt(m)
#endif
          case ('coanth25')
            p_coanth25=ichempt(m)
          case ('bbcod25')
            p_bbcod25=ichempt(m)
          case ('br2')
            p_br2=ichempt(m)
          case default
        if(present(mype))then
            if(mype.eq.0)then
            write(6,*)'missing ',m,chemname(m)
            call flush(6)
            endif
            endif
        end select 
      endif
    enddo
!   if(mype==0)then
!   do m=1,15
!     write(6,*)'idaodgsi',m,idaodgsi(m),chemname(idaodgsi(m))
!   end do
!   endif
!    if(present(mype))then
!    write(6,*)'min ichempt ',minval(ichempt)
!    endif

    if(debug)then
    if(present(mype))then
      if(mype.eq.0)then
!        write(6,*)'p_co',p_co,' p_no2 ',p_no2,' p_no3 ',p_no3,' p_dust1 ',p_dust1
!        call flush(6)
        do m=1,ninput
          write(6,*)'input',m,cheminput(m),' point ',icheminputpt(m)
        enddo
      endif
    endif
    endif
    nsol=0
    if(p_hno3t>0)then
      nsol=nsol+1
      idwetd(nsol)=p_hno3t
    endif
    if(p_h2o2>0)then
      nsol=nsol+1
      idwetd(nsol)=p_h2o2
    endif
    if(p_ch2o>0)then
      nsol=nsol+1
      idwetd(nsol)=p_ch2o
    endif
    if(p_ch3ooh>0)then
      nsol=nsol+1
      idwetd(nsol)=p_ch3ooh
    endif
    if(p_hcl>0)then
      nsol=nsol+1
      idwetd(nsol)=p_hcl
    endif
#ifdef ISOPRENE_PEROX
    if(p_rip>0)then
      nsol=nsol+1
      idwetd(nsol)=p_rip
    endif
    if(p_prdp>0)then
      nsol=nsol+1
      idwetd(nsol)=p_prdp
    endif
    if(p_xmrp>0)then
      nsol=nsol+1
      idwetd(nsol)=p_xmrp
    endif
    if(p_vrp>0)then
      nsol=nsol+1
      idwetd(nsol)=p_vrp
    endif
#else
    if(p_hf>0)then
      nsol=nsol+1
      idwetd(nsol)=p_hf
    endif
#endif
    if(p_hno4>0)then
      nsol=nsol+1
      idwetd(nsol)=p_hno4
    endif
    if(p_hbr>0)then
      nsol=nsol+1
      idwetd(nsol)=p_hbr
    endif
    if(p_ethooh>0)then
      nsol=nsol+1
      idwetd(nsol)=p_ethooh
    endif      
    if(p_rooh>0)then
      nsol=nsol+1
      idwetd(nsol)=p_rooh
    endif
    if(p_mgly>0)then
      nsol=nsol+1
      idwetd(nsol)=p_mgly
    endif
    chemsoluable=.false.
    chemsoluablefull=.false.
    chempointsoluable=0
    chempointsoluablefull=0
    if(mype.eq.0)then
      do msol=1,nsol
        write(6,*)'msol',msol,'idwetd',idwetd(msol),cheminputlist(idwetd(msol))
      end do
    endif
    do msol=1,nsol
      soluablechemname(msol)=chemname(chempoint(idwetd(msol)))
      do nfull=1,nchemfulldim
        if(soluablechemname(msol).eq.chemfull(nfull))then
!          chemsoluable(idwetd(m))=.true.
!          soluablechemname(msol)=chemfull(nfull)
          chemsoluable(nfull)=.true.
          chemsoluablefull(idwetd(msol))=.true.
          if(mype.eq.0)then
            write(6,*)'set chempointsoluable nfull ',nfull,'m',msol,' idwetd ',idwetd(msol)
            call flush(6)
           endif
           chempointsoluable(nfull)=msol
           chempointsoluablefull(idwetd(msol))=msol
          if(mype.eq.0)then
            write(6,*)'chempointsoluable(',nfull,')=',msol
            write(6,*)'chempointsoluablefull(',idwetd(msol),')=',msol
          endif
           exit
         endif
      end do
    end do
    if(mype.eq.0)then
      write(6,*)'nsol',nsol
      do m=1,nsol
        write(6,*)'sol',m,idwetd(m),chemname(chempoint(idwetd(m)))
      end do
      call flush(6)
      do m=1,4+nchemtable
        if(chemsoluable(m))then
          write(6,*)'soluable grell at ',m,' point ',chempointsoluable(m)
          call flush(6)
        endif
        if(chemsoluablefull(m))then
          write(6,*)'soluable ls at ',m,' point ',chempointsoluablefull(m)
        endif
      end do
    endif
    return
    end subroutine setchempointers

end module raqmschem_species_mod
