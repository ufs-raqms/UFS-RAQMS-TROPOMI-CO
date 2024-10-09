#include <options.h>
	subroutine rdstab
      use raqmschem_pmgrid_mod, only : masterproc,iam,prefi,masterproct
      use mpishorthand
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast
      use raqmschem_comm_mod, only : raqmschem_comm_all_bcast_tile
      include 'combined_model_logic'
c
C NAME: rdstab
C PURPOSE: read in radiation source function data table for subsequent
C		interpolation.  Also read table indexing arrays.
C CATEGORY:
C CALLING SEQUENCE: Call rdstab
C INPUTS:
C   304 --  unit number of file containing table index values, ie., pressures,
C		wavelengths, sza, and o3 cols (no3,npr) used to generate stab
C		Values were written at same time as stab.
C   303 --  unit number of file containing s data table.  
C		Data must be on wavelength grid consistent with x-sections.
C		Table is written by wrt_stab.pro
C   nlam -- number of wavelength intervals used (checked against read-in)
C   nsza -- number of zenith angles in table
C   no3s -- number of column o3 values in table
C OPTIONAL INPUT PARAMETERS:
C KEYWORD PARAMETERS:
C OUTPUTS:
C   rlam -- array of wavelength values (nlam) used to generate s table
C   sza_tab -- array of sza values (nsza) used to generate s table
C   o3_tab -- array of col o3 values (no3,npr) used to generate s table
C   sdat -- array of s data (nzens,no3s,nprs,nwavelengths)
C   o2jdat -- array (nzens,no3s,nprs) of J(O2) values from srcfun output
C COMMON BLOCKS:
c   kawa
C SIDE EFFECTS:
C PROCEDURE:
C   Straightforward.
C RESTRICTIONS:
C   
C MODIFICATION HISTORY: 
C   Created 930824 - SR Kawa
C   Modified 950404 to read in J(O2) values from srcfun appended to s table
c   Modified 951128 for LaRC model - Richard Eckman
C   Modified 20000124 to include wavelength dependence in the evaluation
C            of SR band photolysis of O2 -  Murali Natarajan
C            ODF formulation for J(NO) evaluation (Minschwaner and Siskind)
C            included -- Murali Natarajan
C
C   Modified 20000216 --  Read solar flux data from Judith Lean (average data
C            from 1986 and 1989)
C
C   Modified 20020315 --  Read source function for 10% albedo to be used
C                         in the troposphere
C
C   Modified 20020410 --  nx is now 47 due to addition of nmhc species
C            (carbon bond mechanism); also additions for acetone quantum
C            yield -- JRO/jaa

	parameter (nlam=79,nsza=20,no3=12)
        parameter (nx=47,nts=200)
	parameter (npr=24)		!number of model press levels
c
        character*16 xname
      common /kawa/ xtab(nlam,nts,nx),sdat(nsza,no3,npr,nlam),
     1              o2jdat(nsza,no3,npr),rlam(nlam),
     2              o3_tab(no3,npr),pr_tab(npr)
      common /kawa8/sza_tab(nsza)
      real *4 xtab,sdat,o2jdat,rlam,o3_tab,pr_tab
      common /kawa10/sdat10(nsza,no3,npr,nlam)
      real *4 sdat10
c
c       SR  Band data (Opacity Distribution Function)
      common /srbx/ sodf(6,16),wo2(6),sflux(nlam)
      real *4 sodf,wo2,sflux
      integer localrc
c
C Read in dimensions of index arrays and check for consistency
      if(masterproc.or.masterproct)then
        open(304,file=trim(prefi)//'Kawa/indx_larc79.dat',
     1       form='formatted',status='old')
	read (304,*) npr_in,nlam_in,nsza_in,no3_in
c	print *,'s data array dimensions:',nsza_in,no3_in,npr_in,nlam_in
        if((npr_in.NE.npr).OR.(nlam_in.NE.nlam).OR.(nsza_in.NE.nsza).OR.
     1		(no3_in.NE.no3)) then
	   print *,'Array dimensions of table indices do not match',
     1		' those expected',npr_in,npr,nlam_in,nlam,nsza_in,
     2		nsza,no3_in,no3
	   return
	end if

C  Read in indexing values
	read(304,*) pr_tab
	read(304,*) rlam
	read(304,*) sza_tab
	read(304,*) o3_tab
	close(304)

C Read in s function data in binary format.  
        open(303,file=trim(prefi)//'Kawa/stab_larc79.dat',
     1       form='unformatted',status='old')
	read(303) sdat
	read(303) o2jdat
	close(303)
CCCMOD 03/2002
C Read in s function data (10% albedo) in binary format.  
        open(313,file=trim(prefi)//'Kawa/stab10_larc79.dat',
     1       form='unformatted',status='old')
	read(313) sdat10
        close (313)
CCC
C Read solar flux
        open(329,file=trim(prefi)//'Kawa/solarf_79_8689avg.dat',
     1       form='formatted',status='old')
        READ(329,101) xname
101     format (a16)
        READ(329,*)(sflux(i),i=1,nlam)
        CLOSE(329)
c  
c Read  SR Band data  (ODF)
        open(305,file=trim(prefi)//'Kawa/srb_odf.d',
     1      form='unformatted',status='old')
        read (305) sodf
        close (305)
      endif
      call raqmschem_comm_all_bcast_tile(pr_tab,rc=localrc)
      call raqmschem_comm_all_bcast_tile(rlam,rc=localrc)
      call raqmschem_comm_all_bcast_tile(sza_tab,rc=localrc)
      call raqmschem_comm_all_bcast_tile(o3_tab,rc=localrc)
      call raqmschem_comm_all_bcast_tile(sdat,rc=localrc)
      call raqmschem_comm_all_bcast_tile(sdat10,rc=localrc)
      call raqmschem_comm_all_bcast_tile(o2jdat,rc=localrc)
      call raqmschem_comm_all_bcast_tile(sflux,rc=localrc)
      call raqmschem_comm_all_bcast_tile(sodf,rc=localrc)
c
	return
	end
c-----------------------------------------------------------------------

	subroutine jcalc(ik,pfac,sza,colo3,slnto3c,slnto2c,
     1                   slntno,t,dens,rjs,ktf)
c
C+
C NAME: jcalc
C PURPOSE:
C   calculate jrates
C CATEGORY:
C CALLING SEQUENCE:
C   Call jcalc(ik,pfac,sza,colo3,slnto2c,slntno,t,dens,rjs)
C INPUTS:
C   ik -- altitude index (1=bottom, 24=top)
C   sza -- array of zenith angle values
C   colo3 -- overhead O3 value
C   slnto2c -- slant o2 column ( for SR Band evaluation )  
C   t -- CTM temperature value
C   xtab -- table (nlam,nts,nx) of x-section data for range of Ts (150-349 K)
C   sza_tab -- array of sza values (nsza) used to generate s table
C   o3_tab -- array of *vertical* col o3 values (no3,npr) used to generate 
C             s table
C   sdat -- array of s data (nsza,no3s,nprs,nwavelengths)
C OUTPUTS:
C   rjs -- modified array of photolysis rates
C RESTRICTIONS:
C    currently set up for 40-J set (see var nx)
C REQUIRED ROUTINES:
C     interp_s
C MODIFICATION HISTORY: 
C   Created 930826 - SR Kawa
C   931123-SRK-remade xtab to do multiplication by solar flux beforehand and
c   Modified 951128 for LaRC model - Richard Eckman
C     removed inputs:
C     xdat -- cross-section, quantum yield, solar flux data read in once before
C     xnames -- character identifyers for xdat
C	940225 modified to add 3 additional Js, incl N2O
C	950404 modified to handle J(O2) separately
C	950726 modified for LARC model grid
C     20000203 modified to include wavelength dependence in
C              J(O2) evaluation.  J(NO) calculation using ODF (Minschwaner
C              and Siskind JGR).  J values interpolated to proper
C              pressure at each node.
	parameter (npr=24)			!larc model grid
	parameter (nx=47,nlam=79,nsza=20,no3=12)
	parameter (nts=200)

	real colo3,t				!CTM inputs
      common /kawa/ xtab(nlam,nts,nx),sdat(nsza,no3,npr,nlam),
     1              o2jdat(nsza,no3,npr),rlam(nlam),
     2              o3_tab(no3,npr),pr_tab(npr)
      common /kawa8/sza_tab(nsza)
      real *4 xtab,sdat,o2jdat,rlam,o3_tab,pr_tab
c
c
      common /srbx/ sodf(6,16),wo2(6),sflux(nlam)
      real *4 sodf,wo2,sflux
c
	real s(nlam),rjs(nx)
        real pfac,s_l(nlam)               !interp in pressure
        integer ind 
c
        real a220(6),a300(6)
        data a220 /0.00000,5.38571e-20,1.28143e-19,2.02429e-19,
     1	    2.76714e-19,3.51000e-19/
 	data a300 /1.42857e-21,3.21429e-20,6.28571e-20,9.35714e-20,
     1	    1.24286e-19,1.55000e-19/
        save a220,a300
c
CCCCC   WAVELENGTH SUMMATION OVER SPECIFIC WAVELENGTH RANGES FOR EACH SPECIES
c NOTE: Not yet optimized for NMHC species, still summed from 1 to 79
        dimension ivb(nx),ivt(nx)
        data ivb/12,1,1,15,5,6,50,1,8,79,79,1,13,1,1,23,12,36,1,1,
     1       12,1,6,48,48,14,8,17,9,41,15,48,1,31,6,9,9,9,6,6,
     1       1,1,1,1,1,1,1/
        data ivt/31,79,56,73,32,58,67,9,59,79,79,31,64,1,25,43,60,
     1       79,38,32,75,42,22,57,59,64,54,61,78,79,78,63,1,69,67,
     1       38,40,45,28,28,
     1       79,79,79,79,79,79,79/
        save ivb,ivt
c 
        ind = min(max(t-148.5,1.),200.) 
c
c       interpolate radiative flux function values to model conds. 
c
CCCMOD 3/2002
C        call interp_s(ik,sza,colo3,s)
        call interp_s(ik,ktf,sza,colo3,s)
        if (pfac .lt. 1.) then
CCCMOD 3/2002
C          call interp_s(ik-1,sza,colo3,s_l)
          call interp_s(ik-1,ktf,sza,colo3,s_l)
          do ilam=1,nlam
            s(ilam) = pfac*s(ilam) + (1.-pfac)*s_l(ilam)
          enddo
        endif
c
c       calculate J values, multiply cross section.solarflux by s (lam)
        indt = ind 
c
        do ix=1,nx 
         rjs(ix)=0.0
ccc        photolysis rates summation over wavelength.  j(NO) is
ccc        evaluated by a separate call to subroutine 'jnoval'
c
         if (ix .ne. 14) then
cccc          do ilam=1,nlam
           do ilam=ivb(ix),ivt(ix)
            rjs(ix) = rjs(ix) + s(ilam)*xtab(ilam,indt,ix)
           end do
         else
           call jnoval (slnto2c,slnto3c,slntno,dens,rjno)
           rjs(ix)=rjno
         endif
        end do
C
C         Include SR Band photolysis of O2 (first 16 wavelength bands)
c            Index for O2 photolysis = 1
c
        do ilam = 1,16   
         trsx = 0.0
         sgx = 0.0
         do i=1,6
          trsx = trsx + exp(-sodf(i,ilam)*slnto2c)*wo2(i)
          sgx = sgx + sodf(i,ilam)*exp(-sodf(i,ilam)*slnto2c)*wo2(i)
         end do
         if (trsx .gt. 0.0) then
          sgx = sgx/trsx
         else
          sgx = 0.0
         endif
         rjs(1) = rjs(1) + s(ilam)*sgx*sflux(ilam)
        end do
c
c       CH2O quantum yield dependence on m, T, wavel
        do ilam=55,59	   
          ia=ilam-54
       	  atemp=min(a220(ia),max(0.,(t-220.)/80.
     1		*(a300(ia)-a220(ia))+a220(ia)))		!T-dep
          alpha=atemp*7.243e18*pr_tab(ik)/t		!m-dep
          qy=1./(1.+alpha)			!total quantum yield at p & T
          qy300=1./(1.+a300(ia)*2.45e19)	!qy at ats press and room T
          rjs(25) = rjs(25)-s(ilam)*xtab(ilam,indt,25)*(1.-qy/qy300)
          rjs(24) = rjs(24)-s(ilam)*xtab(ilam,indt,24)*(1.-qy/qy300)
        end do
c
c JRO nmhc additions
c     Acetone quantum yield dependence on m, wavel
        do ilam=47,57
         wavl=rlam(ilam)/10.
         aden=pr_tab(ik)/(1.38e-19*t)
         if(wavl.lt.308) then
          xtemp=xtab(ilam,indt,45)
     1     /(-15.696+.05707*wavl+
     2       aden*exp(-88.81+0.15161*wavl))
         else
          xtemp=xtab(ilam,indt,45)
     1     /(-130.2+.42884*wavl+
     2      aden*exp(-55.947+.044913*wavl))
         endif
         rjs(45)=rjs(45)-s(ilam)*xtab(ilam,indt,45)+s(ilam)*xtemp
        end do
c JRO end
c
        return
        end 
c-----------------------------------------------------------------------

	subroutine interp_s(ik,ktf,sza,ovho3,s)

C+
C NAME:
C   interp_s
C PURPOSE:
C   Interpolate s values for each wavelength in table to specified O3
c   col and zenith angles
C CATEGORY:
C CALLING SEQUENCE:
C   Call interp_s(ik,sza,ovho3,s)
C INPUTS:
C   ik --   altitude index (1=bottom, 24=top)
C   nlam -- number of wavelength intervals used
C   sza --  input array zenith angle values for which s
c           is desired
C   ovho3 --input overhead o3 column value 
C   sdat -- input array of values of radiative source function 
C           (nzens,no3,np_ctm,nlam) gridded to ctm p layers
C   sza_tab -- values of sza corresponding to stab table values
C   o3_tab -- array of overhead O3 values at each p-level (no3s,np_ctm)
C                used to index stab
C OPTIONAL INPUT PARAMETERS:
C OUTPUTS:
C   s -- array of s values (nlam) for each wavelength 
C             at model p-level interpolated to ovho3 and sza values
C COMMON BLOCKS:
c   kawa
C SIDE EFFECTS:
C PROCEDURE:
C   bi-linear interpolation, for sza>94 s=0, for o3 out of range use min/max
C RESTRICTIONS:
C REQUIRED ROUTINES:
C	stab and other tabs read in previously with rdstab.f
C MODIFICATION HISTORY: 
C   Created 930825 - SR Kawa
C   950404 modified to handle J(O2) separately
c   Modified 951128 for LaRC model - Richard Eckman
C-
C36
    	parameter (nlam=79,nsza=20,no3=12,npr=24)
       parameter (nx=47,nts=200)
      common /kawa/ xtab(nlam,nts,nx),sdat(nsza,no3,npr,nlam),
     1              o2jdat(nsza,no3,npr),rlam(nlam),
     2              o3_tab(no3,npr),pr_tab(npr)
      common /kawa8/sza_tab(nsza)
      real *4 xtab,sdat,o2jdat,rlam,o3_tab,pr_tab

CCCMOD 3/2002
      common /kawa10/ sdat10(nsza,no3,npr,nlam)
      real *4 sdat10
	real s(nlam),sza,ovho3  
	
        integer locsza
        real t

c       For each input solar zenith angle, find the first element of
c       tabled sza_tab values that is greater than it.  Use this
c       table element and previous table element to determined 
c       interpolated value. 
        do is=1,nsza
          if (sza_tab(is).gt.sza) go to 333 
        end do 
333     locsza = is 
            
        if (locsza.eq.(nsza+1)) then 
c         Point is in darkness.  Set s's to zero. 
          do il=1,nlam
            s(il) = 0. 
          end do 
        else  
          ijj = locsza                       
          t=(sza-sza_tab(ijj-1))/
     1      (sza_tab(ijj)-sza_tab(ijj-1)) 
c         For each input overhead o3 column find the first element
c         of tabled o3_tab values that is > than it.  Use this
c         table element and previous table element to determine
c         interpolated value.             
          do is=1,no3 
            if (o3_tab(is,ik).gt.ovho3) go to 334 
          end do 
 334      ikk = is 
          ikkm = ikk-1 
          if ((ikk.gt.1).and.(ikk.le.no3)) then 
            u=(ovho3-o3_tab(ikkm,ik))/
     1        (o3_tab(ikk,ik)-o3_tab(ikkm,ik))          

c             do bilinear interpolation at ik for each wavelength
c             from numerical recipes, p.96

	    do il=1,nlam       
              if (ktf .le. 0) then
	       s(il)=(1.-t)*(1.-u)*
     1              sdat(ijj-1,ikkm,ik,il)+t*(1.-u)*
     1		    sdat(ijj,ikkm,ik,il)+t*u*
     1              sdat(ijj,ikk,ik,il)+
     1	            (1.-t)*u*sdat(ijj-1,ikk,ik,il)
               else
               s(il)=(1.-t)*(1.-u)*
     1              sdat10(ijj-1,ikkm,ik,il)+t*(1.-u)*
     1              sdat10(ijj,ikkm,ik,il)+t*u*
     1              sdat10(ijj,ikk,ik,il)+
     1              (1.-t)*u*sdat10(ijj-1,ikk,ik,il)
               endif
            end do
                    
          else if (ikk.eq.1) then 
c           print 33,ik,ovho3,o3_tab(1,ik)
c 33        format(' Ovhd o3 col(',i3,') of ',e10.3,
c     1            ' too thin! Lowest tabled val=',e10.3)  
      
            do il=1,nlam
             if (ktf .le. 0) then
	      s(il)=(1.-t)*sdat(ijj-1,1,ik,il)+
     1              t*sdat(ijj,1,ik,il)
              else
              s(il)=(1.-t)*sdat10(ijj-1,1,ik,il)+
     1              t*sdat10(ijj,1,ik,il)
              endif
	    end do        

	  else 
c           print 34,ik,ovho3,o3_tab(no3,ik)
c  34       format(' Ovhd o3 col(',i3,') of ',e10.3,
c     1              ' too thick! Highest tabled val=',e10.3) 
            do il=1,nlam
             if (ktf .le. 0) then
	      s(il)=(1.-t)*sdat(ijj-1,no3,ik,il)+
     1              t*sdat(ijj,no3,ik,il)
              else
              s(il)=(1.-t)*sdat10(ijj-1,no3,ik,il)+
     1              t*sdat10(ijj,no3,ik,il)
              endif 
	    end do 
          end if  
	end if

	return
	end 
c-----------------------------------------------------------------------

        subroutine jnoval (o2c,o3c,sno,dens,rjno)
c
ccccccccc  (Opacity distribution function formulation and parameters from
ccccccccc   Minschwaner and siskind.  
ccccccccc   Note:  Oscillator strength for
ccccccccc   Delta 0-0 band, based on results from
ccccccccc   J.E.Murray et al. J. Chem.Phys, v 101, pp 62-73, July 1994,
ccccccccc   should be increased by about 50%.  This version of the code
ccccccccc   uses the data as given in Minshwaner and Siskind)
c
c        o2c --    o2 slant column
c        o3c --    o3 slant column
c        sno -- NO slant column
c        dens --  total number density
c        rjno --  NO photolysis rate constant (output)
c
      parameter (nlam=79)
      real x1o2(6),x2o2(6),x3o2(6),w11(6),w21(6),w12(6),w22(6),
     1 w13(6),w23(6),x11(6),x21(6),x12(6),x22(6),x13(6),x23(6)
c 
c       SR  Band data
      common /srbx/ sodf(6,16),wo2(6),sflux(nlam)
      real *4 sodf,wo2,sflux
      dimension sgn(3)

c     rayleigh x-section
      data ak1/4.103592e-25/,ak2/4.999e-25/
      save ak1,ak2
c
      data x1o2,x2o2,x3o2/
     1 1.117e-23,2.447e-23,7.188e-23,3.042e-22,1.748e-21,1.112e-20,
     2 1.35e-22,2.991e-22,7.334e-22,3.074e-21,1.689e-20,1.658e-19,
     3 2.968e-22,5.831e-22,2.053e-21,8.192e-21,4.802e-20,2.655e-19/
      save x1o2,x2o2,x3o2
C
      data w11,w21,w12,w22,w13,w23/
     1 0.,0.0512,0.136,0.165,0.141,0.045,
     2 0.,0.00568,0.0152,0.0183,0.0157,0.005,
     3 0.,0.,0.00193,0.0973,0.0975,0.0348,
     4 0.,0.,2.14e-4,1.08e-2,1.08e-2,3.86e-3,
     5 4.5e-2,1.8e-1,2.25e-1,2.25e-1,1.8e-1,4.5e-2,
     6 5.0e-3,0.02,0.025,0.025,0.02,5.0e-3/
      save w11,w21,w12,w22,w13,w23
C
      data x11,x21,x12,x22,x13,x23/
     1 0.,1.32e-18,6.35e-19,7.09e-19,2.18e-19,4.67e-19,
     2 0.,4.41e-17,4.45e-17,4.5e-17,2.94e-17,4.35e-17,
     3 0.,0.,3.05e-21,5.76e-19,2.29e-18,2.21e-18,
     4 0.,0.,3.2e-21,5.71e-17,9.09e-17,6.0e-17,
     5 1.8e-18,1.5e-18,5.01e-19,7.2e-20,6.72e-20,1.49e-21,
     6 1.4e-16,1.52e-16,7.0e-17,2.83e-17,2.73e-17,6.57e-18/
      save x11,x21,x12,x22,x13,x23

      tau1 = 5.055e-19 * o3c + 4.77382*o2c*ak1
      tau2 = 7.46e-19 * o3c + 4.77382*o2c*ak2
C
      fn00 = 1.65e9/(1.65e9+5.1e7+1.17126e-9*dens)
c
      sgn(1)=0.0
      sgn(2)=0.0
      sgn(3)=0.0
c
      do 100 i=1,6
       trm1=o2c*x1o2(i)
       if (trm1 .lt. 36.) sgn(1)=sgn(1)+exp(-trm1)*(w11(i)*x11(i)*
     1  exp(-x11(i)*sno)+w21(i)*x21(i)*exp(-x21(i)*sno))
       trm2 = o2c*x2o2(i)
       if (trm2 .lt. 36.) sgn(2)=sgn(2)+exp(-trm2)*(w12(i)*x12(i)*
     1  exp(-x12(i)*sno)+w22(i)*x22(i)*exp(-x22(i)*sno))
       trm3 = o2c*x3o2(i)
       if (trm3 .lt. 36.) sgn(3)=sgn(3)+exp(-trm3)*(w13(i)*x13(i)*
     1  exp(-x13(i)*sno)+w23(i)*x23(i)*exp(-x23(i)*sno))
100    continue
c
ccccccc  NOTE:  the wavelength intervals used by Minschwaner and Siskind (M&S)
ccccccc         do not exactly match the 3-D model intervals.  Solar flux
ccccccc         in the adjacent intervals are combined with appropriate
ccccccc         weighting factors to get the solar flux corresponding to
ccccccc         M&S's data
c
      rjno = exp(-tau1)*sgn(1)*fn00*(0.1667*sflux(9)+sflux(10)+
     1      0.1053*sflux(11)) + exp(-tau2)*(sgn(2)*(0.235*sflux(5)+
     2      0.647*sflux(6))+sgn(3)*(0.125*sflux(4)+0.765*sflux(5)))
c
      return
      end
c-----------------------------------------------------------------------

      block data photolysisdata
c
      parameter (nlam=79,nsza=20,no3=12)
c
c       SR  Band data (Opacity Distribution Function)
      common /srbx/ sodf(6,16),wo2(6),sflux(nlam)
      real *4 sodf,wo2,sflux
      data wo2/0.05,0.2,0.25,0.25,0.2,0.05/
c
c
      end
