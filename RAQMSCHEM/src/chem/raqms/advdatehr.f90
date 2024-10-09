  integer yy,mon,dd,hr,min,sec,ndate
  real dsec,hrout
  yy=2019
  mon=7
  dd=1
  hr=12
  min=30
  dsec=900.
  call advdatehr(yy,mon,dd,hr,min,sec,dsec,ndate,hrout)
  print *,'ndate',ndate,' hrout ',hrout
  sec=0
  min=45
  call advdatehr(yy,mon,dd,hr,min,sec,dsec,ndate,hrout)
  print *,'ndate',ndate,' hrout ',hrout
  hr=23
  call advdatehr(yy,mon,dd,hr,min,sec,dsec,ndate,hrout)
  print *,'ndate',ndate,' hrout ',hrout
  stop
  end
  call advdatehr(yy,mm,dd,h,m,s,dtsec,ndate,hrout)
  subroutine advdatehr(yyin,monin,ddin,hrin,minin,secin,dsec,ndate,hrout)
  implicit none
  integer yyin,monin,ddin,hrin,minin,secin
  real dsec,hrout
  integer min,sec,minadd
  integer ndate,yy,mm,dd,hr,mdh,dh,dhr
  integer daysmon(12,2),il,ileap,hradd
  data daysmon/31,28,31,30,31,30,31,31,30,31,30,31, &
               31,29,31,30,31,30,31,31,30,31,30,31/ 
  yy=yyin
  il=mod(yy,4)
  ileap=1
  if(il.eq.0)ileap=2
  mm=monin
  dd=ddin
  hr=hrin
  min=minin
  sec=secin
  sec=sec+dsec
  if(sec>=60.0)then
    minadd=sec/60
    min=min+minadd
    sec=sec-minadd*60
  endif
  if(min>=60.)then
    hradd=min/60.
    hr=hr+hradd
    min=min-hradd*60
  endif
  if(hr>23)then
    hr=hr-24
    dd=dd+1
    if(dd>daysmon(mm,ileap))then
      dd=dd-daysmon(mm,ileap)
      mm=mm+1
      if(mm>12)then
        mm=1
        yy=yy+1
      endif
    endif
  elseif(hr<0)then
    hr=hr+24
    dd=dd-1
    if(dd<1)then
      mm=mm-1
      if(mm<1)then
        mm=12
        yy=yy-1
      endif
     ileap=1
     if(il.eq.0)ileap=2
      dd=daysmon(mm,ileap)
    endif
  endif
  ndate=(((yy*100+mm)*100+dd)*100)+hr
  hrout=float(hr)+float(min)/60.+float(sec)/3600.
  return
  end 
