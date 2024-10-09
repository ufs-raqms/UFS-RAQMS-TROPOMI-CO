real*4 zero(384,384)
character *256 file
file='/ships19/aqda/lenzen/CEDS_2019/GEFS-Aerosol_emissions/fengsha/soil/C384/tile1/clay.dat'
zero(:,10)=0.0
!open(10,file='zero.bin',convert='big_endian',form='unformatted')
open(10,file=file,convert='little_endian',form='unformatted')
read(10)zero
ngood=0
do j=1,384
  do i=1,384
    if(.not.isnan(zero(i,j)))then
      write(6,*)i,j,zero(i,j)
      ngood=ngood+1
    endif
   end do
end do
write(6,*)'ngood',ngood
!write(6,*)zero(1,1)
!call flush(6)

close(10)
stop
end
