    character *40 input,output*40
    open(10,file='chemlist',form='formatted')
    open(20,file='outcode',form='formatted')
    open(30,file='outcodefrom',form='formatted')
    do i=1,1000
      read(10,'(a)',end=10)input
      write(6,*)'input',trim(input)
      write(20,*)' ',trim(input(3:))//'grd(i,k)=t3d_in(i,j,k,'//trim(input)//')'
      write(30,*)' ','t3d_out(i,j,k,'//trim(input)//')'//'='//trim(input(3:))//'grd(i,k)'
    enddo 
10    stop
    end


