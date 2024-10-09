pointer :: testp
interface
   subroutine testp
   end subroutine testp
   subroutine testp1
   end subroutine testp1
end interface
testp=>testp1
stop
end
