!======================================================================!
!
       module mathtools
       implicit none
!
       contains
!
!======================================================================!
!
       function factorial(n) result(fact)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)  ::  n
       integer             ::  fact
!
! Local variables
!
       integer             ::  i
!
! Computing the factorial of N
!
       fact = 1
!
       do i = 1, n
         fact = fact*i
       end do
!
       end function factorial
!
!======================================================================!
!
       end module mathtools
!
!======================================================================!
