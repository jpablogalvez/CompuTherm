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
! Adapted from: http://rosettacode.org/wiki/Permutations
!
      function nextp ( n, a )
      logical :: nextp
      integer,intent(in) :: n
      integer,dimension(n),intent(inout) :: a
!
!     local variables:
!
      integer i,j,k,t
!
      i = n-1
   10 if ( a(i) .lt. a(i+1) ) goto 20
      i = i-1
      if ( i .eq. 0 ) goto 20
      goto 10
   20 j = i+1
      k = n
   30 t = a(j)
      a(j) = a(k)
      a(k) = t
      j = j+1
      k = k-1
      if ( j .lt. k ) goto 30
      j = i
      if (j .ne. 0 ) goto 40
!      
      nextp = .false.
!      
      return
!
   40 j = j+1
      if ( a(j) .lt. a(i) ) goto 40
      t = a(i)
      a(i) = a(j)
      a(j) = t
!      
      nextp = .true.
!      
      return
      end function
!
!======================================================================!
!
       end module mathtools
!
!======================================================================!
