!======================================================================!
!
       module inertia
!
       implicit none
!
       private
       public   :: inertia_moments
!
       contains
!
!======================================================================!
!
       subroutine inertia_moments(nat,coord,mass,axes,moment)
!
! This function returns the principal axes AXES(3,3) and principal mo-
!  ments of inertia MOMENTS(3) of a rigid body of NAT atoms defined by
!  their coordinates COORD(3,NAT) and masses MASS(NAT)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,nat),intent(inout)  ::  coord   !  Coordinates
       real(kind=8),dimension(nat),intent(in)       ::  mass    !  Masses
       real(kind=8),dimension(3,3),intent(out)      ::  axes    !  Principal axes
       real(kind=8),dimension(3),intent(out)        ::  moment  !  Principal moments of inertia
       integer,intent(in)                           ::  nat     !  Number of atoms
!
! Declaration of the local variables
!
       real(kind=8),dimension(3,3)                  ::  tensor  !  Inertia tensor
       real(kind=8),dimension(3)                    ::  com     !  Center of mass vector
       real(kind=8)                                 ::  cum     !  Element of the inertia tensor
       integer                                      ::  i,j     !  Indexes
!
! Translating the molecule to the center of mass
!
       com = dcomvec(nat,coord,mass)
!
       do i = 1, nat
         coord(:,i) = coord(:,i) - com(:)
       end do
!
! Calculating the principal moments of inertia
!
       do j = 1, 3
         do i = 1, j
!
           call inertia_ij(i,j,cum,nat,coord,mass)
!
           tensor(i,j) = cum
!
         end do
       end do
!
       call diagonalize(3,tensor,moment,axes)
! 
       return
       end subroutine inertia_moments
!
!======================================================================!
!
       subroutine inertia_ij(i,j,cumij,nat,coord,mass)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  coord  !  Coordinates
       real(kind=8),dimension(nat),intent(in)    ::  mass   !  Masses
       real(kind=8),intent(out)                  ::  cumij  !  
       integer,intent(in)                        ::  nat    !  Number of atoms
       integer,intent(in)                        ::  i,j    !  Indexes
!
! Declaration of the local variables
!
       real(kind=8)                              ::  r2     !
       integer                                   ::  iat    !
!
! Declaration of BLAS functions
!
       real(kind=8)                              ::  DDOT   !
!
! Calculating element ij of the inertia tensor
!
       cumij = 0.0d0
       do iat = 1, nat
         cumij = cumij - mass(iat)*coord(i,iat)*coord(j,iat)
!
         if ( i .eq. j ) then
           r2    = DDOT(3,coord(:,iat),1,coord(:,iat),1)
           cumij = cumij + mass(iat)*r2
         end if               
       end do
! 
       return
       end subroutine inertia_ij
!
!======================================================================!
!
! DCOMVEC - Double precision Center Of Mass VECtor
!
! This function returns the Center of Mass of a molecule
!
       function dcomvec(nat,coord,mass) result(com)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  coord  !  Coordinates
       real(kind=8),dimension(nat),intent(in)    ::  mass   !  Masses
       real(kind=8), dimension(3)                ::  com    !  Center of Mass vector
       integer,intent(in)                        ::  nat    !  Number of atoms
!
! Declaration of the local variables
!
       real(kind=8)                              ::  totm   ! Total mass
       integer                                   ::  i,j    ! Indexes
!
! Calculating the center of mass
!
       totm = 0.0d0
       do i = 1, nat
         totm = totm + mass(i)
       end do
!
       com(:) = 0.0d0
       do i = 1, nat
         do j = 1, 3
           com(j) = com(j) + mass(i)*coord(j,i)
         end do
       end do
!
       com(:) = com(:) / totm
!
       return
       end function dcomvec
!
!======================================================================!
!
       end module inertia
!
!======================================================================!
