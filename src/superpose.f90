!======================================================================!
!
       module superpose
! 
       implicit none
!
       contains
!
!======================================================================!
!
! DSUPERPOSE - Double precision SUPERPOSE structures
!
! This subroutine superpose two input set of points
!
       subroutine dsuperpose(flag,ndim,nat,fcoord,rcoord)
!
       use geometry
       use utils,    only: print_end
!
       implicit none 
!
! Input/output variables
!
       real(kind=8),dimension(ndim,nat),intent(inout)  ::  fcoord  !  Input coordinates
       real(kind=8),dimension(ndim,nat),intent(in)     ::  rcoord  !  Reference coordinates
       integer,intent(in)                              ::  ndim    !  Space dimension
       integer,intent(in)                              ::  nat     !  Number of points
       character(len=*),intent(in)                     ::  flag    !  Calculation type flag
!
! Local variables
!
       real(kind=8),dimension(ndim,nat)                ::  mataux  !  Auxiliary reference coordinates
       real(kind=8),dimension(ndim,ndim)               ::  rota    !  Rotation matrix
       real(kind=8),dimension(ndim)                    ::  rcen    !  Reference structure centroid
!
! Translating both sets of coordinates so that their centroid coincides
!   with the origin of the coordinate system
! 
       mataux(:,:) = fcoord(:,:)
       call translate(ndim,nat,mataux,-1.0d0,dcenvec(ndim,nat,fcoord), &
                      fcoord)
!
       rcen = dcenvec(ndim,nat,rcoord)
       call translate(ndim,nat,rcoord,-1.0d0,rcen,mataux)
!
! Calculating the optimal rotation matrix
!
       select case ( flag )
         case ('KABSCH','KAB','ROTKAB','ROTKABSCH')
           call drotkabsch(ndim,nat,fcoord,mataux,rota)
         case ('QUAT1','KK','QUATKK')
           call drotkk(nat,fcoord,rcoord,rota)
         case ('QUAT2','CTK','QUATCTK')
           call drotctk(nat,fcoord,rcoord,rota) 
         case default
           write(*,'(2X,68("="))')
           write(*,'(3X,A)') 'ERROR:  Unknown rotation type option'
           write(*,*)
           write(*,'(3X,A)') 'Option '//trim(flag)//' not known'
           write(*,'(3X,A)') 'Plase, choose among KABSCH, QUAT1 or'//  &
                                                                ' QUAT2'
           write(*,'(2X,68("="))')
           call print_end()
       end select
!
! Rotating the input set of points
!
       call DGEMM('N','N',ndim,nat,ndim,1.0d0,rota,ndim,               &
                  fcoord,ndim,0.0d0,mataux,ndim)
!
! Translating to the original position of the reference structure
!
       call translate(ndim,nat,mataux,1.0d0,rcen,fcoord)
!
       return
       end subroutine dsuperpose
!
!======================================================================!
!
! DROTKABSCH - Double precision ROTation KABSCH
!
! This subroutine calculates the optimal rotation matrix that minimizes
!  the RMSD between a reference structure RCOORD(NDIM,NAT) and the tar-
!  get structure FCOORD(NDIM,NAT) by the Kabsch algorithm. Structures
!  must be centred on their baricenters before calling this subroutine
!
! References
! ----------
!
!  - https://en.wikipedia.org/wiki/Kabsch_algorithm
!  - Kabsch, W., "A solution for the best rotation to relate two sets of 
!     vectors", Acta Crystallographica Section A, 1976, 32, 922-923. 
!     https://doi.org/10.1107/S0567739476001873
!  - Kabsch, W., "A discussion of the solution for the best rotation to 
!     relate two sets of vectors", Acta Crystallographica Section A, 
!     1978, 34, 827-828. https://doi.org/10.1107/S0567739478001680
!
       subroutine drotkabsch(ndim,nat,fcoord,rcoord,rota)
!
       implicit none 
!
! Input/output variables
!
       real(kind=8),dimension(ndim,nat),intent(in)      ::  fcoord  !  Input coordinates
       real(kind=8),dimension(ndim,nat),intent(in)      ::  rcoord  !  Reference coordinates
       real(kind=8),dimension(ndim,ndim),intent(out)    ::  rota    !  Rotation matrix
       integer,intent(in)                               ::  ndim    !  Space dimension
       integer,intent(in)                               ::  nat     !  Number of points
!
! Local variables
!
       real(kind=8),dimension(ndim,ndim)                ::  cova    !  Covariance matrix
       real(kind=8),dimension(ndim,ndim)                ::  U       !
       real(kind=8),dimension(ndim,ndim)                ::  Vt      !
       real(kind=8),dimension(ndim)                     ::  S       !
       real(kind=8)                                     ::  detU    !
       real(kind=8)                                     ::  detVt   !
       integer                                          ::  i,j     !  Indexes
!
! Calculating the covariance matrix
!
       call DGEMM('N','T',ndim,ndim,nat,1.0d0,fcoord,ndim,             &
                  rcoord,ndim,0.0d0,cova,ndim)
!
! Calculating the SVD of the covariance matrix
!
       call SVD('T',ndim,ndim,cova,U,S,Vt)
!
! Correcting the rotation matrix to ensure a right-handed
!  coordinate system
!
       call determinant(ndim,U,detU)
       call determinant(ndim,Vt,detVt)
!
       if ( detU*detVt .lt. 0.0d0 ) then
         do i = 1, ndim
           U(i,ndim) = -U(i,ndim)
         end do
       end if
!
! Computing the optimal rotation matrix (cova is an auxliary variable)
! 
       call DGEMM('N','N',ndim,ndim,ndim,1.0d0,U,ndim,                 &
                  Vt,ndim,0.0d0,cova,ndim)
!
       do i = 1, ndim
         do j = 1, ndim
           rota(i,j) = cova(j,i)
         end do
       end do
!
       return
       end subroutine drotkabsch
!
!======================================================================!
!
! DROTKK - Double precision ROTation Kearsley and Kneller
!
! This subroutine calculates the optimal rotation matrix that minimizes
!  the RMSD between a reference structure RCOORD(3,NAT) and the target
!  structure FCOORD(3,NAT) by means of a quaternions-based algorithm.
!  Input structures must be centred on their baricenters before calling
!  this subroutine
!
! References
! ----------
!
! - Kearsley, S. K., "On the orthogonal transformation used for structu-
!    ral comparisons", Acta Crystallographica Section A, 1989, 45, 
!    208-210. https://doi.org/10.1107/S0108767388010128
! - Kneller, G. R., "Superposition of Molecular Structures using Quater-
!    nions", Molecular Simulation, Taylor & Francis, 1991, 7, 113-119. 
!    https://doi.org/10.1080/08927029108022453
!
       subroutine drotkk(nat,fcoord,rcoord,rota)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  fcoord  !  Input coordinates
       real(kind=8),dimension(3,nat),intent(in)  ::  rcoord  !  Reference coordinates
       real(kind=8),dimension(3,3),intent(out)   ::  rota    !  Rotation matrix
       integer,intent(in)                        ::  nat     !  Number of points
!
! Local variables
!
       real(kind=8),dimension(3,nat)             ::  pcoord  !  Combination matrix
       real(kind=8),dimension(3,nat)             ::  mcoord  !  Subtraction matrix
       real(kind=8),dimension(4,4)               ::  diag    !  Matrix to be diagonalized
       real(kind=8),dimension(4,4)               ::  eigvec  !  Eigenvectors
       real(kind=8),dimension(4)                 ::  quat    !  Quaternion
       real(kind=8)                              ::  cova12  !  Covariance matrix element
       real(kind=8)                              ::  cova13  !  Covariance matrix element
       real(kind=8)                              ::  cova23  !  Covariance matrix element
       real(kind=8)                              ::  cova21  !  Covariance matrix element
       real(kind=8)                              ::  cova31  !  Covariance matrix element
       real(kind=8)                              ::  cova32  !  Covariance matrix element
       real(kind=8)                              ::  plus11  !  Combination variance matrix element
       real(kind=8)                              ::  plus12  !  Combination variance matrix element
       real(kind=8)                              ::  plus13  !  Combination variance matrix element
       real(kind=8)                              ::  plus22  !  Combination variance matrix element
       real(kind=8)                              ::  plus23  !  Combination variance matrix element
       real(kind=8)                              ::  plus33  !  Combination variance matrix element
       real(kind=8)                              ::  mins11  !  Subtraction variance matrix element 
       real(kind=8)                              ::  mins12  !  Subtraction variance matrix element
       real(kind=8)                              ::  mins13  !  Subtraction variance matrix element
       real(kind=8)                              ::  mins22  !  Subtraction variance matrix element
       real(kind=8)                              ::  mins23  !  Subtraction variance matrix element
       real(kind=8)                              ::  mins33  !  Subtraction variance matrix element
       integer                                   ::  i,j     !  Indexes
!
! Doing the change of variables
!
       do i = 1, nat
         do j = 1, 3
           mcoord(j,i) = rcoord(j,i) - fcoord(j,i)
           pcoord(j,i) = rcoord(j,i) + fcoord(j,i)
         end do
       end do
!
! Calculating the covariance and variance matrix elements by hand
!
       cova12 = 0.0d0
       cova13 = 0.0d0
       cova23 = 0.0d0
       cova21 = 0.0d0
       cova31 = 0.0d0
       cova32 = 0.0d0
!
       plus11 = 0.0d0
       plus12 = 0.0d0
       plus13 = 0.0d0
       plus22 = 0.0d0
       plus23 = 0.0d0
       plus33 = 0.0d0
!
       mins11 = 0.0d0
       mins12 = 0.0d0
       mins13 = 0.0d0
       mins22 = 0.0d0
       mins23 = 0.0d0
       mins33 = 0.0d0
!
       do i = 1, nat
         cova12 = cova12 + mcoord(1,i)*pcoord(2,i)
         cova13 = cova13 + mcoord(1,i)*pcoord(3,i)
         cova23 = cova23 + mcoord(2,i)*pcoord(3,i)
         cova21 = cova21 + mcoord(2,i)*pcoord(1,i)
         cova31 = cova31 + mcoord(3,i)*pcoord(1,i) 
         cova32 = cova32 + mcoord(3,i)*pcoord(2,i)
!
         plus11 = plus11 + pcoord(1,i)*pcoord(1,i)
         plus12 = plus12 + pcoord(1,i)*pcoord(2,i)
         plus13 = plus13 + pcoord(1,i)*pcoord(3,i)
         plus22 = plus22 + pcoord(2,i)*pcoord(2,i)
         plus23 = plus23 + pcoord(2,i)*pcoord(3,i)
         plus33 = plus33 + pcoord(3,i)*pcoord(3,i)
!
         mins11 = mins11 + mcoord(1,i)*mcoord(1,i)
         mins12 = mins12 + mcoord(1,i)*mcoord(2,i)
         mins13 = mins13 + mcoord(1,i)*mcoord(3,i)
         mins22 = mins22 + mcoord(2,i)*mcoord(2,i)
         mins23 = mins23 + mcoord(2,i)*mcoord(3,i)
         mins33 = mins33 + mcoord(3,i)*mcoord(3,i)
       end do
!
! Assemblying the upper triangular part of the symmetric matrix
!  to be diagonalized
!
       diag(1,1) = mins11 + mins22 + mins33
!
       diag(1,2) = cova32 - cova23
       diag(2,2) = mins11 + plus22 + plus33
!
       diag(1,3) = cova13 - cova31
       diag(2,3) = mins12 - plus12
       diag(3,3) = plus11 + mins22 + plus33
!
       diag(1,4) = cova21 - cova12
       diag(2,4) = mins13 - plus13
       diag(3,4) = mins23 - plus23
       diag(4,4) = plus11 + plus22 + mins33
!
! Diagonalizing the symmetric matrix (quat is a dummy variable)
!
       call diagonalize(4,diag,quat,eigvec)
!
! Extracting the quaternion associated to the minimum eigenvalue
!
       quat(:) = eigvec(:,1)
!
! Building the optimal rotation matrix
!
       call quat2mat(quat,rota)
!
       return
       end subroutine drotkk
!
!======================================================================!
!
! DROTCTK - Double precision ROTation Coutsias, Theobald and Kneller
!
! This subroutine calculates the optimal rotation matrix that minimizes
!  the RMSD between a reference structure RCOORD(3,NAT) and the target
!  structure FCOORD(3,NAT) by means of a quaternions-based algorithm.
!  Input structures must be centred on their baricenters before calling
!  this subroutine
!
! References
! ----------
!
! - Coutsias, E. A.; Seok, C. & Dill, K. A. "Using quaternions to cal-
!    culate RMSD", Journal of Computational Chemistry, 2004, 25, 
!    1849-1857. https://doi.org/10.1002/jcc.20110
! - Theobald, D. L. "Rapid calculation of RMSDs using a quaternion-based
!    characteristic polynomial", Acta Crystallographica Section A, 2005, 
!    61, 478-480. https://doi.org/10.1107/S0108767305015266
! - Kneller, G. R. "Comment on "Using quaternions to calculate RMSD" 
!    [J. Comp. Chem. 25, 1849 (2004)]", Journal of Computational Chemis-
!    try, 2005, 26, 1660-1662. https://doi.org/10.1002/jcc.20296
!
       subroutine drotctk(nat,fcoord,rcoord,rota)
!
       implicit none 
!
! Input/output variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  fcoord  !  Input coordinates
       real(kind=8),dimension(3,nat),intent(in)  ::  rcoord  !  Reference coordinates
       real(kind=8),dimension(3,3),intent(out)   ::  rota    !  Rotation matrix
       integer,intent(in)                        ::  nat     !  Number of points
!
! Local variables
!
       real(kind=8),dimension(3,3)               ::  cova    !  Covariance matrix
       real(kind=8),dimension(4,4)               ::  diag    !  Matrix to be diagonalized
       real(kind=8),dimension(4,4)               ::  eigvec  !  Eigenvectors
       real(kind=8),dimension(4)                 ::  quat    !  Quaternion
!
! Calculating the covariance matrix
!
       call DGEMM('N','T',3,3,nat,1.0d0,rcoord,3,fcoord,3,0.0d0,cova,3)
!
! Assemblying the upper triangular part of the symmetric matrix
!  to be diagonalized
!
       diag(1,1) =  cova(1,1) + cova(2,2) + cova(3,3)
!
       diag(1,2) =  cova(2,3) - cova(3,2)
       diag(2,2) =  cova(1,1) - cova(2,2) - cova(3,3)
!
       diag(1,3) =  cova(3,1) - cova(1,3)
       diag(2,3) =  cova(1,2) + cova(2,1)
       diag(3,3) = -cova(1,1) + cova(2,2) - cova(3,3)
!
       diag(1,4) =  cova(1,2) - cova(2,1)
       diag(2,4) =  cova(1,3) + cova(3,1)
       diag(3,4) =  cova(2,3) + cova(3,2)
       diag(4,4) = -cova(1,1) - cova(2,2) + cova(3,3)
!
! Diagonalizing the symmetric matrix (quat contains the eigenvalues)
!
       call diagonalize(4,diag,quat,eigvec)
!
! Extracting the quaternion associated to the maximum eigenvalue
!
       quat(:) = eigvec(:,4)
!
! Building the optimal rotation matrix
!
       call quat2mat(quat,rota)
!
       return
       end subroutine drotctk
!
!======================================================================!
!
! QUAT2MAT - QUATernion to MATrix
!
! This subroutine calculates the rotation matrix ROTA(3,3) associated to
!  the unit quaternion QUAT(4) used as rotation operator
!
       subroutine quat2mat(quat,rota)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(4),intent(in)      ::  quat    !  Quaternion
       real(kind=8),dimension(3,3),intent(out)   ::  rota    !  Rotation matrix
!
! Computing by hand the product of the rotation quaternion by itself
!
       rota(1,1) = quat(1)**2 + quat(2)**2 - quat(3)**2 - quat(4)**2
       rota(2,1) = 2.0d0*(quat(2)*quat(3) - quat(1)*quat(4))
       rota(3,1) = 2.0d0*(quat(2)*quat(4) + quat(1)*quat(3))
!
       rota(1,2) = 2.0d0*(quat(3)*quat(2) + quat(1)*quat(4))
       rota(2,2) = quat(1)**2 - quat(2)**2 + quat(3)**2 - quat(4)**2
       rota(3,2) = 2.0d0*(quat(3)*quat(4) - quat(1)*quat(2))
!
       rota(1,3) = 2.0d0*(quat(4)*quat(2) - quat(1)*quat(3))
       rota(2,3) = 2.0d0*(quat(4)*quat(3) + quat(1)*quat(2))
       rota(3,3) = quat(1)**2 - quat(2)**2 - quat(3)**2 + quat(4)**2
!
       return
       end subroutine quat2mat
!
!======================================================================!
!
! DMCALCRMSD - Double precision Matrix CALCulate Root-Mean-Square Deviation
!
! This subroutine computes the root-mean-square deviation between an in-
!  put set of coordinates FCOORD(NDIM,NAT) with respect to a reference
!  set of coordinates RCOORD(NDIM,NAT)
!
       subroutine dmcalcrmsd(ndim,nat,fcoord,rcoord,rmsd)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(ndim,nat),intent(in)  ::  fcoord  !  Input coordinates
       real(kind=8),dimension(ndim,nat),intent(in)  ::  rcoord  !  Reference coordinates
       real(kind=8),intent(out)                     ::  RMSD    !  Root-mean-square deviation
       integer,intent(in)                           ::  ndim    !  Space dimension
       integer,intent(in)                           ::  nat     !  Number of points
!
! Local variables
!
       integer                                      ::  i,j     !  Indexes
!
! Computing RMSD between input structures
! 
       RMSD = 0.0d0 
       do i = 1, nat
         do j = 1, ndim
           RMSD = RMSD + (rcoord(j,i) - fcoord(j,i))**2
         end do
       end do
!
       RMSD = dsqrt(RMSD/ndim)
!
       return
       end subroutine dmcalcrmsd
!
!======================================================================!
!
! DMCALCMAPE - Double precision Matrix CALCulate Mean Absolute Percentage Error
!
! This subroutine computes the mean absolute percentage error between an
!  input set of coordinates FCOORD(NDIM,NAT) with respect to a referen-
!  ce set of coordinates RCOORD(NDIM,NAT)
!
       subroutine dmcalcmape(ndim,nat,fcoord,rcoord,MAPE,BAPE)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(ndim,nat),intent(in)  ::  fcoord  !  Input coordinates
       real(kind=8),dimension(ndim,nat),intent(in)  ::  rcoord  !  Reference coordinates
       real(kind=8),intent(out)                     ::  MAPE    !  Mean absolute percentage error
       real(kind=8),intent(out)                     ::  BAPE    !  Biggest absolute percentage error
       integer,intent(in)                           ::  ndim    !  Space dimension
       integer,intent(in)                           ::  nat     !  Number of points
!
! Local variables
!
       real(kind=8)                                 ::  daux    !
       integer                                      ::  i,j     !  Indexes
!
! Computing RMSD between input structures
! 
       MAPE = 0.0d0 
       BAPE = 0.0d0
!
       do i = 1, nat
         do j = 1, ndim
           daux = dabs((rcoord(j,i) - fcoord(j,i))/rcoord(j,i))
!
           MAPE = MAPE + daux
           BAPE = max(BAPE,daux*100)
         end do
       end do
!
       MAPE = MAPE/ndim*100
!
       return
       end subroutine dmcalcmape
!
!======================================================================!
!
       end module superpose
!
!======================================================================!
