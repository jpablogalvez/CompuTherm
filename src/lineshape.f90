!======================================================================!
!
       module lineshape
!
       use parameters
!
       implicit none
!
       contains
!
!======================================================================!
!
       function gauss(dw,sig) result(broad)
!
       implicit none
!
! Input/Output variables
!
       real(kind=8),intent(in)  ::  dw     !
       real(kind=8),intent(in)  ::  sig    !
       real(kind=8)             ::  broad  !
!
       broad = exp(-dw**2/2/sig**2)/sig/sqrt(2.0*pi)
!
       end function gauss
!
!======================================================================!
!
       function lor(dw,hwhm) result(broad)
!
       implicit none
!
! Input/Output variables
!
       real(kind=8),intent(in)  ::  dw     !
       real(kind=8),intent(in)  ::  hwhm   !
       real(kind=8)             ::  broad  !
!
       broad = 1.0/( pi*hwhm*(1.0 + (dw/hwhm)**2) )
!
       return
       end function lor
!
!======================================================================!
!
       end module lineshape
!
!======================================================================!
