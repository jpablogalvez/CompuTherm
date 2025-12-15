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
       function eunits(x,ener,hwhm,fbroad) result(y)
!
       implicit none
!
! Input/Output variables
!
       real(kind=8),intent(in)  ::  x     !
       real(kind=8),intent(in)  ::  ener  !
       real(kind=8),intent(in)  ::  hwhm  !
       real(kind=8)             ::  y     !
!
! Input subroutines
!
       real(kind=8),external    ::  fbroad  !
!
       y = x*fbroad(x-ener,hwhm)
!
       end function eunits
!
!======================================================================!
!
       function nmunits(x,ener,hwhm,fbroad) result(y)
!
       implicit none
!
! Input/Output variables
!
       real(kind=8),intent(in)  ::  x     !
       real(kind=8),intent(in)  ::  ener  !
       real(kind=8),intent(in)  ::  hwhm  !
       real(kind=8)             ::  y     !
!
! Input subroutines
!
       real(kind=8),external    ::  fbroad  !
!
       y = (1.0E7/x)*fbroad(1.0E7*(1.0d0/x-1.0d0/ener),hwhm)
!
       end function nmunits
!
!======================================================================!
!
       end module lineshape
!
!======================================================================!
