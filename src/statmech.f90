!======================================================================!
!
       module statmech
!
       use parameters
!
       implicit none
!
       private
       public ::  qtrans,qrot,qvib,qel,                                &
                  Htrans,Hrot,Hvib,                                    &
                  Strans,Srot,Svib,Sel,                                &
                  qivib,Hivib,Sivib,                                   &
                  Gfree,ZPVE,Etrans
!
       contains
!
!======================================================================!
!
       real(kind=8) function qel(T,g0,De)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  T   !  Temperature
       real(kind=8),intent(in)  ::  g0  !  Ground state degeneracy
       real(kind=8),intent(in)  ::  De  !  Potential energy
!
! Computing the electronic partition function
!
       qel = g0*dexp(De*Kb/T)
!
       return
       end function qel
!
!======================================================================!
!
       real(kind=8) function qtrans(T,V,m)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  m  !  Molecular weight
       real(kind=8),intent(in)  ::  T  !  Absolute temperature
       real(kind=8),intent(in)  ::  V  !  Molar volume       
!
! Computing the translational partition function
!
       qtrans = dsqrt(2.0d0*pi*T*(m*uma*Kb/hplanck**2))**3*V/Na
!
       return
       end function qtrans
!
!======================================================================!
!
       real(kind=8) function qrot(T,I)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3),intent(in)  ::  I     !  Moments of inertia       
       real(kind=8),intent(in)               ::  T     !  Absolute temperature
!
! Declaration of the local variables
!
       real(kind=8),dimension(3)             ::  Trot  !  Rotational temperature
!
! Computing the rotational partition function (C1 symmetry)
!
       Trot(:) = ( hplanck**2/(Kb*I(:)*uma*a0**2) )/( 8.0d0*pi**2 )
!
       qrot   = dsqrt(pi)*dsqrt(T/Trot(1))*dsqrt(T/Trot(2))            &
                                                       *dsqrt(T/Trot(3))
!
       return
       end function qrot
!
!======================================================================!
!
       real(kind=8) function qvib(T,dof,freq)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq  !  Moments of inertia       
       real(kind=8),intent(in)                 ::  T     !  Absolute temperature
       integer,intent(in)                      ::  dof   !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8),dimension(Dof)             ::  Tvib  !  Rotational temperature
       integer                                 ::  i     !  Index
!
! Computing the vibrational partition function
!
       Tvib(:)  = hplanck*clight/Kb*freq(:)*100.0d0
!
       qvib = 1.0d0
       do i = 1, dof  
!~          q_vib = q_vib*dexp(-Tvib(i)/(2.0d0*T))/( 1 - dexp(-Tvib(i)/T) )
         qvib = qvib/(1 - dexp(-Tvib(i)/T))
       end do
!
       return
       end function qvib
!
!======================================================================!
!
       real(kind=8) function qivib(T,freq)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  freq  !  Moments of inertia       
       real(kind=8),intent(in)  ::  T     !  Absolute temperature
!
! Declaration of the local variables
!
       real(kind=8)             ::  Tvib  !  Rotational temperature
!
! Computing the vibrational partition function
!
       Tvib  = hplanck*clight/Kb*freq*100.0d0
!
!~        q_vib = q_vib*dexp(-Tvib/(2.0d0*T))/( 1 - dexp(-Tvib/T) )
       qivib = 1.0d0/(1 - dexp(-Tvib/T))
!
       return
       end function qivib
!
!======================================================================!
!
       real(kind=8) function ZPVE(DoF,freq)
!
       implicit none  
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq  !  Moments of inertia       
       integer,intent(in)                      ::  DoF   !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       integer                                 ::  i     !  Index    
!
! Computing the zero point vibrational energy
!
       ZPVE = 0.0d0
!
       do i = 1, dof
         ZPVE = ZPVE + freq(i)
       end do
!
       ZPVE = Na*hplanck*clight*ZPVE*100.0d0/2
!
       return
       end function ZPVE
!
!======================================================================!
!
       real(kind=8) function Etrans(T)
!
       implicit none      
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  T       !  Absolute temperature
!
! Computing the translational contribution to the internal energy
!
       Etrans = 1.50d0*Rjul*T
!
       return
       end function Etrans
!
!======================================================================!
!
       real(kind=8) function Htrans(T)
!
       implicit none     
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  T       !  Absolute temperature  
!
! Computing the translational contribution to the enthalpy
!
       Htrans = 2.50d0*Rjul*T
!
       return
       end function Htrans
!
!======================================================================!
!
       real(kind=8) function Hrot(T,rotdof)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  T       !  Absolute temperature
       integer,intent(in)       ::  rotdof  !  Rotational degrees of freedom 
!
! Computing the rotational contribution to the enthalpy
!
       Hrot = real(rotdof)/2*Rjul*T
!
       return
       end function Hrot
!
!======================================================================!
!
       real(kind=8) function Hvib(T,dof,freq)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq  !  Harmonic frequencies       
       real(kind=8),intent(in)                 ::  T     !  Absolute temperature
       integer,intent(in)                      ::  dof   !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8),dimension(Dof)             ::  Tvib  !  Rotational temperature
       integer                                 ::  i     !  Index
!
! Computing the vibrational contribution to the enthalpy
!
       Tvib(:)  = hplanck*clight/Kb*freq(:)*100.0d0
!
       Hvib = 0.0d0
!
       do i = 1, dof
         Hvib = Hvib + Tvib(i)/(dexp(Tvib(i)/T) - 1.0d0)
       end do
!
       Hvib = Rjul*Hvib
!
       return
       end function Hvib
!
!======================================================================!
!
       real(kind=8) function Hivib(T,freq)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  freq  !  Harmonic frequency    
       real(kind=8),intent(in)  ::  T     !  Absolute temperature
!
! Declaration of the local variables
!
       real(kind=8)             ::  Tvib  !  Rotational temperature
!
! Computing the vibrational contribution to the enthalpy
!
       Tvib   = hplanck*clight/Kb*freq*100.0d0
!
       Hivib = Rjul*Tvib/(dexp(Tvib/T) - 1.0d0)
!
       return
       end function Hivib
!
!======================================================================!
!
       real(kind=8) function Sel(qel)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  qel  !  Electronic partition functions
!
! Computing the translational contribution to the entropy
!
       Sel = Rjul*dlog(qel)
!
       return
       end function Sel
!
!======================================================================!
!
       real(kind=8) function Strans(qtrans)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  qtrans  !  Translational partition functions
!
! Computing the translational contribution to the entropy
!
       Strans = Rjul*(2.5d0 + dlog(qtrans))
!
       return
       end function Strans
!
!======================================================================!
!
       real(kind=8) function Srot(T,qrot,Erot)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  T        !  Temperature
       real(kind=8),intent(in)  ::  qrot     !  Rotational partition function
       real(kind=8),intent(in)  ::  Erot     !  Rotational contribution to the energy
!
! Computing the rotational contribution to the entropy (C1 symmetry)
!
       Srot = Rjul*dlog(qrot) + Erot/T
!
       return
       end function Srot
!
!======================================================================!
!
       real(kind=8) function Svib(T,dof,freq)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq  !  Harmonic frequency 
       real(kind=8),intent(in)                 ::  T     !  Absolute temperature
       integer,intent(in)                      ::  dof   !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8)                            ::  aux   !  Auxiliary double precision number
       integer                                 ::  i     !  Index
!
! Computing the vibrational contribution to the entropy
!
       Svib = 0.0d0
!
       do i = 1, dof
         aux   = hplanck*clight/Kb*freq(i)*100.0d0/T
         Svib = Svib + (aux/(dexp(aux) - 1.0d0) - dlog(1 - dexp(-aux))) 
       end do
!
       Svib = Rjul*Svib
!
       return
       end function Svib
!
!======================================================================!
!
       real(kind=8) function Sivib(T,freq)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  freq  !  Harmonic frequencies
       real(kind=8),intent(in)  ::  T     !  Absolute temperature
!
! Declaration of the local variables
!
       real(kind=8)             ::  aux   !  Auxiliary double precision number
!
! Computing the vibrational contribution to the entropy
!
       aux   = hplanck*clight/Kb*freq*100.0d0/T
!
       Sivib = Rjul*(aux/(dexp(aux) - 1.0d0) - dlog(1 - dexp(-aux))) 
!
       return
       end function Sivib
!
!======================================================================!
!
       real(kind=8) function Gfree(T,q)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  T       !  Temperature
       real(kind=8),intent(in)  ::  q       !  Partition function
!
! Computing the translational contribution to the Gibbs free energy
!
       Gfree = -Rjul*T*dlog(q)
!
       return
       end function Gfree
!
!======================================================================!
!
       end module statmech
!
!======================================================================!
