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
                  qqho,Hqho,Sqho,                                      &
                  qqrrho,Hqrrho,Sqrrho,                                &
                  Gfree,ZPVE,Etrans,ZPVEdamp
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
       real(kind=8),dimension(DoF),intent(in)  ::  freq  !  Vibrational frequencies
       real(kind=8),intent(in)                 ::  T     !  Absolute temperature
       integer,intent(in)                      ::  dof   !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       integer                                 ::  i     !  Index
!
! Computing the vibrational partition function
!
       qvib = 1.0d0
!
       do i = 1, dof  
         qvib = qvib/(1 - dexp(-hplanck*clight/Kb*freq(i)*100.0d0/T))
       end do
!
       return
       end function qvib
!
!======================================================================!
!
! References
! ----------
!
! - Zhao, Y. & Truhlar, D. G., "Computational characterization and mo-
!    deling of buckyball tweezers: density functional study of concave-
!    convex π⋯π interactions", Phys. Chem. Chem. Phys., The Royal So-
!    ciety of Chemistry, 2008, 10, 2813-2818. 
!    http://dx.doi.org/10.1039/B717744E
! - Ribeiro, R. F.; Marenich, A. V.; Cramer, C. J. & Truhlar, D. G., 
!    "Use of Solution-Phase Vibrational Frequencies in Continuum Models 
!    for the Free Energy of Solvation", The Journal of Physical Chemis-
!    try B, 2011, 115, 14556-14562. https://doi.org/10.1021/jp205508z
!
       real(kind=8) function qqho(T,dof,freq,cutoff)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq    !  Moments of inertia       
       real(kind=8),intent(in)                 ::  T       !  Absolute temperature
       real(kind=8),intent(in)                 ::  cutoff  !  Cutoff frequency
       integer,intent(in)                      ::  dof     !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8)                            ::  Tvib    !  Vibrational temperature
       integer                                 ::  i       !  Index
!
! Computing the vibrational partition function (QHO model)
!
       qqho = 1.0d0
!
       do i = 1, dof  
!
         if ( freq(i) .ge. cutoff ) then
           Tvib  = hplanck*clight/Kb*freq(i)*100.0d0
         else
!~            Tvib  = hplanck*clight/Kb*(freq(i)+cutoff)*100.0d0 
           Tvib  = hplanck*clight/Kb*cutoff*100.0d0    ! FLAG: check if it is a sum or not
         end if
!
         qqho = qqho/(1 - dexp(-Tvib/T))
!
       end do
!
       return
       end function qqho
!
!======================================================================!
!
! References
! ----------
!
! - Grimme, S., "Supramolecular Binding Thermodynamics by Dispersion-
!    Corrected Density Functional Theory", Chemistry A European Journal,
!    2012, 18, 9955-9964. https://doi.org/10.1002/chem.201200497
! - Li, Y.-P.; Gomes, J.; Mallikarjun Sharada, S.; Bell, A. T. and 
!    Head-Gordon, M., "Improved Force-Field Parameters for QM/MM Simula-
!    tions of the Energies of Adsorption for Molecules in Zeolites and
!    a Free Rotor Correction to the Rigid Rotor Harmonic Oscillator Mo-
!    del for Adsorption Enthalpies", The Journal of Physical Chemistry 
!    C, 2015, 119, 1840-1850. https://doi.org/10.1021/jp509921r
!
       real(kind=8) function qqrrho(T,dof,freq,cutoff,I)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq    !  Vibrational frequencies
       real(kind=8),dimension(3),intent(in)    ::  I       !  Moments of inertia       
       real(kind=8),intent(in)                 ::  T       !  Absolute temperature
       real(kind=8),intent(in)                 ::  cutoff  !  Cutoff frequency
       integer,intent(in)                      ::  dof     !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8)                            ::  Tvib    !  Vibrational temperature
       real(kind=8)                            ::  Bav     !  Average molecular moment of inertia
       real(kind=8)                            ::  mu      !  Effective moment of inertia
       real(kind=8)                            ::  w       !  Head-Gordon damping function 
       integer                                 ::  idof    !  Index
!
! Computing the vibrational partition function (QRRHO model)
!
       Bav = momentav(I)
!
       qqrrho = 1.0d0
!
       do idof = 1, dof  
         mu = momenteff(freq(idof),Bav)
!
         Tvib = hplanck*clight/Kb*freq(idof)*100.0d0
!
         w      = hgdf(cutoff,freq(idof),4)
!
         qqrrho = qqrrho/(1 - dexp(-Tvib/T))**w                        &
                           *dsqrt(8.0d0*pi**3*mu*Kb*T/hplanck**2)**(1-w)
!
       end do
!
       return
       end function qqrrho
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
       real(kind=8)             ::  Tvib  !  Vibrational temperature
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
       real(kind=8) function ZPVEdamp(DoF,freq,cutoff)
!
       implicit none  
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq  !  Moments of inertia      
       real(kind=8),intent(in)                 ::  cutoff  !  Cutoff frequency 
       integer,intent(in)                      ::  DoF   !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       integer                                 ::  i     !  Index    
!
! Computing the zero point vibrational energy
!
       ZPVEdamp = 0.0d0
!
       do i = 1, dof
         ZPVEdamp = ZPVEdamp + hgdf(cutoff,freq(i),4)*freq(i)
       end do
!
       ZPVEdamp = Na*hplanck*clight*ZPVEdamp*100.0d0/2
!
       return
       end function ZPVEdamp
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
       real(kind=8)                            ::  Tvib  !  Vibrational temperature
       integer                                 ::  i     !  Index
!
! Computing the vibrational contribution to the enthalpy
!
!
       Hvib = 0.0d0
!
       do i = 1, dof
         Tvib = hplanck*clight/Kb*freq(i)*100.0d0
         Hvib = Hvib + Tvib/(dexp(Tvib/T) - 1.0d0)
       end do
!
       Hvib = Rjul*Hvib
!
       return
       end function Hvib
!
!======================================================================!
!
! References
! ----------
!
! - Zhao, Y. & Truhlar, D. G., "Computational characterization and mo-
!    deling of buckyball tweezers: density functional study of concave-
!    convex π⋯π interactions", Phys. Chem. Chem. Phys., The Royal So-
!    ciety of Chemistry, 2008, 10, 2813-2818. 
!    http://dx.doi.org/10.1039/B717744E
! - Ribeiro, R. F.; Marenich, A. V.; Cramer, C. J. & Truhlar, D. G., 
!    "Use of Solution-Phase Vibrational Frequencies in Continuum Models 
!    for the Free Energy of Solvation", The Journal of Physical Chemis-
!    try B, 2011, 115, 14556-14562. https://doi.org/10.1021/jp205508z
!
       real(kind=8) function Hqho(T,dof,freq,cutoff)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq    !  Harmonic frequencies       
       real(kind=8),intent(in)                 ::  T       !  Absolute temperature
       real(kind=8),intent(in)                 ::  cutoff  !  Cutoff frequency
       integer,intent(in)                      ::  dof     !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8)                            ::  Tvib    !  Vibrational temperature of mode i
       integer                                 ::  i       !  Index
!
! Computing the vibrational contribution to the enthalpy (QHO model)
!
       Hqho = 0.0d0
!
       do i = 1, dof
!
         if ( freq(i) .ge. cutoff ) then
           Tvib = hplanck*clight/Kb*freq(i)*100.0d0
         else
!~            Tvib = hplanck*clight/Kb*(freq(i)+cutoff)*100.0d0
           Tvib = hplanck*clight/Kb*cutoff*100.0d0    ! FLAG: check if it is a sum or not
         end if
!
         Hqho = Hqho + Tvib/(dexp(Tvib/T) - 1.0d0)
!
       end do
!
       Hqho = Rjul*Hqho
!
       return
       end function Hqho
!
!======================================================================!
!
! References
! ----------
!
! - Li, Y.-P.; Gomes, J.; Mallikarjun Sharada, S.; Bell, A. T. and 
!    Head-Gordon, M., "Improved Force-Field Parameters for QM/MM Simula-
!    tions of the Energies of Adsorption for Molecules in Zeolites and
!    a Free Rotor Correction to the Rigid Rotor Harmonic Oscillator Mo-
!    del for Adsorption Enthalpies", The Journal of Physical Chemistry 
!    C, 2015, 119, 1840-1850. https://doi.org/10.1021/jp509921r
!
       real(kind=8) function Hqrrho(T,dof,freq,cutoff)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq    !  Harmonic frequencies       
       real(kind=8),intent(in)                 ::  T       !  Absolute temperature
       real(kind=8),intent(in)                 ::  cutoff  !  Cutoff frequency
       integer,intent(in)                      ::  dof     !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8)                            ::  Tvib    !  Vibrational temperature
       real(kind=8)                            ::  w       !  Head-Gordon damping function 
       integer                                 ::  i       !  Index
!
! Computing the vibrational contribution to the enthalpy (QRRHO model)
!
       Hqrrho = 0.0d0
!
       do i = 1, dof
         Tvib   = hplanck*clight/Kb*freq(i)*100.0d0
         w      = hgdf(cutoff,freq(i),4)
         Hqrrho = Hqrrho + w*Tvib/(dexp(Tvib/T) - 1.0d0)               &
                                                       + (1 - w)*0.5d0*T
       end do
!
       Hqrrho = Rjul*Hqrrho
!
       return
       end function Hqrrho
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
       real(kind=8)             ::  Tvib  !  Vibrational temperature
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
       real(kind=8) function Srot(qrot,rotdof)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  qrot    !  Rotational partition function
       integer                  ::  rotdof  !  Rotational degrees of freedom
!
! Computing the rotational contribution to the entropy (C1 symmetry)
!
       Srot = Rjul*(dlog(qrot) + real(rotdof)/2)
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
         aux  = hplanck*clight/Kb*freq(i)*100.0d0/T
         Svib = Svib + aux/(dexp(aux) - 1.0d0) - dlog(1 - dexp(-aux)) 
       end do
!
       Svib = Rjul*Svib
!
       return
       end function Svib
!
!======================================================================!
!
! References
! ----------
!
! - Zhao, Y. & Truhlar, D. G., "Computational characterization and mo-
!    deling of buckyball tweezers: density functional study of concave-
!    convex π⋯π interactions", Phys. Chem. Chem. Phys., The Royal So-
!    ciety of Chemistry, 2008, 10, 2813-2818. 
!    http://dx.doi.org/10.1039/B717744E
! - Ribeiro, R. F.; Marenich, A. V.; Cramer, C. J. & Truhlar, D. G., 
!    "Use of Solution-Phase Vibrational Frequencies in Continuum Models 
!    for the Free Energy of Solvation", The Journal of Physical Chemis-
!    try B, 2011, 115, 14556-14562. https://doi.org/10.1021/jp205508z
!
       real(kind=8) function Sqho(T,dof,freq,cutoff)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq    !  Harmonic frequency 
       real(kind=8),intent(in)                 ::  T       !  Absolute temperature
       real(kind=8),intent(in)                 ::  cutoff  !  Cutoff frequency
       integer,intent(in)                      ::  dof     !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8)                            ::  aux     !  Auxiliary double precision number
       integer                                 ::  i       !  Index
!
! Computing the vibrational contribution to the entropy (QHO model)
!
       Sqho = 0.0d0
!
       do i = 1, dof
!
         if ( freq(i) .ge. cutoff ) then
           aux   = hplanck*clight/Kb*freq(i)*100.0d0/T
         else
!~            aux   = hplanck*clight/Kb*(freq(i)+cutoff)*100.0d0/T 
           aux   = hplanck*clight/Kb*cutoff*100.0d0/T   ! FLAG: check if it is a sum or not
         end if
!
         Sqho = Sqho + (aux/(dexp(aux) - 1.0d0) - dlog(1 - dexp(-aux))) 
!
       end do
!
       Sqho = Rjul*Sqho
!
       return
       end function Sqho
!
!======================================================================!
!
! References
! ----------
!
! - Grimme, S., "Supramolecular Binding Thermodynamics by Dispersion-
!    Corrected Density Functional Theory", Chemistry A European Journal,
!    2012, 18, 9955-9964. https://doi.org/10.1002/chem.201200497
!
       real(kind=8) function Sqrrho(T,dof,freq,cutoff,I)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(DoF),intent(in)  ::  freq    !  Harmonic frequencies   
       real(kind=8),dimension(3),intent(in)    ::  I       !  Moments of inertia           
       real(kind=8),intent(in)                 ::  T       !  Absolute temperature
       real(kind=8),intent(in)                 ::  cutoff  !  Cutoff frequency
       integer,intent(in)                      ::  dof     !  Vibrational degrees of freedom
!
! Declaration of the local variables
!
       real(kind=8)                            ::  Bav     !  Average molecular moment of inertia
       real(kind=8)                            ::  mu      !  Effective moment of inertia
       real(kind=8)                            ::  w       !  Head-Gordon damping function 
       real(kind=8)                            ::  aux     !  Auxiliary double precision number
       integer                                 ::  idof    !  Index
!
! Computing the vibrational contribution to the entropy (QRRHO model)
!
       Bav = momentav(I)
!
       Sqrrho = 0.0d0
!
       do idof = 1, dof
         mu = momenteff(freq(idof),Bav)
!
         aux  = hplanck*clight/Kb*freq(idof)*100.0d0/T
!
         w    = hgdf(cutoff,freq(idof),4)
!
         Sqrrho = Sqrrho + w*(aux/(dexp(aux) - 1.0d0)                  &
                               - dlog(1 - dexp(-aux)))                 &
                  + (1-w)*0.5d0*(1+dlog(8.0d0*pi**3*mu*Kb*T/hplanck**2)) 
       end do
!
       Sqrrho = Rjul*Sqrrho
!
       return
       end function Sqrrho
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
       real(kind=8) function momentav(I)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3),intent(in)  ::  I     !  Principal moments of inertia           
!
! Local variables
!
       integer                               ::  idof  !  Index
!
! Computing the average molecular moment of inertia
!
       momentav = 0.0d0
!
       do idof = 1, 3       ! FLAG: Check if average is over rotational dof
         momentav = momentav + I(idof)*uma*a0**2
       end do
!
       momentav = momentav/3
!
       return
       end function momentav
!
!======================================================================!
!
       real(kind=8) function momenteff(freq,Bav)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),intent(in)  ::  freq  !  Harmonic frequency  
       real(kind=8),intent(in)  ::  Bav   !  Average molecular moment of inertia
!
!  Computing the moment of inertia for a free-rotor with the same
!   frequency as the target normal mode
!
         momenteff = hplanck/(8.0d0*pi**2*freq*clight*100.0d0)
!
!  Calculating the effective moment of inertia
!
         momenteff = momenteff*Bav/(momenteff + Bav)
!
       return
       end function momenteff
!
!======================================================================!
!
! HGDF - Head-Gordon Damping Function 
!
! References
! ----------
!
! - Chai, J.-D. & Head-Gordon, M., "Long-range corrected hybrid density 
!    functionals with damped atom-atom dispersion corrections", Phys. 
!    Chem. Chem. Phys., The Royal Society of Chemistry, 2008, 10, 
!    6615-6620. http://dx.doi.org/10.1039/B810189B
!
       real(kind=8) function hgdf(nu0,nu,alpha)
!
       implicit none
!
! Declaration of the in/out variables
!
       integer,intent(in)       ::  alpha  !  Fitting parameter
       real(kind=8),intent(in)  ::  nu0    !  Threshold frequency
       real(kind=8),intent(in)  ::  nu     !  Target frequency
!
! Computing the Head-Gordon damping function
!
       hgdf = 1.0d0/(1.0d0 + (nu0/nu)**alpha)
!
       return
       end function hgdf
!
!======================================================================!
!
       end module statmech
!
!======================================================================!
