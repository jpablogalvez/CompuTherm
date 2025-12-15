!======================================================================!
!
! Fundamental physical constants (SI units)
!
! Taken from: http://physics.nist.gov/constants
!
       module parameters
!
       implicit none
!
       save
!
! Declaration of numbers
!
       real(kind=8),parameter  ::  zero     =  1.0e-16
       real(kind=8),parameter  ::  one      =  1.0d0
       real(kind=8),parameter  ::  two      =  2.0d0
       real(kind=8),parameter  ::  three    =  3.0d0
!
! Declaration of the physical constants
!
       real(kind=8),parameter  ::  pi       =  4*atan(1.0_8)         !  Pi number
       real(kind=8),parameter  ::  hplanck  =  6.626070150E-34       !  Planck constant (J*s)
       real(kind=8),parameter  ::  hbar     =  hplanck/(2.0d0*pi)    !  Reduced Planck constant (J*s)
       real(kind=8),parameter  ::  clight   =  299792458.0d0         !  Speed of light in vacuum (m/s)
       real(kind=8),parameter  ::  qe       =  1.6021766340E-19      !  Elementary charge (C)
       real(kind=8),parameter  ::  me       =  9.1093837015280E-31   !  Electron mass (kg)
       real(kind=8),parameter  ::  mp       =  1.67262192369510E-27  !  Proton mass (kg)
       real(kind=8),parameter  ::  Na       =  6.022140760E+23       !  Avogadro constant (mol**-1)
       real(kind=8),parameter  ::  Kb       =  1.3806490E-23         !  Boltzmann constant (J/K)
       real(kind=8),parameter  ::  Rjul     =  8.3144626180          !  Universal gas constant (J/mol·K)
       real(kind=8),parameter  ::  eps0     =  8.854187812813E-12    !  Vacuum electric permittivity (s**2·C**2/kg·m**2)
       real(kind=8),parameter  ::  Ratm     =  Rjul/101.325          !  Universal gas constant (L*atm/mol·K)
       real(kind=8),parameter  ::  uma      =  1.6605390666050E-27   !  Atomic mass unit (kg)
       real(kind=8),parameter  ::  a0       =  5.2917721090380E-11   !  Bohr radius (m)
!
!  Magnetic constant or permeability of vacuum [N/A**2]       1.25663706143592E-06
!  Electric constant or permittivity of vacuum [F/m]          8.85418781762039E-12
!  Electron g factor [ ]                                     -2.00231930436220E+00
!  Fine-structure constant                                    7.29735253760000E-03
!  Rydberg constant [1/m]                                     1.09737315685270E+07
!
! Declaration of the conversion factors
!
!  [eV] -> [cm**-1]                            8065.73
       real(kind=8)            ::  ev2cm    =  qe/clight/hplanck/100.0
!  [cm**-1] -> [eV]                            1.23981E-4
       real(kind=8)            ::  cm2ev    =  clight*hplanck/qe*100.0
!  [eV] -> [kJ/mol]                            96.4869
       real(kind=8)            ::  ev2jmol  =  qe*Na/1000.0
!  [u] -> [a.u.]                               1.82288848426455E+03
       real(kind=8),parameter  ::  uma2au   =  1.82288848426455E+03
!  [Angstrom] -> [Bohr] = [a.u.]               1.88972613288564E+00
       real(kind=8),parameter  ::  ang2au   =  1.88972613288564
!  [a.u.] = [Bohr] -> [Angstrom]               5.29177208590000E-01
       real(kind=8),parameter  ::  au2ang   =  5.2917720859E-01
!  [a.u.] -> [s]                               2.41888432650478E-17
       real(kind=8),parameter  ::  au2s     =  2.41888432650478E-17
!  [a.u.] -> [fs]                              2.41888432650478E-02
       real(kind=8),parameter  ::  au2fs    =  2.41888432650478E-02
       real(kind=8),parameter  ::  fs2au    =  41.34137314d0
!  [a.u.] -> [N]                               8.23872205491840E-08
       real(kind=8),parameter  ::  au2N     =  8.23872205491840E-08
!  [a.u.] -> [J]                               4.35974393937059E-18
       real(kind=8),parameter  ::  au2J     =  4.35974393937059E-18
!  [a.u.] -> [kJ/mol]                          2.62549961709828E+03
       real(kind=8),parameter  ::  au2kJ    =  4.35974393937059E-21*Na
!  [a.u.] -> [kcal/mol]                        6.27509468713739E+02
       real(kind=8),parameter  ::  au2kcal  =  6.27509468713739E+02
       real(kind=8),parameter  ::  kcal2kJ  =  4.184d0
!  [a.u.] -> [eV]                              2.72113838565563E+01
       real(kind=8),parameter  ::  au2ev    =  2.72113838565563E+01
!  [a.u.] -> [K]                               3.15774647902944E+05
       real(kind=8),parameter  ::  au2K     =  3.15774647902944E+05
!  [a.u.] -> [Pa]                              2.94210107994716E+13
       real(kind=8),parameter  ::  au2Pa    =  2.94210107994716E+13
!  [a.u.] -> [bar]                             2.94210107994716E+08
       real(kind=8),parameter  ::  au2bar   =  2.94210107994716E+08
!  [a.u.] -> [atm]                             2.90362800883016E+08
       real(kind=8),parameter  ::  au2atm   =  2.90362800883016E+08
       real(kind=8),parameter  ::  cm12kJ   =  hplanck*clight*Na/10.0d0
       real(kind=8),parameter  ::  cm12au   =  4.55633E-6
!  [a.u.] -> [1/cm] (wave numbers)             2.19474631370540E+05
       real(kind=8),parameter  ::  au2cm1   =  219474.63d0
!
!  [a.u.] -> [Hz]                              6.57968392072181E+15
!  [a.u./Bohr**2] -> [1/cm]                    5.14048714338585E+03
!
! Planck's law constants
!
       real(kind=8),parameter  ::  c1  = 2.0*pi*hplanck*clight**2
       real(kind=8),parameter  ::  c2  = hplanck*clight/Kb
       real(kind=8),parameter  ::  cc2 = c2*1.0E9/560.d0
!
! References for the absorption prefactor
! ---------------------------------------
!
! - J. A. Schellman, Chem. Rev. 75, 323 (1975).  [Eq. 18]
! - N. Lin, F. Santoro, A. Rizzo, Y. Luo, X. Zhao, and V. Barone,
!     J. Phys. Chem. A 113, 4198 (2009), PMID: 19253990. [Eq. 5]
! - A. Baiardi, J. Bloino, and V. Barone, J. Chem. Theory Comput.
!    9, 4097 (2013). [Eq 1]
!
       real(kind=8),parameter  ::  Aabs = (2*pi)**2*Na*10.0*(qe*a0)**2/&
                                   (3.0*log(10.0)*hbar*clight*4*pi*eps0)
!
! References for the emission prefactor
! -------------------------------------
!
! - A. Baiardi, J. Bloino, and V. Barone, J. Chem. Theory Comput.
!    9, 4097 (2013). [Eq 1]
!
       real(kind=8),parameter  ::  Aemi = 2.0*Na*(qe*a0)**2/           &
                                                    (3.0*eps0*clight**3)
!
! References for the Einstein coefficients (for frequencies in cm**-1)
! ----------------------------------------
!
! - A. Requena and J. Zúñiga, Espectroscopia.
!
       real(kind=8),parameter  ::  Bnm = 2.0*pi**2*10.0                &
                                               /(3.0*eps0*hplanck**2)  &
                                                             *(a0*qe)**2
       real(kind=8),parameter  ::  Amn = 8.0*pi*hplanck*Bnm            &
                                                       *(a0*qe)**2*1.0E6
       real(kind=8),parameter  ::  fnm = 4.0*eps0*me*hplanck*Bnm/qe**2 &
                                                            *10.0*clight
       real(kind=8),parameter  ::  eps = Bnm*Na*hplanck/               &
                                                      (log(10.0)*clight)
!
       end module parameters
!
!======================================================================!
