!======================================================================!
!
       module datatypes
!
       use utils,    only: leninp,                                     &
                           lename
!
       implicit none
!
       type thermo
         character(len=leninp)                      ::  inp      !  Input file name
         character(len=8)                           ::  wfn      !  Wavefunction type
         real(kind=8),dimension(:,:),allocatable    ::  coord    !  Coordinates
         real(kind=8),dimension(:),allocatable      ::  freq     !  Frequencies
         real(kind=8),dimension(:),allocatable      ::  inten    !  IR Intensities
         real(kind=8),dimension(3)                  ::  moment   !  Princial moments of inertia
         real(kind=8)                               ::  Escf     !  Potential energy
         real(kind=8)                               ::  ZPVE     !  Zero-point energy
         real(kind=8)                               ::  D0       !  Dissociation energy (0 K)
         real(kind=8)                               ::  qtrans   !  Translational partition function
         real(kind=8)                               ::  qrot     !  Rotational partition function
         real(kind=8)                               ::  qvib     !  Vibrational partition function
         real(kind=8)                               ::  qel      !  Electronic partition function
         real(kind=8)                               ::  G        !  Thermal Gibbs free energy
         real(kind=8)                               ::  Gtrans   !  Translational contribution to Gtherm
         real(kind=8)                               ::  Grot     !  Rotational contribution to Gtherm
         real(kind=8)                               ::  Gvib     !  Vibrational contribution to Gtherm
         real(kind=8)                               ::  E        !  Thermal internal energy
         real(kind=8)                               ::  Etrans   !  Translational contribution to Etherm
         real(kind=8)                               ::  Erot     !  Rotational contribution to Etherm
         real(kind=8)                               ::  Evib     !  Vibrational contribution to Etherm
         real(kind=8)                               ::  H        !  Thermal enthalpy
         real(kind=8)                               ::  Htrans   !  Translational contribution to Htherm
         real(kind=8)                               ::  Hrot     !  Rotational contribution to Htherm
         real(kind=8)                               ::  Hvib     !  Vibrational contribution to Htherm
         real(kind=8)                               ::  S        !  Total entropy
         real(kind=8)                               ::  Sel      !  Electronic contribution to Stherm
         real(kind=8)                               ::  Strans   !  Translational contribution to Stherm
         real(kind=8)                               ::  Srot     !  Rotational contribution to Stherm
         real(kind=8)                               ::  Svib     !  Vibrational contribution to Stherm
         integer                                    ::  nequi    !  Number of equivalent conformations
         integer                                    ::  symnum   !  Symmetry number
         integer                                    ::  dof      !  Degrees of freedom 
         integer                                    ::  rotdof   !  Rotational degrees of freedom 
         logical                                    ::  linear   !  Linear molecule flag
         logical                                    ::  chiral   !  Chirality flag
       end type thermo  
!
       type molecule
         type(thermo),dimension(:),allocatable      ::  conf     !  Conformer information
         character(len=5),dimension(:),allocatable  ::  atname   !  Atom names
         character(len=lename)                      ::  molname  !  Molecule name
         character(len=8)                           ::  phase    !  Component phase
         real(kind=8),dimension(:),allocatable      ::  pop      !  Conformers populations
         real(kind=8),dimension(:),allocatable      ::  atmass   !  Atomic masses
         real(kind=8)                               ::  mass     !  Molecule mass
         real(kind=8)                               ::  qtot     !  Multistructural partition function
         real(kind=8)                               ::  Detot    !  Multistructural electronic energy
         real(kind=8)                               ::  D0tot    !  Multistructural zero kelvin energy
         real(kind=8)                               ::  Etot     !  Multistructural energy
         real(kind=8)                               ::  Htot     !  Multistructural enthalpy
         real(kind=8)                               ::  Stot     !  Multistructural entropy
         real(kind=8)                               ::  Gtot     !  Multistructural free energy
         real(kind=8)                               ::  popsum   !  Sum of relative populations
         integer,dimension(:),allocatable           ::  Znum     !  Atomic number
         integer,dimension(:),allocatable           ::  frag     !  Number of indistinguishable fragments
         integer                                    ::  nfrag    !  Number of different fragments
         integer                                    ::  npermu   !  Number of fragments permutations
         integer                                    ::  nat      !  Number of atoms
         integer                                    ::  nconf    !  Number of conformers
       end type molecule
!
       type reaction
         real(kind=8)                               ::  dDe      !  Potential energy change
         real(kind=8)                               ::  dD0      !  Zero kelvin energy change
         real(kind=8)                               ::  dE       !  Internal energy change
         real(kind=8)                               ::  dH       !  Enthalpy change
         real(kind=8)                               ::  dS       !  Entropu energy change
         real(kind=8)                               ::  dG       !  Free energy change
         real(kind=8)                               ::  Keq      !  Equilibrium constant
         integer,dimension(:),allocatable           ::  nu       !  Stoichiometric coefficients
       end type reaction
!
       end module datatypes
!
!======================================================================!
