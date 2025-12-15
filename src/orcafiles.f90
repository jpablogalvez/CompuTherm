!======================================================================!
!
       module orcafiles
!
       use utils,     only: uniinp,leninp,                             &
                            lenarg,lenline
!
       implicit none
!
       private
       public  ::  chk_orca,                                           &
                   read_orca,                                          &
                   read_esorca
!
       contains
!
!======================================================================!
!
       subroutine read_orca(fcalc,inp,auxinp,nat,coord,atname,znum,    &
                            atmass,dof,freq,inten,moment,mass,eldeg,   &
                            Escf,wfn,schm)
!
       use datatypes
       use parameters
       use inertia
       use utils
!
       implicit none
!
! Input/output variables
!
       type(scheme),intent(in)                      ::  schm    !
       character(len=*),intent(in)                  ::  inp     !
       character(len=*),dimension(:),intent(in)     ::  auxinp  !
       character(len=*),intent(in)                  ::  wfn     !  Wavefunction information
       character(len=8),intent(in)                  ::  fcalc   !  Calculation information flag
       character(len=5),dimension(nat),intent(out)  ::  atname  !  Atom names
       real(kind=8),dimension(3,nat),intent(out)    ::  coord   !
       real(kind=8),dimension(nat),intent(out)      ::  atmass  !  Atomic masses
       real(kind=8),dimension(dof),intent(out)      ::  freq    !
       real(kind=8),dimension(dof),intent(out)      ::  inten   !
       real(kind=8),dimension(3),intent(out)        ::  moment  !
       real(kind=8),intent(out)                     ::  Escf    !
       real(kind=8),intent(out)                     ::  mass    !  Molecule mass
       real(kind=8),intent(out)                     ::  eldeg   !
       integer,dimension(nat),intent(out)           ::  znum    !  Atomic number
       integer,intent(in)                           ::  nat     !  Number of atoms
       integer,intent(in)                           ::  dof     !
!
! Local variables
!
       real(kind=8)                                 ::  Ell     !
       real(kind=8)                                 ::  Ehl     !
       real(kind=8),dimension(3,3)                  ::  axes    !  Principal axes
       integer                                      ::  io      !  Input/Output status
!
! Reading information from Gaussian16 output file
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
       if ( io .ne. 0 ) call print_missinp(inp)
!
! Reading electronic degeneracy
!
       call read_eldeg(inp,eldeg)
!
! Reading molecular information
!
       call read_znum(inp,nat,znum,atname,atmass,mass)
!
! Reading converged SCF energy
!
       select case (trim(schm%fschm))
         case ('NORMAL')
!
           call read_energy(uniinp,inp,Escf,wfn)
!
         case ('HL')
!
           open(unit=uniinp+1,file=trim(auxinp(1)),action='read',      &
                status='old',iostat=io)
           if ( io .ne. 0 ) call print_missinp(auxinp(1))
!
           call read_energy(uniinp+1,auxinp(1),Escf,schm%wfn)
!
           close(uniinp+1)
!
         case ('LLSOL')
!
           open(unit=uniinp+1,file=trim(auxinp(1)),action='read',      &
                status='old',iostat=io)
           if ( io .ne. 0 ) call print_missinp(auxinp(1))
!
           open(unit=uniinp+2,file=trim(auxinp(2)),action='read',      &
                status='old',iostat=io)
           if ( io .ne. 0 ) call print_missinp(auxinp(2))
!
           call read_energy(uniinp,inp,Escf,wfn)
           call read_energy(uniinp+1,auxinp(1),Ell,wfn)
           call read_energy(uniinp+2,auxinp(2),Ehl,schm%wfn)
!
           Escf = Escf - Ell + Ehl
!
           close(uniinp+1)
           close(uniinp+2)
!
         case ('HLSOL')
!
           open(unit=uniinp+1,file=trim(auxinp(1)),action='read',      &
                status='old',iostat=io)
           if ( io .ne. 0 ) call print_missinp(auxinp(1))
!
           open(unit=uniinp+2,file=trim(auxinp(2)),action='read',      &
                status='old',iostat=io)
           if ( io .ne. 0 ) call print_missinp(auxinp(2))
!
           open(unit=uniinp+3,file=trim(auxinp(3)),action='read',      &
                status='old',iostat=io)
           if ( io .ne. 0 ) call print_missinp(auxinp(2))
!
           call read_energy(uniinp+1,auxinp(1),Escf,wfn)
           call read_energy(uniinp+2,auxinp(2),Ell,wfn)
           call read_energy(uniinp+3,auxinp(3),Ehl,schm%wfn)
!
           Escf = Escf - Ell + Ehl
!
           close(uniinp+1)
           close(uniinp+2)
           close(uniinp+3)
!
       end select
!
       rewind(uniinp)
!
! Reading optimized geometry
!
       call read_opt(inp,nat,coord)
!
! Computing the principal moments of inertia
!
       coord = coord*ang2au
       call inertia_moments(nat,coord,atmass,axes,moment)
!
       if ( (trim(fcalc).eq.'EONLY') .or. (trim(fcalc).eq.'OPT') ) then
         close(uniinp)
         return
       end if
!
       rewind(uniinp)
!
! Reading vibrational frequencies and IR intensities  ! FIXME: breaks for linear molecules
!
       if ( nat .gt. 2 ) then
         call read_ir(inp,dof,freq,inten)
       end if
!
       close(uniinp)
!
       return
       end subroutine read_orca
!
!======================================================================!
!
       subroutine read_esorca(inp,nband,vener,tdip)
!
       use parameters
       use utils
!
       implicit none
!
! Input/Output variables
!
       character(len=leninp),intent(in)                   ::  inp     !  Input file name
       real(kind=8),dimension(:),allocatable,intent(out)  ::  vener   !  Vertical energies
       real(kind=8),dimension(:),allocatable,intent(out)  ::  tdip    !  Transition dipole moment
       integer,intent(out)                                ::  nband   !  Number of roots
!
! Local variables
!
       character(len=lenline)                             ::  line    !  Line read
       character(len=lenline)                             ::  key
       real(kind=4)                                       ::  s1,s2
       integer                                            ::  i1
       integer                                            ::  keylen  !  Keyword length
       integer                                            ::  io      !  Input/Output status
       integer                                            ::  iost    !  Input/Output status
       integer                                            ::  i       !  Index
!
! Reading vertical transition energies and transition dipole moments
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
                                                 status='old',iostat=io)
!
       if ( io .ne. 0 ) call print_missinp(trim(inp))
!
       call find_key(uniinp,'Number of roots to be determined',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'Number of roots to'//  &
                                                       ' be determined')
!
       io   = scan(line,'.',.TRUE.) 
       line = line(io+1:)
!
       read(line,*) nband
!
       allocate(vener(nband),tdip(nband))
!
       key = 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE M'//  &
                                                                'OMENTS'
       keylen = len_trim(key)
       iost   = 1
!
       do
         read(uniinp,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         line = adjustl(line)
         if ( line(:keylen) .eq. trim(key) ) then
           iost = 0
!
           read(uniinp,'(A)') line
           read(uniinp,'(A)') line
           read(uniinp,'(A)') line
           read(uniinp,'(A)') line
!
           do i = 1, nband
             read(uniinp,*) i1,vener(i),s1,s2,tdip(i)
           end do
         end if
       end do
!
       if ( iost .ne. 0 ) call print_badread(inp,'ABSORPTION SPECT'//  &
                                   'RUM VIA TRANSITION ELECTRIC DIPOLE')
!
       close(uniinp)
!
       return
       end subroutine read_esorca
!
!======================================================================!
!
       subroutine chk_orca(inp,nat,fcalc,fchck)
!
       use utils,     only: print_missinp
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  inp     !  Input file name
       integer,intent(out)               ::  nat     !  Number of atoms
       character(len=*),intent(in)       ::  fcalc   !  Calculation information flag
       logical,intent(in)                ::  fchck
!
! Local variables
!
       character(len=lenline)            ::  line    !  Line read
       integer                           ::  io      !  Input/Output status
!
! Reading information from Gaussian16 output file
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
!
       if ( io .ne. 0 ) call print_missinp(inp)
!
! Reading number of atoms !   TODO: linear molecules not considered
!
       call find_key(uniinp,'Number of atoms',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'Number of atoms')
!
       io   = scan(line,'.',.TRUE.) 
       line = line(io+1:)
!
       read(line,*) nat
!
       if ( fchck ) then
         rewind(uniinp)
       else
         close(uniinp)
         return
       end if
!
! Fatal errors check
!
       call chk_errter(inp)
!
       rewind(uniinp)
!
       select case ( trim(fcalc) )
         case ('OPT')
!
           call chk_opt(inp)
!
         case ('FREQ')
!
           call chk_imagi(inp)
!
         case ('ALL')  ! FIXME: check optimization and frequencies
!
!~            call chk_opt(inp)
!
!~            rewind(uniinp)
!
           call chk_imagi(inp)
!
         case default

       end select
!
       close(uniinp)
!
       return
       end subroutine chk_orca
!
!======================================================================!
!
       subroutine print_badread(inp,key)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)       ::  key     !  Input file name
       character(len=leninp),intent(in)  ::  inp     !  Input file name
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)')      'ERROR:  Missing information in the '//  &
                              'input file'
       write(*,*)
       write(*,'(3X,2(A))')   'Keyword    : ',trim(adjustl(key))
       write(*,'(3X,2(A))')   'Input file : ',trim(inp)
       write(*,'(2X,68("="))')
       write(*,*)
       call exit(0)
!
       return
       end subroutine print_badread
!
!======================================================================!
!
       subroutine find_key(uni,key,line,iost)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key     !  Input file name
       character(len=lenline),intent(out)  ::  line    !  Line read
       integer,intent(in)                  ::  uni     !
       integer,intent(out)                 ::  iost    !  Reading status
!
! Local variables
!
       integer                             ::  keylen  !  Keyword length
       integer                             ::  io      !  Input/Output status
!
! Finding the first line starting with the specified keyword
!
       keylen = len(key)
       iost   = 0
!
       do
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         line = adjustl(line)
         if ( line(:keylen) .eq. trim(key) ) return
       end do
!
       iost = 1
!
       return
       end subroutine find_key
!
!======================================================================!
!
       subroutine find_last(uni,key,line,iost)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key     !  Input file name
       character(len=lenline),intent(out)  ::  line    !  Line read
       integer,intent(in)                  ::  uni     !
       integer,intent(out)                 ::  iost    !  Reading status
!
! Local variables
!
       character(len=lenline)              ::  straux  !  Auxliary string
       integer                             ::  keylen  !  Keyword length
       integer                             ::  io      !  Input/Output status
!
! Finding the last line starting with the specified keyword
!
       keylen = len_trim(key)
       iost   = 1
!
       do
         read(uni,'(A)',iostat=io) straux
         if ( io /= 0 ) exit
         straux = adjustl(straux)
         if ( straux(:keylen) .eq. trim(key) ) then
           iost = 0
           line = straux
         end if
       end do
!
       return
       end subroutine find_last
!
!======================================================================!
!
       subroutine find_coord(nat,coord,iost) 
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(3,nat),intent(out)  ::  coord   !
       integer,intent(in)                         ::  nat     !
       integer,intent(out)                        ::  iost    !  Reading status
!
! Local variables
!
       character(len=lenline)                     ::  key     !
       character(len=lenline)                     ::  line    !  Line read
       character(len=lenarg)                      ::  str1    !  Auxiliary string
       integer                                    ::  keylen  !  Keyword length
       integer                                    ::  iat     !
       integer                                    ::  io      !  Input/Output status
!
! Reading coordinates
!
       key    = 'CARTESIAN COORDINATES (ANGSTROEM)'
       keylen = len_trim(key)
       iost   = 1
!
       do
         read(uniinp,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         line = adjustl(line)
         if ( line(:keylen) .eq. trim(key) ) then
           iost = 0
!
           read(uniinp,*) line
!
           do iat = 1, nat
             read(uniinp,*) str1,coord(:,iat)
           end do
         end if
       end do
!
       return
       end subroutine find_coord
!
!======================================================================!
!
       subroutine chk_errter(inp)
!
       use utils,     only: print_end
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  inp     !
!
! Local variables
!
       character(len=lenline)       ::  line    !  Line read
       integer                      ::  io      !  Input/Output status
!
! Checking if the calculation therminated abnormally
!
       call find_key(uniinp,'****ORCA TERMINATED NORMALLY****',line,io)
!
       if ( io .ne. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  ORCA job terminated abnor'//  &
                                                                 'mally'
         write(*,*)
         write(*,'(3X,A)')      'Please, check the following input'//  &
                                                                ' file:'
         write(*,*)
         write(*,'(4X,A)')       trim(inp)
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
       return
       end subroutine chk_errter
!
!======================================================================!
!
       subroutine chk_opt(inp)
!
       use utils,     only: print_end
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  inp     !
!
! Local variables
!
       character(len=lenline)       ::  line    !  Line read
       integer                      ::  io      !  Input/Output status
!
! Checking if the optimization therminated abnormally
!
       io = 0
       call find_key(uniinp,'***********************HURRAY',line,io)  ! FIXME: breaks if only frequency calculation
!
       if ( io .ne. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Geometry optimization ter'//  &
                                                    'minated abnormally'
         write(*,*)
         write(*,'(3X,A)')      'Please, check the following input'//  &
                                                                ' file:'
         write(*,*)
         write(*,'(4X,A)')       trim(inp)
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
       return
       end subroutine chk_opt
!
!======================================================================!
!
       subroutine chk_imagi(inp)
!
       use utils,     only: print_end
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  inp     !
!
! Local variables
!
       real(kind=8)                 ::  dpaux   !
       character(len=lenline)       ::  line    !  Line read
       integer                      ::  io      !  Input/Output status
!
! Checking if there are imaginary frequencies present
!
       call find_key(uniinp,'VIBRATIONAL FREQUENCIES',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'VIBRATIONAL FREQUE'//  &
                                                                'NCIES')
!
       call find_key(uniinp,'6:',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'6:')
!
       read(line,*) line,dpaux
!
       if ( dpaux .lt. 1.0E-6 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Geometry optimization did'//  &
                                       ' not finished in a true minimum'
         write(*,*)
         write(*,'(3X,A)')      'Please, check the following input'//  &
                                                                ' file:'
         write(*,*)
         write(*,'(4X,A)')       trim(inp)
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
       return
       end subroutine chk_imagi
!
!======================================================================!
!
       subroutine read_energy(uni,inp,Escf,wfn)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  inp     !
       character(len=*),intent(in)  ::  wfn     !  Wavefunction information
       real(kind=8),intent(out)     ::  Escf    !
       integer,intent(in)           ::  uni     !
!
! Local variables
!
       character(len=lenline)       ::  line    !  Line read
       character(len=lenarg)        ::  str1    !  Auxiliary string
       character(len=lenarg)        ::  str2    !  Auxiliary string
       character(len=lenarg)        ::  str3    !  Auxiliary string
       character(len=lenarg)        ::  str4    !  Auxiliary string
       integer                      ::  io      !  Input/Output status
!
! Reading converged SCF energy ! TODO: only supports SCF calculations
!
       if ( trim(wfn) .eq. 'SCF' ) then
!
         call find_last(uni,'FINAL SINGLE POINT ENERGY',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'FINAL SINGLE POI'//  &
                                                            'NT ENERGY')
!
         read(line,*) str1,str2,str3,str4,Escf
!
       else if ( trim(wfn) .eq. 'CC' ) then
!
         call find_last(uni,'FINAL SINGLE POINT ENERGY',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'FINAL SINGLE POI'//  &
                                                            'NT ENERGY')
!
         read(line,*) str1,str2,str3,str4,Escf
!
       else if ( trim(wfn) .eq. 'EXTERNAL' ) then
!
         STOP 'External calculation for ORCA output is not yet imp'//  &
                                                             'lemented!'
!
!~          call find_last(uni,'Energy=',line,io)
!~          if ( io .ne. 0 ) call print_badread(inp,'Energy=')
!~ !
!~          io = index(line,'=',.FALSE.)
!~          line = line(io+1:)
!~ !
!~          read(line,*) Escf,line
!
       else if ( trim(wfn) .eq. 'COUNTERPOISE' ) then
!
         STOP 'Counterpoise correction for ORCA output is not yet '//  &
                                                          'implemented!'
!
!~          call find_last(uni,'Counterpoise corrected energy',line,io)
!~          if ( io .ne. 0 ) call print_badread(inp,'Counterpoise cor'//  &
!~                                              'rected energy')
!~ !
!~          io = index(line,'=',.TRUE.)
!~          line = line(io+1:)
!~ !
!~          read(line,*) Escf
!
       end if
!
       return
       end subroutine read_energy
!
!======================================================================!
!
       subroutine read_opt(inp,nat,coord)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                ::  inp     !
       real(kind=8),dimension(3,nat),intent(out)  ::  coord   !
       integer,intent(in)                         ::  nat     !  Number of atoms
!
! Local variables
!
       integer                                    ::  io      !  Input/Output status
!
! Reading optimized geometry
!
       call find_coord(nat,coord,io)
       if ( io .ne. 0 ) call print_badread(inp,'CARTESIAN COORDINA'//  &
                                                      'TES (ANGSTROEM)')
!
       return
       end subroutine read_opt
!
!======================================================================!
!
       subroutine read_ir(inp,dof,freq,inten)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)              ::  inp     !
       real(kind=8),dimension(dof),intent(out)  ::  freq    !
       real(kind=8),dimension(dof),intent(out)  ::  inten   !
       integer,intent(in)                       ::  dof     !
!
! Local variables
!
       character(len=lenline)                   ::  line    !  Line read
       character(len=lenarg)                    ::  str1    !  Auxiliary string
       character(len=lenarg)                    ::  str2    !  Auxiliary string
       integer                                  ::  io      !  Input/Output status
       integer                                  ::  i       !  Index
!
! Reading vibrational frequencies and IR intensities  ! FIXME: breaks for linear molecules
!
       call find_key(uniinp,'IR SPECTRUM',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'IR SPECTRUM')
!
       do i = 1, 5
         read(uniinp,*)
       end do
!
       do i = 1, dof
!
         read(uniinp,*) str1,freq(i),str2,inten(i)               
!
       end do
!
       inten(:) = 100.0d0/log(10.0)*inten(:)/freq(:)   ! Gaussian formula
!~        inten(:) = inten(:)/freq(:)  !  JC formula without concentration
!
       return
       end subroutine read_ir
!
!======================================================================!
!
       subroutine read_eldeg(inp,eldeg)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  inp     !
       real(kind=8),intent(out)     ::  eldeg   !
!
! Local variables
!
       character(len=lenline)       ::  line    !  Line read
       integer                      ::  io      !  Input/Output status
!
! Reading electronic degeneracy
!
       call find_key(uniinp,'Multiplicity',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'Multiplicity')
!
       io   = scan(line,'.',.TRUE.) 
       line = line(io+1:)
!
       read(line,*) eldeg
!
       return
       end subroutine read_eldeg
!
!======================================================================!
!
       subroutine read_znum(inp,nat,znum,atname,atmass,mass)
!
       use utils,     only : znum2atname,                              &
                             znum2atmass
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                  ::  inp     !
       integer,intent(in)                           ::  nat     !  Number of atoms
       character(len=5),dimension(nat),intent(out)  ::  atname  !  Atom names
       real(kind=8),dimension(nat),intent(out)      ::  atmass  !  Atomic masses
       real(kind=8),intent(out)                     ::  mass    !  Molecule mass
       integer,dimension(nat),intent(out)           ::  znum    !  Atomic number
!
! Local variables
!
       real(kind=8)                                 ::  dpaux   !
       character(len=lenline)                       ::  line    !  Line read
       integer                                      ::  io      !  Input/Output status
       integer                                      ::  iat     !  Index
!
! Reading atomic numbers
!
       call find_key(uniinp,'CARTESIAN COORDINATES (A.U.)',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'CARTESIAN COORDINA'//  &
                                                           'TES (A.U.)')
!
       read(uniinp,*) line
       read(uniinp,*) line
!
       do iat = 1, nat
!
         read(uniinp,*) line,atname(iat),dpaux 
         znum(iat) = int(dpaux)
!
       end do
!
! Setting molecule information
!
       mass = 0.0d0
!
       do iat = 1, nat
         call znum2atmass(znum(iat),atmass(iat))
         mass = mass + atmass(iat)
       end do
!
       return
       end subroutine read_znum
!
!======================================================================!
!
       end module orcafiles
!
!======================================================================!
