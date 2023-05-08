!======================================================================!
!
       module g16files
!
       use utils,     only: uniinp,leninp,                             &
                            lenarg,lenline
!
       implicit none
!
       private
       public  ::  chk_log,                                            &
                   read_log
!
       contains
!
!======================================================================!
!
       subroutine read_log(fcalc,inp,auxinp,nat,coord,atname,znum,     &
                           atmass,dof,freq,inten,moment,mass,eldeg,    &
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
       call inertia_moments(nat,coord*ang2au,atmass,axes,moment)
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
         call read_ir(inp,nat,dof,freq,inten)
       end if
!
       close(uniinp)
!
       return
       end subroutine read_log
!
!======================================================================!
!
       subroutine chk_log(inp,nat,fcalc,fchck)
!
       use utils,     only: print_missinp
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  inp     !  Input file name
       integer,intent(out)               ::  nat     !  Number of atoms
       character(len=8),intent(in)       ::  fcalc   !  Calculation information flag
       logical,intent(in)                ::  fchck
!
! Local variables
!
       character(len=lenline)            ::  line    !  Line read
       character(len=lenarg)             ::  straux  !  Auxiliary string
       integer                           ::  io      !  Input/Output status
!
! Reading information from Gaussian16 output file
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
!
       if ( io .ne. 0 ) then
         call print_missinp(inp)
       end if
! 
! Reading number of atoms !   TODO: linear molecules not considered
!
       call find_key('NAtoms=',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'NAtoms=')
!
       read(line,*) straux,nat,line
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
         case ('ALL')
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
       end subroutine chk_log
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
       subroutine find_key(key,line,iost)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key     !  Input file name
       character(len=lenline),intent(out)  ::  line    !  Line read
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
         read(uniinp,'(A)',iostat=io) line
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
       character(len=leninp)                      ::  key     !  
       character(len=lenline)                     ::  line    !  Line read
       integer                                    ::  keylen  !  Keyword length
       integer                                    ::  iat     !
       integer                                    ::  io      !  Input/Output status
!
! Reading coordinates
!
       key    = 'Standard orientation:'
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
           read(uniinp,*) line
           read(uniinp,*) line
           read(uniinp,*) line
!
           do iat = 1, nat
!
             read(uniinp,'(A)') line
! Skipping Center Number column
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
! Skipping Atomic Number column
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
! Skipping Atomic Type column
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
! Reading coordinates
             read(line,*) coord(:,iat)
!
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
       call find_key('Error termination',line,io)
       if ( io .eq. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Gaussian job terminated a'//  &
                                                             'bnormally'
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
!~        call find_key('-- Stationary point found',line,io)  ! FIXME: breaks if only frequency calculation
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
       character(len=lenline)       ::  line    !  Line read
       integer                      ::  io      !  Input/Output status
!
! Checking if there are imaginary frequencies present
!
       call find_key('****** ',line,io)
       if ( io .eq. 0 ) then
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
       integer                      ::  io      !  Input/Output status
!
! Reading converged SCF energy ! TODO: only supports SCF calculations
!
       if ( trim(wfn) .eq. 'SCF' ) then
!
         call find_last(uni,'SCF Done',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'SCF Done')
!
         io = index(line,'=',.TRUE.)
         line = line(io+1:)
!
         read(line,*) Escf,line
!
       else if ( trim(wfn) .eq. 'EXTERNAL' ) then
!
         call find_last(uni,'Energy=',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'Energy=')
!
         io = index(line,'=',.FALSE.)
         line = line(io+1:)
!
         read(line,*) Escf,line
!
       else if ( trim(wfn) .eq. 'COUNTERPOISE' ) then
!
         call find_last(uni,'Counterpoise corrected energy',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'Counterpoise cor'//  &
                                             'rected energy')
!
         io = index(line,'=',.TRUE.)
         line = line(io+1:)
!
         read(line,*) Escf
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
       if ( io .ne. 0 ) call print_badread(inp,'Standard orientation:')
!
       return
       end subroutine read_opt
!
!======================================================================!
!
       subroutine read_ir(inp,nat,dof,freq,inten)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)              ::  inp     ! 
       real(kind=8),dimension(dof),intent(out)  ::  freq    !  
       real(kind=8),dimension(dof),intent(out)  ::  inten   !
       integer,intent(in)                       ::  dof     !  
       integer,intent(in)                       ::  nat     !  Number of atoms
!
! Local variables
!
       character(len=lenline)                   ::  line    !  Line read
       integer                                  ::  io      !  Input/Output status
       integer                                  ::  i,j     !  Indexes
!
! Reading vibrational frequencies and IR intensities  ! FIXME: breaks for linear molecules
!
       j = 0
       do i = 1, nat - 2
         call find_key('Frequencies',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'Frequencies')
!
         io   = scan(line,'-',.TRUE.)
         line = line(io+1:)
!
         read(line,*) freq(j+1:j+3)
!
         call find_key('IR Inten',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'IR Inten')
!
         io   = scan(line,'-',.TRUE.)
         line = line(io+1:)
!
         read(line,*) inten(j+1:j+3)
!
         j = j + 3
       end do
!
       inten(:) = 100.0d0/log(10.0)*inten(:)/freq(:)
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
       call find_key('Charge',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'Charge')
!
       io   = scan(line,'=',.TRUE.)
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
       character(len=lenline)                       ::  line    !  Line read
       integer                                      ::  io      !  Input/Output status
       integer                                      ::  iat     !  Index
!
! Reading atomic numbers
!
       call find_key('Standard orientation:',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'Standard orientation:')
!
       read(uniinp,*) line
       read(uniinp,*) line
       read(uniinp,*) line
       read(uniinp,*) line
!
       do iat = 1, nat
!
         read(uniinp,'(A)') line
! Skipping Center Number column
         line = adjustl(line)
         io   = scan(line,' ')
         line = line(io+1:)
! Reading Atomic Number column
         line = adjustl(line)
         io   = scan(line,' ')
!
         read(line(:io-1),*) znum(iat)
!
       end do
!
! Setting molecule information
!
       mass = 0.0d0
!
       do iat = 1, nat
         call znum2atname(znum(iat),atname(iat))           
         call znum2atmass(znum(iat),atmass(iat))  
         mass = mass + atmass(iat)
       end do
!
       return
       end subroutine read_znum
!
!======================================================================!
!
       end module g16files
!
!======================================================================!
