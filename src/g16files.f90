!======================================================================!
!
       module g16files
       implicit none
!
       private
       public  ::  chck_log,                                           &
                   read_log
!
       contains
!
!======================================================================!
!
       subroutine print_badread(inp,key)
!
       use utils,  only: leninp
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
       use utils,  only: uniinp,                                       &
                         lenline
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
       subroutine find_last(key,line,iost)
!
       use utils,  only: uniinp,                                       &
                         lenline
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
         read(uniinp,'(A)',iostat=io) straux
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
       use utils,  only: uniinp,                                       &
                         leninp,                                       &
                         lenline
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
!         
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
!
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
!
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
!
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
       subroutine chck_log(inp,nat)
!
       use utils,     only: uniinp,                                    &
                            leninp,                                    &
                            lenarg,                                    &
                            lenline
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  inp     !  Input file name
       integer,intent(out)               ::  nat     !  Number of atoms
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
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Missing input file'
         write(*,*)
         write(*,'(3X,3(A))')   'Input file ',trim(inp),    &
                                ' not found in the current directory'
         write(*,'(2X,68("="))')
         write(*,*)
         call exit(0)
       end if
!
! Fatal errors check
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
         call exit(0)
       end if
       rewind(uniinp)
!
       call find_key('-- Stationary point found',line,io)
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
         call exit(0)
       end if
       rewind(uniinp)
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
         call exit(0)
       end if
       rewind(uniinp)
! 
! Reading number of atoms !   FLAG: linear molecules not considered
!
       call find_key('NAtoms=',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'NAtoms=')
!
       read(line,*) straux,nat,line
!
       close(uniinp)
!
       return
       end subroutine chck_log
!
!======================================================================!
!
       subroutine read_log(inp,nat,coord,atname,znum,atmass,dof,       &
                           freq,inten,moment,mass,eldeg,Escf)
!
       use parameters
       use utils,     only : uniinp,                                   &
                             lenarg,                                   &
                             lenline,                                  &
                             znum2atname
       use inertia
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                  ::  inp     ! 
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
       character(len=lenline)                       ::  line    !  Line read
       character(len=lenarg)                        ::  straux  !  Auxiliary string
       real(kind=8),dimension(3,3)                  ::  axes    !  Principal axes
       integer                                      ::  io      !  Input/Output status
       integer                                      ::  i,j     !  Indexes
!
! Reading information from Gaussian16 output file
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
!
! Reading converged SCF energy ! FLAG: only supports SCF calculations
!
       call find_last('SCF Done',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'SCF Done')
!
       io = index(line,'=',.TRUE.)
       line = line(io+1:)
!
       read(line,*) Escf,line
!
       rewind(uniinp)
!
! Reading optimized geometry
!
       call find_coord(nat,coord,io)
       if ( io .ne. 0 ) call print_badread(inp,'Standard orientation:')
!
       rewind(uniinp)
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
! Reading vibrational frequencies and IR intensities  ! FLAG: breaks for linear molecules
!
       if ( nat .gt. 2 ) then
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
       end if
!
! Reading molecule information
!
       call find_key('- Thermochemistry -',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'- Thermochemistry -')
!
       read(uniinp,'(A)') line
       read(uniinp,'(A)') line
!
       do i = 1, nat
         read(uniinp,'(A)') line
!
         io = index(line,'r',.TRUE.)
         line = line(io+1:)
         read(line,'(I3,(A))') znum(i),line  
!
         call znum2atname(znum(i),atname(i))
!
         io = index(line,'s',.TRUE.)
         line = line(io+1:)
         read(line,*) atmass(i)
       end do
!
! Reading molecular mass
!
       read(uniinp,'(A)') line
!
       io   = scan(line,':')
       line = line(io+1:)
!
       read(line,*) mass,straux
!
! Computing principal moments of inertia
!
       call inertia_moments(nat,coord*ang2au,atmass,axes,moment)
!
       close(uniinp)
!
       return
       end subroutine read_log
!
!======================================================================!
!
       end module g16files
!
!======================================================================!
