!======================================================================!
!
       program g16_computherm
!
       use parameters
       use lineshape
       use datatypes
       use utils,      only:  leninp,lenout,                           &
                              line_str,                                &
                              line_int,                                &
                              line_dp,                                 &
                              line_dvec,                               &
                              line_ivec,                               &
                              line_log,                                &
                              print_title,                             &
                              print_titleint,                          &
                              print_start,                             &
                              print_end
       use input
       use statmech
       use mathtools,  only:  factorial 
!
       implicit none
!
       include 'timings.h'
!
       type(molecule),dimension(:),allocatable  ::  mol      !  Molecules information
       type(reaction),dimension(:),allocatable  ::  reac     !  Reactions information
       type(scheme)                             ::  schm     !
       character(len=leninp)                    ::  inp      !  Input file name
       character(len=lenout)                    ::  outp     !  Output file name
       character(len=8)                         ::  ffree    !  Free energy calculation flag
       character(len=8)                         ::  fentha   !  Enthalpy calculation flag
       character(len=8)                         ::  fentro   !  Entropy calculation flag
       character(len=8)                         ::  fcalc    !  Calculation information flag
       character(len=8)                         ::  forder   !  Conformations order flag
       character(len=8)                         ::  frota    !  Rotation method flag
       real(kind=8)                             ::  temp     !  Temperature (K)
       real(kind=8)                             ::  pres     !  Pressure (atm)
       real(kind=8)                             ::  volu     !  Volume (m**3)
       real(kind=8)                             ::  cutoff   !  Cutoff frequency
       real(kind=8)                             ::  fact     !  Frequencies scaling factor
       real(kind=8)                             ::  hwhm     !  Half width at half maximum
       real(kind=8)                             ::  rmsdmax  !  Maximum value for RMSD
       real(kind=8)                             ::  maemax   !  Maximum value for MAE
       real(kind=8)                             ::  baemax   !  Maximum value for BAE
       real(kind=8)                             ::  dsolv    !
       real(kind=8)                             ::  msolv    !
       integer,dimension(:,:),allocatable       ::  order    !  Energy-based order
       integer                                  ::  nmol     !  Number of chemical species
       integer                                  ::  imol     !  Molecule index
       integer                                  ::  mconf    !  Maximum number of conformers
       integer                                  ::  iconf    !  Conformer index
       integer                                  ::  ifrag    !  Fragment index
       integer                                  ::  nreac    !  Number of reactions
       integer                                  ::  iline    !  Broading function
       integer                                  ::  lin      ! 
       integer                                  ::  lfin     !  
       logical                                  ::  fchck    !  
       logical                                  ::  doir     !  IR calculation flag
       logical                                  ::  fsoln    !  Standard state flag
       logical                                  ::  fenan    !  Enantiomers calculation flag
       logical                                  ::  fpermu   !  Permutations calculation flag
       logical                                  ::  fscreen  !  Screening calculation flag
       logical                                  ::  doconf   !  Conformational analysis flag
       logical                                  ::  doequi   !  Equilibrium calculation flag
       logical                                  ::  debug    !  Debug mode
!
! Declaration of external functions
!
       logical                                  ::  chkenan  !
!
! Printing header
!
       write(*,*)
       write(*,'(5X,82("#"))')
       write(*,'(5X,10("#"),3X,A,3X,13("#"))') 'CompuTherm - Compu'//  &
                                   'tational Thermochemistry Calculator'
       write(*,'(5X,82("#"))')   
       write(*,*)
       write(*,'(9X,A)') 'Welcome to CompuTherm, a very simple pro'//  &
                                           'gram for the calculation of'
       write(*,'(9X,A)') ' thermodynamic properties from electroni'//  &
                                              'c structure calculations'
       write(*,*)
       write(*,'(1X,90("="))')
       write(*,*)    
!
       call system_clock(count_max=count_max,count_rate=count_rate)
       call system_clock(t1)  
!
       call print_start()
!
! Initializing timings
!
       tcpu  = 0.0d0
!
       lin  = 45
       lfin = 90
!
       doconf = .FALSE.
       doequi = .FALSE.
!
! Reading command line options
!
       call command_line(inp,outp,temp,pres,volu,debug)
!
! Reading general input file
!
       write(*,'(8X,54("*"))')
       write(*,'(8X,10("*"),X,("General input processing section")'//  &
                                                          ',X,10("*"))')
       write(*,'(8X,54("*"))')
       write(*,*)
!
       call read_inp(inp,nmol,mol,cutoff,fact,hwhm,iline,doir,ffree,   &
                     fentha,fentro,forder,fcalc,fchck,fenan,fpermu,    &
                     fsoln,nreac,reac,mconf,fscreen,frota,rmsdmax,     &
                     maemax,baemax,dsolv,msolv,schm)
!
       allocate(order(nmol,mconf))
       order(:,:) = -1
!
       if ( fsoln ) volu = 1.0E-3
!
! Printing summary of the general input file information
!
       call print_title(6,1,'Files information','=')
       write(*,*)
!       
       call line_str(6,2,'General input file name',lin,':',            &
                     trim(inp),lfin)
       call line_str(6,2,'Output file name',lin,':',                   &
                     trim(outp),lfin)
!
       if ( debug ) then
         call line_str(6,2,'Debug mode',lin,':','ON',lfin)
       else
         call line_str(6,2,'Debug mode',lin,':','OFF',lfin)
       end if
       write(*,*)
!
       call print_title(6,1,'Calculation type information','=')
       write(*,*)
!
       call line_str(6,2,'Type of calculation requested',lin,':',      &
                     trim(fcalc),lfin)
       call line_str(6,2,'Calculation scheme',lin,':',                 &
                     trim(schm%fschm),lfin)
       if ( trim(schm%fschm) .ne. 'NORMAL' ) call line_str(6,2,        &
                  'High level wavefunction',lin,':',trim(schm%wfn),lfin)
       write(*,*)
!
       if ( (trim(fcalc).eq.'ALL') .or. (trim(fcalc).eq.'FREQ') ) then
         write(*,'(2X,A)') 'Vibrational frequencies (IR spectra)'
         call line_dp(6,3,'Frequencies scaling factor',lin,':',        &
                      'F6.4',fact,lfin)
         if ( doir ) then
           if ( iline .eq. 1 ) then
             call line_str(6,2,'Broadening function',lin,':',          &
                           'Gaussian',lfin)
           else if ( iline .eq. 2 ) then
             call line_str(6,2,'Broadening function',lin,':',          &
                           'Lorentzian',lfin)
           else if ( iline .eq. 3 ) then
             call line_str(6,2,'Broadening function',lin,':',          &
                           'Psudo-Voigt',lfin)
           end if
!
           if ( (iline.eq.1) .or. (iline.eq.2) ) then
             call line_dp(6,2,'Half width at half maximum (cm**-1)',   &
                          lin,':','F12.2',hwhm,lfin)
           else if ( iline .eq. 3 ) then
!~              call line_dp(6,3,'Gaussian HWHM (cm**-1)',lin,':',        &  ! TODO: create ghwhm and lhwhm
!~                           'F12.4',ghwhm,lfin) 
!~              call line_dp(6,3,'Lorentzian HWHM (cm**-1)',lin,':',      &
!~                           'F12.4',lhwhm,lfin) 
           end if
         end if
         write(*,*)
       end if
!
       if ( trim(fcalc) .eq. 'ALL' ) then
         write(*,'(2X,A)') 'Thermodynamic equilibrium properties'
         call line_str(6,3,'Type of free energy calculation',lin,':',  &
                     trim(ffree),lfin)
         call line_str(6,3,'Type of enthalpy calculation',lin,':',     &
                     trim(fentha),lfin)
         call line_str(6,3,'Type of entropy calculation',lin,':',      &
                     trim(fentro),lfin)
         call line_log(6,3,'Enantiomers calculaton',lin,':',fenan,lfin)
         call line_log(6,3,'Indistinguishable aggregates',lin, &
                     ':',fpermu,lfin)
         write(*,*)
!
         write(*,'(2X,A)') 'Standard state'
         call line_dp(6,3,'Temperature (K)',lin,':','F12.2',temp,lfin)
         call line_dp(6,3,'Pressure (atm)',lin,':','F12.4',pres,lfin)
         call line_dp(6,3,'Volume (L)',lin,':','F12.4',volu*1000.0d0,  &
                      lfin)
         call line_log(6,3,'Standard state 1M',lin,':',fsoln,lfin)
         write(*,*) 
       end if
!
       if ( fscreen ) then
         write(*,'(2X,A)') 'Screening information'
         call line_str(6,3,'Superposition method',lin,':',             &
                       trim(frota),lfin)
         call line_dp(6,3,'Maximum value for RMSD',lin,':','D10.4',    &
                      rmsdmax,lfin)
         call line_dp(6,3,'Maximum value for MAE',lin,':','D10.4',     &
                      maemax,lfin)
         call line_dp(6,3,'Maximum value for BAE',lin,':','D10.4',     &
                      baemax,lfin)
         write(*,*)
       end if
!
! Processing Gaussian16 input files
!
       call read_g16(nmol,mol,fcalc,fchck,schm)
!
! Printing summary of the Gaussian16 information
!
       if ( debug ) then
         call print_title(6,1,'Gaussian16 input files information','=')
         write(*,*)
       end if
!
       do imol = 1, nmol
!
         if ( fpermu ) then
           mol(imol)%npermu = 1
           do ifrag = 1, mol(imol)%nfrag
             mol(imol)%npermu = mol(imol)%npermu                       &
                                       *factorial(mol(imol)%frag(ifrag))
           end do
         end if
!
         if ( debug ) then
           call print_titleint(6,3,'Information of molecular system',1,&
                               'I3',imol,'-')
           write(*,*)
           call line_str(6,2,'Molecule name',lin,':',                  &
                         trim(mol(imol)%molname),lfin)
           call line_str(6,2,'Individual calculation scheme',lin,':',  &
                     trim(schm%fschm),lfin)
           if ( trim(mol(imol)%schm%fschm) .ne. 'NORMAL' ) then
             call line_str(6,2,'Individual high level wavefunction',   &
                           lin,':',trim(schm%wfn),lfin)
           end if
           write(*,*)
           call line_dp(6,2,'Molecular mass (g/mol)',lin,':',          &
                        'F12.5',mol(imol)%mass,lfin)
           call line_int(6,2,'Number of atoms',lin,':','I3',           &
                         mol(imol)%nat,lfin)
           call line_str(6,2,'Component phase',lin,':',                &
                         trim(mol(imol)%phase),lfin)
           call line_ivec(6,2,'Indistinguishable fragments',lin,':',   &
                          mol(imol)%nfrag,1,'I3',mol(imol)%frag,lfin)
           call line_int(6,2,'Fragments permutations',lin,':','I3',    &
                         mol(imol)%npermu,lfin)
           call line_int(6,2,'Number of conformers',lin,':','I3',      &
                         mol(imol)%nconf,lfin)
           write(*,*)
         end if
!
         do iconf = 1, mol(imol)%nconf
!
           if ( (trim(fcalc).eq.'ALL') .or.                            &
                                          (trim(fcalc).eq.'FREQ') ) then
             mol(imol)%conf(iconf)%freq = mol(imol)%conf(iconf)%freq   &
                                                                   *fact
           end if
!
           if ( trim(fcalc) .eq. 'ALL' ) then
             mol(imol)%conf(iconf)%qtrans = qtrans(temp,volu,          &
                                                   mol(imol)%mass)
             mol(imol)%conf(iconf)%qrot   = qrot(temp,mol(imol)%       &
                                                 conf(iconf)%moment)/  &
                                                 mol(imol)%conf(iconf)%&
                                                 symnum
             mol(imol)%conf(iconf)%qvib   = qvib(temp,mol(imol)%       &
                                                 conf(iconf)%dof,      &
                                                 mol(imol)%conf(iconf)%&
                                                 freq)
           end if
!
           if ( fpermu ) then
             mol(imol)%conf(iconf)%nequi = mol(imol)%npermu
!~              mol(imol)%conf(iconf)%nequi = 2*mol(imol)%frag(1)*mol(imol)%npermu ! Rotamers
           end if
!
           if ( fenan ) then
             if ( chkenan(mol(imol)%nat,mol(imol)%conf(iconf)%coord) ) &
                                                                    then
               mol(imol)%conf(iconf)%chiral = .TRUE.             
               mol(imol)%conf(iconf)%nequi  = 2*mol(imol)%conf(iconf)% &
                                                                   nequi
             end if
           end if
!
           if ( debug ) then
             call print_titleint(6,5,'Information of conformer',1,     &
                                 'I3',iconf,'.') 
             write(*,*)
             call line_str(6,2,'Input file name',lin,':',              &
                           trim(mol(imol)%conf(iconf)%inp),lfin)
!
             if ( trim(mol(imol)%schm%fschm) .eq. 'HL' ) then
               call line_str(6,2,'Auxiliary input file name',lin,':',  &
                             trim(mol(imol)%conf(iconf)%auxinp(1)),lfin)
             else if ( trim(mol(imol)%schm%fschm) .eq. 'LLSOL' ) then
               call line_str(6,2,'Auxiliary input file name (sol)',    &
                     lin,':',trim(mol(imol)%conf(iconf)%auxinp(1)),lfin)
               call line_str(6,2,'Auxiliary input file name  (gas)',   &
                     lin,':',trim(mol(imol)%conf(iconf)%auxinp(2)),lfin)
             end if
!
             write(*,*)
             call line_dvec(6,2,'Principal moments (au)',lin,':',3,1,  &
                            'F11.5',mol(imol)%conf(iconf)%moment(:),   &
                            lfin)
             call line_int(6,2,'Symmetry number',lin,':','I3',         &
                           mol(imol)%conf(iconf)%symnum,lfin)
             call line_log(6,2,'The structure has an enantiomer',lin,  &
                           ':',mol(imol)%conf(iconf)%chiral,lfin)
             call line_int(6,2,'PES degeneracy',lin,':','I3',          &
                           mol(imol)%conf(iconf)%nequi,lfin)
             call line_dp(6,2,'SCF energy (au)',lin,':','F16.9',       &
                          mol(imol)%conf(iconf)%Escf,lfin)
!
             if ( trim(fcalc) .eq. 'ALL' ) then
               call line_dp(6,2,'Translational partition function',lin,&
                            ':','D16.10',mol(imol)%conf(iconf)%qtrans, &
                            lfin)
               call line_dp(6,2,'Rotational partition function',lin,   &
                            ':','D16.10',mol(imol)%conf(iconf)%qrot,   &
                            lfin)
               call line_dp(6,2,'Vibrational partition function',lin,  &
                            ':','D16.10',mol(imol)%conf(iconf)%qvib,   &
                            lfin)
               call line_dp(6,2,'Electronic partition function',lin,   &
                            ':','F16.4',mol(imol)%conf(iconf)%qel,lfin)
             end if
!
             write(*,*)
           end if

!
         end do 
!
       end do
!
! Computing selected thermodynamic properties for each molecule
!
       if ( trim(fcalc) .eq. 'ALL' ) then
         if ( debug ) then
           write(*,'(4X,65("*"))')
           write(*,'(4X,10("*"),X,("Individual thermodynamic prope'//  &
                                           'rties section"),X,10("*"))')
           write(*,'(4X,65("*"))')
           write(*,*)
         end if
!
         call print_thermo(temp,nmol,mol,ffree,fentha,fentro,cutoff,   &
                           fsoln,dsolv,msolv,debug)
       else if ( trim(fcalc) .ne. 'FREQ' ) then
         do imol = 1, nmol
           do iconf = 1, mol(imol)%nconf
             mol(imol)%conf(iconf)%Escf = mol(imol)%conf(iconf)%Escf   &
                                                             *au2kJ*1000    
           end do
         end do       
       end if
!
! Computing equilibrium properties 
!
       if ( trim(fcalc) .ne. 'FREQ' ) then
         do imol = 1, nmol
           if ( mol(imol)%nconf .gt. 1 ) then
             doconf = .TRUE.
             exit
           end if
         end do
!
         if ( nreac .gt. 0 ) doequi = .TRUE.
!
         if ( doequi .or. doconf ) then
           write(*,'(4X,66("*"))')
           write(*,'(4X,10("*"),X,("Equilibrium thermodynamic prop'//  &
                                          'erties section"),X,10("*"))')
           write(*,'(4X,66("*"))')
           write(*,*)
!
           if ( doconf ) then
             select case ( forder )
               case('GORDER')
                 write(*,'(2X,A)') 'Printing relative values to th'//  &
                                        'e lowest free energy conformer'
                 write(*,*)
               case('EORDER')
                 write(*,'(2X,A)') 'Printing relative values to th'//  &
                                             'e lowest energy conformer'
                 write(*,*)
               case('DORDER')
                 write(*,'(2X,A)') 'Printing relative values to th'//  &
                                       'e lowest energy conformer at 0K'
                 write(*,*)
             end select
           end if
         end if
!
!~          if ( doconf ) call conformations(temp,nmol,mol,mconf,order,   & ! FLAG: create order array first
!~                                           forder,fcalc,fscreen,frota,  &         otherwise if doconf = FALSE
!~                                           rmsdmax,maemax,baemax,outp,  &         the program might break
!~                                           debug)
         call conformations(temp,nmol,mol,mconf,order,forder,fcalc,    &
                            fscreen,frota,rmsdmax,maemax,baemax,outp,  &
                            debug)
! 
         if ( doequi ) call keq(temp,volu,nmol,mol,mconf,order,nreac,  &
                                reac,fcalc,debug)  
!
         if ( doir ) then
           if ( iline .eq. 1 ) then
             call irspec(outp,nmol,mol,hwhm,iline,gauss)
           else if ( iline .eq. 2 ) then
             call irspec(outp,nmol,mol,hwhm,iline,lor)
           end if
         end if
       end if
!
! Deallocating memory
!
       deallocate(mol)
       deallocate(order)
!
       if ( nreac .ne. 0) deallocate(reac)
!
! Printing timings
!
       call system_clock(t2)
!
       tcpu = dble(t2-t1)/dble(count_rate)
!
       write(*,'(1X,90("="))')
       write(*,*)
       write(*,'(1X,A,3(X,I2,X,A))') 'Total CPU time        ',         &
                                      int(tcpu/(60*60)),'hours',       &
                                      mod(int(tcpu/60),60),'minutes',  &
                                      mod(int(tcpu),60),'seconds'  
       write(*,*)
!
! Printing finishing date 
!    
       call print_end()
!
       end program g16_computherm
!
!======================================================================!
!
       subroutine command_line(inp,outp,temp,pres,volu,debug)
!
       use parameters
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(out)  ::  inp     !  Input file name
       character(len=lenout),intent(out)  ::  outp    !  Output file name
       real(kind=8),intent(out)           ::  temp    !  Temperature (K)
       real(kind=8),intent(out)           ::  pres    !  Pressure (atm)
       real(kind=8),intent(out)           ::  volu    !  Volume (m**3)
       logical,intent(out)                ::  debug   !  Debug mode
!
! Local variables
!
       character(len=lencmd)              ::  cmd     !  Command executed
       character(len=lenarg)              ::  code    !  Executable name
       character(len=lenarg)              ::  arg     !  Argument read
       character(len=lenarg)              ::  next    !  Next argument to be read
       integer                            ::  io      !  Status
       integer                            ::  i       !  Index
!
! Setting defaults
!
       inp   = 'thermo.inp'
       outp  = 'thermo.out'
       temp  =  298.15d0
       pres  =  1.0d0
       volu  =  Ratm*temp/1000.0d0/pres!1.0d0/1000.0d0!
       debug = .FALSE.
!
! Reading command line
!
       call get_command_argument(0,code)
       call get_command(cmd)
! Reading command line options
       if ( command_argument_count().eq.0) return
!
       i = 1
       do
         call get_command_argument(i,arg)
         if ( len_trim(arg) == 0 ) exit
         i = i+1
         select case ( arg )
           case ('-f','-file','--file')
             call get_command_argument(i,inp,status=io)
             call check_arg(inp,io,arg,cmd)
             i = i + 1
           case ('-o','-out','-outp','-output','--output')
             call get_command_argument(i,outp,status=io)
             call check_arg(outp,io,arg,cmd)
             i = i + 1
           case ('-t','-T','-temp','--temp','--temperature')
             call get_command_argument(i,next,status=io)
             read(next,*) temp
             i = i + 1
!
             volu = Ratm*temp/1000.0d0/pres
!
           case ('-p','-P','-pres','--pres','--pressure')
             call get_command_argument(i,next,status=io)
             read(next,*) pres
             i = i + 1
!
             volu = Ratm*temp/1000.0d0/pres
!
           case ('-v','--debug','--verbose')
             debug = .TRUE.
             i = i + 1
           case ('-h','-help','--help')
             call print_help()
             call print_end()
           case default
             write(*,*)
             write(*,'(2X,68("="))')
             write(*,'(3X,A)')    'ERROR:  Unknown statements from'//  &
                                  ' command line'
             write(*,*)
             write(*,'(4X,A)')     trim(cmd)
             write(*,*)
             write(*,'(3X,2(A))') 'Unrecognised command-line option'// &
                                  '  :  ', arg
             write(*,'(3X,A)')    'Please, to request help execute'
             write(*,*)
             write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
             write(*,'(2X,68("="))')
             write(*,*)
             call print_end()
         end select
       end do
!
       return
       end subroutine command_line
!
!======================================================================!
!
       subroutine print_help()
!
       implicit none
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'Command-line options:'
       write(*,*)
       write(*,'(5X,A)') '-h,--help             Print usage inform'//  &
                          'tion and exit'
       write(*,'(5X,A)') '-f,--file             Input file name(s)'
       write(*,'(5X,A)') '-o,--output           Output file name'
       write(*,'(5X,A)') '-T,--temperature      Absolute temperatu'//  &
                         're (K)'
       write(*,'(2X,68("="))')
       write(*,*)
!
       return
       end subroutine print_help
!
!======================================================================!
!
! CHKENAN - CHecK ENANtiomer
!
! This function returns
!
       logical function chkenan(nat,coord)
!
       use geometry
       use superpose
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  coord  !  Input coordinates
       integer,intent(in)                        ::  nat    !  Number of particles
!
! Declaration of the local variables
!
       real(kind=8),dimension(3,nat)             ::  dmat   !  Douple precision auxliary matrix
       real(kind=8)                              ::  rmsd   !  Root-mean-square deviation
!
! Checking if the input molecule is chiral
!
       chkenan = .FALSE.
!
! Generating the mirror image of the target molecule
!
       call Sxy(nat,coord,dmat) 
!
! Checking if the two structures are equivalent
!
       call dsuperpose('QUATKK',3,nat,dmat,coord)
!
       call dmcalcrmsd(3,nat,dmat,coord,rmsd)
!
       if ( rmsd .gt. 1.0E-4 ) chkenan = .TRUE.
!
       return
       end function chkenan
!
!======================================================================!
!
       subroutine print_thermo(temp,nmol,mol,ffree,fentha,fentro,      &
                               cutoff,fsoln,dsolv,msolv,debug)
!
       use utils,      only:  print_titleint
       use parameters
       use datatypes
       use statmech
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(nmol),intent(inout)  ::  mol     !  Molecules information
       real(kind=8),intent(in)                       ::  temp    !  Temperature
       real(kind=8),intent(in)                       ::  cutoff  !  Cutoff frequency
       real(kind=8),intent(in)                       ::  dsolv   !
       real(kind=8),intent(in)                       ::  msolv   !
       character(len=8),intent(in)                   ::  ffree   !  Free energy calculation flag
       character(len=8),intent(in)                   ::  fentha  !  Enthalpy calculation flag
       character(len=8),intent(in)                   ::  fentro  !  Entropy calculation flag
       integer,intent(in)                            ::  nmol    !  Number of chemical species
       logical,intent(in)                            ::  fsoln   !  Standard state flag
       logical,intent(in)                            ::  debug   !  Debug mode
!
! Local variables
!
       character(len=64)                             ::  fmt1    !  Format variable
       character(len=64)                             ::  fmt2    !  Format variable
       integer                                       ::  imol    !  Molecule index
       integer                                       ::  iconf   !  Conformer index
! 
! Computing equilibrium thermodynamic properties 
!
       fmt1 = '(1X,A,3(1X,F12.4),2(1X,F16.4))'
!
       fmt2 = '(1X,A,2X,3(1X,F12.8),4X,F12.8,4X,F16.8,1X,F16.8)'
!
       do imol = 1, nmol
         if ( debug ) then
           call print_titleint(6,1,'Starting molecular system',1,      &
                             'I3',imol,'=')
           write(*,*)
         end if
!
         do iconf = 1, mol(imol)%nconf
           if ( debug ) then
             call print_titleint(6,3,'Thermodynamic properties of '//  &
                                 'conformer',1,'I3',iconf,'-')            
             write(*,*)
           end if
!
           mol(imol)%conf(iconf)%Escf   = mol(imol)%conf(iconf)%Escf   &
                                                             *au2kJ*1000
!
           mol(imol)%conf(iconf)%ZPVE   = ZPVE(mol(imol)%conf(iconf)%  &
                                               dof,mol(imol)%          &
                                               conf(iconf)%freq)
           mol(imol)%conf(iconf)%D0     = mol(imol)%conf(iconf)%ZPVE   &
                                        + mol(imol)%conf(iconf)%Escf
! Computing translational and rotational contributions to enthalpy
           mol(imol)%conf(iconf)%Htrans = Htrans(temp)
           mol(imol)%conf(iconf)%Hrot   = Hrot(temp,mol(imol)%         &
                                                conf(iconf)%rotdof)
! Computing translational and rotational contributions to energy
           mol(imol)%conf(iconf)%Etrans = Etrans(temp)
           mol(imol)%conf(iconf)%Erot   = mol(imol)%conf(iconf)%Hrot     
! Computing translational and rotational contributions to entropy
           mol(imol)%conf(iconf)%Strans = Strans(mol(imol)%            &
                                                 conf(iconf)%qtrans)
           mol(imol)%conf(iconf)%Srot   = Srot(mol(imol)%conf(iconf)%  &
                                               qrot,                   &
                                               mol(imol)%conf(iconf)%  &
                                               rotdof)
          mol(imol)%conf(iconf)%Sel     = Sel(mol(imol)%conf(iconf)%qel)
! Computing translational and rotational contributions to free energy
           mol(imol)%conf(iconf)%Gtrans = Gfree(temp,mol(imol)%        &
                                                conf(iconf)%qtrans)
           mol(imol)%conf(iconf)%Grot   = Gfree(temp,mol(imol)%        &
                                                conf(iconf)%qrot)
! Computing vibrational contribution to entropy
           select case ( fentha )
             case ('RRHO')
               mol(imol)%conf(iconf)%Hvib =                            &
                                 Hvib(temp,mol(imol)%conf(iconf)%dof,  &
                                      mol(imol)%conf(iconf)%freq)
             case ('QHO') 
               mol(imol)%conf(iconf)%Hvib =                            &
                                 Hqho(temp,mol(imol)%conf(iconf)%dof,  &
                                      mol(imol)%conf(iconf)%freq,cutoff)
             case ('QRRHO') 
               mol(imol)%conf(iconf)%Hvib =                            &
                               Hqrrho(temp,mol(imol)%conf(iconf)%dof,  &
                                      mol(imol)%conf(iconf)%freq,cutoff)
               mol(imol)%conf(iconf)%ZPVE =                            &
                             ZPVEdamp(mol(imol)%conf(iconf)%dof,       &
                                      mol(imol)%conf(iconf)%freq,cutoff)
           end select
! Computing vibrational contribution to entalphy
           select case ( fentro )
             case ('RRHO')
               mol(imol)%conf(iconf)%Svib =                            &
                                 Svib(temp,mol(imol)%conf(iconf)%dof,  &
                                      mol(imol)%conf(iconf)%freq)
             case ('QHO') 
               mol(imol)%conf(iconf)%Svib =                            &
                                 Sqho(temp,mol(imol)%conf(iconf)%dof,  &
                                      mol(imol)%conf(iconf)%freq,cutoff)
             case ('QRRHO')
               mol(imol)%conf(iconf)%Svib =                            &
                            Sqrrho(temp,mol(imol)%conf(iconf)%dof,     &
                                   mol(imol)%conf(iconf)%freq,cutoff,  &
                                   mol(imol)%conf(iconf)%moment)
           end select
! Computing vibrational contribution to free energy
           select case ( ffree )
             case ('RRHO')
               mol(imol)%conf(iconf)%Gvib =                            &
                                  Gfree(temp,mol(imol)%conf(iconf)%qvib)
             case ('QHO') 
               mol(imol)%conf(iconf)%Gvib =                            &
                      Gfree(temp,qqho(temp,mol(imol)%conf(iconf)%dof,  &
                            mol(imol)%conf(iconf)%freq,cutoff))
             case ('QRRHO')
               mol(imol)%conf(iconf)%Gvib =                            &
                    Gfree(temp,qqrrho(temp,mol(imol)%conf(iconf)%dof,  &
                          mol(imol)%conf(iconf)%freq,cutoff,           &
                          mol(imol)%conf(iconf)%moment))
             case ('MIX')
               mol(imol)%conf(iconf)%Gvib = mol(imol)%conf(iconf)%Hvib &
                                       - temp*mol(imol)%conf(iconf)%Svib
           end select
! Correcting with the zero-point vibrational energy
           mol(imol)%conf(iconf)%Hvib = mol(imol)%conf(iconf)%ZPVE     &
                                            + mol(imol)%conf(iconf)%Hvib 
           mol(imol)%conf(iconf)%Evib = mol(imol)%conf(iconf)%Hvib   
           mol(imol)%conf(iconf)%Gvib = mol(imol)%conf(iconf)%ZPVE     &
                                            + mol(imol)%conf(iconf)%Gvib
! Computing thermal corrections
           mol(imol)%conf(iconf)%H = mol(imol)%conf(iconf)%Htrans      &
                                        + mol(imol)%conf(iconf)%Hrot   &
                                        + mol(imol)%conf(iconf)%Hvib
           mol(imol)%conf(iconf)%E = mol(imol)%conf(iconf)%Etrans      &
                                        + mol(imol)%conf(iconf)%Erot   &
                                        + mol(imol)%conf(iconf)%Evib 
           mol(imol)%conf(iconf)%S = mol(imol)%conf(iconf)%Strans      &
                                        + mol(imol)%conf(iconf)%Srot   &
                                        + mol(imol)%conf(iconf)%Svib   &
                                        + mol(imol)%conf(iconf)%Sel
           mol(imol)%conf(iconf)%G = mol(imol)%conf(iconf)%Gtrans      &
                                        + mol(imol)%conf(iconf)%Grot   &
                                        + mol(imol)%conf(iconf)%Gvib
!
           if ( debug ) then
             write(*,'(5X,A)') '- SI units (kJ/mol and J/K*mol)'
             write(*,*)
             write(*,'(8X,A,3(2X,A),2(6X,A))')                         &
                     'Translation','   Rotation','  Vibration',        &
                     ' Electronic','    Total  '
             write(*,fmt1) 'De ',                                      &
                            zero,zero,zero,                            &
                            mol(imol)%conf(iconf)%Escf/1000,           &
                            mol(imol)%conf(iconf)%Escf/1000
             write(*,fmt1) 'ZPE',                                      &
                            zero,zero,mol(imol)%conf(iconf)%ZPVE/1000, &
                            zero,mol(imol)%conf(iconf)%ZPVE/1000
             write(*,fmt1) 'D0 ',                                      &
                            zero,zero,mol(imol)%conf(iconf)%ZPVE/1000, &
                            mol(imol)%conf(iconf)%Escf/1000,           &                          
                            mol(imol)%conf(iconf)%D0/1000
             write(*,fmt1) 'E  ',                                      &
                            mol(imol)%conf(iconf)%Etrans/1000,         &
                            mol(imol)%conf(iconf)%Erot/1000,           &
                            mol(imol)%conf(iconf)%Evib/1000,           &
                            mol(imol)%conf(iconf)%Escf/1000,           &
                            mol(imol)%conf(iconf)%E/1000               &
                          + mol(imol)%conf(iconf)%Escf/1000
             write(*,fmt1) 'H  ',                                      &
                            mol(imol)%conf(iconf)%Htrans/1000,         &
                            mol(imol)%conf(iconf)%Hrot/1000,           &
                            mol(imol)%conf(iconf)%Hvib/1000,           &
                            mol(imol)%conf(iconf)%Escf/1000,           &
                            mol(imol)%conf(iconf)%H/1000               &
                          + mol(imol)%conf(iconf)%Escf/1000
             write(*,fmt1) 'T*S',                                      &
                            mol(imol)%conf(iconf)%Strans*temp,         &
                            mol(imol)%conf(iconf)%Srot*temp,           &
                            mol(imol)%conf(iconf)%Svib*temp,           &
                            mol(imol)%conf(iconf)%Sel*temp,            &
                            mol(imol)%conf(iconf)%S*temp  
             write(*,fmt1) 'G  ',                                      &
                            mol(imol)%conf(iconf)%Gtrans/1000,         &
                            mol(imol)%conf(iconf)%Grot/1000,           &
                            mol(imol)%conf(iconf)%Gvib/1000,           &
                            mol(imol)%conf(iconf)%Escf/1000,           &
                            mol(imol)%conf(iconf)%G/1000               &
                          + mol(imol)%conf(iconf)%Escf/1000
             write(*,*)
!
             write(*,'(5X,A)') '- Atomic units (Hartree/particle a'//  &
                                                'nd Hartree/K*particle)'
             write(*,*)
             write(*,'(8X,A,2(2X,A),3X,A,2X,A,5X,A)') 'Translation',   &
                   '   Rotation','  Vibration','Thermal correction',   &
                   ' Electronic','      Total'
             write(*,fmt2) 'De ',                                      &
                            zero,zero,zero,zero,                       &
                            mol(imol)%conf(iconf)%Escf/au2kj/1000,     &
                            mol(imol)%conf(iconf)%Escf/au2kj/1000
             write(*,fmt2) 'ZPE',                                      &
                            zero,zero,                                 &
                            mol(imol)%conf(iconf)%ZPVE/au2kj/1000,     &
                            mol(imol)%conf(iconf)%ZPVE/au2kj/1000,     &
                            zero,                                      &
                            mol(imol)%conf(iconf)%ZPVE/au2kj/1000
             write(*,fmt2) 'D0 ',                                      &
                            zero,zero,                                 &
                            mol(imol)%conf(iconf)%ZPVE/au2kj/1000,     &
                            mol(imol)%conf(iconf)%ZPVE/au2kj/1000,     &
                            mol(imol)%conf(iconf)%Escf/au2kj/1000,     &                          
                            mol(imol)%conf(iconf)%D0/au2kj/1000
             write(*,fmt2) 'E  ',                                      &
                            mol(imol)%conf(iconf)%Etrans/au2kj/1000,   &
                            mol(imol)%conf(iconf)%Erot/au2kj/1000,     &
                            mol(imol)%conf(iconf)%Evib/au2kj/1000,     &
                            mol(imol)%conf(iconf)%E/au2kj/1000,        &
                            mol(imol)%conf(iconf)%Escf/au2kj/1000,     &
                            mol(imol)%conf(iconf)%E/au2kj/1000         &
                          + mol(imol)%conf(iconf)%Escf/au2kj/1000
             write(*,fmt2) 'H  ',                                      &
                            mol(imol)%conf(iconf)%Htrans/au2kj/1000,   &
                            mol(imol)%conf(iconf)%Hrot/au2kj/1000,     &
                            mol(imol)%conf(iconf)%Hvib/au2kj/1000,     &
                            mol(imol)%conf(iconf)%H/au2kj/1000,        &
                            mol(imol)%conf(iconf)%Escf/au2kj/1000,     &
                            mol(imol)%conf(iconf)%H/au2kj/1000         &
                          + mol(imol)%conf(iconf)%Escf/au2kj/1000
             write(*,fmt2) 'T*S',                                      &
                         mol(imol)%conf(iconf)%Strans/au2kj/1000*temp, &
                         mol(imol)%conf(iconf)%Srot/au2kj/1000*temp,   &
                         mol(imol)%conf(iconf)%Svib/au2kj/1000*temp,   &
                         mol(imol)%conf(iconf)%S/au2kj/1000*temp,      &
                         mol(imol)%conf(iconf)%Sel/au2kj/1000*temp,    &
                         mol(imol)%conf(iconf)%S/au2kj/1000*temp
             write(*,fmt2) 'G  ',                                      &
                            mol(imol)%conf(iconf)%Gtrans/au2kj/1000,   &
                            mol(imol)%conf(iconf)%Grot/au2kj/1000,     &
                            mol(imol)%conf(iconf)%Gvib/au2kj/1000,     &
                            mol(imol)%conf(iconf)%G/au2kj/1000,        &
                            mol(imol)%conf(iconf)%Escf/au2kj/1000,     &
                            mol(imol)%conf(iconf)%G/au2kj/1000         &
                          + mol(imol)%conf(iconf)%Escf/au2kj/1000         
             write(*,*)
           end if
!
           mol(imol)%conf(iconf)%E = mol(imol)%conf(iconf)%Escf        &
                                               + mol(imol)%conf(iconf)%E 
! 
           mol(imol)%conf(iconf)%H = mol(imol)%conf(iconf)%Escf        &
                                               + mol(imol)%conf(iconf)%H 
!
           mol(imol)%conf(iconf)%G = mol(imol)%conf(iconf)%Escf        &
                                               + mol(imol)%conf(iconf)%G
!
           if ( fsoln ) then
             mol(imol)%conf(iconf)%G = mol(imol)%conf(iconf)%G         &
                                    + Rjul*temp*dlog(msolv/1000.0/dsolv)
             mol(imol)%conf(iconf)%S = mol(imol)%conf(iconf)%S         &
                                         - Rjul*dlog(msolv/1000.0/dsolv)
           end if
!
         end do
       end do
!
       return
       end subroutine print_thermo
!
!======================================================================!
!
       subroutine conformations(temp,nmol,mol,mconf,order,forder,fcalc,&
                                fscreen,frota,rmsdmax,maemax,baemax,   &
                                outp,debug)
!
       use parameters
       use datatypes
       use sorting,    only: dviqsort
       use screening
       use utils,      only: lenout,                                   &
                             print_title,                              &
                             print_titleint,                           &
                             print_end
!
       implicit none
!
! Input/output variables
!
!
       type(molecule),dimension(nmol),intent(inout)  ::  mol      !  Molecules information
       character(len=lenout)                         ::  outp     !  Output file name
       real(kind=8),intent(in)                       ::  temp     !  Temperature
       real(kind=8),intent(in)                       ::  rmsdmax  !  Maximum value for RMSD
       real(kind=8),intent(in)                       ::  maemax   !  Maximum value for MAE
       real(kind=8),intent(in)                       ::  baemax   !  Maximum value for BAE
       integer,dimension(nmol,mconf)                 ::  order    !  Energy-based order
       integer,intent(in)                            ::  nmol     !  Number of chemical species
       integer,intent(in)                            ::  mconf    !  Maximum number of conformers
       character(len=8),intent(in)                   ::  forder   !  Conformations order flag
       character(len=8),intent(in)                   ::  fcalc    !  Calculation information flag
       character(len=8),intent(in)                   ::  frota    !  Rotation method flag
       logical,intent(in)                            ::  fscreen  !  Screening calculation flag
       logical,intent(in)                            ::  debug    !  Debug mode
!
! Local variables
!
       character(len=64)                             ::  fmt1     !  I/O format 
       character(len=64)                             ::  fmt2     !  I/O format 
       character(len=64)                             ::  fmt3     !  I/O format 
       real(kind=8),dimension(mconf)                 ::  Gval     !  Free energy values
       real(kind=8),dimension(mconf)                 ::  pop      !
       real(kind=8)                                  ::  popsum   !
       real(kind=8)                                  ::  daux     !
       integer,dimension(mconf)                      ::  idx      !  Non-redundant indexes
       integer                                       ::  imol     !  Molecule index
       integer                                       ::  iconf    !  Conformer index
       integer                                       ::  isize    !  Largest input file length
!
! Computing thermodynamic properties of the conformational equilibrium
!
       call print_title(6,1,'Analysis of conformations','=')
       write(*,*)
!
       isize = 0
       do imol = 1, nmol
         do iconf = 1, mol(imol)%nconf
           isize = max(isize,len_trim(mol(imol)%conf(iconf)%inp))
         end do
       end do
!
       isize = max(isize,10)
!
       write(fmt1,*) isize
       fmt1 = adjustl(fmt1)
!
       write(fmt2,*) isize - 10 + 1
       fmt2 = adjustl(fmt2)
!
       if ( trim(fcalc) .eq. 'ALL' ) then
         fmt1 = '(4X,I4,4X,A'//trim(fmt1)//',2X,F10.4,1X,F10.4,3X,'//  &
                                               'F10.4,1X,F10.4,4X,F10.6)'
         fmt2 = '(2(1X,A10),'//trim(fmt2)//'X,2(1X,A10,3X,A10),1X,A)'
       else
         fmt1 = '(4X,I4,4X,A'//trim(fmt1)//',1X,3(1X,F10.4))'
         fmt2 = '(2(1X,A10),'//trim(fmt2)//'X,2(1X,A10),1X,A)'
       end if
!
       do imol = 1, nmol
         call print_titleint(6,3,'Conformations of molecular system',1,&
                             'I3',imol,'-')
         write(*,*)
!
         if ( mol(imol)%nconf .eq. 1 ) then
           write(*,'(1X,A,1X,I3,1X,A)') 'Molecular system',imol,       &
                                         'is formed by a single species'
           write(*,'(1X,A)') 'Skipping analysis of conformations'
           write(*,*) 
!
           order(imol,1)    = 1
           mol(imol)%pop(1) = 1.0d0/mol(imol)%conf(1)%nequi
           mol(imol)%popsum = mol(imol)%conf(1)%nequi
         else
!
           if ( fscreen ) then
!
             if ( mol(imol)%nfrag .gt. 1 ) then
               write(*,'(2X,68("="))')
               write(*,'(3X,A)')    'ERROR:  Requested calculation'//  &
                                                  ' not implemented yet'
               write(*,*)
               write(*,'(3X,A)')    'Screening of the conformation'//  &
                                                   's is only available'
               write(*,'(3X,A)')    ' for aggregates with simply 1'//  &
                                                        ' fragment type'
               write(*,'(2X,68("="))')
               write(*,*)
               call print_end()
             end if
!
             call screen(outp,mol(imol),mol(imol)%frag(1),frota,       &
                         rmsdmax,maemax,baemax,idx,debug)
!
           else
             do iconf = 1, mol(imol)%nconf
               idx(iconf) = iconf
             end do
           end if
!
           Gval(:)  = 0.0d0
!
           select case ( forder )
             case('GORDER')
               do iconf = 1, mol(imol)%nconf
                 order(imol,iconf) = idx(iconf)
                 Gval(iconf)       = mol(imol)%conf(idx(iconf))%G
               end do
             case('EORDER')
               do iconf = 1, mol(imol)%nconf
                 order(imol,iconf) = idx(iconf)
                 Gval(iconf)       = mol(imol)%conf(idx(iconf))%Escf
               end do
             case('DORDER')
               do iconf = 1, mol(imol)%nconf
                 order(imol,iconf) = idx(iconf)
                 Gval(iconf)       = mol(imol)%conf(idx(iconf))%D0
               end do
           end select
!
           call dviqsort(mconf,Gval,order(imol,:),1,mol(imol)%nconf)
!
           write(*,'(1X,A,1X,I4)') 'Printing relative values to co'//  &
                                                 'nformer',order(imol,1)
           write(*,*)
!
           if ( (trim(fcalc).eq.'ALL') .and.                           &
                                       (trim(forder).ne.'GORDER') ) then
!
             mol(imol)%popsum = 0.0d0
!
             do iconf = 1, mol(imol)%nconf
               mol(imol)%popsum = mol(imol)%popsum +                   & 
                          mol(imol)%conf(order(imol,iconf))%nequi      &  
                         *dexp((mol(imol)%conf(order(imol,1))%G        &
                               - mol(imol)%conf(order(imol,iconf))%G)  &
                                                             /Rjul/temp)   
             end do
!
             do iconf = 1, mol(imol)%nconf
               mol(imol)%pop(order(imol,iconf)) =                      &
                          dexp((mol(imol)%conf(order(imol,1))%G        &
                               - mol(imol)%conf(order(imol,iconf))%G)  &
                                            /Rjul/temp)/mol(imol)%popsum
             end do
!
           else 
!
             mol(imol)%popsum = 0.0d0
!
             do iconf = 1, mol(imol)%nconf
               mol(imol)%popsum = mol(imol)%popsum +                   & 
                              mol(imol)%conf(order(imol,iconf))%nequi  &  
                             *dexp((Gval(1)-Gval(iconf))/Rjul/temp)   
             end do
!
             do iconf = 1, mol(imol)%nconf
               mol(imol)%pop(order(imol,iconf)) =                      &
                                dexp((Gval(1)-Gval(iconf))/Rjul/temp)  &
                                                       /mol(imol)%popsum          
             end do
! 
           end if
!
           if ( trim(fcalc) .eq. 'ALL' ) then
             popsum = 0.0d0
!
             do iconf = 1, mol(imol)%nconf
               popsum = popsum +                                       & 
                      mol(imol)%conf(order(imol,iconf))%nequi          &  
                     *dexp((mol(imol)%conf(order(imol,1))%Escf         &
                            - mol(imol)%conf(order(imol,iconf))%Escf)  &
                                                             /Rjul/temp)   
             end do
!
             do iconf = 1, mol(imol)%nconf
               pop(order(imol,iconf)) =                                &
                      dexp((mol(imol)%conf(order(imol,1))%Escf         &
                            - mol(imol)%conf(order(imol,iconf))%Escf)  &
                                                      /Rjul/temp)/popsum
             end do           
           end if
!
           call print_title(6,5,'Populations of the conformers','.')
           write(*,*)
           if ( trim(fcalc) .eq. 'ALL' ) then
             write(*,fmt2) 'Conformer','File name','Energy',           &
                             'Population','G','Population','Free energy'
           else
             write(*,fmt2) 'Conformer','File name','Energy','Boltz',   &
                                                            'Population'                                   
           end if
!
           if ( trim(fcalc) .eq. 'ALL' ) then
             daux = 0.0d0
             do iconf = 1, mol(imol)%nconf     
               daux = daux + mol(imol)%conf(order(imol,iconf))%nequi   &
                           *exp(-(mol(imol)%conf(order(imol,iconf))%G  &
                           - mol(imol)%conf(order(imol,1))%G)/Rjul/temp)
               write(*,fmt1) order(imol,iconf),                        &
                        mol(imol)%conf(order(imol,iconf))%inp(:isize), &
                         (mol(imol)%conf(order(imol,iconf))%Escf       &
                          - mol(imol)%conf(order(imol,1))%Escf)/1000,  &
                         mol(imol)%conf(order(imol,iconf))%nequi       &
                                            *pop(order(imol,iconf)),   &
                         (mol(imol)%conf(order(imol,iconf))%G          &
                          - mol(imol)%conf(order(imol,1))%G)/1000,     &
                         mol(imol)%conf(order(imol,iconf))%nequi       &
                                  *mol(imol)%pop(order(imol,iconf)),   &
!~                        (mol(imol)%conf(order(imol,1))%G/1000/au2kj -   & 
!~                                                     Kb*temp*dlog(daux))
                                               Rjul*temp*dlog(daux)/1000
!~ write(*,*) mol(imol)%conf(order(imol,1))%G,Rjul*temp*dlog(daux)/1000
             end do
           else
             do iconf = 1, mol(imol)%nconf   
               write(*,fmt1) order(imol,iconf),                        &
                        mol(imol)%conf(order(imol,iconf))%inp(:isize), &
                        (mol(imol)%conf(order(imol,iconf))%Escf        &
                         - mol(imol)%conf(order(imol,1))%Escf)/1000,   &
                        dexp((Gval(1)-Gval(iconf))/Rjul/temp),         &
                        mol(imol)%conf(order(imol,iconf))%nequi        &
                                       *mol(imol)%pop(order(imol,iconf))
             end do
           end if
!
           write(*,*)
!
           if ( (trim(fcalc).eq.'ALL') .and. debug ) then
             call print_title(6,5,'Interconversion equilibrium pro'//  &
                                                          'perties','.')
             write(*,*)
             write(*,'(8(1X,A10),3X,A)') 'Conformer','D0','De','E',    &
                                                     'H','T*S','G','Keq'
!
             fmt3 = '(1X,I7,3X,6(1X,F10.4),1X,D10.4)'
!
             do iconf = 1, mol(imol)%nconf     
               write(*,fmt3) order(imol,iconf),                        &
                          (mol(imol)%conf(order(imol,iconf))%Escf      &
                           - mol(imol)%conf(order(imol,1))%Escf)/1000, &
                          (mol(imol)%conf(order(imol,iconf))%D0        &
                           - mol(imol)%conf(order(imol,1))%D0)/1000,   &
                          (mol(imol)%conf(order(imol,iconf))%E         &
                           - mol(imol)%conf(order(imol,1))%E)/1000,    &
                          (mol(imol)%conf(order(imol,iconf))%H         &
                           - mol(imol)%conf(order(imol,1))%H)/1000,    &
                          temp*(mol(imol)%conf(order(imol,iconf))%S    &
                           - mol(imol)%conf(order(imol,1))%S)/1000,    &
                          (mol(imol)%conf(order(imol,iconf))%G         &
                           - mol(imol)%conf(order(imol,1))%G)/1000,    &
                          dexp((mol(imol)%conf(order(imol,1))%G        &
                               - mol(imol)%conf(order(imol,iconf))%G)  &
                                                             /Rjul/temp)
             end do
             write(*,*) 
           end if
!
! 
         end if
!
       end do
!
       return
       end subroutine conformations
!
!======================================================================!
!
! References
! ----------
!
! - Jensen, J. H., "Predicting accurate absolute binding energies in 
!    aqueous solution: thermodynamic considerations for electronic 
!    structure methods", Phys. Chem. Chem. Phys., The Royal Society of
!    Chemistry, 2015, 17, 12441. http://dx.doi.org/10.1039/C5CP00628G
!
       subroutine keq(temp,volu,nmol,mol,mconf,order,nreac,reac,fcalc, &
                      debug)
!
       use parameters
       use datatypes
       use utils,      only: lenline,                                  &
                             print_title,                              &
                             print_titleint
       use statmech
!
       implicit none
!
! Input/output variables
!
       type(reaction),dimension(nreac),intent(inout)  ::  reac    !  Reactions information
       type(molecule),dimension(nmol),intent(inout)   ::  mol     !  Molecules information
       real(kind=8),intent(in)                        ::  temp    !  Temperature (K)
       real(kind=8),intent(in)                        ::  volu    !  Volume (L)
       integer,dimension(nmol,mconf),intent(in)       ::  order   !  Energy-based order
       integer,intent(in)                             ::  nmol    !  Number of chemical species
       integer,intent(in)                             ::  nreac   !  Number of reactions
       integer,intent(in)                             ::  mconf   !  Maximum number of conformers
       character(len=8),intent(in)                    ::  fcalc   !  Calculation information flag
       logical,intent(in)                             ::  debug   !  Debug mode
!
! Local variables
!
       character(len=lenline)                         ::  line    !
       integer,dimension(nmol)                        ::  inu     !
       integer                                        ::  ireac   !  Reaction index
       integer                                        ::  iconf   !  Conformation index
       integer                                        ::  imol    !  Molecule index
       integer                                        ::  iaux    !
!
! Computing thermodynamic properties of the reactions
!
       call print_title(6,1,'Thermochemistry of the reactions','=')
       write(*,*)
!
       if ( trim(fcalc) .eq. 'ALL' ) then
         do imol = 1, nmol
           mol(imol)%Gtot = mol(imol)%conf(order(imol,1))%G -          & 
                                        Rjul*temp*dlog(mol(imol)%popsum)
!
           mol(imol)%Detot = 0.0d0
           mol(imol)%D0tot = 0.0d0
           mol(imol)%Htot  = 0.0d0       
           mol(imol)%Htot  = 0.0d0
           mol(imol)%Stot  = 0.0d0
!
           do iconf = 1, mol(imol)%nconf
             mol(imol)%Detot = mol(imol)%Detot + mol(imol)%pop(iconf)  &
                 *mol(imol)%conf(iconf)%nequi*mol(imol)%conf(iconf)%Escf
             mol(imol)%D0tot = mol(imol)%D0tot + mol(imol)%pop(iconf)  &
                   *mol(imol)%conf(iconf)%nequi*mol(imol)%conf(iconf)%D0
             mol(imol)%Etot = mol(imol)%Etot + mol(imol)%pop(iconf)    &
                    *mol(imol)%conf(iconf)%nequi*mol(imol)%conf(iconf)%E
             mol(imol)%Htot = mol(imol)%Htot + mol(imol)%pop(iconf)    &
                    *mol(imol)%conf(iconf)%nequi*mol(imol)%conf(iconf)%H
             mol(imol)%Stot = mol(imol)%Stot +                         &
                       mol(imol)%conf(iconf)%nequi                     &
                       *(mol(imol)%conf(iconf)%S*mol(imol)%pop(iconf)  &
                         - Rjul*mol(imol)%pop(iconf)                   &
                                            *dlog(mol(imol)%pop(iconf)))
           end do
!
         end do 
       else
         do imol = 1, nmol
           mol(imol)%Detot = 0.0d0
!
           do iconf = 1, mol(imol)%nconf
             mol(imol)%Detot = mol(imol)%Detot + mol(imol)%pop(iconf)  &
                 *mol(imol)%conf(iconf)%nequi*mol(imol)%conf(iconf)%Escf
           end do
         end do 
       end if
!
       do ireac = 1, nreac
!
         call print_titleint(6,3,'Starting reaction number',1,         &
                             'I3',ireac,'-')
         write(*,*)
!
         call print_react(nmol,mol,reac(ireac),line)
!
!~          write(*,'(1X,A)') 'Reaction Title :'
         write(*,'(1X,A)') 'Reaction scheme:'
         write(*,'(3X,A)')  trim(line)
         write(*,*)
! 
         iaux = 1
         do imol = 1, nmol
           if ( reac(ireac)%nu(imol) .lt. 0 ) then
             inu(iaux) = imol
             iaux = iaux + 1
           end if
         end do
!
         do imol = 1, nmol
           if ( reac(ireac)%nu(imol) .gt. 0 ) then
             inu(iaux) = imol
             iaux = iaux + 1
           end if
         end do
!
         do imol = 1, nmol
           if ( reac(ireac)%nu(imol) .eq. 0 ) then
             inu(iaux) = imol
             iaux = iaux + 1
           end if
         end do
!
         call print_title(6,5,'Microscopic equilibrium properties','.')
         write(*,*)
!
         call microkeq(temp,nmol,mol,reac(ireac),inu,mconf,order,      &
                       fcalc,debug)
!
         call print_title(6,5,'Macroscopic equilibrium properties','.')
         write(*,*)
!
         call reac_thermo(temp,reac(ireac),nmol,mol,fcalc)
!
         if ( trim(fcalc) .eq. 'ALL' ) then
           write(*,'(10(1X,A10))') 'De','D0','E','H','T*S','G','Keq'
           write(*,'(6(1X,F10.4),1X,D10.4)') reac(ireac)%dDe/1000,     &
                                             reac(ireac)%dD0/1000,     &
                                             reac(ireac)%dE/1000,      &
                                             reac(ireac)%dH/1000,      &
                                        temp*reac(ireac)%dS/1000,      &
                                             reac(ireac)%dG/1000,      &
                                             reac(ireac)%Keq
         else
           write(*,'(1X,A,1X,F10.4)') 'Dissociation energy =',         &
                                                    reac(ireac)%dDe/1000
         end if
         write(*,*)
!
       end do
!
       return
       end subroutine keq
!
!======================================================================!
!
       subroutine reac_thermo(temp,reac,nmol,mol,fcalc)
!
       use parameters
       use datatypes
!
       implicit none
!
! Input/output variables
!
       type(reaction),intent(inout)               ::  reac   !  Reaction information
       type(molecule),dimension(nmol),intent(in)  ::  mol    !  Molecules information
       character(len=8),intent(in)                ::  fcalc  !  Calculation information flag
       real(kind=8),intent(in)                    ::  temp   !  Temperature
       integer,intent(in)                         ::  nmol   !  Number of chemical species
!
! Local variables
!
       integer                                    ::  imol   !  Molecule index
! 
! Computing thermodynamic properties of reaction
!
       if ( trim(fcalc) .eq. 'ALL' ) then
         reac%dDe = 0.0d0
         reac%dD0 = 0.0d0
         reac%dE  = 0.0d0
         reac%dH  = 0.0d0
         reac%dS  = 0.0d0
         reac%dG  = 0.0d0
!
         do imol = 1, nmol
           reac%dDe = reac%dDe + reac%nu(imol)*mol(imol)%Detot
           reac%dD0 = reac%dD0 + reac%nu(imol)*mol(imol)%D0tot
           reac%dE  = reac%dE  + reac%nu(imol)*mol(imol)%Etot
           reac%dH  = reac%dH  + reac%nu(imol)*mol(imol)%Htot
           reac%dS  = reac%dS  + reac%nu(imol)*mol(imol)%Stot
           reac%dG  = reac%dG  + reac%nu(imol)*mol(imol)%Gtot
         end do
!
         reac%Keq = dexp(-reac%dG/Rjul/temp)
       else
         reac%dDe = 0.0d0
!
         do imol = 1, nmol
           reac%dDe = reac%dDe + reac%nu(imol)*mol(imol)%Detot
         end do
       end if
!
       return
       end subroutine reac_thermo
!
!======================================================================!
!
       subroutine microkeq(temp,nmol,mol,reac,inu,mconf,order,fcalc,   &
                           debug)
!
       use datatypes
       use parameters
!
       implicit none
!
! Input/output variables
!
       type(reaction),intent(inout)               ::  reac    !  Reactions information
       type(molecule),dimension(nmol),intent(in)  ::  mol     !  Molecules information
       real(kind=8),intent(in)                    ::  temp    !  Temperature
       integer,dimension(nmol,mconf),intent(in)   ::  order   !  Energy-based order
       integer,dimension(nmol),intent(in)         ::  inu     !  Sorted stoichiometric coefficients
       integer,intent(in)                         ::  nmol    !  Number of chemical species
       integer,intent(in)                         ::  mconf   !  Maximum number of conformers
       character(len=8),intent(in)                ::  fcalc   !  Calculation information flag
       logical,intent(in)                         ::  debug   !  Debug mode
!
! Local variables
!
       character(len=32)                          ::  straux  !  Auxiliary string
       character(len=64)                          ::  fmt1    !  Format variable
       character(len=64)                          ::  fmt2    !  Format variable
       real(kind=8)                               ::  dDe     !  Potential energy change
       real(kind=8)                               ::  dD0     !  Zero kelvin energy change
       real(kind=8)                               ::  dE      !  Internal energy change
       real(kind=8)                               ::  dH      !  Enthalpy change
       real(kind=8)                               ::  dS      !  Entropu energy change
       real(kind=8)                               ::  dG      !  Free energy change
       real(kind=8)                               ::  Keq     !  Equilibrium constant
       integer,dimension(nmol)                    ::  iconf   !  Equilibrium conformations indexes
       integer                                    ::  imol    !  Molecule index
       integer                                    ::  mmol    !  Number of molecules in the reaction
       integer                                    ::  ncombi  !  Number of combinations
       integer                                    ::  mcombi  !
       integer                                    ::  icombi  !  Combination index
       integer                                    ::  iaux1   !  Auxiliary variable
       integer                                    ::  iaux2   !  Auxiliary variable
       integer                                    ::  isize1  !
       integer                                    ::  isize2  !
!
! Computing the microequilibrium thermodynamic properties
!
       ncombi = 1
       mmol   = 0
       do imol = 1, nmol
         if ( reac%nu(imol) .ne. 0 ) then
           mmol = mmol + 1
           ncombi = ncombi*mol(imol)%nconf
         end if
       end do
!
       iaux1 = len('Conformer index') + 1
       iaux2 = 4*mmol
!
       isize1 = max(iaux1,iaux2)
       isize2 = isize1
!
       isize1 = isize1 - iaux1 + 1
       isize2 = isize2 - iaux2 + 1
! 
       write(fmt1,*) isize1
       fmt1 = adjustl(fmt1)
       fmt1 = '(1X,A,'//trim(fmt1)//'X,8(1X,A10))'
!
       write(fmt2,*) mmol
       fmt2 = adjustl(fmt2)

       fmt2 = '('//trim(fmt2)//'(1X,I3),'
!
       write(straux,*) isize2
       straux = adjustl(straux)
       fmt2 = trim(fmt2)//trim(straux)//'X,6(1X,F10.4),1X,D10.4)'
!
       if ( trim(fcalc) .eq. 'ALL' ) then
         write(*,fmt1) 'Conformer index','De','D0','E','H','T*S',    &
                                                               'G','Keq'
       else
         write(*,fmt1) 'Conformer index','De'
       end if
!
       iconf(:) = 1
!
       mcombi = 1
       if ( debug ) mcombi = ncombi
!
       do icombi = 1, mcombi
!
         if ( trim(fcalc) .eq. 'ALL' ) then
           dDe = 0.0d0
           dD0 = 0.0d0
           dE  = 0.0d0
           dH  = 0.0d0
           dS  = 0.0d0
           dG  = 0.0d0
!
           do imol = 1, nmol
             dDe = dDe + reac%nu(inu(imol))*mol(inu(imol))%            &
                            conf(order(inu(imol),iconf(inu(imol))))%Escf
             dD0 = dD0 + reac%nu(inu(imol))*mol(inu(imol))%            &
                              conf(order(inu(imol),iconf(inu(imol))))%D0
             dE  = dE  + reac%nu(inu(imol))*mol(inu(imol))%            &
                               conf(order(inu(imol),iconf(inu(imol))))%E
             dH  = dH  + reac%nu(inu(imol))*mol(inu(imol))%            &
                               conf(order(inu(imol),iconf(inu(imol))))%H
             dS  = dS  + reac%nu(inu(imol))*mol(inu(imol))%            &
                               conf(order(inu(imol),iconf(inu(imol))))%S
             dG  = dG  + reac%nu(inu(imol))*mol(inu(imol))%            &
                               conf(order(inu(imol),iconf(inu(imol))))%G
           end do
!
           Keq = dexp(-dG/Rjul/temp)
!
           write(*,fmt2) (order(inu(imol),iconf(inu(imol))),           &
                                                          imol=1,mmol),&
                          dDe/1000,dD0/1000,dE/1000,dH/1000,           &
                                                temp*dS/1000,dG/1000,Keq
         else
           dDe = 0.0d0
!
           do imol = 1, nmol
             dDe = dDe + reac%nu(inu(imol))*mol(inu(imol))%            &
                            conf(order(inu(imol),iconf(inu(imol))))%Escf
           end do
!
           write(*,fmt2) (order(inu(imol),iconf(inu(imol))),           &
                                                   imol=1,mmol),dDe/1000
         end if
!
         iconf(inu(mmol)) = iconf(inu(mmol)) + 1
         imol = mmol
!
         do while ( (imol.gt.1) .and.                                  &
                       (iconf(inu(imol)).gt.mol(inu(imol))%nconf) )
           iconf(inu(imol)) = 1
           iconf(inu(imol-1)) = iconf(inu(imol-1)) + 1
           imol = imol - 1
         end do
!
       end do
!
       write(*,*)
!
       end subroutine microkeq
!
!======================================================================!
!
       subroutine print_react(nmol,mol,reac,line)
!
       use parameters
       use datatypes
       use utils,      only: lenline,                                  &
                             print_title,                              &
                             print_titleint
!
       implicit none
!
! Input/output variables
!
!
       type(molecule),dimension(nmol),intent(in)  ::  mol     !  Molecules information
       type(reaction),intent(in)                  ::  reac    !  Reactions information
       character(len=lenline),intent(out)         ::  line    !
       integer,intent(in)                         ::  nmol    !  Number of chemical species
!
! Local variables
!
       character(len=80)                          ::  reacts  !  Reactants
       character(len=80)                          ::  prods   !  Products line
       character(len=16)                          ::  straux  !  
       integer                                    ::  imol    !  Molecule index
       integer                                    ::  ireact  !
       integer                                    ::  iprod   !
!
! Printing reaction scheme
!
       reacts = ''
       prods  = ''
!
       ireact = 0
       iprod  = 0
!
       do imol = 1, nmol
         if ( reac%nu(imol) .lt. 0 ) then
           ireact = imol
!
           write(straux,*) abs(reac%nu(imol))
           straux = adjustl(straux)
!
           reacts = trim(straux)//' '//trim(mol(imol)%molname)
! 
           exit
         end if
       end do
!
       if ( ireact .eq. 0 ) return
!
       do imol = 1, nmol
         if ( reac%nu(imol) .gt. 0 ) then
           iprod = imol
!
           write(straux,*) reac%nu(imol)
           straux = adjustl(straux)
!
           prods  = trim(straux)//' '//trim(mol(imol)%molname)
!
           exit
         end if
       end do
!
       do imol = 1, nmol
!
         write(straux,*) abs(reac%nu(imol))
         straux = adjustl(straux)
!
         if ( (reac%nu(imol).lt.0) .and. (imol.gt.ireact) ) then
           reacts = trim(reacts)//' + '//trim(straux)//' '//           &
                                                 trim(mol(imol)%molname)
         else if ( (reac%nu(imol).gt.0) .and. (imol.gt.iprod) ) then
           prods  = trim(prods)//' + '//trim(straux)//' '//            & 
                                                 trim(mol(imol)%molname)
         end if
!
       end do
!
       line = trim(reacts)//' = '//trim(prods)
!
       return
       end subroutine print_react
!
!======================================================================!
!
       subroutine irspec(outp,nmol,mol,hwhm,iline,fbroad)
!
       use datatypes
       use utils,      only: uniout,lenout,leninp
!
       implicit none
!
! Input/output variables
!
       character(len=lenout),intent(in)           ::  outp    !  Output file name
       type(molecule),dimension(nmol),intent(in)  ::  mol     !  Molecules information
       real(kind=8),intent(inout)                 ::  hwhm    !  
       integer,intent(in)                         ::  nmol    !
       integer,intent(in)                         ::  iline   !
       real(kind=8),external                      ::  fbroad  ! 
!
! Local variables
!
       character(len=lenout)                      ::  file1   !
       character(len=leninp)                      ::  file2   !
       character(len=leninp)                      ::  file3   !
       real(kind=8),dimension(:),allocatable      ::  xpts    !
       real(kind=8),dimension(:),allocatable      ::  ymol    !
       real(kind=8),dimension(:),allocatable      ::  ytot    !
       real(kind=8)                               ::  xin     !
       real(kind=8)                               ::  xfin    !
       real(kind=8)                               ::  x       !
       real(kind=8)                               ::  yconf   !
       real(kind=8)                               ::  dx      !
       integer                                    ::  imol    !
       integer                                    ::  npts    !
       integer                                    ::  iconf   !
       integer                                    ::  i,j     !
!
! Printing IR spectrum
!
       xin  = 400.0d0
       xfin = 4000.0d0
       dx   = 0.1d0
!
       npts = int((xfin - xin)/dx) + 1
       dx = (xfin - xin)/(npts - 1)
!
       allocate(xpts(npts),ymol(npts),ytot(npts))
!
       if ( iline .eq. 1 ) hwhm = hwhm/sqrt(2.0*log(2.0))
!
       file1 = outp(:len_trim(outp)-4)//'.dat'
       open(unit=uniout+1,file=trim(file1),action='write')
!
       ytot(:) = 0.0d0
       do imol = 1, nmol
!
         file2 = outp(:len_trim(outp)-4)//'_'//                        &
                                         trim(mol(imol)%molname)//'.dat'
         open(unit=uniout+2,file=trim(file2),action='write')
!
         ymol(:) = 0.0d0
         do iconf = 1, mol(imol)%nconf
!
           file3 = outp(:len_trim(outp)-4)//'_'//                      &
                      mol(imol)%conf(iconf)%inp                        &
                        (:len_trim(mol(imol)%conf(iconf)%inp)-4)//'.dat'
           open(unit=uniout+3,file=trim(file3),action='write')
!
           x = xin
           do i = 1, npts
             xpts(i) = x
!
             yconf = 0.0
             do  j = 1, mol(imol)%conf(iconf)%dof
               yconf = yconf + mol(imol)%conf(iconf)%inten(j)*x        &
                           *fbroad(x-mol(imol)%conf(iconf)%freq(j),hwhm)
             end do
!
             if ( mol(imol)%readw ) then
               yconf = yconf*mol(imol)%conc*mol(imol)%conf(iconf)%weight
             else
               yconf = yconf*mol(imol)%conc                            &
                                         *mol(imol)%conf(iconf)%nequi  &
                                         *mol(imol)%pop(iconf)
             end if
!
             write(uniout+3,*) xpts(i),yconf
!
             ymol(i) = ymol(i) + yconf
             ytot(i) = ytot(i) + yconf
!
             x = x + dx
           end do
!
           close(uniout+3)
!
         end do
!
         do i = 1, npts
           write(uniout+2,*) xpts(i),ymol(i)
         end do
!
         close(uniout+2)
!
       end do
!
       do i = 1, npts
         write(uniout+1,*) xpts(i),ytot(i)
       end do
!
       close(uniout+1)
!
       deallocate(xpts,ymol,ytot)
!
       return
       end subroutine irspec
!
!======================================================================!
