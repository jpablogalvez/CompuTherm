!======================================================================!
!
       program g16_computherm
!
       use parameters
       use datatypes
       use utils,      only: leninp,lenout,                            &
                             line_str,                                 &
                             line_int,                                 &
                             line_dp,                                  &
                             line_dvec,                                &
                             line_ivec,                                &
                             line_log,                                 &
                             print_title,                              &
                             print_titleint,                           &
                             print_start,                              &
                             print_end
       use input
       use g16files
       use statmech
       use mathtools
!
       implicit none
!
       include 'timings.h'
!
       type(molecule),dimension(:),allocatable  ::  mol     !  Molecules information
       type(reaction),dimension(:),allocatable  ::  reac    !  Reactions information
       character(len=leninp)                    ::  inp     !  Input file name
       character(len=lenout)                    ::  outp    !  Output file name
       character(len=8)                         ::  fqvib   !  Qvib calculation flag
       character(len=8)                         ::  ffree   !  Free energy calculation flag
       character(len=64)                        ::  fmt1    !  Format variable
       real(kind=8)                             ::  temp    !  Temperature (K)
       real(kind=8)                             ::  pres    !  Pressure (atm)
       real(kind=8)                             ::  volu    !  Volume (m**3)
       real(kind=8)                             ::  thr     !  Threshold frequency
       real(kind=8)                             ::  fact    !  Frequencies scaling factor
       integer,dimension(:,:),allocatable       ::  order   !  Energy-based order
       integer                                  ::  nmol    !  Number of chemical species
       integer                                  ::  imol    !  Molecule index
       integer                                  ::  mconf   !  Maximum number of conformers
       integer                                  ::  iconf   !  Conformer index
       integer                                  ::  ifrag   !  Fragment index
       integer                                  ::  nreac   !  Number of reactions
       integer                                  ::  lfin    !  
       integer                                  ::  lin     ! 
       logical                                  ::  fsoln   !  Standard state flag
       logical                                  ::  fenan   !  Enantiomers calculation flag
       logical                                  ::  fpermu  !  Permutations calculation flag
       logical                                  ::  debug   !  Debug mode
!
! Declaration of external functions
!
       logical                                  ::  chkenan ! 

!
! Printing header
!
       write(*,*)
       write(*,'(5X,82("#"))')
       write(*,'(5X,10("#"),3X,A,3X,13("#"))') 'CompuTherm - Comput'//  &
                                    'ational Thermochemistry Calculator'
       write(*,'(5X,82("#"))')   
       write(*,*)
       write(*,'(9X,A)') 'Welcome to CompuTherm, a very simple pr'//  &
                                          'ogram for the calculation of'
       write(*,'(9X,A)') ' thermodynamic properties from electroni'//  &
                                              'c structure calculations'
       write(*,*)
       write(*,'(1X,90("-"))')
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
       lin  = 35
       lfin = 80
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
       call read_inp(inp,nmol,mol,thr,fact,fqvib,ffree,                &
                     fenan,fpermu,fsoln,nreac,reac)
!
       mconf = 1
       do imol = 1, nmol
         mconf = max(mconf,mol(imol)%nconf)
       end do
!
       allocate(order(nmol,mconf))
       order(:,:) = -1
!
! Printing summary of the general input file information
!
       call print_title(6,1,'Files information','=')
       write(*,*)
!       
       call line_str(6,2,'General input file name',lin,':',            &
                     trim(inp),lfin)
       call line_str(6,2,'Outputnput file name',lin,':',               &
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
       write(*,'(2X,A)') 'Vibrational partition function'
       call line_str(6,3,'Type of partition function',lin,':',         &
                     trim(fqvib),lfin)
       call line_dp(6,3,'Frequencies scaling factor',lin,':',          &
                    'F6.4',fact,lfin)
       write(*,*)
       write(*,'(2X,A)') 'Gibbs free energy'
       call line_str(6,3,'Type of free energy calculation',lin,':',    &
                     trim(ffree),lfin)
       call line_log(6,3,'Enantiomers calculaton',lin,':',fenan,lfin)
       call line_log(6,3,'Indistinguishable aggregates',lin, &
                     ':',fpermu,lfin)
       call line_log(6,3,'Standard state 1M',lin,':',fsoln,lfin)
       write(*,*)
       write(*,'(2X,A)') 'Standard state'
       call line_dp(6,3,'Temperature (K)',lin,':','F12.2',temp,lfin)
       call line_dp(6,3,'Pressure (atm)',lin,':','F12.4',pres,lfin)
       call line_dp(6,3,'Volume (L)',lin,':','F12.4',volu*1000.0d0,lfin)
       write(*,*) 
!
! Processing Gaussian16 input files
!
       call read_g16(nmol,mol)
!
! Printing summary of the Gaussian16 information
!
       call print_title(6,1,'Gaussian16 input files information','=')
       write(*,*)
!
       do imol = 1, nmol
         if ( fpermu ) then
           mol(imol)%npermu = 1
           do ifrag = 1, mol(imol)%nfrag
             mol(imol)%npermu = mol(imol)%npermu                       &
                                       *factorial(mol(imol)%frag(ifrag))
           end do
         end if
!
         call print_titleint(6,3,'Information of molecular system',1,  &
                             'I3',imol,'-')
         write(*,*)
         call line_str(6,2,'Molecule name',lin,':',                    &
                       trim(mol(imol)%molname),lfin)
         call line_dp(6,2,'Molecular mass (g/mol)',lin,':','F12.5',    &
                      mol(imol)%mass,lfin)
         call line_int(6,2,'Number of atoms',lin,':','I3',             &
                       mol(imol)%nat,lfin)
         call line_str(6,2,'Component phase',lin,':',                  &
                       trim(mol(imol)%phase),lfin)
         call line_ivec(6,2,'Indistinguishable fragments',lin,':',     &
                        mol(imol)%nfrag,1,'I3',mol(imol)%frag,lfin)
         call line_int(6,2,'Fragments permutations',lin,':','I3',      &
                       mol(imol)%npermu,lfin)
         call line_int(6,2,'Number of conformers',lin,':','I3',        &
                       mol(imol)%nconf,lfin)
         write(*,*)
!
         do iconf = 1, mol(imol)%nconf
!
           mol(imol)%conf(iconf)%freq   = mol(imol)%conf(iconf)%freq   &
                                                                   *fact
!
           mol(imol)%conf(iconf)%qtrans = qtrans(temp,volu,mol(imol)%  &
                                                 mass)
           mol(imol)%conf(iconf)%qrot   = qrot(temp,mol(imol)%         &
                                               conf(iconf)%moment)/    &
                                               mol(imol)%conf(iconf)%  &
                                               symnum
           mol(imol)%conf(iconf)%qvib   = qvib(temp,mol(imol)%         &
                                               conf(iconf)%dof,        &
                                               mol(imol)%conf(iconf)%  &
                                               freq)
!
           if ( fpermu ) then
             mol(imol)%conf(iconf)%nequi = mol(imol)%conf(iconf)%nequi &
                                                       *mol(imol)%npermu
           end if
!
           if ( fenan ) then
             if ( chkenan(mol(imol)%nat,mol(imol)%conf(iconf)%coord) ) &
                                                                    then
               mol(imol)%conf(iconf)%nequi = 2*mol(imol)%conf(iconf)%  &
                                                                   nequi
               mol(imol)%conf(iconf)%chiral = .TRUE.             
             end if
           end if
!
           call print_titleint(6,5,'Information of conformer',1,       &
                               'I3',iconf,'.') 
           write(*,*)
           call line_str(6,2,'Input file name',lin,':',                &
                         trim(mol(imol)%conf(iconf)%inp),lfin)

           call line_dvec(6,2,'Principal moments (au)',lin,':',3,1,    &
                          'F11.5',mol(imol)%conf(iconf)%moment(:),lfin)
           call line_int(6,2,'Symmetry number',lin,':','I3',           &
                         mol(imol)%conf(iconf)%symnum,lfin)
           call line_log(6,2,'The structure has an enantiomer',lin,    &
                         ':',mol(imol)%conf(iconf)%chiral,lfin)
           call line_int(6,2,'PES degeneracy',lin,':','I3',            &
                         mol(imol)%conf(iconf)%nequi,lfin)
           call line_dp(6,2,'SCF energy (au)',lin,':','F16.9',         &
                        mol(imol)%conf(iconf)%Escf,lfin)
           call line_dp(6,2,'Translational partition function',lin,    &
                        ':','F16.4',mol(imol)%conf(iconf)%qtrans,lfin)
           call line_dp(6,2,'Rotational partition function',lin,       &
                        ':','F16.4',mol(imol)%conf(iconf)%qrot,lfin)
           call line_dp(6,2,'Vibrational partition function',lin,      &
                        ':','F16.4',mol(imol)%conf(iconf)%qvib,lfin)
           call line_dp(6,2,'Electronic partition function',lin,       &
                        ':','F16.4',mol(imol)%conf(iconf)%qel,lfin)
           write(*,*)
!
!~            if ( fpermu ) then
!~              mol(imol)%conf(iconf)%nequi = mol(imol)%conf(iconf)%nequi &
!~                                            *mol(imol)%conf(iconf)%symnum
!~              mol(imol)%conf(iconf)%qtrans = mol(imol)%conf(iconf)%nequi&
!~                                            *mol(imol)%conf(iconf)%qtrans
!~            end if
         end do 
!
       end do
!
! Computing selected thermodynamic properties for each molecule
!
       write(*,'(4X,65("*"))')
       write(*,'(4X,10("*"),X,("Individual thermodynamic propertie'//  &
                                               's section"),X,10("*"))')
       write(*,'(4X,65("*"))')
       write(*,*)
!
       do imol = 1, nmol
         call print_titleint(6,1,'Starting molecular system',1,        &
                             'I3',imol,'=')
         write(*,*)
!
         do iconf = 1, mol(imol)%nconf
           call print_titleint(6,3,'Thermodynamic properties of co'//  &
                               'nformer',1,'I3',iconf,'-')            
           write(*,*)
!
           mol(imol)%conf(iconf)%Escf   = mol(imol)%conf(iconf)%Escf   &
                                                             *au2kJ*1000
!
           mol(imol)%conf(iconf)%ZPVE   = ZPVE(mol(imol)%conf(iconf)%  &
                                               dof,mol(imol)%          &
                                               conf(iconf)%freq)
           mol(imol)%conf(iconf)%D0     = mol(imol)%conf(iconf)%ZPVE   &
                                        + mol(imol)%conf(iconf)%Escf
!
           mol(imol)%conf(iconf)%Htrans = Htrans(temp)
           mol(imol)%conf(iconf)%Hrot   = Hrot(temp,mol(imol)%         &
                                                conf(iconf)%rotdof)
           mol(imol)%conf(iconf)%Hvib   = mol(imol)%conf(iconf)%ZPVE   &
                                        + Hvib(temp,mol(imol)%         &
                                               conf(iconf)%dof,        &
                                               mol(imol)%conf(iconf)%  &
                                               freq)
!
           mol(imol)%conf(iconf)%Htherm = mol(imol)%conf(iconf)%Htrans &
                                        + mol(imol)%conf(iconf)%Hrot   &
                                        + mol(imol)%conf(iconf)%Hvib
!
           if ( fsoln ) mol(imol)%conf(iconf)%Htherm =                 &
                                          mol(imol)%conf(iconf)%Htherm! &
!~                                         + Rjul*temp*dlog(volu*1000)
!
           mol(imol)%conf(iconf)%Etrans = Etrans(temp)
           mol(imol)%conf(iconf)%Erot   = mol(imol)%conf(iconf)%Hrot
           mol(imol)%conf(iconf)%Evib   = mol(imol)%conf(iconf)%Hvib
!
           mol(imol)%conf(iconf)%Etherm = mol(imol)%conf(iconf)%Etrans &
                                        + mol(imol)%conf(iconf)%Erot   &
                                        + mol(imol)%conf(iconf)%Evib      
!
           mol(imol)%conf(iconf)%Strans = Strans(mol(imol)%            &
                                                 conf(iconf)%qtrans)
           mol(imol)%conf(iconf)%Srot   = Srot(temp,mol(imol)%         &
                                               conf(iconf)%qrot,       &
                                               mol(imol)%conf(iconf)%  &
                                               Hrot)
           mol(imol)%conf(iconf)%Svib   = Svib(temp,mol(imol)%         &
                                               conf(iconf)%dof,        &
                                               mol(imol)%conf(iconf)%  &
                                               freq)
          mol(imol)%conf(iconf)%Sel     = Sel(mol(imol)%conf(iconf)%qel)
!
           mol(imol)%conf(iconf)%Stherm = mol(imol)%conf(iconf)%Strans &
                                        + mol(imol)%conf(iconf)%Srot   &
                                        + mol(imol)%conf(iconf)%Svib   &
                                        + mol(imol)%conf(iconf)%Sel
!
           mol(imol)%conf(iconf)%Gtrans = Gfree(temp,mol(imol)%        &
                                                conf(iconf)%qtrans)
           mol(imol)%conf(iconf)%Grot   = Gfree(temp,mol(imol)%        &
                                                conf(iconf)%qrot)
           mol(imol)%conf(iconf)%Gvib   = mol(imol)%conf(iconf)%ZPVE   &
                                + Gfree(temp,mol(imol)%conf(iconf)%qvib)
!
           mol(imol)%conf(iconf)%Gtherm = mol(imol)%conf(iconf)%Gtrans &
                                        + mol(imol)%conf(iconf)%Grot   &
                                        + mol(imol)%conf(iconf)%Gvib
!
           if ( fsoln ) mol(imol)%conf(iconf)%Gtherm =                 &
                                          mol(imol)%conf(iconf)%Gtherm! &
!~                                         + Rjul*temp*dlog(volu*1000)
!
           fmt1 = '(1X,A,3(1X,F12.4),2(1X,F16.4))'
!
           write(*,'(5X,A)') '- SI units (kJ/mol and J/K*mol)'
           write(*,*)
           write(*,'(8X,A,3(2X,A),2(6X,A))')                           &
                   'Translation','   Rotation','  Vibration',          &
                   ' Electronic','    Total  '
           write(*,fmt1) 'De ',                                        &
                          zero,zero,zero,                              &
                          mol(imol)%conf(iconf)%Escf/1000,             &
                          mol(imol)%conf(iconf)%Escf/1000
           write(*,fmt1) 'ZPE',                                        &
                          zero,zero,mol(imol)%conf(iconf)%ZPVE/1000,   &
                          zero,mol(imol)%conf(iconf)%ZPVE/1000
           write(*,fmt1) 'D0 ',                                        &
                          zero,zero,mol(imol)%conf(iconf)%ZPVE/1000,   &
                          mol(imol)%conf(iconf)%Escf/1000,             &                          
                          mol(imol)%conf(iconf)%D0/1000
           write(*,fmt1) 'E  ',                                        &
                          mol(imol)%conf(iconf)%Etrans/1000,           &
                          mol(imol)%conf(iconf)%Erot/1000,             &
                          mol(imol)%conf(iconf)%Evib/1000,             &
                          mol(imol)%conf(iconf)%Escf/1000,             &
                          mol(imol)%conf(iconf)%Etherm/1000            &
                          + mol(imol)%conf(iconf)%Escf/1000
           write(*,fmt1) 'H  ',                                        &
                          mol(imol)%conf(iconf)%Htrans/1000,           &
                          mol(imol)%conf(iconf)%Hrot/1000,             &
                          mol(imol)%conf(iconf)%Hvib/1000,             &
                          mol(imol)%conf(iconf)%Escf/1000,             &
                          mol(imol)%conf(iconf)%Htherm/1000            &
                          + mol(imol)%conf(iconf)%Escf/1000
           write(*,fmt1) 'S  ',                                        &
                          mol(imol)%conf(iconf)%Strans,                &
                          mol(imol)%conf(iconf)%Srot,                  &
                          mol(imol)%conf(iconf)%Svib,                  &
                          mol(imol)%conf(iconf)%Sel,                   &
                          mol(imol)%conf(iconf)%Stherm
           write(*,fmt1) 'G  ',                                        &
                          mol(imol)%conf(iconf)%Gtrans/1000,           &
                          mol(imol)%conf(iconf)%Grot/1000,             &
                          mol(imol)%conf(iconf)%Gvib/1000,             &
                          mol(imol)%conf(iconf)%Escf/1000,             &
                          mol(imol)%conf(iconf)%Gtherm/1000            &
                          + mol(imol)%conf(iconf)%Escf/1000
           write(*,*)
!
           fmt1 = '(1X,A,2X,3(1X,F12.8),4X,F12.8,4X,F16.8,1X,F16.8)'
!
           write(*,'(5X,A)') '- Atomic units (Hartree/particle and'//  &
                                                  ' Hartree/K*particle)'
           write(*,*)
           write(*,'(8X,A,2(2X,A),3X,A,2X,A,5X,A)')                    &
                        'Translation','   Rotation','  Vibration',     &
                        'Thermal correction',' Electronic','      Total'
           write(*,fmt1) 'De ',                                        &
                          zero,zero,zero,zero,                         &
                          mol(imol)%conf(iconf)%Escf/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Escf/au2kj/1000
           write(*,fmt1) 'ZPE',                                        &
                          zero,zero,                                   &
                          mol(imol)%conf(iconf)%ZPVE/au2kj/1000,       &
                          mol(imol)%conf(iconf)%ZPVE/au2kj/1000,       &
                          zero,                                        &
                          mol(imol)%conf(iconf)%ZPVE/au2kj/1000
           write(*,fmt1) 'D0 ',                                        &
                          zero,zero,                                   &
                          mol(imol)%conf(iconf)%ZPVE/au2kj/1000,       &
                          mol(imol)%conf(iconf)%ZPVE/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Escf/au2kj/1000,       &                          
                          mol(imol)%conf(iconf)%D0/au2kj/1000
           write(*,fmt1) 'E  ',                                        &
                          mol(imol)%conf(iconf)%Etrans/au2kj/1000,     &
                          mol(imol)%conf(iconf)%Erot/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Evib/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Etherm/au2kj/1000,     &
                          mol(imol)%conf(iconf)%Escf/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Etherm/au2kj/1000      &
                          + mol(imol)%conf(iconf)%Escf/au2kj/1000
           write(*,fmt1) 'H  ',                                        &
                          mol(imol)%conf(iconf)%Htrans/au2kj/1000,     &
                          mol(imol)%conf(iconf)%Hrot/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Hvib/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Htherm/au2kj/1000,     &
                          mol(imol)%conf(iconf)%Escf/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Htherm/au2kj/1000      &
                          + mol(imol)%conf(iconf)%Escf/au2kj/1000
           write(*,fmt1) 'S  ',                                        &
                          mol(imol)%conf(iconf)%Strans/au2kj/1000,     &
                          mol(imol)%conf(iconf)%Srot/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Svib/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Stherm/au2kj/1000,     &
                          mol(imol)%conf(iconf)%Sel/au2kj/1000,        &
                          mol(imol)%conf(iconf)%Stherm/au2kj/1000
           write(*,fmt1) 'G  ',                                        &
                          mol(imol)%conf(iconf)%Gtrans/au2kj/1000,     &
                          mol(imol)%conf(iconf)%Grot/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Gvib/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Gtherm/au2kj/1000,     &
                          mol(imol)%conf(iconf)%Escf/au2kj/1000,       &
                          mol(imol)%conf(iconf)%Gtherm/au2kj/1000      &
                          + mol(imol)%conf(iconf)%Escf/au2kj/1000         
           write(*,*)
!
           mol(imol)%conf(iconf)%Etherm = mol(imol)%conf(iconf)%Escf   &
                                        + mol(imol)%conf(iconf)%Etherm 
! 
           mol(imol)%conf(iconf)%Htherm = mol(imol)%conf(iconf)%Escf   &
                                        + mol(imol)%conf(iconf)%Htherm 
!
           mol(imol)%conf(iconf)%Gtherm = mol(imol)%conf(iconf)%Escf   &
                                        + mol(imol)%conf(iconf)%Gtherm
!
         end do
       end do
!
! Computing equilibrium properties
!
       write(*,'(4X,66("*"))')
       write(*,'(4X,10("*"),X,("Equilibrium thermodynamic properti'//  &
                                              'es section"),X,10("*"))')
       write(*,'(4X,66("*"))')
       write(*,*)
!
       call conformations(temp,nmol,mol,mconf,order)
!
       if ( nreac .gt. 0 ) then
         call keq(temp,volu,nmol,mol,mconf,order,nreac,reac,ffree,fsoln)
       end if
!
! Deallocating memory
!
       deallocate(mol)
       deallocate(order)
       if ( nreac .ne. 0) deallocate(reac)
!
! Printing timings
!
       call system_clock(t2)
!
       tcpu = dble(t2-t1)/dble(count_rate)
!
       write(*,'(1X,90("-"))')
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
! Generating the mirror image of the molecule
!
       call Sxy(nat,coord,dmat) 
!
! Checking if the two molecules are the same or different
!
       call dsuperpose('QUATKK',3,nat,dmat,coord)
!
       call dmcalcrmsd(3,nat,dmat,coord,rmsd)
!
       if ( rmsd .gt. 0.0d0 ) chkenan = .TRUE.
!
       return
       end function chkenan
!
!======================================================================!
!
       subroutine conformations(temp,nmol,mol,mconf,order)
!
       use parameters
       use datatypes
       use sorting,    only: dviqsort
       use utils,      only: print_title,                              &
                             print_titleint
!
       implicit none
!
! Input/output variables
!
!
       type(molecule),dimension(nmol),intent(inout)  ::  mol     !  Molecules information
       real(kind=8),intent(in)                       ::  temp    !  Temperature
       integer,dimension(nmol,mconf)                 ::  order   !  Energy-based order
       integer,intent(in)                            ::  nmol    !  Number of chemical species
       integer,intent(in)                            ::  mconf   !  Maximum number of conformers
!
! Local variables
!
       character(len=64)                             ::  fmt1    !  I/O format 
       real(kind=8),dimension(mconf)                 ::  Gval    !  Free energy values
       integer                                       ::  imol    !  Molecule index
       integer                                       ::  iconf   !  Conformer index
!
!
!
       call print_title(6,1,'Analysis of conformations','=')
       write(*,*)
       write(*,'(2X,A)') 'Printing relative values to the lowest f'//  &
                                                  'ree energy conformer'
       write(*,*)
!
       do imol = 1, nmol
         call print_titleint(6,3,'Conformations of molecular system',1,&
                             'I3',imol,'-')
         write(*,*)
!
         if ( mol(imol)%nconf .eq. 1 ) then
           write(*,'(1X,A,1X,I3,1X,A)') 'Molecular system',imol,       &
                                         'is formed by a single species'
           write(*,'(1X,A)') 'Skipping analysis of populations'
           write(*,*) 
!
           order(imol,1)    = 1
           mol(imol)%pop(1) = 1.0d0
           mol(imol)%popsum = 1.0d0
         else
!
           Gval(:)  = 0.0d0
!
           do iconf = 1, mol(imol)%nconf
             order(imol,iconf) = iconf
             Gval(iconf)       = mol(imol)%conf(iconf)%Gtherm
           end do
!
           call dviqsort(mconf,Gval,order(imol,:),1,mol(imol)%nconf)
!
           write(*,'(1X,A,1X,I3)') 'Printing relative values to co'//  &
                                                 'nformer',order(imol,1)
           write(*,*)
!
           write(*,'(15X,7(1X,A10),3X,A)') 'D0','De','E','H','T*S',    &
                                           'G','Keq','Population'
!
           fmt1 = '(1X,A,1X,I3,1X,8(1X,F10.4))'
!
           mol(imol)%popsum = 0.0d0
           do iconf = 1, mol(imol)%nconf
             mol(imol)%popsum = mol(imol)%popsum +                     &
                              mol(imol)%conf(order(imol,iconf))%nequi  &
                              *dexp(-(Gval(iconf)-Gval(1))/Rjul/temp)
           end do
!
           do iconf = 1, mol(imol)%nconf
             mol(imol)%pop(order(imol,iconf)) =                        &
                              mol(imol)%conf(order(imol,iconf))%nequi  & 
                              *dexp(-(Gval(iconf)-Gval(1))/Rjul/temp)  &
                               /mol(imol)%popsum
           end do
!
           do iconf = 1, mol(imol)%nconf     
             write(*,fmt1) 'Conformer',order(imol,iconf),              &
                        (mol(imol)%conf(order(imol,iconf))%Escf        &
                         - mol(imol)%conf(order(imol,1))%Escf)/1000,   &
                        (mol(imol)%conf(order(imol,iconf))%D0          &
                         - mol(imol)%conf(order(imol,1))%D0)/1000,     &
                        (mol(imol)%conf(order(imol,iconf))%Etherm      &
                         - mol(imol)%conf(order(imol,1))%Etherm)/1000, &
                        (mol(imol)%conf(order(imol,iconf))%Htherm      &
                         - mol(imol)%conf(order(imol,1))%Htherm)/1000, &
                        temp*(mol(imol)%conf(order(imol,iconf))%Stherm &
                         - mol(imol)%conf(order(imol,1))%Stherm)/1000, &
                        (mol(imol)%conf(order(imol,iconf))%Gtherm      &
                         - mol(imol)%conf(order(imol,1))%Gtherm)/1000, &
                        dexp(-(Gval(iconf)-Gval(1))/Rjul/temp),        &
                        mol(imol)%pop(order(imol,iconf))
           end do
!
           write(*,*) 
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
       subroutine keq(temp,volu,nmol,mol,mconf,order,nreac,reac,       &
                      ffree,fsoln)
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
       character(len=8),intent(in)                    ::  ffree   !  Free energy calculation flag
       real(kind=8),intent(in)                        ::  temp    !  Temperature (K)
       real(kind=8),intent(in)                        ::  volu    !  Volume (m**3)
       integer,dimension(nmol,mconf),intent(in)       ::  order   !  Energy-based order
       integer,intent(in)                             ::  nmol    !  Number of chemical species
       integer,intent(in)                             ::  nreac   !  Number of reactions
       integer,intent(in)                             ::  mconf   !  Maximum number of conformers
       logical,intent(in)                             ::  fsoln   !  Standard state flag
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
!
!
       call print_title(6,1,'Thermochemistry of the reactions','=')
       write(*,*)
!
       select case ( ffree )
         case ('BOLTZ') 
           do imol = 1, nmol
             mol(imol)%Gtot = mol(imol)%conf(1)%Gtherm -               &
                                        Rjul*temp*dlog(mol(imol)%popsum)
!
             if ( fsoln ) mol(imol)%Gtot = mol(imol)%Gtot +            &
                                               Rjul*temp*dlog(volu*1000)
!
             mol(imol)%Detot = 0.0d0
             mol(imol)%D0tot = 0.0d0
             mol(imol)%Htot  = 0.0d0
             mol(imol)%Htot  = 0.0d0
             mol(imol)%Stot  = 0.0d0
!
             do iconf = 1, mol(imol)%nconf
               mol(imol)%Detot = mol(imol)%Detot +                     &
                     mol(imol)%conf(iconf)%Escf*mol(imol)%pop(iconf)
               mol(imol)%D0tot = mol(imol)%D0tot +                     &
                     mol(imol)%conf(iconf)%D0*mol(imol)%pop(iconf)
               mol(imol)%Etot = mol(imol)%Etot +                       &
                     mol(imol)%conf(iconf)%Etherm*mol(imol)%pop(iconf)
               mol(imol)%Htot = mol(imol)%Htot +                       &
                     mol(imol)%conf(iconf)%Htherm*mol(imol)%pop(iconf)
               mol(imol)%Stot = mol(imol)%Stot +                       &
                     mol(imol)%conf(iconf)%Stherm*mol(imol)%pop(iconf) &
                       - Rjul*mol(imol)%pop(iconf)                     &
                                             *dlog(mol(imol)%pop(iconf))
             end do
!
             if ( fsoln ) then
               mol(imol)%Htot = mol(imol)%Htot +                       &
                                               Rjul*temp*dlog(volu*1000)
               mol(imol)%Etot = mol(imol)%Etot +                       &
                                               Rjul*temp*dlog(volu*1000)
             end if
!
           end do 
         case ('QUAD') 
           do imol = 1, nmol
!
             mol(imol)%qtot = 0.0d0
             do iconf = 1, mol(imol)%nconf
               mol(imol)%qtot = mol(imol)%qtot +                       &
                                 mol(imol)%conf(iconf)%qrot            &
                                *mol(imol)%conf(iconf)%qvib            &
                                *mol(imol)%conf(iconf)%qel             &
                                *dexp((mol(imol)%conf(iconf)%Escf      &
                                  -mol(imol)%conf(order(imol,1))%Escf) &
                                                             /Rjul/temp)
                                    
             end do
!
             mol(imol)%qtot = mol(imol)%qtot*mol(imol)%conf(1)%qtrans  &
                     *mol(imol)%conf(order(imol,1))%qel                &
                     *dexp(mol(imol)%conf(order(imol,1))%Escf/Rjul/temp)   
             mol(imol)%Gtot = Gfree(temp,mol(imol)%qtot)
!
           end do
           
       end select
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
         call microkeq(temp,nmol,mol,reac(ireac),inu,mconf,order)
!
         call print_title(6,5,'Macroscopic equilibrium properties','.')
         write(*,*)
!
         write(*,'(10(1X,A10))') 'De','D0','E','H','T*S','G','Keq'
!
         select case ( ffree )
           case ('INDEP')
!
             call print_thermo(temp,reac(ireac))
!
           case ('BOLTZ','QUAD') 
!
             call reac_thermo(temp,reac(ireac),nmol,mol)
!
             call print_thermo(temp,reac(ireac))
!
         end select
!
       end do
!
       return
       end subroutine keq
!
!======================================================================!
!
       subroutine reac_thermo(temp,reac,nmol,mol)
!
       use parameters
       use datatypes
!
       implicit none
!
! Input/output variables
!
       type(reaction),intent(inout)               ::  reac  !  Reaction information
       type(molecule),dimension(nmol),intent(in)  ::  mol   !  Molecules information
       real(kind=8),intent(in)                    ::  temp  !  Temperature
       integer,intent(in)                         ::  nmol  !  Number of chemical species
!
! Local variables
!
       integer                                    ::  imol  !  Molecule index
! 
! Computing thermodynamic properties of reaction
!
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
!
       return
       end subroutine reac_thermo
!
!======================================================================!
!
       subroutine microkeq(temp,nmol,mol,reac,inu,mconf,order)
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
       integer                                    ::  icombi  !  Combination index
       integer                                    ::  iaux1   !  Auxiliary variable
       integer                                    ::  iaux2   !  Auxiliary variable
       integer                                    ::  isize1  !
       integer                                    ::  isize2  !
!
! Computing the microequilibrium thermodynamic properties
!
       reac%dDe = 0.0d0
       reac%dD0 = 0.0d0
       reac%dE  = 0.0d0
       reac%dH  = 0.0d0
       reac%dS  = 0.0d0
       reac%dG  = 0.0d0
       reac%Keq = 0.0d0
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
!~        fmt2 = trim(fmt2)//trim(straux)//'X,6(1X,F10.4),1X,E10.4)'
       fmt2 = trim(fmt2)//trim(straux)//'X,6(1X,F10.4),1X,D10.4)'

!
!~        write(*,'(1X,A,1X,7(1X,A10))') 'Conformer index','D0','De',   &
!~                                                   'E','H','T*S','G','Keq'
       write(*,fmt1) 'Conformer index','De','D0','E','H','T*S','G','Keq'
!
       iconf(:) = 1
!
       do icombi = 1, ncombi
!
         dDe = 0.0d0
         dD0 = 0.0d0
         dE  = 0.0d0
         dH  = 0.0d0
         dS  = 0.0d0
         dG  = 0.0d0
!
         do imol = 1, nmol
           dDe = dDe + reac%nu(inu(imol))*mol(inu(imol))%              &
                            conf(order(inu(imol),iconf(inu(imol))))%Escf
           dD0 = dD0 + reac%nu(inu(imol))*mol(inu(imol))%              &
                              conf(order(inu(imol),iconf(inu(imol))))%D0
           dE  = dE  + reac%nu(inu(imol))*mol(inu(imol))%              &
                          conf(order(inu(imol),iconf(inu(imol))))%Etherm
           dH  = dH  + reac%nu(inu(imol))*mol(inu(imol))%              &
                          conf(order(inu(imol),iconf(inu(imol))))%Htherm
           dS  = dS  + reac%nu(inu(imol))*mol(inu(imol))%              &
                          conf(order(inu(imol),iconf(inu(imol))))%Stherm
           dG  = dG  + reac%nu(inu(imol))*mol(inu(imol))%              &
                          conf(order(inu(imol),iconf(inu(imol))))%Gtherm
         end do
!
         Keq = dexp(-dG/Rjul/temp)
!
         reac%dDe = reac%dDe + dDe
         reac%dD0 = reac%dD0 + dD0
         reac%dE  = reac%dE  + dE
         reac%dH  = reac%dH  + dH
         reac%dS  = reac%dS  + dS
         reac%dG  = reac%dG  + dG
         reac%Keq = reac%Keq + Keq
!
         write(*,fmt2) (order(inu(imol),iconf(inu(imol))),imol=1,mmol),&
                        dDe/1000,dD0/1000,dE/1000,dH/1000,             &
                                                temp*dS/1000,dG/1000,Keq  ! FLAG: change format keq
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
       subroutine print_thermo(temp,reac)
!
       use datatypes
!
       implicit none
!
       type(reaction),intent(in)  ::  reac    !  Reaction information
       real(kind=8),intent(in)    ::  temp    !  Temperature
!
!~        write(*,'(6(1X,F10.4),1X,E10.4)') reac%dDe/1000,reac%dD0/1000,  &
       write(*,'(6(1X,F10.4),1X,D10.4)') reac%dDe/1000,reac%dD0/1000,  &
                                 reac%dE/1000,reac%dH/1000,            &
                                 temp*reac%dS/1000,reac%dG/1000,reac%Keq  ! FLAG: change format keq
       write(*,*)
!
       return
       end subroutine print_thermo
!
!======================================================================!
