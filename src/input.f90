!======================================================================!
!
       module input
       implicit none
!
       private
       public  ::  read_inp,                                           &
                   read_g16
!
       contains
!
!======================================================================!
!
       subroutine read_g16(nmol,mol,fcalc,fchck,schm)
!
       use datatypes
       use utils,     only: print_end
       use g16files
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(nmol),intent(inout)  ::  mol      !  Input file name
       type(scheme),intent(in)                       ::  schm     !
       integer,intent(in)                            ::  nmol     !  Number of atoms
       character(len=8),intent(in)                   ::  fcalc    !  Calculation information flag
       logical,intent(in)                            ::  fchck
!
! Local variables
!
       character(len=5),dimension(:),allocatable     ::  atname  !  Atom names
       real(kind=8),dimension(:),allocatable         ::  atmass  !  Atomic masses
       real(kind=8)                                  ::  mass    !  Molecule mass
       integer,dimension(:),allocatable              ::  znum    !  Atomic number
       integer                                       ::  nat     !  Number of atoms in the input file
       integer                                       ::  imol    !  Chemical species index
       integer                                       ::  iconf   !  Configuration index
!
! Reading Gaussian16 input files
!
       do imol = 1, nmol
         call chk_log(mol(imol)%conf(1)%inp,mol(imol)%nat,fcalc,fchck)
!
         mol(imol)%conf(1)%chiral = .FALSE.
         mol(imol)%conf(1)%nequi  = 1
         mol(imol)%conf(1)%rotdof = 3                  ! TODO: check if the species are linear
         mol(imol)%conf(1)%dof    = 3*mol(imol)%nat - 3                &
                                              - mol(imol)%conf(1)%rotdof
!
         allocate(mol(imol)%conf(1)%freq(mol(imol)%conf(1)%dof),       &  
                  mol(imol)%conf(1)%inten(mol(imol)%conf(1)%dof),      &
                  mol(imol)%conf(1)%coord(3,mol(imol)%nat),            &
                  mol(imol)%atname(mol(imol)%nat),                     &
                  mol(imol)%atmass(mol(imol)%nat),                     &
                  mol(imol)%znum(mol(imol)%nat))
!
         call read_log(fcalc,mol(imol)%conf(1)%inp,                    &
                       mol(imol)%conf(1)%auxinp,                       &
                       mol(imol)%nat,                                  &
                       mol(imol)%conf(1)%coord,                        &
                       mol(imol)%atname,                               &
                       mol(imol)%znum,                                 &
                       mol(imol)%atmass,                               &
                       mol(imol)%conf(1)%dof,                          &
                       mol(imol)%conf(1)%freq,                         &
                       mol(imol)%conf(1)%inten,                        &
                       mol(imol)%conf(1)%moment,                       &
                       mol(imol)%mass,                                 &
                       mol(imol)%conf(1)%qel,                          &
                       mol(imol)%conf(1)%Escf,                         &
                       mol(imol)%wfn,                                  &
                       mol(imol)%schm)
!
         if ( mol(imol)%nconf .gt. 1 ) then
!
           allocate(atname(mol(imol)%nat),                             &
                    atmass(mol(imol)%nat),                             &
                    znum(mol(imol)%nat))
!
           do iconf = 2, mol(imol)%nconf
!
             call chk_log(mol(imol)%conf(iconf)%inp,nat,fcalc,fchck)
!
             mol(imol)%conf(iconf)%chiral = .FALSE.
             mol(imol)%conf(iconf)%nequi  = 1
             mol(imol)%conf(iconf)%rotdof = 3 ! TODO: check if the species are linear
             mol(imol)%conf(iconf)%dof    = 3*mol(imol)%nat - 3        &
                                          - mol(imol)%conf(iconf)%rotdof
!
             if ( nat .ne. mol(imol)%nat ) then
               write(*,*)
               write(*,'(2X,68("="))')
               write(*,'(3X,A)') 'ERROR:  Number of atoms for diff'//  &
                                                      'erent conformers'
               write(*,'(3X,A)') '        of the same molecule doe'//  &
                                                           's not match'
               write(*,*) 
               write(*,'(3X,A,X,I3)') 'Atoms in the input file '//     &
                                        trim(mol(imol)%conf(1)%inp)//  &
                                                      ' :',mol(imol)%nat
               write(*,'(3X,A,X,I3)') 'Atoms in the input file '//     &
                                    trim(mol(imol)%conf(iconf)%inp)//  &
                                                                ' :',nat

               write(*,'(2X,68("="))')
               write(*,*) 
               call print_end()
             end if
!
             allocate(mol(imol)%conf(iconf)%                           &
                                freq(mol(imol)%conf(iconf)%dof),       &  
                      mol(imol)%conf(iconf)%                           &
                                inten(mol(imol)%conf(iconf)%dof),      &
                      mol(imol)%conf(iconf)%coord(3,mol(imol)%nat))
!
             call read_log(fcalc,mol(imol)%conf(iconf)%inp,            &
                           mol(imol)%conf(iconf)%auxinp,               &
                           mol(imol)%nat,                              &
                           mol(imol)%conf(iconf)%coord,                &
                           atname,znum,atmass,                         &
                           mol(imol)%conf(iconf)%dof,                  &
                           mol(imol)%conf(iconf)%freq,                 &
                           mol(imol)%conf(iconf)%inten,                &
                           mol(imol)%conf(iconf)%moment,               &
                           mass,mol(imol)%conf(iconf)%qel,             &
                           mol(imol)%conf(iconf)%Escf,                 &
                           mol(imol)%wfn,                              &
                           mol(imol)%schm)
!
           end do
!
           deallocate(atname,atmass,znum)
!
         end if
!
       end do
!
       return
       end subroutine read_g16
!
!======================================================================!
!
       subroutine read_inp(inp,nmol,mol,cutoff,fact,hwhm,iline,doir,   &
                           ffree,fentha,fentro,forder,fcalc,fchck,     &
                           fenan,fpermu,fsoln,nreac,reac,mconf,fscreen,&
                           frota,rmsdmax,maemax,baemax,dsolv,msolv,schm)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(:),allocatable,intent(out)  ::  mol      !  Molecules information
       type(reaction),dimension(:),allocatable,intent(out)  ::  reac     !  Reactions information
       type(scheme),intent(out)                             ::  schm     !
       character(len=leninp),intent(in)                     ::  inp      !  General input file name
       character(len=8),intent(out)                         ::  ffree    !  Free energy calculation flag
       character(len=8),intent(out)                         ::  fentha   !  Enthalpy calculation flag
       character(len=8),intent(out)                         ::  fentro   !  Entropy calculation flag
       character(len=8),intent(out)                         ::  fcalc    !  Calculation information flag
       character(len=8),intent(out)                         ::  forder   !  Conformations order flag
       character(len=8),intent(out)                         ::  frota    !  Rotation method flag
       real(kind=8),intent(out)                             ::  cutoff   !  Cutoff frequency
       real(kind=8),intent(out)                             ::  fact     !  Frequencies scaling factor
       real(kind=8),intent(out)                             ::  hwhm     !  Full width at half maximum
       real(kind=8),intent(out)                             ::  rmsdmax  !  Maximum value for RMSD
       real(kind=8),intent(out)                             ::  maemax   !  Maximum value for MAE
       real(kind=8),intent(out)                             ::  baemax   !  Maximum value for BAE
       real(kind=8),intent(out)                             ::  dsolv    !
       real(kind=8),intent(out)                             ::  msolv    !
       integer,intent(out)                                  ::  nmol     !  Number of chemical species
       integer,intent(out)                                  ::  nreac    !  Number of reactions
       integer,intent(out)                                  ::  mconf    !  Maximum number of conformers
       integer,intent(out)                                  ::  iline    !  Broadening function
       logical,intent(out)                                  ::  doir     !  IR calculation flag
       logical,intent(out)                                  ::  fchck    !  
       logical,intent(out)                                  ::  fenan    !  Enantiomers calculation flag
       logical,intent(out)                                  ::  fpermu   !  Permutations calculation flag
       logical,intent(out)                                  ::  fsoln    !  Standard state flag
       logical,intent(out)                                  ::  fscreen  !  Screening calculation flag
!
! Local variables
!
       character(len=lenline)                               ::  line     !  
       character(len=lenline)                               ::  key      !
       real(kind=8)                                         ::  summ     !
       integer                                              ::  imol     !
       integer                                              ::  iconf    !
       integer                                              ::  io       !  Input/Output status
!
! Setting defaults  ! FLAG: missing default file names (breaks if **SYSTEM is not defined)
!
       fact  = 1.0d0
       hwhm  = 1.0d0
       doir  = .FALSE.
       iline = 1      ! 1. Gaussian | 2. Lorentzian | 3. Pseudo-Voigt
!
       frota   = 'QUATCTK'
       rmsdmax = 0.05d0
       maemax  = 0.05d0
       baemax  = 0.5d0
       fscreen = .FALSE.
!
       ffree  = 'MIX'
       fentha = 'RRHO'
       fentro = 'RRHO'
!
       forder = 'GORDER'
       fcalc  = 'ALL'
       fchck  = .TRUE.
!
       cutoff = 100.0d0
!
       fenan  = .FALSE.
       fpermu = .FALSE.
       fsoln  = .FALSE.
!
       dsolv  = 1.59    ! CCl4 data
       msolv  = 153.82  ! CCl4 data
!
       schm%wfn   = 'SCF'
       schm%fschm = 'NORMAL'
!
       nreac  = 0
!
! Reading general input file
! --------------------------
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
!
       if ( io .ne. 0 ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')    'ERROR:  Missing input file'
         write(*,*)
         write(*,'(3X,A)')    'Input file '//trim(inp)//' not foun'//  &
                                            'd in the current directory'
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
! Counting the number of *MOLECULE blocks specified
! 
       nmol   = 0
!
       do
! Reading input file line
         read(uniinp,'(A)',iostat=io) line 
! If the end of the input file is reached exit
         if ( io /= 0 ) exit
! If reads a white line then reads the next line
         if ( len_trim(line) == 0 ) cycle
! Processing the line read 
         call chkcomment(line,key)
! If the line just contains a comment reads the next line
         if ( len_trim(key) == 0 ) cycle
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Reading the different input file blocks       
         select case (key)
           case ('*MOL','*MOLEC','*MOLECULE','*AGG','*AGGREGATE')
             nmol = nmol + 1
         end select  
       end do
!
       if ( nmol .eq. 0 ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')    'ERROR:  Electronic structure inform'//  & 
                                                         'ation missing'
         write(*,*)
         write(*,'(3X,A)')    'Section *MOLECULE of block **SYSTEM'//  &
                                                        ' not specified'
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
       allocate(mol(nmol))
!
       rewind(uniinp)
!
! Reading input file blocks  
!
       do
! Reading input file line
         read(uniinp,'(A)',iostat=io) line 
! If the end of the input file is reached exit
         if ( io /= 0 ) exit
!~          write(*,'(A)') trim(line)  ! FLAG: dump of input data file
! If reads a white line then reads the next line
         if ( len_trim(line) == 0 ) cycle
! Processing the line read 
         call chkcomment(line,key)
! If the line just contains a comment reads the next line
         if ( len_trim(key) == 0 ) cycle
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Reading the different input file blocks       
         select case (key)
           case ('**SYS','**SYSTEM')
!~              write(*,*) 
!~              write(*,*) 'Reading **SYSTEM block'
!~              write(*,*) 
!
             call findline(line,'blck','**SYSTEM')
!
             call read_sys(line,'**SYSTEM',nmol,mol,schm)      
!      
           case ('**PROP','**PROPERTIES')
!~              write(*,*) 
!~              write(*,*) 'Reading **PROPERTIES block'
!~              write(*,*) 
!
             call findline(line,'blck','**PROPERTIES')
!
             call read_prop(line,'**PROPERTIES',cutoff,fact,hwhm,iline,&
                            doir,ffree,fentha,fentro,forder,fcalc,     &
                            fchck,fenan,fpermu,fsoln,fscreen,frota,    &
                            rmsdmax,maemax,baemax,dsolv,msolv) 
!      
           case ('**REACT','**REAC','**REACTOR')
!~              write(*,*) 
!~              write(*,*) 'Reading **REACTOR block'
!~              write(*,*) 
!
             call findline(line,'blck','**REACTOR')
!
             call read_reactor(line,'**REACTOR',nmol,mol,nreac,reac)
!
           case default
             write(*,*)
             write(*,'(2X,68("="))')
             write(*,'(3X,A)') 'ERROR:  Unknown block from input file'
             write(*,*) 
             write(*,'(3X,A)') 'Block '//trim(key)//' not known'
             write(*,'(2X,68("="))')
             write(*,*) 
             call print_end()
         end select  
       end do
! Closing the input file     
       close(uniinp)
!
! General settings and fatal errors check
!
       mconf = 1
       do io = 1, nmol
         mconf = max(mconf,mol(io)%nconf)
       end do
!
       if ( trim(ffree) .ne. 'MIX' ) then
         fentro = ffree
         fentha = ffree
       end if
!
       do imol = 1, nmol
!
         summ = 0.0d0
         do iconf = 1, mol(imol)%nconf
           summ = summ + mol(imol)%conf(iconf)%weight   
         end do
!
         do iconf = 1, mol(imol)%nconf
           mol(imol)%conf(iconf)%weight = mol(imol)%conf(iconf)%weight &
                                                                   /summ
         end do
!
       end do
!
       return
       end subroutine read_inp
!
!======================================================================!
!
       subroutine read_sys(key,blck,nmol,mol,schm)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(nmol),intent(inout)  ::  mol   !  Molecules information
       type(scheme),intent(inout)                    ::  schm
       character(len=lenline),intent(inout)          ::  key   !
       character(len=*),intent(in)                   ::  blck  !  Block name
       integer,intent(in)                            ::  nmol  !
!
! Local variables
!
       integer                                       ::  posi  !
       integer                                       ::  imol  !
!
! Reading SYSTEM block sections 
! -----------------------------
!
       imol = 0
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.TITLE')
!~              write(*,*) '  Reading .TITLE option'
!~              write(*,*)
!
!~              read(uniinp,*) tgrp     ! FLAG: check errors
!~              tgrp = adjustl(tgrp)
!
             call findline(key,'blck',blck)
!
           case ('*MOL','*MOLEC','*MOLECULE')
!
!~              write(*,*) '  Reading *MOLECULE section'
!~              write(*,*)
!
             imol = imol + 1
!
             call findline(key,'sect','*MOLECULE')
!
             call read_mol(key,'*MOLECULE',blck,nmol,mol,imol,schm)
!
           case ('*SCH','*SCHM','*SCHEME')
!
!~              write(*,*) '  Reading *SCHEME section'
!~              write(*,*)
!
             call findline(key,'sect','*SCHEME')
!
             call read_schm(key,'*SCHEME',blck,schm)
!
           case ('**END')
!~              write(*,*) 'Exiting from **SYSTEM block'
!~              write(*,*)
             return
           case default
             call unksect(key,blck)
         end select  
       end do
!
       return
       end subroutine read_sys
!
!======================================================================!
!
       subroutine read_mol(key,sect,blck,nmol,mol,imol,schm)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(nmol),intent(inout)  ::  mol   !  Molecules information
       type(scheme),intent(in)                       ::  schm
       character(len=lenline),intent(inout)          ::  key   !  
       character(len=*),intent(in)                   ::  sect  !  Section name
       character(len=*),intent(in)                   ::  blck  !  Block name
       integer,intent(in)                            ::  imol  !
       integer,intent(in)                            ::  nmol  !
!
! Local variables
!
       character(len=lenline)                        ::  line  !
       character(len=lenline)                        ::  arg   !  
       character(len=lenline)                        ::  str   !  
       integer                                       ::  posi  !
       integer                                       ::  i     !
!
! Reading MOLECULE section keywords 
! ---------------------------------
!
       mol(imol)%nconf  = 1
       mol(imol)%npermu = 1
!
       mol(imol)%wfn    = 'SCF'
       mol(imol)%readw  = .FALSE.
!
       allocate(mol(imol)%frag(1))
       mol(imol)%frag(1) = 1
       mol(imol)%nfrag   = 1
!
       mol(imol)%conc    = 1.0d0
!
       mol(imol)%schm%wfn   = schm%wfn
       mol(imol)%schm%fschm = schm%fschm
!
       write(mol(imol)%molname,'(I5)') imol
       mol(imol)%molname = adjustl(mol(imol)%molname) 
       mol(imol)%molname = 'mol-'//trim(mol(imol)%molname)
!
       mol(imol)%phase   = 'soln'
!  
       line = key
!
       if ( key(1:1) .ne. '.' ) then
         do
           if ( len_trim(line) == 0 ) exit
! Saving the keyword 
           posi = scan(key,'=') 
           if ( posi .ne. 0 ) then 
             line = key(posi+1:)  
             line = adjustl(line)
!
             key  = key(:posi-1)
             key  = lowercase(key)
           else
             call errkey('section',sect)
           end if
! Saving the arguments
           select case (key)
             case ('nconf')
!
               call chkkeyarg(key,line,arg)
!
               read(arg,*) mol(imol)%nconf   
!        
             case ('name')
!
               call chkkeyarg(key,line,arg)
!
               mol(imol)%molname = arg(:len(mol(imol)%molname))
!
             case ('phase')
!
               call chkkeyarg(key,line,arg)
!
               mol(imol)%phase = arg(:len(mol(imol)%phase))
!
             case default
               call unkkeysect(key,sect)  ! FLAG: create UNKINP subroutine
           end select  
!
           key = line
!
         end do
!
         call findline(key,'sect',sect)
       end if
!
       allocate(mol(imol)%conf(mol(imol)%nconf))
!
! Reading MOLECULE section options 
!
       if ( mol(imol)%nconf .gt. 1 ) then
!
         do i = 1, mol(imol)%nconf 
           write(mol(imol)%conf(i)%inp,'(I5)') i
           mol(imol)%conf(i)%inp = adjustl(mol(imol)%conf(i)%inp) 
           mol(imol)%conf(i)%inp = trim(mol(imol)%molname)//'_conf-'// &
                                     trim(mol(imol)%conf(i)%inp)//'.log'
!
           mol(imol)%conf(i)%symnum = 1  ! TODO: change by zero to decide if it must be read from input
           mol(imol)%conf(i)%weight = 0.0d0
         end do
!
       else
!
         mol(imol)%conf(1)%inp = trim(mol(imol)%molname)//'.log'
         mol(imol)%conf(1)%symnum = 1    ! TODO: change by zero to decide if it must be read from input
         mol(imol)%conf(1)%weight = 0.0d0
!
       end if
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.WFNSCHM','.SCHMWFN','.WAVEFUNCTION-SCHEME')
!
!~              write(*,*) '    Reading .WAVEFUNCTION option'
!~              write(*,*)
!
             read(uniinp,*) mol(imol)%schm%wfn 
!
             call select_wfn(mol(imol)%schm%wfn)  
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case ('.SCH','.SCHM','.SCHEME')
!
!~              write(*,*) '    Reading .TYPE option'
!~              write(*,*)
!
             read(uniinp,*) str
!
             call select_schm(str,mol(imol)%schm%fschm)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case ('.FILE','.FILES')
!
!~              write(*,*) '    Reading .FILE option'
!~              write(*,*)
             if ( trim(mol(imol)%schm%fschm) .eq. 'NORMAL' ) then
!
               do i = 1, mol(imol)%nconf        ! TODO: check if input names are introduced correctly
                 read(uniinp,*) mol(imol)%conf(i)%inp
                 mol(imol)%conf(i)%inp = adjustl(mol(imol)%conf(i)%inp) 
               end do
!
             else
!
               do i = 1, mol(imol)%nconf        ! TODO: check if input names are introduced correctly
!
                 if ( trim(mol(imol)%schm%fschm) .eq. 'LLSOL' ) then
                   allocate(mol(imol)%conf(i)%auxinp(2))
                 else if ( trim(mol(imol)%schm%fschm) .eq. 'HL' ) then
                   allocate(mol(imol)%conf(i)%auxinp(1))
                 end if
!
                 read(uniinp,*) mol(imol)%conf(i)%inp,                 &
                                             mol(imol)%conf(i)%auxinp(:)
                 mol(imol)%conf(i)%inp = adjustl(mol(imol)%conf(i)%inp) 
                 mol(imol)%conf(i)%auxinp(:) =                         &
                                    adjustl(mol(imol)%conf(i)%auxinp(:)) 
!~                  if ( trim(schm%fschm) .eq. 'LLSOL' ) then
!~                    mol(imol)%conf(i)%auxinp(1) =                       &
!~                                     adjustl(mol(imol)%conf(i)%auxinp(1)) 
!~                    mol(imol)%conf(i)%auxinp(2) =                       &
!~                                     adjustl(mol(imol)%conf(i)%auxinp(2)) 
!~                  else if ( trim(schm%fschm) .eq. 'HL' ) then
!~                    mol(imol)%conf(i)%auxinp(1) =                       &
!~                                     adjustl(mol(imol)%conf(i)%auxinp(1))
!~                  end if
               end do
!
             end if
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case ('.READ','.READFILE','.READFILES')
!
!~              write(*,*) '    Reading .READ option'
!~              write(*,*)
!
! Reading file name
             read(uniinp,'(A)') line         ! FLAG: check errors
             line = adjustl(line)
! Reading input files names
!~              deallocate(mol(imol)%conf)   ! FLAG:  check what happens if nconf is unknown
!
             call countlines(trim(line),uniaux,i)
!
             if ( i .ne. mol(imol)%nconf ) then 
               write(*,*)
               write(*,'(2X,68("="))')
               write(*,'(3X,A)') 'ERROR:  Number of lines in the i'//  &
                                   'nput file does not match the number'
               write(*,'(3X,A)') '        of conformers specified'
               write(*,*)
               write(*,'(2X,A)') 'Input file name: '//trim(line)
               write(*,*)
               write(*,'(3X,A,X,I4)') 'Number of lines in the inpu'//  &
                                                              't file',i
               write(*,'(3X,A,X,I4)') 'Number of conformers specif'//  &
                                                   'ied',mol(imol)%nconf
               write(*,'(2X,68("="))')
               call print_end()
             end if             
!
             open(unit=uniaux,file=trim(line),action='read',           &
                  status='old',iostat=i)
!
             if ( trim(mol(imol)%schm%fschm) .eq. 'NORMAL' ) then

               do i = 1, mol(imol)%nconf        ! TODO: check if input names are introduced correctly
                 read(uniaux,*) mol(imol)%conf(i)%inp
                 mol(imol)%conf(i)%inp = adjustl(mol(imol)%conf(i)%inp) 
               end do
!
             else 
!
               do i = 1, mol(imol)%nconf        ! TODO: check if input names are introduced correctly
!
                 select case (trim(mol(imol)%schm%fschm))
                   case ('HLSOL')
                     allocate(mol(imol)%conf(i)%auxinp(3))
                   case ('LLSOL')
                     allocate(mol(imol)%conf(i)%auxinp(2))
                   case ('HL')
                     allocate(mol(imol)%conf(i)%auxinp(1))
                 end select
!
                 read(uniaux,*) mol(imol)%conf(i)%inp,                 &
                                             mol(imol)%conf(i)%auxinp(:)
                 mol(imol)%conf(i)%inp = adjustl(mol(imol)%conf(i)%inp) 
                 mol(imol)%conf(i)%auxinp(:) =                         &
                                    adjustl(mol(imol)%conf(i)%auxinp(:)) 
!~                  if ( trim(schm%fschm) .eq. 'LLSOL' ) then
!~                    mol(imol)%conf(i)%auxinp(1) =                       &
!~                                     adjustl(mol(imol)%conf(i)%auxinp(1)) 
!~                    mol(imol)%conf(i)%auxinp(2) =                       &
!~                                     adjustl(mol(imol)%conf(i)%auxinp(2)) 
!~                  else if ( trim(schm%fschm) .eq. 'HL' ) then
!~                    mol(imol)%conf(i)%auxinp(1) =                       &
!~                                     adjustl(mol(imol)%conf(i)%auxinp(1))
!~                  end if
               end do
!
             end if
!
             close(uniaux)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case ('.SYMNUM','.SYMMETRYNUM','.SYMMETRYNUMBER')
!
!~              write(*,*) '    Reading .SYMNUM option'
!~              write(*,*)
!
             read(uniinp,*) mol(imol)%conf(:)%symnum   ! FLAG: check if bad introduced
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case ('.WEIGHT','.WEIGHTS','.POP','.POPULATION',            &
                                                         '.POPULATIONS')
!
!~              write(*,*) '    Reading .SYMNUM option'
!~              write(*,*)
!  
             mol(imol)%readw = .TRUE.
             read(uniinp,*) mol(imol)%conf(:)%weight   ! FLAG: check if bad introduced
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case ('.NFRAG','.NFRAGMENTS','.FRAGMENTS')
!
!~              write(*,*) '    Reading .NFRAG option'
!~              write(*,*)
!
             read(uniinp,'(A)') line
             line = adjustl(line)
!
             key = ''
             mol(imol)%nfrag = 0
!
             do
               if ( len_trim(line) .eq. 0 ) exit
               mol(imol)%nfrag = mol(imol)%nfrag + 1
!
               posi = scan(line,' ')
               key  = trim(key)//' '//line(:posi-1)
               line = line(posi+1:)
               line = adjustl(line)
             end do
!
             deallocate(mol(imol)%frag)
             allocate(mol(imol)%frag(mol(imol)%nfrag))
!
             read(key,*) mol(imol)%frag(:)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case ('.CONC','.CONCENTRATION')
!
!~              write(*,*) '    Reading .CONC option'
!~              write(*,*)
!
             read(uniinp,*) mol(imol)%conc
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case ('.WFN','.WAVEFUNCTION')
!
!~              write(*,*) '    Reading .WFN option'
!~              write(*,*)
!
             read(uniinp,*) mol(imol)%wfn
!
             call select_wfn(mol(imol)%wfn)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) GOTO 1000   ! FLAG: check exit
!
           case default
             call unkopt(key,sect,blck)
         end select  
       end do
!
1000   allocate(mol(imol)%pop(mol(imol)%nconf))
!
       return
       end subroutine read_mol
!
!======================================================================!
!
       subroutine read_prop(key,blck,cutoff,fact,hwhm,iline,doir,      &
                            ffree,fentha,fentro,forder,fcalc,fchck,    &
                            fenan,fpermu,fsoln,fscreen,frota,rmsdmax,  &
                            maemax,baemax,dsolv,msolv)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lenline),intent(inout)  ::  key      !
       character(len=*),intent(in)           ::  blck     !  Block name
       character(len=8),intent(inout)        ::  ffree    !  Free energy calculation flag
       character(len=8),intent(inout)        ::  fentha   !  Enthalpy calculation flag
       character(len=8),intent(inout)        ::  fentro   !  Entropy calculation flag
       character(len=8),intent(inout)        ::  fcalc    !  Calculation information flag
       character(len=8),intent(inout)        ::  forder   !  Conformations order flag
       character(len=8),intent(inout)        ::  frota    !  Rotation method flag
       real(kind=8),intent(inout)            ::  cutoff   !  Cutoff frequency
       real(kind=8),intent(inout)            ::  fact     !  Frequencies scaling factor
       real(kind=8),intent(inout)            ::  hwhm     !  Half width at half maximum
       real(kind=8),intent(inout)            ::  rmsdmax  !  Maximum value for RMSD
       real(kind=8),intent(inout)            ::  maemax   !  Maximum value for MAE
       real(kind=8),intent(inout)            ::  baemax   !  Maximum value for BAE
       real(kind=8),intent(inout)            ::  dsolv    !
       real(kind=8),intent(inout)            ::  msolv    !
       integer,intent(inout)                 ::  iline    !  Broadening function
       logical,intent(inout)                 ::  fchck    !  
       logical,intent(inout)                 ::  doir     !  IR calculation flag
       logical,intent(inout)                 ::  fenan    !  Enantiomers calculation flag
       logical,intent(inout)                 ::  fpermu   !  Permutations calculation flag
       logical,intent(inout)                 ::  fsoln    !  Standard state flag
       logical,intent(inout)                 ::  fscreen  !  Screening calculation flag
!
! Local variables
!
       integer                               ::  posi     !
!
! Reading PROPERTIES block sections 
! ---------------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('*GENERAL','*SETTINGS')
!
!~              write(*,*) '  Reading *GENERAL section'
!~              write(*,*)
!
             call findline(key,'sect','*GENERAL')
!
             call read_general(key,'*GENERAL',blck,forder,fcalc,fchck)
!
!~            case ('.TITLE')
!~              write(*,*) '  Reading .TITLE option'
!~              write(*,*)
!
!~              read(uniinp,*) tgrp
!~              tgrp = adjustl(tgrp)
!
!~              call findline(key,'blck',blck)
!
           case ('*FREQ','*FREQUENCIES','*IRSPECTRA','*IRSPEC')
!
!~              write(*,*) '  Reading *FREQ section'
!~              write(*,*)
!
             call findline(key,'sect','*FREQ')
!
             call read_freq(key,'*FREQ',blck,fact,hwhm,iline,doir)
!
           case ('*SCREEN','*SCREENING')
!
!~              write(*,*) '  Reading *SCREEN section'
!~              write(*,*)
!
             call findline(key,'sect','*SCREEN')
!
             call read_screen(key,'*SCREEN',blck,frota,rmsdmax,        &
                              maemax,baemax)
!
             fscreen = .TRUE.
!
           case ('*THERMO','*FREE','*KEQ','*EQUILIBRIUM')
!
!~              write(*,*) '  Reading *THERMO section'
!~              write(*,*)
!
             call findline(key,'sect','*THERMO')
!
             call read_thermo(key,'*THERMO',blck,cutoff,ffree,fentha,  &
                              fentro,fenan,fpermu,fsoln,dsolv,msolv)
!
           case ('**END')
!~              write(*,*) 'Exiting from **PROPERTIES block'
!~              write(*,*)
             return
           case default
             call unksect(key,blck)
         end select  
       end do
!
       return
       end subroutine read_prop
!
!======================================================================!
!
       subroutine read_general(key,sect,blck,forder,fcalc,fchck)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lenline),intent(inout)  ::  key     !  
       character(len=*),intent(in)           ::  sect    !  Section name
       character(len=*),intent(in)           ::  blck    !  Block name
       character(len=8),intent(inout)        ::  forder  !  Conformations order flag
       character(len=8),intent(inout)        ::  fcalc   !  Calculation information flag
       logical,intent(inout)                 ::  fchck   !  
!
! Local variables
!
       character(len=lenline)                ::  line    !
       integer                               ::  posi    !
!
! Reading FREQ section options 
! ----------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.ORDER')
!
!~              write(*,*) '    Reading .ORDER option'
!~              write(*,*)
!
             read(uniinp,*) forder    ! FLAG: check errors
!
             select case ( forder )
               case ('EORDER','ENERGY','POTENTIAL')
                 forder = 'EORDER' 
               case ('DORDER','VORDER','ZPE','ZPVE')
                 forder = 'DORDER'
               case ('GORDER','FREE','GIBBS')
                 forder = 'GORDER'
               case default        
                 call errkeyoptsectblck(trim(forder),key,sect,blck)  
             end select
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.JOBTYPE','.TYPE','.JOB','.READ','.READONLY')
!
!~              write(*,*) '    Reading .JOBTYPE option'
!~              write(*,*)
!
             read(uniinp,*) line    ! FLAG: check errors
             line = adjustl(line)
!
             select case ( line )
               case ('EONLY','ENERGY','SP')
                 fcalc = 'EONLY'
               case ('OPT','OPTIMIZATION')
                 fcalc = 'OPT'
               case ('FREQ','FREQUENCIES','FREQUENCY')
                 fcalc = 'FREQ'
               case ('ALL','OF','OPTFREQ','FULL','FREE','GIBBS')
                 fcalc = 'ALL'
               case default        
                 call errkeyoptsectblck(trim(line),key,sect,blck)  
             end select
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.CHECK','.CHCK')
!
!~              write(*,*) '    Reading .CHECK option'
!~              write(*,*)
!
             fchck = .TRUE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.NOCHECK','.NOCHCK')
!
!~              write(*,*) '    Reading .CHECK option'
!~              write(*,*)
!
             fchck = .FALSE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect,blck)
         end select  
       end do
!
       return
       end subroutine read_general
!
!======================================================================!
!
       subroutine read_freq(key,sect,blck,fact,hwhm,iline,doir)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lenline),intent(inout)  ::  key    !  
       character(len=*),intent(in)           ::  sect   !  Section name
       character(len=*),intent(in)           ::  blck   !  Block name
       real(kind=8),intent(inout)            ::  fact   !  Frequencies scaling factor
       real(kind=8),intent(inout)            ::  hwhm   !  Full width at half maximum
       integer,intent(inout)                 ::  iline  !  Broadening function
       logical,intent(inout)                 ::  doir   !  IR calculation flag
!
! Local variables
!
       character(len=15)                     ::  next   !
       integer                               ::  posi   !
!
! Reading FREQ section options 
! ----------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.PRINT','.PRINTIR')
!
             doir = .TRUE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.NOPRINT','.NOPRINTIR')
!
             doir = .FALSE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.SCALING','.FACTOR','.SCALINGFACTOR')
!
!~              write(*,*) '    Reading .SCALING option'
!~              write(*,*)
!
             read(uniinp,*) fact    ! FLAG: check errors
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.HWHM','.WIDTH')
!
!~              write(*,*) '    Reading .HWHM option'
!~              write(*,*)
!
             read(uniinp,*) hwhm    ! FLAG: check errors
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.PROFILE','.TYPE')
!
!~              write(*,*) '    Reading .PROFILE option'
!~              write(*,*)
!
             read(uniinp,*) next    ! FLAG: check errors
!
             next = uppercase(next)
             select case (trim(next))
               case('G','GAU','GAUS','GAUSS','GAUSSIAN')
                 iline = 1
               case('L','LOR','LORENT','LORENTZ','LORENTZIAN')
                 iline = 2
               case('V','VOIGT')
                 stop 'VOIGT profile not implemented yet!'
                 iline = 3
               case default
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')    'ERROR:  Invalid value intro'//  &
                                          'duced for .PROFILE option'
                 write(*,*)
                 write(*,'(3X,2(A))') 'Unrecognised value     : ',     &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, choose between : "gaus'//  &
                                        'sian", "lorentzian" or "voigt"'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect,blck)
         end select  
       end do
!
       return
       end subroutine read_freq
!
!======================================================================!
!
       subroutine read_screen(key,sect,blck,frota,rmsd,mae,bae)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lenline),intent(inout)  ::  key     !  
       character(len=*),intent(in)           ::  sect    !  Section name
       character(len=*),intent(in)           ::  blck    !  Block name
       character(len=8),intent(inout)        ::  frota   !  Rotation method flag
       real(kind=8),intent(inout)            ::  rmsd    !  Maximum value for RMSD
       real(kind=8),intent(inout)            ::  mae     !  Maximum value for MAE
       real(kind=8),intent(inout)            ::  bae     !  Maximum value for BAE
!
! Local variables
!
       integer                               ::  posi    !
!
! Reading SCREEN section options 
! ------------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.ROTATION','.ROTA','.METHOD','.TYPE')
!
!~              write(*,*) '    Reading .ROTATION option'
!~              write(*,*)
!
             read(uniinp,*) frota    ! FLAG: check errors
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.RMSDMAX','.RMSD','.MAXRMSD')
!
!~              write(*,*) '    Reading .RMSDMAX option'
!~              write(*,*)
!
             read(uniinp,*) rmsd    ! FLAG: check errors
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.MAEMAX','.MAE','.MAXMAE')
!
!~              write(*,*) '    Reading .MAEMAX option'
!~              write(*,*)
!
             read(uniinp,*) mae    ! FLAG: check errors
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.BAEMAX','.BAE','.MAXBAE')
!
!~              write(*,*) '    Reading .BAEMAX option'
!~              write(*,*)
!
             read(uniinp,*) bae    ! FLAG: check errors
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect,blck)
         end select  
       end do
!
       return
       end subroutine read_screen
!
!======================================================================!
!
       subroutine read_thermo(key,sect,blck,cutoff,ffree,fentha,       &
                              fentro,fenan,fpermu,fsoln,dsolv,msolv)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lenline),intent(inout)  ::  key     !  
       character(len=*),intent(in)           ::  sect    !  Section name
       character(len=*),intent(in)           ::  blck    !  Block name
       character(len=8),intent(inout)        ::  ffree   !  Free energy calculation flag
       character(len=8),intent(inout)        ::  fentha  !  Enthalpy calculation flag
       character(len=8),intent(inout)        ::  fentro  !  Entropy calculation flag
       real(kind=8),intent(inout)            ::  cutoff  !  Cutoff frequency
       real(kind=8),intent(inout)            ::  dsolv   ! 
       real(kind=8),intent(inout)            ::  msolv   !
       logical,intent(inout)                 ::  fenan   !  Enantiomers calculation flag
       logical,intent(inout)                 ::  fpermu  !  Permutations calculation flag
       logical,intent(inout)                 ::  fsoln   !  Standard state flag
!
! Local variables
!
       integer                               ::  posi    !
!
! Reading THERMO section options 
! ------------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.TYPE','.GTYPE','.TYPEG')
!
!~              write(*,*) '    Reading .TYPE option'
!~              write(*,*)
!
             read(uniinp,*) ffree        ! FLAG: check errors
!
             ffree = uppercase(ffree)
!
             ffree = chktype(ffree,key,sect,blck)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.HTYPE','.TYPEH')
!
!~              write(*,*) '    Reading .HTYPE option'
!~              write(*,*)
!
             read(uniinp,*) fentha        ! FLAG: check errors
!
             fentha = uppercase(fentha)
!
             fentha = chktype(fentha,key,sect,blck)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.STYPE','.TYPES')
!
!~              write(*,*) '    Reading .STYPE option'
!~              write(*,*)
!
             read(uniinp,*) fentro        ! FLAG: check errors
!
             fentro = uppercase(fentro)
!
             fentro = chktype(fentro,key,sect,blck)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.CUTOFF','.THR','.THRESH','.THRESHOLD')
!
!~              write(*,*) '    Reading .CUTOFF option'
!~              write(*,*)
!
             read(uniinp,*) cutoff           ! FLAG: check errors
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.ENAN','.ENANTIOMER','.ENANTIOMERS')
!
!~              write(*,*) '    Reading .ENAN option'
!~              write(*,*)
! 
             fenan = .TRUE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.NOENAN','.NOENANTIOMER','.NOENANTIOMERS')
!
!~              write(*,*) '    Reading .NOENAN option'
!~              write(*,*)
! 
             fenan = .FALSE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.INDIS','.INDIST','.INDISTINGUISHABLE','.INDISTI'//  &
                                                          'NGUISHABLES')
!
!~              write(*,*) '    Reading .INDIST option'
!~              write(*,*)
! 
             fpermu = .TRUE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.DISTIN','.DISTINGUISHABLE','.DISTINGUISHABLES')
!
!~              write(*,*) '    Reading .DISTIN option'
!~              write(*,*)
! 
             fpermu = .FALSE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.SOLN','.1MOLAR','.1M')
!
!~              write(*,*) '    Reading .SOLN option'
!~              write(*,*)
! 
             fsoln = .TRUE.
             read(uniinp,*) msolv,dsolv
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.GAS','.1BAR')
!
!~              write(*,*) '    Reading .GAS option'
!~              write(*,*)
! 
             fsoln = .FALSE.
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect,blck)
         end select  
       end do
!
       return
       end subroutine read_thermo
!
!======================================================================!
!
       function chktype(aux,key,sect,blck) result(str)
!
       use utils
!
       implicit none
!
       character(len=*),intent(in)   ::  aux     !  Keyword name
       character(len=*),intent(in)   ::  key     !  Option name  
       character(len=*),intent(in)   ::  sect    !  Section name
       character(len=*),intent(in)   ::  blck    !  Block name
       character(len=len_trim(aux))  ::  str     !  System keyword name
! 
       select case ( aux )
         case ('RRHO','GAS','IDEAL','IDEALGAS')
           str = 'RRHO' 
         case ('QRRHO','ROT','FREEROT','VIBROT') 
           str = 'QRRHO'          
         case ('QHO','RRQHO','QUASIHO') 
           str = 'QHO'  
         case ('MIX','CALC') 
           str = 'MIX'             
         case default        
           call errkeyoptsectblck(aux,key,sect,blck)  
       end select
!     
       return
       end function chktype
!
!======================================================================!
!
       subroutine read_reactor(key,blck,nmol,mol,nreac,reac)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(nmol),intent(in)            ::  mol    !  Molecules information
       type(reaction),dimension(:),allocatable,intent(out)  ::  reac   !  Reactions information
       character(len=lenline),intent(inout)                 ::  key    !
       character(len=*),intent(in)                          ::  blck   !  Block name
       integer,intent(in)                                   ::  nmol   !  Number of molecules
       integer,intent(inout)                                ::  nreac  !  Number of reactions
!
! Local variables
!
       integer                                              ::  posi    !
!
! Reading REACTOR block sections 
! ------------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('*PHASE','*PHASES')
!
!~              write(*,*) '  Reading *PHASES section'
!~              write(*,*)
!
             call findline(key,'sect','*PHASES')
!
!~              call read_phase(key,'*PHASES'l)
!
           case ('*REACTION','*REACTIONS')
!
!~              write(*,*) '  Reading *REACTIONS section'
!~              write(*,*)
!
             call findline(key,'sect','*REACTIONS')
!
             call read_reaction(key,'*REACTIONS',blck,nmol,mol,        &
                                nreac,reac)
!
           case ('**END')
!~              write(*,*) 'Exiting from **REACTOR block'
!~              write(*,*)
             return
           case default
             call unksect(key,blck)
         end select  
       end do
!
       return
       end subroutine read_reactor
!
!======================================================================!
!
       subroutine read_reaction(key,sect,blck,nmol,mol,nreac,reac)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(nmol),intent(in)            ::  mol      !  Molecules information
       type(reaction),dimension(:),allocatable,intent(out)  ::  reac     !  Reactions information
       character(len=*),intent(inout)                       ::  key      !  
       character(len=*),intent(in)                          ::  sect     !  Section name
       character(len=*),intent(in)                          ::  blck     !  Block name
       integer,intent(in)                                   ::  nmol     !  Number of molecules
       integer,intent(out)                                  ::  nreac    !  Number of reactions

!
! Local variables
!
       character(len=lenline)                               ::  line     !
       character(len=lename),dimension(nmol)                ::  molname  !
       integer                                              ::  posi     !
       integer                                              ::  ireac    !
!
! Reading REACTIONS section options 
! ---------------------------------
!
       nreac = 0
!
       do posi = 1, nmol
         molname(posi) = uppercase(mol(posi)%molname)
       end do
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.NREAC','.NREACT','.NREACTION')
!
!~              write(*,*) '    Reading .NREAC option'
!~              write(*,*)
!
             read(uniinp,*) nreac     ! FLAG: check no commentss
!
             allocate(reac(nreac))
!
             do posi = 1, nreac
!
               allocate(reac(posi)%nu(nmol))
               reac(posi)%nu(:) = 0.0d0
!
             end do
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.SCHEME','.SCHEMES')
!
!~              write(*,*) '    Reading .SCHEMES option'
!~              write(*,*)
!
             if ( nreac .eq. 0 ) then
                write(*,*)
                write(*,'(2X,68("="))')
                write(*,'(3X,A)') 'ERROR:  Number of reactions not'//  &
                                                        ' specified yet'
                write(*,*) 
                write(*,'(3X,A)') 'Option .NREAC must be specified'//  &
                                               ' before option .SCHEMES'
                write(*,'(2X,68("="))')
                write(*,*) 
                call print_end()
             end if
!
             do ireac = 1, nreac
               call findline(line,'sect',sect)
!
               line = adjustl(line)
               line = uppercase(line)
!
               posi = scan(line,'=')
!
               call read_scheme(line(:posi-1),-1,      &
                                ireac,reac(ireac),nmol,molname)
!
               call read_scheme(line(posi+1:),1,       &
                                ireac,reac(ireac),nmol,molname)
             end do
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return  
!
           case default
             call unkopt(key,sect,blck)
         end select  
       end do
!
       return
       end subroutine read_reaction
!
!======================================================================!
!
       subroutine read_scheme(inline,nusign,ireac,reac,       &
                              nmol,molname)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(reaction),intent(inout)                      ::  reac     !  Reactions information
       character(len=lename),dimension(nmol),intent(in)  ::  molname  !
       character(len=*),intent(in)                       ::  inline   !  
       integer,intent(in)                                ::  nmol     !  Number of molecules
       integer,intent(in)                                ::  ireac    !  Reaction index
       integer,intent(in)                                ::  nusign   !  Stoichiometric coefficient sign

!
! Local variables
!
       character(len=lenline)                            ::  line     !
       character(len=lename)                             ::  straux   !
       integer                                           ::  nu       !
       integer                                           ::  posi     !
!
! Reading right/left part of the reaction scheme
!
       line = inline
       line = adjustl(line)
!
       do while ( len_trim(line) .gt. 0 )
         straux = ''
         do
           select case ( line(1:1) )
!
             case ( '0':'99999' )
               straux = trim(straux)//line(1:1)
               line   = line(2:)
!
             case default
               if ( trim(straux) .eq. '' ) then
                 nu = 1
               else
                 read(straux,*) nu
               end if
! 
               line = adjustl(line)
               posi = scan(trim(line),' ')
!
               if ( posi .eq. 0 ) then
!
                 straux = line(:len(straux))  
                 line   = ''
!
               else 
                 straux = line(:posi-1)  
                 line   = line(posi+1:)
                 line   = adjustl(line)
!
                 posi = scan(line,' ')
!
                 line   = line(posi+1:)
                 line   = adjustl(line)
               end if
!
               straux = adjustl(straux)
!
               posi = findcv(nmol,molname,trim(straux))
!
               if ( posi .eq. 0 ) then
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A,1X,I2)') 'ERROR:  Unknown species '//  &
                                     'declared in reaction scheme',ireac
                 write(*,*) 
                 write(*,'(3X,A)') 'Chemical species '//trim(straux)// &
                               ' not specified in any section *MOLECULE'
                 write(*,'(3X,A)') ' of block **SYSTEM'
                 write(*,*)
                 write(*,'(1X,20(2X,A))') molname
                 write(*,'(2X,68("="))')
                 write(*,*) 
                 call print_end() 
               end if
!
               reac%nu(posi) = reac%nu(posi) + nusign*nu
!
               exit
!
           end select
         end do
       end do 
!
       return
       end subroutine read_scheme
!
!======================================================================!
!
       subroutine read_schm(key,sect,blck,schm)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(scheme),intent(inout)            ::  schm   !
       character(len=lenline),intent(inout)  ::  key    !  
       character(len=*),intent(in)           ::  sect   !  Section name
       character(len=*),intent(in)           ::  blck   !  Block name
!
! Local variables
!
       character(len=20)                     ::  str    !
       integer                               ::  posi   !
!
! Reading FREQ section options 
! ----------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.WFN','.WAVEFUNCTION')
!
!~              write(*,*) '    Reading .WAVEFUNCTION option'
!~              write(*,*)
!
             read(uniinp,*) schm%wfn 
!
             call select_wfn(schm%wfn)  
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.TYPE')
!
!~              write(*,*) '    Reading .TYPE option'
!~              write(*,*)
!
             read(uniinp,*) str
!
             call select_schm(str,schm%fschm)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect,blck)
         end select  
       end do
!
       return
       end subroutine read_schm
!
!======================================================================!
!
       end module input
!
!======================================================================!
