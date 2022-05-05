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
       subroutine read_g16(nmol,mol)
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
       integer,intent(in)                            ::  nmol     !  Number of atoms
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
         call chck_log(mol(imol)%conf(1)%inp,mol(imol)%nat)
!
         mol(imol)%conf(1)%chiral = .FALSE.
         mol(imol)%conf(1)%nequi  = 1
         mol(imol)%conf(1)%rotdof = 3                  ! FLAG: check if the species are linear
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
         call read_log(mol(imol)%conf(1)%inp,                          &
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
                       mol(imol)%conf(1)%Escf)
!
         if ( mol(imol)%nconf .gt. 1 ) then
!
           allocate(atname(mol(imol)%nat),                             &
                    atmass(mol(imol)%nat),                             &
                    znum(mol(imol)%nat))
!
           do iconf = 2, mol(imol)%nconf
!
             call chck_log(mol(imol)%conf(iconf)%inp,nat)
!
             mol(imol)%conf(iconf)%chiral = .FALSE.
             mol(imol)%conf(iconf)%nequi  = 1
             mol(imol)%conf(iconf)%rotdof = 3                  ! FLAG: check if the species are linear
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
             call read_log(mol(imol)%conf(iconf)%inp,                  &  ! FLAG : check that inputs include the same atoms
                           mol(imol)%nat,                              &
                           mol(imol)%conf(iconf)%coord,                &
                           atname,znum,atmass,                         &
                           mol(imol)%conf(iconf)%dof,                  &
                           mol(imol)%conf(iconf)%freq,                 &
                           mol(imol)%conf(iconf)%inten,                &
                           mol(imol)%conf(iconf)%moment,               &
                           mass,mol(imol)%conf(iconf)%qel,             &
                           mol(imol)%conf(iconf)%Escf)
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
       subroutine read_inp(inp,nmol,mol,thr,fact,fqvib,ffree,          &
                           fenan,fpermu,fsoln,nreac,reac)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(:),allocatable,intent(out)  ::  mol     !  Molecules information
       type(reaction),dimension(:),allocatable,intent(out)  ::  reac    !  Reactions information
       character(len=leninp),intent(in)                     ::  inp     !  General input file name
       character(len=8),intent(out)                         ::  fqvib   !  Qvib calculation flag
       character(len=8),intent(out)                         ::  ffree   !  Free energy calculation flag
       real(kind=8),intent(out)                             ::  thr     !  Threshold frequency
       real(kind=8),intent(out)                             ::  fact    !  Frequencies scaling factor
       integer,intent(out)                                  ::  nmol    !  Number of chemical species
       integer,intent(out)                                  ::  nreac   !  Number of reactions
       logical,intent(out)                                  ::  fenan   !  Enantiomers calculation flag
       logical,intent(out)                                  ::  fpermu  !  Permutations calculation flag
       logical,intent(out)                                  ::  fsoln   !  Standard state flag
!
! Local variables
!
       character(len=lenline)                               ::  line    !  
       character(len=lenline)                               ::  key     !
       integer                                              ::  io      !  Input/Output status
!
! Setting defaults  ! FLAG: missing default file names (breaks if **SYSTEM is not defined)
!
       fact   = 1.0d0
       thr    = 0.0d0
       fqvib  = 'GAS'
       ffree  = 'INDEP'
       fenan  = .FALSE.
       fpermu = .FALSE.
       fsoln  = .FALSE.
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
           case ('*MOL','*MOLECULE')
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
             call read_sys(line,'**SYSTEM',nmol,mol)      
!      
           case ('**PROP','**PROPERTIES')
!~              write(*,*) 
!~              write(*,*) 'Reading **PROPERTIES block'
!~              write(*,*) 
!
             call findline(line,'blck','**PROPERTIES')
!
             call read_prop(line,'**PROPERTIES',thr,fact,fqvib,        &
                            ffree,fenan,fpermu,fsoln)
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

!
       return
       end subroutine read_inp
!
!======================================================================!
!
       subroutine read_sys(key,blck,nmol,mol)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(nmol),intent(inout)  ::  mol   !  Molecules information
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
!~              read(uniinp,*) tgrp
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
             call read_mol(key,'*MOLECULE',nmol,mol,imol)
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
       subroutine read_mol(key,sect,nmol,mol,imol)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       type(molecule),dimension(nmol),intent(inout)  ::  mol   !  Molecules information
       character(len=lenline),intent(inout)          ::  key   !  
       character(len=*),intent(in)                   ::  sect  !  Section name
       integer,intent(in)                            ::  imol  !
       integer,intent(in)                            ::  nmol  !
!
! Local variables
!
       character(len=lenline)                        ::  line  !
       character(len=lenline)                        ::  arg   !  
       integer                                       ::  posi  !
       integer                                       ::  i     !
!
! Reading MOLECULE section keywords 
! ---------------------------------
!
       mol(imol)%nconf = 1
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
               call unkkey(key,sect)  ! FLAG: create UNKINP subroutine
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
       allocate(mol(imol)%pop(mol(imol)%nconf))
!
! Reading MOLECULE section options 
!
       if ( mol(imol)%nconf .gt. 1 ) then
         do i = 1, mol(imol)%nconf 
           write(mol(imol)%conf(i)%inp,'(I5)') i
           mol(imol)%conf(i)%inp = adjustl(mol(imol)%conf(i)%inp) 
           mol(imol)%conf(i)%inp = trim(mol(imol)%molname)//'_conf-'// &
                                     trim(mol(imol)%conf(i)%inp)//'.log'
!
           mol(imol)%conf(i)%symnum = 1  ! FLAG: change by zero to decide if it must be readen from input
         end do
       else
         mol(imol)%conf(1)%inp = trim(mol(imol)%molname)//'.log'
         mol(imol)%conf(1)%symnum = 1    ! FLAG: change by zero to decide if it must be readen from input
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
           case ('.FILE','.FILES')
!
!~              write(*,*) '    Reading .FILE option'
!~              write(*,*)
!
             do i = 1, mol(imol)%nconf        ! FLAG: check if input names are introduced correctly
               read(uniinp,*) mol(imol)%conf(i)%inp
               mol(imol)%conf(i)%inp = adjustl(mol(imol)%conf(i)%inp) 
             end do
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.SYMNUM','.SYMMETRYNUM','.SYMMETRYNUMBER')
!
!~              write(*,*) '    Reading .SYMNUM option'
!~              write(*,*)
!
             read(uniinp,*) mol(imol)%conf(:)%symnum
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
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
             allocate(mol(imol)%frag(mol(imol)%nfrag))
!
             read(key,*) mol(imol)%frag(:)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect)
         end select  
       end do
!
       return
       end subroutine read_mol
!
!======================================================================!
!
       subroutine read_prop(key,blck,thr,fact,fqvib,ffree,             &
                            fenan,fpermu,fsoln)
!
       use datatypes
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lenline),intent(inout)  ::  key     !
       character(len=*),intent(in)           ::  blck    !  Block name
       character(len=8),intent(out)          ::  fqvib   !  Qvib calculation flag
       character(len=8),intent(out)          ::  ffree   !  Free energy calculation flag
       real(kind=8),intent(out)              ::  thr     !  Threshold frequency
       real(kind=8),intent(out)              ::  fact    !  Frequencies scaling factor
       logical,intent(out)                   ::  fenan   !  Enantiomers calculation flag
       logical,intent(out)                   ::  fpermu  !  Permutations calculation flag
       logical,intent(out)                   ::  fsoln   !  Standard state flag
!
! Local variables
!
       integer                               ::  posi    !
!
! Reading PROPERTIES block sections 
! ---------------------------------
!
       fact  = 1.0d0
       thr   = 0.0d0
       fqvib = 'GAS'
!
       ffree  = 'INDEP'
       fenan  = .FALSE.
       fpermu = .FALSE.
       fsoln  = .FALSE.
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
!~              read(uniinp,*) tgrp
!~              tgrp = adjustl(tgrp)
!
             call findline(key,'blck',blck)
!
           case ('*QVIB')
!
!~              write(*,*) '  Reading *QVIB section'
!~              write(*,*)
!
             call findline(key,'sect','*QVIB')
!
             call read_qvib(key,'*QVIB',thr,fact,fqvib)
!
           case ('*FREE')
!
!~              write(*,*) '  Reading *FREE section'
!~              write(*,*)
!
             call findline(key,'sect','*FREE')
!
             call read_free(key,'*FREE',ffree,fenan,fpermu,fsoln)
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
       subroutine read_qvib(key,sect,thr,fact,fqvib)
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
       character(len=8),intent(out)          ::  fqvib  !  Qvib calculation flag
       real(kind=8),intent(out)              ::  fact   !  Frequencies scaling factor
       real(kind=8),intent(out)              ::  thr    !  Threshold frequency
!
! Local variables
!
       integer                               ::  posi   !
!
! Reading QVIB section options 
! ----------------------------
!
       fact  = 1.0d0
       thr   = 0.0d0
       fqvib = 'GAS'
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.SCALING','.FACTOR','.SCALINGFACTOR')
!
!~              write(*,*) '    Reading .SCALING option'
!~              write(*,*)
!
             read(uniinp,*) fact    ! FLAG: check reading error
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.THR','.THRESH','.THRESHOLD')
!
!~              write(*,*) '    Reading .THR option'
!~              write(*,*)
!
             read(uniinp,*) thr
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.TYPE')
!
!~              write(*,*) '    Reading .TYPE option'
!~              write(*,*)
!
             read(uniinp,*) fqvib
!
             fqvib = uppercase(fqvib)
!
             select case ( fqvib )
               case ('GAS','IDEAL','IDEALGAS')
                 fqvib = 'GAS' 
               case ('ROT','VIBROT','GRIMME2012') 
                 fqvib = 'ROT'                   
               case default        
                 call errkeyoptsect(fqvib,key,sect)  
             end select
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect)
         end select  
       end do
!
       return
       end subroutine read_qvib
!
!======================================================================!
!
       subroutine read_free(key,sect,ffree,fenan,fpermu,fsoln)
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
       character(len=8),intent(out)          ::  ffree   !  Qvib calculation flag
       logical,intent(out)                   ::  fenan   !  Enantiomers calculation flag
       logical,intent(out)                   ::  fpermu  !  Permutations calculation flag
       logical,intent(out)                   ::  fsoln   !  Standard state flag
!
! Local variables
!
       integer                               ::  posi    !
!
! Reading FREE section options 
! ----------------------------
!
       ffree  = 'INDEP'
       fenan  = .FALSE.
       fpermu = .FALSE.
       fsoln  = .FALSE.
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.TYPE')
!
!~              write(*,*) '    Reading .TYPE option'
!~              write(*,*)
!
             read(uniinp,*) ffree
!
             ffree = uppercase(ffree)
!
             select case ( ffree )
               case ('INDEP','INDEPENDENT')
                 ffree = 'INDEP' 
               case ('BOLTZ','BOLTZMANN') 
                 ffree = 'BOLTZ' 
               case ('QUAD','QUADRATURE') 
                 ffree = 'QUAD'                   
               case default        
                 call errkeyoptsect(ffree,key,sect)  
             end select
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
             call unkopt(key,sect)
         end select  
       end do
!
       return
       end subroutine read_free
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
       integer,intent(out)                                  ::  nreac  !  Number of reactions
!
! Local variables
!
       integer                                              ::  posi    !
!
! Reading REACTOR block sections 
! ------------------------------
!
       nreac = 0
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
             call read_reaction(key,'*REACTIONS',nmol,mol,nreac,reac)
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
       subroutine read_reaction(key,sect,nmol,mol,nreac,reac)
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
             call unkopt(key,sect)
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
       end module input
!
!======================================================================!