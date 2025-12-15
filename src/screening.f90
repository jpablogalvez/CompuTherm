!======================================================================!
!
       module screening
!
       implicit none
!
       private
       public  ::  screen
!
       contains
!
!======================================================================!
!
       subroutine screen(mol,nfrag,frota,rmsdmax,maemax,baemax,        &
                         idx,debug)
!
       use datatypes
       use mathtools
       use utils,      only: line_str,                                 &
                             line_int,                                 &
                             line_dp,                                  &
                             line_dvec,                                &
                             line_ivec,                                &
                             line_log,                                 &
                             print_title,                              &
                             print_titleint
!
       implicit none
!
! Input/Output variables
!
       type(molecule),intent(inout)              ::  mol       !  Molecules information
       character(len=8),intent(in)               ::  frota     !  Rotation method flag
       real(kind=8),intent(in)                   ::  rmsdmax   !  Maximum value for RMSD
       real(kind=8),intent(in)                   ::  maemax    !  Maximum value for MAE
       real(kind=8),intent(in)                   ::  baemax    !  Maximum value for BAE
       integer,dimension(mol%nconf),intent(out)  ::  idx       !  Non-redundant indexes
       integer,intent(in)                        ::  nfrag     !  Number of fragments
       logical,intent(in)                        ::  debug     !  Debug mode 
!
! Local variables
!
       real(kind=8),dimension(3,mol%nat)         ::  coord     !  Auxiliary coordinates 
       real(kind=8)                              ::  rmsd      !  Root-mean-square deviation
       real(kind=8)                              ::  mae       !  Mean absolute error
       real(kind=8)                              ::  bae       !  Biggest absolute error
       real(kind=8)                              ::  rmsdmin   !  Actual value for RMSD
       real(kind=8)                              ::  maemin    !  Actual value for MAE
       real(kind=8)                              ::  baemin    !  Actual value for BAE
       integer,dimension(nfrag)                  ::  ncoord    !  Atom chunk index
       integer,dimension(nfrag)                  ::  idxpermu  !  Permutation indexes
       integer                                   ::  nat       !  Monomer atoms
       integer                                   ::  ntot      !  Non-redundant structures
       integer                                   ::  iconf     !  Conformer index
       integer                                   ::  iaux      !  Auxiliary integer number
       integer                                   ::  i,j,k     !  Indexes
       logical,dimension(mol%nconf)              ::  new       !
       logical                                   ::  chk       !  Checking variable
!
! Screening conformers of the input molecule
! 
       call print_title(6,5,'Screening of the conformations','.')
       write(*,*)
!
       nat = mol%nat/nfrag
!
       new(:) = .TRUE.     
!
       do iconf = 1, mol%nconf
         idx(iconf) = iconf
       end do
!
       write(*,'(1X,A)') 'Screening input structures'
       write(*,'(1X,A)') 'Please wait, this may take a while...'
       write(*,*)
!
! Setting up the initial chunk of atoms
!
       ncoord(1) = 0
       do i = 2, nfrag
         ncoord(i) = ncoord(i-1) + nat
       end do
!
! Screening input structures
!
       do i = 1, mol%nconf-1
         do j = i+1, mol%nconf 
!
           do k = 1, nfrag
             idxpermu(k) = k
           end do
!
           rmsdmin = rmsdmax
           maemin  = maemax
           baemin  = baemax
!
           chk = .FALSE.
!
           do 
             do k = 1, nfrag
               coord(:,(k-1)*nat+1:k*nat) = mol%conf(i)%               &
                  coord(:,ncoord(idxpermu(k))+1:ncoord(idxpermu(k))+nat)
             end do
!
             call compare(frota,mol%nat,mol%conf(j)%coord,coord,       &
                          rmsd,mae,bae)
!
!~ write(*,*) trim(mol%conf(i)%inp),' and ',trim(mol%conf(j)%inp),rmsd,mae,bae
!
             if ( (rmsd.le.rmsdmax) .and. (mae.le.maemax) .and.      &  ! FLAG: check how to update values
                                                (bae.le.baemax) ) then
!~                if ( (rmsd.lt.rmsdmin) .or. (mae.lt.maemin) .or.      & 
!~                                                 (bae.lt.baemin) ) then
!
!~              if ( rmsd .le. rmsdmax ) then
!
               if ( rmsd .lt. rmsdmin ) then
                 rmsdmin = rmsd
                 maemin  = mae
                 baemin  = bae
               end if
!
               new(j) = .FALSE.
               chk    = .TRUE.
!
             end if
!
             if ( nfrag .gt. 1 ) then
               if ( .NOT. nextp(nfrag,idxpermu) ) exit
             else
               exit
             end if
           end do
!
           if ( chk .and. debug ) then
             write(*,'(1X,2(1X,A,1X,I4))') 'Structure number',j,       &
                                  'is redundant with structure number',i
             write(*,'(2X,A)') 'Structure '//trim(mol%conf(j)%inp)//   &
                               ' is equivalent to structure '//        &
                                                   trim(mol%conf(i)%inp)
             write(*,'(1X,3(1X,A,1X,D10.4))') 'RMSD=',rmsdmin,'MAE=', &
                                                 maemin,'BAE=',baemin
             write(*,'(2X,A,1X,I4,1X,A)') 'Removing structure number', &
                                                        j,'from the set'
             write(*,*)
           end if
!
         end do
       end do
!
! Removing redundant structures from the set of conformers
!
       ntot = 0
       do i = 1, mol%nconf
         if ( new(i) ) ntot = ntot + 1
       end do
!
       write(*,'(1X,A,1X,I4)') 'Number of non-redundant structures =', &
                                                                    ntot
       write(*,*)
!
       ntot  = mol%nconf
       iconf = 1
!
       do while ( iconf .le. ntot )
!
         if ( .NOT. new(iconf) ) then
           do while ( iconf .le. ntot )
!
             if ( new(ntot) ) then
!
               iaux       = idx(iconf) 
               idx(iconf) = idx(ntot)
               idx(ntot)  = iaux
!
               ntot = ntot - 1
!
               exit
!
             end if
!
             ntot = ntot - 1
!
           end do
         end if
!
         iconf = iconf + 1
!
       end do
!
       mol%nconf = ntot
!
! Saving list with non-redundant structures
!
!~        if ( debug ) then
!~          open(unit=uniout,file=trim(outp),action='write')
!~ !
!~          do iconf = 1, mol%nconf
!~ !           write(*,'(4X,A,1X,I4,2X,A,2X,A)') 'Structure',idx(iconf),'-', &
!~ !                                          trim(mol%conf(idx(iconf))%inp)
!~            write(uniout,*) trim(mol%conf(idx(iconf))%inp)
!~          end do
!~ !         write(*,*)
!~ !
!~          close(uniout)
!~        end if  
!
       end subroutine screen
!
!======================================================================!
!
! COMPARE - COMPARE structures
!
! This function returns
!
       subroutine compare(frota,nat,inp,ref,rmsd,mae,bae)
!
       use geometry
       use superpose
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  inp    !  Input coordinates
       real(kind=8),dimension(3,nat),intent(in)  ::  ref    !  Reference coordinates
       real(kind=8),intent(out)                  ::  rmsd   !  Root-mean-square deviation
       real(kind=8),intent(out)                  ::  mae    !  Mean absolute error
       real(kind=8),intent(out)                  ::  bae    !  Biggest absolute error
       character(len=8),intent(in)               ::  frota  !  Rotation method flag
       integer,intent(in)                        ::  nat    !  Number of particles
!
! Local variables
!
       real(kind=8),dimension(3,nat)             ::  aux    !  Auxiliry coordinates
       real(kind=8)                              ::  drmsd  !  Auxiliary RMSD value   
       real(kind=8)                              ::  dmae   !  Auxiliary MAE value   
       real(kind=8)                              ::  dbae   !  Auxiliary BAE value   
!
! Computing RMSD between input structures
!
       aux(:,:) = inp(:,:)
!
       select case ( frota )
         case ('KABSCH','KAB','ROTKAB','ROTKABSCH')
!
           call dsuperpose('KABSCH',3,nat,aux,ref)
!
           call dmcalcrmsd(3,nat,aux,ref,rmsd)
           call dmcalcmae(3,nat,aux,ref,mae,bae)
!
           call Sxy(nat,inp,aux) 
!
           call dsuperpose('KABSCH',3,nat,aux,ref)
!
           call dmcalcrmsd(3,nat,aux,ref,drmsd)
           call dmcalcmae(3,nat,aux,ref,dmae,dbae)
!
           if ( drmsd .lt. rmsd ) then
             rmsd = drmsd
             mae  = dmae
             bae  = dbae
           end if
!
         case ('QUAT1','KK','QUATKK')
!
           call dsuperpose('QUATKK',3,nat,aux,ref)
!
           call dmcalcrmsd(3,nat,aux,ref,rmsd)
           call dmcalcmae(3,nat,aux,ref,mae,bae)
!
           call Sxy(nat,inp,aux) 
!
           call dsuperpose('QUATKK',3,nat,aux,ref)
!
           call dmcalcrmsd(3,nat,aux,ref,drmsd)
           call dmcalcmae(3,nat,aux,ref,dmae,dbae)
!
           if ( drmsd .lt. rmsd ) then
             rmsd = drmsd
             mae = dmae
             bae = dbae
           end if
!
         case ('QUAT2','CTK','QUATCTK')
!
           call dsuperpose('QUATCTK',3,nat,aux,ref)
!
           call dmcalcrmsd(3,nat,aux,ref,rmsd)
           call dmcalcmae(3,nat,aux,ref,mae,bae)
!
       end select
!
       return
       end subroutine compare
!
!======================================================================!
!
       end module screening
!
!======================================================================!
