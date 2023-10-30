program vib_potPGR


!     This program reads an input file for the ground state calculation of Quantum Espresso 
!     (e.g. C10H16_scf.in) which contains a given set of nuclear coordinates, and it generates
!     other QE input files which are identical to the original one except for the fact that 
!     their nuclear positions are equal to those of the original QE input file plus
!     one given eigenvector (U) of the dynamical equation. These U vectors are stored in a
!     file produced by QE (e.g. called C10H16.dynG-cutoff30). The produced files are stored in 
!     the 'frozen+R' folder.
!     This program also produces a file (wQ0.dat) with the phonon frequencies, that is read by Maat.
!
!     A possible example of input file vib_potPGR.in is: 
!
!-------------------------------------------
!
! QEMassesFile = out_scf.out
! InputQEscffile = scf.in
! InputQEpotfile = Veff-unperturb.in 
! FreqsUfile = dynmat.out
! #FreqsUfile = dynG
! BasicDirPath = ./
! Files_for_2nd_derivatives = 1
! DisplacementParameter = 4.000000d0
! 
! %Atoms=
! C  12.0107      10
! H  1.007825035  16
! %
!
!-------------------------------------------
!   
!     Author: Pablo Garcia Risueno (2016-2017).



  implicit none
  logical :: exist
  integer :: i, j, k, l, h, i1, i2, i3, FirstCalculatedMode, LastCalculatedMode, Nunphysfreqs
  integer :: total_num_atoms, min_band_index, max_band_index, dimnu, Ne, Nks, filesfor2ndderivatives 
  real(8) :: auxr1, auxr2, auxr3, DisplParam
  complex(8) :: auxc1, auxc2, auxc3
  real(8),allocatable :: auxrv1(:), xi(:,:), U(:,:), omega(:), omegaorig(:), epsilon(:), mass(:,:), atom_mass(:), atom_posi(:,:)
  character(2), allocatable :: atom_name(:)
  integer, allocatable :: atom_num(:)
  complex(8),allocatable :: auxcv1(:), gDW(:,:,:,:)
  character(100) :: InputQEscf_file_name, InputQEpot_file_name, FreqsU_file_name, prefix, QE_masses_file_name,filename   
  character(250) :: basic_dir_path 

   !! local varibales
   integer  :: ii, jj, kk, nn
   integer  :: Nspecies, n_mode
   integer  :: status, ierror ! I/O status
   logical  :: exist_flag, positivedisplacement



  call read_Nspecies (Nspecies)

  allocate(atom_mass(Nspecies),atom_num(Nspecies),atom_name(Nspecies))
  
  call read_input_file (Nspecies, Nunphysfreqs, basic_dir_path, InputQEscf_file_name, InputQEpot_file_name, &
       & FreqsU_file_name, QE_masses_file_name, prefix, atom_name, atom_mass, atom_num, dimnu, total_num_atoms, &
       & FirstCalculatedMode, LastCalculatedMode, filesfor2ndderivatives, DisplParam)  
 
  allocate(atom_posi(total_num_atoms,3),omegaorig(dimnu),xi(dimnu,dimnu),U(dimnu,dimnu),mass(total_num_atoms,2))  
  
  call read_QE_masses (total_num_atoms,  Nspecies, QE_masses_file_name, mass)   
  
  call read_nuclear_positions (total_num_atoms, InputQEscf_file_name, atom_posi)  ! The nuclear positions (atom_posi) are in Angstrom                

                                                                                          
  call read_omega_U (total_num_atoms,  Nunphysfreqs, dimnu, FreqsU_file_name, mass, omegaorig, xi, U)   

! Here U is in a.u.

  call  write_new_QE_input_single_displ (.true., total_num_atoms, dimnu, DisplParam, basic_dir_path, & 
     & InputQEscf_file_name, InputQEpot_file_name, prefix, atom_posi, U)
  call  write_new_QE_input_single_displ (.false., total_num_atoms, dimnu, DisplParam, basic_dir_path, & 
     & InputQEscf_file_name, InputQEpot_file_name, prefix, atom_posi, U)
  if ( filesfor2ndderivatives .gt. 0) then
    call  write_new_QE_input_double_displ (total_num_atoms, dimnu, DisplParam, basic_dir_path, & 
       & InputQEscf_file_name, InputQEpot_file_name, prefix, atom_posi, U)
  end if     

! Here U is in a.u.
  
  filename='Input_dot'
  call write_U_files (filename,dimnu,total_num_atoms,DisplParam,omegaorig,U)
  filename='Uvec.dat'
  call write_U_files (filename,dimnu,total_num_atoms,1.d0,omegaorig,U)
  
  write(*,*) '             ------ The calculations finished satisfactorily. ------'; write(*,*)

  !! deallocate the arrays
  deallocate(atom_posi,xi,U,omegaorig,atom_name,atom_num,atom_mass)

end program ! End of the main program


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------





! PENG'S CODE
! 
!    
!    read(90,*) Input_dot      !! The input file for dot eignemodes
!    read(90,*) Input_coord    !! The input file for the coorditotal_num_atomses
!    read(90,*) Nspecies          !! number of elements
!    allocate(atom_num(Nspecies))
!    allocate(atom_name(Nspecies))
!    allocate(atom_mass(Nspecies))
!    do ii = 1, Nspecies
!    read(90,*) atom_name(ii)  !! Read the name of element
!    read(90,*) atom_num(ii)  !! Read the number of each type
!    read(90,*) atom_mass(ii)  !! Read the mass of each atom
!    end do
!    read(90,*) n_mode         !! The phonon mode need to be calculated 
!    
!    close(90)
! 
! 
!    write(*,*) 'Input read from file:'
!    write(*,*) '---------------------'
!    write(*,*) 'Input coord file    :', Input_coord
!    write(*,*) 'Input dot eigenmode :', Input_dot
!    write(*,*) 'Type of atoms       :', Nspecies
!    do ii = 1, Nspecies
!    write(*,*) 'Atom type No.       :', ii
!    write(*,*) 'Name of the atom    :', atom_name(ii)
!    write(*,*) 'Number of atom      :', atom_num(ii)
!    write(*,*) 'Mass of the atom    :', atom_mass(ii)
!    end do
!    write(*,*) 'Phonon mode to be calculated:', n_mode
! 
!    total_num_atoms = 0
!    do ii = 1, Nspecies
!       total_num_atoms = total_num_atoms+atom_num(ii)
!       atom_mass(ii) = atom_mass(ii)*P_mass2au
!    end do
! 
!    !! allocate the arrays
!    allocate(atom_posi(total_num_atoms,3))
!    allocate(vib_mode(3*total_num_atoms,3*total_num_atoms))
!    allocate(omegaorig(3*total_num_atoms))
!    
!    !! Read the vibrational modes
!    call dot_eigen(Input_coord, Input_dot, total_num_atoms, atom_posi, vib_mode, omegaorig)
! 
!    !! Calculate the displacement 
!    do ii = 1, 3*total_num_atoms !! For the modes
!       nn = 1
!       do jj = 1, Nspecies
!          do kk = 1, 3*atom_num(jj)
!             vib_mode(ii,nn) = vib_mode(ii,nn)/sqrt(atom_mass(jj))
!             nn = nn+1
!          end do
!       end do
!    end do 
! 
!     open(90, file='frequency.dat')
!     do ii = 1, 3*total_num_atoms
!        write(90,*) 'I=', ii, 'Freq=', omegaorig(ii)
!     end do
!     close(90)
!       
!    !! Transfer to au2A
!    do ii = 1, 3*total_num_atoms
!       vib_mode(ii,:) = vib_mode(ii,:)* P_au2A
!    end do

   !! Output the results
 !  call output_file(total_num_atoms, Nspecies, atom_name, atom_num, vib_mode, atom_posi, n_mode)
   
   
   

 subroutine StripSpaces(string)
  
  ! This subroutine removes blank spaces (' ') from the name of the files
  
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

  end subroutine
  

!-------------------------------------------------------------------------------------------

       
subroutine read_omega_U (total_num_atoms, Nunphysfreqs, dimnu, FreqsU_file_name, mass, omegaorig, xi, U)     
 
 ! This subroutine reads the file which contains the output of QE and gives phonon frequencies and
 ! eigenvectors (xi and U) of the dynamical equation. We assume that the negligible acoustic modes 
 ! correspond to the lowest six freqs.
 ! Note that the output of CPMD gives the xi eigenvectors, while the output of Quantum Espresso
 ! gives the U eigenvectors (U := xi/sqrt(M) ). The original code written by Peng used the U eigenvectors; it read the
 ! xi eigenvectors of CPMD and then it divided them by sqrt(M_I). We don't need to do it, because
 ! we know the U eigenvectors directly from the output of Quantum Espresso. However, we need to normalize the xi 
 ! vectors that we obtain, which changes the value of the U's; this is, the U's that we use in our code
 ! are NOT identical to those given by Quantum Espresso, but scaled ones (each U from us is proportional to the corresponding
 ! U from QE, but the proportionality constant is different for every nu-index)
 

  ! System parameters
  integer, Intent(In) :: total_num_atoms                ! The number of atoms of the system
  integer, Intent(In) :: Nunphysfreqs                   ! The number of unphysical freqs of the system (6 in isolated system, 3 in periodic system)
  integer, Intent(In) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-Nunphysfreqs)
  character(100), Intent(In) :: FreqsU_file_name        ! Name of the file where the phonon frequencies and eigenvectors of the dynamical matrix are stored
  real(8), Intent(In) :: mass(total_num_atoms,2)        ! Masses of atoms
  real(8), Intent(Out) :: omegaorig(dimnu)              ! Phonon frequencies (as they are read from file, with the units provided by the original file)
  real(8), Intent(Out) :: xi(dimnu,dimnu)				! Eigenvectors of the dynamical matrix
  real(8), Intent(Out) :: U(dimnu,dimnu)				! Eigenvectors of the dynamical matrix divided by the square root of the corresponding mass
 
 
   ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT, auxchar1, auxchar2
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist, exist_flag
  integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, atomindex
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, auxr10, P_au2A

 
  write(*,*) '  *** Now reading the output of QE ('&
              & , trim(FreqsU_file_name),') to find ' 
  write(*,*)  '             phonon frequencies and eigenvectors of the dynamical matrix. ***'
     


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!     BLOCK FOR READING PHONON FREQUENCIES AND EIGENVECTORS OF THE DYNAMICAL MATRIX        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  do i=1,total_num_atoms
    if ( mass(i,2) .lt. 0.d0 ) then
      write(*,*)
      write(*,*) ' ERROR: Masses were not properly read (mass of atom ',i,' not read).'
      write(*,*)
      call exit(1)
    end if
  end  do
 
  
 omegaorig(:) = 0.d0
 xi(:,:) = 0.d0
 U(:,:) = 0.d0
 freqcount = 0
 
 open(341,file=FreqsU_file_name,status='unknown')
  
 do m=1,dimnu
  
  do k=1,10*((dimnu+20)**2)

      read(341,'(A)',IOSTAT=stat) inputline
      !write(*,*) inputline
      if (stat .ne. 0) go to 1239

      auxstr1=''; auxstr2=''; auxstr3=''
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( auxstr1 .eq. 'freq') then    ! Example line of the file to read:      freq (    1) =      -4.920305 [THz] =    -164.123717 [cm-1]
          freqcount = freqcount + 1
         ! write(*,*) 'CONSEGUIDO, linea', k, i
          l1 = k       ! l1 is the line where the word 'Frequencies' is found
          go to 1239
        end if
      end do ! i
      
  end do !k   
   
 1239 continue      
      
      
   punto=100
   j=i+1
   do while ( (j .le. 90) )
      ! write(*,*) j, punto
      letter = inputline(j:j)
      if (letter .eq. '=') then
        punto=j
        ! go to 
      end if
      auxstr2 = trim(auxstr2)//letter
      j=j+1
   end do  
   
    
    j = punto+1
    do while ( (j .le. 90) .and. (letter.ne.'[') ) !.and. (j.le.punto+6))
      ! write(*,*) j, punto
      letter = inputline(j:j)
      if (letter .ne. '[') then
        punto=j
        auxstr3 = trim(auxstr3)//letter 
      end if
      j=j+1
   end do   
    
    
   
    ! write(*,*) '    ',trim(auxstr3)
    read(auxstr3,*) auxr1
    omegaorig(freqcount)=auxr1
    
    do i=1,total_num_atoms
      read(341,*) auxchar1, auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxchar2  ! Example line of the file to read:  (  0.000000   0.000000    -0.104934   0.000000     0.104934   0.000000   )
      U((i-1)*3+1,freqcount) = auxr1
      U((i-1)*3+2,freqcount) = auxr3
      U((i-1)*3+3,freqcount) = auxr5
      ! write(*,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6
    end do    

 end do ! m      
  
     
  

!   Input for Maat:  * wQ0.dat: 1st row: number of q-points and the number of phonon frequencies; then: q-index, branch-index, phonon frequency  
!    open(343,file='wQ0.dat',status='unknown')
!    write(343,*) 1, dimnu
!    do m=1,dimnu
!       write(343,*) 1,m,omegaorig(m)
!    end do
!    close(343)
!    
!    ! File defined by Peng
!    open(90, file='frequency.dat')
!    do ii = 1, 3*total_num_atoms
!        write(90,*) 'I=', ii, 'Freq=', omegaorig(ii)
!    end do
!    close(90)   
   
    
   ! We convert the U vectors, given by the QE output, to xi vectors (orthogonal)
   do j=1,dimnu
     do i=1,dimnu
       atomindex = ceiling( dble(i)/3.d0 -0.0001d0 )  ! OBoTP nint( dble(i)/3.d0 ) gives error (e.g. gives 0 for the 1st index)
       xi(i,j) = U(i,j) * sqrt(mass(atomindex,2))
     end do
   end do  
   



   ! We normalize the xi vectors (the fact that the output of QE gives normalized U does not mean that the corresponding xi are normalized)
   ! The orthoNORMALITY of xi is necessary to apply the expression of gDW as a function of gFan (see Risueno/Draxl/Bester 2017)
   do j=1,dimnu 
     auxr1=0.d0
     do i=1,dimnu
       auxr1 = auxr1 + (xi(i,j))**2 
     end do
     auxr1 = sqrt(auxr1)
     do i=1,dimnu
       xi(i,j) = xi(i,j)/auxr1 
     end do
   end do  

      
 
    ! We redefine the U vectors with the new definition of xi
   do j=1,dimnu
     do i=1,dimnu
       atomindex = ceiling( dble(i)/3.d0 -0.0001d0 )
       U(i,j) = xi(i,j) / sqrt(mass(atomindex,2))
     end do
   end do  


! COMMENT (PGR): The eigenvectors of the dynamical matrix do NOT have determined units because
! by definition eigenvectors are always undetermined up to a multiplicative factor. However,
! the eq. (38) of Ponce PRB 90, 214304 (2014) states that if mass is in atomic units
! and \xi, U satisfy the corresponding orthoNORMALITY condition, then \xi and U are in atomic units.
! Since we have just ensured that this orthonormality hols, then in this point of the program
! both \xi and U ARE IN ATOMIC UNITS.

 
 
!  write(*,*) ' CHECKING THE ORTHONORMALITY OF EIGENVECTORS OF THE DYNAMICAL MATRIX'

  auxr4=0.d0  ! summation of the orthogonal terms
  auxr5=0.d0  ! summation of terms equaling 1
  auxr6=0.d0   ! maximum error
  
  do j=Nunphysfreqs+1,dimnu
    do k=Nunphysfreqs+1,dimnu
       
       auxr1=0.d0
       do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
          auxr1=auxr1 + xi(i,j)*xi(i,k)
       end do
       
       if (j.ne.k) then
          auxr1 = abs(auxr1)
          auxr4=auxr4+auxr1
          if (auxr1 .gt. 0.0001d0) then
             write(*,'(A,I4,A,I4,A,f15.8,A)') ' Prod eigv. ',j,' and ',k, '= ',auxr1,' NOT orthogonal'
          end if
       end if
       
       if (j.eq.k) then
         auxr5=auxr5+abs(auxr1)
         ! write(*,*) '  Modulus of vec ',j,' = ',sqrt(auxr1)
       end if
       if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
    
    end do
  end do  
 

  
  auxr4=auxr4/dble((dimnu-Nunphysfreqs)**2 - dimnu+Nunphysfreqs); auxr5=auxr5/dble(dimnu-Nunphysfreqs)
  if ( ( abs(auxr4) .lt. 0.00001d0) .and. (abs(auxr5-1.d0).lt.0.00001) ) then 
    ! write(*,*) ;  write(*,*) '   ORTHONORMALITY  1: OK '
!write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
!write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5 
  else
     write(*,*) '   WARNING: The first condition of orthonormality is not properly satisfied' 
     write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
     write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5     
  end if


  auxr4=0.d0  ! summation of the orthogonal terms
  auxr5=0.d0  ! summation of terms equaling 1
  auxr6=0.d0   ! maximum error
  
  do j=Nunphysfreqs+1,dimnu
    do k=Nunphysfreqs+1,dimnu
       
       auxr1=0.d0
       do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
          auxr1=auxr1 + xi(j,i)*xi(k,i)
       end do
       if (j.ne.k) auxr4=auxr4+abs(auxr1)
       if (j.eq.k) auxr5=auxr5+abs(auxr1)
       if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
    
    end do
  end do  
  auxr4=auxr4/dble((dimnu-Nunphysfreqs)**2 - dimnu+Nunphysfreqs); auxr5=auxr5/dble(dimnu-Nunphysfreqs)
 
  if ( ( abs(auxr4) .lt. 0.00001d0) .and. ( abs(auxr5-1.d0).lt.0.00001) ) then 
   ! write(*,*) '   ORTHONORMALITY  2: OK '
!write(*,*) ' The average orthonormaliz 2 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
!write(*,*) ' The average orthonormaliz 2 is (parallel terms):', auxr5             
  else
     write(*,*) '   WARNING: The second condition of orthonormality is not properly satisfied'  
     write(*,*) ' The average orthonormaliz 2 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
     write(*,*) ' The average orthonormaliz 2 is (parallel terms):', auxr5   
  end if


  auxr4=0.d0  ! summation of the orthogonal terms
  auxr5=0.d0  ! summation of terms equaling 1
  auxr6=0.d0   ! maximum error
  
  
  do j=1,dimnu
    do k=1,dimnu
       
       auxr1=0.d0
       do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
          atomindex = ceiling( dble(i)/3.d0 -0.0001d0 )
          auxr1=auxr1 + U(i,j)*U(i,k)*(mass(atomindex,2))
       end do
       if (j.ne.k) auxr4=auxr4+abs(auxr1)     
       if (j.eq.k) auxr5=auxr5+abs(auxr1-1.d0)
       if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
    
   !       if (j.eq.k) then 
!           write(*,*) 'i=',j, 'j=',k,'cal=',auxr1, 'auxr5=',auxr5
!         end if  
!         if ( (auxr1 .gt. 0.05) .and. (j .ne. k)) stop
!        ! IF ( (J.EQ.K) .and. (abs(auxr1-1.d0) .gt. 0.0005 ) ) stop 
    
    end do
  end do  
  
  auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+Nunphysfreqs)); auxr5=auxr5/(dble(dimnu-n_mode_beg+1))
  if ( ( abs(auxr4) .lt. 0.00001d0) .and. (abs(auxr5).lt.0.0001d0) ) then 
!     write(*,*);  write(*,*) '  U ORTHONORMALITY   1: OK '
!     write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
!     write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5        
  else
     write(*,*) '  U ERROR: The first condition of orthonormality is not properly satisfied' 
     write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
     write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5   
     stop  
  end if
     
     
     
     
     

! WE DON'T WRITE U in Angstrom (hence we comment the text below), that is done later, within write_new_QE_input_single_displ.
!! Transfer atomic units to Angstrom (as done by Peng); this is necessary because the output of QE is expected to be given in Angstrom
!    P_au2A = 0.5291772d0
!    do ii = 1, 3*total_num_atoms
!       U(ii,:) = U(ii,:)* P_au2A       ! PENG:  vib_mode(ii,:) = vib_mode(ii,:)* P_au2A
!    end do


! COMMENT (PGR): In this point of the program, U is in Angstrom but \xi is in atomic units.



!  PENG'S ORIGINAL CODE (it reads xi and it reads U from xi; we don't need it, because QE gives U directly): 
!
!    !! Read the vibrational modes
!    call dot_eigen(Input_coord, Input_dot, total_num_atoms, atom_posi, vib_mode, omegaorig)
! 
!    !! Calculate the displacement 
!    do ii = 1, 3*total_num_atoms !! For the modes
!       nn = 1
!       do jj = 1, Nspecies
!          do kk = 1, 3*atom_num(jj)
!             vib_mode(ii,nn) = vib_mode(ii,nn)/sqrt(atom_mass(jj))
!             nn = nn+1
!          end do
!       end do
!    end do 
! 
!     open(90, file='frequency.dat')
!     do ii = 1, 3*total_num_atoms
!        write(90,*) 'I=', ii, 'Freq=', omegaorig(ii)
!     end do
!     close(90)
!       
!    !! Transfer to au2A
!    do ii = 1, 3*total_num_atoms
!       vib_mode(ii,:) = vib_mode(ii,:)* P_au2A
!    end do
! 
!   


end subroutine read_omega_U




!-------------------------------------------------------------------------------------------


subroutine read_Nspecies (Nspecies)

  ! This subroutine reads the number of species of the system to analyse. They can be read either from
  ! the variable Nspecies or automatically from the block starting with the keyword "%atoms=" and ending
  ! with the character "%"

  integer, Intent(Out) :: Nspecies
   
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,iostat,unitt !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist, exist_flag
  real(8) :: auxr1, auxr2

 
   inquire(file = 'vib_potPGR.in', exist=exist_flag)
   if(exist_flag) then
      open(344, file='vib_potPGR.in')
   else
      write(*,*)
      write(*,'(a)') "  ERROR: An input file with the parameters (vib_potPGR.in) is needed. "
      WRITE(*,*)
      call exit(1)
   endif
   
   
   
   Nspecies = -1
    
   do k=1,500

      read(344,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 1137

      auxstr1=''
      auxstr2=''
      auxstr3=''
   
      do i=1,100
        letter = inputline(i:i)
        if (letter .ne. '=') then
          auxstr1 = trim(auxstr1)//letter 
        else 
          go to 1235
        end if
      end do ! i
 
      1235 continue
  
      do j=i+1,100
        letter = inputline(j:j)
        auxstr2 = trim(auxstr2)//letter 
      end do ! j

          
    if ((auxstr1 .eq. 'SpeciesNumber').or.(auxstr1 .eq. 'speciesnumber'))  then
       Read(auxstr2, '(I4)' ) Nspecies
       close(344)
       return
    end if
     
     
    if ((auxstr1 .eq. '%Atoms').or.(auxstr1 .eq. '%atoms'))  then
        
       do while (letter .ne. '%') 
       
          read(344,'(A)',IOSTAT=stat) letter
          Nspecies=Nspecies+1
          
          if (letter .eq. '%') then
            close(344)
            return
          end if
          
        end do
        
    end if

  end do ! k
   

1137 continue

 close(344)

end subroutine read_Nspecies



!-------------------------------------------------------------------------------------------




subroutine read_input_file (Nspecies, Nunphysfreqs, basic_dir_path, InputQEscf_file_name, InputQEpot_file_name,  &
                   & FreqsU_file_name,  QE_masses_file_name, prefix, atom_name,atom_mass,atom_num, &
                   & dimnu, total_num_atoms, FirstCalculatedMode, LastCalculatedMode, filesfor2ndderivatives , DisplParam)

 ! This subroutine reads the inp.in file, which contains the names of the files where one must read the data
 ! as well as other relevant data. 

  ! System parameters
   Integer, Intent(In) :: Nspecies
   Integer, Intent(Out) :: Nunphysfreqs
   character(250), Intent(Out) :: basic_dir_path         ! Path of the basic directory where the calc_gFan.x program is called
   character(100), Intent(Out) :: InputQEscf_file_name   ! Name of the file where the nuclear positions and other infos are stored
   character(100), Intent(Out) :: InputQEpot_file_name   ! Name of the file where the information for printing the potential
   character(100), Intent(Out) :: FreqsU_file_name      ! Name of the file where the phonon frequencies and eigenvectors of the dynamical matrix are stored
   character(100), Intent(Out) :: QE_masses_file_name    ! Output file of an 'scf' calculation of Quantum Espresso. It contains the masses of the atoms.
   character(100), Intent(Out) :: prefix                ! basis of the QE prefix variable
   Character(2), Intent(Out) :: atom_name(Nspecies)      ! Atomic species name
   Real(8), Intent(Out) :: atom_mass(Nspecies)           ! Atomic mass
   Integer, Intent(Out) :: atom_num(Nspecies)            ! Numerical index for atoms
   Integer, Intent(Out) :: dimnu                         ! Total number of phonon frequencies (and modes)
   Integer, Intent(Out) :: total_num_atoms               ! Total number of atoms
   Integer, Intent(Out) :: FirstCalculatedMode, LastCalculatedMode ! First and last modes whose corresponding QE input files are generated
   Integer, Intent(Out) :: filesfor2ndderivatives        ! If non-zero, the files for calculations of 2nd derivatives (by finite-difference) are generated in the frozenRR folder
   Real(8), Intent(Out) :: DisplParam                    ! Parameter for the displacement; the nuclear positions of the newly generated QE inputs are the original ones plus the ( U^{\nu}_I vectors multiplied by DisplParam) 
   
 
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,iostat,unitt,PeriodicSystem  !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: auxr1, auxr2, P_mass2au


  P_mass2au = 1822.8884842645450d0 !! mass -> a.u.
  filesfor2ndderivatives = 0


  ! real(dp), allocatable :: atom_posi(:,:)
  ! real(dp), allocatable :: vib_mode(:,:) --> xi(:,:)
  ! real(dp), allocatable :: vib_freq(:) --> omegaorig(:)
  

   !parameter (P_mass2au = 1.660538782E-27_DP/9.10938215E-31_DP)  !! mass -> a.u.
  ! parameter (P_au2A = 0.5291772_dp)


  total_num_atoms = 0
  FirstCalculatedMode=1
  LastCalculatedMode=-1
  periodicsystem = 0; Nunphysfreqs=6
  FirstDerivParam=1.d0
  InputQEpot_file_name = 'Veff-unperturb.in'
  InputQEscf_file_name = 'scf.in'

  ! Reading parameters from file
  


   
  
   write(*,*); write(*,*); 
   write(*,*) '  ****************************************************************************' 
   write(*,*) '  ************************  NOW RUNNING vib_potPGR  **************************'   
   write(*,*) '  ****************************************************************************'   
   write(*,*) 
   write(*,*)
   write(*,*) '  ****** FROM THE INPUT FILE OF vib_potPGR (vib_potPGR.in): **************'   
   
   
open(345, file='vib_potPGR.in',status='unknown')

   do k=1,50

      read(345,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
  
      if (stat .ne. 0) go to 1117


      auxstr1=''
      auxstr2=''
      auxstr3=''
   
      do i=1,100
        letter = inputline(i:i)
        if (letter .ne. '=') then
          auxstr1 = trim(auxstr1)//letter 
        else 
          go to 1235
        end if
      end do ! i
 
      1235 continue

  
      do j=i+1,100
        letter = inputline(j:j)
        auxstr2 = trim(auxstr2)//letter 
      end do ! j


     if ((auxstr1 .eq. 'BasicDirPath').or.(auxstr1 .eq. 'Basicdirpath').or. &
        & (auxstr1 .eq. 'basicdirpath'))   then
         basic_dir_path = trim(auxstr2)
         call StripSpaces (basic_dir_path)
     end if    

     if ((auxstr1 .eq. 'InputQEscffile').or.(auxstr1 .eq. 'inputQEscffile').or.(auxstr1 .eq. 'inputqescffile'))   then
         InputQEscf_file_name = trim(auxstr2)
         call StripSpaces (InputQEscf_file_name)
     end if
     
     if ((auxstr1 .eq. 'InputQEpotfile').or.(auxstr1 .eq. 'inputQEpotfile').or.(auxstr1 .eq. 'inputqepotfile'))   then
         InputQEpot_file_name = trim(auxstr2)
         call StripSpaces (InputQEpot_file_name)
     end if     

     if ((auxstr1 .eq. 'QEmassesfile').or.(auxstr1 .eq. 'QEMassesFile').or.(auxstr1 .eq. 'qemassesfile'))   then
         QE_masses_file_name = trim(auxstr2)
         call StripSpaces (QE_masses_file_name)
     end if    


     if ((auxstr1 .eq. 'Prefix').or.(auxstr1 .eq. 'prefix'))   then
         prefix = trim(auxstr2)
         call StripSpaces (prefix)
     end if
     
     if ((auxstr1 .eq. 'FreqsUfile').or.(auxstr1 .eq. 'freqsUfile').or.(auxstr1 .eq. 'freqsufile'))   then
         FreqsU_file_name = trim(auxstr2)
         call StripSpaces (FreqsU_file_name)
     end if
          
     
     if ((auxstr1 .eq. '%Atoms').or.(auxstr1 .eq. '%atoms')) then
       do j=1,Nspecies
          read(345,*) atom_name(j),atom_mass(j),atom_num(j)
          write(*,*) atom_name(j),atom_mass(j),atom_num(j)
       end do  
     end if
     
     if ((auxstr1 .eq. 'PeriodicSystem').or.(auxstr1 .eq. 'Periodicsystem') .or. &
         & (auxstr1 .eq. 'periodicsystem') )  Read(auxstr2, '(I10)' ) periodicsystem
     
     if ((auxstr1 .eq. 'FirstCalculatedMode').or.(auxstr1 .eq. 'Firstcalculatedmode') .or. &
         & (auxstr1 .eq. 'firstcalculatedmode') )  Read(auxstr2, '(I10)' ) FirstCalculatedMode

     if ((auxstr1 .eq. 'LastCalculatedMode').or.(auxstr1 .eq. 'Lastcalculatedmode') .or. &
         & (auxstr1 .eq. 'lastcalculatedmode') )  Read(auxstr2, '(I10)' ) LastCalculatedMode

     if ((auxstr1 .eq. 'DisplacementParameter').or.(auxstr1 .eq. 'Displacementparameter') .or. &
         & (auxstr1 .eq. 'displacementparameter') )  Read(auxstr2, '(G15.8)' ) DisplParam

     if ((auxstr1 .eq. 'Files_for_2nd_derivatives').or.(auxstr1 .eq. 'Files_For_2nd_Derivatives') .or. &
         & (auxstr1 .eq. 'files_for_2nd_derivatives') )  Read(auxstr2, '(I10)' ) filesfor2ndderivatives

  end do ! k
   
1117 continue



 total_num_atoms = 0
 do j=1,Nspecies
   total_num_atoms  = total_num_atoms + atom_num(j)
   atom_mass(j) = atom_mass(j)*P_mass2au
 end do
 dimnu = 3*total_num_atoms 
 
 if (LastCalculatedMode .lt. 0) LastCalculatedMode=dimnu
 
 

  write(*,*) 
  write(*,*) '          Name of the folder where calc_gDW.x is executed:  ',trim(basic_dir_path)
  write(*,*) '             Name of the file that stores the QE-scf input:  ',trim(InputQEscf_file_name)
  write(*,*) 'Name of the file that stores the QE input to print the pot:  ',trim(InputQEpot_file_name)
  write(*,*) ' File of phonon frequencies and eigvecs. of the dyn. matr.:  ',trim(FreqsU_file_name)
  write(*,'(A,I5)') '                                           Number of atoms = ', total_num_atoms
  if (PeriodicSystem .eq. 0) then
    write(*,'(A,I4,A,I4,A)') '                                 Number of phonon branches = ', dimnu, '  (',dimnu-6,' valid)'
  else
    write(*,'(A,I4,A,I4,A)') '                                 Number of phonon branches = ', dimnu, '  (',dimnu-3,' valid)'  
  end if
  write(*,'(A,f12.6)') '                                   Displacement Parameter = ',DisplParam
  if (PeriodicSystem .eq. 0) then
    write(*,*) '  The system is non-periodic (it has 6 unphysical frequencies).'
    Nunphysfreqs=6
  else
    write(*,*) '  The system is periodic (it has 3 unphysical frequencies).'
    Nunphysfreqs=3
  end if  
  write(*,*); write(*,*); 




  inquire(file=InputQEscf_file_name, exist=exist)      
  if (.not. exist) then
      write(*,*) 
      write(*,*) '**********************************************************************************' 
      write(*,*) '    ERROR: The file ',trim(InputQEscf_file_name),' does not exist. Please, provide it.'
      write(*,*) '**********************************************************************************'
      write(*,*)
      call exit(1) 
  end if

  inquire(file=FreqsU_file_name, exist=exist)      ! This checks whether the file gDW.dat is existing.
  if (.not. exist) then
      write(*,*) 
      write(*,*) '**********************************************************************************' 
      write(*,*) '    ERROR: The file ',trim(FreqsU_file_name),' does not exist. Please, provide it.'
      write(*,*) '**********************************************************************************'
      write(*,*)
      call exit(1) 
  end if
  
    
   inquire(file = InputQEpot_file_name, exist=exist)
   if(.not.exist) then
      write(*,*)
      write(*,*) '**********************************************************************************' 
      write(6,'(a)') " ERROR: A Quantum Espresso input file with the instructions to print the potentials (", &
                    & trim(InputQEpot_file_name),") is needed. Please, provide it."
      write(*,*) '**********************************************************************************' 
      write(*,*)
      call exit(1)
   endif
   
   inquire(file = InputQEscf_file_name, exist=exist)
   if(.not.exist) then
       write(*,*)
      write(*,*) '**********************************************************************************' 
      write(6,'(a)') " ERROR: A Quantum Espresso input file with the instructions to run the scf calculation (", &
                    & trim(InputQEscf_file_name),") is needed. Please, provide it."
      write(*,*) '**********************************************************************************' 
      write(*,*)
      call exit(1)
   endif

  
  
  
  
  ! PENGS BLOCK
     
!   read(90,*) Input_dot      !! The input file for dot eignemodes
!   read(90,*) Input_coord    !! The input file for the coorditotal_num_atomses
!   read(90,*) Nspecies          !! number of elements
!   allocate(atom_num(Nspecies))
!    allocate(atom_name(Nspecies))
!    allocate(atom_mass(Nspecies))
!    do ii = 1, Nspecies
!    read(90,*) atom_name(ii)  !! Read the name of element
!    read(90,*) atom_num(ii)  !! Read the number of each type
!    read(90,*) atom_mass(ii)  !! Read the mass of each atom
!    end do
!    read(90,*) n_mode         !! The phonon mode need to be calculated 
!    
!    close(90)
! 
! 
!    write(*,*) 'Input read from file:'
!    write(*,*) '---------------------'
!    write(*,*) 'Input coord file    :', Input_coord
!    write(*,*) 'Input dot eigenmode :', Input_dot
!    write(*,*) 'Type of atoms       :', Nspecies
!    do ii = 1, Nspecies
!    write(*,*) 'Atom type No.       :', ii
!    write(*,*) 'Name of the atom    :', atom_name(ii)
!    write(*,*) 'Number of atom      :', atom_num(ii)
!    write(*,*) 'Mass of the atom    :', atom_mass(ii)
!    end do
!    write(*,*) 'Phonon mode to be calculated:', n_mode
! 
!    total_num_atoms = 0
!    do ii = 1, Nspecies
!       total_num_atoms = total_num_atoms+atom_num(ii)
!       atom_mass(ii) = atom_mass(ii)*P_mass2au
!    end do
! 
!    !! allocate the arrays
!    allocate(atom_posi(total_num_atoms,3))
!    allocate(vib_mode(3*total_num_atoms,3*total_num_atoms))
!    allocate(omegaorig(3*total_num_atoms))
!    

  


end subroutine read_input_file



!-------------------------------------------------------------------------------------------




 subroutine read_nuclear_positions (total_num_atoms, InputQEscf_file_name, atom_posi)
 
 ! This subroutine reads the nuclear unperturbed positions from the file with name given by InputQEscf_file_name.
 
  Integer, Intent(In)        :: total_num_atoms              ! Total number of atoms
  character(100), Intent(In) :: InputQEscf_file_name            ! Name of the file where the nuclear positions and other infos are stored
  Real(8), Intent(Out)       :: atom_posi(total_num_atoms,3) ! Nuclear unperturbed positions of the system
  

  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character(2) :: auxstr4
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,iostat,unitt,positioninangstrom,notspccounter !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: latticeparameter, auxr1, auxr2, auxr3, maxposix, minposix, maxposiy, minposiy, maxposiz, minposiz
   
   positioninangstrom=0
   latticeparameter=-1.d0
   
   write(auxstr3,*) trim(InputQEscf_file_name) 
   call StripSpaces(auxstr3)
   write(*,*) '  *** Now reading the nuclear positions from the ',trim(InputQEscf_file_name),' file ***'
   
   open(838, file=InputQEscf_file_name, status='unknown')


   do k=1,500

      notspccounter=0
      read(838,'(A)',IOSTAT=stat) inputline
  
      if (stat .ne. 0) go to 1133

      auxstr1=''; auxstr2=''; auxstr3=''
   
      do i=1,100
        letter = inputline(i:i)
        if (letter .ne. '') notspccounter=notspccounter+1 
        auxstr1 = trim(auxstr1)//letter
        if ( (auxstr1(notspccounter-1:notspccounter-1) .eq. 'a' ) .and. (auxstr1(notspccounter:notspccounter).eq.'=' ) ) then
          go to 1235
        end if
      end do ! i
      
   end do ! k
   
   write(*,*) '  ERROR: No lattice parameter found. Please write e.g. "a=10.d0" in the'
   write(*,*) '   &system block of the ',trim(InputQEscf_file_name),' file'; write(*,*)
   call exit(1)

1235 continue
      
	do j=i+1,100
	  letter = inputline(j:j)
	  if (letter .eq. ',') cycle
	  auxstr2 = trim(auxstr2)//letter 
	end do ! j
	
    Read(auxstr2, '(G15.8)' ) latticeparameter
    
   ! write(*,*) '     Lattice Parameter = ',latticeparameter
    
1133 continue



   do k=1,500

      read(838,'(A)',IOSTAT=stat) inputline

!write(*,*) inputline
  
      if (stat .ne. 0) go to 9235

      auxstr1='';  auxstr2='';  auxstr3=''
   
!       do i=1,100
!         letter = inputline(i:i)
!         if ( (letter .ne. '}') .and. ((auxstr1 .ne. 'ATOMIC_POSITIONS { angstrom').or. &
!                       & (auxstr1.ne.'atomic_positions { angstrom ') ) ) then
!           auxstr1 = trim(auxstr1)//letter 
!           write(*,*) auxstr1
!         else 
!           positioninangstrom=1
!           write(*,*) 'positioninangstrom',positioninangstrom
!           stop
!           go to 9235
!         end if
!       end do ! i

      do i=1,100
        letter = inputline(i:i)
        if ( (letter .ne. '}') ) then
          auxstr1 = trim(auxstr1)//letter 
          ! write(*,*) auxstr1
          if ((auxstr1 .eq. 'ATOMIC_POSITIONS{angstrom').or. (auxstr1.eq.'atomic_positions{angstrom ') ) then
			positioninangstrom=1
			go to 9235            
          end if
        end if
      end do ! i

  end do ! k
   
9235 continue

  if (positioninangstrom .eq. 0) then
    write(auxstr3,*) trim(InputQEscf_file_name) 
    call StripSpaces(auxstr3)
	write(*,*) '  ERROR: Nuclear positions could not be read, or they are not in Angstrom.'
	write(*,*) '  Please, include them after the '
	write(*,*) '  word << ATOMIC_POSITIONS >> in the ',trim(auxstr3),'file.' 
	write(*,*)
	call exit(1)
  end if

   
   do i=1,total_num_atoms
      read(838,*) auxstr4, auxr1, auxr2, auxr3
      atom_posi(i,1) = auxr1
      atom_posi(i,2) = auxr2
      atom_posi(i,3) = auxr3
      ! write(*,*)  atom_posi(i,1),  atom_posi(i,2),  atom_posi(i,3)
   end do
   
   close(838)

! PGR: We avoid the symmetrization, because after a displacement in the positions, the 
! relaxation is (surprisingly) no longer valid. The scf.in contains the NOT symmetrized
! positions, because its coordinates result from the relaxation calculation, which does not
! symmetrize, and it is used to generate Veff_unperturb.pot. On the other hand, the Veff_modeX.pot
! result from symmetrized coordinates. This means that the symmetrization can be a source of errors.
! If one wants to use the code for calculating the density, the coordinates must be symmetrized.
! In these cases, symmetrize them by hand!!!  
!  
!    ! We now symmetrize the positions
!    maxposix=atom_posi(1,1); minposix=atom_posi(1,1)
!    maxposiy=atom_posi(1,2); minposiy=atom_posi(1,2)
!    maxposiz=atom_posi(1,3); minposiz=atom_posi(1,3)
!    
!    do i=2,total_num_atoms
!       if ( atom_posi(i,1) .gt.  maxposix)    maxposix=atom_posi(i,1)
!       if ( atom_posi(i,1) .lt.  minposix)    minposix=atom_posi(i,1)
!       if ( atom_posi(i,2) .gt.  maxposiy)    maxposiy=atom_posi(i,2)
!       if ( atom_posi(i,2) .lt.  minposiy)    minposiy=atom_posi(i,2)
!       if ( atom_posi(i,3) .gt.  maxposiz)    maxposiz=atom_posi(i,3)
!       if ( atom_posi(i,3) .lt.  minposiz)    minposiz=atom_posi(i,3)
!    end do   
!    
!    auxr1 = latticeparameter/2.d0 - (maxposix+minposix)/2.d0
!    auxr2 = latticeparameter/2.d0 - (maxposiy+minposiy)/2.d0
!    auxr3 = latticeparameter/2.d0 - (maxposiz+minposiz)/2.d0
! 
!    do i=1,total_num_atoms
!       atom_posi(i,1) =  atom_posi(i,1) + auxr1
!       atom_posi(i,2) =  atom_posi(i,2) + auxr2
!       atom_posi(i,3) =  atom_posi(i,3) + auxr3
!    end do
   
  ! Now the atomic positions are symmetric, and the center of the system lies in the center of the cell.

 end subroutine read_nuclear_positions

  
 
!-------------------------------------------------------------------------------------------



 subroutine write_new_QE_input_single_displ (positivedisplacement, total_num_atoms, dimnu, DisplParam, basic_dir_path, &
  & InputQEscf_file_name, InputQEpot_file_name, prefix, atom_posi, U)
 
 ! This subroutine writes the new input file for QE. Its nuclear coordinates are the relaxed ones
 ! (we assume that they are in angstrom) plus the displacement of the normal mode.
 
  logical, Intent(In)        :: positivedisplacement         ! Capaz's frozen-phonon method requires two displacements of opposite directions. This parameter specifies if the displacement is positive or negative.
  Integer, Intent(In)        :: total_num_atoms              ! Total number of atoms
  Integer, Intent(In)        :: dimnu                        ! Total number of phonon frequencies 
  Real(8), Intent(In)        :: DisplParam                   ! Parameter for the displacement; the nuclear positions of the newly generated QE inputs are the original ones plus the ( U^{\nu}_I vectors multiplied by DisplParam) 
  character(250), Intent(In) :: basic_dir_path               ! Path of the basic directory where the calc_gFan.x program is called
  character(100), Intent(In) :: InputQEscf_file_name         ! Name of the file where the nuclear positions and other infos are stored
  character(100), Intent(In) :: InputQEpot_file_name         ! Name of the file where the information for printing the potential
  character(100), Intent(In) :: prefix                       ! basic QE prefix (see prefix variable of pw.x)
  Real(8), Intent(In)        :: atom_posi(total_num_atoms,3) ! Nuclear unperturbed positions of the system
  Real(8), Intent(In)        :: U(dimnu,dimnu)               ! Eigenvectors of the dynamical matrix divided by the square root of the atomic masses (one eigenvector is in one column)
  

  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3, auxstr5, newQEinput,auxstr8
  character(20)  :: Tunits, stringT,auxstr7, auxstr88
  character(50) :: auxstr9
  character(250) :: auxstr6
  character(4) :: auxstr4
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l,nu,iostat,unitt,withinblock,blankcounter !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, atom_new(total_num_atoms,3), displsign
   
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !!!!!!!!!!!! BLOCK 1: WRITING THE QE-SCF FILES  !!!!!!!!!!!!!!  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if ( positivedisplacement ) then 
     write(*,*) '  *** Now writing the new QE-scf input files (see << frozen+R >> folder) ***'
     inquire(file='frozen+R', exist=exist)     
     if ( .not. exist  ) call system('mkdir frozen+R')  
     displsign = 1.d0
   else  
     write(*,*) '  *** Now writing the new QE-scf input files (see << frozen-R >> folder) ***'
     inquire(file='frozen-R', exist=exist)     
     if ( .not. exist  ) call system('mkdir frozen-R') 
     displsign = -1.d0    
   end if  
   
!    k=1
!     do i=1,1!total_num_atoms
!      write(*,*) U((i-1)*3+1,k)  
!      write(*,*) U((i-1)*3+2,k)  
!      write(*,*) U((i-1)*3+3,k) 
!     end do       
!    
!    write(*,*);write(*,*)
   
   do nu=1,dimnu
   
    ! We write the positions of the nuclei for the new calculations of the potentials (for one given phonon mode)  
     do j=1,total_num_atoms
      atom_new(j,1) = atom_posi(j,1) + displsign*DisplParam * U(3*(j-1)+1,nu)  * 0.529177249d0  ! U is in atomic units; we write it in Angstrom
      atom_new(j,2) = atom_posi(j,2) + displsign*DisplParam * U(3*(j-1)+2,nu)  * 0.529177249d0 
      atom_new(j,3) = atom_posi(j,3) + displsign*DisplParam * U(3*(j-1)+3,nu)  * 0.529177249d0 
     end do 

     
!      write(*,*) 'atom_posi(1,1)',atom_posi(1,1)
!      write(*,*) 'U(3*(j-1)+1,nu)',U(1,1)
!      write(*,*) 'scl',DisplParam * U(1,1)  * 0.529177249d0 
!      write(*,*) 'atom_new(j,1) ',atom_new(1,1) 
!      stop
     
    ! The new atomic positions are the original ones plus U (not plus xi). As you one see in Peng's original code (vib_pot.f90):
    ! vib_mode(ii,nn) = vib_mode(ii,nn)/sqrt(atom_mass(jj))  
    !  (...)  
    ! atom_new(ii,jj) = atom_posi(ii,jj) + vib_mode(n_mode,(ii-1)*3+jj)
    
 
    if ( positivedisplacement ) then 
	   ! We define the name of the file that we will write              
	   if ((dimnu .lt. 100) ) then
		 write (newQEinput,"(A,I2,A7)") 'frozen+R/displaced-mode',nu,'_scf.in'
	   else if ((dimnu .lt. 1000)  ) then 
		 write (newQEinput,"(A,I3,A7)") 'frozen+R/displaced-mode',nu,'_scf.in'
	   else if ((dimnu .lt. 10000)  ) then 
		 write (newQEinput,"(A,I4,A7)") 'frozen+R/displaced-mode',nu,'_scf.in'
	   else
		 write (newQEinput,"(A,I6,A7)") 'frozen+R/displaced-mode',nu,'_scf.in'                  
	   end if  
	   call StripSpaces (newQEinput)    
	 else 
	   if ((dimnu .lt. 100) ) then
		 write (newQEinput,"(A,I2,A7)") 'frozen-R/displaced-mode',nu,'_scf.in'
	   else if ((dimnu .lt. 1000)  ) then 
		 write (newQEinput,"(A,I3,A7)") 'frozen-R/displaced-mode',nu,'_scf.in'
	   else if ((dimnu .lt. 10000)  ) then 
		 write (newQEinput,"(A,I4,A7)") 'frozen-R/displaced-mode',nu,'_scf.in'
	   else
		 write (newQEinput,"(A,I6,A7)") 'frozen-R/displaced-mode',nu,'_scf.in'                  
	   end if  
	   call StripSpaces (newQEinput)   
	 end if  
      
   
    open(839, file=InputQEscf_file_name, status='unknown')
    open(840, file=newQEinput, status='unknown')
   
   
    withinblock = 0


    ! Copying: 1st block: we copy the text before the atomic positions
    do k=1,200

      auxstr1='';  auxstr2='';  auxstr3=''; auxstr4=''; auxstr5=''
            
      read(839,'(A)',IOSTAT=iostat) inputline 
      do i=1,100
       letter = inputline(i:i)
       auxstr5 = trim(auxstr5)//letter        
      end do
      
      
      if (auxstr5(1:1) .eq. '&') then
         withinblock=withinblock+1
         write (840,*) trim(inputline)
      else if (auxstr5(1:1) .eq. '/') then 
         withinblock=withinblock-1
         write (840,*) trim(inputline)
         write (840,*) 
      else  
         if (withinblock .eq. 1) then
             call StripSpaces(inputline)
             write (840,*) '   '//trim(inputline)
         else 
            auxstr9=''
            do i=1,50
              letter = inputline(i:i)
              auxstr9(i:i) = letter 
            end do
            blankcounter=0
            do i=51,100
              letter = inputline(i:i)
              if (letter .eq. ' ') blankcounter=blankcounter+1
            end do
            if (blankcounter .eq.50) then
              write (840,*) trim(auxstr9)
            else
              write (840,*) trim(inputline)
            end if
          
         end if   
      end if  
      
      if ( (iostat .ne. 0) ) go to 4561
      
      call StripSpaces(inputline)
      if ((trim(inputline) .eq. '&control')  .or. (trim(inputline)  .eq. '&CONTROL')  ) then
  !         write(840,*) "   prefix='",trim(prefix),"',"
!           write(auxstr7,*) nu,"/',"
!           call StripSpaces(auxstr7)
!           write(auxstr6,*) '    outdir='//"'"//trim(basic_dir_path)//"',"!//'mode'//trim(auxstr7)
!           call StripSpaces(auxstr6)
!           write(840,*) '  '//trim(auxstr6)
!           write(auxstr7,*) nu,"/',"
!           call StripSpaces(auxstr7)
          write(auxstr6,*) 'wfcdir='//"'"//trim(basic_dir_path)//"',"!//'mode'//trim(auxstr7)
          call StripSpaces(auxstr6)
          write(840,*) '   '//trim(auxstr6)
      end if
            
            
      do i=1,100
        letter = inputline(i:i)
        if ( ((auxstr1 .ne. 'ATOMIC_POSITIONS{').and.(auxstr1.ne.'atomic_positions{') ) ) then
          auxstr1 = trim(auxstr1)//letter 
        else 
               ! Copying: 2nd block: we write the new atomic positions
                do j=1,total_num_atoms
                   read(839,*) auxstr4, auxr1, auxr2, auxr3
                   write(840,*) auxstr4, atom_new(j,1), atom_new(j,2), atom_new(j,3)
                end do ! j

                ! Copying: 3rd block: we copy the text after the atomic positions
                do l=1,100
                  read(839,'(A)',IOSTAT=iostat) inputline
                  write (840,*) trim(inputline)
                  if ( (iostat .ne. 0) ) go to 4561  
                end do ! l
                
          go to 4561
        end if
      end do ! i   
      
    end do ! k
    
4561 continue
    
    close(839); close(840)
    
    
    
     
   end do ! nu
     

   
   
   
   
   
   
   
  
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !!!!!!!!!!!! BLOCK 2: WRITING THE QE-POTENTIAL-WRITING FILES !!!!!!!!!!!!!!!!!!!!!!!!!!  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if ( positivedisplacement) then
     write(*,*) '  *** Now writing the new QE-potential-writing input files' 
     write(*,*) '    (see << frozen+R >> folder) ***'
   else  
     write(*,*) '  *** Now writing the new QE-potential-writing input files' 
     write(*,*) '    (see << frozen-R >> folder) ***'   
   end if
   
   
   do nu=1,dimnu
   
    if ( positivedisplacement) then  
	   ! We define the name of the file that we will write              
	   if ((dimnu .lt. 100) ) then
		 write (newQEinput,"(A,I2,A8)") 'frozen+R/displaced-mode',nu,'_Veff.in'
	   else if ((dimnu .lt. 1000)  ) then 
		 write (newQEinput,"(A,I3,A8)") 'frozen+R/displaced-mode',nu,'_Veff.in'
	   else if ((dimnu .lt. 10000)  ) then 
		 write (newQEinput,"(A,I4,A8)") 'frozen+R/displaced-mode',nu,'_Veff.in'
	   else
		 write (newQEinput,"(A,I6,A8)") 'frozen+R/displaced-mode',nu,'_Veff.in'                  
	   end if  
	   call StripSpaces (newQEinput)    
    else
	   ! We define the name of the file that we will write              
	   if ((dimnu .lt. 100) ) then
		 write (newQEinput,"(A,I2,A8)") 'frozen-R/displaced-mode',nu,'_Veff.in'
	   else if ((dimnu .lt. 1000)  ) then 
		 write (newQEinput,"(A,I3,A8)") 'frozen-R/displaced-mode',nu,'_Veff.in'
	   else if ((dimnu .lt. 10000)  ) then 
		 write (newQEinput,"(A,I4,A8)") 'frozen-R/displaced-mode',nu,'_Veff.in'
	   else
		 write (newQEinput,"(A,I6,A8)") 'frozen-R/displaced-mode',nu,'_Veff.in'                  
	   end if  
	   call StripSpaces (newQEinput)        
    end if  
   
    open(839, file=InputQEpot_file_name, status='unknown')
    open(840, file=newQEinput, status='unknown')
   
   
    withinblock = 0


    ! Copying: 1st block: we copy the text before the atomic positions
    do k=1,200

      auxstr1='';  auxstr2='';  auxstr3=''; auxstr4=''; auxstr5=''
            
      read(839,'(A)',IOSTAT=iostat) inputline 
      
      do i=1,100
       letter = inputline(i:i)
       auxstr5 = trim(auxstr5)//letter 
 
!  PGR: I comment this line to include a prefix:       
!        if (((auxstr5 .eq. 'prefix').or.(auxstr5 .eq. 'Prefix'))) then
!            auxstr8='' 
!            do j=1,5
!              letter = inputline(i+j:i+j)
!              auxstr88 = trim(auxstr88)//letter
!            end do
!           if ((auxstr88.eq."=''").or.(auxstr88.eq."='',")) then 
!           write(*,*) 'pasamos'
!             go to 9535
!           end if      
!        end if     
       
      end do
      
      
      if (auxstr5(1:1) .eq. '&') then
         withinblock=withinblock+1
         write (840,*) trim(inputline)       
      else if (auxstr5(1:1) .eq. '/') then 
         withinblock=withinblock-1
         write (840,*) trim(inputline)
         write (840,*) 
      else  
      
        do i=1,90
         letter = inputline(i:i)
         auxstr1 = trim(auxstr1)//letter 
         !if ((auxstr1 .eq. 'prefix').or.(auxstr1 .eq. 'Prefix')) then      
         !  go to 9535
         !end if
         if ((auxstr1 .eq. 'filplot').or. (auxstr1 .eq. 'Filplot')) then      
           go to 9535
         end if
        end do ! i
      
      
         if (withinblock .eq. 1) then
             call StripSpaces(inputline)
             write (840,*) '   '//trim(inputline)
         else 
             write (840,*) inputline
         end if   
      end if  
      
      if ( (iostat .ne. 0) ) go to 4233 
      
      call StripSpaces(inputline)
      if ((trim(inputline) .eq. '&inputPP')  .or. (trim(inputline)  .eq. '&INPUTPP') &
         & .or. (trim(inputline)  .eq. '&inputpp')  ) then
          write(auxstr7,*) nu,".pot',"
          call StripSpaces(auxstr7)
          write(auxstr6,*) '    filplot='//"'"//'Veff_mode'//trim(auxstr7)
          call StripSpaces(auxstr6)
          write(840,*) '  '//trim(auxstr6)
        !  write(840,*) "  prefix='',"
      end if 
 
 9535 continue
      
    end do ! k
    
4233 continue
    
    close(839); close(840)   
     
   end do ! nu
   
   

   

 end subroutine write_new_QE_input_single_displ 
  
 
!-------------------------------------------------------------------------------------------




 subroutine write_new_QE_input_double_displ (total_num_atoms, dimnu, DisplParam, basic_dir_path, &
  & InputQEscf_file_name, InputQEpot_file_name, prefix, atom_posi, U)
 
 ! This subroutine writes the new input file for QE. Its nuclear coordinates are the relaxed ones
 ! (we assume that they are in angstrom) plus the displacement of the normal mode.
 
  Integer, Intent(In)        :: total_num_atoms              ! Total number of atoms
  Integer, Intent(In)        :: dimnu                        ! Total number of phonon frequencies 
  Real(8), Intent(In)        :: DisplParam                   ! Parameter for the displacement; the nuclear positions of the newly generated QE inputs are the original ones plus the ( U^{\nu}_I vectors multiplied by DisplParam) 
  character(250), Intent(In) :: basic_dir_path               ! Path of the basic directory where the calc_gFan.x program is called
  character(100), Intent(In) :: InputQEscf_file_name         ! Name of the file where the nuclear positions and other infos are stored
  character(100), Intent(In) :: InputQEpot_file_name         ! Name of the file where the information for printing the potential
  character(100), Intent(In) :: prefix                       ! basic QE prefix (see prefix variable of pw.x)
  Real(8), Intent(In)        :: atom_posi(total_num_atoms,3) ! Nuclear unperturbed positions of the system
  Real(8), Intent(In)        :: U(dimnu,dimnu)               ! Eigenvectors of the dynamical matrix divided by the square root of the atomic masses (one eigenvector is in one column)
  

  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3, auxstr5, newQEinput,auxstr8
  character(20)  :: Tunits, stringT,auxstr7, auxstr88
  character(50) :: auxstr9
  character(250) :: auxstr6
  character(4) :: auxstr4
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l,isgn1,isgn2,nu,nu1,nu2,iostat,unitt,withinblock,blankcounter !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, atom_new(total_num_atoms,3), displsign1, displsign2
   
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !!!!!!!!!!!! BLOCK 1: WRITING THE QE-SCF FILES  !!!!!!!!!!!!!!  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  
  write(*,*) '  *** Now writing the new QE-scf input files with two displacements of'
  write(*,*) '      (U) eigenvectors of the dynamical matrix (see << frozenRR >> folder) ***'
  inquire(file='frozenRR', exist=exist)     
  if ( .not. exist  ) call system('mkdir frozenRR')    
   
   
do isgn1=-1,1,2

 displsign1 = dble(isgn1)
 
 do isgn2=-1,1,2   
 
   displsign2 = dble(isgn2)
   
   do nu1=1,dimnu
     do nu2=1,dimnu
   
	  ! We write the positions of the nuclei for the new calculations of the potentials (for one given phonon mode)  
	   do j=1,total_num_atoms
		atom_new(j,1) = atom_posi(j,1) + displsign1*DisplParam * U(3*(j-1)+1,nu1)  * 0.529177249d0  ! U is in atomic units; we write it in Angstrom
		atom_new(j,2) = atom_posi(j,2) + displsign1*DisplParam * U(3*(j-1)+2,nu1)  * 0.529177249d0 
		atom_new(j,3) = atom_posi(j,3) + displsign1*DisplParam * U(3*(j-1)+3,nu1)  * 0.529177249d0 
	   end do 

	   do j=1,total_num_atoms
		atom_new(j,1) = atom_new(j,1) + displsign2*DisplParam * U(3*(j-1)+1,nu2)  * 0.529177249d0  ! U is in atomic units; we write it in Angstrom
		atom_new(j,2) = atom_new(j,2) + displsign2*DisplParam * U(3*(j-1)+2,nu2)  * 0.529177249d0 
		atom_new(j,3) = atom_new(j,3) + displsign2*DisplParam * U(3*(j-1)+3,nu2)  * 0.529177249d0 
	   end do 
     
!      write(*,*) 'atom_posi(1,1)',atom_posi(1,1)
!      write(*,*) 'U(3*(j-1)+1,nu)',U(1,1)
!      write(*,*) 'scl',DisplParam * U(1,1)  * 0.529177249d0 
!      write(*,*) 'atom_new(j,1) ',atom_new(1,1) 
!      stop
     
    ! The new atomic positions are the original ones plus U (not plus xi). As you one see in Peng's original code (vib_pot.f90):
    ! vib_mode(ii,nn) = vib_mode(ii,nn)/sqrt(atom_mass(jj))  
    !  (...)  
    ! atom_new(ii,jj) = atom_posi(ii,jj) + vib_mode(n_mode,(ii-1)*3+jj)
    
 
     if ( (isgn1.eq.-1) .and. (isgn2.eq.-1) ) then       
	   write (newQEinput,"(A,I6,A,I6,A7)") 'frozenRR/displacedmode-',nu1,'-',nu2,'_scf.in'  
	 end if   
     if ( (isgn1.eq.1) .and. (isgn2.eq.-1) ) then       
	   write (newQEinput,"(A,I6,A,I6,A7)") 'frozenRR/displacedmode+',nu1,'-',nu2,'_scf.in'  
	 end if 
     if ( (isgn1.eq.-1) .and. (isgn2.eq.1) ) then       
	   write (newQEinput,"(A,I6,A,I6,A7)") 'frozenRR/displacedmode-',nu1,'+',nu2,'_scf.in'  
	 end if 	 	                
     if ( (isgn1.eq.1) .and. (isgn2.eq.1) ) then       
	   write (newQEinput,"(A,I6,A,I6,A7)") 'frozenRR/displacedmode+',nu1,'+',nu2,'_scf.in'  
	 end if 
     call StripSpaces (newQEinput)    
 
 
      
 
 
!     if ( positivedisplacement ) then 
! 	   ! We define the name of the file that we will write              
! 	   if ((dimnu .lt. 100) ) then
! 		 write (newQEinput,"(A,I2,A7)") 'frozen+R/displaced-mode',nu,'_scf.in'
! 	   else if ((dimnu .lt. 1000)  ) then 
! 		 write (newQEinput,"(A,I3,A7)") 'frozen+R/displaced-mode',nu,'_scf.in'
! 	   else if ((dimnu .lt. 10000)  ) then 
! 		 write (newQEinput,"(A,I4,A7)") 'frozen+R/displaced-mode',nu,'_scf.in'
! 	   else
! 		 write (newQEinput,"(A,I6,A,I6,A7)") 'frozen+R/displaced-mode',nu,'_scf.in'                  
! 	   end if  
! 	   call StripSpaces (newQEinput)    
! 	 else 
! 	   if ((dimnu .lt. 100) ) then
! 		 write (newQEinput,"(A,I2,A7)") 'frozen-R/displaced-mode',nu,'_scf.in'
! 	   else if ((dimnu .lt. 1000)  ) then 
! 		 write (newQEinput,"(A,I3,A7)") 'frozen-R/displaced-mode',nu,'_scf.in'
! 	   else if ((dimnu .lt. 10000)  ) then 
! 		 write (newQEinput,"(A,I4,A7)") 'frozen-R/displaced-mode',nu,'_scf.in'
! 	   else
! 		 write (newQEinput,"(A,I6,A,I6,A7)") 'frozen-R/displaced-mode',nu,'_scf.in'                  
! 	   end if  
! 	   call StripSpaces (newQEinput)   
! 	 end if  
!       
   
    open(839, file=InputQEscf_file_name, status='unknown')
    open(840, file=newQEinput, status='unknown')
   
   
    withinblock = 0


    ! Copying: 1st block: we copy the text before the atomic positions
    do k=1,200

      auxstr1='';  auxstr2='';  auxstr3=''; auxstr4=''; auxstr5=''
            
      read(839,'(A)',IOSTAT=iostat) inputline 
      do i=1,100
       letter = inputline(i:i)
       auxstr5 = trim(auxstr5)//letter        
      end do
      
      
      if (auxstr5(1:1) .eq. '&') then
         withinblock=withinblock+1
         write (840,*) trim(inputline)
      else if (auxstr5(1:1) .eq. '/') then 
         withinblock=withinblock-1
         write (840,*) trim(inputline)
         write (840,*) 
      else  
         if (withinblock .eq. 1) then
             call StripSpaces(inputline)
             write (840,*) '   '//trim(inputline)
         else 
            auxstr9=''
            do i=1,50
              letter = inputline(i:i)
              auxstr9(i:i) = letter 
            end do
            blankcounter=0
            do i=51,100
              letter = inputline(i:i)
              if (letter .eq. ' ') blankcounter=blankcounter+1
            end do
            if (blankcounter .eq.50) then
              write (840,*) trim(auxstr9)
            else
              write (840,*) trim(inputline)
            end if
          
         end if   
      end if  
      
      if ( (iostat .ne. 0) ) go to 4561
      
      call StripSpaces(inputline)
      if ((trim(inputline) .eq. '&control')  .or. (trim(inputline)  .eq. '&CONTROL')  ) then
  !         write(840,*) "   prefix='",trim(prefix),"',"
!           write(auxstr7,*) nu,"/',"
!           call StripSpaces(auxstr7)
!           write(auxstr6,*) '    outdir='//"'"//trim(basic_dir_path)//"',"!//'mode'//trim(auxstr7)
!           call StripSpaces(auxstr6)
!           write(840,*) '  '//trim(auxstr6)
!           write(auxstr7,*) nu,"/',"
!           call StripSpaces(auxstr7)
          write(auxstr6,*) 'wfcdir='//"'"//trim(basic_dir_path)//"',"!//'mode'//trim(auxstr7)
          call StripSpaces(auxstr6)
          write(840,*) '   '//trim(auxstr6)
      end if
            
            
      do i=1,100
        letter = inputline(i:i)
        if ( ((auxstr1 .ne. 'ATOMIC_POSITIONS{').and.(auxstr1.ne.'atomic_positions{') ) ) then
          auxstr1 = trim(auxstr1)//letter 
        else 
               ! Copying: 2nd block: we write the new atomic positions
                do j=1,total_num_atoms
                   read(839,*) auxstr4, auxr1, auxr2, auxr3
                   write(840,*) auxstr4, atom_new(j,1), atom_new(j,2), atom_new(j,3)
                end do ! j

                ! Copying: 3rd block: we copy the text after the atomic positions
                do l=1,100
                  read(839,'(A)',IOSTAT=iostat) inputline
                  write (840,*) trim(inputline)
                  if ( (iostat .ne. 0) ) go to 4561  
                end do ! l
                
          go to 4561
        end if
      end do ! i   
      
    end do ! k
    
4561 continue
    
    close(839); close(840)
    
    
    
     end do ! nu2
   end do ! nu1
  end do ! isgn2   
end do ! isgn1
   
   
   
   
   
   
   
  
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !!!!!!!!!!!! BLOCK 2: WRITING THE QE-POTENTIAL-WRITING FILES !!!!!!!!!!!!!!!!!!!!!!!!!!  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   write(*,*) '  *** Now writing the new QE-potential-writing input files' 
   write(*,*) '    (see << frozenRR >> folder) ***'


do isgn1=-1,1,2

 displsign1 = dble(isgn1)
 
 do isgn2=-1,1,2   
 
   displsign2 = dble(isgn2)
   
   do nu1=1,dimnu
     do nu2=1,dimnu
     
     
      if ( (isgn1.eq.-1) .and. (isgn2.eq.-1) ) then       
	   write (newQEinput,"(A,I6,A,I6,A8)") 'frozenRR/displacedmode-',nu1,'-',nu2,'_Veff.in'  
	 end if   
     if ( (isgn1.eq.1) .and. (isgn2.eq.-1) ) then       
	   write (newQEinput,"(A,I6,A,I6,A8)") 'frozenRR/displacedmode+',nu1,'-',nu2,'_Veff.in'  
	 end if 
     if ( (isgn1.eq.-1) .and. (isgn2.eq.1) ) then       
	   write (newQEinput,"(A,I6,A,I6,A8)") 'frozenRR/displacedmode-',nu1,'+',nu2,'_Veff.in'  
	 end if 	 	                
     if ( (isgn1.eq.1) .and. (isgn2.eq.1) ) then       
	   write (newQEinput,"(A,I6,A,I6,A8)") 'frozenRR/displacedmode+',nu1,'+',nu2,'_Veff.in'  
	 end if 
     call StripSpaces (newQEinput)
     
   
    ! if ( positivedisplacement) then  
! 	   ! We define the name of the file that we will write              
! 	   if ((dimnu .lt. 100) ) then
! 		 write (newQEinput,"(A,I2,A8)") 'frozen+R/displaced-mode',nu,'_Veff.in'
! 	   else if ((dimnu .lt. 1000)  ) then 
! 		 write (newQEinput,"(A,I3,A8)") 'frozen+R/displaced-mode',nu,'_Veff.in'
! 	   else if ((dimnu .lt. 10000)  ) then 
! 		 write (newQEinput,"(A,I4,A8)") 'frozen+R/displaced-mode',nu,'_Veff.in'
! 	   else
! 		 write (newQEinput,"(A,I6,A8)") 'frozen+R/displaced-mode',nu,'_Veff.in'                  
! 	   end if  
! 	   call StripSpaces (newQEinput)    
!     else
! 	   ! We define the name of the file that we will write              
! 	   if ((dimnu .lt. 100) ) then
! 		 write (newQEinput,"(A,I2,A8)") 'frozen-R/displaced-mode',nu,'_Veff.in'
! 	   else if ((dimnu .lt. 1000)  ) then 
! 		 write (newQEinput,"(A,I3,A8)") 'frozen-R/displaced-mode',nu,'_Veff.in'
! 	   else if ((dimnu .lt. 10000)  ) then 
! 		 write (newQEinput,"(A,I4,A8)") 'frozen-R/displaced-mode',nu,'_Veff.in'
! 	   else
! 		 write (newQEinput,"(A,I6,A8)") 'frozen-R/displaced-mode',nu,'_Veff.in'                  
! 	   end if  
! 	   call StripSpaces (newQEinput)        
!     end if  
   
    open(839, file=InputQEpot_file_name, status='unknown')
    open(840, file=newQEinput, status='unknown')
   
   
    withinblock = 0


    ! Copying: 1st block: we copy the text before the atomic positions
    do k=1,200

      auxstr1='';  auxstr2='';  auxstr3=''; auxstr4=''; auxstr5=''
            
      read(839,'(A)',IOSTAT=iostat) inputline 
      
      do i=1,100
       letter = inputline(i:i)
       auxstr5 = trim(auxstr5)//letter 
       
!        if (((auxstr5 .eq. 'prefix').or.(auxstr5 .eq. 'Prefix'))) then
!            auxstr8='' 
!            do j=1,5
!              letter = inputline(i+j:i+j)
!              auxstr88 = trim(auxstr88)//letter
!            end do
!           if ((auxstr88.eq."=''").or.(auxstr88.eq."='',")) then 
!           write(*,*) 'pasamos'
!             go to 9535
!           end if      
!        end if     
       
      end do
      
      
      if (auxstr5(1:1) .eq. '&') then
         withinblock=withinblock+1
         write (840,*) trim(inputline)       
      else if (auxstr5(1:1) .eq. '/') then 
         withinblock=withinblock-1
         write (840,*) trim(inputline)
         write (840,*) 
      else  
      
        do i=1,90
         letter = inputline(i:i)
         auxstr1 = trim(auxstr1)//letter 
         if ((auxstr1 .eq. 'filplot').or. (auxstr1 .eq. 'Filplot')) then      
           go to 9535
         end if
    !     if ((auxstr1 .eq. 'prefix').or.(auxstr1 .eq. 'Prefix')) then      
    !       go to 9535
    !     end if
        end do ! i
      
      
         if (withinblock .eq. 1) then
             call StripSpaces(inputline)
             write (840,*) '   '//trim(inputline)
         else 
             write (840,*) inputline
         end if   
      end if  
      
      if ( (iostat .ne. 0) ) go to 4233 
      
      call StripSpaces(inputline)
      if ((trim(inputline) .eq. '&inputPP')  .or. (trim(inputline)  .eq. '&INPUTPP') &
         & .or. (trim(inputline)  .eq. '&inputpp')  ) then
          write(auxstr7,*) ".pot',"
          call StripSpaces(auxstr7)
          write(auxstr6,*) '    filplot='//"'"//'Veff'//trim(auxstr7)
          call StripSpaces(auxstr6)
          write(840,*) '  '//trim(auxstr6)
        !  write(840,*) "  prefix='',"
      end if 
 
 9535 continue
      
    end do ! k
    
4233 continue
    
    close(839); close(840)   
     
   end do ! nu2
  end do ! nu1 
 end do ! isgn2
end do ! isgn1   
   

   

 end subroutine write_new_QE_input_double_displ 
  
 
!-------------------------------------------------------------------------------------------



 subroutine write_U_files (filename,dimnu,total_num_atoms,DisplParam,omegaorig,U)
 
 ! This subroutine writes the Input_dot file. In it, we write the U vectors and the phonon frequencies
 ! in the CPMD format; this is, we write them twice (as required by the elph_pwscf program).
 ! Input_dot contains the U eigenvectors of the dynamical matrix multiplied by the DisplacementParameter,
 ! this is, not normalized. Input_dot is demanded by the elph_pwscf program.
 ! In addition, this subroutine generates the Uvec.dat file, that is like Input_dot, but with DisplacedParameter=1
 ! which is the value that the genmaatinput.x program needs.
 
 
  character(100), Intent(In) :: filename                     ! Name of the file to write (either Input_dot or Uvec.dat)
  Integer, Intent(In)        :: dimnu                        ! Total number of phonon frequencies 
  Integer, Intent(In)        :: total_num_atoms              ! Total number of atoms
  Real(8), Intent(In)        :: DisplParam                   ! Displacement parameter for the finite differences
  Real(8), Intent(In)        :: omegaorig(dimnu)             ! Phonon frequencies
  Real(8), Intent(In)        :: U(dimnu,dimnu)               ! Eigenvectors of the dynamical matrix divided by the square root of the atomic masses (one eigenvector is in one column)
  
  
  ! Local variables
  integer :: stat,auxi,i,j,k,l,nu
  real(8) :: auxr1, auxr2, auxr3
  real(8), allocatable :: neoU(:,:) 
   
   
!     open(439, file='DisplParam.out', status='unknown')  
!       write(439,*) 'DislplParam=',DisplParam
!     close(439)
 
    call StripSpaces(filename)
    write(*,*) '  *** Now writing the ',trim(filename),' file ***'
  
    allocate(neoU(dimnu,dimnu))
    
    do i=1,dimnu
      do nu=1,dimnu
        neoU(i,nu) = DisplParam * U(i,nu)   
      end do
    end do 
    
    open(439, file='frequencies.out', status='unknown')  
    do j=1,dimnu
      write(439,*) 'index ',j,'omega', omegaorig(j)
    end do
    close(439)

    open(439, file='frequenciesbare.out', status='unknown')  
    do j=1,dimnu
      write(439,*) omegaorig(j)
    end do
    close(439)  

!  In Input_dot (a file which is necessary for elph_pwscf) 
!  the U eigenvectors of the dynamical matrix must be in atomic units.
!  As it can be viewed in Peng's original code (vib_pot.f90): 
!
!    !! Read the vibrational modes
!    call dot_eigen(Input_coord, Input_dot, nat, atom_posi, vib_mode, vib_freq)  ! This reads Input_dot and stores the (xi) eigenvectors in the vib_mode variable
!    (...)
!     vib_mode(ii,nn) = vib_mode(ii,nn)/sqrt(atom_mass(jj)) ! This mass is in atomic units: atom_mass(ii) = atom_mass(ii)*P_mass2au ! P_mass2au=1822.84
!    (...) 
!    Later Peng writes the vib_modes in angstrom to rewrite the atomic positions for the input files:
!    !! Transfer to au2A
!    do ii = 1, 3*nat
!       vib_mode(ii,:) = vib_mode(ii,:) * P_au2A
!    end do

    

    
  !   open(439, file='U.out', status='unknown')  
!     do j=1,dimnu
!       write(439,*) 'index ',j,'omega', omegaorig(j)
!       do i=1,dimnu
!         write(439,*) U(i,j)
!       end do
!       write(439,*)
!     end do
!     close(439)
! 
!  
!     open(339, file='neoU.out', status='unknown')  
!     do j=1,dimnu
!       write(339,*) 'index ',j,'omega', omegaorig(j)
!       do i=1,dimnu
!         write(339,*) neoU(i,j)
!       end do
!       write(339,*)
!     end do
!     close(339)  



    j=dimnu/8
    k=mod(dimnu,8)

    open(739, file=filename, status='unknown')
     
    do i=1,j
        write(739,'(I15,A,I15,A,I15,A,I15,A,I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4, '    ', &
          & (i-1)*8+5, '    ',(i-1)*8+6,'    ', (i-1)*8+7, '    ',(i-1)*8+8
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4), '    ',omegaorig((i-1)*8+5), '    ', omegaorig((i-1)*8+6), '    ', &
               &   omegaorig((i-1)*8+7), '    ', omegaorig((i-1)*8+8)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4), '    ', &
               &  neoU(l,(i-1)*8+5), '    ',  neoU(l,(i-1)*8+6), '    ',  neoU(l,(i-1)*8+7), '    ', neoU(l,(i-1)*8+8)   
       end do
                                
    end do ! i
    
    
    if (k.eq.7) then
    
        i=j+1
        write(739,'(I15,A,I15,A,I15,A,I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4, '    ', &
          & (i-1)*8+5, '    ',(i-1)*8+6,'    ', (i-1)*8+7
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4), '    ',omegaorig((i-1)*8+5), '    ', omegaorig((i-1)*8+6), '    ', &
               &   omegaorig((i-1)*8+7)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4), '    ', &
               &  neoU(l,(i-1)*8+5), '    ',  neoU(l,(i-1)*8+6), '    ',  neoU(l,(i-1)*8+7) 
       end do
                                

    
    else if (k.eq.6) then


        i=j+1
        write(739,'(I15,A,I15,A,I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4, '    ', &
          & (i-1)*8+5, '    ',(i-1)*8+6
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4), '    ',omegaorig((i-1)*8+5), '    ', omegaorig((i-1)*8+6)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4), '    ', &
               &  neoU(l,(i-1)*8+5), '    ',  neoU(l,(i-1)*8+6)
       end do
                                

    
    else if (k.eq.5) then

        i=j+1
        write(739,'(I15,A,I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4, '    ', &
          & (i-1)*8+5
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4), '    ',omegaorig((i-1)*8+5)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4), '    ', &
               &  neoU(l,(i-1)*8+5)
       end do ! l
    
    
    else if (k.eq.4) then

       i=j+1
        write(739,'(I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4)
       end do ! l
    
    
    else if (k.eq.3) then

       i=j+1
        write(739,'(I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3
       write(739,'(G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3)
       end do ! l
    
    
    else if (k.eq.2) then
 
        i=j+1
        write(739,'(I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2
       write(739,'(G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2)
       end do ! l
    

    
    else if (k.eq.1) then
 
         i=j+1
        write(739,'(I15)') &
        &  (i-1)*8+1
       write(739,'(G15.7,A)') & 
               & omegaorig((i-1)*8+1)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7)') &
               &  neoU(l,(i-1)*8+1)
       end do ! l
    
    
    end if



    write(739,*) '<<<<<<  NEW DATA  >>>>>>'

   do i=1,j
        write(739,'(I15,A,I15,A,I15,A,I15,A,I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4, '    ', &
          & (i-1)*8+5, '    ',(i-1)*8+6,'    ', (i-1)*8+7, '    ',(i-1)*8+8
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4), '    ',omegaorig((i-1)*8+5), '    ', omegaorig((i-1)*8+6), '    ', &
               &   omegaorig((i-1)*8+7), '    ', omegaorig((i-1)*8+8)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4), '    ', &
               &  neoU(l,(i-1)*8+5), '    ',  neoU(l,(i-1)*8+6), '    ',  neoU(l,(i-1)*8+7), '    ', neoU(l,(i-1)*8+8)   
       end do
                                
    end do ! i
    
 
    
    if (k.eq.7) then
    
        i=j+1
        write(739,'(I15,A,I15,A,I15,A,I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4, '    ', &
          & (i-1)*8+5, '    ',(i-1)*8+6,'    ', (i-1)*8+7
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4), '    ',omegaorig((i-1)*8+5), '    ', omegaorig((i-1)*8+6), '    ', &
               &   omegaorig((i-1)*8+7)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4), '    ', &
               &  neoU(l,(i-1)*8+5), '    ',  neoU(l,(i-1)*8+6), '    ',  neoU(l,(i-1)*8+7) 
       end do
                                

    
    else if (k.eq.6) then


        i=j+1
        write(739,'(I15,A,I15,A,I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4, '    ', &
          & (i-1)*8+5, '    ',(i-1)*8+6
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4), '    ',omegaorig((i-1)*8+5), '    ', omegaorig((i-1)*8+6)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4), '    ', &
               &  neoU(l,(i-1)*8+5), '    ',  neoU(l,(i-1)*8+6)
       end do
                                

    
    else if (k.eq.5) then

        i=j+1
        write(739,'(I15,A,I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4, '    ', &
          & (i-1)*8+5
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4), '    ',omegaorig((i-1)*8+5)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4), '    ', &
               &  neoU(l,(i-1)*8+5)
       end do ! l
    
    
    else if (k.eq.4) then

       i=j+1
        write(739,'(I15,A,I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3,'    ', (i-1)*8+4
       write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3),  &
               & '    ',omegaorig((i-1)*8+4)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3), '    ',  neoU(l,(i-1)*8+4)
       end do ! l
    
    
    else if (k.eq.3) then

       i=j+1
        write(739,'(I15,A,I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2, '    ',(i-1)*8+3
       write(739,'(G15.7,A,G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2),  '    ',omegaorig((i-1)*8+3)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2),  '    ', neoU(l,(i-1)*8+3)
       end do ! l
    
    
    else if (k.eq.2) then
 
        i=j+1
        write(739,'(I15,A,I15)') &
        &  (i-1)*8+1,'    ', (i-1)*8+2
       write(739,'(G15.7,A,G15.7)') & 
               & omegaorig((i-1)*8+1),  '    ',omegaorig((i-1)*8+2)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7,A,G15.7)') &
               &  neoU(l,(i-1)*8+1), '    ', neoU(l,(i-1)*8+2)
       end do ! l
    

    
    else if (k.eq.1) then
 
         i=j+1
        write(739,'(I15)') &
        &  (i-1)*8+1
       write(739,'(G15.7,A)') & 
               & omegaorig((i-1)*8+1)
       write(739,*)
       
       do l=1,dimnu
         write(739,'(G15.7)') &
               &  neoU(l,(i-1)*8+1)
       end do ! l
    
    
    end if

   
    close(739)
    
    deallocate(neoU)
     
    
 end subroutine  write_U_files
   
   
   


! PENG'S SUBROUTINES
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine: output()
!! This subroutine is used to generate
!! the input file for PWscf
!! Begin: 2011.05.23
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine output_file(total_num_atoms, Nspecies, atom_name, atom_num, vib_mode, atom_posi, n_mode)
! implicit none
!     !!I/O variables
!     integer  :: total_num_atoms, n_mode, Nspecies
!     integer  :: atom_num(Nspecies)
!     character(len=2) :: atom_name(Nspecies)
!     real(dp) :: vib_mode(3*total_num_atoms,3*total_num_atoms), atom_posi(total_num_atoms,3)
! 
!     !! local variables
!     integer :: ii, jj, kk
!     real(dp), allocatable :: atom_new(:,:)
! 
!     !! allocate the array
!     allocate(atom_new(total_num_atoms,3))
! 
!     if(n_mode == 0) then
!        atom_new = atom_posi
!     else 
!        do ii = 1, total_num_atoms
!           do jj = 1, 3
!              atom_new(ii,jj) = atom_posi(ii,jj) + vib_mode(n_mode,(ii-1)*3+jj)
!           end do
!        end do
!     end if
! 
!     open(90, file='position_nu.dat')
!     kk = 1
!     do ii = 1, Nspecies
!        do jj = 1, atom_num(ii)
! !          write(90, '(I2, 4X, 3(3XF12.8))') atom_name(ii), atom_new(kk,1), &
! !                      atom_new(kk,2), atom_new(kk,3)
!           write(90, *) atom_name(ii), atom_new(kk,1), &
!                       atom_new(kk,2), atom_new(kk,3)
! 
!          kk = kk+1
!        end do
!     end do
! 
!     deallocate(atom_new)
! 
!     return
! end subroutine output_file
! 
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! This subroutine is used to read and 
! !! normlize the eigen modes of nanocluster
! !! begin : 2011.04.09
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine dot_eigen(Input_coord, Input_dot, total_num_atoms, atom_position, eigenmode_dot, hw_dot)
! !use project_ph_mod
! implicit none
!    !! I/O variables
!    character (len=20) :: Input_coord, Input_dot
!    integer :: total_num_atoms
!    real(dp) :: atom_position(total_num_atoms,3)
!    real(dp) :: eigenmode_dot(3*total_num_atoms, 3*total_num_atoms)
!    real(dp) :: hw_dot(3*total_num_atoms)
! 
!    !! local variables
!   integer :: ii, jj
!   real(dp)  :: temp_re
! 
!   !! read the data
!   call read_data_dot(Input_coord, Input_dot, total_num_atoms, atom_position, eigenmode_dot, hw_dot)
! 
!   !! The eigenmodes read from CPMD don't need to be re-massed!
!   
!   !! Normalization
!   do ii = 1, 3*total_num_atoms !! For different modes
!      temp_re = 0.0
!      do jj = 1, 3*total_num_atoms  !! For different coorditotal_num_atomse 
!         temp_re = temp_re + eigenmode_dot(ii,jj)**2.0
!      end do
!      do jj = 1, 3*total_num_atoms
!         eigenmode_dot(ii,jj) = eigenmode_dot(ii,jj)/sqrt(temp_re)
!      end do
!   end do
! 
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! check  the orthogonalization!
! !  do ii = 1, 3*total_num_atoms
! !     do jj = 1, 3*total_num_atoms
! !        temp_re = dot_product(eigenmode_dot(ii,:), eigenmode_dot(jj,:)) !&
! ! !               /dot_product(eigenmode_dot(ii,:), eigenmode_dot(ii,:))
! !        write(*,*) 'Mode1=', ii, 'Mode2=',jj, 'orthod=', temp_re
! !     end do
! !  end do
! !stop
! 
!   return
! 
! end subroutine dot_eigen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine: Read data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine read_data_dot(Input_coord, Input_eigenmode, N_atom, Atom_posi, xi, Eigen_value)
! implicit none
!     integer :: N_atom
!     character (len=20) :: Input_coord, Input_eigenmode
!     real(dp) :: Atom_posi(N_atom,3)
!     real(dp) :: xi(3*N_atom, 3*N_atom)
!     real(dp) :: Eigen_value(3*N_atom)
! 
!     integer :: ii, jj, kk
!     integer :: N3, temp1, temp2
!     integer :: array_number(8)
!     integer :: status, ierror ! I/O status
!     character (len=2) :: temp_ch
!     character (len=20) :: temp_ch1
!     real(dp) :: array_temp(8)
!     logical :: exist_flag
! 
!     N3=3*N_atom
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! Read data from Coorditotal_num_atomse file
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ! Open the file and check for errors on open
! 
!    OPEN(UNIT=3, FILE=Input_coord, STATUS='OLD', ACTION='READ', IOSTAT=status)
!    IF (status /= 0) THEN
!        STOP
!        write(*,*) 'Error in reading coorditotal_num_atomse file'
!    END IF
!    WRITE(*,*) 'Start to read coorditotal_num_atomse file'
! 
!    do ii = 1, N_atom
!       read(3,*)  temp_ch, Atom_posi(ii,1), Atom_posi(ii,2), Atom_posi(ii,3)
!    end do
! 
!   CLOSE(UNIT = 3)
!    WRITE(*,*) 'File reading done'
! 
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! Read the eigenmode
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! Here, the Input file is in the format of CPMD
!    temp1 = int(real(N3)/8.0)
!    temp2 = mod(N3,8)
! 
!    OPEN(UNIT=3, FILE=Input_eigenmode, STATUS='OLD', ACTION='READ', IOSTAT=status)
!    IF (status /= 0) THEN
!        STOP
!        write(*,*) 'Error in reading coorditotal_num_atomse file'
!    END IF
!    WRITE(*,*) 'Start to read coorditotal_num_atomse file'
!    
!    write(*,*) 'Read the fake data'
!    do ii = 1, temp1
!       do jj = 1, 2
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6), array_temp(7), array_temp(8)      
!       end do
!       read(3,*)
! 
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6), array_temp(7), array_temp(8)
!       end do
!   end do
! !! Judge the mode
!    if(temp2 == 1) then
!       do jj = 1, 2
!          read(3,*)  array_temp(1)
!       end do
!       read(3,*)
! 
!       do jj = 1, N3
!          read(3,*)  array_temp(1)
!       end do
!    else if (temp2 == 2) then
!       do jj = 1, 2
!          read(3,*)  array_temp(1), array_temp(2)
!       end do
!       read(3,*)
! 
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2)
!       end do
!    else if (temp2 == 3) then
!       do jj = 1, 2
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3)
!       end do
!       read(3,*)
! 
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3)
!       end do
!   else if (temp2 == 4) then
!       do jj = 1, 2
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4)
!       end do
!       read(3,*)
! 
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4)
!       end do
!   else if (temp2 == 5) then
!       do jj = 1, 2
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5)
!       end do
!       read(3,*)
! 
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5)
!       end do
!   else if (temp2 == 6) then
!       do jj = 1, 2
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6)
!       end do
!       read(3,*)
! 
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6)
!       end do
!   else if (temp2 == 7) then
!       do jj = 1, 2
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6), array_temp(7)
!       end do
!       read(3,*)
! 
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6), array_temp(7)
!       end do
!   end if
!   
!   read (3,*) temp_ch1
!   write(*,*) temp_ch1
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   write(*,*) 'Read the real data!'
!    do ii = 1, temp1
!       read(3,*)  array_number(1), array_number(2), array_number(3), array_number(4), &
!                  array_number(5), array_number(6), array_number(7), array_number(8)
! !! Read the eigen value
!       read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                  array_temp(5), array_temp(6), array_temp(7), array_temp(8)
!       read(3,*)
! !! Write the eigen values
!       do kk = 1, 8
!          Eigen_value(array_number(kk))=array_temp(kk)
!       end do
! !! Read the eigen modes
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6), array_temp(7), array_temp(8)
! !! Write the eigen modes
!          do kk = 1, 8
!             xi(array_number(kk),jj) = array_temp(kk)
!          end do
!       end do
!   end do
! !! Judge the mode
!    if(temp2 == 1) then
!      !! read the eigen value
!       read(3,*)  array_number(1)
!       read(3,*)  array_temp(1)
!       read(3,*)
!      !! write the eigen value
!          Eigen_value(array_number(1)) = array_temp(1)
!      !! read & write the eigen mode
!       do jj = 1, N3
!          read(3,*)  array_temp(1)
!          Eigen_value(array_number(1)) = array_temp(1)
!       end do
!    else if (temp2 == 2) then
!      !! read the eigen value
!       read(3,*)  array_number(1), array_number(2)      
!       read(3,*)  array_temp(1), array_temp(2)
!       read(3,*)
!      !! write the eigen value
!       do kk = 1, temp2
!          Eigen_value(array_number(kk)) = array_temp(kk)
!       end do
!      !! read & write the eigen mode
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2)
!          do kk = 1, temp2
!             xi(array_number(kk),jj) = array_temp(kk)
!          end do
!       end do
!    else if (temp2 == 3) then
!       !! read the eigen value
!       read(3,*)  array_number(1), array_number(2), array_number(3)
!       read(3,*)  array_temp(1), array_temp(2), array_temp(3)
!       read(3,*)
!       !! write the eigen value
!       do kk = 1, temp2
!          Eigen_value(array_number(kk)) = array_temp(kk)
!       end do
!       !! read & write the eigen mode
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3)
!          do kk = 1, temp2
!             xi(array_number(kk),jj) = array_temp(kk)
!          end do
!       end do
!   else if (temp2 == 4) then
!       !! read the eigen value
!       read(3,*)  array_number(1), array_number(2), array_number(3), array_number(4)
!       read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4)
!       read(3,*)
!       !! write the eigen value
!       do kk = 1, temp2
!          Eigen_value(array_number(kk)) = array_temp(kk)
!       end do
!       !! read & write the eigen mode
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4)
!          do kk = 1, temp2
!             xi(array_number(kk),jj) = array_temp(kk)
!          end do
!       end do
!   else if (temp2 == 5) then
!       !! read the eigen vale
!       read(3,*)  array_number(1), array_number(2), array_number(3), array_number(4), &
!                  array_number(5)
!       read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                  array_temp(5)
!       read(3,*)
!       !! write the eigen value
!       do kk = 1, temp2
!          Eigen_value(array_number(kk)) = array_temp(kk)
!       end do
!       !! read & write the eigen mode
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5)
!          do kk = 1, temp2
!             xi(array_number(kk),jj) = array_temp(kk)
!          end do
!       end do
!   else if (temp2 == 6) then
!       !! read the eigen value
!       read(3,*)  array_number(1), array_number(2), array_number(3), array_number(4), &
!                  array_number(5), array_number(6)
!       read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                  array_temp(5), array_temp(6)
!       read(3,*)
!       !! write the eigen value
!       do kk = 1, temp2
!          Eigen_value(array_number(kk)) = array_temp(kk)
!       end do
!       !! read & write the eigen mode
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6)
!          do kk = 1, temp2
!             xi(array_number(kk),jj) = array_temp(kk)
!          end do
!       end do
!   else if (temp2 == 7) then
!       !! read the eigen value
!       read(3,*)  array_number(1), array_number(2), array_number(3), array_number(4), &
!                  array_number(5), array_number(6), array_number(7)
!       read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                  array_temp(5), array_temp(6), array_temp(7)
!       read(3,*)
!       !! write the eigen value
!       do kk = 1, temp2
!          Eigen_value(array_number(kk)) = array_temp(kk)
!       end do
!       !! read & write the eigen mode
!       do jj = 1, N3
!          read(3,*)  array_temp(1), array_temp(2), array_temp(3), array_temp(4), &
!                     array_temp(5), array_temp(6), array_temp(7)
!          do kk = 1, temp2
!             xi(array_number(kk),jj) = array_temp(kk)
!          end do
!       end do
!   end if 
! 
!   CLOSE(UNIT = 3)
!    WRITE(*,*) 'File reading done'
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   return
! end subroutine read_data_dot


!-------------------------------------------------------------------------------------------    
  
 subroutine read_QE_masses (number_of_atoms,  number_of_species, QE_masses_file_name, mass)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives the atomic masses.
 
 ! Comment: Ponce2014 considers hbar = m_e = e = 1. This means that, since we use the same convention,
 ! we must multiply the masses expressed with the unified atomic mass pattern by 1822.888486192
 ! Also Peng used mass in atomic units; in his original program vib_pot.f90 one can see:
 ! atom_mass(ii) = atom_mass(ii)*P_mass2au ! P_mass2au=1822.84
 

   ! System parameters
  integer, Intent(In) :: number_of_atoms                ! Total number of atoms of the system
  integer, Intent(In) :: number_of_species              ! Total number of different atomic species of the system
  character(100), Intent(In) :: QE_masses_file_name     ! Name of the file where the electronic eigenvalues are stored (it must be an output of an 'scf' calculation of Quantum Espresso)
  real(8), Intent(Out) :: mass(number_of_atoms,2)       ! Masses of atoms
  
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character(12), allocatable :: species_name(:)
  real(8), allocatable :: species_mass(:)
  character(12) :: auxstr4, auxstr5
  character(18) :: auxstr6, auxstr7
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, HOMO, LUMO, fermi_level  
  
  
  inquire(file=QE_masses_file_name, exist=exist) 
  if ( .not. exist) then
     write(*,*)
     write(*,*) '  ERROR: The ', trim(QE_masses_file_name), ' specified in vib_potPGR.in '
     write(*,*) '  does not exist. Please, provide it.'
     write(*,*) 
     call exit(1)
  end if
        
  write(auxstr3,*) trim(QE_masses_file_name) 
  call StripSpaces(auxstr3)      
  write(*,*) '  *** Now reading the output of Quantum Espresso'
  write(*,*) '  (', trim(auxstr3),') to find atomic masses. ***'                    


  allocate(species_name(number_of_species),species_mass(number_of_species))
  
  mass(:,2) = -1.0d0
  
 open(333,file=QE_masses_file_name,status='unknown')

  do k=1,20000

      read(333,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 7456

      auxstr1=''; auxstr2='' 
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( auxstr1 .eq. 'atomicspeciesvalencemass') then      
          go to 2435
        end if
      end do ! i
      
   end do !k   
   
 2435 continue  
  
  
  do i=1,number_of_species
    read(333,*) species_name(i), auxr1, species_mass(i), auxstr1
    call StripSpaces(species_name(i))
  end do

         
   do k=1,20000

      read(333,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 7456

      auxstr1=''; auxstr2='' 
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( auxstr1 .eq. 'siten.atompositions') then  
          go to 1435
        end if
      end do ! i
      
   end do !k   
   
 1435 continue   
 

  do i=1,number_of_atoms
    read(333,*) j, auxstr4
    call StripSpaces(auxstr4)
    do k=1,number_of_species
       if ( auxstr4 .eq. species_name(k)) then
         mass(i,1) = dble(i)       
         mass(i,2) = species_mass(k) * 1822.888486192d0  
       end if
    end do
  end do
  
  
  do i=1,number_of_atoms
     if (mass(i,2) .lt. 0.d0) then
        write(*,*) '  ERROR: Masses are incorrect or they were not properly read.'
        call exit(1)
     end if
  end do
  
  
  close(333)
  
 
 7456 continue
  
  deallocate(species_name)


 end subroutine read_QE_masses
 

!-------------------------------------------------------------------------------------------

