program frozenphononPGR

!     This program (frozenphonon) calculates renormalized electronic eigenvalues a la Capaz (PRL 2005).
!     The input files for this program are the output files from Quantum Espresso (QE); see the end of this file for examples.
!     If you want to calculate the HOMO and the LUMO in paralel, submit 2 runs, the first one with
!     the value of FinalStateFF equal to the highest (maybe degenerate) HOMO index and the second one with
!     the value of InitialStateFF equal to the lowest (maybe degenerate) LUMO index.
!
!     Author:   Pablo Garcia Risueno (2017).
!
!     An example of the input file for this program (frozenphonon.in) is given below:
! 
!      Temperature = 0.0 K
!  
!      Number_of_atoms = 26
!      Number_of_species = 2
! 
! 	   # The variable below is the number multiplied by U (the eigenvector of the dynamical
! 	   # matrix) to perform the displacements used to do the finite-difference.
! 	   Finite_Difference_Parameter = 2.d0
! 
! 	   # The variable below is the name of the output file of a gs calculation of Quantum Espresso
! 	   # (pw.x) storing the electronic eigenvalues of the undisplaced position
! 	   elec_eigval_file_name = out_scf0.out
! 
! 	   Band_occupation = 2
!      
!      # The <Check_crossover> variable below specifies whether or not the possible crossovers in eigenvalues due to
!      # the finite-difference displacements are checked and avoided or not. Doing it makes the
!      # results much more accurate, but it demands the calculation and printing to file of 
!      # wavefunctions, which makes calculations slower and demands much disk-storage capability.
!      # The <Crossovers_from_file> variable reads the information on crossovers from files (generated previously
!      # by making this variable equal to 1). This is useful to avoid the repetition of the slow crossover
!      # calculation when evaluating other temperatures.
!      # The <InitialStateFF> and <FinalStateFF> variables make it possible to specify the range
!      # of wavefunctions to look for the crossovers (these two variables are optional, a default of 5 states 
!      # above the LUMO and below the HOMO is taken as default). 
!      # If you want to calculate the HOMO and the LUMO in paralel, submit 2 runs, the first one with
!      # the value of FinalStateFF equal to the highest (maybe degenerate) HOMO index and the second one with
!      # the value of InitialStateFF equal to the lowest (maybe degenerate) LUMO index.
! 
!      Check_crossover = 0
!      Crossovers_from_file = 0
!      InitialStateFPh = 24
!      FinalStateFPh = 30
!
!
!      # The two variables below are optional. They specify which phonon modes are considered in the calculations
!      # (which can be useful to parallelize the calculations, or to neglect some modes considered not interesting). 
!
!      InitialMode = 7
!      FinalMode = 25
!
! 	   # The variable below is the name of the output file of a calculation of Quantum Espresso
! 	   # (pp.x) to solve the dynamical matrix. It stores phonon frequencies and the U eigenvectors,
! 	   # which are normalized following the eq. (38) of Ponce et al. PRB 90, 214304 (2014).
!      # Check its format at the example file.
!
! 	   phonon_frequency_file_name = Uvec.dat
! 	   phonon_frequency_units=cm-1
! 
! 	   # The variables below are the names of the folders storing the outputs of Quantum Espresso
! 	   # (pw.x) for gs calculations with nuclear positions displaced with respect to the relaxation
! 	   # (displaced + and - U eigenvectors times <<Finite_Difference_Parameter>>). Note that the  
! 	   # names of the files in these folders MUST be "out_scf-modeXXX.out".
!
! 	   folder_output_files_displaced+ = outputsVeffcalc/displaced+
! 	   folder_output_files_displaced- = outputsVeffcalc/displaced-
!
!     # The variables below make it possible to impose degeneracy to calculate the renormalizations by hand; this can be useful for nearly degenerate states.
!       HOMO_imposed_degeneracy = 2
!       LUMO_imposed_degeneracy = 2
!
!  ------------------------------------------------------------------------------------------------------
!
!   To make a sweep of different temperatures, add lines like the following ones to the frozenphonon.in file:
!
!   Initial_Temperature = 0.0001d0 K
!   Final_Temperature = 350.0 K
!   Number_of_temperatures = 4
!
!
!  ======================================================================================================
!  ======================================================================================================
!  ======================================================================================================


  implicit none
  integer :: sizeepsilon, number_of_atoms, number_of_atomic_species, dimnu, Ne, band_occupation,i
  integer :: iband1, iband2, nat, n_atom_type, nfft, InitialMode, FinalMode, number_of_Ts
  integer :: InitialStateFF, FinalStateFF, HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f,HOMO_deg, LUMO_deg
  real(8) :: E_Fermi, displ_parameter, lattparam, prod, T, InitialT, FinalT
  real(8),allocatable ::  U(:,:), omega(:), omegaorig(:) 
  real(8), allocatable :: epsilon0orig(:), epsilon0(:)
  real(8), allocatable ::  epsilon_displ_plus0(:,:),  epsilon_displ_minus0(:,:)
  real(8), allocatable ::  epsilon_displ_plus(:,:),  epsilon_displ_minus(:,:)
  character(100) :: elec_eigval_file_name, xml_eigval_file_name, phon_freq_file_name  
  character(100) :: folder_output_files_displaced_plus, folder_output_files_displaced_minus
  character(len=100) :: Input_wfc_prefix, Input_wfc(2), auxstr1, auxstr2, wfc_directory
  character(40)  :: units_omegaorig, units_epsilon0orig
  logical :: check_crossover, crossoversfromfile, eigvalsfromxml
  
  wfc_directory='Wfcs/';  Input_wfc_prefix='wfc'
     
  call read_input_file (T, InitialT, FinalT, number_of_Ts, &
                         & number_of_atoms, number_of_atomic_species, dimnu, InitialMode, FinalMode, band_occupation, &
                         &  check_crossover, crossoversfromfile, eigvalsfromxml, InitialStateFF, FinalStateFF, displ_parameter,&
                         & HOMO_deg, LUMO_deg,units_omegaorig, elec_eigval_file_name, xml_eigval_file_name, phon_freq_file_name, &
                         &  folder_output_files_displaced_plus, folder_output_files_displaced_minus)
  
  call read_number_of_bands (elec_eigval_file_name, band_occupation, Ne, sizeepsilon)       
                                   
  allocate(omega(dimnu),omegaorig(dimnu),U(dimnu,dimnu))
  allocate(epsilon0(sizeepsilon),epsilon0orig(sizeepsilon))
  allocate(epsilon_displ_plus0(sizeepsilon,dimnu),epsilon_displ_minus0(sizeepsilon,dimnu))
  allocate(epsilon_displ_plus(sizeepsilon,dimnu),epsilon_displ_minus(sizeepsilon,dimnu))
  epsilon0(:)=0.d0; omega(:)=0.d0; epsilon0orig(:)=0.d0;  omegaorig(:)=0.d0
  epsilon_displ_plus0(:,:)=0.d0; epsilon_displ_minus0(:,:)=0.d0
  epsilon_displ_plus(:,:)=0.d0; epsilon_displ_minus(:,:)=0.d0
              
  ! We read the bands of the gs, i.e. of the system with undisplaced positions  
  call find_Fermi_level (eigvalsfromxml, sizeepsilon, Ne, elec_eigval_file_name, xml_eigval_file_name,  E_Fermi)            
  call read_bands_from_QE_output (eigvalsfromxml, sizeepsilon, Ne, elec_eigval_file_name, xml_eigval_file_name, &
                                  &  E_Fermi, epsilon0orig,  units_epsilon0orig)                

  ! We read the output of the dynamical equation (i.e. of the phonon calculation).
  call read_dyneq_output (number_of_atoms, dimnu, phon_freq_file_name, omegaorig, U)   

  ! We read the output files of the DFT calculations with + and - displaced positions of the nuclei.
  call read_QE_displaced_output (eigvalsfromxml, sizeepsilon,Ne,dimnu,folder_output_files_displaced_plus, &
                                & E_Fermi, epsilon_displ_plus0)         
  call read_QE_displaced_output (eigvalsfromxml, sizeepsilon,Ne,dimnu,folder_output_files_displaced_minus, &
                                & E_Fermi, epsilon_displ_minus0)         


  ! We write the read variables in atomic units.
  call omega_epsilon_to_atomic_units (eigvalsfromxml, units_omegaorig,units_epsilon0orig,dimnu,sizeepsilon, &        
                       & omegaorig,epsilon0orig,omega,epsilon0, epsilon_displ_plus0, epsilon_displ_minus0)

 
  ! We find the indices of the HOMO and LUMO states (to calculate their renormalization due to frozen-phonon calculations)
  call find_HOMO_LUMO_indices ( Ne, band_occupation,sizeepsilon,epsilon0, HOMO_deg, LUMO_deg, &
                                  & InitialStateFF, FinalStateFF,  &
                                  &  HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f)
  
  
  ! If required, we check the possible crossovers of eigenvalues due to the displacements. This makes the calculations much
  ! more accurate, but it demands the evaluation of wavefunctions, which is very disk-memory consuming. 
  if ( (check_crossover) .and. (.not.crossoversfromfile)) then
	! We read the variable which gives the size of the wavefunctions (in reciprocal space)                                
	call read_nfft(HOMO_index_f, number_of_atoms, number_of_atomic_species, wfc_directory, Input_wfc_prefix, nfft) 
  end if	                               

  if  (check_crossover) then
  
	! We find the (appropriate eigenvalues with displaced positions) which correspond to the eigenvalues with unperturbed
	! positions. The finite-difference displacement may have resulted into a crossover of them, which leads to completely wrong results.
	call remove_crossover_in_eigenvalues (crossoversfromfile, sizeepsilon, dimnu, number_of_atoms,  &
		  & number_of_atomic_species, nfft, InitialStateFF, FinalStateFF, InitialMode, FinalMode, &
		  & HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f, &
		  & wfc_directory, Input_wfc_prefix, epsilon0, &
		  & epsilon_displ_plus0, epsilon_displ_minus0, &
		  & epsilon_displ_plus, epsilon_displ_minus)
		  
  else	! (check_crossover)	  
  
     epsilon_displ_minus(:,:) = epsilon_displ_minus0(:,:)    
     epsilon_displ_plus(:,:) = epsilon_displ_plus0(:,:)  
     
  end if      
  
  
   
  ! With the collected data, we finally calculate the renormalization of HOMO and LUMO due to electron-phonon interaction.
  call calculate_frozen_phonon_renormalization ( T, InitialT, FinalT, number_of_Ts, &
                                               & sizeepsilon, dimnu, InitialStateFF, FinalStateFF , &
  											   &  InitialMode, FinalMode, &
   											   &  HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f, &	
                                               &  displ_parameter, omega, epsilon0, &
                                               &  epsilon_displ_plus, epsilon_displ_minus)


  deallocate(omega,omegaorig,U,epsilon0,epsilon0orig,epsilon_displ_plus,epsilon_displ_minus)                       
  

end program frozenphononPGR


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------


subroutine read_input_file (T, InitialT, FinalT, number_of_Ts, &
                         & number_of_atoms, number_of_atomic_species, dimnu, InitialMode, FinalMode,  &
						 &  band_occupation, check_crossover, &
                         & crossoversfromfile, eigvalsfromxml, InitialStateFF, FinalStateFF,  displ_parameter, &
                         & HOMO_deg, LUMO_deg, units_omegaorig, &
                         &  elec_eigval_file_name, xml_eigval_file_name, phon_freq_file_name, &
                         & folder_output_files_displaced_plus,folder_output_files_displaced_minus)

 ! This subroutine reads the frozenphonon.in file, which contains the names of the files where one must read the data
 ! as well as other relevant data. 

  ! System parameters
  real(8), Intent(Out) :: T                              ! Temperature of the system
  real(8), Intent(Out) :: InitialT, FinalT               ! Optional variables; if present, they are the initial and final temperatures for a sweep in T
  integer, Intent(Out) :: number_of_Ts                   ! Optional variable; if present, it indicates the number of temperatures considered in the sweep in T
  integer, Intent(Out) :: number_of_atoms                ! The number of atoms of the system
  integer, Intent(Out) :: number_of_atomic_species       ! The number of atomic species (different chemical elements) of the system  
  integer, Intent(Out) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  integer, Intent(Out) :: InitialMode, FinalMode         ! Initial and final phonon modes to be considered in the calculations (they are useful to parallelize the calculations)
  integer, Intent(Out) :: band_occupation                ! Number of electrons per band of the ones given in the output file of Quantum Espresso; it can be 1 or 2
  logical, Intent(Out) :: check_crossover				 ! If true, the possible crossovers of eigenvalues due to the displacements are checked; this makes calculations more accurate, but also more disk-memory consuming.
  logical, Intent(Out) :: crossoversfromfile             ! If true, the crossovers are read from file (which is useful to avoid the calculation of crossovers for every T)
  logical, Intent(Out) :: eigvalsfromxml                 ! If true (default) the eigenvalues are read from .xml files, containing more digits in the output of QE
  integer, Intent(Out) :: InitialStateFF, FinalStateFF   ! Initial and final state indices for the calculation of frozen-phonon where the possible crossover is sought.
  real(8), Intent(Out) :: displ_parameter                ! Number multiplied by U (the eigenvector of the dynamical matrix) to perform the displacements used to do the finite-difference.
  integer, Intent(Out) :: HOMO_deg, LUMO_deg             ! Degeneracies of HOMO and LUMO, imposed by hand by the human user in frozenphonon.in (useful for states nearly-degenerate)
  character(40) :: units_omegaorig						 ! units of the frequencies stored in the phon_freq_file_name file (Input_dot)
  character(100), Intent(Out) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored (apart from other informations)
  character(100), Intent(Out) :: xml_eigval_file_name   ! Name of the xml file where JUST the electronic eigenvalues are stored (this is used just if eigvalsfromxml is true)
  character(100), Intent(Out) :: phon_freq_file_name     ! Name of the file where the phonon frequencies and eigenvectors of the dynamical matrix are stored
  character(100), Intent(Out) ::  folder_output_files_displaced_plus
  character(100), Intent(Out) ::  folder_output_files_displaced_minus 
  
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT, auxstr4, auxstr5
  character :: letter 
  integer :: stat,auxi,i,j,k,l,l2,l3,iostat,unitt !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: auxr1, auxr2


  T = -10.d0; InitialT=-10.d0; FinalT=-10.d0; number_of_Ts=-10
  InitialMode=-10; FinalMode = - 10
  number_of_atoms = 0; number_of_atomic_species = 0; units_omegaorig='cm**-1'
  InitialStateFF=0; FinalStateFF=0
  gFan_i_index_beg=1; gFan_i_index_end=1; gFan_j_index_beg=1; gFan_j_index_end=1
  gDW_i_index_beg=0; gDW_i_index_end=0; band_occupation=2
  folder_output_files_displaced_plus='./';folder_output_files_displaced_minus='./'
  displ_parameter = 1.d0; HOMO_deg=0; LUMO_deg=0
  check_crossover=.false.; crossoversfromfile=.false.; eigvalsfromxml=.true.; l=0; l2=0
  
  
  !gFan_file_name='gFan.dat'

  ! Reading parameters from file
  
   write(*,*); write(*,*); 
   write(*,*) '  ***************************************************************************' 
   write(*,*) '  ***************************  NOW RUNNING   ********************************'   
   write(*,*) '  ************************** frozen-phonon.x  *******************************'       
   write(*,*) '  ***************************************************************************'   
   write(*,*) ;  write(*,*)
   write(*,*) '  ********* FROM THE INPUT FILE OF frozen-phonon.x  (frozenphonon.in): ********'   


   inquire(file='frozenphonon.in', exist=exist)      
   if (.not. exist) then
	   write(*,*) 
	   write(*,*) '  ***************************************************************************' 
	   write(*,*) '    ERROR: The file "frozenphonon.in" does not exist. Please, provide it.'
	   write(*,*) '  ***************************************************************************' 
	   write(*,*)
	   call exit(1) 
   end if


   write(*,*)

   open(345,file='frozenphonon.in',status='unknown')

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




    if ((auxstr1 .eq. 'temperature').or.(auxstr1 .eq. 'Temperature'))   then
       Read(auxstr2, * )  

       do i=1,100
         letter = auxstr2(i:i)
         if ( (letter .eq. 'm') .or. (letter .eq. 'H') .or. (letter .eq. 'h') .or. &
            & (letter .eq. 'K') .or. (letter .eq. 'k') ) then          
           stringT =    auxstr2(1:i-1)
           Tunits = trim(auxstr2(i:i+6))
           go to 9997
         end if  
       end do
           write(*,*) 
           write(*,*) '************************************************************************************'
           write(*,*) '    ERROR: Unrecognised temperature units. Please, write K, Hartree or meV'
           write(*,*) '    after the temperature in your "input.in" file. Note that the units of the '
           write(*,*) '    temperature must be the same as the units of the energy (eigenvalues). '
           write(*,*) '************************************************************************************'
           write(*,*) 
           call exit(1)
  9997 continue

       Read(stringT, '(f10.0)' ) T        
     !  write(*,*) 'T= ',T,'  unitsT= ',Tunits
           
       if ( (Tunits .eq. 'meV') .or. (Tunits .eq. 'mev') .or. (Tunits .eq. 'MEV') ) then
          unitt=1
          write(*,'(A,f9.3,A,f8.2,A)')  '               Temperature = ',T/0.08621738,' K (', T, ' meV )'; write(*,*)
          T=T/(27211.3860217)
       else if ( (Tunits .eq. 'Hartree') .or. (Tunits .eq. 'hartree') .or. (Tunits .eq. 'H') ) then
          unitt=2
          write(*,'(A,f9.3,A,f14.7,A)') '               Temperature = ',T*27211.3860217/0.08621738,' K (', T, ' Hartree )'
          write(*,*)
       else if ( (Tunits .eq. 'k') .or. (Tunits .eq. 'K') .or. (Tunits .eq. 'kelvin').or. (Tunits .eq. 'Kelvin') ) then
          unitt=3
          write(*,'(A,f9.3,A,f8.2,A)')  '                                               Temperature = ',T,' K '
          write(*,*)
          T=T*0.08621738/27211.3860217
       else
          write(*,*) 
          write(*,*) '**********************************************************************************'
          write(*,*) '    ERROR: Unrecognised temperature units. Please, write K, Hartree or meV'
          write(*,*) '    after the temperature in your "input.in" file. Note that the units of the '
          write(*,*) '    temperature must be the same as the units of the energy (eigenvalues). '
          write(*,*) '**********************************************************************************'
          write(*,*) 
          call exit(1)
       end if
 
    end if   ! Reading of temperature
    
    
    if ((auxstr1 .eq. 'Initial_Temperature').or.(auxstr1 .eq. 'Initial_temperature') &
           & .or.(auxstr1 .eq. 'initial_temperature') .or.(auxstr1 .eq. 'initialtemperature') )   then
       Read(auxstr2, * )  

       do i=1,100
         letter = auxstr2(i:i)
         if ( (letter .eq. 'm') .or. (letter .eq. 'H') .or. (letter .eq. 'h') .or. &
            & (letter .eq. 'K') .or. (letter .eq. 'k') ) then          
           stringT =    auxstr2(1:i-1)
           Tunits = trim(auxstr2(i:i+6))
           go to 9998
         end if  
       end do
           write(*,*) 
           write(*,*) '************************************************************************************'
           write(*,*) '    ERROR: Unrecognised temperature units. Please, write K, Hartree or meV'
           write(*,*) '    after the temperature in your "input.in" file. Note that the units of the '
           write(*,*) '    temperature must be the same as the units of the energy (eigenvalues). '
           write(*,*) '************************************************************************************'
           write(*,*) 
           call exit(1)
  9998 continue

       Read(stringT, '(f10.0)' ) InitialT        
           
       if ( (Tunits .eq. 'meV') .or. (Tunits .eq. 'mev') .or. (Tunits .eq. 'MEV') ) then
          unitt=1
          write(*,'(A,f9.3,A,f8.2,A)')  '  Initial T for T-sweeping = ',InitialT/0.08621738,' K (', InitialT , ' meV )'
          InitialT=InitialT/(27211.3860217)
       else if ( (Tunits .eq. 'Hartree') .or. (Tunits .eq. 'hartree') .or. (Tunits .eq. 'H') ) then
          unitt=2
          write(*,'(A,f9.3,A,f14.7,A)') '  Initial T for T-sweeping = ',InitialT *27211.3860217/0.08621738,&
                                                   & ' K (', InitialT , ' Hartree )'
       else if ( (Tunits .eq. 'k') .or. (Tunits .eq. 'K') .or. (Tunits .eq. 'kelvin').or. (Tunits .eq. 'Kelvin') ) then
          unitt=3
          write(*,'(A,f9.3,A,f8.2,A)')  '                                  Initial T for T-sweeping = ',InitialT,' K '
          InitialT=InitialT*0.08621738/27211.3860217
       else
          write(*,*) 
          write(*,*) '**********************************************************************************'
          write(*,*) '    ERROR: Unrecognised temperature units. Please, write K, Hartree or meV'
          write(*,*) '    after the temperature in your "input.in" file. Note that the units of the '
          write(*,*) '    temperature must be the same as the units of the energy (eigenvalues). '
          write(*,*) '**********************************************************************************'
          write(*,*) 
          call exit(1)
       end if
 
    end if   ! Reading of initial temperature   
    
    
   if ((auxstr1 .eq. 'Final_Temperature').or.(auxstr1 .eq. 'Final_temperature') &
           & .or.(auxstr1 .eq. 'Final_temperature') .or.(auxstr1 .eq. 'finaltemperature') )   then
       Read(auxstr2, * )  

       do i=1,100
         letter = auxstr2(i:i)
         if ( (letter .eq. 'm') .or. (letter .eq. 'H') .or. (letter .eq. 'h') .or. &
            & (letter .eq. 'K') .or. (letter .eq. 'k') ) then          
           stringT =    auxstr2(1:i-1)
           Tunits = trim(auxstr2(i:i+6))
           go to 9999
         end if  
       end do
           write(*,*) 
           write(*,*) '************************************************************************************'
           write(*,*) '    ERROR: Unrecognised temperature units. Please, write K, Hartree or meV'
           write(*,*) '    after the temperature in your "input.in" file. Note that the units of the '
           write(*,*) '    temperature must be the same as the units of the energy (eigenvalues). '
           write(*,*) '************************************************************************************'
           write(*,*) 
           call exit(1)
  9999 continue

       Read(stringT, '(f10.0)' ) FinalT        
           
       if ( (Tunits .eq. 'meV') .or. (Tunits .eq. 'mev') .or. (Tunits .eq. 'MEV') ) then
          unitt=1
          write(*,'(A,f9.3,A,f8.2,A)')  '  Final T for T-sweeping = ',FinalT/0.08621738,' K (', FinalT , ' meV )'
          write(*,*)
          FinalT=FinalT/(27211.3860217)
       else if ( (Tunits .eq. 'Hartree') .or. (Tunits .eq. 'hartree') .or. (Tunits .eq. 'H') ) then
          unitt=2
          write(*,'(A,f9.3,A,f14.7,A)') '  Final T for T-sweeping = ',FinalT *27211.3860217/0.08621738,&
                                                   & ' K (', FinalT , ' Hartree )'
          write(*,*)
       else if ( (Tunits .eq. 'k') .or. (Tunits .eq. 'K') .or. (Tunits .eq. 'kelvin').or. (Tunits .eq. 'Kelvin') ) then
          unitt=3
          write(*,'(A,f9.3,A,f8.2,A)')  '                                    Final T for T-sweeping = ',FinalT,' K '
          write(*,*)
          FinalT=FinalT*0.08621738/27211.3860217
       else
          write(*,*) 
          write(*,*) '**********************************************************************************'
          write(*,*) '    ERROR: Unrecognised temperature units. Please, write K, Hartree or meV'
          write(*,*) '    after the temperature in your "input.in" file. Note that the units of the '
          write(*,*) '    temperature must be the same as the units of the energy (eigenvalues). '
          write(*,*) '**********************************************************************************'
          write(*,*) 
          call exit(1)
       end if
 
    end if   ! Reading of final temperature   
     

     if ( T .lt. 0.d0) then                                          ! |  ---------------------------------------------------------------------
         T = 0.0000000001d0                                       !     |  ---------------------------------------------------------------------
       ! write(*,'(A,f9.3,A,f8.2,A)')  '                                               Temperature = 0 K ' !----------------------
      end if                                                      !     |  ---------------------------------------------------------------------
                 

     if ((auxstr1 .eq. 'elec_eigval_file_name').or.(auxstr1 .eq. 'Elec_eigval_file_name'))   then
         elec_eigval_file_name = trim(auxstr2)
         call StripSpaces (elec_eigval_file_name)
     end if

     if ((auxstr1 .eq. 'xml_eigval_file_name').or.(auxstr1 .eq. 'Xml_eigval_file_name'))   then
         xml_eigval_file_name = trim(auxstr2)
         call StripSpaces (xml_eigval_file_name)
     end if

     if ((auxstr1 .eq. 'phonon_frequency_file_name').or.(auxstr1 .eq. 'Phonon_frequency_file_name') &
      & .or.(auxstr1 .eq. 'Phonon_Frequency_File_Name') )   then
         phon_freq_file_name = trim(auxstr2)
         call StripSpaces (phon_freq_file_name)
     end if


     if ((auxstr1 .eq. 'folder_output_files_displaced+').or.(auxstr1 .eq. 'Folder_output_files_displaced+') &
      & .or.(auxstr1 .eq. 'Folder_Output_Files_Displaced+') )   then
         folder_output_files_displaced_plus = trim(auxstr2)
         call StripSpaces (folder_output_files_displaced_plus)
     end if

     if ((auxstr1 .eq. 'folder_output_files_displaced-').or.(auxstr1 .eq. 'Folder_output_files_displaced-') &
      & .or.(auxstr1 .eq. 'Folder_Output_Files_Displaced-') )   then
         folder_output_files_displaced_minus = trim(auxstr2)
         call StripSpaces (folder_output_files_displaced_minus)
     end if 
     
! JUST FOR CPMD INPUT     
!     if ((auxstr1 .eq. 'mass_file_name').or.(auxstr1 .eq. 'Mass_file_name'))   then
!          mass_file_name = trim(auxstr2)
!          call StripSpaces (mass_file_name)
!      end if
     
    if ((auxstr1 .eq. 'Omega_units').or.(auxstr1 .eq. 'omega_units').or.(auxstr1 .eq. 'Omega_Units') .or. &
        & (auxstr1 .eq. 'phonon_frequency_units').or.(auxstr1 .eq. 'Phonon_frequency_units').or. &
        &   (auxstr1 .eq. 'phonon_Frequency_Units')  )   then
         units_omegaorig = trim(auxstr2)
         call StripSpaces (units_omegaorig)
     end if         

   
    if ((auxstr1 .eq. 'Number_of_atoms').or.(auxstr1 .eq. 'number_of_atoms').or. &
             & (auxstr1 .eq. 'Number_atoms'))  Read(auxstr2, '(I3)' ) number_of_atoms 
             
    if ((auxstr1 .eq. 'Number_of_Ts').or.(auxstr1 .eq. 'number_of_Ts').or. &
             & (auxstr1 .eq. 'Number_of_temperatures') .or. (auxstr1 .eq. 'number_of_temperatures')) then
       Read(auxstr2, '(I3)' ) number_of_Ts        
       write(*,'(A,I6)') ' The number of temperatures for T-sweeping is ', number_of_Ts
    end if   
             
    if ((auxstr1 .eq. 'Check_crossover').or.(auxstr1 .eq. 'Check_Crossover').or. &
             & (auxstr1 .eq. 'CheckCrossover') .or.(auxstr1 .eq. 'Check_Crossovers').or. &
             & (auxstr1 .eq. 'Check_crossovers').or.(auxstr1 .eq. 'CheckCrossovers') )  Read(auxstr2, '(I3)' ) l             

    if ((auxstr1 .eq. 'Crossovers_From_File').or.(auxstr1 .eq. 'crossovers_from_file').or. &
             & (auxstr1 .eq. 'Crossovers_from_file') .or.(auxstr1 .eq. 'Crossovers_from_File').or. &
             & (auxstr1 .eq. 'Crossovers_file').or.(auxstr1 .eq. 'Crossoversfromfile') )  Read(auxstr2, '(I3)' ) l2             


    if ((auxstr1 .eq. 'Eigenvalues_from_xml_files').or.(auxstr1 .eq. 'eigenvalues_from_xml_files').or. &
             & (auxstr1 .eq. 'Eigenvalues_from_xml') .or.(auxstr1 .eq. 'Eigenvalues_From_Xml_Files').or. &
             & (auxstr1 .eq. 'Eigenvalues_from_Xml').or.(auxstr1 .eq. 'eigenvalues_from_xml') )  Read(auxstr2, '(I3)' ) l3    
                 
             
    if ((auxstr1 .eq. 'Number_of_atomic_species').or.(auxstr1 .eq. 'number_of_atomic_species').or. &
             & (auxstr1 .eq. 'Number_atomic_species') .or.(auxstr1 .eq. 'number_of_species').or. &
             & (auxstr1 .eq. 'Number_of_species'))  Read(auxstr2, '(I3)' ) number_of_atomic_species 

    if ((auxstr1 .eq. 'InitialFrequency').or.(auxstr1 .eq. 'Initialfrequency').or. &
             & (auxstr1 .eq. 'initialfrequency') .or. (auxstr1 .eq. 'InitialModeuency') &
             & .or.(auxstr1 .eq. 'InitialMode').or.(auxstr1 .eq. 'Initialmode').or. &
             & (auxstr1 .eq. 'initialmode'))  Read(auxstr2, '(I3)' ) InitialMode            
             
    if ((auxstr1 .eq. 'FinalFrequency').or.(auxstr1 .eq. 'Finalfrequency').or. &
             & (auxstr1 .eq. 'finalfrequency') .or. (auxstr1 .eq. 'finalmode') &
             & .or.(auxstr1 .eq. 'FinalMode').or.(auxstr1 .eq. 'Finalmode').or. &
             & (auxstr1 .eq. 'final_mode'))  Read(auxstr2, '(I3)' ) FinalMode   
 
     if ((auxstr1 .eq. 'InitialStateFPh').or.(auxstr1 .eq. 'InitialStatefph').or. &
             & (auxstr1 .eq. 'InitialStateFPh'))  Read(auxstr2, '(I3)' ) InitialStateFF
             
    if ((auxstr1 .eq. 'FinalStatefph').or.(auxstr1 .eq. 'FinalStateFPh').or. &
             & (auxstr1 .eq. 'FinalStatefph'))  Read(auxstr2, '(I3)' ) FinalStateFF                        

    if ((auxstr1 .eq. 'band_occupation').or.(auxstr1 .eq. 'Band_Occupation').or. &
             & (auxstr1 .eq. 'Band_occupation'))  Read(auxstr2, '(I3)' ) band_occupation

    if ((auxstr1 .eq. 'HOMO_degeneracy').or.(auxstr1 .eq. 'HOMO_Imposed_Degeneracy').or. &
             & (auxstr1 .eq. 'HOMO_imposed_degeneracy'))  Read(auxstr2, '(I3)' ) HOMO_deg

    if ((auxstr1 .eq. 'LUMO_degeneracy').or.(auxstr1 .eq. 'LUMO_Imposed_Degeneracy').or. &
             & (auxstr1 .eq. 'LUMO_imposed_degeneracy'))  Read(auxstr2, '(I3)' ) LUMO_deg

    if ((auxstr1 .eq. 'Finite_difference_parameter').or.(auxstr1 .eq. 'finite_difference_parameter').or. &
             & (auxstr1 .eq. 'Finite_Difference_Parameter'))  Read(auxstr2, '(G18.10)' ) displ_parameter             
              
   
  end do ! k
   

1117 continue



   if ( ( (InitialT .ge. 0.d0 ) .and. (FinalT .lt.0.d0) ) .or. ( (InitialT .lt. 0.d0 ) .and. (FinalT .ge.0.d0) ) ) then
	   write(*,*) 
	   write(*,*) '  ***************************************************************************'  
	   write(*,*) '    ERROR: Please, define both the initial and the final temperatures for '
	   write(*,*) '    The temperature-sweeping.'
	   write(*,*) '    (e.g. with:" '
	   write(*,*) '     Initial_Temperature = 10.0 K  '
	   write(*,*) '     Final_Temperature = 300.0 K '
	   write(*,*) '     " in the frozenphonon.in file)  '
	   write(*,*) '  ***************************************************************************' 
	   write(*,*)
	   call exit(1)    
   end if

   if ( ( (InitialT .ge. 0.d0 ) .and. (FinalT .ge.0.d0)  .and. (number_of_Ts .le. 0) ) )  then
	   write(*,*) 
	   write(*,*) '  ***************************************************************************'  
	   write(*,*) '    ERROR: Please, define a number of temperatures for the T-sweeping.'
	   write(*,*) '    (e.g. " Number_of_Temperatures = 10  " in the frozenphonon.in file)  '
	   write(*,*) '  ***************************************************************************' 
	   write(*,*)
	   call exit(1)    
   end if
  
  
  
!    if ( (number_of_Ts .gt. 0) .and. (check_crossovers) .and. (.not.crossoversfromfile)  )  then
! 	   write(*,*) 
! 	   write(*,*) '  ***************************************************************************'  
! 	   write(*,*) '    WARNING: You plan to do a sweeping in temperatures checking crossovers'
! 	   write(*,*) '    but you do not plan to read the crossover information from file. This'
! 	   write(*,*) '    could be extremely slow. We recommend you to generate the files with '
! 	   write(*,*) '    crossover info at ONE given temperature, and then to sweep in temperatures'
! 	   write(*,*) '    using << Crossovers_from_file = 1 >> in the frozenphonon.in file.'
! 	   write(*,*) '  ***************************************************************************' 
! 	   write(*,*)
!    end if  
!   
  
  
   inquire(file=phon_freq_file_name, exist=exist)      
   if (.not. exist) then
	   write(*,*) 
	   write(*,*) '  ***************************************************************************'  
	   write(*,*) '    ERROR: The file ',phon_freq_file_name,' with phonon frequencies does not exist.'
	   write(*,*) '    Please, provide it.'
	   write(*,*) '  ***************************************************************************' 
	   write(*,*)
	   call exit(1) 
   end if
   
   inquire(file=elec_eigval_file_name, exist=exist)      
   if (.not. exist) then
	   write(*,*) 
	   write(*,*) '  ***************************************************************************'  
	   write(*,*) '    ERROR: The file ',elec_eigval_file_name,' with electronic eigenvalues and     '
	   write(*,*) '    atomic masses does not exist. Please, provide it.'
	   write(*,*) '  ***************************************************************************' 
	   write(*,*)
	   call exit(1) 
   end if    

   if ( eigvalsfromxml ) then
	 inquire(file=xml_eigval_file_name, exist=exist)      
	 if (.not. exist) then
		 write(*,*) 
		 write(*,*) '  ***************************************************************************'  
		 write(*,*) '    ERROR: The xml file ',xml_eigval_file_name,' with electronic eigenvalues and     '
		 write(*,*) '    atomic masses does not exist. Please, provide it.'
		 write(*,*) '  ***************************************************************************' 
		 write(*,*)
		 call exit(1) 
	 end if 
   end if 

  if ((band_occupation .ne. 1) .and. (band_occupation .ne. 2)) then
    band_occupation = 2
    write(*,*)
    write(*,*) '  <<<<<<<<<<<< WARNING: The number of electrons per band in ',trim(elec_eigval_file_name)
    write(*,*) '  was assumed to be 2. >>>>>>>>>>>'
    write(*,*)   
  end if
  
   
  dimnu = number_of_atoms*3 ! Strictly speaking it is number_of_atoms*3-6, but we make it bigger to avoid problems of lacks of space (some freqs are 0 and they are anyway read)
 
  if ( InitialMode .le. 0 ) then
     InitialMode = 7
  end if 
  if ( ( FinalMode .lt. 0) .or. ( FinalMode .gt. dimnu ) ) then
     FinalMode = dimnu
  end if

 
  if ( ( trim(units_omegaorig) .eq. 'CM**-1' ) .or. ( trim(units_omegaorig) .eq. 'Cm**-1' ) ) units_omegaorig='cm**-1' 
  if ( ( trim(units_omegaorig) .eq. 'mev' )    .or. ( trim(units_omegaorig) .eq. 'MeV' ) )    units_omegaorig='meV'
  if ( ( trim(units_omegaorig) .eq. 'ev' )    .or. ( trim(units_omegaorig) .eq. 'EV' ) )    units_omegaorig='eV'
  if ( ( trim(units_omegaorig) .eq. 'A.u.' ) .or. ( trim(units_omegaorig) .eq. 'A.u.' ) .or. &
    &  ( trim(units_omegaorig) .eq. 'Au.' )  .or. ( trim(units_omegaorig) .eq. 'AU' ) .or.   &
    & ( trim(units_omegaorig) .eq. 'hartree' ) .or. ( trim(units_omegaorig) .eq. 'HARTREE' )  &
    & .or.  ( trim(units_omegaorig) .eq. 'Hartree' ) )                                          units_omegaorig='a.u.'
  
     
  call StripSpaces(units_omegaorig)
  if ((  ( (trim(units_omegaorig).ne.'cm**-1') .and. (trim(units_omegaorig).ne.'cm-1') ) ) &
   & .and. ( trim(units_omegaorig) .ne. 'meV' ) &
   & .and. ( trim(units_omegaorig) .ne. 'a.u.' ).and. ( trim(units_omegaorig) .ne. 'eV' ) ) then
  write(*,*)
    write(*,*) '  ERROR: Unrecognised units for the phonon frequencies (units_omega variable '
    write(*,*) '         in frozenphonon.in). Accepted values are cm**-1, eV, meV and  '
    write(*,*) '         a.u. (Hartree). Please set the variable to one of these values.  '
    call exit(1)
  end if
 

  if (number_of_atoms .le. 0) then
     write(*,*) ' ERROR: Please, specify a valid number of atoms in frozenphonon.in; '
     write(*,*) '        E.g.: Number_of_atoms =  2 '
     call exit(1)
  end if

  if (number_of_atomic_species .le. 0) then
     write(*,*) ' ERROR: Please, specify a valid number of atomic species in frozenphonon.in; '
     write(*,*) '        E.g.: Number_of_atomic_species =  2 '
     call exit(1)
  end if
  
    
 
  write(*,*) '  File of phonon frequencies and eigvectors of the dynamical matr.:  ',trim(phon_freq_file_name)
  write(*,*) '   with freqs. in ',trim(units_omegaorig)
  write(*,*) '  Name of the file that stores electronic eigenvalues and masses:  ',trim(elec_eigval_file_name)  

  if (band_occupation .eq. 1) then
    write(*,*) '   with up to 1 electron per band.'
  else  if(band_occupation .eq. 2) then
    write(*,*) '   with up to 2 electrons per band.'
  else  
    write(*,*) '  ERROR: Unphysical number of electrons per band (',band_occupation,').'; write(*,*)
    call exit(1)     
  end if     

  if(eigvalsfromxml) write(*,*)' Name of the xml file that stores precise electronic eigenvalues:  ',trim(xml_eigval_file_name)  

  write(*,*) '  Name of the folder that stores displaced (+) outputs:  ',trim(folder_output_files_displaced_plus)
  write(*,*) '  Name of the folder that stores displaced (-) outputs:  ',trim(folder_output_files_displaced_minus)

  write(*,*)
  write(*,'(A,I5)') '                                           Number of atoms = ', number_of_atoms
  write(*,'(A,I5)') '                                  Number of atomic species = ', number_of_atomic_species
  write(auxstr4,*) dimnu-6; call StripSpaces(auxstr4)
  write(auxstr5,*) dimnu; call StripSpaces(auxstr5) 
  write(*,'(A,G18.9)') '                  Finite-difference displacement parameter = ', displ_parameter
  write(*,*)        '                                Number of phonon branches =    ', &
           & trim(auxstr5), ' (',trim(auxstr4),' valid)'
  if (  (InitialMode .gt. 7) .or. ( FinalMode .lt. dimnu ) ) then         
	write(auxstr4,*) InitialMode; call StripSpaces(auxstr4)
	write(auxstr5,*)   FinalMode; call StripSpaces(auxstr5)
	write(*,*)  '         Just the contribution of phonon modes between ',trim(auxstr4),'   and ',trim(auxstr5)
	write(*,*)  '         will be considered.'; write(*,*)
  end if	
  
  write(auxstr4,*) gDW_i_index_beg; call StripSpaces(auxstr4)
  write(auxstr5,*) gDW_i_index_end; call StripSpaces(auxstr5)


  if ( l .gt. 0 ) then
    check_crossover = .true.
    write(*,*) ; write(*,*) '      Possible crossovers in the eigenvalues due to the finite-difference '
    write(*,*) '      displacements will be checked and avoided.' 
  else
    write(*,*) ; write(*,*) '      Possible crossovers in the eigenvalues due to the finite- '
    write(*,*) '      difference displacements will not be checked.'; write(*,*)     
  end if 

  if ( l2 .gt. 0 ) then
    crossoversfromfile = .true.
    write(*,*) '      The possible crossovers in the eigenvalues due to the finite-difference '
    write(*,*) '      displacements will be read from file.'; write(*,*)     
  end if 

  if ( l3 .le. 0 ) then
    eigvalsfromxml  = .false.
    write(*,*) '      The eigenvalues will not be read from the .xml files output of QE, but from '
    write(*,*) '      QE standard output; note that this contains less digits, which can lead to'
    write(*,*) '      lower accuracy.' ; write(*,*)     
  end if 


  inquire(file=phon_freq_file_name, exist=exist)      ! This checks whether the file gDW.dat is existing.
  if (.not. exist) then
      write(*,*) 
      write(*,*) '  ***************************************************************************'  
      write(*,*) '    ERROR: The file ',trim(phon_freq_file_name),' does not exist. Please, provide it.'
      write(*,*) '  ***************************************************************************' 
      write(*,*)
      call exit(1) 
  end if
  
  inquire(file=elec_eigval_file_name, exist=exist)      ! This checks whether the file gDW.dat is existing.
  if (.not. exist) then
      write(*,*) 
      write(*,*) '  ***************************************************************************'  
      write(*,*) '    ERROR: The file ',trim(elec_eigval_file_name),' does not exist. Please, provide it.'
      write(*,*) '  ***************************************************************************' 
      write(*,*)
      call exit(1) 
  end if

 
end subroutine read_input_file


!-------------------------------------------------------------------------------------------
      
 subroutine read_number_of_bands (elec_eigval_file_name, band_occupation, Ne, sizeepsilon)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives the numbers of occ and unocc states.

   ! System parameters
  character(100), Intent(In) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  integer, Intent(In)  :: band_occupation               ! 1 or 2, number of electrons per band in the QE output file
  integer, Intent(Out) :: Ne                            ! The number of electrons of the system (i.e. the number of occupied states)
  integer, Intent(Out) :: sizeepsilon                   ! size of the vector storing the electronic eigenvalues

     
    ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character :: letter 
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount , Nbands
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9  
  
        
  write(*,*); write(*,*) '  ************ FROM THE QUANTUM ESPRESSO OUTPUT (', trim(elec_eigval_file_name),'): *************'
  write(*,*)    



  open(433,file=elec_eigval_file_name,status='unknown')

  do m=1,20000
  
   do k=1,2000

      read(433,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 2117

      auxstr1=''; auxstr2=''; auxstr3=''
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( (auxstr1 .eq. 'numberofelectrons=') .or.  ( auxstr1 .eq. 'Numberofelectrons=')  ) then  
          go to 2235
        end if
      end do ! i
      
   end do !k   
   
 2235 continue      

   do j=i+1,i+18
     letter = inputline(j:j)
     auxstr2 = trim(auxstr2)//letter
   end do
   read(auxstr2,*) auxr1
   Ne=nint(auxr1)
  ! write(*,'(A,I5)') '                                       Number of electrons = ', Ne

  close(433)
  
end do ! m
  
 2117 continue 
  
  close(433)
  
  
  open(451,file=elec_eigval_file_name,status='unknown')

  do m=1,20000
  
   do k=1,2000

      read(451,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 3117

      auxstr1=''; auxstr2=''; auxstr3=''
      
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( (auxstr1 .eq. 'numberofKohn-Shamstates=') .or.  ( auxstr1 .eq. 'NumberofKohn-Shamstates=')  ) then  
          go to 2237
        end if
      end do ! i
      
   end do !k   
   
 2237 continue      
      
   do j=i+1,i+18
     letter = inputline(j:j)
     auxstr2 = trim(auxstr2)//letter
   end do
   
   read(auxstr2,*) Nbands
  write(*,'(A,I5)') '                                       Number of electrons = ', Ne
  write(*,'(A,I5)') '                                 Number of occupied states = ', Ne/band_occupation + mod(Ne,band_occupation)
  write(*,'(A,I5)') '               Total number of considered Kohn-Sham states = ', Nbands
  write(*,*)
      
  close(451)
  
end do ! m
  
  3117 continue
  
  if (band_occupation .eq. 2) then
     sizeepsilon = Nbands   !2*Nbands
  else if (band_occupation .eq. 1) then 
     sizeepsilon = Nbands
  else   
      write(*,*) '  ERROR: Unphysical number of electrons per band'
      call exit(1)
  end if

 end subroutine read_number_of_bands
  

!-------------------------------------------------------------------------------------------    
  
 subroutine read_bands_from_QE_output (eigvalsfromxml, dimband, Ne, elec_eigval_file_name, xml_eigval_file_name, & 
                                      & fermi_level , epsilon_output, units_epsilon_output)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives electronic eigenvalues and atomic masses.
 
 ! Comment: Ponce2014 considers hbar = m_e = e = 1 (atomic units). This means that, since we use the same convention,
 ! we must multiply the masses expressed with the unified atomic mass pattern by 1822.888486192
 
 

   ! System parameters
  logical, Intent(In) :: eigvalsfromxml                 ! True if the eigenvalues are read from the .xml output of QE (which gives eigenvalues in atomic units already) 
  integer, Intent(In) :: dimband                        ! Number of considered electronic eigenvalues
  integer, Intent(In) :: Ne                             ! The number of electrons of the system (i.e. the number of occupied states)
  character(100), Intent(In) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  character(100), Intent(In) :: xml_eigval_file_name    ! Name of the xml file where the electronic eigenvalues are stored
  real(8), Intent(In) ::  fermi_level 
  real(8), Intent(Out) :: epsilon_output(dimband)				! Electronic eigenvalues
  character(40), Intent(Out) :: units_epsilon_output          ! Units of the electronic eigenvalues

  
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character(12), allocatable :: species_name(:)
  real(8), allocatable :: species_mass(:)
  character(12) :: auxstr4, auxstr5
  character(18) :: auxstr6, auxstr7
  character :: letter 
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist
  real(8) :: Nbands,auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, HOMO, LUMO 
  
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																									!
!                             READING EIGENVALUES										!
!																									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We count the number of bands

!number of Kohn-Sham states=


  epsilon_output(:)=0.d0
  
  if ( .not. eigvalsfromxml ) then

	  open(333,file=elec_eigval_file_name,status='unknown')

	  do k=1,20000

		  read(333,'(A)',IOSTAT=stat) inputline

		  if (stat .ne. 0) go to 3456

		  auxstr1=''; auxstr2='' 
	  
		  do i=1,90
			letter = inputline(i:i)
			auxstr1 = trim(auxstr1)//letter 
			if ( auxstr1 .eq. 'Endofself-consistentcalculation') then  
			  go to 2235
			end if
		  end do ! i
	  
	   end do !k   
   
	 2235 continue    
 
 
	  read(333,'(A)',IOSTAT=stat) inputline;  
	  read(333,'(A)',IOSTAT=stat) inputline 
	  auxstr1=''
	  do i=1,90
		letter = inputline(i:i)
		auxstr1 = trim(auxstr1)//letter 
		if ( auxstr1 .eq. 'convergenceNOT') then  
			write(*,*)
			write(*,*) '   ERROR: The calculations summarized in '
			write(*,*) '   ',trim(elec_eigval_file_name)
			write(*,*) '   did not converge. Please, rerun Quantum Espresso (perhaps with lower conv_thr).'; write(*,*)
			call exit(1) 
		end if
	  end do ! i
  
 
	  auxstr1=''
	  do i=1,60
		 letter = inputline(i:i)
		 auxstr1 = trim(auxstr1)//letter 
		 auxstr2=''
		 j=LNBLNK(auxstr1)   ! Last index of the chain that is not a blank space
		 if ( j.ge.5) then
		   do k=1,5
			  letter=inputline(i-6+k:i-6+k)
			  auxstr2 = trim(auxstr2)//letter
		   end do
		 end if
	 
		 if ( auxstr2 .eq. 'bands') then 
			  units_epsilon_output=''
			  do j=i+1,i+7
				  letter=inputline(j:j)
				  if ( (letter .ne. '(') .and.(letter .ne. ')') .and.(letter .ne. ',') .and. &
					& (letter .ne. '.') .and.(letter .ne. ';') .and.(letter .ne. ':')  ) then
					 units_epsilon_output = trim(units_epsilon_output)//letter
				  end if
			  end do
			  if ( ( trim(units_epsilon_output) .eq. 'CM**-1' ) .or. ( trim(units_epsilon_output) .eq. 'Cm**-1' ) ) &
				   &  units_epsilon_output='cm**-1' 
			  if ( ( trim(units_epsilon_output) .eq. 'mev' )    .or. ( trim(units_epsilon_output) .eq. 'MeV' ) )   &
				   &  units_epsilon_output='meV'
			  if ( ( trim(units_epsilon_output) .eq. 'ev' )  .or. ( trim(units_epsilon_output) .eq. 'eV' ) &
			  .or. ( trim(units_epsilon_output) .eq. 'Ev' )     .or. ( trim(units_epsilon_output) .eq. 'EV' ) )   &
				 &   units_epsilon_output='eV'
			  if ( ( trim(units_epsilon_output) .eq. 'A.u.' ) .or. ( trim(units_epsilon_output) .eq. 'A.u.' ) .or. &
			   &  ( trim(units_epsilon_output) .eq. 'Au.' )  .or. ( trim(units_epsilon_output) .eq. 'AU' ) .or.   &
				& ( trim(units_epsilon_output) .eq. 'hartree' ) .or. ( trim(units_epsilon_output) .eq. 'HARTREE' )  &
				& .or.  ( trim(units_epsilon_output) .eq. 'Hartree' ) )                                         &
					 &  units_epsilon_output='a.u.'

		 end if
	  end do ! i       

	  read(333,'(A)',IOSTAT=stat) inputline
		
	  auxstr2='' 
	  i=1             ! i is the band counter
	  do k=1,20000 
		read(333,'(A)',IOSTAT=stat) inputline  
	
		do j=1,100
	 
		 if ((stat .ne. 0)  .or. (i.gt.dimband)) go to 3456
		 letter = inputline(j:j)
		 if (letter .eq. 'h') go to 3456
		 auxstr2 = trim(auxstr2)//letter
		 if ( (auxstr2 .ne. '') .and. (letter .eq. ' ') ) then
			read(auxstr2,*) auxr1
			epsilon_output(i) = auxr1
			i = i+1         
			auxstr2=''
		 end if
	   end do ! j  
	  end do   ! k
	 
		
	3456 continue

	  close(333)
	  
	  

   else !  eigvalsfromxml=.true.	
   
   
	
      open(333,file=xml_eigval_file_name,status='unknown')

	  do k=1,20000

		  read(333,'(A)',IOSTAT=stat) inputline

		  if (stat .ne. 0) go to 6456

		  auxstr1=''; auxstr2='' 
  
		  do i=1,90
			letter = inputline(i:i)
			auxstr1 = trim(auxstr1)//letter 
			if ( auxstr1 .eq. '<EIGENVALUES') then  
			  go to 6235
			end if
		  end do ! i
	  
		end do !k   
   
		 6235 continue    
 

 
         do j=1,dimband
            !read(333,'(A)',IOSTAT=stat) inputline
            !write(*,*) inputline
		    read(333,'(G30.23)') epsilon_output(j)
         end do
	 
	    6456 continue
	    
	    close(333)

   end if ! eigvalsfromxml


  ! Now we make the 0 of energy in the Fermi level; we establish the Fermi level in the medium point of the band gap
  do i=1,dimband
    epsilon_output(i) = epsilon_output(i) - fermi_level
  end do 


 end subroutine read_bands_from_QE_output
 

!-------------------------------------------------------------------------------------------    
  
 subroutine find_Fermi_level (eigvalsfromxml, dimband, Ne, elec_eigval_file_name, xml_eigval_file_name,  fermi_level)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives electronic eigenvalues and atomic masses.
 
 ! Comment: Ponce2014 considers hbar = m_e = e = 1 (atomic units). This means that, since we use the same convention,
 ! we must multiply the masses expressed with the unified atomic mass pattern by 1822.888486192
 
 

   ! System parameters
  logical, Intent(In) :: eigvalsfromxml                 ! True if the eigenvalues are read from the .xml output of QE (which gives eigenvalues in atomic units already) 
  integer, Intent(In) :: dimband                        ! Number of considered electronic eigenvalues
  integer, Intent(In) :: Ne                             ! The number of electrons of the system (i.e. the number of occupied states)
  character(100), Intent(In) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  character(100), Intent(In) :: xml_eigval_file_name    ! Name of the xml file where the electronic eigenvalues are stored (with more digits)
  real(8), Intent(Out) ::  fermi_level 
 
  
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character(12), allocatable :: species_name(:)
  real(8), allocatable :: species_mass(:)
  character(12) :: auxstr4, auxstr5
  character(18) :: auxstr6, auxstr7
  character :: letter 
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist
  real(8) :: Nbands,auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, HOMO, LUMO 
  real(8), allocatable :: epsilon_output(:)
  character(40) :: units_epsilon_output 
  
  allocate(epsilon_output(dimband))
  
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																									!
!                             READING EIGENVALUES										!
!																									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We count the number of bands

!number of Kohn-Sham states=


  epsilon_output(:)=0.d0

  if ( .not. eigvalsfromxml) then
  
      open(333,file=elec_eigval_file_name,status='unknown')
  
	  do k=1,20000

		  read(333,'(A)',IOSTAT=stat) inputline

		  if (stat .ne. 0) go to 3456

		  auxstr1=''; auxstr2='' 
	  
		  do i=1,90
			letter = inputline(i:i)
			auxstr1 = trim(auxstr1)//letter 
			if ( auxstr1 .eq. 'Endofself-consistentcalculation') then  
			  go to 2235
			end if
		  end do ! i
	  
	   end do !k   
   
	 2235 continue    
 
 
	  read(333,'(A)',IOSTAT=stat) inputline;   read(333,'(A)',IOSTAT=stat) inputline 
 
	  auxstr1=''
	  do i=1,60
		 letter = inputline(i:i)
		 auxstr1 = trim(auxstr1)//letter 
		 auxstr2=''
		 j=LNBLNK(auxstr1)   ! Last index of the chain that is not a blank space
		 if ( j.ge.5) then
		   do k=1,5
			  letter=inputline(i-6+k:i-6+k)
			  auxstr2 = trim(auxstr2)//letter
		   end do
		 end if
	 
		 if ( auxstr2 .eq. 'bands') then 
			  units_epsilon_output=''
			  do j=i+1,i+7
				  letter=inputline(j:j)
				  if ( (letter .ne. '(') .and.(letter .ne. ')') .and.(letter .ne. ',') .and. &
					& (letter .ne. '.') .and.(letter .ne. ';') .and.(letter .ne. ':')  ) then
					 units_epsilon_output = trim(units_epsilon_output)//letter
				  end if
			  end do
			  if ( ( trim(units_epsilon_output) .eq. 'CM**-1' ) .or. ( trim(units_epsilon_output) .eq. 'Cm**-1' ) ) &
				   &  units_epsilon_output='cm**-1' 
			  if ( ( trim(units_epsilon_output) .eq. 'mev' )    .or. ( trim(units_epsilon_output) .eq. 'MeV' ) )   &
				   &  units_epsilon_output='meV'
			  if ( ( trim(units_epsilon_output) .eq. 'ev' )  .or. ( trim(units_epsilon_output) .eq. 'eV' ) &
			  .or. ( trim(units_epsilon_output) .eq. 'Ev' )     .or. ( trim(units_epsilon_output) .eq. 'EV' ) )   &
				 &   units_epsilon_output='eV'
			  if ( ( trim(units_epsilon_output) .eq. 'A.u.' ) .or. ( trim(units_epsilon_output) .eq. 'A.u.' ) .or. &
			   &  ( trim(units_epsilon_output) .eq. 'Au.' )  .or. ( trim(units_epsilon_output) .eq. 'AU' ) .or.   &
				& ( trim(units_epsilon_output) .eq. 'hartree' ) .or. ( trim(units_epsilon_output) .eq. 'HARTREE' )  &
				& .or.  ( trim(units_epsilon_output) .eq. 'Hartree' ) )                                         &
					 &  units_epsilon_output='a.u.'

		 end if
	  end do ! i       

	  read(333,'(A)',IOSTAT=stat) inputline
		
	  auxstr2='' 
	  i=1             ! i is the band counter
	  do k=1,20000 
		read(333,'(A)',IOSTAT=stat) inputline  
	
		do j=1,100
	 
		 if ((stat .ne. 0)  .or. (i.gt.dimband)) go to 3456
		 letter = inputline(j:j)
		 if (letter .eq. 'h') go to 3456
		 auxstr2 = trim(auxstr2)//letter
		 if ( (auxstr2 .ne. '') .and. (letter .eq. ' ') ) then
			read(auxstr2,*) auxr1
			epsilon_output(i) = auxr1
			i = i+1         
			auxstr2=''
		 end if
	   end do ! j  
	  end do   ! k
	 
		
	3456 continue
	
	close(333)
	
   else !  eigvalsfromxml=.true.	
	
      open(333,file=xml_eigval_file_name,status='unknown')

	  do k=1,20000

		  read(333,'(A)',IOSTAT=stat) inputline

		  if (stat .ne. 0) go to 6456

		  auxstr1=''; auxstr2='' 
  
		  do i=1,90
			letter = inputline(i:i)
			auxstr1 = trim(auxstr1)//letter 
			if ( auxstr1 .eq. '<EIGENVALUES') then  
			  go to 6235
			end if
		  end do ! i
	  
		end do !k   
   
		 6235 continue    
 

 
         do j=1,dimband
            !read(333,'(A)',IOSTAT=stat) inputline
            !write(*,*) inputline
		    read(333,'(G30.23)') epsilon_output(j)
         end do
	 
	    6456 continue
	    
	    close(333)

   end if ! eigvalsfromxml



  ! Now we make the 0 of energy in the Fermi level; we establish the Fermi level in the medium point of the band gap
  ! Note that if the eigenvalues are read from xml, the Fermi level is in Hartree, but if read from standard output, it will be in meV 
  if ( mod(Ne,2) .eq. 1) then
     HOMO = epsilon_output(Ne);  LUMO = HOMO
  else
     HOMO = epsilon_output(Ne/2);  LUMO = epsilon_output(Ne/2+1)
  end if
  fermi_level = (HOMO+LUMO)/2.d0
  do i=1,dimband
    epsilon_output(i) = epsilon_output(i) - fermi_level
  end do 



  deallocate(epsilon_output)

 end subroutine find_Fermi_level
 
!-------------------------------------------------------------------------------------------    
  
 subroutine read_QE_displaced_output (eigvalsfromxml,sizeepsilon,Ne,dimnu,folder_output_files_displaced, &
                                     & E_Fermi, epsilon_displ)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives electronic eigenvalues
 
 ! Comment: Ponce2014 considers hbar = m_e = e = 1 (atomic units). This means that, since we use the same convention,
 ! we must multiply the masses expressed with the unified atomic mass pattern by 1822.888486192
 

   ! System parameters
  logical, Intent(In) :: eigvalsfromxml                 ! True if the eigenvalues are read from the .xml output of QE (which gives eigenvalues in atomic units already) 
  integer, Intent(In) :: sizeepsilon                        ! Number of considered electronic eigenvalues
  integer, Intent(In) :: Ne                             ! The number of electrons of the system (i.e. the number of occupied states)
  integer, Intent(In) :: dimnu								! Number of phonon branches
  character(100), Intent(In) :: folder_output_files_displaced  ! Name of the folder where the files to read lie
  real(8), Intent(In) :: E_Fermi 
  real(8), Intent(Out) :: epsilon_displ(sizeepsilon,dimnu)  ! Read eigenvalues for the displaced positions
 
   
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT, auxstrA
  character(12), allocatable :: species_name(:)
  real(8), allocatable :: species_mass(:)
  character(12) :: auxstr4, auxstr5
  character(18) :: auxstr6, auxstr7
  character :: letter 
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount , nu
  logical :: exist
  real(8) :: Nbands,auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, HOMO, LUMO, fermi_level  
  character(100) :: elec_eigval_file_name, xml_eigval_file_name 
  character(40)  :: units_epsilon_displ
  real(8), allocatable :: epsilon_aux(:)
  
  allocate(epsilon_aux(sizeepsilon))
   
        
   write(*,*) '  *** Now reading the displaced output of Quantum Espresso ***'                    

   do nu=7,dimnu

      epsilon_aux(:)=0.d0
      
      
	  auxstrA=''; write(auxstrA,*) nu,".out";  call StripSpaces(auxstrA)
	  write(elec_eigval_file_name ,*) trim(folder_output_files_displaced)//'/'//'out_scf-mode'//auxstrA
	  call StripSpaces(elec_eigval_file_name)
	  
	  auxstrA=''; write(auxstrA,*) nu,".xml";  call StripSpaces(auxstrA)
	  write(xml_eigval_file_name ,*) trim(folder_output_files_displaced)//'/'//'eigenval-mode'//auxstrA
	  call StripSpaces(xml_eigval_file_name)
	   
      call read_bands_from_QE_output (eigvalsfromxml, sizeepsilon, Ne, elec_eigval_file_name, &
                      & xml_eigval_file_name, E_Fermi, epsilon_aux, units_epsilon_displ)                

      do i=1,sizeepsilon
        epsilon_displ(i,nu) = epsilon_aux(i)
      end do

   end do
  
  
  deallocate(epsilon_aux)
  
 end subroutine read_QE_displaced_output

!-------------------------------------------------------------------------------------------    

       
subroutine read_dyneq_output (number_of_atoms, dimnu, phon_freq_file_name,  omegaorig,  U)     
 
 ! This subroutine reads the file which contains the phonon frequencies and the eigenvectors of the dynamical equation.
 ! The format of the read file must be that of CPMD; we use Quantum Espresso instead, but our vib_potPGR program transforms
 ! the QE format for these data to CPMD format. Note that the output of QE stores the U vectors, while that of CPMD 
 ! stores the xi vectors (U=xi/sqrt{M}). The CPMD output file is called VIVEIGVEC; for QE it is user-defined.
 ! Since our system of analysis is an isolated system, we assume that the negligible acoustic modes correspond to the 1st six freqs 
 ! appearing in the file to read (if the system was periodic, just the first 3 ones would be negligible).

! IMPORTANT: In the version of PGR of April 2017 of vib_potPGR.f90, the Uvec.dat and the Input_dot files are written 
! in atomic units, NOT in Angstrom (as it was the case of Input_dot of the original version of vib_pot.f90 by Peng).
! Hence this means that we do NOT need to transform U from angstrom to a.u.: i.e., the line
! Uau(i1,i2) =  U(i1,i2)/P_au2A 
! must be replaced by Uau(i1,i2) =  U(i1,i2).


  ! System parameters
  integer, Intent(In) :: number_of_atoms                ! The number of atoms of the system
  integer, Intent(In) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  character(100), Intent(In) :: phon_freq_file_name     ! Name of the file where the phonon frequencies and eigenvectors of the dynamical matrix are stored
 ! real(8), Intent(In)  :: mass(number_of_atoms,2)       ! Masses of the atoms (in the 2nd entry) 
  real(8), Intent(Out) :: omegaorig(dimnu)              ! Phonon frequencies (as they are read from file, with the units provided by the original file)
  real(8), Intent(Out) :: U(dimnu,dimnu)				! Eigenvectors of the dynamical matrix (U)
 
 
   ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character :: letter 
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist
  integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10 
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, auxr10, P_au2A

 
  write(*,*) '  *** Now reading the output of QE (in CPMD format, '&
              & , trim(phon_freq_file_name),' file) to find ' 
  write(*,*)  '             phonon frequencies and eigenvectors of the dynamical matrix. ***'
     


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!     BLOCK FOR READING PHONON FREQUENCIES AND EIGENVECTORS OF THE DYNAMICAL MATRIX        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  omegaorig(:) = 0.d0
  U(:,:)=0.d0

  
  open(341,file=phon_freq_file_name,status='unknown')
  
! In the Uvec.dat there are two blocks with two iterations with the data. We move to the point of the file where the 2nd block begins  
  
  do m=1,dimnu+1
    do k=1,(dimnu+20)**2
      read(341,'(A)',IOSTAT=stat) inputline
      if (stat .ne. 0) go to 1235
      auxstr1=''; auxstr2=''; auxstr3=''
      do i=1,90
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( auxstr1 .eq. '<<<<<<NEWDATA') then  
          go to 1235
        end if
      end do ! i
     end do !k   
   end do ! m  
1235 continue  
  
  
  
  
  do i=1,(ceiling(dble(dimnu)/8.d0)-1)  
    read(341,*) i1, i2, i3, i4, i5, i6, i7, i8 
    read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8
        
    omegaorig(8*(i-1)+1) = auxr1; omegaorig(8*(i-1)+2) = auxr2; omegaorig(8*(i-1)+3) = auxr3;  omegaorig(8*(i-1)+4) = auxr4
    omegaorig(8*(i-1)+5) = auxr5; omegaorig(8*(i-1)+6) = auxr6; omegaorig(8*(i-1)+7) = auxr7;  omegaorig(8*(i-1)+8) = auxr8
     
    read(341,'(A)',IOSTAT=stat) inputline 
    
    
    do j=1,dimnu  ! This way to read means that the U vector for a given I lies in a COLUMN of the U matrix.
        read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8

        U(j,8*(i-1)+1) = auxr1; U(j,8*(i-1)+2) = auxr2; U(j,8*(i-1)+3) = auxr3;  U(j,8*(i-1)+4) = auxr4
        U(j,8*(i-1)+5) = auxr5; U(j,8*(i-1)+6) = auxr6; U(j,8*(i-1)+7) = auxr7;  U(j,8*(i-1)+8) = auxr8
    end do
    
  end do ! i
  

  
  ! Reading of the last block of Uvec.dat (which does not contain 8 vectors, but a lower number)
  
  if ( mod(dimnu,8) .eq. 0) then
  
    read(341,*) i1, i2, i3, i4, i5, i6, i7, i8 
    read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8
    omegaorig(8*(i-1)+1) = auxr1; omegaorig(8*(i-1)+2) = auxr2; omegaorig(8*(i-1)+3) = auxr3;  omegaorig(8*(i-1)+4) = auxr4
    omegaorig(8*(i-1)+5) = auxr5; omegaorig(8*(i-1)+6) = auxr6; omegaorig(8*(i-1)+7) = auxr7;  omegaorig(8*(i-1)+8) = auxr8
    read(341,'(A)',IOSTAT=stat) inputline 
    do j=1,dimnu
        read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8
        U(j,8*(i-1)+1) = auxr1; U(j,8*(i-1)+2) = auxr2; U(j,8*(i-1)+3) = auxr3;  U(j,8*(i-1)+4) = auxr4
        U(j,8*(i-1)+5) = auxr5; U(j,8*(i-1)+6) = auxr6; U(j,8*(i-1)+7) = auxr7;  U(j,8*(i-1)+8) = auxr8
    end do  
    
  else if ( mod(dimnu,8) .eq. 7) then
  
    read(341,*) i1, i2, i3, i4, i5, i6, i7
    read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7
    omegaorig(8*(i-1)+1) = auxr1; omegaorig(8*(i-1)+2) = auxr2; omegaorig(8*(i-1)+3) = auxr3;  omegaorig(8*(i-1)+4) = auxr4
    omegaorig(8*(i-1)+5) = auxr5; omegaorig(8*(i-1)+6) = auxr6; omegaorig(8*(i-1)+7) = auxr7
    read(341,'(A)',IOSTAT=stat) inputline 
    do j=1,dimnu
        read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7
        U(j,8*(i-1)+1) = auxr1; U(j,8*(i-1)+2) = auxr2; U(j,8*(i-1)+3) = auxr3;  U(j,8*(i-1)+4) = auxr4
        U(j,8*(i-1)+5) = auxr5; U(j,8*(i-1)+6) = auxr6; U(j,8*(i-1)+7) = auxr7
    end do 

  else if ( mod(dimnu,8) .eq. 6) then
  
    read(341,*) i1, i2, i3, i4, i5, i6
    read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6
    omegaorig(8*(i-1)+1) = auxr1; omegaorig(8*(i-1)+2) = auxr2; omegaorig(8*(i-1)+3) = auxr3;  omegaorig(8*(i-1)+4) = auxr4
    omegaorig(8*(i-1)+5) = auxr5; omegaorig(8*(i-1)+6) = auxr6
    read(341,'(A)',IOSTAT=stat) inputline 
    do j=1,dimnu
        read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6
        U(j,8*(i-1)+1) = auxr1; U(j,8*(i-1)+2) = auxr2; U(j,8*(i-1)+3) = auxr3;  U(j,8*(i-1)+4) = auxr4
        U(j,8*(i-1)+5) = auxr5; U(j,8*(i-1)+6) = auxr6
    end do 
    
  else if ( mod(dimnu,8) .eq. 5) then
  
    read(341,*) i1, i2, i3, i4, i5 
    read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5 
    omegaorig(8*(i-1)+1) = auxr1; omegaorig(8*(i-1)+2) = auxr2; omegaorig(8*(i-1)+3) = auxr3;  omegaorig(8*(i-1)+4) = auxr4
    omegaorig(8*(i-1)+5) = auxr5 
    read(341,'(A)',IOSTAT=stat) inputline 
    do j=1,dimnu
        read(341,*) auxr1, auxr2, auxr3, auxr4, auxr5
        U(j,8*(i-1)+1) = auxr1; U(j,8*(i-1)+2) = auxr2; U(j,8*(i-1)+3) = auxr3;  U(j,8*(i-1)+4) = auxr4
        U(j,8*(i-1)+5) = auxr5 
    end do     

  else if ( mod(dimnu,8) .eq. 4) then
  
    read(341,*) i1, i2, i3, i4 
    read(341,*) auxr1, auxr2, auxr3, auxr4 
    omegaorig(8*(i-1)+1) = auxr1; omegaorig(8*(i-1)+2) = auxr2; omegaorig(8*(i-1)+3) = auxr3;  omegaorig(8*(i-1)+4) = auxr4
    read(341,'(A)',IOSTAT=stat) inputline 
    do j=1,dimnu
        read(341,*) auxr1, auxr2, auxr3, auxr4
        U(j,8*(i-1)+1) = auxr1; U(j,8*(i-1)+2) = auxr2; U(j,8*(i-1)+3) = auxr3;  U(j,8*(i-1)+4) = auxr4
    end do   

  else if ( mod(dimnu,8) .eq. 3) then
  
    read(341,*) i1, i2, i3 
    read(341,*) auxr1, auxr2, auxr3 
    omegaorig(8*(i-1)+1) = auxr1; omegaorig(8*(i-1)+2) = auxr2; omegaorig(8*(i-1)+3) = auxr3 
    read(341,'(A)',IOSTAT=stat) inputline 
    do j=1,dimnu
        read(341,*) auxr1, auxr2, auxr3 
        U(j,8*(i-1)+1) = auxr1; U(j,8*(i-1)+2) = auxr2; U(j,8*(i-1)+3) = auxr3 
    end do   

  else if ( mod(dimnu,8) .eq. 2) then
  
    read(341,*) i1, i2 
    read(341,*) auxr1, auxr2 
    omegaorig(8*(i-1)+1) = auxr1; omegaorig(8*(i-1)+2) = auxr2 
    read(341,'(A)',IOSTAT=stat) inputline 
    do j=1,dimnu
        read(341,*) auxr1, auxr2 
        U(j,8*(i-1)+1) = auxr1; U(j,8*(i-1)+2) = auxr2 
    end do         
  
  else if ( mod(dimnu,8) .eq. 1) then
  
    read(341,*) i1 
    read(341,*) auxr1 
    omegaorig(8*(i-1)+1) = auxr1 
    read(341,'(A)',IOSTAT=stat) inputline 
    do j=1,dimnu
        read(341,*) auxr1 
        U(j,8*(i-1)+1) = auxr1 
    end do       
    
  end if  
  
  close(341)
  
  
  
!  do i=1,number_of_atoms
!      sqmass(3*(i-1)+1)=sqrt(mass(i,2)); sqmass(3*(i-1)+2)=sqrt(mass(i,2)); sqmass(3*(i-1)+3)=sqrt(mass(i,2))
!  end do
  

  
  !! We re-express the U's in atomic units to check their orthonormality (the U's must be orthonormal in atomic units, but in Uvec.dat they are stored in Angstroms)
   P_au2A = 0.5291772d0
   do i1=1,dimnu
     do i2=1,dimnu
      !!!  En este caso concreto (Si35H36) quitamos el factor P_au2A pq tal factor estaba mal
      !!Uau(i1,i2) = U(i1,i2)!  U(i1,i2)/P_au2A       ! PENG:  vib_mode(ii,:) = vib_mode(ii,:)* P_au2A
      !!U(i1,i2) = Uau(i1,i2)    ! We rescale the U vector to atomic units, because we need everything in atomic units to calculate gDW, and U is an input parameter in the calculate_gDW subroutine.
      ! Quito el cambio comentado arriba
      !Uau(i1,i2) = U(i1,i2)/P_au2A 
     end do 
   end do  
    
  
! 
! 
!   ! Beware of this: Quantum Espresso provides U vectors which are orthogonal but NOT orthonormal. Hence one must renormalize them
!   ! so that they are orthoNORMAL, as in eq. (38) of Ponce et al. PRB 90 214304 (2014).
! 
!    ! We normalize -with masses- the U vectors (the fact that the output of QE gives normalized U does not mean that the corresponding xi are normalized)
!    do j=1,dimnu   
!      auxr1=0.d0
!      do i=1,dimnu
!        auxr1 = auxr1 + ( (U(i,j))**2 ) *sqmass(i)*sqmass(i) 
!      end do
!      auxr1 = sqrt(auxr1)
!      if (abs(auxr1) .lt. 0.000001) then
!         write(*,*) 'ERROR: Singularity of U vector ',j, ', module = ',auxr1
!         call exit(1)
!      end if
!      do i=1,dimnu
!        U(i,j) = U(i,j)/auxr1 
!      end do
!    end do 


  
  ! ORTHONORMALITY CHECK
! 
!   auxr4=0.d0  ! summation of the orthogonal terms
!   auxr5=0.d0  ! summation of terms equaling 1
!   auxr6=0.d0   ! maximum error
!   
!   do j=1,dimnu
!     do k=1,dimnu
!        
!        auxr1=0.d0
!        do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
!           auxr1=auxr1 + Uau(i,j)*Uau(i,k)*sqmass(i)*sqmass(i)
!        end do
!        if (j.ne.k) auxr4=auxr4+abs(auxr1)
!        if (j.eq.k) then
!          ! write(*,*) j, auxr1
!          auxr5=auxr5+abs(auxr1-1.d0)
!        end if  
!        if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
!     
!     end do
!   end do  
!   
!   auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+6)); auxr5=auxr5/(dble(dimnu-n_mode_beg+1))
!   if ( ( abs(auxr4) .lt. 0.00001d0) .and. (abs(auxr5).lt.0.001) ) then 
!      write(*,*);  write(*,*) '    - Orthonormality of U (in atomic units, with mass scaling) checked. '
!   else
!      write(*,*) ' ERROR in U: The first condition of orthonormality is not properly satisfied' 
!      write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
!      write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5   
!      call exit(1)  
!   end if
! 
!   
!   
!   ! We have read the U's; now we get xi:
!    do j=1,dimnu
!      do i=1,dimnu 
!        xi(i,j) = Uau(i,j) * sqmass(i)
!      end do
!    end do  
! 
! 
!    ! We normalize the xi vectors (the fact that the output of QE gives normalized U does not mean that the corresponding xi are normalized)
!    do j=1,dimnu 
!    
!      auxr1=0.d0
!      do i=1,dimnu
!        auxr1 = auxr1 + (xi(i,j))**2 
!      end do
!      auxr1 = sqrt(auxr1)
!      if (abs(auxr1) .lt. 0.000001) then
!         write(*,*) 'ERROR: Singularity of xi vector ',j, ', module = ',auxr1
!         call exit(1)
!      end if
!      do i=1,dimnu
!        xi(i,j) = xi(i,j)/auxr1 
!      end do
!      
!    end do 
!   
!  
! !  write(*,*) ' CHECKING THE ORTHONORMALITY OF EIGENVECTORS OF THE DYNAMICAL MATRIX'
! 
!   auxr4=0.d0  ! summation of the orthogonal terms
!   auxr5=0.d0  ! summation of terms equaling 1
!   auxr6=0.d0   ! maximum error
!   
!   do j=1,dimnu
!     do k=1,dimnu
!        
!        auxr1=0.d0
!        do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
!           auxr1=auxr1 + xi(i,j)*xi(i,k)
!        end do
!        if (j.ne.k) auxr4=auxr4+abs(auxr1)
!        if (j.eq.k) auxr5=auxr5+abs(auxr1-1.d0)
!        if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
!     
!     
!     end do
!   end do  
!   
!  
!   auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+6)); auxr5=auxr5/(dble(dimnu-n_mode_beg+1))
!   if ( ( abs(auxr4) .lt. 0.00001d0) .and. (abs(auxr5).lt.0.00001) ) then 
!      write(*,*) '    - Orthonormality of xi in rows checked. '
!   else
!      write(*,*) '   ERROR: The first condition of orthonormality is not properly satisfied' 
!      write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
!      write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5   
!      call exit(1)  
!   end if
!      
!    
!     
!      
!   auxr4=0.d0  ! summation of the orthogonal terms
!   auxr5=0.d0  ! summation of terms equaling 1
!   auxr6=0.d0   ! maximum error
!   
!   do j=7,dimnu
!     do k=7,dimnu
!        
!        auxr1=0.d0
!        do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
!           auxr1=auxr1 + xi(j,i)*xi(k,i)
!        end do
!        if (j.ne.k) auxr4=auxr4+abs(auxr1)
!        if (j.eq.k) auxr5=auxr5+abs(auxr1)
!        if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
!     
!     end do
!   end do  
!   auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+6)); auxr5=auxr5/(dble(dimnu-n_mode_beg+1))
!  
!   if ( ( abs(auxr4) .lt. 0.00001d0) .and. ( abs(auxr5-1.d0).lt.0.00001) ) then 
!      write(*,*) '    - Orthonormality of xi in columns checked. '; write(*,*)
!   else
!      write(*,*) '   ERROR: The second condition of orthonormality is not properly satisfied'  
!      write(*,*) ' The average orthonormaliz 2 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
!      write(*,*) ' The average orthonormaliz 2 is (parallel terms):', auxr5   
!      call exit(1)  
!   end if

  


end subroutine read_dyneq_output


!-------------------------------------------------------------------------------------------

subroutine omega_epsilon_to_atomic_units (eigvalsfromxml,units_omegaorig,units_epsilon0orig,dimnu,dimband, &        
                                             & omegaorig,epsilon0orig,omega,epsilon0, &
                                             & epsilon_displ_plus,epsilon_displ_minus)
     
  ! We use atomic units as done in Ponce et al. PRB 90, 214304 (2014); this is because our derivation of gDW
  ! depends on the definitions given in that paper.   
  ! This subroutine writes the phonon freqs. in the units of the electronic eigenvectors, multiplying the former
  ! values of the omegas by the 'factor'.   
       
  ! System parameters
  logical, Intent(In) :: eigvalsfromxml                 ! True if the eigenvalues are read from the .xml output of QE (which gives eigenvalues in atomic units already)
  character(40), Intent(In) :: units_omegaorig, units_epsilon0orig
  integer, Intent(In) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  integer, Intent(In) :: dimband                        ! Number of considered electronic eigenvalues
  real(8), Intent(In) :: omegaorig(dimnu)               ! Original omegas, in the original units
  real(8), Intent(In) :: epsilon0orig(dimband)			! Electronic eigenvalues in their original units
  real(8), Intent(Out):: omega(dimnu)                   ! Omegas in atomic units
  real(8), Intent(Out):: epsilon0(dimband)				! Electronic eigenvalues in atomic units
  real(8), Intent(InOut)::  epsilon_displ_plus(dimband,dimnu)  ! The eigenvalues of the system with displaced positions (with positive displacement)
  real(8), Intent(InOut)::  epsilon_displ_minus(dimband,dimnu) ! The eigenvalues of the system with displaced positions (with negative displacement)
   
  ! Local variable
  integer :: i
  real(8) :: factor  

! 
!  if (( trim(units_omegaorig) .ne. 'cm**-1' ) .and. ( trim(units_omegaorig) .ne. 'CM**-1' ).and. &
!   & ( trim(units_omegaorig) .ne. 'Cm**-1' ) &
!    &  .and. ( trim(units_omegaorig) .ne. 'mev' ) .and. ( trim(units_omegaorig) .ne. 'meV' ) &
!     .and. ( trim(units_omegaorig) .ne. 'a.u.' ).and. ( trim(units_omegaorig) .ne. 'A.U.' ) &
!     .and. ( trim(units_omegaorig) .ne. 'eV' ).and. ( trim(units_omegaorig) .ne. 'ev' ) &
!     &.and. ( trim(units_omegaorig) .ne. 'Hartree' ).and. ( trim(units_omegaorig) .ne. 'hartree' ) ) then

  write(*,*) '  *** Now converting the units of the frequency (',trim(units_omegaorig),' to atomic units). ***'
!  write(*,*) '         The units for phonon frequencies are ',units_omegaorig
!                write(*,*) '         The units for electronic eigenvalues are ',units_epsilon0orig; write(*,*)

  if    ( (trim(units_omegaorig).eq.'cm**-1') .or. (trim(units_omegaorig).eq.'cm-1') ) then
    factor=1.d0/219474.631370515d0 ! 0.00000455633d0
  else if (trim(units_omegaorig).eq.'meV') then  
    factor=1.d0/27211.385056d0  ! 0.0000367502d0   
  else if (trim(units_omegaorig).eq.'eV') then  
    factor=1.d0/27.211385056d0 ! 0.0367502d0
  else if (trim(units_omegaorig).eq.'a.u.') then  
    factor=1.d0           
  else
    write(*,*) '  ERROR: Conversion of units not established; please, use "eV", "meV", "cm**-1", or "a.u.", '
    write(*,*) '  or rewrite the <<omega_epsilon_to_atomic_units>> subroutine'
    call exit(1) 
  end if

   
   do i=1,dimnu
     omega(i) = (omegaorig(i))*factor
   end do
   
  do i=1,6
    omega(i) = 1.d0 !! We impose the nonphysical frequencies to be 1 (in a.u.) to avoid singularities when calculating the h factors
  end do
    

  
! We write the phonon frequencies with a format that Maat can read:
! wQ0.dat: 1st row: number of q-points and the number of phonon frequencies; then: q-index, branch-index, phonon frequency
  open(533,file='wQ0.dat',status='unknown')
  write(533,*) '  1  ', dimnu 
  do i=1,dimnu
    write(533,*) '  1  ', i,' ', omega(i)
  end do
  close(533)
  
 
  factor = 0.d0
 
 if ( .not. eigvalsfromxml ) then
   write(*,*) '  *** Now converting the units of the the electronic eigenvalues (',trim(units_epsilon0orig)
   write(*,*) '                                  to atomic units). ***'

  if      (trim(units_epsilon0orig).eq.'cm**-1') then
    factor=1.d0/219474.631370515d0 ! 0.00000455633d0
  else if (trim(units_epsilon0orig).eq.'meV') then  
    factor=1.d0/27211.385056d0  ! 0.0000367502d0 
  else if (trim(units_epsilon0orig).eq.'eV') then  
    factor=1.d0/27.211385056d0 ! 0.0367502d0
  else if (trim(units_epsilon0orig).eq.'a.u.') then  
    factor=1.d0           
  else
    write(*,*) '  ERROR: Conversion of units not established; please, use "eV", "meV", "cm**-1", or "a.u.", '
    write(*,*) '  or rewrite the <<omega_epsilon_to_atomic_units>> subroutine'
    call exit(1)    
  end if


   do i=1,dimband
     epsilon0(i) = (epsilon0orig(i))*factor
   end do
   
   do i=1,dimband
     do j=1,dimnu
        epsilon_displ_plus(i,j) = epsilon_displ_plus(i,j) * factor
     end do
   end do     
   
   do i=1,dimband
     do j=1,dimnu
        epsilon_displ_minus(i,j) = epsilon_displ_minus(i,j) * factor
     end do
   end do      
  
 else 
 
     epsilon0(:) = epsilon0orig(:)
  
 end if ! if ( .not. eigvalsfromxml ) then
  
! We write the electronic eigenvalues with a format that Maat can read:
! epsilon0orig0.dat: 1st row: number of bands, number of k-points, number of spins; then all three indices and the corresponding unrenormalized electronic eigenvalue
  open(633,file='epsilon_undisplaced.dat',status='unknown')
  write(633,*) dimband, '  1    1  ' 
  do i=1,dimband
    write(633,*) i, '  1  ', '  1  ', epsilon0(i)
  end do
  close(633)
        
       
end subroutine omega_epsilon_to_atomic_units   

!-------------------------------------------------------------------------------------------

 subroutine find_HOMO_LUMO_indices ( Ne, band_occupation,dimband,epsilon0, HOMO_deg, LUMO_deg, &
 								  &  InitialStateFF, FinalStateFF, &
                                  &  HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f)

 Integer, Intent(In) :: Ne                                                      ! Total number of electrons
 Integer, Intent(In) :: band_occupation											! Number of electrons per state (usually 2) 
 integer, Intent(In) :: dimband                                                 ! Number of considered electronic eigenvalues
 real(8), Intent(In) :: epsilon0(dimband)			                        	! Electronic eigenvalues in atomic units
 integer, Intent(In) :: HOMO_deg, LUMO_deg                                      ! Degeneracies of HOMO and LUMO, imposed by hand by the human user in frozenphonon.in (useful for states nearly-degenerate)
 integer, Intent(InOut) :: InitialStateFF, FinalStateFF     ! Initial and final state indices for the calculation of frozen-phonon where the possible crossover is sought.       
 Integer, Intent(Out) :: HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f ! Initial and final indices of the states of the HOMO and LUMO (which can be degenerate)

 ! Local variables
  integer :: i 
  character(100):: auxstr1, auxstr2
  


 if ( band_occupation .eq. 1 ) then
    HOMO_index_i = Ne; HOMO_index_f = Ne
    LUMO_index_i = Ne+1; LUMO_index_f = Ne+1
 else   
    HOMO_index_i = Ne/2+mod(Ne,2); HOMO_index_f = Ne/2+mod(Ne,2)
    LUMO_index_i = Ne/2+mod(Ne,2)+1; LUMO_index_f = Ne/2+mod(Ne,2)+1 
 end if   
 
 if (HOMO_deg .gt. 0) then
   HOMO_index_i = max(HOMO_index_f - HOMO_deg + 1,1)
 else
   do i=HOMO_index_f-1,1,-1
	 if ( (abs( epsilon0(HOMO_index_f)- epsilon0(i) ) ) .lt. (1.d0/27211.86d0) ) HOMO_index_i = i
   !  write(*,*) i,epsilon0(i),epsilon0(HOMO_index_f)- epsilon0(i)
   end do
 end if  

 if (LUMO_deg .gt. 0) then
   LUMO_index_f = min(LUMO_index_i + LUMO_deg - 1,dimband)
 else
   do i=LUMO_index_i+1,min(LUMO_index_i+10,dimband)
	 if ( (abs( epsilon0(LUMO_index_i)- epsilon0(i) ) ) .lt. (1.d0/27211.86d0) ) LUMO_index_f = i
   end do 
 end if  
 

 if (InitialStateFF .eq. 0) then
   InitialStateFF=1
 end if
 if (FinalStateFF .gt. dimband) then
   FinalStateFF=dimband
 end if
 
 if ( HOMO_index_i .lt. InitialStateFF ) then
   write(*,*) ' WARNING: Your "InitialStateFPh" variable is not low enough to read all HOMO states. '
 end if
 if ( LUMO_index_f .gt. FinalStateFF ) then
   write(*,*) ' WARNING: Your "FinalStateFPh" variable is not high enough to read all LUMO states. '
 end if 
 
 write(*,*)
 write(auxstr1,*) HOMO_index_i ;    call StripSpaces(auxstr1)
 write(auxstr2,*) HOMO_index_f ;    call StripSpaces(auxstr2) 
 if (HOMO_index_i .eq. HOMO_index_f) then 
   write(*,*) '       The HOMO is the ',trim(auxstr1),'-th state.'
 else   
   write(*,*) '       The HOMO is degenerate (states ',trim(auxstr1),' to ',trim(auxstr2),').'
   if (HOMO_deg .gt. 0) write(*,*) '        (the degeneracy of the HOMO was artificially imposed).'
 end if

 write(auxstr1,*) LUMO_index_i ;    call StripSpaces(auxstr1)
 write(auxstr2,*) LUMO_index_f ;    call StripSpaces(auxstr2) 
 if (LUMO_index_i .eq. LUMO_index_f) then 
   write(*,*) '       The LUMO is the ',trim(auxstr1),'-th state.'
 else   
   write(*,*) '       The LUMO is degenerate (states ',trim(auxstr1),' to ',trim(auxstr2),').'
   if (HOMO_deg .gt. 0) write(*,*) '        (the degeneracy of the LUMO was artificially imposed).'
 end if
 write(*,*)
  
  
  if ( FinalStateFF .lt. LUMO_index_i ) then               ! calculation of just the HOMO
 	write(*,*) '       In this run, just the frozen-phonon renormalization of the HOMO'
 	write(*,*) '       will be calculated.'; write(*,*) 
  else if ( InitialStateFF .gt. HOMO_index_f ) then  	    ! calculation of just the LUMO
 	write(*,*) '       In this run, just the frozen-phonon renormalization of the LUMO'
 	write(*,*) '       will be calculated.'; write(*,*) 	
  end if
 
 
 end subroutine find_HOMO_LUMO_indices

!-------------------------------------------------------------------------------------------


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

 subroutine remove_crossover_in_eigenvalues ( crossoversfromfile, sizeepsilon, dimnu, number_of_atoms, &
        & number_of_atomic_species, nfft, InitialStateFF, FinalStateFF, InitialMode, FinalMode, &
        & HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f, &
        & wfc_directory, Input_wfc_prefix, epsilon0, &
        & epsilon_displ_plus0, epsilon_displ_minus0, &
        & epsilon_displ_plus, epsilon_displ_minus)
        
  ! This subroutine removes the possible crossovers in the eigenvalues due to the finite-difference displacements.
  ! These crossovers result into completely wrong frozen-phonon calculations.
  ! This subroutine is a bit complicated to avoid the repetition of displaced states in case of degeneracy: if there is 
  ! degeneracy of eigenvalues, it is frequent that after doing the displacement, a given displaced eigenvalue has maximum 
  ! coupling with SEVERAL of the undisplaced degenerate HOMO (or LUMO) states. This means that if we use just the criterion
  ! of maximum coupling to establish the correspondence between displaced and undisplaced states we would be taking the same
  ! eigenvalue several times, what is unphysical and would lead to errors in the final result.
   
  
  Logical, Intent(In) :: crossoversfromfile ! True if the crossovers must be read from a given file
  Integer, Intent(In) :: sizeepsilon      ! Number of electronic eigenvalues
  Integer, Intent(In) :: dimnu            ! Total number of phonon frequencies
  Integer, Intent(In) :: number_of_atoms
  Integer, Intent(In) :: number_of_atomic_species
  Integer, Intent(In) :: nfft             ! Dimension in every direction of the wavefunction files from Quantum Espresso
  Integer, Intent(In) :: InitialStateFF, FinalStateFF     ! Initial and final state indices for the calculation of frozen-phonon where the possible crossover is sought.
  integer, Intent(In) :: InitialMode, FinalMode         ! Initial and final phonon modes to be considered in the calculations (they are useful to parallelize the calculations)
  Integer, Intent(In) :: HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f ! Initial and final indices of the states of the HOMO and LUMO (which can be degenerate)
  character(100):: wfc_directory, Input_wfc_prefix 
  real(8), Intent(In) :: epsilon0(sizeepsilon) ! Unrenormalized electronic eigenvalues (in a.u.)
  real(8), Intent(In) :: epsilon_displ_plus0(sizeepsilon,dimnu)  ! electronic eigenvalues after positive displacements (in a.u.), still with crossovers in eigenvalues
  real(8), Intent(In) :: epsilon_displ_minus0(sizeepsilon,dimnu) ! electronic eigenvalues after negative displacements (in a.u.), still with crossovers in eigenvalues
  real(8), Intent(Out) :: epsilon_displ_plus(sizeepsilon,dimnu)   ! electronic eigenvalues after positive displacements (in a.u.), without crossovers in eigenvalues
  real(8), Intent(Out) :: epsilon_displ_minus(sizeepsilon,dimnu)  ! electronic eigenvalues after negative displacements (in a.u.), without crossovers in eigenvalues
       
  ! Local variables
  integer :: nu, i, j, k, l, istateunperturb, istateplus, istateminus, entsprechend, iostat, stat
  real(8) :: auxr1, auxr2, lattparam, prod, prodentsprechend 
  real(8), allocatable :: epsilon_renorm(:)
  real(8),allocatable :: wfc_i(:,:,:), wfc_f(:,:,:)
  integer, allocatable :: taken(:,:), corresponding(:,:), correspwithinHOMOrange(:), correspwithinLUMOrange(:)
  character(100) :: auxstr1,  auxstr2,  auxstr3,  auxstr4,  auxstr5, auxstr6
  character(30) :: auxstr7, auxstr8, auxstr9, auxstr10
  logical :: repeated, exist, exist_flag1, exist_flag2
  
  
 !------------------------------------------------------------------------------------------- 
 !                 BLOCK OF CROSSOVER INFORMATION READ FROM FILE
 !------------------------------------------------------------------------------------------- 
  
 
  if (crossoversfromfile) then
  
	inquire(file = 'crossovers-HOMO.out', exist=exist_flag1)
	inquire(file = 'crossovers-LUMO.out', exist=exist_flag2)
	
    allocate(corresponding(-dimnu:dimnu,LUMO_index_f-HOMO_index_i+1)); corresponding(:,:) = 0 
	
	if ( ( .not. exist_flag1) .or. (.not. exist_flag2 ) ) then
	  write(*,*) '  ***************************************************************************'
	  write(*,*) "    ERROR: Please, make sure that the 'crossovers-HOMO.out' and "
	  write(*,*) "    'crossovers-LUMO.out' exist (or switch to 0 the 'crossoversfromfile' "
	  write(*,*) "     variable in the 'frozenphonon.in' file). "
	  write(*,*) '  ***************************************************************************'
	  write(*,*)
	  call exit(1)
	end if 

    ! We initialize the 'corresponding' variable to its unchanged values
    do istateunperturb=max(HOMO_index_i,InitialStateFF),min(LUMO_index_f,FinalStateFF)
	   do nu=InitialMode,FinalMode
		 corresponding(nu,istateunperturb-HOMO_index_i+1)=istateunperturb
		 corresponding(-nu,istateunperturb-HOMO_index_i+1)=istateunperturb
	   end do 
	 end do  

	open(445,file='crossovers-HOMO.out',status='unknown') 
	read(445,'(A)',IOSTAT=stat) auxstr1
	do while ((stat .eq. 0))
	  read(445,*,IOSTAT=stat) nu, j, k
	  corresponding(nu,j-HOMO_index_i+1) = k
	end do
     
	open(446,file='crossovers-LUMO.out',status='unknown') 
	read(446,'(A)',IOSTAT=stat) auxstr1
	do while ((stat .eq. 0))
	  read(446,*,IOSTAT=stat) nu, j, k
	  corresponding(nu,j-HOMO_index_i+1) = k
	end do
	
 
	do istateunperturb=max(HOMO_index_i,InitialStateFF),min(LUMO_index_f,FinalStateFF)
	   do nu=InitialMode,FinalMode
		 j = corresponding(nu,istateunperturb-HOMO_index_i+1)
		! write(*,*) 'state',istateunperturb,'j=',j
		! write(*,*)  epsilon_displ_plus0(j,nu) 
		 epsilon_displ_plus(istateunperturb,nu) = epsilon_displ_plus0(j,nu) 
		 k = corresponding(-nu,istateunperturb-HOMO_index_i+1)
		! write(*,*) 'state',istateunperturb,'k=',k
		! write(*,*)  epsilon_displ_minus0(k,nu) 
		 epsilon_displ_minus(istateunperturb,nu) = epsilon_displ_minus0(k,nu) 
	   end do 
	 end do  
   
    close(445); close(446)
    
    deallocate(corresponding)
    
	return
  
  end if
  
  
  
 !------------------------------------------------------------------------------------------- 
 !               BLOCK OF CROSSOVER INFORMATION CALCULATED FROM WAVEFUNCTIONS
 !------------------------------------------------------------------------------------------- 
  
 
  allocate(taken(sizeepsilon,dimnu),corresponding(-dimnu:dimnu,LUMO_index_f-HOMO_index_i+1))
  allocate(correspwithinHOMOrange(HOMO_index_f-HOMO_index_i+1),correspwithinLUMOrange(LUMO_index_f-LUMO_index_i+1))
  taken(:,:)=0; corresponding(:,:)=0
 
  do i=1,sizeepsilon
    epsilon_displ_plus(i,:) = epsilon0(i)
    epsilon_displ_minus(i,:)= epsilon0(i)
  end do  
 
  write(auxstr1,*) InitialStateFF ;  call StripSpaces(auxstr1)
  write(auxstr2,*) FinalStateFF ;    call StripSpaces(auxstr2)
  write(*,*) '  *** Now searching for crossovers in the eigenvalues (due to the  '
  write(*,*) '      finite-difference displacements) between the ' 
  write(*,*) '      ',trim(auxstr1), '-th and the ', trim(auxstr2),'-th states. *** '
  write(*,*) 

  write(auxstr1,*) trim(wfc_directory)//'Wfcs-undisplaced/'
  write(auxstr2,*) trim(wfc_directory)//'Wfcs-displaced+/'
  write(auxstr3,*) trim(wfc_directory)//'Wfcs-displaced-/'
    
  allocate(wfc_i(nfft,nfft,nfft),wfc_f(nfft,nfft,nfft))
  
  if ( FinalStateFF .lt. LUMO_index_i ) then               ! calculation of just the HOMO
 	open(444,file='wfc_product_results-HOMO.out',status='unknown')
	open(445,file='crossovers-HOMO.out',status='unknown') 
  else if ( InitialStateFF .gt. HOMO_index_f ) then  	    ! calculation of just the LUMO
	open(444,file='wfc_product_results-LUMO.out',status='unknown')
	open(446,file='crossovers-LUMO.out',status='unknown') 
  else  	                                                 ! calculation of both HOMO and LUMO
	open(444,file='wfc_product_results.out',status='unknown')
	open(445,file='crossovers-HOMO.out',status='unknown') 	
	open(446,file='crossovers-LUMO.out',status='unknown') 
  end if	


! PART OF POSITIVE DISPLACEMENTS
  
  write(444,*) '#        nu      undispl. wf.    displ. wf. +      < | > '
  write(445,*) '# nu | undispl. wf. | corresponding displ. wf.   '
  write(446,*) '# nu | undispl. wf. | corresponding displ. wf.   '
  
  ! Part A: Finding the states which maximum < | > with given state:
  
  do istateunperturb = max(HOMO_index_i,InitialStateFF), min(LUMO_index_f,FinalStateFF)
  
    write(auxstr4,*) 'wfc',istateunperturb,".wfn"; call StripSpaces(auxstr4)
    
    call read_wfc( nfft, Number_of_atoms, Number_of_atomic_species, auxstr1, auxstr4, wfc_i)
    do nu=InitialMode,FinalMode
    
      entsprechend=InitialStateFF
      prodentsprechend=0.d0
      
      do istateplus=InitialStateFF,FinalStateFF !
      
        write(auxstr5,*) 'mode',nu,'-wfc',istateplus,".wfn"; call StripSpaces(auxstr5)
        
        ! We look for the displaced wavefunction which corresponds to the given unperturbed wavefunction
        ! (we call this variable 'entsprechend'). It will be the wavefunction with maximum coupling. 
        call  read_wfc( nfft, Number_of_atoms, Number_of_atomic_species, auxstr2, auxstr5, wfc_f)
        call wfc_product (nfft, wfc_i, wfc_f, prod)
        write(444,*) nu, ' ',istateunperturb,' ', istateplus, '     ', prod
        if ( abs(prod) .gt. prodentsprechend) then
           entsprechend=istateplus
           prodentsprechend=abs(prod)
        end if   
        
      end do ! istateplus
      
      corresponding(nu,istateunperturb-HOMO_index_i+1)=entsprechend

    end do ! nu
    
  end do ! istateunperturb


  
 ! Part B: We check that, in case of degeneracy of states, it does not happen that a given
 !         displaced state is assigned to more than one degenerate undisplaced state.
  

  
  ! Part B.1: HOMO
  if (  (InitialStateFF .le. HOMO_index_i) ) then

	do nu=InitialMode,FinalMode
	 
	   correspwithinHOMOrange(:)=1
	 
	   ! If a given displaced state corresponds to two of the degenerate HOMO states, and the corresponding states lies away of the range of HOMO states, we specify it in the 'correspwithinHOMOrange' variable.    
       repeated = .false.
	   do i=1,HOMO_index_f-HOMO_index_i+1
		 do j=i+1,HOMO_index_f-HOMO_index_i+1
	       
		   if (( corresponding(nu,j) .eq. corresponding(nu,i) ) .and. (corresponding(nu,i) .ge. 0 ) ) then
		      repeated = .true.
		      !write(*,*) 'we have found repetition:',HOMO_index_i-1+i,'and ',HOMO_index_i-1+j
		      !write(*,*) ' unperturbed states correspond to', corresponding(nu,i) 
			  if (  ( corresponding(nu,i) .lt. HOMO_index_i  ) .or. (  corresponding(nu,i) .gt. HOMO_index_f  )  ) then
				correspwithinHOMOrange(i) = 0
			  end if
		   end if
	  

		 end do
	   end do
	   
	  ! If there were repeated states but all them lie within the HOMO range, then we merely consider that no crossover was produced; so we avoid to have repeated displaced eigenvalues used to find the final result.
	  if ( repeated ) then	
	  
		k = 0
		do  i=1,HOMO_index_f-HOMO_index_i+1
		  k = k + correspwithinHOMOrange(i)
		end do 
	  
		if ( k .eq. (HOMO_index_f-HOMO_index_i+1) ) then  ! This is, if all the (perhaps repeated) state indices of the displaced (degenerate) states are in the range of the undisplaced degenerate states.
		
		 ! write(*,*) 'The repetition does not matter, because all the HOMO corresponding indices of mode',nu   
		 ! write(*,*) 'lie within the HOMO range'  
		
		  do i=HOMO_index_i,HOMO_index_f
			 corresponding(nu,i-HOMO_index_i+1) = i
		  end do
		  
		else ! This is, if there are repeated indices of displaced eigenvalues, and these do not lie within the range of degenerate undisplaced states.

          !write(*,*) 'Not in the HOMO range'  

		  do i=1,HOMO_index_f-HOMO_index_i+1
			do j=i+1,HOMO_index_f-HOMO_index_i+1 
			  if ( (corresponding(nu,j) .eq. corresponding(nu,i)) .and. (corresponding(nu,i) .ge. 0 ) ) then

				write(auxstr7,*)nu;  call StripSpaces(auxstr7); write(auxstr8,*) corresponding(nu,i); call StripSpaces(auxstr8)
				write(auxstr9,*)  i+HOMO_index_i-1 ; call StripSpaces(auxstr9)
				write(auxstr10,*) j+HOMO_index_i-1;  call StripSpaces(auxstr10)
				write(*,*)'    WARNING: For the +',trim(auxstr7),'-th phonon mode, the ',trim(auxstr8),'-th eigenvalue is taken'
				write(*,*)'       several times (by the ',trim(auxstr9),'-th and by the ',trim(auxstr10),'-th eigenvalues).'
				write(*,*)

			  end if
			end do
		  end do

        end if ! k .eq. (HOMO_index_f-HOMO_index_i+1) 
		 
	  end if ! repeated
	  
	  ! We write the final corresponding indices
	  write(auxstr7,*)nu;  call StripSpaces(auxstr7);
	  do  i=1,HOMO_index_f-HOMO_index_i+1
		 write(auxstr8,*) corresponding(nu,i); call StripSpaces(auxstr8)
		 write(auxstr9,*)  i+HOMO_index_i-1 ; call StripSpaces(auxstr9)
		 if ( ( corresponding(nu,i) .ge.  HOMO_index_i)  .and. ( corresponding(nu,i) .le.  HOMO_index_f)   ) then
		   write(445,*) '   ',trim(auxstr7), '   ',trim(auxstr9),'   ', trim(auxstr8)
		 else
		   write(445,*) '   ',trim(auxstr7), '   ',trim(auxstr9),'   ', trim(auxstr8), '  *  '
		 end if  
	  end do
	   
	end do  ! nu

  end if ! (  (InitialStateFF .ge. HOMO_index_i) ) 
  
  
      
  ! Part B.2: LUMO
   if (  ( FinalStateFF .ge.  LUMO_index_f) ) then

	do nu=InitialMode,FinalMode
	 
	   correspwithinLUMOrange(:)=1
	 
	   ! If a given displaced state corresponds to two of the degenerate LUMO states, and the corresponding states lies away of the range of LUMO states, we specify it in the 'correspwithinLUMOrange' variable.    
       repeated = .false.
	   do i=LUMO_index_i-HOMO_index_i+1,LUMO_index_f-HOMO_index_i+1
		 do j=i+1,LUMO_index_f-HOMO_index_i+1
	  
		   if ( ( corresponding(nu,j) .eq. corresponding(nu,i) ) .and. (corresponding(nu,i) .ge. 0 ) ) then
		      repeated = .true.
		     ! write(*,*) 'we have found repetition:',HOMO_index_i-1+i,'and ',HOMO_index_i-1+j
		     ! write(*,*) ' unperturbed states correspond to', corresponding(nu,i) 
			  if (  ( corresponding(nu,i) .lt. LUMO_index_i  ) .or. (  corresponding(nu,i) .gt. LUMO_index_f  )  ) then
				correspwithinLUMOrange(i) = 0
			!	write(*,*) 'ATTENTION, the repeated index lies away from the LUMO range'
			  end if
		   end if
	  
		 end do
	   end do
	   
	  ! If there were repeated states but all them lie within the HOMO range, then we merely consider that no crossover was produced; so we avoid to have repeated displaced eigenvalues used to find the final result.
	  if ( repeated ) then	
	  
		k = 0
		do  i=1,LUMO_index_f-LUMO_index_i+1
		  k = k + correspwithinLUMOrange(i)
		end do 
	  
		if ( k .eq. (LUMO_index_f-LUMO_index_i+1) ) then  ! This is, if all the (perhaps repeated) state indices of the displaced (degenerate) states are in the range of the undisplaced degenerate states.
		
		 ! write(*,*) 'The repetition does not matter, because all the HOMO corresponding indices of mode',nu   
		 ! write(*,*) 'lie within the HOMO range'  
		
		  do i=LUMO_index_i,LUMO_index_f
			 corresponding(nu,i-HOMO_index_i+1) = i 
		  end do
		  
		else ! This is, if there are repeated indices of displaced eigenvalues, and these do not lie within the range of degenerate undisplaced states.

         ! write(*,*) 'Not in the LUMO range'  
         			  
		  do i=LUMO_index_i-HOMO_index_i+1,LUMO_index_f-HOMO_index_i+1
			do j=i+1,LUMO_index_f-HOMO_index_i+1
			 if ( (corresponding(nu,j) .eq. corresponding(nu,i)) .and. (corresponding(nu,i) .ge. 0 ) ) then
			 
				write(auxstr7,*)nu;  call StripSpaces(auxstr7); write(auxstr8,*) corresponding(nu,i); call StripSpaces(auxstr8)
				write(auxstr9,*)  i+HOMO_index_i-1 ; call StripSpaces(auxstr9)
				write(auxstr10,*) j+HOMO_index_i-1;  call StripSpaces(auxstr10)
				write(*,*)'    WARNING: For the +',trim(auxstr7),'-th phonon mode, the ',trim(auxstr8),'-th eigenvalue is taken'
				write(*,*)'       several times (by the ',trim(auxstr9),'-th and by the ',trim(auxstr10),'-th eigenvalues).'
				write(*,*)

			  end if
			end do
		  end do

        end if ! k .eq. (LUMO_index_f-LUMO_index_i+1)
		 
	  end if ! repeated



	  ! We write the final corresponding indices
	  write(auxstr7,*)nu;  call StripSpaces(auxstr7);
	  do  i=LUMO_index_i-HOMO_index_i+1,LUMO_index_f-HOMO_index_i+1
		 write(auxstr8,*) corresponding(nu,i); call StripSpaces(auxstr8)
		 write(auxstr9,*)  i+HOMO_index_i-1 ; call StripSpaces(auxstr9)	 
		 if ( ( corresponding(nu,i) .ge.  LUMO_index_i)  .and. ( corresponding(nu,i) .le.  LUMO_index_f)   ) then
		   write(446,*) '   ',trim(auxstr7), '   ',trim(auxstr9),'   ', trim(auxstr8)
		 else
		   write(446,*) '   ',trim(auxstr7), '   ',trim(auxstr9),'   ', trim(auxstr8), '  *  '
		 end if		 
	  end do
	   
	end do  ! nu

  end if !( ( FinalStateFF .ge.  LUMO_index_f) ) 
  

  
  ! Part C: We write the displaced eigenvalues without crossovers
 
   do istateunperturb=max(HOMO_index_i,InitialStateFF),min(LUMO_index_f,FinalStateFF)
     do nu=InitialMode,FinalMode
       j = corresponding(nu,istateunperturb-HOMO_index_i+1)
       epsilon_displ_plus(istateunperturb,nu) = epsilon_displ_plus0(j,nu) 
     end do 
   end do  
  




! PART OF NEGATIVE DISPLACEMENTS
  
  write(444,*) '#        nu      undispl. wf.    displ. wf. -      < | > '
  !write(445,*) '# nu | undispl. wf. | corresponding displ. - wf.   '
  !write(446,*) '# nu | undispl. wf. | corresponding displ. - wf.   '
  
  ! Part A: Finding the states which maximum < | > with given state:
  
  do istateunperturb = max(HOMO_index_i,InitialStateFF), min(LUMO_index_f,FinalStateFF)
  
    write(auxstr4,*) 'wfc',istateunperturb,".wfn"; call StripSpaces(auxstr4)
    
    call read_wfc( nfft, Number_of_atoms, Number_of_atomic_species, auxstr1, auxstr4, wfc_i)
    do nu=InitialMode,FinalMode
    
      entsprechend=InitialStateFF
      prodentsprechend=0.d0
      
      do istateminus=InitialStateFF,FinalStateFF !
      
        write(auxstr5,*) 'mode',nu,'-wfc',istateminus,".wfn"; call StripSpaces(auxstr5)
        
        ! We look for the displaced wavefunction which corresponds to the given unperturbed wavefunction
        ! (we call this variable 'entsprechend'). It will be the wavefunction with maximum coupling. 
        call  read_wfc( nfft, Number_of_atoms, Number_of_atomic_species, auxstr3, auxstr5, wfc_f)
        call wfc_product (nfft, wfc_i, wfc_f, prod)
        write(444,*) -nu, ' ',istateunperturb,' ', istateminus, '     ', prod
        if ( abs(prod) .gt. prodentsprechend) then
           entsprechend=istateminus
           prodentsprechend=abs(prod)
        end if   
        
      end do ! istateminus
      
      corresponding(-nu,istateunperturb-HOMO_index_i+1)=entsprechend

    end do ! nu
    
  end do ! istateunperturb


  
 ! Part B: We check that, in case of degeneracy of states, it does not happen that a given
 !         displaced state is assigned to more than one degenerate undisplaced state.
  
  ! Part B.1: HOMO
  if ( (InitialStateFF .le. HOMO_index_i) ) then

	do nu=InitialMode,FinalMode
	 
	   correspwithinHOMOrange(:)=1
	 
	   ! If a given displaced state corresponds to two of the degenerate HOMO states, and the corresponding states lies away of the range of HOMO states, we specify it in the 'correspwithinHOMOrange' variable.    
       repeated = .false.
	   do i=1,HOMO_index_f-HOMO_index_i+1
		 do j=i+1,HOMO_index_f-HOMO_index_i+1
	       
		   if (( corresponding(-nu,j) .eq. corresponding(-nu,i) ) .and. (corresponding(-nu,i) .ge. 0 ) ) then
		      repeated = .true.
		      !write(*,*) 'we have found repetition:',HOMO_index_i-1+i,'and ',HOMO_index_i-1+j
		      !write(*,*) ' unperturbed states correspond to', corresponding(nu,i) 
			  if (  ( corresponding(-nu,i) .lt. HOMO_index_i  ) .or. (  corresponding(-nu,i) .gt. HOMO_index_f  )  ) then
				correspwithinHOMOrange(i) = 0
			  end if
		   end if
	  

		 end do
	   end do
	   
	  ! If there were repeated states but all them lie within the HOMO range, then we merely consider that no crossover was produced; so we avoid to have repeated displaced eigenvalues used to find the final result.
	  if ( repeated ) then	
	  
		k = 0
		do  i=1,HOMO_index_f-HOMO_index_i+1
		  k = k + correspwithinHOMOrange(i)
		end do 
	  
		if ( k .eq. (HOMO_index_f-HOMO_index_i+1) ) then  ! This is, if all the (perhaps repeated) state indices of the displaced (degenerate) states are in the range of the undisplaced degenerate states.
		
		 ! write(*,*) 'The repetition does not matter, because all the HOMO corresponding indices of mode',nu   
		 ! write(*,*) 'lie within the HOMO range'  
		
		  do i=HOMO_index_i,HOMO_index_f
			 corresponding(-nu,i-HOMO_index_i+1) = i
		  end do
		  
		else ! This is, if there are repeated indices of displaced eigenvalues, and these do not lie within the range of degenerate undisplaced states.

          !write(*,*) 'Not in the HOMO range'  

		  do i=1,HOMO_index_f-HOMO_index_i+1
			do j=i+1,HOMO_index_f-HOMO_index_i+1 
			  if ( (corresponding(-nu,j) .eq. corresponding(-nu,i)) .and. (corresponding(-nu,i) .ge. 0 ) ) then

				write(auxstr7,*)nu;  call StripSpaces(auxstr7); write(auxstr8,*) corresponding(-nu,i); call StripSpaces(auxstr8)
				write(auxstr9,*)  i+HOMO_index_i-1 ; call StripSpaces(auxstr9)
				write(auxstr10,*) j+HOMO_index_i-1;  call StripSpaces(auxstr10)
				write(*,*)'    WARNING: For the -',trim(auxstr7),'-th phonon mode, the ',trim(auxstr8),'-th eigenvalue is taken'
				write(*,*)'       several times (by the ',trim(auxstr9),'-th and by the ',trim(auxstr10),'-th eigenvalues).'
				write(*,*)

			  end if
			end do
		  end do

        end if ! k .eq. (HOMO_index_f-HOMO_index_i+1) 
		 
	  end if ! repeated
	  
	  ! We write the final corresponding indices
	  write(auxstr7,*)-nu;  call StripSpaces(auxstr7);
	  do  i=1,HOMO_index_f-HOMO_index_i+1
		 write(auxstr8,*) corresponding(-nu,i); call StripSpaces(auxstr8)
		 write(auxstr9,*)  i+HOMO_index_i-1 ; call StripSpaces(auxstr9)
		 if ( ( corresponding(-nu,i) .ge.  HOMO_index_i)  .and. ( corresponding(-nu,i) .le.  HOMO_index_f)   ) then
		   write(445,*) '   ',trim(auxstr7), '   ',trim(auxstr9),'   ', trim(auxstr8)
		 else
		   write(445,*) '   ',trim(auxstr7), '   ',trim(auxstr9),'   ', trim(auxstr8), '  *  '
		 end if
	  end do
	   
	end do  ! nu

  end if ! ( HOMO_index_f .ne. HOMO_index_i ) 
  
  
  ! Part B.2: LUMO
   if (  ( FinalStateFF .ge.  LUMO_index_f) ) then

	do nu=InitialMode,FinalMode
	 
	   correspwithinLUMOrange(:)=1
	 
	   ! If a given displaced state corresponds to two of the degenerate LUMO states, and the corresponding states lies away of the range of LUMO states, we specify it in the 'correspwithinLUMOrange' variable.    
       repeated = .false.
	   do i=LUMO_index_i-HOMO_index_i+1,LUMO_index_f-HOMO_index_i+1
		 do j=i+1,LUMO_index_f-HOMO_index_i+1
	  
		   if ( ( corresponding(-nu,j) .eq. corresponding(-nu,i) ) .and. (corresponding(-nu,i) .ge. 0 ) ) then
		      repeated = .true.
		      !write(*,*) 'we have found repetition:',HOMO_index_i-1+i,'and ',HOMO_index_i-1+j
		      !write(*,*) ' unperturbed states correspond to', corresponding(nu,i) 
			  if (  ( corresponding(-nu,i) .lt. LUMO_index_i  ) .or. (  corresponding(-nu,i) .gt. LUMO_index_f  )  ) then
				correspwithinLUMOrange(i) = 0
				!write(*,*) 'ATTENTION, the repeated index lies away from the HOMO range'
			  end if
		   end if
	  
		 end do
	   end do
	   
	  ! If there were repeated states but all them lie within the LUMO range, then we merely consider that no crossover was produced; so we avoid to have repeated displaced eigenvalues used to find the final result.
	  if ( repeated ) then	
	  
		k = 0
		do  i=1,LUMO_index_f-LUMO_index_i+1
		  k = k + correspwithinLUMOrange(i)
		end do 
	  
		if ( k .eq. (LUMO_index_f-LUMO_index_i+1) ) then  ! This is, if all the (perhaps repeated) state indices of the displaced (degenerate) states are in the range of the undisplaced degenerate states.
		
		  !write(*,*) 'The repetition does not matter, because all the HOMO corresponding indices of mode',nu   
		  !write(*,*) 'lie within the HOMO range'  
		
		  do i=LUMO_index_i,LUMO_index_f
			 corresponding(-nu,i-HOMO_index_i+1) = i 
		  end do
		  
		else ! This is, if there are repeated indices of displaced eigenvalues, and these do not lie within the range of degenerate undisplaced states.

          !write(*,*) 'Not in the LUMO range'  
         			  
		  do i=LUMO_index_i-HOMO_index_i+1,LUMO_index_f-HOMO_index_i+1
			do j=i+1,LUMO_index_f-HOMO_index_i+1
			 if ( (corresponding(-nu,j) .eq. corresponding(-nu,i)) .and. (corresponding(-nu,i) .ge. 0 ) ) then
			 
				write(auxstr7,*)nu;  call StripSpaces(auxstr7); write(auxstr8,*) corresponding(-nu,i); call StripSpaces(auxstr8)
				write(auxstr9,*)  i+HOMO_index_i-1 ; call StripSpaces(auxstr9)
				write(auxstr10,*) j+HOMO_index_i-1;  call StripSpaces(auxstr10)
				write(*,*)'    WARNING: For the -',trim(auxstr7),'-th phonon mode, the ',trim(auxstr8),'-th eigenvalue is taken'
				write(*,*)'       several times (by the ',trim(auxstr9),'-th and by the ',trim(auxstr10),'-th eigenvalues).'
				write(*,*)

			  end if
			end do
		  end do

        end if ! k .eq. (LUMO_index_f-LUMO_index_i+1)
		 
	  end if ! repeated

	  ! We write the final corresponding indices
	  write(auxstr7,*)-nu;  call StripSpaces(auxstr7);
	  do  i=LUMO_index_i-HOMO_index_i+1,LUMO_index_f-HOMO_index_i+1
		 write(auxstr8,*) corresponding(-nu,i); call StripSpaces(auxstr8)
		 write(auxstr9,*)  i+HOMO_index_i-1 ; call StripSpaces(auxstr9)
		 if ( ( corresponding(-nu,i) .ge.  LUMO_index_i)  .and. ( corresponding(-nu,i) .le.  LUMO_index_f)   ) then
		   write(446,*) '   ',trim(auxstr7), '   ',trim(auxstr9),'   ', trim(auxstr8)
		 else
		   write(446,*) '   ',trim(auxstr7), '   ',trim(auxstr9),'   ', trim(auxstr8), '  *  '
		 end if
	  end do
	   
	end do  ! nu

  end if ! 
  

  
  ! Part C: We write the displaced eigenvalues without crossovers
 
   do istateunperturb=max(HOMO_index_i,InitialStateFF),min(LUMO_index_f,FinalStateFF)
     do nu=InitialMode,FinalMode
       j = corresponding(-nu,istateunperturb-HOMO_index_i+1)
       epsilon_displ_minus(istateunperturb,nu) = epsilon_displ_minus0(j,nu) 
     end do 
   end do  
  

  write(*,*); write(*,*) '  *** The crossovers due to the finite-difference displacement were avoided '
  write(*,*) '  by checking the similarity between wavefunctions. Results can be checked'
  if ( FinalStateFF .lt. LUMO_index_i ) then               ! calculation of just the HOMO
 	write(*,*) '   at the "wfc_product_results-HOMO.out" file. *** '; 
  else if ( InitialStateFF .gt. HOMO_index_f ) then  	    ! calculation of just the LUMO
	write(*,*) '   at the "wfc_product_results-LUMO.out" file. *** '; 
  else  	                                                 ! calculation of both HOMO and LUMO
	write(*,*) '   at the "wfc_product_results.out" file. *** '; 
  end if 
  write(*,*)
  
  close(444); close(445); close(446)
  deallocate(wfc_i,wfc_f,taken,correspwithinHOMOrange,correspwithinLUMOrange)
  

  
  end  subroutine remove_crossover_in_eigenvalues
  
  
!-------------------------------------------------------------------------------------------


 subroutine calculate_frozen_phonon_renormalization ( T, InitialT, FinalT, number_of_Ts,&
        & sizeepsilon, dimnu, InitialStateFF, FinalStateFF,  &
        & InitialMode, FinalMode, HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f, &
        &  displ_parameter, omega, epsilon0, &
        & epsilon_displ_plus, epsilon_displ_minus)
        
  ! This subroutine calculates the frozen-phonon renormalization to the electronic bands.
  ! Take into account that the size of the displacement must be appropriate for the concrete 
  ! band that you want to renormalize, it must be neither too high nor too small. For example,
  ! for C14H20 (cutoff 30 Ry), displacements of +/- U_{\nu} make the renormalization of the HOMO
  ! to be of the order of hundreds of meV, and the renormalization of the LUMO to be of the order of
  ! tens of meV. But with the same displacement, for example, the band 3 experiences a renormalization
  ! of -15 eV, while the band 49 experiences a renormalization of +381 keV; these are unphysical, 
  ! because the quadratic region is forsaken (in both examples one can appreciate that the renormalizations
  ! due to the + and - displacements are completely different, which does not happen in the HOMO and the LUMO).
  
  Real(8), Intent(In) :: T                ! Temperature of the system (in atomic units)
  real(8), Intent(In) :: InitialT, FinalT ! Optional variables; if present, they are the initial and final temperatures for a sweep in T
  integer, Intent(In) :: number_of_Ts    ! Optional variable; if present, it indicates the number of temperatures considered in the sweep in T
  Integer, Intent(In) :: sizeepsilon      ! Number of electronic eigenvalues
  Integer, Intent(In) :: dimnu            ! Total number of phonon frequencies
  Integer, Intent(In) :: InitialStateFF, FinalStateFF     ! Initial and final state indices for the calculation of frozen-phonon where the possible crossover is sought.
  integer, Intent(In) :: InitialMode, FinalMode         ! Initial and final phonon modes to be considered in the calculations (they are useful to parallelize the calculations)
  Integer, Intent(In) :: HOMO_index_i, HOMO_index_f, LUMO_index_i, LUMO_index_f ! Initial and final indices of the states of the HOMO and LUMO (which can be degenerate)
  real(8), Intent(In) :: displ_parameter  ! Number multiplied by U (the eigenvector of the dynamical matrix) to perform the displacements used to do the finite-difference.
  real(8), Intent(In) :: omega(dimnu)     ! Phonon frequencies (in a.u.)
  real(8), Intent(In) :: epsilon0(sizeepsilon) ! Unrenormalized electronic eigenvalues (in a.u.)
  real(8), Intent(In) :: epsilon_displ_plus(sizeepsilon,dimnu)  ! electronic eigenvalues after positive displacements (in a.u.)
  real(8), Intent(In) :: epsilon_displ_minus(sizeepsilon,dimnu) ! electronic eigenvalues after negative displacements (in a.u.)
       
  ! Local variables
  integer :: nu, i, numTs, iT
  real(8) :: auxr1, auxr2, auxr3, auxr4, Ttemporal, HOMO, LUMO
  real(8), allocatable :: epsilon_renorm(:)
  character(20) :: auxstr1
  logical :: exist, exist_flag
  

  
  allocate(epsilon_renorm(sizeepsilon))
  
 !  write(*,*) 'eigv 49 - undispl',epsilon0(49)*27.211385056d0; write(*,*)
!   do nu=7,dimnu
!     write(*,'(I6,A,G13.6,G12.5,G12.5,G12.5)') nu, '  ',epsilon_displ_plus(49,nu)*27.211, &
!      & (epsilon_displ_plus(49,nu)-epsilon0(49))*27211.385056d0, &
!      & (epsilon_displ_minus(49,nu)-epsilon0(49))*27211.385056d0, &
!      & ((epsilon_displ_plus(49,nu)+epsilon_displ_minus(49,nu))/2.d0-epsilon0(49))*27211.385056d0
!   end do
!   write(*,*)     


 if (number_of_Ts .le. 0) then
   numTs = 1
   Ttemporal = T
 else
   numTs=number_of_Ts 
   inquire(file = 'HOMO_renorm_vs_T.txt', exist=exist_flag)
   if( .not. exist_flag) then
	  OPEN(172, FILE='HOMO_renorm_vs_T.txt',status='unknown')
	  write(172,*) '# T (K)      HOMO renorm. (meV)     '
	  close(172)
   end if 
   inquire(file = 'LUMO_renorm_vs_T.txt', exist=exist_flag)
   if( .not. exist_flag) then
	  OPEN(173, FILE='LUMO_renorm_vs_T.txt',status='unknown')
	  write(173,*) '# T (K)      LUMO renorm. (meV)     '
	  close(173)
   end if   
   open(172, file='HOMO_renorm_vs_T.txt', status="old", position="append", action="write")
   open(173, file='LUMO_renorm_vs_T.txt', status="old", position="append", action="write")
 end if  


 do iT=1,numTs

 
      if (number_of_Ts  .gt. 1 ) then
        Ttemporal = InitialT + (FinalT - InitialT)*dble(iT-1)/(dble(numTs-1))
        write(*,*) '----------------------------------------------------------------------------'
      end if 
      write(*,*);
      write(*,'(A, G15.5, A)') '   Temperature = ', Ttemporal*27211.3860217d0/0.08621738d0, 'K';write(*,*)


   
		 
  
      ! We write the renormalization of HOMO
	  if (  ( InitialStateFF .le.  HOMO_index_i) ) then
	   
	   do i=HOMO_index_i,HOMO_index_f

         if ( HOMO_index_i .eq. HOMO_index_f ) then
           write(*,*) 'HOMO - undisplaced =',epsilon0(HOMO_index_f)*27.211385056d0,' eV'; write(*,*)
         else
           write(auxstr1,*) i; call StripSpaces(auxstr1)
           write(*,*) 'HOMO - undisplaced(',trim(auxstr1),')=',epsilon0(i)*27.211385056d0,' eV'; write(*,*)
         end if  

		 do nu=InitialMode,FinalMode
		   write(*,'(I6,A,G13.6,A,G12.5,A,G12.5,A,G12.5)') nu, '  ', &
			& ( epsilon_displ_plus(i,nu) - epsilon0(i))*27211.385056d0,' ', &
			& (epsilon_displ_minus(i,nu) - epsilon0(i))*27211.385056d0,' ', &
			& ((epsilon_displ_plus(i,nu) + epsilon_displ_minus(i,nu))/2.d0-epsilon0(i))&
																	  & *27211.385056d0 , ' ',&
			& ((epsilon_displ_plus(i,nu) + epsilon_displ_minus(i,nu))/2.d0-epsilon0(i))&
																	  & *27211.385056d0/(2.d0*(omega(nu))*displ_parameter**2)												  													  
		 end do
		 write(*,*)
		 
        end do ! i=HOMO_index_i,HOMO_index_f 
        
	  end if

	  ! We write the renormalization of LUMO
	  if (  ( FinalStateFF .ge.  LUMO_index_f) ) then

	   do i=LUMO_index_i,LUMO_index_f

         if ( LUMO_index_i .eq. LUMO_index_f ) then
           write(*,*) 'LUMO - undisplaced =',epsilon0(LUMO_index_f)*27.211385056d0,' eV'; write(*,*)
         else
           write(auxstr1,*) i; call StripSpaces(auxstr1)
           write(*,*) 'LUMO - undisplaced(',trim(auxstr1),')=',epsilon0(i)*27.211385056d0,' eV'; write(*,*)
         end if  
	
		 do nu=InitialMode,FinalMode
		   write(*,'(I6,A,G13.6,A,G12.5,A,G12.5,A,G12.5)') nu, '  ', &
			& (epsilon_displ_plus(i,nu)-epsilon0(i))*27211.385056d0 ,' ', &
			& (epsilon_displ_minus(i,nu)-epsilon0(i))*27211.385056d0 ,' ', &
			& ((epsilon_displ_plus(i,nu)+epsilon_displ_minus(i,nu))/2.d0-epsilon0(i))&
					 & *27211.385056d0 ,' ', &
			& ((epsilon_displ_plus(i,nu)+epsilon_displ_minus(i,nu))/2.d0-epsilon0(i))&
					 & *27211.385056d0/(2.d0*(omega(nu))*displ_parameter**2)	 
		 end do
 		 write(*,*)
 		
        end do !i=LUMO_index_i,LUMO_index_f
 
	  end if
  
 
	  do i=1,sizeepsilon
		auxr1=0.d0
		do nu=InitialMode,FinalMode
		   auxr2=(epsilon_displ_plus(i,nu)+epsilon_displ_minus(i,nu))/2.d0 - epsilon0(i)
		   auxr2=auxr2/((omega(nu)))
		   auxr1=auxr1+auxr2
		end do
		epsilon_renorm(i) = epsilon0(i) + (auxr1/(2.d0*(displ_parameter)**2))
	  end do


	  if ( FinalStateFF .lt. LUMO_index_i ) then                ! calculation of just the HOMO
		open(445,file='frozen-phonon-HOMO.out',status='unknown')
	  else if ( InitialStateFF .gt. HOMO_index_f ) then  	    ! calculation of just the LUMO
		open(445,file='frozen-phonon-LUMO.out',status='unknown')
	  else  	                                                ! calculation of both HOMO and LUMO
		open(445,file='frozen-phonon.out',status='unknown')
	  end if	  
  
  
	  write(445,*) '# Eigv.   Unrenorm.         Renorm.      Renormalization '
	  write(445,*) '# index   eigv. (eV)       eigv. (eV)        (meV)'
	  do i=1,sizeepsilon
		write(445,'(I7,A,G14.6,A,G14.6,A,G15.5)') i, '  ', epsilon0(i)*27.211385056d0, '   ', epsilon_renorm(i)*27.211385056d0, &
					 &  ' ',(-epsilon0(i)+epsilon_renorm(i))*27211.385056d0
	  end do
	  close(445)

      write(*,*); write(*,*) ' ---------------------------------------------------------------- ';write(*,*)
  
	  ! We calculate the renormalized HOMO
	  if ( InitialStateFF .le. HOMO_index_i ) then 
   
		 auxr1=0.d0
		 do i=HOMO_index_i,HOMO_index_f
		 
		   auxr3=0.d0
		   do nu=InitialMode,FinalMode
	
			  if ( (dble(omega(nu)/Ttemporal) .gt. 30.d0) .or. ( Ttemporal .lt. 0.0000001d0)  ) then ! To avoid overflows
				 auxr4=(0.d0,0.d0)
			  else  
				 auxr4 = ( 1.d0 / ( Exp(omega(nu)/Ttemporal) - 1.d0)  ) 
			  end if    
	
			  auxr2=(epsilon_displ_plus(i,nu)+epsilon_displ_minus(i,nu))/2.d0 - epsilon0(i)
			  auxr2=(auxr2/((omega(nu))))*(1.d0+2.d0*auxr4)
			  auxr3=auxr3+auxr2
			  !auxr1=auxr1+auxr2
	   
		   end do
		   
		   write(auxstr1,*) i; call StripSpaces(auxstr1)
		   if (HOMO_index_i.ne.HOMO_index_f) write(*,'(A,A,A,G15.6,A)') '   The contribution of state ',trim(auxstr1),' is ', &
	         & (auxr3*27211.385056d0/(2.d0*((displ_parameter)**2)))/(dble(1+HOMO_index_f-HOMO_index_i)),' meV'
		   auxr1=auxr1+auxr3
		   
		 end do
		 
		 HOMO=0.d0
		 do i=HOMO_index_i,HOMO_index_f
		   HOMO = HOMO + epsilon0(i)
		 end do 
		 HOMO = HOMO / dble(HOMO_index_f-HOMO_index_i+1)
	 
		 auxr1=((auxr1/(2.d0*((displ_parameter)**2)))/(dble(1+HOMO_index_f-HOMO_index_i)))
	 
		 write(*,*)  ;	write(*,*) ' ========================= HOMO ==================================='
		 write(*,'(A,G13.6,A,G13.6,A)') '    original: ', HOMO*27.211385056d0, &
					 & ' eV  - renormalized: ', (HOMO+auxr1)*27.211385056d0,' eV'
		 write(*,'(A,G16.4,A)') '    difference: ', (auxr1)*27211.385056d0,' meV   <====='
	     write(*,*) ' =================================================================='
		 
		 if (number_of_Ts .gt. 0) write(172,'(G16.7, G16.7)') Ttemporal*27211.3860217d0/0.08621738d0, (auxr1)*27211.385056d0
	 
	  end if
  
      write(*,*)
 
	  ! We calculate the renormalized LUMO
	  if ( FinalStateFF .ge. LUMO_index_f ) then  
		 auxr1=0.d0
		 do i=LUMO_index_i,LUMO_index_f
		   auxr3=0.d0
		   do nu=InitialMode,FinalMode
	
			  if ( (dble(omega(nu)/Ttemporal) .gt. 30.d0) .or. ( Ttemporal .lt. 0.0000001d0)  ) then ! To avoid overflows
				 auxr4=(0.d0,0.d0)
			  else  
				 auxr4 = ( 1.d0 / ( Exp(omega(nu)/Ttemporal) - 1.d0)  ) 
			  end if    
	
			  auxr2=(epsilon_displ_plus(i,nu)+epsilon_displ_minus(i,nu))/2.d0 - epsilon0(i)
			  auxr2=(auxr2/((omega(nu))))*(1.d0+2.d0*auxr4)
			  !auxr1=auxr1+auxr2
			  auxr3=auxr3+auxr2
	   
		   end do
		   
		   write(auxstr1,*) i; call StripSpaces(auxstr1)
		   if (LUMO_index_i.ne.LUMO_index_f)	write(*,'(A,A,A,G15.6,A)') '   The contribution of state ',trim(auxstr1),' is ', &
	         & (auxr3*27211.385056d0/(2.d0*((displ_parameter)**2)))/(dble(1+LUMO_index_f-LUMO_index_i)),' meV '
		   auxr1=auxr1+auxr3	   
		   
		 end do
		 

		 
		 LUMO=0.d0
		 do i=LUMO_index_i,LUMO_index_f
		   LUMO = LUMO + epsilon0(i)
		 end do 
		 LUMO = LUMO / dble(LUMO_index_f-LUMO_index_i+1)		 
		 
		 auxr1=((auxr1/(2.d0*(displ_parameter)**2))/(dble(1+LUMO_index_f-LUMO_index_i))) 
		 
		 write(*,*)  ; write(*,*) ' ========================= LUMO ==================================='
		 write(*,'(A,G13.6,A,G13.6,A)') '    original: ', LUMO*27.211385056d0, &
			 & ' eV  - renormalized: ', (LUMO+auxr1)*27.211385056d0,' eV'
		 write(*,'(A,G16.4,A)') '    difference: ', (auxr1)*27211.385056d0,' meV   <====='
	     write(*,*) ' =================================================================='
		 write(*,*) ; write(*,*)
		 if (number_of_Ts .gt. 0) write(173,'(G16.7, G16.7)')  Ttemporal*27211.3860217d0/0.08621738d0, (auxr1)*27211.385056d0
	  end if
  
  end do !  do iT=1,numTs	 
 
  close(172);close(173)
  
  deallocate(epsilon_renorm)  
    
end  subroutine calculate_frozen_phonon_renormalization


!-------------------------------------------------------------------------------------------


subroutine read_nfft(istate, nat, n_atom_type, wfc_directory, Input_wfc_prefix, nfft)

 ! This subroutine reads the <istate>-th wavefunction and it stores it in the <wfc_i> variable (normalized in reciprocal space).
 ! The wavefunction is stored in <wfc_directory>, and its name is "<Input_wfc_prefix//istate//.wfn>"

!use mytypes
!use elph_pwscf_mod, only: wfc_i, nfft, Input_wfc, nat,n_atom_type, lattparam, wfc_directory

implicit none

    integer, Intent(In) :: istate,  nat, n_atom_type
    character(len=100), Intent(In) :: wfc_directory, Input_wfc_prefix    ! Names of the directory where the potentials and wavefunctions are stored and prefix of the name of the wavefunction files
    integer, Intent(Out) :: nfft                                         ! Dimension in every direction of the wavefunction files from Quantum Espresso


   !! local variables
   integer :: status, ierror ! I/O status
   logical :: exist_flag
   character(50) :: auxstr1
   character(len=100) :: auxstr3, auxstr4

   write(auxstr4,*) istate,".wfn"; call StripSpaces(auxstr4)
   write(auxstr3,*) trim(Input_wfc_prefix)//trim(auxstr4)
   write(auxstr1,*) trim(wfc_directory)//"/Wfcs-undisplaced/"//trim(auxstr3)
   call StripSpaces(auxstr1)
   inquire(file = auxstr1, exist=exist_flag)
   if(exist_flag) then
      OPEN(UNIT=3, FILE=auxstr1, STATUS='OLD', ACTION='READ', IOSTAT=status)
   else
     write(*,*) '  ***************************************************************************'
     write(*,*) "    ERROR: The ",trim(auxstr1),' file does not exist. '
     write(*,*) '           Please, provide it.'
     write(*,*) '  ***************************************************************************'
     write(*,*)
     call exit(1)
   endif     
   
   
   IF (status /= 0) THEN
       call exit(1)
       write(*,*) 'Error while reading wavefunction file.'
   END IF
   
   read(3,*) nfft
   
 end subroutine read_nfft
   
!-------------------------------------------------------------------------------------------


subroutine read_wfc( nfft, nat, n_atom_type, directory_name, wfn_name, wfc_i)

 ! This subroutine reads the <istate>-th wavefunction and it stores it in the <wfc_i> variable (normalized in reciprocal space).
 ! The wavefunction is stored in <wfc_directory>, and its name is "<Input_wfc_prefix//istate//.wfn>"

!use mytypes
!use elph_pwscf_mod, only: wfc_i, nfft, Input_wfc, nat,n_atom_type, lattparam, wfc_directory

implicit none

    integer, Intent(In) :: nfft, nat, n_atom_type
    character(len=100), Intent(In) :: directory_name, wfn_name !wfc_directory, Input_wfc_prefix    !! names of the directories where the potentials and wavefunctions are stored
    real(8), Intent(Out) :: wfc_i(nfft,nfft,nfft) 

   !! local variables
   integer :: ii, jj, kk, ll
   integer :: temp_int, temp1, temp2
   integer :: status, ierror ! I/O status
   logical :: exist_flag
   character(50) :: auxstr1

   real(8):: temp, temp_re 
   real(8):: array_temp(4)

   character(len=40) :: temp_ch1
   character(len=100) :: auxstr2, auxstr3, auxstr4
   character(len=2)  :: temp_ch2

   real(8), allocatable :: wfc_temp(:)

   allocate(wfc_temp(nfft*nfft*nfft))
   
   wfc_i(:,:,:)=0.d0

   write(auxstr1,*) trim(directory_name)//"/"//trim(wfn_name); call StripSpaces(auxstr1)
    
   inquire(file = auxstr1, exist=exist_flag)
   if(exist_flag) then
      OPEN(UNIT=3, FILE=auxstr1, STATUS='OLD', ACTION='READ', IOSTAT=status)
   else
     write(*,*) '  ***************************************************************************'
     write(*,*) "    ERROR: The ",trim(auxstr1),' file does not exist. '
     write(*,*) '           Please, provide it.'
     write(*,*) '  ***************************************************************************'
     write(*,*)
     call exit(1)
   endif     
   
   
   IF (status /= 0) THEN
       call exit(1)
       write(*,*) 'Error in reading wavefunction file'
   END IF
  ! WRITE(*,*) 'Start to read initial state wavefunction file'

   !! Read the head
   do ii = 1, 6+n_atom_type
      read(3,*) temp_ch1
   end do
   do ii = 1, nat
      read(3,*) temp_ch2, array_temp(1), array_temp(2), array_temp(3), ll
   end do

   temp_int = nfft*nfft*nfft

  ! ORIGINAL VERSION: IT GIVES SEGFAULTS
  ! temp1 = int(real(temp_int)/5.0)
  ! temp2 = mod(temp_int,5)
  !
  ! CORRECTION BY PGR 15 AUGUST 2017
   temp2 = mod(temp_int,5)
   temp1 = int(real(temp_int-temp2)/5.0)

   
   do ii = 1, temp1
      read(3,*) (wfc_temp((ii-1)*5+jj), jj=1,5)
!      write(*,*) 'I=', ii, 'N=', temp1
   end do

 
   if(temp2 /= 0) then
      read(3,*) (wfc_temp(temp1*5+jj), jj=1,temp2)
   end if

   ll = 1
   do kk = 1, nfft  !! for z
      do jj = 1, nfft  !! for y
         do ii = 1, nfft   !! for x
            if(wfc_temp(ll) > 0.0) then
               wfc_i(ii,jj,kk) = sqrt(wfc_temp(ll))
            else
               wfc_i(ii,jj,kk) = -sqrt(abs(wfc_temp(ll)))
            end if
            ll = ll+1
         end do
      end do
   end do

   CLOSE(UNIT = 3)

   
   deallocate(wfc_temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Normalization of wfc
   temp = 0.d0
   temp_re = 0.d0
   do kk = 1, nfft
      do jj = 1, nfft
         do ii = 1, nfft
            temp = temp + wfc_i(ii,jj,kk)**2
         end do
      end do
   end do

   ! Beware of this: the normalization condition below means that the wavefunction is normalized in Fourier (reciprocal) space, not in real space.
   temp = sqrt(temp) !OLDER: sqrt(temp*(lattparam/dble(nfft))**3)

   do kk = 1, nfft
      do jj = 1, nfft
         do ii = 1, nfft
            wfc_i(ii,jj,kk) = wfc_i(ii,jj,kk)/temp
         end do
      end do
   end do
   
   
   return

end subroutine read_wfc



!-------------------------------------------------------------------------------------------


subroutine wfc_product (nfft, wfc_i, wfc_f, prod)

 ! This subroutine reads the <istate>-th wavefunction and it stores it in the <wfc_i> variable (normalized).
 ! The wavefunction is stored in <wfc_directory>, and its name is "<Input_wfc_prefix//istate//.wfn>"

!use mytypes
!use elph_pwscf_mod, only: wfc_i, nfft, Input_wfc, nat,n_atom_type, lattparam, wfc_directory

implicit none

    integer, Intent(In) :: nfft
    real(8), Intent(In) :: wfc_i(nfft,nfft,nfft), wfc_f(nfft,nfft,nfft) 
    real(8), Intent(Out) :: prod

   !! local variables
   integer :: ii, jj, kk

   prod=0.d0
   do kk = 1, nfft
      do jj = 1, nfft
         do ii = 1, nfft
            prod = prod + (wfc_i(ii,jj,kk)) * wfc_f(ii,jj,kk) !dcmplx(wfc_i(ii,jj,kk)) * wfc_f(ii,jj,kk)
         end do
      end do
   end do

end subroutine

!-------------------------------------------------------------------------------------------


!  ======================================================================================================
!  ======================================================================================================
!  ======================================================================================================
!   
!      Examples of the input files for Quantum Espresso to generate the necessary files can be viewed below:
!
!
!      QE geometry optimization (pw.x):
! 
! 		   &CONTROL
! 			   calculation = 'relax',
! 			   restart_mode = 'from_scratch',
! 			   prefix='',
! 			   outdir = './',
! 			   pseudo_dir = '/home/risueno/CalcsPengProject/QE_files/PP/',
! 			   forc_conv_thr = 1.0D-11 ,
! 			   etot_conv_thr = 1.0D-12 ,
! 			/
! 		   &system
! 			   ibrav = 0, a=18.0, 
! 			   nat= 26, ntyp= 2,
! 			   ecutwfc = 30d0,
! 			   nbnd = 600, 
! 		   /
! 		   &electrons
! 			   conv_thr = 1.0e-15,
! 			   mixing_beta = 0.7,
! 			   mixing_mode = 'plain',
! 			   diagonalization = 'david'
! 		   /
! 		   &IONS
! 		   /
! 
! 		   ATOMIC_SPECIES
! 		   C   12.0107   C.pz-vbc.UPF
! 		   H   1.007825035  H.pz-vbc.UPF
! 
! 		   ATOMIC_POSITIONS { angstrom }
! 		   C        8.905374526   8.905374526   8.905374526
! 		   C        9.783416934   9.783416934   8.032419367
! 		   C        9.783416934   8.032419367   9.783416934
! 		   C        8.032419367   9.783416934   9.783416934
! 		   C        8.905361195  10.661436254  10.661436254
! 		   C       10.661436254  10.661436254   8.905361195
! 		   C       10.661436254   8.905361195  10.661436254
! 		   C        9.783417452  11.534334298   9.783417452
! 		   C        9.783417452   9.783417452  11.534334298
! 		   C       11.534334298   9.783417452   9.783417452
! 		   H        8.268394663   8.268394663   8.268394663
! 		   H        7.376027810  10.411288572   9.155507252
! 		   H        7.376027810   9.155507252  10.411288572
! 		   H        9.155507252  10.411288572   7.376027810
! 		   H       10.411288572   9.155507252   7.376027810
! 		   H        9.155507252   7.376027810  10.411288572
! 		   H       10.411288572   7.376027810   9.155507252
! 		   H        9.155523979  12.190725092   9.155523979
! 		   H       12.190725092   9.155523979   9.155523979
! 		   H        9.155523979   9.155523979  12.190725092
! 		   H        8.268390473  11.298415367  11.298415367
! 		   H       11.298415367  11.298415367   8.268390473
! 		   H       11.298415367   8.268390473  11.298415367
! 		   H       10.411305333  12.190722482  10.411305333
! 		   H       10.411305333  10.411305333  12.190722482
! 		   H       12.190722482  10.411305333  10.411305333
! 
! 		   CELL_PARAMETERS {cubic}
! 			1.00  0.00  0.00
! 			0.00  1.00  0.00
! 			0.00  0.00  1.00 
! 
! 		   K_POINTS
! 			1
! 			0.0 0.0 0.0  1.0
!
!  ======================================================================================================
!
!      QE undisplaced ground state calculation (pw.x)
! 
! 		  &CONTROL
! 			 calculation = 'scf',
! 			 disk_io='low',
! 			 tprnfor=.true.,
! 			 nstep = 2000,
! 			 restart_mode='from_scratch',
! 			 prefix='C10H16-', 
! 			 outdir='./',
! 			 pseudo_dir='/home/risueno/CalcsPengProject/QE_files/PP/',
! 		   /
!  
! 		  &system
! 			 ibrav=0,a=18.0,
! 			 nat=26,ntyp=2,
! 			 ecutwfc=30d0,
! 			 nbnd=40,
! 		  /
!  
! 		  &electrons
! 			 conv_thr=1.0e-15,
! 			 mixing_beta=0.7,
! 			 mixing_mode='plain',
! 			 diagonalization='david'
! 		  /
!  
!  
! 		  ATOMIC_SPECIES
! 		  C   12.0107   C.pz-vbc.UPF
! 		  H   1.007825035  H.pz-vbc.UPF
!  
! 		  ATOMIC_POSITIONS { angstrom }
! 		  C      8.89842525100000        8.89842525100000        8.89842525100000     
! 		  C      9.78330467200000        9.78330467200000        8.01966606300000     
! 		  C      9.78330467200000        8.01966606300000        9.78330467200000     
! 		  C      8.01966606300000        9.78330467200000        9.78330467200000     
! 		  C      8.89843542000000        10.6681787870000        10.6681787870000     
! 		  C      10.6681787870000        10.6681787870000        8.89843542000000     
! 		  C      10.6681787870000        8.89843542000000        10.6681787870000     
! 		  C      9.78332034500000        11.5470749970000        9.78332034500000     
! 		  C      9.78332034500000        9.78332034500000        11.5470749970000     
! 		  C      11.5470749970000        9.78332034500000        9.78332034500000     
! 		  H      8.25068775200000        8.25068775200000        8.25068775200000     
! 		  H      7.35210328600000        10.4218786750000        9.14510012400000     
! 		  H      7.35210328600000        9.14510012400000        10.4218786750000     
! 		  H      9.14510012400000        10.4218786750000        7.35210328600000     
! 		  H      10.4218786750000        9.14510012400000        7.35210328600000     
! 		  H      9.14510012400000        7.35210328600000        10.4218786750000     
! 		  H      10.4218786750000        7.35210328600000        9.14510012400000     
! 		  H      9.14511714000000        12.2149969440000        9.14511714000000     
! 		  H      12.2149969440000        9.14511714000000        9.14511714000000     
! 		  H      9.14511714000000        9.14511714000000        12.2149969440000     
! 		  H      8.25084989100000        11.3160468160000        11.3160468160000     
! 		  H      11.3160468160000        11.3160468160000        8.25084989100000     
! 		  H      11.3160468160000        8.25084989100000        11.3160468160000     
! 		  H      10.4217698150000        12.2146243630000        10.4217698150000     
! 		  H      10.4217698150000        10.4217698150000        12.2146243630000     
! 		  H      12.2146243630000        10.4217698150000        10.4217698150000     
!  
! 		  CELL_PARAMETERS {cubic}
! 		   1.00  0.00  0.00
! 		   0.00  1.00  0.00
! 		   0.00  0.00  1.00
!  
! 		  K_POINTS
! 		   1
! 		   0.0 0.0 0.0  1.0
! 
! 	
!  ======================================================================================================
!
!   QE solution of the dynamical equation to find phonon frequencies and displacement eigenvectors. (ph.x).
!   The eigenvectors of all eigenvalues (except the 6 lowest ones) are added to the relaxed positions to find
!   the eigenvalues after + and - displacements (generate these QE input files e.g. with vib_potPGR.f90). 
!	 
! 	  Phonons at Gamma
! 	   &inputph
! 		tr2_ph=1.0d-14,
! 		amass(1)=12.0107,
! 		amass(2)=1.007825035,
! 		outdir='./'
! 		fildyn='dynG',
! 		epsil=.true.
! 	   /
! 	   0.0 0.0 0.0
!
!
!  ======================================================================================================
!
!      QE displaced calculation (pw.x)
!
! 		   &CONTROL
! 			   wfcdir='./',
! 			   calculation='scf',
! 			   disk_io='low',
! 			   tprnfor=.true.,
! 			   nstep=2000,
! 			   restart_mode='from_scratch',
! 			   prefix='C10H16-',
! 			   outdir='./',
! 			   pseudo_dir='/home/risueno/CalcsPengProject/QE_files/PP/',
! 			  /
!  
!  
! 			 &system
! 			   ibrav=0,a=18.0,
! 			   nat=26,ntyp=2,
! 			   ecutwfc=30d0,
! 			   nbnd=40,
! 			 /
!  
!  
! 			 &electrons
! 			   conv_thr=1.0e-12,
! 			   mixing_beta=0.7,
! 			   mixing_mode='plain',
! 			   diagonalization='david'
! 			 /
!  
!  
!  
! 			 ATOMIC_SPECIES
! 			 C   12.0107   C.pz-vbc.UPF
! 			 H   1.007825035  H.pz-vbc.UPF
!  
! 			 ATOMIC_POSITIONS { angstrom }
! 			C      8.11407921699733        8.11539414861903        8.11515204238364     
! 			C      8.99975092922971        8.99975751210923        7.23664013493481     
! 			C      8.99975123968866        7.23709845282974        8.99975660947859     
! 			C      7.23460925623545        8.99975613229169        8.99975491920210     
! 			C      8.11409429009884        9.88410950415149        9.88435190934735     
! 			C      9.88542541314396        9.88410884298892        8.11516050385943     
! 			C      9.88542526366373        8.11540112104173        9.88435109870455     
! 			C      8.99976952859275        10.7625418654879        8.99976457274808     
! 			C      8.99977157532210        8.99976549262644        10.7630004191017     
! 			C      10.7650323556612        8.99977454767906        8.99977566303157     
! 			H      7.46709850192498        7.46716315787528        7.46715125119974     
! 			H      6.56707411008172        9.63849900771107        8.36167260578857     
! 			H      6.56706920598021        8.36138477010015        9.63820371525868     
! 			H      8.36233908640736        9.63912253001040        6.56908192283118     
! 			H      9.63751872783712        8.36077576463121        6.56905528085415     
! 			H      8.36220165083098        6.56953012206413        9.63896794445182     
! 			H      9.63765943473092        6.56950838418861        8.36092617049248     
! 			H      8.36221876456806        11.4304678714512        8.36092420825205     
! 			H      11.4329480819837        8.36141411247019        8.36170625433919     
! 			H      8.36236276002700        8.36077605034352        11.4309245335651     
! 			H      7.46726878759771        10.5324722494833        10.5324848863123     
! 			H      10.5325373251281        10.5324686964531        7.46731055582454     
! 			H      10.5325365202346        7.46731998457774        10.5324805283886     
! 			H      9.63755355283711        11.4301192072886        9.63883972676156     
! 			H      9.63741918160552        9.63899561164832        11.4305647101281     
! 			H      11.4325388265833        9.63840019163292        9.63810993551456     
!  
! 			 CELL_PARAMETERS {cubic}
! 			  1.00  0.00  0.00
! 			  0.00  1.00  0.00
! 			  0.00  0.00  1.00
!  
! 			 K_POINTS
! 			  1
! 			  0.0 0.0 0.0  1.0
!
!  ======================================================================================================
!
!      QE printing of wavefunctions calculation (pp.x):
! 
! 			  &InputPP
! 				  filplot = '35.wfn',
! 				  kband=35,
! 				 prefix='C10H16-',
! 				 outdir='./',
! 				 plot_num=7,
! 				 kpoint=1,
! 				 lsign=.true.,
! 				 spin_component=0,
! 			  /
! 			  &plot
! 				 iflag=3,
! 				 output_format=5
! 			  /
!
!
!  ======================================================================================================
!  ======================================================================================================
!  ======================================================================================================

