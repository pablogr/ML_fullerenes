program extract_data_fullerenes

!     This program generates the input files that are necessary to run the Maat program. 
!     This program requires an input file called 'extract_data_fullerenes.in'. An example is:
! 
! # THIS IS THE INPUT FILE FOR THE PROGRAM THAT
! # GENERATES THE INPUT OF Maat (extract_data_fullerenes.x)
! 
! Number_of_atoms = 26
! Number_of_species = 2
! 
! elec_eigval_file_name = out_scf.out
! Band_occupation =2
! 
! Phonon_Frequency_File_Name = Input_dot
! phonon_frequency_units=cm-1
! 
!     
!
!     This program generates:
!    
!     The Maat program requires an input file specifying system parameters (input.in) as well as the following input data files:
!     * wQ0.dat: 1st row: number of q-points and the number of phonon frequencies; then: q-index, branch-index, phonon frequency
!     * epsilonorig0.dat: 1st row: number of bands, number of k-points, number of spins; then all three indices and the corresponding unrenormalized electronic eigenvalue
! 
!     Comment: Ponce2014 considers hbar = m_e = e = 1. This means that, if we use the same convention -atomic units-,
!     we must multiply the masses expressed with the unified atomic mass pattern by 1822.888486192.
!     Our results of renormalized electronic eigenvalues will be given with the same units that were given by 
!     the original files for the electronic eigenvalues (e.g. eV); this means that if the original phonon freqs.
!     were given in another units, a conversion is necessary (such conversion is automatically done by the  change_units_omega subroutine).
!
!     Author:   Pablo Garcia Risueno (2019).



  implicit none
  logical :: exist
  integer :: i, j, k, l, h, i1, i2, i3, sizeepsilon, NstatesforgDW
  integer :: gFan_i_index_beg, gFan_i_index_end, gFan_j_index_beg, gFan_j_index_end,  gDW_i_index_beg, gDW_i_index_end
  integer :: number_of_atoms, number_of_species, min_band_index, max_band_index, dimnu, Ne, band_occupation, n_mode_beg,n_mode_end
  real(8) :: auxr1, auxr2, auxr3
  real(8) :: meanepsilon, varepsilon, skewepsilon, kurtosisepsilon
  real(8) :: meanomega, varomega, skewomega, kurtosisomega
  complex(8) :: auxc1, auxc2, auxc3
  real(8),allocatable :: auxrv1(:), xi(:,:), U(:,:), omega(:), omegaorig(:), epsilonorig(:), epsilon(:), mass(:,:)
  complex(8),allocatable :: auxcv1(:), gFan(:,:,:,:,:,:), gDW(:,:,:,:,:)
  character(100) :: elec_eigval_file_name, phon_freq_file_name, gFan_file_name !, mass_file_name
  character(40) :: units_omegaorig, units_epsilonorig, fo
  character(16) :: cad1, cad2, cad3, cad4, cad5, cad6, cad7, cad8, cad9, cad10, cad11, cad12, cad13, cad14
  character(16) :: cad15, cad16, cad17, cad18, cad19, cad20, cad21, cad22, cad23, cad24, cad25, cad26
  character(416) :: longcad
  
  call read_input_file (number_of_atoms, number_of_species, dimnu, band_occupation, units_omegaorig, &
                         &  elec_eigval_file_name, phon_freq_file_name, gFan_file_name, &
                         &  gFan_i_index_beg, gFan_i_index_end, gFan_j_index_beg, gFan_j_index_end, &
                         &  gDW_i_index_beg, gDW_i_index_end,n_mode_beg,n_mode_end,NstatesforgDW)
        
  call read_QE_occ_unocc_numbers (elec_eigval_file_name, band_occupation, Ne, sizeepsilon)  
                                         
  allocate(omega(dimnu),omegaorig(dimnu),xi(dimnu,dimnu),U(dimnu,dimnu))
  allocate(epsilon(sizeepsilon),epsilonorig(sizeepsilon),mass(number_of_atoms,2))
  allocate(gFan(gFan_i_index_end-gFan_i_index_beg+1, gFan_j_index_end-gFan_j_index_beg+1, 1,1,1, dimnu))
  allocate(gDW(gFan_i_index_end-gFan_i_index_beg+1,1,1,1,dimnu))                       
             
  call read_QE_output (sizeepsilon, Ne, number_of_atoms, number_of_species, band_occupation, &
                      & elec_eigval_file_name, epsilonorig, units_epsilonorig, mass) 

  call read_dyneq_output (number_of_atoms, dimnu,  n_mode_beg, n_mode_end,  phon_freq_file_name, mass, omegaorig, xi, U)   
 
  call omega_epsilon_to_atomic_units (units_omegaorig,units_epsilonorig,dimnu,sizeepsilon, &        
                                             & omegaorig,epsilonorig,omega,epsilon)


  call      mean_omega(dimnu, omega, meanomega)
  call       var_omega(dimnu, omega, meanomega, varomega)
  call  skewness_omega(dimnu, omega, meanomega, varomega, skewomega)
  call  kurtosis_omega(dimnu, omega, meanomega, varomega, skewomega, kurtosisomega)
  write(*,*)
  ! Note that we calculate these quantities just for the occupied states 
  ! (otherwise the result would depend on the arbitrary number of unoccupied states chosen).
  call     mean_epsilon(Ne/2, epsilon, meanepsilon)
  call      var_epsilon(Ne/2, epsilon, meanepsilon, varepsilon)
  call skewness_epsilon(Ne/2, epsilon, meanepsilon, varepsilon, skewepsilon)
  call kurtosis_epsilon(Ne/2, epsilon, meanepsilon, varepsilon, skewepsilon, kurtosisepsilon)
  write(*,*)

  
 ! fo = '(F16.4,a)'
 ! write(*,fo)omega(7),omega(dimnu),meanomega,varomega,skewomega,kurtosisomega,meanepsilon,varepsilon,skewepsilon,kurtosisepsilon

  cad1=''; cad2='';longcad=''
  write(cad1,'(F16.4)') omega(7); call StripSpaces(cad1) 
  write(cad2,'(F16.4)') omega(dimnu); call StripSpaces(cad2) 
  write(cad3,'(F16.4)') meanomega; call StripSpaces(cad3) 
  write(cad4,'(F16.4)') varomega; call StripSpaces(cad4) 
  write(cad5,'(F16.4)') skewomega; call StripSpaces(cad5) 
  write(cad6,'(F16.4)') kurtosisomega; call StripSpaces(cad6) 
  write(cad7,'(F16.4)') meanepsilon; call StripSpaces(cad7) 
  write(cad8,'(F16.4)') varepsilon; call StripSpaces(cad8) 
  write(cad9,'(F16.4)') skewepsilon; call StripSpaces(cad9) 
  write(cad10,'(F16.4)') kurtosisepsilon; call StripSpaces(cad10) 
  write(cad11,'(F16.4)') epsilon(Ne/2)-epsilon(Ne/2-1); call StripSpaces(cad11)
  write(cad12,'(F16.4)') epsilon(Ne/2)-epsilon(Ne/2-2); call StripSpaces(cad12) 
  write(cad13,'(F16.4)') epsilon(Ne/2)-epsilon(Ne/2-3); call StripSpaces(cad13) 
  write(cad14,'(F16.4)') epsilon(Ne/2)-epsilon(Ne/2-4); call StripSpaces(cad14) 
  write(cad15,'(F16.4)') epsilon(Ne/2)-epsilon(Ne/2-5); call StripSpaces(cad15) 
  write(cad16,'(F16.4)') epsilon(Ne/2)-epsilon(Ne/2-6); call StripSpaces(cad16) 
  write(cad17,'(F16.4)') epsilon(Ne/2)-epsilon(Ne/2-7); call StripSpaces(cad17) 
  write(cad18,'(F16.4)') epsilon(Ne/2)-epsilon(Ne/2-8); call StripSpaces(cad18)  
  write(cad19,'(F16.4)') epsilon(Ne/2+2)-epsilon(Ne/2+1); call StripSpaces(cad19)
  write(cad20,'(F16.4)') epsilon(Ne/2+3)-epsilon(Ne/2+1); call StripSpaces(cad20)
  write(cad21,'(F16.4)') epsilon(Ne/2+4)-epsilon(Ne/2+1); call StripSpaces(cad21)
  write(cad22,'(F16.4)') epsilon(Ne/2+5)-epsilon(Ne/2+1); call StripSpaces(cad22)
  write(cad23,'(F16.4)') epsilon(Ne/2+6)-epsilon(Ne/2+1); call StripSpaces(cad23)
  write(cad24,'(F16.4)') epsilon(Ne/2+7)-epsilon(Ne/2+1); call StripSpaces(cad24)
  write(cad25,'(F16.4)') epsilon(Ne/2+8)-epsilon(Ne/2+1); call StripSpaces(cad25)
  write(cad26,'(F16.4)') epsilon(Ne/2+9)-epsilon(Ne/2+1); call StripSpaces(cad26)

  
    
  write(longcad,*) trim(cad1)//'  '//trim(cad2)//'  '//trim(cad3)//'  '//trim(cad4)//'  ' &
     & //trim(cad5)//'  '//trim(cad6)//'  '//trim(cad7)//'  '//trim(cad8)//'  '//trim(cad9)//'  '//trim(cad10)//'  '&
     & //trim(cad11)//'  '//trim(cad12)//'  '//trim(cad13)//'  '//trim(cad14)//'  '//trim(cad15)//'  '//trim(cad16)//'  '&
     & //trim(cad17)//'  '//trim(cad18)//'  '//trim(cad19)//'  '//trim(cad20)//'  '//trim(cad21)//'  '//trim(cad22)//'  '&
     & //trim(cad23)//'  '//trim(cad24)//'  '//trim(cad25)//'  '//trim(cad26) 
 
  write(*,*)longcad

  deallocate(omega,omegaorig,xi,U,epsilon,epsilonorig,mass,gFan,gDW)                       
  
  call exit(0)
 
 
end program extract_data_fullerenes



!-------------------------------------------------------------------------------------------


subroutine mean_epsilon(N, variable, mean)

 integer, Intent(In) ::  N  			! Number of cells of the input vector
 real(8), Intent(In) ::  variable(N)    ! Variable to average
 real(8), Intent(Out) :: mean      	    ! Average of the input variable

 integer :: i


   mean = 0.d0

   do i=1,N
     mean = mean + variable(i) 
   end do
   
   mean = mean / real(N)
   
   write(*,*) "The average epsilon is ", mean

end subroutine mean_epsilon



subroutine mean_omega(N, variable, mean)

 integer, Intent(In) ::  N  			! Number of cells of the input vector
 real(8), Intent(In) ::  variable(N)    ! Variable to average
 real(8), Intent(Out) :: mean      	    ! Average of the input variable

 integer :: i


   mean = 0.d0

   do i=7,N
     mean = mean + variable(i) 
   end do
   
   mean = mean / real(N-6)
   
   write(*,*) "The average omega is ", mean

end subroutine mean_omega

!-------------------------------------------------------------------------------------------
 
 

subroutine var_epsilon(N, variable, mean, variance)

 integer, Intent(In) ::  N  		 ! Number of cells of the input vector
 real(8), Intent(In) ::  variable(N) ! Variable whose variance we calculate
 real(8), Intent(In) ::  mean	     ! Average of the input variable
 real(8), Intent(Out) :: variance    ! Variance of the input variable

 integer :: i


   variance = 0.d0

   do i=1,N
     variance = variance + ( variable(i) - mean )**2
   end do
   
   variance = variance / real(N)
   
   write(*,*) "The variance of epsilon is ", variance

end subroutine var_epsilon



subroutine var_omega(N, variable, mean, variance)

 integer, Intent(In) ::  N  		 ! Number of cells of the input vector
 real(8), Intent(In) ::  variable(N) ! Variable whose variance we calculate
 real(8), Intent(In) ::  mean	     ! Average of the input variable
 real(8), Intent(Out) :: variance    ! Variance of the input variable

 integer :: i


   variance = 0.d0

   do i=7,N
     variance = variance + ( variable(i) - mean )**2
   end do
   
   variance = variance / real(N-6)
   
   write(*,*) "The variance of omega is ", variance

end subroutine var_omega
!-------------------------------------------------------------------------------------------

subroutine skewness_epsilon(N, variable, mean, variance, skewness)

 integer, Intent(In) ::  N  			! Number of cells of the input vector
 real(8), Intent(In) ::  variable(N)    ! Variable to average
 real(8), Intent(In) ::  mean      	    ! Average of the input variable
 real(8), Intent(In) ::  variance       ! Variance of the input variable
 real(8), Intent(Out) :: skewness       ! Skewness of the input variable
  
 integer :: i


   skewness = 0.d0

   do i=1,N
     skewness = skewness + ( variable(i) - mean )**3
   end do
   
   skewness = skewness / ( real(N) * (sqrt(variance))**3 )
   
   write(*,*) "The skewness of epsilon is ", skewness

end subroutine skewness_epsilon


subroutine skewness_omega(N, variable, mean, variance, skewness)

 integer, Intent(In) ::  N  			! Number of cells of the input vector
 real(8), Intent(In) ::  variable(N)    ! Variable to average
 real(8), Intent(In) ::  mean      	    ! Average of the input variable
 real(8), Intent(In) ::  variance       ! Variance of the input variable
 real(8), Intent(Out) :: skewness       ! Skewness of the input variable
  
 integer :: i


   skewness = 0.d0

   do i=7,N
     skewness = skewness + ( variable(i) - mean )**3
   end do
   
   skewness = skewness / ( real(N-6) * (sqrt(variance))**3 )
   
   write(*,*) "The skewness of omega is ", skewness

end subroutine skewness_omega

!-------------------------------------------------------------------------------------------


subroutine kurtosis_epsilon(N, variable, mean, variance, skewness, kurtosis)

 integer, Intent(In) ::  N  			! Number of cells of the input vector
 real(8), Intent(In) ::  variable(N)    ! Variable to average
 real(8), Intent(In) ::  mean      	    ! Average of the input variable
 real(8), Intent(In) ::  variance       ! Variance of the input variable
 real(8), Intent(In) ::  skewness       ! Skewness of the input variable
 real(8), Intent(Out) :: kurtosis       ! Kurtosis of the input variable
  
 integer :: i


   kurtosis = 0.d0

   do i=1,N
     kurtosis = kurtosis + ( variable(i) - mean )**4
   end do
   
   kurtosis = kurtosis / ( real(N) * (variance)**2 )
   
   write(*,*) "The kurtosis of epsilon is ", kurtosis

end subroutine kurtosis_epsilon


subroutine kurtosis_omega(N, variable, mean, variance, skewness, kurtosis)

 integer, Intent(In) ::  N  			! Number of cells of the input vector
 real(8), Intent(In) ::  variable(N)    ! Variable to average
 real(8), Intent(In) ::  mean      	    ! Average of the input variable
 real(8), Intent(In) ::  variance       ! Variance of the input variable
 real(8), Intent(In) ::  skewness       ! Skewness of the input variable
 real(8), Intent(Out) :: kurtosis       ! Kurtosis of the input variable
  
 integer :: i


   kurtosis = 0.d0

   do i=7,N
     kurtosis = kurtosis + ( variable(i) - mean )**4
   end do
   
   kurtosis = kurtosis / ( real(N-6) * (variance)**2 )
   
   write(*,*) "The kurtosis of omega is ", kurtosis

end subroutine kurtosis_omega

!-------------------------------------------------------------------------------------------
 


  
!-------------------------------------------------------------------------------------------
  
  
subroutine read_input_file (number_of_atoms, number_of_species, dimnu, band_occupation, units_omegaorig, &
                         &  elec_eigval_file_name, phon_freq_file_name, gFan_file_name, &
                         & gFan_i_index_beg, gFan_i_index_end, gFan_j_index_beg, gFan_j_index_end, &
                         & gDW_i_index_beg, gDW_i_index_end,n_mode_beg,n_mode_end,NstatesforgDW)

 ! This subroutine reads the extract_data_fullerenes.in file, which contains the names of the files where one must read the data
 ! as well as other relevant data. 

  ! System parameters
  integer, Intent(Out) :: number_of_atoms                ! The number of atoms of the system
  integer, Intent(Out) :: number_of_species                ! The number of different atomic species of the system
  integer, Intent(Out) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  integer, Intent(Out) :: band_occupation                ! number of electrons per band of the ones given in the output file of Quantum Espresso; it can be 1 or 2
  character(40) :: units_omegaorig						 ! units of the frequencies stored in the phon_freq_file_name file (Input_dot)
  character(100), Intent(Out) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  character(100), Intent(Out) :: phon_freq_file_name     ! Name of the file where the phonon frequencies and eigenvectors of the dynamical matrix are stored
  character(100), Intent(Out) :: gFan_file_name          ! Name of the file where the Fan electron-phonon matrix elements are stored
  Integer, Intent(Out)        :: gFan_i_index_beg         ! Initial wavefunction index for < psi_i | for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(Out)        :: gFan_i_index_end         ! Final wavefunction index < psi_i | for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(Out)        :: gFan_j_index_beg         ! Initial wavefunction index for | psi_j >for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(Out)        :: gFan_j_index_end         ! Final wavefunction index  | psi_j >  for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(Out)        :: gDW_i_index_beg         ! Initial index for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(Out)        :: gDW_i_index_end         ! Final index for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(Out)        :: n_mode_beg,n_mode_end    ! Initial and final indices of phonon branch
  Integer, Intent(Out)        :: NstatesforgDW            ! Number of states (bands) considered in the summation to calculate gDW
  
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT, auxstr4, auxstr5
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,iostat,unitt !unitt=1:meV; unitt=2: Hartree; unitt=3:K
  logical :: exist
  real(8) :: auxr1, auxr2


  number_of_atoms = 0; number_of_species = 0; units_omegaorig='cm**-1'
  gFan_i_index_beg=1; gFan_i_index_end=1; gFan_j_index_beg=1; gFan_j_index_end=1
  gDW_i_index_beg=0; gDW_i_index_end=0; band_occupation=2; NstatesforgDW=-100
  
  !gFan_file_name='gFan.dat'

  ! Reading parameters from file
  
   write(*,*); write(*,*); 
   write(*,*) '  ***************************************************************************' 
   write(*,*) '  ***************************  NOW RUNNING   ********************************'   
   write(*,*) '  *************************** extract_data_fullerenes.x ********************************' 
   write(*,*) '  *************** to generate the input for the Maat program ****************'      
   write(*,*) '  ***************************************************************************'   
   write(*,*) ;  write(*,*)
   write(*,*) '  ********* FROM THE INPUT FILE OF extract_data_fullerenes.x (extract_data_fullerenes.in): ********'   


   inquire(file='extract_data_fullerenes.in', exist=exist)      
   if (.not. exist) then
	   write(*,*) 
	   write(*,*) '  ***************************************************************************' 
	   write(*,*) '    ERROR: The file "extract_data_fullerenes.in" does not exist. Please, provide it.'
	   write(*,*) '  ***************************************************************************' 
	   write(*,*)
	   call exit(1) 
   end if




   open(345,file='extract_data_fullerenes.in',status='unknown')

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



     if ((auxstr1 .eq. 'elec_eigval_file_name').or.(auxstr1 .eq. 'Elec_eigval_file_name'))   then
         elec_eigval_file_name = trim(auxstr2)
         call StripSpaces (elec_eigval_file_name)
     end if
     
     if ((auxstr1 .eq. 'phonon_frequency_file_name').or.(auxstr1 .eq. 'Phonon_frequency_file_name') &
      & .or.(auxstr1 .eq. 'Phonon_Frequency_File_Name') )   then
         phon_freq_file_name = trim(auxstr2)
         call StripSpaces (phon_freq_file_name)
     end if
          
     if ((auxstr1 .eq. 'gFan_file_name').or.(auxstr1 .eq. 'gFan_file_name'))   then
         gFan_file_name = trim(auxstr2)
         call StripSpaces (gFan_file_name)
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

    if ((auxstr1 .eq. 'band_occupation').or.(auxstr1 .eq. 'Band_Occupation').or. &
             & (auxstr1 .eq. 'Band_occupation'))  Read(auxstr2, '(I3)' ) band_occupation
             
    if ((auxstr1 .eq. 'Number_of_species').or.(auxstr1 .eq. 'number_of_species').or. &
             & (auxstr1 .eq. 'Number_species'))  Read(auxstr2, '(I3)' ) number_of_species     
             
    if ((auxstr1 .eq. 'Number_of_bands_for_DW').or.(auxstr1 .eq. 'number_of_bands_for_DW').or. &
             & (auxstr1 .eq. 'number_of_bands_for_DW'))  Read(auxstr2, '(I8)' )   NstatesforgDW    
             
    if ((auxstr1 .eq. 'min_band_index').or.(auxstr1 .eq. 'Min_band_index'))  Read(auxstr2, '(I3)' ) min_band_index 
    if ((auxstr1 .eq. 'max_band_index').or.(auxstr1 .eq. 'Max_band_index'))  Read(auxstr2, '(I3)' ) max_band_index
  
!     if ((auxstr1 .eq. 'Initial_gFan_Index_i').or.(auxstr1 .eq. 'initial_gFan_index_i')) &
!                                                                             & Read(auxstr2, '(I6)' ) gFan_i_index_beg
!     if ((auxstr1 .eq. 'Final_gFan_Index_i').or.(auxstr1 .eq. 'final_gFan_index_i')) &
!                                                                             & Read(auxstr2, '(I6)' ) gFan_i_index_end
!     if ((auxstr1 .eq. 'Initial_gFan_Index_j').or.(auxstr1 .eq. 'initial_gFan_index_j')) &
!                                                                             & Read(auxstr2, '(I6)' ) gFan_j_index_beg
!     if ((auxstr1 .eq. 'Final_gFan_Index_j').or.(auxstr1 .eq. 'final_gFan_index_j')) &
!                                                                             & Read(auxstr2, '(I6)' ) gFan_j_index_end 
! 
    if ((auxstr1 .eq. 'Initial_gDW_Index').or.(auxstr1 .eq. 'initial_gdw_index')) &
                                                                            & Read(auxstr2, '(I6)' ) gDW_i_index_beg
    if ((auxstr1 .eq. 'Final_gDW_Index').or.(auxstr1 .eq. 'final_gdw_index')) &
                                                                            & Read(auxstr2, '(I6)' ) gDW_i_index_end
  
  end do ! k
   

1117 continue



   
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




  
  if ((band_occupation .ne. 1) .and. (band_occupation .ne. 2)) then
    band_occupation = 2
    write(*,*)
    write(*,*) '  <<<<<<<<<<<< WARNING: The number of electrons per band in ',trim(elec_eigval_file_name)
    write(*,*) '  was assumed to be 2. >>>>>>>>>>>'
    write(*,*)   
  end if
  
   
  dimnu = number_of_atoms*3 ! Strictly speaking it is number_of_atoms*3-6, but we make it bigger to avoid problems of lacks of space (some freqs are 0 and they are anyway read)
 
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
    write(*,*) '         in extract_data_fullerenes.in). Accepted values are cm**-1, eV, meV and  '
    write(*,*) '         a.u. (Hartree). Please set the variable to one of these values.  '
    call exit(1)
  end if
 

  if (number_of_atoms .le. 0) then
     write(*,*) ' ERROR: Please, specify a valid number of atoms in extract_data_fullerenes.in; '
     write(*,*) '        E.g.: Number_of_atoms =  2 '
     call exit(1)
  end if

  if (number_of_species .le. 0) then
     write(*,*) ' ERROR: Please, specify a valid number of different atomic species in extract_data_fullerenes.in; '
     write(*,*) '        E.g.: Number_of_species =  1 '
     call exit(1)
  end if
  
  write(*,*) 
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
  
  write(*,'(A,I5)') '                        Number of different atomic species = ', number_of_species
  write(*,'(A,I5)') '                                           Number of atoms = ', number_of_atoms
  write(auxstr4,*) dimnu-6; call StripSpaces(auxstr4)
  write(auxstr5,*) dimnu; call StripSpaces(auxstr5) 
  write(*,*)        '                                Number of phonon branches =    ', &
           & trim(auxstr5), ' (',trim(auxstr4),' valid)'
  write(*,*)
  



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
      
 subroutine read_QE_occ_unocc_numbers (elec_eigval_file_name, band_occupation, Ne, sizeepsilon)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives the numbers of occ and unocc states.

   ! System parameters
  character(100), Intent(In) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  integer, Intent(In)  :: band_occupation               ! 1 or 2, number of electrons per band in the QE output file
  integer, Intent(Out) :: Ne                            ! The number of electrons of the system (i.e. the number of occupied states)
  integer, Intent(Out) :: sizeepsilon                   ! size of the vector storing the electronic eigenvalues
  
     
    ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  integer :: Nks     ! Number of Kohn-Sham states
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9  
  
    
  write(*,*)      
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
   write(*,'(A,I5)') '                                       Number of electrons = ', Ne

  close(433)
  
end do ! m
  
 2117 continue 
  
  
  
  
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
   
   read(auxstr2,*) Nks
  write(*,'(A,I5)') '                                Number of Kohn-Sham states = ', Nks
  write(*,*)
      
  close(451)
  
end do ! m
  
  3117 continue
  
  if (band_occupation .eq. 2) then
     sizeepsilon = Nks   !2*Nks
  else if (band_occupation .eq. 1) then 
     sizeepsilon = Nks
  else   
      write(*,*) '  ERROR: Unphysical number of electrons per band'
      call exit(1)
  end if

 end subroutine read_QE_occ_unocc_numbers
  

!-------------------------------------------------------------------------------------------    
  
 subroutine read_QE_output (dimband, Ne, number_of_atoms,  number_of_species, band_occupation, &
           &  elec_eigval_file_name, epsilonorig, units_epsilonorig, mass)         

 ! This subroutine reads the file which contains the output of Quantum Espresso and gives electronic eigenvalues and atomic masses.
 
 ! Comment: Ponce2014 considers hbar = m_e = e = 1 (atomic units). This means that, since we use the same convention,
 ! we must multiply the masses expressed with the unified atomic mass pattern by 1822.888486192
 
 

   ! System parameters
  integer, Intent(In) :: dimband                        ! Number of considered electronic eigenvalues
  integer, Intent(In) :: Ne                             ! The number of electrons of the system (i.e. the number of occupied states)
  integer, Intent(In) :: number_of_atoms                ! Total number of atoms of the system
  integer, Intent(In) :: number_of_species              ! Total number of different atomic species of the system
  integer, Intent(In) :: band_occupation                ! Number of electrons (1 or 2) per band in the file output of QE
  character(100), Intent(In) :: elec_eigval_file_name   ! Name of the file where the electronic eigenvalues are stored
  real(8), Intent(Out) :: epsilonorig(dimband)				! Electronic eigenvalues
  character(40), Intent(Out) :: units_epsilonorig           ! Units of the electronic eigenvalues
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
  
        
   write(*,*) '  *** Now reading the output of Quantum Espresso ('&
              & , trim(elec_eigval_file_name),') to find'
   write(*,*)    '                electronic eigenvalues and atomic masses. ***'                    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																									!
!                               1st BLOCK: READING EIGENVALUES										!
!																									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  epsilonorig(:)=0.d0

  open(333,file=elec_eigval_file_name,status='unknown')

  do k=1,20000

      read(333,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
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
 
  read(333,'(A)',IOSTAT=stat) inputline; read(333,'(A)',IOSTAT=stat) inputline 
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
          units_epsilonorig=''
          do j=i+1,i+7
              letter=inputline(j:j)
              if ( (letter .ne. '(') .and.(letter .ne. ')') .and.(letter .ne. ',') .and. &
                & (letter .ne. '.') .and.(letter .ne. ';') .and.(letter .ne. ':')  ) then
                 units_epsilonorig = trim(units_epsilonorig)//letter
              end if
          end do
          if ( ( trim(units_epsilonorig) .eq. 'CM**-1' ) .or. ( trim(units_epsilonorig) .eq. 'Cm**-1' ) ) units_epsilonorig='cm**-1' 
          if ( ( trim(units_epsilonorig) .eq. 'mev' )    .or. ( trim(units_epsilonorig) .eq. 'MeV' ) )    units_epsilonorig='meV'
          if ( ( trim(units_epsilonorig) .eq. 'ev' )  .or. ( trim(units_epsilonorig) .eq. 'eV' ) &
          .or. ( trim(units_epsilonorig) .eq. 'Ev' )     .or. ( trim(units_epsilonorig) .eq. 'EV' ) )     units_epsilonorig='eV'
          if ( ( trim(units_epsilonorig) .eq. 'A.u.' ) .or. ( trim(units_epsilonorig) .eq. 'A.u.' ) .or. &
           &  ( trim(units_epsilonorig) .eq. 'Au.' )  .or. ( trim(units_epsilonorig) .eq. 'AU' ) .or.   &
            & ( trim(units_epsilonorig) .eq. 'hartree' ) .or. ( trim(units_epsilonorig) .eq. 'HARTREE' )  &
            & .or.  ( trim(units_epsilonorig) .eq. 'Hartree' ) )                                          units_epsilonorig='a.u.'

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
        ! write(*,*) ' auxstr2 = ', auxstr2
        read(auxstr2,*) auxr1
        
!        if (band_occupation .eq. 2) then 
           !epsilonorig(i) = auxr1; epsilonorig(i+1) = auxr1;  i = i+2   
!        else 
           epsilonorig(i) = auxr1
           i = i+1         
!        end if  
        
        auxstr2=''
     end if
   end do ! j  
  end do   ! k
     
        
3456 continue

  close(333)


  ! Now we make the 0 HOMO; we consider this reference for all the eigenvalues, and offset them all
  if ( mod(Ne,2) .eq. 1) then
     HOMO = epsilonorig(Ne) 
  else
     HOMO = epsilonorig(Ne/2)
  end if

  do i=1,dimband
    epsilonorig(i) = epsilonorig(i) - HOMO
  end do 
  
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																									!
!                               2nd BLOCK: READING ATOMIC MASSES									!
!																									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(species_name(number_of_species),species_mass(number_of_species))
  
  mass(:,2) = -1.0d0
  
 open(333,file=elec_eigval_file_name,status='unknown')

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


 end subroutine read_QE_output
 

!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------

       
subroutine read_dyneq_output (number_of_atoms, dimnu, n_mode_beg, n_mode_end, phon_freq_file_name, mass, omegaorig, xi, U)     
 
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
  integer, Intent(In) :: n_mode_beg, n_mode_end
  character(100), Intent(In) :: phon_freq_file_name     ! Name of the file where the phonon frequencies and eigenvectors of the dynamical matrix are stored
  real(8), Intent(In)  :: mass(number_of_atoms,2)       ! Masses of the atoms (in the 2nd entry) 
  real(8), Intent(Out) :: omegaorig(dimnu)              ! Phonon frequencies (as they are read from file, with the units provided by the original file)
  real(8), Intent(Out) :: xi(dimnu,dimnu)				! Eigenvectors of the dynamical matrix (xi)
  real(8), Intent(Out) :: U(dimnu,dimnu)				! Eigenvectors of the dynamical matrix (U)
 
 
   ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist
  integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10 
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, auxr10, P_au2A
  real(8), allocatable :: Uau(:,:),sqmass(:)


 
  write(*,*) '  *** Now reading the output of QE (in CPMD format, '&
              & , trim(phon_freq_file_name),' file) to find ' 
  write(*,*)  '             phonon frequencies and eigenvectors of the dynamical matrix. ***'
     


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!     BLOCK FOR READING PHONON FREQUENCIES AND EIGENVECTORS OF THE DYNAMICAL MATRIX        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  allocate(Uau(dimnu,dimnu),sqmass(dimnu))
  omegaorig(:) = 0.d0
  U(:,:)=0.d0
  xi(:,:) =0.d0
  
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
 !   read(341,'(8(4X,G8.3))') auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8
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
  
  
  
  do i=1,number_of_atoms
      sqmass(3*(i-1)+1)=sqrt(mass(i,2)); sqmass(3*(i-1)+2)=sqrt(mass(i,2)); sqmass(3*(i-1)+3)=sqrt(mass(i,2))
  end do
  

  
  !! We re-express the U's in atomic units to check their orthonormality (the U's must be orthonormal in atomic units, but in Uvec.dat they are stored in Angstroms)
   P_au2A = 0.5291772d0
   do i1=1,dimnu
     do i2=1,dimnu
      ! UUUU: En este caso concreto (Si35H36) quitamos el factor P_au2A pq tal factor estaba mal
      Uau(i1,i2) = U(i1,i2)! UUUU: U(i1,i2)/P_au2A       ! PENG:  vib_mode(ii,:) = vib_mode(ii,:)* P_au2A
      U(i1,i2) = Uau(i1,i2)    ! We rescale the U vector to atomic units, because we need everything in atomic units to calculate gDW, and U is an input parameter in the calculate_gDW subroutine.
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

  auxr4=0.d0  ! summation of the orthogonal terms
  auxr5=0.d0  ! summation of terms equaling 1
  auxr6=0.d0   ! maximum error
  
  do j=1,dimnu
    do k=1,dimnu
       
       auxr1=0.d0
       do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
          auxr1=auxr1 + Uau(i,j)*Uau(i,k)*sqmass(i)*sqmass(i)
       end do
       if (j.ne.k) auxr4=auxr4+abs(auxr1)
       if (j.eq.k) then
         ! write(*,*) j, auxr1
         auxr5=auxr5+abs(auxr1-1.d0)
       end if  
       if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
    
    end do
  end do  
  
  auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+6)); auxr5=auxr5/(dble(dimnu-n_mode_beg+1))
  if ( ( abs(auxr4) .lt. 0.00001d0) .and. (abs(auxr5).lt.0.001) ) then 
     write(*,*);  write(*,*) '    - Orthonormality of U (in atomic units, with mass scaling) checked. '
  else
     write(*,*) ' ERROR in U: The first condition of orthonormality is not properly satisfied' 
     write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
     write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5   
     call exit(1)  
  end if

  
  
  ! We have read the U's; now we get xi:
   do j=1,dimnu
     do i=1,dimnu 
       xi(i,j) = Uau(i,j) * sqmass(i)
     end do
   end do  


   ! We normalize the xi vectors (the fact that the output of QE gives normalized U does not mean that the corresponding xi are normalized)
   do j=1,dimnu 
   
     auxr1=0.d0
     do i=1,dimnu
       auxr1 = auxr1 + (xi(i,j))**2 
     end do
     auxr1 = sqrt(auxr1)
     if (abs(auxr1) .lt. 0.000001) then
        write(*,*) 'ERROR: Singularity of xi vector ',j, ', module = ',auxr1
        call exit(1)
     end if
     do i=1,dimnu
       xi(i,j) = xi(i,j)/auxr1 
     end do
     
   end do 
  
 
!  write(*,*) ' CHECKING THE ORTHONORMALITY OF EIGENVECTORS OF THE DYNAMICAL MATRIX'

  auxr4=0.d0  ! summation of the orthogonal terms
  auxr5=0.d0  ! summation of terms equaling 1
  auxr6=0.d0   ! maximum error
  
  do j=1,dimnu
    do k=1,dimnu
       
       auxr1=0.d0
       do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
          auxr1=auxr1 + xi(i,j)*xi(i,k)
       end do
       if (j.ne.k) auxr4=auxr4+abs(auxr1)
       if (j.eq.k) auxr5=auxr5+abs(auxr1-1.d0)
       if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
    
    
    end do
  end do  
  
 
  auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+6)); auxr5=auxr5/(dble(dimnu-n_mode_beg+1))
  if ( ( abs(auxr4) .lt. 0.00001d0) .and. (abs(auxr5).lt.0.00001) ) then 
     write(*,*) '    - Orthonormality of xi in rows checked. '
  else
     write(*,*) '   ERROR: The first condition of orthonormality is not properly satisfied' 
     write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
     write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5   
     call exit(1)  
  end if
     
   
    
     
  auxr4=0.d0  ! summation of the orthogonal terms
  auxr5=0.d0  ! summation of terms equaling 1
  auxr6=0.d0   ! maximum error
  
  do j=7,dimnu
    do k=7,dimnu
       
       auxr1=0.d0
       do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
          auxr1=auxr1 + xi(j,i)*xi(k,i)
       end do
       if (j.ne.k) auxr4=auxr4+abs(auxr1)
       if (j.eq.k) auxr5=auxr5+abs(auxr1)
       if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
    
    end do
  end do  
  auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+6))
  
  auxr5 = auxr5/(dble(dimnu-7+1))
 
  if ( ( abs(auxr4) .lt. 0.00001d0) .and. ( abs(auxr5-1.d0).lt.0.00001) ) then 
     write(*,*) '    - Orthonormality of xi in columns checked. '; write(*,*)
  else
     write(*,*) '   ERROR: The second condition of orthonormality is not properly satisfied'  
     write(*,*) ' The average orthonormaliz 2 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
     write(*,*) ' The average orthonormaliz 2 is (parallel terms):', auxr5   
     call exit(1)  
  end if

  
 deallocate(sqmass,Uau)

end subroutine read_dyneq_output


!-------------------------------------------------------------------------------------------

       
subroutine read_omega_U_CPMD (total_num_atoms, dimnu, FreqsU_file_name, mass, omegaorig, xi, U)     
 
 ! This subroutine reads the file which contains the output of QE and gives phonon frequencies and
 ! eigenvectors (xi and U) of the dynamical equation. We assume that the negligible acoustic modes 
 ! correspond to the lowest six freqs.
 ! Note that the output of CPMD gives the xi eigenvectors, while the output of Quantum Espresso
 ! gives the U eigenvectors. The original code written by Peng used the U eigenvectors; it read the
 ! xi eigenvectors of CPMD and then it divided them by sqrt(M_I). We don't need to do it, because
 ! we know the U eigenvectors directly from the output of Quantum Espresso.
 

  ! System parameters
  integer, Intent(In) :: total_num_atoms                ! The number of atoms of the system
  integer, Intent(In) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  character(100), Intent(In) :: FreqsU_file_name        ! Name of the file where the phonon frequencies and eigenvectors of the dynamical matrix are stored
  real(8), Intent(In)  :: mass(dimnu,2)                 ! Atomic masses
  real(8), Intent(Out) :: omegaorig(dimnu)              ! Phonon frequencies (as they are read from file, with the units provided by the original file)
  real(8), Intent(Out) :: xi(dimnu,dimnu)				! Eigenvectors of the dynamical matrix
  real(8), Intent(Out) :: U(dimnu,dimnu)				! Eigenvectors of the dynamical matrix divided by the square root of the corresponding mass
 
 
   ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT, auxchar1, auxchar2
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist, exist_flag
  integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10 
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, auxr10, P_au2A
  real(8), allocatable :: sqmass(:)

 
  write(*,*) '  *** Now reading the output of QE ('&
              & , trim(FreqsU_file_name),') to find ' 
  write(*,*)  '             phonon frequencies and (U) eigenvectors of the dynamical matrix. ***'
     


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!     BLOCK FOR READING PHONON FREQUENCIES AND EIGENVECTORS OF THE DYNAMICAL MATRIX        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
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
        if ( auxstr1 .eq. 'freq') then  
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
      read(341,*) auxchar1, auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxchar2
      U((i-1)*3+1,freqcount) = auxr1
      U((i-1)*3+2,freqcount) = auxr3
      U((i-1)*3+3,freqcount) = auxr5
      ! write(*,*) auxr1, auxr2, auxr3, auxr4, auxr5, auxr6
    end do    

 end do ! m      
  

!   Input for Maat:  * wQ0.dat: 1st row: number of q-points and the number of phonon frequencies; then: q-index, branch-index, phonon frequency  
   open(343,file='wQ0.dat',status='unknown')
   write(343,*) 1, dimnu
   do m=1,dimnu
      write(343,*) 1,m,omegaorig(m)
   end do
   close(343)
   
   
   ! File defined by Peng
   open(90, file='frequency.dat')
   do ii = 1, 3*total_num_atoms
       write(90,*) 'I=', ii, 'Freq=', omegaorig(ii)
   end do
   close(90)   
   


  allocate(sqmass(3*number_of_atoms))
  
  do i=1,number_of_atoms
      sqmass(3*(i-1)+1)=sqrt(mass(i,2)); sqmass(3*(i-1)+2)=sqrt(mass(i,2)); sqmass(3*(i-1)+3)=sqrt(mass(i,2))
  end do
  
    do j=1,dimnu
     do i=1,dimnu 
          xi(i,j) = U(i,j) * sqmass(i)
     end do
   end do    


   ! We normalize the xi vectors (the fact that the output of QE gives normalized U does not mean that the corresponding xi are normalized)
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
      
 
 
 
 
!  write(*,*) ' CHECKING THE ORTHONORMALITY OF EIGENVECTORS OF THE DYNAMICAL MATRIX'

  auxr4=0.d0  ! summation of the orthogonal terms
  auxr5=0.d0  ! summation of terms equaling 1
  auxr6=0.d0   ! maximum error
  
  do j=7,dimnu
    do k=7,dimnu
       
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
 

  
  auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+6)); auxr5=auxr5/(dble(dimnu-n_mode_beg+1))
  if ( ( abs(auxr4) .lt. 0.00001d0) .and. (abs(auxr5-1.d0).lt.0.00001) ) then 
      write(*,*) ;  write(*,*) '   ORTHONORMALITY  1: OK '
!write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
!write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5 
  else
     write(*,*) ' <<<<<<<<<<<< WARNING: The first condition of orthonormality is not properly satisfied' 
     write(*,*) ' The average orthonormaliz 1 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
     write(*,*) ' The average orthonormaliz 1 is (parallel terms):', auxr5, ' >>>>>>>>>>>'     
  end if


  auxr4=0.d0  ! summation of the orthogonal terms
  auxr5=0.d0  ! summation of terms equaling 1
  auxr6=0.d0   ! maximum error
  
  do j=7,dimnu
    do k=7,dimnu
       
       auxr1=0.d0
       do i=1,dimnu     ! OJO, NU DEBES TRUNCAR ESTE SUMATORIO, DEBE INCLUIR TODOS LOS SUMANDOS!!!
          auxr1=auxr1 + xi(j,i)*xi(k,i)
       end do
       if (j.ne.k) auxr4=auxr4+abs(auxr1)
       if (j.eq.k) auxr5=auxr5+abs(auxr1)
       if ((abs(auxr1).gt.auxr6).and.(j.ne.k) ) auxr6=abs(auxr1)
    
    end do
  end do  
  auxr4=auxr4/(dble((dimnu-n_mode_beg+1)**2 - dimnu+6)); auxr5=auxr5/(dble(dimnu-n_mode_beg+1))
 
  if ( ( abs(auxr4) .lt. 0.00001d0) .and. ( abs(auxr5-1.d0).lt.0.00001) ) then 
    write(*,*) '   ORTHONORMALITY  2: OK '
!write(*,*) ' The average orthonormaliz 2 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
!write(*,*) ' The average orthonormaliz 2 is (parallel terms):', auxr5             
  else
     write(*,*) ' <<<<<<<<<<<< WARNING: The second condition of orthonormality is not properly satisfied'  
     write(*,*) ' The average orthonormaliz 2 is (orthogonal terms):', auxr4, '(max val=',auxr6,')'
     write(*,*) ' The average orthonormaliz 2 is (parallel terms):', auxr5 , ' >>>>>>>>>>>'  
  end if


!! Transfer atomic units to Angstrom (as done by Peng); this is necessary because the output of QE is expected to be given in Angstrom
   P_au2A = 0.5291772d0
   do ii = 1, 3*total_num_atoms
      U(ii,:) = U(ii,:)* P_au2A       ! PENG:  vib_mode(ii,:) = vib_mode(ii,:)* P_au2A
   end do



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

  deallocate(sqmass)

end subroutine read_omega_U_CPMD






!-------------------------------------------------------------------------------------------

       
subroutine read_mass_CPMD (number_of_atoms, dimnu, Cmass, Hmass, Omass, Nmass, mass_file_name, mass,units_omegaorig)     
 
 ! This subroutine reads the file which contains the output of CPMD and gives the masses and indices of the atoms.  
 ! Comment: Ponce2014 considers hbar = m_e = e = 1. This means that, if we use the same convention,
 ! we must multiply the masses expressed with the unified atomic mass pattern by 1822.888486192
 
  ! System parameters
  integer, Intent(In) :: number_of_atoms                ! The number of atoms of the system
  integer, Intent(In) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  real(8), Intent(In) :: Cmass, Hmass, Omass, Nmass     ! Masses of Carbon, Hydrogen, Oxygen and Nitrogen
  character(100), Intent(In) :: mass_file_name          ! Name of the file where the phonon frequencies and eigenvectors of the dynamical matrix are stored
  real(8), Intent(Out) :: mass (number_of_atoms,2)      ! Variable containing (ordered by atom index) the atomic number and the corresponding mass
  character(40), Intent(Out) :: units_omegaorig             ! The units of the phonon frequencies
 
   ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,k,l1,l2,m,punto,iostat,unitt, freqcount 
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9


 
  write(*,*); write(*,*) '  *** Now reading the output of CPMD ('&
              & , trim(mass_file_name),') to find atomic masses and indices. ***'
    
  
  open(341,file=mass_file_name,status='unknown')
  
  do k=1,(dimnu+20)**2

      read(341,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 2117
      
      auxstr1=''; auxstr2=''
      
      do i=1,50
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
        if ( auxstr1 .eq. 'Harmonicfrequencies' )  then
             units_omegaorig=''
             do j=1,10
                letter=inputline(i+j:i+j)
                if ( (letter .ne. '(') .and.(letter .ne. ')') .and.(letter .ne. ',') .and. &
                  & (letter .ne. '.') .and.(letter .ne. ';') .and.(letter .ne. ':')  ) then
                   units_omegaorig = trim(units_omegaorig)//letter
                end if
             end do
             if (units_omegaorig .eq. 'ev') units_omegaorig='eV'
             if (units_omegaorig .eq. 'mev') units_omegaorig='meV'
             if (units_omegaorig .eq. 'Cm**-1') units_omegaorig='cm**-1'
          go to 4235
        end if
      end do ! i
      
   end do !k   
   
 4235 continue   

  close(341)

!   read(333,'(A)',IOSTAT=stat) inputline; read(333,'(A)',IOSTAT=stat) inputline 
!   auxstr1=''
!   do i=1,60
!      letter = inputline(i:i)
!      auxstr1 = trim(auxstr1)//letter 
!      auxstr2=''
!    !  write(*,*) i,letter, auxstr1, LNBLNK(auxstr1)
!      j=LNBLNK(auxstr1)   ! Last index of the chain that is not a blank space
!      if ( j.ge.19) then
!        do k=1,19
!           letter=inputline(i-19+k:i-19+k)
!           auxstr2 = trim(auxstr2)//letter
!        end do
!      end if
!      
!      if ( auxstr2 .eq. 'Harmonicfrequencies') then 
!           auxstr3=''
!           do j=i+1,i+7
!               letter=inputline(j:j)
!               auxstr3 = trim(auxstr3)//letter
!           end do
!           write(*,*) '         The units for phonon frequencies are ',auxstr3; write(*,*)
!      end if
!   end do ! i        

  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!     BLOCK FOR READING NUCLEAR INDICES AND MASSES       !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  mass(:,:) = 0.d0
  open(341,file=mass_file_name,status='unknown')
  
  do k=1,15!(dimnu+20)**2

      read(341,'(A)',IOSTAT=stat) inputline
      ! write(*,*) inputline
      if (stat .ne. 0) go to 2117

      auxstr1=''; auxstr2=''
      
      do i=1,50
        letter = inputline(i:i)
        auxstr1 = trim(auxstr1)//letter 
     
        if ( auxstr1 .eq. 'NumberNumberTypeX' )   go to 2235
      end do ! i
      
   end do !k   
   
 2235 continue 

  
 read(341,'(A)',IOSTAT=stat) inputline
 
 do l1=1,number_of_atoms
   
    read(341,*) i, j, k, auxr1, auxr2, auxr3
    mass(i,1) = dble(j)
    
    if ( (i.lt.1) .or. (i.gt.number_of_atoms) ) then
       write(*,'(A,I3,A)') '  ERROR: Wrong atomic index (',i,') '; write(*,*)
       call exit(1)
    end if
    
    if      (j .eq. 1) then
       mass(i,2) = Hmass * 1822.888486192d0
    else if (j .eq. 6) then
       mass(i,2) = Cmass * 1822.888486192d0  
    else if (j .eq. 7) then
       mass(i,2) = Nmass * 1822.888486192d0 
    else if (j .eq. 8) then
       mass(i,2) = Omass * 1822.888486192d0         
    else
      write(*,*) '  ERROR: Mass not available for atomic number Z=',j,' (index=',i,') '; write(*,*)
      call exit(1)
    end if
    
  end do ! l1
   
  
 2117 continue  
 

  
!  do i=1,number_of_atoms
!    write(*,*) i, mass(i,1), mass(i,2)
!  end do



! 
! 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !!     BLOCK FOR READING PHONON FREQUENCIES AND EIGENVECTORS OF THE DYNAMICAL MATRIX        !!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
!   
!   omega(:) = 0.d0
!   xi(:,:) = 0.d0
! 
!   
!  freqcount = -2
!   
!  do m=1,dimnu+1
!   
!   do k=1,(dimnu+20)**2
! 
!       read(341,'(A)',IOSTAT=stat) inputline
!       ! write(*,*) inputline
!       if (stat .ne. 0) go to 1117
! 
!       auxstr1=''; auxstr2=''; auxstr3=''
!       
!       do i=1,90
!         letter = inputline(i:i)
!         auxstr1 = trim(auxstr1)//letter 
!         if ( auxstr1 .eq. 'Frequencies--') then  
!           freqcount = freqcount + 3
!           !write(*,*) 'CONSEGUIDO, linea', k, i
!           l1 = k       ! l1 is the line where the word 'Frequencies' is found
!           go to 1235
!         end if
!       end do ! i
!       
!    end do !k   
!    
!  1235 continue      
!       
!    punto=100
!    j=i+1
!    do while ( (j .le. 90) .and. (j.le.punto+6))
!       ! write(*,*) j, punto
!       letter = inputline(j:j)
!       if (letter .eq. '.') punto=j
!       auxstr2 = trim(auxstr2)//letter
!       j=j+1
!    end do  
!    ! write(*,*) '-->',trim(auxstr2)
!    read(auxstr2,*) auxr1
!    omega(freqcount)=auxr1
!    auxstr2=''
! 
!  
!  
!    punto=100
!    j=j+1
!    do while ( (j .le. 90) .and. (j.le.punto+6))
!       ! write(*,*) j, punto
!       letter = inputline(j:j)
!       if (letter .eq. '.') punto=j
!       auxstr2 = trim(auxstr2)//letter
!       j=j+1
!    end do  
!    !write(*,*) '-->',trim(auxstr2)
!    read(auxstr2,*) auxr1
!    omega(freqcount+1)=auxr1
!    auxstr2=''
!       
!    punto=100
!    j=j+1
!    do while ( (j .le. 90) .and. (j.le.punto+6))
!       ! write(*,*) j, punto
!       letter = inputline(j:j)
!       if (letter .eq. '.') punto=j
!       auxstr2 = trim(auxstr2)//letter
!       j=j+1
!    end do  
!    !write(*,*) '-->',trim(auxstr2)
!       read(auxstr2,*) auxr1
!    omega(freqcount+2)=auxr1  
!    auxstr2=''       
!    
!    
!    
! !..........................................................................
!   ! The block below reads three consecutive eigenvectors of the dyn. matr. We cannot move this block to a different subroutine because we need to read file 341
! 
!   do k=l1+1,(dimnu+20)**2
! 
!       read(341,'(A)',IOSTAT=stat) inputline
!       ! write(*,*) inputline
!       if (stat .ne. 0) go to 1118
! 
!       auxstr1=''; auxstr2=''; auxstr3=''
!       
!       do i=1,90
!         letter = inputline(i:i)
!         auxstr1 = trim(auxstr1)//letter 
!         if ( auxstr1 .eq. 'Atom') then  
!          ! write(*,*) 'CONSEGUIDO, linea', k, i
!           l2 = k       ! l is the line where the word 'Frequencies' is found
!           go to 1236
!         end if
!       end do ! i
!       
!    end do !k  
!  
!   1236 continue
!         
!   read(341,IOSTAT=stat) inputline
!       
!   do k=1,number_of_atoms   
!     read(341,*) i, j, auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9
!     ! write(*,*) i, j, auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9
!     xi((3*(k-1)+1),freqcount)  =auxr1;  xi((3*(k-1)+2),freqcount)  =auxr2;  xi((3*(k-1)+3),freqcount)  =auxr3
!     xi((3*(k-1)+1),freqcount+1)=auxr4;  xi((3*(k-1)+2),freqcount+1)=auxr5;  xi((3*(k-1)+3),freqcount+1)=auxr6
!     xi((3*(k-1)+1),freqcount+2)=auxr7;  xi((3*(k-1)+2),freqcount+2)=auxr8;  xi((3*(k-1)+3),freqcount+2)=auxr9
!   end do
!  
! 1118 continue 
! 
! !..........................................................................   
!    
!    
! end do ! m      
! 
!    
! 1117 continue
! 
! close(341)
 
!        write(*,*) 'xi(',1,'-',3,')=',xi(1,3)
!        write(*,*) 'xi(',5,'-',3,')=',xi(5,3)
!  write(*,*) 'xi(',9,'-',3,')=',xi(9,3)
!  
!        write(*,*) 'xi(',1,'-',6,')=',xi(1,6)
!        write(*,*) 'xi(',5,'-',6,')=',xi(5,6)
!  write(*,*) 'xi(',9,'-',6,')=',xi(9,6)
! 
! 
!        write(*,*) 'xi(',1,'-',75,')=',xi(1,75)
!        write(*,*) 'xi(',5,'-',75,')=',xi(5,75)
!        write(*,*) 'xi(',9,'-',75,')=',xi(9,75)
!
!do i=1,dimnu
! write(*,*) i, omega(i)
!end do

 
!  write(*,*) ' CHECK ORTHONORMALITY'
!  j=10; k=21
!  auxr1=0.d0
!  do i=1,78
!    auxr1=auxr1 + xi(j,i)*xi(k,i)
!  end do
!  write(*,*) auxr1
!  call exit(1)

end subroutine read_mass_CPMD


!-------------------------------------------------------------------------------------------



subroutine omega_epsilon_to_atomic_units (units_omegaorig,units_epsilonorig,dimnu,dimband, &        
                                             & omegaorig,epsilonorig,omega,epsilon)
     
  ! We use atomic units as done in Ponce et al. PRB 90, 214304 (2014); this is because our derivation of gDW
  ! depends on the definitions given in that paper.   
  ! This subroutine writes the phonon freqs. in the units of the electronic eigenvectors, multiplying the former
  ! values of the omegas by the 'factor'.   
       
  ! System parameters
  character(40), Intent(In) :: units_omegaorig, units_epsilonorig
  integer, Intent(In) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
  integer, Intent(In) :: dimband                        ! Number of considered electronic eigenvalues
  real(8), Intent(In) :: omegaorig(dimnu)               ! Original omegas, in the original units
  real(8), Intent(In) :: epsilonorig(dimband)			! Electronic eigenvalues in their original units
  real(8), Intent(Out):: omega(dimnu)                   ! Omegas in atomic units
  real(8), Intent(Out):: epsilon(dimband)				! Electronic eigenvalues in atomic units
  
 
  ! Local variable
  integer :: i
  real(8) :: factor  

  factor = 1.d0
  
 ! write(*,*) '  *** Now converting the units of the frequency (',trim(units_omegaorig),' to atomic units). ***'
! 
!   if    ( (trim(units_omegaorig).eq.'cm**-1') .or. (trim(units_omegaorig).eq.'cm-1') ) then
!     factor=1.d0/219474.631370515d0 ! 0.00000455633d0
!   else if (trim(units_omegaorig).eq.'meV') then  
!     factor=1.d0/27211.385056d0  ! 0.0000367502d0   
!   else if (trim(units_omegaorig).eq.'eV') then  
!     factor=1.d0/27.211385056d0 ! 0.0367502d0
!   else if (trim(units_omegaorig).eq.'a.u.') then  
!     factor=1.d0           
!   else
!     write(*,*) '  ERROR: Conversion of units not established; please, use "eV", "meV", "cm**-1", or "a.u.", '
!     write(*,*) '  or rewrite the <<omega_epsilon_to_atomic_units>> subroutine'
!     call exit(1) 
!   end if

   do i=1,dimnu
     omega(i) = (omegaorig(i))*factor
   end do
   
  do i=1,6
    omega(i) = 1.d0 !! We impose the nonphysical frequencies to be 1 (in a.u.) to avoid singularities when calculating the h factors
  end do
    

  
! We write the phonon frequencies with a format that Maat can read:
! wQ0.dat: 1st row: number of q-points and the number of phonon frequencies; then: q-index, branch-index, phonon frequency
  open(533,file='wQ0.dat',status='unknown')
  write(533,*)  dimnu 
  do i=1,dimnu
    write(533,*)  i,' ', omega(i)
  end do
  close(533)
  
 
  factor = 1.d0
 
 
!    write(*,*) '  *** Now converting the units of the the electronic eigenvalues (',trim(units_epsilonorig)
!    write(*,*) '                                  to atomic units). ***'
! 
!   if      (trim(units_epsilonorig).eq.'cm**-1') then
!     factor=1.d0/219474.631370515d0 ! 0.00000455633d0
!   else if (trim(units_epsilonorig).eq.'meV') then  
!     factor=1.d0/27211.385056d0  ! 0.0000367502d0 
!   else if (trim(units_epsilonorig).eq.'eV') then  
!     factor=1.d0/27.211385056d0 ! 0.0367502d0
!   else if (trim(units_epsilonorig).eq.'a.u.') then  
!     factor=1.d0           
!   else
!     write(*,*) '  ERROR: Conversion of units not established; please, use "eV", "meV", "cm**-1", or "a.u.", '
!     write(*,*) '  or rewrite the <<omega_epsilon_to_atomic_units>> subroutine'
!     call exit(1)    
!   end if


   do i=1,dimband
     epsilon(i) = (epsilonorig(i))*factor
   end do
  
! We write the phonon frequencies with a format that Maat can read:
  
! We write the electronic eigenvalues with a format that Maat can read:
! epsilonorig0.dat: 1st row: number of bands, number of k-points, number of spins; then all three indices and the corresponding unrenormalized electronic eigenvalue
  open(633,file='epsilon0.dat',status='unknown')
  write(633,*) dimband
  do i=1,dimband
    write(633,*) i,  epsilon(i)
  end do
  close(633)
  
        
       
end subroutine omega_epsilon_to_atomic_units   


!-------------------------------------------------------------------------------------------


! 
! 
! subroutine change_units_omega (units_omegaorig,units_epsilonorig,dimnu,omegaorig,omega)
!      
!   ! Since the units of the electronic eigenvectors are sometimes different to those of the phonon frequencies,
!   ! this subroutine writes the phonon freqs. in the units of the electronic eigenvectors, multiplying the former
!   ! values of the omegas by the 'factor'.   
!        
!   ! System parameters
!   character(40), Intent(In) :: units_omegaorig, units_epsilonorig
!   integer, Intent(In) :: dimnu                          ! Number of phonon frequencies (i.e. dimension of the dynamical matrix, i.e. 3*N_at-6)
!   real(8), Intent(In) :: omegaorig(dimnu)               ! Original omegas, in the original units
!   real(8), Intent(Out) :: omega(dimnu)                   ! Omegas in the new units (those of epsilonorig)
!   
!  
!   ! Local variable
!   integer :: i
!   real(8) :: factor  
! 
! ! 
! !  if (( trim(units_omegaorig) .ne. 'cm**-1' ) .and. ( trim(units_omegaorig) .ne. 'CM**-1' ).and. &
! !   & ( trim(units_omegaorig) .ne. 'Cm**-1' ) &
! !    &  .and. ( trim(units_omegaorig) .ne. 'mev' ) .and. ( trim(units_omegaorig) .ne. 'meV' ) &
! !     .and. ( trim(units_omegaorig) .ne. 'a.u.' ).and. ( trim(units_omegaorig) .ne. 'A.U.' ) &
! !     .and. ( trim(units_omegaorig) .ne. 'eV' ).and. ( trim(units_omegaorig) .ne. 'ev' ) &
! !     &.and. ( trim(units_omegaorig) .ne. 'Hartree' ).and. ( trim(units_omegaorig) .ne. 'hartree' ) ) then
! 
!   write(*,*) '  *** Now converting the units of the frequency (',trim(units_omegaorig),' to ',trim(units_epsilonorig), '). ***'
! !  write(*,*) '         The units for phonon frequencies are ',units_omegaorig
! !                write(*,*) '         The units for electronic eigenvalues are ',units_epsilonorig; write(*,*)
! 
!   if      ((units_epsilonorig.eq.'eV').and.(units_omegaorig.eq.'cm**-1')) then
!     factor=0.000123984193d0
!   else if ((units_epsilonorig.eq.'eV').and.(units_omegaorig.eq.'meV')) then  
!     factor=0.001d0 
!   else if ((units_epsilonorig.eq.'eV').and.(units_omegaorig.eq.'a.u.')) then  
!     factor=27.2113838668d0
!   else if ((units_epsilonorig.eq.'eV').and.(units_omegaorig.eq.'eV')) then  
!     factor=1.d0   
! !    
!   else if ((units_epsilonorig.eq.'meV').and.(units_omegaorig.eq.'cm**-1')) then
!     factor=0.123984193d0
!   else if ((units_epsilonorig.eq.'meV').and.(units_omegaorig.eq.'eV')) then  
!     factor=1000.d0   
!   else if ((units_epsilonorig.eq.'meV').and.(units_omegaorig.eq.'a.u.')) then  
!     factor=27211.3838668d0
!   else if ((units_epsilonorig.eq.'meV').and.(units_omegaorig.eq.'meV')) then  
!     factor=1.d0
! !  
!   else if ((units_epsilonorig.eq.'cm**-1').and.(units_omegaorig.eq.'eV')) then
!     factor=1.d0/0.000123984193d0
!   else if ((units_epsilonorig.eq.'cm**-1').and.(units_omegaorig.eq.'meV')) then  
!     factor=1.d0/0.123984193d0   
!   else if ((units_epsilonorig.eq.'cm**-1').and.(units_omegaorig.eq.'a.u.')) then  
!     factor=1.d0/0.00000455633d0
!   else if ((units_epsilonorig.eq.'cm**-1').and.(units_omegaorig.eq.'cm**-1')) then  
!     factor=1.d0
! !  
!   else if ((units_epsilonorig.eq.'a.u.').and.(units_omegaorig.eq.'cm**-1')) then
!     factor=0.00000455633d0
!   else if ((units_epsilonorig.eq.'a.u.').and.(units_omegaorig.eq.'meV')) then  
!     factor=0.0000367502d0   
!   else if ((units_epsilonorig.eq.'a.u.').and.(units_omegaorig.eq.'eV')) then  
!     factor=0.0367502d0
!   else if ((units_epsilonorig.eq.'a.u.').and.(units_omegaorig.eq.'a.u.')) then  
!     factor=1.d0    
! !         
!   else
!     write(*,*) '  ERROR: Conversion of units not established; please, rewrite the <<change_units_omega>> subroutine'
!     call exit(1) 
!   end if
! 
! 
! 
!    do i=1,dimnu
!      omega(i) = (omegaorig(i))*factor
!    end do
!    
!    
! 
!   
! ! We write the phonon frequencies with a format that Maat can read:
! ! wQ0.dat: 1st row: number of q-points and the number of phonon frequencies; then: q-index, branch-index, phonon frequency
!   open(533,file='wQ0.dat',status='unknown')
!   write(533,*) '  1  ', dimnu 
!   do i=1,dimnu
!     write(533,*) '  1  ', i,' ', omega(i)
!   end do
!   close(533)
!   
!        
!        
! end subroutine change_units_omega    
! 
! 
! !-------------------------------------------------------------------------------------------


 subroutine read_gFan(sizeepsilon, epsilon, dimnu, omega, gFan_i_index_beg,gFan_i_index_end,&  
                    & gFan_j_index_beg,gFan_j_index_end, units_epsilonorig,gFan_file_name,gFan)
 
 ! This subroutine reads gFan from file and re-expresses it to atomic units; it takes into account the degeneracy of electronic eigenvalues and phonon frequencies,
 ! making average of |g^{ij}_{\nu}|^2 through the degenerate eigenvalues (j) and phonon branches (\nu)
 ! This is valid just for isolated systems. If you want to extend it to be valid to periodic systems you have to change **2 by ||**2, and make epsilonorig dependent on k.
 
  ! System parameters
  
  Integer, Intent(In) :: sizeepsilon
  real(8), Intent(In) :: epsilon(sizeepsilon) ! Initial electronic eigenvectors. They are necessary because, depending on their degeneracy, one must average the corresponding gFan's
  Integer, Intent(In) :: dimnu
  real(8), Intent(In) :: omega(dimnu)             ! Initial phononic frequencies. They are necessary because, depending on their degeneracy, one must average the corresponding gFan's.
  Integer, Intent(In) :: gFan_i_index_beg         ! Initial wavefunction index for < psi_i | for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(In) :: gFan_i_index_end         ! Final wavefunction index < psi_i | for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(In) :: gFan_j_index_beg         ! Initial wavefunction index for | psi_j >for the calculation of (g^{DW})^{ii}_{nu}
  Integer, Intent(In) :: gFan_j_index_end         ! Final wavefunction index  | psi_j >  for the calculation of (g^{DW})^{ii}_{nu}    
  character(40), Intent(In)  :: units_epsilonorig ! The units of the electronic eigenvalues given by Quantum Espresso, which are also the units of the Fan matrix elements by construction; we will reexpress gFan in atomic units
  character(100), Intent(In) :: gFan_file_name    ! Name of the file where the original gFan (generated by elph_pwscf, and not necessarily in a.u.) lie. Not to confuse with the file with gFan in a.u. (gFan.dat).
  complex(8), Intent(Out) :: gFan(gFan_i_index_end-gFan_i_index_beg+1,gFan_j_index_end-gFan_j_index_beg+1,1,1,1,dimnu)
  
  ! Local variables
  character(100) :: inputline, auxstr1, auxstr2, auxstr3
  character(20)  :: Tunits, stringT
  character(25) :: auxstrnum1,auxstrnum2,auxstrnum3,auxstrnum4,auxstrnum5,auxstrnum6,auxstrnum7,auxstrnum8,auxstrnum9
  character :: letter, letter2, letter3, letter4
  integer :: stat,auxi,i,j,nu,k,l1,l2,m,punto,iostat,unitt,freqcount,i1,i2,i3,i4,i5,i6,j1,j2,j3,k1,k2,k3,k4,k5,k6 
  integer :: gFan_nu_index_beg,gFan_nu_index_end, firstj, lastj, dif
  logical :: exist
  real(8) :: auxr1, auxr2, auxr3, auxr4, auxr5, auxr6, auxr7, auxr8, auxr9, factor
  complex(8), allocatable :: gFanaux1(:,:,:,:,:,:), gFanaux2(:,:,:,:,:,:)
  
  allocate(gFanaux1(gFan_i_index_end-gFan_i_index_beg+1,gFan_j_index_end-gFan_j_index_beg+1,1,1,1,dimnu))
  allocate(gFanaux2(gFan_i_index_end-gFan_i_index_beg+1,gFan_j_index_end-gFan_j_index_beg+1,1,1,1,dimnu))

  
  ! gFanaux1 is the value read from file; gFanaux2 is gFanaux1 averaged in degenerate electronic eigenvalues;
  ! gFan is gFanaux2 averaged in degenerate phonon frequencies.
  
  !  
!          write(auxstr1,*)  iband1+wfc_i_index_beg-1; call StripSpaces(auxstr1)
!          write(auxstr2,*)  iband2+wfc_j_index_beg-1; call StripSpaces(auxstr2)
!          do nu=1,n_mode_cal
!            write(auxstr3,*)  nu+n_mode_beg-1;          call StripSpaces(auxstr3)
!            write(155,*) trim(auxstr1),' ',trim(auxstr2),' 1 1 1 ', trim(auxstr3), gFan(iband1,iband2,nu), ' 0.d0'    
!            
 
  write(*,*) '  *** Now reading the Fan electron-phonon matrix elements from the ' 
  write(*,*) '          ',trim(gFan_file_name),' file. We assume atomic units in it.  ***'
  
 
! READING THE gFan CALCULATED BY OUR PROGRAM (elph_pwscf.x) 
  open(241,file=gFan_file_name,status='unknown')

!   
!   gFan(:,:,:,:,:,:)=0.d0
!   read(241,*) k1,k2,k3,k4,k5,k6
!   do i1=1,k1!gFan_i_index_end-gFan_i_index_beg+1
!     do i2=1,k2!gFan_j_index_end-gFan_j_index_beg+1
!       do i3=1,k3
!         do i4=1,k4
!            do i5=1,k5
!              do i6=1,k6
!                read(241,*,IOSTAT=stat) i,j, j1,j2,j3, nu, auxr1, auxr2
!                gFan ( i-gFan_i_index_beg+1, j-gFan_j_index_beg+1, j1,j2,j3, nu-6 ) = auxr1*(1.d0,0.d0) + auxr2*(0.d0,1.d0)
!                !gFan(i-gFan_i_index_beg+1,j-gFan_j_index_beg+1,j1,j2,j3,nu) = auxr1*(1.d0,0.d0) + auxr2*(0.d0,1.d0)
!              end do
!            end do
!          end do
!        end do        
!     end do  
!   end do 
!    

!   if      (units_epsilonorig.eq.'cm**-1') then
!     factor=0.00000455633d0
!   else if (units_epsilonorig.eq.'meV') then  
!     factor=0.0000367502d0   
!   else if (units_epsilonorig.eq.'eV') then  
!     factor=0.0367502d0
!   else if (units_epsilonorig.eq.'a.u.') then  
!     factor=1.d0           
!   else
!     write(*,*) '  ERROR: Conversion of units not established; please, rewrite the <<omega_epsilon_to_atomic_units>> subroutine'
!     call exit(1) 
!   end if



  gFanaux1(:,:,:,:,:,:)=0.d0; gFanaux2(:,:,:,:,:,:)=0.d0; gFan(:,:,:,:,:,:)=0.d0
  read(241,*) k1,i1,k2,i2,k3,k4,k5,gFan_nu_index_beg,gFan_nu_index_end 
  do i1=gFan_i_index_beg,gFan_i_index_end 
    do i2=gFan_j_index_beg,gFan_j_index_end
      do i3=1,k3
        do i4=1,k4
           do i5=1,k5
             do i6=1,gFan_nu_index_end 
               read(241,*,IOSTAT=stat) i,j, j1,j2,j3, nu, auxr1, auxr2
               if (stat .ne. 0) go to 7890
               gFanaux1 ( i-gFan_i_index_beg+1, j-gFan_j_index_beg+1, j1,j2,j3, nu ) = &
                                             &    auxr1*(1.d0,0.d0) + auxr2*(0.d0,1.d0)                             
             end do
           end do
         end do
       end do        
    end do  
  end do 

 7890 continue   
  
  close(241)


! 
!  open(243,file='gFan-av1.dat',status='unknown')
! 
! ! AVERAGE OF gFan FOR DEGENERATE ELECTRONIC STATES AND PHONONIC BRANCHES
! ! Average in electronic degenerate eigenvalues (2nd band index)
!   do i1=1,gFan_i_index_end-gFan_i_index_beg+1 
! 	do i3=1,k3
! 	  do i4=1,k4
! 		do i5=1,k5
! 		  do i6=1,gFan_nu_index_end 
! 		   
! 		     firstj=1; lastj=1
! 		     auxr1 = (gFanaux1(i1,1,i3,i4,i5,i6))**2
! 		     
! 		     do i2=2,gFan_j_index_end-gFan_j_index_beg+1-1 
!              
!                if ( ( abs(epsilon(i2)-epsilon(i2-1)) .gt. 0.00001d0 )  ) then  ! if i2 and i2-1  DO NOT correspond to the same degenerate electronic eigenvalue
!                  auxr1 = ( sqrt(auxr1) ) / (dble(lastj-firstj+1))
!                  do j=firstj,lastj
!                    gFanaux2(i1,j,i3,i4,i5,i6)=dcmplx(auxr1)
!                  end do
!                  auxr1=(gFanaux1(i1,i2,i3,i4,i5,i6))**2
!                  firstj=i2; lastj=i2
!                else                       ! if i2 and i2-1 correspond to the same degenerate electronic eigenvalue
!                  auxr1 = auxr1 + (gFanaux1(i1,i2,i3,i4,i5,i6))**2
!                  lastj = lastj + 1   
!                end if
!               
!              end do ! i2
!              
!              i2=gFan_j_index_end-gFan_j_index_beg+1
! 			 if ( ( abs(epsilon(i2)-epsilon(i2-1)) .gt. 0.00001d0 )  ) then  ! if i2 and i2-1  DO NOT correspond to the same degenerate electronic eigenvalue
!                  auxr1 = ( sqrt(auxr1) ) / (dble(i2-1-firstj+1))
!                  do j=firstj,i2-1
!                    gFanaux2(i1,j,i3,i4,i5,i6) = auxr1
!                  end do
! 				 gFanaux2(i1,i2,i3,i4,i5,i6) = gFanaux1(i1,i2,i3,i4,i5,i6)
! 			 else                       ! if i2 and i2-1 correspond to the same degenerate eigenvalue
! 			   auxr1 = auxr1 + (gFanaux1(i1,i2,i3,i4,i5,i6))**2
! 			   do j=firstj,i2
! 			     gFanaux2(i1,j,i3,i4,i5,i6) = ( sqrt(auxr1) ) / (dble(i2-firstj+1))
! 			   end do  
! 			 end if             
! 			 
! 			 if (i6.eq.1) then 
! 			   do i2=1,gFan_j_index_end-gFan_j_index_beg+1
! 			     write(243,'(I4,A,I4,A,I4,A,f12.8,A,f12.8)')  & 
! 			          & i1,' ',i2,' ',i6,' ',dble(gFanaux2(i1,i2,i3,i4,i5,i6)),' - ',epsilon(i2)
! 			   end do
! 			 end if  
!              
!           end do
!         end do
!       end do        
!     end do  
!   end do 
!   close(243)
! 
! ! 
! ! ! Average in phononic degenerate branches
! !  gFanaux1(:,:,:,:,:,:)=(0.d0,0.d0)
! !  open(244,file='gFan-av2.dat',status='unknown')
! !   do i1=1,gFan_i_index_end-gFan_i_index_beg+1 
! !     do i2=1,gFan_j_index_end-gFan_j_index_beg+1 
! ! 	  do i3=1,k3
! ! 	    do i4=1,k4
! ! 	      do i5=1,k5
! ! 		   
! ! 		     firstj=1; lastj=1
! ! 		     auxr1 = (gFanaux2(i1,i2,i3,i4,i5,1))**2     
! ! 		     
! ! 		     do i6=2,gFan_nu_index_end-1 
! !              
! !                if ( ( abs(omega(i6)-omega(i6-1)) .gt. 0.00000367502d0  )  ) then  ! if i6 and i6-1  DO NOT correspond to the same degenerate freq
! !                  auxr1 = ( sqrt(auxr1) ) / (dble(lastj-firstj+1))
! !                  do j=firstj,lastj
! !                    gFanaux1(i1,i2,i3,i4,i5,j)=dcmplx(auxr1)
! !                  end do
! !                  auxr1=(gFanaux2(i1,i2,i3,i4,i5,i6))**2
! !                  firstj=i6; lastj=i6
! !                !  write(*,*) i6, 'eigval nuevo',epsilon(i6)
! !                else                       ! if i6 and i6-1 correspond to the same degenerate freq
! !                  auxr1 = auxr1 + (gFanaux2(i1,i2,i3,i4,i5,i6))**2
! !                  lastj = lastj + 1   
! !                !  write(*,*) i6, 'eigval repetido',epsilon(i6)
! !                end if
! !               
! !              end do ! i6
! !  
! !              
! !              i6=gFan_nu_index_end 
! ! !			 if ( ( abs(omega(i6)-omega(i6-1)) .gt. 0.00000367502  )  ) then  ! if i2 and i2-1  DO NOT correspond to the same degenerate  phonon freq
! ! 			 if ( ( abs(omega(i6)-omega(i6-1)) .gt. 0.00000367502  )  ) then  ! if i2 and i2-1  DO NOT correspond to the same degenerate  phonon freq
! ! 				 do j=firstj,i6-1
! ! 				   gFanaux1(i1,i2,i3,i4,i5,j) = ( sqrt(auxr1) ) / (dble(i6-1-firstj+1))
! ! 				 end do  
! ! 				 gFanaux1(i1,i2,i3,i4,i5,i6) = gFanaux2(i1,i2,i3,i4,i5,i6)
! ! 			 else                       ! if i2 and i2-1 correspond to the same degenerate phonon freq
! ! 			   auxr1 = auxr1 + (gFanaux2(i1,i2,i3,i4,i5,i6))**2
! ! 			   do j=firstj,i6
! ! 			     gFanaux1(i1,i2,i3,i4,i5,j) = ( sqrt(auxr1) ) / (dble(i6-firstj+1))
! ! 			   end do  
! ! 			 end if             
! ! 			 
! ! 			 if ((i1.eq.1) .and.(i2.lt.3) )then 
! ! 			   do i6=1,dimnu
! ! 			     write(244,'(I4,A,I4,A,I4,A,f12.8,A,f12.8)')  & 
! ! 			          & i1,' ',i2,' ',i6,' ',dble(gFanaux1(i1,i2,i3,i4,i5,i6)),' - ',omega(i6)
! ! 			   end do
! ! 			 end if  
! !              
! !           end do
! !         end do
! !       end do        
! !     end do  
! !   end do 
! !   close(244)  
! 
! gfanaux1(:,:,:,:,:,:)=gfanaux2(:,:,:,:,:,:) 
! 
! ! Average in electronic degenerate eigenvalues (1st band index)
!   gFan(:,:,:,:,:,:)=(0.d0,0.d0)
!   open(243,file='gFan-av3.dat',status='unknown')
!   do i2=1,gFan_j_index_end-gFan_j_index_beg+1
! 	do i3=1,k3
! 	  do i4=1,k4
! 		do i5=1,k5
! 		  do i6=1,gFan_nu_index_end 
! 		   
! 		     firstj=1; lastj=1
! 		     auxr1 = (gFanaux1(1,i2,i3,i4,i5,i6))**2
! 		     
! 		     do i1=2,gFan_i_index_end-gFan_i_index_beg+1-1 
!              
!                if ( ( abs(epsilon(i1)-epsilon(i1-1)) .gt. 0.00001d0 )  ) then  ! if i2 and i2-1  DO NOT correspond to the same degenerate electronic eigenvalue
!                  auxr1 = ( sqrt(auxr1) ) / (dble(lastj-firstj+1))
!                  do j=firstj,lastj
!                    gFan(j,i2,i3,i4,i5,i6)=dcmplx(auxr1)
!                  end do
!                  auxr1=(gFanaux1(i1,i2,i3,i4,i5,i6))**2
!                  firstj=i1; lastj=i1
!                else                       ! if i2 and i2-1 correspond to the same degenerate electronic eigenvalue
!                  auxr1 = auxr1 + (gFanaux1(i1,i2,i3,i4,i5,i6))**2
!                  lastj = lastj + 1   
!                end if
!               
!              end do ! i2
!              
!              i1=gFan_i_index_end-gFan_i_index_beg+1
! 			 if ( ( abs(epsilon(i1)-epsilon(i1-1)) .gt. 0.00001d0 )  ) then  ! if i2 and i2-1  DO NOT correspond to the same degenerate electronic eigenvalue
!                  auxr1 = ( sqrt(auxr1) ) / (dble(i1-1-firstj+1))
!                  do j=firstj,i1-1
!                    gFan(j,i2,i3,i4,i5,i6) = auxr1
!                  end do
! 				 gFan(i1,i2,i3,i4,i5,i6) = gFanaux1(i1,i2,i3,i4,i5,i6)
! 			 else                       ! if i2 and i2-1 correspond to the same degenerate eigenvalue
! 			   auxr1 = auxr1 + (gFanaux1(i1,i2,i3,i4,i5,i6))**2
! 			   do j=firstj,i1
! 			     gFan(j,i2,i3,i4,i5,i6) = ( sqrt(auxr1) ) / (dble(i1-firstj+1))
! 			   end do  
! 			 end if             
! 			 
! 			 if (i6.eq.1) then 
! 			   do i1=1,gFan_i_index_end-gFan_i_index_beg+1
! 			     write(243,'(I4,A,I4,A,I4,A,f12.8,A,f12.8)')  & 
! 			          & i1,' ',i2,' ',i6,' ',dble(gFan(i1,i2,i3,i4,i5,i6)),' - ',epsilon(i1)
! 			   end do
! 			 end if  
!              
!           end do
!         end do
!       end do        
!     end do  
!   end do 
!   close(243)
! 


gfan(:,:,:,:,:,:)=gfanaux1(:,:,:,:,:,:)

! WRITING THE FINAL gFan TO FILE
  open(243,file='gFan.dat',status='unknown')
  
!  write(243,*) gFan_i_index_beg,gFan_i_index_end,gFan_j_index_beg,gFan_j_index_end,k3,k4,k5,gFan_nu_index_beg,gFan_nu_index_end
 write(auxstrnum1,*)  gFan_i_index_beg;      call StripSpaces(auxstrnum1)
 write(auxstrnum2,*)  gFan_i_index_end;      call StripSpaces(auxstrnum2) 
 write(auxstrnum3,*)  gFan_j_index_beg;      call StripSpaces(auxstrnum3)                             
 write(auxstrnum4,*)  gFan_j_index_end;      call StripSpaces(auxstrnum4)                             
 write(auxstrnum5,*)  k3;                    call StripSpaces(auxstrnum5)                             
 write(auxstrnum6,*)  k4;                    call StripSpaces(auxstrnum6)  
 write(auxstrnum7,*)  k5;                    call StripSpaces(auxstrnum7)  
 write(auxstrnum8,*)  gFan_nu_index_beg;     call StripSpaces(auxstrnum8)  
 write(auxstrnum9,*)  gFan_nu_index_end;     call StripSpaces(auxstrnum9)   
 write(243,*,IOSTAT=stat)  trim(auxstrnum1),'  ',trim(auxstrnum2),'  ',trim(auxstrnum3),'  ', &
                         & trim(auxstrnum4),'  ',trim(auxstrnum5),'  ',trim(auxstrnum6),'  ', & 
                         & trim(auxstrnum7),'  ',trim(auxstrnum8),'  ',trim(auxstrnum9)
              
              
                
  do i1=1,gFan_i_index_end-gFan_i_index_beg+1 !k1
    do i2=1,gFan_j_index_end-gFan_j_index_beg+1 !k2
      do i3=1,k3
        do i4=1,k4
           do i5=1,k5
             do i6=1,gFan_nu_index_end! gFan_nu_index_beg,gFan_nu_index_end !1,k6
             
               write(auxstrnum1,*)  i1+gFan_i_index_beg-1; call StripSpaces(auxstrnum1)
               write(auxstrnum2,*)  i2+gFan_j_index_beg-1; call StripSpaces(auxstrnum2) 
               write(auxstrnum3,*)  i3;                    call StripSpaces(auxstrnum3)                             
               write(auxstrnum4,*)  i4;                    call StripSpaces(auxstrnum4)                             
               write(auxstrnum5,*)  i5;                    call StripSpaces(auxstrnum5)                             
               write(auxstrnum6,*)  i6;                    call StripSpaces(auxstrnum6)                                          
               write(243,*,IOSTAT=stat)  trim(auxstrnum1),' ',trim(auxstrnum2),' ',trim(auxstrnum3),' ', &
                                       & trim(auxstrnum4),' ',trim(auxstrnum5),' ',trim(auxstrnum6),' ', & 
!                   &  dble ( gFan ( i1+gFan_i_index_beg-1, i2+gFan_j_index_beg-1, i3,i4,i5,i6-gFan_nu_index_beg+1 ) ), ' ', &
!                   & dimag ( gFan ( i1+gFan_i_index_beg-1, i2+gFan_j_index_beg-1, i3,i4,i5,i6-gFan_nu_index_beg+1 ) )
                   &  dble ( gFan ( i1, i2, i3,i4,i5,i6) ), ' ', &
                   & dimag ( gFan ( i1, i2, i3,i4,i5,i6) )
             end do
           end do
         end do !
       end do        
    end do  
  end do 
  close(243)

  deallocate(gFanaux1,gFanaux2)
  
  
 l1 = gFan_i_index_end
 if ( (l1 .lt. sizeepsilon)) then
   if ( abs(epsilon(l1)-epsilon(l1+1)) .lt. 0.00001d0  ) then
 
      dif=-1
      if ( (l1+1 .lt. sizeepsilon) .and. ( abs(epsilon(l1)-epsilon(l1+2)) .ge. 0.00001d0  ) ) then
         dif=1 
      else if ( (l1+2 .lt. sizeepsilon) .and. ( abs(epsilon(l1)-epsilon(l1+3)) .ge. 0.00001d0  )) then
         dif=2
      else if ( (l1+3 .lt. sizeepsilon) .and. ( abs(epsilon(l1)-epsilon(l1+4)) .ge. 0.00001d0  )) then
         dif=3 
      else if ( (l1+4 .lt. sizeepsilon) .and. ( abs(epsilon(l1)-epsilon(l1+5)) .ge. 0.00001d0  )) then
         dif=4 
      end if                 

      write(auxstr1,*) gFan_i_index_end; call StripSpaces(auxstr1);  write(*,*)  
      write(*,*) '    WARNING: The highest gFan index of the 1st band (',trim(auxstr1),', see first'
      write(*,*) '     line of ',trim(gFan_file_name),') seems to be in the middle of a set of'
      write(*,*) '     degenerate eigenvalues: '
      write(*,*) '       band(',trim(auxstr1),')=',epsilon(gFan_i_index_end),' Hartree, '
      write(auxstr1,*) gFan_i_index_end+1; call StripSpaces(auxstr1) 
      write(*,*) '       band(',trim(auxstr1),')=',epsilon(gFan_i_index_end+1),' Hartree. '
      if ( dif .gt. 0) then
        write(*,*) '     Please, calculate the Fan electron-phonon matrix elements (1st band) '
        write(auxstr1,*) gFan_i_index_end+dif; call StripSpaces(auxstr1) 
        write(*,*) '     ', trim(auxstr1),' up to for higher accuracy.'
      end if
      write(*,*)
    end if  
  end if   


 l1 = gFan_j_index_end
 if ( (l1 .lt. sizeepsilon)) then
   if ( abs(epsilon(l1)-epsilon(l1+1)) .lt. 0.00001d0  ) then
 
      dif=-1
      if ( (l1+1 .lt. sizeepsilon) .and. ( abs(epsilon(l1)-epsilon(l1+2)) .ge. 0.00001d0  ) ) then
         dif=1 
      else if ( (l1+2 .lt. sizeepsilon) .and. ( abs(epsilon(l1)-epsilon(l1+3)) .ge. 0.00001d0  )) then
         dif=2
      else if ( (l1+3 .lt. sizeepsilon) .and. ( abs(epsilon(l1)-epsilon(l1+4)) .ge. 0.00001d0  )) then
         dif=3 
      else if ( (l1+4 .lt. sizeepsilon) .and. ( abs(epsilon(l1)-epsilon(l1+5)) .ge. 0.00001d0  )) then
         dif=4 
      end if                 
      
      write(auxstr1,*) gFan_j_index_end; call StripSpaces(auxstr1);  write(*,*)  
      write(*,*) '    WARNING: The highest gFan index of the 1st band (',trim(auxstr1),', see first'
      write(*,*) '     line of ',trim(gFan_file_name),') seems to be in the middle of a set of'
      write(*,*) '     degenerate eigenvalues: '
      write(*,*) '       band(',trim(auxstr1),')=',epsilon(gFan_j_index_end),' Hartree, '
      write(auxstr1,*) gFan_j_index_end+1; call StripSpaces(auxstr1) 
      write(*,*) '       band(',trim(auxstr1),')=',epsilon(gFan_j_index_end+1),' Hartree. '
      if ( dif .gt. 0) then
        write(*,*) '     If you want to renormalize the highest bands, please calculate the '
        write(*,*) '     Fan electron-phonon matrix elements (2nd band) '
        write(auxstr1,*) gFan_j_index_end+dif; call StripSpaces(auxstr1) 
        write(*,*) '     ', trim(auxstr1),' up to for higher accuracy.'
      end if
      write(*,*)
    end if  
  end if   


 end subroutine read_gFan
 
 
!-------------------------------------------------------------------------------------------


 !subroutine write_arbitrary_gFa n(numbands,dimnu,gFan_file_name)
!  
!    ! System parameters
!   integer, Intent(In) :: numbands
!   integer, Intent(In) :: dimnu
!   character(100), Intent(In) :: gFan_file_name
!  
!   
!   ! Local variables
!   integer ::  i,j,nu
! 
!   
!   open(541,file=gFan_file_name,status='unknown')
!   
!   do i=1,numbands
!     do j=1,numbands
!       do nu=1,dimnu
!         write(541,*) i,j,nu, 1.2345d0, 0.d0
!       end do
!     end do  
!   end do 
!   
!   close(541)
! 
! 
!  end subroutine write_arbitrary_gFan
 
 
!-------------------------------------------------------------------------------------------
  
 subroutine calculate_gDW (sizeepsilon, dimnu, number_of_atoms, band_occupation, gFan_i_index_beg,  &
           & gFan_i_index_end, gFan_j_index_beg, gFan_j_index_end, &
           & gDW_i_index_beg,gDW_i_index_end,NstatesforgDW,n_mode_beg,n_mode_end,epsilon,omega,xi,U,mass, gFan, gDW)
 
  ! This subroutine calculates the Debye-Waller matrix elements using the formula below (a la Giustino, but with a correct formula)
!  
!    \big(g^{DW} \big)^{ii}_{\nu}  = \sum\limits_{\substack{s\alpha  \\s'\alpha'}}
! \frac{ \big( \  U^{\alpha,s}_{\nu} U^{\alpha',s}_{\nu} + U^{\alpha,s'}_{\nu} U^{\alpha',s'}_{\nu}  \ \big)}{4 \ \hbar \ \omega_{\nu}}   \
! \textrm{Re} \left[ \  
! \sum\limits_{\substack{j=1  \\  j \neq i}}^{\infty}
!  \  \frac{h^{j,i;\alpha,s}h^{i,j;\alpha',s'} + h^{j,i;\alpha',s'}h^{i,j;\alpha,s}}{ \epsilon^0_j  - \epsilon^0_i }  \ \right] \ , \\
! with 
!   h^{i,j;\alpha,s} \ \equiv \ \sum_{\nu} \  \sqrt{M_s \omega_{\nu}} \ \xi^{\alpha,s}_{\nu} \ g^{ij}_{\nu} \ ; \qquad  {\bm{U}}^{\textrm{ } s}_{  \nu}  \ \equiv \ {{\bm{\xi}}^{\textrm{ } s}_{\nu}}/{\sqrt{M_s} }  \ ,
! 
!  IMPORTANT: We consider h^{i,j}=h^{j,i} because we consider a nonperiodic system, where gFan are real and g^{ij}=g^{ji}.
 
 
   Integer, Intent(In)  :: sizeepsilon              ! Number of evaluated electronic eigenstates (occupied and unoccupied)
   Integer, Intent(In)  :: dimnu                    ! Number of phonon branches of the system
   Integer, Intent(In)  :: number_of_atoms          ! Number of atoms of the system
   integer, Intent(In)  :: band_occupation          ! Number of electrons (1 or 2) per band in the file output of QE
   Integer, Intent(In)  :: gFan_i_index_beg         ! Initial wavefunction index for < psi_i | for the calculation of (g^{DW})^{ii}_{nu}
   Integer, Intent(In)  :: gFan_i_index_end         ! Final wavefunction index < psi_i | for the calculation of (g^{DW})^{ii}_{nu}
   Integer, Intent(In)  :: gFan_j_index_beg         ! Initial wavefunction index for | psi_j >for the calculation of (g^{DW})^{ii}_{nu}
   Integer, Intent(In)  :: gFan_j_index_end         ! Final wavefunction index  | psi_j >  for the calculation of (g^{DW})^{ii}_{nu}    
   Integer, Intent(In)  :: gDW_i_index_beg          ! Initial index for the calculation of (g^{DW})^{ii}_{nu}
   Integer, Intent(In)  :: gDW_i_index_end          ! Final index for the calculation of (g^{DW})^{ii}_{nu}
   Integer, Intent(In)  :: NstatesforgDW            ! Number of states (bands) considered in the summation to calculate gDW
   Integer, Intent(In)  :: n_mode_beg,n_mode_end    ! Initial and final indices of phonon branch
   real(8), Intent(In)  :: epsilon(sizeepsilon)     ! Electronic eigenvalues
   real(8), Intent(In)  :: omega(dimnu)             ! Phonon frequencies
   real(8), Intent(In)  :: xi(dimnu,dimnu), U(dimnu,dimnu) ! xi and U eigenvectors of the dynamical matrix
   real(8), Intent(In)  :: mass(number_of_atoms,2)  ! Masses of the atoms
   complex(8), Intent(In) :: gFan(gFan_i_index_end-gFan_i_index_beg+1, & 
                               & gFan_j_index_end-gFan_j_index_beg+1, 1,1,1, dimnu)! Fan electron-phonon matrix elements
   complex(8), Intent(Out) :: gDW(gDW_i_index_end-gDW_i_index_beg+1,1,1,1,dimnu)  ! Debye-Waller electron-phonon matrix elements
 

  ! Local variables
  Integer :: i, i1, i2, j, k,l, m, alpha, alpha1, alpha2, s, s1, s2, nu, dif,nup
  real(8) :: auxr1, auxr2, auxr3, auxr4
  complex(8) :: auxc1, auxc2, auxc3, auxc4
  complex(8), allocatable :: h(:,:,:,:)  
  character(20) :: auxstr1, auxstr2
  character(100) :: auxstr3, auxstr4
  

 
   ! This subroutine reads properly: mass; U;  sizeepsilon, dimnu, number_of_atoms; gFan; omega; all six g indices; epsilon  
  
  allocate(h(gFan_i_index_end-gFan_i_index_beg+1, gFan_j_index_end-gFan_j_index_beg+1, 3, number_of_atoms))
  
  write(*,*) '  *** Now calculating the Debye-Waller electron-phonon matrix elements. ***'

  gDW(:,:,:,:,:)=(0.d0,0.d0);   h(:,:,:,:) = (0.d0,0.d0)
    
!        i=1;j=1;nu=18
!        auxc1=0.d0
!        do alpha=1,3
!           do s=1,number_of_atoms
!               write(*,*) s, alpha, sqrt(omega(nu)), dble(xi(3*(s-1)+alpha,nu)), dble(sqrt(dble(mass(s,2))))
!               write(*,*) dble(sqrt(omega(nu)) * xi(3*(s-1)+alpha,nu) * gFan(i,j,1,1,1,nu)* sqrt(mass(s,2)));write(*,*)
!   auxc1=auxc1+dble(sqrt(omega(nu)) * xi(3*(s-1)+alpha,nu) * gFan(i,j,1,1,1,nu)* sqrt(mass(s,2)))
!            end do
!           end do
!           write(*,*)'final',dble(auxc1)

   
  ! Calculation of the h's: h^{i,j;\alpha,s} \ \equiv \ \sum_{\nu} \  \sqrt{M_s \omega_{\nu}} \ \xi^{\alpha,s}_{\nu} \ g^{ij}_{\nu} \ 
  open(623,file='h.out',status='unknown')
  do i=1,gFan_i_index_end-gFan_i_index_beg+1
     do j=1,gFan_j_index_end-gFan_j_index_beg+1
       auxr4=0.d0; nup=0
       do alpha=1,3
          do s=1,number_of_atoms
 
              auxr1=0.d0         
              do nu=1,dimnu ! AQUI NO PUEDES TRUNCAR, si lo haces no se cumple la relacion de ortonormalizacion, que es imprescindible para que las ecuaciones sean validas !!
                 auxr1 = auxr1 + sqrt(omega(nu)) * xi(3*(s-1)+alpha,nu) * gFan(i,j,1,1,1,nu)  ! Beware of this notation: the nu-index of xi and U runs from 1 to dimnu, but that of gFan runs from 1 to dimnu-6; Hence e.g. the nu=7 corresponds to U(I,7) and to gFan(I,J,S,K,q,1).
               !  write(*,*) 'i', i,'j',j  , 'nu', nu 
               !  write(*,*)  'xi', xi(3*(s-1)+alpha,nu)
               !  write(*,*)'gFan(',i,j,nu,')=', dble(gFan(i,j,1,1,1,nu))
               !  write(*,*) 'sqrtOmeg',sqrt(omega(nu))
               !  write(*,*) 'sumaparcial',auxr1; write(*,*) 
              end do

              h(i,j,alpha,s) = auxr1 * sqrt(mass(s,2))
              write(623,'(I5,A,I5,A,I2,A,I5,A,G16.9)')  i+gFan_i_index_beg-1,' ',j+gFan_j_index_beg-1, &
                                                       & ' ',alpha,' ',s,'   ', dble(h(i,j,alpha,s))
          end do
       end do
       
       auxr5=0.d0
       do nu=1,dimnu
          if (abs(gFan(i,j,1,1,1,nu)).gt. auxr5) auxr5=abs(gFan(i,j,1,1,1,nu)) 
       end do                                                     
                                                                                
     end do
  end do   
  close(623)
  
  

 !  
!   ! Calculation of the U´s
!   do s=1,number_of_atoms
!      do alpha=1,3
!         do nu=1,dimnu
!            U(3*(s-1)+alpha,nu) = xi(3*(s-1)+alpha,nu)/sqrt(mass(s,2))
!         end do
!      end do
!   end do   
  
  
  
  ! Calculation of the gDW
  dif = gDW_i_index_beg - gFan_i_index_beg
  if (dif.ne.0) then
    write(*,*)'ERROR: The solution as a function of "dif" is still to be coded'
    call exit(1)
  end if
  

    
  do i=1,gDW_i_index_end-gDW_i_index_beg+1
    do nu=n_mode_beg,n_mode_end
    
      auxc3 = (0.d0,0.d0)
      
      do s1=1,number_of_atoms
        do alpha1=1,3
        
  ! write(*,*)'h(',i+dif,j,alpha1,s1,')=',dble(h(i+dif,j,alpha1,s1))
          
          do s2=1,number_of_atoms
            do alpha2=1,3
    
               auxc1=(0.d0,0.d0)
               
               do j=1,NstatesforgDW-gFan_j_index_beg+1 !UUUU 1,28
               
                 auxr3 = epsilon(j)-epsilon(i+gFan_i_index_beg-1)
                 if ( abs(auxr3) .gt. 0.00001d0 ) then
                   auxc2 = h(i,j,alpha1,s1)*h(i,j,alpha2,s2)  ! (  h(i,j,alpha1,s1)*h(i,j,alpha2,s2) + h(i,j,alpha2,s2)*h(i,j,alpha1,s1) )
                   auxc2 = auxc2 / (dcmplx(auxr3))
                   auxc1  = auxc1 + auxc2      

               auxc4 =  ( U(3*(s1-1)+alpha1,nu) * U(3*(s1-1)+alpha2,nu) + &
                        & U(3*(s2-1)+alpha1,nu) * U(3*(s2-1)+alpha2,nu) )

! if((i.eq.10).and.(nu.eq.18))then
!   write(*,*)'i',i,'j',j,'nu',nu
!   write(*,*) 'alpha1',alpha1,'alpha2',alpha2,'s1',s1,'s2',s2
!    write(*,*) 'hh',h(i,j,alpha1,s1)*h(i,j,alpha2,s2),'h',dble(h(i,j,alpha1,s1)),'h',dble((h(i,j,alpha2,s2)))                                                                          
!     write(*,*) 'UU+UU',  auxc4  
!     write(*,*) U(3*(s1-1)+alpha1,nu),      U(3*(s1-1)+alpha2,nu), U(3*(s2-1)+alpha1,nu), U(3*(s2-1)+alpha2,nu)
!    write(*,*)'eps-eps',auxr3
!    write(*,*) '8w_Q',  8.d0*omega(nu)  
!    write(*,*) 'term parcial',( auxc4 * auxc2 )/ ( 8.d0*omega(nu) )
!    write(*,*)   
!    write(*,*)'den1', mass(s1,2)*omega(nu),' - den2',mass(s2,2)*omega(nu)
!    write(*,*)'xis', xi(3*(s1-1)+alpha1,nu), xi(3*(s1-1)+alpha2,nu), xi(3*(s2-1)+alpha1,nu), xi(3*(s2-1)+alpha2,nu)
!    write(*,*)'xixi',xi(3*(s1-1)+alpha1,nu)*xi(3*(s1-1)+alpha2,nu)+ xi(3*(s2-1)+alpha1,nu)*xi(3*(s2-1)+alpha2,nu)
!     write(*,*) '(UU+UU)/8wQ',  auxc4  /(8.d0*omega(nu) )
!    write(*,*);write(*,*)                                                   
! end if                                

         
                 end if
               end do ! j
               auxc1 = dble(auxc1)   !*band_occupation ! FALSE (spins make both terms orthogonal): If there are 2 electrons per band, then there are two bands on each other, so one must count two bands
           
               auxc4 =  ( U(3*(s1-1)+alpha1,nu) * U(3*(s1-1)+alpha2,nu) + &
                        & U(3*(s2-1)+alpha1,nu) * U(3*(s2-1)+alpha2,nu) )

                                                  
                        
                        
! write(*,*) 'UU+UU',dble(auxc4) 
 !  write(*,*) ' nu',nu,'s',s1,s2,'alpha',alpha1,alpha2                      
!   write(*,*) 'U:',  U(3*(s1-1)+alpha1,nu),U(3*(s1-1)+alpha2,nu) , U(3*(s2-1)+alpha1,nu) , U(3*(s2-1)+alpha2,nu)                      
!   write(*,*) 'UU:',  ( U(3*(s1-1)+alpha1,nu) * U(3*(s1-1)+alpha2,nu) +  U(3*(s2-1)+alpha1,nu) * U(3*(s2-1)+alpha2,nu) ) , &
!    &  'coc:',(U(3*(s1-1)+alpha1,nu) * U(3*(s1-1)+alpha2,nu) +  U(3*(s2-1)+alpha1,nu) * U(3*(s2-1)+alpha2,nu) )/(2.d0*omega(nu) )
!    write(*,*)
               auxc4 = auxc4 / ( 2.d0*omega(nu) )
 
! write(*,*) '(UU+UU)hh/(eps-eps)',dble(auxc4 *auxc1)
               
               auxc3 = auxc3 + ( auxc4 * auxc1 )

! write(*,*) 'partial',dble(auxc3) ; write(*,*)

              end do ! alpha2
            end do ! s2
          end do ! alpha1
        end do ! s1
 
! write(*,*) 'gDW(',i,1,1,1,nu-n_mode_beg+1,')=',auxc3       
        gDW(i,1,1,1,nu-n_mode_beg+1) = auxc3
            
    end do ! nu
  end do ! i

  
 ! We write the gDW to a file. According to the Maat format:
!       * gDW.dat:  ! First row: first band index of the DW matrix elements, last band index of the DW matrix elements, number of spins, number of k-points, number of q-points, first phonon-branch index, last phonon-branch index
!                     then: first 5 columns: band index, spin-index, k-index, q-index and phonon branch-index, respectively. 5th and 6th columns: Real and imaginary parts of the Debye-Waller electron-phonon matrix elements 

  write(*,*) '  *** Now storing Debye-Waller electron-phonon matrix elements in gDW.dat ***'

  open(123,file='gDW.dat',status='unknown')
  write(123,'(I7,A,I7,A,I7,A,I7)') gDW_i_index_beg,' ',gDW_i_index_end,'  1  1  1  ',n_mode_beg,' ',n_mode_end  
  do i=1,gDW_i_index_end-gDW_i_index_beg+1
	do nu=1,dimnu-n_mode_beg+1
	
	  write(auxstr1,*) i+gDW_i_index_beg-1; call StripSpaces(auxstr1)
	  write(auxstr2,*) nu+n_mode_beg-1; call StripSpaces(auxstr2)          
  
      write(123,*) trim(auxstr1), ' 1 1 1 ', trim(auxstr2),' ', &
					& dble(gDW(i,1,1,1,nu)),' ', dimag(gDW(i,1,1,1,nu))
	end do
  end do  
  close(123)
  
  ! We copy the file with the Debye-Waller electron-phonon matrix elements to a file specifying the number of states considered in the calculation of gDW
  write(auxstr3,'(A,I8,A)') 'gDW-',NstatesforgDW,'states.dat'; call StripSpaces(auxstr3)
  write(auxstr4,*) 'cp gDW.dat '//trim(auxstr3)
  call system( auxstr4 )
  

  write(*,*); write(*,*) '                ******* Calculations finished correctly. *******'
  write(*,*) '         ******* The input for the Maat program was generated. *******'
  write(*,*); write(*,*) 

  deallocate (h)

 end  subroutine calculate_gDW 


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
