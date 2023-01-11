!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This module reads the data from the experimental file then 
!! writes the experimental binding energies, errors, protons, neutrons, 
!! theoretical binding energies and theoretical errors.
!!----------------------------------------------------------------------
!! Included subroutines:
!! read_exp_data
!! write_predictions 
!! write_advanced
!!----------------------------------------------------------------------
module read_write

use types
use nuclear_model, only : semi_empirical_mass, semi_empirical_error

implicit none

private
public :: read_exp_data, write_predictions, write_advanced

contains

!-----------------------------------------------------------------------
!! Subroutine: read_exp_data
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Opens the experimental file and reads the experimental data.
!!----------------------------------------------------------------------
!! Output:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
subroutine read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
    implicit none
    integer, intent(out), allocatable :: n_protons(:), n_neutrons(:)
    real(dp), intent(out), allocatable :: exp_values(:), uncertainties(:)

    character(len=128) :: filename 
    integer, allocatable :: dummy_int(:) 
    character, allocatable :: dummy_string(:)
    logical :: file_exists
    integer :: file_unit, unit, M, i 
    real(dp), allocatable :: dummy_real(:)

    print *, 'This program will calculate the binding energies of a given' 
    print *, 'isotope and write the experimental and theoretical results' 
    print *, 'to a file.' 
    print *, 

    print *, 'please provide the file name with the experimental data'
    read(*, '(a)') filename
    
    ! when trying to open a file provided by the user it's good practice to
    ! check if the file exists in the current directory
    inquire(file=trim(filename), exist=file_exists)
    !If the file exists then it will open it and begin reading
    if (file_exists) then
        open(newunit=unit,file=filename, status = 'old') 
        
        !M is the number of data points given at the top of the file
        read(unit, *) M 
        
       
        
        !Now we know how much we need to allocate each array
        allocate(dummy_string(1:M)) 
        allocate(dummy_int(1:M)) 
        allocate(dummy_real(1:M))
        allocate(n_protons(1:M)) 
        allocate(n_neutrons(1:M)) 
        allocate(exp_values(1:M)) 
        allocate(uncertainties(1:M))
        
        !Create dummy variable to skip first 2 lines
        read(unit, *) 
        read(unit, *) 
        
        !Loop to read relevant data while skipping irrelevant ones. I am using dummy to indicate 
        !if it is an irrelevant data point
        do i = 1, M
            read(unit, *) dummy_int(i), dummy_string(i), dummy_int(i), n_neutrons(i), &
            n_protons(i), exp_values(i), dummy_real(i), uncertainties(i) 
        end do
        
       
        
        close(unit)
        
    else
        print *, "File does not exist, please enter name of the correct data file." 
    endif
end subroutine read_exp_data

!-----------------------------------------------------------------------
!! Subroutine: write_predictions
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This subroutine writes the experimental binding energies, error, protons, 
!! neutrons, theoretical binding energies, theoretical error into a file 
!! called results.dat
!!----------------------------------------------------------------------
!! Input:
!!
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! covariance       real        Array containing the elements of the covariance matrix
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------
subroutine write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)
    implicit none
    real(dp), intent(in) :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:) 
    real(dp) :: binding_energy, binding_energy_error
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    integer:: unit, k
    character(len=*), parameter :: file_name = 'results.dat' 
    
    open(newunit=unit,file=file_name) 
    !Header for the file
    write(unit,'(6a23)') 'Protons', 'Neutrons' , 'Exp Binding Energies', 'Exp Error', 'Binding Energy', 'Theoretical Error'
   
    
    !Loop through and write every binding energy and error
    do k = 1, size(exp_values)
       
        binding_energy = semi_empirical_mass(c_parameters, n_protons(k), n_neutrons(k))
        binding_energy_error = semi_empirical_error(covariance, n_protons(k), n_neutrons(k))
    write(unit,'(2i8, 4e28.16, 5X)') n_protons(k), n_neutrons(k), exp_values(k), uncertainties(k), binding_energy, &
    binding_energy_error  
    end do 
    
    close(unit)
    print *,
    print *, 'theoretical binding energies were written in ', file_name 
    print *,  
    print *, 'Now we will calculate the valley of stability and neutron dripline'
end subroutine write_predictions



!-----------------------------------------------------------------------
!! Subroutine: write_advanced
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This subroutine opens a file called results_advanced.dat and writes 
!! the positions of the neutron dripline and valley of stability into a file.
!!----------------------------------------------------------------------
!! Input:
!!
!! stable_neutrons     real        Array containing neutrons that correspond to lowest binding energy per nucleon
!!                                 of a given proton number z. 
!!
!! neutron_dripline    real        Array containing the position of the neutron dripline
!-----------------------------------------------------------------------
subroutine write_advanced(stable_neutrons, neutron_dripline) 
implicit none 
integer :: stable_neutrons(:), z, neutron_dripline(:), unit    
character(len=*), parameter :: file_name = 'results_advanced.dat' 

!Open file
open(newunit = unit, file = file_name) 
!Write the header of the file
write(unit,'(3a28)') 'Protons', 'Stable Isotope Position' , 'Neutron Dripline Position' 
!Loop through and write protons, valley of stability, and position of dripline.
do z = 1, size(stable_neutrons)
    write(unit,'(3i24)') z, stable_neutrons(z), neutron_dripline(z) 
end do 
close(unit)
print *, 
print *, 'Results written in results_advanced.dat'
end subroutine write_advanced
    
end module read_write
