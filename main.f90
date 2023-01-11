! Program: nuclear_energies
! By: Nathan Crawford
!-----------------------------------------------------------------------------
! This program calculates the binding energies of various isotopes using the 
! semi empirical mass formula. In order to do so it has to find the best fit 
! linear parameters for the semi empirical mass formula which it does by solving 
! the matrix equation ax = b. 

! In order to do this, the program uses a technique known as LU decomposition. 
! After the program solves this matrix equation by 
! LU decomposition, finding it's inverse and multiplying it on both sides it will 
! print to screen the best fit parameters of each term in the semi empirical mass 
! formula and the corresponding error of each term. 

! We can now calculate the binding energies which it will write into a file 
! as well as the ones calculated experimentally. 
! It will write both these binding energies and their corresponding 
! errors into a file called results.dat 
! 
! The second part of this project will use the model to make predictions. 
! It will find the valley of stability, meaning it will loop through proton numbers 
! 1-118 and find the most stable isotope of each proton number, the lowest binding 
! energy per nucleon. 
! 
! It will then find the location of the neutron dripline. 
! It will loop through the same number of protons and find the largest neutron 
! number such that the sepration energy for that proton number is positive. 
!-----------------------------------------------------------------------------
program nuclear_energies
use types
use read_write, only : read_exp_data, write_predictions, write_advanced
use nuclear_model, only : find_best_parameters, stable_isotopes, find_dripline
implicit none

integer, allocatable :: n_protons(:), n_neutrons(:), stable_neutrons(:), neutron_dripline(:)
real(dp), allocatable :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:) 
real(dp), allocatable :: binding_energies(:), lowest_energies(:), separation(:)  


call read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)

call find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)

call write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)

call stable_isotopes(binding_energies, lowest_energies, n_neutrons, stable_neutrons, c_parameters)  
call find_dripline(separation, neutron_dripline, c_parameters, n_neutrons) 
call write_advanced(stable_neutrons, neutron_dripline)  


end program nuclear_energies
