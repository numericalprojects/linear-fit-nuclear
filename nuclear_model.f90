!-----------------------------------------------------------------------
!Module: nuclear_model
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This is the physics module that finds the best fit parameters in the 
!! semi empirical mass formula by filling the alpha and beta arrays
!! and calling subroutines from the linear algebra module. 
!! It then prints these best fit parameters to the screen and calculates the 
!! binding energies of the isotopes from the experimental file and calculuates 
!! the corresponding theoretical error.
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! find_best_parameters 
!! stable_isotopes 
!! find_dripline
!!----------------------------------------------------------------------
!! Included functions:
!! semi_empirical_mass 
!! semi_empirical_error
!!----------------------------------------------------------------------
module nuclear_model
use types
use linear_algebra, only : solve_linear_system
implicit none

private

public :: find_best_parameters, semi_empirical_mass, semi_empirical_error, stable_isotopes, find_dripline
contains


!-----------------------------------------------------------------------
!! Subroutine: find_best_parameters
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This subroutine finds the best fit parameters for the semi empirical 
!! mass formula by calling a subroutine to fill the alpha and beta arrays 
!! and then solving the matrix equation by calling another subroutine. 
!! It then calls another subroutine to have the results printed to screen.
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! c_parameters     real        Array containing the semi-empirical mass formula parameters
!! covariance       real        Array containing the covariance matrix of the parameters
!-----------------------------------------------------------------------
subroutine find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out), allocatable ::  c_parameters(:), covariance(:,:)

   
    integer, parameter :: n_parmaters = 5
    !Allocate arrays to the number of parameters in the semi empirical mass 
    !formula
    real(dp) :: alpha(1:n_parmaters,1:n_parmaters), beta(1:n_parmaters)
    allocate(c_parameters(1:n_parmaters)) 
    allocate(covariance(1:n_parmaters, 1:n_parmaters))
   
    
    !Call a subroutine to fill the alpha and beta arrays 
    call construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)

    
    !Call subroutine to solve the matrix equation
    call solve_linear_system(alpha,beta,c_parameters,covariance)
    
    ! Now just print the parameters (with it's uncertainties) to screen
    print*, 'displaying results'
    call print_best_parameters(c_parameters,covariance)
end subroutine find_best_parameters

!-----------------------------------------------------------------------
!! Subroutine: construct_alpha_beta
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Fills the alpha and beta arrays according to the formula in the 
!! least squares technique method. 
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! alpha            real        Array containing the alpha matrix
!! beta             real        Array containing the beta vector
!-----------------------------------------------------------------------
subroutine construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out) :: alpha(:,:), beta(:)

    integer :: n_data, n_parmaters, alpha_shape(1:2), i, j, k
    real(dp):: linear_terms(1:size(beta))
    
    ! Check if the alpha array is a square matrix
    ! Also check that beta has the same number of elements as alpha has rows (or columns) 
    
    !Shape is an intrinsic function that returns the columns and rows of a matrix. 
    !Useful for checking the size of a matrix
    alpha_shape = shape(alpha)
    n_data = size(uncertainties)
    n_parmaters = alpha_shape(1)
    
    if(alpha_shape(1) /= alpha_shape(2)) then 
        print*, "Alpha is not a square matrix. Make sure alpha is square before continuing."
        stop 
    end if 
    if(size(beta) /= alpha_shape(1)) then 
        print*, "Elements in beta is not equal to the number of columns/rows of alpha." 
        print*, "Please make sure beta has the appropriate number of elements."
        stop 
    end if 
    
    alpha = 0._dp
    beta = 0._dp
    
    !According the least squares technique the alpha matrix is 
    !A_ij = sum of k = 1 to # of data points, f_i(Z_k, N_k) * f_j(Z_k, N_k)/o_k^2
    !the f functions are the terms multiplying the best fit parameters. Z and N are 
    !protons and neutrons respectively and the o_k is the experimental uncertainty. 
    
    !Similarly the beta array is  
    !B_i = sum of k = 1 to # of data points f_i(Z_k, N_k) * BE_exp(Z_k, N_k)/o_k^2 
    !Here BE_exp is the experimental binding energy from the data file. 
    
    !We can see that for a single element we have to loop through all the data points 
    !So from outer most to inner most it goes i, j, k where i and j are columns and rows 
    !to indicate an element and k is a data point from the file.
    
    do i = 1, size(alpha, 1) 
        do j = 1, size(alpha, 2)
            do k=1,n_data
        ! The subroutine below should return the f_\alpha(Z_i,N_i) terms defined in
        ! the README file
                call calculate_linear_termns(n_protons(k), n_neutrons(k), linear_terms) 
                alpha(i, j) = alpha(i, j) + linear_terms(i) * linear_terms(j)/(uncertainties(k) ** 2.0_dp) 
            end do 
        end do 
    end do 
    
    do i = 1, size(beta) 
        do k = 1, n_data 
        call calculate_linear_termns(n_protons(k), n_neutrons(k), linear_terms)
        beta(i) = beta(i) + linear_terms(i) * exp_values(k)/(uncertainties(k) ** 2.0_dp)
        end do
       
    enddo
end subroutine construct_alpha_beta

!-----------------------------------------------------------------------
!! Subroutine: calculate_linear_termns
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Calculates the terms multiplying the best fit parameters in the semi 
!! empirical mass formula that are functions of Z protons and N neutrons.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z                integer     number of protons in an isotope
!! N                integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! linear_terms        real        Array containing the linear terms in the semi-empirical mass formula
!-----------------------------------------------------------------------
subroutine calculate_linear_termns(Z, N, linear_terms)
    implicit none
    integer, intent(in) :: Z, N
    real(dp), intent(out) :: linear_terms(:)

    ! We could write down all the formulas for each term here. However, in
    ! order to keep the code readable and easy to understand  we'll  separate
    ! them into different functions
    linear_terms(1) = volume_term(Z,N)
    linear_terms(2) = surface_term(Z,N)
    linear_terms(3) = asymmetry_term(Z,N)
    linear_terms(4) = coulomb_term(Z,N)
    linear_terms(5) = pairing_term(Z,N)
  
end subroutine calculate_linear_termns

!-----------------------------------------------------------------------
!! function: volume_term
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates the volume term
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        volume term
!-----------------------------------------------------------------------
real(dp) function volume_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N 
    !volume term is A in the semi empirical mass formula which is Z + N
    r = Z + N 
end function volume_term

!-----------------------------------------------------------------------
!! function: surface_term
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculate the surface term
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        surface term
!-----------------------------------------------------------------------
real(dp) function surface_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    !Surface term is A^2/3. A is protons + neutrons 
    r = (Z+N)**(2.0_dp/3.0_dp)
end function surface_term

!-----------------------------------------------------------------------
!! function: asymmetry_term
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates the asymmetry term
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        asymmetry term
!-----------------------------------------------------------------------
real(dp) function asymmetry_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N 
    !Asymmetry term is (Neutrons - Protons)^2/(A) 
    !A is protons + neutrons
    r = (N - Z)**(2.0_dp)/(Z + N)
end function asymmetry_term

!-----------------------------------------------------------------------
!! function: coulomb_term
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates the coulomb term
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        coulomb term
!-----------------------------------------------------------------------
real(dp) function coulomb_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    !Formula for coulomb term is Protons(Protons - 1)/A^1/3 
    !A is protons + neutrons
    r = Z*(Z-1)/((Z+N)**(1.0_dp/3.0_dp))
end function coulomb_term

!-----------------------------------------------------------------------
!! function: pairing_term
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates the pairing term using another function I wrote for the 
!! different cases of protons, neutrons and nucleons. 
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        pairing term
!-----------------------------------------------------------------------
real(dp) function pairing_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    !formula for pairing term is A^-3/4 * delta(Protons, Neutrons) 
    !A is protons + neutrons. 
    !Delta is a function that depends on Z and N and is defined below 
    !with information on the different cases.
     
     r = (Z+N)**(-1.0_dp*3.0_dp/4.0_dp) * delta(Z,N)
    
    ! you might wanna define another function for the \delta(Z,N) factor in
    ! this term
end function pairing_term 

!-----------------------------------------------------------------------
!! function: delta
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates the delta function term in the pairing term function above.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        kronnecker delta result
!-----------------------------------------------------------------------
real(dp) function delta(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
     
     !Delta is 1 is Z and N are even, -1 if Z and N are odd and 0 if 
     !Z + N is odd.
     
     !Modulo operator is useful for this because if you divide 2 numbers 
     !by 2 and there's a remainder then it's odd and if there is no remainder 
     !then it's even. This intrinsic function modulo calculates remainder.
     
     if(modulo(Z, 2) == 0 .and. modulo(N,2) == 0) then 
        r = 1.0_dp 
        else if(modulo(Z, 2) > 0 .and. modulo(N,2) > 0) then
            r = -1.0_dp 
     end if 
     
     if(modulo(Z+N,2) > 0) then 
        r = 0.0_dp 
     end if
            
end function delta

!-----------------------------------------------------------------------
!! Subroutine: print_best_parameters
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Prints best fit parameters for semi empirical mass formula and it's 
!! corresponding uncertainties to screen.
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the best fit parameters
!! covariance       real        Array containing covariance matrix
!-----------------------------------------------------------------------
subroutine print_best_parameters(c_parameters, covariance)
    implicit none
    real(dp), intent(in) :: c_parameters(:), covariance(:,:)

  
    ! The covariance matrix is the inverse of the alpha matrix. 
    ! and the uncertainty squared of an individual linear parameter is the corresponding 
    ! diagonal of the covariance matrix. 
    ! So we print the square root of the corresponding diagonals.
    
    print *,
    print *, 'Best fit values:              value                 uncertainty'
    print 1, ' Volume parameter:   ', c_parameters(1) ,          sqrt(covariance(1,1)) 
    print 1, ' Surface parameter:  ', c_parameters(2) ,          sqrt(covariance(2,2))
    print 1, ' Asymmetry parameter:', c_parameters(3) ,          sqrt(covariance(3,3))    
    print 1, ' Coulomb parameter:  ', c_parameters(4) ,          sqrt(covariance(4,4)) 
    print 1, ' Pairing term:       ', c_parameters(5) ,          sqrt(covariance(5,5))

1 format(a,f15.8,e28.16)
end subroutine print_best_parameters



!-----------------------------------------------------------------------
!! function: semi_empirical_mass
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates binding energy for a given proton and neutron number
!!----------------------------------------------------------------------
!! Input:
!!
!! c    real        Array containing the parameters of the semi-empirical mass formula
!! Z    integer     number of protons in an isotope
!! N    integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r    real        Binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_mass(c, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: c(:)
    integer, intent(in) :: Z, N 
    real(dp), allocatable:: linear_terms(:) 
    allocate(linear_terms(1:size(c)))
    
    !Need an array to hold the linear terms multiplying the linear parameters
    call calculate_linear_termns(Z, N, linear_terms) 
   
    !BE = c_vol * A + c_surf * A^2/3 + c_asym * (N - Z)^2/A + c_coul * Z(Z-1)/A^1/3 + c_pair * A^-3/4 * delta(Z, N)
    
    r = c(1) * linear_terms(1) + c(2) * linear_terms(2) + c(3) * linear_terms(3) + c(4) * linear_terms(4) + c(5) * linear_terms(5)
    
    ! You can call the calculate_linear_termns subroutine and use its output
    ! to calculate the binding energy
end function semi_empirical_mass

!-----------------------------------------------------------------------
!! function: semi_empirical_error
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates theoretical error of the binding energies we calculated.
!!----------------------------------------------------------------------
!! Input:
!!
!! covariance   real        2D array containing the parameters' covariance matrix
!! Z            integer     number of protons in an isotope
!! N            integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        statistical uncertainty in the binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_error(covariance, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: covariance(:,:)
    integer, intent(in) :: Z, N
    real(dp), allocatable:: linear_terms(:) 
    integer :: shape_inverse(1:2), i, j  
    shape_inverse = shape(covariance) 
    !We need a linear terms array to hold the terms multiplying the linear parameters.
    allocate(linear_terms(1:shape_inverse(1))) 
    call calculate_linear_termns(Z, N, linear_terms) 
    r = 0.0_dp
    
    !The formula for the uncertainty squared is 
    !sum from i = 1 to number of parameters, sum from j = 1 to number of parameters 
    !partial dg/dx_i * partial dg/dx_j * C_ij 
    
    !C is the covariance matrix or inverse of alpha 
    !g is some function that depends on the parameters x. 
    
    !We can calculate it's derivative and say that the derivative of g 
    !with respect to each parameter are the elements of the linear term 
    !array that we have been using. 
    
    !So we don't have to calculate the derivative using another method 
    !Since it is already coded.
    do i = 1, size(linear_terms) 
        do j = 1, size(linear_terms) 
            r = r + linear_terms(i) * linear_terms(j) * covariance(i,j) 
            end do 
    end do 
    r = sqrt(r)
   
end function semi_empirical_error



!-----------------------------------------------------------------------
!! subroutine: stable_isotopes
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Finds the valley of stability.
!!----------------------------------------------------------------------
!! Input:
!!
!! n_neutrons   real    Array containing number of neutrons from experimental file
!! c            real    Array containing calculated parameters of semi empirical mass formula
!-----------------------------------------------------------------------
!! Output:
!!
!! stable_neutrons   real    Array containing position of valley of stability
!----------------------------------------------------------------------- 
subroutine stable_isotopes(binding_energies, lowest_energies, n_neutrons, stable_neutrons, c)  
implicit none 
real(dp), allocatable :: binding_energies(:), lowest_energies(:) 
real(dp) :: c(:), minimum 
integer, allocatable, intent(out) :: stable_neutrons(:) 
integer :: n_neutrons(:), nucleon, i, z, neutrons


allocate(stable_neutrons(1:118)) 



!A nucleon is defined as the number of protons Z + number of neutrons N. 
!We want to find the binding energy per nucleon so we divide calculated binding 
!energy by the nucleons. 
!The smallest binding energy has the largest negative. 

!Set up a do loop to go from proton numbers 1-118 
!and loop infinitely through neutron values and only
!exit once the binding energy per nucleon is greater than 
!the one previously calculated

!We will save the valley of stability position to an array 
!called stable_neutrons.
do z = 1, size(stable_neutrons) 
    neutrons = 1 
    do  
        nucleon = z + neutrons
        if(semi_empirical_mass(c, z, neutrons)/nucleon > semi_empirical_mass(c, z, neutrons - 1)/nucleon) then 
            exit
        end if 
        neutrons = neutrons + 1
        end do 
    stable_neutrons(z) = neutrons
end do 

end subroutine stable_isotopes 

!-----------------------------------------------------------------------
!! subroutine: find_dripline
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Finds the position of the neutron dripline
!!----------------------------------------------------------------------
!! Input:
!!
!! n_neutrons   real    Array containing number of neutrons from experimental file
!! c            real    Array containing calculated parameters of semi empirical mass formula
!-----------------------------------------------------------------------
!! Output:
!!
!! neutron_dripline   real    Array containing position of neutron dripline
!----------------------------------------------------------------------- 
subroutine find_dripline(separation, neutron_dripline, c, n_neutrons) 
implicit none 
real(dp), allocatable :: separation(:) 
real(dp) :: c(:) 
integer, allocatable, intent(out) :: neutron_dripline(:) 
integer :: n_neutrons(:), z, neutrons 



 
allocate(neutron_dripline(1:118)) 

!We want to find the maximum number of neutrons such that the separation energy stays positive 
!We loop through proton numbers 1-118 and all the neutron values and calculate separation energy 
!using the formula below and check if it's greater than 0 and the current neutron value is greater than 
!the largest neutron value for that proton number. Reset maximum number of neutrons each time we go to a new 
!proton value.

do z = 1, size(neutron_dripline) 
    neutrons = 1
    do 
        !Calculate separation energy for every number of neutrons and if the energy is greater than 0 
        !and if the neutron number is higher than before then that is the position of the neutron dripline 
        !at a given z. 
        
        !Separation Energy Formula: 
        !S_n = BE(Z, N-1) - BE(Z, N)
        
        
        if(semi_empirical_mass(c, z, neutrons - 1) - semi_empirical_mass(c, z, neutrons) < 0) then 
            exit
        end if 
        neutrons = neutrons + 1
       
    end do 
   neutron_dripline(z) = neutrons 
end do 
end subroutine find_dripline


end module nuclear_model
