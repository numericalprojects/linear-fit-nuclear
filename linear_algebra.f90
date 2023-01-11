!-----------------------------------------------------------------------
!Module: linear_algebra
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This is the linear algebra module that solves the matrix equation 
!! ax = b. The arrays are tested to make sure they are of correct size. 
!! Matrix A is inverted using LU decomposition and back substitution 
!! then this A inverse matrix is multiplied by the b vector to get us 
!! the x vector which has the best fit parameters.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_linear_system 
!! test_array_sizes 
!! invert_matrix 
!! matrix_multiply 
!! ludcmp 
!! lubksb
!!----------------------------------------------------------------------
module linear_algebra
use types
implicit none
private
public :: solve_linear_system, invert_matrix
contains

!-----------------------------------------------------------------------
!! Subroutine: solve_linear_system
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! This subroutine will solve the matrix equation Ax = b 
!! This subroutine calls another subroutine to test that the arrays are 
!! the right size for our purposes. It then calls another subroutine to 
!! invert matrix A. Then it multiplies vector B and inverted matrix A 
!! to get us the best fit parameters vector X. 
!!----------------------------------------------------------------------
!! Input:
!!
!! a_matrix        real        2D array containing the $a$ matrix
!! b_vector         real        1D array containing the $b$ vector
!-----------------------------------------------------------------------
!! Output:
!!
!! x_vector         real        1D array with the solution to the system of equations
!! a_inverse        real        2D array with the inverse matrix $a^{-1}$
!-----------------------------------------------------------------------
subroutine solve_linear_system(a_matrix, b_vector, x_vector, a_inverse)
    implicit none
    real(dp), intent(in) :: a_matrix(:,:), b_vector(:)
    real(dp), intent(out) ::  x_vector(:), a_inverse(:,:)
    ! The first thing is to make sure that all the arrays have proper sizes to
    ! make sure that  the matrix and vectors provided actually represent a
    ! system of equations. Otherwise the subroutines below will give errors or
    ! worst won't behave as expected but we won't notice.
    call test_array_sizes(a_matrix, b_vector, a_inverse, x_vector)
    call invert_matrix(a_matrix, a_inverse)

    call matrix_multiply(a_inverse, b_vector, x_vector) 
end subroutine solve_linear_system 

!-----------------------------------------------------------------------
!! Subroutine: matrix_multiply
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Solves for vector x which has the best fit parameters by multiplying 
!! the b vector by the inverse of matrix A using intrinsic function 
!! matmul.
!!----------------------------------------------------------------------
!! Input:
!!
!! a_inverse        real        2D array containing the inverse of alpha matrix
!! b_vector         real        1D array containing the $b$ vector
!-----------------------------------------------------------------------
!! Output:
!!
!! x_vector         real        1D array with the solution to the system of equations
!-----------------------------------------------------------------------
subroutine matrix_multiply (a_inverse, b_vector, x_vector)
    implicit none
    real(dp), intent(in) :: a_inverse(:,:), b_vector(:)
    real(dp), intent(out) ::  x_vector(:) 
    integer:: i, j 
    !We could write a whole subroutine for the matrix multiplication but 
    !there is an intrinsic GNU Fortran function matmul that multiplies 2 matrices
    x_vector = matmul(b_vector, a_inverse)
end subroutine matrix_multiply 

!-----------------------------------------------------------------------
!! Subroutine: test_array_sizes
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Checks if the arrays are the right size and ends the program if it's not.
!!----------------------------------------------------------------------
!! Input:
!!
!! a_matrix        real        2D array containing the $a$ matrix
!! b_vector         real        1D array containing the $b$ vector
!! x_vector         real        1D array with the solution to the system of equations
!! a_inverse        real        2D array with the inverse matrix $a^{-1}$
!-----------------------------------------------------------------------
subroutine test_array_sizes(a_matrix, b_vector, a_inverse, x_vector)
    implicit none
    real(dp), intent(in) :: a_matrix(:,:), b_vector(:), a_inverse(:,:), x_vector(:)

    integer :: shape_a(1:2), b_size, shape_ainv(1:2), x_size
    
    !Variables to hold columns and rows of A and its inverse 
    shape_a = shape(a_matrix)
    shape_ainv = shape(a_inverse)
    !Variables to hold elements in b vector and x vector
    b_size = size(b_vector)
    x_size = size(x_vector)
    
    !Checks if A and it's inverse are square and of the same size 
    !Checks if the number of columns in matrix are equal to the number 
    !of elements in the x vector and if the number of rows in matrix A 
    !are equal to the number of elements in the b vector.
    
    if(shape_a(1) /= shape_a(2) .or. shape_ainv(1) /= shape_ainv(2)) then 
        print *, "Matrix A or it's inverse are not square. Please adjust" 
        print *, "The program so that they are."
        stop
    end if
    if(shape_a(1) /= shape_ainv(1) .or. shape_a(2) /= shape_ainv(2)) then 
            print *, "Matrix A and it's inverse are not the same size." 
            print *, "Please adjust the program so that they are." 
            stop
    end if 
    if(shape_a(1) /= x_size .or. shape_a(2) /= b_size) then
             print *, "Matrix A has more columns than X has elements or" 
             print *, "Matrix A has more rows than B has elements."
             stop 
                
    end if 
   
end subroutine test_array_sizes

!-----------------------------------------------------------------------
!! Subroutine: invert_matrix
!-----------------------------------------------------------------------
!! Nathan Crawford
!!
!! Given a non singular matrix $a$, returns its inverse $a^{-1}$
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real    2D array containing the $a$ matrix
!!----------------------------------------------------------------------
!! Output:
!!
!! a_inv    real    2D array with the $a^{-1}$ matrix
!-----------------------------------------------------------------------
subroutine invert_matrix(a, a_inv)
    implicit none
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: a_inv(:,:)
    real(dp), allocatable :: a_work(:,:)
    integer :: shape_a(1:2), n, i
    real(dp) :: d
    integer, allocatable :: indx(:)
    
    allocate(a_work,mold=a)
    shape_a = shape(a)
    n = shape_a(1)
    allocate(indx(1:n))
    
    ! ludcmp destroys the input matrix a. In order to preserve a we will copy
    ! it into a work array that will be used in ludcmp
    a_work = a
    call ludcmp(a_work,indx,d)
    
    ! We construct a matrix that has orthogonal unit vectors as columns
    a_inv = 0._dp
    do i=1,n
        a_inv(i,i) = 1._dp
    enddo

    ! And then feed each column to the back-substitution routine
    do i = 1,n
        call lubksb(a_work,indx,a_inv(:,i))
    enddo
    ! This results in a_inv being the inverse of a
end subroutine invert_matrix



!-----------------------------------------------------------------------
!! Subroutine: ludcmp
!-----------------------------------------------------------------------
!! Nathan Crawford
!!
!! 
!! Performs LU decomposition on a non singular matrix $a$.
!! The original $a$ matrix is destroyed as the LU decomposition is returned
!! in the same array
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real        2D array containing the $a$ matrix
!!----------------------------------------------------------------------
!! Output:
!!
!! a        real        2D array with LU decomposition of the $a$ matrix
!! indx     integer     1D array that records the row permutation effected by the partial pivoting
!! d        real        +1 or -1 depending on whether the number of row interchanges was even or odd, respectively
!-----------------------------------------------------------------------
subroutine ludcmp(a, indx, d)
    implicit none
    real(dp), intent(inout) :: a(:,:)
    integer, intent(out) :: indx(:)
    real(dp), intent(out) :: d
    integer :: n,i,imax,j,k
    real(dp) aamax,dum,sum
    real(dp), allocatable :: vv(:)
    n = size(indx)
    allocate(vv(1:n))
    d=1._dp
    do i=1,n
        aamax=0._dp
        do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0._dp) then
            print *, 'singular matrix in ludcmp'
            stop
        endif
        vv(i)=1._dp/aamax
    enddo

    do j=1,n
        do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
        enddo
        aamax=0._dp
        do i=j,n
            sum=a(i,j)
            do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
                imax=i
                aamax=dum
            endif
        enddo
        if (j.ne.imax)then
            do k=1,n
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
            enddo
            d=-d
            vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0._dp) a(j,j) = tiny(1._sp)
        if(j.ne.n)then
            dum=1._dp/a(j,j)
            do i=j+1,n
                a(i,j)=a(i,j)*dum
            enddo
        endif
    enddo
end subroutine ludcmp

!-----------------------------------------------------------------------
!! Subroutine: lubksb
!-----------------------------------------------------------------------
!! Nathan Crawford
!!
!!
!! Performs back-substitution after a LU decomposition in order to solve the
!! linear system of equations $a \cdot x = b$. The $b$ vector is given in the b
!! array (which is destroyed) and the solution $x$ is returned in its place
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real        2D array containing the LU decomposition $a$ matrix (as returned by ludecomp)
!! indx     integer     1D array with the record of the row permutation effected by the partial pivoting (as returned by ludecomp)
!! b        real        1D array containing the $b$ vector
!!----------------------------------------------------------------------
!! Output:
!! b        real        1D array containing the $x$ vector
!-----------------------------------------------------------------------
subroutine lubksb(a, indx, b)
    implicit none
    real(dp), intent(in) :: a(:,:)
    integer, intent(in) :: indx(:)
    real(dp), intent(inout) :: b(:)

    integer :: n
    integer :: i,ii,j,ll
    real(dp) :: sum

    n = size(b)
    ii=0
    do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            enddo
        else if (sum.ne.0.) then
            ii=i
        endif
        b(i)=sum
    enddo
    do i=n,1,-1
        sum=b(i)
        do j=i+1,n
            sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
    enddo
end subroutine lubksb

end module linear_algebra
