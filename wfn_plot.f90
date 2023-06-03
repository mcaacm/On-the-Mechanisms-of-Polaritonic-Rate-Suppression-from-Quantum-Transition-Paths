! "On the mechanism of polaritonic rate suppression from quantum transition paths"
! Michelle C Anderson, Esmae J. Woods, Thomas P. Fay, David J. Wales and David T. Limmer
! University of California Berkeley
! Copyright 2023

! This file is part of TPT_SM 
! TPT_SM is free software: you can redistribute it and/or modify it under the terms of the GNU General 
! Public License as published by the Free Software Foundation, either version 3 of the License, 
! or (at your option) any later version.
! TPT_SM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with TPT_SM. 
! If not, see <https://www.gnu.org/licenses/>.

! Reads specifications in params.f90 and plots the wavefunctions of the
! eigenstates requested

PROGRAM plot_wf

USE parameters
USE prequel
USE set_up_H
USE rk_ns_utils

IMPLICIT NONE
INTEGER :: i ! Loop indexes
CHARACTER(LEN=100) :: fname   ! File name
CHARACTER(LEN=20) :: fstring  ! More file name stuff
REAL*8, DIMENSION(sz,sz) :: hf  ! Hermite polynomial factors
REAL*8, DIMENSION(:,:), ALLOCATABLE :: grid  ! Hermite evaluation on grid points
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: sigma, H, evec  ! Denisty matrix, Hamiltonian, energy eigenvectors
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:,:) :: couple_op  ! Coupling operator
REAL*8, DIMENSION(n) :: estate_probs  ! Unneeded setup
COMPLEX*16, DIMENSION(n) :: eval  ! Energy eigenvalues
COMPLEX*16, DIMENSION(m,1) :: wfn_m  ! Wavefunctions in full vs reduced basis
COMPLEX*16, DIMENSION(n,1) :: wfn
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: translated_evecs  ! Unneeded setup
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: plot  ! For wavefunction plotting
REAL*8 :: low, high, lowDVR, highDVR, integral  ! Plot grid limits and an integration temporary
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: out_grid  

IF (bound1 .LT. 0 .OR. bound2 .LT. 0 .OR. bound2 .LT. bound1) THEN
  WRITE(*,*) "Invalid plot bounds:", bound1, bound2
  STOP
END IF
IF (bound3 .LT. 0 .OR. bound4 .LT. 0 .OR. bound4 .LT. bound3) THEN
  WRITE(*,*) "Invalid plot bounds:", bound3, bound4
  STOP
END IF

ALLOCATE(plot(res,res,6))
ALLOCATE(grid(sz,res))
ALLOCATE(out_grid(res,2))

ALLOCATE(sigma(n,n))
ALLOCATE(H(n,n))
ALLOCATE(evec(m,n))
ALLOCATE(couple_op(n,n,nbath))
ALLOCATE(translated_evecs(n,n))
! Allocate global storage found in prequel
ALLOCATE(G_min_global(n,n,nbath))
ALLOCATE(G_pls_global(n,n,nbath))
ALLOCATE(G_global(n,n,nbath))
! Allocate global storage in Setup_H
ALLOCATE(Q_x_g(n,n))
ALLOCATE(Q_y_g(n,n))


! Setup system
WRITE(*,*) "Setting up system Hamiltonian"
CALL setup(H,sigma,evec,eval,translated_evecs,estate_probs,couple_op)

! Fill a grid of Hermite polynomials for the basis
low = -70.0d0
high = 70.0d0
lowDVR = -pi
highDVR = pi
CALL fill_hermite(sz,hf)
CALL fill_wf_grid(sz,res,hf,grid,low,high)
WRITE(*,*) "Filled hermite grid"


! Plot lower energy eigenstates of relevance
DO i = bound1, bound2
  WRITE(*,*) "Plotting state ", i
  wfn = 0.0d0
  wfn(i,1) = 1.0d0
  wfn_m = MATMUL(evec,wfn)
  ! Fill and integrate a grid for each eigenstate
  CALL fill_wf_plot_ho_sinc(sz,res,grid,low,high,lowDVR,highDVR,wfn_m,plot(1:res,1:res,1:3))
  integral = integrate_grid_ho_sinc(res,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR)
  WRITE(*,*) "Integrated density eigen fn ", i, "is", integral

  IF (i .LT. 10) THEN
    fstring = "(A21,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A21,I2,A4)"
  ELSE 
    fstring = "(A21,I3,A4)"
  END IF 
   

  WRITE(fname,fstring) "data/wfn_plots/ef_2d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  CALL print_wf_plot_ho_sinc(43,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR,res)
  CLOSE(43)
END DO


! Check normalization of individual eigenstates
! Plot higher energy eigenstates of relevance
DO i = bound3, bound4
  WRITE(*,*) "Plotting state ", i
  wfn = 0.0d0
  wfn(i,1) = 1.0d0
  wfn_m = MATMUL(evec,wfn)
  CALL fill_wf_plot_ho_sinc(sz,res,grid,low,high,lowDVR,highDVR,wfn_m,plot(1:res,1:res,1:3))
  integral = integrate_grid_ho_sinc(res,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR)
  WRITE(*,*) "Integrated density eigen fn ", i, "is", integral

  IF (i .LT. 10) THEN
    fstring = "(A21,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A21,I2,A4)"
  ELSE 
    fstring = "(A21,I3,A4)"
  END IF 
   
  WRITE(fname,fstring) "data/wfn_plots/ef_2d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  CALL print_wf_plot_ho_sinc(43,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR,res)
  CLOSE(43)

END DO


CLOSE(11)


END PROGRAM

