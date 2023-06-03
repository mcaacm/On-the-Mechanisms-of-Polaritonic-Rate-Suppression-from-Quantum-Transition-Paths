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

! Run a secular redfield calculation starting from each (pure)
! eigenstate and report on the percentage of state B character
! that results after a given period of time in order to assemble
! a Markov state model.
PROGRAM create_msm
USE parameters
USE prequel
USE rk_ns_utils
USE set_up_H
IMPLICIT NONE

! Density matrix, Hamiltonian which will be diagonalized
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: sigma, H
! Energy eigenvector inverse, diagonal and off diagonal matrix parts for secular propogation
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: evec_inv, diag_mat, odiag_mat, diag_sig
! Energy eigenvectors
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: evec
! Piece that is set up for jump dynamics even when that isn't run
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: translated_evecs
! Correlation function fourier transforms, coupling operator
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:,:) :: theta_plus, theta_minus, couple_op
! Right eigenvalues of Hamiltonian 
COMPLEX*16, DIMENSION(n) :: eval, H_d
REAL*8, DIMENSION(n) :: estate_probs, evals
REAL*8 :: normalization
! Loop variables
INTEGER :: i, j, k, l, stop_point
COMPLEX*16, DIMENSION(n) :: dm_evals
COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: dm_evecs, dm_evecs_inv, temp_diag_mat
REAL*8, DIMENSION(30) :: v_elems
! Real eigenvalues
REAL*8, DIMENSION(n) :: revals
! The Markov state model matrix
COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: markov_mat

! Allocate memory
ALLOCATE(sigma(n,n))
ALLOCATE(H(n,n))
ALLOCATE(evec_inv(n,m))
ALLOCATE(evec(m,n))
ALLOCATE(translated_evecs(n,n))
ALLOCATE(theta_plus(n,n,nbath))
ALLOCATE(theta_minus(n,n,nbath))
ALLOCATE(couple_op(n,n,nbath))
! Allocate global storage found in prequel
ALLOCATE(G_min_global(n,n,nbath))
ALLOCATE(G_pls_global(n,n,nbath))
ALLOCATE(G_global(n,n,nbath))
! Allocate global storage in Setup_H
ALLOCATE(Q_x_g(n,n))
ALLOCATE(Q_y_g(n,n))
! These are never deallocated but it doesn't matter;
! They're important until the end
ALLOCATE(markov_mat(n,n))


H = (0.0d0, 0.0d0)
evec = (0.0d0, 0.0d0)
theta_plus = (0.0d0, 0.0d0)
sigma = (0.0d0,0.0d0)
eval = (0.0d0,0.0d0)

WRITE(*,*) "Setting up Hamiltonian."
CALL setup(H,sigma,evec,eval,translated_evecs,estate_probs,couple_op)
! Copy and invert the eigenvector matrix
evec_inv = CONJG(TRANSPOSE(evec))

! Number of propagation steps to run; just one for an MSM probability calculation
stop_point = msm_steps 
DO i = 1, n
  H_d(i) = H(i,i)
END DO

ALLOCATE(diag_mat(n,n))
ALLOCATE(odiag_mat(n,n))
ALLOCATE(diag_sig(n,1))

CALL set_G(couple_op)
CALL get_Theta_plus(theta_plus, H, 0.0d0, 100.0d0)

CALL get_Theta_minus(theta_minus, theta_plus)
CALL set_G_plus(theta_plus)
CALL set_G_minus(theta_minus)

revals = REAL(REAL(eval))
WRITE(*,*) "Setting up secular Redfield machinery."
CALL secular_rf_setup(diag_mat,odiag_mat,theta_plus,theta_minus,H_d)
OPEN(89,FILE="markov_mat.txt")
WRITE(89,*) n, dt*stop_point
evals = REAL(REAL(eval))
evals = evals - evals(1)  ! Set first eval to 0; see if this helps
normalization = SUM(EXP(-beta*evals))
WRITE(89,*) EXP(-beta*evals)/normalization

ALLOCATE(dm_evecs(n,n))
ALLOCATE(dm_evecs_inv(n,n))
ALLOCATE(temp_diag_mat(n,n))
temp_diag_mat = diag_mat
CALL diagonalize_m_g(n,temp_diag_mat,dm_evals,dm_evecs)
dm_evecs_inv = dm_evecs
CALL invert(dm_evecs_inv)

! Quick print of some of the coupling matrix elements
OPEN(90,FILE="overlap_elements_1.txt")
DO i = 1, 30
  DO j = 1, 30
    v_elems(j) = ABS(couple_op(j,i,1)**2.0d0)
  END DO
  WRITE(90,*) v_elems
END DO
CLOSE(90)

WRITE(*,*) "Running."
! For every eigenstate, run one short secular initialization
DO l = 1, n
  sigma = 0.0d0
  sigma(l,l) = 1.0d0  
  DO k = 1, n
    diag_sig(k,1) = sigma(k,k)
  END DO

  DO i = 1, stop_point
    ! Updates the diagonal components
    CALL drho_dt_secular(diag_sig,diag_mat,dt)     
    ! Updates the non-diagonal components
    DO j = 1, n
      DO k = 1, n
        sigma(j,k) = sigma(j,k)*EXP(odiag_mat(j,k)*dt)
      END DO
    END DO
  END DO

  DO j = 1, n
    sigma(j,j) = diag_sig(j,1)
  END DO
  WRITE(89,*) REAL(REAL(diag_sig,8),8)  ! Already normalized; zero probability to transition to uncollapsed state
  WRITE(*,*) "Secular Redfield propagation for eigenstate ", l, "complete"
END DO

WRITE(89,*) evals  ! Useful to have the energy in the MSM file

CLOSE(89)

END PROGRAM create_msm

