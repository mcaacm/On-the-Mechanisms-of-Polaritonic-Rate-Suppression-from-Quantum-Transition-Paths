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

! Sets up the Hamiltonian, initial vector to propagate from in
! Lindblad, coupling operators, etc.

! SM setup follows Li, X., Mandal, A., Huo, P. "Theory of Mode-Selective
! Chemistry through Polaritonic Vibrational Strong Coupling" J. Phys. Chem. Lett.
! 2021, 12, 29, 6974-6982.

MODULE set_up_H
USE parameters
USE prequel
USE rk_ns_utils

IMPLICIT NONE

CONTAINS 

! Output information about a Shin Menitu model
! gives position and energy because diabatic states are not a thing
SUBROUTINE output_info_SM_dmat(time,fd,dmat,H)
  INTEGER, INTENT(IN) :: fd
  COMPLEX*16, INTENT(IN), DIMENSION(n,n) :: dmat  ! Density matrix
  COMPLEX*16, INTENT(IN), DIMENSION(n,n) :: H
  REAL*8, INTENT(IN) :: time
  REAL*8 :: Q_x, Q_y, H_temp
  REAL*8, DIMENSION(2*n) :: output_wfn
  INTEGER :: i

  ! Usual wavefunction format for plotting via wfn_plot, even
  ! though it's a mixed state it should be fine
  DO i = 1, n
    output_wfn(2*i - 1) = dmat(i,i)
    output_wfn(2*i) = 0.0d0
  END DO

  Q_x = trace_mat(n,MATMUL(Q_x_g,dmat),1,n)
  Q_y = trace_mat(n,MATMUL(Q_y_g,dmat),1,n)
  H_temp = trace_mat(n,MATMUL(H,dmat),1,n)
  WRITE(fd,*) time, Q_x, Q_y, H_temp, output_wfn
END SUBROUTINE




! Setup a DVR sinc basis kinetic energy matrix, 1d
! x_i = i \dx*i
SUBROUTINE sinc_ke(dx,dim_coor,mass,kemat)
  REAL*8, INTENT(IN) :: dx  ! delta x interval
  REAL*8, INTENT(IN) :: mass ! mass
  COMPLEX*16, DIMENSION(dim_coor,dim_coor), INTENT(OUT) :: kemat 
  INTEGER, INTENT(IN) :: dim_coor

  INTEGER :: i, j, lim
  REAL*8 :: temp1, temp2
  temp1 = (hbar**2.0d0)/(2.0d0*mass*(dx**2.0d0))
  temp2 = (pi**2.0d0)/3.0d0

  lim = FLOOR(dim_coor/2.0d0)
  kemat = 0.0d0
  DO i = -lim, lim
    DO j = -lim, lim
      IF (i .EQ. j) THEN
        kemat(i+lim+1,j+lim+1) = temp1*temp2*((-1.0d0)**(i-j))
      ELSE
        kemat(i+lim+1,j+lim+1) = temp1*((-1.0d0)**(i-j))*2.0d0/((i-j)**2.0d0)
      END IF
    END DO 
  END DO 

END SUBROUTINE



! Helper function for the DVR sinc basis; just returns the value
! of the full potential at the particular point in question, i*dx = R
! This is the part that will be tensor producted with the identity
! aR^2 + bR^4 + c + omegasc^2/2.0d0*(SQRT(2/hbar*omegasc^3)*xi*(mu(R))^2
COMPLEX*16 FUNCTION sinc_V(i,dx)
  INTEGER, INTENT(IN) :: i
  REAL*8, INTENT(IN) :: dx
  REAL*8 :: R, R2, R4, tanhR  ! position and position squared
  REAL*8 :: prefactor, temp

  R = dx*i
  R2 = R**2.0d0
  R4 = R2**2.0d0
  tanhR = TANH(coeff_y*R)
  temp = (SQRT(2.0d0/(hbar*omegasc**3.0d0))*coeff_xi*(tanhR*coeff_v + coeff_z*R))**2.0d0
  prefactor = (omegasc**2.0d0)/2.0d0
  sinc_V = prefactor*temp + coeff_a*R2 + coeff_b*R4 + coeff_c + coeff_lp*R
END FUNCTION

! Helper function for the DVR sinc basis; returns the value
! of the dipole at the particular point in question
! This is the part that will be tensor producted with q_c
! Includes prefacotrs!
COMPLEX*16 FUNCTION sinc_V_cross(i,dx)
  INTEGER, INTENT(IN) :: i
  REAL*8, INTENT(IN) :: dx
  REAL*8 :: loc

  loc = dx*i
  sinc_V_cross = coeff_z*loc + coeff_v*TANH(coeff_y*loc)
  sinc_V_cross = sinc_V_cross*omegasc**2.0d0*coeff_xi*&
    SQRT(2.0d0/(hbar*omegasc**3.0d0))
END FUNCTION

! Helper function for the DVR sinc basis; returns the value
! of the lamb shift \omegac*eta*R^2/Mpi 
COMPLEX*16 FUNCTION sinc_LS(i,dx)
  INTEGER, INTENT(IN) :: i
  REAL*8, INTENT(IN) :: dx
  REAL*8 :: loc

  loc = dx*i
  sinc_LS = (loc**2.0d0)*omegac(1)*eta(1)/(mass_R*pi)
END FUNCTION



! Call the setup routines for specific types; when multiple different models were supported,
! this had a much bigger role
SUBROUTINE setup(H, sigma, evec, eval, H_trans, estate_probs,couple_op,Qt_full)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: H, sigma  ! Hamiltonian, starting denisty matrix
  COMPLEX*16, DIMENSION(m,n), INTENT(OUT) :: evec  ! Energy eigenvectors
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: couple_op  ! System Bath coupling operator
  COMPLEX*16, DIMENSION(m,m), INTENT(OUT), OPTIONAL :: Qt_full  ! Full basis position operator
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval  ! Energy eigenvalues
  ! Boltzman likelihood to start in oscillator eigenstate; used by thermal init and DEPRECATED
  REAL*8, INTENT(OUT), DIMENSION(n) :: estate_probs
  ! Whose translation to the unshifted oscillator basis
  ! looks like this; each vector is its own column; used by thermal init and DEPRECATED
  COMPLEX*16, INTENT(OUT), DIMENSION(n,n) :: H_trans

  IF (run_type .EQ. SM) THEN 
    IF (PRESENT(Qt_full)) THEN
      CALL setup_H_shin_metiu(H,sigma,evec,eval,H_trans,estate_probs,couple_op,Qt_full)
    ELSE
      CALL setup_H_shin_metiu(H,sigma,evec,eval,H_trans,estate_probs,couple_op)
    END IF
  ELSE
    WRITE(*,*) "Unrecognized runtype"
    STOP
  END IF
END SUBROUTINE


! Set a Shin-Metiu Hamiltonian describing an adiabatic setup
! where a proton/electron pair is confined in an optical cavity
! with two charges resulting in a quartic PES with the dipole
! moment serving as the coupling between the proton and photon modes
SUBROUTINE setup_H_shin_metiu(H,sigma,evec,eval,H_trans,estate_probs,couple_op,Qt_full)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: H, sigma  ! Hamiltonian, density matrix
  COMPLEX*16, DIMENSION(m,n), INTENT(OUT) :: evec  ! H eigenvectors
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: couple_op  ! System-bath coupling operator
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval  ! H evals
  ! Boltzman likelihood to start in oscillator eigenstate
  REAL*8, INTENT(OUT), DIMENSION(n) :: estate_probs  ! Initial probabilities for "vertical" excitation
  COMPLEX*16, OPTIONAL, INTENT(OUT), DIMENSION(m,m) :: Qt_full  ! Full H prior to basis trimming
  ! Translation of initial eigenvectors for "vertical excitation to the unshifted oscillator basis
  ! looks like this; each vector is its own column
  COMPLEX*16, INTENT(OUT), DIMENSION(n,n) :: H_trans
  ! Raising/lowering, q, q^2, p, p^2 operators for the photon mode
  ! Potential, matter-light coupling component for matter, kinetic energy and lambshift for the
  ! proton mode
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: raise_op, lower_op, q_op, q2_op, mom_op, &
    mom2_op, R_V, R_Vc, R_ke, R_ls
  ! Full size hamiltonian, temporaries and identity matrices for performing tensor products 
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: H_m, temp1, temp2, ident_qc, ident_R
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: evec_m, evec_inv
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: evals_m
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: evals_m_2d
  INTEGER :: i, j, sub
  TYPE(wfunc) :: wfn
  ! Full size density matrix and vertical excitation start vectors
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: sigma_m, H_trans_m
  COMPLEX*16 :: temp
  REAL*8 :: entanglement

  ALLOCATE(raise_op(dim_qc,dim_qc))
  ALLOCATE(lower_op(dim_qc,dim_qc))
  ALLOCATE(q_op(dim_qc,dim_qc))
  ALLOCATE(q2_op(dim_qc,dim_qc))
  ALLOCATE(mom_op(dim_qc,dim_qc))
  ALLOCATE(mom2_op(dim_qc,dim_qc))
  ALLOCATE(R_Vc(dim_R,dim_R))
  ALLOCATE(R_V(dim_R,dim_R))
  ALLOCATE(R_ke(dim_R,dim_R))
  ALLOCATE(R_ls(dim_R,dim_R))
  ALLOCATE(H_trans_m(m,m))
  ALLOCATE(H_m(m,m))
  ALLOCATE(sigma_m(m,m))

  CALL fill_raise(dim_qc,raise_op)
  CALL fill_lower(dim_qc,lower_op)

  q_op = SQRT(hbar/(2.0d0*omegasc))*(raise_op + lower_op)
  mom_op = SQRT(hbar*omegasc/2.0d0)*(raise_op - lower_op)*CMPLX(0.0d0,1.0d0)
  q2_op = MATMUL(q_op,q_op)
  mom2_op = MATMUL(mom_op,mom_op)

  ! Set up R coordinate sinc function basis ke and V matrices
  R_V = 0.0d0
  R_Vc = 0.0d0
  R_ls = 0.0d0
  CALL sinc_ke(dx_R,dim_R,mass_R,R_ke)
  sub = CEILING(1.0d0*dim_R/2.0d0)
  DO i = 1, dim_R
    R_V(i,i) = sinc_V(i - sub,dx_R)
    R_Vc(i,i) = sinc_V_cross(i - sub,dx_R)
    R_ls(i,i) = sinc_LS(i - sub,dx_R)
  END DO

  ALLOCATE(evals_m(m))
  ALLOCATE(evals_m_2d(m,1))
  ALLOCATE(evec_m(m,m))
  ALLOCATE(evec_inv(n,m))

  ALLOCATE(temp1(m,m))
  ALLOCATE(temp2(m,m))
  ALLOCATE(ident_qc(1:dim_qc,1:dim_qc))
  ALLOCATE(ident_R(1:dim_R,1:dim_R))  ! m = dim_qc*dim_R


  ident_qc = 0.0d0
  ident_R = 0.0d0
  DO i = 1, dim_qc
    ident_qc(i,i) = 1.0d0
  END DO 
  DO i = 1, dim_R
    ident_R(i,i) = 1.0d0
  END DO

  ! Add in the kinetic and potential term from the R sinc basis
  CALL tensor_product(dim_qc,dim_R,ident_qc,R_V + R_ke + R_ls,temp1)
  H_m = temp1 
  ! Add in the kinetic and potential term from the qc HO basis
  CALL tensor_product(dim_qc,dim_R,mom2_op/2.0d0 + 0.5d0*omegasc**2.0d0*q2_op,ident_R,temp1)
  H_m = H_m + temp1
  ! Add in the coupling term between R and qc
  CALL tensor_product(dim_qc,dim_R,q_op,R_Vc,temp1)
  H_m = H_m + temp1

  IF (RESTART .EQV. .FALSE.) THEN
    CALL diagonalize_m(m,H_m,evals_m,evec_m)
    evals_m_2d(1:m,1) = evals_m(1:m)
    OPEN(22,FILE="restart.dat")
    CALL print_mat_C16(evals_m_2d,22)
    CALL print_mat_C16(evec_m,22)
    CLOSE(22) 
  ELSE
    OPEN(22,FILE="restart.dat")
    i = read_matrix(evals_m_2d,22)
    i = i + read_matrix(evec_m,22)
    evals_m(1:m) = evals_m_2d(1:m,1)
    IF (i .NE. 2) THEN
      WRITE(*,*) "Failed to read in eigenvalues/vectors for restart"
      STOP
    END IF 
    CLOSE(22) 
  END IF
  
  evec = evec_m(1:m,1:n)
  eval = evals_m(1:n)
  H = 0.0d0
  DO i = 1, n
    H(i,i) = eval(i)
  END DO
  evec_inv = CONJG(TRANSPOSE(evec))
  ! The coupling operator is linear in R?
  ! Single bath system I believe
  DO i = 1, dim_R
    R_V(i,i) = dx_R*(i - sub)
  END DO
  CALL tensor_product(dim_qc,dim_R,ident_qc,R_V,temp1)
  IF (PRESENT(Qt_full)) THEN
    Qt_full = temp1
  END IF
  CALL to_basis(temp1,Q_x_g,evec,evec_inv)
  CALL tensor_product(dim_qc,dim_R,q_op,ident_R,temp1)
  CALL to_basis(temp1,Q_y_g,evec,evec_inv)


  ! The specified init_es is presumed to be an eigenstate
  IF (init_stat_therm .EQ. 1) THEN

    sigma = 0.0d0
    H_trans = 0.0d0
    estate_probs = 0.0d0
    sigma(init_es,init_es) = 1.0d0
    estate_probs = 1.0d0
    estate_probs(init_es) = 1.0d0
    H_trans(init_es,1) = 1.0d0

  ELSE  ! Use the init_es as a DVR position index and diagonalize
    estate_probs = 0.0d0
    estate_probs(1) = 1.0d0
    ident_R = 0.0d0  ! Reused as a temporary
    ident_qc = 0.0d0  ! Reused as a temporary
    ident_qc(1,1) = 1.0d0  ! Lowest HO state
    ident_R(init_es,init_es) = 1.0d0  ! Selected DVR center state
    CALL tensor_product(dim_qc,dim_R,ident_qc,ident_R,sigma_m)  ! Full proof but simpler ways exist
    CALL to_basis(sigma_m,sigma,evec,evec_inv)  ! Complete formation of sigma
    H_trans_m = 0.0d0
    H_trans_m(dim_qc*(init_es - 1) + 1, 1) = 1.0d0
    H_trans(1:n,1) = MATMUL(evec_inv,H_trans_m(1:m,1))
    OPEN(32,FILE="init_wfn.txt")
    WRITE(32,*) org_out(H_trans(1:n,1))
    CLOSE(32)
    temp = 0.0d0
    DO i = 1, n
      temp = temp + sigma(i,i)
    END DO
    sigma = sigma/temp  ! Normalize the trace
  END IF


  couple_op(1:n,1:n,1) = Q_x_g

  ! Write out some wavefunction information
  OPEN(63,FILE="plot_ev.txt")
  OPEN(64,FILE="evic.txt")
  ! Write information about the eigenvectors
  CALL init_wfn(wfn)
  DO i = 1, n
    entanglement = 0.0d0
    DO j = 1, m
      entanglement = entanglement - evec(j,i)**2.0d0*LOG(evec(j,i)**2.0d0)
    END DO
    wfn%fe = 0.0d0
    wfn%fe(i,1) = 1.0d0
    WRITE(64,*) i, op_avg(wfn,Q_x_g), op_avg(wfn,Q_y_g), op_avg(wfn,H), entanglement
    WRITE(63,*) i, -100.0, op_avg(wfn,H)
    WRITE(63,*) i, 100.0, op_avg(wfn,H)
    WRITE(63,*) ""
    WRITE(63,*) ""
  END DO
  CALL destroy_wfn(wfn)
  CLOSE(62)
  CLOSE(63)
  CLOSE(64)



  DEALLOCATE(H_m)
  DEALLOCATE(evec_m)
  DEALLOCATE(evals_m)
  DEALLOCATE(evec_inv)
  DEALLOCATE(R_ls)
  DEALLOCATE(raise_op)
  DEALLOCATE(lower_op)
  DEALLOCATE(q_op)
  DEALLOCATE(q2_op)
  DEALLOCATE(mom_op)
  DEALLOCATE(mom2_op)
  DEALLOCATE(R_Vc)
  DEALLOCATE(R_V)
  DEALLOCATE(R_ke)
  DEALLOCATE(H_trans_m)
  DEALLOCATE(sigma_m)

END SUBROUTINE


END MODULE



