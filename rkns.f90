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


! Integration functions and machinery mostly for Redfield
! propagation but Lindblad needs the integrals, too.

! Theta, G and Theta+/Theta- notation in this file follows the derivation in 
! Pollard, W. T. and Frisner, R. A. "Solutions of the Redfield equation for the dissipative quantum dynamics of 
! multilevel systems", J. Chem. Phys. 100, 1994.

! The derivation of the fast, secular Redfield as implemented here is based on
! Balzer, B. and Stock, G. "Modeling of decoherence and dissipation in nonadiabatic
! photoreactions by an effective-scaling nonsecular Redfield algorithm", Chem. Phys.
! 310, 2005.


MODULE rk_ns_utils
USE prequel
USE parameters
IMPLICIT NONE

CONTAINS

! Trace over the matrix sigma
! Start at diagonal index s, end at e, inclusive.
! Return the sum of the real parts traced over
COMPLEX*16 FUNCTION trace_mat_c(sigma, s, e)
  COMPLEX*16, DIMENSION(:,:), INTENT(IN) :: sigma
  INTEGER, INTENT(IN) :: s, e
  INTEGER i
  REAL*8 temp
  temp = 0.0d0
  DO i=s,e
    temp = temp + REAL(REAL(sigma(i,i)))
  END DO
  trace_mat_c = temp
END FUNCTION 

! Trace over the matrix sigma
! Start at diagonal index s, end at e, inclusive.
! return the sum of the real parts traced over
REAL*8 FUNCTION trace_mat(d,sigma, s, e)
  INTEGER, INTENT(IN) :: d
  COMPLEX*16, DIMENSION(d,d), INTENT(IN) :: sigma
  INTEGER, INTENT(IN) :: s, e
  INTEGER i
  REAL*8 temp
  temp = 0.0d0
  DO i=s,e
    temp = temp +  REAL(REAL(sigma(i,i)))
  END DO
  trace_mat = temp
END FUNCTION 

! Get matrix theta plus according to
! rule in paper (integration of correlation 
! function times factor from 0 to infinity)
! i is the index of the coupling operator
! \int_-\infty^\infty e^{-\omega_{i,j}t}*C(t) dt
! Where C is the correlation function
SUBROUTINE get_Theta_plus(th_pl,H,ti,tf)
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: th_pl  ! Matrix to be filled
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: H  ! Hamiltonian
  REAL*8, INTENT(IN) :: ti, tf  ! In the case of time dependent integrals (no longer used)
  REAL*8 :: omegaij  ! Frequency difference between relevant eigenvalues
  INTEGER :: i, j, k
  th_pl = (0.0d0,0.0d0)
  DO k = 1, nbath
    DO j=1,n
      DO i=1,n
        omegaij = REAL(REAL((H(i,i)-H(j,j))/hbar))
        th_pl(i,j,k) = get_thpl_element(omegaij,k,ti,tf)
      END DO
    END DO
  END DO
END SUBROUTINE

! Convenient helper since this stuff needs to be called in a bunch of
! different circumstances; returns a theta+ element given the frequency,
! the bath number (k) and initial/final times. Since time dependence
! is no longer used, ti and tf are no longer important
COMPLEX*16 FUNCTION get_thpl_element(omegaij,k,ti,tf)
  REAL*8, INTENT(IN) :: omegaij, ti, tf
  INTEGER, INTENT(IN) :: k
  ! Time independent: exact value available for real part;
  ! imaginary part from infinite fourier transform
  ! Special things occur at zero temperature; deal with them in another conditional
  IF (beta .GT. 0) THEN 
    get_thpl_element = get_theta_element_nz_ti(omegaij,k)
  ELSE IF (beta .LE. 0) THEN
    WRITE(*,*) "Error: zero temperature bath no longer supported."
    STOP
  END IF
  IF (ISNAN(REAL(REAL(get_thpl_element))) .OR. ISNAN(AIMAG(get_thpl_element))) THEN
    WRITE(*,*) "NaN in thpl ", omegaij, k, ti, tf
  END IF
END FUNCTION

! Calculate a theta+ element when the temperature is not zero
! and the calculation is time independent
COMPLEX*16 FUNCTION get_theta_element_nz_ti(omegaij,i)
  INTEGER, INTENT(IN) :: i
  REAL*8, INTENT(IN) :: omegaij  ! frequency to integrate for
  REAL (KIND=4) :: res, abserr, epsabs, epsrel
  INTEGER (KIND=4) :: neval, ier, sum_ier
  COMPLEX*16 temp
  temp = (0.0d0,0.0d0)
  temp = real_part_plus(omegaij,i)
  sum_ier = 0
  epsabs = 0.000
  epsrel = 0.0005

  CALL set_cbi(i)  ! Tell Quadpack which bath is referenced
  CALL set_omega_cbi(omegaij)  ! Tell quadpack what \omega is

  ! Different approaches must be used for different values
  IF (ABS(omegaij) .LT. 1d-14) THEN
    IF (bath_type(i) .EQ. 0) THEN ! Ohmic bath
      temp = temp + CMPLX(0.0d0,-eta(i)*omegac(i))
    ELSE IF (bath_type(i) .EQ. 1) THEN  ! Debye bath
      temp = temp + CMPLX(0.0d0,-pi*eta(i)/(2.0d0*omegac(i)))
    ELSE
      WRITE(*,*) "Warning, unsupported bath type"
      STOP
    END IF
  ELSE IF (omegaij .GT. 0) THEN
    ! Integrate up to the problematic region
    CALL qags(cpi_1_weighted,0.0,0.9*REAL(omegaij,4),epsabs,epsrel,res,abserr,neval,ier)
    temp = temp + CMPLX(0.0d0,res)
    ! Perform the integral which is actually going to have principal value issues
    CALL qawc(cpi_1,0.9*REAL(omegaij,4),1.1*REAL(omegaij,4),&
              REAL(omegaij,4),epsabs,epsrel,res,abserr,neval,ier)
    temp = temp + CMPLX(0.0d0,res)
    ! Finish off on that integral
    CALL qagi(cpi_1_weighted,1.1*REAL(omegaij,4),1,epsabs,epsrel,res,abserr,neval,ier)
    temp = temp + CMPLX(0.0d0,res)

    ! Perform the integral without principal value issues
    CALL qagi(cpi_2_weighted,0.0,1,epsabs,epsrel,res,abserr,neval,ier)
    temp = temp - CMPLX(0.0d0,res)
  ELSE
    ! As above, the integral that has the singularity has now swapped position
    CALL qagi(cpi_1_weighted,0.0,1,epsabs,epsrel,res,abserr,neval,ier)
    temp = temp + CMPLX(0.0d0,res)
    ! Integrate up to the problematic region
    CALL qags(cpi_2_weighted,0.0,-0.9*REAL(omegaij,4),epsabs,epsrel,res,abserr,neval,ier)
    temp = temp - CMPLX(0.0d0,res)
    ! Integrate the problematic region
    CALL qawc(cpi_2,-0.9*REAL(omegaij,4),-1.1*REAL(omegaij,4),&
              -REAL(omegaij,4),epsabs,epsrel,res,abserr,neval,ier)
    temp = temp - CMPLX(0.0d0,res)
    ! Finish the integral to infinity
    CALL qagi(cpi_2_weighted,-1.1*REAL(omegaij,4),1,epsabs,epsrel,res,abserr,neval,ier) 
    temp = temp - CMPLX(0.0d0,res)
  END IF

  get_theta_element_nz_ti = temp/pi   ! This factor of pi is a real nuisance

  IF (ISNAN(REAL(REAL(temp))) .OR. ISNAN(AIMAG(temp))) THEN
    WRITE(*,*) "NaN in Theta+: omegaij ", omegaij, "final value ", temp
    STOP
  END IF
END FUNCTION

! Since our coupling operator is presumed to be Hermitian,
! this is an acceptable formulation for theta-; if it weren't
! Hermitian the theta+ calculations would need to be redone with
! some tweaks
SUBROUTINE get_Theta_minus(th_mn,th_pl)
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: th_mn, th_pl
  INTEGER :: k
  DO k = 1, nbath
    th_mn(1:n,1:n,k) = TRANSPOSE(CONJG(th_pl(1:n,1:n,k)))
  END DO
END SUBROUTINE

! For the fast version of secular redfield, assemble R_(ii)_(jj) into a matrix
! and also (-\omega_ij + R_(ij)_(ij)) as a matrix as well. The former is
! the diag_mat which deals with eigenvalue population propogation. The later
! deals with the exponential decay of coherences.
SUBROUTINE secular_rf_setup(diag_mat,odiag_mat,theta_plus,theta_minus,H_d)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: diag_mat, odiag_mat  ! Diagonal and off diagonal propogation factors
  COMPLEX*16, DIMENSION(n), INTENT(IN) :: H_d  ! Diagonal hamilotnian (of energy eigenvalues)
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: theta_plus, theta_minus  ! FT of correlation functions
  REAL*8 :: temp
  INTEGER :: i, j
  diag_mat = 0.0d0
  odiag_mat = 0.0d0

  DO i = 1, n
    DO j = 1, n
      diag_mat(i,j) = tensor_element(i,i,j,j,theta_plus,theta_minus,G_global,G_global)
    END DO
  END DO

  DO i = 1, n
    DO j = 1, n
      IF (i .EQ. j) THEN
        odiag_mat(i,j) = 0.0d0
      ELSE
        temp = H_d(j) - H_d(i)
        odiag_mat(i,j) = (CMPLX(0.0d0,temp/hbar) +&
                        tensor_element(i,j,i,j,theta_plus,&
                        theta_minus,G_global,G_global))
        ! WRITE(*,*) "off diagonal matrix element is ", odiag_mat
        !WRITE(*,*) "Phase factor ", CMPLX(0.0d0,temp/hbar)
        ! WRITE(*,*) "Tensor element is ", tensor_element(i,j,i,j,theta_plus,&
        !                theta_minus,G_global,G_global)


      END IF
    END DO
  END DO
END SUBROUTINE


! Find the real part of the fourier transform from 0 to inf
! of the correlation function. This can be done exactly.
! For use when C(tau) is required
! Do NOT call this routine for zero temperature.
REAL*8 FUNCTION real_part_plus(omega_in,i)
  REAL*8, INTENT(IN) :: omega_in
  INTEGER, INTENT(IN) :: i
  REAL*8 :: omega, temp
  omega = omega_in
  real_part_plus = 0.0d0
  ! Deal with the case of almost exactly zero
  IF (ABS(omega) .LT. 1.0d-15) THEN
    IF (bath_type(i) .EQ. 0) THEN
      real_part_plus = eta(i)/beta
    ELSE IF (bath_type(i) .EQ. 1) THEN
      real_part_plus = eta(i)/(beta*omegac(i)**2)
    ELSE  ! Approximation for a bath type not implemented
      omega = 1.0d-14
      real_part_plus = EXP(-beta*hbar*omega)*spectral_density(omega,i)*hbar /&
                     (1.0d0 - EXP(-beta*hbar*omega))
    END IF
  ELSE  ! Nonzero
    temp = EXP(-beta*hbar*omega)
    IF (ABS(temp) .GT. HUGE(0.0d0)) THEN  ! It is infinity or negative infinity
      ! Limit as omega -> - inf is -1
      IF (omega .LT. 0.0d0) THEN
        real_part_plus = spectral_density(-omega,i)*hbar
      ! Limit as omega -> inf is 0
      ELSE
        real_part_plus = 0.0d0
      END IF
    ELSE  ! No problems with overflows. Proceed normally
      IF (omega .LT. 0.0d0) THEN
         real_part_plus = -EXP(-beta*hbar*omega)*spectral_density(-omega,i)*hbar /&
                       (1.0d0 - EXP(-beta*hbar*omega))
      ELSE
         real_part_plus = EXP(-beta*hbar*omega)*spectral_density(omega,i)*hbar /&
                       (1.0d0 - EXP(-beta*hbar*omega))
      END IF
    END IF
  END IF
  real_part_plus = real_part_plus*pi  ! Because everything is divided by pi later
  IF (ISNAN(real_part_plus)) THEN
    WRITE(*,*) "real_part_plus is NAN", omega_in, EXP(-beta*hbar*omega),spectral_density(-omega,i)*hbar,&
                (1.0d0 - EXP(-beta*hbar*omega)), -beta*hbar*omega

  END IF
END FUNCTION 

! Rk4 integrator for the revised secular equation propogation
! ds should be the vector of the diagonal components (populations) and
! diag_mat_ij = R_ii_jj contains the relevant Redfield
! factors
! Could have been done more efficiently but kept running
! into weird problems when the basis was changed to use
! exact exponentiation
SUBROUTINE drho_dt_secular(ds,diag_mat,dt)
  COMPLEX*16, DIMENSION(n,1), INTENT(INOUT) :: ds ! Diagonals to propogate
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: diag_mat  ! R_ii_jj
  REAL*8, INTENT(IN) :: dt  ! Time interval
  COMPLEX*16, DIMENSION(n,1) :: tempi, drdt, slope
  slope = MATMUL(diag_mat,ds)
  tempi = 0.5d0*slope*dt + ds
  drdt = (1.0d0/6.0d0)*slope*dt
  slope = MATMUL(diag_mat,tempi)
  tempi = ds + slope*dt/2.0d0
  drdt = drdt + (1.0d0/3.0d0)*slope*dt
  slope = MATMUL(diag_mat,tempi)
  drdt = drdt + (1.0d0/3.0d0)*slope*dt
  tempi = ds + slope*dt
  slope = MATMUL(diag_mat,tempi)
  drdt = drdt + (1.0d0/6.0d0)*slope*dt
  ds = ds + drdt
END SUBROUTINE

! Exponential integrator of the diagonal terms of the secular
! Redfield equation
SUBROUTINE advance_secular(ds,dt,evals,evec,evec_inv)
  COMPLEX*16, DIMENSION(n,1), INTENT(INOUT) :: ds ! Diagonals to propogate
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: evec, evec_inv
  COMPLEX*16, DIMENSION(n), INTENT(IN) :: evals
  REAL*8, INTENT(IN) :: dt  ! Time interval
  INTEGER :: i

  ds = MATMUL(evec_inv,ds)  
  DO i = 1, n
    ds(i,1) = EXP(evals(i)*dt)*ds(i,1)
  END DO
  ds = MATMUL(evec,ds)

END SUBROUTINE

! Build the redfield tensor, 4d matrix, according to Ga, Gb,
! theta minus etc. 
SUBROUTINE build_tensor(tensor,theta_plus,theta_minus,Ga,Gb)
  COMPLEX*16, DIMENSION(n,n,n,n), INTENT(OUT) :: tensor
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: theta_plus
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: theta_minus
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: Ga, Gb
  INTEGER :: i, j, k, l
  INTEGER :: td_switch
  td_switch = 0
  ! We already calculated this once; just fill it in
  IF (td_switch .EQ. 0 .AND. tensor_flag > 0) THEN
    tensor = tensor_global 
  ! Time dependent version; calculate from scratch
  ELSE
    DO i=1,n
      DO j=1,n
        DO k=1,n
          DO l=1,n
            tensor(i,j,k,l) = tensor_element(i,j,k,l,&
                              theta_plus,theta_minus,Ga,Gb)
            WRITE(*,*) "R ", i, j, k, l, " is ", tensor(i,j,k,l)
          END DO
        END DO
      END DO
    END DO
  END IF
  ! Set global tensor if it is only calculated once
  IF (td_switch .EQ. 0 .AND. tensor_flag .EQ. 0) THEN
    ALLOCATE(tensor_global(n,n,n,n))
    tensor_global = tensor
    tensor_flag = 1
  END IF
END SUBROUTINE

! Get the tensor element; helper function for building the redfield
! tensor
COMPLEX*16 FUNCTION tensor_element(i,j,k,l,theta_plus,theta_minus,Ga,Gb)
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: theta_plus
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: theta_minus
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: Ga, Gb
  INTEGER, INTENT(IN) :: i, j, k, l
  INTEGER :: r, bath_index 
  COMPLEX*16 :: temp
 
  temp = 0.0d0 
  DO bath_index = 1, nbath

    temp = temp + gamma_plus_element(l,j,i,k,theta_plus(1:n,1:n,bath_index),&
           Ga(1:n,1:n,bath_index),Gb(1:n,1:n,bath_index)) + &
           gamma_minus_element(l,j,i,k,theta_minus(1:n,1:n,bath_index),&
           Ga(1:n,1:n,bath_index),Gb(1:n,1:n,bath_index))
    IF (l .EQ. j) THEN
      DO r=1,n
        temp = temp - gamma_plus_element(i,r,r,k,theta_plus(1:n,1:n,bath_index),&
               Ga(1:n,1:n,bath_index),Gb(1:n,1:n,bath_index))
      END DO
    END IF
    IF (i .EQ. k) THEN
      DO r=1,n
        temp = temp - gamma_minus_element(l,r,r,j,theta_minus(1:n,1:n,bath_index),&
               Ga(1:n,1:n,bath_index),Gb(1:n,1:n,bath_index))
      END DO
    END IF

  END DO
!  WRITE(*,*) "Tensor element: ", tensor_element
  tensor_element = temp
END FUNCTION



! Calculate the redfield factor using the equations of
! Pollard and Friesner
! May want to switch to LAPACK routines at some point
SUBROUTINE redfield(sigma,factors,theta_plus,theta_minus,Ga,Gb,G_pls_sig,temp)
  ! redfield factors
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: factors
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: sigma
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: theta_plus, theta_minus
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: Ga, Gb
  COMPLEX*16, INTENT(INOUT), DIMENSION(n,n) :: G_pls_sig, temp  ! Both are temporaries
  COMPLEX*16, DIMENSION(:,:,:,:), ALLOCATABLE :: tensor
  INTEGER :: i, j, k, l

  ! Build the redfield tensor factor by factor and
  ! eliminate those required to make the secular approximation
  factors = (0.0d0,0.0d0)
!  WRITE(*,*) "CHECK 1 FACTORS ", factors
!  WRITE(*,*) "CHECK 1 SIGMA ", sigma
  IF (secular .EQ. 1) THEN
    ALLOCATE(tensor(n,n,n,n))
    tensor = (0.0d0,0.0d0)
    CALL build_tensor(tensor,theta_plus,theta_minus,Ga,Gb) 
!    WRITE(*,*) "Tensor ", tensor
    DO i=1,n
      DO j=1,n
        DO k=1,n
          DO l=1,n
            IF ((i .EQ. k .AND. j .EQ. l) .OR. (i .EQ. j .AND. k .EQ. l)) THEN
!              WRITE(*,*) "CHECK 2.5 factors ", factors(i,j), tensor(i,j,k,l),&
!                         sigma(k,l)
              factors(i,j) = factors(i,j) + tensor(i,j,k,l)*sigma(k,l)
            END IF
          END DO
        END DO
      END DO
    END DO 
    DEALLOCATE(tensor)
  ! Otherwise, perform cheaper matrix multiplication
!    WRITE(*,*) "CHECK 2 FACTORS ", factors
  ELSE
    !WRITE(*,*) "CHECK 1 FACTORS ", factors
    IF (SUM(g_min_flag) + SUM(g_flag) + SUM(g_pls_flag) .LT. 3*nbath) THEN
      WRITE(*,*) "Error. Unset global matrices."
      STOP
    END IF
    factors = 0.0d0
    DO i = 1, nbath
      G_pls_sig = MATMUL(G_pls_global(1:n,1:n,i),sigma) 
      temp = MATMUL(G_pls_sig,G_global(1:n,1:n,i)) - MATMUL(G_global(1:n,1:n,i),G_pls_sig)
      factors = factors + temp + CONJG(TRANSPOSE(temp))
      !WRITE(*,*) "CHECK 2 FACTORS ", factors
    END DO
  END IF
  !WRITE(*,*) "Redfield factor equal to: ", factors
END SUBROUTINE


! Perform rk4 on the sigma matrix; don't worry about
! efficiency for now.
SUBROUTINE rk4(sigma,H,theta_plus,theta_minus,Ga,Gb,G_pls_sig,temp_mat) 
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: sigma
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: G_pls_sig, temp_mat
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(INOUT) :: theta_plus, theta_minus
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: Ga, Gb
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: H

  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: k1, k2, k3, k4, temp, redfield_factor
  INTEGER :: i, j

  ALLOCATE(k1(n,n))
  ALLOCATE(k2(n,n))
  ALLOCATE(k3(n,n))
  ALLOCATE(k4(n,n))
  ALLOCATE(temp(n,n))
  ALLOCATE(redfield_factor(n,n))

  CALL redfield(sigma, redfield_factor, theta_plus, theta_minus, Ga, Gb, G_pls_sig, temp_mat)
  ! Calculate k1=dt*f(sigma)
  DO i=1,n
    DO j=1,n
      k1(j,i) = dt*row_func(j, i, sigma, H, redfield_factor)
    END DO
  END DO  
  ! Calculate k2=dt*f(sigma + k1/2)
  temp = sigma + 0.5d0*k1
  CALL redfield(temp, redfield_factor, theta_plus, theta_minus, Ga, Gb, G_pls_sig, temp_mat)
  DO i = 1,n
    DO j = 1,n
      k2(j,i) = dt*row_func(j, i, temp, H, redfield_factor)
    END DO
  END DO
  ! Calculate k3=dt*f(sigma + k2/2)
  temp = sigma + 0.5d0*k2
  ! As we calculate at middle points t + dt/2, update theta
  CALL redfield(temp, redfield_factor, theta_plus, theta_minus, Ga, Gb, G_pls_sig, temp_mat)
  DO i = 1,n
    DO j = 1,n
      k3(j,i) = dt*row_func(j, i, temp, H, redfield_factor)
    END DO
  END DO
  ! Calculate k4=dt*f(sigma + k3)
  temp = sigma + k3
  ! As we calculate at the end poit, t + dt, update theta
  CALL redfield(temp, redfield_factor, theta_plus, theta_minus, Ga, Gb, G_pls_sig, temp_mat)
  DO i = 1, n
    DO j = 1,n
      k4(j,i) = dt*row_func(j, i, temp, H, redfield_factor)
    END DO
  END DO

  ! Calculate forward step sigma + 1/6(k1 + 2(k2+k3) + k4)
  temp = sigma + (1.0d0/6.0d0)*(k1 + 2.0d0*(k2 + k3) + k4)
  !WRITE(*,*) "Current sigma: "
  !CALL print_mat(sigma)
  !WRITE(*,*) "Advancing by: "
  !CALL print_mat((1.0d0/6.0d0)*(k1 + 2.0d0*(k2 + k3) + k4))
  sigma = temp

  DEALLOCATE(k1)
  DEALLOCATE(k2)
  DEALLOCATE(k3)
  DEALLOCATE(k4)
  DEALLOCATE(temp)
  DEALLOCATE(redfield_factor)

END SUBROUTINE

! Helper function for differential; determines the value for
! the derivative of sigma at any point; requires the redfield
! factors for sigma to be computed already and passed in
COMPLEX*16 FUNCTION row_func(row, col, sigma, H, redfield_factor)
  ! The one we are after
  INTEGER, INTENT(IN) :: row, col
  ! Density matrix, hamiltonian, redfield factors
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: sigma, H, redfield_factor
  row_func = (0.0d0,0.0d0)
  row_func = -CMPLX(0.0d0,1.0d0,8)*sigma(row,col)*(H(row,row) - H(col,col))&
             /hbar + redfield_factor(row, col)
END FUNCTION row_func

END MODULE rk_ns_utils


