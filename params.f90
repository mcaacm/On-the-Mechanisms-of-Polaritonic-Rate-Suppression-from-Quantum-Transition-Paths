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


! All parameters should be speicifed in atomic units to avoid any
! complications.These are all parameters to specify how the run proceeds.


MODULE parameters
IMPLICIT NONE

! In all cases save where otherwise noted, 1 is a 'true' or 'yes' on a binary
! switch and '0' is 'no' or 'false' but usually any number except 1 is ignored.

! Option to print extra information to deal with problems
LOGICAL, PARAMETER :: DEBUG = .FALSE.
! Read eigenvector/value info from restart.dat rather than diagonalizing H
LOGICAL, PARAMETER :: RESTART = .FALSE.

! This used to have lots of numbers for different model systems
! For the case of a shin_metiu model
INTEGER, PARAMETER :: SM = 1
INTEGER, PARAMETER :: secular = 1  ! Do secular Redfield

! What kind of run to perform (SM is the only setup included in this code, so leave unchanged)
INTEGER, PARAMETER :: run_type = SM

! When creating a markov state model, how many timesteps to run (tau = dt*msm_steps)
INTEGER, PARAMETER :: msm_steps = 500

! Which eigenstates LB propogation should end with
INTEGER, PARAMETER :: flag_estate_1 = 1  ! First eigenstate to stop propogation (should be 1)
INTEGER, PARAMETER :: flag_estate_2 = 4  ! Second eigenstate to stop propogation (should not be 1)

! Number of trajectories for a lindblad run to carry out
INTEGER, PARAMETER :: ntraj = 10000
! Duration for the Redfield propagation (au)
REAL*8, PARAMETER :: duration = 10000.0d0
! Time step for rk4 routine
REAL*8, PARAMETER :: dt = 5.000d0 !duration / (ntimes - 1)
! How many iterations of redfield to do between printings of position/diabatic state
INTEGER, PARAMETER :: print_num = 1

! Constants and conversions
REAL*8, PARAMETER :: pi = 3.1415926d0
REAL*8, PARAMETER :: hbar = 1.0d0 

! Bath parameters
! Number of baths which will be coupled to some coordinate
INTEGER, PARAMETER :: nbath = 1
! Cutoff parameter in spectral density for each bath
REAL*8, PARAMETER, DIMENSION(nbath) :: omegac = (/0.006269431d0/)
! Strength parameter in spectral density
REAL*8, PARAMETER, DIMENSION(nbath) :: eta = (/0.0018228d0/)
! Debye bath is 1, Ohmic is 0
! ohmic is the default if no valid type is given
INTEGER, PARAMETER, DIMENSION(nbath) :: bath_type = (/0/)
REAL*8, PARAMETER :: beta = 1052.584412992859d0

REAL*8, PARAMETER :: factor = 0.99
! SM parameters; Initialization WILL ALWAYS USE
! A SINGLE EIGENSTATE (thermal initialization is deprecated)
INTEGER, PARAMETER :: dim_qc = 70 !! Harmonic coordinate dimension; must be even
INTEGER, PARAMETER :: dim_R = 101 ! DVR coordiante dimension; must be odd
REAL*8, PARAMETER :: dx_R = 0.10  ! DVR delta x placement information
REAL*8, PARAMETER :: dx_Qc = 0.10d0  ! DVR delta x placement information
REAL*8, PARAMETER :: omegasc = 0.006269431d0*factor  ! Harmonic coordinate spring constant
REAL*8, PARAMETER :: mass_R = 1836.0d0  ! DVR coordinate mass
! Parameters for the potential energy surface and H
! for the adiabatic QED Hamiltonian:
! H = P^2/2M + E(R) + H_vib + pc^2/2 + wc^2/2 *
! (qc + sqrt(2/(hbar wc^3))xi * mu(R))^2 + coeff_lp*R + H_ls
! mu(R) = coeff_v*tanh(y*R) + coeff_z*R
! E(R) = a*R^2 + b*R^4 + c
REAL*8, PARAMETER :: coeff_a = -0.573838d0/27.211396d0
REAL*8, PARAMETER :: coeff_b = 0.0900909d0/27.211396d0
REAL*8, PARAMETER :: coeff_c = 0.902345d0/27.211396d0
REAL*8, PARAMETER :: coeff_v = -1.90249d0
REAL*8, PARAMETER :: coeff_y = 1.26426d0
REAL*8, PARAMETER :: coeff_z = 0.37044d0
REAL*8, PARAMETER :: coeff_xi = 0.00234562*factor !0.1d0/27.211396d0
REAL*8, PARAMETER :: coeff_lp = 1.0*0.0001d0/27.211396d0  

! Size parameters
! Dimension of density matrix; must be an even number
INTEGER, PARAMETER :: n = 200   ! Truncated basis size, must be less than or equal to m
INTEGER, PARAMETER :: m = dim_qc*dim_R ! Basis size for initial diagonalization (dim_nc*dim_nt*2 for a 2-surface CI)

! Populate a single eigenstate
INTEGER, PARAMETER :: init_stat_therm = 1  ! This must remain set to 1 for this system
! Eigenstate to initialize all density in if init_stat_therm = 1
INTEGER, PARAMETER :: init_es = 17

! Plotting bounds; plots eigenstate wavefunctions as specified between each pair of bounds
INTEGER, PARAMETER :: bound1 = 1
INTEGER, PARAMETER :: bound2 = 10
INTEGER, PARAMETER :: bound3 = 11
INTEGER, PARAMETER :: bound4 = 20

END MODULE parameters
