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

! Run an analysis to try to find useful committor
! and trajectory information from a MSM which is 
! read in from a file whose first entry is the dimension, then time
! then the rest is the matrix (doubles)
PROGRAM analyze_markov
 
USE analyze_markov_utils

IMPLICIT NONE 

REAL*8, DIMENSION(:,:), ALLOCATABLE :: mmat  ! The markov rate matrix
REAL*8 :: dt
! File number, number of entries in intermediate states, reactant states, product states, dimension, loop indices
INTEGER :: fnum, ni, na, nb, d, i, j
! Committors (forward, backward), equilibrium probs, non-equilibrium probabilities (d+1), collapse probs, eigenvalues
! neq_pis and cprobs are depreated; they were bad ideas
REAL*8, DIMENSION(:), ALLOCATABLE :: comm_f_to_1, comm_f_to_B, comm_b_to_1, comm_b_to_B, pis, neq_pis, cprobs, evals  
! Index definitions for reactant, product, and intermediate states
INTEGER, DIMENSION(:), ALLOCATABLE :: A_def, B_def, I_def
! Eigenvector energy, coupling coordiante, tuning coordinate, probability to be in diabatic state 1 and 2
! These may or may not be present and the routine can continue without them
REAL*8, DIMENSION(:),  ALLOCATABLE :: evec_en, evec_pc, evec_pt, evec_p1, evec_p2
REAL*8 :: factor
INTEGER :: B  ! The bottom of well 2; product. Eigenstate 1 is reactant

! Read in the information from the markov state file
fnum = 15
OPEN(15,FILE="markov_mat.txt")
  CALL read_msm(mmat,pis,neq_pis,cprobs,evals,d,dt,fnum)
CLOSE(15)

WRITE(*,*) "MSM read"

ALLOCATE(A_def(d))
ALLOCATE(B_def(d))
ALLOCATE(I_def(d))
ALLOCATE(comm_f_to_1(d))
ALLOCATE(comm_f_to_B(d))
ALLOCATE(comm_b_to_1(d))
ALLOCATE(comm_b_to_B(d))

! Allocate holders for eigenvector information
ALLOCATE(evec_en(d))
ALLOCATE(evec_pc(d))
ALLOCATE(evec_pt(d))
ALLOCATE(evec_p1(d))
ALLOCATE(evec_p2(d))

OPEN(16,FILE="evic.txt")
CALL read_evi(d,16,evec_en,evec_pc,evec_pt,evec_p1,evec_p2)
evals = evec_en  ! This is a more reliable way to read it in
ClOSE(16)

B = get_arg()

A_def = 0
B_def = 0
I_def = 0
comm_f_to_1 = 0.0d0
comm_f_to_B = 0.0d0
comm_b_to_1 = 0.0d0
comm_b_to_B = 0.0d0
factor = 1.0d0

! Thermal barrier crossing: bottom of well 1 to bottom of well 2 1---->B
na = 1
nb = 1
ni = d - 2
j = 1
DO i = 1, d
  IF (i .NE. 1 .AND. i .NE. B) THEN
    I_def(j) = i
    j = j + 1
  END IF
END DO
A_def(1) = 1
B_def(1) = B

OPEN(21,FILE="thermal_comm_1toB.txt")
OPEN(22,FILE="flux.txt")
OPEN(23,FILE="thermal_paths_1toB.txt")
OPEN(24,FILE="thermal_paths_1toB_plot.txt")
WRITE(22,*) "FLUX thermal 1 to B"
CALL perform_analysis(d,ni,na,nb,A_def,B_def,I_def,dt,evals,factor,1,B,21,22,23,24,mmat,pis,&
                      comm_f_to_B,comm_b_to_B,evec_pc,evec_pt,evec_p1,evec_p2)
CLOSE(21)
CLOSE(23)
CLOSE(24)

! Thermal barrier crossing: bottom of well 2 to bottom of well 1 B---->1
A_def(1) = B
B_def(1) = 1

OPEN(21,FILE="thermal_comm_Bto1.txt")
OPEN(23,FILE="thermal_paths_Bto1.txt")
OPEN(24,FILE="thermal_paths_Bto1_plot.txt")
WRITE(22,*) "FLUX thermal B to 1"
CALL perform_analysis(d,ni,na,nb,A_def,B_def,I_def,dt,evals,factor,B,1,21,22,23,24,mmat,pis,&
                      comm_f_to_1,comm_b_to_1,evec_pc,evec_pt,evec_p1,evec_p2)
CLOSE(21)
CLOSE(22)
CLOSE(23)
CLOSE(24)


END PROGRAM


