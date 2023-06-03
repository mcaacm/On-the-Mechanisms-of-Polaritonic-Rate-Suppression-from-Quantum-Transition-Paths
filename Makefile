# "On the mechanism of polaritonic rate suppression from quantum transition paths"
# Michelle C Anderson, Esmae J. Woods, Thomas P. Fay, David J. Wales and David T. Limmer
# University of California Berkeley
# Copyright 2023

# This file is part of TPT_SM 
# TPT_SM is free software: you can redistribute it and/or modify it under the terms of the GNU General 
# Public License as published by the Free Software Foundation, either version 3 of the License, 
# or (at your option) any later version.
# TPT_SM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with TPT_SM. 
# If not, see <https://www.gnu.org/licenses/>.


objs = quadpack.o prequel.o rkns.o setup.o parameters.o markov_analysis_utils.o 

flags = -Wall -Wno-conversion -g -fcheck=all # LAPACK/BLAS/QUADPACK calls may produce INF/NAN and appear to handle such situations appropriately. Do not compile to catch IEE FPEs

link_flags = -llapack -lblas

comp = gfortran

all : $(objs)
	$(comp) $(flags) create_msm.f90 -o create_msm setup.o rkns.o prequel.o parameters.o quadpack.o $(link_flags)
	$(comp) $(flags) wfn_plot.f90 -o wfn_plot setup.o rkns.o prequel.o parameters.o quadpack.o $(link_flags)
	$(comp) $(flags) markov_analysis.f90 markov_analysis_utils.o -o markov_analysis $(link_flags)
	$(comp) $(flags) main.f90 -o main  setup.o rkns.o prequel.o parameters.o quadpack.o $(link_flags)
quadpack.o: quadpack.f90
	$(comp) $(flags) -c quadpack.f90 -o quadpack.o 

parameters.o: params.f90
	$(comp) $(flags) -c params.f90 -o parameters.o

prequel.o: prequel.f90 parameters.o
	$(comp) $(flags) -c prequel.f90 -o prequel.o 

setup.o: setup_H.f90 parameters.o
	$(comp) $(flags) -c setup_H.f90 -o setup.o

rkns.o: prequel.o quadpack.o rkns.f90 parameters.o
	$(comp) $(flags) -c rkns.f90 -o rkns.o

markov_analysis_utils.o: markov_analysis_utils.f90 
	$(comp) $(flags) -c markov_analysis_utils.f90 -o markov_analysis_utils.o $(link_flags)

