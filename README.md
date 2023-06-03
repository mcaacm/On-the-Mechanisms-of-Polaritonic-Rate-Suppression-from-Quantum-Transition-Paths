The following contains code and data accompany 'On the Mechanism of Polaritonic
Rate Suppression from Quantum Transition Paths' by
Michelle Anderson, Esmae Woods, Thomas Fay, David Wales and David Limmer, 2023.


PAPER DATA

The data and gnuplot script necessary to regenerate the figures in the paper are found
in the folder title 'figures' with the exception of the disconnectivity plots which 
are generated by a seperate program set. Some postscript disconnectivity plots are found
in 'figures/disconnectivity.'  Example input files and their associated
Markov state models are located in 'inputs_and_msms.' Mean first passage time data are located
in 'figures/FPT_data.'

PLEASE NOTE THAT THE VALUE OF THE PARAMETER OMEGA_S, the default frequency of the 
proton which other frequencies are measured in comparison to, WAS CHANGED after
data was generated. All omega_c/omega_s values are multiplied by a value, 'ratio=0.9251,'
to accomodate for this. This multiplication is done by the plotting script.
Files labels have undergone multiplication by 'ratio' but file contents have not.
This is generally mentioned as a reminder in file headers. Keep this in mind when
looking at the example input files which are labaled by their adjusted values
but contain input parameters using the original omega_s value.

To generate the plots in 'figures' run:
'gnuplot plot.gs'

Output will be generated in files 'figure_1.ps', 'figure_2.ps', 'figure_3.ps', 
'heom_compare.ps', 'mfpt.ps'

The full 'figure_1' includes a disconnectivity plot which is appended in. The first
figure in the supporting information is a concatenation of three disconnectivity plots.

--------------------------------------------------------------------------------------------

PROGRAMS

The programs necessary to carry out secular Redfield calculations, 
Markov state model analysis, and plotting of eigenstate wavefunctions
are main, create_msm, markov_analysis, and wfn_plot.

The file params.f90 contains specifications to control the
programs, with a few controls being present in other program files, most
of these referring to very specific functions such as whether to impose
detailed balance when calculating committors in markov_analysis.

The programs can be compiled from the given Makefile with

'make'

Most programs require Lapack and BLAS. Quadpack is also required
but the necessary file is included and compiled by make. Quadpack
is opensource software or licensed as specified in 'quadpack.f90'.

--------------------------------------------------------------------------------------------

In order to perform a full run including generating a Markov state
model and producing analysis, run: 

make  
./create_msm  # Generate a Markov state model called markov_mat.txt
./markov_analysis # Run analysis on the model found in markov_mat.txt 

--------------------------------------------------------------------------------------------

Programs in General:

Programs were constructed and tested on Mac OS Ventura 13.3.1 and compiled with 
gfortran 8.2.0. The makefile settings assume that gfortran, Lapack and BLAS are available.
To change the compiler or library linking settings, edit 'Makefile'.

Precision and details of the LAPACK and BLAS installations may result in small 
changes, including inversion of the relative sign of wavefunctions.

Note, all programs that call a setup function will write a file 'evic.txt'
which includes information about the average characteristics of energy eigenvectors in the format
'(R coordinate) (q_c coordinate) (energy) (deprecated measurement for historical reasons)'.

All programs also produce a file, 'plot_ev.txt' which merely contains 
'(eigenvalue number) (arbitrary plot coordinate) (eigenvector energy)' which can be used to plot
the energies of the eigenvectors.

Some programs will produce a file, 'overlap_elements_1.txt' which contains the values of
|R_ij|^2 for the first thirty eigenstates in a matrix format.

It is recommended to use atomic units for all calculations. Although the value of reduced
Planck's constant can be changed in the parameters file, this option was not thoroughly 
tested.

--------------------------------------------------------------------------------------------

main

Depending on parameters specified in 'params.f90', main  will run either
a secular or non-secular Redfield calculation for a fixed period of time starting 
from the specified eigenstate (init_es).
Parameters to specify the number of baths as well as their cut off frequencies,
coupling strengths and temperatures are available. The basis size 'm' must
always be (basis size in R)x(basis size in q_c) but 'n', the truncated
basis size, can be altered to be any number less than or equal to m. The default is 200.

The file 'state_prop.txt' will include (iteration) (time) and
the populations of the first twenty-five eigenstates if 'n' is larger than this.

--------------------------------------------------------------------------------------------

create_msm

Runs a short secular calculation starting from each eigenstate in turn for the 
specified system and assembles the populations for each run into a Markov state model
printed to 'markov_mat.txt' with the following format:
'(n=matrix dimension) (t=duration of propogation)
Equilibrium populations of all energy eigenstates in the system
Populations of eigenstates at t following initialization in energy eigenstate 1
Populations of eigenstates at t following initialization in energy eigenstate 2
...
Populations of eigenstates at t following initialization in energy eigenstate n
Eigenstate energies'

--------------------------------------------------------------------------------------------

markov_analysis

Reads in 'markov_mat.txt' and 'evic.txt.' Committors from reactant (energy eigenstate 1) to product
(energy eigenstate labeled by B, default 2) will be written in 'thermal_comm_Bto1.txt in the 
format '(eigenstate) (committor to 1) (committor to B)'. The format for 'thermal_comm_1toB.txt'
is the same with the destination and origin reversed. The reactive pathway information is
in 'thermal_paths_1toB.txt' and 'thermal_paths_1toB_plot.txt' with the latter being designed
for plotting and the former being human readable and containing extra information.
Eigestate B is specified by a parameter 'B_estate' in 'markov_analysis_utils.f90'.
If 'evic.txt' is not read in, specifications about 'q_c' and 'R' in the analysis will
be garbage.

Note that it is possible to change the committor eigenstate on the command line by providing
an argument, but it's better to do it in 'markov_analysis_utils.f90.'

A file that can be used to plot some of the pathway information is 
'thermal_comm_1toB_plot.txt' but it was not used for this particualr investigation.

'thermal_paths_1toB.txt' prints the total reactive flux then the reactive paths in order of
largest flux. The format for this is
'(pathway flux) (fraction total flux on this pathway) (energy eigenstates on pathway in order)'

In the next section, the same information is then printed again, but this time the eigenstates 
are replaced by the energy of the eigenstates for each entry on the path.

This is followed by some information about the top few paths in terms of cummulative flux
including the transition state energy, ground state energy (lowest eigenstate energy),
portion of the rate accounted for by that path, the change in q_c (Qc) and R (called Qt for
historical reasons) as well as two deprecated parameters (the chance in P1 and P2) then
the change in committor value during the committor jump and the change in energy during the
committor jump. 

Path transmission coefficients (\kappa in the paper) are then presented along with
the change in energy between the transition eigenstates and the change in committor
probability. The average (not weighted by pathway flux) and the average (weighted
by pathway flux) for the transmission coefficients concludes this section.

After pathway information, the total flux is again printed followed by specific information
about the jump at which each given pathway exceeds a committor probability of 50%, i.e. each
committor jump is printed along with the total flux that passes through that jump and the
percentage of total flux that passes through that jump.

The average over all paths of energy, coordinates q_c (QC) and R (QT) as well as deprecated
P1 and P2 (once used in non-adiabatic systems) and the forward committor value (COMM_F) for 
the committor eigenstates is presented along with the difference between pre and
post-committor eigenstates.

The next section of the file is a running sum of the total and total percent reactive flux accounted
for as a function of reactive pathways. The pathways are in order of decreasing flux. This
section concludes by printing the path entropy for the ensemble.

The file 'thermal_paths_Bto1.txt' contains the same information but for the reverse reaction.
The results found here should be very similar but not necessarily identical for numerical
reasons and because of the slight asymmetry of the system. The small difference in rates is due
to the factor \Pi_a in transition path theory rate calculations.

The file 'flux.txt' shows the reactive flux and rate for thermal passage from reactants to products
and from products to reactants. As this system is nearly symmetric, these should be quite
close in value but not perfectly identical.

--------------------------------------------------------------------------------------------

wfn_plot

This plots all wavefunctions specified in 'params.f90' between bound1 and bound2 and 
between bound3 and bound4. The wavefunctions are written to files 'data/wfn_plots/ef_2d_i.txt'.
The containing folders have to be made before the program will run, otherwise it will terminate.

When the program runs, it will print out the integrated wavefunction density of the eigenstates
it plots. These should be approximately 1.

'data/wfn_plots/ef_2d_i.txt' has the format for eigenstate wavefunction 'i':
'(q_c coordinate) (R coordinate) (wfn density) (real wfn component) (imaginary wfn component)

