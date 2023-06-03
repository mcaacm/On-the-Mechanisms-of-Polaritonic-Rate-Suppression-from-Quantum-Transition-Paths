set terminal postscript enhanced portrait colour font "Times-Roman,18"
ratio = 0.9251

# VERY IMPORTANT TO AVOID CONFUSION: The value of \omega_s, the well frequency, was changed after the
# calculations were run from 0.00627 to 0.00677 due to an initial oversight when modifying the model.
# This requires a rescaling ratio of 0.9251 on all the x-axis that show 
# \omega_c/\omega_s. This is why "ratio" exists. It is used to adjust the x-axis.
# It is applied universally.





# Plot the fitted potential energy surfaces alongside the data to which they are fitted
# Note that the disconnectivity plot included in this figure must be appended in manually.
# The tree in question is found in disconnectivity/tree_default_0.8.ps
set output "figure_1.ps"
set multiplot

# Plot the 3d polariton model surface

unset cbtics
unset colorbox
set size 0.9,0.47
set origin 0.05,0.5
set label "a)" at screen 0.05,0.85
set label "b)" at screen 0.05,0.53
set label "c)" at screen 0.05,0.29
set label "d)" at screen 0.57,0.53
set surface
set xlabel "{/Times-Italic R}/au"
set ylabel "{/Times-Italic q_c}/au"

factor = 1.0  # NOTE THAT RATIO HAS NOT BEEN APPLIED TO THIS
Q0 = 1
coeff_a = -0.573838/27.211396
coeff_b = 0.0900909/27.211396
coeff_c = 0.902345/27.211396
coeff_v = -1.90249
coeff_y = 1.26426
coeff_z = 0.37044
coeff_xi = factor*0.00234562 
coeff_lp = 0.0001/27.211396
omegac = 0.006269431
eta = 0.001822766
omegasc = factor*0.006269431
hbar = 1.0
pi = 3.14159
mass_R = 1836

mu(x) = coeff_v*tanh(coeff_y*x) + coeff_z*x
H_ls(x) = (x**2.0)*omegac*eta/(mass_R*pi)
f(x,y) = coeff_a*x**2.0 + coeff_b*x**4.0 + coeff_c + 0.5*(omegasc**2.0)*(y + sqrt(2.0/(hbar*omegasc**3.0))*mu(x)*coeff_xi)**2.0 + H_ls(x) + x*coeff_lp 
u(x) = coeff_a*x**2 + coeff_b*x**4 + coeff_c + coeff_lp*x + H_ls(x)

set xyplane at 0
set yrange [-28:40]
set xrange [-2.8:2.3]
set ytics -20,20,40
set zrange [0:1.2]
set cbrange [0:1.2]
set xlabel "{/Times-Italic R}/au"
set ylabel "{/Times-Italic q_c}/au"
set zlabel "{/Times-Italic E}/eV"

set isosamples 1500,1500
set samples 1000
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot f(x,y)*27.211 with pm3d notitle

set xrange [-5:5]
set ylabel "{/Times-Italic U(R)}/eV"
set origin 0.04,0.30
set size 0.55,0.26
set ytics 0,1,3
unset xlabel
set yrange [-0.05:3.0]
set xrange [-3:3] 
set xtics ("" -3, "" -2, "" -1, "" 0, "" 1, "" 2, "" 3)
plot u(x)*27.211  lw 7 lc "black" notitle  # Convert to electron volts
set ylabel "{/Symbol-Oblique m}{/Times-Italic (R)}/au"
set size 0.55,0.33
set origin 0.04,0.00
set xtics -3,1,3
set ytics -1.0,1.0,1.0
set xlabel "{/Times-Italic R}/au"
set yrange [-1.5:1.5]
plot mu(x) lw 7 lc "black" notitle

unset multiplot





# Plot the rates, transmission coefficients, and path entropy
reset
set terminal postscript enhanced portrait colour font "Times-Roman,20"
set output "figure_2.ps"
set multiplot
# Plot the rates relative to reference values
set label "a)" at screen 0.04,0.92
set label "b)" at screen 0.04,0.62
set label "c)" at screen 0.04,0.32
set size 1.1,0.36
set origin  0.0,0.60
set key bottom right
set xrange [ratio*0.5:ratio*2.0]
set xlabel "{/Symbol-Oblique w}_{/Times-Italic c}/{/Symbol-Oblique w}_{/Times-Italic s}"
set ylabel "{/Times-Italic k}/{/Times-Italic k_{0}}({/Symbol-Oblique h}_{/Times-Italic c})"
set yrange [0:1.1]
set ytics 0.2,0.2,1.1
unset xlabel
unset xtics
set xtics ("" 0.5, "" 0.7, "" 0.9, "" 1.1, "" 1.3, "" 1.5, "" 1.7, "" 1.9)
a = 7.643E-017 
b = 9.782E-017
c = 1.042E-016
plot "ratedefault.txt" u ($1*ratio):($5/a) with lines lw 5 lc "black" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=1.0", "ratehalf.txt" u ($1*ratio):($5/b) with lines lw 5 lc "blue" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.5", "ratequarter.txt" u ($1*ratio):($5/c) with lines lw 5 lc "red" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.25", 


# Plot the kappa values for the dominant reactive pathways relative to references and the entropy as well
unset xtics
set xtics ("" 0.5, "" 0.7, "" 0.9, "" 1.1, "" 1.3, "" 1.5, "" 1.7, "" 1.9)
unset xlabel
unset key
set ytics 0.2,0.2,1.0
set ylabel "{/Symbol-Oblique k}/{/Symbol-Oblique k}_{/Times-Italic 0}({/Symbol-Oblique h}_{/Times-Italic c})"
a = 1.6667E-008
b = 2.1801E-008
c = 2.3307E-008  
set origin 0.0,0.30
set size 1.1,0.36
plot "kappadefault.txt" u ($1*ratio):($2/a) with lines lw 5 lc "black" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=1.0", "kappahalf.txt" u ($1*ratio):($2/b) with lines lw 5 lc "blue" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.5", "kappaquarter.txt" u ($1*ratio):($2/c) with lines lw 5 lc "red" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.25"
set size 1.055,0.42
set origin 0.045,-0.06
unset key
set yrange [1:3]
set ytics 1.0,1.0,3
set ylabel "{/Times-Italic S}"
set xlabel "{/Times-Italic {/Symbol-Oblique w}_c/{/Symbol-Oblique w}_s}"
set xtics 0.5,0.2,2.0
plot "entropydefault.txt" u ($1*ratio):2 with lines lw 5 lc "black" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=1.0", "entropyhalf.txt" u ($1*ratio):2 with lines lw 5 lc "blue" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.5", "entropyquarter.txt" u ($1*ratio):2 with lines lw 5 lc "red" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.25"

unset multiplot






# Plot the overlap matrix element for the committor eigenstates as well
# as the committor eigenstates themselves as a heat map
set terminal postscript enhanced portrait colour font "Times-Roman,22"
reset
set label "a)" at screen -0.02,0.99
set label "b)" at screen -0.02,0.6
set output "figure_3.ps"
set multiplot
unset key
set ylabel "|{/Times-Italic R_{i,j}}|^2x1000"
set yrange [1.50:7.00]
set ytics 1.5,1.0,7.5
set xrange [ratio*0.5:ratio*2.0]
set size 1.20,0.48
set origin -0.15,0.57
set xtics 0.6,0.4,2.0
set xlabel "{/Times-Italic {/Symbol-Oblique w}_c/{/Symbol-Oblique w}_s}"
plot "coupling_elements.txt" u ($1*ratio):($2*1000) with lines lw 5 lc "black" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.25"

reset
set ytics -40,40,40
set cbtics -0.25,0.25,0.25
set colorbox
set xtics ("" -3, "" -2, "" -1, "" 0, "" 1, "" 2, "" 3)
unset ylabel

factor = 0.8  # NOTE THAT RATIO HAS NOT BEEN APPLIED TO THIS
coeff_xi = factor*0.00234562 
omegasc = factor*0.006269431

set yrange [-50:50]
set xrange [-3:3]
set cbrange [-0.3:0.3]
unset cbtics
unset colorbox
set ylabel "{/Times-Italic q_c}/au"
set origin 0.06,0.31
set size 0.58,0.33

set surface
unset contour
set view map
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot "default_0.74_19.txt" u ($2/Q0):($1/Q0):4 with pm3d at b notitle 
unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.01,50.0
splot f(x,y) with lines lc "black" notitle
set origin 0.48,0.31
set cbtics
set colorbox
unset ylabel
unset ytics
set ytics ("" -40, "" 0, "" 40)
set surface
unset contour
set view map
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot "default_0.74_20.txt" u ($2/Q0):($1/Q0):4 with pm3d at b notitle 
unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.01,50.0
splot f(x,y) with lines lc "black" notitle

factor = 0.99  # NOTE THAT RATIO HAS NOT BEEN APPLIED TO THIS
coeff_xi = factor*0.00234562 
omegasc = factor*0.006269431
set origin 0.06,0.13
unset cbtics
unset colorbox
set ylabel "{/Times-Italic q_c}/au"
set ytics -40,40,40

set surface
unset contour
set view map
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot "default_0.916_13.txt" u ($2/Q0):($1/Q0):4 with pm3d at b notitle 
unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.01,50.0
splot f(x,y) with lines lc "black" notitle
set origin 0.48,0.13
set cbtics -0.25,0.25,0.25
set colorbox
unset ylabel
unset ytics
set ytics ("" -40, "" 0, "" 40)
set surface
unset contour
set view map
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot "default_0.916_14.txt" u ($2/Q0):($1/Q0):4 with pm3d at b notitle 
unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.01,50.0
splot f(x,y) with lines lc "black" notitle

factor = 1.1  # NOTE THAT RATIO HAS NOT BEEN APPLIED TO THIS
coeff_xi = factor*0.00234562 
omegasc = factor*0.006269431
set xtics -3,1,3
set origin 0.06,-0.05
set xtics
unset cbtics 
unset colorbox
set ylabel "{/Times-Italic q_c}/au"
set xlabel "{/Times-Italic R}/au"
set ytics -40,40,40

set surface
unset contour
set view map
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot "default_1.02_13.txt" u ($2/Q0):($1/Q0):4 with pm3d at b notitle 
unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.01,50.0
splot f(x,y) with lines lc "black" notitle
set origin 0.48,-0.05
set cbtics -0.25,0.25,0.25
set colorbox
unset ylabel
unset ytics
set ytics ("" -40, "" 0, "" 40)
set surface
unset contour
set view map
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot "default_1.02_14.txt" u ($2/Q0):($1/Q0):4 with pm3d at b notitle 
unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.01,50.0
splot f(x,y) with lines lc "black" notitle

unset multiplot




# Plot mean first passage times
reset
set terminal postscript enhanced landscape colour font "Times-Roman,20"
set xrange [28:42]
set yrange [0:0.40]
set xtics 30,10,40
set ytics 0,0.1,0.4
set label "a)" at screen 0.01,0.55
set label "b)" at screen 0.375,0.55
set label "c)" at screen 0.70,0.55
set output "mfpt.ps"
set multiplot
set xlabel "ln{/Times-Italic (t/}au{/Times-Italic )}"
set ylabel "{/Palatino-Italic P}"
unset key
set origin -0.1,0.1
set size 0.50,0.5
plot "FPT_data/default/lntpdf_0.92_1to2" with lines lw 3 lc "black" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=1.0", "FPT_data/default/lntpdf_0.74_1to2" with lines lw 3 lc "black" dt "-" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.5"
unset ylabel
set size 0.35,0.5
set origin 0.375,0.1
set ytics ("" 0.1, "" 0.2, "" 0.3,  "" 0.4)
plot "FPT_data/half/lntpdf_0.89_1to2" with lines lw 3 lc "blue" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.5", "FPT_data/half/lntpdf_0.83_1to2" with lines lw 3 lc "blue" dt "-" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.5"

set origin 0.70,0.1
plot "FPT_data/quarter/lntpdf_0.86_1to2" with lines lw 3 lc "red" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.25", "FPT_data/quarter/lntpdf_0.83_1to2" with lines lw 3 lc "red" dt "-" title "{/Symbol-Oblique h}_{/Times-Italic c}/{/Symbol-Oblique h}_{/Times-Italic s}=0.25" 


unset multiplot



# Plot the heom/non-secular and secular Redfield short-time rate comparisons
set terminal postscript enhanced colour font "Times-Roman,22"
reset
set output "heom_compare.ps"
set multiplot
set size 0.40,0.65
k0 = 1e-9
set origin 0.0,0.1
set xrange [0.85*ratio:1.05*ratio]
set xtics 0.85,0.1,1.05
set key at screen 1.05,0.5
set xlabel "{/Symbol-Oblique w}/{/Symbol-Oblique w}_{/Times-Italic s}"
set ylabel "{/Times-Italic k_{10000}/k_s}"
plot "heom_rates.txt" u ($1*ratio):($2/(k0*10000)) ps 2 pt 5 lc "black" title "HEOM", "redfield_rates.txt" u ($1*ratio):($2/(k0*10000)) ps 2 pt 7 lc "red" title "Redfield", "secular_rates.txt" u ($1*ratio):($2/(k0*10000)) ps 2 pt 9 lc "blue" title "Secular Redfield"
set ylabel "Linear Short-Time Committor Hop Rate Estimate (au^{-1}) From 2500 au Runs"
set ylabel "{/Times-Italic k_{2500}/k_s}"
set origin 0.35,0.1
unset key
plot "heom_rates.txt" u ($1*ratio):($3/(k0*2500)) ps 2 pt 5 lc "black" title "HEOM", "redfield_rates.txt" u ($1*ratio):($3/(k0*2500)) ps 2 pt 7 lc "red" title "Redfield", "secular_rates.txt" u ($1*ratio):($3/(k0*2500)) ps 2 pt 9 lc "blue" title "Secular Redfield"

unset multiplot



