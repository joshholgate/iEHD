set terminal postscript eps enhanced colour font ",20"
set output "profiles.eps"

set multiplot layout 1,5 \
              margins 0.1,1,0.1,1 \
              spacing 0
set size ratio 2.333333333333333
set xrange [-1.5:1.5]
set yrange [1:7]
unset key
set tics front
set ticscale 0.5
set xtics -1,1
set ytics 1,1

set contour base
set cntrparam levels disc 0
unset surface
set table 'cont.dat'
splot 'E1.2/output/t0_gridvalues.dat' u 1:2:6, 'E1.2/output/t0_gridvalues.dat' u ($1)*(-1):2:6
unset table

set contour base
set cntrparam levels disc 0
unset surface
set table 'cont1.dat'
splot 'E1.2/output/t5_gridvalues.dat' u 1:2:6, 'E1.2/output/t5_gridvalues.dat' u ($1)*(-1):2:6
unset table

set contour base
set cntrparam levels disc 0
unset surface
set table 'cont2.dat'
splot 'E1.2/output/t15_gridvalues.dat' u 1:2:6, 'E1.2/output/t15_gridvalues.dat' u ($1)*(-1):2:6
unset table

set contour base
set cntrparam levels disc 0
unset surface
set table 'cont3.dat'
splot 'E1.2/output/t25_gridvalues.dat' u 1:2:6, 'E1.2/output/t25_gridvalues.dat' u ($1)*(-1):2:6
unset table

set contour base
set cntrparam levels disc 0
unset surface
set table 'cont4.dat'
splot 'E1.2/output/t35_gridvalues.dat' u 1:2:6, 'E1.2/output/t35_gridvalues.dat' u ($1)*(-1):2:6
unset table

set palette defined (-1 "white", 1 "gray") 
set palette maxcolors 2
set cbrange [-0.1:0.1]
unset colorbox


set title 't/{/Symbol t} = 0' offset 0,-0.5
plot 'E1.2/output/t0_gridvalues.dat' u 1:2:6 w image,\
     'E1.2/output/t0_gridvalues.dat' u ($1)*(-1):2:6 w image,\
     'cont.dat' w l lt -1 lw 2
set ytics format ""
set title 't/{/Symbol t} = 2.5' offset 0,-0.5
plot 'E1.2/output/t5_gridvalues.dat' u 1:2:6 w image,\
     'E1.2/output/t5_gridvalues.dat' u ($1)*(-1):2:6 w image,\
     'cont1.dat' w l lt -1 lw 2
set title 't/{/Symbol t} = 7.5' offset 0,-0.5
plot 'E1.2/output/t15_gridvalues.dat' u 1:2:6 w image,\
     'E1.2/output/t15_gridvalues.dat' u ($1)*(-1):2:6 w image,\
     'cont2.dat' w l lt -1 lw 2
set title 't/{/Symbol t} = 12.5' offset 0,-0.5
plot 'E1.2/output/t25_gridvalues.dat' u 1:2:6 w image,\
     'E1.2/output/t25_gridvalues.dat' u ($1)*(-1):2:6 w image,\
     'cont3.dat' w l lt -1 lw 2
set title 't/{/Symbol t} = 17.5' offset 0,-0.5
plot 'E1.2/output/t35_gridvalues.dat' u 1:2:6 w image,\
     'E1.2/output/t35_gridvalues.dat' u ($1)*(-1):2:6 w image,\
     'cont4.dat' w l lt -1 lw 2
