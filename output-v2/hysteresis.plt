set terminal postscript eps enhanced colour font ',20'
set output 'hysteresis.eps'

set size square
set xlabel '{/Symbol W}' font ',25'
set label '{/Symbol a}' at graph -0.18, graph 0.52 font ',25'
set xrange [0:0.5]
set yrange [0.99:2.25]
set key title 'x values'
set key top left

plot 'x0.dat' every 200 u 4:5 w l lw 2 dt 1 lc rgb "#0072bd" smooth bezier t '0',\
     'x0-spindown.dat' every 200 u ($4>0.25 ? $4 : 1/0):5 w l lw 2 dt 2 lc rgb "#0072bd" smooth bezier t '',\
     'x0.2.dat' every 200 u 4:5 w l lw 2 dt 1 lc rgb "#a2142f" smooth bezier t '0.2',\
     'x0.2-spindown.dat' every 200 u 4:5 w l lw 2 dt 2 lc rgb "#a2142f" smooth bezier t '',\
     'x0.4.dat' every 200 u 4:5 w l lw 2 dt 1 lc rgb "#77ac30" smooth bezier t '0.4',\
     'x0.4-spindown.dat' every 200 u 4:5 w l lw 2 dt 2 lc rgb "#77ac30" smooth bezier t '',\
     'x0.6.dat' every 400 u 4:5 w l lw 2 dt 1 lc rgb "#7e2f8e" smooth bezier t '0.6',\
     'x0.6-spindown.dat' every 400 u 4:5 w l lw 2 dt 2 lc rgb "#7e2f8e" smooth bezier t ''
