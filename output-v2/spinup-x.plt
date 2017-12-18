set terminal postscript eps enhanced colour font ',20'
set output 'spinup-x.eps'
set xlabel '{/Symbol W}' font ',25'
set label '{/Symbol a}' at graph -0.18, graph 0.52 font ',25'
set yrange [0.98:2.3]
set xtic 0,0.1
set size ratio 1

set key title 'x values'
set key at 0.14,1.6

#set arrow from 0.35,1.4 to 0.48,1.41
#set label "A" at 0.33,1.4
#set arrow from 0.35,1.6 to 0.465,1.65
#set label "B" at 0.33,1.6
#set arrow from 0.3,1.9 to 0.345,2.07
#set label "C" at 0.285,1.87
#set arrow from 0.25,2 to 0.29,2.17
#set label "D" at 0.235,1.97

plot 'x0.dat' every 200 u 4:5 w l lw 3 dt 1 lc rgb "#4dbeee" smooth bezier t '0',\
     'x0.2.dat' every 200 u 4:5 w l lw 3 dt 2 lc rgb "#edb120" smooth bezier t '0.2',\
     'x0.4.dat' every 200 u 4:5 w l lw 3 dt 3 lc rgb "#7e2f8e" smooth bezier t '0.4',\
     'x0.6.dat' every 400 u 4:5 w l lw 3 dt 4 lc rgb "#77ac30" smooth bezier t '0.6',\
     'x0.8.dat' every 400 u 4:5 w l lw 3 dt 5 lc rgb "#a2142f" smooth bezier t '0.8'
