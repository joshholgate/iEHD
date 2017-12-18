set terminal postscript eps enhanced colour font ',20'
set output "LvsX.eps"

set size square
set xlabel 'x' font ',25' offset 0,0.5
set label 'L' at graph -0.18, graph 0.52 font ',25'

set label "Oblate" at 0.2,0.2 front font ',25'
set label "Prolate" at 0.2,0.7 front font ',25'
set label "Unstable" at 0.6,0.9 front font ',25'

set label 'Spheroid model' at 0.5,0.05 front font ',20' 
set label 'FDLS' at 0.8,0.5 front font ',20' 
set label 'simulations' at 0.75,0.45 front font ',20'
set label 'L^{max} = 1.15 - 0.59x - 0.56x^2' at 0.1,1.15 front font ',20'
set label '  _{II}' at 0.1,1.14 front font ',20'

set arrow from 0.7,0.08 to 0.75,0.24 front 
set arrow from 0.85,0.42 to 0.8,0.32 front
set arrow from 0.85,0.42 to 0.8,0.2 front  
set arrow from 0.3,1.11 to 0.25,0.96 front  

f(x) = a*x*x + b*x + c
#fit f(x) 'crit_points.dat' u 1:3 via a,b,c
a = -0.56
b = -0.59
c = 1.15
set yrange [0:1.2]

plot 'crit_points.dat' u 1:2 smooth bezier w filledcurves x1 lc rgb "#87CEFA" t '',\
     f(x) w filledcurves x2 lc rgb "#ffaaaa" t '',\
     '' u 1:2 w l dt 2 lw 2 lc -1 smooth bezier t '',\
     '' u 1:2 w p pt 4 lc -1 t '',\
     f(x) dt 1 lw 2 lc -1 t '',\
     '' u 1:3 w p pt 5 lc -1 t '',\
     'spheroid.dat' u ($1)*($1)*0.25:2 w l dt 3 lw 2 lc -1 t ''
     



