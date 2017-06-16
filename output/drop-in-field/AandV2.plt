set terminal postscript eps enhanced colour font ",20"
set output "AandV2.eps"
set xrange [0:20]
set yrange [1:4]
set y2range [1:1.0045]
set ytics nomirror 1,0.5
set y2tics nomirror 1,0.0005
set size square
set key bottom right
set xlabel "t / {/Symbol t}"
set ylabel "A = a/b" offset 2
set y2label "(V+{/Symbol D}V) / V" offset -2

plot 'fluid_params_F1.2.dat' u ($1)/2:5 axes x1y1 w lp lt 7 t 'Aspect ratio',\
     '' u ($1)/2:($2)/4.18879 axes x1y2 w lp lt 6 t 'Volume fraction'
