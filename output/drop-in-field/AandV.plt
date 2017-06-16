set terminal postscript eps enhanced colour font ",20"
set output "AandV.eps"
set xrange [0:100]
set yrange [1:1.35]
set y2range [0.9994:1.0008]
set ytics nomirror 1,0.05
set y2tics nomirror 0.9994,0.0002
set size square
set xlabel "t / {/Symbol t}"
set ylabel "A = a/b" offset 2
set y2label "(V+{/Symbol D}V) / V" offset -2

plot 'fluid_params_F0.6.dat' u ($1)/2:5 axes x1y1 w lp lt 7 t 'Aspect ratio',\
     '' u ($1)/2:($2)/4.18879 axes x1y2 w lp lt 6 t 'Volume fraction'
