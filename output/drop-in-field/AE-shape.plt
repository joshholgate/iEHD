set terminal postscript eps enhanced colour font ",20"
set output "AE-shape.eps"
set xrange [0:1.7]
set yrange [1:2]
set size square
set xlabel "E(4{/Symbol p}{/Symbol e}_0r/{/Symbol g})^{1/2}"
set ylabel "A = a/b" offset 1
set key off

plot 'Avalues.dat' u(1.625*sqrt($1)):2:3 w e pt 7 lc 'black',\
     'AE-shape-spheroids.dat' u (1.625*sqrt($1)):2 w l lt 1 lw 2 lc 'black'
