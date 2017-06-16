set terminal postscript eps enhanced colour font ",20"
set output "ERTrates.eps"
set xrange [0:1.1]
set size square
set xlabel "k{/Symbol g} / {/Symbol e}_0 E^2"
set ylabel "Growth rate, N = (n^2{/Symbol r}/k^2{/Symbol e}_0E^2)^{1/2}"
set key off

plot 'ERTrates.dat' w p lt 7 lc 'black' t 'Simulation',\
     7.6632*sqrt(1-x) w l lt 1 lw 2 lc rgb "black" t 'Invisid PT'
