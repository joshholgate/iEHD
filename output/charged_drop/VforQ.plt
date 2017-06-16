set terminal postscript eps enhanced colour font ",20"
set output "VforQ.eps"
set xrange [0:150]
#set yrange [0.9985:1.0015]
set size square
set xlabel "t / {/Symbol t}"
set ylabel "(V+{/Symbol D}V) / V" offset 2
set key bottom left title "Q/Q_c"

plot 'fluid_params_Q1.dat' u 1:($2)/4.18876 every 3 w lp lw 1.5 lc rgb "#d95319" lt 6 t '1.00',\
     'fluid_params_Q1.01.dat' u 1:($2)/4.18876 every 3 w lp lw 1.5 lc rgb "#edb120" lt 7 t '1.01',\
     'fluid_params_Q1.02.dat' u 1:($2)/4.18876 every 3 w lp lw 1.5 lc rgb "#7e2f8e" pt 4 t '1.02',\
     'fluid_params_Q1.03.dat' u 1:($2)/4.18876 every 3 w lp lw 1.5 lc rgb "#77ac30" pt 5 t '1.03',\
     'fluid_params_Q1.04.dat' u 1:($2)/4.18876 every 3 w lp lw 1.5 lc rgb "#4dbeee" pt 8 t '1.04',\
     'fluid_params_Q1.05.dat' u 1:($2)/4.18876 every 3 w lp lw 1.5 lc rgb "#a2142f" pt 9 t '1.05'
