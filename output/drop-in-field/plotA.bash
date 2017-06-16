#!/bin/bash
> Avalues.dat
for E in 0 0.1 0.2 0.3 0.4 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.825 0.85 0.875 0.9 0.925 0.935 0.95 1
do
	DATA="fluid_params_F""$E"".dat"
	gnuplot -persist <<-EOFMarker
		stats "$DATA" u (\$1>100 ? \$5 : 1/0) name "A"
		set print "Avalues.dat" append
		print "$E ", A_mean, A_stddev
	EOFMarker
done
