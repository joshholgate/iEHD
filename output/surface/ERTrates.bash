#!/bin/bash
> ERTrates.dat
for E in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.96 0.97 0.98
do
	DATA="fluid_params_F""$E"".dat"
	gnuplot -persist <<-EOFMarker
		f(x) = a + b*x
		fit f(x) "$DATA" u (\$1)/100:(\$3 >=0.03 && \$3 <=0.08  ? (log(\$3)) : 1/0) via a,b
		set print "ERTrates.dat" append
		print "$E ", b
	EOFMarker
done
