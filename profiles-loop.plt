set terminal png size 450,600
t = t + 1
datafile = sprintf('t%01.0f_gridvalues.dat',t)
outfile = sprintf('t%01.0f_gridvalues.png',t)
heading = sprintf('Time = %d',t)

set contour base
set cntrparam levels disc 0
unset surface
set table 'cont.dat'
splot datafile u 1:2:6, datafile u ($1)*(-1):2:6
unset table

set output outfile
#set title heading
set size ratio 2

set xrange [-2:2]
set yrange [0:8]
set key off
set border lw 1.5
set palette defined (-1 "white", 1 "gray") 
set palette maxcolors 2
unset colorbox

set cbrange [-0.5:0.5]
plot datafile u 1:2:6 w image,\
     datafile u ($1)*(-1):2:6 w image,\
     'cont.dat' w l lt -1 lw 1.5
if(t<end_time) reread;



# The follow lines of code might be found useful for plotting velocity, pressure and electric field profiles
#set palette maxcolors 5
#set palette defined (-1 "#00008B", -0.8 "#0000CD", -0.6 "#4169E1", -0.4 "#1E90FF", -0.2 "#87CEEB", 0 "white")
#set palette defined (-1 "white", -0.8 "#87CEEB", -0.6 "#1E90FF", -0.4 "#4169E1", -0.2 "#0000CD", -0.001 "#00008B", 0 "white")
#plot datafile u 1:2:(($6 < 0)  ? (sqrt(($8)*($8)+($9)*($9))) : 1/0) w image,\
#     datafile u ($1)*(-1):2:(($6 < 0)  ? (sqrt(($8)*($8)+($9)*($9))) : 1/0) w image,\
#     datafile u 1:2:(($6 > 0)  ? ($3)*10.0 : 1/0):(($6 > 0)  ? ($4)*10.0 : 1/0) every 5:5:1 w vectors head filled lt 1,\
#     datafile u ($1)*(-1):2:(($6 > 0)  ? ($3)*(-10.0) : 1/0):(($6 > 0)  ? ($4)*10.0 : 1/0) every 5:5::4 w vectors head filled lt 1,\
#     'cont.dat' w l lt -1 lw 1
