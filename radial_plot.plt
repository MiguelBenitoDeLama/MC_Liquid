 set nokey
 set grid
 set xlabel "r"
 set ylabel "g(r)"
 set title 'Radial Distribution function g(r)'
 plot "radial.dat" u 1:2 w l
