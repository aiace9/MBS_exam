set term postscript portrait colour
set size square
set out 'Energy_evolution.ps'
set grid
show grid

set title 'Energy'
set xlabel "step"
set ylabel "E"

plot 'kin_pot_tot.dat' u 1 : 4 w l t ''

set term postscript portrait colour
set size square
set out 'Temperature_evolution.ps'
set grid
show grid

set title 'Temperature'
set xlabel "step"
set ylabel "T*"

plot 'T_P.dat' u 1 : 3 w l t ''

set term postscript portrait colour
set size square
set out 'gr.ps'
set grid
show grid

set title 'Pair correlation function'
set xlabel "g(r)"
set ylabel "E"

plot 'ist.dat' u 1 : ($2*1000) w l t ''
