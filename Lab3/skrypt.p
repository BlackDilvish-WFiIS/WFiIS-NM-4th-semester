set term jpeg size 800,600

set xlabel 'k'
set ylabel 'x'
set out 'Temp.jpg'

plot 'Double1.dat' u 1:4 w l lt -1 lw 2 t 'x'
