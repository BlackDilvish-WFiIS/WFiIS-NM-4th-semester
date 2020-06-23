set term jpeg size 800,600

set xlabel 'x'
set ylabel 'v'

set out 'V50.jpg'

plot 'dane.dat' u 1:2 w l lt -1 lw 2 t 'V wyznaczone',\
       ''       u 1:3 w l lt  2 lw 2 t 'V teoretyczne'
