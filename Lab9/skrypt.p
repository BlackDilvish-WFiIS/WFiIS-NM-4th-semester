set term jpeg size 800,600

set xlabel "x"
set ylabel "y"

set out "MN3.jpg"

plot "MN3.dat" u 1:2 w l t 'f(x)',\
              '' u 1:3 w l t 'R(x)'
