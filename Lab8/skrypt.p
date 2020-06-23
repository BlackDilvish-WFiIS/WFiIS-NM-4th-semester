set term jpeg size 800,600

set xlabel "x"
set ylabel "y"

set out "f2n21.jpg"

plot "f2n21.dat" u 1:2 w l t 'f2(x)',\
              '' u 1:3 w l t 's(x)'
