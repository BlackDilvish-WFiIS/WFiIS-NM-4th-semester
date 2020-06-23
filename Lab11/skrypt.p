set term jpeg size 800,600

set xlabel "x"
set ylabel "y"

set out "k8cz2.jpg"

plot "k8cz2.dat" u 1:2 w l t '|FFT|',\
              '' u 1:3 w l t 'prog dyskryminacji'
