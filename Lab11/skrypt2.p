set term jpeg size 800,600

set xlabel "x"
set ylabel "y"

set out "k12cz4.jpg"
set pointsize 1

plot "k12cz3.dat" u 1:3 w p pt 5 ps 1 t 'niezaburzony',\
              '' u 1:4 w l t 'odszumiony'
