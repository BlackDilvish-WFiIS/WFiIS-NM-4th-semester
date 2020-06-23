set term jpeg size 800,600

set xlabel "alfa"
set ylabel "Wartosc wlasna"

set out "najmniejsze_wartosci.jpg"

set xrange [0:100]
set yrange [0:0.3]

plot "najmniejsze_wartosci.dat" u 1:2 w l t 'Wartosc wlasna 1',\
                                '' u 1:3 w l t 'Wartosc wlasna 2',\
                                '' u 1:4 w l t 'Wartosc wlasna 3',\
                                '' u 1:5 w l t 'Wartosc wlasna 4',\
                                '' u 1:6 w l t 'Wartosc wlasna 5',\
                                '' u 1:7 w l t 'Wartosc wlasna 6'