set term jpeg size 800,600

set xlabel "nr iteracji"
set ylabel "modul roznicy rozwiazania dokladnego i przyblizonego"

set logscale y

set out "g.jpg"

plot "zlotyG.dat" u 1:2 w l t 'metoda zlotego podzialu',\
     "rowne3G.dat" u 1:2 w l t 'podzial na 3 rowne odcinki'
