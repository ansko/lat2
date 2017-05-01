set terminal png size 2048, 768
set output 'graph-soft.png'
#set xtics 0.25 nomirror
set xtics 1 nomirror
#set yrange [0 to 1000000]
#set xrange [-0.5 to 1]
#set x2range [0 to 10000]
set ytics nomirror
set pointsize 3
#const=1
#set trange [1:4]

set ylabel "Увеличение модуля Юнга композита" font "Times-Roman,20"
set xlabel "Доля интеркалированных частиц" font "Times-Roman,20"



plot 'log' u 1:2 w l lc rgb "red" 
