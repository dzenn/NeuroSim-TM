set title "Network activity (avg time 2 ms)"
set yrange [0:1]
set xrange [0:10000]
set xlabel "time, ms"
set ylabel "net activity, %"
plot "activity.txt" with line lw 1 lt 1