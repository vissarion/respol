# gnuplot script

set xlabel "Number of points" font "Helvetica 30"
set ylabel "time"
set key right bottom
set logscale y 
#set yrange [0:10]

plot "results_convex_hulls/random_d4.txt" using 7:10\
 title "cgal" smooth unique with points, \
    "results_convex_hulls/random_d4.txt" using 7:17\
 title "bb" smooth unique with points, \
    "results_convex_hulls/random_d4.txt" using 7:21\
 title "cdd" smooth unique with points, \
    "results_convex_hulls/random_d4.txt" using 7:25\
 title "lrs" smooth unique with points, \
  "results_convex_hulls/random_d4.txt" using 7:11\
 title "cgal_off" smooth unique with points, \
 "results_convex_hulls/random_d4_2.txt" using 7:10\
 notitle smooth unique with points pt 1 lc 1, \
    "results_convex_hulls/random_d4_2.txt" using 7:17\
 notitle smooth unique with points pt 2 lc 2, \
    "results_convex_hulls/random_d4_2.txt" using 7:21\
 notitle smooth unique with points pt 4 lc 4, \
  "results_convex_hulls/random_d4_2.txt" using 7:11\
 notitle smooth unique with points pt 5 lc 5;

set terminal postscript eps color lw 1 "Helvetica" 22
#terminal postscript eps color lw 15 "Helvetica" 20
set output "my-plot.ps"
replot
