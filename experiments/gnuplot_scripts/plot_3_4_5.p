# gnuplot script

set xlabel "Number of input points" font "Helvetica 30"
set ylabel "cpu time (sec)"
set key left top
set logscale y 
set xrange [10:45]

plot "new_results_W/random_d3_sparse_10_45.txt" using 3:7\
 title "n=2" smooth unique with lp, \
"new_results_W/random_d4_sparse.txt" using 3:7\
 title "n=3" smooth unique with lp, \
"new_results_W/random_d5_sparse.txt" using 3:7\
 title "n=4" smooth unique with lp;

 
set terminal postscript eps lw 1 "Helvetica" 22
#terminal postscript eps color lw 15 "Helvetica" 20
set output "my-plot.ps"
replot
