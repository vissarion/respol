# gnuplot script

set xlabel "Number of input points" font "Helvetica 30"
set ylabel "Time (sec)"
set key left top
#set logscale y 
#set xrange [0:100]

plot "../new_results_W/random_d3_dense.txt" using 3:7\
 title "impl-dense" smooth unique with points, \
"../new_results_W/random_d3_ures_dense.txt" using 3:7\
 title "ures-dense" smooth unique with points, \
"../new_results_W/random_d3_sparse.txt" using 3:7\
 title "impl-sparse" smooth unique with points

 
set terminal postscript eps color lw 1 "Helvetica" 22
#terminal postscript eps color lw 15 "Helvetica" 20
set output "my-plot.ps"
replot
