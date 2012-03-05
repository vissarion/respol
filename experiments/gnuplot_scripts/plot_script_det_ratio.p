# gnuplot script

set xlabel "Number of input points" font "Helvetica 30"
set ylabel "Det time (sec) / Overall time (sec)"
set key right top
#set logscale y 
set xrange [0:100]

plot "aug_new_results/random_d3_sparse.txt" using 3:($8/$7)\
 title "n=2" smooth unique with points, \
"aug_new_results/random_d4_sparse.txt" using 3:($8/$7)\
 title "n=3" smooth unique with points, \
"aug_new_results/random_d5_sparse.txt" using 3:($8/$7)\
 title "n=4" smooth unique with points
 
set terminal postscript eps color lw 1 "Helvetica" 22
#terminal postscript eps color lw 15 "Helvetica" 20
set output "my-plot.ps"
replot
