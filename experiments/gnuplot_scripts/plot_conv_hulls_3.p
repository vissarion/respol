# gnuplot script

set xlabel "Number of points" font "Helvetica 30"
set ylabel "time"
set key left top
set logscale y 
#set yrange [0:.2]

f(x) = a1*x**b1            # define the function to be fit
a1 = 0.5; b1 = 0.2;            # initial guess for a1 and b1
fit f(x) "results_convex_hulls/random_d3.txt" using 7:10 via a1, b1

# This command fits your fitting function to the 
# x and y values contained in columns 1 and 2 
# of the data file "row-data.dat". The fitting
# parameters are listed after the "via" keyword.
#fit f(x) "results_convex_hulls/random_d3.txt" using 7:10 via a,b,c,x0

plot "results_convex_hulls/random_d3.txt" using 7:10\
 title "cgal" smooth unique with points, \
    "results_convex_hulls/random_d3.txt" using 7:17\
 title "bb" smooth unique with points, \
    "results_convex_hulls/random_d3.txt" using 7:21\
 title "cdd" smooth unique with points, \
    "results_convex_hulls/random_d3.txt" using 7:25\
 title "lrs" smooth unique with points,\
 "results_convex_hulls/random_d3_2.txt" using 7:10\
 notitle smooth unique with points pt 1 lc 1, \
    "results_convex_hulls/random_d3_2.txt" using 7:17\
 notitle smooth unique with points pt 2 lc 2, \
    "results_convex_hulls/random_d3_2.txt" using 7:21\
 notitle smooth unique with points pt 3 lc 3, \
    "results_convex_hulls/random_d3_2.txt" using 7:25\
 notitle smooth unique with points  pt 4 lc 4, \
 "results_convex_hulls/random_d3.txt" using 7:11\
 title "cgal_off" smooth unique with points  pt 5 lc 5,\
 "results_convex_hulls/random_d3_2.txt" using 7:11\
 notitle smooth unique with points pt 5 lc 5;
 #f(x) lc 1;

set terminal postscript eps color lw 1 "Helvetica" 22
#terminal postscript eps color lw 15 "Helvetica" 20
set output "my-plot.ps"
replot
