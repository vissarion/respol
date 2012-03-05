# gnuplot script

set xlabel "Number of input points" font "Helvetica 30"
set ylabel "time (sec)"
set key right bottom
set logscale y 
set yrange [0.01:1000]

plot "./triang_experiments/random_cube_hash.txt" using 2:3\
 title "hash"  with points, \
 "/workspace/Triangulation-higher_dimensions-odevil_shornus/\
Triangulation/examples/Triangulation/experiments/random_cube_horn.txt"\
 using 2:3 title "no-hash" with points

set terminal postscript eps color lw 1 "Helvetica" 22
#terminal postscript eps color lw 15 "Helvetica" 20
set output "my-plot.ps"
replot
