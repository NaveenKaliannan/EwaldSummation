set terminal png enhanced
set terminal postscript enhanced
set terminal postscript eps enhanced
set term postscript enhanced color
set output "colorindex.ps"
set border lw 4
set grid
set lmargin at screen 0.19
set style line 1 lt 0 lc rgb "red" lw 3
set style line 2 lt 2 lw 3
set style line 3 lt 3 lw 3
plot "../data/MadelungConstant_rutile.dat" using 1:2  title "" with linespoints ls 1 pt 1 ps 2 
set title "TiO_2 (rutile)" font "arial,24"
set xtic font "arial,15"
set xlabel "k (Å)"  font "arial,18"
set ylabel "Madelung constant" font "arial,18"
set xrange [0:3]
set yrange [1.5:2.8]
set terminal png font arial 20 size 1024,768
set output "../graph/rutile.png"
replot