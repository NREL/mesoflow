set terminal png
set output "cpcompare.png"
set termoption font "Helvetica,20"
set key spacing 1.2
#set lmargin 8
#set rmargin 3
#set bmargin 4
set xlabel "angle (deg.)"
set ylabel "Cp"
set yrange [-1.5:1.8]
set xrange [0:180]
set xtics (0,40,80,120,160)
plot 'cpcyl40.dat' u ($1*180.0/3.141959265):2 w lp ps 1 pt 7 title "Kou et al., JCP, 448, 2022",\
'presdata_000' u ($1*180.0/3.14159265):3 w lp ps 1 pt 5 title "current work"
