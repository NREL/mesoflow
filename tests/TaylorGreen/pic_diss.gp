set terminal pngcairo dashed size 800,500
set output "diss_compare.png"
set termoption font "Helvetica,18"
set key spacing 1.2
set lmargin 11
set rmargin 5
set bmargin 3.5
set xlabel "Time"
set ylabel "dissipation rate"

plot 'tg_diss.dat' u ($1/3.18):(-$2*3.18/0.01) w l lw 3 title "512^3",\
'spectral_Re1600_512.gdiag' u 1:3 w lp lw 1 pt 6 ps 1 pi 50 lc 7 title "HO CFD workshop"
