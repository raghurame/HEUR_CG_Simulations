set term qt font "DejaVu Math TeX Gyre, 20" size 1000,900

set colorsequence podo
set xtics nomirror
set ytics nomirror
set xlabel "Time"
set ylabel "coordination number"
set key spacing 1.2 right bottom

plot \
'./m0.1/all.lammpstrj.coordination' u ($1*2*0.1):2 w l lw 3 title "m_{core} = 10^{-1}", \
'./m0.01/all.lammpstrj.coordination' u ($1*2*0.1):2 w l lw 3 title "m_{core} = 10^{-2}", \
'./m0.001/all.lammpstrj.coordination' u ($1*2*0.1):2 w l lw 3 title "m_{core} = 10^{-3}", \
'./m0.0001/all.lammpstrj.coordination' u ($1*2*0.1):2 w l lw 3 title "m_{core} = 10^{-4}", \
'./m0.00001/all.lammpstrj.coordination' u ($1*2*0.1):2 w l lw 3 title "m_{core} = 10^{-5}"
