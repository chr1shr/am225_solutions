set size ratio -1
xmin=XMIN
xmax=XMAX
ymin=YMIN
ymax=YMAX
del=2*pi/3.
rgb(x)=int(255*(0.4+0.4*sin(x))+0.5)+256*int(255*(0.4+0.4*sin(x+del)))+65536*int(255*(0.4+0.4*sin(x-del))+0.5)
unset key
set xlabel 'x' tc rgb 'white'
set ylabel 'y' tc rgb 'white'
set border lc rgb 'white'
set term pngcairo size 870,870 background rgb 'black' lw 3 font 'Helvetica-Bold'
set lmargin at screen 0.09
set rmargin at screen 0.98
set tmargin at screen 0.975
set bmargin at screen 0.07
set output 'OUTFILE'
#set label 3 at screen 0.08,0.025 "Frame FRAME" tc rgbcolor "#880088"
set label 3 at screen 0.08,0.0125 "T=TIME" tc rgbcolor "#5c98ff" font 'Helvetica-Bold'
plot [xmin:xmax] [ymin:ymax] 'INFILE' u 1:2:(rgb($3)) w p lc rgb variable pt 7 ps 0.5
unset label 3
set output
