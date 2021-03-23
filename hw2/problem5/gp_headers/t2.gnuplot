set size ratio -1
xmin=-1.5
xmax=1.5
ymin=-1.5
ymax=1.5
del=2*pi/3.
rgb(x)=int(255*(0.4+0.4*sin(x))+0.5)+256*int(255*(0.4+0.4*sin(x+del)))+65536*int(255*(0.4+0.4*sin(x-del))+0.5)
unset key
set xlabel '$x$' offset 0,0.75
set ylabel '$y$' offset 3.25 rotate by 0
set term epslatex size 3in,3in header "\\scriptsize" color
if (!exists("ifile")) ifile='sw.0/fr.0'
if (!exists("ofile")) ofile='out.tex'
if (!exists("time")) time='0.0'
set output ofile
set label 3 at -1.4,1.3 "\\textcolor{magenta!50!blue}{T = ".time."}" tc rgbcolor "#0000ff"
set rmargin 0.2
set tmargin 0.2
plot [xmin:xmax] [ymin:ymax] ifile u 1:2:(rgb($3)) w p lc rgb variable pt 7 ps 0.5
unset label 3
set output
