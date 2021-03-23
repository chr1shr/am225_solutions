#!/usr/bin/env bash

for mod in {0..2}
do
	for fr in 20 40 100 400
	do
		gnuplot -e "ifile='sw.$mod/fr.$fr'" -e "ofile='snaps/sn_$mod""_$fr.tex'" -e "time='$(($fr/2)).00'" gp_headers/t2.gnuplot
	done
done
