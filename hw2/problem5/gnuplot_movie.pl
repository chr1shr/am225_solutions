#!/usr/bin/perl
use Getopt::Std;
getopts("bdfg:hn:p:t");

# Print help message if -h flag is given
if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename>\n";
	print "\nOptions:\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-f       (Just make the frames, don't make a movie)\n";
	print "-g <n>   (Use an alternative gnuplot header)\n";
	print "-h       (Print this information)\n";
	print "-n <n>   (Render up to this number of files)\n";
	print "-p <n>   (Override default number of processes to use)\n";
	print "-t       (Don't add the Voronoi tessellation)\n";
	exit 0;
}

# Check command line options
die "Need either two arguments" unless $#ARGV==0;

# Determine the number of processors
if(!defined $opt_p) {
	$uname=`uname`;
	if($uname=~/Linux/) {
		$nodes=`awk '/^processor/ {++n} END {print n+1}' /proc/cpuinfo`;
		chomp $nodes;
	} elsif($uname=~/Darwin/) {
		$nodes=`sysctl -n hw.ncpu`;
		chomp $nodes;
	} else {
		$nodes=4;
	}
} else {
	$nodes=$opt_p;
}

# Prepare output directory
$e=$ARGV[0];
$o=$e.".frames";
$o=~s/\.odr//;
mkdir $o unless -e $o;

# Loop over available files, create a gnuplot script, and render them
$a=0;
$P=0;$queue=$nodes==1?1:0;

# Declare region size
$ax=-1.5;$bx=1.5;$ay=-1.5;$by=1.5;

# Read in the gnuplot header
$maxcol=0;
$gph="gp_headers/t".($opt_g?$opt_g:($opt_b?"3":"1")).".gnuplot";
open A,"$gph" or die "Can't find gnuplot header\n";
$gpn=0;
while(<A>) {
	s/XMIN/$ax/;
	s/XMAX/$bx/;
	s/YMIN/$ay/;
	s/YMAX/$by/;
	$gp[$gpn++]=$_;
}
close A;$gpn--;

while(1) {

	# Check for presence of required file
	$infile="$e/fr.$a";
	last unless -e $infile;

	# Check for termination condition
	last if defined $opt_n && $a>$opt_n;
	$of=sprintf "$o\/fr_%04d.png",$a;

	# If the -d flag is given, then skip making the image if one already
	# exists
	if ($opt_d && -e $of && -M $infile > -M $of ) {
		print "$a (skipped)\n";
		$a++;$infile="$e\/pts.$a";
		next;
	}

	open D,"<$infile";
	$ts=<D>;
	$ts=~s/# Time: //;
	$ts=~s/\n//;
	$ts=sprintf "%0.4f",1.*$ts;
	close $D;

	# Assemble the gnuplot script for this file
	print "Frame $a (thread $P)\n";
	open B,">$e/temp$P.gnuplot";
	foreach $i (0..$gpn) {
		$_=$gp[$i];
		s/OUTFILE/$of/;
		s/FRAME/$a/;
		s/INFILE/$infile/g;
		s/TIME/$ts/g;
		print B;
	}
	close B;

	# Send the Gnuplot file to a forked process
	exec $ex."gnuplot $e/temp$P.gnuplot 2>/dev/null" if ($pid[$P]=fork)==0;

	# Wait for one of the forked jobs to finish
	if ($queue) {
		$piddone=wait;$P=0;
		$P++ while $piddone!=$pid[$P] && $P<$nodes;
		die "PID return error!" if $P>=$nodes;
	} else {
		$P++;$queue=1 if $P>=$nodes-1;
	}

	# Increment file counter and set new filename
	$a++;
}

# Wait for all remaining forked jobs to complete
wait foreach 1..($queue?$nodes:$h-1);

# Make a movie of the output
unless ($opt_f) {
	$mn="${e}_pts";
	$mn=~s/\.odr_pts//;
	$frate=30;
	$uname=`uname`;
	if ($uname=~/Linux/) {
		unlink "$mn.mpg";
		system "ffmpeg -r $frate -y -i $o/fr_%4d.png -vb 20M $mn.mpg"
	} elsif($uname=~/Darwin/) {
		unlink "$mn.mov";
		system "qt_export --sequencerate=$frate $o/fr_0001.png --loadsettings=../../misc/qtprefs/qt --replacefile $mn.mov";
	}
	#system "rm -rf $o";
}
