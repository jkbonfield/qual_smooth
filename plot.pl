#!/usr/bin/perl -w
use strict;

if (scalar(@ARGV) != 1) {
    print "Usage: plot.pl train.stats\n";
    exit 1;
}

my $fn = shift(@ARGV);
my $title = $fn;
$title =~ s#.*/##;
$title =~ s/\.stat$//;
$title =~ s/_/-/g;

#-----------------------------------------------------------------------------
# Kmer plot; overlapping
open(GP, "|gnuplot");
print GP << "EOF"
set terminal "pngcairo"
set output "$fn.kmer.png"
set title "KMERS: $title" font "Helvetica,12"
set xlabel "Measured Phred score"
set ylabel "Predicted Phred score"
set yrange [0:50]
a=rand(-1)
plot x lt 0 notitle, \\
"<grep MATCH2 $fn" using (\$6+rand(0)):(\$8+rand(0)):2 with labels hypertext point lt 7 ps 0.3 lc 1 title "Mismatch", \\
"<grep OVER2  $fn" using (\$6+rand(0)):(\$8+rand(0)):2 with labels hypertext point lt 7 ps 0.3 lc 8 title "Overcall", \\
"<grep UNDER2 $fn" using (\$6+rand(0)):(\$8+rand(0)):2 with labels hypertext point lt 7 ps 0.2 lc 20 title "Undercall", \\
"<grep ALL2   $fn" using (\$6+rand(0)):(\$8+rand(0)):2 with labels hypertext point lt 7 ps 0.1 lc 11 title "Combined", \\
"<grep QUAL   $fn | grep -v -- -2147" using 5:4 with linespoints lw 2 lc 4 title "Qual"
EOF
;
close(GP);

#-----------------------------------------------------------------------------
# Kmer multi-plot
open(GP, "|gnuplot");
print GP << "EOFm"
set terminal pngcairo size 800,600
set output "$fn.kmer4.png"

yl="Predicted Phred score"
xl="Measured Phred score"
set yrange [0:\`awk '/(ALL|UNDER|MATCH)2/ {v=v<\$8?\$8:v} END {print 5*(int((v+9.99)/5))}' $fn\`]
set xrange [0:\`awk '/(ALL|UNDER|MATCH)2/ {v=v<\$6?\$6:v} END {print 5*(int((v+5.99)/5))}' $fn\`]
set xtics offset 0,0.3

#set xrange [0:30]
#set yrange [0:25]

a=rand(-1)

set grid
set multiplot
unset title

# Top left: mismatch
set label "KMERS: $title" font "Helvetica,12" at graph 1.1,1.2 centre
unset xlabel
set ylabel yl offset 1.5
set size 0.53,0.46
set origin 0,0.46
plot "<grep MATCH2 $fn" using (\$6+rand(0)):(\$8+rand(0)):2 with labels hypertext point lt 7 ps 0.2 lc 1 title "Mismatch"
unset ylabel
unset label

# Top right: overcalls
set size 0.5,0.46
set origin 0.5,0.46
plot "<grep OVER2  $fn" using (\$6+rand(0)):(\$8+rand(0)):2 with labels hypertext point lt 7 ps 0.2 lc 8 title "Overcall"

# Bottom left: undercalls
set size 0.53,0.48
set xlabel xl offset 0,0.8
set ylabel yl offset 1.5
set origin 0,0.01
plot "<grep UNDER2 $fn" using (\$6+rand(0)):(\$8+rand(0)):2 with labels hypertext point lt 7 ps 0.2 lc 20 title "Undercall"
unset ylabel

# Bottom right: Combined
set size 0.5,0.48
set origin 0.5,0.01
plot "<grep ALL2   $fn" using (\$6+rand(0)):(\$8+rand(0)):2 with labels hypertext point lt 7 ps 0.2 lc 6 title "Combined"
EOFm
;
close(GP);

#-----------------------------------------------------------------------------
# ACGT* -> ACGT* heatmap
# Figures are a percentage of the total base count.
foreach my $sec (qw/SUBST_HQ SUBST_ALL/) {
    open(my $fh, "<$fn") || die "$fn";
    while(<$fh>) {
	next unless /^$sec+\s+#\s+A\s+C/;
	last;
    }
    my $mat = "$_";
    $mat =~ s/^$sec\s+#/x/;
    my @m;
    my $x = 0;
    my $tot = 0.1;
    while(<$fh>) {
	last unless /^SUBST/;
	chomp($_);
	my @F = split(/\s+/,$_);
	for (my $y=0; $y<=5; $y++) {
	    if ($F[1] ne "N") {
		$tot += $F[$y+2];
		$m[$x][$y]=$F[$y+2];
	    } else {
		# Skip N as ref-N isn't a real thing
		$m[$x][$y] = 0;
	    }
	}
	$x++;
    }
    my @bases = qw/A C G T N */;
    for (my $x=0; $x<=5; $x++) {
	next if ($x == 4); # skip N row
	$mat .= "$bases[$x]\t";
	for (my $y=0; $y<=5; $y++) {
	    my $f = $x==$y ? 0 : 100 * $m[$x][$y] / $tot;
	    $mat .= "\t$f";
	}
	$mat .= "\n";
    }
    close($fh);

    my $subst = lc($sec);
    my $subst2 = $subst; $subst2 =~ s/_/ /g;

    open(GP, "|gnuplot");
    print GP << "EOF2"
\$mat << EOD
$mat
EOD
set terminal "pngcairo"
set output "$fn.$subst.png"
set title "$subst2 error \%age matrix: $title" offset 0,0.5 font "Helvetica,12"
set palette rgb 34,35,36
#set palette rgb 7,3,23
set xlabel "call"
set ylabel "cons" offset -2
set size 0.95,1
set view map
#set cbtics format "%g"
set xrange [-0.5:5.5]
set yrange [-0.5:4.5]
splot '\$mat' matrix rowheaders columnheaders using 1:2:3 with image notitle
EOF2
;
    close(GP);
}
