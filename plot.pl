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

# Kmer plot
print << "EOF"
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
"<grep QUAL   $fn" using 4:5 with linespoints lw 2 lc 4 title "Qual"
EOF
;

# ACGT* -> ACGT* heatmap
# Figures are a percentage of the total base count.
open(my $fh, "<$fn") || die "$fn";
while(<$fh>) {
    next unless /^#\s+A\s+C/; # FIXME: output SUBST here too would help
    last;
}
my $mat = "$_";
$mat =~ s/#/x/;
my @m;
my $x = 0;
my $tot = 0.1;
while(<$fh>) {
    last unless /^SUBST/;
    chomp($_);
    my @F = split(/\s+/,$_);
    for (my $y=0; $y<=5; $y++) {
	$tot += $F[$y+2];
	$m[$x][$y]=$F[$y+2];
    }
    $x++;
}
my @bases = qw/A C G T N */;
for (my $x=0; $x<=5; $x++) {
    $mat .= "$bases[$x]\t";
    for (my $y=0; $y<=5; $y++) {
	my $f = $x==$y ? 0 : 100 * $m[$x][$y] / $tot;
	$mat .= "\t$f";
    }
    $mat .= "\n";
}
#print STDERR "$mat\n";
close($fh);

print << "EOF2"
\$mat << EOD
$mat
EOD
set terminal "pngcairo"
set output "$fn.subst.png"
set title "Error \%age matrix: $title" offset 0,0.5 font "Helvetica,12"
set xlabel "call"
set ylabel "real" offset -2
set size 0.95,1
set view map
#set cbtics format "%g"
set xrange [-0.5:5.5]
set yrange [-0.5:5.5]
splot '\$mat' matrix rowheaders columnheaders using 1:2:3 with image notitle
EOF2
;
