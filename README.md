Building
--------

Edit the Makefile to specify the location of htslib and htscodecs
directories.

Note, these are currently coded up to point to built source trees, but
they could also point to official installs in e.g. /usr/lib and /usr/include.

Then build with:

    make

Running
-------

    ./qual_smooth -O cram,reference=$HREF38 -P5 -B10 PB.bam PB_qs.cram


Options
-------

Multiple options can be combined, so for example we can smooth using
P-block for low quality values while R-block can modify the much
larger values which are being shifted by more than the P-block
parameter permits.  Finally the smoothed values could be quantised
with -B or -b binning.

-P INT   Use the P-block algorithm to move all runs of quality values
         within +/- INT to their mid point.

         This can sometimes help avoid flipping qualities around when
         they're close to a bin boundary.

-R FLOAT Use the R-block algorithm to move all runs of quality values
         within +/- Max/FLOAT and FLOAT/Min where Max and Min are the
         largest and smallest quality values within the window.

-B INT   Break into discrete bins of width INT, using the mid points.
         So -B 10 will set 1-10 as 5, 11-20 as 15, etc. (Or maybe 0-9,
         I forget now, but see below for a better alternative...)

-b LIST  Specify bins explicitly as a list in the form L-H=M[,L-H=M]*.
         Eg the colord profile would be 0-6=3,7-13=10,14-25=18,26-92=35

         I dislike changing qual 0 as this has a special meaning too,
         so I'd probably go with 1-6.  Anything not listed isn't
         changed, so Q93 is still Q93.

-k LIST  Keep quality INT as-is. We never create this quality nor
         modify this quality (unless explicitly requested via -b).

         We can specify it multiple times.  (It probably ought to be a
         list, but for now it's one at a time.)
         By default we preserve 0 and 93.

         LIST may be single digits ("10"), comma-separated lists
         ("10,11,12,20"), and include ranges ("10-12,20").

-X profile
         Use a set of predefined standard options.  These are:

         pbccs      -P5 -R1.5 -B10
         illuimna   -P3 -R1.2 -B8
         hiseq      -P0 -R0 -b0-0=0,1-1=1,2-9=6,10-19=15,20-24=22,25-29=27,\
                        30-34=33,35-39=37,40-99=40
         novaseq    -P0 -R0 -b0-6=2,7-13=12,14-24=18,25-99=36

-t INT   Use a pool of INT threads for the decoding and encoding.

-I FMT
-O FMT   Specifies input and output format, along with format specific
         options.  This is the same as samtools view.  Eg:
         "-O bam", or "-O cram,version=3.1,small,level=7,reference=$HREF38".

         Note CRAM 3.1 is draft so don't use for archival, but it's a
         suitable thing to evaluate when surveying the landscape of
         tools as I'm not expecting major changes. (It's stuck waiting
         for Broad to pull their finger out on the 2nd implementation
         before we can ratify it as an official GA4GH standard.)

TODO: Consider preservation around key sequence motifs as well
(similar to crumble, but purely kmer based).  For example homopolymers
are often error prone, so we may wish to be more cautious on changing
those values.
