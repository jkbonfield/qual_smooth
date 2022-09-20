#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>

#ifndef DEBUG
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#endif

#include "entropy.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define ABS(a)   ((a)>0?(a):-(a))

static int qpreserve[256] = {0};
static uint8_t qmap[256] = {0};

static void gen_map(int blevel) {
    int i, bin = blevel;
    for (i = 0; i < 256; i++) {
	if (qpreserve[i]) {
	    qmap[i] = i;
	    continue;
	}
	int j, k;
	for (j = 0; j < bin && i+j < 256 && !qpreserve[i+j]; j++);
	// bin i..j-1 together
	for (k = i; k < i+j; k++)
	    qmap[k] = i+j/2;
	i += j-1;
    }
}

// Format is e.g. 0-6=3,7-13=10,14-25=18,26-92=35
// It starts with a 1:1 mapping, so anything absent is left as-is
static int parse_map(char *str) {
    for (int i = 0; i < 255; i++)
	qmap[i] = i;

    for(;;) {
	while (*str == ',')
	    str++;
	if (*str && !isdigit(*str)) {
	    fprintf(stderr, "Incorrect format '%s': should be low-high=val\n",
		    str);
	    return 1;
	}

	int low, high, mid, n, ret;
	ret = sscanf(str, "%d-%d=%d%n", &low, &high, &mid, &n);
	if (ret <= 0) {
	    break;
	} else if (ret != 3) {
	    fprintf(stderr, "Incorrect format '%s': should be low-high=val\n",
		    str);
	    return 1;
	}
	str += n;

	for (int i = MAX(0,low); i <= MIN(255, high); i++)
	    qmap[i] = mid;
    }

//    for (int i = 0; i <= 255; i++)
//	printf("%d\t%d\n", i, qmap[i]);

    return 0;
}

static inline void pblock(uint8_t *qual, int len, int level, uint8_t *qmap) {
    int i, j, qmin = INT_MAX, qmax = INT_MIN;
    int last_qmin = 0, last_qmax = 0, mid;

    level *= 2;

    for (i = j = 0; i < len; i++) {
	if (qmin > qual[i])
	    qmin = qual[i];
	if (qmax < qual[i])
	    qmax = qual[i];
	if (qmax - qmin > level || qpreserve[qual[i]]) {
	    mid = (last_qmin + last_qmax) / 2;
	    memset(qual+j, qmap ? qmap[mid] : mid, i-j);
	    while (i < len && qpreserve[qual[i]])
		i++;
	    qmin = qmax = qual[i];
	    j = i;
	}
	last_qmin = qmin;
	last_qmax = qmax;
    }

    mid = (last_qmin + last_qmax) / 2;
    memset(qual+j, qmap ? qmap[mid] : mid, i-j);
}

static inline void rblock(uint8_t *qual, int len, double level, uint8_t *qmap)
{
    int i, j, qmin = INT_MAX, qmax = INT_MIN;
    int last_qmin = 0, last_qmax = 0;
    double last_r = 0;

    for (i = j = 0; i < len; i++) {
	if (qmin > qual[i])
	    qmin = qual[i];
	if (qmax < qual[i])
	    qmax = qual[i];

	// Pick representation r such that max/r = x and r/min = x,
	// and then see if max/r or r/min <= level.
	//
	//    (max/r) / (r/min) = x/x.
	// => (max*min) / (r*r) = 1
	// => r = sqrt(max*min)
	double r;
	if ((r = sqrt(qmax * qmin)) >= level*qmin || qpreserve[qual[i]]) {
	    int ri = (int)last_r+0.5;
	    memset(qual+j, qmap ? qmap[ri] : ri, i-j);
	    while (i < len && qpreserve[qual[i]])
		i++;
	    qmin = qmax = qual[i];
	    r = qual[i];
	    j = i;
	}
	last_qmin = qmin;
	last_qmax = qmax;
	last_r = r;
    }

    int ri = sqrt(last_qmax * last_qmin)+0.5;
    memset(qual+j, qmap ? qmap[ri] : ri, i-j);
}

static inline void qbin(uint8_t *qual, int len, uint8_t *map) {
    for (int i = 0; i < len; i++)
	qual[i] = map[qual[i]];
}

static void sim(uint8_t *q1, uint8_t *q2, int len, long *d1, long *d2) {
    int i;
    long diff1 = 0, diff2 = 0;
    for (i = 0; i < len; i++) {
	diff1 += ABS(q1[i]-q2[i]);
	diff2 += (q1[i]-q2[i]) * (q1[i]-q2[i]);
    }
    *d1 = diff1;
    *d2 = diff2;
}

/* Taylor (deg 3) implementation of the log */
static inline double fast_log2(double val) {
    // FP representation is exponent & mantissa, where
    // value = 2^E * M.
    // Hence log2(value) = log2(2^E * M)
    //                   = log2(2^E)+ log2(M)
    //                   =        E + log2(M)
    union { double d; uint64_t x; } u = {val};
    const int E = ((u.x >> 52) & 2047) - 1024; // exponent E
    // Initial log2(M) based on mantissa
    u.x &= ~(2047LL << 52);
    u.x +=   1023LL << 52;

    val = ((-1/3.) * u.d + 2) * u.d - 2/3.;

    return E + val;
}

// Fix for platforms that call homopolymers as a single thing, with quals
// that decay over time (454, PB CCS).
void homopoly_fix(bam1_t *b) {
    static double ph2err[256] = {0};
    if (!ph2err[0]) {
        for (int i = 0; i < 256; i++)
            ph2err[i] = pow(10, i/-10.0);
    }

    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);
    for (int i = 0; i < b->core.l_qseq; i++) {
        int s = i; // start of homopoly
        int base = bam_seqi(seq, i);
        while (i+1 < b->core.l_qseq && bam_seqi(seq, i+1) == base)
            i++;
        // s..i inclusive is now homopolymer
        if (s == i)
            continue;

        // Simplest:  reverse if end_qual < start_qual
        // Next:      average outer-most two, then next two, etc
        // Best:      fully redistribute so start/end lower qual than centre
        
        // Middle route of averaging outer pairs is sufficient?
        for (int j = s, k = i; j < k; j++,k--) {
            double e = ph2err[qual[j]] + ph2err[qual[k]];
            qual[j] = qual[k] = -fast_log2(e/2)*3.0104+.49;
        }
    }
}

#ifdef DEBUG
int main(int argc, char **argv) {
    qpreserve['Z'] = 1;
    //pblock(argv[1], strlen(argv[1]), atoi(argv[2]), NULL);
    int i, l = strlen(argv[1]);
    for (i = 0; i < l; i++) argv[1][i] -= 33;
    //pblock(argv[1], l, atoi(argv[2]), NULL);
    rblock(argv[1], l, atof(argv[2]), NULL);
    for (i = 0; i < l; i++) argv[1][i] += 33;
    //gen_map(atoi(argv[3]));
    //qbin(argv[1], strlen(argv[1]), qmap);
    puts(argv[1]);
    return 0;
}

#else
#define BLEVEL 10
#define PLEVEL 5
#define RLEVEL 1.5

void usage(void) {
    printf("Usage: qual_smooth [opts] [in_file [out_file]]\n\n");
    printf("Options are:\n");
    printf("  -I FMT[,OPTS]  Specify input format and options\n");
    printf("  -O FMT[,OPTS]  Specify output format and options\n");
    printf("  -k INT[,INT]   Keep quality INT range intact\n");
    printf("  -P INT         P-block smoothing within +/- INT [%d]\n", PLEVEL);
    printf("  -R FLOAT       R-block smoothing within +/- FLOAT [%f]\n", RLEVEL);
    printf("  -B INT         Bin qualities with bin width INT [%d]\n", BLEVEL);
    printf("  -b [A-B=C],... Map qualities A to B inclusive to C; a list\n");
    printf("  -t INT         Use a pool of INT threads for decoding and encoding\n");
    printf("  -H             Redistribute low quals in homopolymers to both ends\n");
    printf("  -X profile     Custom option sets.  Profiles:\n");
    printf("      illumina   -P2 -R1.2 -B4\n");
    printf("      pbccs      -P5 -R1.5 -B10 -H\n");
    printf("      hiseq      -P0 -R0 -b0-0=0,1-1=1,2-9=6,10-19=15,\\\n");
    printf("                     20-24=2225-29=27,30-34=33,35-39=37,40-99=40\n");
    printf("      novaseq    -P0 -R0 -b0-6=2,7-13=12,14-24=18,25-99=36\n");
    printf("  -v             Increase verbosity\n");
    exit(0);
}

#define ENT_BLK 10000000
uint8_t edat[ENT_BLK];
uint8_t Edat[ENT_BLK];

int main(int argc, char **argv) {
    int r;
    bam1_t *b = bam_init1();
    samFile *fp_in = NULL, *fp_out = NULL;
    sam_hdr_t *hdr = NULL;
    htsFormat in_fmt = {0}, out_fmt = {0};
    int blevel = BLEVEL;
    int plevel = PLEVEL;
    int homopoly = 0;
    double rlevel = RLEVEL;
    hts_tpool *tpool = NULL;
    int verbose = 0;

    if (argc == 1 && isatty(0))
	usage();

    // Specify which qualities we must preserve as intact
    qpreserve[0] = 1;
    qpreserve[93] = 1;
    gen_map(blevel); // default

    int c;
    while ((c = getopt(argc, argv, "P:B:O:I:b:t:R:vk:X:H")) >= 0) {
	switch(c) {
	case 'I': hts_parse_format(&in_fmt, optarg); break;
	case 'O': hts_parse_format(&out_fmt, optarg); break;

	case 'k': {
	    char *endp = optarg;
	    do {
		long q1 = strtol(endp, &endp, 10), q2 = q1;
		if (*endp == '-')
		    q2 = strtol(endp+1, &endp, 10);
		do {
		    qpreserve[MAX(0,MIN(255,q1))] = 1 + (c == 'K');
		} while (++q1 <= q2);
	    } while (*endp++ == ',');
	    break;
	}

	case 'P': plevel = MAX(0, atoi(optarg)); break;
	case 'R': rlevel = MAX(1, atof(optarg)); break;
	case 'v': verbose++; break;
	case 'b':
	    if (parse_map(optarg) < 0)
		goto err;
	    blevel = 99;
	    break;

	case 'B':
	    blevel = MAX(1, atoi(optarg));
	    // Set the binning map; blevel==1 is a no-op.
	    gen_map(blevel);
	    break;

	case 't':
	    if (!(tpool = hts_tpool_init(atoi(optarg)))) {
		fprintf(stderr, "Error initialising threads\n");
		goto err;
	    }
	    break;

	case 'X':
	    if (strcmp(optarg, "pbccs") == 0) {
		// -P5 -R1.5 -B10
		plevel = 5;
		rlevel = 1.5;
		gen_map(blevel = 10);
		homopoly = 1;
	    } else if (strcmp(optarg, "illumina") == 0) {
		// -P3 -R1.2 -B8
		plevel = 3;
		rlevel = 1.1;
		gen_map(blevel = 8);
	    } else if (strcmp(optarg, "hiseq") == 0) {
		plevel = 0;
		rlevel = 0;
		parse_map("0-0=0,1-1=1,2-9=6,10-19=15,20-24=22,25-29=27,"
			  "30-34=33,35-39=37,40-99=40");
		blevel = 99;
	    } else if (strcmp(optarg, "novaseq") == 0) {
		plevel = 0;
		rlevel = 0;
		parse_map("0-6=2,7-13=12,14-24=18,25-99=36");
		blevel = 99;
	    } else {
		fprintf(stderr, "Unknown profile \"%s\"\n", optarg);
		exit(1);
	    }
	    break;

	case 'H':
	    homopoly = 1;
	    break;

	case '?':
	default:
	    usage();
	}
    }


    char *fn_in = optind < argc ? argv[optind] : "-";
    if (!(fp_in = sam_open_format(fn_in, "r", &in_fmt)))
	goto err;
    optind++;
    
    char mode[5] = "w";
    char *fn_out = optind < argc ? argv[optind] : "-";
    sam_open_mode(mode+1, fn_out, NULL);
    if (!(fp_out = sam_open_format(fn_out, mode, &out_fmt)))
	goto err;

    if (tpool) {
	htsThreadPool p = {tpool, 0};
	hts_set_thread_pool(fp_in,  &p);
	hts_set_thread_pool(fp_out, &p);
    }

    if (!(hdr = sam_hdr_read(fp_in)))
	goto err;
    if ((r = sam_hdr_write(fp_out, hdr)) < 0)
	goto err;

    size_t longest_read = 0;
    uint8_t *tmp_qual = 0;
    long diff1 = 0, diff2 = 0;
    double e0 = 0, e1 = 0, E0 = 0, E1 = 0;
    int esize = 0;

    while ((r = sam_read1(fp_in, hdr, b)) >= 0) {
        if (verbose && b->core.l_qseq > longest_read) {
	    longest_read = b->core.l_qseq*2;
	    tmp_qual = realloc(tmp_qual, longest_read);
	    if (!tmp_qual) {
		perror("realloc");
		goto err;
	    }
	}
	if (verbose)
	    memcpy(tmp_qual, bam_get_qual(b), b->core.l_qseq);
	if (homopoly)
	    homopoly_fix(b);
        if (plevel > 0)
	    pblock(bam_get_qual(b), b->core.l_qseq, plevel, NULL/*qmap*/);
	if (rlevel > 1)
	    rblock(bam_get_qual(b), b->core.l_qseq, rlevel, NULL/*qmap*/);
	if (blevel > 1)
	    qbin(bam_get_qual(b), b->core.l_qseq, qmap);
	if ((r = sam_write1(fp_out, hdr, b)) < 0)
	    break;

	if (verbose)
	    sim(tmp_qual, bam_get_qual(b), b->core.l_qseq, &diff1, &diff2);
	if (verbose > 1) {
	    if (esize + b->core.l_qseq > ENT_BLK) {
		e0 += entropy0(edat, esize);
		e1 += entropy1(edat, esize);
		E0 += entropy0(Edat, esize);
		E1 += entropy1(Edat, esize);
		esize = 0;
	    }
	    memcpy(edat+esize, tmp_qual,        MIN(ENT_BLK, b->core.l_qseq));
	    memcpy(Edat+esize, bam_get_qual(b), MIN(ENT_BLK, b->core.l_qseq));
	    esize += MIN(ENT_BLK, b->core.l_qseq);
	}
    }

    if (verbose)
	fprintf(stderr, "Diff %ld\nMSE  %ld\n", diff1, diff2);
    if (verbose > 1) {
	e0 += entropy0(edat, esize);
	e1 += entropy1(edat, esize);
	E0 += entropy0(Edat, esize);
	E1 += entropy1(Edat, esize);
	fprintf(stderr, "O(0) %10.0f %10.0f ratio %f\n", e0, E0, E0/e0);
	fprintf(stderr, "O(1) %10.0f %10.0f ratio %f\n", e1, E1, E1/e1);
    }

    if (r != -1) // EOF
	goto err;
    if (sam_close(fp_in) < 0 || sam_close(fp_out) < 0)
	goto err;

    bam_destroy1(b);
    sam_hdr_destroy(hdr);

    if (tpool)
	hts_tpool_destroy(tpool);
    free(tmp_qual);

    hts_opt_free(in_fmt.specific);
    hts_opt_free(out_fmt.specific);

    return 0;

 err:
    fprintf(stderr, "Error\n");
    return 1;
}
#endif
