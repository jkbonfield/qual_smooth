// TODO: also sum average called score, so we can compute real vs expected.

/*
Notes:

Use a self-built consensus (eg samtools consensus, crumble) instead of
a truth set.

Eg GIAB HG002 chr1 claims:

chr1    4376297 .       G       T       50      PASS    ...
chr1    4376333 .       A       G       50      PASS    ...
chr1    4376342 .       T       G       50      PASS    ...

These are all visible in the data set and due to long reads are all
phased together.  However chr1 4376225 also has a SNP which is phased
with the other 3 above and not called.  Consequentially it's labelled
as 3 mismatches.  (It is however missing from the bed file, but it
casts doubt on the completeness of some calls.)

We may be better producing an overly optimistic het consensus call
with confidence, so confident not-het is evaluated, confident is-het
is evaluated, and anything inbetween is skipped.  This means we could
self-train on a denovo data set using the errors in "obvious" regions
to tell us something about the error rates in less obvious regions
too.
*/

/*
 * Given a corrected reference (eg "bcftools consensus" on a VCF truth
 * set), this compares alignments in BAM et al to group sequence and quality
 * into correct and incorrect, permitting kmer based recalibration of
 * quality values.
 *
 * Start with:
 * samtools faidx $HREF38 chr1 \|
 * bcftools consensus -I --mark-ins lc --mark-del "*" truth.vcf > chr1.fa
 * ./qual_train -b truth.bed -f $HREF38 -r chr1 truth.bam
 *
 * Then "*" is del, [acgt] is ins, and RYMSW etc het.
 * => can map ref coords to sample coords.
 *
 * TODO:
 * - Add bed filter regions too.
 *
 * - Reverse complement too
 *   - Either aggregating by a canonical orientation
 *   - (DONE) Or use the REVERSE flag to store in original BAM orientation
 *
 * - Auto-build a hom/het consensus instead of relying on a truth set.
 *   Like Crumble, we could ignore unknown data and concentrate on clear
 *   hom/het/indel positions.
 *
 * - How to handle heterozygous insertions.  Bcftools consensus has no
 *   ambiguity code to represent this?  Samtools consensus uses lowercase,
 *   but then we need an alternative mechanism for cons->ref coord mapping.
 */

#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <ctype.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/thread_pool.h>
#include <htslib/regidx.h>

uint32_t (*kmer_count)[3][2] = NULL;
#define WIN_LEN 5
#define WIN_SHIFT 0 // +ve => move right
#define WIN_MASK ((1<<(3*WIN_LEN))-1)
#define KMER_INIT (05555555555 & WIN_MASK);

// Put kmer mismatch, ins and del into own stats
#define KMAT 0
#define KINS 1
#define KDEL 2

// Accumulate on orig seq strand only.
//
// This may not be as useful as it sounds as gaps during alignment get
// left justified regardless of orientations.  It works OK for mismatches
// though.
#define DO_REVERSE

// Ignore indel cases, so match-pileups only
//#define MATCH_ONLY

// // Discount columns within NEIGHBOUR distance of a previous one
// // TODO
// #define NEIGHBOUR 0

/*
 * Builds a mapping table of reference positions to consensus
 * positions:  map[ref_pos] = cons_pos.
 *
 * Lowercase = insertion
 * * = deletion
 */
hts_pos_t *build_ref_map(const uint8_t *ref, hts_pos_t len) {
    hts_pos_t *map = malloc((len+1) * sizeof(*map));
    if (!map)
	return NULL;

    // As we have deletions shown, can guarantee ref_len <= cons_len.
    hts_pos_t i, r;
    for (i = r = 0; i < len; i++) {
	if (ref[i] == '*' || isupper(ref[i]))
	    map[r++] = i;
    }
    map[r] = i; // so we can compare map[pos] with map[pos+1] to find ins
    return map;
}

// If ref is an ambiguity code and query is compataible, then return
// query, otherwise return N
uint8_t ambig(uint8_t ref, uint8_t query) {
    static uint8_t ambig_base[256] = {0};
    if (!ambig_base[0]) {
	memset(ambig_base, 0xff, 256); // NN
	ambig_base['A'] = ambig_base['a'] = 0x11; // A
	ambig_base['C'] = ambig_base['c'] = 0x22; // C
	ambig_base['G'] = ambig_base['g'] = 0x44; // G
	ambig_base['T'] = ambig_base['t'] = 0x88; // T
	ambig_base['M'] = ambig_base['m'] = 0x12; // AC
	ambig_base['R'] = ambig_base['r'] = 0x14; // AG
	ambig_base['W'] = ambig_base['w'] = 0x18; // AT
	ambig_base['S'] = ambig_base['s'] = 0x24; // CG
	ambig_base['Y'] = ambig_base['y'] = 0x28; // CT
	ambig_base['K'] = ambig_base['k'] = 0x48; // GT
	// any triplets, eg B (not-A) or quads (N) are considered mismatches
    }

    if (ref == query)
	return query;

    if (seq_nt16_str[(ambig_base[ref] & 15)] == query ||
	seq_nt16_str[(ambig_base[ref] >> 4)] == query)
	return query;
    else
	return 'N';
}

/*
 * We iterate over consensus kmers, anchoring them to sequence bases
 * so we can see if they match or mismatch, and then accummulate.
 * We count for seq kmers as that's what we'll be recalibrating later.
 *
 * NB we don't care about differences caused by location of indels.
 * So qry ATTT*G and A*TTTG are the same kmer and both differ
 * to ref ATTTTG.
 *
 * When the seq kmer isn't long enough, at start / end, we compare vs
 * ref kmer with Ns and accumulate in an N-based kmer bin.
 * This also permits us to aggregate many longer kmers into smaller kmers
 * if it gives a better discrimination in quality scores.
 */
static char *context_s(bam1_t *b, int pos) {
    static char ctx[WIN_LEN+1] = {0};
    uint8_t *seq = bam_get_seq(b);
    int len = b->core.l_qseq;

    for (int i = 0, j = pos-WIN_LEN/2+WIN_SHIFT; i < WIN_LEN; i++, j++)
	ctx[i] = j >= 0 && j < len
	    ? seq_nt16_str[bam_seqi(seq, j)]
	    : 'N';

    return ctx;
}

// ACGTN -> {0..4}, 3 bits at a time.
// NOTE: this could be more efficient via a rolling kmer rather than
// recomputing it every time.  This is here for simplicity currently.
static uint32_t context_i(bam1_t *b, int pos) {
    uint8_t *seq = bam_get_seq(b);
    int len = b->core.l_qseq;
    uint32_t ctx = 0;

    // FIXME: switch to two versions for reverse or orig rather than
    // one vers and a second loop to reverse.

#ifdef DO_REVERSE
    int shift = (b->core.flag & BAM_FREVERSE) ? -WIN_SHIFT : WIN_SHIFT;
#else
    int shift = WIN_SHIFT;
#endif
    for (int i = 0, j = pos-WIN_LEN/2+shift; i < WIN_LEN; i++, j++)
	ctx = (ctx<<3) | 
	    (j >= 0 && j < len
	     ? seq_nt16_int[bam_seqi(seq, j)]
	     : 4);

#ifdef DO_REVERSE
    if (b->core.flag & BAM_FREVERSE) {
	uint32_t ctx2 = 0;
	//fprintf(stderr, "%05o ", ctx);
	for (int i = 0; i < WIN_LEN; i++) {
	    ctx2 = (ctx2<<3) | ((ctx & 7)^3);
	    ctx >>= 3;
	}
	//fprintf(stderr, "%05o\n", ctx2);
	ctx = ctx2;
    }
#endif

    return ctx;
}

int in_bed(regitr_t *bed_itr, hts_pos_t pos) {
    if (!bed_itr)
	return 1;

    while (pos > bed_itr->end) {
	if (!regitr_loop(bed_itr)) {
	    bed_itr->beg = HTS_POS_MAX; // next region is never
	    return 0;
	}
    }
    return pos >= bed_itr->beg && pos <= bed_itr->end ? 1 : 0;
}

void accumulate_kmers(sam_hdr_t *hdr, const uint8_t *ref, hts_pos_t *map,
		      bam1_t *b, regidx_t *bed, regitr_t *bed_itr) { 
    const int V=0; // DEBUG only

    if (bed_itr)
	regidx_overlap(bed, sam_hdr_tid2name(hdr, b->core.tid),
		       b->core.pos, bam_endpos(b), bed_itr);

    static uint8_t L[256], L_done = 0;
    if (!L_done) {
	memset(L, 4, 256);
	L['A'] = L['a'] = 0;
	L['C'] = L['c'] = 1;
	L['G'] = L['g'] = 2;
	L['T'] = L['t'] = 3;
	L_done = 1;

	kmer_count = calloc(WIN_MASK+1, sizeof(*kmer_count));
    }

    hts_pos_t rpos = b->core.pos; // ref pos; cons pos = map[rpos]
    int qpos = 0;                 // query pos
    uint32_t *cig = bam_get_cigar(b);
    uint8_t *qseq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);
    int ncig = b->core.n_cigar;

    if (V) printf("SEQ %s\n", bam_get_qname(b));
    for (int i = 0; i < ncig; i++) {
	int cig_op = bam_cigar_op(cig[i]);
	int cig_len = bam_cigar_oplen(cig[i]);
	int nth = 0;

	if (V>1) printf("> %d%c\n", cig_len, BAM_CIGAR_STR[cig_op]);

	if (bam_cigar_type(cig_op) == 0)
	    // Consumes neither seq or ref
	    continue;

	if (cig_op == BAM_CINS)
	    rpos--, nth = 1; // ins after last last, not before next

	for (int j = 0; j < cig_len; j++) {
	    uint8_t qbase = seq_nt16_str[bam_seqi(qseq, qpos)];
	    uint8_t rbase = ref[map[rpos]+nth];

	    uint32_t kmer = context_i(b, qpos);
	    if (V>1) printf("%.5s\t%05o\t", context_s(b, qpos), kmer);

	    switch (cig_op) {
	    case BAM_CSOFT_CLIP:
		// Seq only, just ignore
		if (V>1) printf("S %ld\t. %c\n", rpos, qbase);
		qpos++;
		break;

	    case BAM_CMATCH:
	    case BAM_CDIFF:
	    case BAM_CEQUAL:
		// Both seq and ref continue in sync.
		if (islower(rbase)) {
		    if (V) printf("dm %ld\t%c %c *\n", rpos, rbase, qbase);
#ifndef MATCH_ONLY
		    if (in_bed(bed_itr, rpos)) // base vs cons-ins
			kmer_count[kmer][KDEL][0]++;
#endif
		} else {
		    if (qbase == ambig(rbase, qbase)) {
			if (V>1)
			    printf("M  %ld\t%c %c\n", rpos, rbase, qbase);
			if (in_bed(bed_itr, rpos))
			    kmer_count[kmer][KMAT][1]++;
		    } else if (rbase == '*') {
			if (V)
			    printf("im %ld\t%c %c *\n", rpos, rbase, qbase);
#ifndef MATCH_ONLY
			if (in_bed(bed_itr, rpos)) // ins vs cons-del
			    kmer_count[kmer][KINS][0]++;
#endif
		    } else if (rbase != 'N') {
			if (V) {
			    if (V<2)
				printf("%.5s\t%05o\t", context_s(b, qpos), kmer);
			    printf("X  %ld\t%c %c %d *\n", rpos, rbase, qbase, qual[qpos]);
			}
			if (in_bed(bed_itr, rpos)) // mismatch
			    kmer_count[kmer][KMAT][0]++;
		    }
		}

		// Check if map[rpos] and map[rpos+1] are more than 1 apart,
		// indicating insertion in cons which is absent in seq.
		if (map[rpos+1] - map[rpos] != 1) {
		    if (V)
			printf("mD %ld\t. %c *\n", rpos,  qbase);
#ifndef MATCH_ONLY
		    if (in_bed(bed_itr, rpos)) // ins in cons, absent here
			kmer_count[kmer][KDEL][0]++;
#endif
		}

		qpos++;
		rpos++;
		break;

	    case BAM_CINS:
		if (islower(rbase)) {
		    // TODO: also check base matches
		    if (qbase == ambig(rbase, qbase)) {
			if (V>1) printf("mi %ld\t%c %c\n", rpos, rbase, qbase);
#ifndef MATCH_ONLY
			// Actually [KMAT][1]
			if (in_bed(bed_itr, rpos)) // ins matching cons-ins
			    kmer_count[kmer][KINS][1]++;
#endif
		    } else {
			if (V>1) printf("xi %ld\t%c %c\n", rpos, rbase, qbase);
#ifndef MATCH_ONLY
			// Maybe [KMAT][0]?
			if (in_bed(bed_itr, rpos)) // mismatching cons-ins
			    kmer_count[kmer][KINS][0]++;
#endif
		    }
		} else {
		    if (V) printf("I  %ld\t. %c *\n", rpos, qbase);
#ifndef MATCH_ONLY
		    if (in_bed(bed_itr, rpos))
			kmer_count[kmer][KINS][0]++; // ins vs cons
#endif
		}
		qpos++;
		nth++;
		break;

	    case BAM_CDEL:
		if (rbase == '*') {
		    if (V>1) printf("md %ld\t%c .\n", rpos, rbase);
#ifndef MATCH_ONLY
		    if (in_bed(bed_itr, rpos))
			kmer_count[kmer][KDEL][1]++; // del matching cons-del
#endif
		} else {
		    if (V) printf("D  %ld\t%c . *\n", rpos, rbase);
#ifndef MATCH_ONLY
		    if (in_bed(bed_itr, rpos))
			kmer_count[kmer][KDEL][0]++; // del vs cons-base
#endif
		}
		rpos++;
		break;
	    }
	}

	if (cig_op == BAM_CINS)
	    rpos++; // undo pos-- above
    }
}

void dump_kmers(void) {
    int i, j, k;
    puts("=== kmers ===\n");
    for (i = 0; i <= WIN_MASK; i++) {
	// We record INS here vs INS cons as a valid INS,
	// but logically this is really a "match" (or mismatch) state.
	// Hence we use the total count across all types rather than just
	// the specific one for this type.  Hence we could possibly
	// combine all "pass" states into a single count during analysis
	// and only have "fail" states for the separate types, but this is
	// left to future exploration.
	int cnt =
	    kmer_count[i][0][0]+kmer_count[i][0][1] +
	    kmer_count[i][1][0]+kmer_count[i][1][1] +
	    kmer_count[i][2][0]+kmer_count[i][2][1];
	if (!cnt)
	    continue;

	int all_mis = 0;
	for (k = 0; k < 3; k++) {
	    for (j = WIN_LEN-1; j >= 0; j--)
		putchar("ACGTNNNN"[(i>>(j*3))&7]);
	    double err = cnt ? (double)kmer_count[i][k][0] / cnt : 0;
	    int qval = err ? -10*log10(err)+.5 : 99;
	    printf("\t%c\t%d\t%d\t%d\n",
		   "MID"[k], kmer_count[i][k][0], cnt, qval);

	    all_mis += kmer_count[i][k][0];
	}

	// Combined stats
	for (j = WIN_LEN-1; j >= 0; j--)
	    putchar("ACGTNNNN"[(i>>(j*3))&7]);
	double err = all_mis ? (double)all_mis / cnt : 0;
	int qval = err ? -10*log10(err)+.5 : 99;
	printf("\t?\t%d\t%d\t%d\n", all_mis, cnt, qval);
    }
}

int main(int argc, char **argv) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    htsFormat in_fmt = {0};
    faidx_t *fai = NULL;
    hts_pos_t *map = NULL;
    bam1_t *b = bam_init1();
    uint8_t *ref = NULL;
    regidx_t *bed = NULL;
    regitr_t *bed_itr = NULL;
    char *reg = NULL;
    hts_itr_t *itr = NULL;

    int c;
    while ((c = getopt(argc, argv, "I:f:b:r:")) >= 0) {
	switch(c) {
	case 'I': hts_parse_format(&in_fmt, optarg); break;

	case 'f': {
	    if (!(fai = fai_load(optarg))) {
		fprintf(stderr, "Failed to load consensus fasta\n");
		goto err;
	    }
	    break;
	}
	    
	case 'b':
	    bed = regidx_init(optarg, NULL, NULL, 0, NULL);
	    if (!bed) {
		fprintf(stderr, "Failed to load bed file\n");
		goto err;
	    }
	    break;

	case 'r':
	    reg = optarg;
	    break;

	case '?':
	default:
	    fprintf(stderr, "Usage: [blah]\n");
	    exit(0);
	}
    }

    if (!fai) {
	fprintf(stderr, "Use -f cons.fa to specify the consensus file\n");
	return 1;
    }

    // Open files
    char *fn = optind < argc ? argv[optind] : "-";
    if (!(fp = sam_open_format(fn, "r", &in_fmt))) {
	perror(fn);
	goto err;
    }
    if (!(hdr = sam_hdr_read(fp)))
	goto err;

    bed_itr = bed ? regitr_init(bed) : NULL;

    if (reg) {
	hts_idx_t *idx = sam_index_load(fp, fn);
	if (!idx) {
	    fprintf(stderr, "Unable to load index\n");
	    goto err;
	}
	
	if (!(itr = sam_itr_querys(idx, hdr, reg))) {
	    fprintf(stderr, "Unable to query region '%s'\n", reg);
	    goto err;
	}

	hts_idx_destroy(idx);
    }

    // Process alignments
    int r;
    int last_tid = -1;
    hts_pos_t ref_len;
    while ((r = (reg
		 ? sam_itr_next(fp, itr, b)
		 : sam_read1(fp, hdr, b))) >= 0) {
	if (b->core.tid < 0)
	    continue;
	if (b->core.tid != last_tid) {
	    if (ref)
		free(ref);
	    ref = (uint8_t *)
		faidx_fetch_seq64(fai, sam_hdr_tid2name(hdr, b->core.tid),
				  0, HTS_POS_MAX, &ref_len);
	    if (!ref)
		goto err;
	    if (map)
		free(map);
	    if (!(map = build_ref_map(ref, ref_len)))
		goto err;
	    last_tid = b->core.tid;
	}
	//puts(bam_get_qname(b));

	accumulate_kmers(hdr, ref, map, b, bed, bed_itr);
    }
    if (r != -1)
	goto err;

    if (sam_close(fp) < 0)
	goto err;
    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    fai_destroy(fai);
    if (ref) free(ref);
    if (map) free(map);
    if (bed) regidx_destroy(bed);
    if (bed_itr) regitr_destroy(bed_itr);
    if (itr) hts_itr_destroy(itr);

    dump_kmers();

    return 0;

 err:
    bam_destroy1(b);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp) sam_close(fp);
    if (fai) fai_destroy(fai);
    if (ref) free(ref);
    if (map) free(map);
    return 1;
}
