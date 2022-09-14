/*
 * Given a corrected reference (eg "bcftools consensus" on a VCF truth
 * set), this compares alignments in BAM et al to group sequence and quality
 * into correct and incorrect, permitting kmer based recalibration of
 * quality values.
 *
 * Start with:
 * samtools faidx $HREF38 chr1 \|
 * bcftools consensus -I --mark-ins lc --mark-del "*" truth.vcf > chr1.fa
 *
 * Then "*" is del, [acgt] is ins, and RYMSW etc het.
 * => can map ref coords to sample coords.
 *
 * TODO:
 * - Add bed filter regions too.
 * - Consider both strands
 */

#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <ctype.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/thread_pool.h>

uint32_t (*kmer_count)[2] = NULL;
#define WIN_LEN 5 // FIXME: needs to be centred
#define WIN_MASK ((1<<(3*WIN_LEN))-1)
#define KMER_INIT (05555555555 & WIN_MASK);

/*
 * Builds a mapping table of reference positions to consensus
 * positions:  map[ref_pos] = cons_pos.
 *
 * Lowercase = insertion
 * * = deletion
 */
hts_pos_t *build_ref_map(const uint8_t *ref, hts_pos_t len) {
    hts_pos_t *map = malloc(len * sizeof(*map));
    if (!map)
	return NULL;

    // As we have deletions shown, can guarantee ref_len <= cons_len.
    hts_pos_t i, r;
    for (i = r = 0; i < len; i++) {
	if (ref[i] == '*' || isupper(ref[i]))
	    map[r++] = i;
    }
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
void accumulate_kmers(const uint8_t *ref, hts_pos_t *map, bam1_t *b) {
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

    hts_pos_t pos = b->core.pos;
    uint32_t *cig = bam_get_cigar(b);
    uint8_t *qseq = bam_get_seq(b);
    int ncig = b->core.n_cigar, i = 0;
    int cig_op = 0, cig_len = 0, cig_idx = 0;
    int nth = 0; // nth base in an insertion
    int diff_type = 0;

    //printf("New seq at %ld\n", pos);
    uint32_t qkmer = 0, rkmer = 0, diff_in = 0;
    int kmer_skip = WIN_LEN;
    for (i = 0; i < b->core.l_qseq; i++) {
	uint8_t qbase = seq_nt16_str[bam_seqi(qseq, i)];
	qkmer = qkmer*8 + seq_nt16_int[bam_seqi(qseq, i)];
	qkmer = qkmer & WIN_MASK;
	kmer_skip -= kmer_skip>0;
	int valid = 1;

    next_op:
	if (cig_len == 0) {
	    if (cig_idx >= ncig) {
		fprintf(stderr, "Malformed CIGAR string\n");
		return;
	    }
	    cig_len = bam_cigar_oplen(cig[cig_idx]);
	    cig_op  = bam_cigar_op(cig[cig_idx]);
	    cig_idx++;
	    printf("CIG: %d %c\n", cig_len, BAM_CIGAR_STR[cig_op]);
	    if (cig_op == BAM_CINS)
		pos--, nth = 1; // ins after last base, not before this
	    else
		nth = 0;
	}

	printf("%d %c%c%c\t", i, qbase, ref[map[pos]+nth],
	    ambig(ref[map[pos]+nth], qbase));
	if (bam_cigar_type(cig_op) == 2) {
	    // Consumes ref only
	    do {
		uint8_t rbase = ambig(ref[map[pos++]], qbase);
		if (rbase == '*') {
		    // deletion in consensus too, so ignore
		    break; // ?
		} else {
		    rkmer = rkmer*8 + L[rbase];
		    rkmer = rkmer & WIN_MASK;
		}
		printf("(D)\t%05o\t%05o\t%d\n", qkmer, rkmer, kmer_skip);
		// Even if [qr]kner matches, it's something to log
		valid = 0;
	    } while (--cig_len > 0);
	    cig_len++;
	    //goto next_op;
	} else if (bam_cigar_type(cig_op) & 2) {
	    // Consumes both ref and query
	    uint8_t rbase = ambig(ref[map[pos++]], qbase);
	    rkmer = rkmer*8 + L[rbase];
	    rkmer = rkmer & WIN_MASK;
	    printf("M\t%05o\t%05o\t%d\n", qkmer, rkmer, kmer_skip);
	} else if (bam_cigar_type(cig_op) & 1) {
	    // Consumes seq only
	    // If soft-clip, not an cons diff
	    // If insertion, it's a difference unless cons has lowercase
	    if (cig_op == BAM_CSOFT_CLIP) {
		printf("S\t%05o\t%05o\t%d\n", qkmer, rkmer, kmer_skip);
		i+=cig_len-1;
		cig_len -= cig_len-1;
		kmer_skip = WIN_LEN;
	    } else {
		if (islower(ref[map[pos]+nth])) {
		    uint8_t rbase = ref[map[pos]+nth++];
		    rkmer = rkmer*8 + L[rbase];
		    rkmer = rkmer & WIN_MASK;
		    printf("i\t%05o\t%05o\t%d\n", qkmer, rkmer, kmer_skip);
		} else {
		    printf("I\t%05o\t%05o\t%d\n", qkmer, rkmer, kmer_skip);
		    valid = 0;
		    //kmer_skip++;
		}
		if (cig_len == 1) pos++; // correct for pos-- above
	    }
	} else {
	    // Consumes neither - not aligned so just skip
	    // TODO: may want to clear qkmer to have Ns here?
	}

	// If we have a single bp difference, do we want to only
	// count it as an error when the diff is the middle base (and
	// not count as valid for the others). Ie
	//
	// AACCG AACCG Y
	// ACCGT ACCGT Y
	// CCGTA CCGTT ignore
	// CGTAA CGTTA ignore
	// GTAAC GTTAC N
	// TAACC TTACC ignore
	// AACCG TACCG ignore
	// ACCGA ACCGA Y
	// CCGAG CCGAG Y
	if (!kmer_skip) {
	    //kmer_count[qkmer][qkmer==rkmer]++;
	    if (qkmer == rkmer && valid) {
		kmer_count[qkmer][1]++;
	    } else {
		if (!diff_type) diff_type = BAM_CIGAR_STR[cig_op];
		if (diff_in > 1) {
		    diff_in--;
		    printf("Delay1\n");
		} else if (diff_in == 1) {
//		if (((qkmer >> 3*(WIN_LEN/2)) & 7) !=
//		     ((rkmer >> 3*(WIN_LEN/2)) & 7)) {
		    printf("Diff %c %05o\n", diff_type, qkmer);
		    kmer_count[qkmer][0]++;
		    kmer_skip += WIN_LEN/2+1;
		    diff_in = 0;
		    diff_type = 0;
		} else {
		    printf("Delay0\n");
		    diff_in = WIN_LEN/2;
		}
	    }
	} else {
	    printf("Skip %d\n", kmer_skip);
	}
//	if (qkmer != rkmer)
//	    printf("Diff\n");
	cig_len--;

/*
	switch (cig_op) {
	case BAM_CMATCH:
	case BAM_CEQUAL:
	case BAM_CBACK:
	    
	}
*/

/*
	if (bam_cigar_type(cig_op) & 1) {
	    // consumes query
	}
	if (bam_cigar_type(cig_op) & 2) {
	    // consumes reference
	}
	// if neither loop for next cigar?
*/
    }
}

void dump_kmers(void) {
    int i, j;
    puts("=== kmers ===\n");
    for (i = 0; i <= WIN_MASK; i++) {
//	int valid = 1;
//	for (j = 0; j < WIN_LEN; j++)
//	    if (((i>>(j*3))&7) > 4)
//		valid = 0;
//	if (!valid)
//	    continue;

	int cnt = kmer_count[i][0]+kmer_count[i][1];
	if (!cnt)
	    continue;

	for (j = WIN_LEN-1; j >= 0; j--)
	    putchar("ACGTN---"[(i>>(j*3))&7]);
	double err = cnt ? (double)kmer_count[i][0] / cnt : 0;
	int qval = err ? -10*log10(err)+.5 : 99;
	printf("\t%d\t%d\t%d\t%d\n", kmer_count[i][0]+kmer_count[i][1],
	       kmer_count[i][0], kmer_count[i][1], qval);
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

    int c;
    while ((c = getopt(argc, argv, "I:f:")) >= 0) {
	switch(c) {
	case 'I': hts_parse_format(&in_fmt, optarg); break;
	case 'f': {
	    if (!(fai = fai_load(optarg))) {
		fprintf(stderr, "Failed to load consensus fasta");
		goto err;
	    }
	    break;
	}
	    
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

    // Process alignments
    int r;
    int last_tid = -1;
    hts_pos_t ref_len;
    while ((r = sam_read1(fp, hdr, b)) >= 0) {
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

	accumulate_kmers(ref, map, b);
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
