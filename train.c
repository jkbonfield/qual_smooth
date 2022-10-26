/*
 * Tracking of kmers for correctly identified in consensus vs incorrect,
 * based on the central kmer base.  E.g.
 *
 * CONS:  AGCTTGAGCTAGG
 *              x
 * SEQ:   AGCTTGTGCTAGG
 *
 * 3 types of state per kmer
 * - KMER is correct, cented    (AGCTT GCTTG GCTAG CTAGG)   TRUE
 * - KMER is incorrect, centred (TGTGC)                     FALSE
 * - KMER overlaps error        (CTTGT TTGTG GTGCT TGCTA)   IGNORE
 *
 * To do this, we could compute consensus kmers and seq kmers and
 * compare, but consensus has two alleles and they're unphased.
 * So instead we march through seq and ref one base at a time, comparing
 * that central position, and mark our kmers as we go.  These need to
 * have a history before fully updating counts, as we may need to edit
 * preious <correct> kmers and revise them to <ignore> when we see a
 * <incorrect> status (plus subsequent few few too, if not erroneous).
 *
 * Q: for deletions, should we also include the consensus kmers?
 * Eg if CONS is AGTTTAGG and we have seq AG(3D)AGG then we're only
 * recording one deletion at the D cigar op, but maybe all the kmers with
 * T in are erroneous, so when we report them as seen elsewhere we know
 * their probability of being missed.
 * Currently we don't.  Our kmer is purely from seq.  Se given we've seen
 * it in the seq, what is the probability of X, I and D errors.
 */

/*
 * Types are classified in the context of the consensus, not ref.
 * So consider a homozygous insertion to the reference.
 *
 * If our seq has the insertion, it is a match (recorded as correct ins).
 * If it doesn't it is a mismatch (even if it matches ref) and is recorded
 * as incorrect ins.
 *
 * Similarly for deletions and matches.
 *
 * While arguably seq matching cons perfect is "match", we separate by
 * consensus type as this allows us to judge the impact of reference bias.
 * In theory match_false/match_true should be the same ratio as
 * ins_false/ins_true, but where not this is because we aligned preferentially
 * to reference and not the real sample consensus.  (A realigner should
 * solve that, so we have a way to evaluate those too.)
 *
 * Hence KMAT, KINS and KDEL are all "matching consensus" categories and
 * allow us to detect the substitution errors from the instrument.
 * The "true" count for KMAT/KINS/KDEL are true cons matches.
 * The "false" count for KMAT is substitution errors.
 *
 * We also have overcalls and undercalls categories, where we have additional
 * insertions and deletions in the sequence with respect to the consensus.
 * These are recorded in the "false" count for KINS/KDEL respectively.
 *
 * FIXME: what about correct INS base counts but subtitution errors?
 * For now we're not recording those.
 */

/*
Input: "bcftools consensus" format:
- 1 upper-case base per 1 ref base (so coords map)
  - Ambiguity codes (M, R, Y etc) for known SNPs
- lower-case cons base = inserted base
- * = deletion

Problem: cannot represent ambiguity in insertion, so reads lacking the
insertion are labelled as having a deletion vs the sample.

Input: "samtools consensus --mark-ins" format:
- Standard uppercase bases are in ref coordinates
- * is both-allele deletion
- Lowercase is single-allele deletion
- Underscore means next base is insertion (overules 1st statement)

Thus insertions can themselves be uppercase (both alleles with insertion,
maybe ambiguity code if they differ) or lowercase (one allele only with
insertion).

Quality values, if in FASTQ format, also have the "_" symbol.


Both input formats internally get converted to a third different format:
- Uppercase 7-bit ASCII = aligned to a ref coord
- * or lowercase means a deletion.
- Top bit set means insertion (x80+upper or x80+lower for hom/het).

This gives us easier matching of nth base of insertion vs nth pos in ref.
*/


/* TODO:
- Additional filtering via halo around lower-case and within N bases
  of STRs, etc.  This ensures we only calibrate on regions with
  certainty in the consensus.
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
#include <float.h>
#include <ctype.h>
#include <assert.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/thread_pool.h>
#include <htslib/regidx.h>

// -- good ones
#define K_MAT_M 0 // match in M op
#define K_MAT_I 1 // match in I op
#define K_MAT_D 2 // match in D op

// -- bad ones
#define K_WRONG 3 // wrong below this point
#define K_MIS_M 3 // mismatch in M op
#define K_MIS_I 4 // mismatch in I op
#define K_OVER  5 // overcall;  rogue I op
#define K_UNDER 6 // undercall; rogue D op
#define K_NCAT  (K_UNDER+1)
#define K_CAT "MIDxiou"
uint32_t (*k_count)[K_NCAT] = NULL;
double   (*k_qual) [K_NCAT] = NULL;
double   (*k_qual2)[K_NCAT][99] = NULL;

#define ST_HALO 50 // distance from indel
#define ST_NEAR_INS 1
#define ST_NEAR_DEL 2
#define ST_NEAR_STR 4



// [4<<WIN_LEN][TYPE][IS_CORRECT]  (type being match, ins and del)
uint32_t (*kmer_count)[3][2] = NULL; // counts of true and false bases
double   (*kmer_qual)[3] = NULL;     // sum of estimated errors from qual
int      (*kmer_qual2)[3][99] = NULL;     // sum of estimated errors from qual
#define WIN_LEN 5
#define WIN_SHIFT 0 // +ve => move right
#define WIN_MASK ((1<<(3*WIN_LEN))-1)
#define KMER_INIT (05555555555 & WIN_MASK);

#define MIN(a,b) ((a)<(b)?(a):(b))


double Perr[256];

// Put kmer mismatch, ins and del into own stats
#define KMAT   0 // matching base
#define KINS   1 // matching insertion
#define KDEL   2 // matching deletion

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
 * If mark_ins is zero then we assume this is a bcftools consensus file:
 *    * = deletion
 *    Lowercase = insertion
 *
 * If mark_ins is non-zero then we assume it to be a samtools consensus
 * file produced with -a --show-del yes --mark-ins.
 *    * = homozygous deletion
 *    _ = next base is insertion (_ as + is probelmatic at start of a line)
 *    Lowercase = heterozygous indel
 *
 * Internally we change both consensus sequences to the same third form.
 *    * = homozygous deletion
 *    lowercase = heterozygous deletion (1 allele deleted)
 *    Top bit set = insertion.
 *    Top-bit "lowercase" = heterozygous insertion (1 allele inserted)
 */
hts_pos_t *build_ref_map(uint8_t *ref, hts_pos_t *len_p, uint8_t *stat,
			 int mark_ins) {
    hts_pos_t len = *len_p;
    hts_pos_t *map = malloc((len+1) * sizeof(*map));
    if (!map)
	return NULL;

    // As we have deletions shown, can guarantee ref_len <= cons_len.
    hts_pos_t i, r;
    if (mark_ins == 0) {
	for (i = r = 0; i < len; i++) {
	    if (ref[i] == '*' || isupper(ref[i]))
		map[r++] = i;
	    else if (islower(ref[i]))
		ref[i] = toupper(ref[i]) | 0x80; // internal marker for ins

	    if (ref[i] == '*') {
		for (int z = -ST_HALO; z <= ST_HALO; z++)
		    if (i+z >= 0 && i+z < len)
			stat[i+z] |= ST_NEAR_DEL;
	    }
	    if (ref[i] & 0x80) {
		for (int z = -ST_HALO; z <= ST_HALO; z++)
		    if (i+z >= 0 && i+z < len)
			stat[i+z] |= ST_NEAR_INS;
	    }
	}
    } else {
	// TODO use qual and have a filter for ref bases to skip.
	int k;
	for (i = k = r = 0; i < len; i++) {
	    if (ref[i] == '_') {
		ref[k++] = ref[++i] | 0x80; // mark insertion
		for (int z = -ST_HALO; z <= ST_HALO; z++)
		    if (k+z >= 0 && k+z < len)
			stat[k+z] |= ST_NEAR_INS;
	    } else {
		map[r++] = k;
		ref[k++] = ref[i];

		if (islower(ref[i])) {
		for (int z = -ST_HALO; z <= ST_HALO; z++)
		    if (k+z >= 0 && k+z < len)
			stat[k+z] |= ST_NEAR_DEL;
		}
	    }
	}
	*len_p = k;
    }

    //for (i = 0; i < len; i++) printf("%ld\t%ld\t%d\n", i, map[i], ref[i]);

    map[r] = i; // so we can compare map[pos] with map[pos+1] to find ins
    return map;
}

// If ref is an ambiguity code and query is compataible, then return
// query, otherwise return N
uint8_t ambig(uint8_t ref, uint8_t query) {
    ref &= ~0x80;
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

int kmer_skip = 0;
int kmer_hist_kmer[WIN_LEN/2];
int kmer_hist_type[WIN_LEN/2];
int kmer_hist_ok[WIN_LEN/2];
double kmer_hist_qual[WIN_LEN/2];
int kmer_hist_qual2[WIN_LEN/2];
int64_t kmer_num = 0;

static inline void incr_kmer(regitr_t *bed_itr, hts_pos_t rpos,
			     uint32_t kmer, int type, int ok, int qual) {
//    printf("Incr %05o %c %d %d\n", kmer, "MID"[type], ok, qual);

    if (in_bed(bed_itr, rpos)) {
	// Keep list of last WIN_LEN/2 kmers added (if OK), so we can
	// undo them if we add a not-OK one.  Similarly track the numeber
	// of subsequent kmers to skip (if OK).
	if (kmer_skip && ok) {
	    kmer_skip--;
	    //kmer_hist_ok[(kmer_num++) % (WIN_LEN/2)] = 0;
	    //printf("Skip %05o\n", kmer);
	} else {
	    kmer_count[kmer][type][ok]++;
	    kmer_qual[kmer][type] += Perr[qual];
	    kmer_qual2[kmer][type][qual]++;
	}

	// Revert to ignore
	if (!ok) {
	    for (int i = 0; i < WIN_LEN/2; i++) {
		int h_ok = kmer_hist_ok[i];
		if (!h_ok) {
		    //printf("Kmer %05o wasn't OK anyway\n", kmer_hist_kmer[i]);
		    continue;
		}		    
		kmer_count[kmer_hist_kmer[i]][kmer_hist_type[i]][h_ok]--;
		kmer_qual [kmer_hist_kmer[i]][kmer_hist_type[i]]
		    -= kmer_hist_qual[i];
		kmer_qual2[kmer_hist_kmer[i]][kmer_hist_type[i]][kmer_hist_qual2[i]]--;
		kmer_hist_ok[i] = 0; // prevent another decr on next-err
		//printf("Reset %05o to ignore\n", kmer_hist_kmer[i]);
	    }
	    kmer_skip = WIN_LEN/2;
	}

	// cache
	int idx = kmer_num % (WIN_LEN/2); // could also round and AND
	kmer_hist_kmer[idx] = kmer;
	kmer_hist_type[idx] = type;
	kmer_hist_ok  [idx] = ok;
	kmer_hist_qual[idx] = Perr[qual];
	kmer_hist_qual2[idx] = qual;
	kmer_num++;
    }
}

static inline void incr_kmer2(regitr_t *bed_itr, uint8_t stat, hts_pos_t rpos,
			      uint32_t kmer, int type, int ok, int qual) {
    int ok2 = type<K_WRONG;
    assert(ok == ok2);
    
//    if (stat) // not near an indel
//	return;

//    printf("Incr %05o %c %d %d %d\n", kmer, K_CAT[type], ok, qual, type);
    if (in_bed(bed_itr, rpos)) {
	// Keep list of last WIN_LEN/2 kmers added (if OK), so we can
	// undo them if we add a not-OK one.  Similarly track the numeber
	// of subsequent kmers to skip (if OK).
	int skipped;
	if (kmer_skip && ok) {
	    kmer_skip--;
	    skipped=1;
//	    printf("Skip %05o\n", kmer);
	} else {
	    k_count[kmer][type]++;
	    k_qual[kmer][type] += Perr[qual];
	    k_qual2[kmer][type][qual]++;
//	    if (kmer == 031022 && type == K_OVER)
//		printf("31022: overcall with qual %d\n", qual);
	    skipped=0;
	}

	// Revert to ignore
	if (!ok) {
	    for (int i = 0; i < WIN_LEN/2; i++) {
		int h_ok = kmer_hist_type[i] < K_WRONG;
//		printf("%s %05o to ignore, count %d type %d\n",
//		       h_ok ? "Reset" : "RDone",
//		       kmer_hist_kmer[i],
//		       k_count[kmer_hist_kmer[i]][kmer_hist_type[i]],
//		       kmer_hist_type[i]);
		if (h_ok) {
		    k_count[kmer_hist_kmer[i]][kmer_hist_type[i]]--;
		    k_qual [kmer_hist_kmer[i]][kmer_hist_type[i]]
			-= kmer_hist_qual[i];
		    k_qual2[kmer_hist_kmer[i]][kmer_hist_type[i]][kmer_hist_qual2[i]]--;
		    if (k_count[kmer_hist_kmer[i]][kmer_hist_type[i]] > 2000000000) {
			fflush(stdout);
			fflush(stderr);
			exit(1);
		    }
		}
		kmer_hist_type[i] = 99; // prevent double decr is more err
	    }
	    kmer_skip = WIN_LEN/2;
	}

	// cache
	int idx = kmer_num % (WIN_LEN/2); // could also round and AND
	kmer_hist_kmer[idx] = kmer;
	kmer_hist_type[idx] = skipped ? 99 : type;
	//printf("type[%05o,%d] = %d / %d\n", kmer, idx, type, kmer_hist_type[idx]);
	kmer_hist_qual[idx] = Perr[qual];
	kmer_hist_qual2[idx] = qual;
	kmer_num++;
    } else {
	// TODO: set kmer_hist_type to 99.  Optimise this as only
	// need it on in-bed to out-bed transitions
	for (int i = 0; i < WIN_LEN/2; i++)
	    kmer_hist_type[i] = 99;
    }
}

void accumulate_kmers(sam_hdr_t *hdr, const uint8_t *ref,
		      const uint8_t *stat, hts_pos_t *map,
		      bam1_t *b, regidx_t *bed, regitr_t *bed_itr,
		      int mark_ins) { 
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
	kmer_qual  = malloc((WIN_MASK+1) * sizeof(*kmer_qual));
	kmer_qual2 = calloc((WIN_MASK+1), sizeof(*kmer_qual2));
	int i, j, k;
	for (k = 0; k < WIN_MASK; k++)
	    for (j = 0; j < 3; j++)
		kmer_qual[k][j] = DBL_MIN;

	k_count = calloc(WIN_MASK+1, sizeof(*k_count));
	k_qual  = calloc(WIN_MASK+1, sizeof(*k_qual));
	k_qual2 = calloc(WIN_MASK+1, sizeof(*k_qual2));
	for (k = 0; k < WIN_MASK; k++)
	    for (j = 0; j < K_NCAT; j++)
		kmer_qual[k][j] = DBL_MIN;

	Perr[0] = 1;
	for (i = 1; i < 256; i++)
	    Perr[i] += pow(10, -i/10.0);
    }

    hts_pos_t rpos = b->core.pos; // ref pos; cons pos = map[rpos]
    int qpos = 0;                 // query pos
    uint32_t *cig = bam_get_cigar(b);
    uint8_t *qseq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);
    int ncig = b->core.n_cigar;
    int ndel = 0;

// 0x80 = ins                                   80-F0 1xxx ....
// 0x80+uppercase = HOM ins (both alleles)      C0-D0 110x ....
// 0x80+lowercase = HET ins (one allele only)   E0-F0 111x ....
#define isins(c) ((c) & 0x80)
#define isins_hom(c) (((c) & 0xe0) == 0xc0)
#define isins_het(c) (((c) & 0xe0) == 0xe0)
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
	    uint8_t qqual = qual[qpos];
	    uint8_t rbase = ref[map[rpos]+nth];
	    uint8_t rst   = stat[map[rpos]+nth];

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
		if (isins_hom(rbase)) {
		    if (V) printf("dm %ld\t%c %c *\n", rpos, rbase, qbase);
#ifndef MATCH_ONLY
		    // TYPE: undercall (cons has an insertion, we did not)
		    //incr_kmer(bed_itr, rpos, kmer, KDEL, 0, qqual);
		    incr_kmer2(bed_itr, rst, rpos, kmer, K_UNDER, 0, qqual);
#endif
		} else {
		    if (qbase == ambig(rbase, qbase)) {
			if (V>1)
			    printf("M  %ld\t%c %c\n", rpos, rbase, qbase);
			// TYPE: match-match
			//incr_kmer(bed_itr, rpos, kmer, KMAT, 1, qqual);
			incr_kmer2(bed_itr, rst, rpos, kmer, K_MAT_M, 1, qqual);
		    } else if (rbase == '*') {
			if (V)
			    printf("im %ld\t%c %c *\n", rpos, rbase, qbase);
#ifndef MATCH_ONLY
			// TYPE: overcall
			//incr_kmer(bed_itr, rpos, kmer, KINS, 0, qqual);
			incr_kmer2(bed_itr, rst, rpos, kmer, K_OVER, 0, qqual);
#endif
		    } else if (rbase != 'N') {
			if (V) {
			    if (V<2)
				printf("%.5s\t%05o\t", context_s(b, qpos), kmer);
			    printf("X  %ld\t%c %c %d *\n", rpos, rbase, qbase, qual[qpos]);
			}

			// TYPE: substitution
			//incr_kmer(bed_itr, rpos, kmer, KMAT, 0, qqual);
			incr_kmer2(bed_itr, rst, rpos, kmer, K_MIS_M, 0, qqual);
		    }
		}

		// Check if map[rpos] and map[rpos+1] are more than 1 apart,
		// indicating insertion in cons which is absent in seq.
		//
		// Further more if we have a 5bp insertion in cons and it's
		// absent here, then we have 5 deletion kmers to mark up,
		// not just the one.
		// TODO: (we should have started earlier then)

		// TODO ref pos 9918843 on ONT test.
		// Here cons is G(*/T)TTTTTTC, so 1 het ins and
		// map[rpos+1]-map[rpos]-1 = 1 (1bp extra).
		// BUT: our next cig op is 1D, so we should have +0 or +1
		// base insertion, but have 1 del.
		//
		// Ie if this was HOM INS, it'd be ndel=2, despite
		// map[rpos+1]-map[rpos]-1 being 1, as we have to add
		// the D in too.  Maybe that's covered by CDEL below?
		if (ndel || (map[rpos+1] - map[rpos] != 1 && i < ncig &&
			     !(bam_cigar_op(cig[i+1]) == BAM_CINS ||
			       bam_cigar_op(cig[i+1]) == BAM_CPAD))) {
		    // Nominal deletion length
		    if (!ndel)
			ndel = map[rpos+1] - map[rpos] - 1 - WIN_LEN/2;
		    ndel-=(ndel>0);

		    // 2nd alternative length based on heterozyugous ins
		    // to consensus
		    int ndel2 = 0, k = rpos+1;
		    while (ref[k] && isins(ref[k]))
			ndel2 += isins_hom(ref[k++]);
		    if (V)
			printf("Ins of = %d,%d in ref, absent in seq\n",
			       ndel, ndel2);
		    // Correct?
		    // May wish to check the one that closest matches the
		    // next cigar op.
		    ndel = MIN(ndel,ndel2);
		    if (ndel < 0) ndel = 0; // seq-D next to cons-I

		    if (ndel) {
			if (V)
			    printf("mD %ld\t. %c *\n", rpos,  qbase);
#ifndef MATCH_ONLY
			// TYPE: undercall (cons had insertion, we did not)
			//incr_kmer(bed_itr, rpos, kmer, KDEL, 0, qqual);
			incr_kmer2(bed_itr, rst, rpos, kmer, K_UNDER, 0, qqual);
#endif
		    }
		}

		qpos++;
		rpos++;
		break;

	    case BAM_CINS:
		if (isins(rbase)) {
		    if (qbase == ambig(rbase, qbase)) {
			if (V>1) printf("mi %ld\t%c %c\n", rpos, rbase, qbase);
#ifndef MATCH_ONLY
			// TYPE: ins-match
			//incr_kmer(bed_itr, rpos, kmer, KINS, 1, qqual);
			incr_kmer2(bed_itr, rst, rpos, kmer, K_MAT_I, 1, qqual);
#endif
		    } else {
			if (V>1) printf("xi %ld\t%c %c\n", rpos, rbase, qbase);
#ifndef MATCH_ONLY
			// TYPE: ins-substitution
			//incr_kmer(bed_itr, rpos, kmer, KINS, 0, qqual);
			incr_kmer2(bed_itr, rst, rpos, kmer, K_MIS_I, 0, qqual);
#endif
		    }
		} else {
		    if (V) printf("I  %ld\t. %c *\n", rpos, qbase);
#ifndef MATCH_ONLY
		    // TYPE: overcall (no insertion in cons)
		    //incr_kmer(bed_itr, rpos, kmer, KINS, 0, qqual);
		    incr_kmer2(bed_itr, rst, rpos, kmer, K_OVER, 0, qqual);
#endif
		}
		qpos++;
		nth++;
		break;

	    case BAM_CDEL:
		if (rbase == '*' || islower(rbase)) {
		    if (V>1) printf("md %ld\t%c .\n", rpos, rbase);
#ifndef MATCH_ONLY
		    // TYPE: del-match
		    //incr_kmer(bed_itr, rpos, kmer, KDEL, 1, qqual);
		    incr_kmer2(bed_itr, rst, rpos, kmer, K_MAT_D, 1, qqual);
#endif
		} else {
		    if (V) printf("D  %ld\t%c . *\n", rpos, rbase);
#ifndef MATCH_ONLY
		    // TYPE: undercall
		    //incr_kmer(bed_itr, rpos, kmer, KDEL, 0, qqual);
		    incr_kmer2(bed_itr, rst, rpos, kmer, K_UNDER, 0, qqual);
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
	//
	// Q: should we still do this however for I and D states?
	int cnt =
	    kmer_count[i][0][0]+kmer_count[i][0][1] +
	    kmer_count[i][1][0]+kmer_count[i][1][1] +
	    kmer_count[i][2][0]+kmer_count[i][2][1];
	if (!cnt)
	    continue;

	int all_mis = 0;
	double all_err = 0;
	for (k = 0; k < 3; k++) {
	    if (!(kmer_count[i][k][0]+kmer_count[i][k][1]))
		continue;

	    for (j = WIN_LEN-1; j >= 0; j--)
		putchar("ACGTNNNN"[(i>>(j*3))&7]);
	    double err = cnt ? (double)kmer_count[i][k][0] / cnt : 0;
	    int qval = err ? -10*log10(err)+.5 : 99;
	    printf("\t%c\t%12d\t%12d\t%d\t%d\n",
		   "MID"[k], kmer_count[i][k][0], cnt, qval,
		   (int)(-4.343*log(kmer_qual[i][k]/cnt)+.5));

	    all_mis += kmer_count[i][k][0];
	    all_err += kmer_qual[i][k];
	}

	// Combined stats
	for (j = WIN_LEN-1; j >= 0; j--)
	    putchar("ACGTNNNN"[(i>>(j*3))&7]);
	double err = all_mis ? (double)all_mis / cnt : 0;
	int qval = err ? -10*log10(err)+.5 : 99;
	printf("\t?\t%12d\t%12d\t%d\t%d\n", all_mis, cnt, qval,
	       (int)(-4.343*log(all_err/cnt)+.5));
    }
}

// TODO: compare k_count[i][M_MAT_M] with K_MAT_I/D to see if the ratios
// of false to true agree.  The discrepancy represents possible alignment
// artefacts (ref bias, errors) or consensus call inaccuracies.
void dump_kmers2(void) {
    int i, j, k;
    puts("=== kmers ===\n");
    for (i = 0; i <= WIN_MASK; i++) {
	// For now skip kmers containing N.  We gathered stats, but I'm
	// not sure what to do with them.
	for (j = 0; j < WIN_LEN; j++)
	    if (((i >> (j*3))&4) == 4)
		break;
	if (j < WIN_LEN)
	    continue;

	// K_MAT_[MID] are all matches, but for M, I and D cigar ops.
	// (Eg cons is I and seq is I and sequences match).
	// So this is the total of "good".
	//
	// K_MIS_[MI] are matching cigar ops (eg both cons/seq is I) but
	// a substitution error.
	//
	// So K_MIS_M/K_MAT_M vs K_MIS_I/K_MAT_I can give us an indication
	// of potential alignment reference bias.  Both sequence data and
	// sample are in agreement, but the reference is the thing changing.

	// K_OVER and K_UNDER represent cigar op issues where cons/seq
	// differ on the number of base calls
	int ntrue = k_count[i][K_MAT_M] 
	          + k_count[i][K_MAT_I]
	          + k_count[i][K_MAT_D];
	int nsubst= k_count[i][K_MIS_M]
	          + k_count[i][K_MIS_I];
	int nunder= k_count[i][K_UNDER];
	int nover = k_count[i][K_OVER];

	int ntot  = ntrue + nsubst + nunder + nover;

	if (!ntot)
	    continue;

	for (k = 0; k < 5; k++) {
	    int nerr, ncnt;
	    double qerr;
	    qerr = k_qual[i][K_MAT_M]
		+ k_qual[i][K_MAT_I]
		+ k_qual[i][K_MAT_D]
		+ k_qual[i][K_MIS_M]
		+ k_qual[i][K_MIS_I];
	    double qerr2 = 0;
	    for (int z = 99; z >= 0; z--) {
		// amortised phred error score
		qerr2 += Perr[z]*k_qual2[i][K_MAT_M][z];
		qerr2 += Perr[z]*k_qual2[i][K_MAT_I][z];
		qerr2 += Perr[z]*k_qual2[i][K_MAT_D][z];
		qerr2 += Perr[z]*k_qual2[i][K_MAT_M][z];
		qerr2 += Perr[z]*k_qual2[i][K_MIS_M][z];
		qerr2 += Perr[z]*k_qual2[i][K_MIS_I][z];
	    }
	    ncnt = ntot;
	    switch (k) {
	    case 0:
		nerr = nsubst;
		break;
	    case 1:
		nerr = nunder;
		qerr += k_qual[i][K_UNDER];
		for (int z = 99; z >= 0; z--)
		    qerr2 += Perr[z]*k_qual2[i][K_UNDER][z];
		break;
	    case 2:
		nerr = nover;
		qerr += k_qual[i][K_OVER];
		for (int z = 99; z >= 0; z--)
		    qerr2 += Perr[z]*k_qual2[i][K_OVER][z];
		break;


	    // Looking for reference bias.
	    // Compare substitution errors in places where sample and cons
	    // are in sync, but both cigar-Ms vs both cigar-Is.
	    // A cigar-I match here is where the consensus is insertion
	    // as well as sample (FIXME: limit to homozygous only?), but
	    // the reference doesn't have this base.
	    // If there is a difference, it's due to alignment artifacts
	    // instead.
	    case 3: // match in M cigar ops only
		qerr = qerr2 = 0;
		qerr = k_qual[i][K_MAT_M] + k_qual[i][K_MIS_M];
		for (int z = 99; z >= 0; z--) {
		    qerr2 += Perr[z]*k_qual2[i][K_MAT_M][z];
		    qerr2 += Perr[z]*k_qual2[i][K_MIS_M][z];
		}
		nerr = k_count[i][K_MIS_M];
		ncnt = k_count[i][K_MIS_M] + k_count[i][K_MAT_M];
		break;
	    case 4: // match in I cigar ops only
		qerr = qerr2 = 0;
		qerr = k_qual[i][K_MAT_I] + k_qual[i][K_MIS_I];
		for (int z = 99; z >= 0; z--) {
		    qerr2 += Perr[z]*k_qual2[i][K_MAT_I][z];
		    qerr2 += Perr[z]*k_qual2[i][K_MIS_I][z];
		}
		nerr = k_count[i][K_MIS_I];
		ncnt = k_count[i][K_MIS_I] + k_count[i][K_MAT_I];
		break;
	    }

	    if (!nerr) continue;

	    double err = (double)nerr / ncnt;
	    int qreal = nerr ? -10*log10(err)+.5 : 99;
	    int qcall = (int)(-4.343*log(qerr / ncnt)+.5);
	    int qcall2 = (int)(-4.343*log(qerr2 / ncnt)+.5);
	    char *s[] = {"MATCH", "UNDER", "OVER", "MAT_M", "MAT_I"};

	    printf("%05o ", i);
	    for (j = WIN_LEN-1; j >= 0; j--)
		putchar("ACGTNNNN"[(i>>(j*3))&7]);
	    printf("\t%s\t%12d\t%12d\t%d\t%d\t%d\n",
		   s[k], nerr, ncnt, qreal, qcall,qcall2);
	    
	}
    }
}

int main(int argc, char **argv) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    htsFormat in_fmt = {0};
    faidx_t *fai = NULL;
    int mark_ins = 0;
    hts_pos_t *map = NULL;
    bam1_t *b = bam_init1();
    uint8_t *ref = NULL;
    uint8_t *stat = NULL; // status bits
    regidx_t *bed = NULL;
    regitr_t *bed_itr = NULL;
    char *reg = NULL;
    hts_itr_t *itr = NULL;

    int c;
    while ((c = getopt(argc, argv, "I:f:b:r:m")) >= 0) {
	switch(c) {
	case 'I': hts_parse_format(&in_fmt, optarg); break;

	case 'f':
	    if (!(fai = fai_load(optarg))) {
		fprintf(stderr, "Failed to load consensus fasta\n");
		goto err;
	    }
	    break;

	case 'm':
	    mark_ins = 1;
	    break;
	    
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
	    if (stat)
		free(stat);
	    stat = (uint8_t *)calloc(ref_len, sizeof(*stat));
	    if (!ref || !stat)
		goto err;
	    if (map)
		free(map);
	    if (!(map = build_ref_map(ref, &ref_len, stat, mark_ins)))
		goto err;
	    last_tid = b->core.tid;
	}
	//puts(bam_get_qname(b));

	accumulate_kmers(hdr, ref, stat, map, b, bed, bed_itr, mark_ins);
    }
    if (r != -1)
	goto err;

    if (sam_close(fp) < 0)
	goto err;
    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    fai_destroy(fai);
    if (ref) free(ref);
    if (stat) free(stat);
    if (map) free(map);
    if (bed) regidx_destroy(bed);
    if (bed_itr) regitr_destroy(bed_itr);
    if (itr) hts_itr_destroy(itr);

    //dump_kmers();
    dump_kmers2();

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
