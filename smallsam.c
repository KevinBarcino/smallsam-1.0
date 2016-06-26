/***************************************************************************************************/
/*                                                                                                 */
/*   Program Name:     smallsam.c                                                                  */
/*                                                                                                 */
/*   Author:           Kevin Chow                                                                  */
/*                                                                                                 */
/*   Date:             March 11, 2015                                                              */
/*                                                                                                 */
/*   Version:          1.0                                                                         */
/*                                                                                                 */
/*   Purpose:          Parse a BAM file, retrieve the data from the position, length and           */
/*                     reference columns and save them in csv format.                              */
/*                                                                                                 */
/*                     Functions in this program were compiled from the following samtools files:  */
/*                     (1) ./bamtk.c                                                               */
/*                     (2) ./sam_view.c                                                            */
/*                     (3) ./htslib-1.2.1/bgzf.c                                                   */
/*                     (4) ./htslib-1.2.1/hfile.c                                                  */
/*                     (5) ./htslib-1.2.1/hts.c                                                    */
/*                     (6) ./htslib-1.2.1/kstring.c                                                */
/*                     (7) ./htslib-1.2.1/knetfile.c                                               */
/*                     (8) ./htslib-1.2.1/hfile_net.c                                              */
/*                     (9) ./htslib-1.2.1/sam.c                                                    */
/*                                                                                                 */
/*                     Compile example (include the header files listed in the code):              */
/*                     gcc -Ofast -pthread -Wall smallsam.c -o smallsam -lz                       */
/*                                                                                                 */
/*                     Run example:                                                                */
/*                     ./smallsam <in.bam>                                                         */
/*                                                                                                 */
/***************************************************************************************************/

#define SMALLSAM

#include <pthread.h>
#include <errno.h>
#include <assert.h>
#include <inttypes.h> // PRId64
#include <unistd.h>
//#ifndef _WIN32
#include <netdb.h>
//#endif

#include "./htslib-1.2.1/htslib/bgzf.h"
#include "./htslib-1.2.1/htslib/sam.h"
#include "./htslib-1.2.1/htslib/hts.h"
#include "./htslib-1.2.1/htslib/hfile.h"
#include "./htslib-1.2.1/htslib/khash.h"
#include "./htslib-1.2.1/htslib/kstring.h"
#include "./htslib-1.2.1/htslib/kseq.h"
#include "./htslib-1.2.1/htslib/knetfile.h"
#include "./htslib-1.2.1/hfile_internal.h"



static FILE *pFile2 = NULL;
static char alignmentData[1000] = {0};


///////////////////////////////////////////////////////////
//                                                       //
//                  (0)  main() function                 //
//                                                       //
///////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	bam_hdr_t bamheader;
	FILE *pFile1 = NULL;
	int c = 0;
	int count = 0;
	int line = 0;
	int n = 0;
	char *data = NULL;
	char buffer[3000] = {0};
	char file1[1024] = {0};
	char file2[1024] = {0};
	char filenameBuf[1024] = {0};
	int  bytesRead = 0;
	int  bytesWritten = 0;
	long totalBytesRead = 0;
	int  currBuffer = 0;
	char *pos1 = NULL;
	char *pos2 = NULL;
	
	char *msg1 = "Loading file again\n";
	
	if (argc == 1 || argc > 3)
	{
		fprintf(stderr, "\nMissing file name or too many parameters.\nUsage:   smallsam <in.bam>\n\n");
		return 1;
	}

  int ret = 0;	

  //if (strcmp(argv[1], "view") == 0)
    //ret = main_samview(argc-1, argv+1);
  ret = main_samview(argc, argv);

//char *pos = NULL;
//strcpy(filenameBuf, argv[1]);
//strcpy(file1, filenameBuf);


/*
	pFile1 = fopen(file1, "r");
	if (pFile1 == NULL)
	{	
		sprintf(buffer, "Error opening file %s", file1);
		perror(buffer);
		return 1;	
	}

	if (bytesRead = fread(&bamheader, sizeof (bamheader), 1, pFile1))
	{
		printf("Bytes Read: %d\n", bytesRead);
		printf("%d\n", bamheader.n_targets);
		printf("%d\n", bamheader.ignore_sam_err);
		printf("%d\n", sizeof (int32_t));
		printf("%d\n", sizeof (char *));
		printf("%d\n", sizeof (void *));
	}

	fclose(pFile1);
*/
	return ret;
}

///////////////////////////////////////////////////////////
//                                                       //
//                   (1)  bamtk.c                        //
//                                                       //
///////////////////////////////////////////////////////////
void print_error(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "samtools: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

void print_error_errno(const char *format, ...)
{
    int err = errno;
    va_list args;
    va_start(args, format);
    fprintf(stderr, "samtools: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, ": %s\n", strerror(err));
    va_end(args);
}


///////////////////////////////////////////////////////////
//                                                       //
//                 (2)  sam_view.c                       //
//                                                       //
///////////////////////////////////////////////////////////
KHASH_SET_INIT_STR(rg)
typedef khash_t(rg) *rghash_t;
// This structure contains the settings for a samview run
typedef struct samview_settings {
    rghash_t rghash;
    int min_mapQ;
    int flag_on;
    int flag_off;
    int min_qlen;
    int remove_B;
    uint32_t subsam_seed;
    double subsam_frac;
    char* library;
    void* bed;
    size_t remove_aux_len;
    char** remove_aux;
} samview_settings_t;

static int usage(int is_long_help)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   smallsam <in.bam>\n\n");
/*
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]\n\n");
    // output options
    fprintf(stderr, "Options: -b       output BAM\n");
    fprintf(stderr, "         -C       output CRAM (requires -T)\n");
    fprintf(stderr, "         -1       use fast BAM compression (implies -b)\n");
    fprintf(stderr, "         -u       uncompressed BAM output (implies -b)\n");
    fprintf(stderr, "         -h       include header in SAM output\n");
    fprintf(stderr, "         -H       print SAM header only (no alignments)\n");
    fprintf(stderr, "         -c       print only the count of matching records\n");
    fprintf(stderr, "         -o FILE  output file name [stdout]\n");
    fprintf(stderr, "         -U FILE  output reads not selected by filters to FILE [null]\n");
    // extra input
    fprintf(stderr, "         -t FILE  FILE listing reference names and lengths (see long help) [null]\n");
    fprintf(stderr, "         -T FILE  reference sequence FASTA FILE [null]\n");
    // read filters
    fprintf(stderr, "         -L FILE  only include reads overlapping this BED FILE [null]\n");
    fprintf(stderr, "         -r STR   only include reads in read group STR [null]\n");
    fprintf(stderr, "         -R FILE  only include reads with read group listed in FILE [null]\n");
    fprintf(stderr, "         -q INT   only include reads with mapping quality >= INT [0]\n");
    fprintf(stderr, "         -l STR   only include reads in library STR [null]\n");
    fprintf(stderr, "         -m INT   only include reads with number of CIGAR operations\n");
    fprintf(stderr, "                  consuming query sequence >= INT [0]\n");
    fprintf(stderr, "         -f INT   only include reads with all bits set in INT set in FLAG [0]\n");
    fprintf(stderr, "         -F INT   only include reads with none of the bits set in INT\n");
    fprintf(stderr, "                  set in FLAG [0]\n");
    // read processing
    fprintf(stderr, "         -x STR   read tag to strip (repeatable) [null]\n");
    fprintf(stderr, "         -B       collapse the backward CIGAR operation\n");
    fprintf(stderr, "         -s FLOAT integer part sets seed of random number generator [0];\n");
    fprintf(stderr, "                  rest sets fraction of templates to subsample [no subsampling]\n");
    // general options
    fprintf(stderr, "         -@ INT   number of BAM compression threads [0]\n");
    fprintf(stderr, "         -?       print long help, including note about region specification\n");
    fprintf(stderr, "         -S       ignored (input format is auto-detected)\n");
    fprintf(stderr, "\n");
    if (is_long_help)
        fprintf(stderr, "Notes:\n\
\n\
  1. This command now auto-detects the input format (BAM/CRAM/SAM).\n\
\n\
  2. The file supplied with `-t' is SPACE/TAB delimited with the first\n\
     two fields of each line consisting of the reference name and the\n\
     corresponding sequence length. The `.fai' file generated by \n\
     `samtools faidx' is suitable for use as this file. This may be an\n\
     empty file if reads are unaligned.\n\
\n\
  3. SAM->BAM conversion: `samtools view -bT ref.fa in.sam.gz'.\n\
\n\
  4. BAM->SAM conversion: `samtools view -h in.bam'.\n\
\n\
  5. A region should be presented in one of the following formats:\n\
     `chr1', `chr2:1,000' and `chr3:1000-2,000'. When a region is\n\
     specified, the input alignment file must be a sorted and indexed\n\
     alignment (BAM/CRAM) file.\n\
\n\
  6. Option `-u' is preferred over `-b' when the output is piped to\n\
     another samtools command.\n\
\n");
*/
    return 1;
}

static void check_sam_close(samFile *fp, const char *fname, const char *null_fname, int *retp)
{
    int r = sam_close(fp);
    if (r >= 0) return;

    // TODO Need error infrastructure so we can print a message instead of r
    if (fname) print_error("error closing \"%s\": %d", fname, r);
    else print_error("error closing %s: %d", null_fname, r);

    *retp = EXIT_FAILURE;
}

int main_samview(int argc, char *argv[])
{
/************************************************************/  
FILE *pFile1 = NULL;
char file1[1024] = {0};
char file2[1024] = {0};
char filenameBuf[1024] = {0};
char *pos1 = NULL;
char *pos2 = NULL;
/************************************************************/
    //char alignmentData[1000] = {0};
    int c, is_header = 0, is_header_only = 0, ret = 0, compress_level = -1, is_count = 0;
    int is_long_help = 0, n_threads = 0;
    int64_t count = 0;
    samFile *in = 0, *out = 0, *un_out=0;
    bam_hdr_t *header = NULL;
    char out_mode[5], *out_format = "", *fn_out = 0, *fn_list = 0, *fn_ref = 0, *q, *fn_un_out = 0;

    samview_settings_t settings = {
        .rghash = NULL,
        .min_mapQ = 0,
        .flag_on = 0,
        .flag_off = 0,
        .min_qlen = 0,
        .remove_B = 0,
        .subsam_seed = 0,
        .subsam_frac = -1.,
        .library = NULL,
        .bed = NULL,
    };

    /* parse command-line options */
    /* TODO: convert this to getopt_long we're running out of letters */
    strcpy(out_mode, "w");
    while ((c = getopt(argc, argv, "SbBcCt:h1Ho:q:f:F:ul:r:?T:R:L:s:@:m:x:U:")) >= 0) {
        switch (c) {
        case 's':
            if ((settings.subsam_seed = strtol(optarg, &q, 10)) != 0) {
                srand(settings.subsam_seed);
                settings.subsam_seed = rand();
            }
            settings.subsam_frac = strtod(q, &q);
            break;
        case 'm': settings.min_qlen = atoi(optarg); break;
        case 'c': is_count = 1; break;
        case 'S': break;
        case 'b': out_format = "b"; break;
        case 'C': out_format = "c"; break;
        case 't': fn_list = strdup(optarg); break;
        case 'h': is_header = 1; break;
        case 'H': is_header_only = 1; break;
        case 'o': fn_out = strdup(optarg); break;
        case 'U': fn_un_out = strdup(optarg); break;
        case 'f': settings.flag_on |= strtol(optarg, 0, 0); break;
        case 'F': settings.flag_off |= strtol(optarg, 0, 0); break;
        case 'q': settings.min_mapQ = atoi(optarg); break;
        case 'u': compress_level = 0; break;
        case '1': compress_level = 1; break;
        case 'l': settings.library = strdup(optarg); break;
/*        
        case 'L':
            if ((settings.bed = bed_read(optarg)) == NULL) {
                print_error_errno("Could not read file \"%s\"", optarg);
                ret = 1;
                goto view_end;
            }
            break;
*/
/*        
          case 'r':
            if (add_read_group_single(&settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            break;
*/
/*
        case 'R':
            if (add_read_groups_file(&settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            break;
*/
                /* REMOVED as htslib doesn't support this
        //case 'x': out_format = "x"; break;
        //case 'X': out_format = "X"; break;
                 */
        case '?': is_long_help = 1; break;
        case 'T': fn_ref = strdup(optarg); break;
        case 'B': settings.remove_B = 1; break;
        case '@': n_threads = strtol(optarg, 0, 0); break;
        case 'x':
            {
                if (strlen(optarg) != 2) {
                    fprintf(stderr, "main_samview: Error parsing -x auxiliary tags should be exactly two characters long.\n");
                    return usage(is_long_help);
                }
                settings.remove_aux = (char**)realloc(settings.remove_aux, sizeof(char*) * (++settings.remove_aux_len));
                settings.remove_aux[settings.remove_aux_len-1] = optarg;
            }
            break;
        default: return usage(is_long_help);
        }
    }
    if (compress_level >= 0) out_format = "b";
    if (is_header_only) is_header = 1;
    strcat(out_mode, out_format);
    if (compress_level >= 0) {
        char tmp[2];
        tmp[0] = compress_level + '0'; tmp[1] = '\0';
        strcat(out_mode, tmp);
    }
    if (argc == optind) return usage(is_long_help); // potential memory leak...

    // generate the fn_list if necessary
/* commented out
    if (fn_list == 0 && fn_ref) fn_list = samfaipath(fn_ref);
*/
//printf("sam_open ***************************************%s\n", argv[optind]);
strcpy(file1, argv[optind]);
/****************************************/
strcpy(filenameBuf, argv[optind]);
pos1 = strrchr(filenameBuf, '.');
if (pos1 != NULL) {
	*pos1 = '\0';
}	
// extract last file name path
pos2 = strrchr(filenameBuf, '/');
if (pos2 != NULL && pos2+1 != '\0') {
  strcpy(file2, pos2+1);
}
else {
  strcpy(file2, filenameBuf);
}
strcat(file2, ".csv");

printf("[smallsam] Processing: %s\n", file1);
printf("[smallsam] Saving data to: %s\n", file2);
/*
if( access( file2, F_OK ) != -1 )
{
    printf("The csv file '%s' already exists.\n\n", file2);
    ret = 1;
    goto view_end;
}
*/
pFile1 = fopen(file1, "r");
pFile2 = fopen(file2, "w");
/*****************************************/

    // open file handlers
    if ((in = sam_open(argv[optind], "r")) == 0) {

        print_error_errno("failed to open \"%s\" for reading", argv[optind]);
        ret = 1;
//printf("\nview_end1\n");
        goto view_end;
    }
    if (fn_list) hts_set_fai_filename(in, fn_list);
    if ((header = sam_hdr_read(in)) == 0) {
        fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", argv[optind]);
        ret = 1;
//printf("\nview_end2\n");
        goto view_end;
    }
/* commented out to expedite conversion
    if (settings.rghash) { // FIXME: I do not know what "bam_header_t::n_text" is for...
        char *tmp;
        int l;
        tmp = drop_rg(header->text, settings.rghash, &l);
        free(header->text);
        header->text = tmp;
        header->l_text = l;
    }
*/
    if (!is_count) {
        if ((out = sam_open(fn_out? fn_out : "-", out_mode)) == 0) {
            print_error_errno("failed to open \"%s\" for writing", fn_out? fn_out : "standard output");
            ret = 1;
//printf("\nview_end3\n");
            goto view_end;
        }
        if (fn_list) hts_set_fai_filename(out, fn_list);
        if (*out_format || is_header)  {
            if (sam_hdr_write(out, header) != 0) {
                fprintf(stderr, "[main_samview] failed to write the SAM header\n");
                ret = 1;
//printf("\nview_end4\n");
                goto view_end;
            }
        }
        if (fn_un_out) {
            if ((un_out = sam_open(fn_un_out, out_mode)) == 0) {
                print_error_errno("failed to open \"%s\" for writing", fn_un_out);
                ret = 1;
//printf("\nview_end5\n");
                goto view_end;
            }
            if (*out_format || is_header) {
                if (sam_hdr_write(un_out, header) != 0) {
                    fprintf(stderr, "[main_samview] failed to write the SAM header\n");
                    ret = 1;
//printf("\nview_end6\n");
                    goto view_end;
                }
            }
        }
    }
    //commented out to expedite conversion
    if (n_threads > 1) { if (out) hts_set_threads(out, n_threads); }

    if (is_header_only) goto view_end; // no need to print alignments

//    int startCount = 0;

//printf("\nhere 1\n");
    if (argc == optind + 1) { // convert/print the entire file
//printf("\nhere 2\n");
        bam1_t *b = bam_init1();
        int r;

        int tid = 0;
        int pos = 0;
        int len = 0;
        int match = 0;
        while ((r = sam_read1(in, header, b, &tid, &pos, &match, &len)) >= 0) { // read one alignment from `in'
          if (match && tid >= 0)
          {
            sprintf(alignmentData, "%d,%d,%s\n", pos, pos + len, header->target_name[tid]);
            fwrite(alignmentData, 1, strlen(alignmentData), pFile2);
          }
          tid = 0;
          pos = 0;
          len = 0;
          match = 0;
            
/*
            if (!process_aln(header, b, &settings))
            {
                if (!is_count)
                {                       
                    if (check_sam_write1(out, header, b, fn_out, &ret) < 0)
                        break;
                   
                }
                count++;
            }
            else 
            {
                if (un_out) { if (check_sam_write1(un_out, header, b, fn_un_out, &ret) < 0) break; }
            }
*/
        }
        if (r < -1) {
            fprintf(stderr, "[main_samview] truncated file.\n");
            ret = 1;
        }
//printf("\nbam_destroy\n");
        bam_destroy1(b);
    }/*
    else
    { // retrieve alignments in specified regions
//printf("\nhere 5\n");
        int i;
        bam1_t *b;
        hts_idx_t *idx = sam_index_load(in, argv[optind]); // load index
        if (idx == 0) { // index is unavailable
            fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.\n");
            ret = 1;
            goto view_end;
        }
        b = bam_init1();
        for (i = optind + 1; i < argc; ++i) {
            int result;
            hts_itr_t *iter = sam_itr_querys(idx, header, argv[i]); // parse a region in the format like `chr2:100-200'
            if (iter == NULL) { // reference name is not found
                fprintf(stderr, "[main_samview] region \"%s\" specifies an unknown reference name. Continue anyway.\n", argv[i]);
                continue;
            }
            // fetch alignments
            while ((result = sam_itr_next(in, iter, b)) >= 0) {

                /////if (!process_aln(header, b, &settings)) {
                /////    if (!is_count) { if (check_sam_write1(out, header, b, fn_out, &ret) < 0) break; }
                /////    count++;
                /////} else {
                /////    if (un_out) { if (check_sam_write1(un_out, header, b, fn_un_out, &ret) < 0) break; }
                /////}

            }
            hts_itr_destroy(iter);
            if (result < -1) {
                fprintf(stderr, "[main_samview] retrieval of region \"%s\" failed due to truncated file or corrupt BAM index file\n", argv[i]);
                ret = 1;
                break;
            }
        }
        bam_destroy1(b);
        hts_idx_destroy(idx); // destroy the BAM index


    }*/
//printf("\nfclose(pFile1)\n");
	fclose(pFile1);
	fclose(pFile2);
view_end:
//printf("\nview_end\n");
    if (is_count && ret == 0)
        printf("%" PRId64 "\n", count);

    // close files, free and return
    if (in) check_sam_close(in, argv[optind], "standard input", &ret);
    if (out) check_sam_close(out, fn_out, "standard output", &ret);
    if (un_out) check_sam_close(un_out, fn_un_out, "file", &ret);

    free(fn_list); free(fn_ref); free(fn_out); free(settings.library);  free(fn_un_out);
/* commented out to expedite conversion
    if ( header ) bam_hdr_destroy(header);
*/
/* commented out to expedite conversion
    if (settings.bed) bed_destroy(settings.bed);
*/
    if (settings.rghash) {
        khint_t k;
        for (k = 0; k < kh_end(settings.rghash); ++k)
            if (kh_exist(settings.rghash, k)) free((char*)kh_key(settings.rghash, k));
        kh_destroy(rg, settings.rghash);
    }
    if (settings.remove_aux_len) {
        free(settings.remove_aux);
    }
    return ret;
}


///////////////////////////////////////////////////////////
//                                                       //
//            (3)  htslib-1.2.1/bgzf.c                   //
//                                                       //
///////////////////////////////////////////////////////////
#define BGZF_CACHE

#ifdef BGZF_CACHE
typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
#include "./htslib-1.2.1/htslib/khash.h"
KHASH_MAP_INIT_INT64(cache, cache_t)
#endif

static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";
static inline void packInt16(uint8_t *buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

static inline void packInt32(uint8_t *buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

#define KS_BGZF 1
#if KS_BGZF
    // bgzf now supports gzip-compressed files, the gzFile branch can be removed
    KSTREAM_INIT2(, BGZF*, bgzf_read, 65536)
#else
    KSTREAM_INIT2(, gzFile, gzread, 16384)
#endif

KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

#define BGZF_CACHE

/* BGZF/GZIP header (speciallized from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  B   C

  BGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^16 bytes and adds and an extra "BC" field in the gzip header which
  records the size.

*/

typedef struct
{
    uint64_t uaddr;  // offset w.r.t. uncompressed data
    uint64_t caddr;  // offset w.r.t. compressed data
}
bgzidx1_t;

struct __bgzidx_t
{
    int noffs, moffs;       // the size of the index, n:used, m:allocated
    bgzidx1_t *offs;        // offsets
    uint64_t ublock_addr;   // offset of the current block (uncompressed data)
};


typedef struct {
    struct bgzf_mtaux_t *mt;
    void *buf;
    int i, errcode, toproc, compress_level;
} worker_t;

typedef struct bgzf_mtaux_t {
    int n_threads, n_blks, curr, done;
    volatile int proc_cnt;
    void **blk;
    int *len;
    worker_t *w;
    pthread_t *tid;
    pthread_mutex_t lock;
    pthread_cond_t cv;
} mtaux_t;


static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}

// Returns: 0 on success (BGZF header); -1 on non-BGZF GZIP header; -2 on error
static int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}

static int bgzf_compress(void *_dst, int *dlen, void *src, int slen, int level)
{
    uint32_t crc;
    z_stream zs;
    uint8_t *dst = (uint8_t*)_dst;

    // compress the body
    zs.zalloc = NULL; zs.zfree = NULL;
    zs.next_in  = (Bytef*)src;
    zs.avail_in = slen;
    zs.next_out = dst + BLOCK_HEADER_LENGTH;
    zs.avail_out = *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
    if (deflateInit2(&zs, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) return -1; // -15 to disable zlib header/footer
    if (deflate(&zs, Z_FINISH) != Z_STREAM_END) return -1;
    if (deflateEnd(&zs) != Z_OK) return -1;
    *dlen = zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    // write the header
    memcpy(dst, g_magic, BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
    packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes
    // write the footer
    crc = crc32(crc32(0L, NULL, 0L), (Bytef*)src, slen);
    packInt32((uint8_t*)&dst[*dlen - 8], crc);
    packInt32((uint8_t*)&dst[*dlen - 4], slen);
    return 0;
}


static int worker_aux(worker_t *w)
{
    int i, stop = 0;
    // wait for condition: to process or all done
    pthread_mutex_lock(&w->mt->lock);
    while (!w->toproc && !w->mt->done)
        pthread_cond_wait(&w->mt->cv, &w->mt->lock);
    if (w->mt->done) stop = 1;
    w->toproc = 0;
    pthread_mutex_unlock(&w->mt->lock);
    if (stop) return 1; // to quit the thread
    w->errcode = 0;
    for (i = w->i; i < w->mt->curr; i += w->mt->n_threads) {
        int clen = BGZF_MAX_BLOCK_SIZE;
        if (bgzf_compress(w->buf, &clen, w->mt->blk[i], w->mt->len[i], w->compress_level) != 0)
            w->errcode |= BGZF_ERR_ZLIB;
        memcpy(w->mt->blk[i], w->buf, clen);
        w->mt->len[i] = clen;
    }
    __sync_fetch_and_add(&w->mt->proc_cnt, 1);
    return 0;
}
static void *mt_worker(void *data)
{
    while (worker_aux((worker_t*)data) == 0);
    return 0;
}

int bgzf_check_EOF(BGZF *fp)
{
    uint8_t buf[28];
    off_t offset = htell(fp->fp);
    if (hseek(fp->fp, -28, SEEK_END) < 0) {
        if (errno == ESPIPE) { hclearerr(fp->fp); return 2; }
        else return -1;
    }
    if ( hread(fp->fp, buf, 28) != 28 ) return -1;
    if ( hseek(fp->fp, offset, SEEK_SET) < 0 ) return -1;
    return (memcmp("\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", buf, 28) == 0)? 1 : 0;
}

ssize_t bgzf_read(BGZF *fp, void *data, size_t length)
{
//printf ("**** Entering function bgzf_read() in bgzf.c\n");
//printf("length: %d\n", length);
    ssize_t bytes_read = 0;
    uint8_t *output = (uint8_t*)data;
    if (length <= 0) return 0;
    assert(fp->is_write == 0);
    while (bytes_read < length) {
//printf("fp->block_length: %d, fp->block_offset: %d\n", fp->block_length, fp->block_offset);
        int copy_length, available = fp->block_length - fp->block_offset;
        uint8_t *buffer;
        if (available <= 0) {
            if (bgzf_read_block(fp) != 0) return -1;
            available = fp->block_length - fp->block_offset;
            if (available <= 0) break;
        }
        copy_length = length - bytes_read < available? length - bytes_read : available;
        buffer = (uint8_t*)fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;
    }
    if (fp->block_offset == fp->block_length) {
        fp->block_address = htell(fp->fp);
        fp->block_offset = fp->block_length = 0;
    }
    fp->uncompressed_address += bytes_read;
    return bytes_read;
}

static int lazy_flush(BGZF *fp)
{
    if (fp->mt) {
        ////if (fp->block_offset) mt_queue(fp);
        ////return (fp->mt->curr < fp->mt->n_blks)? 0 : mt_flush_queue(fp);
        return 0;
    }
    else return bgzf_flush(fp);
}

ssize_t bgzf_write(BGZF *fp, const void *data, size_t length)
{
    if ( !fp->is_compressed )
        return hwrite(fp->fp, data, length);

    const uint8_t *input = (const uint8_t*)data;
    ssize_t remaining = length;
    assert(fp->is_write);
    while (remaining > 0) {
        uint8_t* buffer = (uint8_t*)fp->uncompressed_block;
        int copy_length = BGZF_BLOCK_SIZE - fp->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (fp->block_offset == BGZF_BLOCK_SIZE) {
            if (lazy_flush(fp) != 0) return -1;
        }
    }
    return length - remaining;
}

int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)
{
    int i;
    mtaux_t *mt;
    pthread_attr_t attr;
    if (!fp->is_write || fp->mt || n_threads <= 1) return -1;
    mt = (mtaux_t*)calloc(1, sizeof(mtaux_t));
    mt->n_threads = n_threads;
    mt->n_blks = n_threads * n_sub_blks;
    mt->len = (int*)calloc(mt->n_blks, sizeof(int));
    mt->blk = (void**)calloc(mt->n_blks, sizeof(void*));
    for (i = 0; i < mt->n_blks; ++i)
        mt->blk[i] = malloc(BGZF_MAX_BLOCK_SIZE);
    mt->tid = (pthread_t*)calloc(mt->n_threads, sizeof(pthread_t)); // tid[0] is not used, as the worker 0 is launched by the master
    mt->w = (worker_t*)calloc(mt->n_threads, sizeof(worker_t));
    for (i = 0; i < mt->n_threads; ++i) {
        mt->w[i].i = i;
        mt->w[i].mt = mt;
        mt->w[i].compress_level = fp->compress_level;
        mt->w[i].buf = malloc(BGZF_MAX_BLOCK_SIZE);
    }
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_mutex_init(&mt->lock, 0);
    pthread_cond_init(&mt->cv, 0);
    for (i = 1; i < mt->n_threads; ++i) // worker 0 is effectively launched by the master thread
        pthread_create(&mt->tid[i], &attr, mt_worker, &mt->w[i]);
    fp->mt = mt;
    return 0;
}


static int inflate_gzip_block(BGZF *fp, int cached)
{
    int ret = Z_OK;
    do
    {
        if ( !cached && fp->gz_stream->avail_out!=0 )
        {
            fp->gz_stream->avail_in = hread(fp->fp, fp->compressed_block, BGZF_BLOCK_SIZE);
            if ( fp->gz_stream->avail_in<=0 ) return fp->gz_stream->avail_in;
            if ( fp->gz_stream->avail_in==0 ) break;
            fp->gz_stream->next_in = fp->compressed_block;
        }
        else cached = 0;
        do
        {
            fp->gz_stream->next_out = (Bytef*)fp->uncompressed_block + fp->block_offset;
            fp->gz_stream->avail_out = BGZF_MAX_BLOCK_SIZE - fp->block_offset;
            ret = inflate(fp->gz_stream, Z_NO_FLUSH);
            if ( ret==Z_BUF_ERROR ) continue;   // non-critical error
            if ( ret<0 ) return -1;
            unsigned int have = BGZF_MAX_BLOCK_SIZE - fp->gz_stream->avail_out;
            if ( have ) return have;
        }
        while ( fp->gz_stream->avail_out == 0 );
    }
    while (ret != Z_STREAM_END);
    return BGZF_MAX_BLOCK_SIZE - fp->gz_stream->avail_out;
}

// (bgzf_read_block)
int bgzf_index_add_block(BGZF *fp)
{
    fp->idx->noffs++;
    if ( fp->idx->noffs > fp->idx->moffs )
    {
        fp->idx->moffs = fp->idx->noffs;
        kroundup32(fp->idx->moffs);
        fp->idx->offs = (bgzidx1_t*) realloc(fp->idx->offs, fp->idx->moffs*sizeof(bgzidx1_t));
        if ( !fp->idx->offs ) return -1;
    }
    fp->idx->offs[ fp->idx->noffs-1 ].uaddr = fp->idx->ublock_addr;
    fp->idx->offs[ fp->idx->noffs-1 ].caddr = fp->block_address;
    return 0;
}

#ifdef BGZF_CACHE
static void free_cache(BGZF *fp)
{
    khint_t k;
    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    if (fp->is_write) return;
    for (k = kh_begin(h); k < kh_end(h); ++k)
        if (kh_exist(h, k)) free(kh_val(h, k).block);
    kh_destroy(cache, h);
}



static int load_block_from_cache(BGZF *fp, int64_t block_address)
{
    khint_t k;
    cache_t *p;
    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    k = kh_get(cache, h, block_address);
    if (k == kh_end(h)) return 0;
    p = &kh_val(h, k);
    if (fp->block_length != 0) fp->block_offset = 0;
    fp->block_address = block_address;
    fp->block_length = p->size;
    memcpy(fp->uncompressed_block, p->block, BGZF_MAX_BLOCK_SIZE);
    if ( hseek(fp->fp, p->end_offset, SEEK_SET) < 0 )
    {
        // todo: move the error up
        fprintf(stderr,"Could not hseek to %"PRId64"\n", p->end_offset);
        exit(1);
    }
    return p->size;
}

static void cache_block(BGZF *fp, int size)
{
    int ret;
    khint_t k;
    cache_t *p;
    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    if (BGZF_MAX_BLOCK_SIZE >= fp->cache_size) return;
    if ((kh_size(h) + 1) * BGZF_MAX_BLOCK_SIZE > (uint32_t)fp->cache_size) {
        /* A better way would be to remove the oldest block in the
         * cache, but here we remove a random one for simplicity. This
         * should not have a big impact on performance. */
        for (k = kh_begin(h); k < kh_end(h); ++k)
            if (kh_exist(h, k)) break;
        if (k < kh_end(h)) {
            free(kh_val(h, k).block);
            kh_del(cache, h, k);
        }
    }
    k = kh_put(cache, h, fp->block_address, &ret);
    if (ret == 0) return; // if this happens, a bug!
    p = &kh_val(h, k);
    p->size = fp->block_length;
    p->end_offset = fp->block_address + size;
    p->block = (uint8_t*)malloc(BGZF_MAX_BLOCK_SIZE);
    memcpy(kh_val(h, k).block, fp->uncompressed_block, BGZF_MAX_BLOCK_SIZE);
}
#else
static void free_cache(BGZF *fp) {}
static int load_block_from_cache(BGZF *fp, int64_t block_address) {return 0;}
static void cache_block(BGZF *fp, int size) {}
#endif

// Inflate the block in fp->compressed_block into fp->uncompressed_block
static int inflate_block(BGZF* fp, int block_length)
{
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = (Bytef*)fp->compressed_block + 18;
    zs.avail_in = block_length - 16;
    zs.next_out = (Bytef*)fp->uncompressed_block;
    zs.avail_out = BGZF_MAX_BLOCK_SIZE;

    if (inflateInit2(&zs, -15) != Z_OK) {
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
        inflateEnd(&zs);
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    if (inflateEnd(&zs) != Z_OK) {
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    return zs.total_out;
}

// get the compress level from the mode string: compress_level==-1 for the default level, -2 plain uncompressed
static int mode2level(const char *__restrict mode)
{
    int i, compress_level = -1;
    for (i = 0; mode[i]; ++i)
        if (mode[i] >= '0' && mode[i] <= '9') break;
    if (mode[i]) compress_level = (int)mode[i] - '0';
    if (strchr(mode, 'u')) compress_level = -2;
    return compress_level;
}

static BGZF *bgzf_read_init(hFILE *hfpr)
{
    BGZF *fp;
    uint8_t magic[18];
    ssize_t n = hpeek(hfpr, magic, 18);
    if (n < 0) return NULL;

    fp = (BGZF*)calloc(1, sizeof(BGZF));
    if (fp == NULL) return NULL;

    fp->is_write = 0;
    fp->is_compressed = (n==2 && magic[0]==0x1f && magic[1]==0x8b);
    fp->uncompressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
    fp->compressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
    fp->is_compressed = (n==18 && magic[0]==0x1f && magic[1]==0x8b) ? 1 : 0;
    fp->is_gzip = ( !fp->is_compressed || ((magic[3]&4) && memcmp(&magic[12], "BC\2\0",4)==0) ) ? 0 : 1;
#ifdef BGZF_CACHE
    fp->cache = kh_init(cache);
#endif
    return fp;
}

static BGZF *bgzf_write_init(const char *mode)
{
    BGZF *fp;
    fp = (BGZF*)calloc(1, sizeof(BGZF));
    fp->is_write = 1;
    int compress_level = mode2level(mode);
    if ( compress_level==-2 )
    {
        fp->is_compressed = 0;
        return fp;
    }
    fp->is_compressed = 1;
    fp->uncompressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
    fp->compressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
    fp->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1
    if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;
    if ( strchr(mode,'g') )
    {
        // gzip output
        fp->is_gzip = 1;
        fp->gz_stream = (z_stream*)calloc(1,sizeof(z_stream));
        fp->gz_stream->zalloc = NULL;
        fp->gz_stream->zfree  = NULL;
        if ( deflateInit2(fp->gz_stream, fp->compress_level, Z_DEFLATED, 15|16, 8, Z_DEFAULT_STRATEGY)!=Z_OK ) return NULL;
    }
    return fp;
}

static int bgzf_gzip_compress(BGZF *fp, void *_dst, int *dlen, void *src, int slen, int level)
{
    uint8_t *dst = (uint8_t*)_dst;
    z_stream *zs = fp->gz_stream;
    int flush = slen ? Z_NO_FLUSH : Z_FINISH;
    zs->next_in   = (Bytef*)src;
    zs->avail_in  = slen;
    zs->next_out  = dst;
    zs->avail_out = *dlen;
    if ( deflate(zs, flush) == Z_STREAM_ERROR ) return -1;
    *dlen = *dlen - zs->avail_out;
    return 0;
}


static int deflate_block(BGZF *fp, int block_length)
{
    int comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret;
    if ( !fp->is_gzip )
        ret = bgzf_compress(fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level);
    else
        ret = bgzf_gzip_compress(fp, fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level);

    if ( ret != 0 )
    {
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    fp->block_offset = 0;
    return comp_size;
}

int bgzf_read_block(BGZF *fp)
{
    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size = 0, block_length, remaining;

    // Reading an uncompressed file
    if ( !fp->is_compressed )
    {
        count = hread(fp->fp, fp->uncompressed_block, BGZF_MAX_BLOCK_SIZE);
        if ( count==0 )
        {
            fp->block_length = 0;
            return 0;
        }
        if (fp->block_length != 0) fp->block_offset = 0;
        fp->block_address += count;
        fp->block_length = count;
        return 0;
    }

    // Reading compressed file
    int64_t block_address;
    block_address = htell(fp->fp);
    if ( fp->is_gzip && fp->gz_stream ) // is this is a initialized gzip stream?
    {
        count = inflate_gzip_block(fp, 0);
        if ( count<0 )
        {
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->block_length = count;
        fp->block_address = block_address;
        return 0;
    }
    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
    count = hread(fp->fp, header, sizeof(header));
    if (count == 0) { // no data read
        fp->block_length = 0;
        return 0;
    }
    int ret;
    if ( count != sizeof(header) || (ret=check_header(header))==-2 )
    {
        fp->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    if ( ret==-1 )
    {
        // GZIP, not BGZF
        uint8_t *cblock = (uint8_t*)fp->compressed_block;
        memcpy(cblock, header, sizeof(header));
        count = hread(fp->fp, cblock+sizeof(header), BGZF_BLOCK_SIZE - sizeof(header)) + sizeof(header);
        int nskip = 10;

        // Check optional fields to skip: FLG.FNAME,FLG.FCOMMENT,FLG.FHCRC,FLG.FEXTRA
        // Note: Some of these fields are untested, I did not have appropriate data available
        if ( header[3] & 0x4 ) // FLG.FEXTRA
        {
            nskip += unpackInt16(&cblock[nskip]) + 2;
        }
        if ( header[3] & 0x8 ) // FLG.FNAME
        {
            while ( nskip<BGZF_BLOCK_SIZE && cblock[nskip] ) nskip++;
            if ( nskip==BGZF_BLOCK_SIZE )
            {
                fp->errcode |= BGZF_ERR_HEADER;
                return -1;
            }
            nskip++;
        }
        if ( header[3] & 0x10 ) // FLG.FCOMMENT
        {
            while ( nskip<BGZF_BLOCK_SIZE && cblock[nskip] ) nskip++;
            if ( nskip==BGZF_BLOCK_SIZE )
            {
                fp->errcode |= BGZF_ERR_HEADER;
                return -1;
            }
            nskip++;
        }
        if ( header[3] & 0x2 ) nskip += 2;  //  FLG.FHCRC

        fp->is_gzip = 1;
        fp->gz_stream = (z_stream*) calloc(1,sizeof(z_stream));
        int ret = inflateInit2(fp->gz_stream, -15);
        if (ret != Z_OK)
        {
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->gz_stream->avail_in = count - nskip;
        fp->gz_stream->next_in  = cblock + nskip;
        count = inflate_gzip_block(fp, 1);
        if ( count<0 )
        {
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->block_length = count;
        fp->block_address = block_address;
        if ( fp->idx_build_otf ) return -1; // cannot build index for gzip
        return 0;
    }
    size = count;
    block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
    compressed_block = (uint8_t*)fp->compressed_block;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
    if (count != remaining) {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }
    size += count;
    if ((count = inflate_block(fp, block_length)) < 0) return -1;
    if (fp->block_length != 0) fp->block_offset = 0; // Do not reset offset if this read follows a seek.
    fp->block_address = block_address;
    fp->block_length = count;
    if ( fp->idx_build_otf )
    {
        bgzf_index_add_block(fp);
        fp->idx->ublock_addr += count;
    }
    cache_block(fp, size);
    return 0;
}

BGZF *bgzf_hopen(hFILE *hfp, const char *mode)
{
    BGZF *fp = NULL;
    assert(compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);
    if (strchr(mode, 'r')) {
        fp = bgzf_read_init(hfp);
        if (fp == NULL) return NULL;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        fp = bgzf_write_init(mode);
    }
    else { errno = EINVAL; return 0; }

    fp->fp = hfp;
    fp->is_be = ed_is_big();
    return fp;
}

int bgzf_flush(BGZF *fp)
{
    if (!fp->is_write) return 0;
#ifdef BGZF_MT
    if (fp->mt) {
        if (fp->block_offset) mt_queue(fp); // guaranteed that assertion does not fail
        return mt_flush_queue(fp);
    }
#endif
    while (fp->block_offset > 0) {
        if ( fp->idx_build_otf )
        {
            bgzf_index_add_block(fp);
            fp->idx->ublock_addr += fp->block_offset;
        }
        int block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) return -1;
        if (hwrite(fp->fp, fp->compressed_block, block_length) != block_length) {
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }
        fp->block_address += block_length;
    }
    return 0;
}

void bgzf_index_destroy(BGZF *fp)
{
    if ( !fp->idx ) return;
    free(fp->idx->offs);
    free(fp->idx);
    fp->idx = NULL;
    fp->idx_build_otf = 0;
}

int bgzf_close(BGZF* fp)
{
    int ret, block_length;
    if (fp == 0) return -1;
    if (fp->is_write && fp->is_compressed) {
        if (bgzf_flush(fp) != 0) return -1;
        fp->compress_level = -1;
        block_length = deflate_block(fp, 0); // write an empty block
        if (hwrite(fp->fp, fp->compressed_block, block_length) < 0
            || hflush(fp->fp) != 0) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
#ifdef BGZF_MT
        if (fp->mt) mt_destroy(fp->mt);
#endif
    }
    if ( fp->is_gzip )
    {
        if (!fp->is_write) (void)inflateEnd(fp->gz_stream);
        else (void)deflateEnd(fp->gz_stream);
        free(fp->gz_stream);
    }
    ret = hclose(fp->fp);
    if (ret != 0) return -1;
    bgzf_index_destroy(fp);
    free(fp->uncompressed_block);
    free(fp->compressed_block);
    free_cache(fp);
    free(fp);
    return 0;
}



///////////////////////////////////////////////////////////
//                                                       //
//             (4)  htslib-1.2.1/hfile.c                 //
//                                                       //
///////////////////////////////////////////////////////////
static inline int writebuffer_is_nonempty(hFILE *fp)
{
    return fp->begin > fp->end;
}

/* Flushes the write buffer, fp->[buffer,begin), out through the backend
   returning 0 on success or negative if an error occurred.  */
static ssize_t flush_buffer(hFILE *fp)
{
    const char *buffer = fp->buffer;
    while (buffer < fp->begin) {
        ssize_t n = fp->backend->write(fp, buffer, fp->begin - buffer);
        if (n < 0) { fp->has_errno = errno; return n; }
        buffer += n;
        fp->offset += n;
    }

    fp->begin = fp->buffer;  // Leave the buffer empty
    return 0;
}

off_t hseek(hFILE *fp, off_t offset, int whence)
{
    off_t pos;

    if (writebuffer_is_nonempty(fp)) {
        int ret = flush_buffer(fp);
        if (ret < 0) return ret;
    }
    else {
        // Convert relative offsets from being relative to the hFILE's stream
        // position (at begin) to being relative to the backend's physical
        // stream position (at end, due to the buffering read-ahead).
        if (whence == SEEK_CUR) offset -= fp->end - fp->begin;
    }

    pos = fp->backend->seek(fp, offset, whence); ////////////////////////////////////////////////////
    if (pos < 0) { fp->has_errno = errno; return pos; }

    // Seeking succeeded, so discard any non-empty read buffer
    fp->begin = fp->end = fp->buffer;
    fp->at_eof = 0;

    fp->offset = pos;
    return pos;
}

/* Refills the read buffer from the backend (once, so may only partially
   fill the buffer), returning the number of additional characters read
   (which might be 0), or negative when an error occurred.  */
static ssize_t refill_buffer(hFILE *fp)
{
    ssize_t n;

    // Move any unread characters to the start of the buffer
    if (fp->begin > fp->buffer) {
        fp->offset += fp->begin - fp->buffer;
        memmove(fp->buffer, fp->begin, fp->end - fp->begin);
        fp->end = &fp->buffer[fp->end - fp->begin];
        fp->begin = fp->buffer;
    }

    // Read into the available buffer space at fp->[end,limit)
    if (fp->at_eof || fp->end == fp->limit) n = 0;
    else {
        n = fp->backend->read(fp, fp->end, fp->limit - fp->end);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
    }

    fp->end += n;
    return n;
}

/* Called only from hread(); when called, our buffer is empty and nread bytes
   have already been placed in the destination buffer.  */
ssize_t hread2(struct hFILE *fp, void *destv, size_t nbytes, size_t nread)
{
    const size_t capacity = fp->limit - fp->buffer;
    char *dest = (char *) destv;
    dest += nread, nbytes -= nread;

    // Read large requests directly into the destination buffer
    while (nbytes * 2 >= capacity && !fp->at_eof) {
        ssize_t n = fp->backend->read(fp, dest, nbytes);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
        fp->offset += n;
        dest += n, nbytes -= n;
        nread += n;
    }

    while (nbytes > 0 && !fp->at_eof) {
        size_t n;
        ssize_t ret = refill_buffer(fp);
        if (ret < 0) return ret;

        n = fp->end - fp->begin;
        if (n > nbytes) n = nbytes;
        memcpy(dest, fp->begin, n);
        fp->begin += n;
        dest += n, nbytes -= n;
        nread += n;
    }

    return nread;
}

// hfile.c
int hclose(hFILE *fp)
{
    int err = fp->has_errno;

    if (writebuffer_is_nonempty(fp) && hflush(fp) < 0) err = fp->has_errno;
    if (fp->backend->close(fp) < 0) err = errno;
    hfile_destroy(fp);

    if (err) {
        errno = err;
        return EOF;
    }
    else return 0;
}

// hfile.c
int hflush(hFILE *fp)
{
    if (flush_buffer(fp) < 0) return EOF;
    if (fp->backend->flush(fp) < 0) { fp->has_errno = errno; return EOF; }
    return 0;
}


// hfile.c
void hfile_destroy(hFILE *fp)
{
    int save = errno;
    if (fp) free(fp->buffer);
    free(fp);
    errno = save;
}

// hfile.c
typedef struct {
    hFILE base;
    const char *buffer;
    size_t length, pos;
} hFILE_mem;


typedef struct {
    hFILE base;
    int fd;
    int is_socket:1;
} hFILE_fd;

static ssize_t fd_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    ssize_t n;
    do {
        n = fp->is_socket? recv(fp->fd, buffer, nbytes, 0)
                         : read(fp->fd, buffer, nbytes);
    } while (n < 0 && errno == EINTR);
    return n;
}

static ssize_t fd_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    ssize_t n;
    do {
        n = fp->is_socket?  send(fp->fd, buffer, nbytes, 0)
                         : write(fp->fd, buffer, nbytes);
    } while (n < 0 && errno == EINTR);
    return n;
}

static off_t fd_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    return lseek(fp->fd, offset, whence);
}

static int fd_flush(hFILE *fpv)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    int ret;
    do {
#ifdef HAVE_FDATASYNC
        ret = fdatasync(fp->fd);
#else
        ret = fsync(fp->fd);
#endif
        // Ignore invalid-for-fsync(2) errors due to being, e.g., a pipe,
        // and operation-not-supported errors (Mac OS X)
        if (ret < 0 && (errno == EINVAL || errno == ENOTSUP)) ret = 0;
    } while (ret < 0 && errno == EINTR);
    return ret;
}

static int fd_close(hFILE *fpv)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    int ret;
    do {
#ifdef HAVE_CLOSESOCKET
        ret = fp->is_socket? closesocket(fp->fd) : close(fp->fd);
#else
        ret = close(fp->fd);
#endif
    } while (ret < 0 && errno == EINTR);
    return ret;
}

static const struct hFILE_backend fd_backend =
{
    fd_read, fd_write, fd_seek, fd_flush, fd_close
};

// hfile.c
static size_t blksize(int fd)
{
    struct stat sbuf;
    if (fstat(fd, &sbuf) != 0) return 0;
    return sbuf.st_blksize;
}


// hfile.c
int hfile_oflags(const char *mode)
{
    int rdwr = 0, flags = 0;
    const char *s;
    for (s = mode; *s; s++)
        switch (*s) {
        case 'r': rdwr = O_RDONLY;  break;
        case 'w': rdwr = O_WRONLY; flags |= O_CREAT | O_TRUNC;  break;
        case 'a': rdwr = O_WRONLY; flags |= O_CREAT | O_APPEND;  break;
        case '+': rdwr = O_RDWR;  break;
        default:  break;
        }

#ifdef O_BINARY
    flags |= O_BINARY;
#endif

    return rdwr | flags;
}

// hfile.c
hFILE *hfile_init(size_t struct_size, const char *mode, size_t capacity)
{
    hFILE *fp = (hFILE *) malloc(struct_size);
    if (fp == NULL) goto error;

    if (capacity == 0) capacity = 32768;
    // FIXME For now, clamp input buffer sizes so mpileup doesn't eat memory
    if (strchr(mode, 'r') && capacity > 32768) capacity = 32768;

    fp->buffer = (char *) malloc(capacity);
    if (fp->buffer == NULL) goto error;

    fp->begin = fp->end = fp->buffer;
    fp->limit = &fp->buffer[capacity];

    fp->offset = 0;
    fp->at_eof = 0;
    fp->has_errno = 0;
    return fp;

error:
    hfile_destroy(fp);
    return NULL;
}



// htslib-1.2.1/hfile.c
static hFILE *hopen_fd(const char *filename, const char *mode)
{
    hFILE_fd *fp = NULL;
    int fd = open(filename, hfile_oflags(mode), 0666);
    if (fd < 0) goto error;

    fp = (hFILE_fd *) hfile_init(sizeof (hFILE_fd), mode, blksize(fd));
    if (fp == NULL) goto error;

    fp->fd = fd;
    fp->is_socket = 0;
    fp->base.backend = &fd_backend;
    return &fp->base;

error:
    if (fd >= 0) { int save = errno; (void) close(fd); errno = save; }
    hfile_destroy((hFILE *) fp);
    return NULL;
}

// htslib-1.2.1/hfile.c
void hclose_abruptly(hFILE *fp)
{
    int save = errno;
    if (fp->backend->close(fp) < 0) { /* Ignore subsequent errors */ }
    hfile_destroy(fp);
    errno = save;
}
// hfile.c
ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes)
{
    size_t n = fp->end - fp->begin;
    while (n < nbytes) {
        ssize_t ret = refill_buffer(fp);
        if (ret < 0) return ret;
        else if (ret == 0) break;
        else n += ret;
    }

    if (n > nbytes) n = nbytes;
    memcpy(buffer, fp->begin, n);
    return n;
}


// hfile.c
/* Called only from hwrite() and hputs2(); when called, our buffer is full and
   ncopied bytes from the source have already been copied to our buffer.  */
ssize_t hwrite2(hFILE *fp, const void *srcv, size_t totalbytes, size_t ncopied)
{
    const char *src = (const char *) srcv;
    ssize_t ret;
    const size_t capacity = fp->limit - fp->buffer;
    size_t remaining = totalbytes - ncopied;
    src += ncopied;

    ret = flush_buffer(fp);
    if (ret < 0) return ret;

    // Write large blocks out directly from the source buffer
    while (remaining * 2 >= capacity) {
        ssize_t n = fp->backend->write(fp, src, remaining);
        if (n < 0) { fp->has_errno = errno; return n; }
        fp->offset += n;
        src += n, remaining -= n;
    }

    // Just buffer any remaining characters
    memcpy(fp->begin, src, remaining);
    fp->begin += remaining;

    return totalbytes;
}



///////////////////////////////////////////////////////////
//                                                       //
//            (5)  htslib-1.2.1/hts.c                    //
//                                                       //
///////////////////////////////////////////////////////////
int hts_verbose = 3;

hFILE *hopen(const char *fname, const char *mode)
{
    if (strncmp(fname, "http://", 7) == 0 ||
        strncmp(fname, "ftp://", 6) == 0) return hopen_net(fname, mode);
/////#ifdef HAVE_IRODS
/////    else if (strncmp(fname, "irods:", 6) == 0) return hopen_irods(fname, mode);
/////#endif
/////    else if (strncmp(fname, "data:", 5) == 0) return hopen_mem(fname + 5, mode);
/////    else if (strcmp(fname, "-") == 0) return hopen_fd_stdinout(mode);
    else return hopen_fd(fname, mode);
}

// Parse "x.y" text, taking care because the string is not NUL-terminated
// and filling in major/minor only when the digits are followed by a delimiter,
// so we don't misread "1.10" as "1.1" due to reaching the end of the buffer.
static void
parse_version(htsFormat *fmt, const unsigned char *u, const unsigned char *ulim)
{
    const char *str  = (const char *) u;
    const char *slim = (const char *) ulim;
    const char *s;

    fmt->version.major = fmt->version.minor = -1;

    for (s = str; s < slim; s++) if (!isdigit(*s)) break;
    if (s < slim) {
        fmt->version.major = atoi(str);
        if (*s == '.') {
            str = &s[1];
            for (s = str; s < slim; s++) if (!isdigit(*s)) break;
            if (s < slim)
                fmt->version.minor = atoi(str);
        }
        else
            fmt->version.minor = 0;
    }
}

// Decompress up to ten or so bytes by peeking at the file, which must be
// positioned at the start of a GZIP block.
static size_t decompress_peek(hFILE *fp, unsigned char *dest, size_t destsize)
{
    // Typically at most a couple of hundred bytes of input are required
    // to get a few bytes of output from inflate(), so hopefully this buffer
    // size suffices in general.
    unsigned char buffer[512];
    z_stream zs;
    ssize_t npeek = hpeek(fp, buffer, sizeof buffer);

    if (npeek < 0) return 0;

    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = buffer;
    zs.avail_in = npeek;
    zs.next_out = dest;
    zs.avail_out = destsize;
    if (inflateInit2(&zs, 31) != Z_OK) return 0;

    while (zs.total_out < destsize)
        if (inflate(&zs, Z_SYNC_FLUSH) != Z_OK) break;

    destsize = zs.total_out;
    inflateEnd(&zs);

    return destsize;
}

int hts_detect_format(hFILE *hfile, htsFormat *fmt)
{
    unsigned char s[21];
    ssize_t len = hpeek(hfile, s, 18);
    if (len < 0) return -1;

    if (len >= 2 && s[0] == 0x1f && s[1] == 0x8b) {
        // The stream is either gzip-compressed or BGZF-compressed.
        // Determine which, and decompress the first few bytes.
        fmt->compression = (len >= 18 && (s[3] & 4) &&
                            memcmp(&s[12], "BC\2\0", 4) == 0)? bgzf : gzip;
        len = decompress_peek(hfile, s, sizeof s);
    }
    else {
        fmt->compression = no_compression;
        len = hpeek(hfile, s, sizeof s);
    }
    if (len < 0) return -1;

    fmt->compression_level = -1;
    fmt->specific = NULL;

    if (len >= 6 && memcmp(s,"CRAM",4) == 0 && s[4]>=1 && s[4]<=3 && s[5]<=1) {
        fmt->category = sequence_data;
        fmt->format = cram;
        fmt->version.major = s[4], fmt->version.minor = s[5];
        fmt->compression = custom;
        return 0;
    }
    else if (len >= 4 && s[3] <= '\4') {
        if (memcmp(s, "BAM\1", 4) == 0) {
            fmt->category = sequence_data;
            fmt->format = bam;
            // TODO Decompress enough to pick version from @HD-VN header
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BAI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = bai;
            fmt->version.major = -1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BCF\4", 4) == 0) {
            fmt->category = variant_data;
            fmt->format = bcf;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BCF\2", 4) == 0) {
            fmt->category = variant_data;
            fmt->format = bcf;
            fmt->version.major = s[3];
            fmt->version.minor = (len >= 5 && s[4] <= 2)? s[4] : 0;
            return 0;
        }
        else if (memcmp(s, "CSI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = csi;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "TBI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = tbi;
            fmt->version.major = -1, fmt->version.minor = -1;
            return 0;
        }
    }
    else if (len >= 16 && memcmp(s, "##fileformat=VCF", 16) == 0) {
        fmt->category = variant_data;
        fmt->format = vcf;
        if (len >= 21 && s[16] == 'v')
            parse_version(fmt, &s[17], &s[len]);
        else
            fmt->version.major = fmt->version.minor = -1;
        return 0;
    }
    else if (len >= 4 && s[0] == '@' &&
             (memcmp(s, "@HD\t", 4) == 0 || memcmp(s, "@SQ\t", 4) == 0 ||
              memcmp(s, "@RG\t", 4) == 0 || memcmp(s, "@PG\t", 4) == 0)) {
        fmt->category = sequence_data;
        fmt->format = sam;
        // @HD-VN is not guaranteed to be the first tag, but then @HD is
        // not guaranteed to be present at all...
        if (len >= 9 && memcmp(s, "@HD\tVN:", 7) == 0)
            parse_version(fmt, &s[7], &s[len]);
        else
            fmt->version.major = 1, fmt->version.minor = -1;
        return 0;
    }
    else {
        // Various possibilities for tab-delimited text:
        // .crai   (gzipped tab-delimited six columns: seqid 5*number)
        // .bed    ([3..12] tab-delimited columns)
        // .bedpe  (>= 10 tab-delimited columns)
        // .sam    (tab-delimited >= 11 columns: seqid number seqid...)
        // FIXME For now, assume it's SAM
        fmt->category = sequence_data;
        fmt->format = sam;
        fmt->version.major = 1, fmt->version.minor = -1;
        return 0;
    }

    fmt->category = unknown_category;
    fmt->format = unknown_format;
    fmt->version.major = fmt->version.minor = -1;
    fmt->compression = no_compression;
    return 0;
}

static enum htsFormatCategory format_category(enum htsExactFormat fmt)
{
    switch (fmt) {
    case bam:
    case sam:
    case cram:
        return sequence_data;

    case vcf:
    case bcf:
        return variant_data;

    case bai:
    case crai:
    case csi:
    case gzi:
    case tbi:
        return index_file;

    case bed:
        return region_list;

    case unknown_format:
    case binary_format:
    case text_format:
    case format_maximum:
        break;
    }

    return unknown_category;
}

htsFile *hts_hopen(struct hFILE *hfile, const char *fn, const char *mode)
{

    htsFile *fp = (htsFile*)calloc(1, sizeof(htsFile));
    if (fp == NULL) goto error;

    fp->fn = strdup(fn);
    fp->is_be = ed_is_big();

    if (strchr(mode, 'r')) {
        if (hts_detect_format(hfile, &fp->format) < 0) goto error;
    }
    else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        htsFormat *fmt = &fp->format;
        fp->is_write = 1;

        if (strchr(mode, 'b')) fmt->format = binary_format;
        else if (strchr(mode, 'c')) fmt->format = cram;
        else fmt->format = text_format;

        if (strchr(mode, 'z')) fmt->compression = bgzf;
        else if (strchr(mode, 'g')) fmt->compression = gzip;
        else if (strchr(mode, 'u')) fmt->compression = no_compression;
        else {
            // No compression mode specified, set to the default for the format
            switch (fmt->format) {
            case binary_format: fmt->compression = bgzf; break;
            case cram: fmt->compression = custom; break;
            case text_format: fmt->compression = no_compression; break;
            default: abort();
            }
        }

        // Fill in category (if determinable; e.g. 'b' could be BAM or BCF)
        fmt->category = format_category(fmt->format);

        fmt->version.major = fmt->version.minor = -1;
        fmt->compression_level = -1;
        fmt->specific = NULL;
    }
    else goto error;

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        fp->fp.bgzf = bgzf_hopen(hfile, mode);
        if (fp->fp.bgzf == NULL) goto error;
        fp->is_bin = 1;
        break;

    /////case cram:
    /////    fp->fp.cram = cram_dopen(hfile, fn, mode);
    /////    if (fp->fp.cram == NULL) goto error;
    /////    if (!fp->is_write)
    /////        cram_set_option(fp->fp.cram, CRAM_OPT_DECODE_MD, 1);
    /////     fp->is_cram = 1;
    /////    break;

    case text_format:
    case sam:
    case vcf:
        if (!fp->is_write) {
        #if KS_BGZF
            BGZF *gzfp = bgzf_hopen(hfile, mode);
        #else
            // TODO Implement gzip hFILE adaptor
            hclose(hfile); // This won't work, especially for stdin
            gzFile gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
        #endif
            if (gzfp) fp->fp.voidp = ks_init(gzfp);
            else goto error;
        }
        else if (fp->format.compression != no_compression) {
            fp->fp.bgzf = bgzf_hopen(hfile, mode);
            if (fp->fp.bgzf == NULL) goto error;
        }
        else
            fp->fp.hfile = hfile;
        break;

    default:
        goto error;
    }

    return fp;

error:
    if (hts_verbose >= 2)
        fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);

    if (fp) {
        free(fp->fn);
        free(fp->fn_aux);
        free(fp);
    }
    return NULL;
}

htsFile *hts_open(const char *fn, const char *mode)
{
    htsFile *fp = NULL;
    hFILE *hfile = hopen(fn, mode);
    if (hfile == NULL) goto error;

    fp = hts_hopen(hfile, fn, mode);
//printf("1 %s\n", fp->fn);
//printf("2 %s\n", fn);


    if (fp == NULL) goto error;

    return fp;

error:
    if (hts_verbose >= 2)
        fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);

    if (hfile)
        hclose_abruptly(hfile);

    return NULL;
}

int hts_close(htsFile *fp)
{
    int ret, save;

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        ret = bgzf_close(fp->fp.bgzf);
        break;

    case cram:
/*
        if (!fp->is_write) {
            switch (cram_eof(fp->fp.cram)) {
            case 0:
                fprintf(stderr, "[E::%s] Failed to decode sequence.\n", __func__);
                return -1;
            case 2:
                fprintf(stderr, "[W::%s] EOF marker is absent. The input is probably truncated.\n", __func__);
                break;
            default: // case 1, expected EOF
                break;
            }
        }
        ret = cram_close(fp->fp.cram);
*/
        ret = -1;
        break;

    case text_format:
    case sam:
    case vcf:
        if (!fp->is_write) {
        #if KS_BGZF
            BGZF *gzfp = ((kstream_t*)fp->fp.voidp)->f;
            ret = bgzf_close(gzfp);
        #else
            gzFile gzfp = ((kstream_t*)fp->fp.voidp)->f;
            ret = gzclose(gzfp);
        #endif
            ks_destroy((kstream_t*)fp->fp.voidp);
        }
        else if (fp->format.compression != no_compression)
            ret = bgzf_close(fp->fp.bgzf);
        else
            ret = hclose(fp->fp.hfile);
        break;

    default:
        ret = -1;
        break;
    }

    save = errno;
    free(fp->fn);
    free(fp->fn_aux);
    free(fp->line.s);
    free(fp);
    errno = save;
    return ret;
}

int hts_set_fai_filename(htsFile *fp, const char *fn_aux)
{
    free(fp->fn_aux);
    if (fn_aux) {
        fp->fn_aux = strdup(fn_aux);
        if (fp->fn_aux == NULL) return -1;
    }
    else fp->fn_aux = NULL;

    ////if (fp->format.format == cram)
    ////    cram_set_option(fp->fp.cram, CRAM_OPT_REFERENCE, fp->fn_aux);

    return 0;
}

int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
{
    int ret, dret;
    ret = ks_getuntil((kstream_t*)fp->fp.voidp, delimiter, str, &dret);
    ++fp->lineno;
    return ret;
}

int hts_set_threads(htsFile *fp, int n)
{
    if (fp->format.compression == bgzf) {
        return bgzf_mt(fp->fp.bgzf, n, 256);
    } /*else if (fp->format.format == cram) {
        return hts_set_opt(fp, CRAM_OPT_NTHREADS, n);
    }*/
    else return 0;
}

///////////////////////////////////////////////////////////
//                                                       //
//             (6)   htslib-1.2.1/kstring.c              //
//                                                       //
///////////////////////////////////////////////////////////
// htslib-1.2.1/kstring.c
int kvsprintf(kstring_t *s, const char *fmt, va_list ap)
{
	va_list args;
	int l;
	va_copy(args, ap);
	l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args); // This line does not work with glibc 2.0. See `man snprintf'.
	va_end(args);
	if (l + 1 > s->m - s->l) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
		va_copy(args, ap);
		l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args);
		va_end(args);
	}
	s->l += l;
	return l;
}

// htslib-1.2.1/kstring.c
int ksprintf(kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	int l;
	va_start(ap, fmt);
	l = kvsprintf(s, fmt, ap);
	va_end(ap);
	return l;
}




///////////////////////////////////////////////////////////
//                                                       //
//            (7)  htslib-1.2.1/knetfile.c               //
//                                                       //
///////////////////////////////////////////////////////////
/* This function tests if the file handler is ready for reading (or
 * writing if is_read==0). */
static int socket_wait(int fd, int is_read)
{
	fd_set fds, *fdr = 0, *fdw = 0;
	struct timeval tv;
	int ret;
	tv.tv_sec = 5; tv.tv_usec = 0; // 5 seconds time out
	FD_ZERO(&fds);
	FD_SET(fd, &fds);
	if (is_read) fdr = &fds;
	else fdw = &fds;
	ret = select(fd+1, fdr, fdw, 0, &tv);
#ifndef _WIN32
	if (ret == -1) perror("select");
#else
	if (ret == 0)
		fprintf(stderr, "select time-out\n");
	else if (ret == SOCKET_ERROR)
		fprintf(stderr, "select: %d\n", WSAGetLastError());
#endif
	return ret;
}

/* This function does not work with Windows due to the lack of
 * getaddrinfo() in winsock. It is addapted from an example in "Beej's
 * Guide to Network Programming" (http://beej.us/guide/bgnet/). */
static int socket_connect(const char *host, const char *port)
{
#define __err_connect(func) do { perror(func); freeaddrinfo(res); return -1; } while (0)

	int ai_err, on = 1, fd;
	struct linger lng = { 0, 0 };
	struct addrinfo hints, *res = 0;
	memset(&hints, 0, sizeof(struct addrinfo));
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	/* In Unix/Mac, getaddrinfo() is the most convenient way to get
	 * server information. */
	if ((ai_err = getaddrinfo(host, port, &hints, &res)) != 0) { fprintf(stderr, "can't resolve %s:%s: %s\n", host, port, gai_strerror(ai_err)); return -1; }
	if ((fd = socket(res->ai_family, res->ai_socktype, res->ai_protocol)) == -1) __err_connect("socket");
	/* The following two setsockopt() are used by ftplib
	 * (http://nbpfaus.net/~pfau/ftplib/). I am not sure if they
	 * necessary. */
	if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on)) == -1) __err_connect("setsockopt");
	if (setsockopt(fd, SOL_SOCKET, SO_LINGER, &lng, sizeof(lng)) == -1) __err_connect("setsockopt");
	if (connect(fd, res->ai_addr, res->ai_addrlen) != 0) __err_connect("connect");
	freeaddrinfo(res);
	return fd;
}

/*************************
 * FTP specific routines *
 *************************/
static int kftp_get_response(knetFile *ftp)
{
#ifndef _WIN32
	unsigned char c;
#else
	char c;
#endif
	int n = 0;
	char *p;
	if (socket_wait(ftp->ctrl_fd, 1) <= 0) return 0;
	while (netread(ftp->ctrl_fd, &c, 1)) { // FIXME: this is *VERY BAD* for unbuffered I/O
		//fputc(c, stderr);
		if (n >= ftp->max_response) {
			ftp->max_response = ftp->max_response? ftp->max_response<<1 : 256;
			ftp->response = (char*)realloc(ftp->response, ftp->max_response);
		}
		ftp->response[n++] = c;
		if (c == '\n') {
			if (n >= 4 && isdigit(ftp->response[0]) && isdigit(ftp->response[1]) && isdigit(ftp->response[2])
				&& ftp->response[3] != '-') break;
			n = 0;
			continue;
		}
	}
	if (n < 2) return -1;
	ftp->response[n-2] = 0;
	return strtol(ftp->response, &p, 0);
}

static int kftp_send_cmd(knetFile *ftp, const char *cmd, int is_get)
{
	if (socket_wait(ftp->ctrl_fd, 0) <= 0) return -1; // socket is not ready for writing
    int len = strlen(cmd);
	if ( netwrite(ftp->ctrl_fd, cmd, len) != len ) return -1;
	return is_get? kftp_get_response(ftp) : 0;
}

static int kftp_pasv_prep(knetFile *ftp)
{
	char *p;
	int v[6];
	kftp_send_cmd(ftp, "PASV\r\n", 1);
	for (p = ftp->response; *p && *p != '('; ++p);
	if (*p != '(') return -1;
	++p;
	sscanf(p, "%d,%d,%d,%d,%d,%d", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]);
	memcpy(ftp->pasv_ip, v, 4 * sizeof(int));
	ftp->pasv_port = (v[4]<<8&0xff00) + v[5];
	return 0;
}


static int kftp_pasv_connect(knetFile *ftp)
{
	char host[80], port[10];
	if (ftp->pasv_port == 0) {
		fprintf(stderr, "[kftp_pasv_connect] kftp_pasv_prep() is not called before hand.\n");
		return -1;
	}
	sprintf(host, "%d.%d.%d.%d", ftp->pasv_ip[0], ftp->pasv_ip[1], ftp->pasv_ip[2], ftp->pasv_ip[3]);
	sprintf(port, "%d", ftp->pasv_port);
	ftp->fd = socket_connect(host, port);
	if (ftp->fd == -1) return -1;
	return 0;
}

int kftp_connect(knetFile *ftp)
{
	ftp->ctrl_fd = socket_connect(ftp->host, ftp->port);
	if (ftp->ctrl_fd == -1) return -1;
	kftp_get_response(ftp);
	kftp_send_cmd(ftp, "USER anonymous\r\n", 1);
	kftp_send_cmd(ftp, "PASS kftp@\r\n", 1);
	kftp_send_cmd(ftp, "TYPE I\r\n", 1);
	return 0;
}

int kftp_reconnect(knetFile *ftp)
{
	if (ftp->ctrl_fd != -1) {
		netclose(ftp->ctrl_fd);
		ftp->ctrl_fd = -1;
	}
	netclose(ftp->fd);
	ftp->fd = -1;
	return kftp_connect(ftp);
}

// initialize ->type, ->host, ->retr and ->size
knetFile *kftp_parse_url(const char *fn, const char *mode)
{
	knetFile *fp;
	char *p;
	int l;
	if (strstr(fn, "ftp://") != fn) return 0;
	for (p = (char*)fn + 6; *p && *p != '/'; ++p);
	if (*p != '/') return 0;
	l = p - fn - 6;
	fp = (knetFile*)calloc(1, sizeof(knetFile));
	fp->type = KNF_TYPE_FTP;
	fp->fd = -1;
	/* the Linux/Mac version of socket_connect() also recognizes a port
	 * like "ftp", but the Windows version does not. */
	fp->port = strdup("21");
	fp->host = (char*)calloc(l + 1, 1);
	if (strchr(mode, 'c')) fp->no_reconnect = 1;
	strncpy(fp->host, fn + 6, l);
	fp->retr = (char*)calloc(strlen(p) + 8, 1);
	sprintf(fp->retr, "RETR %s\r\n", p);
    fp->size_cmd = (char*)calloc(strlen(p) + 8, 1);
    sprintf(fp->size_cmd, "SIZE %s\r\n", p);
	fp->seek_offset = 0;
	return fp;
}

// place ->fd at offset off
int kftp_connect_file(knetFile *fp)
{
	int ret;
	long long file_size;
	if (fp->fd != -1) {
		netclose(fp->fd);
		if (fp->no_reconnect) kftp_get_response(fp);
	}
	kftp_pasv_prep(fp);
    kftp_send_cmd(fp, fp->size_cmd, 1);
#ifndef _WIN32
    // If the file does not exist, the response will be "550 Could not get file
    // size". Be silent on failure, hts_idx_load can be trying the existence of .csi or .tbi.
    if ( sscanf(fp->response,"%*d %lld", &file_size) != 1 ) return -1;
#else
	const char *p = fp->response;
	while (*p != ' ') ++p;
	while (*p < '0' || *p > '9') ++p;
	file_size = strtoint64(p);
#endif
	fp->file_size = file_size;
	if (fp->offset>=0) {
		char tmp[32];
#ifndef _WIN32
		sprintf(tmp, "REST %lld\r\n", (long long)fp->offset);
#else
		strcpy(tmp, "REST ");
		int64tostr(tmp + 5, fp->offset);
		strcat(tmp, "\r\n");
#endif
		kftp_send_cmd(fp, tmp, 1);
	}
	kftp_send_cmd(fp, fp->retr, 0);
	kftp_pasv_connect(fp);
	ret = kftp_get_response(fp);
	if (ret != 150) {
		fprintf(stderr, "[kftp_connect_file] %s\n", fp->response);
		netclose(fp->fd);
		fp->fd = -1;
		return -1;
	}
	fp->is_ready = 1;
	return 0;
}


// htslib-1.2.1/knetfile.c
/**************************
 * HTTP specific routines *
 **************************/
knetFile *khttp_parse_url(const char *fn, const char *mode)
{
	knetFile *fp;
	char *p, *proxy, *q;
	int l;
	if (strstr(fn, "http://") != fn) return 0;
	// set ->http_host
	for (p = (char*)fn + 7; *p && *p != '/'; ++p);
	l = p - fn - 7;
	fp = (knetFile*)calloc(1, sizeof(knetFile));
	fp->http_host = (char*)calloc(l + 1, 1);
	strncpy(fp->http_host, fn + 7, l);
	fp->http_host[l] = 0;
	for (q = fp->http_host; *q && *q != ':'; ++q);
	if (*q == ':') *q++ = 0;
	// get http_proxy
	proxy = getenv("http_proxy");
	// set ->host, ->port and ->path
	if (proxy == 0) {
		fp->host = strdup(fp->http_host); // when there is no proxy, server name is identical to http_host name.
		fp->port = strdup(*q? q : "80");
		fp->path = strdup(*p? p : "/");
	} else {
		fp->host = (strstr(proxy, "http://") == proxy)? strdup(proxy + 7) : strdup(proxy);
		for (q = fp->host; *q && *q != ':'; ++q);
		if (*q == ':') *q++ = 0; 
		fp->port = strdup(*q? q : "80");
		fp->path = strdup(fn);
	}
	fp->type = KNF_TYPE_HTTP;
	fp->ctrl_fd = fp->fd = -1;
	fp->seek_offset = 0;
	return fp;
}


//htslib-1.2.1/knetfile.c
static off_t my_netread(int fd, void *buf, off_t len)
{
	off_t rest = len, curr, l = 0;
	/* recv() and read() may not read the required length of data with
	 * one call. They have to be called repeatedly. */
	while (rest) {
		if (socket_wait(fd, 1) <= 0) break; // socket is not ready for reading
		curr = netread(fd, (void*)((char*)buf + l), rest);
		/* According to the glibc manual, section 13.2, a zero returned
		 * value indicates end-of-file (EOF), which should mean that
		 * read() will not return zero if EOF has not been met but data
		 * are not immediately available. */
		if (curr == 0) break;
		l += curr; rest -= curr;
	}
	return l;
}

// htslib-1.2.1/knetfile.c
int khttp_connect_file(knetFile *fp)
{
	int ret, l = 0;
	char *buf, *p;
	if (fp->fd != -1) netclose(fp->fd);
	fp->fd = socket_connect(fp->host, fp->port);
	buf = (char*)calloc(0x10000, 1); // FIXME: I am lazy... But in principle, 64KB should be large enough.
	l += sprintf(buf + l, "GET %s HTTP/1.0\r\nHost: %s\r\n", fp->path, fp->http_host);
    l += sprintf(buf + l, "Range: bytes=%lld-\r\n", (long long)fp->offset);
	l += sprintf(buf + l, "\r\n");
	if ( netwrite(fp->fd, buf, l) != l ) { free(buf); return -1; }
	l = 0;
	while (netread(fp->fd, buf + l, 1)) { // read HTTP header; FIXME: bad efficiency
		if (buf[l] == '\n' && l >= 3)
			if (strncmp(buf + l - 3, "\r\n\r\n", 4) == 0) break;
		++l;
	}
	buf[l] = 0;
	if (l < 14) { // prematured header
		free(buf);
		netclose(fp->fd);
		fp->fd = -1;
		return -1;
	}
	ret = strtol(buf + 8, &p, 0); // HTTP return code
	if (ret == 200 && fp->offset>0) { // 200 (complete result); then skip beginning of the file
		off_t rest = fp->offset;
		while (rest) {
			off_t l = rest < 0x10000? rest : 0x10000;
			rest -= my_netread(fp->fd, buf, l);
		}
	} else if (ret != 206 && ret != 200) {
		// failed to open file
		free(buf);
		netclose(fp->fd);
		switch (ret) {
		case 401: errno = EPERM; break;
		case 403: errno = EACCES; break;
		case 404: errno = ENOENT; break;
		case 407: errno = EPERM; break;
		case 408: errno = ETIMEDOUT; break;
		case 410: errno = ENOENT; break;
		case 503: errno = EAGAIN; break;
		case 504: errno = ETIMEDOUT; break;
		default:  errno = (ret >= 400 && ret < 500)? EINVAL : EIO; break;
		}
		fp->fd = -1;
		return -1;
	}
	free(buf);
	fp->is_ready = 1;
	return 0;
}

/********************
 * Generic routines *
 ********************/
knetFile *knet_open(const char *fn, const char *mode)
{
	knetFile *fp = 0;
	if (mode[0] != 'r') {
		fprintf(stderr, "[kftp_open] only mode \"r\" is supported.\n");
		return 0;
	}
	if (strstr(fn, "ftp://") == fn) {
		fp = kftp_parse_url(fn, mode);
		if (fp == 0) return 0;
		if (kftp_connect(fp) == -1) {
			knet_close(fp);
			return 0;
		}
		kftp_connect_file(fp);
	} else if (strstr(fn, "http://") == fn) {
		fp = khttp_parse_url(fn, mode);
		if (fp == 0) return 0;
		khttp_connect_file(fp);
	} else { // local file
#ifdef _WIN32
		/* In windows, O_BINARY is necessary. In Linux/Mac, O_BINARY may
		 * be undefined on some systems, although it is defined on my
		 * Mac and the Linux I have tested on. */
		int fd = open(fn, O_RDONLY | O_BINARY);
#else		
		int fd = open(fn, O_RDONLY);
#endif
		if (fd == -1) {
			perror("open");
			return 0;
		}
		fp = (knetFile*)calloc(1, sizeof(knetFile));
		fp->type = KNF_TYPE_LOCAL;
		fp->fd = fd;
		fp->ctrl_fd = -1;
	}
	if (fp && fp->fd == -1) {
		knet_close(fp);
		return 0;
	}
	return fp;
}

ssize_t knet_read(knetFile *fp, void *buf, size_t len)
{
	off_t l = 0;
	if (fp->fd == -1) return 0;
	if (fp->type == KNF_TYPE_FTP) {
		if (fp->is_ready == 0) {
			if (!fp->no_reconnect) kftp_reconnect(fp);
			kftp_connect_file(fp);
		}
	} else if (fp->type == KNF_TYPE_HTTP) {
		if (fp->is_ready == 0)
			khttp_connect_file(fp);
	}
	if (fp->type == KNF_TYPE_LOCAL) { // on Windows, the following block is necessary; not on UNIX
		size_t rest = len;
		ssize_t curr;
		while (rest) {
			do {
				curr = read(fp->fd, (void*)((char*)buf + l), rest);
			} while (curr < 0 && EINTR == errno);
			if (curr < 0) return -1;
			if (curr == 0) break;
			l += curr; rest -= curr;
		}
	} else l = my_netread(fp->fd, buf, len);
	fp->offset += l;
	return l;
}

off_t knet_seek(knetFile *fp, off_t off, int whence)
{
	if (whence == SEEK_SET && off == fp->offset) return 0;
	if (fp->type == KNF_TYPE_LOCAL) {
		/* Be aware that lseek() returns the offset after seeking, while fseek() returns zero on success. */
		off_t offset = lseek(fp->fd, off, whence);
		if (offset == -1) return -1;
		fp->offset = offset;
		return fp->offset;
	} else if (fp->type == KNF_TYPE_FTP) {
		if (whence == SEEK_CUR) fp->offset += off;
		else if (whence == SEEK_SET) fp->offset = off;
		else if (whence == SEEK_END) fp->offset = fp->file_size + off;
		else return -1;
		fp->is_ready = 0;
		return fp->offset;
	} else if (fp->type == KNF_TYPE_HTTP) {
		if (whence == SEEK_END) { // FIXME: can we allow SEEK_END in future?
			fprintf(stderr, "[knet_seek] SEEK_END is not supported for HTTP. Offset is unchanged.\n");
			errno = ESPIPE;
			return -1;
		}
		if (whence == SEEK_CUR) fp->offset += off;
		else if (whence == SEEK_SET) fp->offset = off;
		else return -1;
		fp->is_ready = 0;
		return fp->offset;
	}
	errno = EINVAL;
	fprintf(stderr,"[knet_seek] %s\n", strerror(errno));
	return -1;
}

int knet_close(knetFile *fp)
{
	if (fp == 0) return 0;
	if (fp->ctrl_fd != -1) netclose(fp->ctrl_fd); // FTP specific
	if (fp->fd != -1) {
		/* On Linux/Mac, netclose() is an alias of close(), but on
		 * Windows, it is an alias of closesocket(). */
		if (fp->type == KNF_TYPE_LOCAL) close(fp->fd);
		else netclose(fp->fd);
	}
	free(fp->host); free(fp->port);
	free(fp->response); free(fp->retr); // FTP specific
	free(fp->path); free(fp->http_host); // HTTP specific
	free(fp);
	return 0;
}


///////////////////////////////////////////////////////////
//                                                       //
//            (8) htslib-1.2.1/hfile_net.c               //
//                                                       //
///////////////////////////////////////////////////////////
typedef struct {
    hFILE base;
    knetFile *netfp;
} hFILE_net;

static int net_inited = 0;

static int net_init(void)
{
#ifdef _WIN32
    if (knet_win32_init() != 0) return -1;

    // In the unlikely event atexit() fails, it's better to succeed here and
    // carry on and do the I/O; then eventually when the program exits, we'll
    // merely have failed to clean up properly, as if we had aborted.
    (void) atexit(net_exit);
#endif

    net_inited = 1;
    return 0;
}

static ssize_t net_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_net *fp = (hFILE_net *) fpv;
    return knet_read(fp->netfp, buffer, nbytes);
}

static off_t net_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_net *fp = (hFILE_net *) fpv;
    return knet_seek(fp->netfp, offset, whence);
}

static int net_close(hFILE *fpv)
{
    hFILE_net *fp = (hFILE_net *) fpv;
    return knet_close(fp->netfp);
}

static const struct hFILE_backend net_backend =
{
    net_read, NULL, net_seek, NULL, net_close
};

hFILE *hopen_net(const char *filename, const char *mode)
{
    hFILE_net *fp;

    // Do any networking initialisation if this is the first use.
    if (! net_inited) { if (net_init() < 0) return NULL; }

    fp = (hFILE_net *) hfile_init(sizeof (hFILE_net), mode, 0);
    if (fp == NULL) return NULL;

    fp->netfp = knet_open(filename, mode);
    if (fp->netfp == NULL) { hfile_destroy((hFILE *) fp); return NULL; }

    fp->base.backend = &net_backend;
    return &fp->base;
}




























///////////////////////////////////////////////////////////
//                                                       //
//             (9)  htslib-1.2.1/sam.c                   //
//                                                       //
///////////////////////////////////////////////////////////
typedef khash_t(s2i) sdict_t;

// htslib-1.2.1/sam.c
static bam_hdr_t *hdr_from_dict(sdict_t *d)
{
    bam_hdr_t *h;
    khint_t k;
    h = bam_hdr_init();
    h->sdict = d;
    h->n_targets = kh_size(d);
    h->target_len = (uint32_t*)malloc(sizeof(uint32_t) * h->n_targets);
    h->target_name = (char**)malloc(sizeof(char*) * h->n_targets);
    for (k = kh_begin(d); k != kh_end(d); ++k) {
        if (!kh_exist(d, k)) continue;
        h->target_name[kh_val(d, k)>>32] = (char*)kh_key(d, k);
        h->target_len[kh_val(d, k)>>32]  = kh_val(d, k)<<32>>32;
        kh_val(d, k) >>= 32;
    }
    return h;
}

// htslib-1.2.1/sam.c
bam_hdr_t *sam_hdr_read(htsFile *fp)
{
    switch (fp->format.format) {
    case bam:
        return bam_hdr_read(fp->fp.bgzf);

    case cram:
        abort();
        //return cram_header_to_bam(fp->fp.cram->header);

    case sam: {
        kstring_t str;
        bam_hdr_t *h;
        int has_SQ = 0;
        str.l = str.m = 0; str.s = 0;
        while (hts_getline(fp, KS_SEP_LINE, &fp->line) >= 0) {
            if (fp->line.s[0] != '@') break;
            if (fp->line.l > 3 && strncmp(fp->line.s,"@SQ",3) == 0) has_SQ = 1;
            kputsn(fp->line.s, fp->line.l, &str);
            kputc('\n', &str);
        }
        if (! has_SQ && fp->fn_aux) {
            char line[2048];
            FILE *f = fopen(fp->fn_aux, "r");
            if (f == NULL) return NULL;
            while (fgets(line, sizeof line, f)) {
                const char *name = strtok(line, "\t");
                const char *length = strtok(NULL, "\t");
                ksprintf(&str, "@SQ\tSN:%s\tLN:%s\n", name, length);
            }
            fclose(f);
        }
        if (str.l == 0) kputsn("", 0, &str);
        h = sam_hdr_parse(str.l, str.s);
        h->l_text = str.l; h->text = str.s;
        return h;
        }

    default:
        abort();
    }
}


//////////////
// HEADER
//////////////
// int main(int argc, char *argv[])         [(SAMToolsRootDir)/bamtk.c]            ->
// int main_samview(int argc, char *argv[]) [(SAMToolsRootDir)/sam_view.c]         ->
// bam_hdr_t *sam_hdr_read(htsFile *fp)     [(SAMToolsRootDir)/htslib-1.2.1/sam.c] ->
// bam_hdr_t *bam_hdr_read(BGZF *fp)        [(SAMToolsRootDir)/htslib-1.2.1/sam.c] ->

// htslib-1.2.1/sam.c
bam_hdr_t *bam_hdr_read(BGZF *fp)
{
//printf("*** Entering function bam_hdr_read() in sam.c \n");
    bam_hdr_t *h;
    char buf[4];
    int magic_len, has_EOF;
    int32_t i = 1, name_len;
    // check EOF
    has_EOF = bgzf_check_EOF(fp);
    if (has_EOF < 0) {
        perror("[W::sam_hdr_read] bgzf_check_EOF");
    } else if (has_EOF == 0 && hts_verbose >= 2)
        fprintf(stderr, "[W::%s] EOF marker is absent. The input is probably truncated.\n", __func__);
    // read "BAM1"

    magic_len = bgzf_read(fp, buf, 4);
    if (magic_len != 4 || strncmp(buf, "BAM\1", 4)) {
        if (hts_verbose >= 1) fprintf(stderr, "[E::%s] invalid BAM binary header\n", __func__);
        return 0;
    }
    h = bam_hdr_init();
    // read plain text and the number of reference sequences
    bgzf_read(fp, &h->l_text, 4);
    if (fp->is_be) ed_swap_4p(&h->l_text);
    h->text = (char*)malloc(h->l_text + 1);
    h->text[h->l_text] = 0; // make sure it is NULL terminated
//printf("h->l_text: %d\n", h->l_text);
    bgzf_read(fp, h->text, h->l_text);
    bgzf_read(fp, &h->n_targets, 4);
    if (fp->is_be) ed_swap_4p(&h->n_targets);
    // read reference sequence names and lengths
    h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
    h->target_len = (uint32_t*)calloc(h->n_targets, sizeof(uint32_t));
//printf("h->n_targets: %d\n", h->n_targets);
    for (i = 0; i != h->n_targets; ++i) {
        bgzf_read(fp, &name_len, 4);
        if (fp->is_be) ed_swap_4p(&name_len);
        h->target_name[i] = (char*)calloc(name_len, 1);
        bgzf_read(fp, h->target_name[i], name_len);
        bgzf_read(fp, &h->target_len[i], 4);
//printf("TARGET_NAME: h->target_name[%d] = %s\n", i, h->target_name[i]);
        if (fp->is_be) ed_swap_4p(&h->target_len[i]);
    }
    return h;
}

// htslib-1.2.1/sam.c
int bam_hdr_write(BGZF *fp, const bam_hdr_t *h)
{
    char buf[4];
    int32_t i, name_len, x;
    // write "BAM1"
    strncpy(buf, "BAM\1", 4);
    bgzf_write(fp, buf, 4);
    // write plain text and the number of reference sequences
    if (fp->is_be) {
        x = ed_swap_4(h->l_text);
        bgzf_write(fp, &x, 4);
        if (h->l_text) bgzf_write(fp, h->text, h->l_text);
        x = ed_swap_4(h->n_targets);
        bgzf_write(fp, &x, 4);
    } else {
        bgzf_write(fp, &h->l_text, 4);
        if (h->l_text) bgzf_write(fp, h->text, h->l_text);
        bgzf_write(fp, &h->n_targets, 4);
    }
    // write sequence names and lengths
    for (i = 0; i != h->n_targets; ++i) {
        char *p = h->target_name[i];
        name_len = strlen(p) + 1;
        if (fp->is_be) {
            x = ed_swap_4(name_len);
            bgzf_write(fp, &x, 4);
        } else bgzf_write(fp, &name_len, 4);
        bgzf_write(fp, p, name_len);
        if (fp->is_be) {
            x = ed_swap_4(h->target_len[i]);
            bgzf_write(fp, &x, 4);
        } else bgzf_write(fp, &h->target_len[i], 4);
    }
    bgzf_flush(fp);
    return 0;
}


// htslib-1.2.1/sam.c
bam_hdr_t *bam_hdr_init()
{
    return (bam_hdr_t*)calloc(1, sizeof(bam_hdr_t));
}

// htslib-1.2.1/sam.c
static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

// htslib-1.2.1/sam.c
static void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
    uint8_t *s;
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i, n;
    s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
    while (s < data + l_data) {
        int size;
        s += 2; // skip key
        size = aux_type2size(*s); ++s; // skip type
        switch (size) {
        case 1: ++s; break;
        case 2: ed_swap_2p(s); s += 2; break;
        case 4: ed_swap_4p(s); s += 4; break;
        case 8: ed_swap_8p(s); s += 8; break;
        case 'Z':
        case 'H':
            while (*s) ++s;
            ++s;
            break;
        case 'B':
            size = aux_type2size(*s); ++s;
            if (is_host) memcpy(&n, s, 4), ed_swap_4p(s);
            else ed_swap_4p(s), memcpy(&n, s, 4);
            s += 4;
            switch (size) {
            case 1: s += n; break;
            case 2: for (i = 0; i < n; ++i, s += 2) ed_swap_2p(s); break;
            case 4: for (i = 0; i < n; ++i, s += 4) ed_swap_4p(s); break;
            case 8: for (i = 0; i < n; ++i, s += 8) ed_swap_8p(s); break;
            }
            break;
        }
    }
}


//////////////
// CORE
//////////////

// htslib-1.2.1/sam.c
int bam_read1(BGZF *fp, bam1_t *b, int *tid, int *pos, int *match, int *len)
{
//printf ("*** Entering core function bam_read1() in sam.c \n");
    bam1_core_t *c = &b->core;
    int32_t block_len, ret, i;
    uint32_t x[8];
    if ((ret = bgzf_read(fp, &block_len, 4)) != 4) {
        if (ret == 0) return -1; // normal end-of-file
        else return -2; // truncated
    }
//printf("block length: %d\n", block_len);
    if (bgzf_read(fp, x, 32) != 32) return -3;
    if (fp->is_be) {

        ed_swap_4p(&block_len);

        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
    }
    c->tid = x[0]; c->pos = x[1];

//printf("c->tid: %d\n", c->tid);
    *tid = c->tid;
    *pos = c->pos+1;

    
    c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
    c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
    b->l_data = block_len - 32;
    if (b->l_data < 0 || c->l_qseq < 0) return -4;
    if ((char *)bam_get_aux(b) - (char *)b->data > b->l_data)
        return -4;
    if (b->m_data < b->l_data) {
        b->m_data = b->l_data;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
        if (!b->data)
            return -4;
    }
//printf("l_data: %d\n", b->l_data); 

    if (bgzf_read(fp, b->data, b->l_data) != b->l_data) return -4;
    //b->l_aux = b->l_data - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);

	  if (c->n_cigar) // cigar
	  {
      uint32_t *cigar = bam_get_cigar(b);
		  for (i = 0; i < c->n_cigar; ++i)
		  {
        int oplength = bam_cigar_oplen(cigar[i]);
        char opcharacter = bam_cigar_opchr(cigar[i]);
			  if (opcharacter == 'M' || 
            opcharacter == 'I' ||
            opcharacter == 'S' ||
            opcharacter == 'X' ||
            opcharacter == '=')
			  {
				  *match = 1;
				  *len += oplength;
			  }
      }
    }

    return 4 + block_len;
}

// htslib-1.2.1/sam.c
bam1_t *bam_init1()
{
//printf("bam_init: %d", sizeof(bam1_t));
    return (bam1_t*)calloc(1, sizeof(bam1_t));
}

// htslib-1.2.1/sam.c
void bam_destroy1(bam1_t *b)
{
    if (b == 0) return;
    free(b->data); free(b);
}

// htslib-1.2.1/sam.c
int sam_read1(htsFile *fp, bam_hdr_t *h, bam1_t *b, int *tid, int *pos, int *match, int *len)
{
    switch (fp->format.format) {
    case bam: {
//printf("** Entering function sam_read1() in sam.c - switch: bam\n");
        int r = bam_read1(fp->fp.bgzf, b, tid, pos, match, len);
//printf("%d\n", r);
        if (r >= 0) {
            if (b->core.tid  >= h->n_targets || b->core.tid  < -1 ||
                b->core.mtid >= h->n_targets || b->core.mtid < -1)
                return -3;
        }
//printf("return %d\n", r);
        return r;
        }

    case cram:
        /////return cram_get_bam_seq(fp->fp.cram, &b);
        return -1;

    case sam: {
//printf("\nsam_read1: sam\n");
        int ret;
err_recover:
        if (fp->line.l == 0) {
            ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
            if (ret < 0) return -1;
        }
        ret = sam_parse1(&fp->line, h, b);
        fp->line.l = 0;
        if (ret < 0) {
            if (hts_verbose >= 1)
                fprintf(stderr, "[W::%s] parse error at line %lld\n", __func__, (long long)fp->lineno);
            if (h->ignore_sam_err) goto err_recover;
        }
        return ret;
        }

    default:
        abort();
    }
}

// htslib-1.2.1/sam.c
/**********************
 *** SAM record I/O ***
 **********************/
int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b)
{
#define _read_token(_p) (_p); for (; *(_p) && *(_p) != '\t'; ++(_p)); if (*(_p) != '\t') goto err_ret; *(_p)++ = 0
#define _read_token_aux(_p) (_p); for (; *(_p) && *(_p) != '\t'; ++(_p)); *(_p)++ = 0 // this is different in that it does not test *(_p)=='\t'
#define _get_mem(type_t, _x, _s, _l) ks_resize((_s), (_s)->l + (_l)); *(_x) = (type_t*)((_s)->s + (_s)->l); (_s)->l += (_l)
#define _parse_err(cond, msg) do { if ((cond) && hts_verbose >= 1) { fprintf(stderr, "[E::%s] " msg "\n", __func__); goto err_ret; } } while (0)
#define _parse_warn(cond, msg) if ((cond) && hts_verbose >= 2) fprintf(stderr, "[W::%s] " msg "\n", __func__)

    uint8_t *t;
    char *p = s->s, *q;
    int i;
    kstring_t str;
    bam1_core_t *c = &b->core;

    str.l = b->l_data = 0;
    str.s = (char*)b->data; str.m = b->m_data;
    memset(c, 0, 32);
    if (h->cigar_tab == 0) {
        h->cigar_tab = (int8_t*) malloc(128);
        for (i = 0; i < 128; ++i)
            h->cigar_tab[i] = -1;
        for (i = 0; BAM_CIGAR_STR[i]; ++i)
            h->cigar_tab[(int)BAM_CIGAR_STR[i]] = i;
    }
    // qname
    q = _read_token(p);
    kputsn_(q, p - q, &str);
    c->l_qname = p - q;
    // flag
    c->flag = strtol(p, &p, 0);
    if (*p++ != '\t') goto err_ret; // malformated flag
    // chr
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(h->n_targets == 0, "missing SAM header");
        /////c->tid = bam_name2id(h, q);
        _parse_warn(c->tid < 0, "urecognized reference name; treated as unmapped");
    } else c->tid = -1;
    // pos
    c->pos = strtol(p, &p, 10) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->pos < 0 && c->tid >= 0) {
        _parse_warn(1, "mapped query cannot have zero coordinate; treated as unmapped");
        c->tid = -1;
    }
    if (c->tid < 0) c->flag |= BAM_FUNMAP;
    // mapq
    c->qual = strtol(p, &p, 10);
    if (*p++ != '\t') goto err_ret;
    // cigar
    if (*p != '*') {
        uint32_t *cigar;
        size_t n_cigar = 0;
        for (q = p; *p && *p != '\t'; ++p)
            if (!isdigit(*p)) ++n_cigar;
        if (*p++ != '\t') goto err_ret;
        _parse_err(n_cigar >= 65536, "too many CIGAR operations");
        c->n_cigar = n_cigar;
        _get_mem(uint32_t, &cigar, &str, c->n_cigar<<2);
        for (i = 0; i < c->n_cigar; ++i, ++q) {
            int op;
            cigar[i] = strtol(q, &q, 10)<<BAM_CIGAR_SHIFT;
            op = (uint8_t)*q >= 128? -1 : h->cigar_tab[(int)*q];
            _parse_err(op < 0, "unrecognized CIGAR operator");
            cigar[i] |= op;
        }
        /////i = bam_cigar2rlen(c->n_cigar, cigar);
    } else {
        _parse_warn(!(c->flag&BAM_FUNMAP), "mapped query must have a CIGAR; treated as unmapped");
        c->flag |= BAM_FUNMAP;
        q = _read_token(p);
        i = 1;
    }
    c->bin = hts_reg2bin(c->pos, c->pos + i, 14, 5);
    // mate chr
    q = _read_token(p);
    if (strcmp(q, "=") == 0) c->mtid = c->tid;
    else if (strcmp(q, "*") == 0) c->mtid = -1;
    /////else c->mtid = bam_name2id(h, q);
    // mpos
    c->mpos = strtol(p, &p, 10) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->mpos < 0 && c->mtid >= 0) {
        _parse_warn(1, "mapped mate cannot have zero coordinate; treated as unmapped");
        c->mtid = -1;
    }
    // tlen
    c->isize = strtol(p, &p, 10);
    if (*p++ != '\t') goto err_ret;
    // seq
    q = _read_token(p);
    if (strcmp(q, "*")) {
        c->l_qseq = p - q - 1;
        /////i = bam_cigar2qlen(c->n_cigar, (uint32_t*)(str.s + c->l_qname));
        _parse_err(c->n_cigar && i != c->l_qseq, "CIGAR and query sequence are of different length");
        i = (c->l_qseq + 1) >> 1;
        _get_mem(uint8_t, &t, &str, i);
        memset(t, 0, i);
        /////for (i = 0; i < c->l_qseq; ++i)
        /////    t[i>>1] |= seq_nt16_table[(int)q[i]] << ((~i&1)<<2);
    } else c->l_qseq = 0;
    // qual
    q = _read_token_aux(p);
    _get_mem(uint8_t, &t, &str, c->l_qseq);
    if (strcmp(q, "*")) {
        _parse_err(p - q - 1 != c->l_qseq, "SEQ and QUAL are of different length");
        for (i = 0; i < c->l_qseq; ++i) t[i] = q[i] - 33;
    } else memset(t, 0xff, c->l_qseq);
    // aux
    // Note that (like the bam1_core_t fields) this aux data in b->data is
    // stored in host endianness; so there is no byte swapping needed here.
    while (p < s->s + s->l) {
        uint8_t type;
        q = _read_token_aux(p); // FIXME: can be accelerated for long 'B' arrays
        _parse_err(p - q - 1 < 6, "incomplete aux field");
        kputsn_(q, 2, &str);
        q += 3; type = *q++; ++q; // q points to value
        if (type == 'A' || type == 'a' || type == 'c' || type == 'C') {
            kputc_('A', &str);
            kputc_(*q, &str);
        } else if (type == 'i' || type == 'I') {
            if (*q == '-') {
                long x = strtol(q, &q, 10);
                if (x >= INT8_MIN) {
                    kputc_('c', &str); kputc_(x, &str);
                } else if (x >= INT16_MIN) {
                    int16_t y = x;
                    kputc_('s', &str); kputsn_((char*)&y, 2, &str);
                } else {
                    int32_t y = x;
                    kputc_('i', &str); kputsn_(&y, 4, &str);
                }
            } else {
                unsigned long x = strtoul(q, &q, 10);
                if (x <= UINT8_MAX) {
                    kputc_('C', &str); kputc_(x, &str);
                } else if (x <= UINT16_MAX) {
                    uint16_t y = x;
                    kputc_('S', &str); kputsn_(&y, 2, &str);
                } else {
                    uint32_t y = x;
                    kputc_('I', &str); kputsn_(&y, 4, &str);
                }
            }
        } else if (type == 'f') {
            float x;
            x = strtod(q, &q);
            kputc_('f', &str); kputsn_(&x, 4, &str);
        } else if (type == 'd') {
            double x;
            x = strtod(q, &q);
            kputc_('d', &str); kputsn_(&x, 8, &str);
        } else if (type == 'Z' || type == 'H') {
            kputc_(type, &str);kputsn_(q, p - q, &str); // note that this include the trailing NULL
        } else if (type == 'B') {
            int32_t n;
            char *r;
            _parse_err(p - q - 1 < 3, "incomplete B-typed aux field");
            type = *q++; // q points to the first ',' following the typing byte
            for (r = q, n = 0; *r; ++r)
                if (*r == ',') ++n;
            kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
            // FIXME: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
            if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0); kputc_(x, &str); }
            else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtoul(q + 1, &q, 0); kputc_(x, &str); }
            else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
            else if (type == 'S') while (q + 1 < p) { uint16_t x = strtoul(q + 1, &q, 0); kputsn_(&x, 2, &str); }
            else if (type == 'i') while (q + 1 < p) { int32_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
            else if (type == 'I') while (q + 1 < p) { uint32_t x = strtoul(q + 1, &q, 0); kputsn_(&x, 4, &str); }
            else if (type == 'f') while (q + 1 < p) { float    x = strtod(q + 1, &q);    kputsn_(&x, 4, &str); }
            else _parse_err(1, "unrecognized type");
        } else _parse_err(1, "unrecognized type");
    }
    b->data = (uint8_t*)str.s; b->l_data = str.l; b->m_data = str.m;
    return 0;

#undef _parse_warn
#undef _parse_err
#undef _get_mem
#undef _read_token_aux
#undef _read_token
err_ret:
    b->data = (uint8_t*)str.s; b->l_data = str.l; b->m_data = str.m;
    return -2;
}

// htslib-1.2.1/sam.c
/* Called only from hputs(), when our buffer is already full.  */
int hputs2(const char *text, size_t totalbytes, size_t ncopied, hFILE *fp)
{
    return (hwrite2(fp, text, totalbytes, ncopied) >= 0)? 0 : EOF;
}

// htslib-1.2.1/sam.c
int sam_hdr_write(htsFile *fp, const bam_hdr_t *h)
{
    switch (fp->format.format) {
    case binary_format:
        fp->format.category = sequence_data;
        fp->format.format = bam;
        /* fall-through */
    case bam:
        bam_hdr_write(fp->fp.bgzf, h);
        break;

    case cram: 
        /////{
        /////cram_fd *fd = fp->fp.cram;
        /////if (cram_set_header(fd, bam_header_to_cram((bam_hdr_t *)h)) < 0) return -1;
        /////if (fp->fn_aux)
        /////    cram_load_reference(fd, fp->fn_aux);
        /////if (cram_write_SAM_hdr(fd, fd->header) < 0) return -1;
        /////}
        abort();
        break;

    case text_format:
        fp->format.category = sequence_data;
        fp->format.format = sam;
        /* fall-through */
    case sam: {
        char *p;
        hputs(h->text, fp->fp.hfile);
        p = strstr(h->text, "@SQ\t"); // FIXME: we need a loop to make sure "@SQ\t" does not match something unwanted!!!
        if (p == 0) {
            int i;
            for (i = 0; i < h->n_targets; ++i) {
                fp->line.l = 0;
                kputsn("@SQ\tSN:", 7, &fp->line); kputs(h->target_name[i], &fp->line);
                kputsn("\tLN:", 4, &fp->line); kputw(h->target_len[i], &fp->line); kputc('\n', &fp->line);
                if ( hwrite(fp->fp.hfile, fp->line.s, fp->line.l) != fp->line.l ) return -1;
            }
        }
        if ( hflush(fp->fp.hfile) != 0 ) return -1;
        }
        break;

    default:
        abort();
    }
    return 0;
}

// htslib-1.2.1/sam.c
bam_hdr_t *sam_hdr_parse(int l_text, const char *text)
{
    const char *q, *r, *p;
    khash_t(s2i) *d;
    d = kh_init(s2i);
    for (p = text; *p; ++p) {
        if (strncmp(p, "@SQ", 3) == 0) {
            char *sn = 0;
            int ln = -1;
            for (q = p + 4;; ++q) {
                if (strncmp(q, "SN:", 3) == 0) {
                    q += 3;
                    for (r = q; *r != '\t' && *r != '\n'; ++r);
                    sn = (char*)calloc(r - q + 1, 1);
                    strncpy(sn, q, r - q);
                    q = r;
                } else if (strncmp(q, "LN:", 3) == 0)
                    ln = strtol(q + 3, (char**)&q, 10);
                while (*q != '\t' && *q != '\n') ++q;
                if (*q == '\n') break;
            }
            p = q;
            if (sn && ln >= 0) {
                khint_t k;
                int absent;
                k = kh_put(s2i, d, sn, &absent);
                if (!absent) {
                    if (hts_verbose >= 2)
                        fprintf(stderr, "[W::%s] duplicated sequence '%s'\n", __func__, sn);
                    free(sn);
                } else kh_val(d, k) = (int64_t)(kh_size(d) - 1)<<32 | ln;
            }
        }
        while (*p != '\n') ++p;
    }
    return hdr_from_dict(d);
}

