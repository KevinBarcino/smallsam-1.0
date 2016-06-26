#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "./htslib-1.2.1/htslib/sam.h"
#include "./htslib-1.2.1/htslib/hts.h"
//#include "./htslib-1.2.1/htslib/faidx.h"
//#include "./htslib-1.2.1/htslib/kstring.h"
#include "./htslib-1.2.1/htslib/khash.h"
#include "samtools.h"
#include "kccommon.h"

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

FILE *pFile2 = NULL;
int main(int argc, char *argv[])
{
  /************************************************************/  
  FILE *pFile1 = NULL;
  char file1[1024] = {0};
  char file2[1024] = {0};
  char filenameBuf[1024] = {0};
  char *pos1 = NULL;
  /************************************************************/

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
/*
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
      case 'L':
          if ((settings.bed = bed_read(optarg)) == NULL) {
              print_error_errno("Could not read file \"%s\"", optarg);
              ret = 1;
              goto view_end;
          }
          break;
      case 'r':
          if (add_read_group_single(&settings, optarg) != 0) {
              ret = 1;
              goto view_end;
          }
          break;
      case 'R':
          if (add_read_groups_file(&settings, optarg) != 0) {
              ret = 1;
              goto view_end;
          }
          break;
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
*/
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
  if (fn_list == 0 && fn_ref) fn_list = samfaipath(fn_ref);
strcpy(file1, argv[optind]);
/****************************************/
strcpy(filenameBuf, argv[optind]);
pos1 = strrchr(filenameBuf, '.');
if (pos1 != NULL)
{
*pos1 = '\0';
}	
strcpy(file2, filenameBuf);
strcat(file2, ".csv");
//printf("%s\n", file1);
//printf("%s\n", file2);
pFile1 = fopen(file1, "r");
pFile2 = fopen(file2, "w");
/*****************************************/
  // open file handlers
  if ((in = sam_open(argv[optind], "r")) == 0) {

      print_error_errno("failed to open \"%s\" for reading", argv[optind]);
      ret = 1;
      goto view_end;
  }
  if (fn_list) hts_set_fai_filename(in, fn_list);
  if ((header = sam_hdr_read(in)) == 0) {
      fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", argv[optind]);
      ret = 1;
      goto view_end;
  }

  if (settings.rghash) { // FIXME: I do not know what "bam_header_t::n_text" is for...
      char *tmp;
      int l;
      tmp = drop_rg(header->text, settings.rghash, &l);
      free(header->text);
      header->text = tmp;
      header->l_text = l;
  }
  if (!is_count) {
      if ((out = sam_open(fn_out? fn_out : "-", out_mode)) == 0) {
          print_error_errno("failed to open \"%s\" for writing", fn_out? fn_out : "standard output");
          ret = 1;
          goto view_end;
      }
      if (fn_list) hts_set_fai_filename(out, fn_list);
      if (*out_format || is_header)  {
          if (sam_hdr_write(out, header) != 0) {
              fprintf(stderr, "[main_samview] failed to write the SAM header\n");
              ret = 1;
              goto view_end;
          }
      }
      if (fn_un_out) {
          if ((un_out = sam_open(fn_un_out, out_mode)) == 0) {
              print_error_errno("failed to open \"%s\" for writing", fn_un_out);
              ret = 1;
              goto view_end;
          }
          if (*out_format || is_header) {
              if (sam_hdr_write(un_out, header) != 0) {
                  fprintf(stderr, "[main_samview] failed to write the SAM header\n");
                  ret = 1;
                  goto view_end;
              }
          }
      }
  }
  if (n_threads > 1) { if (out) hts_set_threads(out, n_threads); }
  if (is_header_only) goto view_end; // no need to print alignments

//    int startCount = 0;

//printf("\nhere 1\n");
  if (argc == optind + 1) { // convert/print the entire file
//printf("\nhere 2\n");
      bam1_t *b = bam_init1();
      int r;
      while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
          if (!process_aln(header, b, &settings))
          {
              if (!is_count)
              {   /* printf("\nhere 3b\n"); */
                  
                  if (check_sam_write1(out, header, b, fn_out, &ret) < 0)
                      break;
                  //startCount++;
                  //printf("%d", startCount);
                  
              }
              count++;
          }
          else 
          {
              if (un_out) { /* printf("\nhere 3c\n"); */ if (check_sam_write1(un_out, header, b, fn_un_out, &ret) < 0) break; }
          }
      }
      if (r < -1) {
          fprintf(stderr, "[main_samview] truncated file.\n");
          ret = 1;
      }
      bam_destroy1(b);
  } else { // retrieve alignments in specified regions
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
              if (!process_aln(header, b, &settings)) {
                  if (!is_count) { if (check_sam_write1(out, header, b, fn_out, &ret) < 0) break; }
                  count++;
              } else {
                  if (un_out) { if (check_sam_write1(un_out, header, b, fn_un_out, &ret) < 0) break; }
              }
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

      fclose(pFile1);
      fclose(pFile2);
  }

view_end:
  if (is_count && ret == 0)
      printf("%" PRId64 "\n", count);

  // close files, free and return
  if (in) check_sam_close(in, argv[optind], "standard input", &ret);
  if (out) check_sam_close(out, fn_out, "standard output", &ret);
  if (un_out) check_sam_close(un_out, fn_un_out, "file", &ret);

  free(fn_list); free(fn_ref); free(fn_out); free(settings.library);  free(fn_un_out);
  if ( header ) bam_hdr_destroy(header);
  if (settings.bed) bed_destroy(settings.bed);
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
