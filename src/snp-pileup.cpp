// #include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>
#include "snp-pileup-rev.h"
#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib> // For exit()
#include <cstdio>  // For printf, FILE
#include <sstream>
#include <ctime>

// Define the arguments structure
#ifndef SNP_PILEUP_REV_H_DEFINED
#define SNP_PILEUP_REV_H_DEFINED
struct arguments 
{
  std::vector<std::string> args;
  bool count_orphans = true;
  bool ignore_overlaps = false;
  int min_base_quality = 0;
  int min_map_quality = 0;
  std::vector<int> min_read_counts;
  int max_depth = 4000;
  int rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  int pseudo_snps = 0;
  BGZF *gzippedPointer = nullptr; // Added
};
#endif


// Writes gzipped output to a file 
void gzip_output(arguments arguments, std::string str, BGZF *fp)
{
  if (bgzf_write(arguments.gzippedPointer, str.c_str(), str.length()) < 0)
  {
    Rprintf("Error: Failed to write to bgzf file.");
  }
}

// Checks if a string ends with a certain substring
inline bool ends_with(std::string const &value, std::string const &ending)
{
  if (ending.size() > value.size())
    return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

static int mplp_func(void *data, bam1_t *b)
{
  // it seems to me that this function is run once every read, and lets you control if the read should be skipped or not
  // here is where you would put quality checks and things like that
  //	printf("hello\n");
  // char *ref;
  mplp_aux_t *ma = (mplp_aux_t *)data;
  int ret, skip = 0; // ref_len;
  do
  {
    // int has_ref;
    ret = sam_read1(ma->fp, ma->h, b);
    if (ret < 0)
      break;

    /*if (b->core.flag & READ_FAILED_QUALITY_CHECKS || b->core.flag & READ_SECONDARY_ALIGNMENT || b->core.flag & READ_PCR_OPTICAL_DUPLICATE || b->core.flag & READ_UNMAPPED) {
        skip = 1;
        continue;
    }*/

    if (ma->conf->rflag_filter && ma->conf->rflag_filter & b->core.flag)
    {
      skip = 1;
      continue;
    }

    if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP))
    { // exclude unmapped reads
      skip = 1;
      continue;
    }

    /*if (ma->conf->bed && ma->conf->all == 0) { // test overlap
        skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
        if (skip) continue;
    }*/

    skip = 0;
    /*if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
    if (has_ref && ma->conf->capQ_thres > 10) {
        int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
        if (q < 0) skip = 1;
        else if (b->core.qual > q) b->core.qual = q;
        } */
    if (b->core.qual < ma->conf->min_map_quality)
    {
      skip = 1;
    }
    else if (b->core.qual == 0)
    {
      skip = 1;
    }
    else if (!ma->conf->count_orphans && (b->core.flag & BAM_FPAIRED) && !(b->core.flag & BAM_FPROPER_PAIR))
    {
      skip = 1;
    }
  } while (skip);
  return ret;
}

static int vcf_chr_to_bam(char *chromosome, char **bam_chrs, int32_t n_targets)
{
  // try to remove chr prefix if it exists
  // Rcpp::Rcout << "Debug: remove chr prefix if exists" << std::endl;
  if (!strncmp(chromosome, "chr", 3))
  {
    chromosome += 3;
  }
  bool bam_chr_prefix = (!strncmp(bam_chrs[0], "chr", 3));
  for (int i = 0; i < n_targets; i++)
  {
    if (!bam_chr_prefix)
    {
      if (!strcmp(bam_chrs[i], chromosome))
      {
        return i;
      }
    }
    else
    {
      if (!strcmp(bam_chrs[i] + 3, chromosome))
      {
        return i;
      }
    }
  }
  return -1;
}

uint64_t get_snp_count(char *file)
{

  uint64_t count = 0;
  // Open BGZF compressed vcf file in read mode
  BGZF *bgzf = bgzf_open(file, "r");
    if (!bgzf) {
        REprintf("Failed to open BGZF file: %s\n", file);
        return -1; // return -1 if there is an error
    } else {
        //Rprintf("Opened BGZF file: %s\n", file);
        bgzf_close(bgzf); // Ensure we close the file after opening
        return 0;
    }

  // Buffer for reading lines
  kstring_t line = {0, 0, NULL}; // Initialize kstring_t
  // size_t len = 0;

  // read lines from BGZF file
  while (bgzf_getline(bgzf, '\n', &line) >= 0)
  { // reaD A SINGLE LINE FROM BGZF FILE UNTIL newline character, store line in provided buffer 'line'
    // ignore header lines (start with #)
    if (line.s[0] != '#')
    {
      count++;
    }
  }
  // clean up
  free(line.s);     // free the allocated line buffer
  bgzf_close(bgzf); // close the bgzf handle

  return count;
}

int snp_pileup_main(arguments arguments)
{
  // debugPrint("Debug: enter program main", arguments.debug_mode);
  //std::cout << "In Main yay" << std::endl;
  clock_t start = clock();

  // initialize engine
  int i = 0;
  int n = arguments.args.size() - 2; // n is the number of files. currently hardcoded at 2
  bam_mplp_t iter;
  const bam_pileup1_t **plp = (const bam_pileup1_t **)calloc(n, sizeof(bam_pileup1_t *));
  int *n_plp = (int *)calloc(n, sizeof(int));
  mplp_ref_t mp_ref = MPLP_REF_INIT;
  struct arguments *conf = &arguments;
  hts_verbose = 1;

  // Load VCF file
  // debugPrint("Debug: About to load vcf file", arguments.debug_mode);

  // open file
  BGZF *bgzf = bgzf_open(arguments.args[0].c_str(), "r");
  if (!bgzf)
  {
    Rprintf("Failed to open BGZF file: %s\n", arguments.args[0].c_str());
    return 1;
  }

  std::vector<std::string> header_columns;

  //  Initialize BAM file data
  mplp_aux_t **data = (mplp_aux_t **)calloc(n, sizeof(mplp_aux_t *));
  for (i = 0; i < n; ++i)
  {
    data[i] = (mplp_aux_t *)calloc(1, sizeof(mplp_aux_t));

    // open file
    hFILE *hfp = hopen(arguments.args[i + 2].c_str(), "r");
    if (!hfp)
    {
      REprintf("Failed to open BAM file: %s\n", arguments.args[i + 2].c_str());
      return 1;
    }

    // ommitting detect file format because ... ???

    data[i]->fp = hts_hopen(hfp, arguments.args[i + 2].c_str(), "rb");
    if (!data[i]->fp)
    {
      // Error handling: Couldn't open the file
      REprintf("Couldn't open sequence file: %s (%s)\n", 
         arguments.args[i + 2].c_str(), strerror(errno));
      // hclose(hfp); // Clean up the file handle
      if (hclose(hfp) != 0)
      {
        // Handle the error, e.g., log an error message or throw an exception.
        REprintf("Error: Failed to close file handle.");
      }
      return 1;
    }

    hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0);

    data[i]->conf = conf;
    data[i]->h = sam_hdr_read(data[i]->fp);
    data[i]->ref = &mp_ref;
  }

  // Read BAM header
  bam_hdr_t *hdr = data[0]->h;

  // Start pileup engine
  iter = bam_mplp_init(n, mplp_func, (void **)data);
  if (!arguments.ignore_overlaps)
  {
    bam_mplp_init_overlaps(iter);
  }
  bam_mplp_set_maxcnt(iter, arguments.max_depth);

  // setup output file

  string fname = string(arguments.args[1]);

  if (!ends_with(fname, ".gz"))
  {
    fname += ".gz";
  }
  
  //  DON'T CLOSE test_output HERE because if you are here, test_output is null
  BGZF *output_file = NULL;

  arguments.gzippedPointer = bgzf_open(fname.c_str(), "w+");
  if (!arguments.gzippedPointer)
  {
    return 1;
  }

  ostringstream output;

  // output header to file
  output << "Chromosome,Position,Ref,Alt";
  for (int i = 0; i < n; ++i)
  {
    output << ",File" << (i + 1) << "R";
    output << ",File" << (i + 1) << "A";
    output << ",File" << (i + 1) << "E";
    output << ",File" << (i + 1) << "D";
  }
  output << "\n";

  // output header
  std::string header = output.str();

  // make sure there is a valid bgzf pointer
  if (!arguments.gzippedPointer) {
    Rprintf("Error: BGZF pointer is NULL!");
}

// write the header to the output file 
  gzip_output(arguments, output.str(), output_file); // this wrote the header~
 
  // Process VCF records
  int ret;
  int tid, pos;
  vector<file_info> f_info = vector<file_info>();
  bool first = true;
  kstring_t line = {0, 0, NULL};
  uint64_t current_count = 0;
  //float last_progress = 100.0;
  bool have_snpped = false;
  // ready to loop
  while (bgzf_getline(bgzf, '\n', &line) >= 0)
  {
    // Skip header lines
    // debugPrint("Debug: in the loop, no parse yet", arguments.debug_mode);
    if (line.s[0] == '#')
      continue;

    // Parse the record
    std::istringstream iss(line.s);
    std::vector<std::string> fields;
    std::string field;

    // debugPrint("Debug: did parse", arguments.debug_mode);

    while (std::getline(iss, field, '\t'))
    {
      fields.push_back(field);
      // debugPrint("Debug: in field while", arguments.debug_mode);
      // debugPrint(field, arguments.debug_mode);
    }

    // Validate record
    if (fields.size() < 4)
    {
      Rprintf("Malformed VCF record: %s\n", line.s);
      continue;
    }


    std::string chrom = fields[0];

    int vcf_pos = std::stoi(fields[1]);
    std::string ref = fields[2];
    // debugPrint(ref, arguments.debug_mode);
    std::string alt = fields[3];


    // Filter for SNPs (single nucleotide variants only)
    if (ref.size() != 1 || alt.size() != 1)
    {
      //std::cout << "Skipping non-snp" << std::endl;
      // debugPrint("Skipping non-SNP record: " + std::string(line.s), arguments.debug_mode);
      continue;
    }
    // Map VCF chromosome and position to BAM target indices
    char chrom_buffer[chrom.size() + 1];
    std::strcpy(chrom_buffer, chrom.c_str()); // copy string into mutable buffer

    int vcf_tid = vcf_chr_to_bam(chrom_buffer, hdr->target_name, hdr->n_targets);
    current_count++;

    // plp processing loop
    while ((ret = bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0)
    {
      if (!first && tid > vcf_tid)
        continue; // Skip ahead if pileup is beyond VCF tid
      if (!first && tid == vcf_tid && pos < vcf_pos)
        continue; // Skip if position is earlier

      // moving forward....
      //std::cout << "Checking BAM for SNP at Chrom: " << chrom << " Pos: " << vcf_pos << std::endl;

      if (tid == vcf_tid && pos == vcf_pos)
      {
        // Process the pileup data for each file
        // output data to output file
        bool is_not_zero = false;
        bool fails_min = false;
        for (i = 0; i < n; ++i)
        {
          file_info this_file = {0, 0, 0, 0}; // Initialize counts

          // check that there are sufficient reads for this file
          if (n_plp[i] >= arguments.min_read_counts[i])
          {
            // loop and pile
            for (int j = 0; j < n_plp[i]; ++j)
            {
              const bam_pileup1_t *p = plp[i] + j;
              char base = bam_seqi(bam_get_seq(p->b), p->qpos);
              if (p->is_del)
              {
                this_file.deletions++;
              }
              else if (seq_nt16_str[(unsigned int)base] == ref[0])
              {
                this_file.refs++;
              }
              else if (seq_nt16_str[(unsigned int)base] == alt[0])
              {
                this_file.alts++;
              }
              else
              {
                this_file.errors++;
              }
            }
          }
          else
          {

            fails_min = true;
          }
          if (this_file.refs != 0 || this_file.alts != 0 || this_file.errors != 0)
          {
            is_not_zero = true;
          }
          f_info.push_back(this_file);
        }

        // Write the SNP information to the output
        if (is_not_zero && !fails_min)
        {
          output.str("");
          output.clear();
          output << chrom << "," << vcf_pos << "," << ref << "," << alt;
          for (i = 0; i < n; ++i)
          {
            output << "," << f_info[i].refs
                   << "," << f_info[i].alts
                   << "," << f_info[i].errors
                   << "," << f_info[i].deletions;
          }
          output << "\n";
          if (!arguments.gzippedPointer) {
    Rprintf("ERROR: BGZF pointer is NULL!");
}
          // writey write
          gzip_output(arguments, output.str(), output_file);
          
        }
        f_info.clear();
        break;
      }
      if (arguments.pseudo_snps && !have_snpped && (((pos + 1) % arguments.pseudo_snps) == 0))
      {
        bool is_not_zero = false;
        bool fails_min = false;
        for (i = 0; i < n; i++)
        {
          file_info this_file = file_info();
          this_file.refs = 0;
          this_file.alts = 0;
          this_file.errors = 0;
          this_file.deletions = 0;
          if (n_plp[i] >= arguments.min_read_counts[i])
          {
            for (int j = 0; j < n_plp[i]; ++j)
            {
              const bam_pileup1_t *p = plp[i] + j;
              int c = p->qpos < p->b->core.l_qseq
                          ? bam_get_qual(p->b)[p->qpos]
                          : 0;
              if (c == 0)
              {
                continue; // no
              }
              if (c < arguments.min_base_quality)
              {
                continue; // skip anything with quality below threshold
              }
              this_file.refs++;
            }
          }
          else
          {
            fails_min = true;
          }
          if (this_file.refs != 0)
          {
            is_not_zero = true;
          }
          f_info.push_back(this_file);
        }
        // add a pseudo snp!
        if (is_not_zero && !fails_min)
        {
          output.str("");
          output.clear();
          output << hdr->target_name[tid];
          output << ",";
          output << (pos + 1);
          output << ",.,.";
          for (i = 0; i < n; i++)
          {
            output << ",";
            output << f_info[i].refs;
            output << ",0,0,0";
          }
          output << "\n";
          
          //write to file 
          gzip_output(arguments, output.str(), output_file);
        }
        f_info.clear();
      }
      // Break the loop if the BAM position surpasses the VCF position
      if (tid > vcf_tid || (tid == vcf_tid && pos > vcf_pos))
      {
        break;
      }
    }
  }

  // Cleanup
  free(line.s);
  bgzf_close(bgzf);
  bam_mplp_destroy(iter);

  for (i = 0; i < n; ++i)
  {
    if (data[i]->fp)
      hts_close(data[i]->fp);
    free(data[i]);
  }
  free(data);

  // debugPrint("Debug: Completed snp_pileup_main", arguments.debug_mode);
  //printf("Finished in %f seconds.\n", (clock() - start) / (double)CLOCKS_PER_SEC);

  bgzf_close(arguments.gzippedPointer);

  return 0;
}

//' Run Snp Pileup
//'
//' This function runs snp_pileup
//'
//' @return None.
//' @export
extern "C" SEXP run_snp_pileup_logic(SEXP args_in)
{

  arguments args;

  args.count_orphans = INTEGER(VECTOR_ELT(args_in, 0))[0] != 0;
  //std::cout << "count_orphans: " << args.count_orphans << "\n";
  args.ignore_overlaps = INTEGER(VECTOR_ELT(args_in, 1))[0] != 0;// LOGICAL(VECTOR_ELT(args_in, 1))[0]; // this is of type integer (13)
  //std::cout << "ignore_overlaps: " << args.ignore_overlaps << "\n";
  args.min_base_quality = INTEGER(VECTOR_ELT(args_in, 2))[0]; 
  //std::cout << "min_base_quality: " << args.min_base_quality << "\n";
  args.min_map_quality = INTEGER(VECTOR_ELT(args_in, 3))[0]; 
  //std::cout << "min_map_quality: " << args.min_map_quality << "\n";

  // Convert SEXP to vector<int> for min_read_counts
  SEXP read_counts = VECTOR_ELT(args_in, 4);
  int len = LENGTH(read_counts);
  for (int i = 0; i < len; ++i)
  {
    args.min_read_counts.push_back(INTEGER(read_counts)[i]);
  }

  args.max_depth = INTEGER(VECTOR_ELT(args_in, 5))[0];   // Convert SEXP to int
  args.pseudo_snps = INTEGER(VECTOR_ELT(args_in, 6))[0]; // Convert SEXP to int

    SEXP char_args = VECTOR_ELT(args_in, 7);
  int char_len = LENGTH(char_args);
  for (int i = 0; i < char_len; ++i)
  {
    args.args.push_back(CHAR(STRING_ELT(char_args, i)));
    //std::cout << "args[" << i << "]: " << args.args[i] << "\n";
  }
  args.gzippedPointer = nullptr;

  //std::cout << "parsed the args" << std::endl;
  int status = snp_pileup_main(args);

  if (status != 0)
  {
    // Rcpp::stop("Program terminated with errors.");
  }
  return ScalarInteger(status);
}

//' Check Htslib Version
//'
//' This function prints the version of HTSLIB being used for snp-pileup
//'
//' @return None.
//' @export
extern "C" SEXP htslib_version_cpp()
{
  // Rcpp::Rcout << "Htslib linked successfully - snp-pileup!" << std::endl;
  //std::cout << "HTSlib version:" << hts_version() << std::endl;
  return(mkString(hts_version()));
}

