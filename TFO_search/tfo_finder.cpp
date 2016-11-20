#include "commandline_parse.hpp"
#include "find_tfo.hpp"
#include "seqan/sequence.h"
#include <boost/tokenizer.hpp>
#include <functional>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/stream.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace seqan;
using namespace std;
using namespace boost;

///////////////////////////////////////////////////////////////////////////////
// :Info
///////////////////////////////////////////////////////////////////////////////
// Find triplex forming oligonucleotieds (TFOs) in a given (set of) sequence(s)
// The sequences should be introduced using the fasta format
// There are a couple of modes, mainly.
//
// a) single stranded --> gives back the sequences of the binding partner of
//                        the TFOs it finds. It gives back both the sense
//                        and antisense. Works for DNA and RNA, but using
//                        always the same set of rules. You can also request
//                        using RNAfold to fold the single stranded, this only
//                        makes sense for RNA. Futhermore, if you suply
//                        constraints in the RNAfold format (modified fastq)
//                        these can also be used. Only TFOs found in an unpaired
//                        region will be reported.
//
// b) double stranded --> gives back the TFO sequences that would bind to
//                        the double stranded DNA (e.g. genome). It searches
//                        both the fwd and rev strands. The output fasta file
//                        containts the coordinates of the motif wrt the
//                        + strand in the header. E.g., if it reports |1:10|-|
//                        this means that the binding site for the TFO is
//                        in the reverse strand, and its 5' starts at 10 and
//                        the 3' finishes at 1.
//
// TODO --> better document the output format
//
// The code is organised such that you have the command line parsing
// in one hpp file, the routines to compute the predicted stabilities of TFOs
// in another one, and main class that does the job. The class is called
// LncRna because that was the original intention, but the code grew to
// treat double stranded DNA (tribially). Each object of the lncrna class
// has one or severall Tfo objects, which store the TFO for each input
// sequence. find_tfo.hpp, where the code is, should be well documented.
// The class makes heavy usage of Tags to differentiate between function
// calls that do the same but speciallized (eg. folding vs non-folding).
//
// Quite important: do use -O3 for production, makes the codes way faster
// although the compiler takes a tad longer
///////////////////////////////////////////////////////////////////////////////

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

int main(int argc, char const **argv) {

  // OpenMP setting the maximum number of threads possible
  int nProcessors = omp_get_max_threads();
  omp_set_num_threads(nProcessors);
  // Parse the command line
  findMatchOptions parseOptions;
  ArgumentParser::ParseResult res = parseCommandLine(parseOptions, argc, argv);

  // If parsing did not work, then exit with code 1
  // Otherwise exit with code 0
  if (res != ArgumentParser::PARSE_OK)
    return res == ArgumentParser::PARSE_ERROR;

  CharString sequenceFileName = parseOptions.needlesFileName;
  CharString resultFileName = parseOptions.outFileName;
  bool b_validation_set = parseOptions.b_validation_set;
  bool b_read_rna = parseOptions.b_read_rna;
  bool b_double_strand = parseOptions.b_input_ds;
  bool b_rnafold = parseOptions.b_rnafold;
  bool b_verbose = parseOptions.b_verbose;
  bool b_constrain = parseOptions.b_constrain;
  bool b_silent = parseOptions.b_silent;
  // declarations for fasta inputs
  StringSet<CharString> ids;
  StringSet<Dna5String> dseqs;
  StringSet<Rna5String> rseqs;
  StringSet<CharString> quals;
  SeqFileIn seqFileIn;
  SeqFileIn chrFileIn;

  // to sort out the options
  bool b_has_constraints = false;

  if (b_verbose) {
    cout << "Set " << nProcessors << " OpenMP threads" << endl;
  }

  // check the extension and decide if we go for fasta or fastq (constraints)
  string fn = toCString(sequenceFileName);
  if (fn.substr(fn.find_last_of(".") + 1) == "fq") {
    // we assume right now that the file might be valid
    b_has_constraints = true;
    b_double_strand = false;
  }
  // if we don't read constraints, we can do constrained folding
  if (b_constrain && !b_rnafold) {
    if (!b_silent) {
      cout << "Note: You requested the use of constraints, assuming "
              "you want to use RNAfold "
           << endl;
    }
    b_rnafold = true;
    b_double_strand = false;
  }
  if (b_double_strand) {
    // it does not make any sense to search both ways
    // if we expect the folding of a single stranded RNA
    b_rnafold = false;
  }

  if (!b_has_constraints && b_constrain) {
    cout << "Warning: You requested the use of constraints, but the input does "
            "not have them "
         << endl;
    b_constrain = false;
  }
  // print out options
  if (!b_silent) {
    cout << "--------------------------------------------------" << endl;
    cout << "These are the options we will use" << endl;
    cout << "Validating predictor : " << boolalpha << b_validation_set << endl;
    cout << "RNA : " << boolalpha << b_read_rna << endl;
    cout << "Read fastq : " << boolalpha << b_has_constraints << endl;
    cout << "Search both sense (double stranded) : " << boolalpha
         << b_double_strand << endl;
    cout << "Do folding : " << boolalpha << b_rnafold << endl;
    cout << "Use constraints : " << boolalpha << b_constrain << endl;
    cout << "--------------------------------------------------" << endl;
  }

  // open sequence file
  if (!open(seqFileIn, toCString(sequenceFileName))) {
    cerr << "ERROR: Cound not open input file (sequences).\n";
    return 1;
  }
  try {
    // in case you want to read RNA, should be less frequently
    if (b_read_rna) {

      if (b_has_constraints) {
        readRecords(ids, rseqs, quals, seqFileIn);
      } else {
        readRecords(ids, rseqs, seqFileIn);
      }
      dseqs = rseqs; // this should translate RNA to DNA
      // how can I delete the rseq now?? TODO
    } else {
      if (b_has_constraints) {
        readRecords(ids, dseqs, quals, seqFileIn);
      } else {
        readRecords(ids, dseqs, seqFileIn);
      }
    }
  } catch (Exception const &e) {
    cout << "ERROR: " << e.what() << endl;
    return 1;
  }

  if (b_validation_set) {
    vector<std::string> tokens;
    std::string::size_type sz;
    for (unsigned i = 0; i < length(ids); ++i) {
      // Maybe we need to refactor this later, now I'm keeping it
      string text = toCString(ids[i]);
      split(text, '|', tokens);
      double pH = stod(tokens[2], &sz);
      double conc = stod(tokens[3], &sz);
      tokens.clear();
      auto dg_tm = compute_tfo_tm(dseqs[i], conc, pH);
      // You could also use the Tag_dna_nn() if you wanted to use
      // Roberts&Crothers parameters for DNA-DNA:RNA triplexs
      // auto dg_tm = compute_tfo_tm(dseqs[i], conc, pH, Tag_dna_nn());
      cout << ids[i] << " Tm: " << dg_tm.second + 273.15 << endl;
    }
  } else {

    if (!b_rnafold) {

      if (b_double_strand) {
        // Let's go double stranded, no folding
        vector<LncRna> results = find_tfo_from_fasta(dseqs, ids, Tag_ds());
        if (print_putative_tfos_to_fasta(results, resultFileName, Tag_ds()) ==
            1) {
          cout << "ERROR printing the results" << endl;
          return 1;
        }
      } else {
        // Let's go single stranded, no folding
        vector<LncRna> results = find_tfo_from_fasta(dseqs, ids);
        if (print_putative_tfos_to_fasta(results, resultFileName) == 1) {
          cout << "ERROR printing the results" << endl;
          return 1;
        }
      }
    } else {

      if (!b_constrain) {

        // the trouble is that I can't declare results it without initializing
        // (it will call de constructor), so I have to list the printing after
        // You might want to fix the behaviour with the constructor later
        // If I place the print outside the ifs "results" looses scope
        // (desctructed) and can't be printed

        // Let's go no constraints
        vector<LncRna> results = find_tfo_from_fasta(dseqs, ids, Tag_folding());
        if (print_putative_tfos_to_fasta(results, resultFileName) == 1) {
          cout << "ERROR printing the results" << endl;
        }

      } else {

        // Let's go with constrained folding
        vector<LncRna> results =
            find_tfo_from_fasta(dseqs, ids, quals, Tag_constraints());
        if (print_putative_tfos_to_fasta(results, resultFileName) == 1) {
          cout << "ERROR printing the results" << endl;
          return 1;
        }
      }
    }
  }
  return 0;
  ///////////////////////////////////////////////////////////////////////////////
  // The end of the main code
  ///////////////////////////////////////////////////////////////////////////////
}
