#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

struct findMatchOptions {
  bool b_verbose;
  bool b_validation_set;
  bool b_read_rna;
  bool b_input_ds;
  bool b_rnafold;
  bool b_constrain;
  bool b_silent;
  CharString needlesFileName;
  CharString outFileName;

  // I guess this is how to initialize
  findMatchOptions()
      : b_verbose(false), b_validation_set(false), b_read_rna(false),
        b_input_ds(false), b_rnafold(false), b_constrain(false),
        b_silent(true) {}
};

// ArgumentParser::ParseResult parseCommandLine(findMatchOptions &parseOptions,
//                                             int argc, char const **argv);

ArgumentParser::ParseResult parseCommandLine(findMatchOptions &parseOptions,
                                             int argc, char const **argv) {
  // Setup ArgumentParser.
  ArgumentParser parser("tfo_finder");
  setShortDescription(parser, "Finds TFOs  ");
  addDescription(parser, "Finds TFOs in "
                         "fasta from "
                         "multiple "
                         "input fasta file.");
  addUsageLine(parser, "[\\fIOPTIONS\\fP] ");
  setVersion(parser, "1.0");
  setDate(parser, "October 2016");

  // Define Options
  addOption(parser, ArgParseOption("i", "sequence_file", "A fasta input file "
                                                         "with sequence(s) to "
                                                         "search for.",
                                   ArgParseArgument::INPUT_FILE));
  addOption(parser, ArgParseOption("o", "result_file", "Output file ",
                                   ArgParseArgument::OUTPUT_FILE));
  addOption(parser, ArgParseOption("rna", "input_rna", "Input data is  "
                                                       "RNA."));
  addOption(
      parser,
      ArgParseOption("ds", "double_stranded",
                     "Input data is genome-like. It will "
                     "search for binding partners along"
                     " both fwd and reverse. It does "
                     "report the TFO sequence in the fasta seq, "
                     "along with the strand (+/-) and coordinates "
                     "in the header. Coordinates are defined wrt the fwd "
                     "strand, if the strand is - the start position is "
                     "actually the largest value of the pair of coordinates."));

  addOption(parser, ArgParseOption("fold", "fold_rna", "Fold RNA with RNAfold "
                                                       "with noconstr now."));
  addOption(parser,
            ArgParseOption("const", "constrained_fold", "Fold RNA with RNAfold "
                                                        "using constraints."));
  addOption(parser,
            ArgParseOption("val", "validation_set", "Read the conditions "
                                                    "from the fasta header "
                                                    "and do the "
                                                    "validation."));
  addOption(parser, ArgParseOption("v", "be_verbose", "Be verbose in "
                                                      "what you do."));
  addOption(parser,
            ArgParseOption("silent", "be_silent",
                           "Do not write to "
                           "std output, except for errors and warnings"));
  setDefaultValue(parser, "result_file", "out_predict.txt");
  setDefaultValue(parser, "result_file", "found_tfos.fa");
  setValidValues(parser, "sequence_file", "fna fa fq");
  setValidValues(parser, "result_file", "fa fna");
  setRequired(parser, "sequence_file");

  // Parse command line.
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // Only extract options if the
  // program continues after
  // parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

  // Extract option values.
  getOptionValue(parseOptions.needlesFileName, parser, "sequence_file");
  getOptionValue(parseOptions.outFileName, parser, "result_file");
  parseOptions.b_verbose = isSet(parser, "be_verbose");
  parseOptions.b_validation_set = isSet(parser, "validation_set");
  parseOptions.b_read_rna = isSet(parser, "input_rna");
  parseOptions.b_input_ds = isSet(parser, "double_stranded");
  parseOptions.b_rnafold = isSet(parser, "fold_rna");
  parseOptions.b_constrain = isSet(parser, "constrained_fold");
  parseOptions.b_silent = isSet(parser, "be_silent");

  return ArgumentParser::PARSE_OK;
}
