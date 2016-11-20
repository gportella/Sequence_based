#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/seq_io.h>
#include <sstream>

// Readme TODO -- > this should be easily incorporated
// in a generic program that calls the right function for
// mapping (exact or not) depending on what the user wants
// instead of having two sepparate programs. I don't quite
// the energy at the moment.

using namespace seqan;
using namespace std;

struct findMatchOptions {
  bool bCompressedOutput;
  bool bReverse;
  CharString seqId; // to store an index, e.g. the chromosome number
  int edit_dist;
  CharString needlesFileName;
  CharString haystackFileName;
  CharString outFileName;

  // I guess this is how to initialize
  findMatchOptions() : bCompressedOutput(false), bReverse(true) {}
};

ArgumentParser::ParseResult parseCommandLine(findMatchOptions &parseOptions,
                                             int argc, char const **argv) {
  // Setup ArgumentParser.
  ArgumentParser parser("approxseqmatch");
  setShortDescription(parser, "Find approximate sequence matches ");
  addDescription(parser, "Find approximate matches in fasta from multiple "
                         "input fasta file.");
  addUsageLine(parser, "[\\fIOPTIONS\\fP] ");
  setVersion(parser, "1.0");
  setDate(parser, "October 2016");

  // Define Options
  addOption(parser,
            ArgParseOption("i", "sequence_file",
                           "A fasta input file with sequence(s) to search for.",
                           ArgParseArgument::INPUT_FILE));
  addOption(parser, ArgParseOption("c", "chr_file",
                                   "A fasta input file with a single, longer "
                                   "sequence (e.g. chromosome).",
                                   ArgParseArgument::INPUT_FILE));
  addOption(parser,
            ArgParseOption(
                "o", "result_file",
                "Output file with id record of the input sequences appearing "
                "in the longer sequence, its sequence, a index_id and the "
                "begin,end position in the longer sequence.",
                ArgParseArgument::OUTPUT_FILE));
  setDefaultValue(parser, "result_file", "seqmatch.txt");
  addOption(
      parser,
      ArgParseOption("id", "index_id",
                     "An integer/text referencing the single, longer sequence "
                     "(e.g. chromosome number) ",
                     ArgParseArgument::STRING, "TEXT"));
  addOption(parser,
            ArgParseOption("ed", "edit_distance",
                           "An integer giving the maximum edit distance "
                           "for approximate patter matchin ",
                           ArgParseArgument::INTEGER, "TEXT"));
  setDefaultValue(parser, "edit_distance", "2");
  addOption(parser,
            ArgParseOption(
                "gz", "gzip_compress",
                "Compress the output using Gzip (If result_file becomes "
                "too large you might have memory problems with the current "
                "implementation)."));
  addOption(parser, ArgParseOption("r", "search_reverse",
                                   "Search also the reverse complement "
                                   "of each sequence"));
  setValidValues(parser, "sequence_file", "fna fa");
  setValidValues(parser, "chr_file", "fa");
  setRequired(parser, "sequence_file");
  setRequired(parser, "chr_file");
  setRequired(parser, "index_id");

  // Parse command line.
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // Only extract options if the program continues after parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

  // Extract option values.
  getOptionValue(parseOptions.needlesFileName, parser, "sequence_file");
  getOptionValue(parseOptions.haystackFileName, parser, "chr_file");
  getOptionValue(parseOptions.outFileName, parser, "result_file");
  getOptionValue(parseOptions.seqId, parser, "index_id");
  getOptionValue(parseOptions.edit_dist, parser, "edit_distance");
  parseOptions.bCompressedOutput = isSet(parser, "gzip_compress");
  parseOptions.bReverse = isSet(parser, "search_reverse");

  return ArgumentParser::PARSE_OK;
}

int main(int argc, char const **argv) {
  // Parse the command line
  findMatchOptions parseOptions;
  ArgumentParser::ParseResult res = parseCommandLine(parseOptions, argc, argv);

  // If parsing did not work, then exit with code 1
  // Otherwise exit with code 0
  if (res != ArgumentParser::PARSE_OK)
    return res == ArgumentParser::PARSE_ERROR;

  CharString sequenceFileName = parseOptions.needlesFileName;
  CharString chrFileName = parseOptions.haystackFileName;
  CharString resultFileName = parseOptions.outFileName;
  CharString chId = parseOptions.seqId;
  int edit_dist = parseOptions.edit_dist;
  if (parseOptions.bCompressedOutput)
    resultFileName += ".gz";

  // declarations for fasta inputs
  StringSet<CharString> ids;
  StringSet<Dna5String> seqs;
  CharString idchr;
  Dna5String seqchr;
  SeqFileIn seqFileIn;
  SeqFileIn chrFileIn;

  // open sequence file
  if (!open(seqFileIn, toCString(sequenceFileName))) {
    cerr << "ERROR: Cound not open input file (sequences).\n";
    return 1;
  }
  try {
    readRecords(ids, seqs, seqFileIn);
  } catch (Exception const &e) {
    cout << "ERROR: " << e.what() << endl;
    return 1;
  }
  // open chromosome file
  if (!open(chrFileIn, toCString(chrFileName))) {
    cerr << "ERROR: Cound not open the file.\n";
    return 1;
  }
  try {
    readRecord(idchr, seqchr, chrFileIn);
  } catch (Exception const &e) {
    cout << "ERROR: " << e.what() << endl;
    return 1;
  }

  CharString haystack = seqchr;

  // I think I have to declare it anyway, even if I don't use it
  std::stringstream ss, comp;
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;

  // I branch here, otherwise I will have to evaluate the IF everytime
  // I branch one for bCompressedOutput and once for bReverse
  // It could be made a bit more elegant by having a function do the work
  // but I am ok for the moment
  if (!parseOptions.bCompressedOutput) {
    ofstream outFile(toCString(resultFileName));
    if (parseOptions.bReverse) {
      for (unsigned i = 0; i < length(ids); ++i) {
        Finder<CharString> finder(haystack);
        Pattern<CharString, Myers<>> pattern(seqs[i]);
        while (find(finder, pattern, -edit_dist)) {
          while (findBegin(finder, pattern, getScore(pattern))) {
            outFile << toCString(ids[i]) << "\t" << infix(finder) << "\t";
            outFile << chId << "\t" << beginPosition(finder) << ","
                    << endPosition(finder) << "\t" << getScore(pattern) << "\t"
                    << "+" << std::endl;
          }
        }
        // Now we look for its reverse complement
        Finder<CharString> finder_rev(haystack);
        Dna5StringReverseComplement rev_comp(seqs[i]);
        Pattern<CharString, Myers<>> pattern_rev(rev_comp);
        while (find(finder_rev, pattern_rev, -edit_dist)) {
          while (findBegin(finder_rev, pattern_rev, getScore(pattern_rev))) {
            outFile << toCString(ids[i]) << "\t" << infix(finder_rev) << "\t";
            outFile << chId << "\t" << beginPosition(finder_rev) << ","
                    << endPosition(finder_rev) << "\t" << getScore(pattern_rev)
                    << "\t"
                    << "-" << std::endl;
          }
        }
      }
    } else {
      for (unsigned i = 0; i < length(ids); ++i) {
        Finder<CharString> finder(haystack);
        Pattern<CharString, Myers<>> pattern(seqs[i]);
        while (find(finder, pattern, -edit_dist)) {
          while (findBegin(finder, pattern, getScore(pattern))) {
            outFile << toCString(ids[i]) << "\t" << infix(finder) << "\t";
            outFile << chId << "\t" << beginPosition(finder) << ","
                    << endPosition(finder) << "\t" << getScore(pattern)
                    << std::endl;
          }
        }
      }
    }
  } else {
    if (parseOptions.bReverse) {
      for (unsigned i = 0; i < length(ids); ++i) {
        Finder<CharString> finder(haystack);
        Pattern<CharString, Myers<>> pattern(seqs[i]);
        while (find(finder, pattern, -edit_dist)) {
          while (findBegin(finder, pattern, getScore(pattern))) {
            ss << toCString(ids[i]) << "\t" << infix(finder) << "\t";
            ss << chId << "\t" << beginPosition(finder) << "\t"
               << endPosition(finder) << "\t" << getScore(pattern) << "\t"
               << "+" << std::endl;
          }
        }
        // Now we look for its reverse complement
        Finder<CharString> finder_rev(haystack);
        Dna5StringReverseComplement rev_comp(seqs[i]);
        Pattern<CharString, Myers<>> pattern_rev(rev_comp);
        while (find(finder_rev, pattern_rev, -edit_dist)) {
          while (findBegin(finder_rev, pattern_rev, getScore(pattern_rev))) {
            ss << toCString(ids[i]) << "\t" << infix(finder_rev) << "\t";
            ss << chId << "\t" << beginPosition(finder_rev) << "\t"
               << endPosition(finder_rev) << "\t" << getScore(pattern_rev)
               << "\t"
               << "-" << std::endl;
          }
        }
      }
    } else {
      for (unsigned i = 0; i < length(ids); ++i) {
        Finder<CharString> finder(haystack);
        Pattern<CharString, Myers<>> pattern(seqs[i]);
        while (find(finder, pattern, -edit_dist)) {
          while (findBegin(finder, pattern, getScore(pattern))) {
            ss << toCString(ids[i]) << "\t" << infix(finder) << "\t";
            ss << chId << "\t" << beginPosition(finder) << "\t"
               << endPosition(finder) << "\t" << getScore(pattern) << std::endl;
          }
        }
      }
    }
    ofstream compOutFile(toCString(resultFileName));
    in.push(boost::iostreams::gzip_compressor());
    in.push(ss);
    boost::iostreams::copy(in, comp);
    compOutFile << comp.str();
  }

  return 0;
}
