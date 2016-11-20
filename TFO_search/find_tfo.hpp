#include "nn_tm.hpp"
#include <boost/algorithm/string/replace.hpp>
#include <boost/regex.hpp>
#include <functional>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <vector>

#include "do_rnafold.h"

const double tm_threshold(45); // arbitrary
const double tfo_conc(12);     // mM , I think
const double pH(7.0);
const int max_fold_rna(1200);

using namespace seqan;
using namespace std;

struct Tag_folding {};
struct Tag_constraints {};
struct Tag_complement {};
struct Tag_ds {}; // double stranded search both sense

///////////////////////////////////////////////////////////////////////////////
// :Info:
///////////////////////////////////////////////////////////////////////////////
// Classes and functions to find TFOs in a sequence, filter them according to
// predicted Tm and writting out a fasta file with the findings.
// The fasta file contains the id of the original fasta sequence, the DG
// and Tm of the sequence found and the sequence  of the binding partner of the
// found TFO. I write a second entry containing the - strand, such that
// we can also search for that one if we want later on.
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// ... Class LncRna: for sequence entry it finds the possible TFOs using their
// Tm as limit, and stores their id, the sequence and the postion in the
// original sequence. ...
// Containts private functions to search putative TFOs by sequence (given a
// minimum length of 12 bases), and a prive function that computes the Tm
// for those putative sites and fills up the information.
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Functor for mapping Pietro's notation from DMS+RNA-seq to RnaFold constraint
// notation. This is used by Seqan's  ModifiedString modifier, which allows to
// change string's without copying them. They way I think it works is by
// creating class that points at them, and once called they give back the
// watever you told them to do. You can either call pre-designe functors from
// Seqan or make you own.
//
// This is the mapping to fix quals from Pietro-type files
// 1 = base paired --> | in rnafold
// 0 = non-base paired --> . in rnafold
// n = no info (corresponds to G or U) --> x in rnafold
///////////////////////////////////////////////////////////////////////////////

struct Pietro2RnafoldFunctor : public std::unary_function<char, char> {
  inline char operator()(char x) const {
    if (('0' == x)) {
      return 'x';
    } else if (('1' == x)) {
      return '|';
    } else if (('n' == x)) {
      return '.';
    } else {
      return x;
    }
  }
};
// Above I could add a final condition and break, but I want to keep it
// in case I feed the code the right nomenclature. Anyway RNAfold complaints
// if it can not read the constraints, and computes the energy without it
// anyway.
/* else {
     cout << "Found weird symbol in constraints" << endl;
     cout << x << endl;
     exit(1);
 }*/

///////////////////////////////////////////////////////////////////////////////

class LncRna {

public:
  LncRna(); // default constructor, just for declarations
  LncRna(Dna5String haystack, CharString id, int index);
  // Tag dispatching here...
  LncRna(Dna5String haystack, CharString id, int index, Tag_ds const & /*Tag*/);
  LncRna(Dna5String haystack, CharString id, int index,
         Tag_folding const & /*Tag*/);
  LncRna(Dna5String haystack, CharString id, CharString quals, int index,
         Tag_constraints const & /*Tag*/);
  bool found_tfos;
  bool struct_computed;
  int index_id; // no idea why I need it, but it's in the python code
  CharString id_lncrna;

  class Tfos {
  public:
    // We might not need it, actually
    tuple<int, int> position; // hold the position in the genome
    string strand = "+";      // + or -, initialized to +
    Dna5String seq;           // TODO: hold the whole sequence ?
    struct thermodynamics {   // want to write the names explicitelty, hence
                              // I use a struct and not a pair
      double dg = 0;
      double tm = 0;
    } thermo;
  };

  // make an array of tfos to hold all the records
  vector<Tfos> lnctfos; // what about the memory?

private:
  // this is quite horrible...
  template <typename TText1>
  vector<tuple<int, int>> search_by_seq(TText1 const &sequence,
                                        boost::regex const &re);

  template <typename TText1, typename TText2, typename TText3>
  void find_tfo(TText1 const &local_record, TText2 id_lncrna, TText3 index);

  template <typename TText1, typename TText2, typename TText3>
  void find_tfo(TText1 const &local_record, TText2 id_lncrna, TText3 index,
                Tag_ds const & /*Tag*/);

  template <typename TText1, typename TText2, typename TText3>
  void find_tfo(TText1 const &local_record, TText2 local_id_lncrna,
                TText3 index, Tag_folding const & /*Tag*/);

  template <typename TText1, typename TText2, typename TText3, typename TText4>
  void find_tfo(TText1 const &local_record, TText2 local_id_lncrna,
                TText3 &quals, TText4 index, Tag_constraints const & /*Tag*/);

  template <typename TText1, typename TText2>
  void do_fold(TText1 const &seq, TText2 &results);

  template <typename TText1, typename TText2, typename TText3>
  void do_fold(TText1 const &seq, TText2 const &quals, TText3 &results);

  template <typename TText1, typename TText2>
  int paired_tfo(TText1 sequence, TText2 structure);

  template <typename TText1, typename TText2, typename TText3>
  Dna5String refilter_records(TText1 const &structure,
                              TText2 const &pairs_lim_tfo, TText3 const &seq);
  template <typename TText1, typename TText2, typename TText3>
  void search_limits(TText1 const &pair_limits, TText2 const &local_record,
                     TText3 &found_tfo);
  template <typename TText1, typename TText2, typename TText3>
  void search_limits(TText1 const &pair_limits, TText2 const &local_record,
                     TText3 &found_tfo, Tag_ds const &);
  template <typename TText1, typename TText2, typename TText3>
  void search_limits_revc(TText1 const &pair_limits, TText2 const &local_record,
                          TText3 &found_tfo, Tag_ds const &);
};

// Default constructor of LncRna() does nothing
LncRna::LncRna() {}
LncRna::LncRna(Dna5String haystack, CharString id, int index) {
  found_tfos = false;
  index_id = index;
  id_lncrna = id;
  struct_computed = false;
  find_tfo(haystack, id, index);
}
LncRna::LncRna(Dna5String haystack, CharString id, int index,
               Tag_folding const & /*Tag*/) {
  found_tfos = false;
  index_id = index;
  id_lncrna = id;
  struct_computed = false;
  find_tfo(haystack, id, index, Tag_folding());
}
LncRna::LncRna(Dna5String haystack, CharString id, int index,
               Tag_ds const & /*Tag*/) {
  found_tfos = false;
  index_id = index;
  id_lncrna = id;
  struct_computed = false;
  find_tfo(haystack, id, index, Tag_ds());
}
LncRna::LncRna(Dna5String haystack, CharString id, CharString quals, int index,
               Tag_constraints const & /*Tag*/) {
  found_tfos = false;
  index_id = index;
  id_lncrna = id;
  struct_computed = false;
  find_tfo(haystack, id, quals, index, Tag_constraints());
}
// LncRna::~LncRna() {} // destructor left for the compiler (!!)

///////////////////////////////////////////////////////////////////////////////
/// Possible TODO ==> write a common interface for finding patterns,
///                   pass the pattern in the function call and that you can
///                   can get rid two redundant functions
/// Note that to makes things clear, I have a different mechanism to call a
/// two functions that do the same (search for the same pattern). This is done
/// like that to facilitate things.
///////////////////////////////////////////////////////////////////////////////

template <typename TText1>
vector<tuple<int, int>> LncRna::search_by_seq(TText1 const &sequence,
                                              boost::regex const &re) {
  // keep the index of possible tfos in a tuple list
  vector<tuple<int, int>> putative_tfo_limits;

  // 1) not sure if ther is a better way to get a string from the sequence
  // record, at the moment I write it to a stream and from there to a string
  stringstream support_stream;
  support_stream << sequence;
  string str_sequence = support_stream.str();
  try {
    // notice that the pattern re is defined in the calling function
    // search, with the regular expression
    boost::sregex_iterator next(str_sequence.begin(), str_sequence.end(), re);
    boost::sregex_iterator end; // on default constructor to check against end
    while (next != end) {
      boost::smatch match = *next; // dereference match and write it as string
                                   //    cout << match.str() << endl;
      int init = match.position();
      int endit = match.position() + match.length();
      putative_tfo_limits.emplace_back(init, endit);
      next++;
    }
  } catch (boost::regex_error &e) {
    // Syntax error in the regular expression
    cout << "Something went wrong setting up the regular_expression" << endl;
  }
  return putative_tfo_limits;
}
///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2, typename TText3>
void LncRna::search_limits(TText1 const &pair_limits,
                           TText2 const &local_record, TText3 &found_tfo) {

  int tfo_init = get<0>(pair_limits);
  int tfo_end = get<1>(pair_limits);
  found_tfo.position = make_tuple(tfo_init, tfo_end);
  found_tfo.seq = (infix(local_record, tfo_init, tfo_end));
  auto results_tm = compute_tfo_tm(found_tfo.seq, tfo_conc, pH);
  found_tfo.thermo.dg = results_tm.first;
  found_tfo.thermo.tm = results_tm.second;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2, typename TText3>
void LncRna::search_limits(TText1 const &pair_limits,
                           TText2 const &local_record, TText3 &found_tfo,
                           Tag_ds const & /*Tag*/) {

  int tfo_init = get<0>(pair_limits);
  int tfo_end = get<1>(pair_limits);
  found_tfo.position = make_tuple(tfo_init, tfo_end);
  Dna5String local_seq = infix(local_record, tfo_init, tfo_end);
  Dna5StringComplement complement(local_seq);
  found_tfo.seq = complement;
  auto results_tm = compute_tfo_tm(found_tfo.seq, tfo_conc, pH);
  found_tfo.thermo.dg = results_tm.first;
  found_tfo.thermo.tm = results_tm.second;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2, typename TText3>
void LncRna::search_limits_revc(TText1 const &pair_limits,
                                TText2 const &local_record, TText3 &found_tfo,
                                Tag_ds const & /*Tag*/) {

  int tfo_init = get<0>(pair_limits);
  int tfo_end = get<1>(pair_limits);
  found_tfo.position = make_tuple(tfo_init, tfo_end);
  Dna5String local_seq = infix(local_record, tfo_init, tfo_end);
  Dna5StringReverse rev(local_seq);
  found_tfo.seq = rev;
  found_tfo.strand = "-";
  auto results_tm = compute_tfo_tm(found_tfo.seq, tfo_conc, pH);
  found_tfo.thermo.dg = results_tm.first;
  found_tfo.thermo.tm = results_tm.second;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2, typename TText3>
void LncRna::find_tfo(TText1 const &local_record, TText2 local_id_lncrna,
                      TText3 index) {

  Tfos found_tfo;
  static const boost::regex re("[U|T|C]{12,}");
  vector<tuple<int, int>> tfos_limits = this->search_by_seq(local_record, re);

  this->id_lncrna = local_id_lncrna;
  for (auto pairs_lim_tfo : tfos_limits) {
    search_limits(pairs_lim_tfo, local_record, found_tfo);
    if (found_tfo.thermo.tm > tm_threshold) {
      this->found_tfos = true;
      this->lnctfos.push_back(found_tfo);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2, typename TText3>
void LncRna::find_tfo(TText1 const &local_record, TText2 local_id_lncrna,
                      TText3 index, Tag_ds const & /*Tag*/) {

  // In this case we search and the loop for each of the strands
  // The goal is to report the TFO and the coords of site it binds to
  // If we find {A,G} in the fwd strand, we report its
  // complement (not reversed!) and annotate a strand="+" and the coordinates
  // If we find a {T,C} in the fwd strand, we report its reverse
  // and annotate a strand="-" along with the coordinates.
  // Please note that in this case search_limits searches for {A,G} patterns
  // and search_limits_revc searches for {T,C} patterns. I qualify the function
  // with a Tag_ds() for consistency, although it does not really need it
  // as it is the only case where the function is being called.

  Tfos found_tfo;
  static const boost::regex re("[A|G]{12,}");
  vector<tuple<int, int>> tfos_limits = this->search_by_seq(local_record, re);
  static const boost::regex re2("[T|C]{12,}");
  vector<tuple<int, int>> tfos_limits_revc =
      this->search_by_seq(local_record, re2);

  this->id_lncrna = local_id_lncrna;
  for (auto pairs_lim_tfo : tfos_limits) {
    search_limits(pairs_lim_tfo, local_record, found_tfo, Tag_ds());
    if (found_tfo.thermo.tm > tm_threshold) {
      this->found_tfos = true;
      this->lnctfos.push_back(found_tfo);
    }
  }
  for (auto pairs_lim_tfo : tfos_limits_revc) {
    search_limits_revc(pairs_lim_tfo, local_record, found_tfo, Tag_ds());
    if (found_tfo.thermo.tm > tm_threshold) {
      found_tfo.strand = "-";
      this->found_tfos = true;
      this->lnctfos.push_back(found_tfo);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2, typename TText3>
void LncRna::find_tfo(TText1 const &local_record, TText2 local_id_lncrna,
                      TText3 index, Tag_folding const & /*Tag*/) {

  Tfos found_tfo;
  // declared in do_rnafold.h
  res_rnafold fold;
  this->id_lncrna = local_id_lncrna;

  static const boost::regex re("[U|T|C]{12,}");
  vector<tuple<int, int>> tfos_limits = this->search_by_seq(local_record, re);

  for (auto pairs_lim_tfo : tfos_limits) {
    search_limits(pairs_lim_tfo, local_record, found_tfo);
    if (found_tfo.thermo.tm > tm_threshold) {
      // if the sequnce is too long it does not make sense
      // to compute the fold, we keep it for the moment
      if (length(local_record) < max_fold_rna) {
        // we avoid computing the folding more than once
        if (!this->struct_computed) {
          // Folding, remember to free() !!
          fold.structure =
              new char[(sizeof(char) * (length(local_record) + 1))];
          fold.fe = 0;
          // A big mess to change the sequence from "SeqAn" Alphabet
          // to a const char * needed in RNAfold
          stringstream support_stream;
          support_stream << local_record;
          const string &tmp = support_stream.str();
          const char *lncrna_seq = tmp.c_str();
          do_fold(lncrna_seq, fold);
          support_stream.str(std::string()); // clear it
          this->struct_computed = true;
        }
        // if we find a nice region we store it
        if (!paired_tfo(pairs_lim_tfo, fold.structure)) {
          this->found_tfos = true;
          this->lnctfos.push_back(found_tfo);
        } else {
          // search for smaller segments of unpaired regions
          Dna5String refilt_record = this->refilter_records(
              fold.structure, pairs_lim_tfo, local_record);
          vector<tuple<int, int>> refilt_tfos_limits =
              this->search_by_seq(refilt_record, re);
          for (auto pairs_lim_tfo : refilt_tfos_limits) {
            search_limits(pairs_lim_tfo, local_record, found_tfo);
            if (found_tfo.thermo.tm > tm_threshold) {
              this->found_tfos = true;
              this->lnctfos.push_back(found_tfo);
            }
          }
        }
      }
    } // else  we don't store it
  }
  // now we can free the structure if we have actually computed it
  if (this->struct_computed) {
    free(fold.structure);
  }
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2, typename TText3, typename TText4>
void LncRna::find_tfo(TText1 const &local_record, TText2 local_id_lncrna,
                      TText3 &quals, TText4 index,
                      Tag_constraints const & /*Tag*/) {

  Tfos found_tfo;
  // declared in do_rnafold.h
  res_rnafold fold;
  static const boost::regex re("[U|T|C]{12,}");
  this->id_lncrna = local_id_lncrna;
  vector<tuple<int, int>> tfos_limits = this->search_by_seq(local_record, re);

  for (auto pairs_lim_tfo : tfos_limits) {
    search_limits(pairs_lim_tfo, local_record, found_tfo);
    if (found_tfo.thermo.tm > tm_threshold) {
      // if the sequnce is too long it does not make sense
      // to compute the fold, we keep it for the moment
      if (length(local_record) < max_fold_rna) {
        // we avoid computing the folding more than once
        if (!this->struct_computed) {
          // Folding, remember to free() !!
          fold.structure =
              new char[(sizeof(char) * (length(local_record) + 1))];
          fold.fe = 0;

          // A big mess to change the sequence from "SeqAn" Alphabet
          // to a const char * needed in RNAfold
          stringstream support_stream;
          support_stream << local_record;
          const string &tmp = support_stream.str();
          const char *lncrna_seq = tmp.c_str();
          // Now change the notation from Pietro's to Rnafold (see comments
          // above on what's going on)
          ModifiedString<String<char>, ModView<Pietro2RnafoldFunctor>>
              pietro2rnafold_quals(quals);
          stringstream support_stream2;
          support_stream2 << pietro2rnafold_quals;
          const string &tmp2 = support_stream2.str();
          const char *lncrna_quals = tmp2.c_str();
          do_fold(lncrna_seq, lncrna_quals, fold);
          support_stream2.str(std::string()); // clear it
          support_stream.str(std::string());  // clear it

          this->struct_computed = true;
        }
        // if we find a nice region we store it
        if (!paired_tfo(pairs_lim_tfo, fold.structure)) {
          this->found_tfos = true;
          this->lnctfos.push_back(found_tfo);
        } else {
          // search for smaller segments of unpaired regions
          Dna5String refilt_record = this->refilter_records(
              fold.structure, pairs_lim_tfo, local_record);
          vector<tuple<int, int>> refilt_tfos_limits =
              this->search_by_seq(refilt_record, re);
          for (auto pairs_lim_tfo : refilt_tfos_limits) {
            search_limits(pairs_lim_tfo, local_record, found_tfo);
            if (found_tfo.thermo.tm > tm_threshold) {
              this->found_tfos = true;
              this->lnctfos.push_back(found_tfo);
            }
          }
        }
      }
      // ignore too large by now
      /*else if (length(local_record) > max_fold_rna) {
         this->found_tfos = true;
         this->lnctfos.push_back(found_tfo);
       }
       */
    } // else  we don't store it
  }
  // now we can free the structure if we have actually computed it
  if (this->struct_computed) {
    free(fold.structure);
  }
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2, typename TText3>
Dna5String LncRna::refilter_records(TText1 const &structure,
                                    TText2 const &pairs_lim_tfo,
                                    TText3 const &seq) {

  Dna5String filtered;
  assign(filtered, seq);
  int tfo_init = get<0>(pairs_lim_tfo);
  int tfo_end = get<1>(pairs_lim_tfo);
  // not using iterators
  // would be great to avoid using ifs, but I can't think of a faster way
  for (unsigned i = tfo_init; i <= tfo_end; ++i) {
    if (structure[i] == '.') {
      filtered[i] = seq[i];
    } else {
      filtered[i] = 'N';
    }
  }
  return filtered;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
int LncRna::paired_tfo(TText1 pairs_lim_tfo, TText2 structure) {
  // I have decided to do it like this to avoid more ifs
  // the program will only return false (int=0) if all the chars are '.'
  int tfo_init = get<0>(pairs_lim_tfo);
  int tfo_end = get<1>(pairs_lim_tfo);
  unsigned int i = 0;
  // the return value will be 0 (true) if all characters are '.'
  int return_value = (tfo_end - tfo_init) + 1;
  while (structure[tfo_init + i] != '\0' && structure[tfo_init + i] == '.' &&
         (tfo_init + i) <= tfo_end) {
    return_value--;
    ++i; // Never ever touch that or you will die in here
  }
  return return_value;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
void LncRna::do_fold(TText1 const &seq, TText2 &results) {
  do_rnafold_nc(seq, &results);
}

template <typename TText1, typename TText2, typename TText3>
void LncRna::do_fold(TText1 const &seq, TText2 const &quals, TText3 &results) {
  do_rnafold(seq, quals, &results);
}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// End of class LncRna and its associated functions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  Functions
///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
vector<LncRna> find_tfo_from_fasta(TText1 const &seqs, TText2 const &ids) {

  vector<LncRna> results;
  LncRna *hack = new LncRna[length(ids)];
// this hack is to get an array of objects such that we can have random
// access. I still wanted to return the vector of objects, so I have to
// copy the solutions by hand, a bit ugly but not a big deal
#pragma omp parallel for
  for (unsigned i = 0; i < length(ids); ++i) {
    hack[i] = LncRna(seqs[i], ids[i], i);
    // results.push_back(LncRna(seqs[i], ids[i], i));
  }
  for (unsigned i = 0; i < length(ids); ++i) {
    results.push_back(hack[i]);
  }
  delete[] hack;
  return results;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
vector<LncRna> find_tfo_from_fasta(TText1 const &seqs, TText2 const &ids,
                                   Tag_folding const & /*Tag*/) {

  vector<LncRna> results;
  LncRna *hack = new LncRna[length(ids)];
// this hack is to get an array of objects such that we can have random
// access. I still wanted to return the vector of objects, so I have to
// copy the solutions by hand, a bit ugly but not a big deal
#pragma omp parallel for
  for (unsigned i = 0; i < length(ids); ++i) {
    hack[i] = LncRna(seqs[i], ids[i], i, Tag_folding());
    // results.push_back(LncRna(seqs[i], ids[i], i, Tag_folding()));
  }
  for (unsigned i = 0; i < length(ids); ++i) {
    results.push_back(hack[i]);
  }
  delete[] hack;
  return results;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
vector<LncRna> find_tfo_from_fasta(TText1 const &seqs, TText2 const &ids,
                                   Tag_ds const & /*Tag*/) {

  vector<LncRna> results;
  LncRna *hack = new LncRna[length(ids)];
// this hack is to get an array of objects such that we can have random
// access. I still wanted to return the vector of objects, so I have to
// copy the solutions by hand, a bit ugly but not a big deal
#pragma omp parallel for
  for (unsigned i = 0; i < length(ids); ++i) {
    hack[i] = LncRna(seqs[i], ids[i], i, Tag_ds());
    // results.push_back(LncRna(seqs[i], ids[i], i, Tag_folding()));
  }
  for (unsigned i = 0; i < length(ids); ++i) {
    results.push_back(hack[i]);
  }
  delete[] hack;
  return results;
}

///////////////////////////////////////////////////////////////////////////////
template <typename TText1, typename TText2, typename TText3>
vector<LncRna> find_tfo_from_fasta(TText1 const &seqs, TText2 const &ids,
                                   TText3 &quals,
                                   Tag_constraints const & /*Tag*/) {

  vector<LncRna> results;
  LncRna *hack = new LncRna[length(ids)];
// this hack is to get an array of objects such that we can have random
// access. I still wanted to return the vector of objects, so I have to
// copy the solutions by hand, a bit ugly but not a big deal
#pragma omp parallel for
  for (unsigned i = 0; i < length(ids); ++i) {
    hack[i] = LncRna(seqs[i], ids[i], quals[i], i, Tag_constraints());
    //    results.push_back(LncRna(seqs[i], ids[i], quals[i], i,
    //    Tag_constraints()));
  }
  for (unsigned i = 0; i < length(ids); ++i) {
    results.push_back(hack[i]);
  }
  delete[] hack;
  return results;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
CharString compose_id_string(TText1 &it, TText2 const &it2) {
  CharString id_info = "";
  ostringstream sstrs;
  id_info += it->id_lncrna;
  id_info += "|sense|";
  sstrs << setprecision(4) << it2->thermo.dg;
  string support_str = sstrs.str();
  sstrs.str(string()); // clear the contents
  id_info += support_str;
  id_info += "|";
  sstrs << setprecision(4) << it2->thermo.tm;
  support_str = sstrs.str();
  sstrs.str(string()); // clear the contents
  id_info += support_str;
  return id_info;
}

///////////////////////////////////////////////////////////////////////////////
template <typename TText1, typename TText2, typename TText3>
CharString compose_id_string(TText1 &it, TText2 const &it2, TText3 &counter,
                             Tag_ds const & /*Tag*/) {
  CharString id_info = "";
  ostringstream sstrs;
  id_info += it->id_lncrna;
  id_info += "|";
  sstrs << get<0>(it2->position) << ":" << get<1>(it2->position);
  string support_str = sstrs.str();
  id_info += support_str;
  sstrs.str(string()); // clear the contents
  id_info += "|";
  id_info += it2->strand;
  id_info += "|";
  sstrs << setprecision(4) << it2->thermo.dg;
  support_str = sstrs.str();
  sstrs.str(string()); // clear the contents
  id_info += support_str;
  id_info += "|";
  sstrs << setprecision(4) << it2->thermo.tm;
  support_str = sstrs.str();
  sstrs.str(string()); // clear the contents
  id_info += support_str;
  id_info += "|";
  id_info += "Id=";
  sstrs << counter;
  support_str = sstrs.str();
  sstrs.str(string()); // clear the contents
  id_info += support_str;

  return id_info;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
CharString compose_id_string(TText1 &it, TText2 const &it2,
                             Tag_complement const &) {
  CharString id_info = "";
  ostringstream sstrs;
  id_info += it->id_lncrna;
  id_info += "|antisense|";
  sstrs << setprecision(4) << it2->thermo.dg;
  string support_str = sstrs.str();
  sstrs.str(string()); // clear the contents
  id_info += support_str;
  id_info += "|";
  sstrs << setprecision(4) << it2->thermo.tm;
  support_str = sstrs.str();
  sstrs.str(string()); // clear the contents
  id_info += support_str;
  return id_info;
}

///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
int print_putative_tfos_to_fasta(TText1 &results, TText2 const &outname,
                                 Tag_ds const & /*Tag*/) {

  SeqFileOut seqFileOut;
  StringSet<CharString> ids;
  StringSet<Dna5String> seqs;

  // In this case, we will have separated already fwd and revcomplement
  // strands in the finding routines, so just go ahead and print
  // (this is as opposed to the un-tagged one, without Tag_ds)
  vector<LncRna>::iterator it = results.begin();
  int found_hits = 0;
  for (; it != results.end(); ++it) {
    if (it->found_tfos) {
      vector<LncRna::Tfos> local_tfo = it->lnctfos;
      vector<LncRna::Tfos>::iterator it2 = local_tfo.begin();
      for (; it2 != local_tfo.end(); ++it2) {
        CharString id_info = compose_id_string(it, it2, found_hits, Tag_ds());
        appendValue(ids, id_info);
        Dna5String out_transcribed = it2->seq;
        appendValue(seqs, out_transcribed);
        found_hits++;
      }
    }
  }
  // open sequence file
  if (found_hits > 0) {
    if (!open(seqFileOut, toCString(outname))) {
      cerr << "ERROR: Cound not open output file (sequences).\n";
      return 1;
    }
    try {
      writeRecords(seqFileOut, ids, seqs);
    } catch (Exception const &e) {
      cout << "ERROR: " << e.what() << endl;
      return 1;
    }
  } else {
    cout << "Not writting records, nothing was found" << endl;
  }
  return 0;
}
///////////////////////////////////////////////////////////////////////////////

template <typename TText1, typename TText2>
int print_putative_tfos_to_fasta(TText1 &results, TText2 const &outname) {

  SeqFileOut seqFileOut;
  StringSet<CharString> ids;
  // We will write the transcribed DNA, as we want to find
  // the binding partner. Luckily, the rules of TFO formation for this
  // type of triplex are the same as WC, e.g. C with G and T with A.
  // In casting the Rna5String as Dna5String, seqan automatically
  // translates U into T, so we don't have to transcribe it
  // In any case, I have written nn_tm to deal with Dna5Strings, instead
  // of Rna because the data I'm anlysing is actually as genes.
  StringSet<Dna5String> seqs;

  vector<LncRna>::iterator it = results.begin();
  int found_hits = 0;
  for (; it != results.end(); ++it) {
    if (it->found_tfos) {
      vector<LncRna::Tfos> local_tfo = it->lnctfos;
      vector<LncRna::Tfos>::iterator it2 = local_tfo.begin();
      for (; it2 != local_tfo.end(); ++it2) {
        CharString id_info = compose_id_string(it, it2);
        appendValue(ids, id_info);
        Dna5String out_transcribed = it2->seq;
        // This is actually a modifier pre-defined by Seqan
        // This works because the pairing rules for our TFOs are the same
        // as for WC
        Dna5StringComplement transc_complement(out_transcribed);
        appendValue(seqs, transc_complement);
        // now we also add the complement of the complement (the genome
        // might only report on the sense and we might want target
        // the antisense). This is, we simply write the out_transcribed/
        // and we are done
        CharString id_info_complement =
            compose_id_string(it, it2, Tag_complement());
        appendValue(ids, id_info_complement);
        appendValue(seqs, out_transcribed);
        found_hits++;
      }
    }
  }
  // open sequence file
  if (found_hits > 0) {
    if (!open(seqFileOut, toCString(outname))) {
      cerr << "ERROR: Cound not open output file (sequences).\n";
      return 1;
    }
    try {
      writeRecords(seqFileOut, ids, seqs);
    } catch (Exception const &e) {
      cout << "ERROR: " << e.what() << endl;
      return 1;
    }
  } else {
    cout << "Not writting records, nothing was found" << endl;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// The end
///////////////////////////////////////////////////////////////////////////////
