#include "constants_nn.hpp"
#include <functional>
#include <iostream>
#include <math.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <typeinfo>
#include <utility>
#include <vector>
using namespace seqan;
using namespace std;
/////////////////////////////////////////////////////////////////////////////
//
// !!!WARNING!!!! ALTHOUGH THE MODEL IS FOR RNA, I DECLARE AND
// THE SEQUENCES AS DNA5STRINGS, AND SEARCH FOR CT/TC ETC...
// This way I can also use the DNA model from Crothers et al without
// specially dedicated functions
// This is because the data in lncipedia is of genes, and hence in dna, and
// transcribing it takes time (specially if the input is large).
// Be aware of this.
//
// In any case, I have implemented the bool b_read_rna in case you read RNA,
// which basically transcribes the RNA into DNA, as this should be less
// frequent, in case you do something else.
//
/////////////////////////////////////////////////////////////////////////////

template <typename TText> vector<double> nncount(TText &haystack) {
  vector<double> nn(3);
  int counter = 0;
  bool bVerbose = false;
  // See if I can make them constants as well!!
  // trouble is the odd thing you do to add CT/TC
  // might end up unexplained
  String<Dna5String> needles_dh;
  // T vs U is cosmetic here, as Dna5String converts U to C anyhow
  // do not change this order or else grouping base on TC/CT won't work
  appendValue(needles_dh, "CC");
  appendValue(needles_dh, "CT");
  appendValue(needles_dh, "TC");
  appendValue(needles_dh, "TT");

  // you could do it with an iterator over needles, but you'd need
  // to keep a counter to fiddle with nn, or use it - vec.begin()
  // according to stackoverflow (http://goo.gl/2CFXw8)
  for (unsigned i = 0; i < (nn.size()) + 1; ++i) {
    // Finder has to be reset each iteration, or else
    // it will continue searching after the last match
    Finder<Dna5String> finder(haystack);
    Pattern<Dna5String, Horspool> pattern(needles_dh[i]);
    if (bVerbose) {
      cout << needles_dh[i] << '\t';
    }
    while (find(finder, pattern)) {
      counter++;
      if (bVerbose) {
        cout << '[' << beginPosition(finder) << ',';
        cout << endPosition(finder) << ") ";
      }
    }
    if (bVerbose) {
      cout << endl;
    }
    // we add CT and TC together
    if (i < 2) {
      nn[i] += counter;
    } else {
      nn[i - 1] += counter;
    }
    counter = 0;
  }
  return nn;
}

struct Tag_nn_dg {};
template <typename TText>
vector<double> nncount(TText &haystack, Tag_nn_dg const & /*Tag*/) {
  vector<double> nn(4);
  for (unsigned j = 0; j < nn.size(); ++j) {
    nn[j] = 0;
  }
  int counter = 0;
  bool bVerbose = false;
  String<Dna5String> needles_dg;
  // do not change this order or else grouping base on TC/CT won't work
  appendValue(needles_dg, "C");
  appendValue(needles_dg, "T");
  appendValue(needles_dg, "CC");

  for (unsigned i = 0; i < nn.size() - 1; ++i) {
    // Finder has to be reset each iteration, or else
    // it will continue searching after the last match
    Finder<Dna5String> finder(haystack);
    Pattern<Dna5String, Horspool> pattern(needles_dg[i]);
    if (bVerbose) {
      cout << needles_dg[i] << '\t';
    }
    while (find(finder, pattern)) {
      counter++;
      if (bVerbose) {
        cout << '[' << beginPosition(finder) << ',' << endPosition(finder)
             << ") ";
      }
    }
    if (bVerbose) {
      cout << endl;
    }
    nn[i] += counter;
    counter = 0;
  }
  nn[nn.size() - 1] = 1.0;
  return nn;
}

struct Tag_dna_nn {};
template <typename TText1> double compute_dh(TText1 &nn_vect) {
  return inner_product(nn_vect.begin(), nn_vect.end(),
                       nn_rna_params.dh_coeff.begin(), 0.0);
}

template <typename TText1>
double compute_dh(TText1 &nn_vect, Tag_dna_nn const &) {
  return inner_product(nn_vect.begin(), nn_vect.end(),
                       nn_dna_params.dh_coeff.begin(), 0.0);
}
template <typename TText1, typename TText2>
double compute_dg(TText1 &nn_vect, TText2 &pH) {
  double dg = inner_product(nn_vect.begin(), nn_vect.end(),
                            nn_rna_params.dg_coeff.begin(), 0.0);
  return dg +
         nn_vect[0] * (pH - 5.6) * (nn_rna_params.pH_coeff[0] +
                                    nn_rna_params.pH_coeff[1] * nn_vect[1]);
}

template <typename TText1, typename TText2>
double compute_dg(TText1 &nn_vect, TText2 &pH, Tag_dna_nn const & /*Tag*/) {
  double dg = inner_product(nn_vect.begin(), nn_vect.end(),
                            nn_dna_params.dg_coeff.begin(), 0.0);
  return dg +
         nn_vect[0] * (pH - 5.6) * (nn_dna_params.pH_coeff[0] +
                                    nn_dna_params.pH_coeff[1] * nn_vect[1]);
}
template <typename TText1, typename Tconc, typename TpH>
pair<double, double> compute_tfo_tm(TText1 &genomeFragment, Tconc const &conc,
                                    TpH const &pH) {
  // TODO --> Apparently genomeFragment can no be passed as const because
  // it is used by Finder later on, and he does not want it to be const
  // perhaps I could do a local copy ... and prevent nasty things happenning
  // to my string . <-- FIXME??
  vector<double> nn_dh_vect = nncount(genomeFragment);
  vector<double> nn_dg_vect = nncount(genomeFragment, Tag_nn_dg());
  double dh = compute_dh(nn_dh_vect);
  double dg = compute_dg(nn_dg_vect, pH);
  // prediction of Tm base on dh and dg, in C

  pair<double, double> dg_tm_pair;
  dg_tm_pair.first = dg;
  dg_tm_pair.second =
      ((298.0 * dh) / (dh - dg - 298.0 * R * log(4.0 / conc / 1e-6))) - 273.15;
  return dg_tm_pair;
}

template <typename TText1, typename Tconc, typename TpH>
pair<double, double> compute_tfo_tm(TText1 &genomeFragment, Tconc const &conc,
                                    TpH const &pH, Tag_dna_nn const & /*Tag*/) {
  // TODO --> Apparently genomeFragment can no be passed as const because
  // it is used by Finder later on, and he does not want it to be const
  // perhaps I could do a local copy ... and prevent nasty things happenning
  // to my string . <-- FIXME??
  vector<double> nn_dh_vect = nncount(genomeFragment);
  vector<double> nn_dg_vect = nncount(genomeFragment, Tag_nn_dg());
  double dh = compute_dh(nn_dh_vect, Tag_dna_nn());
  double dg = compute_dg(nn_dg_vect, pH, Tag_dna_nn());
  // prediction of Tm base on dh and dg, in C

  pair<double, double> dg_tm_pair;
  dg_tm_pair.first = dg;
  dg_tm_pair.second =
      ((298.0 * dh) / (dh - dg - 298.0 * R * log(4.0 / conc / 1e-6))) - 273.15;
  return dg_tm_pair;
}
