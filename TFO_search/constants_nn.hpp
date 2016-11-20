#include <vector>
using namespace std;

static const double R = 0.0019858;

// From my parameterization
vector<const double> rna_dh = {-10.95, -5.73, -6.44};
vector<const double> rna_dg = {-1.891, -0.758, -0.331, 2.646};
vector<const double> rna_ph = {0.893, -0.005};
// From Roberts & Crothers  PNAS 1996, vol 93 pp 4320
vector<const double> dna_dh = {-4.99, -8.9, -7.4};
vector<const double> dna_dg = {-3.00, -0.65, -1.65, 6.0};
vector<const double> dna_ph = {1.26, -0.08};

typedef struct my_nn_rna_mode {
  vector<const double> dh_coeff;
  vector<const double> dg_coeff;
  vector<const double> pH_coeff;
} nn_model_params_t;

nn_model_params_t nn_rna_params = {rna_dh, rna_dg, rna_ph};
nn_model_params_t nn_dna_params = {dna_dh, dna_dg, dna_ph};

// This type of initializations is in principle also possible but
// I think the other one is clearer
//
// my_rna_params_t nn_rna_params = {
//    {-10.95, -5.73, -6.44},          /* DH  */
//    {-1.891, -0.758, -0.331, 2.646}, /* DG */
//    {0.893, -0.005}                  /* pH */
// };
