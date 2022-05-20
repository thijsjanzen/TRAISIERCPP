#include <Rcpp.h>

#include "TRAISIE_rates.h"

//' function to test get_immigration_rate
//' @param gam gam
//' @param A A
//' @param num_spec num species
//' @param K K value
//' @param mainland_n num mainland species trait 1
//' @param mainland_n2 num mainland species trait 2
//' @param immig_rate2 gam2
//' @return two rates
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector test_get_immig_rate(double gam,
                                        double A,
                                        int num_spec,
                                        double K,
                                        double mainland_n,
                                        double mainland_n2,
                                        double immig_rate2) {
  auto answer = get_immig_rate(gam, A, num_spec, K, mainland_n,
                               mainland_n2, immig_rate2);
  Rcpp::NumericVector out = {answer.rate1, answer.rate2};
  return out;
}

//' function to test get_ext_rate
//' @param mu mu
//' @param num_spec num species
//' @param A A
//' @param ext_rate2 ext_rate2
//' @param num_spec_trait1 num_spec trait 1
//' @param num_spec_trait2 num_spec trait 2
//' @return two rates
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector test_get_ext_rate(double mu,
                                      int num_spec,
                                      double A,
                                      double ext_rate2,
                                      int num_spec_trait1,
                                      int num_spec_trait2) {
  auto answer = get_ext_rate(mu, num_spec,
                             A, ext_rate2,
                             num_spec_trait1, num_spec_trait2);
  Rcpp::NumericVector out = {answer.rate1, answer.rate2};
  return out;
}

//' function to test get_ana_rate
//' @param laa anagenesis rate trait 1
//' @param num_immigrants number of immigrants
//' @param ana_rate2 ana_rate trait 2
//' @param num_spec_trait1 num_spec trait 1
//' @param num_spec_trait2 num_spec trait 2
//' @return two rates
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector test_get_ana_rate(double laa,
                                      int num_immigrants,
                                      double ana_rate2,
                                      size_t num_spec_trait1,
                                      size_t num_spec_trait2) {

  auto answer = get_ana_rate(laa, num_immigrants, ana_rate2,
                             num_spec_trait1, num_spec_trait2);
  Rcpp::NumericVector out = {answer.rate1, answer.rate2};
  return out;
}

//' function to test_get_clado_rate
//' @param lac clado rate trait 1
//' @param num_spec number of species
//' @param K K
//' @param A A
//' @param clado_rate2 clado rate trait 2
//' @param num_spec_trait1 num_spec trait 1
//' @param num_spec_trait2 num_spec trait 2
//' @return two rates
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector test_get_clado_rate(double lac,
                                        double num_spec,
                                        double K,
                                        double A,
                                        double clado_rate2,
                                        size_t num_spec_trait1,
                                        size_t num_spec_trait2) {

  auto answer = get_clado_rate(lac, num_spec, K, A, clado_rate2,
                               num_spec_trait1, num_spec_trait2);
  Rcpp::NumericVector out = {answer.rate1, answer.rate2};
  return out;
}

//' function to test get_trans_rate
//' @param trans_rate
//' @param trans_rate2
//' @param num_spec_trait1 num_spec trait 1
//' @param num_spec_trait2 num_spec trait 2
//' @return two rates
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector test_get_trans_rate(double trans_rate,
                                        double trans_rate2,
                                        size_t num_spec_trait1,
                                        size_t num_spec_trait2) {

  auto answer = get_trans_rate(trans_rate, trans_rate2,
                               num_spec_trait1, num_spec_trait2);
  Rcpp::NumericVector out = {answer.rate1, answer.rate2};
  return out;
}

