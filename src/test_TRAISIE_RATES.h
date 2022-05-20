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
//' @param trans_rate tr1
//' @param trans_rate2 tr2
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


//' function to test get_trans_rate
//' @param timeval timeval
//' @param total_time total time
//' @param gam gam
//' @param laa laa
//' @param lac lac
//' @param mu mu
//' @param K K
//' @param num_spec num_spec
//' @param num_immigrants num_immigrants
//' @param mainland_n mainland_n
//' @param island_spec island_spec
//' @param trait_pars trait pars
//' @return two rates
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector test_update_rates(double timeval,
                                      double total_time,
                                      double gam,
                                      double laa,
                                      double lac,
                                      double mu,
                                      double K,
                                      double num_spec,
                                      double num_immigrants,
                                      double mainland_n,
                                      const Rcpp::NumericMatrix& island_spec_R,
                                      const Rcpp::List& trait_pars_R) {

  island_spec island_spec_;
  for (size_t i = 0; i < island_spec_R.nrow(); ++i) {
    double colonist = island_spec_R(i, 0);
    double timeval  = island_spec_R(i, 1);
    double trait_val = island_spec_R(i, 7);

    island_spec_row to_add(colonist, timeval, species_type::I, trait_val);
    island_spec_.push_back(to_add);
  }

  rates trait_pars(gam, laa, lac, mu, trait_pars_R);

  auto answer = update_rates(timeval,
                             total_time,
                             gam,
                             laa,
                             lac,
                             mu,
                             K,
                             num_spec,
                             num_immigrants,
                             mainland_n,
                             island_spec_,
                             trait_pars);

  Rcpp::NumericVector out = {answer.immig_rate,
                             answer.ext_rate,
                             answer.ana_rate,
                             answer.clado_rate,
                             answer.trans_rate,
                             answer.immig_rate2,
                             answer.ext_rate2,
                             answer.ana_rate2,
                             answer.clado_rate2,
                             answer.trans_rate2,
                             answer.M2};
  return out;
}








