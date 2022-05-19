#ifndef DAISIE_RATES_H
#define DAISIE_RATES_H

#include <Rcpp.h>
#include <vector>
#include <array>
#include <string>

#include "TRAISIE_island_spec.h"

struct area_pars { // not really used though
  double total_island_age;
  double current_area;
  double max_area;
  double proportional_peak_t;
  double sea_level_amplitude;
  double sea_level_frequency;
  double island_gradient_angle;

  area_pars(const Rcpp::List& from_R) {
    total_island_age = from_R["total_island_age"];
    current_area = from_R["current_area"];
    max_area = from_R["max_area"];
    proportional_peak_t = from_R["proportional_peak_t"];
    sea_level_amplitude = from_R["sea_level_amplitude"];
    sea_level_frequency = from_R["sea_level_frequency"];
    island_gradient_angle = from_R["island_gradient_angle"];
  }
};

struct two_rates {
  double rate1;
  double rate2;
};

struct rates {
  double immig_rate;
  double ext_rate;
  double ana_rate;
  double clado_rate;
  double trans_rate;
  double immig_rate2;
  double ext_rate2;
  double ana_rate2;
  double clado_rate2;
  double trans_rate2;
  double M2;

  double sum() const {
    return immig_rate + immig_rate2 +
           ext_rate + ext_rate2 +
           ana_rate + ana_rate2 +
           clado_rate + clado_rate2 +
           trans_rate + trans_rate2;
  }

  int sample_event() {

  }

  rates(two_rates immig,
        two_rates ext,
        two_rates ana,
        two_rates clado,
        two_rates trans,
        double m2) :
    immig_rate(immig.rate1),
    immig_rate2(immig.rate2),
    ext_rate(ext.rate1),
    ext_rate2(ext.rate2),
    ana_rate(ana.rate1),
    ana_rate2(ana.rate2),
    clado_rate(clado.rate1),
    clado_rate2(clado.rate2),
    trans_rate(trans.rate1),
    trans_rate2(trans.rate2),
    M2(m2)
  {
  }

  rates(double gam,
        double laa,
        double lac,
        double mu,
        const Rcpp::List& trait_pars_from_R) {
    immig_rate = gam;
    immig_rate2 = trait_pars_from_R["immig_rate2"];
    ext_rate =  mu;
    ext_rate2 = trait_pars_from_R["ext_rate2"];
    ana_rate = laa;
    ana_rate2 = trait_pars_from_R["ana_rate2"];
    clado_rate = lac;
    clado_rate2 = trait_pars_from_R["clado_rate2"];
    trans_rate =  trait_pars_from_R["trans_rate"];
    trans_rate2 = trait_pars_from_R["trans_rate2"];
    M2 = trait_pars_from_R["M2"];
  }
};

//' function to draw event.
//' @param event_prob
//' @return event
//' @export
// [[Rcpp::export]]
int sample_event(const std::vector<double>& event_prob) {
/*  std::vector<double> event_prob = {immig_rate, ext_rate,
                                    ana_rate, clado_rate,
                                    trans_rate, immig_rate2,
                                    ext_rate2, ana_rate2,
                                    clado_rate2,
                                    trans_rate2};
  */
  double s = std::accumulate(event_prob.begin(), event_prob.end(), 0.0);
  double r = R::runif(0.0, s);
  int index = 0;

  for( ; index < event_prob.size(); ++index) {
    r -= event_prob[index];
    if (r <= 0.0) {
      break;
    }
  }
  return index + 1;
}

two_rates get_immig_rate(double gam,
                         double A,
                         int num_spec,
                         double K,
                         double mainland_n,
                         double mainland_n2,
                         double gam2) {

 // auto mainland_n2 = trait_pars.M2;
//  auto gam2 = trait_pars.immig_rate2;

  two_rates immig_rate;
  immig_rate.rate1 = std::max(mainland_n * gam * (1.0 - (num_spec / (A * K))),
                              0.0);
  immig_rate.rate2 = std::max(mainland_n2 * gam2 * (1.0 - (num_spec / (A * K))),
                              0.0);

  return immig_rate;
}

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




two_rates get_ext_rate(double mu,
                       int num_spec,
                       double A,
                       const rates& trait_pars,
                       size_t num_spec_trait1,
                       size_t num_spec_trait2) {
  two_rates ext_rate;
  ext_rate.rate1 = mu * num_spec_trait1;
  ext_rate.rate2 = trait_pars.ext_rate2 * num_spec_trait2;
  return ext_rate;
}

two_rates get_ana_rate(double laa,
                       int num_immigrants,
                       const rates& trait_pars,
                       const island_spec& island_spec_) {

  size_t num_I_1 = 0;
  size_t num_I_2 = 0;

  for (size_t i = 0; i < island_spec_.data_.size(); ++i) {
    auto s1 = island_spec_.data_[i].type_species;
    auto s2 = island_spec_.data_[i].trait;

    if (s1 == species_type::I && s2 == 1) num_I_1++;
    if (s1 == species_type::I && s2 == 2) num_I_2++;
  }

  two_rates ana_rate;
  ana_rate.rate1 = laa * num_I_1;
  ana_rate.rate2 = trait_pars.ana_rate2* num_I_2;
  return ana_rate;
}

two_rates get_clado_rate(double lac,
                         double num_spec,
                         double K,
                         double A,
                         const rates& trait_pars,
                         size_t num_spec_trait1,
                         size_t num_spec_trait2) {

  two_rates clado_rate;
  clado_rate.rate1 = std::max(0.0,
                              lac * num_spec_trait1 * (1.0 - num_spec / K) );

  clado_rate.rate2 = std::max(0.0,
                              trait_pars.clado_rate2 * num_spec_trait2 * (1.0 - num_spec / K));
  return clado_rate;
}

two_rates get_trans_rate(const rates& trait_pars,
                         size_t num_spec_trait1,
                         size_t num_spec_trait2) {
  two_rates trans_rate;

  trans_rate.rate1 = trait_pars.trans_rate * num_spec_trait1;
  trans_rate.rate2 = trait_pars.trans_rate2 * num_spec_trait2;

  return trans_rate;
}


rates update_rates(double timeval,
                   double total_time,
                   double gam,
                   double laa,
                   double lac,
                   double mu,
                   const area_pars& ap,
                   double K,
                   double num_spec,
                   double num_immigrants,
                   double mainland_n,
                   const island_spec& island_spec_,
                   const rates& trait_pars) {

  double A = 1.0; // island ontogeny dynamics are currently not needed.

  size_t num_spec_trait1 = 0;
  size_t num_spec_trait2 = 0;
  if (!island_spec_.data_.empty()) {
   for (const auto& i : island_spec_.data_) {
     if (i.trait == 1) num_spec_trait1++;
     if (i.trait == 2) num_spec_trait2++;
   }
  }

  two_rates immig_rate = get_immig_rate(gam, A, num_spec, K, mainland_n,
                                        trait_pars.M2,
                                        trait_pars.immig_rate2);

  two_rates ext_rate = get_ext_rate(mu, num_spec, A,
                                    trait_pars,
                                    num_spec_trait1,
                                    num_spec_trait2);

  two_rates ana_rate = get_ana_rate(laa, num_immigrants, trait_pars, island_spec_);
  two_rates clado_rate = get_clado_rate(lac, num_spec, K, A,
                                        trait_pars,
                                        num_spec_trait1,
                                        num_spec_trait2);
  two_rates trans_rate = get_trans_rate(trait_pars,
                                        num_spec_trait1,
                                        num_spec_trait2);


  rates r(immig_rate, ext_rate, ana_rate, clado_rate, trans_rate, trait_pars.M2);

  return(r);
}


#endif /* DAISIE_RATES_H */
