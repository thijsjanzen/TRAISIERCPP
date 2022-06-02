#ifndef TRAISIE_RATES_H
#define TRAISIE_RATES_H

#include <vector>
#include <array>
#include <string>

#include "TRAISIE_structs.h"

two_rates get_immig_rate(double gam,
                         double A,
                         int num_spec,
                         double K,
                         double mainland_n,
                         double mainland_n2,
                         double gam2) {

  two_rates immig_rate;
  immig_rate.rate1 = std::max(mainland_n  * gam  * (1.0 - (num_spec / (A * K))),
                              0.0);
  immig_rate.rate2 = std::max(mainland_n2 * gam2 * (1.0 - (num_spec / (A * K))),
                              0.0);

  return immig_rate;
}


two_rates get_ext_rate(double mu,
                       int num_spec,
                       double A,
                       double ext_rate2,
                       size_t num_spec_trait1,
                       size_t num_spec_trait2) {
  two_rates ext_rate;
  ext_rate.rate1 = mu        * num_spec_trait1;
  ext_rate.rate2 = ext_rate2 * num_spec_trait2;
  return ext_rate;
}

two_rates get_ana_rate(double laa,
                       int num_immigrants,
                       double ana_rate2,
                       size_t num_immig_trait1,
                       size_t num_immig_trait2) {
  two_rates ana_rate;
  ana_rate.rate1 = laa       * num_immig_trait1;
  ana_rate.rate2 = ana_rate2 * num_immig_trait2;
  return ana_rate;
}

two_rates get_clado_rate(double lac,
                         double num_spec,
                         double K,
                         double A,
                         double clado_rate2,
                         size_t num_spec_trait1,
                         size_t num_spec_trait2) {

  two_rates clado_rate;
  clado_rate.rate1 = std::max(0.0,
                              lac * num_spec_trait1 *
                              (1.0 - num_spec / K) );

  clado_rate.rate2 = std::max(0.0,
                              clado_rate2 * num_spec_trait2 *
                              (1.0 - num_spec / K));
  return clado_rate;
}

two_rates get_trans_rate(double trans_rate1,
                         double trans_rate2,
                         size_t num_spec_trait1,
                         size_t num_spec_trait2) {
  two_rates trans_rate;

  trans_rate.rate1 = trans_rate1 * num_spec_trait1;
  trans_rate.rate2 = trans_rate2 * num_spec_trait2;

  return trans_rate;
}


rates update_rates(double timeval,
                   double total_time,
                   double gam,
                   double laa,
                   double lac,
                   double mu,
                   double K,
                   double num_spec,
                   double num_immigrants,
                   double mainland_n,
                   const island_spec& island_spec_,
                   const rates& trait_pars) {

  double A = 1.0; // island ontogeny dynamics are currently not needed.

  size_t num_spec_trait1 = 0;
  size_t num_spec_trait2 = 0;

  size_t num_immig_trait1 = 0;
  size_t num_immig_trait2 = 0;
  if (!island_spec_.data_.empty()) {
   for (const auto& i : island_spec_.data_) {
     if (i.trait == 1) num_spec_trait1++;
     if (i.trait == 2) num_spec_trait2++;
     if (i.type_species == species_type::I && i.trait == 1) num_immig_trait1++;
     if (i.type_species == species_type::I && i.trait == 2) num_immig_trait2++;
   }
  }

  two_rates immig_rate = get_immig_rate(gam, A, num_spec, K, mainland_n,
                                        trait_pars.M2,
                                        trait_pars.immig_rate2);

  two_rates ext_rate = get_ext_rate(mu, num_spec, A,
                                    trait_pars.ext_rate2,
                                    num_spec_trait1,
                                    num_spec_trait2);

  two_rates ana_rate = get_ana_rate(laa, num_immigrants,
                                    trait_pars.ana_rate2,
                                    num_immig_trait1,
                                    num_immig_trait2);
  two_rates clado_rate = get_clado_rate(lac, num_spec, K, A,
                                        trait_pars.clado_rate2,
                                        num_spec_trait1,
                                        num_spec_trait2);
  two_rates trans_rate = get_trans_rate(trait_pars.trans_rate,
                                        trait_pars.trans_rate2,
                                        num_spec_trait1,
                                        num_spec_trait2);


  rates r(immig_rate, ext_rate, ana_rate, clado_rate, trans_rate, trait_pars.M2);

  return(r);
}


#endif /* DAISIE_RATES_H */
