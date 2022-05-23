#ifndef TEST_TRAISIE_PICKEVENT_H
#define TEST_TRAISIE_PICKEVENT_H

#include <Rcpp.h>
#include "TRAISIE_pickevent.h"
#include "TRAISIE_util.h"

//' test immigration
//' @param timeval timeval
//' @param mainland_spec vector
//' @param island_spec_R matrix
//' @return stringmatrix island_spec
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix test_immigration(double timeval,
                      Rcpp::NumericVector mainland_spec_R,
                      const Rcpp::StringVector& island_spec_R) {

  island_spec island_spec_;
    std::string s_colonist = Rcpp::as<std::string>(island_spec_R(0));
    std::string s_timeval  = Rcpp::as<std::string>(island_spec_R(1));
    std::string s_trait_val = Rcpp::as<std::string>(island_spec_R(7));

    double colonist = std::stod(s_colonist);
    double local_timeval  = std::stod(s_timeval);
    double trait_val = std::stod(s_trait_val);

    island_spec_row to_add(colonist, local_timeval, species_type::I, trait_val);
    island_spec_.push_back(to_add);


  std::vector<double> mainland_spec(mainland_spec_R.begin(),
                                    mainland_spec_R.end());

  immigration(timeval, mainland_spec,
              island_spec_);

  auto out = make_island_spec_for_R(island_spec_);

  return out;
}




#endif
