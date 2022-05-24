#ifndef TEST_TRAISIE_PICKEVENT_H
#define TEST_TRAISIE_PICKEVENT_H

#include <Rcpp.h>
#include "TRAISIE_pickevent.h"
#include "TRAISIE_util.h"

//' test draw prop
//' @param Rcpp::NumericVector probs
//' @return drawn index
//' @export
// [[Rcpp::export]]
int test_draw_prop(const Rcpp::NumericVector& v) {
  std::vector<double> v2(v.begin(), v.end());
  auto x = draw_prop(v2);
  return x;
}

island_spec create_island_spec(const Rcpp::StringMatrix& island_spec_R) {
  island_spec island_spec_;
  for (size_t i = 0; i < island_spec_R.nrow(); ++i) {

    std::string s_c = Rcpp::as<std::string>(island_spec_R(i, 0));
    int colonist = std::stoi(s_c);

    std::string s_t = Rcpp::as<std::string>(island_spec_R(i, 2));
    double timeval = std::stod(s_t);

    std::string s_tv = Rcpp::as<std::string>(island_spec_R(i, 7));
    int trait_val = std::stoi(s_tv);

    std::string spec_type = Rcpp::as<std::string>(island_spec_R(i, 3));
    species_type local_spec_type;
    if (spec_type == "I") local_spec_type = species_type::I;
    if (spec_type == "C") local_spec_type = species_type::C;
    if (spec_type == "A") local_spec_type = species_type::A;

    island_spec_row to_add(colonist, timeval, local_spec_type, trait_val);

    std::string motif = Rcpp::as<std::string>(island_spec_R(i, 4));
    for (char& c : motif) {
      if (c == 'A') to_add.anc_type.push_back(species::A);
      if (c == 'B') to_add.anc_type.push_back(species::B);
    }

    island_spec_.push_back(to_add);
  }
  return island_spec_;
}



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

//' test extinction
//' @param island_spec_R island spec table
//' @param focal_trait focal trait
//' @return island spec table
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix test_extinction(const Rcpp::StringMatrix& island_spec_R,
                                   int focal_trait) {

  island_spec island_spec_ = create_island_spec(island_spec_R);

  extinction(island_spec_, focal_trait);

  Rcpp::StringMatrix  out = make_island_spec_for_R(island_spec_);
  return out;
}

//' test execute extinction without rng
//' @param island_spec_R island spec table
//' @param index index
//' @return island spec table
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix test_execute_extinction(const Rcpp::StringMatrix& island_spec_R,
                                           int index) {
  index = index - 1;
  island_spec island_spec_ = create_island_spec(island_spec_R);

  execute_extinction(island_spec_, index);

  Rcpp::StringMatrix  out = make_island_spec_for_R(island_spec_);
  return out;
}

//' test execute anagenesis
//' @param island_spec_R island spec table
//' @param maxspecID maxspecid
//' @param focal_trait focaltrait
//' @return island spec table
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix test_anagenesis(Rcpp::StringMatrix& island_spec_R,
                                    int maxspecID,
                                    int focal_trait) {

  island_spec island_spec_ = create_island_spec(island_spec_R);

  anagenesis(island_spec_, maxspecID, focal_trait);

  Rcpp::StringMatrix  out = make_island_spec_for_R(island_spec_);
  return out;
}

//' test execute cladogenesis
//' @param island_spec_R island spec table
//' @param timeval timeval
//' @param maxspecID maxspecid
//' @param focal_trait focaltrait
//' @return island spec table
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix test_cladogenesis(Rcpp::StringMatrix& island_spec_R,
                                     double timeval,
                                     int maxspecID,
                                     int focal_trait) {
  island_spec island_spec_ = create_island_spec(island_spec_R);
  cladogenesis(island_spec_, timeval, maxspecID, focal_trait);
  Rcpp::StringMatrix  out = make_island_spec_for_R(island_spec_);
  return out;
}

//' test transition
//' @param island_spec_R island spec table
//' @param focal_trait focaltrait
//' @return island spec table
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix test_transition(Rcpp::StringMatrix& island_spec_R,
                                    int focal_trait) {
  island_spec island_spec_ = create_island_spec(island_spec_R);
  transition(island_spec_, focal_trait);
  Rcpp::StringMatrix  out = make_island_spec_for_R(island_spec_);
  return out;
}

//' test sample spec
//' @param mainland_spec mainland spec table
//' @param M2 number species of trait 2 on mainland
//' @return species index
//' @export
// [[Rcpp::export]]
int test_sample_spec(const Rcpp::NumericVector& mainland_spec_R,
                     int M2) {
  std::vector<double> mainland_spec(mainland_spec_R.begin(),
                                    mainland_spec_R.end());
  return sample_spec(mainland_spec, M2);
}

//' test immigration2
//' @param island_spec_R island spec table
//' @param mainland_spec_R mainland species vector
//' @param timeval timeval
//' @param M2 M2
//' @return island spec table
//' @export
// [[Rcpp::export]]
Rcpp::StringMatrix test_immigration2(Rcpp::StringMatrix& island_spec_R,
                                     const Rcpp::NumericVector& mainland_spec_R,
                                     double timeval,
                                     int M2) {
  island_spec island_spec_ = create_island_spec(island_spec_R);
  std::vector<double> mainland_spec(mainland_spec_R.begin(),
                                    mainland_spec_R.end());
  immigration_state2(island_spec_,
                     mainland_spec,
                     timeval,
                     M2);

  Rcpp::StringMatrix  out = make_island_spec_for_R(island_spec_);
  return out;
}





#endif
