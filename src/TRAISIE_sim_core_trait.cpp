#include <Rcpp.h>

#include <vector>

#include <thread>
#include <chrono>

#include "TRAISIE_rates.h"
#include "TRAISIE_pickevent.h"
#include "TRAISIE_island_spec.h"
#include "TRAISIE_util.h"

double calc_next_time_eval(const rates& r,
                           double timeval) {
  double s = r.sum();
  double dt = R::rexp(s);

  return timeval + dt;
}

void output(std::string focal_text) {
  Rcpp::Rcout << focal_text << "\n";
 // std::this_thread::sleep_for(std::chrono::milliseconds(3));
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}

int update_num_immigrants(const island_spec& is) {
  int num_immigrants = 0;
  for (const auto& i : is.data_) {
    if (i.type_species == species_type::I) num_immigrants++;
  }
  return num_immigrants;
}


//' CPP function to execute expensive loop
//' @param timeval timeval
//' @param total_time total time
//' @param gam immigration rate
//' @param laa anagenetic rate
//' @param lac cladogenetic rate
//' @param mu extinction rate
//' @param area_pars_from_R area pars
//' @param K K
//' @param num_spec number of species
//' @param num_immigrants number of immigrants
//' @param mainland_n number of species on mainland
//' @param maxspecID maxpsecid
//' @param trait_pars_R trait pars
//' @param mainland_spec_R mainland spec
//' @return list with stt_table and island_spec
//' @export
// [[Rcpp::export]]
Rcpp::List execute_time_loop(double timeval,
                      double total_time,
                      double gam,
                      double laa,
                      double lac,
                      double mu,
                      Rcpp::List area_pars_from_R,
                      double K,
                      double mainland_n,
                      Rcpp::List trait_pars_R) {

  area_pars ap(area_pars_from_R);
 // output("area_pars loaded");
  rates trait_pars(gam, laa, lac, mu, trait_pars_R);
 // output("trait_pars loaded");


  auto mainland_ntotal = mainland_n + trait_pars.M2;
  std::vector<double> mainland_spec;
  if (mainland_n != 0) {
    for (double i = 1.0; i <= mainland_n; ++i) {
      mainland_spec.push_back(i);
    }
  }
  int maxspecID = static_cast<int>(mainland_ntotal);

  std::vector< std::array< double, 7> > stt_table;
  stt_table.push_back({total_time, 0, 0, 0, 0, 0, 0});
  island_spec island_spec_;

  int num_spec = island_spec_.size();
  int num_immigrants = update_num_immigrants(island_spec_);

 // output("starting loop");



  while (timeval < total_time) {
    auto all_rates = update_rates(
      timeval,
      total_time,
      gam,
      laa,
      lac,
      mu,
      ap,
      K,
      num_spec,
      num_immigrants,
      mainland_n,
      island_spec_,
      trait_pars);


    timeval = calc_next_time_eval(all_rates, timeval);
  //  Rcpp::Rcout << timeval << "\n";

    if (timeval < total_time) {
      auto possible_event = all_rates.sample_event();
   //   Rcpp::Rcout << possible_event << "\n";

      DAISIE_sim_update_state_trait_dep(timeval,
                                        total_time,
                                        possible_event,
                                        maxspecID,
                                        mainland_spec,
                                        island_spec_,
                                        trait_pars,
                                        stt_table);

      num_spec = island_spec_.size();
      num_immigrants = update_num_immigrants(island_spec_);
    }
  }
  // finalize stt_table
  //output("done loop, finalizing");

  std::array<double, 7> add = stt_table.back();
  add[0] = 0.0;
  stt_table.push_back(add);

//  output("starting conversion stt");
  Rcpp::NumericMatrix stt_table_for_R = make_stt_table_for_R(stt_table);
//  output("done conversion stt");
 // output("starting conversion island spec");
  Rcpp::StringMatrix  island_spec_for_R = make_island_spec_for_R(island_spec_);
 // output("done conversion island spec");

  Rcpp::List output = Rcpp::List::create(
          Rcpp::Named("stt_table") = stt_table_for_R,
          Rcpp::Named("island_spec") = island_spec_for_R);
  return output;
}
