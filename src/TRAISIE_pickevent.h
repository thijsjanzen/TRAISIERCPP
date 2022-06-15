#ifndef PICK_EVENT_H
#define PICK_EVENT_H

#include <Rcpp.h>
#include <vector>
#include "TRAISIE_rates.h"

#include <chrono>
#include <thread>

void force_output() {
  std::this_thread::sleep_for(std::chrono::milliseconds(30));
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}

void remove_species(island_spec& is,
                    int index) {
  is.data_.erase(is.data_.begin() + index);
}

bool match_motif(const std::vector< species >& anc,
                 const std::vector< species >& motif) {

  if (motif.empty()) return false;

  for (size_t i = 0; i < anc.size(); ++i) {
    if (anc[i] == motif[0]) {
      int num_matches = 0;
      for (size_t j = 0; j < motif.size(); ++j) {
        if (i + j < anc.size()) {
          if (anc[i + j] == motif[j]) num_matches++;
        }
      }
      if (num_matches == motif.size()) return true;
    }
  }
  return false;
}

void remove_cladogenetic(island_spec& is,
                         int extinct) {

 //  std::cerr << "welcome to remove cladogenetic\n"; force_output();

  std::vector<size_t> sisters;
  int num_sisters = 0;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].parent == is[extinct].parent &&
        is[i].colonisation_time == is[extinct].colonisation_time) {
      if (i != extinct) sisters.push_back(i);
      num_sisters++;
    }
  }

 // std::cerr << "num sisters:" << num_sisters << "\n"; force_output();

  if (num_sisters == 2) {
    auto survivors = sisters.front(); // using survivors to match R code, but this is only a single survivor at a time

    is[survivors].type_species = species_type::A;
    is[survivors].ext_type = anagenesis_type::clado_extinct;
    is[survivors].anc_type.clear();

    remove_species(is, extinct);
  }

  if (num_sisters > 2) {
    std::sort(sisters.begin(), sisters.end()); force_output();

    int number_of_splits = static_cast<int>(is[extinct].anc_type.size());
    auto most_recent_split = is[extinct].anc_type.back();

    auto sister_most_recent_split = species::B;

    if (most_recent_split == species::B) {
      sister_most_recent_split = species::A;
    }
    if (most_recent_split == species::A) { // for completeness and matching R.
      sister_most_recent_split = species::B;
    }
//    std::cerr << "starting looking for motif\n"; force_output();
 //   std::cerr << is.size() << " " << extinct << "\n"; force_output();
 //   std::cerr << is[extinct].anc_type.size() << "\n"; force_output();

    std::vector< species > motif_to_find = is[extinct].anc_type;

 //   std::cerr << "starting pop back\n"; force_output();
    if (!motif_to_find.empty()) motif_to_find.pop_back();

 //   std::cerr << "adding sister\n"; force_output();
    motif_to_find.push_back(sister_most_recent_split);

 //   std::cerr << "motif: ";
  //  for (const auto& m : motif_to_find) {
//      if( m == species::A) std::cerr << "A";
//      if( m == species::B) std::cerr << "B";
  //  } std::cerr << "\n";

    std::vector<int> possible_sister;

 //   std::cerr << "starting motif matching\n"; force_output();
    for (size_t i = 0; i < sisters.size(); ++i) {
      size_t survivor = sisters[i];
      if (match_motif(is[survivor].anc_type, motif_to_find)) {
        possible_sister.push_back(survivor);
      }
    }
 //   std::cerr << "possible sisters: "; force_output();
 //   for (auto p : possible_sister) {
 //     std::cerr << p << " ";
//    } std::cerr << "\n"; force_output();

    if (most_recent_split == species::A) {
      size_t first_sister = possible_sister.front();
      is[first_sister].branching_time = is[extinct].branching_time;
    }
    // remove the offending A/B

 //   std::cerr << "resulting motifs: ";
    for (auto ps : possible_sister) {
      is[ps].anc_type.erase(is[ps].anc_type.begin() + number_of_splits - 1);


 //     for (const auto& m : is[ps].anc_type) {
   //     if( m == species::A) std::cerr << "A";
  //      if( m == species::B) std::cerr << "B";
  //    } std::cerr << "\n"; force_output();
    }
    remove_species(is, extinct);
  }
  return;
}

template <typename T>
int draw_prop(const std::vector<T>& v) {
  T s = std::accumulate(v.begin(), v.end(), 0.0);
  double r = R::runif(0, s);
  int index = 0;
  for( ; index < v.size(); ++index) {
    r -= v[index];
    if (r <= 0.0) {
      break;
    }
  }
  return index;
}

template <typename T>
T draw_from(const std::vector<T>& v) {
  double r = R::runif(0, v.size());
  int index = static_cast<size_t>(r);
  return v[index];
}

void execute_immigration(island_spec& is,
                         double timeval,
                         int colonist,
                         int focal_trait) {
  int is_it_there = -1;
  if (!is.empty()) {
    for (size_t i = 0; i < is.size(); ++i) {
      if (is[i].id == colonist) {
        is_it_there = i;
        break;
      }
    }
  }

  island_spec_row add(colonist, timeval, species_type::I, focal_trait);

  if (is_it_there == -1) { // it is not there
    is.push_back(add);
  } else {
    is[is_it_there] = add;
  }
  return;
}



void immigration(double timeval,
                 const std::vector<double>& mainland_spec,
                 island_spec& is) {

  int colonist = draw_prop(mainland_spec);
  execute_immigration(is, timeval, colonist, 1);
  return;
}

void execute_extinction(island_spec& is,
                        int extinct) {

  auto type_of_species = is[extinct].type_species;

  if (type_of_species == species_type::I) { // 0 == "I
    remove_species(is, extinct); // island_spec = island_spec[-extinct,]
  }
  if (type_of_species == species_type::A) { // 1 == "A"
    remove_species(is, extinct); // island_spec = island_spec[-extinct,]
  }
  if (type_of_species == species_type::C) { // 2 == "C"
    remove_cladogenetic(is, extinct);
  }
}


void extinction(island_spec& is,
                int focal_trait) {

  if (is.empty()) return;

  std::vector<size_t> island_spec_state;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].trait == focal_trait) island_spec_state.push_back(i);
  }

  int extinct = draw_from(island_spec_state);
  execute_extinction(is, extinct);

  return;
}

void anagenesis(island_spec& is,
                int& maxspecID,
                int focal_trait) {

  std::vector<size_t> immi_specs;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].type_species == species_type::I && is[i].trait == focal_trait) {
      immi_specs.push_back(i);
    }
  }

  size_t anagenesis = immi_specs[0];

  if (immi_specs.size() > 1) {
    int index = static_cast<int>(R::runif(0, immi_specs.size()));
    anagenesis = immi_specs[index];
  }
   maxspecID++;

  is[anagenesis].type_species = species_type::A;
  is[anagenesis].id = maxspecID;
  is[anagenesis].ext_type = anagenesis_type::immig_parent;
  is[anagenesis].trait = focal_trait;

  return;
}

void cladogenesis(island_spec& is,
                  double timeval,
                  int maxspecID,
                  int focal_trait) {

  std::vector<size_t> island_spec_state;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].trait == focal_trait) island_spec_state.push_back(i);
  }

  int index = static_cast<int>(R::runif(0, island_spec_state.size()));
  auto to_split = island_spec_state[index];

  if (is[to_split].type_species == species_type::C) {
    // for daughter A
    is[to_split].type_species = species_type::C; // redundant, but following R
    is[to_split].id = maxspecID + 1;
    auto oldstatus = is[to_split].anc_type;
    is[to_split].anc_type.push_back(species::A);
    is[to_split].ext_type = anagenesis_type::NA;
    is[to_split].trait = focal_trait;

    // for daughter B
    island_spec_row add;
    add.id = maxspecID + 2;
    add.parent = is[to_split].parent;
    add.colonisation_time = is[to_split].colonisation_time;
    add.type_species = species_type::C;
    add.anc_type = oldstatus;
    add.anc_type.push_back(species::B);
    add.branching_time = timeval;
    add.ext_type = anagenesis_type::NA;
    add.trait = focal_trait;
    is.push_back(add);

    maxspecID += 2;
  } else {
    // not cladogenetic
    // for daughter A
    is[to_split].type_species = species_type::C;
    is[to_split].id = maxspecID + 1;
    is[to_split].anc_type = {species::A};
    is[to_split].branching_time = is[to_split].colonisation_time;
    is[to_split].ext_type = anagenesis_type::NA;
    is[to_split].trait = focal_trait;

    island_spec_row add;
    add.id = maxspecID + 2;
    add.parent = is[to_split].parent;
    add.colonisation_time = is[to_split].colonisation_time;
    add.type_species = species_type::C;
    add.anc_type = {species::B};
    add.branching_time = timeval;
    add.ext_type = anagenesis_type::NA;
    add.trait = focal_trait;
    is.push_back(add);
    maxspecID += 2;
  }
  return;
}

void transition(island_spec& is,
                int focal_trait) {
  std::vector<size_t> island_spec_state;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].trait == focal_trait) island_spec_state.push_back(i);
  }

  int index = static_cast<int>(R::runif(0, island_spec_state.size()));
  auto totrans = island_spec_state[index];
  if (focal_trait == 1) {
    is[totrans].trait = 2;
  } else {
    is[totrans].trait = 1;
  }

  return;
}



int sample_spec(const std::vector<double>& mainland_spec,
                int mainland2) {
  auto mainland1 = mainland_spec.size();
  //  auto mainland_total = mainland1 + mainland2;
  int index = static_cast<int>(R::runif(0, mainland2));
  auto colonist = index + mainland1 + 1;
  return colonist;
}

void immigration_state2(island_spec& is,
                        std::vector<double> mainland_spec,
                        double timeval,
                        const int M2) {

  int colonist = sample_spec(mainland_spec, M2);
  execute_immigration(is, timeval, colonist, 2);

  return;
}

void DAISIE_sim_update_state_trait_dep(double timeval,
                                       double total_time,
                                       int event,
                                       int& maxspecID,
                                       const std::vector<double>& mainland_spec,
                                       island_spec& is,
                                       const rates& trait_pars,
                                       std::vector< std::array<double, 7>>& stt_table) {

  switch(event) {
    case 1: { immigration(timeval, mainland_spec, is); break;}
    case 2: { extinction(is, 1); break;}
    case 3: { anagenesis(is, maxspecID, 1); break;}
    case 4: { cladogenesis(is, timeval, maxspecID, 1); break;}
    case 5: { transition(is, 1); break;}
    case 6: { immigration_state2(is, mainland_spec, timeval, trait_pars.M2); break;}
    case 7: { extinction(is, 2); break;}
    case 8: { anagenesis(is, maxspecID, 2); break;}
    case 9: { cladogenesis(is, timeval, maxspecID, 2); break;}
    case 10:{ transition(is, 2); break;}
  }

  if (total_time >= timeval) {

    double nI1 = 0, nA1 = 0, nC1 = 0, nI2 = 0, nA2 = 0, nC2 = 0;
    for (const auto& i : is.data_) {
      if (i.type_species == species_type::I && i.trait == 1) nI1++;
      if (i.type_species == species_type::I && i.trait == 2) nI2++;
      if (i.type_species == species_type::A && i.trait == 1) nA1++;
      if (i.type_species == species_type::A && i.trait == 2) nA2++;
      if (i.type_species == species_type::C && i.trait == 1) nC1++;
      if (i.type_species == species_type::C && i.trait == 2) nC2++;
    }

    std::array<double, 7> add =
                          { total_time - timeval, nI1, nA1, nC1, nI2, nA2, nC2};
    stt_table.push_back(add);
  }
  return;
}




#endif
