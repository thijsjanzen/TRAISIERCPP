#ifndef PICK_EVENT_H
#define PICK_EVENT_H

#include <Rcpp.h>
#include <vector>
#include <random>
#include "TRAISIE_rates.h"
#include "TRAISIE_island_spec.h"

void remove_species(island_spec& is,
                    int index) {
  is.data_.erase(is.data_.begin() + index);
}

bool match_motif(const std::vector< species >& anc,
                 const std::vector< species >& motif) {
  for (size_t i = 0; i < motif.size(); ++i) {
    if (i > anc.size()) return false;
    if (anc[i] != motif[i]) return false;
  }
  return true;
}

void remove_cladogenetic(island_spec& is,
                         int extinct) {

  std::vector<size_t> sisters;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].parent2 == is[extinct].parent2 &&
        is[i].colonisation_time == is[extinct].colonisation_time) {
      if (i != extinct) sisters.push_back(i);
    }
  }

  if (sisters.size() == 2) {
    for (size_t i = 0; i < sisters.size(); ++i) {
      auto survivors = sisters[i]; // using survivors to match R code, but this is only a single survivor at a time
      is[survivors].type_species = species_type::A;
      is[survivors].ext_type = extinction_type::clado_extinct;
    }
    is.data_.erase(is.data_.begin() + extinct);
  }
  if (sisters.size() > 2) {
    auto number_of_splits = is[extinct].anc_type.size();
    auto most_recent_split = is[extinct].anc_type.back();

    auto sister_most_recent_split = species::B;

    if (most_recent_split == species::B) {
      sister_most_recent_split = species::A;
    }
    if (most_recent_split == species::A) { // for completeness and matching R.
      sister_most_recent_split = species::B;
    }

    std::vector< species > motif_to_find;
    for (size_t i = 0; i < number_of_splits; ++i) {
      motif_to_find.push_back(is[extinct].anc_type[i]);
    }
    motif_to_find.push_back(sister_most_recent_split);

    auto possible_sister = 0;

    for (size_t i = 0; i < sisters.size(); ++i) {
      size_t survivor = sisters[i];
      if (match_motif(is[survivor].anc_type, motif_to_find)) {
        possible_sister = survivor; break;
      }
    }

    if (most_recent_split == species::A) {
      is[possible_sister].extinction_time = is[extinct].extinction_time;
    }
    // remove the offending A/B
    is[possible_sister].anc_type.erase(is[possible_sister].anc_type.begin() + number_of_splits);
    remove_species(is, extinct);
  }
  return;
}


void immigration(double timeval,
                 const std::vector<double>& mainland_spec,
                 island_spec& is,
                 std::mt19937& rndgen) {

  std::discrete_distribution<int> dist(mainland_spec.begin(), mainland_spec.end());
  int colonist = dist(rndgen);

  int is_it_there = -1;
  if (!is.empty()) {
      for (size_t i = 0; i < is.size(); ++i) {
        if (is[i].parent == colonist) {
          is_it_there = i;
          break;
        }
      }
  }

  island_spec_row add(colonist, timeval, species_type::I, 1);

  if (is_it_there == -1) { // it is not there
     is.push_back(add);
  } else {
    is[is_it_there] = add;
  }
  return;
}

void extinction(island_spec& is,
                std::mt19937& rndgen,
                int focal_trait) {

  if (is.empty()) return;

  std::vector<size_t> island_spec_state;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].trait == focal_trait) island_spec_state.push_back(i);
  }
  std::discrete_distribution<int> dist(island_spec_state.begin(),
                                       island_spec_state.end());

  int extinct = dist(rndgen);

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
  return;
}

void anagenesis(island_spec& is,
                int& maxspecID,
                std::mt19937& rndgen,
                int focal_trait) {

   std::vector<size_t> immi_specs;
   for (size_t i = 0; i < is.size(); ++i) {
     if (is[i].type_species == species_type::I && is[i].trait == focal_trait) {
       immi_specs.push_back(i);
     }
   }

   auto anagenesis = immi_specs[0];

   if (immi_specs.size() > 1) {
     std::uniform_int_distribution<int> dist(0, immi_specs.size() - 1);
     auto index = dist(rndgen);
     anagenesis = immi_specs[index];
   }

   maxspecID++;

   is[anagenesis].type_species = species_type::A;
   is[anagenesis].parent = maxspecID;
   is[anagenesis].ext_type = extinction_type::immig_parent;
   is[anagenesis].trait = focal_trait;

   return;
}

void cladogenesis(island_spec& is,
                  double timeval,
                  int maxspecID,
                  std::mt19937& rndgen,
                  int focal_trait) {

  std::vector<size_t> island_spec_state;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].trait == focal_trait) island_spec_state.push_back(i);
  }

  std::uniform_int_distribution<int> dist(0, island_spec_state.size() - 1);
  auto index = dist(rndgen);
  auto to_split = island_spec_state[index];

  if (is[to_split].type_species == species_type::C) {
      // for daughter A
      is[to_split].type_species = species_type::C; // redunandant, but following R
      is[to_split].parent = maxspecID + 1;
      auto oldstatus = is[to_split].anc_type;
      is[to_split].anc_type.push_back(species::A);
      is[to_split].trait = focal_trait;

      // for daughter B
      island_spec_row add;
      add.parent = maxspecID + 2;
      add.parent2 = is[to_split].parent2;
      add.colonisation_time = is[to_split].colonisation_time;
      add.type_species = species_type::C;
      add.anc_type = oldstatus; add.anc_type.push_back(species::B);
      add.extinction_time = timeval;
      add.trait = focal_trait;
      is.push_back(add);

      maxspecID += 2;
  } else {
    // not cladogenetic
    // for daughter A
    is[to_split].type_species = species_type::C;
    is[to_split].parent = maxspecID + 1;
    is[to_split].anc_type = {species::A};
    is[to_split].extinction_time = is[to_split].colonisation_time;
    is[to_split].trait = focal_trait;

    island_spec_row add;
    add.parent = maxspecID + 2;
    add.parent2 = is[to_split].parent2;
    add.colonisation_time = is[to_split].colonisation_time;
    add.type_species = species_type::C;
    add.anc_type = {species::B};
    add.extinction_time = timeval;
    add.trait = focal_trait;
    is.push_back(add);
    maxspecID += 2;
  }
  return;
}

void transition(island_spec& is,
                std::mt19937& rndgen,
                int focal_trait) {
  std::vector<size_t> island_spec_state;
  for (size_t i = 0; i < is.size(); ++i) {
    if (is[i].trait == focal_trait) island_spec_state.push_back(i);
  }
  std::uniform_int_distribution<int> dist(0, island_spec_state.size() - 1);
  auto index = dist(rndgen);
  auto totrans = island_spec_state[index];
  if (focal_trait == 1) {
    is[totrans].trait = 2;
  } else {
    is[totrans].trait = 1;
  }

  return;
}

void immigration_state2(island_spec& is,
                        std::vector<double> mainland_spec,
                        double timeval,
                        const rates& trait_pars,
                        std::mt19937& rndgen) {

  auto mainland1 = mainland_spec.size();
  auto mainland2 = trait_pars.M2;
//  auto mainland_total = mainland1 + mainland2;
  std::uniform_int_distribution<int> dist(0, mainland2);
  auto index = dist(rndgen);
  auto colonist = index + mainland1 + 1;

  int is_it_there = -1;
  if (!is.empty()) {
    for (size_t i = 0; i < is.size(); ++i) {
      if (is[i].parent == colonist) {
        is_it_there = i;
        break;
      }
    }
  }

  island_spec_row add(colonist, timeval, species_type::I, 2);

  if (is_it_there == -1) { // it is not there
    is.push_back(add);
  } else {
    is[is_it_there] = add;
  }
  return;
}

void DAISIE_sim_update_state_trait_dep(double timeval,
                                       double total_time,
                                       int event,
                                       int& maxspecID,
                                       const std::vector<double>& mainland_spec,
                                       island_spec& is,
                                       const rates& trait_pars,
                                       std::vector< std::array<double, 7>>& stt_table,
                                       std::mt19937& rndgen) {

  switch(event) {
    case 1: { immigration(timeval, mainland_spec, is, rndgen); break;}
    case 2: { extinction(is, rndgen, 1); break;}
    case 3: { anagenesis(is, maxspecID, rndgen, 1); break;}
    case 4: { cladogenesis(is, maxspecID, timeval, rndgen, 1); break;}
    case 5: { transition(is, rndgen, 1); break;}
    case 6: { immigration_state2(is, mainland_spec, timeval, trait_pars, rndgen); break;}
    case 7: { extinction(is, rndgen, 2); break;}
    case 8: { anagenesis(is, maxspecID, rndgen, 2); break;}
    case 9: { cladogenesis(is, maxspecID, timeval, rndgen, 2); break;}
    case 10:{ transition(is, rndgen, 2); break;}
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
