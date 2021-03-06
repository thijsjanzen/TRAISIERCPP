#ifndef TRAISIE_UTIL_H
#define TRAISIE_UTIL_H

#include <vector>
#include "TRAISIE_structs.h"

Rcpp::NumericMatrix make_stt_table_for_R(const std::vector< std::array< double, 7 >>& stt_table) {
  int num_rows = stt_table.size();
  int num_cols = 7;
  Rcpp::NumericMatrix out(num_rows, num_cols);
  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < num_cols; ++j) {
      out(i, j) = stt_table[i][j];
    }
  }
  return out;
}

std::string get_string(const extinction_type& ex) {
  if (ex == extinction_type::clado_extinct) {
    return "Clado_extinct";
  }
  if (ex == extinction_type::immig_parent) {
    return "Immig_parent";
  }
  return "Incorrect_entry";
}

std::string get_string(const species_type& st) {
  if (st == species_type::I) return "I";
  if (st == species_type::A) return "A";
  if (st == species_type::C) return "C";
  return "NA";
}

std::string get_string(const std::vector<species>& anc) {
  std::string out;
  for (const auto& i : anc) {
    if (i == species::A) out.append("A");
    if (i == species::B) out.append("B");
  }
  return out;
}

Rcpp::StringMatrix make_island_spec_for_R(const island_spec& is) {
  int num_rows = is.size();
  int num_cols = 8;
  Rcpp::StringMatrix out(num_rows, num_cols);
  for (size_t i = 0; i < is.size(); ++i) {
    out(i, 0) = std::to_string(is[i].id);
    out(i, 1) = std::to_string(is[i].parent);
    out(i, 2) = std::to_string(is[i].colonisation_time);
    out(i, 3) = get_string(is[i].type_species);
    out(i, 4) = get_string(is[i].anc_type);
    out(i, 5) = std::to_string(is[i].extinction_time);
    out(i, 6) = get_string(is[i].ext_type);
    out(i, 7) = std::to_string(is[i].trait);
  }
  return out;
}

#endif
