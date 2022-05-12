#ifndef ISLAND_SPEC_H
#define ISLAND_SPEC_H

enum class species_type {I, A, C}; // I, A, C
enum extinction_type {clado_extinct, immig_parent};
enum class species {A, B};

struct island_spec_row {
  // 1
  double parent = -1.0;
  // 2
  double parent2 = -1.0;
  // 3
  double colonisation_time = -1.0;
  // 4
  species_type type_species;
  // 5
  std::vector< species > anc_type;
  // 6
  double extinction_time = -1.0;
  // 7
  extinction_type ext_type;
  // 8
  int trait;

  island_spec_row() {
    parent = -1; parent2 = -1;
    colonisation_time = -1.0;
    extinction_time = -1.0;
  }

  island_spec_row(double colonist, double timeval, species_type st, int trait_val) {
    parent = colonist;
    parent2 = colonist;
    colonisation_time = timeval;
    type_species = st;
    trait = trait_val;
  }

};

struct island_spec {
  std::vector< island_spec_row > data_;

  bool empty() const {
    return data_.empty();
  }

  size_t size() const {
    return data_.size();
  }

  island_spec_row& operator[](int index) {
    return data_[index];
  }

  island_spec_row operator[](int index) const {
    return data_[index];
  }

  void push_back(const island_spec_row& add) {
    data_.push_back(add);
  }
};


#endif
