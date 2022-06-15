#ifndef ISLAND_SPEC_H
#define ISLAND_SPEC_H

enum class species_type {I, A, C}; // I, A, C
enum anagenesis_type {NA, clado_extinct, immig_parent};
enum class species {A, B};

struct island_spec_row {
  // 1
  int id = -1.0;
  // 2
  int parent = -1.0;
  // 3
  double colonisation_time = -1.0;
  // 4
  species_type type_species;
  // 5
  std::vector< species > anc_type;
  // 6
  double branching_time = -1.0;
  // 7
  anagenesis_type ext_type;
  // 8
  int trait;

  island_spec_row() {
    id = -1;
    parent = -1;
    colonisation_time = -1.0;
    branching_time = -1.0;
  }

  island_spec_row(int colonist, double timeval, species_type st, int trait_val) {
    id = colonist;
    parent = colonist;
    colonisation_time = timeval;
    type_species = st;
    trait = trait_val;
  }

  island_spec_row(int id_,
                  int parent_,
                  double timeval, species_type st, int trait_val) {
    id = id_;
    parent = parent_;
    branching_time = timeval;
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

  rates(two_rates immig,
        two_rates ext,
        two_rates ana,
        two_rates clado,
        two_rates trans,
        double m2) :
    immig_rate(immig.rate1),
    ext_rate(ext.rate1),
    ana_rate(ana.rate1),
    clado_rate(clado.rate1),
    trans_rate(trans.rate1),
    immig_rate2(immig.rate2),
    ext_rate2(ext.rate2),
    ana_rate2(ana.rate2),
    clado_rate2(clado.rate2),
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


#endif
