// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sample_event
int sample_event(const std::vector<double>& event_prob);
RcppExport SEXP _TRAISIERCPP_sample_event(SEXP event_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type event_prob(event_probSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_event(event_prob));
    return rcpp_result_gen;
END_RCPP
}
// execute_time_loop
Rcpp::List execute_time_loop(double timeval, double total_time, double gam, double laa, double lac, double mu, double K, double mainland_n, Rcpp::List trait_pars_R);
RcppExport SEXP _TRAISIERCPP_execute_time_loop(SEXP timevalSEXP, SEXP total_timeSEXP, SEXP gamSEXP, SEXP laaSEXP, SEXP lacSEXP, SEXP muSEXP, SEXP KSEXP, SEXP mainland_nSEXP, SEXP trait_pars_RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type timeval(timevalSEXP);
    Rcpp::traits::input_parameter< double >::type total_time(total_timeSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< double >::type laa(laaSEXP);
    Rcpp::traits::input_parameter< double >::type lac(lacSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type mainland_n(mainland_nSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type trait_pars_R(trait_pars_RSEXP);
    rcpp_result_gen = Rcpp::wrap(execute_time_loop(timeval, total_time, gam, laa, lac, mu, K, mainland_n, trait_pars_R));
    return rcpp_result_gen;
END_RCPP
}
// test_get_immig_rate
Rcpp::NumericVector test_get_immig_rate(double gam, double A, int num_spec, double K, double mainland_n, double mainland_n2, double immig_rate2);
RcppExport SEXP _TRAISIERCPP_test_get_immig_rate(SEXP gamSEXP, SEXP ASEXP, SEXP num_specSEXP, SEXP KSEXP, SEXP mainland_nSEXP, SEXP mainland_n2SEXP, SEXP immig_rate2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type num_spec(num_specSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type mainland_n(mainland_nSEXP);
    Rcpp::traits::input_parameter< double >::type mainland_n2(mainland_n2SEXP);
    Rcpp::traits::input_parameter< double >::type immig_rate2(immig_rate2SEXP);
    rcpp_result_gen = Rcpp::wrap(test_get_immig_rate(gam, A, num_spec, K, mainland_n, mainland_n2, immig_rate2));
    return rcpp_result_gen;
END_RCPP
}
// test_get_ext_rate
Rcpp::NumericVector test_get_ext_rate(double mu, int num_spec, double A, double ext_rate2, int num_spec_trait1, int num_spec_trait2);
RcppExport SEXP _TRAISIERCPP_test_get_ext_rate(SEXP muSEXP, SEXP num_specSEXP, SEXP ASEXP, SEXP ext_rate2SEXP, SEXP num_spec_trait1SEXP, SEXP num_spec_trait2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type num_spec(num_specSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type ext_rate2(ext_rate2SEXP);
    Rcpp::traits::input_parameter< int >::type num_spec_trait1(num_spec_trait1SEXP);
    Rcpp::traits::input_parameter< int >::type num_spec_trait2(num_spec_trait2SEXP);
    rcpp_result_gen = Rcpp::wrap(test_get_ext_rate(mu, num_spec, A, ext_rate2, num_spec_trait1, num_spec_trait2));
    return rcpp_result_gen;
END_RCPP
}
// test_get_ana_rate
Rcpp::NumericVector test_get_ana_rate(double laa, int num_immigrants, double ana_rate2, size_t num_spec_trait1, size_t num_spec_trait2);
RcppExport SEXP _TRAISIERCPP_test_get_ana_rate(SEXP laaSEXP, SEXP num_immigrantsSEXP, SEXP ana_rate2SEXP, SEXP num_spec_trait1SEXP, SEXP num_spec_trait2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type laa(laaSEXP);
    Rcpp::traits::input_parameter< int >::type num_immigrants(num_immigrantsSEXP);
    Rcpp::traits::input_parameter< double >::type ana_rate2(ana_rate2SEXP);
    Rcpp::traits::input_parameter< size_t >::type num_spec_trait1(num_spec_trait1SEXP);
    Rcpp::traits::input_parameter< size_t >::type num_spec_trait2(num_spec_trait2SEXP);
    rcpp_result_gen = Rcpp::wrap(test_get_ana_rate(laa, num_immigrants, ana_rate2, num_spec_trait1, num_spec_trait2));
    return rcpp_result_gen;
END_RCPP
}
// test_get_clado_rate
Rcpp::NumericVector test_get_clado_rate(double lac, double num_spec, double K, double A, double clado_rate2, size_t num_spec_trait1, size_t num_spec_trait2);
RcppExport SEXP _TRAISIERCPP_test_get_clado_rate(SEXP lacSEXP, SEXP num_specSEXP, SEXP KSEXP, SEXP ASEXP, SEXP clado_rate2SEXP, SEXP num_spec_trait1SEXP, SEXP num_spec_trait2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lac(lacSEXP);
    Rcpp::traits::input_parameter< double >::type num_spec(num_specSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type clado_rate2(clado_rate2SEXP);
    Rcpp::traits::input_parameter< size_t >::type num_spec_trait1(num_spec_trait1SEXP);
    Rcpp::traits::input_parameter< size_t >::type num_spec_trait2(num_spec_trait2SEXP);
    rcpp_result_gen = Rcpp::wrap(test_get_clado_rate(lac, num_spec, K, A, clado_rate2, num_spec_trait1, num_spec_trait2));
    return rcpp_result_gen;
END_RCPP
}
// test_get_trans_rate
Rcpp::NumericVector test_get_trans_rate(double trans_rate, double trans_rate2, size_t num_spec_trait1, size_t num_spec_trait2);
RcppExport SEXP _TRAISIERCPP_test_get_trans_rate(SEXP trans_rateSEXP, SEXP trans_rate2SEXP, SEXP num_spec_trait1SEXP, SEXP num_spec_trait2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type trans_rate(trans_rateSEXP);
    Rcpp::traits::input_parameter< double >::type trans_rate2(trans_rate2SEXP);
    Rcpp::traits::input_parameter< size_t >::type num_spec_trait1(num_spec_trait1SEXP);
    Rcpp::traits::input_parameter< size_t >::type num_spec_trait2(num_spec_trait2SEXP);
    rcpp_result_gen = Rcpp::wrap(test_get_trans_rate(trans_rate, trans_rate2, num_spec_trait1, num_spec_trait2));
    return rcpp_result_gen;
END_RCPP
}
// test_update_rates
Rcpp::NumericVector test_update_rates(double timeval, double total_time, double gam, double laa, double lac, double mu, double K, double num_spec, double num_immigrants, double mainland_n, const Rcpp::NumericMatrix& island_spec_R, const Rcpp::List& trait_pars_R);
RcppExport SEXP _TRAISIERCPP_test_update_rates(SEXP timevalSEXP, SEXP total_timeSEXP, SEXP gamSEXP, SEXP laaSEXP, SEXP lacSEXP, SEXP muSEXP, SEXP KSEXP, SEXP num_specSEXP, SEXP num_immigrantsSEXP, SEXP mainland_nSEXP, SEXP island_spec_RSEXP, SEXP trait_pars_RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type timeval(timevalSEXP);
    Rcpp::traits::input_parameter< double >::type total_time(total_timeSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< double >::type laa(laaSEXP);
    Rcpp::traits::input_parameter< double >::type lac(lacSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type num_spec(num_specSEXP);
    Rcpp::traits::input_parameter< double >::type num_immigrants(num_immigrantsSEXP);
    Rcpp::traits::input_parameter< double >::type mainland_n(mainland_nSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type island_spec_R(island_spec_RSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type trait_pars_R(trait_pars_RSEXP);
    rcpp_result_gen = Rcpp::wrap(test_update_rates(timeval, total_time, gam, laa, lac, mu, K, num_spec, num_immigrants, mainland_n, island_spec_R, trait_pars_R));
    return rcpp_result_gen;
END_RCPP
}
// test_draw_prop
int test_draw_prop(const Rcpp::NumericVector& v);
RcppExport SEXP _TRAISIERCPP_test_draw_prop(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(test_draw_prop(v));
    return rcpp_result_gen;
END_RCPP
}
// test_immigration
Rcpp::StringMatrix test_immigration(double timeval, Rcpp::NumericVector mainland_spec_R, const Rcpp::StringVector& island_spec_R);
RcppExport SEXP _TRAISIERCPP_test_immigration(SEXP timevalSEXP, SEXP mainland_spec_RSEXP, SEXP island_spec_RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type timeval(timevalSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mainland_spec_R(mainland_spec_RSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type island_spec_R(island_spec_RSEXP);
    rcpp_result_gen = Rcpp::wrap(test_immigration(timeval, mainland_spec_R, island_spec_R));
    return rcpp_result_gen;
END_RCPP
}
// test_extinction
Rcpp::StringMatrix test_extinction(const Rcpp::StringMatrix& island_spec_R, int focal_trait);
RcppExport SEXP _TRAISIERCPP_test_extinction(SEXP island_spec_RSEXP, SEXP focal_traitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringMatrix& >::type island_spec_R(island_spec_RSEXP);
    Rcpp::traits::input_parameter< int >::type focal_trait(focal_traitSEXP);
    rcpp_result_gen = Rcpp::wrap(test_extinction(island_spec_R, focal_trait));
    return rcpp_result_gen;
END_RCPP
}
// test_execute_extinction
Rcpp::StringMatrix test_execute_extinction(const Rcpp::StringMatrix& island_spec_R, int index);
RcppExport SEXP _TRAISIERCPP_test_execute_extinction(SEXP island_spec_RSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringMatrix& >::type island_spec_R(island_spec_RSEXP);
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(test_execute_extinction(island_spec_R, index));
    return rcpp_result_gen;
END_RCPP
}
// test_anagenesis
Rcpp::StringMatrix test_anagenesis(Rcpp::StringMatrix& island_spec_R, int maxspecID, int focal_trait);
RcppExport SEXP _TRAISIERCPP_test_anagenesis(SEXP island_spec_RSEXP, SEXP maxspecIDSEXP, SEXP focal_traitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type island_spec_R(island_spec_RSEXP);
    Rcpp::traits::input_parameter< int >::type maxspecID(maxspecIDSEXP);
    Rcpp::traits::input_parameter< int >::type focal_trait(focal_traitSEXP);
    rcpp_result_gen = Rcpp::wrap(test_anagenesis(island_spec_R, maxspecID, focal_trait));
    return rcpp_result_gen;
END_RCPP
}
// test_cladogenesis
Rcpp::StringMatrix test_cladogenesis(Rcpp::StringMatrix& island_spec_R, double timeval, int maxspecID, int focal_trait);
RcppExport SEXP _TRAISIERCPP_test_cladogenesis(SEXP island_spec_RSEXP, SEXP timevalSEXP, SEXP maxspecIDSEXP, SEXP focal_traitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type island_spec_R(island_spec_RSEXP);
    Rcpp::traits::input_parameter< double >::type timeval(timevalSEXP);
    Rcpp::traits::input_parameter< int >::type maxspecID(maxspecIDSEXP);
    Rcpp::traits::input_parameter< int >::type focal_trait(focal_traitSEXP);
    rcpp_result_gen = Rcpp::wrap(test_cladogenesis(island_spec_R, timeval, maxspecID, focal_trait));
    return rcpp_result_gen;
END_RCPP
}
// test_transition
Rcpp::StringMatrix test_transition(Rcpp::StringMatrix& island_spec_R, int focal_trait);
RcppExport SEXP _TRAISIERCPP_test_transition(SEXP island_spec_RSEXP, SEXP focal_traitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type island_spec_R(island_spec_RSEXP);
    Rcpp::traits::input_parameter< int >::type focal_trait(focal_traitSEXP);
    rcpp_result_gen = Rcpp::wrap(test_transition(island_spec_R, focal_trait));
    return rcpp_result_gen;
END_RCPP
}
// test_sample_spec
int test_sample_spec(const Rcpp::NumericVector& mainland_spec_R, int M2);
RcppExport SEXP _TRAISIERCPP_test_sample_spec(SEXP mainland_spec_RSEXP, SEXP M2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mainland_spec_R(mainland_spec_RSEXP);
    Rcpp::traits::input_parameter< int >::type M2(M2SEXP);
    rcpp_result_gen = Rcpp::wrap(test_sample_spec(mainland_spec_R, M2));
    return rcpp_result_gen;
END_RCPP
}
// test_immigration2
Rcpp::StringMatrix test_immigration2(Rcpp::StringMatrix& island_spec_R, const Rcpp::NumericVector& mainland_spec_R, double timeval, int M2);
RcppExport SEXP _TRAISIERCPP_test_immigration2(SEXP island_spec_RSEXP, SEXP mainland_spec_RSEXP, SEXP timevalSEXP, SEXP M2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type island_spec_R(island_spec_RSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mainland_spec_R(mainland_spec_RSEXP);
    Rcpp::traits::input_parameter< double >::type timeval(timevalSEXP);
    Rcpp::traits::input_parameter< int >::type M2(M2SEXP);
    rcpp_result_gen = Rcpp::wrap(test_immigration2(island_spec_R, mainland_spec_R, timeval, M2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TRAISIERCPP_sample_event", (DL_FUNC) &_TRAISIERCPP_sample_event, 1},
    {"_TRAISIERCPP_execute_time_loop", (DL_FUNC) &_TRAISIERCPP_execute_time_loop, 9},
    {"_TRAISIERCPP_test_get_immig_rate", (DL_FUNC) &_TRAISIERCPP_test_get_immig_rate, 7},
    {"_TRAISIERCPP_test_get_ext_rate", (DL_FUNC) &_TRAISIERCPP_test_get_ext_rate, 6},
    {"_TRAISIERCPP_test_get_ana_rate", (DL_FUNC) &_TRAISIERCPP_test_get_ana_rate, 5},
    {"_TRAISIERCPP_test_get_clado_rate", (DL_FUNC) &_TRAISIERCPP_test_get_clado_rate, 7},
    {"_TRAISIERCPP_test_get_trans_rate", (DL_FUNC) &_TRAISIERCPP_test_get_trans_rate, 4},
    {"_TRAISIERCPP_test_update_rates", (DL_FUNC) &_TRAISIERCPP_test_update_rates, 12},
    {"_TRAISIERCPP_test_draw_prop", (DL_FUNC) &_TRAISIERCPP_test_draw_prop, 1},
    {"_TRAISIERCPP_test_immigration", (DL_FUNC) &_TRAISIERCPP_test_immigration, 3},
    {"_TRAISIERCPP_test_extinction", (DL_FUNC) &_TRAISIERCPP_test_extinction, 2},
    {"_TRAISIERCPP_test_execute_extinction", (DL_FUNC) &_TRAISIERCPP_test_execute_extinction, 2},
    {"_TRAISIERCPP_test_anagenesis", (DL_FUNC) &_TRAISIERCPP_test_anagenesis, 3},
    {"_TRAISIERCPP_test_cladogenesis", (DL_FUNC) &_TRAISIERCPP_test_cladogenesis, 4},
    {"_TRAISIERCPP_test_transition", (DL_FUNC) &_TRAISIERCPP_test_transition, 2},
    {"_TRAISIERCPP_test_sample_spec", (DL_FUNC) &_TRAISIERCPP_test_sample_spec, 2},
    {"_TRAISIERCPP_test_immigration2", (DL_FUNC) &_TRAISIERCPP_test_immigration2, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_TRAISIERCPP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
