library(DAISIE)
library(TraisieABC)

obs_sim_pars <- c(2.0,  # lac   // cladogenesis
                  0.0,  # mu    // extinction
                  0.001, # gam   // immigration rate 1
                  0.0,  # laa   // anagenesis
                  0.3,  # clado rate 2
                  0.0,  # ext_rate2
                  0.001,  # immig_rate2
                  0.0,   # ana_rate2
                  0.1,   # trans_rate
                  0.1)   # trans_rate2
obs_sim <- TraisieABC::get_TraiSIE_sim(parameters = obs_sim_pars,
                           K = Inf,
                           replicates = 1,
                           verbose = FALSE)

library(TRAISIERCPP)
#profvis::profvis({
obs_sim2 <- TRAISIERCPP::get_TraiSIE_sim_cpp(parameters = obs_sim_pars,
                                             K = Inf,
                                             replicates = 1)
#})

microbenchmark::microbenchmark(TraisieABC::get_TraiSIE_sim(parameters = obs_sim_pars,
                                           K = Inf,
                                           replicates = 1),
               TRAISIERCPP::get_TraiSIE_sim_cpp(parameters = obs_sim_pars,
                                                K = Inf,
                                                replicates = 1),
               time = 10L)


parameters = obs_sim_pars
K = Inf
replicates = 1
time = 2
M = 500
pars = c(parameters[1], parameters[2], K, parameters[3], parameters[4])
replicates = 1
sample_freq  = Inf
plot_sims = TRUE
cond = 1
verbose = FALSE
trait_pars = DAISIE::create_trait_pars(clado_rate2 = parameters[5],
                                       ext_rate2 = parameters[6],
                                       immig_rate2 = parameters[7],
                                       ana_rate2 = parameters[8],
                                       trans_rate = parameters[9],
                                       trans_rate2 = parameters[10],
                                       M2 = 500)

divdepmodel = "CS"
sample_freq = 25
plot_sims = FALSE
island_ontogeny = "const"
sea_level = "const"
hyper_pars = TRAISIERCPP::create_hyper_pars(d = 0, x = 0)
area_pars = TRAISIERCPP::create_area_pars(
  max_area = 1,
  current_area = 1,
  proportional_peak_t = 0,
  total_island_age = 0,
  sea_level_amplitude = 0,
  sea_level_frequency = 0,
  island_gradient_angle = 0)
extcutoff = 1000
trait_pars_onecolonize <- create_trait_pars(
  trans_rate = trait_pars$trans_rate,
  immig_rate2 = trait_pars$immig_rate2,
  ext_rate2 = trait_pars$ext_rate2,
  ana_rate2 = trait_pars$ana_rate2,
  clado_rate2 = trait_pars$clado_rate2,
  trans_rate2 = trait_pars$trans_rate2,
  M2 = 0)
total_time <- time
time = total_time
  mainland_n = 1
  pars = pars
  island_ontogeny = island_ontogeny
  sea_level = sea_level
  hyper_pars = hyper_pars
  area_pars = area_pars
  extcutoff = extcutoff
  trait_pars = trait_pars_onecolonize
