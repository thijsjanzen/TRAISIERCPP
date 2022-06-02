#' Simulation function to create simulated data as observed data in ABC.
#'
#' @param parameters A vector for CES rates.
#' @param K Carrying capacity, Inf for diverdity-independent models.
#' @param replicates The number of replicates(islands) for TraiSIE simulation.
#'
#' @return A list contains simulated islands
#' @author Shu Xie
#' @export
#' @rawNamespace useDynLib(TRAISIERCPP)
#' @rawNamespace import(Rcpp)
#' @rawNamespace import(RcppArmadillo)
get_TraiSIE_sim_cpp <- function(parameters, K, replicates){
  sim <- list()
  for (j in seq_len(replicates)) {
    sim[[j]] <- TRAISIERCPP::DAISIE_sim_trait_dep_cpp(
      time = 2,
      M = 500,
      pars = c(parameters[1], parameters[2], K, parameters[3], parameters[4]),
      replicates = 1,
      sample_freq  = Inf,
      plot_sims = FALSE,
      cond = 1,
      verbose = FALSE,
      trait_pars = DAISIE::create_trait_pars(clado_rate2 = parameters[5],
                                             ext_rate2 = parameters[6],
                                             immig_rate2 = parameters[7],
                                             ana_rate2 = parameters[8],
                                             trans_rate = parameters[9],
                                             trans_rate2 = parameters[10],
                                             M2 = 500)
    )
  }
  return(sim)
}
