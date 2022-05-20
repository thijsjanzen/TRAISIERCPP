#' @keywords internal
check_param <- function(param, param_name) {
  if (is.null(param)) {
    msg <- paste(param_name, "was null")
    warning(msg)
  }
}


#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
#' @export
DAISIE_sim_core_trait_dep_cpp <- function(
    time,
    mainland_n,
    pars,
    island_ontogeny = 0,
    sea_level = 0,
    hyper_pars,
    area_pars,
    extcutoff = 1000,
    trait_pars = NULL
) {

  #### Initialization ####
  timeval <- 0
  total_time <- time
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  sea_level <- translate_sea_level(sea_level)

  if (is.null(trait_pars)) {
    stop("A second set of rates should be contain considering two trait states.
         If only one state,run DAISIE_sim_cr instead.")
  }

  if (pars[4] == 0 && trait_pars$immig_rate2 == 0) {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }

  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]

  results <- execute_time_loop(timeval,
                               total_time,
                               gam,
                               laa,
                               lac,
                               mu,
                               K,
                               mainland_n,
                               trait_pars)

  stt_table = results$stt_table
  colnames(stt_table) <- c("Time","nI","nA","nC","nI2","nA2","nC2")
  island_spec = results$island_spec

  island <- DAISIE_create_island(
    stt_table = stt_table,
    total_time = total_time,
    island_spec = island_spec,
    mainland_n = mainland_n,
    trait_pars = trait_pars)
  return(island)
}
