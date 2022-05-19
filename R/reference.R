#' Calculates algorithm rates
#' @description Internal function that updates the all the rates and
#' max extinction horizon at time t.
#' @family rate calculations
#'
#' @inheritParams default_params_doc
#'
#' @seealso \code{\link{update_max_rates}()}
#' @export
#' @return a named list with the updated effective rates.
DAISIE_update_rates_trait <- function(timeval,
                               total_time,
                               gam,
                               laa,
                               lac,
                               mu,
                               hyper_pars = hyper_pars,
                               area_pars = NULL,
                               peak = NULL,
                               island_ontogeny = NULL,
                               sea_level = NULL,
                               extcutoff,
                               K,
                               num_spec,
                               num_immigrants,
                               mainland_n,
                               trait_pars = NULL,
                               island_spec = NULL) {
  # Function to calculate rates at time = timeval. Returns list with each rate.

  A <- 1

  immig_rate <- DAISIE_get_immig_rate(
    gam = gam,
    A = A,
    num_spec = num_spec,
    K = K,
    mainland_n = mainland_n,
    trait_pars = trait_pars,
    island_spec = island_spec
  )

  ext_rate <- DAISIE_get_ext_rate(
    mu = mu,
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    A = A,
    trait_pars = trait_pars,
    island_spec = island_spec
  )

  ana_rate <- DAISIE_get_ana_rate(
    laa = laa,
    num_immigrants = num_immigrants,
    trait_pars = trait_pars,
    island_spec = island_spec
  )
  clado_rate <- DAISIE_get_clado_rate(
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    K = K,
    A = A,
    trait_pars = trait_pars,
    island_spec = island_spec
  )



  trans_rate <- DAISIE_get_trans_rate(trait_pars = trait_pars,
                               island_spec = island_spec)


  rates <- list(
    immig_rate = immig_rate$immig_rate1,
    ext_rate = ext_rate$ext_rate1,
    ana_rate = ana_rate$ana_rate1,
    clado_rate = clado_rate$clado_rate1,
    trans_rate = trans_rate$trans_rate1,
    immig_rate2 = immig_rate$immig_rate2,
    ext_rate2 = ext_rate$ext_rate2,
    ana_rate2 = ana_rate$ana_rate2,
    clado_rate2 = clado_rate$clado_rate2,
    trans_rate2 = trans_rate$trans_rate2,
    M2 = trait_pars$M2)

  return(rates)
}


#' Function to describe changes in extinction rate through time.
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @family rate calculations
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784
#' (2014): 20133227.
#' @author Pedro Neves, Joshua Lambert, Shu Xie
DAISIE_get_ext_rate <- function(mu,
                         hyper_pars,
                         extcutoff = 1000,
                         num_spec,
                         A,
                         trait_pars = NULL,
                         island_spec = NULL) {

  x <- hyper_pars$x
  if (is.null(trait_pars)) {
    ext_rate <- max(0, mu * (A ^ -x) * num_spec, na.rm = TRUE)
    ext_rate <- min(ext_rate, extcutoff, na.rm = TRUE)
    # testit::assert(ext_rate >= 0)
    return(ext_rate)
  } else {   ##species have two states
    if (is.matrix(island_spec) || is.null(island_spec)) {
      num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
      num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    }
    ext_rate1 <- mu * num_spec_trait1
    ext_rate2 <- trait_pars$ext_rate2 * num_spec_trait2
    ext_list <- list(ext_rate1 = ext_rate1,
                     ext_rate2 = ext_rate2)
    return(ext_list)
  }
}

#' Calculate anagenesis rate
#' @description Internal function.
#' Calculates the anagenesis rate given the current number of
#' immigrant species and the per capita rate.
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert, Shu Xie
DAISIE_get_ana_rate <- function(laa,
                         num_immigrants,
                         island_spec = NULL,
                         trait_pars = NULL) {

  if (is.null(trait_pars)) {
    ana_rate <- laa * num_immigrants
    return(ana_rate)
  }else{
    ana_rate1 = laa * length(intersect(which(island_spec[,4] == "I"),
                                       which(island_spec[,8] == "1")))
    ana_rate2 = trait_pars$ana_rate2 * length(intersect(which(island_spec[,4] == "I"),
                                                        which(island_spec[,8] == "2")))
    ana_list <- list(ana_rate1 = ana_rate1,
                     ana_rate2 = ana_rate2)
    return(ana_list)
  }
}

#' Calculate cladogenesis rate
#' @description Internal function.
#' Calculates the cladogenesis rate given the current number of
#' species in the system, the carrying capacity and the per capita cladogenesis
#' rate
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @author Pedro Neves, Joshua Lambert, Shu Xie
DAISIE_get_clado_rate <- function(lac,
                           hyper_pars,
                           num_spec,
                           K,
                           A,
                           trait_pars = NULL,
                           island_spec = NULL) {
  # testit::assert(are_hyper_pars(hyper_pars))

  d <- hyper_pars$d
  if (is.null(trait_pars)) {
    clado_rate <- max(
      0, lac * num_spec * (A ^ d) * (1 - num_spec / (K * A)), na.rm = TRUE
    )
    # testit::assert(clado_rate >= 0)
    # testit::assert(is.numeric(clado_rate))
    return(clado_rate)
  }else{
    num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
    num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    clado_rate1 <- max(
      0, lac * num_spec_trait1 * (1 - num_spec / K),
      na.rm = TRUE)
    clado_rate2 <- max(
      0, trait_pars$clado_rate2 * num_spec_trait2 * (1 - num_spec / K),
      na.rm = TRUE
    )
    # testit::assert(clado_rate1 >= 0)
    # testit::assert(clado_rate2 >= 0)
    # testit::assert(is.numeric(clado_rate1))
    # testit::assert(is.numeric(clado_rate2))
    clado_list <- list(clado_rate1 = clado_rate1,
                       clado_rate2 = clado_rate2)
    return(clado_list)
  }
}

#' Calculate immigration rate
#' @description Internal function.
#' Calculates the immigration rate given the current number of
#' species in the system, the carrying capacity
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
DAISIE_get_immig_rate <- function(gam,
                           A,
                           num_spec,
                           K,
                           mainland_n,
                           trait_pars = NULL,
                           island_spec = NULL) {

  if (is.null(trait_pars)) {
    immig_rate <- max(c(mainland_n * gam * (1 - (num_spec / (A * K))),
                        0), na.rm = TRUE)
    return(immig_rate)
  } else {
    mainland_n2 <- trait_pars$M2
    gam2 <- trait_pars$immig_rate2
    immig_rate1 <- max(c(mainland_n * gam * (1 - (num_spec / (A * K))),
                         0), na.rm = TRUE)
    immig_rate2 <- max(c(mainland_n2 * gam2 * (1 - (num_spec / (A * K))),
                         0), na.rm = TRUE)
    immig_list <- list(immig_rate1 = immig_rate1,
                       immig_rate2 = immig_rate2)
    return(immig_list)
  }
}

#' Calculate transition rate
#' @description Internal function.
#' Calculates the transition rate given the current number of
#' immigrant species and the per capita rate.
#' @param trait_pars A named list containing diversification rates considering
#' two trait states created by \code{\link{create_trait_pars}}:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param island_spec Matrix with current state of simulation containing number
#' of species.
#' @export
#' @family rates calculation
DAISIE_get_trans_rate <- function(trait_pars,
                           island_spec){

  # Make function accept island_spec matrix or numeric
  if (is.matrix(island_spec) || is.null(island_spec)) {
    num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
    num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
  }
  trans_rate1 <- trait_pars$trans_rate * num_spec_trait1
  trans_rate2 <- trait_pars$trans_rate2 * num_spec_trait2
  # testit::assert(is.numeric(trans_rate1))
  # testit::assert(trans_rate1 >= 0)
  # testit::assert(is.numeric(trans_rate2))
  # testit::assert(trans_rate2 >= 0)
  trans_list <- list(trans_rate1 = trans_rate1,
                     trans_rate2 = trans_rate2)
  return(trans_list)

}

#' Calculates when the next timestep will be.
#'
#' @param timeval current time of simulation
#' @param max_rates named list of max rates as returned by
#' \code{\link{update_rates}}.
#'
#' @return named list with numeric vector containing the time of the next
#' timestep and the change in time.
#'
#' @export
#'
#' @author Joshua Lambert, Pedro Neves, Shu Xie
DAISIE_calc_next_timeval <- function(max_rates, timeval) {
  # testit::assert(timeval >= 0)

  if (length(max_rates) == 4) {   ## no considering about two trait states
    totalrate <- max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]]
  } else {
    totalrate <- max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]] +
      max_rates[[5]] + max_rates[[6]] + max_rates[[7]] + max_rates[[8]] +
      max_rates[[9]] + max_rates[[10]]
  }
  dt <- stats::rexp(1, totalrate)
  timeval <- timeval + dt
  return(list(timeval = timeval, dt = dt))
}

