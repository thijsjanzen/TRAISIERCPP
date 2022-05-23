#' Create list of hyperparameters
#'
#' @inheritParams default_params_doc
#'
#' @return Named list with hyperparameters
#' @export
#' @author Pedro Neves, Joshua Lambert
create_hyper_pars <- function(d, x) {
  testit::assert(d >= 0.0)
  testit::assert(is.numeric(x))
  list(
    d = d,
    x = x
  )
}

#' Create named list of area parameters
#'
#' @inheritParams default_params_doc
#'
#' @return list of numerical values containing area and sea level parameters
#' for island ontogeny simulation
#' @export
#' @author Richel J.C Bilderbeek, Joshua Lambert, Pedro Neves
#'
#'
create_area_pars <- function(max_area = 1,
                             current_area = 1,
                             proportional_peak_t = 0.5,
                             total_island_age = 5,
                             sea_level_amplitude = 1,
                             sea_level_frequency = 1,
                             island_gradient_angle = 0.45) {
  testit::assert(max_area > 0.0)
  testit::assert(current_area > 0.0)
  testit::assert(proportional_peak_t >= 0.0)
  testit::assert(proportional_peak_t <= 1.0)
  testit::assert(total_island_age >= 0.0)
  testit::assert(sea_level_amplitude >= 0.0)
  testit::assert(sea_level_frequency >= 0.0)
  testit::assert(island_gradient_angle >= 0)
  testit::assert(island_gradient_angle <= 90)
  list(max_area = max_area,
       current_area = current_area,
       proportional_peak_t = proportional_peak_t,
       total_island_age = total_island_age,
       sea_level_amplitude = sea_level_amplitude,
       sea_level_frequency = sea_level_frequency,
       island_gradient_angle = island_gradient_angle)
}

#' Converts simulation output into island output
#'
#' @inheritParams default_params_doc
#'
#' @return list with the island information, composed stt table,
#' branching times of extant species, status of species on
#' the island and number of missing species.
#' @keywords internal
DAISIE_create_island <- function(stt_table,
                                 total_time,
                                 island_spec,
                                 mainland_n,
                                 trait_pars = NULL) {

  if (!is.null(trait_pars)) {
    return(
      DAISIE_create_island_trait(
        stt_table = stt_table,
        total_time = total_time,
        island_spec = island_spec,
        mainland_n = mainland_n,
        trait_pars = trait_pars
      )
    )
  }
  ### if there are no species on the island branching_times = island_age,
  ### stac = 0, missing_species = 0
  if (length(island_spec[, 1]) == 0) {
    island <- list(stt_table = stt_table,
                   branching_times = total_time,
                   stac = 0,
                   missing_species = 0)
  } else {
    cnames <- c("Species",
                "Mainland Ancestor",
                "Colonisation time (BP)",
                "Species type",
                "branch_code",
                "branching time (BP)",
                "Anagenetic_origin")
    colnames(island_spec) <- cnames
    ### set ages as counting backwards from present
    island_spec[, "branching time (BP)"] <- total_time -
      as.numeric(island_spec[, "branching time (BP)"])
    island_spec[, "Colonisation time (BP)"] <- total_time -
      as.numeric(island_spec[, "Colonisation time (BP)"])
    if (mainland_n == 1) {
      island <- DAISIE_ONEcolonist(total_time,
                                   island_spec,
                                   stt_table)
    } else if (mainland_n > 1) {
      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(
        island_spec[, "Mainland Ancestor"])))
      number_colonists_present <- length(colonists_present)
      island_clades_info <- list()
      for (i in 1:number_colonists_present) {
        subset_island <- island_spec[which(island_spec[, "Mainland Ancestor"] ==
                                             colonists_present[i]), ]
        if (!is.matrix(subset_island)) {
          subset_island <- rbind(subset_island[1:7])
          colnames(subset_island) <- cnames
        }
        island_clades_info[[i]] <- DAISIE_ONEcolonist(
          total_time,
          island_spec = subset_island,
          stt_table = NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
      island <- list(stt_table = stt_table,
                     taxon_list = island_clades_info)
    }
  }
  return(island)
}

DAISIE_create_island_trait <- function(stt_table,
                                       total_time,
                                       island_spec,
                                       mainland_n,
                                       trait_pars){

  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0
  if (length(island_spec[,1]) == 0) {
    island <- list(stt_table = stt_table,
                   branching_times = total_time,
                   stac = 0,
                   missing_species = 0)

  } else {
    cnames <- c("Species",
                "Mainland Ancestor",
                "Colonisation time (BP)",
                "Species type",
                "branch_code",
                "branching time (BP)",
                "Anagenetic_origin",
                "trait_state")

    colnames(island_spec) <- cnames

    ### set ages as counting backwards from present
    island_spec[, "branching time (BP)"] <- total_time - as.numeric(island_spec[, "branching time (BP)"])
    island_spec[, "Colonisation time (BP)"] <- total_time - as.numeric(island_spec[, "Colonisation time (BP)"])

    mainland_ntotal = mainland_n + trait_pars$M2

    if (mainland_ntotal == 1) {
      island <- DAISIE_ONEcolonist(total_time,
                                   island_spec,
                                   stt_table)


    } else if (mainland_ntotal > 1) {

      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(island_spec[, 'Mainland Ancestor'])))
      number_colonists_present <- length(colonists_present)

      island_clades_info <- list()
      for (i in 1:number_colonists_present) {
        subset_island <- island_spec[which(island_spec[, "Mainland Ancestor"] ==
                                             colonists_present[i]), ]
        if (!is.matrix(subset_island)) {
          subset_island <- rbind(subset_island[1:8])
          colnames(subset_island) <- cnames
        }
        island_clades_info[[i]] <- DAISIE_ONEcolonist(
          total_time,
          island_spec = subset_island,
          stt_table = NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
      island <- list(stt_table = stt_table,
                     taxon_list = island_clades_info)
    }
  }
  return(island)
}

#' Convert intermediate output to final simulation output
#'
#' @inheritParams default_params_doc
#'
#' @return a list with these elements:
#' \itemize{
#'   \item{[1]: \code{stt_table}, the same stt_table as put in.}
#'   \item{[2]: \code{branching_times}, a sorted numeric vector, as required
#'     by the ML estimation functions. The first element always refers to
#'     the island age. Subsequent elements refer to colonisation, speciation and
#'     recolonisation times. The most recent recolonisation time, if any is
#'     always omitted to approximate simulation results to the mathematical
#'     formulation of the likelihood functions used for MLE.}
#'   \item{[3]: \code{stac}, status of colonist. In this function it can be
#'     returned as either 2, 4 or 3. If \code{stac} is 2, then there is only one
#'     independent colonisation present on the island and the extant species are
#'     endemic. If stac is 4, then only a singleton endemic is present at the
#'     present. If stac is 3, then recolonisation occurred, and more than one
#'     colonising lineage.}
#'   \item{[4]: \code{missing_species}, a numeric value with the number of
#'     missing species, that is, species not sampled in the phylogeny but
#'     present on the island. As this code only runs for simulation models,
#'     here \code{missing_species} is always set to 0.}
#'   \item{[5]:
#'   \code{all_colonisations}, on recolonising lineages only. It is comprised of
#'     \code{$event_times} and \code{$species_type}:
#'     \describe{
#'       \item{\code{$event_times}}{ordered numeric vectors containing all
#'       events for each extant recolonising lineage. This includes all
#'       colonisation and branching times. Each vector pertains to one
#'       colonising lineage.}
#'       \item{\code{$species_type}}{a string. Can be \code{"A"}, \code{"C"} or
#'       \code{"I"} depending on whether the extant clade is of anagenetic,
#'       cladogenetic or immigrant origin, respectively.}
#'     }
#'   }
#' }
#' @keywords internal
DAISIE_ONEcolonist <- function(time,
                               island_spec,
                               stt_table) {
  ### number of independent colonisations
  uniquecolonisation <- as.numeric(unique(
    island_spec[, "Colonisation time (BP)"]))
  number_colonisations <- length(uniquecolonisation)
  ### if there is only one independent colonisation - anagenetic and
  ### cladogenetic species are classed as stac=2; immigrant classed as stac=4:
  if (number_colonisations == 1) {
    if (island_spec[1, "Species type"] == "I") {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(
                            time,
                            as.numeric(island_spec[1, "Colonisation time (BP)"])
                          ),
                          stac = 4,
                          missing_species = 0)
    }
    if (island_spec[1, "Species type"] == "A") {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(
                            time,
                            as.numeric(island_spec[1, "Colonisation time (BP)"])
                          ),
                          stac = 2,
                          missing_species = 0)
    }
    if (island_spec[1, "Species type"] == "C") {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(
                            time,
                            sort(
                              as.numeric(island_spec[, "branching time (BP)"]),
                              decreasing = TRUE
                            )
                          ),
                          stac = 2,
                          missing_species = 0)
    }
  }

  ### if there are two or more independent colonisations, all species are
  ### classed as stac=3 and put within same list item:
  if (number_colonisations > 1) {
    descendants <- list(stt_table = stt_table,
                        branching_times = NA,
                        stac = 3,
                        missing_species = 0,
                        all_colonisations = list())

    # Get branching and colonisation times
    btimes_all_clado_desc <- rev(
      sort(as.numeric(island_spec[, "branching time (BP)"]))
    )
    col_times <- sort(
      unique(as.numeric(island_spec[, "Colonisation time (BP)"])),
      decreasing = TRUE
    )

    # If there are endemic descendants find youngest col time
    if (length(btimes_all_clado_desc) != 0) {
      # Ensure all col_times are in b_times at this point.
      # Covers cases of one recolonization followed by cladogenesis and
      # potential extinction
      if (any(!(col_times %in% btimes_all_clado_desc))) {
        miss_col_time <- which(!(col_times %in% btimes_all_clado_desc))
        btimes_all_clado_desc <- sort(
          c(btimes_all_clado_desc, col_times[miss_col_time]),
          decreasing = TRUE
        )
      }
      youngest_col_time <- min(col_times)
      i_youngest_col_btimes <- which(btimes_all_clado_desc == youngest_col_time)

      # Remove youngest col time in branching times
      testit::assert(youngest_col_time %in% btimes_all_clado_desc)
      btimes_all_clado_desc <- btimes_all_clado_desc[-i_youngest_col_btimes]

      descendants$branching_times <- c(time, btimes_all_clado_desc)
      testit::assert(!(youngest_col_time %in% btimes_all_clado_desc))

      # If no cladogenetic species is present, remove the youngest col time
    } else if (length(btimes_all_clado_desc) == 0) {
      youngest_col_time <- min(col_times)
      i_youngest_col_time <- which(col_times == youngest_col_time)
      col_times <- col_times[-i_youngest_col_time]

      descendants$branching_times <- c(time, col_times)
    }


    # all_colonisations section
    uniquecol <- sort(as.numeric(
      unique(island_spec[, "Colonisation time (BP)"])), decreasing = TRUE
    )
    for (i in seq_along(uniquecol)) {
      descendants$all_colonisations[[i]] <- list(
        event_times = NA,
        species_type = NA
      )

      samecolonisation <- which(as.numeric(
        island_spec[, "Colonisation time (BP)"]) == uniquecol[i]
      )

      if (island_spec[samecolonisation[1], "Species type"] == "I") {
        descendants$all_colonisations[[i]]$event_times <- as.numeric(
          c(time,island_spec[samecolonisation, "Colonisation time (BP)"])
        )
        descendants$all_colonisations[[i]]$species_type <- "I"
      }

      if (island_spec[samecolonisation[1], "Species type"] == "A") {
        descendants$all_colonisations[[i]]$event_times <- as.numeric(
          c(time, island_spec[samecolonisation, "Colonisation time (BP)"])
        )
        descendants$all_colonisations[[i]]$species_type <- "A"
      }

      if (island_spec[samecolonisation[1], "Species type"] == "C") {
        descendants$all_colonisations[[i]]$event_times <-
          sort(c(time, as.numeric(
            island_spec[samecolonisation, "branching time (BP)"]
          )), decreasing = TRUE)
        descendants$all_colonisations[[i]]$species_type <- "C"
      }
    }
  }
  return(descendants)
}


#' Create named list of trait state parameters
#'
#' @param trans_rate   A numeric with the per capita transition rate with state1
#' @param immig_rate2  A numeric with the per capita immigration rate with state2
#' @param ext_rate2    A numeric with the per capita extinction rate with state2
#' @param ana_rate2    A numeric with the per capita anagenesis rate with state2
#' @param clado_rate2  A numeric with the per capita cladogenesis rate with state2
#' @param trans_rate2  A numeric with the per capita transition rate with state2
#' @param M2           A numeric with the number of species with trait state 2 on mainland
#'
#' @return list of numerical values containing trait state parameters
#' @export
#'

create_trait_pars <- function(trans_rate,
                              immig_rate2,
                              ext_rate2,
                              ana_rate2,
                              clado_rate2,
                              trans_rate2,
                              M2) {
  list(trans_rate = trans_rate,
       immig_rate2 = immig_rate2,
       ext_rate2 = ext_rate2,
       ana_rate2 = ana_rate2,
       clado_rate2 = clado_rate2,
       trans_rate2 = trans_rate2,
       M2 = M2)
}

#' Wrapper function around for the DAISIE_format_CS_full_stt and
#' DAISIE_format_CS_sampled_stt
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
#' @return List with CS DAISIE simulation output
DAISIE_format_CS <- function(island_replicates,
                             time,
                             M,
                             sample_freq = 25,
                             verbose = TRUE,
                             trait_pars = NULL) {
  total_time <- time
  testit::assert(
    !is.na(sample_freq) && !is.null(sample_freq) && sample_freq >= 1
  )
  if (is.infinite(sample_freq)) {
    several_islands <- DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = total_time,
      M = M,
      verbose = verbose,
      trait_pars = trait_pars
    )
  } else {
    several_islands <- DAISIE_format_CS_sampled_stt(
      island_replicates = island_replicates,
      time = total_time,
      sample_freq = sample_freq,
      M = M,
      verbose = verbose,
      trait_pars = trait_pars
    )
  }
  return(several_islands)
}

#' Formats clade-specific simulation output into standard
#' DAISIE list output with complete STT table
#'
#' @inheritParams default_params_doc
#'
#' @return List with CS DAISIE simulation output
DAISIE_format_CS_full_stt <- function(island_replicates,
                                      time,
                                      M,
                                      verbose = TRUE,
                                      trait_pars = NULL
) {
  total_time <- time
  several_islands <- list()
  for (rep in seq_along(island_replicates)) {
    full_list <- island_replicates[[rep]]
    stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
    number_not_present <- length(which(stac_vec == 0))
    present <- which(stac_vec != 0)
    number_present <- length(present)
    type_vec <- unlist(full_list)[which(names(unlist(full_list)) == "type1or2")]
    prop_type2_pool <- length(which(type_vec == 2)) / M
    number_type2_cols <- length(which(match(which(stac_vec != 0),
                                            which(type_vec == 2)) > 0))
    number_type1_cols <- number_present - number_type2_cols
    island_list <- list()
    for (i in 1:(number_present + 1)) {
      island_list[[i]] <- list()
    }
    ### all species
    stt_list <- list()
    if (is.null(trait_pars)) {
      for (i in 1:M) {
        stt_list[[i]] <- full_list[[i]]$stt_table
      }
    } else {
      for (i in 1:(M + trait_pars$M2)) {
        stt_list[[i]] <- full_list[[i]]$stt_table
      }
    }

    #### Keep full STT ####
    stt_all <- create_full_CS_stt(
      stt_list = stt_list,
      stac_vec = stac_vec,
      total_time = total_time,
      trait_pars = trait_pars
    )
    # ####  two trait states####
    if(!is.null(trait_pars)){
      stt_combined <- stt_all$stt_all
      stt_two_states <- stt_all$stt_two_states
      stt_all <- stt_combined

      immig_spec <- c()
      ana_spec <- c()
      immig_spec2 <- c()
      ana_spec2 <- c()
      for (i in 1:(M + trait_pars$M2)) {
        immig_spec[i] <- sum(full_list[[i]]$stt_table[1, 2])
        ana_spec[i] <- sum(full_list[[i]]$stt_table[1, 3])
        immig_spec2[i] <- sum(full_list[[i]]$stt_table[1, 5])
        ana_spec2[i] <- sum(full_list[[i]]$stt_table[1, 6])
      }
      immig_spec <- sum(immig_spec)
      ana_spec <- sum(ana_spec)
      immig_spec2 <- sum(immig_spec2)
      ana_spec2 <- sum(ana_spec2)
      init_present <- immig_spec + ana_spec + immig_spec2 + ana_spec2
      stt_two_states[1, 2:8] <- c(immig_spec, ana_spec, 0, immig_spec2, ana_spec2, 0, init_present)
      stt_all[1, 2:5] <- c(immig_spec + immig_spec2, ana_spec + ana_spec2, 0, init_present)
    }else{
      #### Oceanic vs nonoceanic ####

      immig_spec <- c()
      ana_spec <- c()
      for (i in 1:M) {
        immig_spec[i] <- sum(full_list[[i]]$stt_table[1, 2])
        ana_spec[i] <- sum(full_list[[i]]$stt_table[1, 3])
      }
      immig_spec <- sum(immig_spec)
      ana_spec <- sum(ana_spec)
      init_present <- immig_spec + ana_spec
      stt_all[1, 2:5] <- c(immig_spec, ana_spec, 0, init_present)
    }


    #### 2 type ####
    if (number_type2_cols > 0) {
      # Type 1
      stt_list_type1 <- list()
      for (i in 1:max(which(type_vec == 1))) {
        stt_list_type1[[i]] <- full_list[[i]]$stt_table
      }
      stt_type1 <- create_full_CS_stt(
        stt_list = stt_list_type1,
        stac_vec = stac_vec,
        total_time = total_time,
        trait_pars = trait_pars
      )


      ######################################################### list type2
      type2len <- length(which(type_vec == 2))
      stt_list_type2 <- list()
      for (i in 1:type2len) {
        stt_list_type2[[i]] <- full_list[[which(type_vec == 2)[i]]]$stt_table
      }

      stt_type2 <- create_full_CS_stt(
        stt_list = stt_list_type2,
        stac_vec = stac_vec,
        total_time = total_time,
        trait_pars = trait_pars
      )

      island_list[[1]] <- list(island_age = total_time,
                               not_present_type1 = DDD::roundn(
                                 M * (1 - prop_type2_pool)) -
                                 (number_type1_cols),
                               not_present_type2 = DDD::roundn(
                                 M * prop_type2_pool) - number_type2_cols,
                               stt_all = stt_all,
                               stt_type1 = stt_type1,
                               stt_type2 = stt_type2)
    } else if(!is.null(trait_pars)){
      island_list[[1]] <- list(island_age = total_time,
                               not_present = number_not_present,
                               stt_all = stt_all,
                               stt_two_states = stt_two_states)
    }else {
      island_list[[1]] <- list(island_age = total_time,
                               not_present = number_not_present,
                               stt_all = stt_all)
    }
    if (number_present > 0) {
      for (i in 1:number_present) {
        island_list[[1 + i]] <- full_list[[present[i]]]
        island_list[[1 + i]]$stt_table <- NULL
      }
    }
    if (number_present == 0) {
      island_list <- list()
      island_list[[1]] <- list(island_age = total_time,
                               not_present = M,
                               stt_all = stt_all)
    }
    several_islands[[rep]] <- island_list
    if (verbose == TRUE) {
      message(paste0("Island being formatted: ",
                     rep,
                     "/",
                     length(island_replicates)))
    }
  }
  return(several_islands)
}

#' Unsampled CS full STT
#'
#' @param stt_list List of full stt tables as
#' returned by DAISIE_sim_core functions
#' @param total_time Numeric double with total time of simulation.
#' @param stac_vec Vector with status of species on island.
#' @param trait_pars A named list containing diversification rates considering
#' two trait states created by \code{\link{create_trait_pars}}:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on
#'    mainland}
#' }
#'
#' @return 1 complete, unsampled STT table from all clades in an island of a
#' CS model as generated by DAISIE_sim_core functions.
#' @keywords internal
#' @author Pedro Neves, Joshua Lambert, Shu Xie, Giovanni Laudanno
create_full_CS_stt <- function(stt_list,
                               stac_vec,
                               total_time,
                               trait_pars = NULL) {

  if(!is.null(trait_pars)){
    return(
      create_full_CS_stt_trait(
        stt_list = stt_list,
        stac_vec = stac_vec,
        total_time = total_time,
        trait_pars = trait_pars
      )
    )
  }
  # Return empty island, if empty
  present <- which(stac_vec != 0)

  # Checks if stt has only 2 rows and is empty at present (nothing happened)
  second_line_stts <- lapply(stt_list, "[", 2,)
  zeros_second_line <- sapply(second_line_stts, sum) == 0


  filled_stt_lists <- stt_list[!zeros_second_line]

  # Calculate 'present' and append to filled_stt_list
  num_indep_colonists <- list()
  for (i in seq_along(filled_stt_lists)) {
    num_indep_colonists[[i]] <- filled_stt_lists[[i]][, 2] +
      filled_stt_lists[[i]][, 3] +
      filled_stt_lists[[i]][, 4]

    num_indep_colonists[[i]][which(num_indep_colonists[[i]] > 0)] <- 1
    filled_stt_lists[[i]] <- cbind(
      filled_stt_lists[[i]],
      present = num_indep_colonists[[i]]
    )
  }


  # If no colonization ever happened, just return 0 values
  if (length(filled_stt_lists) == 0) {
    times <- c(total_time, 0)
    nI <- c(0, 0)
    nA <- c(0, 0)
    nC <- c(0, 0)
    diff_present <- c(0, 0)
  } else {

    deltas_matrix <- lapply(filled_stt_lists, FUN = diff)
    for (i in seq_along(filled_stt_lists)) {
      if (any(filled_stt_lists[[i]][1, ] !=
              c("Time" = total_time, "nI" = 0, "nA" = 0, "nC" = 0, "present" = 0))) {
        deltas_matrix[[i]] <- rbind(
          filled_stt_lists[[i]][1, ],
          deltas_matrix[[i]]
        )
      } else {
        deltas_matrix[[i]] <- rbind(
          c("Time" = total_time, "nI" = 0, "nA" = 0, "nC" = 0, "present" = 0),
          deltas_matrix[[i]]
        )
      }
    }

    times_list <- lapply(filled_stt_lists, "[", , 1) # nolint
    times <- unlist(times_list)

    nI_list <- lapply(deltas_matrix, "[", , 2) # nolint
    nA_list <- lapply(deltas_matrix, "[", , 3) # nolint
    nC_list <- lapply(deltas_matrix, "[", , 4) # nolint
    present_list <- lapply(deltas_matrix, "[", , 5) # nolint

    nI <- unlist(nI_list)
    nA <- unlist(nA_list)
    nC <- unlist(nC_list)
    diff_present <- unlist(present_list)
  }

  full_stt <- data.frame(
    times = times,
    nI = nI,
    nA = nA,
    nC = nC,
    present = diff_present
  )
  ordered_diffs <- full_stt[order(full_stt$times, decreasing = TRUE), ]

  complete_stt_table <- mapply(ordered_diffs[2:5], FUN = cumsum)
  complete_stt_table <- cbind(ordered_diffs$times, complete_stt_table)
  colnames(complete_stt_table) <- c("Time", "nI", "nA", "nC", "present")

  while (complete_stt_table[1, 1] == complete_stt_table[2, 1]) {
    complete_stt_table <- complete_stt_table[-1, ]
  }

  stt <- complete_stt_table
  # Remove final duplicate lines, if any
  while (
    all(stt[nrow(stt) - 1, ] == stt[nrow(stt), ])
  ) {
    stt <- stt[1:(nrow(stt) - 1), ]
  }
  return(stt)
}

create_full_CS_stt_trait <- function(stt_list, stac_vec, total_time, trait_pars) {
  # Return empty island, if empty
  present <- which(stac_vec != 0)

  # Checks if stt has only 2 rows and is empty at present (nothing happened)
  second_line_stts <- lapply(stt_list, "[", 2,)
  zeros_second_line <- sapply(second_line_stts, sum) == 0
  filled_stt_lists <- stt_list[!zeros_second_line]

  # Calculate 'present' and append to filled_stt_list
  # no_time_stts <- lapply(filled_stt_lists, "[", , 2:4)
  num_indep_colonists <- list()
  for (i in seq_along(filled_stt_lists)) {
    num_indep_colonists[[i]] <- filled_stt_lists[[i]][, 2] +
      filled_stt_lists[[i]][, 3] +
      filled_stt_lists[[i]][, 4] +
      filled_stt_lists[[i]][, 5] +
      filled_stt_lists[[i]][, 6] +
      filled_stt_lists[[i]][, 7]

    num_indep_colonists[[i]][which(num_indep_colonists[[i]] > 0)] <- 1
    filled_stt_lists[[i]] <- cbind(
      filled_stt_lists[[i]],
      present = num_indep_colonists[[i]]
    )
  }
  # If no colonization ever happened, just return 0 values
  if (length(filled_stt_lists) == 0) {
    times <- c(total_time, 0)
    nI <- c(0, 0)
    nA <- c(0, 0)
    nC <- c(0, 0)
    nI2 <- c(0, 0)
    nA2 <- c(0, 0)
    nC2 <- c(0, 0)
    diff_present <- c(0, 0)
  } else {

    deltas_matrix <- lapply(filled_stt_lists, FUN = diff)
    for (i in seq_along(filled_stt_lists)) {
      if (any(filled_stt_lists[[i]][1, ] != c("Time" = total_time,
                                              "nI" = 0,
                                              "nA" = 0,
                                              "nC" = 0,
                                              "nI2" = 0,
                                              "nA2" = 0,
                                              "nC2" = 0,
                                              "present" = 0))) {

        deltas_matrix[[i]] <- rbind(
          filled_stt_lists[[i]][1, ],
          deltas_matrix[[i]]
        )
      } else {
        deltas_matrix[[i]] <- rbind(
          c("Time" = total_time,
            "nI" = 0,
            "nA" = 0,
            "nC" = 0,
            "nI2" = 0,
            "nA2" = 0,
            "nC2" = 0,
            "present" = 0),
          deltas_matrix[[i]]
        )
      }
    }

    times_list <- lapply(filled_stt_lists, "[", , 1) # nolint
    all_times <- unlist(times_list)
    times <- all_times

    nI_list <- lapply(deltas_matrix, "[", , 2) # nolint
    nA_list <- lapply(deltas_matrix, "[", , 3) # nolint
    nC_list <- lapply(deltas_matrix, "[", , 4) # nolint
    nI2_list <- lapply(deltas_matrix, "[", , 5) # nolint
    nA2_list <- lapply(deltas_matrix, "[", , 6) # nolint
    nC2_list <- lapply(deltas_matrix, "[", , 7) # nolint
    present_list <- lapply(deltas_matrix, "[", , 8) #nolint

    nI <- unlist(nI_list)
    nA <- unlist(nA_list)
    nC <- unlist(nC_list)
    nI2 <- unlist(nI2_list)
    nA2 <- unlist(nA2_list)
    nC2 <- unlist(nC2_list)
    diff_present <- unlist(present_list)
  }

  ###  create full_list separate traits(8 columns stt table)

  full_stt_complete <- data.frame(
    times = times,
    nI = nI,
    nA = nA,
    nC = nC,
    nI2 = nI2,
    nA2 = nA2,
    nC2 = nC2,
    present = diff_present
  )
  ordered_diffs_complete <- full_stt_complete[order(full_stt_complete$times, decreasing = TRUE), ]

  complete_stt_table_complete <- mapply(ordered_diffs_complete[2:8], FUN = cumsum)
  complete_stt_table_complete <- cbind(ordered_diffs_complete$times, complete_stt_table_complete)
  colnames(complete_stt_table_complete) <- c("Time", "nI", "nA", "nC", "nI2", "nA2", "nC2", "present")
  while (complete_stt_table_complete[1, 1] == complete_stt_table_complete[2, 1]) {
    complete_stt_table_complete <- complete_stt_table_complete[-1, ]
  }

  ###  create full_list combine nI, nA, nC from two states into one(5 columns stt table)
  full_stt_combined <- data.frame(
    times = times,
    nI = nI + nI2,
    nA = nA + nA2,
    nC = nC + nC2,
    present = diff_present
  )
  ordered_diffs_combined <- full_stt_combined[order(full_stt_combined$times, decreasing = TRUE), ]

  complete_stt_table_combined <- mapply(ordered_diffs_combined[2:5], FUN = cumsum)
  complete_stt_table_combined <- cbind(ordered_diffs_combined$times, complete_stt_table_combined)
  colnames(complete_stt_table_combined) <- c("Time", "nI", "nA", "nC", "present")

  while (complete_stt_table_combined[1, 1] == complete_stt_table_combined[2, 1]) {
    complete_stt_table_combined <- complete_stt_table_combined[-1, ]
  }

  stt_all <- complete_stt_table_combined
  stt_two_states <- complete_stt_table_complete
  # Remove final duplicate lines, if any
  while (
    all(stt_all[nrow(stt_all) - 1, ] == stt_all[nrow(stt_all), ])
  ) {
    stt_all <- stt_all[1:(nrow(stt_all) - 1), ]
  }
  while (
    all(stt_two_states[nrow(stt_two_states) - 1, ] == stt_two_states[nrow(stt_two_states), ])
  ) {
    stt_two_states <- stt_two_states[1:(nrow(stt_two_states) - 1), ]
  }
  stt <- list(stt_all = stt_all,
              stt_two_states = stt_two_states)
  return(stt)
}

#' Formats clade-specific simulation output into standard
#' DAISIE list output
#'
#' @inheritParams default_params_doc
#'
#' @return List with CS DAISIE simulation output
#' @keywords internal
DAISIE_format_CS_sampled_stt <- function(island_replicates,
                                         time,
                                         M,
                                         sample_freq,
                                         verbose = TRUE,
                                         trait_pars = NULL) {

  if (!is.null(trait_pars)) {
    return(
      DAISIE_format_CS_trait(
        island_replicates = island_replicates,
        time = time,
        M = M,
        sample_freq = sample_freq,
        verbose = verbose,
        trait_pars = trait_pars
      )
    )
  }
  total_time <- time
  several_islands <- list()
  for (rep in seq_along(island_replicates)) {
    full_list <- island_replicates[[rep]]
    stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
    number_not_present <- length(which(stac_vec == 0))
    present <- which(stac_vec != 0)
    number_present <- length(present)
    type_vec <- unlist(full_list)[which(names(unlist(full_list)) == "type1or2")]
    prop_type2_pool <- length(which(type_vec == 2)) / M
    number_type2_cols <- length(which(match(which(stac_vec != 0),
                                            which(type_vec == 2)) > 0))
    number_type1_cols <- number_present - number_type2_cols
    island_list <- list()
    for (i in 1:(number_present + 1)) {
      island_list[[i]] <- list()
    }
    ### all species
    stt_list <- list()
    for (i in 1:M) {
      stt_list[[i]] <- full_list[[i]]$stt_table
    }
    stt_all <- matrix(ncol = 5, nrow = sample_freq + 1)
    colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
    stt_all[, "Time"] <- rev(seq(from = 0,
                                 to = total_time,
                                 length.out = sample_freq + 1))

    immig_spec <- c()
    ana_spec <- c()
    for (i in 1:M) {
      immig_spec[i] <- sum(full_list[[i]]$stt_table[1, 2])
      ana_spec[i] <- sum(full_list[[i]]$stt_table[1, 3])
    }
    immig_spec <- sum(immig_spec)
    ana_spec <- sum(ana_spec)
    init_present <- immig_spec + ana_spec
    stt_all[1, 2:5] <- c(immig_spec, ana_spec, 0, init_present)

    for (i in 2:nrow(stt_all)) {
      the_age <- stt_all[i, "Time"]
      store_richness_time_slice <- matrix(nrow = M, ncol = 3)
      colnames(store_richness_time_slice) <- c("I", "A", "C")
      for (x in 1:M) {
        row_index <- max(which(stt_list[[x]][, "Time"] >= the_age))
        store_richness_time_slice[x, ] <- stt_list[[x]][row_index, 2:4]
      }
      count_time_slice <- store_richness_time_slice[, 1] +
        store_richness_time_slice[, 2] +
        store_richness_time_slice[, 3]
      present_time_slice <- rep(0, M)
      present_time_slice[which(count_time_slice > 0)] <- 1
      store_richness_time_slice <- cbind(store_richness_time_slice,
                                         present_time_slice)
      stt_all[i, c(2:5)] <- apply(store_richness_time_slice, 2, sum)
    }
    if (number_type2_cols > 0) {
      ######################################################### list type1
      stt_list_type1 <- list()
      for (i in 1:max(which(type_vec == 1))) {
        stt_list_type1[[i]] <- full_list[[i]]$stt_table
      }
      stt_type1 <- matrix(ncol = 5, nrow = sample_freq + 1)
      colnames(stt_type1) <- c("Time", "nI", "nA", "nC", "present")
      stt_type1[, "Time"] <- rev(seq(from = 0,
                                     to = total_time,
                                     length.out = sample_freq + 1))
      stt_type1[1, 2:5] <- c(0, 0, 0, 0)
      for (i in 2:nrow(stt_type1)) {
        the_age <- stt_type1[i, "Time"]
        store_richness_time_slice <- matrix(nrow = max(which(type_vec == 1)),
                                            ncol = 3)
        colnames(store_richness_time_slice) <- c("I", "A", "C")
        for (x in 1:max(which(type_vec == 1))) {
          store_richness_time_slice[x, ] <- stt_list_type1[[x]][max(
            which(stt_list_type1[[x]][, "Time"] >= the_age)), 2:4]
        }
        count_time_slice <- store_richness_time_slice[, 1] +
          store_richness_time_slice[, 2] +
          store_richness_time_slice[, 3]
        present_time_slice <- rep(0, max(which(type_vec == 1)))
        present_time_slice[which(count_time_slice > 0)] <- 1
        store_richness_time_slice <- cbind(store_richness_time_slice,
                                           present_time_slice)
        stt_type1[i, c(2:5)] <- apply(store_richness_time_slice, 2, sum)
      }
      ######################################################### list type2
      type2len <- length(which(type_vec == 2))
      stt_list_type2 <- list()
      for (i in 1:type2len) {
        stt_list_type2[[i]] <- full_list[[which(type_vec == 2)[i]]]$stt_table
      }
      stt_type2 <- matrix(ncol = 5, nrow = sample_freq + 1)
      colnames(stt_type2) <- c("Time", "nI", "nA", "nC", "present")
      stt_type2[, "Time"] <- rev(seq(from = 0,
                                     to = total_time,
                                     length.out = sample_freq + 1))
      stt_type2[1, 2:5] <- c(0, 0, 0, 0)
      for (i in 2:nrow(stt_type2)) {
        the_age <- stt_type2[i, "Time"]
        store_richness_time_slice <- matrix(nrow = type2len, ncol = 3)
        colnames(store_richness_time_slice) <- c("I", "A", "C")
        for (x in 1:type2len) {
          store_richness_time_slice[x, ] <- stt_list_type2[[x]][max(
            which(stt_list_type2[[x]][, "Time"] >= the_age)), 2:4]
        }
        count_time_slice <- store_richness_time_slice[, 1] +
          store_richness_time_slice[, 2] +
          store_richness_time_slice[, 3]
        present_time_slice <- rep(0, prop_type2_pool * M)
        present_time_slice[which(count_time_slice > 0)] <- 1
        store_richness_time_slice <- cbind(store_richness_time_slice,
                                           present_time_slice)
        stt_type2[i, c(2:5)] <- apply(store_richness_time_slice, 2, sum)
      }
      island_list[[1]] <- list(island_age = total_time,
                               not_present_type1 = DDD::roundn(
                                 M * (1 - prop_type2_pool)) -
                                 (number_type1_cols),
                               not_present_type2 = DDD::roundn(
                                 M * prop_type2_pool) - number_type2_cols,
                               stt_all = stt_all,
                               stt_type1 = stt_type1,
                               stt_type2 = stt_type2)
    } else {
      island_list[[1]] <- list(island_age = total_time,
                               not_present = number_not_present,
                               stt_all = stt_all)
    }
    if (number_present > 0) {
      for (i in 1:number_present) {
        island_list[[1 + i]] <- full_list[[present[i]]]
        island_list[[1 + i]]$stt_table <- NULL
      }
    }
    if (number_present == 0) {
      island_list <- list()
      island_list[[1]] <- list(island_age = total_time,
                               not_present = M,
                               stt_all = stt_all)
    }
    several_islands[[rep]] <- island_list
    if (verbose == TRUE) {
      message(
        "Island being formatted: ", rep, "/", length(island_replicates)
      )
    }
  }
  return(several_islands)
}


#' @title Plot island species-through-time (STT) plots
#' @description Produces STT plots. If only one type of species is present in
#' the simulated islands, STT is plotted for all species. If two types are
#' present, three plots are produced: STT for all, STT for type 1 and STT
#' for type 2.
#'
#' R plots with number of total, endemic and non-endemic STTs for different
#' types of species for the entire time span the islands were simulated.
#' 2.5-97.5th percentiles are plotted in light grey, 25-75th percentiles
#' plotted in dark grey.
#'
#' @inheritParams default_params_doc
#'
#' @return R plot.
#' @author Luis Valente
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @export
DAISIE_plot_sims <- function(
    island_replicates,
    plot_plus_one = TRUE,
    type = "all_species",
    sample_freq = 25,
    trait_pars = NULL
) {
  time <- max(island_replicates[[1]][[1]]$stt_all[, 1])
  if (sample_freq != Inf) {
    # Prepare dataset
    plot_lists <- DAISIE_convert_to_classic_plot(island_replicates,
                                                 trait_pars = trait_pars)
  } else {
    stop("Plotting STT with sample_freq = Inf not yet available. \n")
  }
  if (type == "all") {
    types <- names(plot_lists)
  } else {
    types <- type
  }
  num_plots <- sum(!sapply(plot_lists[types], FUN = is.null))
  graphics::par(mfrow = c(1, num_plots))
  for (type_here in types) {
    DAISIE_plot_stt(
      plot_plus_one = plot_plus_one,
      time = time,
      plot_lists = plot_lists,
      type = type_here)
  }
}


#' Translate user-friendly ontogeny codes to numerics
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric, 0 for null-ontogeny, 1 for beta function
#' @keywords internal
#' @examples translated_ontogeny <- DAISIE:::translate_island_ontogeny("const")
translate_island_ontogeny <- function(island_ontogeny) {

  if (island_ontogeny == "const" || island_ontogeny == 0) {
    island_ontogeny <- 0
  }
  if (island_ontogeny == "beta" || island_ontogeny == 1) {
    island_ontogeny <- 1
  }
  return(island_ontogeny)
}

#' Translate user-friendly sea-level codes to numerics
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric, 0 for null-sea-level, 1 for sine function
#' @keywords internal
#' @examples translated_sea_level <- DAISIE:::translate_sea_level("const")
translate_sea_level <- function(sea_level) {

  if (sea_level == "const" || sea_level == 0) {
    sea_level <- 0
  }

  if (sea_level == "sine" || sea_level == 1) {
    sea_level <- 1
  }
  return(sea_level)
}



#' Prepare input for DAISIE_stt
#'
#' @inheritParams default_params_doc
#'
#' @seealso \code{\link{DAISIE_plot_stt}}, \code{\link{DAISIE_plot_sims}}
#' @return a list with wrangled data to be used for plotting STT plots with
#' DAISIE_plot_stt
#' @keywords internal
DAISIE_convert_to_classic_plot <- function(simulation_outputs,
                                           trait_pars = NULL) {
  if (!is_simulation_outputs(simulation_outputs)) {
    stop(
      "'simulation_outputs' should be a set of simulation outputs. \n",
      "Actual value: ", simulation_outputs
    )
  }
  replicates <- length(simulation_outputs)
  ### STT ALL species
  s_freq <- length(simulation_outputs[[1]][[1]]$stt_all[, 1])
  complete_arr <- array(dim = c(s_freq, 6, replicates))
  for (x in 1:replicates) {
    if(is.null(trait_pars)){
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC"]
      total <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI"]
      complete_arr[, , x] <- cbind(simulation_outputs[[x]][[1]]$stt_all[, c("Time", "nI", "nA", "nC")],
                                   sum_endemics,
                                   total)
    }else{
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nA2"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC2"]
      total <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nA2"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC2"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI2"]
      nI <- simulation_outputs[[x]][[1]]$stt_all[, "nI"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI2"]
      nA <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nA2"]
      nC <- simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC2"]
      complete_arr[,,x]<-cbind(simulation_outputs[[x]][[1]]$stt_all[, 'Time'],
                               nI,
                               nA,
                               nC,
                               sum_endemics,
                               total)
    }
  }
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  testit::assert(stt_average_all == DAISIE_extract_stt_median(
    island_replicates = simulation_outputs,
    trait_pars = trait_pars
  ))
  stt_q0.025_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
  stt_q0.25_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
  stt_q0.75_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
  stt_q0.975_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
  colnames(stt_average_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.025_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.25_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.75_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.975_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  all_species <- list(
    stt_average = stt_average_all,
    stt_q0.025 = stt_q0.025_all,
    stt_q0.25 = stt_q0.25_all,
    stt_q0.75 = stt_q0.75_all,
    stt_q0.975 = stt_q0.975_all
  )
  if (is.null(simulation_outputs[[1]][[1]]$stt_type1) == FALSE) {
    ### STT TYPE1
    s_freq <- length(simulation_outputs[[1]][[1]]$stt_type1[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    for (x in 1:replicates) {
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_type1[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_type1[, "nC"]
      total <- simulation_outputs[[x]][[1]]$stt_type1[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_type1[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_type1[, "nI"]
      complete_arr[, , x] <- cbind(simulation_outputs[[x]][[1]]$stt_type1,
                                   sum_endemics,
                                   total)
    }
    stt_average_type1 <- apply(complete_arr, c(1, 2), stats::median)
    stt_q0.025_type1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
    stt_q0.25_type1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
    stt_q0.75_type1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
    stt_q0.975_type1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
    colnames(stt_average_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.025_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.25_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.75_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.975_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    type1_species <- list(
      stt_average = stt_average_type1,
      stt_q0.025 = stt_q0.025_type1,
      stt_q0.25 = stt_q0.25_type1,
      stt_q0.75 = stt_q0.75_type1,
      stt_q0.975 = stt_q0.975_type1
    )
    ### STT TYPE2
    s_freq <- length(simulation_outputs[[1]][[1]]$stt_type2[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    for (x in 1:replicates) {
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_type2[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_type2[, "nC"]
      total <- simulation_outputs[[x]][[1]]$stt_type2[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_type2[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_type2[, "nI"]
      complete_arr[, , x] <- cbind(
        simulation_outputs[[x]][[1]]$stt_type2,
        sum_endemics,
        total
      )
    }
    stt_average_type2 <- apply(complete_arr, c(1, 2), stats::median)
    stt_q0.025_type2 <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
    stt_q0.25_type2 <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
    stt_q0.75_type2 <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
    stt_q0.975_type2 <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
    colnames(stt_average_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.025_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.25_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.75_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.975_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    type2_species <- list(
      stt_average = stt_average_type2,
      stt_q0.025 = stt_q0.025_type2,
      stt_q0.25 = stt_q0.25_type2,
      stt_q0.75 = stt_q0.75_type2,
      stt_q0.975 = stt_q0.975_type2
    )
    return(list(
      all_species = all_species,
      type1_species = type1_species,
      type2_species = type2_species
    )
    )
  } else {
    return(list(
      all_species = all_species,
      type1_species = NULL,
      type2_species = NULL)
    )
  }
}

#' Create the Species-Through-Time plot. This is used to visualize
#' the output of DAISIE_sim functions
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
DAISIE_plot_stt <- function(
    plot_plus_one = TRUE,
    time,
    plot_lists = plot_lists,
    type = type
) {
  # Plot the y axis iff plus one
  y_axis_type <- "n"
  y_axis_label <- "No of species"
  if (plot_plus_one == TRUE) {
    y_axis_type <- "s"
    y_axis_label <- "No of species + 1"
  }
  stt <- plot_lists[[type]]
  if (is.null(stt)) {
    return()
  }
  suppressWarnings(
    graphics::plot(
      NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt$stt_q0.975)),
      ylab = y_axis_label,
      bty = "l", xaxs = "i", xlab = "Time before present",
      main = "Species-through-time - All species",
      log = "y", cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2,
      yaxt = y_axis_type
    )
  )
  graphics::polygon(c(stt$stt_average[, "Time"], rev(stt$stt_average[, "Time"])), c(stt$stt_q0.025[, "Total"] +
                                                                                      1, rev(stt$stt_q0.975[, "Total"] + 1)), col = "light grey", border = NA)
  graphics::polygon(c(stt$stt_average[, "Time"], rev(stt$stt_average[, "Time"])), c(stt$stt_q0.25[, "Total"] +
                                                                                      1, rev(stt$stt_q0.75[, "Total"] + 1)), col = "dark grey", border = NA)
  graphics::lines(stt$stt_average[, "Time"], stt$stt_average[, "Total"] + 1, lwd = 2)
  graphics::lines(stt$stt_average[, "Time"], stt$stt_average[, "nI"] + 1, lwd = 2, col = "cyan3")
  graphics::lines(stt$stt_average[, "Time"], stt$stt_average[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")
  legend_names <- c("Total", "Non-endemic", "Endemic")
  legend_colors <- c("black", "cyan3", "dodgerblue1")
  graphics::legend(
    time, max(stt$stt_q0.975), legend_names, lty = 1, lwd = 2,
    col = legend_colors, cex = 1.2, border = NA, bty = "n"
  )
  if (plot_plus_one == FALSE) {
    y_axis_values <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
    graphics::axis(2, at = y_axis_values, labels = y_axis_values - 1)
  }
}

#' Extract the STT median from the output of DAISIE_sim functions
#'
#' @inheritParams default_params_doc
#'
#' @return a matrix (?)
#' @keywords internal
DAISIE_extract_stt_median <- function(
    island_replicates,
    trait_pars = NULL) {

  replicates <- length(island_replicates)
  time <- max(island_replicates[[1]][[1]]$stt_all[, 1])
  ### STT ALL species
  s_freq <- length(island_replicates[[1]][[1]]$stt_all[, 1])
  complete_arr <- array(dim = c(s_freq, 6, replicates))
  for (x in 1:replicates) {

    if(is.null(trait_pars)){
      sum_endemics <- island_replicates[[x]][[1]]$stt_all[, "nA"] +
        island_replicates[[x]][[1]]$stt_all[, "nC"]
      total <- island_replicates[[x]][[1]]$stt_all[, "nA"] +
        island_replicates[[x]][[1]]$stt_all[, "nC"] +
        island_replicates[[x]][[1]]$stt_all[, "nI"]
      complete_arr[, , x] <- cbind(
        island_replicates[[x]][[1]]$stt_all[, c("Time", "nI", "nA", "nC")],
        sum_endemics,
        total)
    }else{
      sum_endemics <- island_replicates[[x]][[1]]$stt_all[, "nA"] +
        island_replicates[[x]][[1]]$stt_all[, "nC"] +
        island_replicates[[x]][[1]]$stt_all[, "nA2"] +
        island_replicates[[x]][[1]]$stt_all[, "nC2"]
      total <- island_replicates[[x]][[1]]$stt_all[, "nA"] +
        island_replicates[[x]][[1]]$stt_all[, "nC"] +
        island_replicates[[x]][[1]]$stt_all[, "nI"] +
        island_replicates[[x]][[1]]$stt_all[, "nA2"] +
        island_replicates[[x]][[1]]$stt_all[, "nC2"] +
        island_replicates[[x]][[1]]$stt_all[, "nI2"]
      nI <- island_replicates[[x]][[1]]$stt_all[, "nI"] +
        island_replicates[[x]][[1]]$stt_all[, "nI2"]
      nA <- island_replicates[[x]][[1]]$stt_all[, "nA"] +
        island_replicates[[x]][[1]]$stt_all[, "nA2"]
      nC <- island_replicates[[x]][[1]]$stt_all[, "nC"] +
        island_replicates[[x]][[1]]$stt_all[, "nC2"]
      complete_arr[,,x]<-cbind(island_replicates[[x]][[1]]$stt_all[, 'Time'],
                               nI,
                               nA,
                               nC,
                               sum_endemics,
                               total)
    }
  }
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  colnames(stt_average_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  return(stt_average_all)
}

DAISIE_format_CS_trait <- function(island_replicates,
                                   time,
                                   M,
                                   sample_freq,
                                   verbose = TRUE,
                                   trait_pars = NULL)
{
  total_time <- time
  several_islands <- list()

  for(rep in 1:length(island_replicates))
  {
    full_list <- island_replicates[[rep]]
    stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
    number_not_present <- length(which(stac_vec == 0))
    present <- which(stac_vec!=0)
    number_present <- length(present)
    type_vec <- unlist(full_list)[which(names(unlist(full_list)) == "type1or2")]
    prop_type2_pool <- length(which(type_vec == 2)) / M

    number_type2_cols <- length(which(match(which(stac_vec != 0),which(type_vec == 2)) > 0))
    number_type1_cols <- number_present-number_type2_cols

    island_list <- list()
    for(i in 1:(number_present + 1))
    {
      island_list[[i]] = list()
    }

    ### all species
    stt_list = list()
    for(i in 1:(M + trait_pars$M2))
    {
      stt_list[[i]] = full_list[[i]]$stt_table
    }
    stt_all = matrix(ncol = 8,nrow = sample_freq + 1)

    colnames(stt_all) = c("Time","nI","nA","nC","nI2","nA2","nC2","present")
    stt_all[,"Time"] = rev(seq(from = 0,to = total_time,length.out = sample_freq + 1))

    ####
    immig_spec <- c()
    ana_spec <- c()
    immig_spec2 <- c()
    ana_spec2 <- c()
    for (i in 1:(M + trait_pars$M2)) {
      immig_spec[i] <- sum(full_list[[i]]$stt_table[1, 2])
      ana_spec[i] <- sum(full_list[[i]]$stt_table[1, 3])
      immig_spec2[i] <- sum(full_list[[i]]$stt_table[1, 5])
      ana_spec2[i] <- sum(full_list[[i]]$stt_table[1, 6])
    }
    immig_spec <- sum(immig_spec)
    ana_spec <- sum(ana_spec)
    immig_spec2 <- sum(immig_spec2)
    ana_spec2 <- sum(ana_spec2)
    init_present <- immig_spec + ana_spec + immig_spec2 + ana_spec2
    stt_all[1, 2:8] <- c(immig_spec, ana_spec, 0, immig_spec2, ana_spec2, 0, init_present)

    ####
    for(i in 2:nrow(stt_all))
    {
      the_age = stt_all[i,"Time"]
      store_richness_time_slice = matrix(nrow = M + trait_pars$M2,ncol = 6)
      colnames(store_richness_time_slice) = c("I","A","C","I2","A2","C2")
      for(x in 1:(M + trait_pars$M2))
      {
        # testit::assert(x >= 1)
        # testit::assert(x <= length(stt_list))
        # testit::assert("Time" %in% colnames(stt_list[[x]]))
        store_richness_time_slice[x,] = stt_list[[x]][max(which(stt_list[[x]][,"Time"] >= the_age)),2:7]
      }
      count_time_slice = store_richness_time_slice[,1] +
        store_richness_time_slice[,2] +
        store_richness_time_slice[,3] +
        store_richness_time_slice[,4] +
        store_richness_time_slice[,5] +
        store_richness_time_slice[,6]
      present_time_slice = rep(0, M + trait_pars$M2)
      present_time_slice[which(count_time_slice>0)] = 1
      store_richness_time_slice = cbind(store_richness_time_slice,present_time_slice)
      stt_all[i,c(2:8)] = apply(store_richness_time_slice,2,sum)
    }

    island_list[[1]] = list(
      island_age = total_time,
      not_present = number_not_present,
      stt_all = stt_all
    )


    if(number_present > 0)
    {
      for(i in 1:number_present)
      {
        island_list[[1 + i]] = full_list[[present[i]]]
        island_list[[1 + i]]$stt_table = NULL
      }
    }

    if(number_present == 0)
    {
      island_list = list()
      island_list[[1]] = list(island_age = total_time,not_present = M, stt_all = stt_all)

    }

    several_islands[[rep]] = island_list
    if (verbose == TRUE) {
      message(
        "Island being formatted: ", rep, "/", length(island_replicates)
      )
    }
  }
  return(several_islands)
}


#' Measures if the input is a valid collection of simulation
#' outputs.
#'
#' @inheritParams default_params_doc
#'
#' @return TRUE if the input is a valid collection of simulation
#' outputs.
#' @author Richel J.C Bilderbeek, Pedro Neves
#' @keywords internal
is_simulation_outputs <- function(simulation_outputs) {
  for (n_replicate in seq_along(simulation_outputs)) {
    if (!"island_age" %in% names(simulation_outputs[[n_replicate]][[1]]))
      return(FALSE)
    if (!(names(simulation_outputs[[n_replicate]][[1]])[2] %in%
          c("not_present","not_present_type1"))) {
      return(FALSE)
    }
    if (!"stt_all" %in% names(simulation_outputs[[n_replicate]][[1]]))
      return(FALSE)
    # TODO: Figure out how to test this?
    # if (!"branching_times" %in% names(simulation_outputs)) return(FALSE)
    # if (!"stac" %in% names(simulation_outputs)) return(FALSE)
    # if (!"missing_species" %in% names(simulation_outputs)) return(FALSE)
  }
  if (is.list(simulation_outputs) && length(simulation_outputs) >= 1) {
    return(TRUE)
  }
}

#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
#' @keywords internal
DAISIE_sim_core_trait_dep <- function(
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

  if (is.null(trait_pars)){
    stop("A second set of rates should be contain considering two trait states.
         If only one state,run DAISIE_sim_cr instead.")
  }

  if (pars[4] == 0 && trait_pars$immig_rate2 == 0) {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }

  mainland_n2 <- trait_pars$M2
  mainland_ntotal <- mainland_n + mainland_n2
  if (mainland_n != 0){
    mainland_spec <- seq(1, mainland_n, 1)
  }else{
    mainland_spec <- c()
  }
  maxspecID <- mainland_ntotal

  island_spec <- c()
  stt_table <- matrix(ncol = 7)
  colnames(stt_table) <- c("Time","nI","nA","nC","nI2","nA2","nC2")
  stt_table[1,] <- c(total_time,0,0,0,0,0,0)
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]

  num_spec <- length(island_spec[, 1])
  num_immigrants <- length(which(island_spec[, 4] == "I"))

  #### Start Monte Carlo iterations ####
  while (timeval < total_time) {
    rates <- update_rates(
      timeval = timeval,
      total_time = total_time,
      gam = gam,
      laa = laa,
      lac = lac,
      mu = mu,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
      K = K,
      num_spec = num_spec,
      num_immigrants = num_immigrants,
      mainland_n = mainland_n,
      extcutoff = extcutoff,
      island_ontogeny = 0,
      sea_level = 0,
      island_spec = island_spec,
      trait_pars = trait_pars
    )
    timeval_and_dt <- calc_next_timeval(
      max_rates = rates,
      timeval = timeval
    )
    timeval <- timeval_and_dt$timeval

    if (timeval < total_time) {
      possible_event <- DAISIE_sample_event_trait_dep(
        rates = rates
      )

      updated_state <- DAISIE_sim_update_state_trait_dep(
        timeval = timeval,
        total_time = total_time,
        possible_event = possible_event,
        maxspecID = maxspecID,
        mainland_spec = mainland_spec,
        island_spec = island_spec,
        stt_table = stt_table,
        trait_pars = trait_pars
      )

      island_spec <- updated_state$island_spec
      maxspecID <- updated_state$maxspecID
      stt_table <- updated_state$stt_table
      num_spec <- length(island_spec[, 1])
      num_immigrants <- length(which(island_spec[, 4] == "I"))
    }
  }
  #### Finalize STT ####
  stt_table <- rbind(
    stt_table,
    c(
      0,
      stt_table[nrow(stt_table), 2],
      stt_table[nrow(stt_table), 3],
      stt_table[nrow(stt_table), 4],
      stt_table[nrow(stt_table), 5],
      stt_table[nrow(stt_table), 6],
      stt_table[nrow(stt_table), 7]
    )
  )

  return(island_spec)
}



#' create island spec helper function
#' @param time time
#' @param mainland_n mainland n
#' @param immig_rate igr
#' @param ext_rate er
#' @param ana_rate ar
#' @param clado_rate cr
#' @param immig_rate2 igr2
#' @param ext_rate2 er2
#' @param ana_rate2 ar2
#' @param clado_rate2 cr2
#' @param trans_rate tr
#' @param trans_rate2 tr2
#' @param M2 m2
#' @export
#' @return stochastic island spec table
create_island_spec <- function(time,
                               mainland_n,
                               K,
                               immig_rate,
                               ext_rate,
                               ana_rate,
                               clado_rate,
                               immig_rate2,
                               ext_rate2,
                               ana_rate2,
                               clado_rate2,
                               trans_rate,
                               trans_rate2,
                               M2) {

  local_trait_pars <- TRAISIERCPP::create_trait_pars(trans_rate = trans_rate,
                                                     immig_rate2 = immig_rate2, ext_rate2 = ext_rate2,
                                                     ana_rate2 = ana_rate2,
                                                     clado_rate2 = clado_rate2,
                                                     trans_rate2 = trans_rate2,
                                                     M2 = M2)


  answer <- DAISIE_sim_core_trait_dep(time = time,
                                            mainland_n = mainland_n,
                                            pars = c(clado_rate,
                                                     ext_rate,
                                                     K,
                                                     immig_rate,
                                                     ana_rate),
                                            island_ontogeny = 0,
                                            sea_level = 0,
                                            hyper_pars = TRAISIERCPP::create_hyper_pars(d = 0, x = 0),
                                            area_pars = TRAISIERCPP::create_area_pars(),
                                            extcutoff = 1000,
                                            trait_pars = local_trait_pars)

 return(answer)
}


#' Updates state of island given sampled event with two trait states.
#'
#' Makes the event happen by updating island species matrix and species IDs.
#' What event happens is determined by the sampling in the algorithm.
#'
#' @inheritParams default_params_doc
#'
#' @return The updated state of the system, which is a list with the
#' \code{island_spec} matrix, an integer \code{maxspecID} with the most recent
#' ID of species and the \code{stt_table}, a matrix with the current species
#' through time table.
#'
#' @keywords internal
#'
#' @seealso \link{DAISIE_sim_core_trait_dep}
DAISIE_sim_update_state_trait_dep <- function(timeval,
                                              total_time,
                                              possible_event,
                                              maxspecID,
                                              mainland_spec,
                                              island_spec,
                                              stt_table,
                                              trait_pars)
{
  if (possible_event > 10) {
    # Nothing happens
  }

  ##########################################
  #IMMIGRATION
  if (possible_event == 1)
  {
    colonist = DDD::sample2(mainland_spec,1)

    if (length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }

    if (length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA,1))
    }

    if (length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA,1)
    }
  }

  ##########################################
  #EXTINCTION
  if (possible_event == 2)
  {
    island_spec_state1 = which(island_spec[,8] == "1")
    extinct = DDD::sample2(island_spec_state1, 1)

    #this chooses the row of species data to remove

    typeofspecies = island_spec[extinct,4]

    if(typeofspecies == "I")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove immigrant

    if(typeofspecies == "A")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove anagenetic

    if(typeofspecies == "C")
    {
      #remove cladogenetic
      #first find species with same ancestor AND arrival total_time
      sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]),which(island_spec[,3] == island_spec[extinct,3]))
      survivors = sisters[which(sisters != extinct)]

      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic
        island_spec[survivors,4] = "A"
        island_spec[survivors,c(5,6)] = c(NA,NA)
        island_spec[survivors,7] = "Clado_extinct"
        island_spec = island_spec[-extinct,]
      }

      if(length(sisters) >= 3)
      {
        numberofsplits = nchar(island_spec[extinct, 5])

        mostrecentspl = substring(island_spec[extinct,5],numberofsplits)

        if(mostrecentspl=="B")
        {
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }

        motiftofind = paste(substring(island_spec[extinct,5],1,numberofsplits-1),sistermostrecentspl,sep = "")

        possiblesister = survivors[which(substring(island_spec[survivors,5],1,numberofsplits) == motiftofind)]

        #different rules depending on whether a B or A is removed. B going extinct is simpler because it only
        #carries a record of the most recent speciation
        if(mostrecentspl == "A")
        {
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec[possiblesister,6] == min(as.numeric(island_spec[possiblesister,6])))]
          island_spec[tochange,6] = island_spec[extinct,6]
        }

        #remove the offending A/B from these species
        island_spec[possiblesister,5] = paste(
          substring(island_spec[possiblesister,5],1,numberofsplits - 1),
          substring(island_spec[possiblesister,5],numberofsplits + 1, nchar(island_spec[possiblesister,5])),sep = "")
        island_spec = island_spec[-extinct,]
      }
    }
    island_spec = rbind(island_spec)
  }

  ##########################################
  #ANAGENESIS
  if(possible_event == 3)
  {
    immi_specs = intersect(which(island_spec[,4] == "I"), which(island_spec[,8] == "1"))
    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1)
    {
      anagenesis = immi_specs
    }
    if(length(immi_specs) > 1)
    {
      anagenesis = DDD::sample2(immi_specs,1)
    }

    maxspecID = maxspecID + 1
    island_spec[anagenesis,4] = "A"
    island_spec[anagenesis,1] = maxspecID
    island_spec[anagenesis,7] = "Immig_parent"
    if(!is.null(trait_pars)){
      island_spec[anagenesis,8] = "1"
    }
  }

  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive
  if(possible_event == 4)
  {
    island_spec_state1 = which(island_spec[,8] == "1")
    tosplit = DDD::sample2(island_spec_state1,1)

    #if the species that speciates is cladogenetic
    if(island_spec[tosplit,4] == "C")
    {
      #for daughter A

      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      oldstatus = island_spec[tosplit,5]
      island_spec[tosplit,5] = paste(oldstatus,"A",sep = "")
      #island_spec[tosplit,6] = timeval
      island_spec[tosplit,7] = NA
      island_spec[tosplit,8] = "1"

      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                        "C",paste(oldstatus,"B",sep = ""),timeval,NA,1))
      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic

      #for daughter A

      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,5] = "A"
      island_spec[tosplit,6] = island_spec[tosplit,3]
      island_spec[tosplit,7] = NA
      island_spec[tosplit,8] = "1"

      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA,1))
      maxspecID = maxspecID + 2
    }
  }

  ##########################
  ##transition from state1 to state2
  if(possible_event == 5){
    ##select a species with trait state1
    island_spec_state1 = which(island_spec[,8] == "1")
    totrans = DDD::sample2(island_spec_state1,1)
    island_spec[totrans,8] = "2"
  }

  ##########################
  ##immigration with state2
  if (possible_event == 6)
  {
    mainland1 = length(mainland_spec)
    mainland2 = trait_pars$M2
    mainland_total = mainland1 + mainland2
    colonist = DDD::sample2((mainland1 + 1):mainland_total,1)

    if (length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }

    if (length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA,2))
    }

    if (length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA,2)
    }
  }

  ##########################################
  #EXTINCTION
  if (possible_event == 7)
  {
    island_spec_state2 = which(island_spec[,8] == "2")
    extinct = DDD::sample2(island_spec_state2,1)
    #this chooses the row of species data with state2 to remove

    typeofspecies = island_spec[extinct,4]

    if(typeofspecies == "I")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove immigrant

    if(typeofspecies == "A")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove anagenetic

    if(typeofspecies == "C")
    {
      #remove cladogenetic
      #first find species with same ancestor AND arrival total_time
      sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]), which(island_spec[,3] == island_spec[extinct,3]))
      survivors = sisters[which(sisters != extinct)]

      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic
        island_spec[survivors,4] = "A"
        island_spec[survivors,c(5,6)] = c(NA,NA)
        island_spec[survivors,7] = "Clado_extinct"
        island_spec = island_spec[-extinct,]
      }

      if(length(sisters) >= 3)
      {
        numberofsplits = nchar(island_spec[extinct,5])

        mostrecentspl = substring(island_spec[extinct,5],numberofsplits)

        if(mostrecentspl=="B")
        {
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }

        motiftofind = paste(substring(island_spec[extinct,5],1,numberofsplits-1),sistermostrecentspl,sep = "")

        possiblesister = survivors[which(substring(island_spec[survivors,5],1,numberofsplits) == motiftofind)]

        #different rules depending on whether a B or A is removed. B going extinct is simpler because it only
        #carries a record of the most recent speciation
        if(mostrecentspl == "A")
        {
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec[possiblesister,6] == min(as.numeric(island_spec[possiblesister,6])))]
          island_spec[tochange,6] = island_spec[extinct,6]
        }

        #remove the offending A/B from these species
        island_spec[possiblesister,5] = paste(substring(island_spec[possiblesister,5],1,numberofsplits - 1),
                                              substring(island_spec[possiblesister,5],numberofsplits + 1,
                                                        nchar(island_spec[possiblesister,5])),sep = "")
        island_spec = island_spec[-extinct,]
      }
    }
    island_spec = rbind(island_spec)
  }

  ##########################################
  #ANAGENESIS
  if(possible_event == 8)
  {
    immi_specs = intersect(which(island_spec[,4] == "I"), which(island_spec[,8] == "2"))

    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1)
    {
      anagenesis = immi_specs
    }
    if(length(immi_specs) > 1)
    {
      anagenesis = DDD::sample2(immi_specs,1)
    }

    maxspecID = maxspecID + 1
    island_spec[anagenesis,4] = "A"
    island_spec[anagenesis,1] = maxspecID
    island_spec[anagenesis,7] = "Immig_parent"
    if(!is.null(trait_pars)){
      island_spec[anagenesis,8] = "2"
    }
  }

  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive
  if(possible_event == 9)
  {

    island_spec_state1 = which(island_spec[,8] == "2")
    tosplit = DDD::sample2(island_spec_state1,1)
    #if the species that speciates is cladogenetic
    if(island_spec[tosplit,4] == "C")
    {
      #for daughter A

      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      oldstatus = island_spec[tosplit,5]
      island_spec[tosplit,5] = paste(oldstatus,"A",sep = "")
      #island_spec[tosplit,6] = timeval
      island_spec[tosplit,7] = NA
      if(!is.null(trait_pars)){
        island_spec[tosplit,8] = "2"
      }
      #for daughter B
      if(is.null(trait_pars)){
        island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                          "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      }else{
        island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                          "C",paste(oldstatus,"B",sep = ""),timeval,NA,2))
      }


      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic

      #for daughter A

      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,5] = "A"
      island_spec[tosplit,6] = island_spec[tosplit,3]
      island_spec[tosplit,7] = NA
      if(!is.null(trait_pars)){
        island_spec[tosplit,8] = "2"
      }
      #for daughter B
      if(is.null(trait_pars)){
        island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA))
      }else{
        island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA,2))
      }
      maxspecID = maxspecID + 2
    }
  }


  ##########################
  ##transition from state2 to state1
  if(possible_event == 10){
    ##select a species with trait state1
    island_spec_state1 = which(island_spec[,8] == "2")
    totrans = DDD::sample2(island_spec_state1,1)
    island_spec[totrans,8] = "1"
  }



  if (possible_event <= 10 && total_time >= timeval) {
    stt_table <- rbind(stt_table,
                       c(total_time - timeval,
                         length(intersect(which(island_spec[,4] == "I"),which(island_spec[,8] == "1"))),    #nI1
                         length(intersect(which(island_spec[,4] == "A"),which(island_spec[,8] == "1"))),    #nA1
                         length(intersect(which(island_spec[,4] == "C"),which(island_spec[,8] == "1"))),    #nC1
                         length(intersect(which(island_spec[,4] == "I"),which(island_spec[,8] == "2"))),    #nI2
                         length(intersect(which(island_spec[,4] == "A"),which(island_spec[,8] == "2"))),    #nA2
                         length(intersect(which(island_spec[,4] == "C"),which(island_spec[,8] == "2")))))   #nC2
  }

  updated_state <- list(island_spec = island_spec,
                        maxspecID = maxspecID,
                        stt_table = stt_table)
  return(updated_state)
}

#' Calculates algorithm rates
#' @description Internal function that updates the all the rates and
#' max extinction horizon at time t.
#' @family rate calculations
#'
#' @inheritParams default_params_doc
#'
#' @seealso \code{\link{update_max_rates}()}
#' @keywords internal
#' @return a named list with the updated effective rates.
update_rates <- function(timeval,
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
  # Function to calculate rates at time = timeval. Returns list with each rate
  if (!is.null(trait_pars)) {
    return(
      update_rates_trait(
        timeval = timeval,
        total_time = total_time,
        gam = gam,
        mu = mu,
        laa = laa,
        lac = lac,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        extcutoff = extcutoff,
        K = K,
        mainland_n = mainland_n,
        num_spec = num_spec,
        num_immigrants = num_immigrants,
        trait_pars = trait_pars,
        island_spec = island_spec
      )
    )
  }

  A <- island_area(
    timeval = timeval,
    total_time = total_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  immig_rate <- get_immig_rate(
    gam = gam,
    A = A,
    num_spec = num_spec,
    K = K,
    mainland_n = mainland_n
  )
  # testit::assert(is.numeric(immig_rate))
  ext_rate <- get_ext_rate(
    mu = mu,
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    A = A
  )
  # testit::assert(is.numeric(ext_rate))
  ana_rate <- get_ana_rate(
    laa = laa,
    num_immigrants = num_immigrants
  )
  # testit::assert(is.numeric(ana_rate))
  clado_rate <- get_clado_rate(
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    K = K,
    A = A
  )
  # testit::assert(is.numeric(clado_rate))

  rates <- list(
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate
  )
  return(rates)
}

update_rates_trait <- function(timeval,
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

  A <- island_area(
    timeval = timeval,
    total_time = total_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  immig_rate <- get_immig_rate(
    gam = gam,
    A = A,
    num_spec = num_spec,
    K = K,
    mainland_n = mainland_n,
    trait_pars = trait_pars,
    island_spec = island_spec
  )

  ext_rate <- get_ext_rate(
    mu = mu,
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    A = A,
    trait_pars = trait_pars,
    island_spec = island_spec
  )

  ana_rate <- get_ana_rate(
    laa = laa,
    num_immigrants = num_immigrants,
    trait_pars = trait_pars,
    island_spec = island_spec
  )
  clado_rate <- get_clado_rate(
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    K = K,
    A = A,
    trait_pars = trait_pars,
    island_spec = island_spec
  )



  trans_rate <- get_trans_rate(trait_pars = trait_pars,
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
#' Function to describe changes in area through time.
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
#' @references
#' Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological
#' Sciences 281.1784 (2014): 20133227.
island_area <- function(timeval,
                        total_time,
                        area_pars,
                        peak,
                        island_ontogeny,
                        sea_level) {
  # testit::assert(are_area_pars(area_pars))
  Tmax <- area_pars$total_island_age
  Amax <- area_pars$max_area
  Acurr <- area_pars$current_area
  proptime_max <- area_pars$proportional_peak_t
  ampl <- area_pars$sea_level_amplitude
  freq <- area_pars$sea_level_frequency
  theta <- area_pars$island_gradient_angle
  proptime <- timeval / Tmax
  proptime_curr <- total_time / Tmax
  theta <- theta * (pi / 180)
  # Constant ontogeny and sea-level
  if (island_ontogeny == 0 & sea_level == 0) {
    if (Amax != 1 || is.null(Amax)) {
      warning("Constant island area requires a maximum area of 1.")
    }
    return(1)
  }

  # Beta function ontogeny and constant sea-level
  if (island_ontogeny == 1 & sea_level == 0) {
    At <- calc_Abeta(proptime = proptime,
                     proptime_max = proptime_max,
                     peak = peak,
                     Amax = Amax)
    return(At)
  }

  if (island_ontogeny == 0 & sea_level == 1) {
    angular_freq <- 2 * pi * freq
    delta_sl <- ampl * cos((proptime_curr - proptime) * angular_freq)
    r_curr <- sqrt((Acurr) / pi)
    h_curr <- tan(theta) * r_curr
    h_delta <- max(0, h_curr - ampl + delta_sl)
    At <- pi * (h_delta / tan(theta)) ^ 2
    return(At)
  }
  if (island_ontogeny == 1 && sea_level == 1) {
    A_beta <- calc_Abeta(proptime,
                         proptime_max,
                         peak,
                         Amax)
    angular_freq <- 2 * pi * freq
    delta_sl <- ampl * cos((proptime_curr - proptime) * angular_freq)
    r_curr <- sqrt(A_beta / pi)
    h_curr <- tan(theta) * r_curr
    h_delta <- max(0, h_curr - ampl + delta_sl)
    At <- pi * (h_delta / tan(theta)) ^ 2
    return(At)
  }
}

#' Function to describe changes in extinction rate through time.
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#' @family rate calculations
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784
#' (2014): 20133227.
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_ext_rate <- function(mu,
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
#' @keywords internal
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_ana_rate <- function(laa,
                         num_immigrants,
                         island_spec = NULL,
                         trait_pars = NULL) {

  if(is.null(trait_pars)){
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
#' @keywords internal
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_clado_rate <- function(lac,
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
#' @keywords internal
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_immig_rate <- function(gam,
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
#' @keywords internal
#' @family rates calculation
get_trans_rate <- function(trait_pars,
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
#' @keywords internal
#'
#' @author Joshua Lambert, Pedro Neves, Shu Xie
calc_next_timeval <- function(max_rates, timeval) {
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


#' Calculates when the next timestep will be, and if a shift has occured.
#'
#' @param timeval current time of simulation
#' @param max_rates named list of max rates as returned by
#' \code{\link{update_rates}}.
#' @param dynamic_shift_times numeric vector of times of rate shifts.
#'
#' @return named list with numeric vector containing the time of the next
#' timestep and the change in time.
#' @keywords internal
#'
#' @author Joshua Lambert, Pedro Neves, Shu Xie
calc_next_timeval_shift <- function(max_rates,
                                    timeval,
                                    dynamic_shift_times) {
  # testit::assert(timeval >= 0)
  totalrate <- max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]]
  dt <- stats::rexp(1, totalrate)
  timeval <- timeval + dt
  rate_shift <- FALSE

  if (timeval >= dynamic_shift_times[1]) {
    timeval <- dynamic_shift_times[1]
    dynamic_shift_times <- dynamic_shift_times[-1]
    rate_shift <- TRUE
  }

  out <- list(
    timeval = timeval,
    dt = dt,
    dynamic_shift_times = dynamic_shift_times,
    rate_shift = rate_shift
  )
  return(out)
}

#' Calculates the area at a point in time from a beta function
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
#' @author Joshua Lambert, Pedro Neves, Shu Xie
#'
#' @return Numeric
calc_Abeta <- function(proptime,
                       proptime_max,
                       peak,
                       Amax) {
  f <- proptime_max / (1 - proptime_max)
  a <- f * peak / (1 + f)
  b <- peak / (1 + f)
  At <- Amax * proptime ^ a *
    (1 - proptime) ^ b / ((a / (a + b)) ^ a * (b / (a + b)) ^ b)
  return(At)
}

#' Calculates the peak of ontogeny curve (beta function)
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
#' @return numeric
calc_peak <- function(total_time,
                      area_pars) {
  Amax <- area_pars$max_area
  Acurr <- area_pars$current_area
  proptime_max <- area_pars$proportional_peak_t
  proptime_curr <- total_time / area_pars$total_island_age
  # testit::assert(Acurr <= Amax)
  # testit::assert(proptime_max <= 1 & proptime_max >= 0)
  # testit::assert(proptime_curr <= 1 & proptime_curr >= 0)

  Abeta2 <- function(x) {
    calc_Abeta(proptime_curr, proptime_max, x, Amax) - Acurr
  }
  peak <- stats::uniroot(Abeta2, c(0.01, 1000))$root
  # testit::assert(is.numeric(peak))
  # testit::assert(is.finite(peak))
  return(peak)
}


#' Samples what event to happen next
#'
#' @inheritParams default_params_doc
#'
#' @return numeric indicating what event will happen, or a supposed event that
#' would happen in some timesteps of the ontogeny algorithm.
#' \itemize{
#'   \item{[1]: immigration event with trait1}
#'   \item{[2]: extinction event with trait1}
#'   \item{[3]: cladogenesis event with trait1}
#'   \item{[4]: anagenesis event with trait1}
#'   \item{[5]: transition event with trait1}
#'   \item{[6]: immigration event with trait2}
#'   \item{[7]: extinction event with trait2}
#'   \item{[8]: cladogenesis event with trait2}
#'   \item{[9]: anagenesis event with trait2}
#'   \item{[10]: transition event with trait2}
#' }
#' @author Shu Xie
#' @keywords internal
DAISIE_sample_event_trait_dep <- function(rates) {
  # testit::assert(are_rates(rates))
  possible_event <- sample(x = 1:10,
                           size = 1,
                           replace = FALSE,
                           prob = c(rates$immig_rate,
                                    rates$ext_rate,
                                    rates$ana_rate,
                                    rates$clado_rate,
                                    rates$trans_rate,
                                    rates$immig_rate2,
                                    rates$ext_rate2,
                                    rates$ana_rate2,
                                    rates$clado_rate2,
                                    rates$trans_rate2)
  )
  # testit::assert(is.numeric(possible_event))
  # testit::assert(possible_event >= 1)
  return(possible_event)
}

#' testing function
#' @param island_spec island_spec
#' @param extinct index
#' @return island spec
#' @export
DAISIE_test_execute_extinction <- function(island_spec,
                                           extinct) {


  typeofspecies = island_spec[extinct,4]

  if(typeofspecies == "I")
  {
    island_spec = island_spec[-extinct,]
  }
  #remove immigrant

  if(typeofspecies == "A")
  {
    island_spec = island_spec[-extinct,]
  }
  #remove anagenetic

  if(typeofspecies == "C")
  {
    #remove cladogenetic
    #first find species with same ancestor AND arrival total_time
    sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]),which(island_spec[,3] == island_spec[extinct,3]))
    survivors = sisters[which(sisters != extinct)]

    if(length(sisters) == 2)
    {
      #survivors status becomes anagenetic
      island_spec[survivors, 4] = "A"
      island_spec[survivors, c(5,6)] = c(NA,NA)
      island_spec[survivors,7] = "Clado_extinct"
      island_spec = island_spec[-extinct,]
    }

    if(length(sisters) >= 3)
    {
      numberofsplits = nchar(island_spec[extinct, 5])

      mostrecentspl = substring(island_spec[extinct,5],numberofsplits)

      if(mostrecentspl=="B")
      {
        sistermostrecentspl = "A"
      }
      if(mostrecentspl=="A")
      {
        sistermostrecentspl = "B"
      }

      motiftofind = paste(substring(island_spec[extinct,5],1,numberofsplits-1),sistermostrecentspl,sep = "")

      possiblesister = survivors[which(substring(island_spec[survivors,5],1,numberofsplits) == motiftofind)]

      #different rules depending on whether a B or A is removed. B going extinct is simpler because it only
      #carries a record of the most recent speciation
      if(mostrecentspl == "A") {
        #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
        tochange = possiblesister[which(island_spec[possiblesister,6] == min(as.numeric(island_spec[possiblesister,6])))]
        island_spec[tochange,6] = island_spec[extinct,6]
      }

      #remove the offending A/B from these species
      island_spec[possiblesister,5] = paste(
        substring(island_spec[possiblesister,5], 1, numberofsplits - 1),
        substring(island_spec[possiblesister,5], numberofsplits + 1, nchar(island_spec[possiblesister,5])),sep = "")
      island_spec = island_spec[-extinct,]
    }
  }
  island_spec = rbind(island_spec)
  return(island_spec)
}
