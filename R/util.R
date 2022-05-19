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
create_area_pars <- function(max_area,
                             current_area,
                             proportional_peak_t,
                             total_island_age,
                             sea_level_amplitude,
                             sea_level_frequency,
                             island_gradient_angle) {
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

#' Translate user-friendly ontogeny codes to numerics
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric, 0 for null-ontogeny, 1 for beta function
#' @keywords internal
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
translate_sea_level <- function(sea_level) {

  if (sea_level == "const" || sea_level == 0) {
    sea_level <- 0
  }

  if (sea_level == "sine" || sea_level == 1) {
    sea_level <- 1
  }
  return(sea_level)
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
#' @seealso \code{\link{DAISIE_sim_cr}},
#' \code{\link{DAISIE_sim_time_dep}},
#' \code{\link{DAISIE_sim_cr_shift}}, \code{\link{DAISIE_format_CS}}
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



