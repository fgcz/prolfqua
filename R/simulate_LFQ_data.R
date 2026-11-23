#' Generate random identifier strings of a fixed length
#' @param N number of identifiers to generate
#' @param str_length length of each generated identifier
#' @noRd
.generate_random_string <- function(N, str_length) {
  random_string <- vector(mode = "character", length = N)
  digits <- 0:9
  for (i in seq_len(N)) {
    random_string[i] <- paste(sample(digits, str_length, replace = TRUE), collapse = "")
  }
  random_string
}


#' Simulate protein-level data with two groups
#'
#' Draws protein abundances from a log-normal distribution and multiplies them by group-specific fold changes for
#' down-, up- and not-regulated proteins; with `peptide = TRUE`, each protein gets a geometric number of peptides.
#' The `_config` variants add an [AnalysisConfiguration] and missing values.
#' @return A data frame containing simulated LFQ data.
#' @export
#' @param n_proteins one positive whole number of proteins
#' @param n_per_group one positive whole number of samples per group
#' @param fc named list with one named numeric vector per group: fold changes of the down (D), up (U) and not (N)
#'   regulated proteins
#' @param prop named list with one named numeric vector per group: percentages of down (D), up (U) and not (N)
#'   regulated proteins
#' @param mean_prot one positive number, the mean protein abundance on the natural scale
#' @param sdlog one non-negative number, the standard deviation of the log-normal abundance distribution
#' @param probability_of_success one number in (0, 1], the probability of success of the geometric distribution of
#'   peptide counts
#' @param peptide one TRUE/FALSE: TRUE simulates peptide-level data, FALSE protein-level data
#' @param Nprot,N,PEPTIDE deprecated; use `n_proteins`, `n_per_group` and `peptide`
#' @examples
#'
#' res <- sim_lfq_data(n_proteins = 10)
#' res <- sim_lfq_data(n_proteins = 10, peptide = TRUE)
#' res <- sim_lfq_data(n_proteins = 10, n_per_group = 4)
#'
sim_lfq_data <- function(
  n_proteins = 20,
  n_per_group = 4,
  fc = list(A = c(D = -2, U = 2, N = 0), B = c(D = 1, U = -4)),
  prop = list(A = c(D = 10, U = 10), B = c(D = 5, U = 20)),
  mean_prot = 20,
  sdlog = log(1.2),
  probability_of_success = 0.3,
  peptide = FALSE,
  Nprot,
  N,
  PEPTIDE
) {
  if (!missing(Nprot)) {
    n_proteins <- .deprecate_arg("Nprot", "n_proteins", Nprot)
  }
  if (!missing(N)) {
    n_per_group <- .deprecate_arg("N", "n_per_group", N)
  }
  if (!missing(PEPTIDE)) {
    peptide <- .deprecate_arg("PEPTIDE", "peptide", PEPTIDE)
  }
  n_proteins <- .assert_count(n_proteins)
  n_per_group <- .assert_count(n_per_group)
  .assert_number(mean_prot, "be one positive number", function(v) v > 0)
  .assert_number(sdlog, "be one non-negative number", function(v) v >= 0)
  .assert_number(probability_of_success, "be one number in (0, 1]", function(v) v > 0 && v <= 1)
  .assert_flag(peptide)
  proteins <- stringi::stri_rand_strings(n_proteins, 6)
  idtype2 <- .generate_random_string(n_proteins, 4)
  # simulate number of peptides per protein
  nrpeptides <- rgeom(n = n_proteins, probability_of_success) + 1

  prot <- data.frame(
    proteinID = proteins,
    idtype2 = idtype2,
    nrPeptides = nrpeptides,
    average_prot_abundance = rlnorm(n_proteins, log(mean_prot), sdlog = sdlog),
    mean_Ctrl = 0,
    N_Ctrl = n_per_group,
    sd = 1
  )

  for (i in seq_along(fc)) {
    name <- names(fc)[i]
    fcx <- fc[[i]]
    propx <- prop[[i]]
    if (!"N" %in% names(fcx)) {
      fcx["N"] <- 0
    }
    if (!"N" %in% names(propx)) {
      propx["N"] <- 100 - sum(propx)
    }

    FC <- rep(fcx, ceiling(propx / 100 * n_proteins)) |> head(n = n_proteins)

    # add regulation to group A.
    group_mean <- paste0("mean_", name)
    group_size <- paste0("N_", name)
    prot <- prot |> dplyr::mutate(!!group_mean := FC, !!group_size := n_per_group)
  }

  if (peptide) {
    # add row for each protein
    peptide_df <- prot |> tidyr::uncount(nrPeptides)
    # create peptide ids
    peptide_df$peptideID <- stringi::stri_rand_strings(sum(prot$nrPeptides), 8)
  } else {
    peptide_df <- prot
  }

  # transform into long format
  peptide_df2 <- peptide_df |>
    tidyr::pivot_longer(cols = dplyr::starts_with(c("mean", "N_")), names_to = "group", values_to = "mean")
  peptide_df2 <- peptide_df2 |> tidyr::separate_wider_delim(group, delim = "_", names = c("what", "group"))
  peptide_df2 <- peptide_df2 |> tidyr::pivot_wider(names_from = "what", values_from = mean)

  nrpep <- nrow(peptide_df2)
  sampled_data <- matrix(nrow = nrpep, ncol = n_per_group)
  colnames(sampled_data) <- paste0("V", seq_len(ncol(sampled_data)))

  peptide_df2$average_prot_abundance <- peptide_df2$average_prot_abundance - peptide_df2$mean
  means <- peptide_df2$average_prot_abundance
  if (peptide) {
    peptide_df2$avg_peptide_abd <- rlnorm(nrpep, meanlog = log(means), sdlog = sdlog)
    means <- peptide_df2$avg_peptide_abd
  }
  for (i in seq_len(nrpep)) {
    sampled_data[i, ] <- rnorm(peptide_df2$N[i], means[i], peptide_df2$sd[i])
  }

  x <- dplyr::bind_cols(peptide_df2, sampled_data)

  peptide_abundances <- x |>
    tidyr::pivot_longer(
      dplyr::starts_with("V"),
      names_to = "Replicate",
      values_to = "abundance"
    )
  peptide_abundances <- peptide_abundances |>
    tidyr::unite("sample", group, Replicate, remove = FALSE)
  peptide_abundances
}


#' Draw missing values for an intensity vector based on its values
#'
#' Each value goes missing at random with a probability that falls with its intensity (through its normal
#' quantile), so low intensities go missing more often, as in LFQ data.
#' @return A logical vector indicating which values are missing.
#' @export
#' @param x numeric vector of intensities without NA
#' @param weight_missing one number in \[0, 1\]; greater weight means more missing values
#' @examples
#' which_missing(2**rnorm(10,2,0.4))
#'
which_missing <- function(x, weight_missing = 0.2) {
  if (!is.numeric(x) || anyNA(x)) {
    abort_bad_argument("x", "be a numeric vector without NA")
  }
  .assert_number(weight_missing, "be one number in [0, 1]", function(v) v >= 0 && v <= 1)
  missing_prop <- pnorm(x, mean = mean(x), sd = sd(x))
  # sample TRUE or FALSE with propability in missing_prop
  samplemiss <- function(missing_prop) {
    mp <- c((1 - missing_prop) * weight_missing, missing_prop * 3)
    mp <- mp / sum(mp)
    sample(c(TRUE, FALSE), size = 1, replace = TRUE, prob = mp)
  }

  missing_values <- vapply(missing_prop, samplemiss, logical(1))
  missing_values
}

#' Evaluates code with a fixed seed, or with the current RNG state when seed is NULL.
#' @param seed random seed or NULL
#' @param code expression evaluated with the selected random-number state
#' @noRd
.with_sim_seed <- function(seed, code) {
  if (is.null(seed)) code else withr::with_seed(seed, code)
}

#' Shared tail of the sim_lfq_data_*_config() simulators: configuration of the simulated columns and setup_analysis().
#' @param data long-format measurement data frame
#' @param factors named list of explanatory-variable columns
#' @param peptide whether the input contains peptide-level measurements
#' @param nr_children column tracking measured child-feature counts
#' @param factor_depth number of configured factors included in summaries and plots
#' @noRd
.sim_setup <- function(data, factors, peptide, nr_children = "nr_children", factor_depth = 1) {
  data$qValue <- 0
  config <- AnalysisConfiguration$new()
  config$file_name <- "sample"
  config$nr_children <- nr_children
  config$factors <- factors
  config$factor_depth <- factor_depth
  config$hierarchy[["protein_Id"]] <- c("proteinID", "idtype2")
  if (peptide) {
    config$hierarchy[["peptide_Id"]] <- "peptideID"
  }
  config$set_response("abundance")
  list(data = setup_analysis(data, config), config = config)
}


#' Simulate peptide-level data with its configuration
#'
#' Simulated peptide intensities of proteins in three groups (A, B, Ctrl), with missing values that are more likely
#' at low intensity (see [which_missing()]), returned with a matching [AnalysisConfiguration], ready for
#' `LFQData$new()`.
#' @param n_proteins one positive whole number of proteins to simulate
#' @param with_missing one TRUE/FALSE: add missing values (default TRUE)
#' @param weight_missing one number in \[0, 1\]; greater weight means more missingness
#' @param seed NULL or one whole number; NULL sets no seed
#' @param n_per_group one positive whole number of replicates per group
#' @param Nprot,N deprecated; use `n_proteins` and `n_per_group`
#' @return A list containing simulated data and an analysis configuration.
#' @export
#' @examples
#'
#' x <- sim_lfq_data_peptide_config()
#' stopifnot("data.frame" %in% class(x$data))
#' stopifnot("AnalysisConfiguration" %in% class(x$config))
sim_lfq_data_peptide_config <- function(
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  seed = 1234,
  n_per_group = 4,
  Nprot,
  N
) {
  if (!missing(Nprot)) {
    n_proteins <- .deprecate_arg("Nprot", "n_proteins", Nprot)
  }
  if (!missing(N)) {
    n_per_group <- .deprecate_arg("N", "n_per_group", N)
  }
  .assert_simulation_args(n_proteins, weight_missing)
  .assert_flag(with_missing)
  .assert_seed(seed)
  n_per_group <- .assert_count(n_per_group)
  .with_sim_seed(seed, {
    data <- sim_lfq_data(n_proteins = n_proteins, peptide = TRUE, n_per_group = n_per_group)

    data$nr_children <- 1
    if (with_missing) {
      data <- data[!which_missing(data$abundance, weight_missing = weight_missing), ]
    }
    .sim_setup(data, list(group_ = "group"), peptide = TRUE)
  })
}
#' Simulate protein-level data with its configuration
#'
#' Simulated protein intensities in three groups (A, B, Ctrl), with missing values that are more likely at low
#' intensity (see [which_missing()]), returned with a matching [AnalysisConfiguration], ready for `LFQData$new()`.
#' @param n_proteins one positive whole number of proteins to simulate
#' @param with_missing one TRUE/FALSE: add missing values (default TRUE)
#' @param weight_missing one number in \[0, 1\]; greater weight means more missingness
#' @param seed NULL or one whole number; NULL sets no seed
#' @param paired one TRUE/FALSE: add a paired subject factor to the design
#' @param Nprot deprecated; use `n_proteins`
#' @return A list containing simulated data and an analysis configuration.
#' @export
#' @examples
#'
#' x <- sim_lfq_data_protein_config()
#' stopifnot("data.frame" %in% class(x$data))
#' stopifnot("AnalysisConfiguration" %in% class(x$config))
#' x <- sim_lfq_data_protein_config(with_missing = FALSE)
#'
#' stopifnot(sum(is.na(x$data$abundance)) == 0)
#' xp <- sim_lfq_data_protein_config(with_missing = FALSE, paired = TRUE)
#' stopifnot(length(xp$config$factors) == 2)
#' stopifnot(nrow(xp$data) == nrow(x$data))
#' @srrstats {G2.0a} The help pages state the expected length of every vector input of the simulators,
#'   `LFQData$new()`, `null_split_check()`, the sample-size methods and the strategies ("one TRUE/FALSE", "one positive
#'   whole number", "numeric vector").
#' @srrstats {G2.4a} Counts (`n_proteins`, `n_per_group`, `n_splits`, `n_features`) must be whole numbers and are
#'   converted with `as.integer()` (`.assert_count()`, test-input-types.R).
#' @srrstats {G5.1} Documented simulators and bundled datasets are publicly available for tests and examples.
#' @srrstats {G5.5} Simulators use explicit reproducible seeds and restore the calling RNG state.
sim_lfq_data_protein_config <- function(
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  seed = 1234,
  paired = FALSE,
  Nprot
) {
  if (!missing(Nprot)) {
    n_proteins <- .deprecate_arg("Nprot", "n_proteins", Nprot)
  }
  .assert_simulation_args(n_proteins, weight_missing)
  .assert_flag(with_missing)
  .assert_seed(seed)
  .assert_flag(paired)
  .with_sim_seed(seed, {
    data <- sim_lfq_data(n_proteins = n_proteins, peptide = FALSE)
    if (paired) {
      annot <- data |> select(sample, group) |> distinct()
      annot <- annot |>
        mutate(subject = paste0("P", row_number()), .by = "group")
      data <- inner_join(annot, data, by = c("group", "sample"), relationship = "one-to-many")
    }
    if (with_missing) {
      data <- data[!which_missing(data$abundance, weight_missing = weight_missing), ]
    }
    factors <- if (paired) list(group_ = "group", subject_ = "subject") else list(group_ = "group")
    .sim_setup(data, factors, peptide = FALSE, nr_children = "nrPeptides")
  })
}


#' Simulate a two-factor design with its configuration
#'
#' Simulated protein (or peptide) intensities for a 2 x 2 design of `Treatment` and `Background` with four samples
#' per cell, returned with a matching [AnalysisConfiguration], ready for `LFQData$new()`.
#' @param n_proteins one positive whole number of proteins to simulate
#' @param with_missing one TRUE/FALSE: add missing values (default TRUE)
#' @param weight_missing one number in \[0, 1\]; greater weight means more missingness
#' @param peptide one TRUE/FALSE: TRUE simulates peptide-level data, FALSE protein-level data
#' @param seed NULL or one whole number; NULL sets no seed
#' @param two_factor one TRUE/FALSE: TRUE models `Treatment` and `Background` as two factors, FALSE as one combined
#'   factor
#' @param Nprot,PEPTIDE,TWO deprecated; use `n_proteins`, `peptide` and `two_factor`
#' @return A list containing simulated data and an analysis configuration.
#' @export
#' @examples
#' x <- sim_lfq_data_2factor_config(peptide = FALSE)
#' dim(x$data)
#' stopifnot("data.frame" %in% class(x$data))
#' stopifnot("AnalysisConfiguration" %in% class(x$config))
#' x <- sim_lfq_data_2factor_config(peptide = TRUE)
#'
#' head(x$data)
#' x <- sim_lfq_data_2factor_config(peptide = TRUE, two_factor = TRUE)
#' names(x)
#' nrow(x$data) > 10
#' x$data$Treatment |> table()
sim_lfq_data_2factor_config <- function(
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  peptide = FALSE,
  seed = 1234,
  two_factor = TRUE,
  Nprot,
  PEPTIDE,
  TWO
) {
  if (!missing(Nprot)) {
    n_proteins <- .deprecate_arg("Nprot", "n_proteins", Nprot)
  }
  if (!missing(PEPTIDE)) {
    peptide <- .deprecate_arg("PEPTIDE", "peptide", PEPTIDE)
  }
  if (!missing(TWO)) {
    two_factor <- .deprecate_arg("TWO", "two_factor", TWO)
  }
  .assert_simulation_args(n_proteins, weight_missing)
  .assert_flag(with_missing)
  .assert_flag(peptide)
  .assert_seed(seed)
  .assert_flag(two_factor)
  .with_sim_seed(seed, {
    res <- sim_lfq_data(
      n_proteins = n_proteins,
      peptide = peptide,
      fc = list(A = c(D = -2, U = 2, N = 0), B = c(D = 1, U = -4), C = c(D = -1, U = -4)),
      prop = list(A = c(D = 10, U = 10), B = c(D = 5, U = 20), C = c(D = 15, U = 25))
    )
    res <- res |> mutate(Treatment = case_when(group %in% c("Ctrl", "A") ~ "A", TRUE ~ "B"))
    data <- res |> mutate(Background = case_when(group %in% c("Ctrl", "C") ~ "Z", TRUE ~ "X"))

    if (!"nrPeptides" %in% colnames(data)) {
      data$nrPeptides <- 1
    }
    if (with_missing) {
      data <- data[!which_missing(data$abundance, weight_missing = weight_missing), ]
    }
    if (two_factor) {
      factors <- list(Treatment = "Treatment", Background = "Background")
    } else {
      data <- data |> tidyr::unite(Group, c("Treatment", "Background"))
      factors <- list(Group = "Group")
    }
    .sim_setup(data, factors, peptide = peptide, nr_children = "nrPeptides", factor_depth = if (two_factor) 2 else 1)
  })
}


#' Shared body of sim_build_models_lm() and sim_build_models_lmer(). mixed = TRUE simulates peptides, adds the random
#'
#' effects (1|peptide_Id) + (1|sampleName) and fits with strategy_lmer(); its "parallel3" data ignores n_proteins,
#' with_missing and weight_missing.
#' @param model model or model specification used by the operation
#' @param n_proteins number of simulated proteins
#' @param with_missing include simulated missing measurements
#' @param weight_missing weight controlling simulated missingness
#' @param mixed fit a nested mixed model rather than an aggregated linear model
#' @noRd
.sim_build_models <- function(model, n_proteins, with_missing, weight_missing, mixed) {
  if (model != "parallel3") {
    istar <- prolfqua::sim_lfq_data_2factor_config(
      n_proteins = n_proteins,
      with_missing = with_missing,
      peptide = mixed,
      weight_missing = weight_missing
    )
  } else if (mixed) {
    istar <- prolfqua::sim_lfq_data_peptide_config()
  } else {
    istar <- prolfqua::sim_lfq_data_protein_config(
      n_proteins = n_proteins,
      with_missing = with_missing,
      weight_missing = weight_missing
    )
  }
  istar <- prolfqua::LFQData$new(istar$data, istar$config)
  rhs <- c(
    factors = "~ Treatment + Background",
    interaction = "~ Treatment * Background",
    parallel2 = "~ Treatment",
    parallel3 = "~ group_"
  )[[model]]
  if (mixed) {
    rhs <- paste0(rhs, " + (1|peptide_Id) + (1|sampleName)")
  }
  strategy <- if (mixed) strategy_lmer else strategy_lm
  build_model(istar, strategy(paste0(istar$response(), rhs)))
}

#' Build linear models on simulated data for testing
#' @param model design: `"parallel2"`, `"parallel3"`, `"factors"` or `"interaction"`
#' @param n_proteins number of simulated proteins
#' @param with_missing include simulated missing measurements
#' @param weight_missing weight controlling simulated missingness
#' @param Nprot deprecated; use `n_proteins`
#' @return a [Model] object with one linear model per simulated protein, fitted with [build_model()]
#' @noRd
.sim_build_models_lm <- function(
  model = c("parallel2", "parallel3", "factors", "interaction"),
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 1,
  Nprot
) {
  if (!missing(Nprot)) {
    n_proteins <- .deprecate_arg("Nprot", "n_proteins", Nprot)
  }
  .sim_build_models(match.arg(model), n_proteins, with_missing, weight_missing, mixed = FALSE)
}

#' Build mixed-effects models on simulated data for testing
#' @return a [Model] object with one mixed-effects model per simulated protein, fitted with [build_model()]
#' @noRd
.sim_build_models_lmer <- function(
  model = c("parallel2", "parallel3", "factors", "interaction"),
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 1,
  Nprot
) {
  if (!missing(Nprot)) {
    n_proteins <- .deprecate_arg("Nprot", "n_proteins", Nprot)
  }
  .sim_build_models(match.arg(model), n_proteins, with_missing, weight_missing, mixed = TRUE)
}


#' Fit a linear model to simulated data for examples
#' @param model design: `"parallel2"`, `"parallel3"`, `"factors"` or `"interaction"`
#' @return an `lm` object fitted to one simulated protein
#' @noRd
.sim_make_model_lm <- function(model = c("parallel2", "parallel3", "factors", "interaction")) {
  .sim_build_models_lm(model = match.arg(model), n_proteins = 1, with_missing = FALSE)$model_df$linear_model[[1]]
}


#' Fit a mixed-effects model to simulated data for examples
#' @param model design: `"parallel2"`, `"parallel3"`, `"factors"` or `"interaction"`
#' @return an `lmerModLmerTest` object, the first successful fit among ten simulated proteins
#' @noRd
.sim_make_model_lmer <- function(model = c("parallel2", "parallel3", "factors", "interaction")) {
  model_df <- .sim_build_models_lmer(model = match.arg(model), n_proteins = 10, with_missing = FALSE)$model_df
  # first protein with a successful fit (isSingular is NA where the fit failed)
  model_df$linear_model[!is.na(model_df$isSingular)][[1]]
}
