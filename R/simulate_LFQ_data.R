.generate_random_string <- function(N, str_length) {
  random_string <- vector(mode = "character", length = N)
  digits <- 0:9
  for (i in seq_len(N)) {
    random_string[i] <- paste0(sample(digits, str_length, replace = TRUE), collapse = "")
  }
  return(random_string)
}


#' simulate protein level data with two groups
#' @export
#' @param Nprot number of porteins
#' @param N group size
#' @param fc D - down regulation N - matrix,  U -  regulation
#' @param prop proportion of down (D), up (U) and not regulated (N)
#' @examples
#'
#' res <- sim_lfq_data(Nprot = 10)
#' res <- sim_lfq_data(Nprot = 10, PEPTIDE = TRUE)
#'

sim_lfq_data <- function(
    Nprot = 20,
    N = 4,
    fc = list(A = c(D = -2,  U = 2, N = 0), B = c(D = 1, U = -4)),
    prop = list(A = c(D = 10, U = 10), B = c(D = 5, U = 20)),
    mean_prot = 20,
    sdlog = log(1.2),
    probability_of_success = 0.3,
    PEPTIDE = FALSE
) {


  proteins <- stringi::stri_rand_strings(Nprot, 6)
  idtype2 <- .generate_random_string(Nprot, 4)
  # simulate number of peptides per protein
  nrpeptides <- rgeom(n = Nprot, probability_of_success) + 1


  prot <- data.frame(
    proteinID = proteins,
    idtype2 = idtype2,
    nr_peptides = nrpeptides,
    average_prot_abundance = rlnorm(Nprot,log(20),sdlog = sdlog),
    mean_Ctrl = 0,
    N_Ctrl = N,
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

    FC = rep(fcx, ceiling(propx / 100 * Nprot)) |> head(n = Nprot)

    # add regulation to group A.
    groupMean <- paste0("mean_", name)
    groupSize <- paste0("N_", name)
    prot <- prot |> dplyr::mutate(!!groupMean := FC, !!groupSize := N)
  }

  if (PEPTIDE) {

    # add row for each protein
    peptide_df <- prot |> tidyr::uncount( nr_peptides )
    # create peptide ids
    peptide_df$peptideID <- stringi::stri_rand_strings(sum(prot$nr_peptides), 8)
  } else {
    peptide_df <- prot
  }


  # transform into long format
  peptide_df2 <- peptide_df |> tidyr::pivot_longer(cols = tidyselect::starts_with(c("mean", "N_")),
                                                   names_to = "group" , values_to = "mean")
  peptide_df2 <-  peptide_df2 |> tidyr::separate(group, c("what", "group"))
  peptide_df2 <- peptide_df2 |> tidyr::pivot_wider(names_from = "what", values_from = mean)

  sample_from_normal <- function(mean, sd, n) {
    rnorm(n, mean, sd)
  }
  nrpep <- nrow(peptide_df2)
  sampled_data <- matrix(nrow = nrpep, ncol = N)
  colnames(sampled_data) <- paste0("V", 1:ncol(sampled_data))

  peptide_df2$average_prot_abundance <- peptide_df2$average_prot_abundance - peptide_df2$mean

  if (PEPTIDE) {
    peptide_df2$avg_peptide_abd <-
      with(peptide_df2,
           rlnorm(nrow(peptide_df2),
                  meanlog = log(average_prot_abundance),
                  sdlog = sdlog))
    for (i in seq_len(nrpep)) {
      sampled_data[i,] <- sample_from_normal(peptide_df2$avg_peptide_abd[i], peptide_df2$sd[1], peptide_df2$N[i])
    }

  } else {
    for (i in seq_len(nrpep)) {
      sampled_data[i,] <- sample_from_normal(peptide_df2$average_prot_abundance[i], peptide_df2$sd[1], peptide_df2$N[i])
    }
  }

  x <- dplyr::bind_cols(peptide_df2,sampled_data)

  peptideAbudances <- x |>
    tidyr::pivot_longer(
      tidyselect::starts_with("V"),
      names_to = "Replicate",
      values_to = "abundance")
  peptideAbundances <- peptideAbudances |>
    tidyr::unite("sample", group, Replicate, remove =  FALSE)
  return(peptideAbundances)
}


#' add missing values to x vector based on the values of x
#' @export
#' @param x vector of intensities
#' @param weight_missing greater weight more missing
#' @examples
#' which_missing(2**rnorm(10,2,0.4))
#'
which_missing <- function(x, weight_missing = 0.2){
  missing_prop <- pnorm(x, mean = mean(x), sd = sd(x))
  # sample TRUE or FALSE with propability in missing_prop
  samplemiss <- function(missing_prop) {
    mp <- c((1 - missing_prop)*weight_missing, missing_prop*3)
    mp <- mp / sum(mp)
    sample(c(TRUE, FALSE), size = 1, replace = TRUE, prob = mp)
  }

  missing_values <- sapply(missing_prop, samplemiss)
  # Introduce missing values into the vector x
  #x[missing_values] <- NA
  return(missing_values)
}


#' Simulate data, protein and peptide, with config
#' @param description Nprot number of proteins
#' @param with_missing add missing values, default TRUE
#' @param seed seed for reproducibility, if NULL no seed is set.
#' @export
#' @examples
#'
#' x <- sim_lfq_data_peptide_config()
#' stopifnot("data.frame" %in% class(x$data))
#' stopifnot("AnalysisConfiguration" %in% class(x$config))
sim_lfq_data_peptide_config <- function(
    Nprot = 10,
    with_missing = TRUE,
    weight_missing = 0.2,
    seed = 1234){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  data <- sim_lfq_data(Nprot = Nprot, PEPTIDE = TRUE)

  not_missing <- !which_missing(data$abundance, weight_missing = weight_missing)
  # data <- data[not_missing,]
  data$nr_children <- as.numeric(not_missing)
  if (with_missing) {
    data <- data[data$nr_children > 0,]
  }
  data$isotopeLabel <- "light"
  data$qValue <- 0

  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "sample"

  atable$factors["group_"] = "group"
  atable$hierarchy[["protein_Id"]] = c("proteinID", "idtype2")
  atable$hierarchy[["peptide_Id"]] = "peptideID"
  atable$set_response("abundance")

  config <- AnalysisConfiguration$new(atable)
  adata <- setup_analysis(data, config)
  return(list(data = adata, config = config))
}
#' Simulate data, protein, with config
#' @param description Nprot number of proteins
#' @param with_missing add missing values, default TRUE
#' @param seed seed for reproducibility, if NULL no seed is set.
#' @export
#' @examples
#'
#' x <- sim_lfq_data_protein_config()
#' stopifnot("data.frame" %in% class(x$data))
#' stopifnot("AnalysisConfiguration" %in% class(x$config))
#' x <- sim_lfq_data_protein_config(with_missing = FALSE)
#' stopifnot(sum(is.na(x$data$abundance)) == 0)
sim_lfq_data_protein_config <- function(Nprot = 10, with_missing = TRUE, weight_missing = 0.2, seed = 1234){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  data <- sim_lfq_data(Nprot = Nprot, PEPTIDE = FALSE)

  data$nr_peptides[which_missing(data$abundance,weight_missing = weight_missing)] <- 0
  if (with_missing) {
    data <- data[data$nr_peptides > 0,]
  }

  data$isotopeLabel <- "light"
  data$qValue <- 0

  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "sample"
  atable$nr_children = "nr_peptides"
  atable$factors["group_"] = "group"
  atable$hierarchy[["protein_Id"]] = c("proteinID", "idtype2")
  atable$set_response("abundance")

  config <- AnalysisConfiguration$new(atable)
  adata <- setup_analysis(data, config)
  return(list(data = adata, config = config))
}


#' Simulate data, protein, with config with 2 factros Treatment and Background
#' @param description Nprot number of proteins
#' @param with_missing add missing values, default TRUE
#' @param seed seed for reproducibility, if NULL no seed is set.
#' @param TWO use two factors for modellin
#' @export
#' @examples
#' x <- sim_lfq_data_2Factor_config(PEPTIDE= FALSE)
#' dim(x$data)
#' stopifnot("data.frame" %in% class(x$data))
#' stopifnot("AnalysisConfiguration" %in% class(x$config))
#' x <- sim_lfq_data_2Factor_config(PEPTIDE = TRUE)
#'
#' head(x$data)
#' x <- sim_lfq_data_2Factor_config(PEPTIDE = TRUE, TWO = FALSE)
#' x$data$Group |> table()
sim_lfq_data_2Factor_config <- function(Nprot = 10,
                                        with_missing = TRUE,
                                        weight_missing = 0.2,
                                        PEPTIDE = FALSE,
                                        seed = 1234,
                                        TWO = TRUE
){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  res <- sim_lfq_data(Nprot = Nprot, PEPTIDE = PEPTIDE,
                      fc = list(A = c(D = -2,  U = 2, N = 0), B = c(D = 1, U = -4), C = c(D = -1, U = -4)),
                      prop = list(A = c(D = 10, U = 10), B = c(D = 5, U = 20), C = c(D = 15, U = 25)))
  res <- res |> mutate(Treatment = case_when(group %in% c("Ctrl", "A") ~ "A", TRUE ~ "B"))
  data <- res |> mutate(Background = case_when(group %in% c("Ctrl", "C") ~ "Z", TRUE ~ "X"))

  if (is.null(data$nr_peptides)) {
    data$nr_peptides <- 1
  }
  data$nr_peptides[which_missing(data$abundance,weight_missing = weight_missing)] <- 0
  if (with_missing) {
    data <- data[data$nr_peptides > 0,]
  }

  data$isotopeLabel <- "light"
  data$qValue <- 0
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "sample"
  atable$nr_children = "nr_peptides"

  if (TWO) {
    atable$factors["Treatment"] = "Treatment"
    atable$factors["Background"] = "Background"
    atable$factorDepth <- 2
  } else {
    data <- data |> tidyr::unite(Group, c("Treatment", "Background"))
    atable$factors["Group"] = "Group"
    atable$factorDepth <- 1
  }
  atable$hierarchy[["protein_Id"]] = c("proteinID", "idtype2")
  if (PEPTIDE) {
    atable$hierarchy[["peptide_Id"]] = c("peptideID")
  }
  atable$set_response("abundance")

  config <- AnalysisConfiguration$new(atable)
  adata <- setup_analysis(data, config)
  return(list(data = adata, config = config))
}

#' build dataframe with models for testing
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#' modi <- sim_build_models_lm(model = "interaction", weight_missing = 1)
#' stopifnot(dim(modi$modelDF) == c(10,9))
#' mod2 <- sim_build_models_lm(model = "parallel2", weight_missing = 1)
#' mod2$modelDF$linear_model[[1]]
#' mod3 <- sim_build_models_lm(model = "parallel3", weight_missing = 1)
#' modf <- sim_build_models_lm(model = "factors", weight_missing = 1)
#'
sim_build_models_lm <- function(model = c("parallel2","parallel3","factors", "interaction"),
                                Nprot = 10,
                                with_missing = TRUE,
                                weight_missing = 1) {
  model <- match.arg(model)
  if (model != "parallel3") {
    istar <- prolfqua::sim_lfq_data_2Factor_config(
      Nprot = Nprot,
      with_missing = with_missing,
      weight_missing = weight_missing)
  } else {
    istar <- prolfqua::sim_lfq_data_protein_config()
  }
  istar <- prolfqua::LFQData$new(istar$data,istar$config)

  model <- if (model == "factors") {
    "~ Treatment + Background"
  } else if (model == "interaction") {
    "~ Treatment * Background"
  } else if (model == "parallel2") {
    "~ Treatment"
  } else if (model == "parallel3") {
    "~ group_"
  } else {NULL}
  modelFunction <- strategy_lm(paste0(istar$response(), model))
  mod <- build_model(
    istar,
    modelFunction)
  return(mod)
}

#' build lmer model from simulated data
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#' undebug(sim_build_models_lmer)
#' modi <- sim_build_models_lmer(model = "interaction", weight_missing = 1)
#' stopifnot(sum(modi$modelDF$exists_lmer) == 6)
#' mod2 <- sim_build_models_lmer(model = "parallel2", weight_missing = 1)
#' stopifnot(sum(mod2$modelDF$exists_lmer) == 6)
#' mod4 <- sim_build_models_lmer(model = "parallel3", weight_missing = 1)
#' stopifnot(sum(mod4$modelDF$exists_lmer) == 6)
#' modf <- sim_build_models_lmer(model = "factors", weight_missing = 1)
#' stopifnot(sum(modf$modelDF$exists_lmer) == 6)
#'
sim_build_models_lmer <- function(model = c("parallel2", "parallel3","factors", "interaction"),
                                  Nprot = 10,
                                  with_missing = TRUE,
                                  weight_missing = 1) {
  model <- match.arg(model)
  if (model != "parallel3") {
    istar <- prolfqua::sim_lfq_data_2Factor_config(
      Nprot = Nprot,
      with_missing = with_missing,
      PEPTIDE = TRUE,
      weight_missing = weight_missing)
  } else {
    istar <- prolfqua::sim_lfq_data_peptide_config()
  }
  istar <- prolfqua::LFQData$new(istar$data,istar$config)

  model <- if (model == "factors") {
    "~ Treatment + Background + (1|peptide_Id) + (1|sampleName)"
  } else if (model == "interaction") {
    "~ Treatment * Background + (1|peptide_Id) + (1|sampleName)"
  } else if (model == "parallel2") {
    "~ Treatment + (1|peptide_Id) + (1|sampleName)"
  } else if (model == "parallel3") {
    "~ group_ + (1|peptide_Id) + (1|sampleName)"
  } else {NULL}
  modelFunction <- strategy_lmer(paste0(istar$response(), model))
  mod <- build_model(
    istar,
    modelFunction)
  return(mod)
}


#' make interaction model for examples
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#' m <- sim_make_model_lm()
#' mi <- sim_make_model_lm("interaction")
#' stopifnot(length(coefficients(summary(mi))[,"Estimate"]) == 4)
#' mf <- sim_make_model_lmer("factors")
#' m2 <- sim_make_model_lmer("parallel2")
#' m3 <- sim_make_model_lmer("parallel3")
sim_make_model_lm <- function(model = c("parallel2", "parallel3","factors", "interaction")){
  model <- match.arg(model)
  mod <- sim_build_models_lm(model = model, Nprot = 1, with_missing = FALSE)
  return(mod$modelDF$linear_model[[1]])
}


#' make interaction model for examples
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#' mf <- sim_make_model_lmer("factors")
#' mi <- sim_make_model_lmer("interaction")
#'
sim_make_model_lmer <- function(model = c("parallel2", "parallel3","factors", "interaction"),
                                singular = FALSE){
  model <- match.arg(model)
  mod <- sim_build_models_lmer(model = model, Nprot = 10, with_missing = FALSE)
  m <- mod$modelDF |> dplyr::filter(isSingular == isSingular) |> dplyr::pull(linear_model)
  return(m[[1]])
}
