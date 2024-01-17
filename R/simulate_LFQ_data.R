#' simulate peptide level data with two groups
#' @export
#' @param Nprot number of porteins
#' @param N group size
#' @param fc D - down regulation N - matrix,  U -  regulation
#' @param prop proportion of down (D), up (U) and not regulated (N)
#' @examples
#' res <- sim_lfq_data()
#' sim_lfq_data(Nprot = 10)
#'
#' dd <- sim_lfq_data()
#' dd$abundance <- add_missing(dd$abundance)
#'
#' sim_lfq_data_config()
#' sim_lfq_data_config(with_missing = FALSE)

sim_lfq_data <- function(
    Nprot = 20,
    N = 4,
    fc = list(A = c(D = -2,  U = 2, N = 0), B = c(D = 1, U = -4)),
    prop = list(A = c(D = 10, U = 10), B = c(D = 5, U = 20)),
    mean_prot = 20,
    sd = 1.2,
    probability_of_success = 0.3
) {

  proteins <- stringi::stri_rand_strings(Nprot, 6)
  # simulate number of peptides per protein
  nrpeptides <- rgeom(n = Nprot, probability_of_success) + 1


  prot <- data.frame(
    proteinID = proteins,
    nrPeptides = nrpeptides,
    average_prot_abundance = rlnorm(Nprot,log(20),sdlog = log(sd)),
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

  # add row for each protein
  peptide_df <- prot |> tidyr::uncount( nrPeptides )
  # create peptide ids
  peptide_df$peptideID <- stringi::stri_rand_strings(sum(prot$nrPeptides), 8)

  # transform into long format
  peptide_df2 <- peptide_df |> tidyr::pivot_longer(cols = tidyselect::starts_with(c("mean", "N_")),
                                                   names_to = "group" , values_to = "mean")
  peptide_df2 <-  peptide_df2 |> tidyr::separate(group, c("what", "group"))
  peptide_df2 <- peptide_df2 |> tidyr::pivot_wider(names_from = "what", values_from = mean)

  peptide_df2$avg_peptide_abd <-
    with(peptide_df2,
         rlnorm(nrow(peptide_df2),
                meanlog = log(average_prot_abundance - mean),
                sdlog = log(sd)))

  sample_from_normal <- function(mean, sd, n) {
    rnorm(n, mean, sd)
  }

  nrpep <- nrow(peptide_df2)
  sampled_data <- matrix(nrow = nrpep, ncol = N)
  colnames(sampled_data) <- c("V1","V2","V3","V4")
  for (i in seq_len(nrpep)) {
    sampled_data[i,] <- sample_from_normal(peptide_df2$avg_peptide_abd[i], peptide_df2$sd[1], peptide_df2$N[i])
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
#'
#'
add_missing <- function(x){
  missing_prop <- pnorm(x, mean = mean(x), sd = sd(x))
  # sample TRUE or FALSE with propability in missing_prop
  samplemiss <- function(missing_prop) {
    mp <- c((1 - missing_prop)*0.2, missing_prop*3)
    mp <- mp / sum(mp)
    sample(c(TRUE, FALSE), size = 1, replace = TRUE, prob = mp)
  }

  missing_values <- sapply(missing_prop, samplemiss)
  # Introduce missing values into the vector x
  x[missing_values] <- NA
  return(x)
}


#' Simulate data with config
#' @export
#' @examples
#'
sim_lfq_data_config <- function(Nprot = 10, with_missing = TRUE){
  data <- sim_lfq_data(Nprot = Nprot)
  if (with_missing) {
    data$abundance <- add_missing(data$abundance)
  }
  data$isotopeLabel <- "light"
  data$qValue <- 0

  atable <- AnalysisTableAnnotation$new()
  atable$sampleName = "sample"
  atable$factors["group_"] = "group"
  atable$hierarchy[["protein_Id"]] = "proteinID"
  atable$hierarchy[["peptide_Id"]] = "peptideID"
  atable$set_response("abundance")

  config <- AnalysisConfiguration$new(atable)
  adata <- setup_analysis(data, config)
  return(list(data = adata, config = config))
}
