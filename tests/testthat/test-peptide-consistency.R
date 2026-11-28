# Peptide consistency of aggregated proteins and its use in the facades' default weights

# Reference implementation, the prolfquabenchmark helper (workflow/scripts/peptide_consistency.R) with its column
# names: median polish of each protein's observed peptides, residual variances per protein and sample, two-level
# limma::squeezeVar(). Returns one row per protein and sample with an observed peptide.
reference_consistency <- function(peptides) {
  response <- peptides$response()
  sname <- peptides$sample_name()
  data <- as.data.frame(peptides$data_long())[, c("protein_Id", "peptide_Id", sname, response)]
  data <- data[!is.na(data[[response]]), ]
  observations <- dplyr::bind_rows(lapply(split(data, data$protein_Id), function(d) {
    features <- unique(d$peptide_Id)
    samples <- unique(d[[sname]])
    i <- match(d$peptide_Id, features)
    j <- match(d[[sname]], samples)
    mat <- matrix(NA_real_, length(features), length(samples))
    mat[cbind(i, j)] <- d[[response]]
    X <- suppressWarnings(stats::medpolish(mat, na.rm = TRUE, trace.iter = FALSE, maxiter = 10, eps = 0.01))
    data.frame(
      protein_Id = d$protein_Id[1],
      sample = d[[sname]],
      residual = X$residuals[cbind(i, j)],
      n_peptides = length(features),
      n_samples = length(samples)
    )
  }))
  cells <- dplyr::summarise(
    observations,
    n = dplyr::n(),
    rss = sum(.data$residual^2),
    .by = c("protein_Id", "sample", "n_peptides", "n_samples")
  )
  protein <- dplyr::summarise(
    cells,
    N = sum(.data$n),
    rss = sum(.data$rss),
    .by = c("protein_Id", "n_peptides", "n_samples")
  )
  protein$d <- pmax(protein$N - protein$n_peptides - protein$n_samples + 1, 0)
  protein$s2 <- ifelse(protein$d > 0, protein$rss / protein$d, NA)
  protein$sigma2 <- limma::squeezeVar(ifelse(protein$d > 0, protein$s2, 0), protein$d)$var.post
  cells <- dplyr::left_join(cells, protein[, c("protein_Id", "N", "d", "sigma2")], by = "protein_Id")
  cells$d_cell <- ifelse(cells$N > cells$n_samples, (cells$n - 1) * cells$d / (cells$N - cells$n_samples), 0)
  cells$ratio_raw <- ifelse(cells$d_cell > 0, cells$rss / cells$d_cell / cells$sigma2, NA)
  cells$factor <- 1 / limma::squeezeVar(ifelse(cells$d_cell > 0, cells$ratio_raw, 0), cells$d_cell)$var.post
  cells <- dplyr::mutate(cells, factor_rel = .data$factor / mean(.data$factor), .by = "protein_Id")
  dplyr::rename(cells, !!sname := "sample")
}

# the first `n_proteins` proteins of the IonStar peptides, raw intensities, and their configuration
ionstar_subset <- function(n_proteins) {
  ionstar <- prolfqua_data("data_ionstar")$filtered()
  proteins <- head(unique(ionstar$data$protein_Id), n_proteins)
  list(data = ionstar$data[ionstar$data$protein_Id %in% proteins, ], config = ionstar$config)
}

ionstar_peptides <- function(n_proteins = 80) {
  ionstar <- ionstar_subset(n_proteins)
  suppressMessages(LFQData$new(ionstar$data, ionstar$config)$get_Transformer()$log2()$lfq)
}

aggregated <- function(lfq, method = "medpolish") {
  agg <- lfq$get_Aggregator(method)
  suppressWarnings(suppressMessages(agg$aggregate()))
  agg
}

contrast_a <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

test_that("medpolish aggregation reports the consistency columns of the reference implementation (IonStar)", {
  lfq <- ionstar_peptides()
  agg <- aggregated(lfq)
  result <- agg$lfq_agg$data_long()
  config <- agg$lfq_agg$get_config()
  reference <- reference_consistency(lfq)
  both <- dplyr::inner_join(result, reference, by = c("protein_Id", config$sample_name))
  expect_equal(nrow(both), nrow(reference))
  expect_equal(nrow(both), sum(!is.na(result[[config$get_response()]])))
  expect_lt(max(abs(both$consistency_sigma2 - both$sigma2)), 1e-10)
  expect_lt(max(abs(both$consistency_factor - both$factor)), 1e-10)
  expect_lt(max(abs(both$consistency_factor_rel - both$factor_rel)), 1e-10)
  n <- both[[config$nr_children]]
  confidence <- both[[config$identification_confidence]]
  expect_equal(both$weight_count, n)
  expect_equal(both$weight_count_conf, n * confidence)
  expect_equal(both$weight_protein_ivw, n / both$sigma2)
  expect_equal(both$weight_sample_ivw, n * both$factor / both$sigma2)
  expect_equal(both$weight_sample_ivw_conf, n * both$factor / both$sigma2 * confidence)
  expect_equal(both$weight_sample_rel, n * both$factor_rel)
  expect_equal(both$weight_sample_rel_conf, n * both$factor_rel * confidence)
  expect_named(agg$consistency_priors, c("protein_df_prior", "protein_var_prior", "cell_df_prior", "cell_var_prior"))
})

test_that("the aggregation names its consistency columns in the configuration", {
  lfq <- ionstar_peptides(30)
  for (method in c("medpolish", "rlm")) {
    config <- aggregated(lfq, method)$lfq_agg$get_config()
    expect_equal(config$consistency_factor, "consistency_factor_rel", info = method)
    expect_true(all(config$weight_columns %in% config$value_vars()), info = method)
    expect_setequal(
      config$weight_columns,
      c(
        "consistency_sigma2",
        "consistency_factor",
        "consistency_factor_rel",
        "weight_count",
        "weight_count_conf",
        "weight_protein_ivw",
        "weight_sample_ivw",
        "weight_sample_ivw_conf",
        "weight_sample_rel",
        "weight_sample_rel_conf"
      )
    )
  }
})

test_that("an outlying sample gets a smaller factor, one-peptide proteins factor 1, factors average 1", {
  ionstar <- ionstar_subset(80)
  data <- ionstar$data
  response <- ionstar$config$get_response()
  observed <- data[!is.na(data[[response]]), ]
  # a protein with at least six peptides in every sample; its peptides disagree in one sample
  counts <- dplyr::summarise(observed, n = dplyr::n(), .by = c("protein_Id", "sampleName"))
  target <- counts$protein_Id[counts$n >= 6][1]
  sample <- counts$sampleName[counts$protein_Id == target][1]
  outlier <- which(data$protein_Id == target & data$sampleName == sample & !is.na(data[[response]]))
  data[[response]][outlier] <- data[[response]][outlier] * 2^(4 * rep(c(-1, 1), length.out = length(outlier)))
  single <- counts$protein_Id[!counts$protein_Id %in% target][1]
  first_peptide <- data$peptide_Id[data$protein_Id == single][1]
  data <- data[data$protein_Id != single | data$peptide_Id == first_peptide, ]
  lfq <- suppressMessages(LFQData$new(data, ionstar$config)$get_Transformer()$log2()$lfq)
  for (method in c("medpolish", "rlm")) {
    result <- aggregated(lfq, method)$lfq_agg$data_long()
    observed <- result[!is.na(result$consistency_factor), ]
    target_rows <- observed[observed$protein_Id == target, ]
    expect_lt(target_rows$consistency_factor_rel[target_rows$sampleName == sample], 0.5, label = method)
    single_rows <- observed[observed$protein_Id == single, ]
    expect_gt(nrow(single_rows), 0)
    expect_equal(single_rows$consistency_factor_rel, rep(1, nrow(single_rows)), info = method)
    means <- as.vector(tapply(observed$consistency_factor_rel, observed$protein_Id, mean))
    expect_equal(means, rep(1, length(means)), info = method)
  }
})

test_that("aggregation works when too few proteins have residual df to estimate the prior", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 10, seed = 2)
  data <- istar$data |>
    dplyr::filter(.data$peptide_Id == min(.data$peptide_Id), .by = "protein_Id")
  lfq <- suppressMessages(LFQData$new(data, istar$config)$get_Transformer()$log2()$lfq)
  agg <- aggregated(lfq)
  result <- agg$lfq_agg$data_long()
  observed <- result[!is.na(result$consistency_factor), ]
  expect_true(all(observed$consistency_factor == 1))
  expect_true(all(observed$consistency_factor_rel == 1))
  expect_true(all(is.na(observed$consistency_sigma2)))
  expect_true(all(is.na(unlist(agg$consistency_priors))))
})

test_that("topN and limpa aggregations write no consistency columns", {
  lfq <- ionstar_peptides(20)
  config <- aggregated(lfq, "medpolish")$lfq_agg$get_config()
  expect_length(config$weight_columns, 10)
  raw <- prolfqua_data("data_ionstar")$filtered()
  lfq_raw <- LFQData$new(raw$data[raw$data$protein_Id %in% unique(raw$data$protein_Id)[1:20], ], raw$config)
  topn <- lfq_raw$get_Aggregator("topN")
  suppressMessages(topn$aggregate())
  expect_length(topn$lfq_agg$get_config()$weight_columns, 0)
  expect_null(topn$lfq_agg$get_config()$consistency_factor)
  expect_false("consistency_factor_rel" %in% names(topn$lfq_agg$data_long()))
})

test_that("the default weights multiply in the consistency factor, rlm included", {
  lfq <- ionstar_peptides(60)
  proteins <- aggregated(lfq)$lfq_agg
  config <- proteins$get_config()
  explicit <- c(config$nr_children, config$identification_confidence, "consistency_factor_rel")
  contrasts <- c("b_vs_a" = "dilution.b - dilution.a")
  for (method in c("lm", "limma", "deqms", "rlm")) {
    default <- suppressWarnings(suppressMessages(build_contrast_analysis(
      proteins,
      "~ dilution.",
      contrasts,
      method = method
    )))
    weighted <- suppressWarnings(suppressMessages(build_contrast_analysis(
      proteins,
      "~ dilution.",
      contrasts,
      method = method,
      weights = explicit
    )))
    expect_equal(default$get_contrasts()$p.value, weighted$get_contrasts()$p.value, info = method)
  }
})

test_that("use_peptide_consistency = FALSE gives the default weights without the factor", {
  lfq <- ionstar_peptides(60)
  proteins <- aggregated(lfq)$lfq_agg
  config_off <- proteins$get_config()$clone(deep = TRUE)
  config_off$use_peptide_consistency <- FALSE
  off <- LFQData$new(proteins$data_long(), config_off)
  counts <- c(config_off$nr_children, config_off$identification_confidence)
  contrasts <- c("b_vs_a" = "dilution.b - dilution.a")
  for (method in c("lm", "limma", "deqms", "rlm")) {
    default <- suppressWarnings(suppressMessages(build_contrast_analysis(
      off,
      "~ dilution.",
      contrasts,
      method = method
    )))
    weighted <- suppressWarnings(suppressMessages(build_contrast_analysis(
      off,
      "~ dilution.",
      contrasts,
      method = method,
      weights = counts
    )))
    expect_equal(default$get_contrasts()$p.value, weighted$get_contrasts()$p.value, info = method)
  }
})

test_that("only an explicit use_peptide_consistency = FALSE turns the consistency factor off", {
  config <- AnalysisConfiguration$new()
  expect_null(.used_consistency(config))
  config$consistency_factor <- "consistency_factor_rel"
  expect_equal(.used_consistency(config), "consistency_factor_rel")
  config$use_peptide_consistency <- NULL
  expect_equal(.used_consistency(config), "consistency_factor_rel")
  config$use_peptide_consistency <- FALSE
  expect_null(.used_consistency(config))
})

test_that("the consistency fields survive the configuration round trip", {
  config <- aggregated(ionstar_peptides(20))$lfq_agg$get_config()
  config$use_peptide_consistency <- FALSE
  values <- r6_extract_values(config)
  restored <- r6_extract_values(list_to_analysis_configuration(values))
  expect_equal(restored[sort(names(restored))], values[sort(names(values))])
  expect_equal(values$consistency_factor, "consistency_factor_rel")
  expect_false(values$use_peptide_consistency)
})
