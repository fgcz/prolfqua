# Identification confidence: setup checks, aggregation, and its use as model weight

# give observed rows a random identification confidence and recompute the feature-level column
with_random_confidence <- function(data, keys, seed = 1) {
  withr::local_seed(seed)
  data$identification_confidence <- ifelse(is.na(data$abundance), NA, stats::runif(nrow(data), 0.3, 1))
  dplyr::mutate(
    data,
    identification_confidence_feature = max(.data$identification_confidence, na.rm = TRUE),
    .by = dplyr::all_of(keys)
  )
}

contrast_a <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

test_that("setup_analysis rejects confidences outside [0, 1] and negative child counts", {
  data <- sim_lfq_data(n_proteins = 3, n_per_group = 2, peptide = TRUE)
  data$nr_children <- 1
  config <- sim_lfq_data_peptide_config(n_proteins = 3)$config
  data$identification_confidence <- 1.2
  expect_error(suppressMessages(setup_analysis(data, config)), class = "prolfqua_error_bad_argument")
  data$identification_confidence <- 1
  data$nr_children <- -1
  expect_error(suppressMessages(setup_analysis(data, config)), class = "prolfqua_error_bad_argument")
})

test_that("setup_analysis derives the feature-level columns as the maximum over samples and reports variation", {
  data <- sim_lfq_data(n_proteins = 3, n_per_group = 2, peptide = TRUE)
  data$nr_children <- 1
  withr::local_seed(2)
  data$identification_confidence <- stats::runif(nrow(data), 0.5, 1)
  config <- sim_lfq_data_peptide_config(n_proteins = 3)$config
  withr::local_options(prolfqua.verbose = TRUE)
  messages <- capture_messages(res <- setup_analysis(data, config))
  expect_true(any(grepl("identification_confidence varies across samples within", messages)))
  expect_true(any(grepl("derived identification_confidence_feature, nr_children_feature", messages)))
  expected <- dplyr::summarise(
    res,
    m = max(.data$identification_confidence, na.rm = TRUE),
    f = dplyr::n_distinct(.data$identification_confidence_feature),
    v = unique(.data$identification_confidence_feature),
    .by = c("protein_Id", "peptide_Id")
  )
  expect_true(all(expected$f == 1))
  expect_equal(expected$v, expected$m)
  expect_true(all(res$nr_children_feature == 1))
})

test_that("setup_analysis rejects a supplied feature-level column with several values per feature", {
  data <- sim_lfq_data(n_proteins = 3, n_per_group = 2, peptide = TRUE)
  data$nr_children <- 1
  data$identification_confidence <- 1
  data$identification_confidence_feature <- seq_len(nrow(data)) / nrow(data)
  config <- sim_lfq_data_peptide_config(n_proteins = 3)$config
  expect_error(suppressMessages(setup_analysis(data, config)), class = "prolfqua_error_bad_argument")
})

test_that("complete_cases gives added rows their feature-level values", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 5, with_missing = FALSE)
  data <- with_random_confidence(istar$data, c("protein_Id", "peptide_Id"))
  lfq <- LFQData$new(data[-c(1, 2), ], istar$config)
  lfq$complete_cases()
  res <- lfq$data_long()
  expect_equal(nrow(res), nrow(data))
  per_feature <- dplyr::summarise(
    res,
    n = dplyr::n_distinct(.data$identification_confidence_feature),
    na = sum(is.na(.data$identification_confidence_feature)),
    .by = c("protein_Id", "peptide_Id")
  )
  expect_true(all(per_feature$n == 1 & per_feature$na == 0))
})

test_that("aggregation keeps the maximum confidence of the observed children and the feature-level columns", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 8, seed = 3)
  data <- with_random_confidence(istar$data, c("protein_Id", "peptide_Id"))
  lfq <- LFQData$new(data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  expected <- data |>
    dplyr::filter(!is.na(.data$abundance)) |>
    dplyr::summarise(expected = max(.data$identification_confidence), .by = c("protein_Id", "sample"))
  aggregators <- list(
    medpolish = function(x) x$get_Aggregator("medpolish")$aggregate(),
    rlm = function(x) x$get_Aggregator("rlm")$aggregate(),
    topN = function(x) x$get_Aggregator("topN")$aggregate()
  )
  for (name in names(aggregators)) {
    agg <- suppressMessages(suppressWarnings(aggregators[[name]](lfq)))
    config <- agg$get_config()
    res <- dplyr::inner_join(agg$data_long(), expected, by = c("protein_Id", "sample"))
    expect_equal(res[[config$identification_confidence]], res$expected, info = name)
    per_protein <- dplyr::summarise(
      agg$data_long(),
      n = dplyr::n_distinct(.data[[config$identification_confidence_feature]]),
      m = max(.data[[config$identification_confidence]], na.rm = TRUE),
      f = unique(.data[[config$identification_confidence_feature]]),
      nc = dplyr::n_distinct(.data[[config$nr_children_feature]]),
      .by = "protein_Id"
    )
    expect_true(all(per_protein$n == 1 & per_protein$nc == 1), info = name)
    expect_equal(per_protein$f, per_protein$m, info = name)
  }
})

test_that("with confidence 1 the default weights reproduce the previous weighting", {
  istar <- sim_lfq_data_protein_config(n_proteins = 30, seed = 4)
  lfq <- LFQData$new(istar$data, istar$config)
  lfq$rename_response("transformedIntensity")
  for (method in c("lm", "limma", "lm_impute")) {
    default <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = method)$get_contrasts()
    counts <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = method, weights = "nrPeptides")
    expect_equal(default$diff, counts$get_contrasts()$diff, info = method)
    expect_equal(default$p.value, counts$get_contrasts()$p.value, info = method)
  }
  for (method in c("rlm", "deqms_voom")) {
    # rlm may not converge for one simulated protein; that warning is not under test
    default <- suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrast_a, method = method))$get_contrasts()
    unweighted <- suppressWarnings(build_contrast_analysis(
      lfq,
      "~ group_",
      contrast_a,
      method = method,
      weights = NULL
    ))
    expect_equal(default$diff, unweighted$get_contrasts()$diff, info = method)
    expect_equal(default$p.value, unweighted$get_contrasts()$p.value, info = method)
  }
})

test_that("with confidence 1 the nested facades equal their unweighted fits", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 15, seed = 5)
  lfq <- LFQData$new(istar$data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  methods <- c("ropeca_nested", "ropecaV2_nested")
  if (requireNamespace("lme4", quietly = TRUE) && requireNamespace("lmerTest", quietly = TRUE)) {
    methods <- c(methods, "lmer_nested")
  }
  for (method in methods) {
    default <- suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrast_a, method = method))
    unweighted <- suppressWarnings(build_contrast_analysis(
      lfq,
      "~ group_",
      contrast_a,
      method = method,
      weights = NULL
    ))
    expect_equal(default$get_contrasts()$diff, unweighted$get_contrasts()$diff, info = method)
    expect_equal(default$get_contrasts()$p.value, unweighted$get_contrasts()$p.value, info = method)
  }
})

test_that("the lm facade weights by child count times identification confidence", {
  istar <- sim_lfq_data_protein_config(n_proteins = 20, seed = 6)
  data <- with_random_confidence(istar$data, "protein_Id")
  data$product <- data$nrPeptides * data$identification_confidence
  lfq <- LFQData$new(data, istar$config)
  lfq$rename_response("transformedIntensity")
  fa <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = "lm")
  direct <- build_model(lfq, strategy_lm("transformedIntensity ~ group_", weights = "product"))
  unweighted <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = "lm", weights = "nrPeptides")
  expect_equal(fa$model$get_coefficients(), direct$get_coefficients())
  expect_false(isTRUE(all.equal(fa$get_contrasts()$diff, unweighted$get_contrasts()$diff)))
})

test_that("lm_impute fills every weight column on imputed rows, so the LOD refit keeps them", {
  istar <- sim_lfq_data_protein_config(n_proteins = 40, weight_missing = 0.5, seed = 7)
  data <- with_random_confidence(istar$data, "protein_Id")
  lfq <- LFQData$new(data, istar$config)
  fa <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = "lm_impute")
  res <- fa$get_contrasts()
  imputed <- res[res$estimate_type == "lod_imputed", ]
  expect_gt(nrow(imputed), 0)
  expect_true(all(is.finite(imputed$diff)))
})

test_that("limma weights from several columns are the product of their feature x sample matrices", {
  istar <- sim_lfq_data_protein_config(n_proteins = 20, seed = 7)
  data <- with_random_confidence(istar$data, "protein_Id")
  lfq <- LFQData$new(data, istar$config)
  lfq$rename_response("transformedIntensity")
  strat <- strategy_limma("transformedIntensity ~ group_", weights = c("nrPeptides", "identification_confidence"))
  setup <- .lfqdata_to_elist(lfq, strat$formula)
  w <- .resolve_weights(lfq, strat, setup$annotation)
  expected <- lfq$data_wide(as.matrix = TRUE, value = "nrPeptides")$data *
    lfq$data_wide(as.matrix = TRUE, value = "identification_confidence")$data
  expect_equal(w, expected)
  strat_bad <- strategy_limma("transformedIntensity ~ group_", weights = c("nrPeptides", "sample"))
  expect_error(.resolve_weights(lfq, strat_bad, setup$annotation), class = "prolfqua_error_bad_argument")
})

test_that("rlm and lmer strategies pass the weight product to their backends", {
  withr::local_seed(8)
  x <- data.frame(y = stats::rnorm(12), g = rep(c("a", "b"), each = 6), w1 = 2, w2 = stats::runif(12, 0.2, 1))
  fit <- strategy_rlm("y ~ g", weights = c("w1", "w2"))$model_fun(x)
  expect_equal(coef(fit), coef(MASS::rlm(y ~ g, data = x, weights = w1 * w2, method = "M", maxit = 20, acc = 1e-4)))
  unweighted <- strategy_rlm("y ~ g")$model_fun(x)
  expect_false(isTRUE(all.equal(coef(fit), coef(unweighted))))
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  x$pep <- rep(c("p1", "p2", "p3"), 4)
  fit_lmer <- strategy_lmer("y ~ g + (1 | pep)", weights = c("w1", "w2"))$model_fun(x)
  direct <- suppressMessages(lme4::lmer(y ~ g + (1 | pep), data = x, weights = w1 * w2))
  expect_equal(lme4::fixef(fit_lmer), lme4::fixef(direct), tolerance = 1e-6)
})

test_that("sigma of an rlm fit treats prior weights as inverse variances", {
  withr::local_seed(9)
  x <- data.frame(y = stats::rnorm(10), g = rep(c("a", "b"), each = 5))
  plain <- MASS::rlm(y ~ g, data = x)
  doubled <- MASS::rlm(y ~ g, data = x, weights = rep(4, 10))
  # weights of 4 mean a variance of sigma^2 / 4 for the same residuals: the unit-weight scale doubles, as in lm()
  expect_equal(sigma(doubled), 2 * sigma(plain), tolerance = 1e-6)
})

test_that("weighted median equals the median for equal weights and follows the weights otherwise", {
  for (n in 1:7) {
    x <- stats::rnorm(n)
    expect_equal(.weighted_median(x, rep(1, n)), stats::median(x))
  }
  expect_equal(.weighted_median(c(1, 2, 3), c(1, 1, 5)), 3)
  expect_equal(.weighted_median(c(1, 2, NA), c(1, 1, 1)), 1.5)
  expect_equal(.weighted_median(c(1, 2, 3), c(0, 1, 1)), 2.5)
  expect_true(is.na(.weighted_median(c(1, 2), c(0, 0))))
  expect_equal(.effective_n(rep(1, 4)), 4)
  expect_equal(.effective_n(c(1, 0.5)), 1.5^2 / 1.25)
  # equal weights below 1 keep their count when rounded down
  for (n in 1:12) {
    expect_identical(.effective_n_floor(rep(0.9, n)), as.numeric(n))
    expect_identical(.effective_n_floor(rep(1 / 3, n)), as.numeric(n))
  }
  expect_equal(.effective_n_floor(c(1, 0.3)), 1)
})

test_that("the null of the weighted median matches simulated weighted medians of normals", {
  withr::local_seed(3)
  for (w in list(c(1, 0.3), c(1, 0.6, 0.2), stats::runif(7, 0.2, 1), stats::runif(25, 0.2, 1))) {
    z <- matrix(stats::rnorm(4e4 * length(w)), ncol = length(w))
    q <- stats::quantile(apply(z, 1, .weighted_median, w = w), c(0.95, 0.99))
    expect_true(all(abs(.weighted_median_normal_sf(q, w) - c(0.05, 0.01)) < 0.006), info = length(w))
  }
})

test_that("ROPECA V2 p-values of equally weighted proteins equal the unweighted null", {
  median_z <- c(0.5, 1.5, -2)
  expect_equal(
    .ropeca_v2_weighted_pvalue(median_z, list(rep(1, 3), rep(0.7, 4), 1), 0.2),
    .ropeca_v2_pvalue(median_z, c(3, 4, 1), 0.2)
  )
})

test_that("ROPECA facades weight the peptide fits and combine the peptides unweighted", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 15, seed = 10)
  data <- with_random_confidence(istar$data, c("protein_Id", "peptide_Id"))
  lfq <- LFQData$new(data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  for (method in c("ropeca_nested", "ropecaV2_nested")) {
    weighted <- suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrast_a, method = method))
    unweighted <- suppressWarnings(build_contrast_analysis(
      lfq,
      "~ group_",
      contrast_a,
      method = method,
      weights = NULL
    ))
    expect_null(weighted$contrast$child_weights, info = method)
    res <- weighted$get_contrasts()
    expect_true(all(res$p.value >= 0 & res$p.value <= 1, na.rm = TRUE), info = method)
    expect_false(isTRUE(all.equal(res$diff, unweighted$get_contrasts()$diff)), info = method)
  }
  fa <- suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrast_a, method = "ropeca_nested"))
  peptide_median <- dplyr::summarise(
    fa$contrast$Contrast$get_contrasts(),
    median_diff = stats::median(.data$diff, na.rm = TRUE),
    .by = c("protein_Id", "contrast")
  )
  res <- dplyr::inner_join(fa$get_contrasts(), peptide_median, by = c("protein_Id", "contrast"))
  expect_equal(res$diff, res$median_diff)
})

test_that("ContrastsROPECA and ContrastsROPECAV2 combine peptides with child_weights", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 15, seed = 10)
  data <- with_random_confidence(istar$data, c("protein_Id", "peptide_Id"))
  lfq <- LFQData$new(data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  child_weights <- dplyr::distinct(lfq$data_long()[, c(
    "protein_Id",
    "peptide_Id",
    "identification_confidence_feature"
  )])
  child_weights <- dplyr::rename(child_weights, weight = "identification_confidence_feature")
  model <- build_model(
    lfq,
    strategy_lm(paste(lfq$response(), "~ group_"), weights = "identification_confidence"),
    subject_id = lfq$hierarchy_keys()
  )
  plain <- ContrastsROPECA$new(Contrasts$new(model, contrast_a))$get_contrasts()
  weighted <- ContrastsROPECA$new(Contrasts$new(model, contrast_a), child_weights = child_weights)$get_contrasts()
  expect_false(isTRUE(all.equal(weighted$diff, plain$diff)))
  peptide_contrasts <- ContrastsModerated$new(Contrasts$new(model, contrast_a))
  plain_v2 <- ContrastsROPECAV2$new(peptide_contrasts, rho = 0.1)$get_contrasts()
  weighted_v2 <- ContrastsROPECAV2$new(peptide_contrasts, rho = 0.1, child_weights = child_weights)$get_contrasts()
  expect_false(isTRUE(all.equal(weighted_v2$p.value, plain_v2$p.value)))
  expect_true(all(weighted_v2$p.value >= 0 & weighted_v2$p.value <= 1, na.rm = TRUE))
})

test_that("limpa_nested passes the protein confidence as prior weights", {
  skip_if_not_installed("limpa")
  istar <- sim_lfq_data_peptide_config(n_proteins = 20, seed = 13)
  lfq <- LFQData$new(istar$data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  default <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = "limpa_nested")$get_contrasts()
  unweighted <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = "limpa_nested", weights = NULL)
  expect_equal(default$p.value, unweighted$get_contrasts()$p.value)
  data <- with_random_confidence(istar$data, c("protein_Id", "peptide_Id"))
  lfq_w <- LFQData$new(data, istar$config)$get_Transformer()$log2()$lfq
  weighted <- build_contrast_analysis(lfq_w, "~ group_", contrast_a, method = "limpa_nested")$get_contrasts()
  expect_false(isTRUE(all.equal(weighted$p.value, default$p.value)))
})

test_that("detection counts are expected true detections", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 6, seed = 11, weight_missing = 0.5)
  data <- with_random_confidence(istar$data, c("protein_Id", "peptide_Id"))
  lfq <- LFQData$new(data, istar$config)
  counts <- .summarize_detection_counts(.prepare_detection_lfqdata(lfq))
  expected <- data |>
    dplyr::summarise(
      detected = sum(.data$identification_confidence[!is.na(.data$abundance)]),
      available = sum(.data$identification_confidence_feature),
      .by = c("protein_Id", "sampleName")
    )
  res <- dplyr::inner_join(counts, expected, by = c("protein_Id", "sampleName"), suffix = c("", ".expected"))
  expect_equal(nrow(res), nrow(counts))
  expect_equal(res$detected, res$detected.expected)
  expect_equal(res$available, res$available.expected)
})

test_that("firth detection rows split observed rows by confidence", {
  data <- data.frame(
    id = c("a", "b", "c"),
    bin_resp = c(1L, 1L, 0L),
    conf = c(1, 0.4, NA),
    conf_f = c(1, 0.9, 0.8)
  )
  rows <- .detection_weight_rows(data, "bin_resp", "conf", "conf_f")
  expect_equal(nrow(rows), 4)
  expect_equal(sort(rows$.detection_weight[rows$id == "b"]), c(0.4, 0.5))
  expect_equal(rows$.detection_weight[rows$id == "a"], 1)
  expect_equal(rows$bin_resp[rows$id == "c"], 0L)
  expect_equal(rows$.detection_weight[rows$id == "c"], 0.8)
})

test_that("rfit says that it ignores the identification confidence", {
  skip_if_not_installed("Rfit")
  istar <- sim_lfq_data_protein_config(n_proteins = 10, seed = 12)
  data <- with_random_confidence(istar$data, "protein_Id")
  lfq <- LFQData$new(data, istar$config)
  lfq$rename_response("transformedIntensity")
  withr::local_options(prolfqua.verbose = TRUE)
  expect_message(
    build_contrast_analysis(lfq, "~ group_", contrast_a, method = "rfit"),
    "Rfit takes no observation weights"
  )
})

test_that("use_identification_confidence = FALSE drops the confidence from the default weights", {
  istar <- sim_lfq_data_protein_config(n_proteins = 20, seed = 6)
  data <- with_random_confidence(istar$data, "protein_Id")
  istar$config$use_identification_confidence <- FALSE
  lfq <- LFQData$new(data, istar$config)
  lfq$rename_response("transformedIntensity")
  for (method in c("lm", "limma", "deqms")) {
    off <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = method)$get_contrasts()
    counts <- build_contrast_analysis(lfq, "~ group_", contrast_a, method = method, weights = "nrPeptides")
    expect_equal(off$diff, counts$get_contrasts()$diff, info = method)
    expect_equal(off$p.value, counts$get_contrasts()$p.value, info = method)
  }
  methods <- c("rlm", "deqms_voom", if (requireNamespace("limpa", quietly = TRUE)) "limpa")
  for (method in methods) {
    off <- suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrast_a, method = method))$get_contrasts()
    unweighted <- suppressWarnings(build_contrast_analysis(
      lfq,
      "~ group_",
      contrast_a,
      method = method,
      weights = NULL
    ))
    expect_equal(off$diff, unweighted$get_contrasts()$diff, info = method)
    expect_equal(off$p.value, unweighted$get_contrasts()$p.value, info = method)
  }
})

test_that("use_identification_confidence = FALSE fits the nested facades unweighted", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 15, seed = 10)
  data <- with_random_confidence(istar$data, c("protein_Id", "peptide_Id"))
  istar$config$use_identification_confidence <- FALSE
  lfq <- LFQData$new(data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  methods <- c("ropeca_nested", "ropecaV2_nested")
  if (requireNamespace("lme4", quietly = TRUE) && requireNamespace("lmerTest", quietly = TRUE)) {
    methods <- c(methods, "lmer_nested")
  }
  for (method in methods) {
    off <- suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrast_a, method = method))
    unweighted <- suppressWarnings(build_contrast_analysis(
      lfq,
      "~ group_",
      contrast_a,
      method = method,
      weights = NULL
    ))
    expect_equal(off$get_contrasts()$diff, unweighted$get_contrasts()$diff, info = method)
    expect_equal(off$get_contrasts()$p.value, unweighted$get_contrasts()$p.value, info = method)
  }
})

test_that("use_identification_confidence = FALSE lets the detection models count every observation once", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 6, seed = 11, weight_missing = 0.5)
  data <- with_random_confidence(istar$data, c("protein_Id", "peptide_Id"))
  ones <- data
  ones$identification_confidence <- ifelse(is.na(ones$abundance), NA, 1)
  ones$identification_confidence_feature <- 1
  config_off <- istar$config$clone(deep = TRUE)
  config_off$use_identification_confidence <- FALSE
  off <- LFQData$new(data, config_off)
  counts_off <- .summarize_detection_counts(.prepare_detection_lfqdata(off))
  counts_one <- .summarize_detection_counts(.prepare_detection_lfqdata(LFQData$new(ones, istar$config)))
  expect_equal(counts_off, counts_one)
  expect_true(all(counts_off$detected == round(counts_off$detected)))

  rows <- .detection_weight_rows(data.frame(id = c("a", "b"), bin_resp = c(1L, 0L)), "bin_resp", NULL, NULL)
  expect_equal(rows$.detection_weight, c(1, 1))
  expect_equal(rows$bin_resp, c(1L, 0L))

  off_fit <- suppressWarnings(build_contrast_analysis(off, "~ group_", contrast_a, method = "binomial_nested"))
  one_fit <- suppressWarnings(build_contrast_analysis(
    LFQData$new(ones, istar$config),
    "~ group_",
    contrast_a,
    method = "binomial_nested"
  ))
  expect_equal(off_fit$get_contrasts()$diff, one_fit$get_contrasts()$diff)
  expect_equal(off_fit$get_contrasts()$p.value, one_fit$get_contrasts()$p.value)
})

test_that("use_identification_confidence = FALSE silences the rfit message", {
  skip_if_not_installed("Rfit")
  istar <- sim_lfq_data_protein_config(n_proteins = 10, seed = 12)
  data <- with_random_confidence(istar$data, "protein_Id")
  istar$config$use_identification_confidence <- FALSE
  lfq <- LFQData$new(data, istar$config)
  lfq$rename_response("transformedIntensity")
  withr::local_options(prolfqua.verbose = TRUE)
  expect_no_message(
    build_contrast_analysis(lfq, "~ group_", contrast_a, method = "rfit"),
    message = "Rfit takes no observation weights"
  )
})

test_that("use_identification_confidence = FALSE fits firth as with confidence 1", {
  skip_if_not_installed("logistf")
  istar <- sim_lfq_data_protein_config(n_proteins = 10, seed = 13, weight_missing = 0.5)
  data <- with_random_confidence(istar$data, "protein_Id")
  ones <- data
  ones$identification_confidence <- ifelse(is.na(ones$abundance), NA, 1)
  ones$identification_confidence_feature <- 1
  config_off <- istar$config$clone(deep = TRUE)
  config_off$use_identification_confidence <- FALSE
  fit <- function(data, config) {
    analysis <- suppressWarnings(build_contrast_analysis(LFQData$new(data, config), "~ group_", contrast_a, "firth"))
    dplyr::arrange(analysis$get_contrasts(), .data$protein_Id)
  }
  off <- fit(data, config_off)
  one <- fit(ones, istar$config)
  expect_equal(off$protein_Id, one$protein_Id)
  expect_equal(off$diff, one$diff)
  expect_equal(off$p.value, one$p.value)
})

test_that("only an explicit use_identification_confidence = FALSE turns the confidence off", {
  config <- AnalysisConfiguration$new()
  expect_equal(.used_confidence(config), "identification_confidence")
  config$use_identification_confidence <- NULL
  expect_equal(.used_confidence(config), "identification_confidence")
  config$use_identification_confidence <- FALSE
  expect_null(.used_confidence(config))
})
