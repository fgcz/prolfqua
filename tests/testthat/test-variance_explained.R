# summarize_variance_explained() and pooled_variance_explained()

.ve_lfq <- function(data, config) {
  lfq <- LFQData$new(data, config)
  lfq$is_transformed(TRUE)
  lfq
}

.ve_r2_lm <- function(lfq, ve, rhs) {
  d <- lfq$data_long()
  key <- lfq$hierarchy_keys()
  ids <- unique(ve[[key[length(key)]]])
  vapply(
    ids,
    function(id) {
      s <- d[d[[key[length(key)]]] == id, , drop = FALSE]
      summary(stats::lm(stats::as.formula(paste(lfq$response(), rhs)), data = s))$r.squared
    },
    numeric(1)
  )
}

test_that("one factor: R2 equals lm R-squared per feature, with missing values", {
  bb <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30)
  lfq <- .ve_lfq(bb$data, bb$config)
  ve <- suppressMessages(summarize_variance_explained(lfq))
  expect_setequal(unique(ve$factor), c("group_", "Residuals"))
  r2 <- ve$R2[match(unique(ve$protein_Id), ve$protein_Id)]
  expect_equal(r2, unname(.ve_r2_lm(lfq, ve, "~ group_")), tolerance = 1e-10)
})

test_that("two factors: default formula is the factorial; R2 equals the one-way eta2 over cells", {
  bb <- prolfqua::sim_lfq_data_2factor_config(n_proteins = 30)
  bb$config$factor_depth <- 2
  lfq <- .ve_lfq(bb$data, bb$config)
  ve <- suppressMessages(summarize_variance_explained(lfq))
  fk <- lfq$relevant_factor_keys()
  expect_identical(unique(ve$factor), c(fk, paste(fk, collapse = ":"), "Residuals"))
  r2 <- ve$R2[match(unique(ve$protein_Id), ve$protein_Id)]
  cells <- paste0("~ interaction(", paste(fk, collapse = ", "), ")")
  expect_equal(r2, unname(.ve_r2_lm(lfq, ve, cells)), tolerance = 1e-10)
})

test_that("sums of squares add up and shares lie in [0, 1]", {
  bb <- prolfqua::sim_lfq_data_2factor_config(n_proteins = 30)
  bb$config$factor_depth <- 2
  ve <- suppressMessages(summarize_variance_explained(.ve_lfq(bb$data, bb$config)))
  total <- tapply(ve$Sum.Sq, ve$protein_Id, sum)
  expect_equal(unname(total), unname(tapply(ve$SS_total, ve$protein_Id, unique)), tolerance = 1e-10)
  expect_true(all(ve$eta2 >= 0 & ve$eta2 <= 1))
  expect_true(all(ve$R2 >= 0 & ve$R2 <= 1))
})

test_that("pooled R2 equals PERMANOVA R2 with Euclidean distances on complete features", {
  bb <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, with_missing = FALSE)
  lfq <- .ve_lfq(bb$data, bb$config)
  pooled <- pooled_variance_explained(summarize_variance_explained(lfq))
  r2_pooled <- 1 - pooled$eta2[pooled$factor == "Residuals"]

  wide <- tidyr::pivot_wider(
    lfq$data_long(),
    id_cols = dplyr::all_of(c(lfq$sample_name(), "group_")),
    names_from = "protein_Id",
    values_from = dplyr::all_of(lfq$response())
  )
  x <- as.matrix(wide[, -(1:2)])
  d2 <- as.matrix(stats::dist(x))^2
  n <- nrow(x)
  ss_total <- sum(d2[upper.tri(d2)]) / n
  ss_within <- sum(vapply(
    split(seq_len(n), wide$group_),
    function(i) {
      s <- d2[i, i, drop = FALSE]
      sum(s[upper.tri(s)]) / length(i)
    },
    numeric(1)
  ))
  expect_equal(r2_pooled, 1 - ss_within / ss_total, tolerance = 1e-10)
})

test_that("peptide-level data is fitted per peptide, not pooled per protein", {
  bb <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 5)
  lfq <- .ve_lfq(bb$data, bb$config)
  ve <- suppressMessages(summarize_variance_explained(lfq))
  expect_true(all(lfq$hierarchy_keys() %in% names(ve)))
  observed <- dplyr::count(dplyr::filter(lfq$data_long(), !is.na(.data[[lfq$response()]])), .data$peptide_Id)
  df_total <- tapply(ve$Df, ve$peptide_Id, sum)
  expect_equal(as.vector(df_total) + 1, observed$n[match(names(df_total), observed$peptide_Id)])
  r2 <- ve$R2[match(unique(ve$peptide_Id), ve$peptide_Id)]
  expect_equal(r2, unname(.ve_r2_lm(lfq, ve, "~ group_")), tolerance = 1e-10)
})

test_that("a custom modelstr is used", {
  bb <- prolfqua::sim_lfq_data_2factor_config(n_proteins = 10)
  bb$config$factor_depth <- 2
  lfq <- .ve_lfq(bb$data, bb$config)
  fk <- lfq$relevant_factor_keys()
  ve <- suppressMessages(summarize_variance_explained(lfq, paste(lfq$response(), "~", fk[1])))
  expect_identical(unique(ve$factor), c(fk[1], "Residuals"))
})

test_that("features without a usable fit are left out with a message", {
  bb <- prolfqua::sim_lfq_data_protein_config(n_proteins = 10, with_missing = FALSE)
  first <- bb$data$protein_Id == unique(bb$data$protein_Id)[1] & bb$data$group_ == "A"
  bb$data$abundance[first] <- NA
  lfq <- .ve_lfq(bb$data, bb$config)
  expect_message(ve <- summarize_variance_explained(lfq), "1 of 10 features left out", class = "prolfqua_message")
  expect_length(unique(ve$protein_Id), 9)
})

test_that("untransformed data and saturated models stop", {
  bb <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5, with_missing = FALSE)
  expect_error(
    summarize_variance_explained(LFQData$new(bb$data, bb$config)),
    class = "prolfqua_error_bad_argument"
  )
  lfq <- .ve_lfq(bb$data, bb$config)
  expect_error(
    summarize_variance_explained(lfq, paste(lfq$response(), "~", lfq$sample_name())),
    "No feature has a full-rank fit",
    class = "prolfqua_error"
  )
})

test_that("pooled_variance_explained requires the columns of summarize_variance_explained", {
  expect_error(pooled_variance_explained(data.frame(factor = "a")), class = "prolfqua_error_missing_column")
})

test_that("R2 rises after removing per-sample shifts", {
  bb <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, with_missing = FALSE)
  set.seed(1)
  samples <- unique(bb$data$sampleName)
  shift <- stats::setNames(stats::rnorm(length(samples), sd = 1), samples)
  bb$data$abundance <- bb$data$abundance + shift[bb$data$sampleName]
  lfq <- .ve_lfq(bb$data, bb$config)
  normalized <- lfq$get_Transformer()$robscale()$lfq
  r2 <- function(x) {
    ve <- summarize_variance_explained(x)
    stats::median(ve$R2[match(unique(ve$protein_Id), ve$protein_Id)])
  }
  expect_gt(r2(normalized), r2(lfq))
})
