#' Null-split calibration check of normalization and model
#'
#' Splits the samples within every stratum at random into two pseudo-groups, A
#' and B, and tests B against A with [build_contrast_analysis()]. No feature
#' differs between the pseudo-groups, so a calibrated pipeline gives uniform
#' p-values: about 5% below 0.05, 1% below 0.01, and no feature with FDR < 0.05.
#' The check uses no labels beyond the strata, so it compares normalizations
#' without favouring any of them by construction.
#'
#' @param lfqdata the model input: a normalized [LFQData] object, aggregated or
#'   not as `method` requires.
#' @param method facade name passed to [build_contrast_analysis()].
#' @param n_splits number of distinct random splits.
#' @param strata factor keys within whose levels the samples are split; the
#'   model is `~ <strata> + pseudo_group_`. `character(0)` splits all samples
#'   as one group; subset `lfqdata` to one condition first.
#' @param n_features number of randomly chosen top-level features
#'   (`hierarchy_keys()[1]`) to test; `NULL` (default) tests all.
#' @param seed `NULL` (default) uses the current random number generator state;
#'   otherwise the draws run under [withr::with_seed()].
#' @param ... further arguments to [build_contrast_analysis()].
#' @return A tibble with one row per split: `split`, `n` (tests with a
#'   p-value), `p_below_0.05`, `p_below_0.01` (shares of p-values),
#'   `n_FDR_below_0.05`, `ks_D` (Kolmogorov-Smirnov distance of the p-values
#'   from the uniform distribution), and the list-columns `assignment` (sample
#'   name and `pseudo_group_` per sample) and `contrasts` (subject keys,
#'   `diff`, `p.value` and `FDR`).
#' @details
#' Within each stratum the samples are split into halves; an odd number of
#' samples gives halves that differ by one. Splits that differ only by swapping
#' A and B count as one. When the design has fewer distinct splits than
#' `n_splits`, all of them are returned with a message.
#'
#' Read the median and range over splits. Splits share samples, so they are not
#' independent replicates. A split that follows a batch or the run order puts a
#' real technical effect into the contrast: compare `assignment` with the batch
#' of each sample, and persistent hits point to a batch effect the
#' normalization left in. With few samples per pseudo-group, variance
#' moderation and imputation behave differently than in the real comparison;
#' calibration carries over, power does not.
#' @seealso [build_contrast_analysis()]
#' @family stats
#' @export
#' @examples
#' bb <- sim_lfq_data_protein_config(n_proteins = 50, with_missing = FALSE)
#' lfqdata <- LFQData$new(bb$data, bb$config)
#' lfqdata$is_transformed(TRUE)
#' res <- null_split_check(lfqdata, method = "limma", n_splits = 5, seed = 1)
#' res[, c("split", "n", "p_below_0.05", "n_FDR_below_0.05", "ks_D")]
null_split_check <- function(
  lfqdata,
  method = "lm",
  n_splits = 20,
  strata = lfqdata$relevant_factor_keys(),
  n_features = NULL,
  seed = NULL,
  ...
) {
  if (!inherits(lfqdata, "LFQData")) {
    abort_bad_argument("lfqdata", "be an LFQData object")
  }
  if (!is.numeric(n_splits) || length(n_splits) != 1L || n_splits < 1) {
    abort_bad_argument("n_splits", "be a single number >= 1", not = toString(n_splits))
  }
  unknown <- setdiff(strata, lfqdata$factor_keys())
  if (length(unknown) > 0L) {
    abort_bad_argument(
      "strata",
      "name factor keys of `lfqdata`",
      not = paste(toString(unknown), "; factor keys:", toString(lfqdata$factor_keys()))
    )
  }
  sample_name <- lfqdata$sample_name()
  samples <- dplyr::distinct(lfqdata$data_long(), dplyr::across(dplyr::all_of(c(sample_name, strata))))
  samples <- samples[order(samples[[sample_name]]), , drop = FALSE]
  stratum <- if (length(strata) == 0L) {
    rep("all", nrow(samples))
  } else {
    do.call(paste, c(as.list(samples[strata]), sep = "~"))
  }
  stratum_size <- table(stratum)
  if (any(stratum_size < 2L)) {
    abort_bad_argument(
      "strata",
      "define strata with at least two samples each",
      not = paste("strata with one sample:", toString(names(stratum_size)[stratum_size < 2L]))
    )
  }

  draw <- function() {
    data <- lfqdata$data_long()
    if (!is.null(n_features)) {
      top <- lfqdata$hierarchy_keys()[1]
      ids <- unique(data[[top]])
      data <- data[data[[top]] %in% sample(ids, min(n_features, length(ids))), , drop = FALSE]
    }
    list(data = data, splits = .draw_null_splits(stratum, n_splits))
  }
  drawn <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  if (length(drawn$splits) < n_splits) {
    inform_prolfqua(
      "null_split_check: the design has only ",
      length(drawn$splits),
      " distinct splits; returning all of them."
    )
  }

  config <- lfqdata$get_config()$clone(deep = TRUE)
  config$factors[["pseudo_group_"]] <- "pseudo_group_"
  modelstr <- paste("~", paste(c(strata, "pseudo_group_"), collapse = " + "))
  contrast <- c(B_vs_A = "pseudo_group_B - pseudo_group_A")

  rows <- lapply(seq_along(drawn$splits), function(i) {
    assignment <- tibble::tibble(samples[sample_name], pseudo_group_ = drawn$splits[[i]])
    data <- drawn$data
    data$pseudo_group_ <- assignment$pseudo_group_[match(data[[sample_name]], assignment[[sample_name]])]
    lfq <- LFQData$new(data, config)
    lfq$is_transformed(lfqdata$is_transformed())
    fa <- build_contrast_analysis(lfq, modelstr, contrast, method = method, ...)
    ct <- dplyr::ungroup(fa$get_contrasts())
    keys <- intersect(lfq$hierarchy_keys(), colnames(ct))
    ct <- ct[, c(keys, "diff", "p.value", "FDR")]
    p <- ct$p.value[!is.na(ct$p.value)]
    tibble::tibble(
      split = i,
      n = length(p),
      p_below_0.05 = mean(p < 0.05),
      p_below_0.01 = mean(p < 0.01),
      n_FDR_below_0.05 = sum(ct$FDR < 0.05, na.rm = TRUE),
      ks_D = unname(suppressWarnings(stats::ks.test(p, "punif"))$statistic),
      assignment = list(assignment),
      contrasts = list(ct)
    )
  })
  dplyr::bind_rows(rows)
}

#' Draw distinct random half-splits within strata
#'
#' @param stratum stratum label per sample
#' @param n_splits number of splits wanted
#' @return list of character vectors ("A"/"B" per sample); fewer than
#'   `n_splits` when the design has fewer distinct splits. Splits equal up to
#'   swapping A and B are duplicates.
#' @keywords internal
#' @noRd
.draw_null_splits <- function(stratum, n_splits) {
  splits <- list()
  keys <- character(0)
  for (attempt in seq_len(100L * n_splits)) {
    split <- character(length(stratum))
    for (s in unique(stratum)) {
      idx <- which(stratum == s)
      split[idx] <- sample(rep(c("A", "B"), length.out = length(idx)))
    }
    canonical <- if (split[1] == "A") split else c(A = "B", B = "A")[split]
    key <- paste(canonical, collapse = "")
    if (!key %in% keys) {
      keys <- c(keys, key)
      splits[[length(splits) + 1L]] <- unname(split)
      if (length(splits) == n_splits) {
        break
      }
    }
  }
  splits
}
