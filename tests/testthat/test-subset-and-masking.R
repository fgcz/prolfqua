# get_subset() needs a subject-id column; data columns must not shadow function arguments inside filter().

.masking_lfq <- function() {
  istar <- sim_lfq_data_peptide_config(Nprot = 10, with_missing = FALSE, seed = 1)
  LFQData$new(istar$data, istar$config)
}

.with_column <- function(lfq, name, value) {
  copy <- lfq$get_copy()
  data <- copy$data_long()
  data[[name]] <- value
  copy$set_data(data)
  copy
}

test_that("get_subset() errors when x has no subject-id column and keeps partial keys", {
  lfq <- .masking_lfq()
  expect_error(lfq$get_subset(data.frame(foo = 1)), class = "prolfqua_error_missing_column")
  protein <- lfq$data_long()$protein_Id[1]
  sub <- lfq$get_subset(data.frame(protein_Id = protein))
  expect_setequal(unique(sub$data_long()$protein_Id), protein)
  expect_equal(nrow(sub$data_long()), sum(lfq$data_long()$protein_Id == protein))
})

test_that("a data column named threshold does not replace the threshold argument", {
  lfq <- .masking_lfq()
  expected <- lfq$get_copy()
  expected$remove_small_intensities(4)
  shadowed <- .with_column(lfq, "threshold", 1e6)
  shadowed$remove_small_intensities(4)
  expect_equal(nrow(shadowed$data_long()), nrow(expected$data_long()))

  data <- lfq$data_long()
  plain <- filter_proteins_by_peptide_count(data, lfq$get_config(), threshold = 2)
  data$threshold <- 1e6
  masked <- filter_proteins_by_peptide_count(data, lfq$get_config(), threshold = 2)
  expect_equal(nrow(masked$data), nrow(plain$data))
})

test_that("a data column named N does not replace the top-N argument", {
  lfq <- .masking_lfq()
  sum_na <- function(x, name = FALSE) if (name) "sum" else sum(x, na.rm = TRUE)
  top_n <- function(l) {
    ranked <- rank_peptide_by_intensity(l$data_long(), l$response(), l$hierarchy_keys())
    aggregate_intensity_top_n(ranked, l, .func = sum_na, N = 3)$data
  }
  expect_equal(nrow(top_n(.with_column(lfq, "N", 0))), nrow(top_n(lfq)))
})
