test_that("built-in facades are registered with package prolfqua", {
  entry <- lookup_facade("lm")
  expect_equal(entry$class, "ContrastsLMFacade")
  expect_equal(entry$needs, "aggregated")
  expect_equal(entry$package, "prolfqua")
  expect_false(entry$needs_saint_annotation)
})

test_that("lookup_facade returns NULL for unknown names", {
  expect_null(lookup_facade("does-not-exist"))
  expect_null(lookup_facade(""))
  expect_null(lookup_facade(NA_character_))
})

test_that("register_facade and unregister_facade roundtrip", {
  on.exit(unregister_facade("test_fake"), add = TRUE)
  invisible(register_facade(
    "test_fake",
    class = "FakeFacade",
    needs = "aggregated",
    package = "fakepkg",
    needs_saint_annotation = TRUE,
    extra = "hello"
  ))
  entry <- lookup_facade("test_fake")
  expect_equal(entry$class, "FakeFacade")
  expect_equal(entry$package, "fakepkg")
  expect_true(entry$needs_saint_annotation)
  expect_equal(entry$extra, "hello")

  expect_true(unregister_facade("test_fake"))
  expect_null(lookup_facade("test_fake"))
  expect_false(unregister_facade("test_fake"))
})

test_that("list_facades reflects registered entries", {
  on.exit(unregister_facade("test_listed"), add = TRUE)
  register_facade("test_listed", class = "Whatever", needs = "aggregated")
  all_facades <- list_facades()
  expect_true("test_listed" %in% names(all_facades))
  expect_true("lm" %in% names(all_facades))
  expect_equal(all_facades$test_listed$class, "Whatever")
})

test_that("FACADE_REGISTRY snapshot lists built-in entries", {
  expect_true("lm" %in% names(FACADE_REGISTRY))
  expect_true("limma" %in% names(FACADE_REGISTRY))
  expect_equal(FACADE_REGISTRY$lm$class, "ContrastsLMFacade")
  expect_equal(FACADE_REGISTRY$lm$package, "prolfqua")
})

test_that("register_facade accepts only aggregated or nested input", {
  on.exit(unregister_facade("test_needs"), add = TRUE)
  expect_error(
    register_facade("test_needs", class = "Whatever", needs = "same"),
    class = "prolfqua_error_bad_argument"
  )
  expect_null(lookup_facade("test_needs"))
  register_facade("test_needs", class = "Whatever", needs = "nested")
  expect_equal(lookup_facade("test_needs")$needs, "nested")
})

test_that("register_facade refuses the names of built-in methods", {
  expect_error(
    register_facade("lm", class = "Whatever", needs = "nested", package = "other"),
    class = "prolfqua_error_bad_argument"
  )
  expect_identical(lookup_facade("lm")$class, "ContrastsLMFacade")
  expect_identical(lookup_facade("lm")$needs, "aggregated")
})

test_that("built-in facades check their input shape themselves, not through the registry", {
  protein <- sim_lfq_data_protein_config(n_proteins = 10, with_missing = FALSE)
  lfq_protein <- LFQData$new(protein$data, protein$config)
  lfq_protein$is_transformed(TRUE)
  peptide <- sim_lfq_data_peptide_config(n_proteins = 10, with_missing = FALSE)
  lfq_peptide <- LFQData$new(peptide$data, peptide$config)
  lfq_peptide$is_transformed(TRUE)
  contrasts <- c(A_vs_Ctrl = "group_A - group_Ctrl")

  # a registry entry changed behind the class's back no longer changes what the class accepts
  on.exit(prolfqua:::.seed_facade_registry(), add = TRUE)
  entry <- lookup_facade("lm")
  entry$needs <- "nested"
  assign("lm", entry, envir = prolfqua:::.facade_registry_env)
  expect_s3_class(ContrastsLMFacade$new(lfq_protein, "~ group_", contrasts)$get_contrasts(), "data.frame")

  expect_error(ContrastsLMFacade$new(lfq_peptide, "~ group_", contrasts), "ContrastsLMFacade requires aggregated")
  expect_error(
    ContrastsROPECANestedFacade$new(lfq_protein, "~ group_", contrasts),
    "ContrastsROPECANestedFacade requires LFQData with additional hierarchy"
  )
})

test_that("the registry's needs equals the .needs field of each built-in facade class", {
  class_needs <- function(generator) {
    while (!is.null(generator)) {
      if (".needs" %in% names(generator$private_fields)) {
        return(generator$private_fields$.needs)
      }
      generator <- generator$get_inherit()
    }
    NA_character_
  }
  builtin <- Filter(function(facade) facade$package == "prolfqua", list_facades())
  expect_length(builtin, 22L)
  for (name in names(builtin)) {
    generator <- get(builtin[[name]]$class, envir = asNamespace("prolfqua"))
    expect_equal(class_needs(generator), builtin[[name]]$needs, label = name)
  }
})

test_that("facade_styles.tsv gives every built-in method one distinct style", {
  styles <- utils::read.delim(
    system.file("extdata", "facade_styles.tsv", package = "prolfqua", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
  expect_setequal(styles$method, names(.builtin_facade_entries()))
  expect_false(anyDuplicated(styles$method) > 0)
  expect_false(anyDuplicated(styles[c("colour", "shape")]) > 0)
})
