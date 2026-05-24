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
