group_stats <- function(nrMeasured, meanAbundance, nrReplicates = 3) {
  data.frame(nrMeasured = nrMeasured, nrReplicates = nrReplicates, meanAbundance = meanAbundance)
}

test_that("the LOD is read from groups seen once when there are any", {
  stats <- group_stats(c(1, 1, 2, 3, 0), c(10, 12, 5, 20, NA))
  expect_equal(.lod_from_group_stats(stats, 0.5), 11)
})

test_that("without groups seen once, the partly observed groups with the fewest values give the LOD", {
  stats <- group_stats(c(2, 2, 3, 4, 0), c(8, 10, 5, 20, NA), nrReplicates = 4)
  expect_equal(.lod_from_group_stats(stats, 0.5), 9)
})

test_that("when every group is complete or empty, the lowest group mean is the LOD", {
  stats <- group_stats(c(3, 3, 0), c(15, 12, NA))
  expect_equal(.lod_from_group_stats(stats, 0.5), 12)
})

test_that("data without any observation has no LOD", {
  expect_true(is.na(.lod_from_group_stats(group_stats(c(0, 0), c(NA, NA)), 0.5)))
})
