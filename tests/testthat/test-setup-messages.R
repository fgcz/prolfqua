# setup_analysis() announces every transformation it applies (srr RE2.0, G2.9).

setup_messages <- function(data, config) {
  messages <- character(0)
  withCallingHandlers(
    setup_analysis(data, config),
    prolfqua_message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning")
  )
  messages
}

raw_setup_input <- function() {
  data <- data.frame(
    run = rep(c("s1", "s2", "s3", "s4"), each = 3),
    protein = rep(c("P1", "P1", "P2"), 4),
    source = "db",
    peptide = rep(c("a", "b", "c"), 4),
    group = factor(rep(c("A", "A", "B", "B"), each = 3), levels = c("B", "A", "C"), ordered = TRUE),
    extra = 1,
    intensity = c(10, 11, 12, 10, 11, NA, 13, 14, 15, 13, 14, 15)
  )
  data <- data[-6, ]
  config <- AnalysisConfiguration$new()
  config$file_name <- "run"
  config$hierarchy[["protein_Id"]] <- c("protein", "source")
  config$hierarchy[["peptide_Id"]] <- "peptide"
  config$factors[["group_"]] <- "group"
  config$set_response("intensity")
  list(data = data, config = config)
}

test_that("a factor named by several columns holds their united labels", {
  withr::local_options(prolfqua.verbose = TRUE)
  input <- raw_setup_input()
  input$data$group <- as.character(input$data$group)
  input$data$batch <- c(s1 = "x", s2 = "y", s3 = "x", s4 = "y")[input$data$run]
  input$config$factors[["group_"]] <- c("group", "batch")
  messages <- setup_messages(input$data, input$config)
  expect_true(any(grepl("`group_`: columns `group`, `batch` united", messages)))
  data <- suppressWarnings(suppressMessages(setup_analysis(input$data, input$config)))
  samples <- dplyr::distinct(data, .data$run, .data$group_)
  expected <- paste(c("A", "A", "B", "B"), c("x", "y", "x", "y"), sep = input$config$sep)
  expect_equal(samples$group_[match(c("s1", "s2", "s3", "s4"), samples$run)], expected)
})

test_that("each setup transformation is announced", {
  withr::local_options(prolfqua.verbose = TRUE)
  input <- raw_setup_input()
  messages <- setup_messages(input$data, input$config)
  expect_true(any(grepl("`protein_Id`: columns `protein`, `source` united", messages)))
  expect_true(any(grepl(
    "Factor `group_`: ordered column `group` kept as an unordered factor with levels B, A",
    messages
  )))
  expect_true(any(grepl("creating sampleName", messages)))
  expect_true(any(grepl("dropped [0-9]+ unconfigured columns: .*extra", messages)))
  expect_true(any(grepl("added 1 rows with NA response", messages)))
  input$data$group <- rep(c(1L, 1L, 2L, 2L), each = 3)[-6]
  messages <- setup_messages(input$data, input$config)
  expect_true(any(grepl("Factor `group_`: integer column `group` converted to character labels", messages)))
})

test_that("input already in the configured form gets no transformation messages", {
  withr::local_options(prolfqua.verbose = TRUE)
  input <- raw_setup_input()
  data <- input$data
  data$group <- as.character(data$group)
  data$extra <- NULL
  data$source <- NULL
  data <- rbind(data, data.frame(run = "s2", protein = "P2", peptide = "c", group = "A", intensity = 9))
  data$sampleName <- data$run
  data$qValue <- 0
  data$nr_children <- 1
  config <- input$config
  config$hierarchy[["protein_Id"]] <- "protein"
  expect_length(setup_messages(data, config), 0L)
})
