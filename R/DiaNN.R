#' read DiaNN diann-output.tsv file
#'
#' filter for 2 peptides per protein, and for Q.Value < 0.01 (default)
#'
#' @export
#'
read_diann_output <- function(path, nrPeptides = 2, Q.Value = 0.01){
  add_nr_pep <- function(report){
    nrPEP <- report |>
      dplyr::select(rlang::.data$Protein.Group, rlang::.data$Stripped.Sequence) |>
      dplyr::distinct() |>
      dplyr::group_by(rlang::.data$Protein.Group) |>
      dplyr::summarize(nrPeptides = dplyr::n())
    report <- dplyr::inner_join(report, nrPEP)
    return(list(nrPEP = nrPEP, report = report))
  }

  select_PG <- function(report){
    columns  <- c("File.Name",
                  "Protein.Group",
                  #Protein.Ids,
                  "Protein.Names",
                  "PG.Quantity",
                  "PG.MaxLFQ",
                  "PG.Q.Value",
                  "Global.PG.Q.Value",
                  "Lib.PG.Q.Value",
                  "nrPeptides")
    columns %in% names(report)

    PG <- report |> dplyr::select_at( columns
    ) |> dplyr::distinct() |> dplyr::ungroup()
    return(PG)
  }

  filter_PG <- function(PG, nrPeptides = 2, Q.Value = 0.01){
    PG <- PG |> dplyr::filter(nrPeptides >= nrPeptides)
    PG <- PG |> dplyr::filter(rlang::.data$Lib.PG.Q.Value < Q.Value)
    PG <- PG |> dplyr::filter(rlang::.data$PG.Q.Value < Q.Value)
  }

  report <- read.table(path, sep = "\t", header = TRUE)
  rNR <- add_nr_pep(report)
  PG <- select_PG(rNR$report)
  PG2 <- filter_PG(PG, nrPeptides = nrPeptides, Q.Value = Q.Value)
  report2 <- inner_join(PG2, report)
  return(report2)
}
