% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{matrix_to_tibble}
\alias{matrix_to_tibble}
\title{matrix or data.frame to tibble (taken from tidyquant)}
\usage{
matrix_to_tibble(x, preserve_row_names = "row.names", ...)
}
\arguments{
\item{x}{a matrix}

\item{preserve_row_names}{give name to rownames column, if NULL discard rownames}

\item{...}{further parameters passed to as_tibble}
}
\description{
matrix or data.frame to tibble (taken from tidyquant)
}
\examples{
x <- matrix(rnorm(20), ncol=4)
rownames(x) <- LETTERS[seq_len(nrow(x))]
matrix_to_tibble(x)
!(is.matrix(x) || is.data.frame(x))

}
\seealso{
Other utilities: 
\code{\link{INTERNAL_FUNCTIONS_BY_FAMILY}},
\code{\link{get_UniprotID_from_fasta_header}()},
\code{\link{multigroup_volcano}()},
\code{\link{names_to_matrix}()},
\code{\link{pairs_smooth}()},
\code{\link{pairs_w_abline}()},
\code{\link{panel_cor}()},
\code{\link{panel_hist}()},
\code{\link{remove_NA_rows}()},
\code{\link{table_facade}()}
}
\concept{utilities}
\keyword{internal}
