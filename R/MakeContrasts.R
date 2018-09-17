#' create contrast matrix to simplify work with multicomp. (taken from MSqRob)
#' @export
#' @examples
#' make_contrast(c("a- b", "c -d "), levels = c("a","b","c","d") )
#'
make_contrast <- function(contrasts, levels)
{
  n <- length(levels)
  if (n < 1)
    stop("No levels to construct contrasts from.")
  if (is.factor(levels))
    levels <- levels(levels)
  if (!is.character(levels))
    levels <- colnames(levels)

  indicator <- function(i,n) {
    out <- rep(0,n)
    out[i] <- 1
    out
  }

  if (!is.null(contrasts)) {

    #In order to remove invalid level values which are not syntactically valid variable names in R and, for example, do not begin with a letter.
    from <- levels
    to <- make.names(levels)

    gsub2 <- function(pattern, replacement, x, ...) {
      for(i in 1:length(pattern))
        x <- gsub(pattern[i], replacement[i], x, ...)
      x
    }

    e <- gsub2(from, to, contrasts)

    levelsenv <- new.env()
    for (i in 1:n) assign(to[i], indicator(i, n), pos = levelsenv)

    ne <- length(contrasts)
    L <- matrix(0, nrow=length(levels), ncol=length(contrasts), dimnames=list(Levels = levels, Contrasts = contrasts))

    for (j in 1:ne) {
      ej <- parse(text = e[j])
      L[, j] <- eval(ej, envir = levelsenv)
    }

  }

  return(L)
}
