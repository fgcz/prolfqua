# ModelInterface -----
#' R6 interface class representing modelling result
ModelInterface <- R6::R6Class(
  "ModelInterface",
  public = list(
    #' @description
    #' return model coefficients
    #' @return data.frame
    get_coefficients = function(){ stop("get_coefficients not implmeneted")},
    #' @description
    #' perform ANOVA analysis
    #' @return data.frame
    get_anova = function(){stop("get_anova not implmeneted")},
    #' @description
    #' plot histogram of coefficients
    #' @return ggplot
    coef_histogram = function(){stop("coef_histogram not implmeneted")},
    #' @description
    #' plot volcano of coefficients
    #' @return ggplot
    coef_volcano = function(){stop("coef_volcano not implmeneted")},
    #' @description
    #' pairs plot of coefficients
    #' @return ggplot
    coef_pairs = function(){stop("coef_pairs not implmeneted")},
    #' @description
    #' histogram of p-values and FDR for anova results
    #' @return ggplot
    anova_histogram = function(){stop("anova_histogram not implmeneted")}
  )
)
