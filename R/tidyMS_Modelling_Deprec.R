#' compute group averages
#'
#' used in p2621, p2109
#'
#' @export
#' @keywords internal
#' @examples
#'
#' modelSummary_A <- LFQServiceData::modellingResult_A
#' m <- get_complete_model_fit(modelSummary_A$modelProtein)
#'
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#' factor_contrasts
#'
#' factor_levelContrasts <- contrasts_linfct( m,
#'         factor_contrasts,
#'         subject_Id = "Compound",
#'         contrastfun = LFQService::my_contrast_V2)
#'
#' #usethis::use_data(factor_levelContrasts, overwrite = TRUE)
#'
#' models_interaction <- LFQServiceData::models_interaction
#'
#' m <- get_complete_model_fit(models_interaction$modelProtein)
#' m$linear_model[[1]]
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#' factor_levelContrasts <- contrasts_linfct_deprec( m,
#'                            factor_contrasts,
#'                        subject_Id = "protein_Id")
#' head(factor_levelContrasts)
#' m$linear_model[[1]]
#' my_contest(m$linear_model[[1]], factor_contrasts )
#'
#' plot(factor_levelContrasts$df, factor_levelContrasts$df.residual.model )
#' abline(c(0,1))
#' plot(factor_levelContrasts$df.residual.model , factor_levelContrasts$df - factor_levelContrasts$df.residual.model )
#'
contrasts_linfct_deprec <- function(models,
                             linfct,
                             subject_Id = "protein_Id" ,
                             contrastfun = LFQService::my_contest){
  #computeGroupAverages
  modelcol <- "linear_model"
  models <- models %>% dplyr::filter(exists_lmer == TRUE)


  interaction_model_matrix <- models %>%
    dplyr::mutate(contrast = map(!!sym(modelcol) , contrastfun , linfct = linfct ))


  mclass <- function(x){
    class(x)[1]
  }

  interaction_model_matrix <-  interaction_model_matrix %>%
    dplyr::mutate(classC = purrr::map_chr(contrast,mclass)) %>%
    dplyr::filter(classC != "logical")

  contrasts <- interaction_model_matrix %>%
    dplyr::select_at( c(subject_Id, "contrast") ) %>%
    tidyr::unnest_legacy()

  # take sigma and df from somewhere else.
  modelInfos <- models %>%
    dplyr::select_at(c(subject_Id,
                       "isSingular",
                       "sigma.model" = "sigma",
                       "df.residual.model" = "df.residual" )) %>%

    dplyr::distinct()
  contrasts <- dplyr::inner_join(contrasts, modelInfos, by = subject_Id)

  # adjust
  contrasts <- contrasts %>% group_by_at("lhs") %>%
    dplyr::mutate(p.value.adjusted = p.adjust(p.value, method = "BH")) %>%
    dplyr::ungroup()

  return(contrasts)
}

