# ContrastsPlotter ----
#' plot contrasts
#' @export
#' @family modelling
#' @family plotting
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$normalized()
#' istar$config <- old2new(istar$config )
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelName <- "Model"
#' modelFunction <-
#'   strategy_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' pepIntensity <- istar_data
#' config <- istar$config
#' config$table$hierarchy_keys_depth()
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  modelName = modelName,
#'  subject_Id = config$table$hierarchy_keys_depth())
#'
#'  #mod$get_coefficients()
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.b - dilution.a",
#'   "dil.e_vs_a" = "dilution.e - dilution.a"
#'   ,"dil.e_vs_b" = "dilution.e - dilution.b",
#'   "dil.c_vs_b" = "dilution.c - dilution.b"
#'  )
#' contrast <- prolfqua::Contrasts$new(mod,
#'   Contr)
#' tmp <- contrast$get_contrasts()
#'
#' cp <- ContrastsPlotter$new(tmp ,
#'  contrast$subject_Id,
#' volcano = list(list(score = "FDR", thresh = 0.1)),
#' histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
#'                  list(score = "FDR", xlim = c(0,1,0.05))),
#' score =list(list(score = "statistic",  thresh = 5)))
#' cp$volcano_plotly()
#'
#' cp <- ContrastsPlotter$new(tmp ,
#'  contrast$subject_Id,
#' volcano = list(list(score = "FDR", thresh = 0.1)),
#' histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
#'                  list(score = "FDR", xlim = c(0,1,0.05))),
#'                  fcthresh = NULL,
#' score =list(list(score = "statistic", thresh = 5)))
#' cp$fcthresh
#' cp$volcano_plotly()
#' p <- cp$score_plot(legend=FALSE)
#' cp$score_plotly()
#' p <- cp$histogram()
#' p <- cp$histogram_estimate()
#'
#' res <- cp$volcano()
#' respltly <- cp$volcano_plotly()
#'
#' length(respltly)
#' cp$ma_plot()
#' cp$ma_plotly(rank=TRUE)
#' res  <- cp$barplot_threshold()
#' names(res)
#' cp$histogram_diff()
#' cp$volcano()
ContrastsPlotter <- R6::R6Class(
  "ContrastsPlotter",
  public = list(
    #' @field contrastDF data frame with contrasts
    contrastDF = NULL,
    #' @field modelName of column with model name
    modelName =  character(),
    #' @field subject_Id hierarchy key columns
    subject_Id = character(),
    #' @field prefix default Contrasts - used to generate file names
    prefix = "Contrasts",
    #' @field diff column with fold change differences
    diff = "diff",
    #' @field contrast column with contrasts names, default "contrast"
    contrast = "contrast",
    #' @field volcano_spec volcano plot specification
    volcano_spec = NULL,
    #' @field score_spec score plot specification
    score_spec = NULL,
    #' @field histogram_spec plot specification
    histogram_spec = NULL,
    #' @field fcthresh fold change threshold
    fcthresh = 1,
    #' @field avg.abundance name of column containing avg abundance values.
    avg.abundance = character(),
    #' @description
    #' create Crontrast_Plotter
    #' @param contrastDF frame with contrast data
    #' @param subject_Id columns containing subject Identifier
    #' @param volcano which score to plot and which ablines to add.
    #' @param histogram which scores to plot and which range (x) should be shown.
    #' @param score score parameters
    #' @param fcthresh default 1 (log2 FC threshold)
    #' @param modelName name of column with model names
    #' @param diff fold change (difference) diff column
    #' @param contrast contrast column
    #' @param avg.abundance name of column with average abundance
    initialize = function(contrastDF,
                          subject_Id,
                          volcano = list(list(score = "FDR", thresh = 0.1)),
                          histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                                           list(score = "FDR", xlim = c(0,1,0.05))),
                          score = list(list(score = "statistic",  thresh = NULL)),
                          fcthresh = 1,
                          modelName = "modelName",
                          diff = "diff",
                          contrast = "contrast",
                          avg.abundance = "avgAbd"
    ){
      self$contrastDF <- tidyr::unite(
        contrastDF,
        "subject_Id", subject_Id, sep = "~", remove = FALSE)

      self$modelName  = modelName
      self$subject_Id = subject_Id
      self$diff = diff
      self$volcano_spec = volcano
      self$score_spec = score
      self$histogram_spec = histogram
      self$fcthresh = fcthresh
      self$contrast = contrast
      self$avg.abundance = avg.abundance
    },
    #' @description
    #' plot histogram of selected scores (e.g. p-value, FDR, t-statistics)
    histogram = function(){
      res <- list()
      if (length(self$histogram_spec) > 0) {
        for (spec in self$histogram_spec) {
          fig <- private$.histogram(score = spec )
          res[[spec$score]] <- fig
        }
        return(res)
      }
    },
    #' @description
    #' plot histogram of effect size - difference between groups
    #' @param binwidth with of bin in histogram
    histogram_estimate = function(binwidth = 0.05){
      re <- range(self$contrastDF[[self$diff]], na.rm = TRUE)
      re[1] <- floor(re[1])
      re[2] <- ceiling(re[2])

      fig <- ggplot(self$contrastDF, aes(x = !!sym(self$diff))) +
        geom_histogram(breaks = seq(from = re[1], to = re[2], by = binwidth)) +
        geom_vline(aes(xintercept = median(!!sym( self$diff ), na.rm = T)),   # Ignore NA values for mean
                   color = "red", linetype = "dashed", size = 1) +
        geom_vline(xintercept = 0, col = "green" , size = 1) +
        facet_wrap(vars(!!sym(self$contrast)))

      return(fig)
    },
    #' @description
    #' plot histogram of differences (diff) fold change
    #' @param binwidth with of bin in histogram
    histogram_diff = function(binwidth = 0.05){
      self$histogram_estimate(binwidth = binwidth)
    },
    #' @description
    #' volcano plots (fold change vs FDR)
    #' @param colour column name with color information default modelName
    #' @param legend default TRUE
    #' @param scales default fixed \code{\link{facet_wrap}}, scales argument
    volcano = function(colour,
                       legend = TRUE,
                       scales = c("fixed","free","free_x","free_y")){
      if(missing(colour)){
        colour <- self$modelName
      }
      scales <- match.arg(scales)
      fig <- private$.volcano(self$contrastDF,
                              self$volcano_spec,
                              colour = colour,
                              legend = legend,
                              scales = scales )
      return(fig)
    },
    #' @description
    #' plotly volcano plots
    #' @param colour column in contrast matrix with colour coding
    #' @return list of ggplots
    #' @param legend default TRUE
    #' @param scales default fixed \code{\link{facet_wrap}}, scales argument
    volcano_plotly = function(colour,
                              legend = TRUE,
                              scales = c("fixed","free","free_x","free_y")){
      if(missing(colour)){
        colour <- self$modelName
      }
      scales <- match.arg(scales)
      res <- private$.volcano(self$contrastDF,
                              self$volcano_spec,
                              colour = colour,
                              legend = legend,
                              scales = scales,
                              plotly = TRUE)
      return(res)
    },
    #' @description
    #' ma plot
    #'
    #' MA plot displays the effect size estimate as a function
    #' of the mean protein intensity across groups.
    #' Each dot represents an observed protein.
    #' Red horizontal lines represent the fold-change threshold.
    #'
    #' Sometimes measured effects sizes (differences between samples groups)
    #' are biased by the signal intensity (here protein abundance).
    #' Such systematic effects can be explored using MA-plots.
    #'
    #' @param fc fold change abline
    #' @param colour column in contrast matrix with colour coding
    #' @param legend enable legend default TRUE
    #' @param rank default FALSE, if TRUE then rank of avgAbd is used.
    #' @return ggplot
    ma_plot = function(fc, colour, legend = TRUE, rank = TRUE){
      if ( missing(fc))
        fc <- self$fcthresh
      if (missing(colour)) {
        colour <- self$modelName
      }
      contrastDF <- self$contrastDF
      if (!is.null(contrastDF[[self$avg.abundance]])) {
        # pdf version
        if (rank) {
          rankcol <- paste0("rank_", self$avg.abundance)

          contrastDF <- contrastDF |>
            dplyr::group_by(!!sym(self$contrast)) |>
            mutate(!!rankcol := rank(!!sym(self$avg.abundance)))

          #contrastDF[[ rankcol ]] <- rank(contrastDF[[self$avg.abundance]])
          fig <- private$.ma_plot(
            contrastDF,
            rankcol,
            self$diff,
            self$contrast,
            fc,
            colour = colour,
            legend = legend )
        }else{
          fig <- private$.ma_plot(
            contrastDF,
            self$avg.abundance,
            self$diff,
            self$contrast,
            fc,
            colour = colour,
            legend = legend)
        }
      }else{
        warning("no group_1 group_2 columns can't generate MA")
        fig <- NULL
      }
      return(fig)
    },
    #' @description
    #' ma plotly
    #' @param fc horizontal lines
    #' @param colour column in contrast matrix with colour coding.
    #' @param legend enable legend default TRUE
    #' @param rank default FALSE, if TRUE then rank of avgAbd is used.
    #' @return list of ggplots
    ma_plotly = function(fc, colour, legend = TRUE, rank = FALSE){
      # html version
      if (missing(fc))
        fc <- self$fcthresh
      if(missing(colour))
        colour <- self$modelName
      contrastDF  <- self$contrastDF
      if (!is.null(contrastDF[[self$avg.abundance]])) {
        if(rank){
          rankcol <- paste0("rank_", self$avg.abundance)
          contrastDF <- contrastDF |>
            dplyr::group_by(!!sym(self$contrast)) |>
            mutate(!!rankcol := rank(!!sym(self$avg.abundance)))

          fig <- private$.ma_plot(
            contrastDF,
            rankcol,
            self$diff,
            self$contrast,
            fc,
            colour = colour,
            legend = legend
          )
        } else {
          contrastDF  <- contrastDF |>
            plotly::highlight_key(~subject_Id)
          fig <- private$.ma_plot(
            contrastDF,
            self$avg.abundance,
            self$diff,
            self$contrast,
            fc,
            colour = colour,
            legend = legend
          )
        }

        fig_plotly <- fig |>
          plotly::ggplotly(tooltip = "subject_Id")

        return(fig_plotly)
      }else{
        return(NULL)
      }
    },
    #' @description
    #' plot a score against the log2 fc e.g. t-statistic
    #' @param scorespec list(score="statistics", fcthres = 2, thresh = 5)
    #' @param colour column with colour coding
    #' @param legend enable legend default TRUE
    #' @return list of ggplots
    score_plot = function(scorespec,  colour, legend = TRUE ){
      if (!missing(scorespec)) {
        self$score_spec[[scorespec$score]] <- scorespec
      }
      if(missing(colour))
        colour <- self$modelName
      res <- list()
      if (length(self$score_spec) > 0) {
        res <- private$.score_plot(
          self$contrastDF,
          self$score_spec,
          colour = colour,
          legend = legend )
      }
      return(res)
    },
    #' @description
    #' plot a score against the log2 fc e.g. t-statistic
    #' @param scorespec list(score="statistics", fcthres = 2, thresh = 5)
    #' @param colour column with colour coding
    #' @param legend enable legend default TRUE
    #' @return list of ggplots
    score_plotly = function(scorespec,
                            colour,
                            legend = TRUE ){
      if (!missing(scorespec)) {
        self$score_spec[[scorespec$score]] <- scorespec
      }
      if(missing(colour))
        colour <- self$modelName
      contrastDF <- self$contrastDF |> plotly::highlight_key( ~subject_Id )
      res <- private$.score_plot(
        contrastDF,
        self$score_spec,
        colour = colour,
        legend = legend )

      for (i in seq_along(res)) {
        res[[i]] <- plotly::ggplotly(res[[i]], tooltip = "subject_Id")
      }
      return(res)
    },
    #' @description
    #' shows the number of significant proteins per contrasts
    #' @return list which contains ggplots and summary tables
    barplot_threshold = function(){
      resBar <- list()
      for (i in seq_along(self$volcano_spec) ) {
        scN <- self$volcano_spec[[i]]$score
        scT <- self$volcano_spec[[i]]$thresh
        filt <- dplyr::filter(
          self$contrastDF ,
          !is.na(!!sym(scN)) & !!sym(scN)  < scT)
        if (is.numeric(self$fcthresh)) {
          filt <-  dplyr::filter(filt, abs(!!sym(self$diff)) > self$fcthresh)
        }
        sumC <- group_by(filt, !!sym(self$contrast), !!sym(self$modelName)) |>
          dplyr::summarize(n = n())
        p <- ggplot(sumC,
                    aes(x = !!sym(self$contrast), y = n, fill = !!sym(self$modelName))) +
          geom_bar(position = "stack", stat = "identity")
        resBar[[scN]] <- list(plot = p, summary = sumC)
      }
      return(resBar)
    }
  ),
  private = list(
    .volcano = function(contrasts,
                        scores,
                        colour = NULL,
                        legend = TRUE,
                        scales = "free_y",
                        plotly = FALSE,
                        min_score = 0.0001 ){
      fig <- list()
      for (score in scores) {
        column <- score$score
        contrasts2 <- contrasts |>
          dplyr::filter(!is.na(!!sym(self$diff))) |>
          dplyr::filter(!is.na(!!sym(column))) |>
          dplyr::mutate(!! column := case_when(!!sym(column) < min_score ~ min_score, TRUE ~ !!sym(column)))
        if (plotly) {
          contrasts2 <- contrasts2 |> plotly::highlight_key(~subject_Id)
        }
        p <- prolfqua:::.multigroup_volcano(
          contrasts2,
          effect = self$diff,
          significance = column,
          contrast = self$contrast,
          text = "subject_Id",
          xintercept = if (is.numeric(self$fcthresh)) { c(-self$fcthresh, self$fcthresh) } else {NULL},
          yintercept = score$thresh,
          colour = colour,
          scales = scales)
        if (!legend) {
          p <- p + guides(colour = "none")
        }
        if (plotly) {
          p <-  plotly::ggplotly(p, tooltip = "subject_Id")
        }
        fig[[column]] <- p
      }
      return(fig)
    },
    .histogram  = function(score){
      xlim = score$xlim
      score = score$score
      plot <- self$contrastDF |> ggplot(aes(x = !!sym(score))) +
        geom_histogram(breaks = seq(from = xlim[1], to = xlim[2], by = xlim[3])) +
        facet_wrap(vars(!!sym(self$contrast)))
      return(plot)
    },
    .ma_plot = function(x, avg.abundance, diff, contrast, fc, colour = NULL, legend = TRUE){
      p <- ggplot(x , aes(x = !!sym(avg.abundance),
                          y = !!sym(diff),
                          text = !!sym("subject_Id"),
                          colour = !!sym(colour))) +
        geom_point(alpha = 0.5) +
        scale_colour_manual(values = c("black", "green")) +
        facet_wrap(vars(!!sym(contrast)))
      if(FALSE){ ylab("log fold change (M)") + xlab("mean log intensities (A)") } else { NULL }
      if ( is.numeric(fc) ) {
        p <- p + geom_hline(yintercept = c(-fc, fc), linetype = "dashed", colour = "red")
      }

      if (!legend) {
        p <- p + guides(colour = "none")
      }
      return(p)
    },
    .score_plot = function(x, scores, colour, legend = TRUE){
      fig <- list()
      for (score in scores) {
        xlim = self$fcthresh
        ylim = score$thresh
        score = score$score
        scoreVal <- if ("data.frame" %in% class(x)) {
          x[[score]]
        } else {
          x$data()[[score]]
        }

        ylims <- c( sign(min(scoreVal, na.rm = TRUE)) * ylim, sign( max(scoreVal, na.rm = TRUE)) * ylim)
        p <- ggplot(x, aes(x = !!sym(self$diff),
                           y = !!sym(score),
                           text = !!sym("subject_Id"),
                           colour = !!sym(colour))) +
          scale_colour_manual(values = c("black", "green")) +
          geom_point(alpha = 0.5) +
          facet_wrap(vars(!!sym(self$contrast))) +
          geom_hline(yintercept = c(0), colour = 1) +
          geom_vline(xintercept = c(0), colour = 1 ) +
          geom_hline(yintercept = ylims, colour = 2, linetype = "dashed")

        if (is.numeric(xlim)) {
          p <- p + geom_vline(xintercept = c(-xlim, xlim), colour = 2, linetype = "dashed" )
        }

        if (!legend) {
          p <- p + guides(colour = "none")
        }
        fig[[score]] <- p
      }
      return(fig)
    }
  )
)
