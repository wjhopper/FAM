#' @import dplyr
#' @importFrom lazyeval lazy_dots
#'
#' @title Hierarchical Summary
#' @description Summarise a data frame at cascading levels of granularity
#'
#' @param collapse
#'  A character vector of variables to group the data frame by,
#'  which will be succesively collapsed over.
#' @param hold
#'  A character vector of variables to consistently group the data frame by,
#'  even while other variables are collapsed
#' @param rawData
#'  The data frame to summmarise iteratively
#' @param ...
#'  Name-value pairs of summary functions like mean(), sd() etc.
#'  See details for warnings about naming conventions
#' @return
#'  A list of tbl_df's, the same length as the character vector passed as the
#'  \code{collapse} argument
#' @export
#'
#' @examples
#' derp <- heirarchicalSummary(collapse = c("subject","group"),
#'                             hold = c("practice","other_type",
#'                                      "prac_score", "other_prac_acc"),
#'                             rawData=ungroup(cbind(LB4L_allSs,n=1)),
#'                             weighted.mean(final_score,n), mean(final_score),
#'                             n = sum(n))
#' str(derp)
#'
heirarchicalSummary <- function(collapse, hold,
                                rawData, ...) {

  summaryFcn <- function(x,  ...) {
    group_by_(x, .dots = lapply(...,as.symbol)) %>%
      summarise_(.dots=fcns)
  }

  fcns <- lazyeval::lazy_dots(...)
  unnamed_fncs <- sapply(names(fcns),identical,"")
  if (any(unnamed_fncs)) {
    fun_arg_pairs <- lapply(lapply(fcns,
                                   `[[`, 'expr'),
                            as.character)
    proposed_names <- unlist(lapply(fun_arg_pairs,`[[`, 2))
    dup_names <- duplicated(proposed_names)
    proposed_names[dup_names] <- ''
    proposed_names[!unnamed_fncs] <- names(fcns)[!unnamed_fncs]
    fcns <- setNames(fcns,proposed_names)
  }
  fcns <- c(fcns, lazyeval::lazy_dots(groupsize=n()))
  index <- mapply(seq, 1:length(collapse), length(collapse))
  gvars <- lapply(index,function(x) c(collapse[x], hold))
  out <- Reduce(summaryFcn, x = c(list(rawData),gvars), accumulate = TRUE)
  out <- setNames(out[-1], collapse)
  return(out)
}



#' @title Performance Bins
#' @description Summarise a variables by groups and bin the observations
#'
#' @param data
#'  data frame to work with
#' @param bin_by
#'  Character vector defining the variables to group the data by
#' @param cutpoints
#'  numeric vector of value to bin the measured variable at. Values not within the range
#'  of this vector will be grouped into the NA category. By default, data is binned
#'  into quartiles.
#' @param ...
#'  Name-value pairs of summary functions like mean(), sd() etc.
#' @return
#'  An ungrouped table df summmarising the observed count and frequency in each bin
#' @export
#'
#' @examples
#' filter(LB4L_allSs, list==1) %>%
#' performanceBins(bin_by = "subject",
#'                 cutpoints=c(0,.25, .75,1),
#'                 acc=mean(final_score))
#'
performanceBins <- function(data, bin_by,
                            cutpoints = c(0,.25,.5,.75,1), ...) {

  fcns <- lazyeval::lazy_dots(...)
  binned <- data %>%
    group_by_(.dots = bin_by) %>%
    summarise_(.dots=fcns) %>%
    mutate_each_(funs(cut(., breaks = cutpoints, include.lowest=TRUE)),
                 vars = names(fcns)) %>%
    group_by_(.dots=names(fcns)) %>%
    summarise(count=n()) %>%
    ungroup() %>%
    mutate(percentage = count/sum(count))
  return(binned)

}
#
#
#   tm <- as.data.frame(xtabs( ~ bin,tailoring_means)) %>%
#     mutate(Freq = Freq/sum(Freq))

#   tailor_plot <- ggplot(tm, aes(x=bin,y=Freq)) +
#     geom_bar(stat='identity') +
#     scale_y_continuous("Proportion of Subjects", limits = c(0,max(tm$Freq+.025))) +
#     geom_text(aes(label =c("6 Secs","5 secs","4 Secs"), y= Freq+.025)) +
#     xlab("Practice List Performance") +
#     theme_larger() +
#     ggtitle('Percent of Subjects with Specific Performance Levels')
#
#   if (return.plot & return.data) {
#     return(list(data = tm, figure = tailor_plot))
#   } else if (return.plot & !return.data) {
#     return(tailor_plot)
#   } else if  (!return.plot & return.data) {
#     return(tm)
#   } else {
#     return(invisible())
#   }
# }
