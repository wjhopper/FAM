#' Summarizing the raw LB4L2 dataset by the experimentally manipulated factors
#'
#' The LB4L2 experiment data is stored in 4 different data sets, each corresponding
#' to a separate phase of the LB4L2 experiment (see LB4L2 datasets for more info).
#' This is a method for the S3 generic function \code{\link{summary}} which reports
#' either the subject-level or condition-levels means for each experimental condition.
#'
#' @param data LB4L or LB4L2 data set from the FAM package.
#' @param level A character vector of length one, used for determing the level
#' of the summary. If the value is \code{"subject"} (the default) then the means
#' for each subject are returned, and if the value is \code{"group"} then the means
#' of each experimental condition is returned, aggregated over individual subjects.
#' @param DV A character vector of length one, used for determing the which dependent
#' variable should be summarised in each experimental condition. If the value is
#' \code{"accuracy"} (the default) then the percent of items correctly recalled
#' is returned for each condition. If the value is \code{"RT"}, then the first-press
#' latency in each experimental condition is returned, conditioned on accuracy of
#' the response.
#'
#' @note Trials where no keys were pressed by the participant are excluded from
#' the RT summary.
#'
#' @export
summary.LB4L_IV <- function(data, level = "subject", DV = "accuracy") {

  level <- tolower(level)
  if ( length(level) > 1 || !(level %in% c("subject", "group")) ) {
    stop('level argument must be either "subject" or "group"')
    group_level_summary <- level == 'group'
  }

  if ( length(DV) > 1 || !(DV %in% c("RT", "accuracy"))) {
    stop('DV argument must be either "RT" or "accuracy"')
  }

  if (DV == "RT") {
    DVfun <- IV_RT
  } else {
    DVfun <- IV_acc
  }

  data %>%
    DVfun(group_level = (level == 'group')) %>%
    ungroup() %>%
    mutate(cond = interaction(practice, OCpractice)) %>%
    as.LB4L_IV_summary() %>%
    return()

}

# Internal method for summarizing accuracy data from LB4L experiments
IV_acc <- function(data, group_level) {

  summarized_data <- data %>%
    group_by(group, practice, OCpractice, subject) %>%
    summarise(avgAcc = mean(acc),
              nObs = n())
  if (group_level) {
    summarized_data <- summarized_data %>%
      summarise(sdAcc = sd(avgAcc),
                avgAcc = mean(avgAcc),
                nObs = sum(nObs),
                N = length(unique(subject))) %>%
      mutate(sem = sdAcc/sqrt(N)) %>%
      select(group, practice, OCpractice, avgAcc, sdAcc, sem, nObs, N)
  }

  return(summarized_data)
}

# Internal method for summarizing RT data from LB4L experiments
IV_RT <- function(data, group_level) {

  data <- filter(data, is.finite(RT))
  summarized_data <- data %>%
    group_by(group, practice, OCpractice, acc, subject) %>%
    summarise(medianRT = median(RT),
              nObs = n())
  if (group_level) {
    summarized_data <- data %>%
      group_by(group, practice, OCpractice, acc) %>%
      summarise(medianRT = median(RT),
                MAD = mad(RT, constant = 1),
                nObs = n(),
                N = length(unique(subject))) %>%
      select(group, practice, OCpractice, acc, medianRT, MAD, nObs, N)
  }
  return(summarized_data)
}


#' Plot Accuracy or RT for each experimental condition of the LB4L2 dataset
#'
#' @param data An LB4L2_IV_summary data frame from the FAM package.
#' @export
autoplot.LB4L_IV_summary <- function(data) {

    if ("medianRT" %in% names(data)) {
      acc_labeller <- function(labels, ...) { list(c("Incorrect", "Correct")[labels$acc+1]) }
      p <- ggplot(data = as.data.frame(data),
                  aes_string(x = "group", y = "medianRT", color = "cond", group = "cond")) +
        geom_point(size=3, shape = 2) +
        geom_line(size=1) +
        facet_grid(~acc, labeller = acc_labeller) +
        scale_color_discrete("Practice\nCondition",
                             breaks = c("C.C","C.S","C.T","S.C","T.C"),
                             labels = c(C.C = "Baseline", C.S ="Other Cue Study", S.C  = "Restudy",
                                        C.T = "Other Cue Test", T.C = "Test Same Cue")) +
        scale_x_discrete("Group",expand=c(0,.25), limits = c("immediate", "delay"),
                         labels=c("Immediate","Delay")) +
        ylab("Median First-Press Latency") +
        theme(legend.key.height = unit(2,"line")) +
        ggtitle('Final Test')

  } else if ("avgAcc" %in% names(data)) {
    p <- ggplot(data = as.data.frame(data),
                aes_string(x = "group", y = "avgAcc", color = "cond", group = "cond")) +
      geom_point(size=3) +
      geom_line(size=1) +
      geom_errorbar(aes(ymax = avgAcc + sem,
                        ymin = avgAcc - sem),
                    width = .025) +
      scale_color_discrete("Practice\nCondition",
                           breaks = c("C.C","C.S","C.T","S.C","T.C"),
                           labels = c(C.C = "Baseline", C.S ="Other Cue Study", S.C  = "Restudy",
                                      C.T = "Other Cue Test", T.C = "Test Same Cue")) +
      scale_x_discrete("Group",expand=c(0,.25), limits = c("immediate", "delay"),
                       labels=c("Immediate","Delay")) +
      ylab("Accuracy") +
      theme(legend.key.height = unit(2,"line")) +
      ggtitle('Final Test')
  } else {
    stop("Did not recognize any DV columns in this LB4L_IV_summary object")
  }

  return(p)


}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_IV"
#' @export
as.LB4L_IV <- function(data, ...) {
  data <- as.data.frame(data)
  class(data) <- c("LB4L_IV", class(data))
  return(data)
}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_IV_summary"
#' @export
as.LB4L_IV_summary <- function(data, ...) {
  data <- as.data.frame(data)
  class(data) <- c("LB4L_IV_summary", class(data))
  return(data)
}

#' @describeIn FAM_classes Checks if an R object has \code{LB4L_IV} as it's
#' first value of the \code{class} attribute.
#' @export
is.LB4L2_IV <- function(x) {
  return("LB4L_IV" %in% class(x)[1])
}


#'@describeIn FAM_classes Checks if an R object has \code{LB4L_IV_summary} as it's
#' first value of the \code{class} attribute.
#' #' @export
is.LB4L2_IV_summary <- function(x) {
  return("LB4L_IV_summary" %in% class(x)[1])
}


#' Summarizing the raw LB4L2 dataset given practice test performance
#'
#' The LB4L2 experiment data is stored in 4 different data sets, each corresponding
#' to a separate phase of the LB4L2 experiment (see LB4L2 datasets for more info). This
#' is an S3 method for the generic function \code{\link{summary}} which returns a data frame reporting
#' either the subject-level or condition-levels data from the final test in each experimental condition,
#' given recall accuracy on one or both of the practice tests.
#'
#'
#' @param data
#' @param level
#'
#' @return
#' @export
#' @examples
# summary.LB4L2_CD<- function(data, level) {
#
# }
