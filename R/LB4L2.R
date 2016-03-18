#' Summarizing the raw LB4L2 dataset by the experimentally manipulated factors
#'
#' The LB4L2 experiment data is stored in 4 different data sets, each corresponding
#' to a separate phase of the LB4L2 experiment (see LB4L2 datasets for more info).
#' This is a method for the S3 generic function \code{\link{summary}} which reports
#' either the subject-level or condition-levels means for each experimental condition.
#'
#' @param data LB4L or LB4L2 data set from the FAM package.
#' @param DV A character vector of length one, used for determing the which dependent
#' variable should be summarised in each experimental condition. If the value is
#' \code{"accuracy"} (the default) then the percent of items correctly recalled
#' is returned for each condition. If the value is \code{"RT"}, then the first-press
#' latency in each experimental condition is returned, conditioned on accuracy of
#' the response.
#' @param level A character vector of length one, used for determing the level
#' of the summary. If the value is \code{"subject"} (the default) then the means
#' for each subject are returned, and if the value is \code{"group"} then the means
#' of each experimental condition is returned, aggregated over individual subjects.
#' @param given_practice A numeric vector listing which practice tests should be
#' used as grouping variables for calculating conditional final test averages. A logical
#' vector of length one (i.e., \code{TRUE} or \code{FALSE} can be used to indicate
#' all practice tests should be used, or no practice test should be used.
#' @param given_final A logical vector of length one, indicating whether or not
#' group observations for specific targets according to the accuracy outcome on the final test.
#' Only used when the value of the \code{DV} argument is \code{"RT"}.
#'
#' @note Trials where no keys were pressed by the participant are excluded from
#' the RT summary.
#'
#' @importFrom Hmisc wtd.var
#' @export
summary.LB4L <- function(data, DV = "accuracy", level = "subject",
                         given_practice = FALSE, given_final = (DV == "RT")) {

  level <- tolower(level)
  if ( length(level) > 1 || !(level %in% c("subject", "group")) ) {
    stop('level argument must be either "subject" or "group"')
  }

  if ( length(DV) > 1 || !(DV %in% c("RT", "accuracy"))) {
    stop('DV argument must be either "RT" or "accuracy"')
  }

  if (any(given_practice > 0) && level == "subject" && DV == "accuracy") {
    stop("Cannot calculate conditional accuracy averages at the subject level")
  }

  if (given_final && DV == "accuracy") {
    stop("'given_final' must be FALSE if the DV is accuracy")
  }

  if (DV == "RT") {
    if (!("RT" %in% names(data))) {
      stop(paste(DV, "variable not found in this data set"))
    }
    DVfun <- RT

  } else if (DV == "accuracy") {
    if (!("acc" %in% names(data))) {
      stop(paste(DV, "variable not found in this data set"))
    }
    DVfun <- accuracy
  }

  if (given_final || any(given_practice > 0)) {
    newClass <- as.LB4L_CD_summary
  } else {
    newClass <- as.LB4L_IV_summary
  }

  if (given_final && !all(given_practice > 0)) {
    given <- "final_acc"

  } else if (all(given_practice > 0)) {
    given <- vector(mode = "list",length = 2L)
    names(given) <- c("tp_data", "vars")

    tp <- get(paste0(attr(data,"experiment"),"_tp"))
    if (is.numeric(given_practice) && "round" %in% names(tp)) {
      tp <- filter(tp, round %in% given_practice)
    } else if (identical(given_practice, TRUE) && "round" %in% names(tp)) {
      given_practice <- sort(unique(tp$round))
    }

    given[[1]] <- tp
    given[[2]] <- c(paste0("practice", as.numeric(given_practice),"acc")[given_practice>0],
                  "final_acc"[given_final])
  } else {
    given <- NULL
  }

  data %>%
    DVfun(group_level = (level == 'group'), given = given) %>%
    ungroup() %>% newClass() %>%
    return()
}

# Internal function for de-normalizing the test practice data to make a column
# for each round of test practice
tp_reshaping <- function(tp_data, DV) {
  spread_names <- as.character(unique(LB4L2_tp$round))
  prefixed_names <- paste0("practice", spread_names, DV)
  LB4L2_tp %>%
    select(subject:target, test, round, acc) %>%
    spread(round, acc) %>%
    rename_(.dots=setNames(sapply(spread_names,as.name),
                           prefixed_names)) %>%
    return()
}

# Internal funcion for summarizing accuracy data from LB4L experiments
accuracy <- function(data, group_level, given) {

  if (!is.null(given)) {
    grouping_vars <- c("subject", "group", "sameCue", given$vars)

    if (sum(grepl("practice*", given$vars)) == 1) {
      given$tp_data <- rename(given$tp_data, practice1_acc = acc)

    } else if (sum(grepl("practice*", given$vars)) > 1) {
      given$tp_data <- tp_reshaping(given$tp_data, DV = "acc")
    }

    data <- given$tp_data %>%
      left_join(select(data, subject:target, acc),
                by = c("subject", "group", "list", "target")) %>%
      select(subject:list, pracCue = cue.x, finalCue = cue.y,
             sameCue = test, target, starts_with("practice"), final_acc = acc) %>%
      mutate(sameCue = factor(sameCue, labels = c('no','yes')))

  } else {
    grouping_vars <-  c("subject","group", "practice", "OCpractice")
    data <- rename(data, final_acc=acc)
  }

  summarized_data <- data %>%
    group_by_(.dots = grouping_vars) %>%
    summarise(avg_final_acc = mean(final_acc), nObs = n())

  if (group_level) {
    summarized_data <- summarized_data %>%
      group_by_(.dots = grouping_vars[grouping_vars != "subject"]) %>%
      summarise(sd_final_acc = sd(avg_final_acc),
                w_sd_final_acc = sqrt(wtd.var(avg_final_acc, nObs)),
                w_avg_final_acc = weighted.mean(avg_final_acc, nObs),
                avg_final_acc = mean(avg_final_acc),
                nObs = sum(nObs),
                N = length(unique(subject))) %>%
      mutate(sem_final_acc = sd_final_acc/sqrt(N)) %>%
      select_(.dots = c(grouping_vars[grouping_vars != "subject"],"avg_final_acc",
                        "sd_final_acc", "nObs", "N", "w_avg_final_acc", "w_sd_final_acc"))

  }

  return(summarized_data)
}


# Internal method for summarizing RT data from LB4L experiments
RT <- function(data, group_level, given = NULL) {

  data <- filter(data, is.finite(RT)) %>%
    rename(final_acc=acc, final_RT = RT)

  if (!is.null(given) && is.list(given)) {
    grouping_vars <- c("subject", "group", "sameCue", given$vars)
    given$tp_data <- given$given$tp_data %>%
      filter(is.finite(RT))

    if (sum(grepl("practice*", given$vars)) == 1) {
      given$tp_data <- rename(given$tp_data, practice1_acc = acc)

    } else if (sum(grepl("practice*", given$vars)) > 1) {
      given$tp_data <- full_join(tp_reshaping(given$tp_data, DV = "acc"),
                                 tp_reshaping(given$tp_data, DV = "RT"))
    }

    data <- given$tp_data %>%
      left_join(select(data, subject:target, acc, RT),
                by = c("subject", "group", "list", "target")) %>%
      select(subject:list, pracCue = cue.x, finalCue = cue.y, ameCue = test,
             target, starts_with("practice"), final_acc, final_RT) %>%
      mutate(sameCue = factor(sameCue, labels = c('no','yes')))

  } else if (!is.null(given) && !is.list(given)) {
    grouping_vars <- c("subject", "group", "practice", "OCpractice", given)
  } else {
    grouping_vars <-  c("subject","group", "practice", "OCpractice")
  }

  if (group_level) {
    summarized_data <- data %>%
      group_by_(.dots = grouping_vars[grouping_vars != "subject"]) %>%
      summarise(medianRT = median(final_RT),
                MAD = mad(final_RT, constant = 1),
                nObs = n(),
                N = length(unique(subject))) %>%
      select_(.dots = c(grouping_vars[grouping_vars != "subject"],
                        "medianRT"," MAD", "nObs"," N"))
  } else {
    summarized_data <- data %>%
      group_by_(.dots = grouping_vars) %>%
      summarise(medianRT = median(final_RT),
                nObs = n())
  }

  return(summarized_data)
}


#' Plot Accuracy or RT for each experimental condition of the LB4L2 dataset
#'
#' @param data An LB4L2_IV_summary data frame from the FAM package.
#' @export
autoplot.LB4L_IV_summary <- function(data) {

  if ("medianRT" %in% names(data)) {
    p <- ggplot(data = as.data.frame(data),
                aes_string(x = "group", y = "medianRT", color = "cond", group = "cond")) +
      geom_point(size=3, shape = 2) +
      ylab("Median First-Press Latency")

  } else if ("avg_final_acc" %in% names(data)) {
    p <- ggplot(data = as.data.frame(data),
                aes_string(x = "group", y = "avgAcc", color = "cond", group = "cond")) +
      geom_point(size=3) +
      geom_errorbar(aes(ymax = avgAcc + sem, ymin = avgAcc - sem),
                    width = .025) +
      ylab("Accuracy")

  } else {
    stop("Did not recognize any DV columns in this LB4L_IV_summary object")
  }

  p <- p +  geom_line(size=1) +
    scale_color_discrete("Practice\nCondition",
                         breaks = c("C.C","C.S","C.T","S.C","T.C"),
                         labels = c(C.C = "Baseline", C.S ="Other Cue Study", S.C  = "Restudy",
                                    C.T = "Other Cue Test", T.C = "Test Same Cue")) +
    scale_x_discrete("Group",expand=c(0,.25),
                     limits = c("immediate", "delay"),
                     labels=c("Immediate","Delay")) +
    theme(legend.key.height = unit(2,"line")) +
    ggtitle('Final Test')

  return(p)

}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L"
#' @export
as.LB4L <- function(data, ...) {
  data <- as.data.frame(data)
  class(data) <- c("LB4L", class(data))
  return(data)
}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_IV_summary"
#' @export
as.LB4L_IV_summary <- function(data, ...) {
  data <- as.data.frame(data)
  class(data) <- c("LB4L_IV_summary", class(data))
  return(data)
}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_CD_summary"
#' @export
as.LB4L_CD_summary <- function(data, ...) {
  data <- as.data.frame(data)
  class(data) <- c("LB4L_CD_summary", class(data))
  return(data)
}

#' @describeIn FAM_classes Checks if an R object has \code{LB4L_IV} as it's
#' first value of the \code{class} attribute.
#' @export
is.LB4L2 <- function(x) {
  return("LB4L" %in% class(x)[1])
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
