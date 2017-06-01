#' Summarizing the raw LB4L2 dataset by the experimentally manipulated factors
#'
#' The LB4L2 experiment data is stored in 4 different data sets, each corresponding
#' to a separate phase of the LB4L2 experiment (see LB4L2 datasets for more info).
#' This is a method for the S3 generic function \code{\link{summary}} which reports
#' either the subject-level or condition-levels means for each experimental condition.
#'
#' @param data An LB4L data set from the FAM package.
#' @param level A character vector of length one, used for determing the level
#' of the summary. If the value is \code{"subject"} (the default) then the summary statistics
#' for each subject are returned. If the value is \code{"group"} then the means
#' of each experimental condition is returned, pooled across individual subjects.

#' @note Trials where no keys were pressed by the participant are excluded from
#' the RT summary.
#'
#' @importFrom tidyr replace_na
#' @export
summary.LB4L <- function(data, level = "subject") {

  level <- match.arg(tolower(level), c("subject"))

  if (!all(c("acc","RT") %in% names(data))) {
    stop(paste("Variables 'acc' and 'RT' not found in this data set"))
  }

  data$logRT <- log(data$RT)

  orig_groups <- sapply(groups(data), as.character)
  all_groups <-  c(orig_groups, "acc")
  fill_vars <- grep("\\w*practice$|group", all_groups, value=TRUE, invert=TRUE)

  possible_n <- group_by_(data, .dots = orig_groups) %>%
    summarise(cell_max = n())
  fill_frame <- do.call(expand.grid, c(lapply(data[,fill_vars],  unique),
                                       list(stringsAsFactors = FALSE))) %>%
    left_join(distinct(ungroup(possible_n), subject, group),
              by="subject") %>%
    left_join(possible_n, by = intersect(names(.), names(possible_n)))

  summarized_data <- data %>%
    group_by_(.dots = all_groups, add = TRUE) %>%
    summarise(RT_median = median(RT, na.rm = TRUE),
              RT_mean = mean(RT, na.rm = TRUE),
              logRT_mean = mean(logRT, na.rm = TRUE),
              N_acc = n(),
              N_RT = sum(!is.na(RT))
    ) %>%
    right_join(y = fill_frame, by = intersect(names(.), names(fill_frame))) %>%
    ungroup() %>%
    replace_na(list(N_acc = 0, N_RT = 0)) %>%
    mutate(p = N_acc/cell_max) %>%
    arrange(subject)

  return(as_LB4L_summary(summarized_data))
}

#' @export
summary.LB4L_summary <- function(data) {

  data <- filter(data, acc == 1) %>%
    WISEsummary(dependentvars = c("p", "RT_mean", "logRT_mean", "RT_median"),
                betweenvars = "group",
                withinvars = c("practice","OCpractice"),
                idvar = "subject",
                na.rm = TRUE) %>%
    rename_(.dots=setNames(names(.),
                           gsub("_mean$", "", names(.))
                           )
            ) %>%
    select(-RT_median_n, -logRT_mean_n) %>%
    rename(RT_n = RT_mean_n)
}

#' @export
recode_conditions <- function(data) {

  x <- paste0(data$practice,data$OCpractice)
  data$Condition <- "N"
  data$Condition[grepl("S", x)] <- "S"
  data$Condition[grepl("T", x)] <- "T"

  data$cue_type <- NA
  same <- data$OCpractice =="N" & data$practice != "N"
  data$cue_type[same] <- "Same Cue"
  other <- data$OCpractice !="N" & data$practice == "N"
  data$cue_type[other] <- "Other Cue"

  data <- select(data, -practice, -OCpractice)
  data <- rename(data, practice = Condition)
  return(data)
}


# Internal function for de-normalizing the test practice data to make a column
# for each round of test practice
#' @importFrom tidyr spread_
reshaping <- function(long_data) {

  long_data <- select(long_data,subject:sameCue, round:acc)
  spread_names <- unique(long_data$round)

  if (length(spread_names) == 1) {
    oldnames <- grep("acc|RT", names(long_data),value = TRUE)
    newnames <- setNames(sapply(oldnames, as.name),
                         paste0("practice", spread_names, oldnames))
    wide_data <- rename_(long_data, .dots = newnames)

  } else if (length(spread_names) > 1) {

    acc <- select(long_data, -RT) %>%
      spread(round, acc) %>%
      rename_(.dots=setNames(sapply(spread_names, as.name), paste0("practice", spread_names,"acc")))
    RT <- select(long_data, -acc) %>%
      spread(round, RT) %>%
      rename_(.dots=setNames(sapply(spread_names, as.name), paste0("practice", spread_names,"RT")))
    wide_data <- full_join(acc, RT, by=c("subject","group","list","cue","target","sameCue"))
  }

  return(wide_data)
}

#' Plot Accuracy or RT for each experimental condition of the LB4L2 dataset
#'
#' @param data An LB4L2_summary data frame from the FAM package.
#' @param DV A character vector specifying which DV to plot, \code{"accuracy"} or
#' \code{"RT"}
#' @importFrom lazyeval interp
#' @export
autoplot.LB4L_summary <- function(data, DV = "accuracy", ...) {

  if (DV == c("accuracy")) {
    data %<>% filter(final_acc == 1)%>% select(-final_acc)
  }

  y_vars <- list(accuracy = list(name = "mean_p",  ylabel = "Recall Accuracy",
                                 lim =c(min(data$mean_p)-.05,max(data$mean_p)+.05)),
                 RT = list(name = "mean_RT", ylabel = "Average Median RT",
                           lim = c(min(data$mean_RT)-.2,max(data$mean_RT)+.2)))


  p <- base_plot(data, aes_string(x = "group", y = y_vars[[DV]]$name,
                                  color = "interaction(practice, OCpractice)",
                                  group = "interaction(practice, OCpractice)"),
                 ...) +
    scale_color_discrete("Practice\nCondition",
                         breaks = c("N.N","N.S","N.T","S.N","T.N"),
                         labels = c(N.N = "Baseline", N.S ="Other Cue Study",
                                    S.N = "Same Cue Study", N.T = "Other Cue Test",
                                    T.N = "Same Cue Test")) +
    scale_y_continuous(y_vars[[DV]]$ylabel, limits = y_vars[[DV]]$lim) +
    ggtitle(paste('Final Test', DV))

  if ("final_acc" %in% names(data)){
    p <- p + facet_grid(~final_acc, labeller = labels_fun)
  }

  p <- error_bars(p, data, y_vars[[DV]]$name, ...)
  return(p)
}

#' Plot Accuracy or RT for each experimental condition of the LB4L2 dataset
#'
#' @param data An LB4L2_CD_summary data frame from the FAM package.
#' @export
autoplot.LB4L_CD_summary <- function(data, DV = "accuracy", ...) {

  var <- intersect(c("conditional_p", "joint_p"), names(data))
  if (var == c("conditional_p")) {
    data %<>% filter(final_acc == 1)%>% select(-final_acc)
  }

  y_vars <- list(joint_p = list(accuracy = list(name = "joint_p",
                                                y = scale_y_continuous("Joint Accuracy",
                                                                       limits = c(-.05,.7))),
                                RT = list(name = "mean_RT",
                                          y = scale_y_continuous("Joint RT"))),
                 conditional_p = list(accuracy = list(name = "conditional_p",
                                                      y = scale_y_continuous("Conditional Accuracy")),
                                      RT = list(name = "mean_RT",
                                                y = scale_y_continuous("Conditional RT"))))


  p <-  p <- base_plot(data, aes_string(x = "group", y = y_vars[[var]][[DV]]$name,
                                        color = "sameCue", group = "sameCue"),
                       ...)  +
    scale_color_discrete("Same Cue",
                         breaks = c("yes", "no"),
                         labels = c("yes" = "Yes", "no" = "No")) +
    y_vars[[var]][[DV]]$y +
    ggtitle(paste("Practice & Final Test", DV))

  facet_vars <- paste(grep("^practice[0-9]acc$", names(data), value = TRUE),
                      intersect("final_acc", names(data)),
                      collapse = "+", sep = "~")

  if (facet_vars != "") {
    if (grepl("~$", facet_vars)) {
      facet_vars <- paste0("~", gsub("~$", "", facet_vars))
    }
    p <- p + facet_grid(facet_vars, labeller = labels_fun)
  }

  p <- error_bars(plot = p, data = data, DV = y_vars[[var]][[DV]]$name, ...)

  return(p)
}

#'
labels_fun <- labeller(practice1acc = setNames(paste("Practice Test", c("Incorrect","Correct")),
                                               c("0","1")),
                       practice2acc = setNames(paste("Practice Test 2",c("Incorrect","Correct")),
                                               c("0","1")),
                       final_acc = setNames(paste("Final Test", c("Incorrect","Correct")),
                                            c("0","1")))
#'
base_plot <- function(data, aesthetics, size = 2) {
  p <- ggplot(data, aesthetics) +
    geom_point(size = size) +
    geom_line(size= .375*size) +
    scale_x_discrete("Group",expand=c(0,.15),
                     limits = c("immediate", "delay"),
                     labels=c("Immediate","Delay")) +
    theme(legend.key.height = unit(2,"line"))
  return(p)
}

#'
error_bars <- function(plot, data, DV, size = 2, prefix = "sem_") {
  if (!("subject" %in% names(data))) {
    err_val <- paste0("sem_",strsplit(DV, "_")[[1]][2])
    vals <- list(x = as.name(DV), y = as.name(err_val))
    data %<>% mutate_(upper = interp(~ x + y, .values = vals),
                      lower = interp(~ x - y, .values = vals))
    plot <- plot + geom_errorbar(aes(ymax = upper, ymin = lower),
                                 data = data, width = .05, size = .125*size)
  }
  return(plot)
}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L"
#'
#' @param data An R object that is, or can be coerced to, a data frame.
#' @param tables A list of named character vectors which give the name of the datasets
#' holding observations from each experimental phase. The datasets can take any names,
#' but the names of each item in the list must be one of \code{"study", "study_practice",
#' "test_practice"} or \code{"final"}.
#' @param IDvars A Character vector specifying variables that can be used to group
#' the data (see \code{\link{group_by}}). Defaults to the variables given by any
#' present grouping.
#' @param experiment A character vector giving the name of the experiment. Used to name
#' figures and any prefix any saved Rda/Rds objects.
#' @param practice_tests A numeric scalar giving the number of practice tests
#' taken.
#' @export
as_LB4L <- function(data, tables = attr(data, "tables"),
                    IDvars = vapply(groups(data), as.character, character(1L)),
                    experiment = attr(data, "experiment"),
                    practice_tests = attr(data, "practice_tests")) {
  mostattributes(data) <- c(attributes(data),
                            list(tables = tables, experiment = experiment,
                                 IDvars = IDvars, practice_tests = practice_tests))
  class(data) <- c("LB4L", class(data))
  return(data)
}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_summary"
#' @export
as_LB4L_summary <- function(data, ...) {
  class(data) <- c("LB4L_summary", class(data))
  return(data)
}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_CD_summary"
#' @export
as_LB4L_CD_summary <- function(data, ...) {
  class(data) <- c("LB4L_CD_summary", class(data))
  return(data)
}

#' @describeIn FAM_classes Checks if an R object has \code{LB4L} as it's
#' first value of the \code{class} attribute.
#' @export
is.LB4L <- function(x) {
  return("LB4L" %in% class(x)[1])
}

#'@describeIn FAM_classes Checks if an R object has \code{LB4L_summary} as it's
#' first value of the \code{class} attribute.
#' #' @export
is.LB4L_summary <- function(x) {
  return("LB4L_summary" %in% class(x)[1])
}

#'@describeIn FAM_classes Checks if an R object has \code{LB4L_CD_summary} as it's
#' first value of the \code{class} attribute.
#' #' @export
is.LB4L_CD_summary <- function(x) {
  return("LB4L_CD_summary" %in% class(x)[1])
}

