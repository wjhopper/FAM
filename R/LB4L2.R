#' Summarizing the raw LB4L2 dataset by the experimentally manipulated factors
#'
#' The LB4L2 experiment data is stored in 4 different data sets, each corresponding
#' to a separate phase of the LB4L2 experiment (see LB4L2 datasets for more info).
#' This is a method for the S3 generic function \code{\link{summary}} which reports
#' either the subject-level or condition-levels means for each experimental condition.
#'
#' @param data An LB4L data set from the FAM package.
#' @param level A character vector of length one, used for determing the level
#' of the summary. If the value is \code{"subject"} (the default) then the means
#' for each subject are returned. If the value is \code{"group"} then the means
#' of each experimental condition is returned, aggregated over individual subjects.
#' If the value is "raw", then the raw data is returned without any aggregation.
#' This options is useful along with \code{given_practice = TRUE} to retreive the
#' raw joint data.
#' @param given_practice A numeric vector listing which practice tests should be
#' used as grouping variables for calculating conditional final test averages. A logical
#' vector of length one (i.e., \code{TRUE} or \code{FALSE} can be used to indicate
#' all practice tests should be used, or no practice test should be used.
#' @param conditional A logical scalar indicating whether or not analysis should
#' condition on practice test performance. If \code{FALSE}, and given_practice is
#' not \code{FALSE}, a joint analysis of practice and final test outcomes is performed.
#' @param given_data The dataset to use for conditional or joint analysis. If not
#' supplied, the datset specified in the "test_practice" field of the "tables" list
#' attribute of the dataset given as the \code{data} argument.
#'
#' @note Trials where no keys were pressed by the participant are excluded from
#' the RT summary.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr expand_ replace_na
#' @export
summary.LB4L <- function(data, level = c("subject", "group", "raw"),
                         given_practice = FALSE,
                         conditional = FALSE,
                         given_data = NULL) {

  level <- match.arg(tolower(level), c("subject", "group", "raw"))

  if (!all(c("acc","RT") %in% names(data))) {
    stop(paste("Variables 'acc' and 'RT' not found in this data set"))
  }

  if (!is.numeric(given_practice) && !is.logical(given_practice)) {
    stop("'given_practice' must be logical scalars or numeric vectors")
  }

  if (conditional && !any(given_practice > 0)) {
    stop('"given_practice" can not be FALSE if "conditional" is TRUE')
  }

  data <- rename(data, final_acc = acc, final_RT = RT)

  if (any(given_practice > 0)) {

    if (is.null(given_data)) {
      given_data <- get(attr(data,"tables")$test_practice)
    }

    if (is.numeric(given_practice)) {

      if ("round" %in% names(given_data)) {
        rounds <- unique(given_data$round)
      } else {
        given_data$round <- rounds <- 1
      }

      if (!setequal(given_practice, rounds)) {
        given_data <- filter(given_data, round %in% given_practice)
      }
    }
  }

  include <- character(0)
  if (!is.null(given_data)) {
    data <- ungroup(data) %>%
      select(subject:target, final_acc, final_RT) %>%
      left_join(x = reshaping(given_data),
                y = ., by = c("subject", "group", "list", "target")) %>%
      select(subject:list, pracCue = cue.x, finalCue = cue.y,
             target, sameCue, starts_with("practice"), final_acc, final_RT)

    if (conditional) {
      DVname <- "conditional_p"
      include <- grep("^practice[0-9]acc", names(data), value = TRUE)
    }  else {
      DVname <- "joint_p"
    }
  } else {
    DVname <- "mean_p"
  }

  if (level == "raw") {
    return(data)
  }

  orig_groups <- sapply(groups(data), as.character)
  new_groups <- grep("acc$", names(data), value = TRUE)
  all_groups <- n_grouping <- c(orig_groups, new_groups)
  fill_vars <- grep("\\w*practice$|group", all_groups, value=TRUE, invert=TRUE)

  possible_n <- group_by_(data, .dots = c(orig_groups, include)) %>%
    summarise(cell_max = n())
  fill_frame <- do.call(expand.grid, c(lapply(data[,fill_vars],  unique),
                                       list(stringsAsFactors = FALSE))) %>%
    left_join(distinct(select(ungroup(possible_n),subject,group)), by="subject") %>%
    left_join(possible_n, by = names(.)[names(.) %in% names(possible_n)])

  prac_RTcols <- grep("practice[0-9]+RT", names(data), value = TRUE)
  prac_RT_funs <- lapply(prac_RTcols, function(x) {
    as.formula(paste0("~median(", x, ", na.rm = TRUE)"))
    }) %>%
    setNames(prac_RTcols)

  summarized_data <- data %>%
    group_by_(.dots = new_groups, add = TRUE) %>%
    summarise_(freq = ~n(),
               RT = ~median(final_RT, na.rm = TRUE),
               nRT = ~sum(!is.na(final_RT)),
               .dots = prac_RT_funs) %>%
    right_join(y = fill_frame, by = intersect(names(.), names(fill_frame))) %>%
    ungroup() %>%
    mutate(freq = ifelse(is.na(freq) & !is.na(cell_max), 0, freq),
           p = ifelse(is.na(freq), NA_real_, freq/cell_max))

  if (level == "group") {
    summarized_data %<>% group_by_(.dots = all_groups[all_groups != "subject"])
    N <- summarise(summarized_data,
                   nObs = sum(freq, na.rm=TRUE),
                   nRT = sum(nRT, na.rm=TRUE),
                   N = n())
    summarized_data %<>% group_by_(.dots = all_groups[all_groups != "subject"]) %>%
      summarise_each(funs(mean(., na.rm = TRUE), sd(., na.rm = TRUE)), p, RT) %>%
      left_join(N, by = intersect(names(.), names(N))) %>%
      mutate(sem_p = sqrt(p_sd)/sqrt(N),
             sem_RT = sqrt(RT_sd)/sqrt(N)) %>%
      ungroup()
  }

  DVcols <- grep("^p(?:_)", names(summarized_data), value = TRUE, perl = TRUE)
  newDVcols <- sapply(strsplit(DVcols, "_"), function(x) paste(rev(x), collapse = "_"))
  newDVcols[grepl("^mean|^p$", newDVcols)] <- DVname
  RTcols <- grep("^RT", names(summarized_data), value = TRUE)
  newRTcols <- sapply(strsplit(RTcols, "_"), function(x) paste(rev(x), collapse = "_"))

  summarized_data %<>% rename_(.dots = setNames(c(DVcols, RTcols), c(newDVcols,  newRTcols)))
  conditions <- sum(grepl("^practice[0-9]+", names(summarized_data)))
  if (conditions >= 1) {
    summarized_data %<>% as_LB4L_CD_summary()
  } else {
    summarized_data %<>%  as_LB4L_IV_summary()
  }

  return(summarized_data)
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
#' @param data An LB4L2_IV_summary data frame from the FAM package.
#' @param DV A character vector specifying which DV to plot, \code{"accuracy"} or
#' \code{"RT"}
#' @importFrom lazyeval interp
#' @export
autoplot.LB4L_IV_summary <- function(data, DV = "accuracy") {

  y_vars <- list(accuracy = list(name = "mean_p",  ylabel = "Recall Accuracy", lim = 0:1),
                 RT = list(name = "mean_RT", ylabel = "Average Median RT", lim = c(0,6)))
  p <- LB4L_plotbuilder(data, y_vars[[DV]]$name) +
    scale_color_discrete("Practice\nCondition",
                         breaks = c("N.N","N.S","N.T","S.N","T.N"),
                         labels = c(N.N = "Baseline", N.S ="Other Cue Study",
                                    S.N = "Same Cue Study", N.T = "Other Cue Test",
                                    T.N = "Same Cue Test")) +
    scale_y_continuous(y_vars[[DV]]$ylabel, limits = y_vars[[DV]]$lim) +
    ggtitle(paste('Final Test', DV))
  return(p)
}

#' Plot Accuracy or RT for each experimental condition of the LB4L2 dataset
#'
#' @param data An LB4L2_IV_summary data frame from the FAM package.
#' @export
autoplot.LB4L_CD_summary <- function(data, DV = "accuracy") {

  y_vars <- list(joint_p = list(accuracy = list(name = "joint_p",
                                                y = scale_y_continuous("Joint Accuracy", limits = c(-.05,1))),
                                RT = list(name = "mean_RT",
                                          y = scale_y_continuous("Joint RT"))),
                 conditional_p = list(accuracy = list(name = "conditional_p",
                                                      y = scale_y_continuous("Conditional Accuracy")),
                                      RT = list(name = "mean_RT",
                                                y = scale_y_continuous("Conditional RT"))))
  var <- intersect(c("conditional_p", "joint_p"), names(data))

  p <- LB4L_plotbuilder(data, y_vars[[var]][[DV]]$name) +
    scale_color_discrete("Same Cue",
                         breaks = c("yes", "no"),
                         labels = c("yes" = "Yes", "no" = "No")) +
    y_vars[[var]][[DV]]$y +
    ggtitle(paste("Practice & Final Test", DV))
  return(p)
}

#' @export
autoplot.LB4L_joint_summary <- function(data, DV = "accuracy") {

}

LB4L_plotbuilder <- function(data, DV) {

  if (DV %in% c("mean_p", "conditional_p")) {
    data %<>% filter(final_acc == 1)  %>% select(-final_acc)
  }
  vars <- names(data)
  lookup_table <- c("0" = "Incorrect", "1" = "Correct")
  labels_fun <- labeller(practice1acc = setNames(paste("Practice Test", lookup_table), names(lookup_table)),
                         practice2acc = setNames(paste("Practice Test 2", lookup_table), names(lookup_table)),
                         final_acc = setNames(paste("Final Test", lookup_table), names(lookup_table)))

  color_vars <- intersect(c("practice","OCpractice","sameCue"), vars)
  if (length(color_vars) > 1) {
    color_vars <- paste0("interaction(", paste0(color_vars, collapse = ","), ")",
                        collapse = "")
  }

  facet_vars <- paste(grep("^practice[0-9]acc$", names(data), value = TRUE),
                      intersect("final_acc", names(data)),
                      collapse = "+", sep = "~")

  p <- ggplot(data, aes_string(x = "group", y = DV, color = color_vars, group = color_vars)) +
    geom_point(size=2) +
    geom_line(size=.75) +
    scale_x_discrete("Group",expand=c(0,.35),
                     limits = c("immediate", "delay"),
                     labels=c("Immediate","Delay")) +
    theme(legend.key.height = unit(2,"line"))

  if (facet_vars != "") {
    if (grepl("~$", facet_vars)) {
      facet_vars <- paste0("~", gsub("~$", "", facet_vars))
    }
    p <- p + facet_grid(facet_vars, labeller = labels_fun)
  }

  if (!("subject" %in% vars)) {
    err_val <- paste0("sem_",strsplit(DV, "_")[[1]][2])
    vals <- list(x = as.name(DV), y = as.name(err_val))
    data %<>% mutate_(upper = interp(~ x + y, .values = vals),
                      lower = interp(~ x - y, .values = vals))
    p <- p + geom_errorbar(aes(ymax = upper, ymin = lower),
                           data = data, width = .025)
  }
}

#'
grep_for_error <- function(DV, variables) {
  table <- c("probability" = "acc", "RT" = "RT")
  grep(paste0('^sem_', table[DV]), variables,
       value = TRUE, perl = TRUE, ignore.case = TRUE)
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

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_IV_summary"
#' @export
as_LB4L_IV_summary <- function(data, ...) {
  class(data) <- c("LB4L_IV_summary", class(data))
  return(data)
}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_CD_summary"
#' @export
as_LB4L_CD_summary <- function(data, ...) {
  class(data) <- c("LB4L_CD_summary", class(data))
  return(data)
}

#' @describeIn FAM_classes Creates a data frame with the primary class "LB4L_joint_summary"
#' @export
as_LB4L_joint_summary <- function(data, ...) {
  class(data) <- c("LB4L_joint_summary", class(data))
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

#'@describeIn FAM_classes Checks if an R object has \code{LB4L_IV_summary} as it's
#' first value of the \code{class} attribute.
#' #' @export
is.LB4L2_CD_summary <- function(x) {
  return("LB4L_CD_summary" %in% class(x)[1])
}

