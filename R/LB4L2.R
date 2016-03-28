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
#' for each subject are returned. If the value is \code{"group"} then the means
#' of each experimental condition is returned, aggregated over individual subjects.
#' If the value is "raw", then the raw data is returned without any aggregation.
#' This options is useful along with \code{given_practice = TRUE} to retreive the
#' raw joint data.
#' @param given_practice A numeric vector listing which practice tests should be
#' used as grouping variables for calculating conditional final test averages. A logical
#' vector of length one (i.e., \code{TRUE} or \code{FALSE} can be used to indicate
#' all practice tests should be used, or no practice test should be used.
#' @param given_final A logical vector of length one, indicating whether or not
#' group observations for specific targets according to the accuracy outcome on the final test.
#' If this argument is set to \code{TRUE}, given_practice is set to \code{TRUE},
#' and the DV is set to \code{"accuracy"}, then the joint distribution of practice and final
#' test outcomes is returned
#'
#' @note Trials where no keys were pressed by the participant are excluded from
#' the RT summary.
#'
#' @importFrom Hmisc wtd.var
#' @export
summary.LB4L <- function(data, DV = "accuracy", level = "subject",
                         given_practice = FALSE, given_final = (DV == "RT")) {

  level <- tolower(level)
  if ( length(level) > 1 || !(level %in% c("subject", "group", "raw")) ) {
    stop('level argument must be either "subject", "group", or "raw"')
  }

  if ( length(DV) > 1 || !(DV %in% c("RT", "accuracy"))) {
    stop('DV argument must be either "RT" or "accuracy"')
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

  if (!all(is.numeric(given_practice) || is.logical(given_practice),
           is.numeric(given_final) || is.logical(given_final))) {
    stop("'given_practice' and 'given_final' arguments must be logical scalars or numeric vectors")
  }

  data <- rename(data, final_acc = acc, final_RT = RT)

  if (!any(given_practice > 0) && given_final) {
    data <- DVfun(data, level = level,
                  conditional_vars = "final_acc")

  } else if (is.numeric(given_practice) && all(given_practice > 0)) {
    rounds <- 1:attr(data,"practice_tests")
    given_data <- get(attr(data,"tables")$test_practice)

    if (!setequal(given_practice, rounds)) {
      given_data <- filter(given_data, round %in% given_practice)
    }

    if (given_final) {
      conditional_vars <- c(paste0("practice", given_practice, "acc"), "final_acc")
    } else {
      conditional_vars <- paste0("practice", given_practice, "acc")
    }

    data <- DVfun(data, level = level,
                  given_data = given_data,
                  conditional_vars = conditional_vars)

  } else if (identical(given_practice, TRUE)) {
    given_practice <- 1:attr(data,"practice_tests")
    given_data <- get(attr(data,"tables")$test_practice)

    if (given_final) {
      conditional_vars <- c(paste0("practice", given_practice,"acc"), "final_acc")
    } else {
      conditional_vars <- paste0("practice", given_practice,"acc")
    }

    data <- DVfun(data, level = level,
                  given_data = given_data,
                  conditional_vars = conditional_vars)

  } else {
    data <- DVfun(data, level = level)
  }

  return(data)
}

# Internal function for de-normalizing the test practice data to make a column
# for each round of test practice
#' @importFrom tidyr spread_
reshaping <- function(given_data, conditional_vars, DV) {

  size <- unique(given_data$round)
  spread_names <- as.character(size)
  prefixed_names <- paste0("practice", spread_names, DV)

  if (length(size) == 1) {
    given_data <- rename_(given_data, .dots = setNames(DV, prefixed_names))

  } else if (length(size) > 1) {
    spread_names <- as.character(unique(given_data$round))
    prefixed_names <- paste0("practice", spread_names, DV)
    given_data <- given_data %>%
      select_("subject:sameCue", "round", .dots = DV) %>%
      spread_("round", DV) %>%
      rename_(.dots=setNames(sapply(spread_names,as.name),
                             prefixed_names))
  }

  return(given_data)
}


# Internal funcion for summarizing accuracy data from LB4L experiments
#' @importFrom tidyr expand_ replace_na
accuracy <- function(data, level, given_data = NULL, conditional_vars = NULL) {

  if (any(grepl("practice", conditional_vars))) {
    data <- ungroup(data) %>%
      select(subject:target, final_acc) %>%
      left_join(x = reshaping(given_data, conditional_vars, DV = "acc"),
                y = ., by = c("subject", "group", "list", "target")) %>%
      select(subject:list, pracCue = cue.x, finalCue = cue.y,
             target, sameCue, starts_with("practice"), final_acc)
  }

  if (level == "raw") {
    return(data)
  }

  unconditional_vars <- vapply(groups(data), as.character, character(1L))
  data <- group_by_(data, .dots = conditional_vars, add = TRUE)

  if (any(grepl("final", conditional_vars))) {

    fill_vars <- c(unconditional_vars[unconditional_vars != "group"],
                   conditional_vars)
    possible_n <- data %>%
      group_by_(.dots = unconditional_vars) %>%
      summarise(cell_max = n())

    summarized_data <- data %>%
      summarise(freq = n()) %>%
      left_join(possible_n) %>%
      mutate(probability = freq/cell_max) %>%
      select(-freq) %>%
      left_join(x = expand_(ungroup(.), fill_vars), y = .) %>%
      replace_na(list(probability = 0)) %>%
      group_by_(.dots = unconditional_vars[unconditional_vars != "group"]) %>%
      mutate_each(funs = funs(replace(., is.na(.), first(.[!is.na(.)]))),
                  group, cell_max)

    if (level == "group") {
      # grouping_vars <- vapply(groups(data), as.character, character(1L))
      possible_n <- distinct(select(ungroup(possible_n), -subject))
      summarized_data <- summarized_data %>%
        group_by_(.dots = c(unconditional_vars[unconditional_vars != "subject"],
                            conditional_vars)) %>%
        left_join(possible_n) %>%
        summarise(probability = mean(probability))
    }

    summarized_data <- as_LB4L_joint_summary(summarized_data)

  } else {
    summarized_data <- data %>%
      summarise(avg_final_acc = mean(final_acc),
                nObs = n())

    if (level == "group") {
      grouping_vars <- vapply(groups(data), as.character, character(1L))
      summarized_data <- summarized_data %>%
        group_by_(.dots = grouping_vars[grouping_vars != "subject"]) %>%
        summarise(sd_final_acc = sd(avg_final_acc),
                  w_sd_final_acc = sqrt(wtd.var(avg_final_acc, nObs)),
                  w_avg_final_acc = weighted.mean(avg_final_acc, nObs),
                  avg_final_acc = mean(avg_final_acc),
                  nObs = sum(nObs),
                  N = length(unique(subject))) %>%
        mutate(sem_final_acc = sd_final_acc/sqrt(N)) %>%
        select_(.dots = c(grouping_vars[grouping_vars != "subject"],
                          "avg_final_acc", "sd_final_acc", "sem_final_acc",
                          "nObs", "N", "w_avg_final_acc", "w_sd_final_acc"))
    }

    if (any(grepl("practice[0-9]+", conditional_vars))) {
      summarized_data <- as_LB4L_CD_summary(summarized_data)
    } else {
      summarized_data <- as_LB4L_IV_summary(summarized_data)
    }
  }

  return(summarized_data)
}


# Internal method for summarizing RT data from LB4L experiments
RT <-function(data, level, given_data = NULL, conditional_vars = NULL) {

  if (!is.null(given_data) && any(grepl("practice*", conditional_vars))) {
    given_data <- given_data %>% filter(is.finite(RT))
    data <- ungroup(data) %>%
      select(subject:target, final_acc, final_RT) %>%
      left_join(x = reshaping(given_data, conditional_vars, DV = "acc"),
               y = .,
               by = c("subject", "group", "list", "target")) %>%
      left_join(x = reshaping(given_data, conditional_vars, DV = "RT"),
                y = .,
                by = c("subject", "group", "list", "target", "sameCue")) %>%
      select(subject:list, pracCue = cue.x, finalCue = cue.y,
             target, sameCue, starts_with("practice"), final_acc, final_RT)
  }

  if (level == "raw") {
    return(data)
  }

  grouping_vars <- c(vapply(groups(data), as.character, character(1L)), conditional_vars)
  unconditional_vars <- vapply(groups(data), as.character, character(1L))
  data <- ungroup(data) %>%
    filter( !(Reduce('|', lapply(data[,grepl("RT", names(data))], is.na))) ) %>%
    group_by_(.dots = grouping_vars)

  summarized_data <- data %>%
    summarise_each(~median, contains("RT")) %>%
    left_join(summarise(data, nObs = n()))

  if (level == "group") {
    counts <- summarized_data %>%
      group_by_(.dots = grouping_vars[grouping_vars != "subject"]) %>%
      summarise(nObs = sum(nObs),
                N = length(unique(subject)))
    summarized_data <- summarized_data %>%
      group_by_(.dots = grouping_vars[grouping_vars != "subject"]) %>%
      summarise_each(funs = funs('median_RT' = median, mad(., constant = 1)), contains("RT")) %>%
      left_join(counts, by = grouping_vars[grouping_vars != "subject"])
  }

  if (!is.null(conditional_vars)) {
    summarized_data <- as_LB4L_CD_summary(summarized_data)
  } else {
    summarized_data <- as_LB4L_IV_summary(summarized_data)
  }

  return(summarized_data)
}

#' Plot Accuracy or RT for each experimental condition of the LB4L2 dataset
#'
#' @param data An LB4L2_IV_summary data frame from the FAM package.
#' @export
autoplot.LB4L_IV_summary <- function(data) {

  data$cond <- interaction(data$practice, data$OCpractice)
  if ("median_RT" %in% names(data)) {
    p <- ggplot(data = data,
                aes_string(x = "group", y = "median_RT", color = "cond", group = "cond")) +
      geom_point(size=3, shape = 2) +
      ylab("Median First-Press Latency")

  } else if ("avg_final_acc" %in% names(data)) {
    p <- ggplot(data = data,
                aes_string(x = "group", y = "avg_final_acc", color = "cond", group = "cond")) +
      geom_point(size=3) +
      geom_errorbar(aes(ymax = avg_final_acc + sem_final_acc,
                        ymin = avg_final_acc - sem_final_acc),
                    width = .025) +
      ylab("Accuracy")

  } else {
    stop("Did not recognize any DV columns in this LB4L_IV_summary object")
  }

  p <- p +  geom_line(size=1) +
    scale_color_discrete("Practice\nCondition",
                         breaks = c("N.N","N.S","N.T","S.N","T.N"),
                         labels = c(N.N = "Baseline", N.S ="Other Cue Study", S.N  = "Restudy",
                                    N.T = "Other Cue Test", T.N = "Test Same Cue")) +
    scale_x_discrete("Group",expand=c(0,.25),
                     limits = c("immediate", "delay"),
                     labels=c("Immediate","Delay")) +
    theme(legend.key.height = unit(2,"line")) +
    ggtitle('Final Test')

  return(p)

}

#' Plot Accuracy or RT for each experimental condition of the LB4L2 dataset
#'
#' @param data An LB4L2_IV_summary data frame from the FAM package.
#' @export
autoplot.LB4L_CD_summary <- function(data) {

  lookup_table <- c("0" = "Incorect", "1" = "Correct")
  labels_fun <- labeller(practice1acc = setNames(paste("Practice Test 1", lookup_table), names(lookup_table)),
                         practice2acc = setNames(paste("Practice Test 2", lookup_table), names(lookup_table)),
                         final_acc = setNames(paste("Final Test", lookup_table), names(lookup_table)))

  if ("median_RT" %in% names(data)) {
    data$cond <- interaction(data$practice, data$OCpractice)
    p <- ggplot(data = data,
                aes_string(x = "group", y = "median_RT", color = "cond", group = "cond")) +
      geom_point(size=3, shape = 2) +
      facet_grid(~final_acc, labeller = labels_fun) +
      ylab("Median First-Press Latency") +
      scale_color_discrete("Practice\nCondition",
                           breaks = c("N.N","N.S","N.T","S.N","T.N"),
                           labels = c(N.N = "Baseline", N.S ="Other Cue Study", S.N  = "Restudy",
                                      N.T = "Other Cue Test", T.N = "Test Same Cue"))

  } else if ("avg_final_acc" %in% names(data)) {
    givens <- grep("practice", names(data), value = TRUE)


    p <- ggplot(data = data,
                aes_string(x = "group", y = "avg_final_acc", color = "sameCue", group = "sameCue")) +
      geom_point(size=3) +
      facet_grid(~practice1acc, labeller = labels_fun) +
      geom_errorbar(aes(ymax = avg_final_acc + sem_final_acc,
                        ymin = avg_final_acc - sem_final_acc),
                    width = .025) +
      scale_color_discrete("Cue Used", breaks = c("no","yes"),
                           labels = c(no = "Other Cue", yes ="Same Cue")) +
      ylab("Condtional Accuracy")

  } else {
    stop("Did not recognize any DV columns in this LB4L_IV_summary object")
  }

  p <- p +  geom_line(size=1) +
    scale_x_discrete("Group",expand=c(0,.25),
                     limits = c("immediate", "delay"),
                     labels=c("Immediate","Delay")) +
    theme(legend.key.height = unit(2,"line")) +
    ggtitle('Final Test Accuracy given Practice Test Accuracy')

  return(p)

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
