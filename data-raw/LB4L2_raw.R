# LB4L2 data munging
options(stringsAsFactors = FALSE)
library(dplyr)
library(whoppeR)
library(tidyr)
source(file.path("R", "LB4L2.R"))

manual_scoring <- function(data) {

  scoring_loop <- function(data_to_score) {
    has_been_scored <- data_to_score$acc %in% c(0,1)
    while (!all(has_been_scored)) {
      scored <- edit(data_to_score[!has_been_scored, ])
      has_been_scored <- scored$acc %in% c(0,1)
    }
  }

  prefix <- as.character(as.list(match.call())$data)
  clean_path <- file.path("data-raw", paste("LB4L2", prefix, "cleaned.Rdata", sep = "_"))
  scoring_cols <- c("subject","target", "response")

  # if the cleaned file exists, load it
  if (file.exists(clean_path)) {
    f <- load(clean_path)
    clean_data <- get(f)
    rm(list = f)
    clean_data$subject <- as.integer(clean_data$subject)

    if (!identical(as.data.frame(data[!data$acc %in% c(0,1), scoring_cols]),
                   clean_data[, scoring_cols])) {
      to_score <- anti_join(data, clean_data, by = scoring_cols)
      clean_data <- rbind(clean_data,
                      scoring_loop(to_score[, c(scoring_cols,"acc")]))
      assign(f, clean_data)
      save(list = f, file = clean_path)
    }

  } else {
    clean_data <- scoring_loop(data[!data$acc %in% c(0,1), c(scoring_cols,"acc")])
    assign(f, clean_data)
    save(list = f, file = clean_path)
  }

  data[data$acc %in% 2, c(scoring_cols,"acc")] <- clean_data
  stopifnot(all(data$acc %in% c(0,1)))
  return(data)
}

## Step 1: Find the datafiles ####

# Have to escape the \ with an \ so the second \ can escape the . in .csv =)
# yo dawg I head you like \ ...
sFList <- list.files("data-raw", pattern="LB4L2+.+Study\\.csv", full.names = TRUE)
spFList <- list.files("data-raw", pattern="LB4L2+.+StudyPractice\\.csv", full.names = TRUE)
tpFList <- list.files("data-raw", pattern="LB4L2+.+TestPractice\\.csv", full.names = TRUE)
fFList<- list.files("data-raw", pattern="LB4L2+.+Final\\.csv", full.names = TRUE)

# Check that no file is in multiple lists, if so the regex failed to partition them
overlap <- Map(intersect,
               list(sFList,  sFList,  sFList, spFList, spFList, tpFList),
               list(spFList, tpFList, fFList, tpFList, fFList,  fFList))
stopifnot(all(vapply(overlap, identical, logical(1), character(0))))

# Check that each group of files has the same number of subjects
# If not we're missing files or matching the wrong set of files somewhere
stopifnot(all(vapply(list(spFList,tpFList,fFList),
                     function(x) length(x) == length(sFList),
                     logical(1))))

## Step 2: Read the data files in each group into their own table #####

## Start with study data
study <- data.frame(do.call(rbind,lapply(sFList, read.csv)), phase = "study") %>%
  replace_na(list(test = FALSE)) %>%
  mutate(practice = replace(practice, practice == '', 'N')) %>%
  group_by(subject) %>%
  mutate(trial = 1:n()) %>%
  group_by(subject, target) %>%
  mutate(OCpractice = rev(practice)) %>%
  ungroup()

# Now study practice data
sp <- data.frame(do.call(rbind,lapply(spFList, read.csv)), phase = "practice")

# Now test practice data
tp <- data.frame(do.call(rbind,lapply(tpFList, read.csv)), phase = "practice") %>%
  mutate(practice = 'T', test = replace(test, !is.finite(test), 0),
         RT = firstPress-onset) %>%
  group_by(subject, list, cue, round) %>%
  mutate(acc = score(target,response))
tp <- ungroup(manual_scoring(tp))

# Now final test data
final <- data.frame(do.call(rbind,lapply(fFList, read.csv)), phase = "final") %>%
  mutate(practice = replace(practice, practice=='', 'N'),
         RT = firstPress - onset) %>%
  group_by(subject, list, cue) %>%
  mutate(acc = score(target,response))
final <- ungroup(manual_scoring(final))

# Join the other-cue practice and condition columns from the study
# table onto the final test table
final  <- study %>%
  select(subject:target, OCpractice) %>%
  right_join(select(final,-test),
             by = c("subject", "group", "list", "cue","target")) %>%
  mutate(RT = firstPress-onset) %>%
  select(subject:target, practice, OCpractice, response, onset, firstPress, lastPress, acc, RT)


overwrite_flag <- FALSE

LB4L2_final <- as.LB4L(final)
attr(LB4L2_final,"experiment") <- "LB4L2"
devtools::use_data(LB4L2_final,overwrite = overwrite_flag)

LB4L2_study <- as.LB4L(study)
attr(LB4L2_study,"experiment") <- "LB4L2"
devtools::use_data(LB4L2_study,overwrite = overwrite_flag)

LB4L2_tp <- as.LB4L(tp)
attr(LB4L2_tp ,"experiment") <- "LB4L2"
devtools::use_data(LB4L2_tp,overwrite = overwrite_flag)

LB4L2_sp <- as.LB4L(sp)
attr(LB4L2_sp ,"experiment") <- "LB4L2"
devtools::use_data(LB4L2_sp,overwrite = overwrite_flag)




