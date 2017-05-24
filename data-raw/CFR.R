library(dplyr)
library(whoppeR)

encodeClass <- function(x) {
  if (x[1] == 'S') {
    return('sp')
  } else if (x[1] == 'T')  {
    return('tp')
  } else {
    return('np')}
}

columnNames <- c("subject","list","practice", "target","abs_order","final_order",
                 "onset","FP","LP", "resp","score")
# Find all the files from the 3 phases
final_flist <- list.files("data-raw", pattern="CFR+.+Final_Data.csv",
                          full.names = TRUE)
study_flist <- list.files("data-raw", pattern="CFR+.+Study_Data.csv",
                          full.names = TRUE)
practice_flist <- list.files("data-raw", pattern="CFR+.+Practice_Data.csv",
                             full.names = TRUE)

# Create data frame of Final Phase Data
final_data <- data.frame(do.call(rbind,
                                 lapply(final_flist, read.csv, header=T, sep=',',
                                        col.names=columnNames)),
                         phase="final") %>%
  group_by(subject,list) %>%
  mutate(class = encodeClass(practice)) %>%
  ungroup() %>%
  select(-final_order) %>%
  mutate(subject = replace(subject, subject==55, 20))

# Create data frame of Practice Phase Data
prac_data <- data.frame(do.call(rbind,
                                lapply(practice_flist, read.csv, header=T, sep=',',
                                       col.names=columnNames)),
                        phase= "prac",
                        class= 'np') %>%
  select(-final_order) %>%
  mutate(subject = replace(subject, subject==55, 20))

# Create study phase data frame
study_data <- data.frame(do.call(rbind,
                                lapply(study_flist, read.csv, header=T, sep=',')),
                         FP = NA,
                         LP=NA,
                         resp = NA,
                         score= NA,
                         phase= "study",
                         class = 'np') %>%
  select(subject = sub_num, list,practice,target,abs_order = study_order,
         onset=study_onset, FP, LP, resp, score, phase, class) %>%
  mutate(subject = replace(subject, subject==55, 20))

test_data <- bind_rows(filter(prac_data, practice != 'S', list != 0),
                       filter(final_data, list != 0)
                       )%>%
  arrange(subject, list, phase) %>%
  mutate(score = NA_real_,
         onset = replace(onset, onset == 0, NA_real_),
         FP = replace(FP, !is.finite(FP) | FP == 0, NA),
         LP = replace(LP, !is.finite(LP) | LP == 0, NA),
         resp = replace(resp, resp=="", NA_character_),
         resp_dur = LP-FP) %>%
  group_by(subject, list, phase) %>%
  mutate(practice = replace(practice, practice == "", practice[1])) %>%
  ungroup()

target_stimuli <- select(.data = test_data,
                         subject, list, class, practice, phase, target) %>%
  group_by(subject, list, phase) %>%
  mutate(order = 1:n()) %>%
  ungroup()

test_responses <- filter(.data = test_data, !is.na(resp)) %>%
  select(-target, -abs_order) %>%
  group_by(subject, list, phase) %>%
  mutate(order = 1:n()) %>%
  ungroup()

test_data <- full_join(target_stimuli, test_responses,
                       by = c("subject", "list", "practice", "class", "phase", "order")
                       ) %>%
  group_by(subject, list, phase) %>%
  mutate(RT = c(FP[1]-onset[1], diff(FP)-resp_dur[1:n()-1]),
         CRT = cumsum(RT),
         score = if (phase[1]=='final'| (phase[1]=='prac' & practice[1]=='T')) {
                   whoppeR::score(target,resp,checkDuplicates = TRUE,
                         ignoreVals=NA_real_)
                 } else {
                   NA_real_
                 },
         where = if (phase[1]=='final'| (phase[1]=='prac' & practice[1]=='T')) {
           score(target,resp,checkDuplicates = TRUE,
                 ignoreVals=NA,index.return = TRUE)$position
         } else {
           NA_real_
         }) %>%
  ungroup() %>%
  arrange(subject, list, desc(phase), order)

if (file.exists("data-raw/CFRcleaned.Rdata")) {
  load(file="data-raw/CFRcleaned.Rdata")
  CFRcleaned <- filter(CFRcleaned, list != 0)
  cleaned_items = nrow(CFRcleaned)
} else {
  cleaned_items = NA
}

is_fuzzy_match <- test_data$score %in% 2

if (is.na(cleaned_items) || cleaned_items != sum(is_fuzzy_match)) {
  CFRcleaned <- filter(test_data3, any(score == 2)) %>%
    select(subject,list, phase, target, resp, score, where) %>%
    mutate(target = target[where]) %>%
    filter(score %in% 2) %>%
    arrange(subject, list, phase) %>%
    edit()

  valid_answer <- FALSE
  while (!valid_answer) {
    answer <- pmatch(tolower(readline("Save as CFRcleaned.Rdata (yes/no)? ")),
                     c("yes", "no"),
                     nomatch = 0
                     )
    if (answer == 1) {
      save(CFRcleaned, file = "data-raw/CFRcleaned.Rdata")
    }

    if (answer %in% c(1, 2)) {
      valid_answer <- TRUE
    }

  }

}
# Note: I know from visual inspection that 'tale' (row 24) is correct
#   it matches 'tail' not cottage, set it to 1, and its position to 1
# Note: I know from visual inspection that "flooe" (row 45) is a duplicate
#   of "floor" in the same list", set it to 0
# Note: I know from visual inspection that "canon" (row 73) is a duplicate
#   of "cannon" in the same list", set it to 0
# Note: I know "staek" in the test_data data frame (row 1358) is a correct
#   answer, dunno how it got missed in fuzzy matching
# Note: I know from visual inspection that "bridge" (row 82) should not
# be considered a match to "ridge", as ridge is recalled later.

test_data$score[test_data$resp %in% 'staek']<-1
test_data$score[test_data$resp %in% 'flooe' &
                  test_data$subject == 13 &
                  test_data$order == 13 &
                  test_data$list == 3 &
                  test_data$phase == "final"] <- 2

# Replace the fuzzy matches with the manual decisions
test_data[is_fuzzy_match, c("score", "where")] <- NA
test_data <- left_join(test_data,
                       select(CFRcleaned, -target),
                       by = c("subject", "list", "phase", "resp")) %>%
  mutate(score = coalesce(score.x, score.y),
         where = coalesce(where.x, where.y)) %>%
  # filter(subject == 13, list == 3, phase=="final") %>%
  group_by(subject, list, phase) %>%
  mutate(resp = replace(resp,
                        is.na(score.x) & score == 1,
                        target[where[is.na(score.x) & score == 1]])
         ) %>%
  select(-ends_with(".y"), -ends_with(".x")) %>%
  mutate(repeats = as.numeric(duplicated(resp,incomparables = c(NA_character_,''))),
         score = replace(score, repeats == 1, 0),
         where = replace(where, score == 0, NA),
         intrusions  = as.numeric(score==0 & repeats ==0 & resp != ''),
         target = replace(target, target == "", NA_character_)) %>%
  group_by(subject,phase, practice) %>%
  mutate(cond_list =  rep(1:length(unique(list)),as.vector(table(list)))) %>%
  ungroup() %>%
  arrange(subject, list, desc(phase))

duplicate_test <- group_by(test_data, subject, class, cond_list) %>%
  mutate(dup = duplicated(where, incomparables = NA_real_)) %>%
  filter(where %in% where[dup]) %>%
  arrange(subject, list, phase, where)

# CFR <- test_data
# save(CFR, file="data/CFR.rda")

