# LB4L data munging
library(dplyr)
library(whoppeR)
library(tidyr)
library(ggplot2)
library(grid)

## Step 1: Find the datafiles ####
options(stringsAsFactors = FALSE)
# Have to escape the \ with an \ so the second \ can escape the . in .csv =)
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

## Start with study data and build a 'condition' index column on the study table
# since it has the full design for each each list
study <- data.frame(do.call(rbind,lapply(sFList, read.csv)), phase = "study") %>%
  group_by(subject) %>%
  mutate(trial = 1:n()) %>%
  ungroup() %>%
  replace_na(list(test = FALSE)) %>%
  mutate(pracFactor  = factor(replace(practice, practice=='', 'C'))) %>%
  group_by(subject, target) %>%
  mutate(OCpractice = rev(pracFactor),
         condition = interaction(pracFactor, OCpractice)) %>%
  select(-pracFactor)

# Now study practice data
sp <- data.frame(do.call(rbind,lapply(spFList, read.csv)), phase = "practice")

# Now test practice data
tp <- data.frame(do.call(rbind,lapply(tpFList, read.csv)), phase = "practice") %>%
  mutate(practice = 'T', test = replace(test, !is.finite(test), 0))
final <- data.frame(do.call(rbind,lapply(fFList, read.csv)), phase = "final")

## Step 1: Score the data ####

final <- final %>% group_by(subject, list, cue) %>%
  mutate(acc = score(target,response))
tp <- tp %>% group_by(subject, list, cue, round) %>%
  mutate(acc = score(target,response))

## Use the interactive editor to make decisions about fuzzy matches
# First, the final test data
final_cleaned_path  = file.path("data-raw","LB4L2_final_cleaned.Rdata")
# if the cleaned file exists, load it
if (file.exists(final_cleaned_path)) {
  load(final_cleaned_path)
  # If the number of rows don't match, we probably added more data and need to score it again
  # This check should really by smarter, i.e., make you clean only new rows in the data files
  # that are not in the clean file, or not rely on the row number heuristic
  # But it ought to save us from clobbering our data with the wrong stuff from disk
  if (sum(final$acc %in% 2) != nrow(final_cleaned)) {
    final_cleaned <- edit(final[final$acc %in% 2, c("subject","target", "response","acc")])
    save(final_cleaned, file = file.path("data-raw","LB4L2_final_cleaned.Rdata"))
  }
} else {
  final_cleaned <- edit(final[final$acc %in% 2, c("subject","target", "response","acc")])
  save(final_cleaned, file = file.path("data-raw","LB4L2_final_cleaned.Rdata"))
}
stopifnot(all(final_cleaned$acc %in% c(0,1)))
final[final$acc %in% 2, c("subject","target", "response","acc")] <- final_cleaned


# Next, the final test data
tp_cleaned_path  = file.path("data-raw","LB4L2_tp_cleaned.Rdata")
if (file.exists(tp_cleaned_path)) {
  load(tp_cleaned_path)
  # If the number of rows don't match, we probably added more data and need to score it again
  # This check should really by smarter, i.e., make you clean only new rows in the data files
  # that are not in the clean file, or not rely on the row number heuristic
  # But it ought to save us from clobbering our data with the wrong stuff from disk
  if (sum(tp$acc %in% 2) != nrow(tp_cleaned)) {
    tp_cleaned <- edit(tp[tp$acc %in% 2, c("subject","target", "response","acc")])
    save(tp_cleaned, file = file.path("data-raw","LB4L2_tp_cleaned.Rdata"))
  }
} else {
  tp_cleaned <- edit(tp[tp$acc %in% 2, c("subject","target", "response","acc")])
  save(tp_cleaned, file = file.path("data-raw","LB4L2_tp_cleaned.Rdata"))
}
stopifnot(all(tp_cleaned$acc %in% c(0,1)))
tp[tp$acc %in% 2, c("subject","target", "response","acc")] <- tp_cleaned


# Join the other-cue practice and condition columns from the study
# table onto the final test table
final <- study %>%
  select(subject:target, OCpractice, condition) %>%
  right_join(select(final,-test),
             by = c("subject", "group", "list", "cue","target")) %>%
  select(subject:target, practice, OCpractice, condition,onset:acc)


IVsummary <- final %>%
  group_by(group, practice, OCpractice, condition) %>%
  summarise(avgAcc = mean(acc), sdAcc = sd(acc), groupN = length(unique(subject))) %>%
  mutate(sem = sdAcc/sqrt(groupN))

ggplot(data =IVsummary, aes(x = group, y = avgAcc, color = condition,
                            group = condition)) +
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

# Join the test practice and final test data for test practiced items
# for conditional accuracy analysis

tested <- tp %>%
  replace_na(list(test = FALSE)) %>%
  select(subject:target, test, round, acc) %>%
  spread(round, acc) %>%
  rename(round1 = `1`, round2= `2`) %>%
  left_join(select(final, subject:target, acc),
            by = c("subject", "group", "list", "target")) %>%
  select(subject:list, pracCue = cue.x, finalCue = cue.y, target:round2, finalAcc=acc) %>%
  mutate(test = factor(test, labels = c('other','same'), exclude = NULL))

F_given_R1acc <- tested %>%
  group_by(group, test, round1) %>%
  summarise(condAcc = mean(finalAcc), sdAcc = sd(finalAcc), groupN = length(unique(subject))) %>%
              mutate(sem = sdAcc/sqrt(groupN))
F_given_R2acc <- tested %>%
  group_by(group, test, round2) %>%
  summarise(condAcc = mean(finalAcc), sdAcc = sd(finalAcc), groupN = length(unique(subject))) %>%
  mutate(sem = sdAcc/sqrt(groupN))
R2_given_R1 <- tested %>%
  group_by(group, test, round1) %>%
  summarise(condAcc = mean(round2), sdAcc = sd(round2), groupN = length(unique(subject))) %>%
  mutate(sem = sdAcc/sqrt(groupN)) %>%
  ungroup() %>%
  complete(group, test, round1)
R1_R2_joint  = tested %>%
  group_by(group, test, round1, round2) %>%
  summarise(n = n()) %>%
  group_by(group, test) %>%
  mutate(jointAcc = n/sum(n))
F_given_R1R2_joint <- tested %>%
  group_by(group, test, round1, round2) %>%
  summarise(condAcc = mean(finalAcc), sdAcc = sd(finalAcc), groupN = length(unique(subject))) %>%
  mutate(sem = sdAcc/sqrt(groupN)) %>%
  ungroup() %>%
  complete(group, test, round1)

ggplot(F_given_R1R2_joint,
       aes(x = interaction(round1,round2), y = condAcc)) +
  geom_point(size = 3) +
  facet_grid(test ~ group) +
  scale_x_discrete("Joint Practice Accuracy")


## RT Analysis ####
tp_RT <- tp %>%
  filter(is.finite(firstPress), is.finite(lastPress)) %>%
  mutate(RT = firstPress-onset) %>%
  group_by(group, test) %>%
  summarise(medianRT = mean(RT), MAD = mad(RT, constant = 1))

final_RT <- final %>%
  filter(is.finite(firstPress), is.finite(lastPress)) %>%
  mutate(RT = firstPress-onset) %>%
  group_by(group, condition) %>%
  summarise(medianRT = median(RT), MAD = mad(RT, constant = 1))


ggplot(data = final_RT,
       aes(x = group, y = medianRT, color = condition, group = condition)) +
  geom_point(size=3, shape = 2) +
  geom_line(size=1) +
  scale_color_discrete("Practice\nCondition",
                       breaks = c("C.C","C.S","C.T","S.C","T.C"),
                       labels = c(C.C = "Baseline", C.S ="Other Cue Study", S.C  = "Restudy",
                                  C.T = "Other Cue Test", T.C = "Test Same Cue")) +
  scale_x_discrete("Group",expand=c(0,.25), limits = c("immediate", "delay"),
                   labels=c("Immediate","Delay")) +
  ylab("Median First-Press Latency") +
  theme(legend.key.height = unit(2,"line")) +
  ggtitle('Final Test')
