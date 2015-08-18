# LB4L data munging
library(dplyr)
library(whoppeR)

flist <- list.files("data-raw", pattern="LB4L+.+csv", full.names = TRUE)
rawData <-lapply(flist, read.csv, header=T, sep=',',
           col.names=c('subject','group','version','datetime','exp_onset',
                       'list','repeating','practice','cue','target',
                       'abs_order','study_onset', 'abs_prac_order', 'prac_onset',
                       'prac_FP','prac_LP','prac_resp','prac_score',
                       'final_order', 'final_onset','final_FP','final_LP',
                       'final_resp','final_score'))
rawData <- do.call(rbind,rawData)

# Compute accuracy on cued-recall practice test
rawData <- rawData %>%
  group_by(subject,group,list,cue) %>%
  mutate(prac_score = score(target,prac_resp),
         final_score = score(target,final_resp))

# There never was a test for these so fix them at NA
rawData$prac_score[rawData$practice %in% c("C","S") | rawData$list ==1]  <- NA

# Use the interactive editor to make decisions about fuzzy matches
cleaned1 <- edit(rawData[rawData$prac_score %in% 2, c("subject","target",
                                              "prac_resp","prac_score")])
cleaned2 <- edit(rawData[rawData$final_score %in% 2, c("subject","target",
                                               "final_resp","final_score")])

# Replace the fuzzy matches with the manual decisions
rawData$prac_score[rawData$prac_score %in% 2] <- cleaned1$prac_score
rawData$final_score[rawData$final_score %in% 2] <- cleaned2$final_score

rawData <- rawData %>%group_by(subject,target) %>%
  # Add identifier and score for the other half of the repeating target 'pair'
  mutate(other_prac_acc  = if (n() > 1) { rev(prac_score) } else { NA_real_},
         other_type  =  if (n() > 1) { rev(practice) } else { NA_character_})

rawData$final_score[rawData$practice %in% "T" & rawData$other_type %in% "C"]  <- NA


# set data types so we know what to expect elsewhere
cols <- colnames(LB4L_allSs) %in% c("subject","list","practice","repeating","target")
LB4L_allSs[cols] <- lapply(LB4L_allSs[cols], factor)
LB4L_allSs$group <- factor(LB4L_allSs$group, levels = c("immediate", "delay"))
LB4L_allSs$other_type <- factor(LB4L_allSs$other_type,exclude = NULL)

# Save as binary rda
# save(rawData, file = file.path("data","LB4L_allSs.rda"))

