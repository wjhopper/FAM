library(dplyr)
library(whoppeR)

flist <- list.files("./data-raw", pattern="^SCKT+.+csv$", full.names = TRUE)
rawData <-lapply(flist, read.csv, header=T, sep=',',
                 col.names=c('subject','group','version','dataetime','list',
                             'half','repeating','practice','condition','cue',
                             'target', 'prac_order','prac_resp','prac_score',
                             'final_order','final_resp','final_score'))

# rbind to data frame and group
rawData <- do.call(rbind,rawData) %>%
  group_by(subject,list,cue) %>%
  # Compute accuracy on cued-recall practice test
  mutate(prac_score = score(target,prac_resp),
         final_score = score(target,final_resp))

# There never was a practice test for these so fix them to NA
rawData$prac_score[rawData$practice %in% c("C","S")]  <- NA

SCKTprac_clean <- edit(rawData[rawData$prac_score %in% 2, c("subject","target",
                                                      "prac_resp","prac_score")])
SCKTfinal_clean <- edit(rawData[rawData$final_score %in% 2, c("subject","target",
                                                       "final_resp","final_score")])
rawData$prac_score[rawData$prac_score %in% 2] <- SCKTprac_clean$prac_score
rawData$final_score[rawData$final_score %in% 2] <- SCKTfinal_clean$final_score

rawData <- rawData %>% group_by(subject,target) %>%
  # Add identifier and score for the other half of the repeating target 'pair'
  mutate(other_prac = if (n() > 1) { rev(prac_score) } else { NA_real_},
         other_type =  if (n() > 1) { rev(practice) } else { NA_character_},
         other_final = if (n() > 1) { rev(final_score) } else { NA_real_})

cols<- colnames(rawData) %in% c("subject","group", "half", "list", "practice",
                             "repeating","condition", "cond", "target",
                             "other_prac","other_final","other_type")
rawData[cols] <- lapply(rawData[cols], factor)

# Save as binary rda
# SCKT_allSs <- ungroup(rawData)
# save(SCKT_allSs, file = file.path("data","SCKT_allSs.rda"))
