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
rawData$prac_score[rawData$pratice %in% c("C","S")]  <- NA

cleaned1 <- edit(rawData[rawData$prac_score %in% 2, c("subject","target",
                                                      "prac_resp","prac_score")])
cleaned2 <- edit(rawData[rawData$final_score %in% 2, c("subject","target",
                                                       "final_resp","final_score")])
rawData$prac_score[rawData$prac_score %in% 2] <- cleaned1$prac_score
rawData$final_score[rawData$final_score %in% 2] <- cleaned2$final_score

rawData <- rawData %>% group_by(subject,target) %>%
  # Add identifier and score for the other half of the repeating target 'pair'
  mutate(other_prac = if (n() > 1) { rev(prac_score) } else { NA_real_},
         other_type =  if (n() > 1) { rev(practice) } else { NA_character_},
         other_final = if (n() > 1) { rev(final_score) } else { NA_real_})


# Save as binary rda
# save(rawData, file = file.path("data","SCKT_allSs.rda"))
