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
  mutate(class = encodeClass(practice))

# Create data frame of Practice Phase Data
prac_data <- data.frame(do.call(rbind,
                                lapply(practice_flist, read.csv, header=T, sep=',',
                                       col.names=columnNames)),
                        phase= "prac",
                        class= 'np')

# Create study phase data frame
study_data <- data.frame(do.call(rbind,
                                lapply(study_flist, read.csv, header=T, sep=',')),
                         final_order=NA,
                         FP = NA,
                         LP=NA,
                         resp = NA,
                         score= NA,
                         phase= "study",
                         class = 'np') %>%
  select(subject = sub_num, list,practice,target,abs_order = study_order,
         final_order, onset=study_onset, FP, LP, resp, score, phase, class)

rawData <- rbind(study_data, prac_data,final_data)
rawData$subject[rawData$subject==55] <- 20
rawData <- arrange(rawData, subject) %>%
  group_by(subject,list,phase) %>%
  mutate(practice = practice[1],
         order = 1:n()) %>%
  ungroup() %>%
  mutate(onset = replace(onset, !is.finite(onset) | onset ==0, NA),
         FP = replace(FP, !is.finite(FP) | FP ==0, NA),
         LP = replace(LP, !is.finite(LP) | LP ==0, NA),
         resp_dur = LP-FP) %>%
  group_by(subject,list,phase) %>%
  mutate(RT = c(FP[1]-onset[1],diff(FP)-resp_dur[1:n()-1]),
         CRT = cumsum(RT),
         score = if (phase[1]=='study'| (phase[1]=='prac' & practice[1]=='S')) {
                   return(NA_real_)
                 } else {
                   as.numeric(score(target,resp,checkDuplicates = TRUE,
                         ignoreVals=NA))
                 },
         where = if (phase[1]=='study'| (phase[1]=='prac' & practice[1]=='S')) {
           return(NA_real_)
         } else {
           score(target,resp,checkDuplicates = TRUE,
                 ignoreVals=NA,index.return = TRUE)$position
         })
CFRcleaned <- filter(rawData,any(score == 2)) %>%
  select(subject,target,resp, score, where) %>%
  mutate(target = target[where]) %>%
  filter(score %in% 2) %>%
  edit()
# Note: I know from visual inspection that 'tale' (row 24) is correct
#   it matches 'tail' not cottage, set it to 1
# Note: I know from visual inspection that "flooe" (row 45) is a duplicate
#   of "floor" in the same list", set it to 0
# Note: I know from visual inspection that "canon" (row 73) is a duplicate
#   of "cannon" in the same list", set it to 0


# Replace the fuzzy matches with the manual decisions
rawData$score[rawData$score %in% 2] <- CFRcleaned$score

rawData <- rawData %>%
  group_by(subject,list,phase) %>%
  mutate(repeats = if (phase[1]=='study'| (phase[1]=='prac' & practice[1]=='S')) {
           return(NA_real_)
         } else {
           as.numeric(duplicated(resp,incomparables = NA))
         },
         intrusions  = as.numeric(score==0 & repeats ==0 & response != '')) %>%
  group_by(subject,phase,practice) %>%
  mutate(cond_list =  rep(1:length(unique(list)),as.vector(table(list))))

#CFR_allSs <- rawData
# save(CFR_allSs, file="data/CFR_allSs.rda")

