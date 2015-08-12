# LB4L data munging
library(dplyr)

score <- function(target,resp,s) {
  for (k in which(grepl('.', resp))) {
    # Try hard string matching first
    acc <- resp[k] == target[k]
    if (acc) {
      s[k] <- 1 # response is exact match to target string
    } else { 
      # try fuzzy matching 
      inds <- agrepl(resp[k], target[k])
      if (inds) {
        s[k] <- 2 # A fuzzy substring match was found
      } else  {
        s[k] <- 0 # hard and fuzzy matching failed
        
      }
    }
  }

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
tmp <- data %>% filter(list !=1) %>% 
  group_by(subject,group,list) %>% 
  mutate(prac_score = score(target,prac_resp,prac_score),
         final_score = score(target,final_resp,final_score))

# Use the interactive editor to make decisions about fuzzy matches 
cleaned1 <- edit(tmp[tmp$prac_score %in% 2, c("subject","target","prac_resp","prac_score")])
cleaned2 <- edit(tmp[tmp$final_score %in% 2, c("subject","target","final_resp","final_score")])

# Replace the fuzzy matches with the manual decisions
tmp$prac_score[tmp$prac_score %in% 2] <- cleaned1$prac_score
tmp$final_score[tmp$final_score %in% 2] <- cleaned2$final_score

tmp$final_score[is.nan(tmp$final_score)] <- 0
tmp$prac_score[is.nan(tmp$prac_score) & tmp$practice =='T'] <- 0
data <- rbind(filter(rawData, list==1), tmp)
save(data, file.path("data","LB4L.Rdata"))

