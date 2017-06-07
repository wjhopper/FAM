#' @import dplyr
#' @import ggplot2
NULL

#' @export
LB4L_joint <- function(data) {

  columns <- c("prac_score", "other_prac_acc")
  data[,columns] <- lapply(data[,columns],factor,exclude=NULL)

  jointAcc <- data  %>% filter(list != 1)  %>%
    group_by(subject, group, practice, other_type,
             prac_score,other_prac_acc,final_score) %>%
    summarise(cell_count = n())

  # I know from visual inspection that subject 20 has complete cases in each cell
  # So, use subject 20's data as a template to join by
  joinerFrame <- ungroup(filter(jointAcc,subject==20)) %>%
    select(practice,other_type,prac_score,other_prac_acc,final_score)
  jointAcc <- left_join(cbind(subject = rep(unique(data$subject),
                                            each=nrow(joinerFrame)),
                              joinerFrame),
                        jointAcc) %>%
    group_by(subject) %>%
    mutate(group = group[!is.na(group)][1],
           cell_count = replace(cell_count, is.na(cell_count),0)) %>%
    mutate(acc = cell_count/15,
           merged_prac_score = ifelse(is.na(levels(prac_score)[prac_score]),
                                      as.numeric(levels(other_prac_acc))[other_prac_acc],
                                      as.numeric(levels(prac_score))[prac_score]))
  jointAcc$merged_prac_score <- factor(jointAcc$merged_prac_score)

  jointAcc_grouped <- jointAcc %>%
    group_by(group, practice, other_type, prac_score,
             other_prac_acc,final_score) %>%
    summarise(avgAcc = mean(acc))

  return(list(subject = jointAcc, groups = jointAcc_grouped))
}

#   left_join(select(ungroup(conds_by_ss),
#                    subject,group,practice,other_type,grp_size=n),
#             by =c("subject","group","practice","other_type")) %>%
#    mutate(acc = cell_count/grp_size)

#   bigJoinerFrame <- cbind(subject = rep(unique(data$subject),
#                                         each=nrow(joinerFrame)),
#                           joinerFrame)

#' @export
LB4L_conditional <- function(data) {
  columns <- c("prac_score", "other_prac_acc")
  data[,columns] <- lapply(data[,columns],factor,exclude=NULL)

  conAcc <- data  %>% filter(list != 1)  %>%
    group_by(subject, group, practice, other_type, prac_score,other_prac_acc) %>%
    summarise(count = n(), # of observations contributing to mean
              final_acc= mean(final_score),# mean of conditional group for each subject, dropping missings
              final_n = final_acc*count)

  # To deal with missing cells dues to practice performance, join the grouped
  # data frame with another data frame that has all combinations of factor levels that could go missing
  # including subject number is important
  # the joiner frame was a pain to build by hand but I don't know a better way!
  nsubs <- length(levels(conAcc$subject))
  conAcc <- left_join(data.frame(subject = rep(factor(levels(conAcc$subject)),each=8),
                                 practice=rep(factor(c("C","C","T","T","T","T","C","S")),nsubs),
                                 other_type=rep(factor(c("T","T","C","C",NA,NA,NA,NA),exclude=NULL),nsubs),
                                 prac_score=rep(factor(c(NA,NA,0,1,0,1,NA,NA),exclude=NULL),nsubs),
                                 other_prac_acc=rep(factor(c(0,1,NA,NA,NA,NA,NA,NA),exclude = NULL),nsubs)),
                      conAcc) %>%
    group_by(subject) %>%
    mutate(count = replace(count, is.na(count), 0),
           final_acc = replace(final_acc,
                               is.na(final_acc) & count != 0 &
                                 !(practice=='T' & other_type=='C'),
                               0),
           # Correct the group identifier
           group = ifelse(is.na(group), group[!is.na(group)][1],group),
           merged_prac_score = ifelse(is.na(levels(prac_score)[prac_score]),
                                      as.numeric(levels(other_prac_acc))[other_prac_acc],
                                      as.numeric(levels(prac_score))[prac_score]))
  conAcc$merged_prac_score <- factor(conAcc$merged_prac_score)
  conAcc$group <- factor(conAcc$group)


  conAcc_grouped <- conAcc %>%
    group_by(group,practice, other_type, prac_score,other_prac_acc,merged_prac_score) %>%
    summarise(weighted_final_acc=weighted.mean(final_acc, count),
              final_acc= mean(final_acc,na.rm=TRUE),
              n_obs = sum(count),
              missing=length(count[count==0])/n(), #length(count),
              upper = qbinom(sqrt(.025), n_obs,final_acc,lower.tail=F)/n_obs,
              lower = qbinom(sqrt(.025), n_obs,final_acc,lower.tail=T)/n_obs)
  #             final_acc2=weighted.mean(final_acc[count != 0], w=count[count != 0]),
  #             final_acc3 = sum(final_n[count != 0])/sum(count[count != 0]))

  return(list(subject= conAcc, groups= conAcc_grouped))

}
