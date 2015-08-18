#' @import dplyr
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
