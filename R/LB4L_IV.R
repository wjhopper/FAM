#' @import dplyr
#' @export
LB4L_IV <- function(data) {

  cols <- colnames(data) %in% c("subject","list","practice","repeating","target")
  data[cols] <- lapply(data[cols], factor)
  data$group <- factor(data$group, levels = c("immediate", "delay"))
  data$other_type <- factor(data$other_type,exclude = NULL)

  conds_by_ss <- data %>%   filter(list != 1) %>%
    group_by(subject, group, practice,other_type)  %>%
    summarize(prac_acc = mean(prac_score,na.rm=TRUE),
              final_acc=mean(final_score),
              n=n()) %>%
    group_by(subject) %>%
    mutate(chain=factor(1:n()))

  conds_collapsed <- data %>%  filter(list != 1) %>%
    group_by(group, practice,other_type)  %>%
    summarize(prac_acc = mean(prac_score,na.rm=TRUE),
              final_acc=mean(final_score),
              n=n()) %>%
    group_by(group) %>%
    mutate(chain=factor(1:n()))

  return(list(subject= conds_by_ss, groups = conds_collapsed))
}
