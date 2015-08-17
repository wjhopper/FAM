#' @export
badSubs <- function(data) {

  sub_means <- data %>% filter(list != 1) %>%
    group_by(subject,group) %>%
    summarise(final_acc=mean(final_score, na.rm=T),
              prac_acc = mean(prac_score, na.rm=T))

  # Find the terrible performers  and remove
  bad <- sub_means$subject[sub_means$final_acc < .05 | sub_means$prac_acc <.05]
  sub_means <- sub_means[!(sub_means$subject %in% bad),] %>% droplevels()

  return(list(removed = bad, means = sub_means))
}
