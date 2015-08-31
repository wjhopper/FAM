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

#' @export
calc_odds <- function(final_acc,prac_score) {
  # if we take the log of (0/1), we get -Inf
  # if we take the log of (1/0), we get Inf
  # if we take the log of (0/0), we get NaN
  #   if ( length(prac_score[prac_score == 1]) != length(prac_score[prac_score == 0]) ) {
  #     stop("invalid input vector, lengths of groups are not equal)

  #   }
  if (!all(prac_score %in% 0:1)) {
    stop("Not all values in practice score vector are zero or one")
  }
  odds <- log(final_acc[prac_score == 1]/final_acc[prac_score == 0])
  return(odds )

}
