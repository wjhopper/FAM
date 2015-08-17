#' @import dplyr
#' @import ggplot2
#' @export
#'
tailoringPhase <- function(data, return.data = TRUE, return.plot = TRUE) {

  tailoring_means <- data %>% filter(list == 1) %>%
    group_by(subject) %>%
    summarise(acc=mean(final_score, na.rm=T)) %>%
    group_by(bin = as.factor(cut(acc,c(0,.28, .76,1),
                                 labels =  c("< 25%", "25-75%","75% <"),
                                 include.lowest = TRUE)))
  tm <- as.data.frame(xtabs( ~ bin,tailoring_means)) %>%
    mutate(Freq = Freq/sum(Freq))

  tailor_plot <- ggplot(tm, aes(x=bin,y=Freq)) +
    geom_bar(stat='identity') +
    scale_y_continuous("Proportion of Subjects", limits = c(0,max(tm$Freq+.025))) +
    geom_text(aes(label =c("6 Secs","5 secs","4 Secs"), y= Freq+.025)) +
    xlab("Practice List Performance") +
    theme_larger() +
    ggtitle('Percent of Subjects with Specific Performance Levels')

  if (return.plot & return.data) {
    return(list(data = tm, figure = tailor_plot))
  } else if (return.plot & !return.data) {
    return(tailor_plot)
  } else if  (!return.plot & return.data) {
    return(tm)
  } else {
    return(invisible())
  }
}

