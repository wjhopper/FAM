#' @importFrom tidyr complete

#' @export
findSubjectClusters <- function(data, cut) {

  condsWide <- data %>% group_by(subject,half,practice,other_type) %>%
    summarize(final_score=sum(final_score)) %>%
    dcast(subject~half*practice*other_type,value.var = 'final_score')
  row.names(condsWide) <- condsWide$subject
  subjectDist <- dist(select(condsWide,-subject))
  tree <- hclust(subjectDist)
  return(list(tree = tree, cluster = data.frame(subject=unique(data$subject),
                                          cluster = cutree(tree,k=cut))))

}

#'IVsummary
#' @export
#'
#' @examples
#' IVdata_grouped <- IVsummary(SCKT_allSs,
#'                             grouping.vars =c("subject","group",
#'                              "half","practice","other_type"))
IVsummary <- function(data,
                      grouping.vars = c("subject","group","practice","other_type"),
                      measure.vars = c("prac_score","final_score"),
                      fn = "mean") {
  summary <- data %>%
    group_by_(.dots = grouping.vars)  %>%
    summarise_each_(fn, measure.vars)
  return(summary)
}

#' @export
SCKT_condSummary <- function(data) {
  CDxSS <-SCKT_allSs %>%
    group_by(subject, cluster, half, practice, other_type,prac_score) %>%
    summarise(count = n(),
              final_acc=mean(final_score)) %>%
    ungroup() %>%
    complete(c(subject,cluster), half,c(practice,other_type,prac_score))

  CDxCluster <-SCKT_allSs %>%
    group_by(cluster, half, practice, other_type,prac_score) %>%
    summarise(count = n(),
              final_acc=mean(final_score)) %>%
    ungroup() %>%
    complete(cluster, half,c(practice,other_type,prac_score))

  CD <- SCKT_allSs %>%
    group_by(half, practice, other_type,prac_score) %>%
    summarise(count = n(),
              final_acc=mean(final_score)) %>%
    ungroup() %>%
    complete(half,c(practice,other_type,prac_score))

  return(list(subject=CDxSS, clusters = CDxCluster, group = CD))
}


#' @export
SCKT_jointSummary <- function(data) {

  JDxSS <- SCKT_allSs %>% # filter(practice=='T') %>%
    group_by(subject, cluster, half, practice, other_type,prac_score,final_score) %>%
    summarise(cell_count = n()) %>% # of observations )) %>%
    ungroup() %>%
    complete(c(subject,cluster),half,c(practice,other_type,prac_score,final_score),
             fill=list(cell_count = 0)) %>%
    left_join(SCKT_allSs %>%
              group_by(subject, cluster, half, practice, other_type) %>%
              summarise(maxCount = n())) %>%
    mutate(cell_percent = cell_count/maxCount)

  JDxCluster <- SCKT_allSs %>% # filter(practice=='T') %>%
    group_by(cluster, half, practice, other_type,prac_score,final_score) %>%
    summarise(cell_count = n()) %>% # of observations )) %>%
    ungroup() %>%
    complete(cluster,half,c(practice,other_type,prac_score,final_score),
             fill=list(cell_count = 0)) %>%
    left_join(SCKT_allSs %>%
                group_by(cluster, half, practice, other_type) %>%
                summarise(maxCount = n())) %>%
    mutate(cell_percent = cell_count/maxCount)

  JD <- SCKT_allSs %>% # filter(practice=='T') %>%
    group_by(half, practice, other_type,prac_score,final_score) %>%
    summarise(cell_count = n()) %>% # of observations contributing to mean)) %>%
    ungroup() %>%
    complete(half,c(practice,other_type,prac_score,final_score),
             fill=list(cell_count = 0)) %>%
    left_join(SCKT_allSs %>%
                group_by(half, practice, other_type) %>%
                summarise(maxCount = n())) %>%
    mutate(cell_percent = cell_count/maxCount)

  return(list(subject= JDxSS, clusters = JDxCluster, group=JD))
}
