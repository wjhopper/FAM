#' LBL42 PCR
#'
#' Implementation of the PCR model of LB4L2 experiment
#' @param PCR_model_obj
#'
#' @return A data frame
#'
#' @import PCR
#' @export
#'
#' @examples
#' model_params <- initPCR(params = list(ER=.6, LR=.2, FR=.1, TR=.15, TV = .05),
#'                    distribution = "beta", nItems = 20, nSims = 1000,
#'                    nFeatures = 100, time = 10)
#' results <- LB4L2_PCR(model_params)
#'
LB4L2_PCR <- function(PCR_model_obj) {

  control <- study(PCR_model_obj, nCues = 1, tests_per_cue = 1) %>%
    forget(cue = 1) %>%
    cuedRecall(cue = 1, increment = FALSE)

  SC_study <- study(PCR_model_obj, nCues = 2, tests_per_cue = c(1,0)) %>%
    study(cue = 1) %>%
    forget(cue = 1) %>%
    cuedRecall(cue = 1, increment = FALSE)

  OC_study <- study(PCR_model_obj, nCues = 2, tests_per_cue = c(1,0)) %>%
    study(cue = 2) %>%
    forget(cue = 1) %>%
    cuedRecall(cue = 1, increment = FALSE)

  SC_test <- study(PCR_model_obj, nCues = 2, tests_per_cue = c(2,0)) %>%
    cuedRecall(cue = 1) %>%
    forget(cue = 1) %>%
    cuedRecall(cue = 1, increment = FALSE)

  OC_test <- study(PCR_model_obj, nCues = 2, tests_per_cue = c(1,1)) %>%
    cuedRecall(cue = 2) %>%
    forget(cue = 1) %>%
    cuedRecall(cue = 1, increment = FALSE)

  results <- list("SC_test" = SC_test, "SC_study" = SC_study,
                  "OC_test" = OC_test, "OC_study" = OC_study,
                  "control" = control)
  class(results) <- c("LB4L2_PCR", class(results))
  return(results)

}

#' @export
optim_wrapper <- function(parameters, PCR_model_obj, ...) {
  PCR_model_obj <- updatePCRparams(as.list(parameters))
  return(LB4L2_PCR(PCR_model_obj, ...))
}

#' @export
#' @importFrom tidyr separate
#' @importFrom tidyr gather
summary.LB4L2_PCR <- function(x, DV = "recalled") {

  IVresults <- bind_rows(lapply(x, summary), .id = "condition") %>%
    mutate(practice = toupper(substr(practice,1,1)))

  practice_test <- filter(IVresults, practice == "T", test == 1) %>%
    mutate(sameCue = ifelse(grepl("OC", condition), "no", "yes"))

  final_test <- anti_join(IVresults,practice_test,
                          by = c("condition","practice","test")) %>%
    group_by(condition) %>%
    mutate(OCpractice = switch(first(condition),
                               OC_test = "T",
                               OC_study = "S",
                               "N")) %>%
    ungroup()

  # Subset all the data from items tested twice.
  SC_test <- x[["SC_test"]]$recalled[[1]]
  SC_joint <- mapply(`&`,
                     list(cor_inc = SC_test[,,1], inc_inc = !SC_test[,,1],
                          cor_cor = SC_test[,,1], inc_cor = !SC_test[,,1]),
                     list(!SC_test[,,2], !SC_test[,,2], SC_test[,,2], SC_test[,,2]))

  OC_tp <- x[["OC_test"]]$recalled[[2]]
  OC_final <- x[["OC_test"]]$recalled[[1]]
  OC_joint <- mapply(`&`,
                    list(cor_inc = OC_tp, inc_inc = !OC_tp,
                         cor_cor = OC_tp, inc_cor = !OC_tp),
                    list(!OC_final, !OC_final, OC_final, OC_final))

  ## Joint Accuracy Results
  SC_joint_acc <- colMeans(SC_joint)
  OC_joint_acc <- colMeans(OC_joint)

  ## Joint RT results
  raw_joint_RT <- list(OC_joint_RT = sapply(x[['OC_test']]$RT, `[`, OC_joint[,"cor_cor"]),
                       SC_joint_RT = apply(x[['SC_test']]$RT[[1]], 3, `[`, SC_joint[,"cor_cor"])) %>%
    lapply(`dimnames<-`, (list(NULL, c("practice", "final")))) %>%
    mutate(tbl_df(data.frame(sameCue = c("no","yes"), practice = 1, final =1)),
           raw_joint_RTs = .)

  joint_results <- cbind(as.data.frame(SC_joint_acc), as.data.frame(OC_joint_acc)) %>%
    mutate(condition = rownames(.)) %>%
    separate(condition, c("practice", "final")) %>%
    mutate(practice = ifelse(practice == "cor", 1, 0),
           final = ifelse(final == "cor", 1, 0)) %>%
    gather(sameCue, probability, SC_joint_acc, OC_joint_acc) %>%
    mutate(sameCue = ifelse(sameCue == "SC_joint_acc", "yes", "no")) %>%
    left_join(y = raw_joint_RT, by = c("sameCue", "practice", "final")) %>%
    tbl_df()

  # Conditional Accuracy
  SC_CD <- SC_joint_acc[c("cor_cor", "inc_cor")]/abs((c(0,1) - practice_test$accuracy[practice_test$sameCue=="yes"]))
  OC_CD <- OC_joint_acc[c("cor_cor", "inc_cor")]/abs((c(0,1) - practice_test$accuracy[practice_test$sameCue=="no"]))

  CD_results <- cbind(as.data.frame(SC_CD), as.data.frame(OC_CD)) %>%
    mutate(practice = gsub("_cor", "", rownames(.))) %>%
    gather(sameCue, accuracy, SC_CD, OC_CD) %>%
    mutate(sameCue = ifelse(sameCue == "SC_CD", "yes", "no")) %>%
    mutate(sameCue, practice, accuracy)

  return(list(practice= practice_test, final = final_test,
              conditional = CD_results, joint = joint_results))
}
