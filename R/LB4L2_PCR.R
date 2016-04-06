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
#' model_params <- initPCRparams(params = list(ER=.6, LR=.2, FR=.1, TR=.15, TV = .05),
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

  if (DV == "recalled") {
    fun <- colMeans
  } else if (DV == "RT"){
    fun = I
  }

    IVresults <- bind_rows(lapply(x, summary), .id = "condition")

    practice_test <- filter(IVresults, practice == "T", test == 1) %>%
      mutate(sameCue = ifelse(grepl("OC", condition), "no", "yes")) %>%
      select(sameCue, practice, accuracy, RT)

    final_test <- filter(IVresults, !(practice == "T" & test == 1)) %>%
      group_by(condition) %>%
      mutate(OCpractice = switch(condition, OC_test = "T", OC_study = "S", "C"))

    x <- lapply(x, function(y){
                    y$practice[is.na(y$practice)] <- "control"
                    return(y)
                   })

    SC_test <- x[["SC_test"]][[DV]][[1]]
    SC_joint <- fun(mapply(`&`,
                         list(cor_inc = SC_test[,,1], inc_inc = !SC_test[,,1],
                              cor_cor = SC_test[,,1], inc_cor = !SC_test[,,1]),
                         list(!SC_test[,,2], !SC_test[,,2], SC_test[,,2], SC_test[,,2])))
    SC_CD <- SC_joint[c("cor_cor", "inc_cor")]/abs((c(0,1) - practice_test$accuracy[practice_test$sameCue=="yes"]))

    OC_tp <- x[["OC_test"]][[DV]][[2]]
    OC_final <- x[["OC_test"]][[DV]][[1]]
    OC_joint <- fun(mapply(`&`,
                                list(cor_inc = OC_tp, inc_inc = !OC_tp,
                                     cor_cor = OC_tp, inc_cor = !OC_tp),
                                list(!OC_final, !OC_final, OC_final, OC_final)))
    OC_CD <- OC_joint[c("cor_cor", "inc_cor")]/abs((c(0,1) - practice_test$accuracy[practice_test$sameCue=="no"]))

    joint_results <- cbind(as.data.frame(SC_joint), as.data.frame(OC_joint)) %>%
      mutate(condition = rownames(.)) %>%
      separate(condition, c("practice", "final")) %>%
      mutate(practice = ifelse(practice == "cor", 1, 0),
             final = ifelse(final == "cor", 1, 0)) %>%
      gather(sameCue, probability, SC_joint, OC_joint) %>%
      mutate(sameCue = ifelse(sameCue == "SC_joint", "yes", "no")) %>%
      mutate(sameCue, practice, final, probability)

    CD_results <- cbind(as.data.frame(SC_CD), as.data.frame(OC_CD)) %>%
      mutate(practice = gsub("_cor", "", rownames(.))) %>%
      gather(sameCue, accuracy, SC_CD, OC_CD) %>%
      mutate(sameCue = ifelse(sameCue == "SC_CD", "yes", "no")) %>%
      mutate(sameCue, practice, accuracy)

    return(list(practice= practice_test, final = final_test,
                conditional = CD_results, joint = joint_results))

}
