#' @import optimx
#' @import PCR
#' @import dplyr
#' @import whoppeR
NULL


#' Title
#'
#' @param model
#' @param inpar
#' @param debugLevel
#'
#' @return A Model list, updated with a list of  the fitting results and model predictions
#' @export
#'
#' @examples
#'  model <- list(par = c(ER=.52,LR=.2,TR =.03, F1=.05,space=.03),
#'                fn = PCRss, method="Nelder-Mead", itnmax = 1000,
#'                fix= c(theta=.5,nFeat=100,nSim=1000,nList=15),
#'                IVdata= data, jointData = data, condData = data,
#'                fitting = TRUE, name = 'std_ss'))
#'  results <- fitLB4L(model,inpar = TRUE)
#'
fitLB4L <- function(model,inpar=FALSE,debugLevel = 0) {

  if (inpar) {
    packs <- requireNamespace(c('foreach', 'doParallel'), quietly = TRUE)
    if (all(packs)) {
      library(foreach)
      library(doParallel)
      cl <- makeCluster(detectCores()-1)
      registerDoParallel(cl)
      #     library(doRNG)
      #     registerDoRNG(456)
    } else {
      warning("Packages foreach and doParallel not found to do parallel processing, falling back")
      inpar = FALSE
    }
  }

  checkRequiredParams(c(model$par,model$fixed), model$fn)
#   if (debugLevel[1]>0){
#     trace(model$fn, browser, at=debugLevel[2])
# #     setBreakpoint('PCR.R', debugLevel[2], envir = environment())
#   }
  if (inpar) {
    logfile <- tempfile("parlog", fileext = ".txt")
    writeLines(paste('[',Sys.time(),']',"INIT parlog"),
               con=logfile,sep='\n')
    clusterExport(cl,"logfile",envir = environment())
    model$results <- foreach(j =unique(model$IVdata$subject), .verbose=T,
                             .packages=c("optimx","PCR","whoppeR")) %dopar% {
        sink(logfile, append=TRUE)
        cat(paste("Fitting subject", j,"\n"))
        sink()
        m  <- model
        m$IVdata <- model$IVdata[model$IVdata$subject == j,]
        m$jointData <- model$jointData[model$jointData$subject == j,]
        fit <- optimx(par=m$par,fn=m$fn, method = m$method, itnmax = m$itnmax,
                      fixed=m$fixed,IVdata=m$IVdata, jointData = m$jointData,
                      fitting= m$fitting)
        best <- as.vector(coef(fit))
        names(best) <- colnames(coef(fit))
        preds <- m$fn(free=best,fixed=m$fixed,IVdata=m$IVdata, jointData = m$jointData,
                      fitting=FALSE)$preds
        end <- list(fit=fit, preds = preds )
      }

    stopCluster(cl)
    cat(paste('[',Sys.time(),']',"Finished Fitting, Goodbye"),
        con=logfile,sep='\n')

  } else {
    model$results <- vector(mode = 'list', length = length(unique(model$IVdata$subject)))
    for (j in unique(model$IVdata$subject)) {

      message(paste("Fitting subject", j))
      m <- model
      m$IVdata <- model$IVdata[model$IVdata$subject == j,]
      m$jointData <- model$jointData[model$jointData$subject == j,]
      fit <- optimx(par=m$par,fn=m$fn, method = m$method, itnmax = m$itnmax,
                    fixed=m$fixed,IVdata=m$IVdata, jointData = m$jointData,
                    fitting= m$fitting)
      best <- as.vector(coef(fit))
      names(best) <- colnames(coef(fit))
      preds <- m$fn(free=best,fixed=m$fixed,IVdata=m$IVdata, jointData = m$jointData,
                     fitting=FALSE)$preds
      m$results[[j]] <- list(fit=fit, preds = preds )
    }
  }

  if (debugLevel[1]>0){
    untrace(PCR)
  }

#   resultsFile <- paste(model$objective_fcn, model$name, 'results',sep='_')
#   assign(resultsFile, model[c("results","objective_fcn","fixed","par")])
#   save(list = resultsFile,
#        file = file.path("data",paste(model$objective_fcn,
#                                      model$name, 'results.rda',
#                                      sep='_')))
  return(model)
}



#' Title
#'
#' @param free Named vector of free parameters. See Details for possible names.
#' @param fixed Named vector of free parameters. See Details for possible names.
#' @param IVdata
#' @param jointdata
#' @param fitting Controls  the return value.  See Value for details.
#'
#' @return If \code{fitting = TRUE}, scalar vector is returned,
#'  useful for optimization routine. If \code{fittting = FALSE}, list of data frames with predicted joint,
#'  conditional, and IV condition-level predictions is returned.
#'
#' @examples
#' LB4L_PCRss(free= c(ER=.53,LR=.15,TR =.1, F1=.1,space=.03),
#'            fixed = c(theta=.5,nFeat=100,nSim=1000,nList=15,
#'                     Tmin=NA, Tmax=NA, lambda=NA,Time=NA),
#'                     fitting=FALSE)
#' @export
LB4L_PCRss <- function(free= c(ER=.52,LR=.2,TR =.03, F1=.05,space=.03),
                       fixed = c(theta=.5,nFeat=100,nSim=1000,
                                 nList=15,Tmin=NA, Tmax=NA, lambda=NA,Time=NA),
                       ...) {

  p <- c(free,fixed)
  if (!paramBounds(p)) {
    return(1000000)
  }

  set.seed(456)

  mxn <-  p['nSim']*p['nList'] #dimensions precalculation

  # Initial Learning
  # Cue 1
  init_mem_C1 <- matrix(rbinom(mxn,p['nFeat'], p['ER']),nrow=p['nSim'],ncol=p['nList'])
  # Cue 2
  init_mem_C2 <- matrix(rbinom(mxn,p['nFeat'], p['ER']),nrow=p['nSim'],ncol=p['nList'])
  # Common Thresholds for all items
  init_thresh <- matrix(rbinom(mxn,p['nFeat'], p['theta']),nrow=p['nSim'],ncol=p['nList'])

  # Practice test
  prac <- cuedRecall(init_mem_C1,init_thresh, p['space'])


  #control no practice
  controlStrengths <- init_mem_C1 - rbinom(mxn, init_mem_C1, p['F1'])
  controlAcc <- cuedRecall(controlStrengths, init_thresh, p['space'])

  # study practice effects
  restudyStrengths <- study(mem = init_mem_C1, nFeatures=p['nFeat'],
                            LR = p['LR'], FR = p['F1'])
  # Final test on study practiced items
  restudyAcc<-cuedRecall(restudyStrengths, init_thresh, p['space'])

  # test practice effects
  testStrengths <- test(mem = init_mem_C1, nFeatures=p['nFeat'],
                        thresh = init_thresh, acc = prac, LR = p['LR'],
                        TR = p['TR'], FR = p['F1'])
  # Final test on test practiced items
  testAcc <- cuedRecall(testStrengths$mem, testStrengths$thresh, p['space'])

  # no practice, other cue test practice
  testOCStrengths <- init_mem_C2 - rbinom(mxn, init_mem_C2, p['F1'])
  testOCAcc<- cuedRecall(testOCStrengths, testStrengths$thresh, p['space'])

  IV <- data.frame(practice = c("T","T","C","C","S"),
                   other_type = c(NA,"C",NA,"T",NA),
                   prac_acc = c(mean(prac),mean(prac), NA, NA, NA),
                   other_prac_acc = c(NA,NA,mean(prac),mean(prac),NA),
                   final_acc = unlist(lapply(list(testAcc, NA , controlAcc,
                                                  testOCAcc, restudyAcc),
                                      mean)))
  CD <- data.frame(practice = c("T","T","C","C"),
                   other_type = c(NA, NA, "T","T"),
                   prac_score = c(1,0,1,0),
                   final_acc = unlist(lapply(list(testAcc[prac], testAcc[!prac],
                                                  testOCAcc[prac], testOCAcc[!prac]),
                                      mean)))
  JD <- data.frame(practice = c("T","T","T","T","C","C","C","C"),
                   other_type = c(NA,NA,NA,NA,"T","T","T","T"),
                   prac_acc = c(1,0,1,0,1,0,1,0),
                   final_acc =c(1,1,0,0,1,1,0,0),
                   acc = unlist(lapply(list((prac & testAcc), (!prac & testAcc),
                                            (prac & !testAcc),(!prac & !testAcc),
                                            (prac & testOCAcc), (!prac & testOCAcc),
                                            (prac & !testOCAcc),(!prac & !testOCAcc)),
                                       mean)))

  return(list(IV = IV, CD = CD, JD = JD))
}
# Average over simulations
# avgs <- lapply(list(prac=prac, # done
#                     C = controlAcc, # done
#                     S= restudyAcc, # done
#                     CT = testOCAcc, # done
#                     CTplus = testOCAcc[prac], # done
#                     CTneg = testOCAcc[!prac], #done
#                     CT_not_prac_not_final = (!prac & !testOCAcc),
#                     CT_not_prac_final = (!prac & testOCAcc),
#                     CT_prac_not_final = (prac & !testOCAcc),
#                     CT_prac_final = (prac & testOCAcc),
#                     `T` = testAcc, # done
#                     Tplus = testAcc[prac], # done
#                     Tneg = testAcc[!prac], # done
#                     T_not_prac_not_final = (!prac & !testAcc), # done
#                     T_not_prac_final = (!prac & testAcc), # done
#                     T_prac_not_final = (prac & !testAcc), # done
#                     T_prac_final = (prac & testAcc)), # done
#                mean)
#   preds <- unlist(avgs)
#   obs <- c(IVdata$final_acc[IVdata$practice %in% c('C','S') & IVdata$other_type!='T'],
#            jointData$acc[(jointData$prac_score %in% 0:1 |
#                            jointData$other_prac_acc %in% 0:1) &
#                          jointData$final_score %in% 0:1])
#   N <- c(IVdata$n[IVdata$practice %in% c('C','S') & IVdata$other_type!='T'],
#          jointData$cell_count[(jointData$prac_score %in% 0:1 |
#                        jointData$other_prac_acc %in% 0:1) &
#                      jointData$final_score %in% 0:1])
#
#   err <- binomialLL(obs=obs[1:2],pred=preds[c("C","S")],N=p['nList']) +
#     multinomialLL(obs = obs[3:10],pred = preds[grepl('T_', names(preds))],
#                   N = p['nList'])
#
#   if (fitting) {
#     return(err)
#   } else {
#     return(list(error=err, preds =avgs))
#   }
# }
