#' CFR_PCL
#' An implementation of the PCL model to the FAM_CFR data set
#'
#' @param free
#'  A Named vector of numeric values.
#' @param fixed
#'  A Named vector of numeric values.
#' @param returnDist
#'  Logical, controls whether a density estimate of the simulated RT distribution
#'  is returned along with accuracy esitmates
#'
#' @importFrom KernSmooth bkde
#' @import reshape2
#' @import dplyr
#' @import PCL
#' @export
#'
CFR_PCL <- function(free = c(ER=.53,LR=.3,Ta =49.5,TR = .1, FR=.1,Tmin=1,Tmax=30,lambda=.8),
                    fixed = c(theta=.5,nFeat=100,nSim=1000,nList=15,Time=90),
                    summarised = TRUE) {

  p <- c(free,fixed)
  if (!paramBounds(p)) {
    stop("Parameter out of allowed range")
  }

  set.seed(456)
  mxn <-  p['nSim']*p['nList'] #dimensions precalculation

  # we want to sample from a beta distribution with the same
  # mean and sd as binomial(N = nFeat,p=ER) distribution.
  # The problem is that the binomial is [0,nFeat], and the equations we
  # to solve for beta distributions alpha and beta from mean and sd
  # are for beta bounded between zero and 1 aka the 2 parameter beta.
  # So we need to divide the mean and varianace of the binomial by nFeat
  # and nFeat^2. Then we can use our equations, and multiply by nFeat after sampling

  binomVAR <- p['nFeat']*p['ER']*(1-p['ER'])
  binomM <- p['nFeat']*p['ER']
  beta_pars = betaABfromMeanSD(mean = binomM/p['nFeat'],
                               sd = sqrt(binomVAR/p['nFeat']^2))
  mem <- matrix(rbeta(mxn, beta_pars$a, beta_pars$b),
                nrow=p['nSim'],ncol=p['nList']) * p['nFeat']
  thresh <- matrix(rbeta(mxn, p['Ta'],  p['Ta']),
                   nrow=p['nSim'],ncol=p['nList']) * p['nFeat']

  # practice test
  prac <- freeRecall(mem,thresh, Tmin = p['Tmin'], Tmax = p['Tmax'],
                     Time = p['Time'], lambda=p['lambda'])

  # study practice
  restudyStrengths <- study_beta(mem=mem, nFeatures=p['nFeat'],
                            LR = p['LR'], FR = p['FR'])
  restudy<-freeRecall(restudyStrengths, thresh, Tmin = p['Tmin'], Tmax = p['Tmax'],
                      Time = p['Time'], lambda=p['lambda'])

  # test practice
  testStrengths <- test_beta(mem=mem,  nFeatures=p['nFeat'],
                        thresh = thresh, acc = prac$Acc, LR = p['LR'],
                        TR = p['TR'], FR = p['FR'])
  tested <- freeRecall(testStrengths$mem, testStrengths$thresh,
                     Tmin = p['Tmin'], Tmax = p['Tmax'],
                     Time = p['Time'], lambda=p['lambda'])

  # Putting the output together
  order <- rbind(prac$order,restudy$order,tested$order)
  RT <-rbind(prac$RT,restudy$RT,tested$RT)
  RTcor <- rbind(prac$RTcor,restudy$RTcor,tested$RTcor)
  rec <- rbind(prac$recoverable,restudy$recoverable,tested$recoverable)
  acc <- rbind(prac$Acc,restudy$Acc,tested$Acc)

  # Sorting the output
  for (x in 1:p['nSim']) {
    RT[x,] <- RT[x,order[x,]]
    RTcor[x,] <- RTcor[x,order[x,]]
    rec[x,] <- rec[x,order[x,]]
    acc[x,] <- acc[x,order[x,]]
  }

  # Reshaping the output
  acc <-melt(acc, varnames=c("class","memOrder"),value.name = "acc")
  RT <- melt(RT, varnames=c("class","memOrder"),value.name = "memRT")
  RTcor <- melt(RTcor, varnames=c("class","memOrder"),value.name = "obsRT")
  rec <- melt(rec, varnames=c("class","memOrder"),value.name = "rec")

  preds <- Reduce(function(x,y) left_join(x,y, by = c("class", "memOrder")),
                  x=list(acc,RT,RTcor,rec)) %>%
    mutate(class =  rep(rep(c("np","sp","tp"), each = nrow(prac$Acc)),
                        ncol(prac$Acc)),
           sim  = rep(1:p['nSim'], times = p['nList']*3),
           unrec = !rec & !acc,
           timeout = !acc & rec) %>%
    group_by(class,sim,acc) %>%
    mutate(obsOrder = 1:n(),
           obsOrder = replace(obsOrder, acc==FALSE, NA)) %>%
    select(-rec)

  # Check summarise switch
  if (summarised) {
    preds <- preds %>%
      # group_by(class, obsOrder) %>%
      summarise_CFR_PCL()
  }
  return(preds)
}


#' RTdist: Normalized Density Estimates
#'
#' Estimated Gaussian kernel density at each 1/10 second of recall test period,
#' using a bandwidth of 1. Density results are normalized so that the sum
#' of the densities at each output positions is 1.
#'
#' @param RT
#' @param time Totall test duration measured in seconds. Density is estimated
#'  at every tenth of a second between .1 and \code{time}
#'
#' @return A data frame
#' @export
#' @examples
#'  preds <- CFR_PCL(summarised=FALSE)
#'  dist <- preds %>% group_by(class,order) %>% RTdist()

RTdist <- function(RT, time=90) {
  if (all(group_size(RT) == 1)) {
    stop("Cannot estimate densities using summarised model predictions")
  }
  acc  <- RT %>% group_by(class,memOrder) %>%
    summarise(acc=mean(acc))
  safe_density <- failwith(list(x=seq(.1,time,.1), y = 0),
                           density)
  dist <- RT %>%
    filter(!is.na(obsOrder)) %>%
    do(d = safe_density(.$obsRT, bw=1,n=time*10,from=.1,to=time)) %>%
    do(data_frame(class = .$class, obsOrder = .$obsOrder, RT = .$d$x,y=.$d$y)) %>%
#     mutate(y = replace(y,y<=0, .Machine$double.xmin)) %>%
    group_by_(.dots = as.character(groups(RT))) %>%
    mutate(y = (y/sum(y))*acc$acc[acc$class==class[1] & acc$memOrder==obsOrder[1]]) %>%
    group_by(class) %>%
    mutate(y =y/sum(y)) %>%
    rename(order=obsOrder)
  return(dist)
}

#' @export
summarise_CFR_PCL <- function(preds) {
  mem_summary <- preds %>%
    group_by(class,memOrder) %>%
    summarise(unrec = mean(unrec),
              timeout = mean(timeout),
              RT = median(memRT),
              acc = mean(acc)) %>%
    rename(order = memOrder)
  obs_summary <- preds %>%
    filter(!is.na(obsOrder)) %>%
    group_by(class,obsOrder) %>%
    summarise(RTcor = median(obsRT)) %>%
    rename(order = obsOrder)
  return(left_join(mem_summary,obs_summary, by = c("class","order")))
}
