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
CFR_PCL <- function(free= c(ER=.53,LR=.3,TR =.3, FR=.1,Tmin=1, Tmax=30,
                            lambda=.8, alpha = .02),
                    fixed = c(theta=.5,nFeat=100,nSim=1000,nList=15,Time=90),
                    summarised = TRUE) {

  p <- c(free,fixed)
  if (!paramBounds(p)) {
    stop("Parameter out of allowed range")
  }

  set.seed(456)
  mxn <-  p['nSim']*p['nList'] #dimensions precalculation

  if (!is.na(p['alpha'])) {
    p['alpha'] <- 1/p['alpha']
    p['beta'] <- (p['alpha']/p['ER']) - p['alpha']

    # get the list-wise encoding rates from beta hyper distribution
    ERates <- matrix(rbeta(n = p['nSim'],p['alpha'],p['beta']),
                     nrow=p['nSim'],ncol=p['nList'])
  } else {
    ERates <- p['ER']
    # initial learning
  }

  mem <- matrix(rbinom(mxn,p['nFeat'], ERates),
                nrow=p['nSim'],ncol=p['nList'])
  thresh <- matrix(rbinom(mxn,p['nFeat'], p['theta']),
                   nrow=p['nSim'],ncol=p['nList'])

  #practice test
  prac <- freeRecall(mem,thresh, Tmin = p['Tmin'], Tmax = p['Tmax'],
                     Time = p['Time'], lambda=p['lambda'])

  # study practice
  restudyStrengths <- study(mem, nFeatures=p['nFeat'],
                            LR = p['LR'], FR = p['FR'])
  restudy<-freeRecall(restudyStrengths, thresh, Tmin = p['Tmin'], Tmax = p['Tmax'],
                      Time = p['Time'], lambda=p['lambda'])

  # test practice
  testStrengths <- test(mem,  nFeatures=p['nFeat'],
                        thresh = thresh, acc = prac$Acc, LR = p['LR'],
                        TR = p['TR'], FR = p['FR'])
  tested <- freeRecall(testStrengths$mem, testStrengths$thresh,
                     Tmin = p['Tmin'], Tmax = p['Tmax'],
                     Time = p['Time'], lambda=p['lambda'])

  # Putting the output together
  order <- rbind(prac$order,restudy$order,tested$order)
  RT <-rbind(prac$RT,restudy$RT,tested$RT)
  rec <- rbind(prac$recoverable,restudy$recoverable,tested$recoverable)
  acc <- rbind(prac$Acc,restudy$Acc,tested$Acc)

  # Sorting the output
  RT <- t(sapply(seq(nrow(RT)),function(x) RT[x,order[x,]]))
  rec <-t(sapply(seq(nrow(rec)), function(x) rec[x, order[x,]]))
  acc <-t(sapply(seq(nrow(acc)), function(x) acc[x, order[x,]]))

  # Reshaping the output
  acc <-melt(acc, varnames=c("class","order"),value.name = "acc")
  RT <- melt(RT, varnames=c("class","order"),value.name = "RT")
  rec <- melt(rec, varnames=c("class","order"),value.name = "rec")

  preds <- Reduce(function(x,y) left_join(x,y, by = c("class", "order")),
                  x=list(acc,RT,rec)) %>%
    mutate(class =  rep(rep(c("np","sp","tp"), each = nrow(prac$Acc)),
                        ncol(prac$Acc)),
           unrec = !rec & !acc,
           timeout = !acc & rec) %>%
    select(-rec) %>%
    group_by(class,order)

  # Check summarise switch
  if (summarised) {
    preds <- preds %>%
      summarise(unrec = mean(unrec),
                timeout = mean(timeout),
                RT = median(RT[acc]),
                acc = mean(acc))
  }

  return(preds)
  # Check if we want a density estimate of RTs
  # For correct items only
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
  acc  <- RT %>%
    summarise(acc=mean(acc))
  dist <- RT %>%
    do(d = density(.$RT[.$acc], bw=1,n=time*10,from=.1,to=time)) %>%
    do(data_frame(class = .$class, order = .$order, RT = .$d$x,y=.$d$y)) %>%
    mutate(y = replace(y,y<=0, .Machine$double.xmin)) %>%
    group_by_(.dots = as.character(groups(RT))) %>%
    mutate(y = (y/sum(y))*acc$acc[acc$class==class[1] & acc$order==order[1]]) %>%
    group_by(class) %>%
    mutate(y =y/sum(y))
  return(dist)
}

#' @export
summarise_CFR_PCL <- function(preds) {
    summarise(preds,
              unrec = mean(unrec),
              timeout = mean(timeout),
              RT = median(RT[acc]),
              acc = mean(acc))
}
