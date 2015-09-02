#' CFR_PCL
#'
#' @param free
#' @param fixed
#' @param yoke
#' @param data
#' @param fittingRT
#' @param fittingAcc
#'
#' @importFrom KernSmooth bkde
#' @import reshape2
#' @import dplyr
#' @export
#'
CFR_PCL <- function(free= c(ER=.53,LR=.3,TR =.3, FR=.1,Tmin=2, Tmax=10,
                            lambda=.8, alpha = 13),
                    fixed = c(theta=.5,nFeat=100,nSim=1000,nList=15,Time=90),
                    return_dist=TRUE) {

  p <- c(free,fixed)
  if (!paramBounds(p)) {
    stop("Parameter out of allowed range")
  }

  set.seed(456)
  mxn <-  p['nSim']*p['nList'] #dimensions precalculation

  if (!is.na(p['alpha'])) {
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

  RT <-rbind(prac$RT,restudy$RT,tested$RT)
  order <- rbind(prac$order,restudy$order,tested$order)
  RT <-t(sapply(seq(nrow(RT)),
                function(x) c(RT[x, order[x,][!is.na(RT[x,order[x,]])]],
                              RT[x, order[x,][is.na(RT[x,order[x,]])]])))
  acc <-melt(!is.na(RT), varnames=c("class","order"),value.name = "acc")
  RT <- melt(RT, varnames=c("class","order"),value.name = "RT")
  preds <- left_join(acc,RT, by = c("class", "order")) %>%
    mutate(class =  rep(rep(c("np","sp","tp"),
                            each = nrow(prac$Acc)),
                        ncol(prac$Acc))) %>%
    group_by(class,order) %>%
    summarise(acc = mean(acc),
              RT = median(RT,na.rm = TRUE))
  if (!anyNA(p[c("Tmin","Tmax","theta","Time")]) && return_dist)   {
    dist <- data.frame(class = rep(c('np','sp','tp'),each=13500), # 13500 = 15 items * 900 points
                       rbind(RTdis(prac$RT, prac$order, p['Time']),
                             RTdis(restudy$RT, restudy$order, p['Time']),
                             RTdis(tested$RT, tested$order, p['Time']),
                       row.names = NULL))
    return(list(preds = preds,distribution = dist))
  } else {
    return(preds = preds)
  }
}


RTdis <- function(RT = NULL, order = NULL, Time= NULL) {

  # for each rank ordering recall position (Nlist),
  # the KS density at each 1/10 second of recall test period (TestTime)

  RTdist=matrix(0,nrow=10*Time,ncol=ncol(RT))
  RT <-t(sapply(seq(nrow(RT)),
                function(x) c(RT[x, order[x,][!is.na(RT[x,order[x,]])]],
                              RT[x, order[x,][is.na(RT[x,order[x,]])]])
  ))
  for (n in 1:ncol(RT)) {
    RTs<- RT[!is.na(RT[,n]),n]
    if (length(RTs)>2) {
#       D <- bkde(RTs,bandwidth=1,gridsize=900,range.x=c(.1,90))
      D <- density(RTs,bw=1,n=900,from=.1,to=90)
      D$y[D$y <= 0] <-  (0.5)*.Machine$double.xmin
      height <- D$y/sum(D$y)
      RTdist[,n] <- (length(RTs)/nrow(RT))*height
    }
  }
  RTdist <- melt(RTdist/sum(RTdist),varnames=c("RTrounded","order"),value.name = "RTdist")
  RTdist$RTrounded <- RTdist$RTrounded/10
  return(RTdist)
}

LL <- function(obs,pred) {
  likelihoods <- inner_join(obs,pred)
  #   likelihoods$RTdist[likelihoods$RTdist == 0] <- (0.5)*.Machine$double.xmin
  err <- -sum(log(likelihoods$RTdist))
  return(err)
}


SSE <- function(obs,pred) {
  obs <- obs %>% group_by(class,order) %>%
    summarise(acc = sum(score)/4)
  data <- pred %>% group_by(class,order) %>%
    summarise(pred_acc = mean(pred_acc)) %>%
    left_join(obs)
  data$acc[is.na(data$acc)] <- 0 # set na's for accuracy for zeros
  err <- sum((data$pred_acc-data$acc)^2)
  return(err)
}
