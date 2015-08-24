
#' @export
SCKT_PCLss <- function(free= c(ER=.52,LR=.2,TR =.03, F1=.05,space=.03),
                       fixed = c(theta=.5,nFeat=100,nSim=1000,
                                 nList=48,Tmin=NA, Tmax=NA, lambda=NA,Time=NA)) {

  p <- c(free,fixed)
  if (!paramBounds(p)) {
    return(1000000)
  }
  set.seed(456)
  mxn <-  p['nSim']*p['nList'] #dimensions precalculation

  # C_ _C  S_ _S T_ _T
  # CC TT SS
  # CT TC ST TS ST TS

  # Initial Learning
  # Cue 1
  memC1 <- matrix(rbinom(mxn,p['nFeat'], p['ER']),nrow=p['nSim'],ncol=p['nList'])
  # Cue 2
  memC2 <- matrix(rbinom(mxn,p['nFeat'], p['ER']),nrow=p['nSim'],ncol=p['nList'])
  # Common Thresholds for all items
  thresh <- matrix(rbinom(mxn,p['nFeat'], p['theta']),nrow=p['nSim'],ncol=p['nList'])

  # Practice test C1 and C2
  pracC1 <- cuedRecall(memC1,thresh, p['space'])
  pracC2 <- cuedRecall(memC2,thresh, p['space'])


  #control, no practice
  C1strengths <- memC1
  C2strengths <- memC2 - rbinom(mxn, memC2, p['F1'])

  # study practice effects
  S1strengths <- study(memC1, nFeatures=p['nFeat'], LR = p['LR'])
  S2strengths <- study(memC2, nFeatures=p['nFeat'],
                        LR = p['LR'], FR = p['F1'])
  # Final test on study practiced items

  # test practice effects
  T1strengths <- test(mem = memC1, nFeatures=p['nFeat'],thresh = thresh,
                      acc = pracC1, LR = p['LR'], TR = p['TR'])
  T2strengths<- test(mem = memC2, nFeatures=p['nFeat'], thresh = thresh,
                     acc = pracC2, LR = p['LR'], TR = p['TR'], FR = p['F1'])
  T1 <- cuedRecall(T1strengths$mem, T1strengths$thresh, p['space'])
  T2 <- cuedRecall(T2strengths$mem, T2strengths$thresh, p['space'])
  TTstrengths <- test(mem = T1strengths$mem, nFeatures=p['nFeat'],thresh = T1strengths$thresh,
                      acc = T1, LR = p['LR'], TR = p['TR'])

  # Conditions 1 and 2 C_ and _C
  C1 <- pracC1
  C2 <- cuedRecall(C2strengths,thresh, p['space'])

  # Conditions 3 and 4 S_ and _S
  S1 <- cuedRecall(S1strengths,thresh, p['space'])
  S2 <- cuedRecall(S2strengths,thresh, p['space'])


  # Conditions 5 and 6 T_ and _T
  T1 <- T1
  T2 <- T2

  # Conditions 7 and 8 C^C and CC^
  CC1 <- pracC1
  CC2 <- C2

  # Conditions 9 and 10 C^S and CS^
  CS1 <- pracC1
  CS2 <- S2

  # Conditions 11 and 12 C^T and CT^
  CT1 <- cuedRecall(mem = memC1,thresh = T1strengths$thresh,p['space'])
  CT2 <- T2

  # Conditions 13 and 14 S^C and SC^
  SC1 <- S1
  SC2 <- C2

  # Conditions 15 and 16 S^S and SS^
  SS1 <- S1
  SS2 <- S2

  # Conditions 17 and 18 S^T and ST^
  ST1 <- cuedRecall(mem = S1strengths,thresh = T2strengths$thresh,p['space'])
  ST2 <- T2

  # Conditions 19 and 20 T^C and TC^
  TC1 <- T1
  TC2 <- cuedRecall(mem = memC2,thresh = T1strengths$thresh,p['space'])

  # Conditions 21 and 22 T^S and TS^
  TS1 <- T1
  TS2 <- cuedRecall(mem = S2strengths,thresh = T1strengths$thresh,p['space'])

  # Conditions 23 and 24 T^T and TT^
  TT1 <- cuedRecall(TTstrengths$mem, TTstrengths$thresh, p['space'])
  TT2 <- cuedRecall(TTstrengths$mem, TTstrengths$thresh, p['space'])


  # Average over simulations
  avgs <- lapply(list(pracC1=pracC1, pracC2 = pracC2,
                      C1 = C1, C2 = C2,
                      S1= S1, S2 = S2,
                      T1 = T1, #T1plus = T1[pracC1], T1neg = T1[!pracC1],
                      T1_p_f = (pracC1 & T1), T1_p_nf =  (pracC1 & !T1),
                      T1_np_f = (!pracC1 & T1), T1_np_nf = (!pracC1 & !T1),
                      T2 = T2, #T2plus = T2[pracC2], T2neg = T2[!pracC2],
                      T2_p_f = (pracC2 & T2), T2_p_nf =  (pracC2 & !T2),
                      T2_np_f = (!pracC2 & T2), T2_np_nf = (!pracC2 & !T2),
                      CC1 = CC1, CC2 = CC2,
                      CS1 = CS1, CS2 = CS2,
                      CT1 = CT1, CT2 = CT2,
                      CT2plus = CT2[pracC2], CT2neg = CT2[!pracC2],
                      CT2_p_f = (pracC2 & CT2), CT2_p_nf =  (pracC2 & !CT2),
                      CT2_np_f = (!pracC2 & CT2), CT2_np_nf = (!pracC2 & !CT2),
                      SC1 = SC1, SC2 = SC2,
                      SS1 = SS1, SS2 = SS2,
                      ST1 = ST1, ST2 = ST2,
                      ST2plus = ST2[pracC2], ST2neg = ST2[!pracC2],
                      ST2_p_f = (pracC2 & ST2), ST2_p_nf =  (pracC2 & !ST2),
                      ST2_np_f = (!pracC2 & ST2), ST2_np_nf = (!pracC2 & !ST2),
                      TC1 = TC1, TC2 = TC2,
                      TC1plus = TC1[pracC1], TC1neg = TC1[!pracC1],
                      TC1_p_f = (pracC1 & TC1), TC1_p_nf =  (pracC1 & !TC1),
                      TC1_np_f = (!pracC1 & TC1), TC1_np_nf = (!pracC1 & !TC1),
                      TS1 = TC1, TS2 = TC2,
                      TS1plus = TS1[pracC1], TS1neg = TS1[!pracC1],
                      TS1_p_f = (pracC1 & TS1), TS1_p_nf =  (pracC1 & !TS1),
                      TS1_np_f = (!pracC1 & TS1), TS1_np_nf = (!pracC1 & !TS1),
                      TT1 = TT1, TT2 = TT2,
                      TT1plus = TT1[pracC1], TT1neg = TT1[!pracC1],
                      TT1_p_f = (pracC1 & TT1), TT1_p_nf =  (pracC1 & !TT1),
                      TT1_np_f = (!pracC1 & TT1), TT1_np_nf = (!pracC1 & !TT1),
                      TT2plus = TT2[pracC1], TT2neg = TT2[!pracC2],
                      TT2_p_f = (pracC2 & TT2), TT2_p_nf =  (pracC2 & !TT2),
                      TT2_np_f = (!pracC2 & TT2), TT2_np_nf = (!pracC2 & !TT2)),
                 mean)
  return(avgs)

}
