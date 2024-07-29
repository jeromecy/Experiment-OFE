source("Initialization.R")

b0 <- 65
b1 <- 0.05
b2 <- -0.0003

x1 <- seq(0,140,length=5)
# x <- c(0,39,50.6,75.4,99.8,124.6)
y1 <- b0 + b1*x1+b2*x1^2

dat1 <- data.frame(x=x1,y=y1)

x2 <- seq(0,140,length=100)
# x <- c(0,39,50.6,75.4,99.8,124.6)
y2 <- b0 + b1*x2+b2*x2^2

dat2 <- data.frame(x=x2,y=y2)

ggplot(data=dat1) + geom_point(aes(x,y),size=3) + 
  scale_x_continuous(labels = x1,breaks = x1) +
  geom_line(data=dat2,aes(x,y)) + thm1 +
  xlab("Nitrogen Level") + ylab("Yield (quintals/ha)")


y3 <- b0 + b1*x1
dat3 <- data.frame(x=x1,y=y3)

ggplot(data=dat3) + geom_point(aes(x,y),size=3) + 
  scale_x_continuous(labels = x1,breaks = x1) +
  geom_line(aes(x,y)) + thm1 +
  xlab("Nitrogen Level") + ylab("Yield (quintals/ha)")






y <- b0 + b1*x

plot(x,y)


####### parameter preparation #########

M <- 100 ## simulation iterations

paralist <- list(
  Nitro  = seq(0,140,length=5),
  nRep   = 4,
  nRow   = 93,
  b_0    = 65,
  b_1    = 0.05,
  b_2    = -0.0003,
  sig_u0 = 5,
  sig_u1 = 0.01,
  sig_u2 = 0.0001,
  rho_range = 0.15,
  rho_row   = 0.5
)
lev    <- 2  # linear trend
nRange <- 20
nRow   <- paralist$nRow


##### simulation study data generating functions #######
SimuGenDat_Lin <- function(paras,covmat,eta){
  
  nitro  <- sample(paras$Nitro)
  
  nRep   <- paras$nRep
  nRange <- nRep * length(nitro)
  nRow   <- paras$nRow
  N      <- nRange * nRow
  
  SystLevel <- rep(nitro,nRep)
  RandLevel <- c(replicate(nRep,sample(nitro,length(nitro))))
  
  b_0 <- paras$b_0
  b_1 <- paras$b_1
  
  LKJMat <- rlkjcorr(1,2,eta)
  sig_u0 <- paras$sig_u0
  sig_u1 <- paras$sig_u1
  
  Sigma_u <- diag(c(sig_u0,sig_u1)) %*% LKJMat %*% diag(c(sig_u0,sig_u1))
  
  if(covmat == "NS"){
    Sigma_S = diag(N)
  }else if(covmat == "AR1"){
    Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))
  }else if(covmat == "Matern"){
    Sigma_S = cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                     sp.type = "matern", 
                     sp.par = c(1,1),
                     smoothness = 1.5)$V
  }
  
  SigmaAll <- kronecker(Sigma_S,Sigma_u)
  
  Us_mat <- mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)
  Us_mat <- matrix(Us_mat,2,N)
  u_0 <- Us_mat[1,]
  u_1 <- Us_mat[2,]
  
  LinearDat <- data.frame(Range = rep(1:nRange,each=nRow),
                          Row   = rep(1:nRow,nRange),
                          SystNitro = rep(SystLevel,each=nRow),
                          RandNitro = rep(RandLevel,each=nRow),
                          Rep   = rep(1:4,each = 5*nRow),
                          gridId = 1:N,
                          b0 = b_0,
                          b1 = b_1,
                          u0 = u_0,
                          u1 = u_1)
  
  ## systematic, linear y = x * a + b + e
  LinearDat$SYield<- LinearDat$b0 + LinearDat$u0 + (LinearDat$b1 + LinearDat$u1) * LinearDat$SystNitro + rnorm(N)
  LinearDat$RYield<- LinearDat$b0 + LinearDat$u0 + (LinearDat$b1 + LinearDat$u1) * LinearDat$RandNitro + rnorm(N)
  
  return(LinearDat)
}
SimuGenDat_Qua <- function(paras,covmat,eta){
  
  nitro  <- sample(paras$Nitro)
  nRep   <- paras$nRep
  nRange <- nRep * length(nitro)
  nRow   <- paras$nRow
  N      <- nRange * nRow
  
  SystLevel <- rep(nitro,nRep)
  RandLevel <- c(replicate(nRep,sample(nitro,length(nitro))))
  
  b_0 <- paras$b_0
  b_1 <- paras$b_1
  b_2 <- paras$b_2
  
  LKJMat <- rlkjcorr(1,3,eta)
  sig_u0 <- paras$sig_u0
  sig_u1 <- paras$sig_u1
  sig_u2 <- paras$sig_u2
  
  Sigma_u <- diag(c(sig_u0,sig_u1,sig_u2)) %*% LKJMat %*% diag(c(sig_u0,sig_u1,sig_u2))
  
  if(covmat == "NS"){
    Sigma_S = diag(N)
  }else if(covmat == "AR1"){
    Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))
  }else if(covmat == "Matern"){
    Sigma_S = cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                     sp.type = "matern", 
                     sp.par = c(1,1),
                     smoothness = 1.5)$V
  }
  
  SigmaAll <- kronecker(Sigma_S,Sigma_u)
  
  Us_mat <- mvrnorm(1, mu= rep(c(0,0,0),N), Sigma = SigmaAll)
  Us_mat <- matrix(Us_mat,3,N)
  u_0 <- Us_mat[1,]
  u_1 <- Us_mat[2,]
  u_2 <- Us_mat[3,]
  
  QuadraDat <-data.frame(Range = rep(1:nRange,each=nRow),
                         Row = rep(1:nRow,nRange),
                         SystNitro = rep(SystLevel,each=nRow),
                         RandNitro = rep(RandLevel,each=nRow),
                         Rep = rep(1:4,each = 5*nRow),
                         gridId = 1:N,
                         b0 = b_0, b1 = b_1, b2 = b_2,
                         u0 = u_0, u1 = u_1, u2 = u_2)
  QuadraDat$SystNitro2 <- QuadraDat$SystNitro^2
  QuadraDat$RandNitro2 <- QuadraDat$RandNitro^2
  
  ## Quadratic  y = a * x2 + b * x + c + e

  QuadraDat$SYield <- QuadraDat$b0 + QuadraDat$u0 + (QuadraDat$b1 + QuadraDat$u1) * QuadraDat$SystNitro + 
    (QuadraDat$b2 + QuadraDat$u2) * QuadraDat$SystNitro2 + rnorm(N)
  QuadraDat$RYield <- QuadraDat$b0 + QuadraDat$u0 + (QuadraDat$b1 + QuadraDat$u1) * QuadraDat$RandNitro + 
    (QuadraDat$b2 + QuadraDat$u2) * QuadraDat$RandNitro2 + rnorm(N)
  
  return(QuadraDat)
}

################## simulation weak correlation eta = 1 ##################
###### no spatial simulation (NS) ######
MSE_ns_LS5 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear systematic, b0, b1
MSE_ns_LS9 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_ns_LSbw <- array(0,dim=c(M,lev,nRange*nRow))

BW_ns_LS <- numeric(M)

MSE_ns_LR5 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear randomized
MSE_ns_LR9 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_ns_LRbw <- array(0,dim=c(M,lev,nRange*nRow))

BW_ns_LR <- numeric(M)


MSE_ns_QS5 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic systematic
MSE_ns_QS9 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_ns_QSbw <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_ns_QS <- numeric(M)

MSE_ns_QR5 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic randomized
MSE_ns_QR9 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_ns_QRbw <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_ns_QR <- numeric(M)

set.seed(2021)
######## basic gwr NS Linear simulation, eta = 1  ##########
eta <- 1  ## weak correlation, 0.1 high correlation
for(i in 1:M){
  LinearDat <- SimuGenDat_Lin(paralist,"NS",eta)  ## higher correlation

  sp_LinearDat <- SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                         LinearDat)
  
  gwr_syst5 <- gwr.basic(SYield ~ SystNitro,
                         data = sp_LinearDat, bw=5,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_LS5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_syst5$SDF$Intercept
  MSE_ns_LS5[i,2,] <- (LinearDat$b1 + LinearDat$u1) -gwr_syst5$SDF$SystNitro
  
  
  gwr_syst9 <- gwr.basic(SYield ~ SystNitro,
                         data = sp_LinearDat, bw=9,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_LS9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_syst9$SDF$Intercept
  MSE_ns_LS9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_syst9$SDF$SystNitro
  
  
  bwLinSyst <- bw.gwr(SYield ~ SystNitro, data = sp_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_syst.bw <- gwr.basic(SYield ~ SystNitro,
                           data = sp_LinearDat, bw=bwLinSyst,
                           kernel = "gaussian")  ## boxcar kernel
  
  BW_ns_LS[i] <- bwLinSyst
  
  MSE_ns_LSbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_syst.bw$SDF$Intercept
  MSE_ns_LSbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_syst.bw$SDF$SystNitro
  
  
  
  gwr_rand5 <- gwr.basic(RYield ~ RandNitro,
                         data = sp_LinearDat, bw = 5,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_LR5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_rand5$SDF$Intercept
  MSE_ns_LR5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_rand5$SDF$RandNitro
  
  
  
  gwr_rand9 <- gwr.basic(RYield ~ RandNitro,
                         data = sp_LinearDat, bw=9, ## bw = 5
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_LR9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_rand9$SDF$Intercept
  MSE_ns_LR9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_rand9$SDF$RandNitro
  
  
  bwLinRand <- bw.gwr(RYield ~ RandNitro, data = sp_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_rand.bw <- gwr.basic(RYield ~ RandNitro,
                           data = sp_LinearDat, bw=bwLinRand, ## bw = 5
                           kernel = "gaussian")  ## boxcar kernel
  
  BW_ns_LR[i] <- bwLinRand
  
  MSE_ns_LRbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_rand.bw$SDF$Intercept
  MSE_ns_LRbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_rand.bw$SDF$RandNitro
  
  cat(i)
}

######## basic gwr NS Quadratic simulation,  eta = 1  ###########
eta <- 1 ## weak correlation
for(i in 1:M){
  QuadraDat <- SimuGenDat_Qua(paralist,"NS",eta)
  
  sp_QuadraDat <- SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),
                                         QuadraDat)
  
  gwr_qua.syst5 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                             data = sp_QuadraDat, bw=5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_QS5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.syst5$SDF$Intercept
  MSE_ns_QS5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.syst5$SDF$SystNitro
  MSE_ns_QS5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.syst5$SDF$SystNitro2
  
  gwr_qua.syst9 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                             data = sp_QuadraDat, bw=9,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_QS9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.syst9$SDF$Intercept
  MSE_ns_QS9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.syst9$SDF$SystNitro
  MSE_ns_QS9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.syst9$SDF$SystNitro2
  
  
  quad.syst.bw <- bw.gwr(SYield ~ SystNitro + SystNitro2, 
                         data = sp_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_qua.syst.bw <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                               data = sp_QuadraDat, bw=quad.syst.bw,
                               kernel = "gaussian")  ## boxcar kernel
  BW_ns_QS[i] <- quad.syst.bw
  
  MSE_ns_QSbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.syst.bw$SDF$Intercept
  MSE_ns_QSbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.syst.bw$SDF$SystNitro
  MSE_ns_QSbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.syst.bw$SDF$SystNitro2
  
  
  gwr_qua.rand5 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                             data = sp_QuadraDat, bw = 5, 
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_QR5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.rand5$SDF$Intercept
  MSE_ns_QR5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.rand5$SDF$RandNitro
  MSE_ns_QR5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.rand5$SDF$RandNitro2
  
  
  gwr_qua.rand9 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                             data = sp_QuadraDat, bw = 9, 
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_QR9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.rand9$SDF$Intercept
  MSE_ns_QR9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.rand9$SDF$RandNitro
  MSE_ns_QR9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.rand9$SDF$RandNitro2
  
  
  quad.rand.bw <- bw.gwr(RYield ~ RandNitro + RandNitro2, 
                         data = sp_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_qua.rand.bw <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                               data = sp_QuadraDat, bw = quad.rand.bw, 
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_ns_QR[i] <- quad.rand.bw
  
  MSE_ns_QRbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.rand.bw$SDF$Intercept
  MSE_ns_QRbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.rand.bw$SDF$RandNitro
  MSE_ns_QRbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.rand.bw$SDF$RandNitro2
  
  cat(i)
  
}

####### ar1 x ar1 Sigma simulation #######
set.seed(1257)

MSE_ar1_LS5 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear systematic, b0, b1
MSE_ar1_LS9 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_ar1_LSbw <- array(0,dim=c(M,lev,nRange*nRow))

BW_ar1_LS <- numeric(M)

MSE_ar1_LR5 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear randomized
MSE_ar1_LR9 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_ar1_LRbw <- array(0,dim=c(M,lev,nRange*nRow))

BW_ar1_LR <- numeric(M)

MSE_ar1_QS5 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic systematic
MSE_ar1_QS9 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_ar1_QSbw <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_ar1_QS <- numeric(M)

MSE_ar1_QR5 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic randomized
MSE_ar1_QR9 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_ar1_QRbw <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_ar1_QR <- numeric(M)


######## basic gwr linear data ##########
eta <- 1
for(i in 1:M){
  LinearDat <- SimuGenDat_Lin(paralist,"AR1",eta)  
  sp_ar1_LinearDat <- SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                         LinearDat)
  gwr_ar1_syst5 <- gwr.basic(SYield ~ SystNitro,
                         data = sp_ar1_LinearDat, bw=5,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_LS5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst5$SDF$Intercept
  MSE_ar1_LS5[i,2,] <- (LinearDat$b1 + LinearDat$u1) -gwr_ar1_syst5$SDF$SystNitro
  
  
  gwr_ar1_syst9 <- gwr.basic(SYield ~ SystNitro,
                         data = sp_ar1_LinearDat, bw=9,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_LS9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst9$SDF$Intercept
  MSE_ar1_LS9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst9$SDF$SystNitro
  
  bwLinSyst <- bw.gwr(SYield ~ SystNitro, data = sp_ar1_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_ar1_syst.bw <- gwr.basic(SYield ~ SystNitro,
                           data = sp_ar1_LinearDat, bw=bwLinSyst,
                           kernel = "gaussian")  ## boxcar kernel
  
  BW_ar1_LS[i] <- bwLinSyst
  
  MSE_ar1_LSbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst.bw$SDF$Intercept
  MSE_ar1_LSbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst.bw$SDF$SystNitro
  
  gwr_ar1_rand5 <- gwr.basic(RYield ~ RandNitro,
                         data = sp_ar1_LinearDat, bw = 5,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_LR5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand5$SDF$Intercept
  MSE_ar1_LR5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand5$SDF$RandNitro
  
  gwr_ar1_rand9 <- gwr.basic(RYield ~ RandNitro,
                         data = sp_ar1_LinearDat, bw=9, ## bw = 5
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_LR9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand9$SDF$Intercept
  MSE_ar1_LR9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand9$SDF$RandNitro
  
  bwLinRand <- bw.gwr(RYield ~ RandNitro, data = sp_ar1_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_ar1_rand.bw <- gwr.basic(RYield ~ RandNitro,
                           data = sp_ar1_LinearDat, bw=bwLinRand, ## bw = 5
                           kernel = "gaussian")  ## boxcar kernel
  
  BW_ar1_LR[i] <- bwLinRand
  
  MSE_ar1_LRbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand.bw$SDF$Intercept
  MSE_ar1_LRbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand.bw$SDF$RandNitro
  
  cat(i)
}

######## basic gwr Quadratic data ##########
eta <- 1
for(i in 1:M){
  QuadraDat <- SimuGenDat_Qua(paralist,"AR1",eta)
  
  sp_ar1_QuadraDat <- SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
  
  gwr_ar1_qua.syst5 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                             data = sp_ar1_QuadraDat, bw=5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_QS5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst5$SDF$Intercept
  MSE_ar1_QS5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst5$SDF$SystNitro
  MSE_ar1_QS5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst5$SDF$SystNitro2
  
  
  gwr_ar1_qua.syst9 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                             data = sp_ar1_QuadraDat, bw=9,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_QS9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst9$SDF$Intercept
  MSE_ar1_QS9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst9$SDF$SystNitro
  MSE_ar1_QS9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst9$SDF$SystNitro2
  
  
  quad.syst.bw <- bw.gwr(SYield ~ SystNitro + SystNitro2, 
                         data = sp_ar1_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_ar1_qua.syst.bw <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                               data = sp_ar1_QuadraDat, bw=quad.syst.bw,
                               kernel = "gaussian")  ## boxcar kernel
  BW_ar1_QS[i] <- quad.syst.bw
  
  MSE_ar1_QSbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst.bw$SDF$Intercept
  MSE_ar1_QSbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst.bw$SDF$SystNitro
  MSE_ar1_QSbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst.bw$SDF$SystNitro2
  
  
  gwr_ar1_qua.rand5 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                             data = sp_ar1_QuadraDat, bw = 5, 
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_QR5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand5$SDF$Intercept
  MSE_ar1_QR5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand5$SDF$RandNitro
  MSE_ar1_QR5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand5$SDF$RandNitro2
  
  
  gwr_ar1_qua.rand9 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                             data = sp_ar1_QuadraDat, bw = 9, 
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_QR9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand9$SDF$Intercept
  MSE_ar1_QR9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand9$SDF$RandNitro
  MSE_ar1_QR9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand9$SDF$RandNitro2
  
  
  quad.rand.bw <- bw.gwr(RYield ~ RandNitro + RandNitro2, 
                         data = sp_ar1_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_ar1_qua.rand.bw <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                               data = sp_ar1_QuadraDat, bw = quad.rand.bw, 
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_ar1_QR[i] <- quad.rand.bw
  
  MSE_ar1_QRbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand.bw$SDF$Intercept
  MSE_ar1_QRbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand.bw$SDF$RandNitro
  MSE_ar1_QRbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand.bw$SDF$RandNitro2
  
  cat(i)
  
}


###### Matern simulation #####

set.seed(8515)

MSE_mat_LS5 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear systematic, b0, b1
MSE_mat_LS9 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_mat_LSbw <- array(0,dim=c(M,lev,nRange*nRow))

BW_mat_LS <- numeric(M)

MSE_mat_LR5 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear randomized
MSE_mat_LR9 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_mat_LRbw <- array(0,dim=c(M,lev,nRange*nRow))

BW_mat_LR <- numeric(M)

MSE_mat_QS5 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic systematic
MSE_mat_QS9 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_mat_QSbw <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_mat_QS <- numeric(M)

MSE_mat_QR5 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic randomized
MSE_mat_QR9 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_mat_QRbw <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_mat_QR <- numeric(M)


######## basic gwr linear data ##########
eta <- 1
for(i in 1:M){
  LinearDat <- SimuGenDat_Lin(paralist,"Matern",eta)  
  sp_mat_LinearDat <- SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                             LinearDat)
  gwr_mat_syst5 <- gwr.basic(SYield ~ SystNitro,
                             data = sp_mat_LinearDat, bw=5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_LS5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst5$SDF$Intercept
  MSE_mat_LS5[i,2,] <- (LinearDat$b1 + LinearDat$u1) -gwr_mat_syst5$SDF$SystNitro
  
  
  gwr_mat_syst9 <- gwr.basic(SYield ~ SystNitro,
                             data = sp_mat_LinearDat, bw=9,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_LS9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst9$SDF$Intercept
  MSE_mat_LS9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst9$SDF$SystNitro
  
  bwLinSyst <- bw.gwr(SYield ~ SystNitro, data = sp_mat_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_mat_syst.bw <- gwr.basic(SYield ~ SystNitro,
                               data = sp_mat_LinearDat, bw=bwLinSyst,
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_mat_LS[i] <- bwLinSyst
  
  MSE_mat_LSbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst.bw$SDF$Intercept
  MSE_mat_LSbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst.bw$SDF$SystNitro
  
  gwr_mat_rand5 <- gwr.basic(RYield ~ RandNitro,
                             data = sp_mat_LinearDat, bw = 5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_LR5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand5$SDF$Intercept
  MSE_mat_LR5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand5$SDF$RandNitro
  
  gwr_mat_rand9 <- gwr.basic(RYield ~ RandNitro,
                             data = sp_mat_LinearDat, bw=9, ## bw = 5
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_LR9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand9$SDF$Intercept
  MSE_mat_LR9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand9$SDF$RandNitro
  
  bwLinRand <- bw.gwr(RYield ~ RandNitro, data = sp_mat_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_mat_rand.bw <- gwr.basic(RYield ~ RandNitro,
                               data = sp_mat_LinearDat, bw=bwLinRand, ## bw = 5
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_mat_LR[i] <- bwLinRand
  
  MSE_mat_LRbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand.bw$SDF$Intercept
  MSE_mat_LRbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand.bw$SDF$RandNitro
  
  cat(i)
}

######## basic gwr Quadratic data ##########
eta <- 1
for(i in 1:M){
  QuadraDat <- SimuGenDat_Qua(paralist,"Matern",eta)
  
  sp_mat_QuadraDat <- SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
  
  gwr_mat_qua.syst5 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                 data = sp_mat_QuadraDat, bw=5,
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_QS5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst5$SDF$Intercept
  MSE_mat_QS5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst5$SDF$SystNitro
  MSE_mat_QS5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst5$SDF$SystNitro2
  
  
  gwr_mat_qua.syst9 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                 data = sp_mat_QuadraDat, bw=9,
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_QS9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst9$SDF$Intercept
  MSE_mat_QS9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst9$SDF$SystNitro
  MSE_mat_QS9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst9$SDF$SystNitro2
  
  
  quad.syst.bw <- bw.gwr(SYield ~ SystNitro + SystNitro2, 
                         data = sp_mat_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_mat_qua.syst.bw <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                   data = sp_mat_QuadraDat, bw=quad.syst.bw,
                                   kernel = "gaussian")  ## boxcar kernel
  BW_mat_QS[i] <- quad.syst.bw
  
  MSE_mat_QSbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst.bw$SDF$Intercept
  MSE_mat_QSbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst.bw$SDF$SystNitro
  MSE_mat_QSbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst.bw$SDF$SystNitro2
  
  
  gwr_mat_qua.rand5 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                 data = sp_mat_QuadraDat, bw = 5, 
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_QR5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand5$SDF$Intercept
  MSE_mat_QR5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand5$SDF$RandNitro
  MSE_mat_QR5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand5$SDF$RandNitro2
  
  
  gwr_mat_qua.rand9 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                 data = sp_mat_QuadraDat, bw = 9, 
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_QR9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand9$SDF$Intercept
  MSE_mat_QR9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand9$SDF$RandNitro
  MSE_mat_QR9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand9$SDF$RandNitro2
  
  
  quad.rand.bw <- bw.gwr(RYield ~ RandNitro + RandNitro2, 
                         data = sp_mat_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_mat_qua.rand.bw <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                   data = sp_mat_QuadraDat, bw = quad.rand.bw, 
                                   kernel = "gaussian")  ## boxcar kernel
  
  BW_mat_QR[i] <- quad.rand.bw
  
  MSE_mat_QRbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand.bw$SDF$Intercept
  MSE_mat_QRbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand.bw$SDF$RandNitro
  MSE_mat_QRbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand.bw$SDF$RandNitro2
  
  cat(i)
  
}


################## simulation high correlation eta =0.1 ##################
lev    <- 2  # linear trend
nRange <- 20
nRow   <- paralist$nRow
###### no spatial simulation (NS) ######
MSE_ns_LS5_eta01 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear systematic, b0, b1
MSE_ns_LS9_eta01 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_ns_LSbw_eta01 <- array(0,dim=c(M,lev,nRange*nRow))

BW_ns_LS_eta01 <- numeric(M)

MSE_ns_LR5_eta01 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear randomized
MSE_ns_LR9_eta01 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_ns_LRbw_eta01 <- array(0,dim=c(M,lev,nRange*nRow))

BW_ns_LR_eta01 <- numeric(M)


MSE_ns_QS5_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic systematic
MSE_ns_QS9_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_ns_QSbw_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_ns_QS_eta01 <- numeric(M)

MSE_ns_QR5_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic randomized
MSE_ns_QR9_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_ns_QRbw_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_ns_QR_eta01 <- numeric(M)

set.seed(2114)
######## basic gwr NS Linear simulation ##########
eta <- 0.1  ## weak correlation, 0.1 high correlation
for(i in 1:M){
  LinearDat <- SimuGenDat_Lin(paralist,"NS",eta)  ## higher correlation
  
  sp_LinearDat <- SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                         LinearDat)
  
  gwr_syst5 <- gwr.basic(SYield ~ SystNitro,
                         data = sp_LinearDat, bw=5,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_LS5_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_syst5$SDF$Intercept
  MSE_ns_LS5_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) -gwr_syst5$SDF$SystNitro
  
  
  gwr_syst9 <- gwr.basic(SYield ~ SystNitro,
                         data = sp_LinearDat, bw=9,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_LS9_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_syst9$SDF$Intercept
  MSE_ns_LS9_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_syst9$SDF$SystNitro
  
  
  bwLinSyst <- bw.gwr(SYield ~ SystNitro, data = sp_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_syst.bw <- gwr.basic(SYield ~ SystNitro,
                           data = sp_LinearDat, bw=bwLinSyst,
                           kernel = "gaussian")  ## boxcar kernel
  
  BW_ns_LS_eta01[i] <- bwLinSyst
  
  MSE_ns_LSbw_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_syst.bw$SDF$Intercept
  MSE_ns_LSbw_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_syst.bw$SDF$SystNitro
  
  
  
  gwr_rand5 <- gwr.basic(RYield ~ RandNitro,
                         data = sp_LinearDat, bw = 5,
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_LR5_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_rand5$SDF$Intercept
  MSE_ns_LR5_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_rand5$SDF$RandNitro
  
  
  
  gwr_rand9 <- gwr.basic(RYield ~ RandNitro,
                         data = sp_LinearDat, bw=9, ## bw = 5
                         kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_LR9_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_rand9$SDF$Intercept
  MSE_ns_LR9_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_rand9$SDF$RandNitro
  
  
  bwLinRand <- bw.gwr(RYield ~ RandNitro, data = sp_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_rand.bw <- gwr.basic(RYield ~ RandNitro,
                           data = sp_LinearDat, bw=bwLinRand, ## bw = 5
                           kernel = "gaussian")  ## boxcar kernel
  
  BW_ns_LR_eta01[i] <- bwLinRand
  
  MSE_ns_LRbw_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_rand.bw$SDF$Intercept
  MSE_ns_LRbw_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_rand.bw$SDF$RandNitro
  
  cat(i)
}

######## basic gwr NS Quadratic simulation  ###########
eta <- 0.1 ## weak correlation
for(i in 1:M){
  QuadraDat <- SimuGenDat_Qua(paralist,"NS",eta)
  
  sp_QuadraDat <- SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),
                                         QuadraDat)
  
  gwr_qua.syst5 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                             data = sp_QuadraDat, bw=5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_QS5_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.syst5$SDF$Intercept
  MSE_ns_QS5_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.syst5$SDF$SystNitro
  MSE_ns_QS5_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.syst5$SDF$SystNitro2
  
  gwr_qua.syst9 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                             data = sp_QuadraDat, bw=9,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_QS9_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.syst9$SDF$Intercept
  MSE_ns_QS9_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.syst9$SDF$SystNitro
  MSE_ns_QS9_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.syst9$SDF$SystNitro2
  
  
  quad.syst.bw <- bw.gwr(SYield ~ SystNitro + SystNitro2, 
                         data = sp_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_qua.syst.bw <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                               data = sp_QuadraDat, bw=quad.syst.bw,
                               kernel = "gaussian")  ## boxcar kernel
  BW_ns_QS_eta01[i] <- quad.syst.bw
  
  MSE_ns_QSbw_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.syst.bw$SDF$Intercept
  MSE_ns_QSbw_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.syst.bw$SDF$SystNitro
  MSE_ns_QSbw_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.syst.bw$SDF$SystNitro2
  
  
  gwr_qua.rand5 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                             data = sp_QuadraDat, bw = 5, 
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_QR5_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.rand5$SDF$Intercept
  MSE_ns_QR5_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.rand5$SDF$RandNitro
  MSE_ns_QR5_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.rand5$SDF$RandNitro2
  
  
  gwr_qua.rand9 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                             data = sp_QuadraDat, bw = 9, 
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ns_QR9_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.rand9$SDF$Intercept
  MSE_ns_QR9_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.rand9$SDF$RandNitro
  MSE_ns_QR9_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.rand9$SDF$RandNitro2
  
  
  quad.rand.bw <- bw.gwr(RYield ~ RandNitro + RandNitro2, 
                         data = sp_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_qua.rand.bw <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                               data = sp_QuadraDat, bw = quad.rand.bw, 
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_ns_QR_eta01[i] <- quad.rand.bw
  
  MSE_ns_QRbw_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_qua.rand.bw$SDF$Intercept
  MSE_ns_QRbw_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_qua.rand.bw$SDF$RandNitro
  MSE_ns_QRbw_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_qua.rand.bw$SDF$RandNitro2
  
  cat(i)
  
}

####### ar1 x ar1 Sigma simulation #######
set.seed(3217882)

MSE_ar1_LS5_eta01 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear systematic, b0, b1
MSE_ar1_LS9_eta01 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_ar1_LSbw_eta01 <- array(0,dim=c(M,lev,nRange*nRow))

BW_ar1_LS_eta01 <- numeric(M)

MSE_ar1_LR5_eta01 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear randomized
MSE_ar1_LR9_eta01 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_ar1_LRbw_eta01 <- array(0,dim=c(M,lev,nRange*nRow))

BW_ar1_LR_eta01 <- numeric(M)

MSE_ar1_QS5_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic systematic
MSE_ar1_QS9_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_ar1_QSbw_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_ar1_QS_eta01 <- numeric(M)

MSE_ar1_QR5_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic randomized
MSE_ar1_QR9_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_ar1_QRbw_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_ar1_QR_eta01 <- numeric(M)


######## basic gwr linear data ##########
eta <- 0.1
for(i in 1:M){
  LinearDat <- SimuGenDat_Lin(paralist,"AR1",eta)  
  sp_ar1_LinearDat <- SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                             LinearDat)
  gwr_ar1_syst5 <- gwr.basic(SYield ~ SystNitro,
                             data = sp_ar1_LinearDat, bw=5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_LS5_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst5$SDF$Intercept
  MSE_ar1_LS5_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) -gwr_ar1_syst5$SDF$SystNitro
  
  
  gwr_ar1_syst9 <- gwr.basic(SYield ~ SystNitro,
                             data = sp_ar1_LinearDat, bw=9,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_LS9_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst9$SDF$Intercept
  MSE_ar1_LS9_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst9$SDF$SystNitro
  
  bwLinSyst <- bw.gwr(SYield ~ SystNitro, data = sp_ar1_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_ar1_syst.bw <- gwr.basic(SYield ~ SystNitro,
                               data = sp_ar1_LinearDat, bw=bwLinSyst,
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_ar1_LS_eta01[i] <- bwLinSyst
  
  MSE_ar1_LSbw_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst.bw$SDF$Intercept
  MSE_ar1_LSbw_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst.bw$SDF$SystNitro
  
  gwr_ar1_rand5 <- gwr.basic(RYield ~ RandNitro,
                             data = sp_ar1_LinearDat, bw = 5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_LR5_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand5$SDF$Intercept
  MSE_ar1_LR5_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand5$SDF$RandNitro
  
  gwr_ar1_rand9 <- gwr.basic(RYield ~ RandNitro,
                             data = sp_ar1_LinearDat, bw=9, ## bw = 5
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_LR9_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand9$SDF$Intercept
  MSE_ar1_LR9_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand9$SDF$RandNitro
  
  bwLinRand <- bw.gwr(RYield ~ RandNitro, data = sp_ar1_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_ar1_rand.bw <- gwr.basic(RYield ~ RandNitro,
                               data = sp_ar1_LinearDat, bw=bwLinRand, ## bw = 5
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_ar1_LR_eta01[i] <- bwLinRand
  
  MSE_ar1_LRbw_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand.bw$SDF$Intercept
  MSE_ar1_LRbw_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand.bw$SDF$RandNitro
  
  cat(i)
}

######## basic gwr Quadratic data ##########
eta <- 0.1
for(i in 1:M){
  QuadraDat <- SimuGenDat_Qua(paralist,"AR1",eta)
  
  sp_ar1_QuadraDat <- SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
  
  gwr_ar1_qua.syst5 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                 data = sp_ar1_QuadraDat, bw=5,
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_QS5_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst5$SDF$Intercept
  MSE_ar1_QS5_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst5$SDF$SystNitro
  MSE_ar1_QS5_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst5$SDF$SystNitro2
  
  
  gwr_ar1_qua.syst9 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                 data = sp_ar1_QuadraDat, bw=9,
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_QS9_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst9$SDF$Intercept
  MSE_ar1_QS9_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst9$SDF$SystNitro
  MSE_ar1_QS9_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst9$SDF$SystNitro2
  
  
  quad.syst.bw <- bw.gwr(SYield ~ SystNitro + SystNitro2, 
                         data = sp_ar1_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_ar1_qua.syst.bw <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                   data = sp_ar1_QuadraDat, bw=quad.syst.bw,
                                   kernel = "gaussian")  ## boxcar kernel
  BW_ar1_QS_eta01[i] <- quad.syst.bw
  
  MSE_ar1_QSbw_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst.bw$SDF$Intercept
  MSE_ar1_QSbw_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst.bw$SDF$SystNitro
  MSE_ar1_QSbw_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst.bw$SDF$SystNitro2
  
  
  gwr_ar1_qua.rand5 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                 data = sp_ar1_QuadraDat, bw = 5, 
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_QR5_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand5$SDF$Intercept
  MSE_ar1_QR5_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand5$SDF$RandNitro
  MSE_ar1_QR5_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand5$SDF$RandNitro2
  
  
  gwr_ar1_qua.rand9 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                 data = sp_ar1_QuadraDat, bw = 9, 
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_ar1_QR9_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand9$SDF$Intercept
  MSE_ar1_QR9_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand9$SDF$RandNitro
  MSE_ar1_QR9_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand9$SDF$RandNitro2
  
  
  quad.rand.bw <- bw.gwr(RYield ~ RandNitro + RandNitro2, 
                         data = sp_ar1_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_ar1_qua.rand.bw <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                   data = sp_ar1_QuadraDat, bw = quad.rand.bw, 
                                   kernel = "gaussian")  ## boxcar kernel
  
  BW_ar1_QR_eta01[i] <- quad.rand.bw
  
  MSE_ar1_QRbw_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand.bw$SDF$Intercept
  MSE_ar1_QRbw_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand.bw$SDF$RandNitro
  MSE_ar1_QRbw_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand.bw$SDF$RandNitro2
  
  cat(i)
  
}


###### Matern simulation #####

set.seed(77128)

MSE_mat_LS5_eta01 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear systematic, b0, b1
MSE_mat_LS9_eta01 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_mat_LSbw_eta01 <- array(0,dim=c(M,lev,nRange*nRow))

BW_mat_LS_eta01 <- numeric(M)

MSE_mat_LR5_eta01 <- array(0,dim=c(M,lev,nRange*nRow))  ## linear randomized
MSE_mat_LR9_eta01 <- array(0,dim=c(M,lev,nRange*nRow))
MSE_mat_LRbw_eta01 <- array(0,dim=c(M,lev,nRange*nRow))

BW_mat_LR_eta01 <- numeric(M)

MSE_mat_QS5_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic systematic
MSE_mat_QS9_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_mat_QSbw_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_mat_QS_eta01 <- numeric(M)

MSE_mat_QR5_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))  ## quadratic randomized
MSE_mat_QR9_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))
MSE_mat_QRbw_eta01 <- array(0,dim=c(M,lev+1,nRange*nRow))

BW_mat_QR_eta01 <- numeric(M)


######## basic gwr linear data ##########
eta <- 0.1
for(i in 1:M){
  LinearDat <- SimuGenDat_Lin(paralist,"Matern",eta)  
  sp_mat_LinearDat <- SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                             LinearDat)
  gwr_mat_syst5 <- gwr.basic(SYield ~ SystNitro,
                             data = sp_mat_LinearDat, bw=5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_LS5_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst5$SDF$Intercept
  MSE_mat_LS5_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) -gwr_mat_syst5$SDF$SystNitro
  
  
  gwr_mat_syst9 <- gwr.basic(SYield ~ SystNitro,
                             data = sp_mat_LinearDat, bw=9,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_LS9_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst9$SDF$Intercept
  MSE_mat_LS9_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst9$SDF$SystNitro
  
  bwLinSyst <- bw.gwr(SYield ~ SystNitro, data = sp_mat_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_mat_syst.bw <- gwr.basic(SYield ~ SystNitro,
                               data = sp_mat_LinearDat, bw=bwLinSyst,
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_mat_LS_eta01[i] <- bwLinSyst
  
  MSE_mat_LSbw_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst.bw$SDF$Intercept
  MSE_mat_LSbw_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst.bw$SDF$SystNitro
  
  gwr_mat_rand5 <- gwr.basic(RYield ~ RandNitro,
                             data = sp_mat_LinearDat, bw = 5,
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_LR5_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand5$SDF$Intercept
  MSE_mat_LR5_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand5$SDF$RandNitro
  
  gwr_mat_rand9 <- gwr.basic(RYield ~ RandNitro,
                             data = sp_mat_LinearDat, bw=9, ## bw = 5
                             kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_LR9_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand9$SDF$Intercept
  MSE_mat_LR9_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand9$SDF$RandNitro
  
  bwLinRand <- bw.gwr(RYield ~ RandNitro, data = sp_mat_LinearDat,
                      kernel = "gaussian",
                      approach = "AICc",
                      adaptive = F)
  
  gwr_mat_rand.bw <- gwr.basic(RYield ~ RandNitro,
                               data = sp_mat_LinearDat, bw=bwLinRand, ## bw = 5
                               kernel = "gaussian")  ## boxcar kernel
  
  BW_mat_LR_eta01[i] <- bwLinRand
  
  MSE_mat_LRbw_eta01[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand.bw$SDF$Intercept
  MSE_mat_LRbw_eta01[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand.bw$SDF$RandNitro
  
  cat(i)
}

######## basic gwr Quadratic data ##########
eta <- 0.1
for(i in 1:M){
  QuadraDat <- SimuGenDat_Qua(paralist,"Matern",eta)
  
  sp_mat_QuadraDat <- SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
  
  gwr_mat_qua.syst5 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                 data = sp_mat_QuadraDat, bw=5,
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_QS5_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst5$SDF$Intercept
  MSE_mat_QS5_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst5$SDF$SystNitro
  MSE_mat_QS5_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst5$SDF$SystNitro2
  
  
  gwr_mat_qua.syst9 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                 data = sp_mat_QuadraDat, bw=9,
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_QS9_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst9$SDF$Intercept
  MSE_mat_QS9_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst9$SDF$SystNitro
  MSE_mat_QS9_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst9$SDF$SystNitro2
  
  
  quad.syst.bw <- bw.gwr(SYield ~ SystNitro + SystNitro2, 
                         data = sp_mat_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_mat_qua.syst.bw <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                                   data = sp_mat_QuadraDat, bw=quad.syst.bw,
                                   kernel = "gaussian")  ## boxcar kernel
  BW_mat_QS_eta01[i] <- quad.syst.bw
  
  MSE_mat_QSbw_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst.bw$SDF$Intercept
  MSE_mat_QSbw_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst.bw$SDF$SystNitro
  MSE_mat_QSbw_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst.bw$SDF$SystNitro2
  
  
  gwr_mat_qua.rand5 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                 data = sp_mat_QuadraDat, bw = 5, 
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_QR5_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand5$SDF$Intercept
  MSE_mat_QR5_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand5$SDF$RandNitro
  MSE_mat_QR5_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand5$SDF$RandNitro2
  
  
  gwr_mat_qua.rand9 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                 data = sp_mat_QuadraDat, bw = 9, 
                                 kernel = "gaussian")  ## boxcar kernel
  
  MSE_mat_QR9_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand9$SDF$Intercept
  MSE_mat_QR9_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand9$SDF$RandNitro
  MSE_mat_QR9_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand9$SDF$RandNitro2
  
  
  quad.rand.bw <- bw.gwr(RYield ~ RandNitro + RandNitro2, 
                         data = sp_mat_QuadraDat,
                         kernel = "gaussian",
                         approach = "AICc",
                         adaptive = F)
  
  gwr_mat_qua.rand.bw <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                                   data = sp_mat_QuadraDat, bw = quad.rand.bw, 
                                   kernel = "gaussian")  ## boxcar kernel
  
  BW_mat_QR_eta01[i] <- quad.rand.bw
  
  MSE_mat_QRbw_eta01[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand.bw$SDF$Intercept
  MSE_mat_QRbw_eta01[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand.bw$SDF$RandNitro
  MSE_mat_QRbw_eta01[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand.bw$SDF$RandNitro2
  
  cat(i)
  
}


# save.image("~/Curtin/Research/BayesTechOFE/ExptOut1Ksimu.RData")

