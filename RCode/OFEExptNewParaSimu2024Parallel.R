m0 <- function (mu, Sigma, tol = 1e-06) {
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE, EISPACK = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  list(eSv=eS$vector,ev=ev, dd=diag(sqrt(pmax(ev, 0)), p))
}

m1 <- function(n,mu,eL,fancy=FALSE) {
  p <- length(mu)
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eL$eSv %*% eL$dd %*% t(X)
  if (fancy) {
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
      nm <- dn[[1L]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1) 
      drop(X)
    else t(X)
  } else t(X)
}


fastMN <- function(n,mu,sigma){
  e <-  m0(mu,sigma)
  x2 <- m1(n,mu,eL=e,fancy=TRUE)
}

mu <- rep(2,20)
Sigma <- matrix(0.1,nrow=20,ncol=20)
diag(Sigma) <- 4

set.seed(1001)
x1 <- mvrnorm(1000,mu=mu,Sigma=Sigma)
set.seed(1001)
e <- m0(mu,Sigma); 
x2 <- m1(1000,mu=mu,eL=e,fancy=TRUE)
identical(x1,x2)

system.time(replicate(1000,mvrnorm(1000,mu=mu,Sigma=Sigma)))
system.time({e <- m0(mu,Sigma); replicate(1000,m1(1000,mu=mu,eL=e))})
system.time({e <- m0(mu,Sigma); replicate(1000,m1(1000,mu=mu,eL=e,fancy=FALSE))})




paralist <- list(
  Nitro  = seq(0,140,length=5),
  nRep   = 3,
  nRow   = 40,
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
nRange <- 2
nRow   <- paralist$nRow


##### simulation study data generating functions #######
paras <- paralist

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

Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))

SigmaAll <- kronecker(Sigma_S,Sigma_u)


system.time({
  replicate(10,mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll))
})

system.time({
  e <- m0(rep(c(0,0),N),SigmaAll); 
  replicate(10,m1(1,mu=rep(c(0,0),N),eL=e))
})



set.seed(1001)
x1 <- mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)

set.seed(1001)
e <-  m0(rep(c(0,0),N),SigmaAll); 
x2 <- m1(1,mu= rep(c(0,0),N),eL=e,fancy=TRUE)

set.seed(1001)
x3<- fastMN(1, mu= rep(c(0,0),N), sigma = SigmaAll)

identical(x1,x3)


# For mvrnorm
start_time <- Sys.time()
for(i in 1:10){
  set.seed(1001)
  x1 <- mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)
}
end_time <- Sys.time()
mvrnorm_time <- end_time - start_time
print(paste("Time taken by mvrnorm: ", mvrnorm_time))

# For rmvnorm
start_time <- Sys.time()
for(i in 1:10){
  set.seed(1001)
  x3<- fastMN(1, mu= rep(c(0,0),N), sigma = SigmaAll)
}
end_time <- Sys.time()
rmvnorm_time <- end_time - start_time
print(paste("Time taken by rmvnorm: ", rmvnorm_time))


identical(x1,x3)


numCores <- detectCores()
cl <- makeCluster(numCores)
simu_func <- function(i) {
  AR1CorMat <- function(rho,n) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  paralist <- list(
    Nitro  = seq(0,140,length=5),
    nRep   = 3,
    nRow   = 40,
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
  nRange <- 2
  nRow   <- paralist$nRow
  
  eta = 1
  ##### simulation study data generating functions #######
  paras <- paralist
  
  nitro  <- sample(paras$Nitro)
  nRep   <- paras$nRep
  nRange <- nRep * length(nitro)
  nRow   <- paras$nRow
  N      <- nRange * nRow
  
  SystLevel <- rep(nitro,nRep)
  RandLevel <- c(replicate(nRep,sample(nitro,length(nitro))))
  
  b_0 <- paras$b_0
  b_1 <- paras$b_1
  
  LKJMat <- rethinking::rlkjcorr(1,2,eta)
  sig_u0 <- paras$sig_u0
  sig_u1 <- paras$sig_u1
  
  Sigma_u <- diag(c(sig_u0,sig_u1)) %*% LKJMat %*% diag(c(sig_u0,sig_u1))
  
  Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))
  
  SigmaAll <- kronecker(Sigma_S,Sigma_u)
  
  
  set.seed(1001)
  x1 <- MASS::mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)
  return(x1)
}

start_time <- Sys.time()
result <- parLapply(cl, 1:10, simu_func)
end_time <- Sys.time()
parrellel_time <- end_time - start_time
print(paste("Time taken by parrellel: ", parrellel_time))


identical(x1,x3,result[[1]],result[[2]])
identical(result[[1]],result[[2]])





Us_mat <- matrix(Us_mat,2,N)
u_0 <- Us_mat[1,]
u_1 <- Us_mat[2,]

LinearDat <- data.frame(Range = rep(1:nRange,each=nRow),
                        Row   = rep(1:nRow,nRange),
                        SystNitro = rep(SystLevel,each=nRow),
                        RandNitro = rep(RandLevel,each=nRow),
                        Rep   = rep(1:nRep,each = length(nitro)*nRow),
                        gridId = 1:N,
                        b0 = b_0,
                        b1 = b_1,
                        u0 = u_0,
                        u1 = u_1)

## systematic, linear y = x * a + b + e
LinearDat$SYield<- LinearDat$b0 + LinearDat$u0 + (LinearDat$b1 + LinearDat$u1) * LinearDat$SystNitro + rnorm(N)
LinearDat$RYield<- LinearDat$b0 + LinearDat$u0 + (LinearDat$b1 + LinearDat$u1) * LinearDat$RandNitro + rnorm(N)




######### parallel computing ###############

library(parallel)

######## AR1 x AR1 ############

numCores <- detectCores()
cl <- makeCluster(numCores)
ar1_simu_func <- function(k) {
  AR1CorMat <- function(rho,n) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
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
    
    LKJMat <- rethinking::rlkjcorr(1,2,eta)
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
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,2,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    
    LinearDat <- data.frame(Range = rep(1:nRange,each=nRow),
                            Row   = rep(1:nRow,nRange),
                            SystNitro = rep(SystLevel,each=nRow),
                            RandNitro = rep(RandLevel,each=nRow),
                            # Rep   = rep(1:4,each = 5*nRow),
                            Rep   = rep(1:nRep,each = length(nitro)*nRow),
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
    
    LKJMat <- rethinking::rlkjcorr(1,3,eta)
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
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,3,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    u_2 <- Us_mat[3,]
    
    QuadraDat <-data.frame(Range = rep(1:nRange,each=nRow),
                           Row = rep(1:nRow,nRange),
                           SystNitro = rep(SystLevel,each=nRow),
                           RandNitro = rep(RandLevel,each=nRow),
                           # Rep = rep(1:4,each = 5*nRow),
                           Rep = rep(1:nRep,each = length(nitro)*nRow),
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
  nRange <- paralist$nRep * length(paralist$Nitro)
  nRow   <- paralist$nRow
  
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
  
  eta <- 1
  
  for(i in 1:M){
    
    
    LinearDat <- SimuGenDat_Lin(paralist,"AR1",eta)  
    sp_ar1_LinearDat <- sp::SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                               LinearDat)
    gwr_ar1_syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                               data = sp_ar1_LinearDat, bw=5,
                               kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_LS5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst5$SDF$Intercept
    MSE_ar1_LS5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst5$SDF$SystNitro
    
    
    gwr_ar1_syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                               data = sp_ar1_LinearDat, bw=9,
                               kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_LS9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst9$SDF$Intercept
    MSE_ar1_LS9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst9$SDF$SystNitro
    
    bwLinSyst <- GWmodel::bw.gwr(SYield ~ SystNitro, data = sp_ar1_LinearDat,
                        kernel = "gaussian",
                        approach = "AICc",
                        adaptive = F)
    
    gwr_ar1_syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                 data = sp_ar1_LinearDat, bw=bwLinSyst,
                                 kernel = "gaussian")  ## boxcar kernel
    
    BW_ar1_LS[i] <- bwLinSyst
    
    MSE_ar1_LSbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst.bw$SDF$Intercept
    MSE_ar1_LSbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst.bw$SDF$SystNitro
    
    gwr_ar1_rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                               data = sp_ar1_LinearDat, bw = 5,
                               kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_LR5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand5$SDF$Intercept
    MSE_ar1_LR5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand5$SDF$RandNitro
    
    gwr_ar1_rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                               data = sp_ar1_LinearDat, bw=9, ## bw = 5
                               kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_LR9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand9$SDF$Intercept
    MSE_ar1_LR9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand9$SDF$RandNitro
    
    bwLinRand <- GWmodel::bw.gwr(RYield ~ RandNitro, data = sp_ar1_LinearDat,
                        kernel = "gaussian",
                        approach = "AICc",
                        adaptive = F)
    
    gwr_ar1_rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                 data = sp_ar1_LinearDat, bw=bwLinRand, ## bw = 5
                                 kernel = "gaussian")  ## boxcar kernel
    
    BW_ar1_LR[i] <- bwLinRand
    
    MSE_ar1_LRbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand.bw$SDF$Intercept
    MSE_ar1_LRbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand.bw$SDF$RandNitro
    
    cat(i)
    
    QuadraDat <- SimuGenDat_Qua(paralist,"AR1",eta)
    
    sp_ar1_QuadraDat <- sp::SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
    
    gwr_ar1_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                   data = sp_ar1_QuadraDat, bw=5,
                                   kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_QS5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst5$SDF$Intercept
    MSE_ar1_QS5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst5$SDF$SystNitro
    MSE_ar1_QS5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst5$SDF$SystNitro2
    
    
    gwr_ar1_qua.syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                   data = sp_ar1_QuadraDat, bw=9,
                                   kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_QS9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst9$SDF$Intercept
    MSE_ar1_QS9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst9$SDF$SystNitro
    MSE_ar1_QS9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst9$SDF$SystNitro2
    
    
    quad.syst.bw <- GWmodel::bw.gwr(SYield ~ SystNitro + SystNitro2, 
                           data = sp_ar1_QuadraDat,
                           kernel = "gaussian",
                           approach = "AICc",
                           adaptive = F)
    
    gwr_ar1_qua.syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                     data = sp_ar1_QuadraDat, bw=quad.syst.bw,
                                     kernel = "gaussian")  ## boxcar kernel
    BW_ar1_QS[i] <- quad.syst.bw
    
    MSE_ar1_QSbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst.bw$SDF$Intercept
    MSE_ar1_QSbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst.bw$SDF$SystNitro
    MSE_ar1_QSbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst.bw$SDF$SystNitro2
    
    
    gwr_ar1_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                   data = sp_ar1_QuadraDat, bw = 5, 
                                   kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_QR5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand5$SDF$Intercept
    MSE_ar1_QR5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand5$SDF$RandNitro
    MSE_ar1_QR5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand5$SDF$RandNitro2
    
    
    gwr_ar1_qua.rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                   data = sp_ar1_QuadraDat, bw = 9, 
                                   kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_QR9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand9$SDF$Intercept
    MSE_ar1_QR9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand9$SDF$RandNitro
    MSE_ar1_QR9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand9$SDF$RandNitro2
    
    
    quad.rand.bw <- GWmodel::bw.gwr(RYield ~ RandNitro + RandNitro2, 
                           data = sp_ar1_QuadraDat,
                           kernel = "gaussian",
                           approach = "AICc",
                           adaptive = F)
    
    gwr_ar1_qua.rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                     data = sp_ar1_QuadraDat, bw = quad.rand.bw, 
                                     kernel = "gaussian")  ## boxcar kernel
    
    BW_ar1_QR[i] <- quad.rand.bw
    
    MSE_ar1_QRbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand.bw$SDF$Intercept
    MSE_ar1_QRbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand.bw$SDF$RandNitro
    MSE_ar1_QRbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand.bw$SDF$RandNitro2
    
    cat(i)
  }
  return(list(MSE_ar1_LS5 = MSE_ar1_LS5, MSE_ar1_LS9 = MSE_ar1_LS9, MSE_ar1_LSbw = MSE_ar1_LSbw, 
              BW_ar1_LS = BW_ar1_LS, MSE_ar1_LR5 = MSE_ar1_LR5, MSE_ar1_LR9 = MSE_ar1_LR9, 
              MSE_ar1_LRbw = MSE_ar1_LRbw, BW_ar1_LR = BW_ar1_LR, MSE_ar1_QS5 = MSE_ar1_QS5, 
              MSE_ar1_QS9 = MSE_ar1_QS9, MSE_ar1_QSbw = MSE_ar1_QSbw, BW_ar1_QS = BW_ar1_QS, 
              MSE_ar1_QR5 = MSE_ar1_QR5, MSE_ar1_QR9 = MSE_ar1_QR9, MSE_ar1_QRbw = MSE_ar1_QRbw, 
              BW_ar1_QR = BW_ar1_QR))
  
}

parallel_AR1_results <- parLapply(cl, 1:10, ar1_simu_func)

saveRDS(parallel_AR1_results,'parallel_AR1_results_lap.rds')


######## Matern ############

numCores <- detectCores()
cl <- makeCluster(numCores)
set.seed(29183)
Mat_simu_func <- function(j) {
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
    
    LKJMat <- rethinking::rlkjcorr(1,2,eta)
    sig_u0 <- paras$sig_u0
    sig_u1 <- paras$sig_u1
    
    Sigma_u <- diag(c(sig_u0,sig_u1)) %*% LKJMat %*% diag(c(sig_u0,sig_u1))
    
    if(covmat == "NS"){
      Sigma_S = diag(N)
    }else if(covmat == "AR1"){
      Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))
    }else if(covmat == "Matern"){
      Sigma_S = SpatialTools::cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                       sp.type = "matern", 
                       sp.par = c(1,1),
                       smoothness = 1.5)$V
    }
    
    SigmaAll <- kronecker(Sigma_S,Sigma_u)
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,2,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    
    LinearDat <- data.frame(Range = rep(1:nRange,each=nRow),
                            Row   = rep(1:nRow,nRange),
                            SystNitro = rep(SystLevel,each=nRow),
                            RandNitro = rep(RandLevel,each=nRow),
                            # Rep   = rep(1:4,each = 5*nRow),
                            Rep   = rep(1:nRep,each = length(nitro)*nRow),
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
    
    LKJMat <- rethinking::rlkjcorr(1,3,eta)
    sig_u0 <- paras$sig_u0
    sig_u1 <- paras$sig_u1
    sig_u2 <- paras$sig_u2
    
    Sigma_u <- diag(c(sig_u0,sig_u1,sig_u2)) %*% LKJMat %*% diag(c(sig_u0,sig_u1,sig_u2))
    
    if(covmat == "NS"){
      Sigma_S = diag(N)
    }else if(covmat == "AR1"){
      Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))
    }else if(covmat == "Matern"){
      Sigma_S = SpatialTools::cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                       sp.type = "matern", 
                       sp.par = c(1,1),
                       smoothness = 1.5)$V
    }
    
    SigmaAll <- kronecker(Sigma_S,Sigma_u)
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,3,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    u_2 <- Us_mat[3,]
    
    QuadraDat <-data.frame(Range = rep(1:nRange,each=nRow),
                           Row = rep(1:nRow,nRange),
                           SystNitro = rep(SystLevel,each=nRow),
                           RandNitro = rep(RandLevel,each=nRow),
                           # Rep = rep(1:4,each = 5*nRow),
                           Rep = rep(1:nRep,each = length(nitro)*nRow),
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
  
  M <- 200 ## simulation iterations
  
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
  nRange <- paralist$nRep * length(paralist$Nitro)
  nRow   <- paralist$nRow
  
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
  
  eta <- 1
  
  for(i in 1:M){
    
    LinearDat <- SimuGenDat_Lin(paralist,"Matern",eta)  
    sp_mat_LinearDat <- sp::SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                                   LinearDat)
    gwr_mat_syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                        data = sp_mat_LinearDat, bw=5,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_LS5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst5$SDF$Intercept
    MSE_mat_LS5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst5$SDF$SystNitro
    
    
    gwr_mat_syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                        data = sp_mat_LinearDat, bw=9,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_LS9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst9$SDF$Intercept
    MSE_mat_LS9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst9$SDF$SystNitro
    
    bwLinSyst <- GWmodel::bw.gwr(SYield ~ SystNitro, data = sp_mat_LinearDat,
                                 kernel = "gaussian",
                                 approach = "AICc",
                                 adaptive = F)
    
    gwr_mat_syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                          data = sp_mat_LinearDat, bw=bwLinSyst,
                                          kernel = "gaussian")  ## boxcar kernel
    
    BW_mat_LS[i] <- bwLinSyst
    
    MSE_mat_LSbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst.bw$SDF$Intercept
    MSE_mat_LSbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst.bw$SDF$SystNitro
    
    gwr_mat_rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                        data = sp_mat_LinearDat, bw = 5,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_LR5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand5$SDF$Intercept
    MSE_mat_LR5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand5$SDF$RandNitro
    
    gwr_mat_rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                        data = sp_mat_LinearDat, bw=9, ## bw = 5
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_LR9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand9$SDF$Intercept
    MSE_mat_LR9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand9$SDF$RandNitro
    
    bwLinRand <- GWmodel::bw.gwr(RYield ~ RandNitro, data = sp_mat_LinearDat,
                                 kernel = "gaussian",
                                 approach = "AICc",
                                 adaptive = F)
    
    gwr_mat_rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                          data = sp_mat_LinearDat, bw=bwLinRand, ## bw = 5
                                          kernel = "gaussian")  ## boxcar kernel
    
    BW_mat_LR[i] <- bwLinRand
    
    MSE_mat_LRbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand.bw$SDF$Intercept
    MSE_mat_LRbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand.bw$SDF$RandNitro
    
    cat(i)
    
    QuadraDat <- SimuGenDat_Qua(paralist,"Matern",eta)
    
    sp_mat_QuadraDat <- sp::SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
    
    gwr_mat_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_mat_QuadraDat, bw=5,
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_QS5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst5$SDF$Intercept
    MSE_mat_QS5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst5$SDF$SystNitro
    MSE_mat_QS5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst5$SDF$SystNitro2
    
    
    gwr_mat_qua.syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_mat_QuadraDat, bw=9,
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_QS9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst9$SDF$Intercept
    MSE_mat_QS9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst9$SDF$SystNitro
    MSE_mat_QS9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst9$SDF$SystNitro2
    
    
    quad.syst.bw <- GWmodel::bw.gwr(SYield ~ SystNitro + SystNitro2, 
                                    data = sp_mat_QuadraDat,
                                    kernel = "gaussian",
                                    approach = "AICc",
                                    adaptive = F)
    
    gwr_mat_qua.syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                              data = sp_mat_QuadraDat, bw=quad.syst.bw,
                                              kernel = "gaussian")  ## boxcar kernel
    BW_mat_QS[i] <- quad.syst.bw
    
    MSE_mat_QSbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst.bw$SDF$Intercept
    MSE_mat_QSbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst.bw$SDF$SystNitro
    MSE_mat_QSbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst.bw$SDF$SystNitro2
    
    
    gwr_mat_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_mat_QuadraDat, bw = 5, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_QR5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand5$SDF$Intercept
    MSE_mat_QR5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand5$SDF$RandNitro
    MSE_mat_QR5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand5$SDF$RandNitro2
    
    
    gwr_mat_qua.rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_mat_QuadraDat, bw = 9, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_QR9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand9$SDF$Intercept
    MSE_mat_QR9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand9$SDF$RandNitro
    MSE_mat_QR9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand9$SDF$RandNitro2
    
    
    quad.rand.bw <- GWmodel::bw.gwr(RYield ~ RandNitro + RandNitro2, 
                                    data = sp_mat_QuadraDat,
                                    kernel = "gaussian",
                                    approach = "AICc",
                                    adaptive = F)
    
    gwr_mat_qua.rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                              data = sp_mat_QuadraDat, bw = quad.rand.bw, 
                                              kernel = "gaussian")  ## boxcar kernel
    
    BW_mat_QR[i] <- quad.rand.bw
    
    MSE_mat_QRbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand.bw$SDF$Intercept
    MSE_mat_QRbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand.bw$SDF$RandNitro
    MSE_mat_QRbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand.bw$SDF$RandNitro2
    
    cat(i)
  }
  
  return(list(MSE_mat_LS5 = MSE_mat_LS5, MSE_mat_LS9 = MSE_mat_LS9, MSE_mat_LSbw = MSE_mat_LSbw, 
              BW_mat_LS = BW_mat_LS, MSE_mat_LR5 = MSE_mat_LR5, MSE_mat_LR9 = MSE_mat_LR9, 
              MSE_mat_LRbw = MSE_mat_LRbw, BW_mat_LR = BW_mat_LR, MSE_mat_QS5 = MSE_mat_QS5, 
              MSE_mat_QS9 = MSE_mat_QS9, MSE_mat_QSbw = MSE_mat_QSbw, BW_mat_QS = BW_mat_QS, 
              MSE_mat_QR5 = MSE_mat_QR5, MSE_mat_QR9 = MSE_mat_QR9, MSE_mat_QRbw = MSE_mat_QRbw, 
              BW_mat_QR = BW_mat_QR))
  
}
st<- Sys.time()
parallel_mat_results <- parLapply(cl, 1:5, Mat_simu_func)
ed<- Sys.time()

ed-st

# saveRDS(parallel_mat_results,"parallel_mat_results_mac.rds")
# saveRDS(parallel_mat_results,"parallel_mat_results_pawsey.rds")





#################  eta = 0.1 ###############


library(parallel)
############## NS #####################

numCores <- detectCores()
cl <- makeCluster(numCores)
ns_simu_func <- function(k) {
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
    
    LKJMat <- rethinking::rlkjcorr(1,2,eta)
    sig_u0 <- paras$sig_u0
    sig_u1 <- paras$sig_u1
    
    Sigma_u <- diag(c(sig_u0,sig_u1)) %*% LKJMat %*% diag(c(sig_u0,sig_u1))
    
    if(covmat == "NS"){
      Sigma_S = diag(N)
    }else if(covmat == "NS"){
      Sigma_S = kronecker(NSCorMat(paras$rho_range,nRange),NSCorMat(paras$rho_row,nRow))
    }else if(covmat == "Matern"){
      Sigma_S = cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                       sp.type = "matern", 
                       sp.par = c(1,1),
                       smoothness = 1.5)$V
    }
    
    SigmaAll <- kronecker(Sigma_S,Sigma_u)
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,2,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    
    LinearDat <- data.frame(Range = rep(1:nRange,each=nRow),
                            Row   = rep(1:nRow,nRange),
                            SystNitro = rep(SystLevel,each=nRow),
                            RandNitro = rep(RandLevel,each=nRow),
                            # Rep   = rep(1:4,each = 5*nRow),
                            Rep   = rep(1:nRep,each = length(nitro)*nRow),
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
    
    LKJMat <- rethinking::rlkjcorr(1,3,eta)
    sig_u0 <- paras$sig_u0
    sig_u1 <- paras$sig_u1
    sig_u2 <- paras$sig_u2
    
    Sigma_u <- diag(c(sig_u0,sig_u1,sig_u2)) %*% LKJMat %*% diag(c(sig_u0,sig_u1,sig_u2))
    
    if(covmat == "NS"){
      Sigma_S = diag(N)
    }else if(covmat == "NS"){
      Sigma_S = kronecker(NSCorMat(paras$rho_range,nRange),NSCorMat(paras$rho_row,nRow))
    }else if(covmat == "Matern"){
      Sigma_S = cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                       sp.type = "matern", 
                       sp.par = c(1,1),
                       smoothness = 1.5)$V
    }
    
    SigmaAll <- kronecker(Sigma_S,Sigma_u)
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,3,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    u_2 <- Us_mat[3,]
    
    QuadraDat <-data.frame(Range = rep(1:nRange,each=nRow),
                           Row = rep(1:nRow,nRange),
                           SystNitro = rep(SystLevel,each=nRow),
                           RandNitro = rep(RandLevel,each=nRow),
                           # Rep = rep(1:4,each = 5*nRow),
                           Rep = rep(1:nRep,each = length(nitro)*nRow),
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
  nRange <- paralist$nRep * length(paralist$Nitro)
  nRow   <- paralist$nRow
  
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
  
  eta <- 0.1
  
  for(i in 1:M){
    
    
    LinearDat <- SimuGenDat_Lin(paralist,"NS",eta)  
    sp_ns_LinearDat <- sp::SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                                   LinearDat)
    gwr_ns_syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                        data = sp_ns_LinearDat, bw=5,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_ns_LS5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ns_syst5$SDF$Intercept
    MSE_ns_LS5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ns_syst5$SDF$SystNitro
    
    
    gwr_ns_syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                        data = sp_ns_LinearDat, bw=9,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_ns_LS9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ns_syst9$SDF$Intercept
    MSE_ns_LS9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ns_syst9$SDF$SystNitro
    
    bwLinSyst <- GWmodel::bw.gwr(SYield ~ SystNitro, data = sp_ns_LinearDat,
                                 kernel = "gaussian",
                                 approach = "AICc",
                                 adaptive = F)
    
    gwr_ns_syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                          data = sp_ns_LinearDat, bw=bwLinSyst,
                                          kernel = "gaussian")  ## boxcar kernel
    
    BW_ns_LS[i] <- bwLinSyst
    
    MSE_ns_LSbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ns_syst.bw$SDF$Intercept
    MSE_ns_LSbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ns_syst.bw$SDF$SystNitro
    
    gwr_ns_rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                        data = sp_ns_LinearDat, bw = 5,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_ns_LR5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ns_rand5$SDF$Intercept
    MSE_ns_LR5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ns_rand5$SDF$RandNitro
    
    gwr_ns_rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                        data = sp_ns_LinearDat, bw=9, ## bw = 5
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_ns_LR9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ns_rand9$SDF$Intercept
    MSE_ns_LR9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ns_rand9$SDF$RandNitro
    
    bwLinRand <- GWmodel::bw.gwr(RYield ~ RandNitro, data = sp_ns_LinearDat,
                                 kernel = "gaussian",
                                 approach = "AICc",
                                 adaptive = F)
    
    gwr_ns_rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                          data = sp_ns_LinearDat, bw=bwLinRand, ## bw = 5
                                          kernel = "gaussian")  ## boxcar kernel
    
    BW_ns_LR[i] <- bwLinRand
    
    MSE_ns_LRbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ns_rand.bw$SDF$Intercept
    MSE_ns_LRbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ns_rand.bw$SDF$RandNitro
    
    cat(i)
    
    QuadraDat <- SimuGenDat_Qua(paralist,"NS",eta)
    
    sp_ns_QuadraDat <- sp::SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
    
    gwr_ns_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_ns_QuadraDat, bw=5,
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_ns_QS5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ns_qua.syst5$SDF$Intercept
    MSE_ns_QS5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ns_qua.syst5$SDF$SystNitro
    MSE_ns_QS5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ns_qua.syst5$SDF$SystNitro2
    
    
    gwr_ns_qua.syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_ns_QuadraDat, bw=9,
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_ns_QS9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ns_qua.syst9$SDF$Intercept
    MSE_ns_QS9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ns_qua.syst9$SDF$SystNitro
    MSE_ns_QS9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ns_qua.syst9$SDF$SystNitro2
    
    
    quad.syst.bw <- GWmodel::bw.gwr(SYield ~ SystNitro + SystNitro2, 
                                    data = sp_ns_QuadraDat,
                                    kernel = "gaussian",
                                    approach = "AICc",
                                    adaptive = F)
    
    gwr_ns_qua.syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                              data = sp_ns_QuadraDat, bw=quad.syst.bw,
                                              kernel = "gaussian")  ## boxcar kernel
    BW_ns_QS[i] <- quad.syst.bw
    
    MSE_ns_QSbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ns_qua.syst.bw$SDF$Intercept
    MSE_ns_QSbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ns_qua.syst.bw$SDF$SystNitro
    MSE_ns_QSbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ns_qua.syst.bw$SDF$SystNitro2
    
    
    gwr_ns_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_ns_QuadraDat, bw = 5, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_ns_QR5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ns_qua.rand5$SDF$Intercept
    MSE_ns_QR5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ns_qua.rand5$SDF$RandNitro
    MSE_ns_QR5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ns_qua.rand5$SDF$RandNitro2
    
    
    gwr_ns_qua.rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_ns_QuadraDat, bw = 9, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_ns_QR9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ns_qua.rand9$SDF$Intercept
    MSE_ns_QR9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ns_qua.rand9$SDF$RandNitro
    MSE_ns_QR9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ns_qua.rand9$SDF$RandNitro2
    
    
    quad.rand.bw <- GWmodel::bw.gwr(RYield ~ RandNitro + RandNitro2, 
                                    data = sp_ns_QuadraDat,
                                    kernel = "gaussian",
                                    approach = "AICc",
                                    adaptive = F)
    
    gwr_ns_qua.rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                              data = sp_ns_QuadraDat, bw = quad.rand.bw, 
                                              kernel = "gaussian")  ## boxcar kernel
    
    BW_ns_QR[i] <- quad.rand.bw
    
    MSE_ns_QRbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ns_qua.rand.bw$SDF$Intercept
    MSE_ns_QRbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ns_qua.rand.bw$SDF$RandNitro
    MSE_ns_QRbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ns_qua.rand.bw$SDF$RandNitro2
    
    cat(i)
  }
  return(list(MSE_ns_LS5 = MSE_ns_LS5, MSE_ns_LS9 = MSE_ns_LS9, MSE_ns_LSbw = MSE_ns_LSbw, 
              BW_ns_LS = BW_ns_LS, MSE_ns_LR5 = MSE_ns_LR5, MSE_ns_LR9 = MSE_ns_LR9, 
              MSE_ns_LRbw = MSE_ns_LRbw, BW_ns_LR = BW_ns_LR, MSE_ns_QS5 = MSE_ns_QS5, 
              MSE_ns_QS9 = MSE_ns_QS9, MSE_ns_QSbw = MSE_ns_QSbw, BW_ns_QS = BW_ns_QS, 
              MSE_ns_QR5 = MSE_ns_QR5, MSE_ns_QR9 = MSE_ns_QR9, MSE_ns_QRbw = MSE_ns_QRbw, 
              BW_ns_QR = BW_ns_QR))
  
}
parallel_NS_results_eta01 <- parLapply(cl, 1:10, ns_simu_func)

saveRDS(parallel_NS_results_eta01,'parallel_NS_results_lap_eta01.rds')



######## AR1 x AR1 ############

numCores <- detectCores()
cl <- makeCluster(numCores)

ar1_simu_func <- function(k) {
  AR1CorMat <- function(rho,n) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
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
    
    LKJMat <- rethinking::rlkjcorr(1,2,eta)
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
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,2,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    
    LinearDat <- data.frame(Range = rep(1:nRange,each=nRow),
                            Row   = rep(1:nRow,nRange),
                            SystNitro = rep(SystLevel,each=nRow),
                            RandNitro = rep(RandLevel,each=nRow),
                            # Rep   = rep(1:4,each = 5*nRow),
                            Rep   = rep(1:nRep,each = length(nitro)*nRow),
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
    
    LKJMat <- rethinking::rlkjcorr(1,3,eta)
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
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,3,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    u_2 <- Us_mat[3,]
    
    QuadraDat <-data.frame(Range = rep(1:nRange,each=nRow),
                           Row = rep(1:nRow,nRange),
                           SystNitro = rep(SystLevel,each=nRow),
                           RandNitro = rep(RandLevel,each=nRow),
                           # Rep = rep(1:4,each = 5*nRow),
                           Rep = rep(1:nRep,each = length(nitro)*nRow),
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
  nRange <- paralist$nRep * length(paralist$Nitro)
  nRow   <- paralist$nRow
  
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
  
  eta <- 0.1
  
  for(i in 1:M){
    LinearDat <- SimuGenDat_Lin(paralist,"AR1",eta)  
    sp_ar1_LinearDat <- sp::SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                                   LinearDat)
    gwr_ar1_syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                        data = sp_ar1_LinearDat, bw=5,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_LS5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst5$SDF$Intercept
    MSE_ar1_LS5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst5$SDF$SystNitro
    
    
    gwr_ar1_syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                        data = sp_ar1_LinearDat, bw=9,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_LS9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst9$SDF$Intercept
    MSE_ar1_LS9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst9$SDF$SystNitro
    
    bwLinSyst <- GWmodel::bw.gwr(SYield ~ SystNitro, data = sp_ar1_LinearDat,
                                 kernel = "gaussian",
                                 approach = "AICc",
                                 adaptive = F)
    
    gwr_ar1_syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                          data = sp_ar1_LinearDat, bw=bwLinSyst,
                                          kernel = "gaussian")  ## boxcar kernel
    
    BW_ar1_LS[i] <- bwLinSyst
    
    MSE_ar1_LSbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_syst.bw$SDF$Intercept
    MSE_ar1_LSbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_syst.bw$SDF$SystNitro
    
    gwr_ar1_rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                        data = sp_ar1_LinearDat, bw = 5,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_LR5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand5$SDF$Intercept
    MSE_ar1_LR5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand5$SDF$RandNitro
    
    gwr_ar1_rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                        data = sp_ar1_LinearDat, bw=9, ## bw = 5
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_LR9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand9$SDF$Intercept
    MSE_ar1_LR9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand9$SDF$RandNitro
    
    bwLinRand <- GWmodel::bw.gwr(RYield ~ RandNitro, data = sp_ar1_LinearDat,
                                 kernel = "gaussian",
                                 approach = "AICc",
                                 adaptive = F)
    
    gwr_ar1_rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                          data = sp_ar1_LinearDat, bw=bwLinRand, ## bw = 5
                                          kernel = "gaussian")  ## boxcar kernel
    
    BW_ar1_LR[i] <- bwLinRand
    
    MSE_ar1_LRbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_ar1_rand.bw$SDF$Intercept
    MSE_ar1_LRbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_ar1_rand.bw$SDF$RandNitro
    
    cat(i)
    
    QuadraDat <- SimuGenDat_Qua(paralist,"AR1",eta)
    
    sp_ar1_QuadraDat <- sp::SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
    
    gwr_ar1_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_ar1_QuadraDat, bw=5,
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_QS5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst5$SDF$Intercept
    MSE_ar1_QS5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst5$SDF$SystNitro
    MSE_ar1_QS5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst5$SDF$SystNitro2
    
    
    gwr_ar1_qua.syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_ar1_QuadraDat, bw=9,
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_QS9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst9$SDF$Intercept
    MSE_ar1_QS9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst9$SDF$SystNitro
    MSE_ar1_QS9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst9$SDF$SystNitro2
    
    
    quad.syst.bw <- GWmodel::bw.gwr(SYield ~ SystNitro + SystNitro2, 
                                    data = sp_ar1_QuadraDat,
                                    kernel = "gaussian",
                                    approach = "AICc",
                                    adaptive = F)
    
    gwr_ar1_qua.syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                              data = sp_ar1_QuadraDat, bw=quad.syst.bw,
                                              kernel = "gaussian")  ## boxcar kernel
    BW_ar1_QS[i] <- quad.syst.bw
    
    MSE_ar1_QSbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.syst.bw$SDF$Intercept
    MSE_ar1_QSbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.syst.bw$SDF$SystNitro
    MSE_ar1_QSbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.syst.bw$SDF$SystNitro2
    
    
    gwr_ar1_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_ar1_QuadraDat, bw = 5, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_QR5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand5$SDF$Intercept
    MSE_ar1_QR5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand5$SDF$RandNitro
    MSE_ar1_QR5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand5$SDF$RandNitro2
    
    
    gwr_ar1_qua.rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_ar1_QuadraDat, bw = 9, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_ar1_QR9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand9$SDF$Intercept
    MSE_ar1_QR9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand9$SDF$RandNitro
    MSE_ar1_QR9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand9$SDF$RandNitro2
    
    
    quad.rand.bw <- GWmodel::bw.gwr(RYield ~ RandNitro + RandNitro2, 
                                    data = sp_ar1_QuadraDat,
                                    kernel = "gaussian",
                                    approach = "AICc",
                                    adaptive = F)
    
    gwr_ar1_qua.rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                              data = sp_ar1_QuadraDat, bw = quad.rand.bw, 
                                              kernel = "gaussian")  ## boxcar kernel
    
    BW_ar1_QR[i] <- quad.rand.bw
    
    MSE_ar1_QRbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_ar1_qua.rand.bw$SDF$Intercept
    MSE_ar1_QRbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_ar1_qua.rand.bw$SDF$RandNitro
    MSE_ar1_QRbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_ar1_qua.rand.bw$SDF$RandNitro2
    
    cat(i)
  }
  return(list(MSE_ar1_LS5 = MSE_ar1_LS5, MSE_ar1_LS9 = MSE_ar1_LS9, MSE_ar1_LSbw = MSE_ar1_LSbw, 
              BW_ar1_LS = BW_ar1_LS, MSE_ar1_LR5 = MSE_ar1_LR5, MSE_ar1_LR9 = MSE_ar1_LR9, 
              MSE_ar1_LRbw = MSE_ar1_LRbw, BW_ar1_LR = BW_ar1_LR, MSE_ar1_QS5 = MSE_ar1_QS5, 
              MSE_ar1_QS9 = MSE_ar1_QS9, MSE_ar1_QSbw = MSE_ar1_QSbw, BW_ar1_QS = BW_ar1_QS, 
              MSE_ar1_QR5 = MSE_ar1_QR5, MSE_ar1_QR9 = MSE_ar1_QR9, MSE_ar1_QRbw = MSE_ar1_QRbw, 
              BW_ar1_QR = BW_ar1_QR))
  
}
parallel_AR1_results_eta01 <- parLapply(cl, 1:10, ar1_simu_func)

saveRDS(parallel_AR1_results_eta01,'parallel_AR1_results_lap_eta01.rds')


######## Matern, eta = 0.1 ############

numCores <- detectCores()
cl <- makeCluster(numCores)
set.seed(48123)
Mat_simu_func <- function(j) {
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
    
    LKJMat <- rethinking::rlkjcorr(1,2,eta)
    sig_u0 <- paras$sig_u0
    sig_u1 <- paras$sig_u1
    
    Sigma_u <- diag(c(sig_u0,sig_u1)) %*% LKJMat %*% diag(c(sig_u0,sig_u1))
    
    if(covmat == "NS"){
      Sigma_S = diag(N)
    }else if(covmat == "AR1"){
      Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))
    }else if(covmat == "Matern"){
      Sigma_S = SpatialTools::cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                                     sp.type = "matern", 
                                     sp.par = c(1,1),
                                     smoothness = 1.5)$V
    }
    
    SigmaAll <- kronecker(Sigma_S,Sigma_u)
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,2,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    
    LinearDat <- data.frame(Range = rep(1:nRange,each=nRow),
                            Row   = rep(1:nRow,nRange),
                            SystNitro = rep(SystLevel,each=nRow),
                            RandNitro = rep(RandLevel,each=nRow),
                            # Rep   = rep(1:4,each = 5*nRow),
                            Rep   = rep(1:nRep,each = length(nitro)*nRow),
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
    
    LKJMat <- rethinking::rlkjcorr(1,3,eta)
    sig_u0 <- paras$sig_u0
    sig_u1 <- paras$sig_u1
    sig_u2 <- paras$sig_u2
    
    Sigma_u <- diag(c(sig_u0,sig_u1,sig_u2)) %*% LKJMat %*% diag(c(sig_u0,sig_u1,sig_u2))
    
    if(covmat == "NS"){
      Sigma_S = diag(N)
    }else if(covmat == "AR1"){
      Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))
    }else if(covmat == "Matern"){
      Sigma_S = SpatialTools::cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                                     sp.type = "matern", 
                                     sp.par = c(1,1),
                                     smoothness = 1.5)$V
    }
    
    SigmaAll <- kronecker(Sigma_S,Sigma_u)
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,3,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    u_2 <- Us_mat[3,]
    
    QuadraDat <-data.frame(Range = rep(1:nRange,each=nRow),
                           Row = rep(1:nRow,nRange),
                           SystNitro = rep(SystLevel,each=nRow),
                           RandNitro = rep(RandLevel,each=nRow),
                           # Rep = rep(1:4,each = 5*nRow),
                           Rep = rep(1:nRep,each = length(nitro)*nRow),
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
  
  M <- 200 ## simulation iterations
  
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
  nRange <- paralist$nRep * length(paralist$Nitro)
  nRow   <- paralist$nRow
  
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
  
  eta <- 0.1
  
  for(i in 1:M){
    
    LinearDat <- SimuGenDat_Lin(paralist,"Matern",eta)  
    sp_mat_LinearDat <- sp::SpatialPointsDataFrame(cbind(LinearDat$Range,LinearDat$Row),
                                                   LinearDat)
    gwr_mat_syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                        data = sp_mat_LinearDat, bw=5,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_LS5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst5$SDF$Intercept
    MSE_mat_LS5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst5$SDF$SystNitro
    
    
    gwr_mat_syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                        data = sp_mat_LinearDat, bw=9,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_LS9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst9$SDF$Intercept
    MSE_mat_LS9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst9$SDF$SystNitro
    
    bwLinSyst <- GWmodel::bw.gwr(SYield ~ SystNitro, data = sp_mat_LinearDat,
                                 kernel = "gaussian",
                                 approach = "AICc",
                                 adaptive = F)
    
    gwr_mat_syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro,
                                          data = sp_mat_LinearDat, bw=bwLinSyst,
                                          kernel = "gaussian")  ## boxcar kernel
    
    BW_mat_LS[i] <- bwLinSyst
    
    MSE_mat_LSbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_syst.bw$SDF$Intercept
    MSE_mat_LSbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_syst.bw$SDF$SystNitro
    
    gwr_mat_rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                        data = sp_mat_LinearDat, bw = 5,
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_LR5[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand5$SDF$Intercept
    MSE_mat_LR5[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand5$SDF$RandNitro
    
    gwr_mat_rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                        data = sp_mat_LinearDat, bw=9, ## bw = 5
                                        kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_LR9[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand9$SDF$Intercept
    MSE_mat_LR9[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand9$SDF$RandNitro
    
    bwLinRand <- GWmodel::bw.gwr(RYield ~ RandNitro, data = sp_mat_LinearDat,
                                 kernel = "gaussian",
                                 approach = "AICc",
                                 adaptive = F)
    
    gwr_mat_rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro,
                                          data = sp_mat_LinearDat, bw=bwLinRand, ## bw = 5
                                          kernel = "gaussian")  ## boxcar kernel
    
    BW_mat_LR[i] <- bwLinRand
    
    MSE_mat_LRbw[i,1,] <- (LinearDat$b0 + LinearDat$u0) - gwr_mat_rand.bw$SDF$Intercept
    MSE_mat_LRbw[i,2,] <- (LinearDat$b1 + LinearDat$u1) - gwr_mat_rand.bw$SDF$RandNitro
    
    cat(i)
    
    QuadraDat <- SimuGenDat_Qua(paralist,"Matern",eta)
    
    sp_mat_QuadraDat <- sp::SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
    
    gwr_mat_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_mat_QuadraDat, bw=5,
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_QS5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst5$SDF$Intercept
    MSE_mat_QS5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst5$SDF$SystNitro
    MSE_mat_QS5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst5$SDF$SystNitro2
    
    
    gwr_mat_qua.syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_mat_QuadraDat, bw=9,
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_QS9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst9$SDF$Intercept
    MSE_mat_QS9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst9$SDF$SystNitro
    MSE_mat_QS9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst9$SDF$SystNitro2
    
    
    quad.syst.bw <- GWmodel::bw.gwr(SYield ~ SystNitro + SystNitro2, 
                                    data = sp_mat_QuadraDat,
                                    kernel = "gaussian",
                                    approach = "AICc",
                                    adaptive = F)
    
    gwr_mat_qua.syst.bw <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                              data = sp_mat_QuadraDat, bw=quad.syst.bw,
                                              kernel = "gaussian")  ## boxcar kernel
    BW_mat_QS[i] <- quad.syst.bw
    
    MSE_mat_QSbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.syst.bw$SDF$Intercept
    MSE_mat_QSbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.syst.bw$SDF$SystNitro
    MSE_mat_QSbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.syst.bw$SDF$SystNitro2
    
    
    gwr_mat_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_mat_QuadraDat, bw = 5, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_QR5[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand5$SDF$Intercept
    MSE_mat_QR5[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand5$SDF$RandNitro
    MSE_mat_QR5[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand5$SDF$RandNitro2
    
    
    gwr_mat_qua.rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_mat_QuadraDat, bw = 9, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    MSE_mat_QR9[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand9$SDF$Intercept
    MSE_mat_QR9[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand9$SDF$RandNitro
    MSE_mat_QR9[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand9$SDF$RandNitro2
    
    
    quad.rand.bw <- GWmodel::bw.gwr(RYield ~ RandNitro + RandNitro2, 
                                    data = sp_mat_QuadraDat,
                                    kernel = "gaussian",
                                    approach = "AICc",
                                    adaptive = F)
    
    gwr_mat_qua.rand.bw <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                              data = sp_mat_QuadraDat, bw = quad.rand.bw, 
                                              kernel = "gaussian")  ## boxcar kernel
    
    BW_mat_QR[i] <- quad.rand.bw
    
    MSE_mat_QRbw[i,1,] <- (QuadraDat$b0+QuadraDat$u0) - gwr_mat_qua.rand.bw$SDF$Intercept
    MSE_mat_QRbw[i,2,] <- (QuadraDat$b1+QuadraDat$u1) - gwr_mat_qua.rand.bw$SDF$RandNitro
    MSE_mat_QRbw[i,3,] <- (QuadraDat$b2+QuadraDat$u2) - gwr_mat_qua.rand.bw$SDF$RandNitro2
    
    cat(i)
  }
  
  return(list(MSE_mat_LS5 = MSE_mat_LS5, MSE_mat_LS9 = MSE_mat_LS9, MSE_mat_LSbw = MSE_mat_LSbw, 
              BW_mat_LS = BW_mat_LS, MSE_mat_LR5 = MSE_mat_LR5, MSE_mat_LR9 = MSE_mat_LR9, 
              MSE_mat_LRbw = MSE_mat_LRbw, BW_mat_LR = BW_mat_LR, MSE_mat_QS5 = MSE_mat_QS5, 
              MSE_mat_QS9 = MSE_mat_QS9, MSE_mat_QSbw = MSE_mat_QSbw, BW_mat_QS = BW_mat_QS, 
              MSE_mat_QR5 = MSE_mat_QR5, MSE_mat_QR9 = MSE_mat_QR9, MSE_mat_QRbw = MSE_mat_QRbw, 
              BW_mat_QR = BW_mat_QR))
  
}
st<- Sys.time()
parallel_mat_results_eta01 <- parLapply(cl, 1:5, Mat_simu_func)
ed<- Sys.time()

ed-st

saveRDS(parallel_mat_results_eta01,"parallel_mat_results_mac_eta01.rds")
##saveRDS(parallel_mat_results_eta01,"parallel_mat_results_pawsey_eta01.rds")



##### parallel test #######

numCores <- detectCores()
cl <- makeCluster(numCores)
set.seed(113132)
paratest_func <- function(j) {
  a <- rnorm(3)
  return(a)
}

paratest_results <- parLapply(cl, 1:50, paratest_func)



######## for coefficients only ########

library(parallel)
numCores <- detectCores()
cl <- makeCluster(numCores)
ar1_simu_coef_func <- function(k) {
  AR1CorMat <- function(rho,n) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  ##### simulation study data generating functions #######
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
    
    LKJMat <- rethinking::rlkjcorr(1,3,eta)
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
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,3,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    u_2 <- Us_mat[3,]
    
    QuadraDat <-data.frame(Range = rep(1:nRange,each=nRow),
                           Row = rep(1:nRow,nRange),
                           SystNitro = rep(SystLevel,each=nRow),
                           RandNitro = rep(RandLevel,each=nRow),
                           # Rep = rep(1:4,each = 5*nRow),
                           Rep = rep(1:nRep,each = length(nitro)*nRow),
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
  nRange <- paralist$nRep * length(paralist$Nitro)
  nRow   <- paralist$nRow
  
  eta <- 1
  
  data_frames_list <- list()
  
  for(i in 1:M){
    
    QuadraDat <- SimuGenDat_Qua(paralist,"AR1",eta)
    
    sp_ar1_QuadraDat <- sp::SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
    
    gwr_ar1_qua.syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_ar1_QuadraDat, bw=9,
                                            kernel = "gaussian")  ## boxcar kernel
    
    QuadraDat$SystGWRB0_bw9 <- gwr_ar1_qua.syst9$SDF$Intercept
    QuadraDat$SystGWRB1_bw9 <- gwr_ar1_qua.syst9$SDF$SystNitro
    QuadraDat$SystGWRB2_bw9 <- gwr_ar1_qua.syst9$SDF$SystNitro2
    
    gwr_ar1_qua.rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_ar1_QuadraDat, bw = 9, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    QuadraDat$RandGWRB0_bw9 <- gwr_ar1_qua.rand9$SDF$Intercept
    QuadraDat$RandGWRB1_bw9 <- gwr_ar1_qua.rand9$SDF$RandNitro
    QuadraDat$RandGWRB2_bw9 <- gwr_ar1_qua.rand9$SDF$RandNitro2
    
    gwr_ar1_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_ar1_QuadraDat, bw=5,
                                            kernel = "gaussian")  ## boxcar kernel
    
    QuadraDat$SystGWRB0_bw5 <- gwr_ar1_qua.syst5$SDF$Intercept
    QuadraDat$SystGWRB1_bw5 <- gwr_ar1_qua.syst5$SDF$SystNitro
    QuadraDat$SystGWRB2_bw5 <- gwr_ar1_qua.syst5$SDF$SystNitro2
    
    gwr_ar1_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_ar1_QuadraDat, bw = 5, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    QuadraDat$RandGWRB0_bw5 <- gwr_ar1_qua.rand5$SDF$Intercept
    QuadraDat$RandGWRB1_bw5 <- gwr_ar1_qua.rand5$SDF$RandNitro
    QuadraDat$RandGWRB2_bw5 <- gwr_ar1_qua.rand5$SDF$RandNitro2
    
    data_frames_list[[i]] <- QuadraDat
  }
  return(data_frames_list)
}

## 1k
set.seed(202432)
st <- Sys.time()
parallel_AR1_coef_results <- parLapply(cl, 1:10, ar1_simu_coef_func)
saveRDS(parallel_AR1_coef_results,'parallel_AR1_coef_results_mac.rds')
ed <- Sys.time()
ed - st





## 1k
numCores <- detectCores()
cl <- makeCluster(numCores)
mat_simu_coef_func <- function(k) {
  AR1CorMat <- function(rho,n) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  ##### simulation study data generating functions #######
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
    
    LKJMat <- rethinking::rlkjcorr(1,3,eta)
    sig_u0 <- paras$sig_u0
    sig_u1 <- paras$sig_u1
    sig_u2 <- paras$sig_u2
    
    Sigma_u <- diag(c(sig_u0,sig_u1,sig_u2)) %*% LKJMat %*% diag(c(sig_u0,sig_u1,sig_u2))
    
    if(covmat == "NS"){
      Sigma_S = diag(N)
    }else if(covmat == "AR1"){
      Sigma_S = kronecker(AR1CorMat(paras$rho_range,nRange),AR1CorMat(paras$rho_row,nRow))
    }else if(covmat == "Matern"){
      Sigma_S = SpatialTools::cov.sp(coords = matrix(c(rep(1:nRange,each=nRow),rep(1:nRow,nRange)),ncol=2), 
                                     sp.type = "matern", 
                                     sp.par = c(1,1),
                                     smoothness = 1.5)$V
    }
    
    SigmaAll <- kronecker(Sigma_S,Sigma_u)
    
    Us_mat <- MASS::mvrnorm(1, mu= rep(c(0,0,0),N), Sigma = SigmaAll)
    Us_mat <- matrix(Us_mat,3,N)
    u_0 <- Us_mat[1,]
    u_1 <- Us_mat[2,]
    u_2 <- Us_mat[3,]
    
    QuadraDat <-data.frame(Range = rep(1:nRange,each=nRow),
                           Row = rep(1:nRow,nRange),
                           SystNitro = rep(SystLevel,each=nRow),
                           RandNitro = rep(RandLevel,each=nRow),
                           # Rep = rep(1:4,each = 5*nRow),
                           Rep = rep(1:nRep,each = length(nitro)*nRow),
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
  nRange <- paralist$nRep * length(paralist$Nitro)
  nRow   <- paralist$nRow
  
  eta <- 1
  data_frames_list <- list()
  
  for(i in 1:M){
    
    QuadraDat <- SimuGenDat_Qua(paralist,"Matern",eta)
    
    sp_mat_QuadraDat <- sp::SpatialPointsDataFrame(cbind(QuadraDat$Range,QuadraDat$Row),QuadraDat)
    
    gwr_mat_qua.syst9 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_mat_QuadraDat, bw=9,
                                            kernel = "gaussian")  ## boxcar kernel
    
    QuadraDat$SystGWRB0_bw9 <- gwr_mat_qua.syst9$SDF$Intercept
    QuadraDat$SystGWRB1_bw9 <- gwr_mat_qua.syst9$SDF$SystNitro
    QuadraDat$SystGWRB2_bw9 <- gwr_mat_qua.syst9$SDF$SystNitro2
    
    gwr_mat_qua.rand9 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_mat_QuadraDat, bw = 9, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    QuadraDat$RandGWRB0_bw9 <- gwr_mat_qua.rand9$SDF$Intercept
    QuadraDat$RandGWRB1_bw9 <- gwr_mat_qua.rand9$SDF$RandNitro
    QuadraDat$RandGWRB2_bw9 <- gwr_mat_qua.rand9$SDF$RandNitro2
    

    gwr_mat_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_mat_QuadraDat, bw=5,
                                            kernel = "gaussian")  ## boxcar kernel
    
    QuadraDat$SystGWRB0_bw5 <- gwr_mat_qua.syst5$SDF$Intercept
    QuadraDat$SystGWRB1_bw5 <- gwr_mat_qua.syst5$SDF$SystNitro
    QuadraDat$SystGWRB2_bw5 <- gwr_mat_qua.syst5$SDF$SystNitro2
    
    gwr_mat_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_mat_QuadraDat, bw = 5, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    QuadraDat$RandGWRB0_bw5 <- gwr_mat_qua.rand5$SDF$Intercept
    QuadraDat$RandGWRB1_bw5 <- gwr_mat_qua.rand5$SDF$RandNitro
    QuadraDat$RandGWRB2_bw5 <- gwr_mat_qua.rand5$SDF$RandNitro2
    
    data_frames_list[[i]] <- QuadraDat
    
    cat(i)
  }
  return(data_frames_list)
}

set.seed(649)
st <- Sys.time()
parallel_Mat_coef_1k_results <- parLapply(cl, 1:10, mat_simu_coef_func)
saveRDS(parallel_Mat_coef_1k_results,'parallel_Mat_coef_results_1k_mac.rds')
ed <- Sys.time()
ed - st



############## fit bw= 5 ####################

st<- Sys.time()
for(i in 1:10){
  for(j in 1:100){
    temp <- parallel_Mat_coef_results[[i]][[j]]
    
    sp_mat_QuadraDat <- sp::SpatialPointsDataFrame(cbind(temp$Range,temp$Row),
                                                   temp)
    gwr_mat_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_mat_QuadraDat, bw=5,
                                            kernel = "gaussian")  ## boxcar kernel
    
    temp$SystGWRB0_bw5 <- gwr_mat_qua.syst5$SDF$Intercept
    temp$SystGWRB1_bw5 <- gwr_mat_qua.syst5$SDF$SystNitro
    temp$SystGWRB2_bw5 <- gwr_mat_qua.syst5$SDF$SystNitro2
    
    gwr_mat_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_mat_QuadraDat, bw = 5, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    temp$RandGWRB0_bw5 <- gwr_mat_qua.rand5$SDF$Intercept
    temp$RandGWRB1_bw5 <- gwr_mat_qua.rand5$SDF$RandNitro
    temp$RandGWRB2_bw5 <- gwr_mat_qua.rand5$SDF$RandNitro2
    
    parallel_Mat_coef_results[[i]][[j]] <- temp
  }
  cat(i)
}
ed<- Sys.time()



st<- Sys.time()
for(i in 1:10){
  for(j in 1:100){
    temp <- parallel_AR1_coef_results[[i]][[j]]
    
    sp_ar1_QuadraDat <- sp::SpatialPointsDataFrame(cbind(temp$Range,temp$Row),
                                                   temp)
    gwr_ar1_qua.syst5 <- GWmodel::gwr.basic(SYield ~ SystNitro + SystNitro2,
                                            data = sp_ar1_QuadraDat, bw=5,
                                            kernel = "gaussian")  ## boxcar kernel
    
    temp$SystGWRB0_bw5 <- gwr_ar1_qua.syst5$SDF$Intercept
    temp$SystGWRB1_bw5 <- gwr_ar1_qua.syst5$SDF$SystNitro
    temp$SystGWRB2_bw5 <- gwr_ar1_qua.syst5$SDF$SystNitro2
    
    gwr_ar1_qua.rand5 <- GWmodel::gwr.basic(RYield ~ RandNitro + RandNitro2,
                                            data = sp_ar1_QuadraDat, bw = 5, 
                                            kernel = "gaussian")  ## boxcar kernel
    
    temp$RandGWRB0_bw5 <- gwr_ar1_qua.rand5$SDF$Intercept
    temp$RandGWRB1_bw5 <- gwr_ar1_qua.rand5$SDF$RandNitro
    temp$RandGWRB2_bw5 <- gwr_ar1_qua.rand5$SDF$RandNitro2
    
    parallel_AR1_coef_results[[i]][[j]] <- temp
  }
  cat(i)
}
ed<- Sys.time()



# saveRDS(parallel_Mat_coef_results,'Outcome/ExptNewPara/Simu1K2024/parallel_Mat_coef_results_1k_pc_withbd5.rds')





