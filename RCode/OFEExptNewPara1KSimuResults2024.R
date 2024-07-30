######### load data ############
rm(list = ls())
load("~/Curtin/Research/BayesTechOFE/ExptOut1Ksimu.RData")

### 1 k simulations

### low correlation eta = 1
#### results , low correlation eta = 1 ####
# NS 
MSE_ns_LR5_sq  <- apply(MSE_ns_LR5, c(1,2), function(a) mean(a^2))
MSE_ns_LR9_sq  <- apply(MSE_ns_LR9, c(1,2), function(a) mean(a^2))
MSE_ns_LRbw_sq <- apply(MSE_ns_LRbw, c(1,2), function(a) mean(a^2))

MSE_ns_LS5_sq <- apply(MSE_ns_LS5, c(1,2), function(a) mean(a^2))
MSE_ns_LS9_sq <- apply(MSE_ns_LS9, c(1,2), function(a) mean(a^2))
MSE_ns_LSbw_sq <- apply(MSE_ns_LSbw, c(1,2), function(a) mean(a^2))

# ar1
MSE_ar1_LR5_sq <- apply(MSE_ar1_LR5, c(1,2), function(a) mean(a^2))
MSE_ar1_LR9_sq <- apply(MSE_ar1_LR9, c(1,2), function(a) mean(a^2))
MSE_ar1_LRbw_sq <- apply(MSE_ar1_LRbw, c(1,2), function(a) mean(a^2))

MSE_ar1_LS5_sq <- apply(MSE_ar1_LS5, c(1,2), function(a) mean(a^2))
MSE_ar1_LS9_sq <- apply(MSE_ar1_LS9, c(1,2), function(a) mean(a^2))
MSE_ar1_LSbw_sq <- apply(MSE_ar1_LSbw, c(1,2), function(a) mean(a^2))

# matern
MSE_mat_LR5_sq <- apply(MSE_mat_LR5, c(1,2), function(a) mean(a^2))
MSE_mat_LR9_sq <- apply(MSE_mat_LR9, c(1,2), function(a) mean(a^2))
MSE_mat_LRbw_sq <- apply(MSE_mat_LRbw, c(1,2), function(a) mean(a^2))

MSE_mat_LS5_sq <- apply(MSE_mat_LS5, c(1,2), function(a) mean(a^2))
MSE_mat_LS9_sq <- apply(MSE_mat_LS9, c(1,2), function(a) mean(a^2))
MSE_mat_LSbw_sq <- apply(MSE_mat_LSbw, c(1,2), function(a) mean(a^2))


# M <- 1000

LinBeta0_ns <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_LR5_sq[,1],MSE_ns_LR9_sq[,1],MSE_ns_LRbw_sq[,1],
            MSE_ns_LS5_sq[,1],MSE_ns_LS9_sq[,1],MSE_ns_LSbw_sq[,1])
)

LinBeta1_ns <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_LR5_sq[,2],MSE_ns_LR9_sq[,2],MSE_ns_LRbw_sq[,2],
            MSE_ns_LS5_sq[,2],MSE_ns_LS9_sq[,2],MSE_ns_LSbw_sq[,2])
)


LinBeta0_ar1 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_LR5_sq[,1],MSE_ar1_LR9_sq[,1],MSE_ar1_LRbw_sq[,1],
            MSE_ar1_LS5_sq[,1],MSE_ar1_LS9_sq[,1],MSE_ar1_LSbw_sq[,1])
)

LinBeta1_ar1 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_LR5_sq[,2],MSE_ar1_LR9_sq[,2],MSE_ar1_LRbw_sq[,2],
            MSE_ar1_LS5_sq[,2],MSE_ar1_LS9_sq[,2],MSE_ar1_LSbw_sq[,2])
)


LinBeta0_mat <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_LR5_sq[,1],MSE_mat_LR9_sq[,1],MSE_mat_LRbw_sq[,1],
            MSE_mat_LS5_sq[,1],MSE_mat_LS9_sq[,1],MSE_mat_LSbw_sq[,1])
)


LinBeta1_mat <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_LR5_sq[,2],MSE_mat_LR9_sq[,2],MSE_mat_LRbw_sq[,2],
            MSE_mat_LS5_sq[,2],MSE_mat_LS9_sq[,2],MSE_mat_LSbw_sq[,2])
)


## quadratic 

MSE_ns_QR5_sq <- apply(MSE_ns_QR5, c(1,2), function(a) mean(a^2))
MSE_ns_QR9_sq <- apply(MSE_ns_QR9, c(1,2), function(a) mean(a^2))
MSE_ns_QRbw_sq <- apply(MSE_ns_QRbw, c(1,2), function(a) mean(a^2))


MSE_ns_QS5_sq <- apply(MSE_ns_QS5, c(1,2), function(a) mean(a^2))
MSE_ns_QS9_sq <- apply(MSE_ns_QS9, c(1,2), function(a) mean(a^2))
MSE_ns_QSbw_sq <- apply(MSE_ns_QSbw, c(1,2), function(a) mean(a^2))



MSE_ar1_QR5_sq <- apply(MSE_ar1_QR5, c(1,2), function(a) mean(a^2))
MSE_ar1_QR9_sq <- apply(MSE_ar1_QR9, c(1,2), function(a) mean(a^2))
MSE_ar1_QRbw_sq <- apply(MSE_ar1_QRbw, c(1,2), function(a) mean(a^2))


MSE_ar1_QS5_sq <- apply(MSE_ar1_QS5, c(1,2), function(a) mean(a^2))
MSE_ar1_QS9_sq <- apply(MSE_ar1_QS9, c(1,2), function(a) mean(a^2))
MSE_ar1_QSbw_sq <- apply(MSE_ar1_QSbw, c(1,2), function(a) mean(a^2))


MSE_mat_QR5_sq <- apply(MSE_mat_QR5, c(1,2), function(a) mean(a^2))
MSE_mat_QR9_sq <- apply(MSE_mat_QR9, c(1,2), function(a) mean(a^2))
MSE_mat_QRbw_sq <- apply(MSE_mat_QRbw, c(1,2), function(a) mean(a^2))


MSE_mat_QS5_sq <- apply(MSE_mat_QS5, c(1,2), function(a) mean(a^2))
MSE_mat_QS9_sq <- apply(MSE_mat_QS9, c(1,2), function(a) mean(a^2))
MSE_mat_QSbw_sq <- apply(MSE_mat_QSbw, c(1,2), function(a) mean(a^2))



QuaBeta0_ns <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_QR5_sq[,1],MSE_ns_QR9_sq[,1],MSE_ns_QRbw_sq[,1],
            MSE_ns_QS5_sq[,1],MSE_ns_QS9_sq[,1],MSE_ns_QSbw_sq[,1])
)


QuaBeta1_ns <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_QR5_sq[,2],MSE_ns_QR9_sq[,2],MSE_ns_QRbw_sq[,2],
            MSE_ns_QS5_sq[,2],MSE_ns_QS9_sq[,2],MSE_ns_QSbw_sq[,2])
)

QuaBeta2_ns <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_QR5_sq[,3],MSE_ns_QR9_sq[,3],MSE_ns_QRbw_sq[,3],
            MSE_ns_QS5_sq[,3],MSE_ns_QS9_sq[,3],MSE_ns_QSbw_sq[,3])
)


QuaBeta0_ar1 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_QR5_sq[,1],MSE_ar1_QR9_sq[,1],MSE_ar1_QRbw_sq[,1],
            MSE_ar1_QS5_sq[,1],MSE_ar1_QS9_sq[,1],MSE_ar1_QSbw_sq[,1])
)

QuaBeta1_ar1 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_QR5_sq[,2],MSE_ar1_QR9_sq[,2],MSE_ar1_QRbw_sq[,2],
            MSE_ar1_QS5_sq[,2],MSE_ar1_QS9_sq[,2],MSE_ar1_QSbw_sq[,2])
)

QuaBeta2_ar1 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_QR5_sq[,3],MSE_ar1_QR9_sq[,3],MSE_ar1_QRbw_sq[,3],
            MSE_ar1_QS5_sq[,3],MSE_ar1_QS9_sq[,3],MSE_ar1_QSbw_sq[,3])
)


QuaBeta0_mat <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_QR5_sq[,1],MSE_mat_QR9_sq[,1],MSE_mat_QRbw_sq[,1],
            MSE_mat_QS5_sq[,1],MSE_mat_QS9_sq[,1],MSE_mat_QSbw_sq[,1])
)

QuaBeta1_mat <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_QR5_sq[,2],MSE_mat_QR9_sq[,2],MSE_mat_QRbw_sq[,2],
            MSE_mat_QS5_sq[,2],MSE_mat_QS9_sq[,2],MSE_mat_QSbw_sq[,2])
)

QuaBeta2_mat <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_QR5_sq[,3],MSE_mat_QR9_sq[,3],MSE_mat_QRbw_sq[,3],
            MSE_mat_QS5_sq[,3],MSE_mat_QS9_sq[,3],MSE_mat_QSbw_sq[,3])
)


###### high correlation eta = 0.1 ######### 
#### results ####
# NS 
MSE_ns_LR5_sq_eta01 <- apply(MSE_ns_LR5_eta01, c(1,2), function(a) mean(a^2))
MSE_ns_LR9_sq_eta01 <- apply(MSE_ns_LR9_eta01, c(1,2), function(a) mean(a^2))
MSE_ns_LRbw_sq_eta01 <- apply(MSE_ns_LRbw_eta01, c(1,2), function(a) mean(a^2))

MSE_ns_LS5_sq_eta01 <- apply(MSE_ns_LS5_eta01, c(1,2), function(a) mean(a^2))
MSE_ns_LS9_sq_eta01 <- apply(MSE_ns_LS9_eta01, c(1,2), function(a) mean(a^2))
MSE_ns_LSbw_sq_eta01 <- apply(MSE_ns_LSbw_eta01, c(1,2), function(a) mean(a^2))

# ar1
MSE_ar1_LR5_sq_eta01 <- apply(MSE_ar1_LR5_eta01, c(1,2), function(a) mean(a^2))
MSE_ar1_LR9_sq_eta01 <- apply(MSE_ar1_LR9_eta01, c(1,2), function(a) mean(a^2))
MSE_ar1_LRbw_sq_eta01 <- apply(MSE_ar1_LRbw_eta01, c(1,2), function(a) mean(a^2))

MSE_ar1_LS5_sq_eta01 <- apply(MSE_ar1_LS5_eta01, c(1,2), function(a) mean(a^2))
MSE_ar1_LS9_sq_eta01 <- apply(MSE_ar1_LS9_eta01, c(1,2), function(a) mean(a^2))
MSE_ar1_LSbw_sq_eta01 <- apply(MSE_ar1_LSbw_eta01, c(1,2), function(a) mean(a^2))

# matern
MSE_mat_LR5_sq_eta01 <- apply(MSE_mat_LR5_eta01, c(1,2), function(a) mean(a^2))
MSE_mat_LR9_sq_eta01 <- apply(MSE_mat_LR9_eta01, c(1,2), function(a) mean(a^2))
MSE_mat_LRbw_sq_eta01 <- apply(MSE_mat_LRbw_eta01, c(1,2), function(a) mean(a^2))

MSE_mat_LS5_sq_eta01 <- apply(MSE_mat_LS5_eta01, c(1,2), function(a) mean(a^2))
MSE_mat_LS9_sq_eta01 <- apply(MSE_mat_LS9_eta01, c(1,2), function(a) mean(a^2))
MSE_mat_LSbw_sq_eta01 <- apply(MSE_mat_LSbw_eta01, c(1,2), function(a) mean(a^2))


LinBeta0_ns_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_LR5_sq_eta01[,1],MSE_ns_LR9_sq_eta01[,1],MSE_ns_LRbw_sq_eta01[,1],
            MSE_ns_LS5_sq_eta01[,1],MSE_ns_LS9_sq_eta01[,1],MSE_ns_LSbw_sq_eta01[,1])
)

LinBeta1_ns_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_LR5_sq_eta01[,2],MSE_ns_LR9_sq_eta01[,2],MSE_ns_LRbw_sq_eta01[,2],
            MSE_ns_LS5_sq_eta01[,2],MSE_ns_LS9_sq_eta01[,2],MSE_ns_LSbw_sq_eta01[,2])
)

LinBeta0_ar1_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_LR5_sq_eta01[,1],MSE_ar1_LR9_sq_eta01[,1],MSE_ar1_LRbw_sq_eta01[,1],
            MSE_ar1_LS5_sq_eta01[,1],MSE_ar1_LS9_sq_eta01[,1],MSE_ar1_LSbw_sq_eta01[,1])
)

LinBeta1_ar1_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_LR5_sq_eta01[,2],MSE_ar1_LR9_sq_eta01[,2],MSE_ar1_LRbw_sq_eta01[,2],
            MSE_ar1_LS5_sq_eta01[,2],MSE_ar1_LS9_sq_eta01[,2],MSE_ar1_LSbw_sq_eta01[,2])
)

LinBeta0_mat_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_LR5_sq_eta01[,1],MSE_mat_LR9_sq_eta01[,1],MSE_mat_LRbw_sq_eta01[,1],
            MSE_mat_LS5_sq_eta01[,1],MSE_mat_LS9_sq_eta01[,1],MSE_mat_LSbw_sq_eta01[,1])
)


LinBeta1_mat_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_LR5_sq_eta01[,2],MSE_mat_LR9_sq_eta01[,2],MSE_mat_LRbw_sq_eta01[,2],
            MSE_mat_LS5_sq_eta01[,2],MSE_mat_LS9_sq_eta01[,2],MSE_mat_LSbw_sq_eta01[,2])
)


## quadratic 
MSE_ns_QR5_sq_eta01 <- apply(MSE_ns_QR5_eta01, c(1,2), function(a) mean(a^2))
MSE_ns_QR9_sq_eta01 <- apply(MSE_ns_QR9_eta01, c(1,2), function(a) mean(a^2))
MSE_ns_QRbw_sq_eta01 <- apply(MSE_ns_QRbw_eta01, c(1,2), function(a) mean(a^2))


MSE_ns_QS5_sq_eta01 <- apply(MSE_ns_QS5_eta01, c(1,2), function(a) mean(a^2))
MSE_ns_QS9_sq_eta01 <- apply(MSE_ns_QS9_eta01, c(1,2), function(a) mean(a^2))
MSE_ns_QSbw_sq_eta01 <- apply(MSE_ns_QSbw_eta01, c(1,2), function(a) mean(a^2))


MSE_ar1_QR5_sq_eta01 <- apply(MSE_ar1_QR5_eta01, c(1,2), function(a) mean(a^2))
MSE_ar1_QR9_sq_eta01 <- apply(MSE_ar1_QR9_eta01, c(1,2), function(a) mean(a^2))
MSE_ar1_QRbw_sq_eta01 <- apply(MSE_ar1_QRbw_eta01, c(1,2), function(a) mean(a^2))


MSE_ar1_QS5_sq_eta01 <- apply(MSE_ar1_QS5_eta01, c(1,2), function(a) mean(a^2))
MSE_ar1_QS9_sq_eta01 <- apply(MSE_ar1_QS9_eta01, c(1,2), function(a) mean(a^2))
MSE_ar1_QSbw_sq_eta01 <- apply(MSE_ar1_QSbw_eta01, c(1,2), function(a) mean(a^2))


MSE_mat_QR5_sq_eta01 <- apply(MSE_mat_QR5_eta01, c(1,2), function(a) mean(a^2))
MSE_mat_QR9_sq_eta01 <- apply(MSE_mat_QR9_eta01, c(1,2), function(a) mean(a^2))
MSE_mat_QRbw_sq_eta01 <- apply(MSE_mat_QRbw_eta01, c(1,2), function(a) mean(a^2))

MSE_mat_QS5_sq_eta01 <- apply(MSE_mat_QS5_eta01, c(1,2), function(a) mean(a^2))
MSE_mat_QS9_sq_eta01 <- apply(MSE_mat_QS9_eta01, c(1,2), function(a) mean(a^2))
MSE_mat_QSbw_sq_eta01 <- apply(MSE_mat_QSbw_eta01, c(1,2), function(a) mean(a^2))



QuaBeta0_ns_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_QR5_sq_eta01[,1],MSE_ns_QR9_sq_eta01[,1],MSE_ns_QRbw_sq_eta01[,1],
            MSE_ns_QS5_sq_eta01[,1],MSE_ns_QS9_sq_eta01[,1],MSE_ns_QSbw_sq_eta01[,1])
)

QuaBeta1_ns_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_QR5_sq_eta01[,2],MSE_ns_QR9_sq_eta01[,2],MSE_ns_QRbw_sq_eta01[,2],
            MSE_ns_QS5_sq_eta01[,2],MSE_ns_QS9_sq_eta01[,2],MSE_ns_QSbw_sq_eta01[,2])
)


QuaBeta2_ns_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ns_QR5_sq_eta01[,3],MSE_ns_QR9_sq_eta01[,3],MSE_ns_QRbw_sq_eta01[,3],
            MSE_ns_QS5_sq_eta01[,3],MSE_ns_QS9_sq_eta01[,3],MSE_ns_QSbw_sq_eta01[,3])
)


QuaBeta0_ar1_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_QR5_sq_eta01[,1],MSE_ar1_QR9_sq_eta01[,1],MSE_ar1_QRbw_sq_eta01[,1],
            MSE_ar1_QS5_sq_eta01[,1],MSE_ar1_QS9_sq_eta01[,1],MSE_ar1_QSbw_sq_eta01[,1])
)

QuaBeta1_ar1_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_QR5_sq_eta01[,2],MSE_ar1_QR9_sq_eta01[,2],MSE_ar1_QRbw_sq_eta01[,2],
            MSE_ar1_QS5_sq_eta01[,2],MSE_ar1_QS9_sq_eta01[,2],MSE_ar1_QSbw_sq_eta01[,2])
)


QuaBeta2_ar1_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_ar1_QR5_sq_eta01[,3],MSE_ar1_QR9_sq_eta01[,3],MSE_ar1_QRbw_sq_eta01[,3],
            MSE_ar1_QS5_sq_eta01[,3],MSE_ar1_QS9_sq_eta01[,3],MSE_ar1_QSbw_sq_eta01[,3])
)


QuaBeta0_mat_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_QR5_sq_eta01[,1],MSE_mat_QR9_sq_eta01[,1],MSE_mat_QRbw_sq_eta01[,1],
            MSE_mat_QS5_sq_eta01[,1],MSE_mat_QS9_sq_eta01[,1],MSE_mat_QSbw_sq_eta01[,1])
)


QuaBeta1_mat_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_QR5_sq_eta01[,2],MSE_mat_QR9_sq_eta01[,2],MSE_mat_QRbw_sq_eta01[,2],
            MSE_mat_QS5_sq_eta01[,2],MSE_mat_QS9_sq_eta01[,2],MSE_mat_QSbw_sq_eta01[,2])
)


QuaBeta2_mat_eta01 <- data.frame(
  BandWidth = c(rep("BW5",M),rep("BW9",M),rep("BWAICc",M),
           rep("BW5",M),rep("BW9",M),rep("BWAICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  Value = c(MSE_mat_QR5_sq_eta01[,3],MSE_mat_QR9_sq_eta01[,3],MSE_mat_QRbw_sq_eta01[,3],
            MSE_mat_QS5_sq_eta01[,3],MSE_mat_QS9_sq_eta01[,3],MSE_mat_QSbw_sq_eta01[,3])
)



########## further figures #########

LinBeta0_ns$Cov <- "NS"
LinBeta0_ar1$Cov <- "AR1"
LinBeta0_mat$Cov <- "Matern"

LinBeta0_ns$Beta <- 0
LinBeta0_ar1$Beta <- 0
LinBeta0_mat$Beta <- 0

LinBeta1_ns$Cov <- "NS"
LinBeta1_ar1$Cov <- "AR1"
LinBeta1_mat$Cov <- "Matern"

LinBeta1_ns$Beta <- 1
LinBeta1_ar1$Beta <- 1
LinBeta1_mat$Beta <- 1


LinComb <- rbind(LinBeta0_ns,LinBeta0_ar1,LinBeta0_mat,
                 LinBeta1_ns,LinBeta1_ar1,LinBeta1_mat)

LinComb$Cov <- factor(LinComb$Cov,levels = c("NS","AR1","Matern"))
LinComb$FacBeta <- factor(LinComb$Beta)


QuaBeta0_ns$Cov <- "NS"
QuaBeta0_ar1$Cov <- "AR1"
QuaBeta0_mat$Cov <- "Matern"

QuaBeta0_ns$Beta <- 0
QuaBeta0_ar1$Beta <- 0
QuaBeta0_mat$Beta <- 0

QuaBeta1_ns$Cov <- "NS"
QuaBeta1_ar1$Cov <- "AR1"
QuaBeta1_mat$Cov <- "Matern"

QuaBeta1_ns$Beta <- 1
QuaBeta1_ar1$Beta <- 1
QuaBeta1_mat$Beta <- 1

QuaBeta2_ns$Cov <- "NS"
QuaBeta2_ar1$Cov <- "AR1"
QuaBeta2_mat$Cov <- "Matern"

QuaBeta2_ns$Beta <- 2
QuaBeta2_ar1$Beta <- 2
QuaBeta2_mat$Beta <- 2



QuaComb <- rbind(QuaBeta0_ns,QuaBeta0_ar1,QuaBeta0_mat,
                 QuaBeta1_ns,QuaBeta1_ar1,QuaBeta1_mat,
                 QuaBeta2_ns,QuaBeta2_ar1,QuaBeta2_mat)

QuaComb$Cov <- factor(QuaComb$Cov,levels = c("NS","AR1","Matern"))
QuaComb$FacBeta <- factor(QuaComb$Beta)



LinBeta0_ns_eta01$Cov <- "NS"
LinBeta0_ar1_eta01$Cov <- "AR1"
LinBeta0_mat_eta01$Cov <- "Matern"

LinBeta0_ns_eta01$Beta <- 0
LinBeta0_ar1_eta01$Beta <- 0
LinBeta0_mat_eta01$Beta <- 0

LinBeta1_ns_eta01$Cov <- "NS"
LinBeta1_ar1_eta01$Cov <- "AR1"
LinBeta1_mat_eta01$Cov <- "Matern"

LinBeta1_ns_eta01$Beta <- 1
LinBeta1_ar1_eta01$Beta <- 1
LinBeta1_mat_eta01$Beta <- 1

LinComb_eta01 <- rbind(LinBeta0_ns_eta01,LinBeta0_ar1_eta01,LinBeta0_mat_eta01,
                       LinBeta1_ns_eta01,LinBeta1_ar1_eta01,LinBeta1_mat_eta01)

LinComb_eta01$Cov <- factor(LinComb_eta01$Cov,levels = c("NS","AR1","Matern"))
LinComb_eta01$FacBeta <- factor(LinComb_eta01$Beta)


QuaBeta0_ns_eta01$Cov <- "NS"
QuaBeta0_ar1_eta01$Cov <- "AR1"
QuaBeta0_mat_eta01$Cov <- "Matern"

QuaBeta0_ns_eta01$Beta <- 0
QuaBeta0_ar1_eta01$Beta <- 0
QuaBeta0_mat_eta01$Beta <- 0

QuaBeta1_ns_eta01$Cov <- "NS"
QuaBeta1_ar1_eta01$Cov <- "AR1"
QuaBeta1_mat_eta01$Cov <- "Matern"

QuaBeta1_ns_eta01$Beta <- 1
QuaBeta1_ar1_eta01$Beta <- 1
QuaBeta1_mat_eta01$Beta <- 1

QuaBeta2_ns_eta01$Cov <- "NS"
QuaBeta2_ar1_eta01$Cov <- "AR1"
QuaBeta2_mat_eta01$Cov <- "Matern"

QuaBeta2_ns_eta01$Beta <- 2
QuaBeta2_ar1_eta01$Beta <- 2
QuaBeta2_mat_eta01$Beta <- 2


QuaComb_eta01 <- rbind(QuaBeta0_ns_eta01,QuaBeta0_ar1_eta01,QuaBeta0_mat_eta01,
                       QuaBeta1_ns_eta01,QuaBeta1_ar1_eta01,QuaBeta1_mat_eta01,
                       QuaBeta2_ns_eta01,QuaBeta2_ar1_eta01,QuaBeta2_mat_eta01)

QuaComb_eta01$Cov <- factor(QuaComb_eta01$Cov,levels = c("NS","AR1","Matern"))
QuaComb_eta01$FacBeta <- factor(QuaComb_eta01$Beta)



########### anova test #############

library(asreml)

LinComb$Eta <- 1
LinComb_eta01$Eta <- 0.1
LinAll <- rbind(LinComb,LinComb_eta01) 

QuaComb$Eta <- 1
QuaComb_eta01$Eta <- 0.1
QuaAll <- rbind(QuaComb,QuaComb_eta01) 


nms <- c("BandWidth","Design","Eta")

LinAll[nms] <- lapply(LinAll[nms], as.factor)
QuaAll[nms] <- lapply(QuaAll[nms], as.factor)


# linall.as <- asreml(Value ~ Design + BandWidth + Cov +  Eta, ## had + FacBeta
#                   data = LinAll)
# wald(linall.as,denDF = "default",ssType = "conditional")
# 
# 
# quaall.as <- asreml(Value ~ Design + BandWidth + Cov + Eta, ## had + FacBeta
#                     data = QuaAll)
# wald(quaall.as,denDF = "default",ssType = "conditional")
# 


linall.aov <- aov(Value ~ Design + BandWidth + Cov  + Eta +  ## had + FacBeta
                    Design:BandWidth + Design:Cov + Design:Eta +
                    BandWidth:Cov + BandWidth:Eta + Cov:Eta,
                  data = LinAll)
summary(linall.aov)


quaall.aov <- aov(Value ~ Design + BandWidth + Cov  + Eta +  ## had + FacBeta
                    Design:BandWidth + Design:Cov + Design:Eta +
                    BandWidth:Cov + BandWidth:Eta + Cov:Eta,
                  data = QuaAll)
summary(quaall.aov)


combined_anova_table <- cbind(summary(linall.aov)[[1]], 
                              summary(quaall.aov)[[1]])

latex_anova_table <- xtable(combined_anova_table,digits = 3)
# print(latex_anova_table, type = "latex")



library(tidyverse)
### Median tables 
# linear eta = 1
LinTable <- c(apply(MSE_ns_LR5_sq,2, median), ## random
              apply(MSE_ar1_LR5_sq,2, median),
              apply(MSE_mat_LR5_sq,2, median),
              apply(MSE_ns_LR9_sq,2, median),
              apply(MSE_ar1_LR9_sq,2, median),
              apply(MSE_mat_LR9_sq,2, median),
              apply(MSE_ns_LRbw_sq,2, median),
              apply(MSE_ar1_LRbw_sq,2, median),
              apply(MSE_mat_LRbw_sq,2, median),
              ## system
              apply(MSE_ns_LS5_sq,2, median),
              apply(MSE_ar1_LS5_sq,2, median),
              apply(MSE_mat_LS5_sq,2, median),
              apply(MSE_ns_LS9_sq,2, median),
              apply(MSE_ar1_LS9_sq,2, median),
              apply(MSE_mat_LS9_sq,2, median),
              apply(MSE_ns_LSbw_sq,2, median),
              apply(MSE_ar1_LSbw_sq,2, median),
              apply(MSE_mat_LSbw_sq,2, median)
)

LinTable_df <- data.frame(matrix(LinTable, ncol = 6))
LinTable_df <- LinTable_df %>% 
  mutate_all(~ifelse(row_number() %% 2 == 0, . * 1e4, .))

xtable(LinTable_df, digits = 3)



# quadratic eta = 1
QuaTable <- c(apply(MSE_ns_QR5_sq,2, median), ## random
              apply(MSE_ar1_QR5_sq,2, median),
              apply(MSE_mat_QR5_sq,2, median),
              apply(MSE_ns_QR9_sq,2, median),
              apply(MSE_ar1_QR9_sq,2, median),
              apply(MSE_mat_QR9_sq,2, median),
              apply(MSE_ns_QRbw_sq,2, median),
              apply(MSE_ar1_QRbw_sq,2, median),
              apply(MSE_mat_QRbw_sq,2, median),
              ## system
              apply(MSE_ns_QS5_sq,2, median),
              apply(MSE_ar1_QS5_sq,2, median),
              apply(MSE_mat_QS5_sq,2, median),
              apply(MSE_ns_QS9_sq,2, median),
              apply(MSE_ar1_QS9_sq,2, median),
              apply(MSE_mat_QS9_sq,2, median),
              apply(MSE_ns_QSbw_sq,2, median),
              apply(MSE_ar1_QSbw_sq,2, median),
              apply(MSE_mat_QSbw_sq,2, median)
)

QuaTable_df <- data.frame(matrix(QuaTable, ncol = 6))

QuaTable_df <- QuaTable_df %>% 
  mutate_all(~ifelse(row_number() %% 3 == 2, . * 1e4, 
                     ifelse(row_number() %% 3 == 0, . * 1e8, .)))

latex_table <- xtable(QuaTable_df, digits = 3)

xtable(QuaTable_df, digits = 3)



### Median tables 
# linear eta = 0.1, high correlation
LinTable_eta01 <- c(apply(MSE_ns_LR5_sq_eta01,2, median), ## random
              apply(MSE_ar1_LR5_sq_eta01,2, median),
              apply(MSE_mat_LR5_sq_eta01,2, median),
              apply(MSE_ns_LR9_sq_eta01,2, median),
              apply(MSE_ar1_LR9_sq_eta01,2, median),
              apply(MSE_mat_LR9_sq_eta01,2, median),
              apply(MSE_ns_LRbw_sq_eta01,2, median),
              apply(MSE_ar1_LRbw_sq_eta01,2, median),
              apply(MSE_mat_LRbw_sq_eta01,2, median),
              ## system
              apply(MSE_ns_LS5_sq_eta01,2, median),
              apply(MSE_ar1_LS5_sq_eta01,2, median),
              apply(MSE_mat_LS5_sq_eta01,2, median),
              apply(MSE_ns_LS9_sq_eta01,2, median),
              apply(MSE_ar1_LS9_sq_eta01,2, median),
              apply(MSE_mat_LS9_sq_eta01,2, median),
              apply(MSE_ns_LSbw_sq_eta01,2, median),
              apply(MSE_ar1_LSbw_sq_eta01,2, median),
              apply(MSE_mat_LSbw_sq_eta01,2, median)
)

LinTable_df_eta01 <- data.frame(matrix(LinTable_eta01, ncol = 6))
LinTable_df_eta01 <- LinTable_df_eta01 %>% 
  mutate_all(~ifelse(row_number() %% 2 == 0, . * 1e4, .))

xtable(LinTable_df_eta01, digits = 3)


# quadratic eta = 0.1, high correlation
QuaTable_eta01 <- c(apply(MSE_ns_QR5_sq_eta01,2, median), ## random
              apply(MSE_ar1_QR5_sq_eta01,2, median),
              apply(MSE_mat_QR5_sq_eta01,2, median),
              apply(MSE_ns_QR9_sq_eta01,2, median),
              apply(MSE_ar1_QR9_sq_eta01,2, median),
              apply(MSE_mat_QR9_sq_eta01,2, median),
              apply(MSE_ns_QRbw_sq_eta01,2, median),
              apply(MSE_ar1_QRbw_sq_eta01,2, median),
              apply(MSE_mat_QRbw_sq_eta01,2, median),
              ## system
              apply(MSE_ns_QS5_sq_eta01,2, median),
              apply(MSE_ar1_QS5_sq_eta01,2, median),
              apply(MSE_mat_QS5_sq_eta01,2, median),
              apply(MSE_ns_QS9_sq_eta01,2, median),
              apply(MSE_ar1_QS9_sq_eta01,2, median),
              apply(MSE_mat_QS9_sq_eta01,2, median),
              apply(MSE_ns_QSbw_sq_eta01,2, median),
              apply(MSE_ar1_QSbw_sq_eta01,2, median),
              apply(MSE_mat_QSbw_sq_eta01,2, median)
)

QuaTable_df_eta01 <- data.frame(matrix(QuaTable_eta01, ncol = 6))

QuaTable_df_eta01 <- QuaTable_df_eta01 %>% 
  mutate_all(~ifelse(row_number() %% 3 == 2, . * 1e4, 
                     ifelse(row_number() %% 3 == 0, . * 1e8, .)))

latex_table_eta01 <- xtable(QuaTable_df_eta01, digits = 3)

xtable(QuaTable_df_eta01, digits = 3)




##### figures ##### 

library(ggplot2)
library(ggpubr)
library(viridis)

## "NS" = "#FFB3BA", 
# custom_colors <- c("Randomised" = "#BAE1FF", "Systematic" = "#BFFCC6")

custom_colors <- c("Randomised" = "gray70", "Systematic" = "white")

generate_plots <- function(data, custom_colors) {
  data$Cov <- factor(data$Cov, levels = c("NS","AR1","Matern" ))
  
  plot_list <- lapply(unique(data$Beta), function(beta, Eta = unique(data$Eta)) {
    lapply(sort(unique(data$Cov)), function(cov) {
      p <- ggplot(subset(data, Beta == beta & Cov == as.character(cov)), 
                  aes(BandWidth, Value)) +
        geom_boxplot(aes(fill = Design)) +
        ylab("") + xlab("") +
        thm1 +
        scale_x_discrete(labels = c("5", "9", "AICc")) +
        scale_fill_manual(values = custom_colors) +
        theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
      
      if (beta == 1) {
        p <- p +       
          labs(title = bquote(beta[.(beta)] ~ .("x") ~10^4 ~ .(", ") ~ 
                                .(as.character(cov)) ~ .(", ") ~ epsilon ~"=" ~.(Eta)))+
          scale_y_log10(labels = function(x) x * 1e4)
      } else if (beta == 2) {
        p <- p +       
          labs(title = bquote(beta[.(beta)] ~ .("x") ~10^8 ~ .(", ") ~ 
                                .(as.character(cov)) ~ .(", ") ~ epsilon ~ "=" ~.(Eta))) +
          scale_y_log10(labels = function(x) x * 1e8)
      } else {
        p <- p +
          labs(title = bquote(beta[.(beta)] ~ ", " ~.(as.character(cov)) ~ .(", ") ~ epsilon ~ "=" ~.(Eta))) +
          scale_y_log10()
      }
      
      return(p)
    })
  })
  
  # Flatten the list of lists
  plot_list <- do.call(c, plot_list)
  
  return(plot_list)
}


plot_list_LinComb <- generate_plots(LinComb, custom_colors)
combined_plot_LinComb <- ggarrange(plotlist = plot_list_LinComb, 
                                   ncol = 3, nrow = 2,
                                   align = "hv", 
                                   # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                   common.legend = TRUE, legend = "top")
combined_plot_LinComb <- annotate_figure(combined_plot_LinComb,
                                         bottom = text_grob("Bandwidth", size = 20),
                                         left = text_grob("MSE", size = 20, rot = 90))
print(combined_plot_LinComb)
## Col_LinCombMSE_newpar_V3


plot_list_QuaComb <- generate_plots(QuaComb, custom_colors)
combined_plot_QuaComb <- ggarrange(plotlist = plot_list_QuaComb, 
                                   ncol = 3, nrow = 3,
                                   align = "hv", 
                                   # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                   common.legend = TRUE, legend = "top")
combined_plot_QuaComb <- annotate_figure(combined_plot_QuaComb,
                                         bottom = text_grob("Bandwidth", size = 20),
                                         left = text_grob("MSE", size = 20, rot = 90))
print(combined_plot_QuaComb)
## Col_QuaCombMSE_newpar_V3


plot_list_LinComb_eta01 <- generate_plots(LinComb_eta01, custom_colors)
combined_plot_LinComb_eta01 <- ggarrange(plotlist = plot_list_LinComb_eta01, 
                                         ncol = 3, nrow = 2,
                                         align = "hv", 
                                         # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                         common.legend = TRUE, legend = "top")
combined_plot_LinComb_eta01 <- annotate_figure(combined_plot_LinComb_eta01,
                                               bottom = text_grob("Bandwidth", size = 20),
                                               left = text_grob("MSE", size = 20, rot = 90))
print(combined_plot_LinComb_eta01)
## Col_LinCombMSE_newpar_eta01_V3



plot_list_QuaComb_eta01 <- generate_plots(QuaComb_eta01, custom_colors)
combined_plot_QuaComb_eta01 <- ggarrange(plotlist = plot_list_QuaComb_eta01, 
                                         ncol = 3, nrow = 3,
                                         align = "hv", 
                                         # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                         common.legend = TRUE, legend = "top")
combined_plot_QuaComb_eta01 <- annotate_figure(combined_plot_QuaComb_eta01,
                                               bottom = text_grob("Bandwidth", size = 20),
                                               left = text_grob("MSE", size = 20, rot = 90))
print(combined_plot_QuaComb_eta01)
## Col_QuaCombMSE_newpar_eta01_V3



















######### archive ############

levels(LinComb$Cov)
LinComb$Cov <- factor(LinComb$Cov, levels = c("NS","AR1","Matern" ))

plot_list_LinComb <- lapply(unique(LinComb$Beta), function(beta, Eta = unique(LinComb$Eta)) {
  lapply(sort(unique(LinComb$Cov)), function(cov) {
    p <- ggplot(subset(LinComb, Beta == beta & Cov == as.character(cov)), 
                aes(BandWidth, Value)) +
      geom_boxplot(aes(fill = Design)) +
      ylab("") + xlab("") +
      thm1 +
      scale_x_discrete(labels = c("5", "9", "AICc")) +
      scale_fill_manual(values = custom_colors) +
      theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
    
    if (beta == 1) {
      p <- p +       
        labs(title = bquote(beta[.(beta)] ~ .("x") ~10^4 ~ .(", ") ~ 
                              .(as.character(cov)) ~ .(", ") ~ epsilon ~"=" ~.(Eta)))+
      scale_y_log10(labels = function(x) x * 1e4)
    } else {
      p <- p +
        labs(title = bquote(beta[.(beta)] ~ ", " ~.(as.character(cov)) ~ 
                              .(", ") ~ epsilon ~ "=" ~.(Eta))) +
        scale_y_log10()
    }
    
    return(p)
  })
})

# Flatten the list of lists
plot_list_LinComb <- do.call(c, plot_list_LinComb)

combined_plot_LinComb <- ggarrange(plotlist = plot_list_LinComb, 
                           ncol = 3, nrow = 2,
                           align = "hv", 
                           # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                           common.legend = TRUE, legend = "top")

# Add common x and y labels
combined_plot_LinComb <- annotate_figure(combined_plot_LinComb,
                                 bottom = text_grob("Bandwidth", size = 20),
                                 left = text_grob("MSE", size = 20, rot = 90))

# Print the combined plot
print(combined_plot_LinComb)
## Col_LinCombMSE_newpar_V3


QuaComb$Cov <- factor(QuaComb$Cov, levels = c("NS","AR1","Matern" ))
# Create individual plots without x and y labels
plot_list_QuaComb <- lapply(unique(QuaComb$Beta), function(beta,Eta = unique(QuaComb$Eta)) {
  lapply(sort(unique(QuaComb$Cov)), function(cov) {
    p <- ggplot(subset(QuaComb, Beta == beta & Cov == as.character(cov)), aes(BandWidth, Value)) +
      geom_boxplot(aes(fill = Design)) +
      ylab("") + xlab("") +
      thm1 +
      scale_x_discrete(labels = c("5", "9", "AICc")) +
      # scale_fill_viridis_d() + # Apply viridis color palette
      scale_fill_manual(values = custom_colors) + # Apply custom pastel colors
      theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
    
    # If the plot is the 4th to 6th in the list, modify the y-axis labels
    if (beta == 2) {
      p <- p +       
        labs(title = bquote(beta[.(beta)] ~ .("x") ~10^8 ~ .(", ") ~ 
                              .(as.character(cov)) ~ .(", ") ~ epsilon ~ "=" ~.(Eta))) +
        scale_y_log10(labels = function(x) x * 1e8)
    }else if (beta == 1) {
      p <- p +       
        labs(title = bquote(beta[.(beta)] ~ .("x") ~10^4 ~ .(", ") ~ 
                              .(as.character(cov)) ~ .(", ") ~ epsilon ~ "=" ~.(Eta))) +
        scale_y_log10(labels = function(x) x * 1e4)
    } else {
      p <- p +
        labs(title = bquote(beta[.(beta)] ~ ", " ~.(as.character(cov)) ~ .(", ") ~ epsilon ~ "=" ~.(Eta))) +
        scale_y_log10()
    }
  })
})

# Flatten the list of lists
plot_list_QuaComb <- do.call(c, plot_list_QuaComb)

# Arrange plots using ggarrange
combined_plot_QuaComb <- ggarrange(plotlist = plot_list_QuaComb, 
                           ncol = 3, nrow = 3,
                           align = "hv", 
                           # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                           common.legend = TRUE, legend = "top")

# Add common x and y labels
combined_plot_QuaComb <- annotate_figure(combined_plot_QuaComb,
                                 bottom = text_grob("Bandwidth", size = 20),
                                 left = text_grob("MSE", size = 20, rot = 90))
# Print the combined plot
print(combined_plot_QuaComb)
## Col_QuaCombMSE_newpar




LinComb_eta01$Cov <- factor(LinComb_eta01$Cov, levels = c("NS","AR1","Matern" ))

plot_list_LinComb_eta01 <- lapply(unique(LinComb_eta01$Beta), 
                                  function(beta, Eta = unique(LinComb_eta01$Eta)) {
  lapply(sort(unique(LinComb_eta01$Cov)), function(cov) {
    p <- ggplot(subset(LinComb_eta01, Beta == beta & Cov == as.character(cov)), 
                aes(BandWidth, Value)) +
      geom_boxplot(aes(fill = Design)) +
      ylab("") + xlab("") +
      thm1 +
      scale_x_discrete(labels = c("5", "9", "AICc")) +
      scale_fill_manual(values = custom_colors) +
      theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
    
    if (beta == 1) {
      p <- p +       
        labs(title = bquote(beta[.(beta)] ~ .("x") ~10^4 ~ .(", ") ~ 
                              .(as.character(cov)) ~ .(", ") ~ epsilon ~"=" ~.(Eta)))+
        scale_y_log10(labels = function(x) x * 1e4)
    } else {
      p <- p +
        labs(title = bquote(beta[.(beta)] ~ ", " ~.(as.character(cov)) ~ 
                              .(", ") ~ epsilon ~ "=" ~.(Eta))) +
        scale_y_log10()
    }
    
    return(p)
  })
})

# Flatten the list of lists
plot_list_LinComb_eta01 <- do.call(c, plot_list_LinComb_eta01)

combined_plot_LinComb_eta01 <- ggarrange(plotlist = plot_list_LinComb_eta01, 
                                   ncol = 3, nrow = 2,
                                   align = "hv", 
                                   # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                   common.legend = TRUE, legend = "top")

# Add common x and y labels
combined_plot_LinComb_eta01 <- annotate_figure(combined_plot_LinComb_eta01,
                                         bottom = text_grob("Bandwidth", size = 20),
                                         left = text_grob("MSE", size = 20, rot = 90))

# Print the combined plot
print(combined_plot_LinComb_eta01)
## Col_LinCombMSE_newpar_eta01


QuaComb_eta01$Cov <- factor(QuaComb_eta01$Cov, levels = c("NS","AR1","Matern" ))
# Create individual plots without x and y labels
plot_list_QuaComb_eta01 <- lapply(unique(QuaComb_eta01$Beta), 
                                  function(beta,Eta = unique(QuaComb_eta01$Eta)) {
  lapply(sort(unique(QuaComb_eta01$Cov)), function(cov) {
    p <- ggplot(subset(QuaComb_eta01, Beta == beta & Cov == as.character(cov)), aes(BandWidth, Value)) +
      geom_boxplot(aes(fill = Design)) +
      ylab("") + xlab("") +
      thm1 +
      scale_x_discrete(labels = c("5", "9", "AICc")) +
      # scale_fill_viridis_d() + # Apply viridis color palette
      scale_fill_manual(values = custom_colors) + # Apply custom pastel colors
      theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
    
    # If the plot is the 4th to 6th in the list, modify the y-axis labels
    if (beta == 2) {
      p <- p +       
        labs(title = bquote(beta[.(beta)] ~ .("x") ~10^8 ~ .(", ") ~ 
                              .(as.character(cov)) ~ .(", ") ~ epsilon ~ "=" ~.(Eta))) +
        scale_y_log10(labels = function(x) x * 1e8)
    }else if (beta == 1) {
      p <- p +       
        labs(title = bquote(beta[.(beta)] ~ .("x") ~10^4 ~ .(", ") ~ 
                              .(as.character(cov)) ~ .(", ") ~ epsilon ~ "=" ~.(Eta))) +
        scale_y_log10(labels = function(x) x * 1e4)
    } else {
      p <- p +
        labs(title = bquote(beta[.(beta)] ~ ", " ~.(as.character(cov)) ~ .(", ") ~ epsilon ~ "=" ~.(Eta))) +
        scale_y_log10()
    }
  })
})

# Flatten the list of lists
plot_list_QuaComb_eta01 <- do.call(c, plot_list_QuaComb_eta01)

# Arrange plots using ggarrange
combined_plot_QuaComb_eta01 <- ggarrange(plotlist = plot_list_QuaComb_eta01, 
                                   ncol = 3, nrow = 3,
                                   align = "hv", 
                                   # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                   common.legend = TRUE, legend = "top")

# Add common x and y labels
combined_plot_QuaComb_eta01 <- annotate_figure(combined_plot_QuaComb_eta01,
                                         bottom = text_grob("Bandwidth", size = 20),
                                         left = text_grob("MSE", size = 20, rot = 90))
# Print the combined plot
print(combined_plot_QuaComb_eta01)
## Col_QuaCombMSE_newpar_eta01







