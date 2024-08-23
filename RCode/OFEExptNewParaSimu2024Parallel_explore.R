rm(list = ls())
library(asreml)
library(ggplot2)
library(tidyverse)

thm1 <- theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.key.width = unit(1,"cm"),
        legend.position = "top",
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size=20))


######### NS , eta = 1 ###############

MSE_ns_LR5_PC  <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_LR5_PC2024.rds")
MSE_ns_LR9_PC  <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_LR9_PC2024.rds")
MSE_ns_LRbw_PC <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_LRbw_PC2024.rds")

MSE_ns_LS5_PC  <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_LS5_PC2024.rds")
MSE_ns_LS9_PC  <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_LS9_PC2024.rds")
MSE_ns_LSbw_PC <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_LSbw_PC2024.rds")


MSE_ns_LR5_sq  <- apply(MSE_ns_LR5_PC, c(1,2), function(a) mean(a^2))
MSE_ns_LR9_sq  <- apply(MSE_ns_LR9_PC, c(1,2), function(a) mean(a^2))
MSE_ns_LRbw_sq <- apply(MSE_ns_LRbw_PC, c(1,2), function(a) mean(a^2))

MSE_ns_LS5_sq <- apply(MSE_ns_LS5_PC, c(1,2), function(a) mean(a^2))
MSE_ns_LS9_sq <- apply(MSE_ns_LS9_PC, c(1,2), function(a) mean(a^2))
MSE_ns_LSbw_sq <- apply(MSE_ns_LSbw_PC, c(1,2), function(a) mean(a^2))

M <- 1000

LinBeta0_ns <- data.frame(
  Bandwidth = c(rep("5",M),rep("9",M),rep("AICc",M),
                rep("5",M),rep("9",M),rep("AICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  MSE    = c(MSE_ns_LR5_sq[,1],MSE_ns_LR9_sq[,1],MSE_ns_LRbw_sq[,1],
             MSE_ns_LS5_sq[,1],MSE_ns_LS9_sq[,1],MSE_ns_LSbw_sq[,1])
)

LinBeta1_ns <- data.frame(
  Bandwidth = c(rep("5",M),rep("9",M),rep("AICc",M),
                rep("5",M),rep("9",M),rep("AICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  MSE    = c(MSE_ns_LR5_sq[,2],MSE_ns_LR9_sq[,2],MSE_ns_LRbw_sq[,2],
             MSE_ns_LS5_sq[,2],MSE_ns_LS9_sq[,2],MSE_ns_LSbw_sq[,2])
)


LinBeta0_ns$Cov <- "NS"
LinBeta0_ns$Beta <- 0
LinBeta1_ns$Cov <- "NS"
LinBeta1_ns$Beta <- 1

LinComb_ns <- bind_rows(LinBeta0_ns, LinBeta1_ns)
LinComb_ns <- LinComb_ns %>% mutate(Cov = factor(Cov, levels = c("NS", "AR1", "Matern")))


ggplot(LinComb_ns) + 
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")



MSE_ns_QR5_PC <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_QR5_PC2024.rds")
MSE_ns_QR9_PC <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_QR9_PC2024.rds")
MSE_ns_QRbw_PC <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_QRbw_PC2024.rds")

MSE_ns_QS5_PC <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_QS5_PC2024.rds")
MSE_ns_QS9_PC <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_QS9_PC2024.rds")
MSE_ns_QSbw_PC <-readRDS("Outcome/ExptNewPara/Simu1K2024/MSE_ns_QSbw_PC2024.rds")


MSE_ns_QR5_sq  <- apply(MSE_ns_QR5_PC, c(1,2), function(a) mean(a^2))
MSE_ns_QR9_sq  <- apply(MSE_ns_QR9_PC, c(1,2), function(a) mean(a^2))
MSE_ns_QRbw_sq <- apply(MSE_ns_QRbw_PC, c(1,2), function(a) mean(a^2))

MSE_ns_QS5_sq <- apply(MSE_ns_QS5_PC, c(1,2), function(a) mean(a^2))
MSE_ns_QS9_sq <- apply(MSE_ns_QS9_PC, c(1,2), function(a) mean(a^2))
MSE_ns_QSbw_sq <- apply(MSE_ns_QSbw_PC, c(1,2), function(a) mean(a^2))


QuaBeta0_ns <- data.frame(
  Bandwidth = c(rep("5",M),rep("9",M),rep("AICc",M),
                rep("5",M),rep("9",M),rep("AICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  MSE = c(MSE_ns_QR5_sq[,1],MSE_ns_QR9_sq[,1],MSE_ns_QRbw_sq[,1],
          MSE_ns_QS5_sq[,1],MSE_ns_QS9_sq[,1],MSE_ns_QSbw_sq[,1])
)

QuaBeta1_ns <- data.frame(
  Bandwidth = c(rep("5",M),rep("9",M),rep("AICc",M),
                rep("5",M),rep("9",M),rep("AICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  MSE = c(MSE_ns_QR5_sq[,2],MSE_ns_QR9_sq[,2],MSE_ns_QRbw_sq[,2],
          MSE_ns_QS5_sq[,2],MSE_ns_QS9_sq[,2],MSE_ns_QSbw_sq[,2])
)

QuaBeta2_ns <- data.frame(
  Bandwidth = c(rep("5",M),rep("9",M),rep("AICc",M),
                rep("5",M),rep("9",M),rep("AICc",M)),
  Design = rep(c("Randomised","Systematic"),each=M*3),
  MSE = c(MSE_ns_QR5_sq[,3],MSE_ns_QR9_sq[,3],MSE_ns_QRbw_sq[,3],
          MSE_ns_QS5_sq[,3],MSE_ns_QS9_sq[,3],MSE_ns_QSbw_sq[,3])
)


QuaBeta0_ns$Cov <- "NS"
QuaBeta0_ns$Beta <- 0
QuaBeta1_ns$Cov <- "NS"
QuaBeta1_ns$Beta <- 1
QuaBeta2_ns$Cov <- "NS"
QuaBeta2_ns$Beta <- 2

QuaComb_ns <- bind_rows(QuaBeta0_ns,QuaBeta1_ns,QuaBeta2_ns)
QuaComb_ns <- QuaComb_ns %>% mutate(Cov = factor(Cov, levels = c("NS", "AR1", "Matern")))


ggplot(QuaComb_ns) + 
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


LinComb_ns <- LinComb_ns %>% 
  mutate(Design = as.factor(Design),
         Bandwidth = as.factor(Bandwidth))

QuaComb_ns <- QuaComb_ns %>% 
  mutate(Design = as.factor(Design),
         Bandwidth = as.factor(Bandwidth))


summary(aov(MSE ~ Design * Bandwidth, LinComb_ns))
summary(aov(MSE ~ Design * Bandwidth, QuaComb_ns))



######### ar1, eta = 1 #################

parallel_ar1_pc <- readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_AR1_results_lap.rds")

length(parallel_ar1_pc)

names(parallel_ar1_pc[[1]])

dim(parallel_ar1_pc[[1]]$MSE_ar1_LS5)

ar1_lin_var_names <- c("MSE_ar1_LS5", "MSE_ar1_LS9", "MSE_ar1_LSbw", 
                       "MSE_ar1_LR5", "MSE_ar1_LR9", "MSE_ar1_LRbw")

# Initialize an empty list to store the results
ar1_lin_results <- list()
# Loop over each variable name
for (var in ar1_lin_var_names) {
  # Initialize two empty vectors for each variable
  B0 <- NULL
  B1 <- NULL
  # Loop over each element in parallel_ar1_mac
  for(i in 1:length(parallel_ar1_pc)){
    for(j in 1:dim(parallel_ar1_pc[[1]]$MSE_ar1_LS5)[1]){
      B0 <- c(B0, mean((parallel_ar1_pc[[i]][[var]][j,1,])^2))
      B1 <- c(B1, mean((parallel_ar1_pc[[i]][[var]][j,2,])^2))
    }
  }
  
  # Store the results in the list
  ar1_lin_results[[paste0(var, "_B0")]] <- B0
  ar1_lin_results[[paste0(var, "_B1")]] <- B1
}


names(ar1_lin_results)

length(ar1_lin_results$MSE_ar1_LS5_B0)


# Function to create a data frame
create_df <- function(var_name, results) {
  design <- substr(var_name, 10, 10)
  bandwidth <- substr(var_name, 11, ifelse(nchar(var_name) == 14, 11, 12))
  beta <- substr(var_name, nchar(var_name)-1, nchar(var_name))
  
  df <- data.frame(MSE = results[[var_name]],
                   Design = design,
                   Bandwidth = ifelse(bandwidth == "bw", "bw", as.numeric(bandwidth)),
                   Beta = beta)
  return(df)
}


# List of variable names
ar1_lin_var_names <- names(ar1_lin_results)

ar1_lin_df_list <- list()

# Loop over each variable name
for (var in ar1_lin_var_names) {
  ar1_lin_df_list[[var]] <- create_df(var, ar1_lin_results)
}

# Combine all data frames into one
ar1_lin_df_list <- do.call(rbind, ar1_lin_df_list)
ar1_lin_df_list <- ar1_lin_df_list %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "AR1"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )


str(ar1_lin_df_list)


ggplot(ar1_lin_df_list) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


summary(aov(MSE ~ Design * Bandwidth, ar1_lin_df_list))



ar1_qua_var_names <- c("MSE_ar1_QS5", "MSE_ar1_QS9", "MSE_ar1_QSbw", 
                       "MSE_ar1_QR5", "MSE_ar1_QR9", "MSE_ar1_QRbw")

# Initialize an empty list to store the results
ar1_qua_results <- list()
for (var in ar1_qua_var_names) {
  B0 <- NULL
  B1 <- NULL
  B2 <- NULL

  for(i in 1:length(parallel_ar1_pc)){
    for(j in 1:dim(parallel_ar1_pc[[1]]$MSE_ar1_QS5)[1]){
      B0 <- c(B0, mean((parallel_ar1_pc[[i]][[var]][j,1,])^2))
      B1 <- c(B1, mean((parallel_ar1_pc[[i]][[var]][j,2,])^2))
      B2 <- c(B2, mean((parallel_ar1_pc[[i]][[var]][j,3,])^2))
    }
  }
  
  ar1_qua_results[[paste0(var, "_B0")]] <- B0
  ar1_qua_results[[paste0(var, "_B1")]] <- B1
  ar1_qua_results[[paste0(var, "_B2")]] <- B2
}

names(ar1_qua_results)

length(ar1_qua_results$MSE_ar1_QS5_B0)

ar1_qua_var_names <- names(ar1_qua_results)

ar1_qua_df_list <- list()
for (var in ar1_qua_var_names) {
  ar1_qua_df_list[[var]] <- create_df(var, ar1_qua_results)
}

ar1_qua_df_list <- do.call(rbind, ar1_qua_df_list)
ar1_qua_df_list <- ar1_qua_df_list %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "AR1"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )


str(ar1_qua_df_list)

ggplot(ar1_qua_df_list) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


summary(aov(MSE ~ Design * Bandwidth, ar1_qua_df_list))



########### mat mac, eta =1 ############

parallel_mat_mac <-readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_mat_results_mac.rds")

length(parallel_mat_mac)

dim(parallel_mat_mac[[1]]$MSE_mat_LS5[1,,])

names(parallel_mat_mac[[1]])
# [1] "MSE_mat_LS5"  "MSE_mat_LS9"  "MSE_mat_LSbw" "BW_mat_LS"    "MSE_mat_LR5" 
# [6] "MSE_mat_LR9"  "MSE_mat_LRbw" "BW_mat_LR"    "MSE_mat_QS5"  "MSE_mat_QS9" 
# [11] "MSE_mat_QSbw" "BW_mat_QS"    "MSE_mat_QR5"  "MSE_mat_QR9"  "MSE_mat_QRbw"
# [16] "BW_mat_QR"   

mat_lin_var_names <- c("MSE_mat_LS5", "MSE_mat_LS9", "MSE_mat_LSbw", 
                       "MSE_mat_LR5", "MSE_mat_LR9", "MSE_mat_LRbw")

mat_lin_results <- list()
for (var in mat_lin_var_names) {
  B0 <- NULL
  B1 <- NULL
  for(i in 1:1000){
    B0 <- c(B0, mean((parallel_mat_mac[[i]][[var]][1,1,])^2))
    B1 <- c(B1, mean((parallel_mat_mac[[i]][[var]][1,2,])^2))
  }
  mat_lin_results[[paste0(var, "_B0")]] <- B0
  mat_lin_results[[paste0(var, "_B1")]] <- B1
}

mat_qua_var_names <- c("MSE_mat_QS5", "MSE_mat_QS9", "MSE_mat_QSbw", 
                       "MSE_mat_QR5", "MSE_mat_QR9", "MSE_mat_QRbw")
mat_qua_results <- list()

for (var in mat_qua_var_names) {
  B0 <- NULL
  B1 <- NULL
  B2 <- NULL
  for(i in 1:1000){
    B0 <- c(B0, mean((parallel_mat_mac[[i]][[var]][1,1,])^2))
    B1 <- c(B1, mean((parallel_mat_mac[[i]][[var]][1,2,])^2))
    B2 <- c(B2, mean((parallel_mat_mac[[i]][[var]][1,3,])^2))
  }
  mat_qua_results[[paste0(var, "_B0")]] <- B0
  mat_qua_results[[paste0(var, "_B1")]] <- B1
  mat_qua_results[[paste0(var, "_B2")]] <- B2
}

names(mat_qua_results)
names(mat_lin_results)

create_df <- function(var_name, results) {
  design <- substr(var_name, 10, 10)
  bandwidth <- substr(var_name, 11, ifelse(nchar(var_name) == 14, 11, 12))
  beta <- substr(var_name, nchar(var_name)-1, nchar(var_name))
  df <- data.frame(MSE = results[[var_name]],
                   Design = design,
                   Bandwidth = ifelse(bandwidth == "bw", "bw", as.numeric(bandwidth)),
                   Beta = beta)
  return(df)
}

mat_lin_var_names <- c("MSE_mat_LS5_B0", "MSE_mat_LS5_B1", "MSE_mat_LS9_B0", "MSE_mat_LS9_B1", 
               "MSE_mat_LSbw_B0", "MSE_mat_LSbw_B1", "MSE_mat_LR5_B0", "MSE_mat_LR5_B1", 
               "MSE_mat_LR9_B0", "MSE_mat_LR9_B1", "MSE_mat_LRbw_B0", "MSE_mat_LRbw_B1")
mat_lin_df_list <- list()
for (var in mat_lin_var_names) {
  mat_lin_df_list[[var]] <- create_df(var, mat_lin_results)
}
mat_lin_df_list <- do.call(rbind, mat_lin_df_list)

mat_lin_df_list <- mat_lin_df_list %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "Matern"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )

str(mat_lin_df_list)

ggplot(mat_lin_df_list) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


names(mat_qua_results)

# [1] "MSE_mat_QS5_B0"  "MSE_mat_QS5_B1"  "MSE_mat_QS5_B2"  "MSE_mat_QS9_B0" 
# [5] "MSE_mat_QS9_B1"  "MSE_mat_QS9_B2"  "MSE_mat_QSbw_B0" "MSE_mat_QSbw_B1"
# [9] "MSE_mat_QSbw_B2" "MSE_mat_QR5_B0"  "MSE_mat_QR5_B1"  "MSE_mat_QR5_B2" 
# [13] "MSE_mat_QR9_B0"  "MSE_mat_QR9_B1"  "MSE_mat_QR9_B2"  "MSE_mat_QRbw_B0"
# [17] "MSE_mat_QRbw_B1" "MSE_mat_QRbw_B2"

mat_qua_df_list <- list()
for (var in names(mat_qua_results)) {
  mat_qua_df_list[[var]] <- create_df(var, mat_qua_results)
}
mat_qua_df_list <- do.call(rbind, mat_qua_df_list)
mat_qua_df_list <- mat_qua_df_list %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "Matern"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )

str(mat_qua_df_list)

ggplot(mat_qua_df_list) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


summary(aov(MSE ~ Design * Bandwidth,mat_lin_df_list))
summary(aov(MSE ~ Design * Bandwidth,mat_qua_df_list))



####### combind eta = 1 all data frames, output ##############

combined_lin_df <- bind_rows(LinComb_ns, ar1_lin_df_list, mat_lin_df_list)
str(combined_lin_df)

summary(aov(MSE ~ Design * Bandwidth * Cov, combined_lin_df))

combined_qua_df <- bind_rows(QuaComb_ns, ar1_qua_df_list, mat_qua_df_list)
str(combined_qua_df)

summary(aov(MSE ~ Design * Bandwidth * Cov,combined_qua_df))

ggplot(combined_lin_df) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(Beta~Cov, scales = "free_y")

ggplot(combined_qua_df) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(Beta~Cov, scales = "free_y")


############# plots for eta = 1 #############

library(ggpubr)
library(viridis)

custom_colors <- c("Randomised" = "gray70", "Systematic" = "white")

plot_lin_list <- lapply(unique(combined_lin_df$Beta), function(beta) {
  lapply(sort(unique(combined_lin_df$Cov)), function(cov) {
    data_subset <- subset(combined_lin_df, Beta == beta & Cov == as.character(cov))
    
    # Adjust MSE values if Beta == 1
    if (beta == 1) {
      data_subset$MSE <- data_subset$MSE * 10^4
      beta_label <- bquote(epsilon == 1 * ", " * beta[.(beta)] * " x " * 10^4)
    } else {
      beta_label <- bquote(epsilon == 1 * ", " * beta[.(beta)])
    }
    
    ggplot(data_subset, aes(Bandwidth, MSE)) +
      geom_boxplot(aes(fill = Design)) +
      scale_y_log10() +
      ylab("") + xlab("") +
      thm1 +
      scale_x_discrete(labels = c("5", "9", "AICc")) +
      labs(title = bquote(.(beta_label) * ", " * .(as.character(cov)))) +
      scale_fill_manual(values = custom_colors) + # Apply custom pastel colors
      theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
  })
})


plot_lin_list <- do.call(c, plot_lin_list)
combined_lin_plot <- ggarrange(plotlist = plot_lin_list, 
                               ncol = 3, nrow = 2,
                               align = "hv", 
                               # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                               common.legend = TRUE, legend = "top")
combined_lin_plot <- annotate_figure(combined_lin_plot,
                                     bottom = text_grob("Bandwidth", size = 20),
                                     left = text_grob("MSE", size = 20, rot = 90))

print(combined_lin_plot)


plot_qua_list <- lapply(unique(combined_qua_df$Beta), function(beta) {
  lapply(sort(unique(combined_qua_df$Cov)), function(cov) {
    data_subset <- subset(combined_qua_df, Beta == beta & Cov == as.character(cov))
    
    # Adjust MSE values if Beta == 1
    if (beta == 2) {
      data_subset$MSE <- data_subset$MSE * 10^8
      beta_label <- bquote(epsilon == 1 * ", " * beta[.(beta)] * " x " * 10^8)
    } else if (beta == 1) {
      data_subset$MSE <- data_subset$MSE * 10^4
      beta_label <- bquote(epsilon == 1 * ", " * beta[.(beta)] * " x " * 10^4)
    } else {
      beta_label <- bquote(epsilon == 1 * ", " * beta[.(beta)])
    }
    
    ggplot(data_subset, aes(Bandwidth, MSE)) +
      geom_boxplot(aes(fill = Design)) +
      scale_y_log10() +
      ylab("") + xlab("") +
      thm1 +
      scale_x_discrete(labels = c("5", "9", "AICc")) +
      labs(title = bquote(.(beta_label) * ", " * .(as.character(cov)))) +
      scale_fill_manual(values = custom_colors) + # Apply custom pastel colors
      theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
  })
})

plot_qua_list <- do.call(c, plot_qua_list)
combined_qua_plot <- ggarrange(plotlist = plot_qua_list, 
                           ncol = 3, nrow = 3,
                           align = "hv", 
                           # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                           common.legend = TRUE, legend = "top")
combined_qua_plot <- annotate_figure(combined_qua_plot,
                                 bottom = text_grob("Bandwidth", size = 20),
                                 left = text_grob("MSE", size = 20, rot = 90))

print(combined_qua_plot)



summarized_lin_df <- combined_lin_df %>%
  group_by(Design, Cov, Beta, Bandwidth) %>%
  summarise(median_MSE = median(MSE), .groups = 'drop')

adjusted_lin_df <- summarized_lin_df %>%
  mutate(
    median_MSE = case_when(
      Beta == 1 ~ median_MSE * 10^4,
      TRUE ~ median_MSE
    ),
    median_MSE = median_MSE
  )

pivot_lin_df <- adjusted_lin_df %>%
  pivot_wider(
    names_from = c(Design, Bandwidth),
    values_from = median_MSE
  )

pivot_lin_df <- pivot_lin_df %>%
  arrange(Cov, Beta)

print(pivot_lin_df)
xtable::xtable(pivot_lin_df,digits=3)


summarized_qua_df <- combined_qua_df %>%
  group_by(Design, Cov, Beta, Bandwidth) %>%
  summarise(median_MSE = median(MSE), .groups = 'drop')
adjusted_qua_df <- summarized_qua_df %>%
  mutate(
    median_MSE = case_when(
      Beta == 1 ~ median_MSE * 10^4,
      Beta == 2 ~ median_MSE * 10^8,
      TRUE ~ median_MSE
    ),
    median_MSE = median_MSE
  )

pivot_qua_df <- adjusted_qua_df %>%
  pivot_wider(
    names_from = c(Design, Bandwidth),
    values_from = median_MSE
  )

pivot_qua_df <- pivot_qua_df %>%
  arrange(Cov, Beta)

print(pivot_qua_df)

xtable::xtable(pivot_qua_df,digits=3)
# xtable::xtable(pivot_qua_df,digits=4)





######### NS , eta = 0.1 ###############

parallel_ns_pawsey_eta01 <- readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_NS_results_eta01_pawsey.rds")

length(parallel_ns_pawsey_eta01)

names(parallel_ns_pawsey_eta01[[1]])

dim(parallel_ns_pawsey_eta01[[1]]$MSE_ns_LS5)

ns_lin_var_names <- c("MSE_ns_LS5", "MSE_ns_LS9", "MSE_ns_LSbw", 
                       "MSE_ns_LR5", "MSE_ns_LR9", "MSE_ns_LRbw")

ns_lin_results_eta01 <- list()
for (var in ns_lin_var_names) {
  # Initialize two empty vectors for each variable
  B0 <- NULL
  B1 <- NULL
  # Loop over each element in parallel_ns_mac
  for(i in 1:length(parallel_ns_pawsey_eta01)){
    for(j in 1:dim(parallel_ns_pawsey_eta01[[1]]$MSE_ns_LS5)[1]){
      B0 <- c(B0, mean((parallel_ns_pawsey_eta01[[i]][[var]][j,1,])^2))
      B1 <- c(B1, mean((parallel_ns_pawsey_eta01[[i]][[var]][j,2,])^2))
    }
  }
  
  # Store the results in the list
  ns_lin_results_eta01[[paste0(var, "_B0")]] <- B0
  ns_lin_results_eta01[[paste0(var, "_B1")]] <- B1
}


names(ns_lin_results_eta01)

length(ns_lin_results_eta01$MSE_ns_LS5_B0)

create_df_ns <- function(var_name, results) {
  # Extract the design, bandwidth, and beta from the variable name
  design <- substr(var_name, 9, 9)
  bandwidth <- substr(var_name, 10, ifelse(nchar(var_name) == 13, 10, 11))
  beta <- substr(var_name, nchar(var_name)-1, nchar(var_name))
  
  # Create the data frame
  df <- data.frame(MSE = results[[var_name]],
                   Design = design,
                   Bandwidth = ifelse(bandwidth == "bw", "bw", as.numeric(bandwidth)),
                   Beta = beta)
  
  return(df)
}


ns_lin_var_names <- names(ns_lin_results_eta01)
ns_lin_df_list_eta01 <- list()

for (var in ns_lin_var_names) {
  ns_lin_df_list_eta01[[var]] <- create_df_ns(var, ns_lin_results_eta01)
}

ns_lin_df_list_eta01 <- do.call(rbind, ns_lin_df_list_eta01)
ns_lin_df_list_eta01 <- ns_lin_df_list_eta01 %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "NS"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )

str(ns_lin_df_list_eta01)

ggplot(ns_lin_df_list_eta01) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")

summary(aov(MSE ~ Design * Bandwidth, ns_lin_df_list_eta01))




ns_qua_var_names <- c("MSE_ns_QS5", "MSE_ns_QS9", "MSE_ns_QSbw", 
                       "MSE_ns_QR5", "MSE_ns_QR9", "MSE_ns_QRbw")

ns_qua_results_eta01 <- list()
for (var in ns_qua_var_names) {
  # Initialize two empty vectors for each variable
  B0 <- NULL
  B1 <- NULL
  B2 <- NULL
  
  # Loop over each element in parallel_ns_mac
  for(i in 1:length(parallel_ns_pawsey_eta01)){
    for(j in 1:dim(parallel_ns_pawsey_eta01[[1]]$MSE_ns_QS5)[1]){
      B0 <- c(B0, mean((parallel_ns_pawsey_eta01[[i]][[var]][j,1,])^2))
      B1 <- c(B1, mean((parallel_ns_pawsey_eta01[[i]][[var]][j,2,])^2))
      B2 <- c(B2, mean((parallel_ns_pawsey_eta01[[i]][[var]][j,3,])^2))
    }
  }
  ns_qua_results_eta01[[paste0(var, "_B0")]] <- B0
  ns_qua_results_eta01[[paste0(var, "_B1")]] <- B1
  ns_qua_results_eta01[[paste0(var, "_B2")]] <- B2
}

names(ns_qua_results_eta01)

length(ns_qua_results_eta01$MSE_ns_QS5_B0)

ns_qua_var_names <- names(ns_qua_results_eta01)

ns_qua_df_list_eta01 <- list()

for (var in ns_qua_var_names) {
  ns_qua_df_list_eta01[[var]] <- create_df_ns(var, ns_qua_results_eta01)
}

ns_qua_df_list_eta01 <- do.call(rbind, ns_qua_df_list_eta01)

ns_qua_df_list_eta01 <- ns_qua_df_list_eta01 %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "NS"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )


str(ns_qua_df_list_eta01)


ggplot(ns_qua_df_list_eta01) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


summary(aov(MSE ~ Design * Bandwidth, ns_qua_df_list_eta01))


######### ar1, eta = 0.1 #################

parallel_ar1_mac_eta01 <- readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_AR1_results_eta01_mac.rds")
length(parallel_ar1_mac_eta01)

names(parallel_ar1_mac_eta01[[1]])

dim(parallel_ar1_mac_eta01[[1]]$MSE_ar1_LS5)

ar1_lin_var_names <- c("MSE_ar1_LS5", "MSE_ar1_LS9", "MSE_ar1_LSbw", 
                       "MSE_ar1_LR5", "MSE_ar1_LR9", "MSE_ar1_LRbw")

ar1_lin_results_eta01 <- list()
for (var in ar1_lin_var_names) {
  B0 <- NULL
  B1 <- NULL
  for(i in 1:length(parallel_ar1_mac_eta01)){
    for(j in 1:dim(parallel_ar1_mac_eta01[[1]]$MSE_ar1_LS5)[1]){
      B0 <- c(B0, mean((parallel_ar1_mac_eta01[[i]][[var]][j,1,])^2))
      B1 <- c(B1, mean((parallel_ar1_mac_eta01[[i]][[var]][j,2,])^2))
    }
  }
  ar1_lin_results_eta01[[paste0(var, "_B0")]] <- B0
  ar1_lin_results_eta01[[paste0(var, "_B1")]] <- B1
}


names(ar1_lin_results_eta01)

length(ar1_lin_results_eta01$MSE_ar1_LS5_B0)

create_df <- function(var_name, results) {
  # Extract the design, bandwidth, and beta from the variable name
  design <- substr(var_name, 10, 10)
  bandwidth <- substr(var_name, 11, ifelse(nchar(var_name) == 14, 11, 12))
  beta <- substr(var_name, nchar(var_name)-1, nchar(var_name))
  
  # Create the data frame
  df <- data.frame(MSE = results[[var_name]],
                   Design = design,
                   Bandwidth = ifelse(bandwidth == "bw", "bw", as.numeric(bandwidth)),
                   Beta = beta)
  
  return(df)
}

ar1_lin_var_names <- names(ar1_lin_results_eta01)
ar1_lin_df_list_eta01 <- list()
for (var in ar1_lin_var_names) {
  ar1_lin_df_list_eta01[[var]] <- create_df(var, ar1_lin_results_eta01)
}

ar1_lin_df_list_eta01 <- do.call(rbind, ar1_lin_df_list_eta01)
ar1_lin_df_list_eta01 <- ar1_lin_df_list_eta01 %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "AR1"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )

str(ar1_lin_df_list_eta01)

ggplot(ar1_lin_df_list_eta01) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


summary(aov(MSE ~ Design * Bandwidth, ar1_lin_df_list_eta01))


ar1_qua_var_names <- c("MSE_ar1_QS5", "MSE_ar1_QS9", "MSE_ar1_QSbw", 
                       "MSE_ar1_QR5", "MSE_ar1_QR9", "MSE_ar1_QRbw")

ar1_qua_results_eta01 <- list()
for (var in ar1_qua_var_names) {
  B0 <- NULL
  B1 <- NULL
  B2 <- NULL
  for(i in 1:length(parallel_ar1_mac_eta01)){
    for(j in 1:dim(parallel_ar1_mac_eta01[[1]]$MSE_ar1_QS5)[1]){
      B0 <- c(B0, mean((parallel_ar1_mac_eta01[[i]][[var]][j,1,])^2))
      B1 <- c(B1, mean((parallel_ar1_mac_eta01[[i]][[var]][j,2,])^2))
      B2 <- c(B2, mean((parallel_ar1_mac_eta01[[i]][[var]][j,3,])^2))
    }
  }
  ar1_qua_results_eta01[[paste0(var, "_B0")]] <- B0
  ar1_qua_results_eta01[[paste0(var, "_B1")]] <- B1
  ar1_qua_results_eta01[[paste0(var, "_B2")]] <- B2
}

names(ar1_qua_results_eta01)

length(ar1_qua_results_eta01$MSE_ar1_QS5_B0)

ar1_qua_var_names <- names(ar1_qua_results_eta01)
ar1_qua_df_list_eta01 <- list()

for (var in ar1_qua_var_names) {
  ar1_qua_df_list_eta01[[var]] <- create_df(var, ar1_qua_results_eta01)
}

ar1_qua_df_list_eta01 <- do.call(rbind, ar1_qua_df_list_eta01)
ar1_qua_df_list_eta01 <- ar1_qua_df_list_eta01 %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "AR1"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )


str(ar1_qua_df_list_eta01)
ggplot(ar1_qua_df_list_eta01) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


summary(aov(MSE ~ Design * Bandwidth, ar1_qua_df_list_eta01))


########### mat mac eta = 0.1 ############

parallel_mat_mac_eta01 <-readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_mat_results_mac_eta01.rds")

length(parallel_mat_mac_eta01)

names(parallel_mat_mac_eta01[[1]])
# [1] "MSE_mat_LS5"  "MSE_mat_LS9"  "MSE_mat_LSbw" "BW_mat_LS"    "MSE_mat_LR5" 
# [6] "MSE_mat_LR9"  "MSE_mat_LRbw" "BW_mat_LR"    "MSE_mat_QS5"  "MSE_mat_QS9" 
# [11] "MSE_mat_QSbw" "BW_mat_QS"    "MSE_mat_QR5"  "MSE_mat_QR9"  "MSE_mat_QRbw"
# [16] "BW_mat_QR"   

# List of variable names
mat_lin_var_names <- c("MSE_mat_LS5", "MSE_mat_LS9", "MSE_mat_LSbw", 
                       "MSE_mat_LR5", "MSE_mat_LR9", "MSE_mat_LRbw")

# Initialize an empty list to store the results
mat_lin_results_eta01 <- list()
# Loop over each variable name
for (var in mat_lin_var_names) {
  # Initialize two empty vectors for each variable
  B0 <- NULL
  B1 <- NULL
  # Loop over each element in parallel_mat_mac
  for(i in 1:length(parallel_mat_mac_eta01)){
    for(j in 1:dim(parallel_mat_mac_eta01[[1]]$MSE_mat_LS5)[1]){
      B0 <- c(B0, mean((parallel_mat_mac_eta01[[i]][[var]][j,1,])^2))
      B1 <- c(B1, mean((parallel_mat_mac_eta01[[i]][[var]][j,2,])^2))
    }
  }
  
  # Store the results in the list
  mat_lin_results_eta01[[paste0(var, "_B0")]] <- B0
  mat_lin_results_eta01[[paste0(var, "_B1")]] <- B1
}

mat_qua_var_names <- c("MSE_mat_QS5", "MSE_mat_QS9", "MSE_mat_QSbw", 
                       "MSE_mat_QR5", "MSE_mat_QR9", "MSE_mat_QRbw")
mat_qua_results_eta01 <- list()

for (var in mat_qua_var_names) {
  B0 <- NULL
  B1 <- NULL
  B2 <- NULL
  for(i in 1:length(parallel_mat_mac_eta01)){
    for(j in 1:dim(parallel_mat_mac_eta01[[1]]$MSE_mat_LS5)[1]){
      B0 <- c(B0, mean((parallel_mat_mac_eta01[[i]][[var]][1,1,])^2))
      B1 <- c(B1, mean((parallel_mat_mac_eta01[[i]][[var]][1,2,])^2))
      B2 <- c(B2, mean((parallel_mat_mac_eta01[[i]][[var]][1,3,])^2))
    }
  }
  mat_qua_results_eta01[[paste0(var, "_B0")]] <- B0
  mat_qua_results_eta01[[paste0(var, "_B1")]] <- B1
  mat_qua_results_eta01[[paste0(var, "_B2")]] <- B2
}

create_df <- function(var_name, results) {
  # Extract the design, bandwidth, and beta from the variable name
  design <- substr(var_name, 10, 10)
  bandwidth <- substr(var_name, 11, ifelse(nchar(var_name) == 14, 11, 12))
  beta <- substr(var_name, nchar(var_name)-1, nchar(var_name))
  
  # Create the data frame
  df <- data.frame(MSE = results[[var_name]],
                   Design = design,
                   Bandwidth = ifelse(bandwidth == "bw", "bw", as.numeric(bandwidth)),
                   Beta = beta)
  
  return(df)
}

mat_lin_var_names <- c("MSE_mat_LS5_B0", "MSE_mat_LS5_B1", "MSE_mat_LS9_B0", "MSE_mat_LS9_B1", 
                       "MSE_mat_LSbw_B0", "MSE_mat_LSbw_B1", "MSE_mat_LR5_B0", "MSE_mat_LR5_B1", 
                       "MSE_mat_LR9_B0", "MSE_mat_LR9_B1", "MSE_mat_LRbw_B0", "MSE_mat_LRbw_B1")

mat_lin_df_list_eta01 <- list()

for (var in mat_lin_var_names) {
  mat_lin_df_list_eta01[[var]] <- create_df(var, mat_lin_results_eta01)
}

mat_lin_df_list_eta01 <- do.call(rbind, mat_lin_df_list_eta01)

mat_lin_df_list_eta01 <- mat_lin_df_list_eta01 %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "Matern"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )

str(mat_lin_df_list_eta01)

ggplot(mat_lin_df_list_eta01) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


mat_qua_df_list_eta01 <- list()

for (var in names(mat_qua_results_eta01)) {
  mat_qua_df_list_eta01[[var]] <- create_df(var, mat_qua_results_eta01)
}
mat_qua_df_list_eta01 <- do.call(rbind, mat_qua_df_list_eta01)

mat_qua_df_list_eta01 <- mat_qua_df_list_eta01 %>%
  mutate(
    Design = recode(Design, "S" = "Systematic", "R" = "Randomised"),
    Bandwidth = recode(Bandwidth, "bw" = "AICc"),
    Beta = str_remove(Beta, "B"),
    Cov = "Matern"
  ) %>% 
  mutate(
    Design = as.factor(Design),
    Bandwidth = as.factor(Bandwidth),
    Beta = as.numeric(Beta),
    Cov = as.factor(Cov)
  )

str(mat_qua_df_list_eta01)

ggplot(mat_qua_df_list_eta01) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(~Beta, scales = "free_y")


summary(aov(MSE ~ Design * Bandwidth,mat_lin_df_list_eta01))
summary(aov(MSE ~ Design * Bandwidth,mat_qua_df_list_eta01))



####### combind eta = 0.1 all data frames, output ##############

combined_lin_df_eta01 <- bind_rows(ns_lin_df_list_eta01, 
                                   ar1_lin_df_list_eta01, 
                                   mat_lin_df_list_eta01)

str(combined_lin_df_eta01)

summary(aov(MSE ~ Design * Bandwidth * Cov, combined_lin_df_eta01))


combined_qua_df_eta01 <- bind_rows(ns_qua_df_list_eta01, 
                                   ar1_qua_df_list_eta01, 
                                   mat_qua_df_list_eta01)
str(combined_qua_df_eta01)

summary(aov(MSE ~ Design * Bandwidth * Cov,combined_qua_df_eta01))


ggplot(combined_lin_df_eta01) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(Beta~Cov, scales = "free_y")


ggplot(combined_qua_df_eta01) +
  geom_boxplot(aes(Bandwidth,MSE,fill=Design)) +
  scale_y_log10() + 
  facet_wrap(Beta~Cov, scales = "free_y")



library(ggpubr)
library(viridis)

custom_colors <- c("Randomised" = "gray70", "Systematic" = "white")

plot_lin_list_eta01 <- lapply(unique(combined_lin_df_eta01$Beta), function(beta) {
  lapply(sort(unique(combined_lin_df_eta01$Cov)), function(cov) {
    data_subset <- subset(combined_lin_df_eta01, Beta == beta & Cov == as.character(cov))
    
    # Adjust MSE values if Beta == 1
    if (beta == 1) {
      data_subset$MSE <- data_subset$MSE * 10^4
      beta_label <- bquote(epsilon == 0.1 * ", " * beta[.(beta)] * " x " * 10^4)
    } else {
      beta_label <- bquote(epsilon == 0.1 * ", " * beta[.(beta)])
    }
    
    ggplot(data_subset, aes(Bandwidth, MSE)) +
      geom_boxplot(aes(fill = Design)) +
      scale_y_log10() +
      ylab("") + xlab("") +
      thm1 +
      scale_x_discrete(labels = c("5", "9", "AICc")) +
      labs(title = bquote(.(beta_label) * ", " * .(as.character(cov)))) +
      scale_fill_manual(values = custom_colors) + # Apply custom pastel colors
      theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
  })
})


plot_lin_list_eta01 <- do.call(c, plot_lin_list_eta01)
combined_lin_plot_eta01 <- ggarrange(plotlist = plot_lin_list_eta01, 
                               ncol = 3, nrow = 2,
                               align = "hv", 
                               # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                               common.legend = TRUE, legend = "top")
combined_lin_plot_eta01 <- annotate_figure(combined_lin_plot_eta01,
                                     bottom = text_grob("Bandwidth", size = 20),
                                     left = text_grob("MSE", size = 20, rot = 90))

print(combined_lin_plot_eta01)


plot_qua_list_eta01 <- lapply(unique(combined_qua_df_eta01$Beta), function(beta) {
  lapply(sort(unique(combined_qua_df_eta01$Cov)), function(cov) {
    data_subset <- subset(combined_qua_df_eta01, Beta == beta & Cov == as.character(cov))
    
    # Adjust MSE values if Beta == 1
    if (beta == 2) {
      data_subset$MSE <- data_subset$MSE * 10^8
      beta_label <- bquote(epsilon == 0.1 * ", " * beta[.(beta)] * " x " * 10^8)
    } else if (beta == 1) {
      data_subset$MSE <- data_subset$MSE * 10^4
      beta_label <- bquote(epsilon == 0.1 * ", " * beta[.(beta)] * " x " * 10^4)
    } else {
      beta_label <- bquote(epsilon == 0.1 * ", " * beta[.(beta)])
    }
    
    ggplot(data_subset, aes(Bandwidth, MSE)) +
      geom_boxplot(aes(fill = Design)) +
      scale_y_log10() +
      ylab("") + xlab("") +
      thm1 +
      scale_x_discrete(labels = c("5", "9", "AICc")) +
      labs(title = bquote(.(beta_label) * ", " * .(as.character(cov)))) +
      scale_fill_manual(values = custom_colors) + # Apply custom pastel colors
      theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) 
  })
})

plot_qua_list_eta01 <- do.call(c, plot_qua_list_eta01)
combined_qua_plot_eta01 <- ggarrange(plotlist = plot_qua_list_eta01, 
                               ncol = 3, nrow = 3,
                               align = "hv", 
                               # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                               common.legend = TRUE, legend = "top")
combined_qua_plot_eta01 <- annotate_figure(combined_qua_plot_eta01,
                                     bottom = text_grob("Bandwidth", size = 20),
                                     left = text_grob("MSE", size = 20, rot = 90))

print(combined_qua_plot_eta01)



summarized_lin_df_eta01 <- combined_lin_df_eta01 %>%
  group_by(Design, Cov, Beta, Bandwidth) %>%
  summarise(median_MSE = median(MSE), .groups = 'drop')

adjusted_lin_df_eta01 <- summarized_lin_df_eta01 %>%
  mutate(
    median_MSE = case_when(
      Beta == 1 ~ median_MSE * 10^4,
      TRUE ~ median_MSE
    ),
    median_MSE = median_MSE
  )

pivot_lin_df_eta01 <- adjusted_lin_df_eta01 %>%
  pivot_wider(
    names_from = c(Design, Bandwidth),
    values_from = median_MSE
  )

pivot_lin_df_eta01 <- pivot_lin_df_eta01 %>%
  arrange(Cov, Beta)

print(pivot_lin_df_eta01)
xtable::xtable(pivot_lin_df_eta01,digits=3)




summarized_qua_df_eta01 <- combined_qua_df_eta01 %>%
  group_by(Design, Cov, Beta, Bandwidth) %>%
  summarise(median_MSE = median(MSE), .groups = 'drop')

adjusted_qua_df_eta01 <- summarized_qua_df_eta01 %>%
  mutate(
    median_MSE = case_when(
      Beta == 1 ~ median_MSE * 10^4,
      Beta == 2 ~ median_MSE * 10^8,
      TRUE ~ median_MSE
    ),
    median_MSE = median_MSE
  )

pivot_qua_df_eta01 <- adjusted_qua_df_eta01 %>%
  pivot_wider(
    names_from = c(Design, Bandwidth),
    values_from = median_MSE
  )

pivot_qua_df_eta01 <- pivot_qua_df_eta01 %>%
  arrange(Cov, Beta)

print(pivot_qua_df_eta01)

xtable::xtable(pivot_qua_df_eta01,digits=3)







pivot_lin_df %>%
  rowwise() %>% 
  summarise(Min_Value = format(round(min(c_across(Randomised_5:Systematic_AICc)), 3), 
                               nsmall = 3))


pivot_lin_df_eta01 %>%
  rowwise() %>% 
  summarise(Min_Value = format(round(min(c_across(Randomised_5:Systematic_AICc)), 3), 
                               nsmall = 3))



pivot_qua_df %>%
  rowwise() %>% 
  summarise(Min_Value = format(round(min(c_across(Randomised_5:Systematic_AICc)), 3), 
                               nsmall = 3))


pivot_qua_df_eta01 %>%
  rowwise() %>% 
  summarise(Min_Value = format(round(min(c_across(Randomised_5:Systematic_AICc)), 3), 
                               nsmall = 3))




###### compare eta ##############


names(combined_lin_df)
names(combined_lin_df_eta01)

combined_lin_df <- combined_lin_df %>%
  mutate(Eta = "low")

combined_lin_df_eta01 <- combined_lin_df_eta01 %>%
  mutate(Eta = "high")

Lin_df <- bind_rows(combined_lin_df, combined_lin_df_eta01)

Lin_df <- Lin_df %>%
  mutate(Eta = as.factor(Eta),
         FacBeta = as.factor(Beta))

str(Lin_df)


summary(aov(MSE ~ Design * Bandwidth * Cov * Eta, Lin_df))
summary(aov(MSE ~ Design * Bandwidth * Cov * Eta * FacBeta, Lin_df))


combined_qua_df <- combined_qua_df %>%
  mutate(Eta = "low")
combined_qua_df_eta01 <- combined_qua_df_eta01 %>%
  mutate(Eta = "high")

Qua_df <- bind_rows(combined_qua_df, combined_qua_df_eta01)

Qua_df <- Qua_df %>%
  mutate(Eta = as.factor(Eta),
         FacBeta = as.factor(Beta))

str(Qua_df)


summary(aov(MSE ~ Design * Bandwidth * Cov * Eta, Qua_df))
summary(aov(MSE ~ Design * Bandwidth * Cov * Eta * FacBeta, Qua_df))

summary(aov(MSE ~ Design * Bandwidth * Cov , Lin_df))
summary(aov(MSE ~ Design * Bandwidth * Cov , Qua_df))


Lin_df <- read.table('Outcome/ExptNewPara/Simu1K2024/Lin_df.txt',header = TRUE)
Qua_df <- read.table('Outcome/ExptNewPara/Simu1K2024/Qua_df.txt',header = TRUE)


summary(aov(MSE ~ Design * Bandwidth * Cov * Eta, Lin_df))
summary(aov(MSE ~ Design * Bandwidth * Cov * Eta, Qua_df))


library(tidyve)

summary(aov(MSE ~ Design * Bandwidth * Cov * Eta, subset(Lin_df,Beta==1)))
summary(aov(MSE ~ Design * Bandwidth * Cov * Eta, subset(Qua_df,Beta==1)))
summary(aov(MSE ~ Design * Bandwidth * Cov * Eta, subset(Qua_df,Beta==2)))


summary(aov(MSE ~ Design + Bandwidth + Cov + Eta +
              Design:Bandwidth + Design:Cov + Design:Eta +
              Bandwidth:Cov + Bandwidth:Eta + Cov:Eta, 
            subset(Lin_df,Beta==1)))

summary(aov(MSE ~ Design + Bandwidth + Cov + Eta +
              Design:Bandwidth + Design:Cov + Design:Eta +
              Bandwidth:Cov + Bandwidth:Eta + Cov:Eta, 
            subset(Qua_df,Beta==1)))

summary(aov(MSE ~ Design + Bandwidth + Cov + Eta +
              Design:Bandwidth + Design:Cov + Design:Eta +
              Bandwidth:Cov + Bandwidth:Eta + Cov:Eta, 
            subset(Qua_df,Beta==2)))


library(effectsize)
# options(es.use_symbols = TRUE) # get nice symbols when printing! (On Windows, requires R >= 4.2.0)

cohens_d(MSE ~ Design + Bandwidth + Cov + Eta +
           Design:Bandwidth + Design:Cov + Design:Eta +
           Bandwidth:Cov + Bandwidth:Eta + Cov:Eta, 
         data=subset(Lin_df,Beta==1))


# 
# 
# ########## out put tables ############
# write.table(Lin_df,'Outcome/ExptNewPara/Simu1K2024/Lin_df.txt',quote = FALSE,row.names = FALSE)
# write.table(Qua_df,'Outcome/ExptNewPara/Simu1K2024/Qua_df.txt',quote = FALSE,row.names = FALSE)
# 
# 


library(GWmodel)
library(rethinking)
library(MASS)
library(ggplot2)
library(rethinking)
library(SpatialTools)

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

eta <- 1

set.seed(2024)
# QuadraDat_ar1 <- SimuGenDat_Qua(paralist,"AR1",eta)
QuadraDat_mat <- SimuGenDat_Qua(paralist,"Matern",eta)

QuadraDat

# QuadraDat <- parallel_Mat_coef_results[[1]][[1]]

QuadraDat$x <- QuadraDat$Row
QuadraDat$y <- QuadraDat$Range

sp_ar1_QuadraDat <- SpatialPointsDataFrame(cbind(QuadraDat$x,QuadraDat$y),QuadraDat)

gwr_ar1_qua.syst9 <- gwr.basic(SYield ~ SystNitro + SystNitro2,
                               data = sp_ar1_QuadraDat, bw=9,
                               kernel = "gaussian")  ## boxcar kernel

# 
# a0<- QuadraDat$b0+QuadraDat$u0
# b0<- QuadraDat$b1+QuadraDat$u1
# c0<- QuadraDat$b1+QuadraDat$u2
# 
# N0 <- pmax(0,pmin(140,-b0/(2*c0)))

hist(QuadraDat$b0+QuadraDat$u0)
hist(QuadraDat$b1+QuadraDat$u1)
hist(QuadraDat$b2+QuadraDat$u2)


a1<- gwr_ar1_qua.syst9$SDF$Intercept
b1<- gwr_ar1_qua.syst9$SDF$SystNitro
c1<- gwr_ar1_qua.syst9$SDF$SystNitro2

N1 <- pmax(0,pmin(140,-b1/(2*c1)))




gwr_ar1_qua.rand9 <- gwr.basic(RYield ~ RandNitro + RandNitro2,
                               data = sp_ar1_QuadraDat, bw = 9, 
                               kernel = "gaussian")  ## boxcar kernel

a2 <- gwr_ar1_qua.rand9$SDF$Intercept
b2 <- gwr_ar1_qua.rand9$SDF$RandNitro
c2 <- gwr_ar1_qua.rand9$SDF$RandNitro2

N2 <- pmax(0,pmin(140,-b2/(2*c2)))

QuadraDat$N1 <- N1
QuadraDat$N2 <- N2

ggplot(QuadraDat, aes(x=Row, y=Range, fill=N1)) +
  geom_tile() +
  labs(y="Column",fill="Optimal Nitrogen rate (kg/ha)") + 
  scale_fill_gradientn(colors=rev(terrain.colors(100)),limit=c(0,150)) +
  scale_x_continuous(breaks=c(1,93)) + 
  scale_y_continuous(breaks=c(1,20)) + 
  thm1

ggplot(QuadraDat, aes(x=Row, y=Range, fill=N2)) +
  geom_tile() +
  labs(y="Column",fill="Optimal Nitrogen rate (kg/ha)") + 
  scale_fill_gradientn(colors=rev(terrain.colors(100)),limit=c(0,150)) +
  scale_x_continuous(breaks=c(1,93)) + 
  scale_y_continuous(breaks=c(1,20)) + 
  thm1



ggplot(QuadraDat, aes(x=Row, y=Range, fill=SYield)) +
  geom_tile() +
  labs(y="Column",fill="Syst yield") + 
  scale_fill_gradientn(colors=rev(terrain.colors(100))) +
  scale_x_continuous(breaks=c(1,93)) + 
  scale_y_continuous(breaks=c(1,20)) + 
  thm1

ggplot(QuadraDat, aes(x=Row, y=Range, fill=RYield)) +
  geom_tile() +
  labs(y="Column",fill="Rand yield") + 
  scale_fill_gradientn(colors=rev(terrain.colors(100))) +
  scale_x_continuous(breaks=c(1,93)) + 
  scale_y_continuous(breaks=c(1,20)) + 
  thm1



dat.lab <- aggregate(cbind(SystNitro,RandNitro) ~ Range,
                     QuadraDat,mean)
dat.lab$Row <- 50

shade <- data.frame(x1=c(0,  0,   0,   0),
                    y1=c(0,  5.5, 10.5,15.5),
                    # x2=c(102,102,102,102),
                    x2=c(94,94,94,94),
                    y2=c(5.5,10.5,15.5,21))



ggplot(QuadraDat) + 
  geom_tile(aes(Row,Range,fill=as.factor(SystNitro)),
            size=0.25, alpha=0.7) + 
  #scale_fill_manual(values = c("gray10","gray30","gray50",
  #                             "gray70", "gray90")) +
  # geom_rect(data = shade, aes(xmin = x1, ymin = y1, xmax = x2, ymax = y2),
  #           color = "blue", alpha = 0, size =1) +
  geom_label(data= dat.lab, aes(Row,Range,label=paste(SystNitro,"kg/ha")),
             size= 6,fontface = "bold") + 
  annotate("text", x = 98, y = c(3,8,13,18), 
           label = paste("Rep", 1:4),size=7) + 
  ylab("Column") +
  scale_x_continuous(breaks=c(1,93)) + 
  scale_y_continuous(breaks=c(1,20)) + 
  scale_fill_brewer(palette="Paired") + 
  guides(fill=guide_legend(title="Nitrogen rate (kg/ha)")) + thm1 +
  theme(legend.position = "na")

# mysave("Col_SystNitro_V2")


ggplot(QuadraDat) + 
  geom_tile(aes(Row,Range,fill=as.factor(RandNitro)),
            size=0.25,alpha=0.7) + 
  #scale_fill_manual(values = c("gray10","gray30","gray50",
  #                             "gray70", "gray90")) +
  geom_rect(data = shade, aes(xmin = x1, ymin = y1, xmax = x2, ymax = y2),
            color = "blue", alpha = 0, size =1) +
  geom_label(data= dat.lab, aes(Row,Range,label=paste(RandNitro,"kg/ha")),
             size= 6,fontface = "bold") + 
  annotate("text", x = 98, y = c(3,8,13,18), 
           label = paste("Rep", 1:4),size=7) + 
  ylab("Column") +
  scale_x_continuous(breaks=c(1,93)) + 
  scale_y_continuous(breaks=c(1,20)) + 
  scale_fill_brewer(palette="Paired") + 
  guides(fill=guide_legend(title="Nitrogen rate (kg/ha)")) + thm1 +
  theme(legend.position = "na")








library(sp)
library(rgdal)
library(raster)
library(maptools)
library(spatstat)
library(tidyverse)

sp::coordinates(QuadraDat) <- c("Row","Range")
spplot(QuadraDat, 'RYield')

QuadraDat_ra_yield<- as(QuadraDat["RYield"], "ppp")
sm_QuadraDat_ra_yield <- Smooth(QuadraDat_ra_yield)

plot(sm_QuadraDat_ra_yield, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="")
contour(sm_QuadraDat_ra_yield,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)


QuadraDat_st_yield<- as(QuadraDat["SYield"], "ppp")
sm_QuadraDat_st_yield <- Smooth(QuadraDat_st_yield)

plot(sm_QuadraDat_st_yield, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="")
contour(sm_QuadraDat_st_yield,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)



QuadraDat$Diff <- QuadraDat$SYield-QuadraDat$RYield

QuadraDat_diff_yield<- as(QuadraDat["Diff"], "ppp")
sm_QuadraDat_diff_yield <- Smooth(QuadraDat_diff_yield)

plot(sm_QuadraDat_diff_yield,  ribsep=0.01, 
     ribside="bottom", main="")
contour(sm_QuadraDat_diff_yield,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)




a0 = QuadraDat$b0+QuadraDat$u0
b0 = QuadraDat$b1+QuadraDat$u1
c0 = QuadraDat$b2+QuadraDat$u2

QuadraDat$N0 <- pmax(0,pmin(140,-b0/(2*c0)))

QuadraDat_true_N<- as(QuadraDat["N0"], "ppp")
sm_QuadraDat_true_N <- Smooth(QuadraDat_true_N)

plot(sm_QuadraDat_true_N, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="", zlim=c(50, 140))
contour(sm_QuadraDat_true_N,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)


QuadraDat_st_N<- as(QuadraDat["N1"], "ppp")
sm_QuadraDat_st_N <- Smooth(QuadraDat_st_N)

plot(sm_QuadraDat_st_N, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="", zlim=c(50, 140))
contour(sm_QuadraDat_st_N,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)


QuadraDat_ra_N<- as(QuadraDat["N2"], "ppp")
sm_QuadraDat_ra_N <- Smooth(QuadraDat_ra_N)

plot(sm_QuadraDat_ra_N, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="", zlim=c(50, 140))
contour(sm_QuadraDat_ra_N,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)



###### explore density ###########

library(scales)
library(ggridges)
library(ggplot2)

# parallel_AR1_coef_results <- readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_AR1_coef_results_mac.rds")
# parallel_Mat_coef_results <- readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_Mat_coef_results_mac.rds")

parallel_AR1_coef_results <- readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_AR1_coef_results_1K_mac.rds")
parallel_Mat_coef_results <- readRDS("Outcome/ExptNewPara/Simu1K2024/parallel_Mat_coef_results_1K_mac.rds")


length(parallel_AR1_coef_results)
dim(parallel_AR1_coef_results[[1]][[1]])

ar1_b0 <- NULL
ar1_b1 <- NULL
ar1_b2 <- NULL

ar1_SystGWRB0 <- NULL
ar1_SystGWRB1 <- NULL
ar1_SystGWRB2 <- NULL

ar1_RandGWRB0 <- NULL
ar1_RandGWRB1 <- NULL
ar1_RandGWRB2 <- NULL

for(i in 1:10){
  for(j in 1:100){
    temp <- parallel_AR1_coef_results[[i]][[j]]
    
    ar1_b0 <- c(ar1_b0, temp$b0+temp$u0)
    ar1_b1 <- c(ar1_b1, temp$b1+temp$u1)
    ar1_b2 <- c(ar1_b2, temp$b2+temp$u2)
    
    ar1_SystGWRB0 <- c(ar1_SystGWRB0,temp$SystGWRB0)
    ar1_SystGWRB1 <- c(ar1_SystGWRB1,temp$SystGWRB1)
    ar1_SystGWRB2 <- c(ar1_SystGWRB2,temp$SystGWRB2)
    
    ar1_RandGWRB0 <- c(ar1_RandGWRB0,temp$RandGWRB0)
    ar1_RandGWRB1 <- c(ar1_RandGWRB1,temp$RandGWRB1)
    ar1_RandGWRB2 <- c(ar1_RandGWRB2,temp$RandGWRB2)
  }
  cat(i)
}


ar1_df_b0 <- data.frame(value = c(ar1_b0,ar1_SystGWRB0,ar1_RandGWRB0),
                        type = rep(c("True","Systematic","Randomised"),each=length(ar1_b0)))

ggplot(ar1_df_b0) +
  geom_density(aes(value, y = ..density.., fill = type), alpha = 0.3) +
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[0]),
       x = "Value",
       y = "Density") + 
  # scale_y_continuous(breaks = c(0.1, 0.2, 0.3,0.4,0.5),
  #                    labels = c(0.1, 0.2, 0.3,0.4,0.5)) + 
  thm1 

  
ar1_df_b1 <- data.frame(value = c(ar1_b1,ar1_SystGWRB1,ar1_RandGWRB1),
                        type = rep(c("True","Systematic","Randomised"),each=length(ar1_b1)))
ggplot(ar1_df_b1) +
  geom_density(aes(value,y=..density..,fill=type),alpha=0.3) +
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[1]),
       x = "Value",
       y = "Density") + 
  thm1 


ar1_df_b2 <- data.frame(value = c(ar1_b2,ar1_SystGWRB2,ar1_RandGWRB2),
                        type = rep(c("True","Systematic","Randomised"),each=length(ar1_b2)))
ggplot(ar1_df_b2) +
  geom_density(aes(value,y=..density..,fill=type),alpha=0.3) +
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[1]),
       x = "Value",
       y = "Density") + 
  thm1 


mat_b0 <- NULL
mat_b1 <- NULL
mat_b2 <- NULL

mat_SystGWRB0 <- NULL
mat_SystGWRB1 <- NULL
mat_SystGWRB2 <- NULL

mat_RandGWRB0 <- NULL
mat_RandGWRB1 <- NULL
mat_RandGWRB2 <- NULL

for(i in 1:10){
  for(j in 1:100){
    temp <- parallel_Mat_coef_results[[i]][[j]]
    
    mat_b0 <- c(mat_b0,temp$b0+temp$u0)
    mat_b1 <- c(mat_b1,temp$b1+temp$u1)
    mat_b2 <- c(mat_b2,temp$b2+temp$u2)
    
    mat_SystGWRB0 <- c(mat_SystGWRB0,temp$SystGWRB0)
    mat_SystGWRB1 <- c(mat_SystGWRB1,temp$SystGWRB1)
    mat_SystGWRB2 <- c(mat_SystGWRB2,temp$SystGWRB2)
    
    mat_RandGWRB0 <- c(mat_RandGWRB0,temp$RandGWRB0)
    mat_RandGWRB1 <- c(mat_RandGWRB1,temp$RandGWRB1)
    mat_RandGWRB2 <- c(mat_RandGWRB2,temp$RandGWRB2)
  }
  cat(i)
}


mat_df_b0 <- data.frame(value = c(mat_b0,mat_SystGWRB0,mat_RandGWRB0),
                        type = rep(c("True","Systematic","Randomised"),each=length(mat_b0)))

ggplot(mat_df_b0) +
  geom_density(aes(value,y=..density..,fill=type),alpha=0.3) +
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[0]),
       x = "Value",
       y = "Density") + 
  thm1 


mat_df_b1 <- data.frame(value = c(mat_b1,mat_SystGWRB1,mat_RandGWRB1),
                        type = rep(c("True","Systematic","Randomised"),each=length(mat_b1)))
ggplot(mat_df_b1) +
  geom_density(aes(value,y=..density..,fill=type),alpha=0.3)+
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[1]),
       x = "Value",
       y = "Density") + 
  thm1 


mat_df_b2 <- data.frame(value = c(mat_b2,mat_SystGWRB2,mat_RandGWRB2),
                        type = rep(c("True","Systematic","Randomised"),each=length(mat_b2)))
ggplot(mat_df_b2) +
  geom_density(aes(value,y=..density..,fill=type),alpha=0.3) +
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[2]),
       x = "Value",
       y = "Density") + 
  thm1 




mat_b0 <- NULL
mat_b1 <- NULL
mat_b2 <- NULL

mat_SystGWRB0_bw5 <- NULL
mat_SystGWRB1_bw5 <- NULL
mat_SystGWRB2_bw5 <- NULL

mat_RandGWRB0_bw5 <- NULL
mat_RandGWRB1_bw5 <- NULL
mat_RandGWRB2_bw5 <- NULL

for(i in 1:10){
  for(j in 1:100){
    temp <- parallel_Mat_coef_results[[i]][[j]]
    
    mat_b0 <- c(mat_b0,temp$b0+temp$u0)
    mat_b1 <- c(mat_b1,temp$b1+temp$u1)
    mat_b2 <- c(mat_b2,temp$b2+temp$u2)
    
    mat_SystGWRB0_bw5 <- c(mat_SystGWRB0_bw5,temp$SystGWRB0_bw5)
    mat_SystGWRB1_bw5 <- c(mat_SystGWRB1_bw5,temp$SystGWRB1_bw5)
    mat_SystGWRB2_bw5 <- c(mat_SystGWRB2_bw5,temp$SystGWRB2_bw5)
    
    mat_RandGWRB0_bw5 <- c(mat_RandGWRB0_bw5,temp$RandGWRB0_bw5)
    mat_RandGWRB1_bw5 <- c(mat_RandGWRB1_bw5,temp$RandGWRB1_bw5)
    mat_RandGWRB2_bw5 <- c(mat_RandGWRB2_bw5,temp$RandGWRB2_bw5)
  }
  cat(i)
}


mat_df_b0_bw5 <- data.frame(value = c(mat_b0,mat_SystGWRB0_bw5,mat_RandGWRB0_bw5),
                        type = rep(c("True","Systematic","Randomised"),each=length(mat_b0)))

ggplot(mat_df_b0_bw5) +
  geom_density(aes(value,y=..density..,fill=type),alpha=0.3) +
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[0]),
       x = "Value",
       y = "Density") + 
  thm1 


mat_df_b1_bw5 <- data.frame(value = c(mat_b1,mat_SystGWRB1_bw5,mat_RandGWRB1_bw5),
                        type = rep(c("True","Systematic","Randomised"),each=length(mat_b1)))
ggplot(mat_df_b1_bw5) +
  geom_density(aes(value,y=..density..,fill=type),alpha=0.3)+
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[1]),
       x = "Value",
       y = "Density") + 
  thm1 


mat_df_b2_bw5 <- data.frame(value = c(mat_b2,mat_SystGWRB2_bw5,mat_RandGWRB2_bw5),
                        type = rep(c("True","Systematic","Randomised"),each=length(mat_b2)))
ggplot(mat_df_b2_bw5) +
  geom_density(aes(value,y=..density..,fill=type),alpha=0.3) +
  scale_fill_manual(values = c("True" = "white", "Systematic" = "gold", "Randomised" = "navy")) +
  labs(fill = expression(beta[2]),
       x = "Value",
       y = "Density") + 
  thm1 





plot(mat_b2,mat_SystGWRB2_bw5)
abline(0,1)

plot(mat_b2,mat_RandGWRB2_bw5)
abline(0,1)












