# Import libraries
{
  library(collapse)
  library(magrittr)
  library(tseries)
  library(lmtest)
  library(sandwich)
  library(jtools)
  library(xts)
  library(ggfortify)
  library(lubridate)
  library(MTE)
  library(quantreg)
  source("lsa.r")
  library(meboot)
  library(foreach)
  library(doParallel)
  library(future)
  library(furrr)
  library(purrr)
  library(rqPen)
}

# Initialize the asynchronous notification on my phone 
initialize_info <- function() {
  if (!require(pacman)) install.packages('pacman'); library(pacman)
  pacman::p_load(devtools)
  if (!require(RPushbullet)) devtools::install_github("eddelbuettel/rpushbullet"); library(RPushbullet)
  library(RPushbullet)
  pbSetup(apikey = "o.WfW4TDpXDx9Szb6Wk0JvwRHNW462ULP3", defdev = 0)
}

# Data reading
# data.path <- file.choose("APRIL_MIN_3.csv")
data <- read.csv("APRIL_MIN_3.csv")

# Date Construction 
fixDate <- function(x, year){
  n_row = sum(data$DAY == data$DAY[1])
  date.df <- matrix(NA, ncol = length(unique(x)), nrow= n_row)
  index <- 1
  for (d in unique(data$DAY)){
    date.df[,index] <- seq(as.POSIXct(paste0(year,"-04-",d," 00:00:00"),tz="Europe/Athens"), by="3 min", length.out=n_row)
    index <- index+1
  }
  return(c(as.POSIXct((date.df))))
}

colnames(data) <- c("L101_volume",    "L102_volume",    "L103_volume",    "L104_volume",    "L106_volume",    "L107_volume",   
                    "L108_volume",    "L101_occupancy", "L102_occupancy", "L103_occupancy", "L104_occupancy", "L106_occupancy",
                    "L107_occupancy", "L108_occupancy", "DAY",            "TIME")
n_row = sum(data$DAY == data$DAY[1])*5
settfm(data, DATE = fixDate(data$DAY,year=2000))
# tsdata <- data %$% xts(cbind(L101_volume, L101_occupancy, L102_volume, L102_occupancy,
#                              L103_volume, L103_occupancy, L104_volume, L104_occupancy,
#                              L106_volume, L106_occupancy, L107_volume, L107_occupancy,
#                              L108_volume, L108_occupancy),order.by = DATE, Frequency = n_row)
# plot_list <- list()
# for (loc in Filter(function(x) grepl("L10._volume",x), colnames(data))){
#   df <- data[c(loc,"DATE")]
#   names(df) <- c("Y1","DATE")
#   plot_list[[loc]] <- ggplot(df, aes(x = DATE, y = Y1)) +
#                           geom_line(color = "navyblue") +
#                           labs(title = paste0("Plot of Traffic Volume of ",strsplit(loc,"_")[[1]][1]), x = "Date", y = "Volume") +
#                           theme_minimal()
# }
# library(gridExtra)
# num_plots <- length(plot_list)
# num_cols <- 2  
# num_rows <- ceiling(num_plots / num_cols)
# dev.new()
# grid.arrange(grobs = plot_list, nrow = num_rows, ncol = num_cols)

# Construct the regressors respecting Fourier Harmonics. We need the first 500. 
# High harmonic frequencies seems to have low impact.  
regressors <- c()
threshold <- 500 #n_row/2 - 1
for (i in 1:(threshold)){
  regressors <- c(regressors, paste("sin(2 * pi *",as.character(i),"*TIME/",as.character(n_row),")"))
}
for (i in 1:(threshold)){
  regressors <- c(regressors, paste("cos(2 * pi *",as.character(i),"*TIME/",as.character(n_row),")"))
}
f <- paste("~",paste(regressors,collapse="+"))

# AIC and BIC respecting L1 and L2 norm in correspondence.
robustIC <- function(X,x,y,measure="AIC"){length(y)*log(mad(y-X%*%x))+sum(x!=0)*ifelse(measure=="AIC",2,log(length(y)))}
MyIC <- function(X,x,y,measure="AIC"){length(y)*log(mean((y-X%*%x)^2))+sum(x!=0)*ifelse(measure=="AIC",2,log(length(y)))}


# Collect all routes
routes <- Filter(function(x) grepl("_volume",x), colnames(data))
train_rows <- 1:(n_row*2)

# Create MEBoot samples and store them locally.
# setwd(choose.dir())
genData <- function(Niter=20){
  for (selRoute in routes){
    cat("FOR",selRoute,"\n")
    setwd(strsplit(selRoute,"_")[[1]][1])
    MEB = meboot(data[[selRoute]],reps=Niter)$ensemble
    save(MEB,file= "synth_dataset.RData")
  }
}
# genData(20)

doProcess <- function(){
  # Set up the communication with my phone and all the proper configurations for the L1/L2 LSA/MTE procedures on
  # MEBoot samples, extract coefficients using AIC/BIC criterion and store them locally.
  initialize_info()
  coeff_names <- c("(Intercept)", regressors)
  num_coeffs <- length(coeff_names)
  na_result_vec <- numeric(num_coeffs)
  names(na_result_vec) <- coeff_names
  n_cores <- 15
  Niter <- 20
  
  # Choose which procedure would you like to apply

  input = readline("Order of parallel procedures: l2/l1 LSA l2/l1 MTE.\n Select with T/F which procedures you want to perform.")
  procedures <- sapply(strsplit(toupper(input),"")[[1]], as.logical)

  for (selRoute in routes){
    setwd(strsplit(selRoute,"_")[[1]][1])
    env = new.env()
    load("synth_dataset.RData",envir=env)
    MEB_DATA <- t(env$MEB)
    rm(env)
    if (procedures[1]){
      beepr::beep_on_error({
        cl <- makeCluster(n_cores)
        registerDoParallel(cl)
        clusterExport(cl, c("MEB_DATA", "data", "f", "train_rows", "na_result_vec", "num_coeffs", "coeff_names", "Niter")) # Add others if needed
        cat(sprintf("Registered parallel backend with %d cores.\n", n_cores))
        
        coef_LS_LSA <- foreach(
          row = 1:Niter,
          .combine = 'list',
          .multicombine = TRUE,
          .errorhandling = 'pass',
          .packages = c("quantreg", "lars")
        ) %dopar% {
          cat(sprintf("Number of Iteration:", row,"\n"))
          iter_results <- data.frame(
            lm_aic  = na_result_vec,
            lm_bic  = na_result_vec
          )
          current_volume <- MEB_DATA[row, ]
          dt <- data.frame("volume" = current_volume, "TIME" = data$TIME)[train_rows, ]
          
          PROGRESS <- paste0("STARTING POINT FOR MEBOOT SAMPLE: ",row,"   ", Sys.time(),"\n"); write.csv(PROGRESS,file= "progress.txt")
          lm.reg  <- lm(paste0("volume", f), data = dt)
          X <- as.matrix(cbind("(Intercept)"=1,lm.reg$model[,-1]))
          y <- lm.reg$model[,1]
          N <- nrow(X)
          P <- ncol(X)
          Sigma0.lm  <- vcov(lm.reg)
          b0_lm  <- coef(lm.reg)
          
          PROGRESS <- paste0(PROGRESS, "STARTING LSA LM: ", Sys.time(),"\n"); write.csv(PROGRESS,file= "progress.txt")
          LSA.lm <- lars.lsa(Sigma0 = Sigma0.lm, b0 = b0_lm, intercept = FALSE, n = N, type = "lasso")
          PROGRESS <- paste0(PROGRESS, "FINISH LSA ESTIMATION:",Sys.time()," \n"); write.csv(PROGRESS,file= "progress.txt")
          
          iter_results$lm_aic <- LSA.lm$beta[which(LSA.lm$AIC==min(LSA.lm$AIC))[1], ]
          iter_results$lm_bic <- LSA.lm$beta[which(LSA.lm$BIC==min(LSA.lm$BIC))[1], ]
          
          do.call(rbind, iter_results)
          
          return(iter_results)
        }
        L2LSA_AIC <- parallel::mclapply(coef_LS_LSA,FUN = function(x) x$lm_aic)
        L2LSA_BIC <- parallel::mclapply(coef_LS_LSA,FUN = function(x) x$lm_bic)
        coef_L2LSA <- list(
          AIC = matrix(unlist(L2LSA_AIC), nrow=num_coeffs, ncol = Niter), 
          BIC = matrix(unlist(L2LSA_BIC), nrow=num_coeffs, ncol = Niter)
        )
        rownames(coef_L2LSA$AIC) <- rownames(coef_L2LSA$BIC) <- coeff_names
        save(coef_L2LSA,file = "coef_L2LSA.RData")
        stopCluster(cl)
        registerDoSEQ()
      },7)
    }
    if (procedures[2]){
      beepr::beep_on_error({
        cl <- makeCluster(n_cores)
        registerDoParallel(cl)
        clusterExport(cl, c("MEB_DATA", "data", "f", "train_rows",  "na_result_vec", "num_coeffs", "coeff_names", "Niter")) # Add others if needed
        cat(sprintf("Registered parallel backend with %d cores.\n", n_cores))
        
        coef_LAD_LSA <- foreach(
          row = 1:Niter,
          .combine = 'list',
          .multicombine = TRUE,
          .errorhandling = 'pass',
          .packages = c("quantreg", "lars")
        ) %dopar% {
          cat(sprintf("Number of Iteration:", row,"\n"))
          iter_results <- data.frame(
            lad_aic  = na_result_vec,
            lad_bic  = na_result_vec
          )
          current_volume <- MEB_DATA[row, ]
          dt <- data.frame("volume" = current_volume, "TIME" = data$TIME)[train_rows, ]
          
          PROGRESS <- paste0("STARTING POINT FOR MEBOOT SAMPLE: ",row,"   ", Sys.time(),"\n"); write.csv(PROGRESS,file= "progress.txt")
          lad.reg  <- rq(paste0("volume", f), data = dt,method = "lasso", lambda = 0.001)
          
          X <- as.matrix(cbind("(Intercept)"=1,lad.reg$model[,-1]))
          y <- lad.reg$model[,1]
          N <- nrow(X)
          P <- ncol(X)
          Sigma0.lad  <- solve(t(X)%*%X)/density(lad.reg$residuals,n=1,from=0,to=0)$y^2 #vcov.rq(lad.reg)
          b0_lad  <- coef(lad.reg)
          
          PROGRESS <- paste0(PROGRESS, "STARTING LSA LAD: ", Sys.time(),"\n"); write.csv(PROGRESS,file= "progress.txt")
          LSA.lad <- lars.lsa(Sigma0 = Sigma0.lad, b0 = b0_lad, intercept = FALSE, n = N, type = "lasso")
          PROGRESS <- paste0(PROGRESS, "FINISH LSA ESTIMATION:",Sys.time()," \n"); write.csv(PROGRESS,file= "progress.txt")
          
          iter_results$lad_aic <- LSA.lad$beta[which(LSA.lad$AIC==min(LSA.lad$AIC))[1], ]
          iter_results$lad_bic <- LSA.lad$beta[which(LSA.lad$BIC==min(LSA.lad$BIC))[1], ]
          
          do.call(rbind, iter_results)
          
          return(iter_results)
        }
        L1LSA_AIC <- parallel::mclapply(coef_LAD_LSA,FUN = function(x) x$lad_aic)
        L1LSA_BIC <- parallel::mclapply(coef_LAD_LSA,FUN = function(x) x$lad_bic)
        coef_L1LSA <- list(
          AIC = matrix(unlist(L1LSA_AIC), nrow=num_coeffs, ncol = Niter), 
          BIC = matrix(unlist(L1LSA_BIC), nrow=num_coeffs, ncol = Niter)
        )
        rownames(coef_L1LSA$AIC) <- rownames(coef_L1LSA$BIC) <- coeff_names
        save(coef_L1LSA,file = "coef_L1LSA.RData")
        stopCluster(cl)
        registerDoSEQ()
      },7 )
    }
    if (procedures[3]){
      beepr::beep_on_error({
        cl <- makeCluster(n_cores)
        registerDoParallel(cl)
        clusterExport(cl, c("MEB_DATA", "data", "f", "train_rows",  "na_result_vec","MyIC", "num_coeffs", "coeff_names", "Niter")) # Add others if needed
        cat(sprintf("Registered parallel backend with %d cores.\n", n_cores))
        
        coef_LS_MTE <- foreach(
          row = 1:Niter,                  # Iterate over each parameter combination
          .combine = 'list',              # Collect results as a list
          .multicombine = TRUE,           # Efficient list combining
          .packages = c("RPushbullet"),      # Load necessary packages on workers
          .errorhandling = 'pass'         # Continue if one task fails
        ) %dopar% {
          
          lambda_sequence <- seq(0.02, 0.1, 0.02)
          n_lambda <- length(lambda_sequence)
          iter_results <- data.frame(
            lm_aic  = na_result_vec,
            lm_bic  = na_result_vec
          )
          dt <- data.frame(volume = MEB_DATA[row,], TIME= data$TIME)[train_rows, ]
          
          # PROGRESS <- paste0("Starting LS estimation: ", Sys.time()," \n"); pbPost("note", title="LS MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
          lm.reg  <- lm(paste0("volume", f), data = dt)
          X <- as.matrix(cbind("(Intercept)"=1,lm.reg$model[,-1]))
          y <- lm.reg$model[,1]
          
          # PROGRESS <- paste0(PROGRESS, "STARTING MTE PROCESS:",Sys.time()," \n"); pbPost("note", title="LS MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
          MTE_COEFFS <- matrix(NA, ncol=n_lambda,nrow=num_coeffs)
          for (lambda_iter in 1:length(lambda_sequence)){
            mte_result <- MTE::MTElasso(X, y, lm.reg$coefficients, intercept = FALSE, lambda = lambda_sequence[lambda_iter])
            MTE_COEFFS[, lambda_iter] <- mte_result$beta
          }
          # PROGRESS <- paste0(PROGRESS, "FINISH MTE ESTIMATION:",Sys.time()," \n"); pbPost("note", title="LS MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
          
          AICs <- apply(MTE_COEFFS,2, FUN = function(x) MyIC(X,x,y,"AIC"))
          BICs <- apply(MTE_COEFFS,2, FUN = function(x) MyIC(X,x,y,"BIC"))
          iter_results$lm_aic <- MTE_COEFFS[,which(AICs==min(AICs))[1]] 
          iter_results$lm_bic <- MTE_COEFFS[,which(BICs==min(BICs))[1]] 
          names(iter_results$lm_aic) <- names(iter_results$lm_bic) <- coeff_names
          return(iter_results)
        } 
        L2MTE_AIC <- parallel::mclapply(coef_LS_MTE,FUN = function(x) x$lm_aic)
        L2MTE_BIC <- parallel::mclapply(coef_LS_MTE,FUN = function(x) x$lm_bic)
        coef_L2MTE <- list(
          AIC = matrix(unlist(L2MTE_AIC), nrow=num_coeffs, ncol = Niter), 
          BIC = matrix(unlist(L2MTE_BIC), nrow=num_coeffs, ncol = Niter)
        )
        rownames(coef_L2MTE$AIC) <- rownames(coef_L2MTE$BIC) <- coeff_names
        save(coef_L2MTE,file = "coef_L2MTE.RData")
        stopCluster(cl)
        registerDoSEQ()
      },7)
    }
    if (procedures[4]){
      beepr::beep_on_error({
        cl <- makeCluster(n_cores)
        registerDoParallel(cl)
        clusterExport(cl, c("MEB_DATA", "data", "f", "train_rows", "robustIC", "na_result_vec", "num_coeffs", "coeff_names", "Niter")) # Add others if needed
        cat(sprintf("Registered parallel backend with %d cores.\n", n_cores))
        
        coef_LAD_MTE <- foreach(
          row = 1:Niter,                  # Iterate over each parameter combination
          .combine = 'list',              # Collect results as a list
          .multicombine = TRUE,           # Efficient list combining
          .packages = c("quantreg", "RPushbullet"),      # Load necessary packages on workers
          .errorhandling = 'pass'         # Continue if one task fails
        ) %dopar% {
          
          lambda_sequence <- seq(0.02, 0.1, 0.02)
          n_lambda <- length(lambda_sequence)
          iter_results <- data.frame(
            lad_aic  = na_result_vec,
            lad_bic  = na_result_vec
          )
          dt <- data.frame(volume = MEB_DATA[row,], TIME= data$TIME)[train_rows, ]
          # PROGRESS <- paste0("Starting LAD estimation: ", Sys.time()," \n"); pbPost("note", title="LAD MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
          lad.reg <- rq(as.formula(paste0("volume", f)), data = dt)
          X <- lad.reg$x
          y <- dt$volume 
          # PROGRESS <- paste0(PROGRESS, "STARTING MTE PROCESS:",Sys.time()," \n"); pbPost("note", title="LAD MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
          MTE_COEFFS <- matrix(NA, ncol=n_lambda,nrow=num_coeffs)
          for (lambda_iter in 1:length(lambda_sequence)){
            mte_result <- MTE::MTElasso(X, y, lad.reg$coefficients, intercept = FALSE, lambda = lambda_sequence[lambda_iter])
            MTE_COEFFS[, lambda_iter] <- mte_result$beta
          }
          # PROGRESS <- paste0(PROGRESS, "FINISH MTE ESTIMATION:",Sys.time()," \n"); pbPost("note", title="LAD MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
          
          AICs <- apply(MTE_COEFFS,2, FUN = function(x) robustIC(X,x,y,"AIC"))
          BICs <- apply(MTE_COEFFS,2, FUN = function(x) robustIC(X,x,y,"BIC"))
          iter_results$lad_aic <- MTE_COEFFS[,which(AICs==min(AICs))[1]] 
          iter_results$lad_bic <- MTE_COEFFS[,which(BICs==min(BICs))[1]] 
          names(iter_results$lad_aic) <- names(iter_results$lad_bic) <- coeff_names 
          return(iter_results)
        } 
        L1MTE_AIC <- parallel::mclapply(coef_LAD_MTE,FUN = function(x) x$lad_aic)
        L1MTE_BIC <- parallel::mclapply(coef_LAD_MTE,FUN = function(x) x$lad_bic)
        coef_L1MTE <- list(
          AIC = matrix(unlist(L1MTE_AIC), nrow=num_coeffs, ncol = Niter), 
          BIC = matrix(unlist(L1MTE_BIC), nrow=num_coeffs, ncol = Niter)
        )
        rownames(coef_L1MTE$AIC) <- rownames(coef_L1MTE$BIC) <- coeff_names
        save(coef_L1MTE,file = "coef_L1MTE.RData") 
        stopCluster(cl)
        registerDoSEQ()
      },7)
    }
    setwd("../")
  }

  # Choose which procedures and on which locations  would you like to apply
  inputMethod = readline("Which coefficients would you like to check? \n Order of  procedures: l2/l1 LSA l2/l1 MTE.\n Select with T/F which methods you want to check")
  procedures <- sapply(strsplit(toupper(inputMethod),"")[[1]], as.logical)
  inputRoute = readline("Which locations would you like to check? \n Order of  procedures: L101, L102, L103, L104, L106, L107, L108.\n Select with T/F which locations you want to check.")
  routesToCheck <- sapply(strsplit(toupper(inputRoute),"")[[1]], as.logical)
  n_cores = 7
  
  # L1/L2 LSA/MTE Coefficients on initial data. 
  if (procedures[1]){
    beepr::beep_on_error({
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      clusterExport(cl, c("data", "f", "train_rows", "na_result_vec", "num_coeffs", "coeff_names")) # Add others if needed
      cat(sprintf("Registered parallel backend with %d cores.\n", n_cores))
      
      coef_LS_LSA <- foreach(
        selRoute = routes[routesToCheck],
        .combine = 'list',
        .multicombine = TRUE,
        .errorhandling = 'pass',
        .packages = c("lars")
      ) %dopar% {
        cat(sprintf("Number of Iteration:", row,"\n"))
        iter_results <- data.frame(
          lm_aic  = na_result_vec,
          lm_bic  = na_result_vec
        )
        dt <- data.frame("volume" = data[[selRoute]], "TIME" = data$TIME)[train_rows, ]
        
        PROGRESS <- paste0("STARTING POINT FOR ORIGINAL SAMPLE OF LOCATION: ",selRoute,"   ", Sys.time(),"\n"); write.csv(PROGRESS,file= "progress.txt")
        lm.reg  <- lm(paste0("volume", f), data = dt)
        X <- as.matrix(cbind("(Intercept)"=1,lm.reg$model[,-1]))
        y <- lm.reg$model[,1]
        N <- nrow(X)
        P <- ncol(X)
        Sigma0.lm  <- vcov(lm.reg)
        b0_lm  <- coef(lm.reg)
        
        PROGRESS <- paste0(PROGRESS, "STARTING LSA LM: ", Sys.time(),"\n"); write.csv(PROGRESS,file= "progress.txt")
        LSA.lm <- lars.lsa(Sigma0 = Sigma0.lm, b0 = b0_lm, intercept = FALSE, n = N, type = "lasso")
        PROGRESS <- paste0(PROGRESS, "FINISH LSA ESTIMATION:",Sys.time()," \n"); write.csv(PROGRESS,file= "progress.txt")
        
        iter_results$lm_aic <- LSA.lm$beta[which(LSA.lm$AIC==min(LSA.lm$AIC))[1], ]
        iter_results$lm_bic <- LSA.lm$beta[which(LSA.lm$BIC==min(LSA.lm$BIC))[1], ]
        
        do.call(rbind, iter_results)
        
        return(iter_results)
      }
      L2LSA_AIC <- parallel::mclapply(coef_LS_LSA,FUN = function(x) x$lm_aic)
      L2LSA_BIC <- parallel::mclapply(coef_LS_LSA,FUN = function(x) x$lm_bic)
      coef_L2LSA_Original <- list(
        AIC = matrix(unlist(L2LSA_AIC), nrow=num_coeffs, ncol = length(routes)), 
        BIC = matrix(unlist(L2LSA_BIC), nrow=num_coeffs, ncol = length(routes))
      )
      rownames(coef_L2LSA_Original$AIC) <- rownames(coef_L2LSA_Original$BIC) <- coeff_names
      colnames(coef_L2LSA_Original$AIC) <- colnames(coef_L2LSA_Original$BIC) <- sapply(routes,FUN = function(x) strsplit(x,"_")[[1]][1])
      
      save(coef_L2LSA_Original,file = "coef_L2LSA_Original.RData")
      stopCluster(cl)
      registerDoSEQ()
    },7)
  }
  if (procedures[2]){
    beepr::beep_on_error({
        cl <- makeCluster(n_cores)
        registerDoParallel(cl)
        clusterExport(cl, c("data", "f", "train_rows", "na_result_vec", "num_coeffs", "coeff_names"))
        cat(sprintf("Registered parallel backend with %d cores.\n", n_cores))
        
        coef_LAD_LSA <- foreach(
          selRoute = routes[routesToCheck],
          .combine = 'list',
          .multicombine = TRUE,
          .errorhandling = 'pass',
          .packages = c("quantreg", "lars")
        ) %dopar% {
        cat(sprintf("Number of Iteration:", row,"\n"))
        iter_results <- data.frame(
          lad_aic  = na_result_vec,
          lad_bic  = na_result_vec
        )
        dt <- data.frame("volume" = data[[selRoute]], "TIME" = data$TIME)[train_rows, ]
        
        PROGRESS <- paste0("STARTING POINT FOR ORIGINAL SAMPLE OF LOCATION: ",selRoute,"   ", Sys.time(),"\n"); write.csv(PROGRESS,file= "progress.txt")
        lad.reg  <- rq(paste0("volume", f), data = dt,method = "lasso", lambda = 0.001)
        
        X <- as.matrix(cbind("(Intercept)"=1,lad.reg$model[,-1]))
        y <- lad.reg$model[,1]
        N <- nrow(X)
        P <- ncol(X)
        Sigma0.lad  <- solve(t(X)%*%X)/density(lad.reg$residuals,n=1,from=0,to=0)$y^2 #vcov.rq(lad.reg)
        b0_lad  <- coef(lad.reg)
        
        PROGRESS <- paste0(PROGRESS, "STARTING LSA LAD: ", Sys.time(),"\n"); write.csv(PROGRESS,file= "progress.txt")
        LSA.lad <- lars.lsa(Sigma0 = Sigma0.lad, b0 = b0_lad, intercept = FALSE, n = N, type = "lasso")
        PROGRESS <- paste0(PROGRESS, "FINISH LSA ESTIMATION:",Sys.time()," \n"); write.csv(PROGRESS,file= "progress.txt")
        
        iter_results$lad_aic <- LSA.lad$beta[which(LSA.lad$AIC==min(LSA.lad$AIC))[1], ]
        iter_results$lad_bic <- LSA.lad$beta[which(LSA.lad$BIC==min(LSA.lad$BIC))[1], ]
        
        do.call(rbind, iter_results)
        
        return(iter_results)
      }
      L1LSA_AIC <- parallel::mclapply(coef_LAD_LSA,FUN = function(x) x$lad_aic)
      L1LSA_BIC <- parallel::mclapply(coef_LAD_LSA,FUN = function(x) x$lad_bic)
      coef_L1LSA_Original <- list(
        AIC = matrix(unlist(L1LSA_AIC), nrow=num_coeffs, ncol = length(routes)), 
        BIC = matrix(unlist(L1LSA_BIC), nrow=num_coeffs, ncol = length(routes))
      )
      rownames(coef_L1LSA_Original$AIC) <- rownames(coef_L1LSA_Original$BIC) <- coeff_names
      colnames(coef_L1LSA_Original$AIC) <- colnames(coef_L1LSA_Original$BIC) <- sapply(routes,FUN = function(x) strsplit(x,"_")[[1]][1])
      
      save(coef_L1LSA_Original,file = "coef_L1LSA_Original.RData")
      stopCluster(cl)
      registerDoSEQ()
    },7 )
  }
  if (procedures[3]){
    beepr::beep_on_error({
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      clusterExport(cl, c("MEB_DATA", "data", "f", "train_rows",  "na_result_vec","MyIC", "num_coeffs", "coeff_names", "Niter")) # Add others if needed
      cat(sprintf("Registered parallel backend with %d cores.\n", n_cores))
      
      coef_LS_MTE <- foreach(
        selRoute = routes[routesToCheck],
        .combine = 'list',
        .multicombine = TRUE,
        .packages = c("RPushbullet"),
        .errorhandling = 'pass'
      ) %dopar% {
        
        lambda_sequence <- seq(0.02, 0.1, 0.02)
        n_lambda <- length(lambda_sequence)
        iter_results <- data.frame(
          lm_aic  = na_result_vec,
          lm_bic  = na_result_vec
        )
        dt <- data.frame("volume" = data[[selRoute]], "TIME" = data$TIME)[train_rows, ]
        
        # PROGRESS <- paste0("Starting LS estimation: ", Sys.time()," \n"); pbPost("note", title="LS MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
        lm.reg  <- lm(paste0("volume", f), data = dt)
        X <- as.matrix(cbind("(Intercept)"=1,lm.reg$model[,-1]))
        y <- lm.reg$model[,1]
        
        # PROGRESS <- paste0(PROGRESS, "STARTING MTE PROCESS:",Sys.time()," \n"); pbPost("note", title="LS MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
        MTE_COEFFS <- matrix(NA, ncol=n_lambda,nrow=num_coeffs)
        for (lambda_iter in 1:length(lambda_sequence)){
          mte_result <- MTE::MTElasso(X, y, lm.reg$coefficients, intercept = FALSE, lambda = lambda_sequence[lambda_iter])
          MTE_COEFFS[, lambda_iter] <- mte_result$beta
        }
        # PROGRESS <- paste0(PROGRESS, "FINISH MTE ESTIMATION:",Sys.time()," \n"); pbPost("note", title="LS MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
        
        AICs <- apply(MTE_COEFFS,2, FUN = function(x) MyIC(X,x,y,"AIC"))
        BICs <- apply(MTE_COEFFS,2, FUN = function(x) MyIC(X,x,y,"BIC"))
        iter_results$lm_aic <- MTE_COEFFS[,which(AICs==min(AICs))[1]] 
        iter_results$lm_bic <- MTE_COEFFS[,which(BICs==min(BICs))[1]] 
        names(iter_results$lm_aic) <- names(iter_results$lm_bic) <- coeff_names
        return(iter_results)
      } 
      L2MTE_AIC <- parallel::mclapply(coef_LS_MTE,FUN = function(x) x$lm_aic)
      L2MTE_BIC <- parallel::mclapply(coef_LS_MTE,FUN = function(x) x$lm_bic)
      coef_L2MTE_Original <- list(
        AIC = matrix(unlist(L2MTE_AIC), nrow=num_coeffs, ncol = length(routes)), 
        BIC = matrix(unlist(L2MTE_BIC), nrow=num_coeffs, ncol = length(routes))
      )
      rownames(coef_L2MTE_Original$AIC) <- rownames(coef_L2MTE_Original$BIC) <- coeff_names
      colnames(coef_L2MTE_Original$AIC) <- colnames(coef_L2MTE_Original$BIC) <- sapply(routes,FUN = function(x) strsplit(x,"_")[[1]][1])
      save(coef_L2MTE_Original,file = "coef_L2MTE_Original.RData")
      stopCluster(cl)
      registerDoSEQ()
    },7)
  }
  if (procedures[4]){
    beepr::beep_on_error({
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      clusterExport(cl, c("MEB_DATA", "data", "f", "train_rows",  "na_result_vec","MyIC", "num_coeffs", "coeff_names", "Niter")) # Add others if needed
      cat(sprintf("Registered parallel backend with %d cores.\n", n_cores))
      
      coef_LAD_MTE <- foreach(
        selRoute = routes[routesToCheck],
        .combine = 'list',
        .multicombine = TRUE,
        .packages = c("quantreg", "RPushbullet"),
        .errorhandling = 'pass'
      ) %dopar% {
        
        lambda_sequence <- seq(0.02, 0.1, 0.02)
        n_lambda <- length(lambda_sequence)
        iter_results <- data.frame(
          lm_aic  = na_result_vec,
          lm_bic  = na_result_vec
        )
        dt <- data.frame("volume" = data[[selRoute]], "TIME" = data$TIME)[train_rows, ]
        # PROGRESS <- paste0("Starting LAD estimation: ", Sys.time()," \n"); pbPost("note", title="LAD MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
        lad.reg <- rq(as.formula(paste0("volume", f)), data = dt)
        X <- lad.reg$x
        y <- dt$volume 
        # PROGRESS <- paste0(PROGRESS, "STARTING MTE PROCESS:",Sys.time()," \n"); pbPost("note", title="LAD MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
        MTE_COEFFS <- matrix(NA, ncol=n_lambda,nrow=num_coeffs)
        for (lambda_iter in 1:length(lambda_sequence)){
          mte_result <- MTE::MTElasso(X, y, lad.reg$coefficients, intercept = FALSE, lambda = lambda_sequence[lambda_iter])
          MTE_COEFFS[, lambda_iter] <- mte_result$beta
        }
        # PROGRESS <- paste0(PROGRESS, "FINISH MTE ESTIMATION:",Sys.time()," \n"); pbPost("note", title="LAD MTE PROGRESS", PROGRESS) #write.csv(PROGRESS,file= "progress.txt")
        
        AICs <- apply(MTE_COEFFS,2, FUN = function(x) robustIC(X,x,y,"AIC"))
        BICs <- apply(MTE_COEFFS,2, FUN = function(x) robustIC(X,x,y,"BIC"))
        iter_results$lad_aic <- MTE_COEFFS[,which(AICs==min(AICs))[1]] 
        iter_results$lad_bic <- MTE_COEFFS[,which(BICs==min(BICs))[1]] 
        names(iter_results$lad_aic) <- names(iter_results$lad_bic) <- coeff_names 
        return(iter_results)
      } 
      L1MTE_AIC <- parallel::mclapply(coef_LAD_MTE,FUN = function(x) x$lad_aic)
      L1MTE_BIC <- parallel::mclapply(coef_LAD_MTE,FUN = function(x) x$lad_bic)
      coef_L1MTE_Original <- list(
        AIC = matrix(unlist(L1MTE_AIC), nrow=num_coeffs, ncol = length(routes)), 
        BIC = matrix(unlist(L1MTE_BIC), nrow=num_coeffs, ncol = length(routes))
      )
      rownames(coef_L1MTE_Original$AIC) <- rownames(coef_L1MTE_Original$BIC) <- coeff_names
      colnames(coef_L1MTE_Original$AIC) <- colnames(coef_L1MTE_Original$BIC) <- sapply(routes,FUN = function(x) strsplit(x,"_")[[1]][1])
      save(coef_L1MTE_Original,file = "coef_L1MTE_Original.RData")
      stopCluster(cl)
      registerDoSEQ()
    },7)
  }
  
  
  # 1) Construct VIPs and apply Adaptive LAD for each MEBoot Sample.
  # L1/L2 LSA/MTE Coefficients on initial data. 
  min_crit_fun <- function(X,x,y,M,C){
    if (grepl("L1",M)){
      if (grepl("A",C)){ return(robustIC(X,x,y,measure="AIC")) }else{ return(robustIC(X,x,y,measure="BIC")) }
    } else {
      if (grepl("A",C)){ return(MyIC(X,x,y,measure="AIC")) }else{ return(MyIC(X,x,y,measure="BIC")) }
    }
  }
  lambda_sequence <- seq(0.02, 0.1, 0.02)
  n_cores <- length(lambda_sequence)
  crits <- c("AIC","BIC")
  n_lambda <- length(lambda_sequence)
  X <- model.matrix(lm(paste0("L101_volume",f),data=data))
  methods <- Filter(function(x) !grepl("AdL1Lasso",x),list.files("L101",pattern="^coef.*\\.RData$"))
  
  plan(multisession)  
  {# AD_LAD_LASSO <- function(sel_route, methods, X, data, coeff_names,
  #                                    train_rows, lambda_sequence, crits, min_crit_fun) {
  #   y <- data[[sel_route]]
  #   str_sel_route <- strsplit(sel_route, "_")[[1]][1]
  #   
  #   walk(methods, function(meth) {
  #     env <- new.env()
  #     load(paste0(str_sel_route, "/", meth), envir = env)
  #     COEFFS <- env[[ls(env)]]
  #     rm(env)
  #     
  #     results <- map(crits, function(crit) {
  #       coeffs <- COEFFS[[crit]]
  #       penalty_factor <- rowMeans(coeffs != 0)
  #       
  #       vip <- penalty_factor
  #       penalty_factor[penalty_factor != 0] <- 1 / penalty_factor[penalty_factor != 0]
  #       penalty_factor <- penalty_factor[-1]
  #       
  #       if (all(penalty_factor == 0)) {
  #         penalty_factor <- rep(1, length(penalty_factor))
  #       }
  #       
  #       model <- rqPen::rq.pen(
  #         x = X[train_rows, -1],
  #         y = y[train_rows],
  #         tau = 0.5,
  #         penalty = "aLASSO",
  #         lambda = lambda_sequence,
  #         penalty.factor = penalty_factor
  #       )
  #       
  #       betas <- model$models$tau0.5a1$coefficients
  #       crit_vals <- apply(betas, 2, function(beta) min_crit_fun(X, beta, y, meth, crit))
  #       best_beta <- betas[, which(crit_vals==min(crit_vals))]
  #       
  #       list(beta = best_beta, vip = vip)
  #     })
  #     ad_l1_lasso <- do.call(cbind, map(results, "beta"))
  #     VIPs <- do.call(cbind, map(results, "vip"))
  #     colnames(ad_l1_lasso) <- colnames(VIPs) <- crits
  #     rownames(ad_l1_lasso) <- rownames(VIPs) <- coeff_names
  #     
  #     filename_ending <- strsplit(meth, "_")[[1]][2]
  #     save(VIPs, file = paste0(str_sel_route, "/VIPs_", filename_ending))
  #     save(ad_l1_lasso, file = paste0(str_sel_route, "/coef_AdL1Lasso_", filename_ending))
  #   })
  # }
  # 
  # beepr::beep_on_error(expr ={
  #   future_walk(
  #     .x = routes,
  #     .f = AD_LAD_LASSO,
  #     methods = methods,
  #     X = X,
  #     data = data,
  #     coeff_names = coeff_names,
  #     train_rows = train_rows,
  #     lambda_sequence = lambda_sequence,
  #     crits = crits,
  #     min_crit_fun = min_crit_fun
  #   )
  # },7 )
    }
  
  AD_LAD_LASSO_CV <- function(sel_route, methods, X, data, coeff_names,
                              train_rows, lambda_sequence, crits, min_crit_fun,
                              num_blocks = 7) {
    y <- data[[sel_route]]
    str_sel_route <- strsplit(sel_route, "_")[[1]][1]
    
    walk(methods, function(meth) {
      env <- new.env()
      load(paste0(str_sel_route, "/", meth), envir = env)
      COEFFS <- env[[ls(env)]]
      rm(env)
      
      results <- map(crits, function(crit) {
        coeffs <- COEFFS[[crit]]
        penalty_factor <- rowMeans(coeffs != 0)
        vip <- penalty_factor
        
        penalty_factor[penalty_factor != 0] <- 1 / penalty_factor[penalty_factor != 0]
        penalty_factor <- penalty_factor[-1]
        
        if (all(penalty_factor == 0)) {
          penalty_factor <- rep(1, length(penalty_factor))
        }
        
        n_train <- length(train_rows)
        fold_ids <- rep(1:num_blocks, length.out = n_train)
        lambda_perf <- matrix(NA, nrow = num_blocks, ncol = length(lambda_sequence))
        
        for (k in 1:num_blocks) {
          val_idx <- which(fold_ids == k)
          train_idx <- setdiff(seq_len(n_train), val_idx)
          
          X_train <- X[train_rows[train_idx], -1, drop=FALSE]
          y_train <- y[train_rows[train_idx]]
          X_val <- X[train_rows[val_idx], , drop=FALSE]
          y_val <- y[train_rows[val_idx]]
          
          model <- tryCatch({
            rqPen::rq.pen(
              x = X_train,
              y = y_train,
              tau = 0.5,
              penalty = "aLASSO",
              lambda = lambda_sequence,
              penalty.factor = penalty_factor
            )
          }, error = function(e) {
            warning(sprintf("rq.pen failed: %s", e$message))
            return(NULL)
          })
          
          if (!is.null(model) && !is.null(model$models$tau0.5a1)) {
            betas <- model$models$tau0.5a1$coefficients
            for (j in seq_along(lambda_sequence)) {
              beta <- betas[, j]
              lambda_perf[k, j] <- min_crit_fun(X_val, beta, y_val, meth, crit)
            }
          } else {
            lambda_perf[k, ] <- Inf
          }
        }
        
        mean_perf <- colMeans(lambda_perf, na.rm = TRUE)
        best_lambda_idx <- which.min(mean_perf)
        best_lambda <- lambda_sequence[best_lambda_idx]
        
        # Fit final model on all training data
        final_model <- rqPen::rq.pen(
          x = X[train_rows, -1],
          y = y[train_rows],
          tau = 0.5,
          penalty = "aLASSO",
          lambda = c(0.2,best_lambda),
          penalty.factor = penalty_factor
        )
        best_beta <- final_model$models$tau0.5a1$coefficients[, 2]
        
        list(beta = best_beta, vip = vip)
      })
      
      ad_l1_lasso <- do.call(cbind, map(results, "beta"))
      VIPs <- do.call(cbind, map(results, "vip"))
      colnames(ad_l1_lasso) <- colnames(VIPs) <- crits
      rownames(ad_l1_lasso) <- rownames(VIPs) <- coeff_names
      
      filename_ending <- strsplit(meth, "_")[[1]][2]
      save(VIPs, file = paste0(str_sel_route, "/VIPs_CV_", filename_ending))
      save(ad_l1_lasso, file = paste0(str_sel_route, "/coef_AdL1Lasso_CV_", filename_ending))
    })
  }
  
  beepr::beep_on_error(expr = {
    future_walk(
      .x = routes,
      .f = AD_LAD_LASSO_CV,
      methods = methods,
      X = X,
      data = data,
      coeff_names = coeff_names,
      train_rows = train_rows,
      lambda_sequence = lambda_sequence,
      crits = crits,
      min_crit_fun = min_crit_fun,
      num_blocks = 7
    )
  }, 7)
  
  beepr::beep(4)
}

#doProcess()

lambda_sequence <- seq(0.02, 0.1, 0.02)
n_cores <- length(lambda_sequence)
crits <- c("AIC","BIC")
n_lambda <- length(lambda_sequence)
X <- model.matrix(lm(paste0("L101_volume",f),data=data))
methods <- Filter(function(x) !grepl("AdL1Lasso",x),list.files("L101",pattern="^coef.*\\.RData$"))
selRoute <- location <- "L101_volume"
method <- methods[1]
crit <- "AIC"
method <- c(T,T)
meb_ <- F
coef_ <- T

huber_loss <- function(y_true, y_pred, delta = 1) {
  error_ <- y_true - y_pred
  is_small_error <- abs(error_) <= delta
  squared_loss <- 0.5 * error_^2
  linear_loss <- delta * (abs(error_) - 0.5 * delta)
  mean(ifelse(is_small_error, squared_loss, linear_loss))
}

mape <- function(y_pred, y_true) {
  mean(abs(y_pred / y_true)) * 100
}

mdape <- function(y_pred, y_true) {
  median(abs(y_pred / y_true)) * 100
}

smape <- function(y_pred, y_true) {
  numerator <- abs(y_pred - y_true)
  denominator <- (abs(y_true) + abs(y_pred)) / 2
  mean((numerator / denominator)* 100, na.rm = TRUE)
}

smdape <- function(y_pred, y_true) {
  numerator <- abs(y_pred - y_true)
  denominator <- (abs(y_true) + abs(y_pred)) / 2
  median((numerator / denominator)* 100, na.rm = TRUE)
}


plotData_meb <- function(location="L101_volume", method=c(T,T), crit="AIC",show_= c("Coefficients", "Fitted Values", "Residuals", "Predict Values","Prediction Errors")){
  methods <- list.files("L101",pattern="\\.RData$")
  filter_files <- "coef_"
  files <- methods[sapply(methods, FUN=function(x) (grepl(filter_files,x) & !grepl("AdL1Lasso",x)))]
  file2load <- Filter(function(x) grepl(paste0(ifelse(method[1], "L1","L2"),ifelse(method[2], "MTE","LSA")),x),files)
  split_loc <- strsplit(location,"_")[[1]][1]
  env <- new.env() 
  load(paste0(split_loc,"/",file2load),envir = env)
  COEFS <- env[[ls(env)]]
  load(paste0(split_loc,"/synth_dataset.RData"),envir = env); meb_data <- env$MEB
  rm(env)
  Coef <-  COEFS[[crit]]
  y <- meb_data
  train_rows <- 1:floor(nrow(X)*0.5)
  
  test_rows <- (floor(nrow(X)*0.5)+1):floor(nrow(X)*0.75)
  X_train <- X[train_rows,]
  X_test  <- X[test_rows,]
  y_train <- y[train_rows,]
  y_test  <- y[test_rows,]
  
  dt <- switch(show_,
               "Coefficients" = data.frame(Values = Coef, Predictor = colnames(X)),
               "Fitted Values" = data.frame(Values = X_train%*%Coef, Predictor = data$DATE[train_rows]),
               "Residuals" = data.frame(Values = y_train-X_train%*%Coef, Predictor = data$DATE[train_rows]),
               "Predict Values" = data.frame(Values = X_test%*%Coef, Predictor = data$DATE[test_rows]),
               "Prediction Errors" = data.frame(Values = y_test-X_test%*%Coef, Predictor = data$DATE[test_rows]),
               stop("Invalid `show_` input. Please, select one of the following:\n 'Coefficients', 'Fitted Values', 'Residuals', 'Predict Values', or 'Prediction Errors'.")
               )
  colnames(dt) <- c(paste0("MEB_",1:20), "Predictor")
  return(dt)
}
  
plotData <- function(location, method, crit, meb_, show_=c("Coefficients", "Fitted Values", "Residuals", "Predict Values","Prediction Errors")){
  methods_meb <- list.files("L101",pattern="\\.RData$")
  methods <- list.files(pattern="\\.RData$")
  methods <- c(methods_meb, methods)
  filter_files <- ifelse(meb_,"^coef_AdL1Lasso_","^coef_.*_Original\\.RData$")
  files <- methods[sapply(methods, FUN=function(x) (grepl(filter_files,x) ))]
  file2load <- Filter(function(x) grepl(paste0(ifelse(method[1], "L1","L2"),ifelse(method[2], "MTE","LSA")),x),files)
  filepath2load <- ifelse(meb_,paste0(strsplit(location,"_")[[1]][1],"/",file2load),file2load)
  env <- new.env() 
  load(filepath2load,envir = env)
  COEFS <- env[[ls(env)]]
  rm(env)
  if (meb_){
    Coef <-  COEFS[,crit]
  }else{
    Coef <-  COEFS[[crit]]
  }
  y <- data[[location]]
  train_rows <- 1:floor(nrow(X)*0.5)
  test_rows <- (floor(nrow(X)*0.5)+1):floor(nrow(X)*0.75)
  X_train <- X[train_rows,]
  X_test  <- X[test_rows,]
  y_train <- y[train_rows]
  y_test  <- y[test_rows]
  locations <- sapply(routes, function(l) strsplit(l,"_")[[1]][1])
  dt <- switch(show_,
               "Coefficients" = data.frame(Values = Coef, Predictor = colnames(X)), 
               "Fitted Values" = data.frame(Values = X_train%*%Coef, Predictor = data$DATE[train_rows]),
               "Residuals" = data.frame(Values = y_train-X_train%*%Coef, Predictor = data$DATE[train_rows]),
               "Predict Values" = data.frame(Values = X_test%*%Coef, Predictor = data$DATE[test_rows]),
               "Prediction Errors" = data.frame(Values = y_test-X_test%*%Coef, Predictor = data$DATE[test_rows]),
               stop("Invalid `show_` input. Please, select one of the following:\n 'Coefficients', 'Fitted Values', 'Residuals', 'Predict Values', or 'Prediction Errors'.")
  )
  if (ncol(dt) >2) {colnames(dt) <- c(locations, "Predictor")}
  return(dt)
}

plotData_VIPs <- function(location, method, crit){
  files <- list.files("L101",pattern="^VIPs_.*\\.RData$")
  file2load <- Filter(function(x) grepl(paste0(ifelse(method[1], "L1","L2"),ifelse(method[2], "MTE","LSA")),x),files) #file2load <- files[(method%*%2^(0:1))+1]
  filepath2load <- paste0(strsplit(location,"_")[[1]][1],"/",file2load)
  env <- new.env() 
  load(filepath2load,envir = env)
  VIM <- env[[ls(env)]]
  rm(env)
  VIPs <-  VIM[,crit]
  dt <- data.frame(Values = VIPs, Predictor = colnames(X))
  return(dt)
}

