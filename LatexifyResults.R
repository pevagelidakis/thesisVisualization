library(shiny)
library(shinydashboard)
library(plotly)
library(leaflet)
library(reshape2)
library(leaflet.extras)
library(shinyWidgets)
library(Polychrome)
library(dplyr)
library(tseries)
library(DT)
library(purrr)
setwd(choose.dir())
source("../../Scripts/L1L2_LSA_MTE_Estimations.R")

METRICS <- function(location, method, crit){
  split_loc <- strsplit(location,"_")[[1]][1]
  env <- new.env() 
  load(paste0(split_loc,"/synth_dataset.RData"),envir = env); meb_data <- env$MEB
  rm(env)
  train_rows <- 1:(dim(data)[1]/2)
  test_rows <- (dim(data)[1]/2 + 1):(dim(data)[1]*3/4)
  
  y_pred_tr <- cbind(plotData_meb(location, method, crit, show_="Fitted Values")[,-21], plotData(location,method,crit, F,show_="Fitted Values")[,split_loc], plotData(location,method,crit, T,show_="Fitted Values")[,1])
  y_pred_te <- cbind(plotData_meb(location, method, crit, show_="Predict Values")[,-21], plotData(location,method,crit, F,show_="Predict Values")[,split_loc], plotData(location,method,crit, T,show_="Predict Values")[,1])
  y_true_tr <- cbind(meb_data[train_rows,], data[train_rows,location], data[train_rows,location])
  y_true_te <- cbind(meb_data[test_rows,], data[test_rows,location], data[test_rows,location])
  
  mae_tr <- colMeans(abs(y_pred_tr-y_true_tr))
  rmse_tr <- sqrt(colMeans((y_pred_tr-y_true_tr)^2))
  huber_loss_tr <- mapply(
    function(col_true, col_pred) huber_loss(col_true, col_pred),
    as.data.frame(y_true_tr),
    as.data.frame(y_pred_tr)
  )
  mape_tr <- mapply(
    function(col_true, col_pred) mape(col_true, col_pred),
    as.data.frame(y_true_tr),
    as.data.frame(y_pred_tr)
  )
  mdape_tr <- mapply(
    function(col_true, col_pred) mdape(col_true, col_pred),
    as.data.frame(y_true_tr),
    as.data.frame(y_pred_tr)
  )
  smape_tr <- mapply(
    function(col_true, col_pred) smape(col_true, col_pred),
    as.data.frame(y_true_tr),
    as.data.frame(y_pred_tr)
  )
  smdape_tr <- mapply(
    function(col_true, col_pred) smdape(col_true, col_pred),
    as.data.frame(y_true_tr),
    as.data.frame(y_pred_tr)
  )
  
  
  mae_te <- colMeans(abs(y_pred_te-y_true_te))
  rmse_te <- sqrt(colMeans((y_pred_tr-y_true_tr)^2))
  huber_loss_te <- mapply(
    function(col_true, col_pred) huber_loss(col_true, col_pred),
    as.data.frame(y_true_te),
    as.data.frame(y_pred_te)
  )
  mape_te <- mapply(
    function(col_true, col_pred) mape(col_true, col_pred),
    as.data.frame(y_true_te),
    as.data.frame(y_pred_te)
  )
  mdape_te <- mapply(
    function(col_true, col_pred) mdape(col_true, col_pred),
    as.data.frame(y_true_te),
    as.data.frame(y_pred_te)
  )
  smape_te <- mapply(
    function(col_true, col_pred) smape(col_true, col_pred),
    as.data.frame(y_true_te),
    as.data.frame(y_pred_te)
  )
  smdape_te <- mapply(
    function(col_true, col_pred) smdape(col_true, col_pred),
    as.data.frame(y_true_te),
    as.data.frame(y_pred_te)
  )
  list(
    data = data.frame(
      Variable = c(rep("MAE_TRAIN",20),rep("MAE_TEST",20), rep("RMSE_TRAIN",20),rep("RMSE_TEST",20),
                   rep("HUBER_LOSS_TRAIN",20),rep("HUBER_LOSS_TEST",20), rep("MAPE_TRAIN",20),rep("MAPE_TEST",20),
                   rep("MdAPE_TRAIN",20),rep("MdAPE_TEST",20), rep("SMAPE_TRAIN",20),rep("SMAPE_TEST",20),
                   rep("SMdAPE_TRAIN",20),rep("SMdAPE_TEST",20)),
      Value = c(mae_tr[-(21:22)], mae_te[-(21:22)], rmse_tr[-(21:22)], rmse_te[-(21:22)],
                huber_loss_tr[-(21:22)], huber_loss_te[-(21:22)], mape_tr[-(21:22)],mape_tr[-(21:22)],
                mdape_tr[-(21:22)], mdape_te[-(21:22)],smape_tr[-(21:22)], smape_tr[-(21:22)],
                smdape_tr[-(21:22)], smdape_te[-(21:22)])),
    artifacts = data.frame(
      Actual = c(mae_tr[21], mae_te[21], rmse_tr[21], rmse_te[21],
                 huber_loss_tr[21], huber_loss_te[21], mape_tr[21],mape_te[21],
                 mdape_tr[21], mdape_te[21],smape_tr[21], smape_te[21],
                 smdape_tr[21], smdape_te[21]),
      ADL1LASSO = c(mae_tr[22], mae_te[22], rmse_tr[22], rmse_te[22],
                    huber_loss_tr[22], huber_loss_te[22], mape_tr[22],mape_te[22],
                    mdape_tr[22], mdape_te[22],smape_tr[22], smape_te[22],
                    smdape_tr[22], smdape_te[22])
    ) 
  )
}
RESULTS_COLLECTION <- function(location, method, crit){
  dist_data <- METRICS(location, method, crit)
  train_rows <- 1:(dim(data)[1]/2)
  test_rows <- (dim(data)[1]/2 + 1):(dim(data)[1]*3/4)
  meb_coef <- rowMeans(plotData_meb(location, method, crit, show_="Coefficients")[,-21])
  meb_fit <- meb_coef%*%t(X)
  meb_pred <- meb_coef%*%t(X)
  y_train <-data[train_rows,location]
  y_test <-data[test_rows,location]
  meb_res <- data[train_rows,location] - meb_fit
  meb_err <- data[test_rows,location]  - meb_pred
  
  plot_data <- dist_data$data
  CI95_mae_tr         <- round(quantile(plot_data[plot_data$Variable=="MAE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
  CI95_mae_te         <- round(quantile(plot_data[plot_data$Variable=="MAE_TEST","Value"],c(0.05,0.5,0.95)),3)
  CI95_rmse_tr        <- round(quantile(plot_data[plot_data$Variable=="RMSE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
  CI95_rmse_te        <- round(quantile(plot_data[plot_data$Variable=="RMSE_TEST","Value"],c(0.05,0.5,0.95)),3)
  CI95_huber_loss_tr  <- round(quantile(plot_data[plot_data$Variable=="HUBER_LOSS_TRAIN","Value"],c(0.05,0.5,0.95)),3)
  CI95_huber_loss_te  <- round(quantile(plot_data[plot_data$Variable=="HUBER_LOSS_TEST","Value"],c(0.05,0.5,0.95)),3)
  CI95_mape_tr        <- round(quantile(plot_data[plot_data$Variable=="MAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
  CI95_mape_te        <- round(quantile(plot_data[plot_data$Variable=="MAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
  CI95_mdape_tr       <- round(quantile(plot_data[plot_data$Variable=="MdAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
  CI95_mdape_te       <- round(quantile(plot_data[plot_data$Variable=="MdAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
  CI95_smape_tr       <- round(quantile(plot_data[plot_data$Variable=="SMAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
  CI95_smape_te       <- round(quantile(plot_data[plot_data$Variable=="SMAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
  CI95_smdape_tr      <- round(quantile(plot_data[plot_data$Variable=="SMdAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
  CI95_smdape_te      <- round(quantile(plot_data[plot_data$Variable=="SMdAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
  
  artifacts <- data.frame(
    "MAE_TRAIN"         = c(CI95_mae_tr, mean(plot_data[plot_data$Variable=="MAE_TRAIN","Value"]),               round(mean(abs(meb_res)),3),          round(dist_data$artifacts$Actual[1],3),  round(dist_data$artifacts$ADL1LASSO[1],3)),
    "MAE_TEST"          = c(CI95_mae_te, mean(plot_data[plot_data$Variable=="MAE_TEST","Value"]),                round(mean(abs(meb_err)),3),          round(dist_data$artifacts$Actual[2],3),  round(dist_data$artifacts$ADL1LASSO[2],3)),
    "RMSE_TRAIN"        = c(CI95_rmse_tr, mean(plot_data[plot_data$Variable=="RMSE_TRAIN","Value"]),             round(sqrt(mean(meb_res^2)),3),       round(dist_data$artifacts$Actual[3],3),  round(dist_data$artifacts$ADL1LASSO[3],3)),
    "RMSE_TEST"         = c(CI95_rmse_te, mean(plot_data[plot_data$Variable=="RMSE_TEST","Value"]),              round(sqrt(mean(meb_err^2)),3),       round(dist_data$artifacts$Actual[4],3),  round(dist_data$artifacts$ADL1LASSO[4],3)),
    "HUBER_LOSS_TRAIN"  = c(CI95_huber_loss_tr, mean(plot_data[plot_data$Variable=="HUBER_LOSS_TRAIN","Value"]), round(huber_loss(y_train,meb_fit),3), round(dist_data$artifacts$Actual[5],3),  round(dist_data$artifacts$ADL1LASSO[5],3)),
    "HUBER_LOSS_TEST"   = c(CI95_huber_loss_te, mean(plot_data[plot_data$Variable=="HUBER_LOSS_TEST","Value"]),  round(huber_loss(y_test,meb_pred),3), round(dist_data$artifacts$Actual[6],3),  round(dist_data$artifacts$ADL1LASSO[6],3)),
    "MAPE_TRAIN"        = c(CI95_mape_tr, mean(plot_data[plot_data$Variable=="MAPE_TRAIN","Value"]),             round(mape(y_train,meb_fit),3),       round(dist_data$artifacts$Actual[7],3),  round(dist_data$artifacts$ADL1LASSO[7],3)),
    "MAPE_TEST"         = c(CI95_mape_te, mean(plot_data[plot_data$Variable=="MAPE_TEST","Value"]),              round(mape(y_test,meb_pred),3),       round(dist_data$artifacts$Actual[8],3),  round(dist_data$artifacts$ADL1LASSO[8],3)),
    "MdAPE_TRAIN"       = c(CI95_mdape_tr, mean(plot_data[plot_data$Variable=="MdAPE_TRAIN","Value"]),           round(mdape(y_train,meb_fit),3),      round(dist_data$artifacts$Actual[9],3),  round(dist_data$artifacts$ADL1LASSO[9],3)),
    "MdAPE_TEST"        = c(CI95_mdape_te, mean(plot_data[plot_data$Variable=="MdAPE_TEST","Value"]),            round(mdape(y_test,meb_pred),3),      round(dist_data$artifacts$Actual[10],3), round(dist_data$artifacts$ADL1LASSO[10],3)),
    "SMAPE_TRAIN"       = c(CI95_smape_tr, mean(plot_data[plot_data$Variable=="SMAPE_TRAIN","Value"]),           round(smape(y_train,meb_fit),3),      round(dist_data$artifacts$Actual[11],3), round(dist_data$artifacts$ADL1LASSO[11],3)),
    "SMAPE_TEST"        = c(CI95_smape_te, mean(plot_data[plot_data$Variable=="SMAPE_TEST","Value"]),            round(smape(y_test,meb_pred),3),      round(dist_data$artifacts$Actual[12],3), round(dist_data$artifacts$ADL1LASSO[12],3)),
    "SMdAPE_TRAIN"      = c(CI95_smdape_tr, mean(plot_data[plot_data$Variable=="SMdAPE_TRAIN","Value"]),         round(smdape(y_train,meb_fit),3),     round(dist_data$artifacts$Actual[13],3), round(dist_data$artifacts$ADL1LASSO[13],3)),
    "SMdAPE_TEST"       = c(CI95_smdape_te, mean(plot_data[plot_data$Variable=="SMdAPE_TEST","Value"]),          round(smdape(y_test,meb_pred),3),     round(dist_data$artifacts$Actual[14],3), round(dist_data$artifacts$ADL1LASSO[14],3))
  )
  
  artifacts_long <- artifacts %>%
    t() %>%
    as.data.frame() %>%
    mutate(Variable = rownames(.)) %>%
    rename(
      `5%` = V1,
      `50%` = V2,
      `95%` = V3,
      mean = V4,
      modelAvg= V5,
      Actual = V6,
      ADL1LASSO = V7
    )
  var_levels <- unique(plot_data$Variable)
  var_map <- setNames(seq_along(var_levels), var_levels)
  plot_data <- plot_data %>%
    mutate(x = var_map[Variable])
  artifacts_long <- artifacts_long %>%
    mutate(x = var_map[Variable])
  
  return(artifacts_long)
}

storeThis <- ""
for (i in c(TRUE, FALSE)) {
  for (j in c(TRUE, FALSE)) {
    method <- c(i,j)
    storeThis <- c(storeThis, paste0("\nFor ", paste0(ifelse(i, "L2", "L1"), ifelse(j, "MTE", "LSA")), " method:\n"))
    for (crit in crits) {
      storeThis <- c(storeThis, paste0("\tFor ", crit, "\n"))
      for (location in routes) {
        storeThis <- c(storeThis, paste0("\t\tFor ", location, "\n"))
        result <- RESULTS_COLLECTION(location, method, crit)
        if (!is.character(result)) {
          result <- capture.output(print(result))
        }
        storeThis <- c(storeThis, paste0("\t\t\t", result, collapse = "\n"))
      }
    }
  }
}
final_output <- paste(storeThis, collapse = "\n")

write.csv(final_output,file= "storeThis.txt")


L1MTE <- L2MTE <- L1LSA <- L2LSA <- list(AIC=list(),BIC = list())
for (i in c(TRUE, FALSE)) {
  for (j in c(TRUE, FALSE)) {
    method <- c(i,j)
    if ((i) && (j)){
      for (crit in crits) {
        for (location in routes) {
          L2MTE[[crit]][[location]] <- RESULTS_COLLECTION(location, method, crit)
        }
      }
    }
    if ((i) && (!j)){
      for (crit in crits) {
        for (location in routes) {
          L2LSA[[crit]][[location]] <- RESULTS_COLLECTION(location, method, crit)
        }
      }
    }
    if ((!i) && (j)){
      for (crit in crits) {
        for (location in routes) {
          L1MTE[[crit]][[location]] <- RESULTS_COLLECTION(location, method, crit)
        }
      }
    }
    if ((!i) && (!j)){
      for (crit in crits) {
        for (location in routes) {
          L1LSA[[crit]][[location]] <- RESULTS_COLLECTION(location, method, crit)
        }
      }
    }
  }
}
train_metrics <- Filter(function(x) grepl("TRAIN",x), L1LSA[[crit]][[location]]$Variable)
test_metrics <- Filter(function(x) grepl("TEST",x), L1LSA[[crit]][[location]]$Variable)
metrics_dashboard = list(L101_volume = list(TRAIN= list(),TEST=list()),
                         L102_volume = list(TRAIN= list(),TEST=list()),
                         L103_volume = list(TRAIN= list(),TEST=list()),
                         L104_volume = list(TRAIN= list(),TEST=list()),
                         L106_volume = list(TRAIN= list(),TEST=list()),
                         L107_volume = list(TRAIN= list(),TEST=list()),
                         L108_volume = list(TRAIN= list(),TEST=list())
                    )
for (location in routes){
  meb_l1lsa_aic <- L1LSA[["AIC"]][[location]]$modelAvg[which(L1LSA[["AIC"]][[location]]$Variable %in%train_metrics)]
  act_l1lsa_aic <- L1LSA[["AIC"]][[location]]$Actual[which(L1LSA[["AIC"]][[location]]$Variable %in%train_metrics)]
  all_l1lsa_aic <- L1LSA[["AIC"]][[location]]$ADL1LASSO[which(L1LSA[["AIC"]][[location]]$Variable %in%train_metrics)]
  meb_l2lsa_aic <- L2LSA[["AIC"]][[location]]$modelAvg[which(L2LSA[["AIC"]][[location]]$Variable %in%train_metrics)]
  act_l2lsa_aic <- L2LSA[["AIC"]][[location]]$Actual[which(L2LSA[["AIC"]][[location]]$Variable %in%train_metrics)]
  all_l2lsa_aic <- L2LSA[["AIC"]][[location]]$ADL1LASSO[which(L2LSA[["AIC"]][[location]]$Variable %in%train_metrics)]
  meb_l1mte_aic <- L1MTE[["AIC"]][[location]]$modelAvg[which(L1MTE[["AIC"]][[location]]$Variable %in%train_metrics)]
  act_l1mte_aic <- L1MTE[["AIC"]][[location]]$Actual[which(L1MTE[["AIC"]][[location]]$Variable %in%train_metrics)]
  all_l1mte_aic <- L1MTE[["AIC"]][[location]]$ADL1LASSO[which(L1MTE[["AIC"]][[location]]$Variable %in%train_metrics)]
  meb_l2mte_aic <- L2MTE[["AIC"]][[location]]$modelAvg[which(L2MTE[["AIC"]][[location]]$Variable %in%train_metrics)]
  act_l2mte_aic <- L2MTE[["AIC"]][[location]]$Actual[which(L2MTE[["AIC"]][[location]]$Variable %in%train_metrics)]
  all_l2mte_aic <- L2MTE[["AIC"]][[location]]$ADL1LASSO[which(L2MTE[["AIC"]][[location]]$Variable %in%train_metrics)]
  
  meb_l1lsa_bic <- L1LSA[["BIC"]][[location]]$modelAvg[which(L1LSA[["BIC"]][[location]]$Variable %in%train_metrics)]
  act_l1lsa_bic <- L1LSA[["BIC"]][[location]]$Actual[which(L1LSA[["BIC"]][[location]]$Variable %in%train_metrics)]
  all_l1lsa_bic <- L1LSA[["BIC"]][[location]]$ADL1LASSO[which(L1LSA[["BIC"]][[location]]$Variable %in%train_metrics)]
  meb_l2lsa_bic <- L2LSA[["BIC"]][[location]]$modelAvg[which(L2LSA[["BIC"]][[location]]$Variable %in%train_metrics)]
  act_l2lsa_bic <- L2LSA[["BIC"]][[location]]$Actual[which(L2LSA[["BIC"]][[location]]$Variable %in%train_metrics)]
  all_l2lsa_bic <- L2LSA[["BIC"]][[location]]$ADL1LASSO[which(L2LSA[["BIC"]][[location]]$Variable %in%train_metrics)]
  meb_l1mte_bic <- L1MTE[["BIC"]][[location]]$modelAvg[which(L1MTE[["BIC"]][[location]]$Variable %in%train_metrics)]
  act_l1mte_bic <- L1MTE[["BIC"]][[location]]$Actual[which(L1MTE[["BIC"]][[location]]$Variable %in%train_metrics)]
  all_l1mte_bic <- L1MTE[["BIC"]][[location]]$ADL1LASSO[which(L1MTE[["BIC"]][[location]]$Variable %in%train_metrics)]
  meb_l2mte_bic <- L2MTE[["BIC"]][[location]]$modelAvg[which(L2MTE[["BIC"]][[location]]$Variable %in%train_metrics)]
  act_l2mte_bic <- L2MTE[["BIC"]][[location]]$Actual[which(L2MTE[["BIC"]][[location]]$Variable %in%train_metrics)]
  all_l2mte_bic <- L2MTE[["BIC"]][[location]]$ADL1LASSO[which(L2MTE[["BIC"]][[location]]$Variable %in%train_metrics)]
  
  
  TRAIN = data.frame(rbind(meb_l1lsa_aic, act_l1lsa_aic, all_l1lsa_aic, meb_l2lsa_aic, act_l2lsa_aic, all_l2lsa_aic, meb_l1mte_aic, act_l1mte_aic, all_l1mte_aic, meb_l2mte_aic, act_l2mte_aic, all_l2mte_aic, meb_l1lsa_bic, act_l1lsa_bic, all_l1lsa_bic, meb_l2lsa_bic, act_l2lsa_bic, all_l2lsa_bic, meb_l1mte_bic, act_l1mte_bic, all_l1mte_bic, meb_l2mte_bic, act_l2mte_bic, all_l2mte_bic))
  colnames(TRAIN) = train_metrics
  finalTrain <- data.frame("HL" = TRAIN$HUBER_LOSS_TRAIN,
                           "MAE" = TRAIN$MAE_TRAIN,
                           "RMSE" = TRAIN$RMSE_TRAIN,
                           "MAPE" = TRAIN$MAPE_TRAIN,
                           "MdAPE" = TRAIN$MdAPE_TRAIN,
                           "SMAPE" = TRAIN$SMAPE_TRAIN,
                           "SMdAPE" = TRAIN$SMdAPE_TRAIN)
  rownames(finalTrain) = c("MEB_L1LSA_AIC", "Actual_L1LSA_AIC", "ADL1LASSO_L1LSA_AIC",
                           "MEB_L2LSA_AIC", "Actual_L2LSA_AIC", "ADL1LASSO_L2LSA_AIC",
                           "MEB_L1MTE_AIC", "Actual_L1MTE_AIC", "ADL1LASSO_L1MTE_AIC",
                           "MEB_L2MTE_AIC", "Actual_L2MTE_AIC", "ADL1LASSO_L2MTE_AIC",
                           "MEB_L1LSA_BIC", "Actual_L1LSA_BIC", "ADL1LASSO_L1LSA_BIC",
                           "MEB_L2LSA_BIC", "Actual_L2LSA_BIC", "ADL1LASSO_L2LSA_BIC",
                           "MEB_L1MTE_BIC", "Actual_L1MTE_BIC", "ADL1LASSO_L1MTE_BIC",
                           "MEB_L2MTE_BIC", "Actual_L2MTE_BIC", "ADL1LASSO_L2MTE_BIC")
  

  meb_l1lsa_aic <- L1LSA[["AIC"]][[location]]$modelAvg[which(L1LSA[["AIC"]][[location]]$Variable %in%test_metrics)]
  act_l1lsa_aic <- L1LSA[["AIC"]][[location]]$Actual[which(L1LSA[["AIC"]][[location]]$Variable %in%test_metrics)]
  all_l1lsa_aic <- L1LSA[["AIC"]][[location]]$ADL1LASSO[which(L1LSA[["AIC"]][[location]]$Variable %in%test_metrics)]
  meb_l2lsa_aic <- L2LSA[["AIC"]][[location]]$modelAvg[which(L2LSA[["AIC"]][[location]]$Variable %in%test_metrics)]
  act_l2lsa_aic <- L2LSA[["AIC"]][[location]]$Actual[which(L2LSA[["AIC"]][[location]]$Variable %in%test_metrics)]
  all_l2lsa_aic <- L2LSA[["AIC"]][[location]]$ADL1LASSO[which(L2LSA[["AIC"]][[location]]$Variable %in%test_metrics)]
  meb_l1mte_aic <- L1MTE[["AIC"]][[location]]$modelAvg[which(L1MTE[["AIC"]][[location]]$Variable %in%test_metrics)]
  act_l1mte_aic <- L1MTE[["AIC"]][[location]]$Actual[which(L1MTE[["AIC"]][[location]]$Variable %in%test_metrics)]
  all_l1mte_aic <- L1MTE[["AIC"]][[location]]$ADL1LASSO[which(L1MTE[["AIC"]][[location]]$Variable %in%test_metrics)]
  meb_l2mte_aic <- L2MTE[["AIC"]][[location]]$modelAvg[which(L2MTE[["AIC"]][[location]]$Variable %in%test_metrics)]
  act_l2mte_aic <- L2MTE[["AIC"]][[location]]$Actual[which(L2MTE[["AIC"]][[location]]$Variable %in%test_metrics)]
  all_l2mte_aic <- L2MTE[["AIC"]][[location]]$ADL1LASSO[which(L2MTE[["AIC"]][[location]]$Variable %in%test_metrics)]
  
  meb_l1lsa_bic <- L1LSA[["BIC"]][[location]]$modelAvg[which(L1LSA[["BIC"]][[location]]$Variable %in%test_metrics)]
  act_l1lsa_bic <- L1LSA[["BIC"]][[location]]$Actual[which(L1LSA[["BIC"]][[location]]$Variable %in%test_metrics)]
  all_l1lsa_bic <- L1LSA[["BIC"]][[location]]$ADL1LASSO[which(L1LSA[["BIC"]][[location]]$Variable %in%test_metrics)]
  meb_l2lsa_bic <- L2LSA[["BIC"]][[location]]$modelAvg[which(L2LSA[["BIC"]][[location]]$Variable %in%test_metrics)]
  act_l2lsa_bic <- L2LSA[["BIC"]][[location]]$Actual[which(L2LSA[["BIC"]][[location]]$Variable %in%test_metrics)]
  all_l2lsa_bic <- L2LSA[["BIC"]][[location]]$ADL1LASSO[which(L2LSA[["BIC"]][[location]]$Variable %in%test_metrics)]
  meb_l1mte_bic <- L1MTE[["BIC"]][[location]]$modelAvg[which(L1MTE[["BIC"]][[location]]$Variable %in%test_metrics)]
  act_l1mte_bic <- L1MTE[["BIC"]][[location]]$Actual[which(L1MTE[["BIC"]][[location]]$Variable %in%test_metrics)]
  all_l1mte_bic <- L1MTE[["BIC"]][[location]]$ADL1LASSO[which(L1MTE[["BIC"]][[location]]$Variable %in%test_metrics)]
  meb_l2mte_bic <- L2MTE[["BIC"]][[location]]$modelAvg[which(L2MTE[["BIC"]][[location]]$Variable %in%test_metrics)]
  act_l2mte_bic <- L2MTE[["BIC"]][[location]]$Actual[which(L2MTE[["BIC"]][[location]]$Variable %in%test_metrics)]
  all_l2mte_bic <- L2MTE[["BIC"]][[location]]$ADL1LASSO[which(L2MTE[["BIC"]][[location]]$Variable %in%test_metrics)]
  
  TEST = data.frame(rbind(meb_l1lsa_aic, act_l1lsa_aic, all_l1lsa_aic, meb_l2lsa_aic, act_l2lsa_aic, all_l2lsa_aic, meb_l1mte_aic, act_l1mte_aic, all_l1mte_aic, meb_l2mte_aic, act_l2mte_aic, all_l2mte_aic, meb_l1lsa_bic, act_l1lsa_bic, all_l1lsa_bic, meb_l2lsa_bic, act_l2lsa_bic, all_l2lsa_bic, meb_l1mte_bic, act_l1mte_bic, all_l1mte_bic, meb_l2mte_bic, act_l2mte_bic, all_l2mte_bic))
  
  colnames(TEST) = test_metrics
  finalTest <- data.frame( "HL" = TEST$HUBER_LOSS_TEST,
                           "MAE" = TEST$MAE_TEST,
                           "RMSE" = TEST$RMSE_TEST,
                           "MAPE" = TEST$MAPE_TEST,
                           "MdAPE" = TEST$MdAPE_TEST,
                           "SMAPE" = TEST$SMAPE_TEST,
                           "SMdAPE" = TEST$SMdAPE_TEST)
  rownames(finalTest) = c("MEB_L1LSA_AIC", "Actual_L1LSA_AIC", "ADL1LASSO_L1LSA_AIC",
                          "MEB_L2LSA_AIC", "Actual_L2LSA_AIC", "ADL1LASSO_L2LSA_AIC",
                          "MEB_L1MTE_AIC", "Actual_L1MTE_AIC", "ADL1LASSO_L1MTE_AIC",
                          "MEB_L2MTE_AIC", "Actual_L2MTE_AIC", "ADL1LASSO_L2MTE_AIC",
                          "MEB_L1LSA_BIC", "Actual_L1LSA_BIC", "ADL1LASSO_L1LSA_BIC",
                          "MEB_L2LSA_BIC", "Actual_L2LSA_BIC", "ADL1LASSO_L2LSA_BIC",
                          "MEB_L1MTE_BIC", "Actual_L1MTE_BIC", "ADL1LASSO_L1MTE_BIC",
                          "MEB_L2MTE_BIC", "Actual_L2MTE_BIC", "ADL1LASSO_L2MTE_BIC")
  metrics_dashboard[[location]]$TEST <- finalTest
  metrics_dashboard[[location]]$TRAIN <- finalTrain
}

boldMin <- function(df,mins){
  edited <- matrix(NA, nrow=dim(df)[1],ncol=dim(df)[2])
  for (column in seq_len(ncol(df))){
    edited[,column] <- sapply(df[,column], FUN = function(l) ifelse(!(l %in% mins[,column]), sprintf("%.2f", as.numeric(l)), paste0("\\textbf{",sprintf("%.2f", as.numeric(l)),"}")))
  }
  edited
}

generate_latex_table <- function(metrics_dashboard, location) {
  train_df <- metrics_dashboard[[location]]$TRAIN
  test_df <- metrics_dashboard[[location]]$TEST
  
  get_row <- function(name) {
    min_tr <- rbind(apply(train_df[seq_len(nrow(train_df)) %% 3 == 1, ], 2, min), apply(train_df[seq_len(nrow(train_df)) %% 3 == 2, ], 2, min), apply(train_df[seq_len(nrow(train_df)) %% 3 == 0, ], 2, min))
    min_te <- rbind(apply(test_df[seq_len(nrow(test_df)) %% 3 == 1, ], 2, min),apply(test_df[seq_len(nrow(test_df)) %% 3 == 2, ], 2, min),apply(test_df[seq_len(nrow(test_df)) %% 3 == 0, ], 2, min))
    
    train_vals <- boldMin(train_df[name, ], min_tr)
    test_vals <- boldMin(test_df[name, ], min_te)
    paste("  \\ \\ \\  \\ \\ \\ \\textbf{", gsub("_.*","",name), "} &", 
          paste(c(train_vals, test_vals), collapse = " & "), "\\\\")
  }
  
  # Order of models under AIC/BIC
  methods <- list(
    AIC = c("MEB_L1LSA_AIC", "Actual_L1LSA_AIC", "ADL1LASSO_L1LSA_AIC",
            "MEB_L2LSA_AIC", "Actual_L2LSA_AIC", "ADL1LASSO_L2LSA_AIC",
            "MEB_L1MTE_AIC", "Actual_L1MTE_AIC", "ADL1LASSO_L1MTE_AIC",
            "MEB_L2MTE_AIC", "Actual_L2MTE_AIC", "ADL1LASSO_L2MTE_AIC"),
    BIC = c("MEB_L1LSA_BIC", "Actual_L1LSA_BIC", "ADL1LASSO_L1LSA_BIC",
            "MEB_L2LSA_BIC", "Actual_L2LSA_BIC", "ADL1LASSO_L2LSA_BIC",
            "MEB_L1MTE_BIC", "Actual_L1MTE_BIC", "ADL1LASSO_L1MTE_BIC",
            "MEB_L2MTE_BIC", "Actual_L2MTE_BIC", "ADL1LASSO_L2MTE_BIC")
  )
  
  model_names <- list(
    L1MTE = c("MEB", "Actual", "ADL1LASSO"),
    L2MTE = c("MEB", "Actual", "ADL1LASSO"),
    L1LSA = c("MEB", "Actual", "ADL1LASSO"),
    L2LSA = c("MEB", "Actual", "ADL1LASSO")
  )
  
  latex <- c(paste("\\begin{table}[H]
\\centering
\\resizebox{\\textwidth}{!}{%
\\begin{tabular}{l||l*{6}{c}|*{7}{c}}
\\toprule
\\multicolumn{2}{l}{\\multirow{3}{*}{\\textbf{\\huge{Methods}}}} & 
\\multicolumn{13}{c}{\\textbf{Location:", strsplit(location,"_")[[1]][1],"}} \\\\
\\cmidrule(lr){2-15}
& \\multicolumn{7}{c}{\\textbf{TRAIN}} & \\multicolumn{7}{c}{\\textbf{TEST}} \\\\
\\cmidrule(lr){2-8} \\cmidrule(lr){9-15}
&  \\textbf{HL} & \\textbf{MAE} & \\textbf{RMSE} & \\textbf{MAPE} & \\textbf{MdAPE} & \\textbf{SMAPE} & \\textbf{SMdAPE} &  \\textbf{HL} & \\textbf{MAE} & \\textbf{RMSE} & \\textbf{MAPE} & \\textbf{MdAPE} & \\textbf{SMAPE} & \\textbf{SMdAPE} \\\\
\\midrule"))
  
  for (criterion in names(methods)) {
    latex <- c(latex, "\\noalign{\\vskip 0.1cm}")
    latex <- c(latex, paste0("\\multirow{4}{*}{\\rotatebox{90}{\\parbox{5.3cm}{\\centering \\textbf{\\large{", criterion, "}}}}}"))
    models <- unique(gsub(".*_(L[12][A-Z]+)_.*", "\\1", methods[[criterion]]))
    for (model in models) {
      rows <- methods[[criterion]][grepl(model, methods[[criterion]])]
      prefix <- ifelse(model!="L1LSA","\\ \\ \\ ","")
      latex <- c(latex, paste0(prefix,"\\multirow{4}{*}{\\rotatebox{90}{\\parbox{0.6cm}{\\centering \\textbf{", model, "}}}}"))
      for (row in rows) {
        latex <- c(latex,get_row(row))
      }
      latex <- c(latex, "\\noalign{\\vskip 0.1cm}","\\cmidrule{2-15}","\\noalign{\\vskip 0.1cm}")
    }
    latex <- c(latex, "\\midrule")
  }
  
  latex <- c(latex, "\\bottomrule", "\\end{tabular}}", "\\caption{\\textit{Comparison of seven performance metrics (HL, MAE, RMSE, MAPE, MdAPE, SMAPE, SMdAPE) for L1LSA, L2LSA, L1MTE, and L2MTE methods, evaluated using model averaging from MEBoot samples (\\texttt{MEB}), the corresponding models applied on the original traffic data (\\texttt{Actual}), and post-hoc variable selectoin algorithm \\texttt{ADL1LASSO} models under AIC and BIC minimization criteria, with bold values indicating the best-performing metrics within each model for location ", strsplit(location,"_")[[1]][1], "}}", paste0("\\label{tab:metrics",strsplit(location,"_")[[1]][1],"}"),"\\end{table}")
  return(paste(latex, collapse = "\n"))
}


generate_latex_tableCorrected <- function(metrics_dashboard, location) {
  test_df <- metrics_dashboard[[location]]$TEST
  test_df <- test_df[,(colnames(test_df)!="MAPE")&(colnames(test_df)!="MdAPE")]
  
  get_row <- function(name) {
    min_te <- rbind(apply(test_df[seq_len(nrow(test_df)) %% 3 == 1, ], 2, min),apply(test_df[seq_len(nrow(test_df)) %% 3 == 2, ], 2, min),apply(test_df[seq_len(nrow(test_df)) %% 3 == 0, ], 2, min))
    test_vals <- boldMin(test_df[name, ], min_te)
    paste("  \\ \\ \\  \\ \\ \\ \\textbf{", gsub("_.*","",name), "} &", 
          paste(c(test_vals), collapse = " & "), "\\\\")
  }
  
  # Order of models under AIC/BIC
  methods <- list(
    AIC = c("MEB_L1LSA_AIC", "Actual_L1LSA_AIC", "ADL1LASSO_L1LSA_AIC",
            "MEB_L2LSA_AIC", "Actual_L2LSA_AIC", "ADL1LASSO_L2LSA_AIC",
            "MEB_L1MTE_AIC", "Actual_L1MTE_AIC", "ADL1LASSO_L1MTE_AIC",
            "MEB_L2MTE_AIC", "Actual_L2MTE_AIC", "ADL1LASSO_L2MTE_AIC"),
    BIC = c("MEB_L1LSA_BIC", "Actual_L1LSA_BIC", "ADL1LASSO_L1LSA_BIC",
            "MEB_L2LSA_BIC", "Actual_L2LSA_BIC", "ADL1LASSO_L2LSA_BIC",
            "MEB_L1MTE_BIC", "Actual_L1MTE_BIC", "ADL1LASSO_L1MTE_BIC",
            "MEB_L2MTE_BIC", "Actual_L2MTE_BIC", "ADL1LASSO_L2MTE_BIC")
  )
  
  model_names <- list(
    L1MTE = c("MEB", "Actual", "ADL1LASSO"),
    L2MTE = c("MEB", "Actual", "ADL1LASSO"),
    L1LSA = c("MEB", "Actual", "ADL1LASSO"),
    L2LSA = c("MEB", "Actual", "ADL1LASSO")
  )
  
  latex <- c(paste("\\begin{table}[H]
\\centering
\\resizebox{0.5\\textwidth}{!}{%
\\centering

\\begin{tabular}{l||l*{4}{c}}
\\toprule
\\multicolumn{2}{l}{\\multirow{3}{*}{\\textbf{\\huge{Methods}}}} & 
\\multicolumn{4}{c}{\\textbf{Location:", strsplit(location,"_")[[1]][1],"}} \\\\
\\cmidrule(lr){2-6} \\\\
&  \\textbf{HL} & \\textbf{MAE} & \\textbf{RMSE} & \\textbf{SMAPE} & \\textbf{SMdAPE} \\\\
\\midrule"))
  
  for (criterion in names(methods)) {
    latex <- c(latex, "\\noalign{\\vskip 0.1cm}")
    latex <- c(latex, paste0("\\multirow{4}{*}{\\rotatebox{90}{\\parbox{5.3cm}{\\centering \\textbf{\\large{", criterion, "}}}}}"))
    models <- unique(gsub(".*_(L[12][A-Z]+)_.*", "\\1", methods[[criterion]]))
    for (model in models) {
      rows <- methods[[criterion]][grepl(model, methods[[criterion]])]
      prefix <- ifelse(model!="L1LSA","\\ \\ \\ ","")
      latex <- c(latex, paste0(prefix,"\\multirow{4}{*}{\\rotatebox{90}{\\parbox{0.6cm}{\\centering \\textbf{", model, "}}}}"))
      for (row in rows) {
        latex <- c(latex,get_row(row))
      }
      latex <- c(latex, "\\noalign{\\vskip 0.1cm}","\\cmidrule{2-6}","\\noalign{\\vskip 0.1cm}")
    }
    latex <- c(latex, "\\midrule")
  }
  
  latex <- c(latex, "\\bottomrule", "\\end{tabular}}", "\\caption{\\textit{Comparison of five performance metrics (HL, MAE, RMSE, SMAPE, SMdAPE) for L1LSA, L2LSA, L1MTE, and L2MTE methods, evaluated using model averaging from MEBoot samples (\\texttt{MEB}), the corresponding models applied on the original traffic data (\\texttt{Actual}), and post-hoc variable selectoin algorithm \\texttt{ADL1LASSO} models under AIC and BIC minimization criteria, with bold values indicating the best-performing metrics within each model for location ", strsplit(location,"_")[[1]][1], "}}", paste0("\\label{tab:metrics",strsplit(location,"_")[[1]][1],"}"),"\\end{table}")
  return(paste(latex, collapse = "\n"))
}
copyThis <- c()
for (location in routes) copyThis <- c(copyThis, generate_latex_tableCorrected(metrics_dashboard, location))
writeClipboard(paste(copyThis, collapse= "\n \n"))

coeff_names <- c("(Intercept)",regressors)
VIP_list <- lapply(routes,function(loc) {  
  res <- data.frame(
    apply(
      cbind(
        "L1MTE_AIC" = c(T,T,"AIC"),
        "L1MTE_BIC" = c(T,T,"BIC"),
        "L1LSA_AIC" = c(T,F,"AIC"),
        "L1LSA_BIC" = c(T,F,"BIC"),
        "L2MTE_AIC" = c(F,T,"AIC"),
        "L2MTE_BIC" = c(F,T,"BIC"),
        "L2LSA_AIC" = c(F,F,"AIC"),
        "L2LSA_BIC" = c(F,F,"BIC")), 2,
      function(x) plotData_VIPs(loc,
                                c(as.logical(x[1]),as.logical(x[2])),
                                x[3])[,1]),
    "Predictor" = coeff_names)
  res <- res %>%
    arrange(across(1:8, ~ desc(.)))
  res[rowSums(res[,-9])>0,]  
  })
names(VIP_list) <- routes

get_name <- function(x) gsub(
  "\\)",
  " \\\\right)",
  gsub(
    "TIME/ 2400",
    " \\\\frac{t}{2400}",
    gsub(
      "pi",
      " \\\\pi \\\\cdot",
      gsub(
        " (\\*|\\* )",
        "",
        gsub(
          "\\(",
          "\\\\left(",
          x)))))

latexify_results <- function(location, object_res,rounding =2) {
  res <- object_res[[location]] 
  latex <- c(paste("\\begin{table}[H]
    \\centering
    \\resizebox{\\textwidth}{!}{%
    \\begin{tabular}{l||*{8}{c}}
    \\toprule
    \\multicolumn{1}{c||}{
    \\begin{tikzpicture}[baseline=(current bounding box.center)]
        \\draw[black, line width=0.3pt] (0,0.8) -- (2.3,0);
        \\node[rotate=-20, font=\\small\\bfseries] at (0.7,0.3) {Features};
        \\node[rotate=-20, font=\\small\\bfseries] at (1.6,0.45) {Methods};
    \\end{tikzpicture}
} & ",
    paste0(sapply(colnames(res[,-length(res)]), function(x) paste0("$\\mathbf{",paste0(paste0(strsplit(x,"_")[[1]],collapse="_{"),"}}$"))),collapse=" & "),"\\\\ 
    \\midrule"))
  for (r in 1:dim(res)[1]){
    latex <- c(latex, paste(c(paste0("$",get_name(res[r,length(res)]),"$"),sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(res[r,-9]))),collapse = " & "),collapse = "\\\\")
  }
  latex <- c(latex, "\\bottomrule", "\\end{tabular}}", "\\caption{\\textit{Variable Inclusion Probabilities for location ", strsplit(location,"_")[[1]][1], "}}", paste0("\\label{tab:VIPs",strsplit(location,"_")[[1]][1],"}"),"\\end{table}")
  paste(latex,collapse = "\n")
  
}

writeClipboard(paste0(sapply(routes, function(loc) latexify_results(loc,VIP_list)),collapse = " \n \n"))


COEF_avg_meb_list <- function(extract = c("Mean", "2.5q", "97.5q","Median")){
  func <- function(x, extr) {
    switch(
      extr,
      "Median" = quantile(x,0.5),
      "2.5q" = quantile(x,0.025),
      "97.5q" = quantile(x,0.975),
      mean(x)
    )
  }
  RES <- lapply(routes,function(loc) {
    res <- data.frame(
      apply(
        cbind(
          "L1MTE_AIC" = c(T,T,"AIC"),
          "L1MTE_BIC" = c(T,T,"BIC"),
          "L1LSA_AIC" = c(T,F,"AIC"),
          "L1LSA_BIC" = c(T,F,"BIC"),
          "L2MTE_AIC" = c(F,T,"AIC"),
          "L2MTE_BIC" = c(F,T,"BIC"),
          "L2LSA_AIC" = c(F,F,"AIC"),
          "L2LSA_BIC" = c(F,F,"BIC")), 2,
        function(x) apply(plotData_meb(loc,c(as.logical(x[1]),as.logical(x[2])),x[3],show_ ="Coefficients")[,-21],1,function(z) func(z,extract))),
      "Predictor" = coeff_names)
    res <- res %>%
      arrange(across(1:(length(res)-1), ~ desc(.)))
    res
  })
  names(RES) <- routes
  RES
}
# # Mean
# writeClipboard(paste0(sapply(routes, function(loc) latexify_results(loc,COEF_avg_meb_list("Mean"))),collapse = " \n \n"))
# # 2.5q
# writeClipboard(paste0(sapply(routes, function(loc) latexify_results(loc,COEF_avg_meb_list("2.5q"),7)),collapse = " \n \n"))
# # Median
# writeClipboard(paste0(sapply(routes, function(loc) latexify_results(loc,COEF_avg_meb_list("97.5q"))),collapse = " \n \n"))
# # 97.5q
# writeClipboard(paste0(sapply(routes, function(loc) latexify_results(loc,COEF_avg_meb_list("Median"))),collapse = " \n \n"))

latexify_tiny_nums <- function(pmean_vec,p2.5_vec,pmed_vec,p97.5_vec){sapply(1:length(pmean_vec),  function(id) paste0("$\\stackrel{",pmean_vec[id],"}{\\text{ \\tiny (\\textcolor{coral}{",p2.5_vec[id],"}, \\textcolor{navyblue}{",pmed_vec[id],"} , \\textcolor{coral}{",p97.5_vec[id],"})}}$"))}

latexify_mean_N_CI <- function(location, rounding=2){
  pmean <- COEF_avg_meb_list("Mean")[[location]]
  pmed <-  COEF_avg_meb_list("Median")[[location]]
  p2.5 <-  COEF_avg_meb_list("2.5q")[[location]]
  p97.5 <- COEF_avg_meb_list("97.5q")[[location]]
  
  pmean <- pmean[rowSums(abs(pmean[,-length(pmean)]))>0,]
  pmed  <- pmed[pmed$Predictor %in% pmean$Predictor,]
  p2.5  <- p2.5[p2.5$Predictor %in% pmean$Predictor,]
  p97.5 <- p97.5[p97.5$Predictor %in% pmean$Predictor,]
  columnNames <- colnames(pmean)
  numCols <-  length(columnNames)-1

  latex <- c(paste("\\begin{table}[H]
  \\centering
  \\resizebox{\\textwidth}{!}{%
  \\begin{tabular}{l||*{8}{c}}
  \\toprule
  \\textbf{Frequencies} & ",
                   paste0(sapply(colnames(pmean[,-length(pmean)]), function(x) paste0("$\\mathbf{",paste0(paste0(strsplit(x,"_")[[1]],collapse="_{"),"}}$"))),collapse=" & "),"\\\\ 
  \\midrule"))
  lastCol <- dim(pmean)[2]
  for (r in 1:dim(pmean)[1]){
    round_pmean  <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(pmean[r,-lastCol]))
    round_pmed   <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(pmed[r,-lastCol]))
    round_p2.5   <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(p2.5[r,-lastCol]))
    round_p97.5  <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(p97.5[r,-lastCol]))
    actual       <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(p97.5[r,-lastCol]))
    adL1Lasso    <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(p97.5[r,-lastCol]))
    latex <- c(latex, paste(c(paste0("$",get_name(pmean[r,lastCol]),"$"),latexify_tiny_nums(round_pmean,round_p2.5,round_pmed,round_p97.5)),collapse = " & "))
  }
  latex <- c(paste(latex,collapse = "\\\\ \n"), "\\\\ \n \\bottomrule", "\\end{tabular}}", "\\caption{Location  ", strsplit(location,"_")[[1]][1],
             paste0(" location. The \\textcolor{coral}{coral} tiny numbers indicate the 2.5th and 97.5th quantiles,",
                    " the \\textcolor{navyblue}{navy-blue} tiny numbers represent the median values ",
                    "and normal-sized numbers refer to the mean value of the estimated coefficients model-wise}"),
             paste0("\\label{tab:Coef",strsplit(location,"_")[[1]][1],"}"),
             "\\end{table}")
  paste(latex,collapse = "\n")
}

writeClipboard(paste0(sapply(routes, function(loc) {print(loc); latexify_mean_N_CI(loc)}),collapse = " \n \n"))

outerJoin <- function(.list){
  merged <- do.call(rbind,.list)
  concl <- as.data.frame(do.call(rbind,lapply(unique(merged$Predictor),FUN=function(pred) colMeans(merged[merged$Predictor==pred,-length(merged)]))))
  concl$Overall <- rowMeans(concl[,-length(merged)])
  concl$Predictors <- unique(merged$Predictor)
  return(concl)
}

concl <- outerJoin(VIP_list)
concl <- concl %>%
  arrange(across(1:(length(concl)-1), ~ desc(.)))
concl <- concl[rowSums(abs(concl[,-length(concl)]))>0,]  
concl$Overall <- rowMeans(concl[,-length(concl)])
latex <- c(paste("\\begin{table}[H]
    \\centering
    \\resizebox{\\textwidth}{!}{%
    \\begin{tabular}{l||*{9}{c}}
    \\toprule
    \\multicolumn{1}{c||}{
    \\begin{tikzpicture}[baseline=(current bounding box.center)]
        \\draw[black, line width=0.3pt] (0,0.8) -- (2.3,0);
        \\node[rotate=-20, font=\\small\\bfseries] at (0.7,0.3) {Features};
        \\node[rotate=-20, font=\\small\\bfseries] at (1.6,0.45) {Methods};
    \\end{tikzpicture}
} & ",
                 paste0(sapply(colnames(concl[,-((length(concl)-1):length(concl))]), function(x) paste0("$\\mathbf{",paste0(paste0(strsplit(x,"_")[[1]],collapse="_{"),"}}$"))),collapse=" & ")," & $\\mathbf{",colnames(concl)[(length(concl)-1)],"} \\\\ 
    \\midrule"))
for (r in 1:dim(concl)[1]){
  latex <- c(latex, paste(c(paste0("$",get_name(concl[r,length(concl)]),"$"),sprintf(paste0("%.",as.character(2),"f"), as.numeric(concl[r,-length(concl)]))),collapse = " & "),collapse = "\\\\")
}
latex <- c(latex, "\\bottomrule", "\\end{tabular}}", "\\caption{\\textit{Variable Inclusion Probabilities for according to all locations}}", paste0("\\label{tab:VIPs_total}"),"\\end{table}")
latexVIPsTotal <- paste(latex,collapse = "\n")
writeClipboard(latexVIPsTotal)


plot_ly(
  x = concl[,9],
  y = colnames(concl[,-9]),
  z = t(concl),  # transpose to flip axes
  type = "heatmap",
  colors = colorRamp(c("navyblue", "white", "coral"))
)



pmean <- outerJoin(COEF_avg_meb_list("Mean"))
pmed <-  outerJoin(COEF_avg_meb_list("Median"))
p2.5 <-  outerJoin(COEF_avg_meb_list("2.5q"))
p97.5 <- outerJoin(COEF_avg_meb_list("97.5q"))

columnNames <- colnames(pmean)
numCols <-  length(columnNames)-1
dfs <- list(pmean, p2.5, pmed, p97.5)
joined_dt <- reduce(dfs, dplyr::full_join, by = "Predictor")
joined_dt[is.na(joined_dt)] <- 0
preds <- joined_dt$Predictor
joined_dt <- joined_dt[,!(names(joined_dt) %in% "Predictor")]

pmean <- cbind(joined_dt[,1:numCols], preds)
pmed  <- cbind(joined_dt[,(numCols+1):(2*numCols)], preds)
p2.5  <- cbind(joined_dt[,(2*numCols+1):(3*numCols)], preds)
p97.5 <- cbind(joined_dt[,(3*numCols+1):(4*numCols)], preds)
colnames(pmean) <- colnames(pmed) <- colnames(p2.5) <- colnames(p97.5) <- columnNames

latex <- c(paste("\\begin{table}[H]
  \\centering
  \\resizebox{\\textwidth}{!}{%
  \\begin{tabular}{l||*{8}{c}}
  \\toprule
  \\textbf{Frequencies} & ",
                 paste0(sapply(colnames(pmean[,-length(pmean)]), function(x) paste0("$\\mathbf{",paste0(paste0(strsplit(x,"_")[[1]],collapse="_{"),"}}$"))),collapse=" & "),"\\\\ 
  \\midrule"))
lastCol <- dim(pmean)[2]
for (r in 1:dim(pmean)[1]){
  round_pmean  <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(pmean[r,-lastCol]))
  round_pmed   <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(pmed[r,-lastCol]))
  round_p2.5   <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(p2.5[r,-lastCol]))
  round_p97.5  <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(p97.5[r,-lastCol]))
  actual       <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(p97.5[r,-lastCol]))
  adL1Lasso    <- sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(p97.5[r,-lastCol]))
  latex <- c(latex, paste(c(paste0("$",get_name(pmean[r,lastCol]),"$"),latexify_tiny_nums(round_pmean,round_p2.5,round_pmed,round_p97.5)),collapse = " & "))
}
latex <- c(paste(latex,collapse = "\\\\ \n"), "\\\\ \n \\bottomrule", "\\end{tabular}}", "\\caption{Location  ", strsplit(location,"_")[[1]][1],
           paste0(" location. The \\textcolor{coral}{coral} tiny numbers indicate the 2.5th and 97.5th quantiles,",
                  " the \\textcolor{navyblue}{navy-blue} tiny numbers represent the median values ",
                  "and normal-sized numbers refer to the mean value of the estimated coefficients model-wise}"),
           paste0("\\label{tab:Coef",strsplit(location,"_")[[1]][1],"}"),
           "\\end{table}")
paste(latex,collapse = "\n")



METRMERGEDtr <- rbind(metrics_dashboard[[1]]$TRAIN,
                      metrics_dashboard[[2]]$TRAIN,
                      metrics_dashboard[[3]]$TRAIN,
                      metrics_dashboard[[4]]$TRAIN,
                      metrics_dashboard[[5]]$TRAIN,
                      metrics_dashboard[[6]]$TRAIN,
                      metrics_dashboard[[7]]$TRAIN)

min_METRMERGEDtr <- apply(METRMERGEDtr,2,min)
METRMERGEDte <- rbind(metrics_dashboard[[1]]$TEST,
                      metrics_dashboard[[2]]$TEST,
                      metrics_dashboard[[3]]$TEST,
                      metrics_dashboard[[4]]$TEST,
                      metrics_dashboard[[5]]$TEST,
                      metrics_dashboard[[6]]$TEST,
                      metrics_dashboard[[7]]$TEST)
min_METRMERGEDte <- apply(METRMERGEDte,2,min)
meanIncrease <-  matrix(NA, nrow= 8, ncol = 7 )
k <- 0 
for (i in c("L1", "L2")) for (j in c("AIC", "BIC")) {
  k <- k +1 
  cat(i, j, ":\n")
  A <- METRMERGEDtr[Filter(function(x) grepl(paste0("^MEB_",i,"MTE_",j),x),  rownames(METRMERGEDtr)), ]
  B <- METRMERGEDtr[Filter(function(x) grepl(paste0("^MEB_",i,"LSA_",j),x),  rownames(METRMERGEDtr)), ]
  print((B-A)/B*100)
  meanIncrease[k,] <- colMeans((B-A)/B*100)
}
for (i in c("L1", "L2")) for (j in c("AIC", "BIC")) {
  k <- k +1 
  cat(i, j, ":\n")
  A <- METRMERGEDte[Filter(function(x) grepl(paste0("^MEB_",i,"MTE_",j),x),  rownames(METRMERGEDte)), ]
  B <- METRMERGEDte[Filter(function(x) grepl(paste0("^MEB_",i,"LSA_",j),x),  rownames(METRMERGEDte)), ]
  print((B-A)/B*100)
  meanIncrease[k,] <- colMeans((B-A)/B*100)
}
apply(meanIncrease[1:4,],2,min)
apply(meanIncrease[1:4,],2,max)
colMeans(meanIncrease[1:4,])

apply(meanIncrease[5:8,],2,min)
apply(meanIncrease[5:8,],2,max)
colMeans(meanIncrease[5:8,])


METRMERGEDtr <- METRMERGEDtr %>%
  as.data.frame() %>%
  arrange(.[[1]], .[[2]], .[[3]])

METRMERGEDte <- METRMERGEDte %>%
  as.data.frame() %>%
  arrange(.[[1]], .[[2]], .[[3]])

min_rows_tr <- apply(METRMERGEDtr, 2, function(col) {
  rownames(METRMERGEDtr)[which.min(col)]
})

# Row names for test metrics
min_rows_te <- apply(METRMERGEDte, 2, function(col) {
  rownames(METRMERGEDte)[which.min(col)]
})




latexify_results_locationwise <- function(methpd, object_res,rounding =2) {
  res <- object_res[[location]] 
  latex <- c(paste("\\begin{table}[H]
    \\centering
    \\resizebox{\\textwidth}{!}{%
    \\begin{tabular}{l||*{8}{c}}
    \\toprule
    \\multicolumn{1}{c||}{
    \\begin{tikzpicture}[baseline=(current bounding box.center)]
        \\draw[black, line width=0.3pt] (0,0.8) -- (2.3,0);
        \\node[rotate=-20, font=\\small\\bfseries] at (0.7,0.3) {Features};
        \\node[rotate=-20, font=\\small\\bfseries] at (1.6,0.45) {Methods};
    \\end{tikzpicture}
} & ",
                   paste0(sapply(colnames(res[,-length(res)]), function(x) paste0("$\\mathbf{",paste0(paste0(strsplit(x,"_")[[1]],collapse="_{"),"}}$"))),collapse=" & "),"\\\\ 
    \\midrule"))
  for (r in 1:dim(res)[1]){
    latex <- c(latex, paste(c(paste0("$",get_name(res[r,length(res)]),"$"),sprintf(paste0("%.",as.character(rounding),"f"), as.numeric(res[r,-9]))),collapse = " & "),collapse = "\\\\")
  }
  latex <- c(latex, "\\bottomrule", "\\end{tabular}}", "\\caption{\\textit{Variable Inclusion Probabilities for location ", strsplit(location,"_")[[1]][1], "}}", paste0("\\label{tab:VIPs",strsplit(location,"_")[[1]][1],"}"),"\\end{table}")
  paste(latex,collapse = "\n")
  
}

refine_dataframe <- function(metrics_dashboard, method){
  refine_df <- list()
  ref_test_df <- ref_train_df <- matrix(NA, ncol=ncol(metrics_dashboard[[1]]$TRAIN))
  colnames(ref_test_df) <- colnames(ref_train_df) <- colnames(metrics_dashboard[[1]]$TRAIN)
  for (location in routes){
    ref_train_df <- rbind(ref_train_df, metrics_dashboard[[location]]$TRAIN[grepl(paste0(method,"_AIC"),rownames(metrics_dashboard[[location]]$TRAIN)),])
    ref_test_df <- rbind(ref_test_df, metrics_dashboard[[location]]$TEST[grepl(paste0(method,"_AIC"),rownames(metrics_dashboard[[location]]$TEST)),])
  }
  ref_train_df <- ref_train_df[-1,]; ref_test_df <- ref_test_df[-1,]; 
  for (location in routes){
    ref_train_df <- rbind(ref_train_df, metrics_dashboard[[location]]$TRAIN[grepl(paste0(method,"_BIC"),rownames(metrics_dashboard[[location]]$TRAIN)),])
    ref_test_df <- rbind(ref_test_df, metrics_dashboard[[location]]$TEST[grepl(paste0(method,"_BIC"),rownames(metrics_dashboard[[location]]$TEST)),])
  }
  rownames(ref_test_df) <- rownames(ref_train_df) <- c(as.vector(sapply(routes, function(x) c(paste0(strsplit(x,"_")[[1]][1],"_MEB_AIC"),paste0(strsplit(x,"_")[[1]][1],"_Actual_AIC"),paste0(strsplit(x,"_")[[1]][1],"_ADL1LASSO_AIC")))),as.vector(sapply(routes, function(x) c(paste0(strsplit(x,"_")[[1]][1],"_MEB_BIC"),paste0(strsplit(x,"_")[[1]][1],"_Actual_BIC"),paste0(strsplit(x,"_")[[1]][1],"_ADL1LASSO_BIC")))))
  refine_df$TRAIN <- ref_train_df; refine_df$TEST <- ref_test_df
  refine_df
  
}

generate_latex_table_locationwise <- function(metrics_dashboard, method) {
  boldMin <- function(df,mins){
    edited <- matrix(NA, nrow=dim(df)[1],ncol=dim(df)[2])
    for (column in seq_len(ncol(df))){
      edited[,column] <- sapply(df[,column], FUN = function(l) ifelse(!(l %in% mins[,column]), sprintf("%.2f", as.numeric(l)), paste0("\\textbf{",sprintf("%.2f", as.numeric(l)),"}")))
    }
    edited
  }
  refine_df <- refine_dataframe(metrics_dashboard, method= "L1MTE")
  train_df <- refine_df$TRAIN
  test_df <- refine_df$TEST
  get_row <- function(name) {
    min_tr <- rbind(apply(train_df[seq_len(nrow(train_df)) %% 3 == 1, ], 2, min), apply(train_df[seq_len(nrow(train_df)) %% 3 == 2, ], 2, min), apply(train_df[seq_len(nrow(train_df)) %% 3 == 0, ], 2, min))
    min_te <- rbind(apply(test_df[seq_len(nrow(test_df)) %% 3 == 1, ], 2, min),apply(test_df[seq_len(nrow(test_df)) %% 3 == 2, ], 2, min),apply(test_df[seq_len(nrow(test_df)) %% 3 == 0, ], 2, min))
    dt_tr <- train_df[grepl(name,rownames(train_df)), ]
    dt_te <- test_df[grepl(name,rownames(test_df)), ]
    final_names <- sapply(rownames(dt_tr),function(x){
      parts <- strsplit(x, "_")[[1]]
      paste0("\\textbf{",parts[2],"}")
      #paste0("$\\mathbf{", parts[1], "_{", parts[2], "_{", parts[3], "}}}$")
    } )
    train_vals <- cbind(final_names,boldMin(dt_tr, min_tr))
    test_vals <- boldMin(dt_te, min_te)
    parts <- strsplit(name, "_")[[1]]
    paste("  \\ \\ \\  \\ \\ \\ \\ ", 
          paste(c(train_vals, test_vals), collapse = " & "), "\\\\")
  }
  
  # Order of models under AIC/BIC
  methods <- list(
    AIC = rownames(train_df)[grepl("AIC",rownames(train_df))], 
    BIC = rownames(train_df)[grepl("BIC",rownames(train_df))]
  )

  # \\multicolumn{2}{l}{\\multirow{3}{*}{\\textbf{\\large{Locationwise Methods}}}}
  latex <- c(paste("\\begin{table}[H]
\\centering
\\resizebox{\\textwidth}{!}{%
\\begin{tabular}{l||l*{6}{c}|*{7}{c}}
\\toprule
\\multicolumn{2}{l}{\\multirow{3}{*}{\\begin{tikzpicture}[baseline=(current bounding box.center)]
        \\draw[black, line width=0.3pt] (0,0.8) -- (3.4,-0.5);
        \\node[rotate=-20, font=\\large\\bfseries] at (1.4,-0.05) {Methods};
        \\node[rotate=-20, font=\\large\\bfseries] at (1.6,0.55) {Metrics};
    \\end{tikzpicture}}}
 & 
\\multicolumn{13}{c}{\\textbf{Method:} $\\mathbf{", paste0("\\ell_",substr(method, 2, nchar(method))), "}$} \\\\
\\cmidrule(lr){2-15}
& \\multicolumn{7}{c}{\\textbf{TRAIN}} & \\multicolumn{7}{c}{\\textbf{TEST}} \\\\
\\cmidrule(lr){2-8} \\cmidrule(lr){9-15}
&  \\textbf{HL} & \\textbf{MAE} & \\textbf{RMSE} & \\textbf{MAPE} & \\textbf{MdAPE} & \\textbf{SMAPE} & \\textbf{SMdAPE} &  \\textbf{HL} & \\textbf{MAE} & \\textbf{RMSE} & \\textbf{MAPE} & \\textbf{MdAPE} & \\textbf{SMAPE} & \\textbf{SMdAPE} \\\\
\\midrule"))
  
  for (criterion in names(methods)) {
    latex <- c(latex, "\\noalign{\\vskip 0.1cm}")
    latex <- c(latex, paste0("\\multirow{25}{*}{\\rotatebox{90}{\\parbox{5.3cm}{\\centering \\textbf{\\large{", criterion, "}}}}}"))
    # models <- unique(gsub(".*_(L[12][A-Z]+)_.*", "\\1", methods[[criterion]]))
    for (route in routes) {
      rows <- methods[[criterion]][grepl(strsplit(route,"_")[[1]][1], methods[[criterion]])]
      prefix <- ifelse(route!="L101_volume","\\ \\ \\ ","")
      latex <- c(latex, paste0(prefix,"\\multirow{4}{*}{\\rotatebox{90}{\\parbox{0.6cm}{\\centering \\textbf{", strsplit(route,"_")[[1]][1], "}}}}"))
      for (row in rows) {
        latex <- c(latex,get_row(row))
      }
      latex <- c(latex, "\\noalign{\\vskip 0.1cm}","\\cmidrule{2-15}","\\noalign{\\vskip 0.1cm}")
    }
    latex <- c(latex, "\\midrule ")
  }
  
  latex <- c(latex, "\\bottomrule", "\\end{tabular}}", "\\caption{\\textit{Comparison of seven performance metrics (HL, MAE, RMSE, MAPE, MdAPE, SMAPE, SMdAPE) for $", paste0("\\ell_",substr(method, 2, nchar(method))), "$ method, evaluated using model averaging from MEBoot samples (\\texttt{MEB}), the corresponding models applied on the original traffic data (\\texttt{Actual}), and post-hoc variable selection algorithm \\texttt{ADL1LASSO} models under AIC and BIC minimization criteria, with bold values indicating the best-performing metrics within each model.}}", paste0("\\label{tab:metrics",method,"_locationwise}"),"\\end{table}")
  return(paste(latex, collapse = "\n"))
}

writeClipboard(paste0(sapply(c("L1MTE","L2MTE","L1LSA","L2LSA"), function(x) generate_latex_table_locationwise(metrics_dashboard,x)),"\n \n "))







refine_dataframe2 <- function(metrics_dashboard, method){
  refine_df <- list() 
  metrics <- colnames(metrics_dashboard$L101_volume$TRAIN)
  for (metric in metrics){
    train_dataframe <- test_dataframe <- matrix(NA, ncol = 24,nrow = 7)
    colnames(train_dataframe) <- colnames(test_dataframe) <- rownames(metrics_dashboard$L101_volume$TRAIN)
    for (method in colnames(test_dataframe)){
      train_dataframe[,method] <- as.numeric(sapply(routes, FUN = function(route) metrics_dashboard[[route]]$TRAIN[method,metric]))
      test_dataframe[,method] <- as.numeric(sapply(routes, FUN = function(route) metrics_dashboard[[route]]$TEST[method,metric]))
    }
    rownames(train_dataframe) <- rownames(test_dataframe) <- sapply(routes,function(x) strsplit(x,"_")[[1]][1])
    refine_df[[metric]]$TRAIN <- train_dataframe
    refine_df[[metric]]$TEST <- test_dataframe
  }
  refine_df
}

generate_latex_table_metricwise <- function(metrics_dashboard, metric) {
  boldMin <- function(df,mins){
    edited <- matrix(NA, nrow=dim(df)[1],ncol=dim(df)[2])
    for (row in seq_len(nrow(df))){
      edited[row,] <- sapply(df[row,], FUN = function(l) ifelse(!(l == min_tr[row]), sprintf("%.2f", as.numeric(l)), paste0("\\textbf{",sprintf("%.2f", as.numeric(l)),"}")))
    }
    colnames(edited) <- colnames(df)
    edited
  }
  refine_df <- refine_dataframe2(metrics_dashboard, method)[[metric]]
  train_df <- refine_df$TRAIN
  test_df <- refine_df$TEST
  get_row <- function() {
    min_tr <- apply(train_df, 1, min)
    min_te <- apply(test_df, 1, min)
    dt_tr <- train_df
    dt_te <- test_df
    final_names <- sapply(rownames(dt_tr),function(x){ paste0("\\textbf{",x,"}") })
    train_vals <- cbind(final_names,boldMin(dt_tr, min_tr))
    test_vals <- cbind(final_names,boldMin(dt_te, min_te))
    res <-  paste0("\\multirow{2}{*}{\\rotatebox{90}{\\parbox{3cm}{\\centering \\textbf{\\large{TRAIN SET}}}}} \n \\ \\ \\ ",paste0(apply(train_vals,1, function(x) paste(x,collapse=" & ")),collapse="  \\\\ \n \\ \\ \\ \\ \\ \\ \\  "),"\\\\ \n \\midrule")
    paste0(res,"\\multirow{2}{*}{\\rotatebox{90}{\\parbox{3cm}{\\centering \\textbf{\\large{TEST SET}}}}} \n \\ \\ \\ ",paste0(apply(test_vals,1, function(x) paste(x,collapse=" & ")),collapse="  \\\\ \n \\ \\ \\ \\ \\ \\ \\  "))
  }# cbind(train_vals, test_vals)
  
  # Order of models under AIC/BIC
  metrics <- list("HL" = "Huber Loss",
                  "MAE" = "Mean Absolute Error (MAE)",
                  "RMSE" = "Root Mean Square Error (RMSE)",
                  "MAPE" = "Mean Absolute Percentage Error (MAPE)",
                  "MdAPE" = "Median Absolute Percentage Error (MdAPE)",
                  "SMAPE" = "Symmetric Mean Absolute Percentage Error (SMAPE)",
                  "SMdAPE"= "Symmetric Median Absolute Percentage Error (SMdAPE)")
  
  # \\multicolumn{2}{l}{\\multirow{3}{*}{\\textbf{\\large{Locationwise Methods}}}}
  latex <- c(paste("\\begin{table}[H]
\\centering
\\resizebox{\\textwidth}{!}{%
\\begin{tabular}{l||*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}}
\\toprule
\\multicolumn{2}{l}{\\multirow{3}{*}{\\begin{tikzpicture}[baseline=(current bounding box.center)]
        \\draw[black, line width=0.3pt] (-0.8,1.) -- (0.5,-0.7);
        \\node[rotate=-53, font=\\large\\bfseries] at (-0.3,0) {Locations};
        \\node[rotate=-53, font=\\large\\bfseries] at (0,0.55) {Methods};
    \\end{tikzpicture}}}
 & 
 \\multicolumn{23}{c}{\\hspace{0.5cm}\\textbf{Metric:", metrics[[metric]],"}} \\\\
\\cmidrule(lr){2-25}\\\\
& \\multicolumn{3}{c}{$\\mathbf{\\ell_1LSA_{AIC}}$} \n & \\multicolumn{3}{c}{$\\mathbf{\\ell_2LSA_{AIC}}$} \n &  \\multicolumn{3}{c}{$\\mathbf{\\ell_1MTE_{AIC}}$} \n &  \\multicolumn{3}{c}{$\\mathbf{\\ell_2MTE_{AIC}}$} \n &   \\multicolumn{3}{c}{$\\mathbf{\\ell_1LSA_{BIC}}$} \n &   \\multicolumn{3}{c}{$\\mathbf{\\ell_2LSA_{BIC}}$} \n & \\multicolumn{3}{c}{$\\mathbf{\\ell_1MTE_{BIC}}$} \n & \\multicolumn{3}{c}{$\\mathbf{\\ell_2MTE_{BIC}}$} \\\\ 
\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}  \\cmidrule(lr){11-13} \\cmidrule(lr){14-16} \\cmidrule(lr){17-19} \\cmidrule(lr){20-22} \\cmidrule(lr){23-25} \\\\
& \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO} \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO} \n  & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \\\\
\\midrule
"))
  latex<- c(latex,get_row())
  latex <- c(latex, "\\\\ \n \\bottomrule", "\\end{tabular}}", "\\caption{\\textit{Comparison of all the methods ($\\ell_1MTE$,$\\ell_2MTE$, $\\ell_1LSA$, $\\ell_2LSA$ under AIC and BIC minimization) under",metrics[[metric]],"metric , evaluated using model averaging from MEBoot samples (\\texttt{MEB}), the corresponding models applied on the original traffic data (\\texttt{Actual}), and post-hoc variable selection algorithm \\texttt{ADL1LASSO} models under AIC and BIC minimization criteria, with bold values indicating the best-performing metrics within each model.}}", paste0("\\label{tab:methods",metric,"_metricwise}"),"\\end{table}")
  return(paste(latex, collapse = "\n"))
}



generate_latex_table_metricwiseCorrected <- function(metrics_dashboard, metric) {
  
  refine_df <- refine_dataframe2(metrics_dashboard, method)[[metric]]
  test_df <- refine_df$TEST
  get_row2 <- function() {
    min_te <- apply(test_df, 1, min)
    dt_te <- test_df
    final_names <- sapply(rownames(dt_te),function(x){ paste0("\\textbf{",x,"}") })
    boldMin <- function(df,mins){
      edited <- matrix(NA, nrow=dim(df)[1],ncol=dim(df)[2])
      for (row in seq_len(nrow(df))){
        edited[row,] <- sapply(df[row,], FUN = function(l) ifelse(!(l == min_te[row]), sprintf("%.2f", as.numeric(l)), paste0("\\textbf{",sprintf("%.2f", as.numeric(l)),"}")))
      }
      colnames(edited) <- colnames(df)
      edited
    }
    test_vals <- cbind(final_names,boldMin(dt_te, min_te))
    paste0("\\ \\ \\ \\ \\ \\ \\ ",paste0(apply(test_vals,1, function(x) paste(x,collapse=" & ")),collapse="  \\\\ \n \\ \\ \\ \\ \\ \\ \\  "))
  }# cbind(train_vals, test_vals)
  
  # Order of models under AIC/BIC
  metrics <- list("HL" = "Huber Loss",
                  "MAE" = "Mean Absolute Error (MAE)",
                  "RMSE" = "Root Mean Square Error (RMSE)",
                  "SMAPE" = "Symmetric Mean Absolute Percentage Error (SMAPE)",
                  "SMdAPE"= "Symmetric Median Absolute Percentage Error (SMdAPE)")
  
  # \\multicolumn{2}{l}{\\multirow{3}{*}{\\textbf{\\large{Locationwise Methods}}}}
  latex <- c(paste("\\begin{table}[H]
\\centering
\\resizebox{\\textwidth}{!}{%
\\begin{tabular}{l||*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}|*{3}{c}}
\\toprule
\\multicolumn{2}{l}{\\multirow{3}{*}{\\begin{tikzpicture}[baseline=(current bounding box.center)]
        \\draw[black, line width=0.3pt] (-0.8,1.) -- (0.5,-0.7);
        \\node[rotate=-53, font=\\large\\bfseries] at (-0.3,0) {Locations};
        \\node[rotate=-53, font=\\large\\bfseries] at (0,0.55) {Methods};
    \\end{tikzpicture}}}
 & 
 \\multicolumn{23}{c}{\\hspace{0.5cm}\\textbf{Metric:", metrics[[metric]],"}} \\\\
\\cmidrule(lr){2-25}\\\\
& \\multicolumn{3}{c}{$\\mathbf{\\ell_1LSA_{AIC}}$} \n & \\multicolumn{3}{c}{$\\mathbf{\\ell_2LSA_{AIC}}$} \n &  \\multicolumn{3}{c}{$\\mathbf{\\ell_1MTE_{AIC}}$} \n &  \\multicolumn{3}{c}{$\\mathbf{\\ell_2MTE_{AIC}}$} \n &   \\multicolumn{3}{c}{$\\mathbf{\\ell_1LSA_{BIC}}$} \n &   \\multicolumn{3}{c}{$\\mathbf{\\ell_2LSA_{BIC}}$} \n & \\multicolumn{3}{c}{$\\mathbf{\\ell_1MTE_{BIC}}$} \n & \\multicolumn{3}{c}{$\\mathbf{\\ell_2MTE_{BIC}}$} \\\\ 
\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}  \\cmidrule(lr){11-13} \\cmidrule(lr){14-16} \\cmidrule(lr){17-19} \\cmidrule(lr){20-22} \\cmidrule(lr){23-25} \\\\
& \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO} \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \n & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO} \n  & \\textbf{MEB} & \\textbf{Actual} & \\textbf{ADL1LASSO}  \\\\
\\midrule
"))
  latex<- c(latex,get_row2())
  latex <- c(latex, "\\\\ \n \\bottomrule", "\\end{tabular}}", "\\caption{\\textit{Comparison of all the methods ($\\ell_1MTE$,$\\ell_2MTE$, $\\ell_1LSA$, $\\ell_2LSA$ under AIC and BIC minimization) under",metrics[[metric]],"metric , evaluated using model averaging from MEBoot samples (\\texttt{MEB}), the corresponding models applied on the original traffic data (\\texttt{Actual}), and post-hoc variable selection algorithm \\texttt{ADL1LASSO} models under AIC and BIC minimization criteria, with bold values indicating the best-performing metrics within each model.}}", paste0("\\label{tab:methods_",metric,"_metricwise}"),"\\end{table}")
  return(paste(latex, collapse = "\n"))
}
metrics_ <- c("HL","MAE","RMSE","SMAPE","SMdAPE")
writeClipboard(paste0(sapply(metrics_, function(x) generate_latex_table_metricwiseCorrected(metrics_dashboard,x)),"\n \n "))

