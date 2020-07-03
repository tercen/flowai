library(tercen)
library(tidyverse)
library(flowCore)
library(FlowSOM)
library(flowAI)

matrix2flowset <- function(a_matrix){ 
  
  minRange<- matrixStats::colMins(a_matrix)
  maxRange<- matrixStats::colMaxs(a_matrix)
  range<- maxRange-minRange
  
  df_params <- data.frame(name=colnames(a_matrix), desc=colnames(a_matrix), range=range, minRange=minRange, maxRange=maxRange)
  params <- Biobase::AnnotatedDataFrame()
  Biobase::pData(params) <- df_params
  Biobase::varMetadata(params) <- data.frame(labelDescription=c("Name of Parameter", "Description of Parameter","Range of Parameter","Minimum Parameter Value after Transformation","Maximum Parameter Value after Transformation"))
  keyval <-list()
  keyval$FILENAME <-"data"
  keyval$`$PAR` <- ncol(a_matrix)
  keyval$`$P1R` <- range[1]
  flowset <- flowCore::flowFrame(a_matrix, params, keyval)
}

ctx = tercenCtx()
data = ctx$as.matrix() 

stopifnot("Time" %in% ctx$cnames)
time <- ctx$cselect("Time")

data = t(data)
data <- as.matrix(cbind(data,time))

fc_frame <- matrix2flowset(data)

qc_frame <- suppressWarnings(flowAI::flow_auto_qc(fc_frame, output = 2, html_report = FALSE, mini_report = FALSE, fcs_QC = FALSE, folder_results = FALSE))

qc_df <- as.data.frame(exprs(qc_frame))
flag <- ifelse(qc_df[["QCvector"]] >= 10000, "fail", "pass")
qc_df <- cbind(qc_df["QCvector"], flag, .ci = (1:nrow(qc_df)))

result <- ctx$addNamespace(qc_df)
ctx$save(result)