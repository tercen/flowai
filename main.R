library(tercen)
library(dplyr)
library(tidyverse)
library(flowCore)
library(FlowSOM)
library(flowAI)

matrix2flowset <- function(a_matrix){ 
  
  minRange <- matrixStats::colMins(a_matrix)
  maxRange <- matrixStats::colMaxs(a_matrix)
  rnge <- maxRange - minRange
  
  df_params <- data.frame(
    name = colnames(a_matrix),
    desc = colnames(a_matrix),
    range = rnge,
    minRange = minRange,
    maxRange = maxRange
  )
  
  params <- Biobase::AnnotatedDataFrame()
  Biobase::pData(params) <- df_params
  Biobase::varMetadata(params) <- data.frame(
    labelDescription = c("Name of Parameter",
                         "Description of Parameter",
                         "Range of Parameter",
                         "Minimum Parameter Value after Transformation",
                         "Maximum Parameter Value after Transformation")
  )
  
  keyval <- list()
  keyval$FILENAME <-"data"
  keyval$`$PAR` <- ncol(a_matrix)
  sq <- seq_len(ncol(a_matrix))
  eval(parse(text = paste0("keyval$`$P", sq, "R` <- rnge[", sq, "];")))
  flowset <- flowCore::flowFrame(a_matrix, params, keyval)
  
  return(flowset)
}

ctx <- tercenCtx()

input.pars <- list(
  second_fractionFR = ifelse(is.null(ctx$op.value('second_fractionFR')), 0.1, ctx$op.value('second_fractionFR')),
  alphaFR = ifelse(is.null(ctx$op.value('alphaFR')), 0.01, ctx$op.value('alphaFR')),
  decompFR = ifelse(is.null(ctx$op.value('decompFR')), TRUE, ctx$op.value('decompFR')),
  outlier_binsFS = ifelse(is.null(ctx$op.value('outlier_binsFS')), FALSE, ctx$op.value('outlier_binsFS')),
  pen_valueFS = ifelse(is.null(ctx$op.value('pen_valueFS')), 500, ctx$op.value('pen_valueFS')),
  max_cptFS = ifelse(is.null(ctx$op.value('max_cptFS')), 3, ctx$op.value('max_cptFS')),
  sideFM = ifelse(is.null(ctx$op.value('sideFM')), "both", ctx$op.value('sideFM')),
  neg_valuesFM = ifelse(is.null(ctx$op.value('neg_valuesFM')), 1, as.double(ctx$op.value('neg_valuesFM')))
)

time <- ctx$cselect(ctx$cnames[[1]])

data <- ctx$as.matrix() %>% t()
data <- as.matrix(cbind(data, time))

fc_frame <- matrix2flowset(data)

qc_frame <- suppressWarnings(flowAI::flow_auto_qc(
  fcsfiles = fc_frame,
  output = 2,
  timeCh = NULL,
  second_fractionFR = input.pars$second_fractionFR,
  alphaFR = input.pars$alphaFR,
  decompFR = input.pars$decompFR,
  outlier_binsFS = input.pars$outlier_binsFS, 
  pen_valueFS = input.pars$pen_valueFS,
  max_cptFS = input.pars$max_cptFS,
  sideFM = input.pars$sideFM,
  neg_valuesFM = input.pars$neg_valuesFM,
  html_report = FALSE,
  mini_report = FALSE,
  fcs_QC = FALSE,
  folder_results = FALSE
))

qc_df <- as.data.frame(exprs(qc_frame))
flag <- ifelse(qc_df[["QCvector"]] >= 10000, "fail", "pass")
qc_df <- cbind(qc_df["QCvector"], flag, .ci = (1:nrow(qc_df)))

result <- ctx$addNamespace(qc_df)
ctx$save(result)
