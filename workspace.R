library(tercen)
library(dplyr)
library(tibble)
library(flowCore)
library(flowAI)

convert_to_seconds <- function(t) {
  t<- t[[1, drop = TRUE]]
  if (length(t) > 2) {
    t2 <-(t[2])
    unit_size <- ceiling(log10(t2))
    t <- t/(10^unit_size)
    return(enframe(t, name = NULL))
  } else
    return(t)
}


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

ctx <- tercenCtx(workflowId = "2b5291026ae8514a1a77d7d9190002cc",
                 stepId = "e260a25e-5ae5-4f58-9743-e66e7093becb")


input.pars <- list(
  second_fractionFR = ifelse(is.null(ctx$op.value('second_fractionFR')), 0.1, as.double(ctx$op.value('second_fractionFR'))),
  alphaFR = ifelse(is.null(ctx$op.value('alphaFR')), 0.01, as.double(ctx$op.value('alphaFR'))),
  decompFR = ifelse(is.null(ctx$op.value('decompFR')), TRUE, as.logical(ctx$op.value('decompFR'))),
  outlier_binsFS = ifelse(is.null(ctx$op.value('outlier_binsFS')), FALSE, as.logical(ctx$op.value('outlier_binsFS'))),
  pen_valueFS = ifelse(is.null(ctx$op.value('pen_valueFS')), 500, as.double(ctx$op.value('pen_valueFS'))),
  max_cptFS = ifelse(is.null(ctx$op.value('max_cptFS')), 3, as.double(ctx$op.value('max_cptFS'))),
  sideFM = ifelse(is.null(ctx$op.value('sideFM')), "both", ctx$op.value('sideFM')),
  neg_valuesFM = ifelse(is.null(ctx$op.value('neg_valuesFM')), 1, as.double(ctx$op.value('neg_valuesFM')))
)

time <- ctx$cselect(ctx$cnames[[1]])
time <- convert_to_seconds(time)

data <- ctx$as.matrix() %>% t()
data <- as.matrix(cbind(data, time))

fc_frame <- matrix2flowset(data)

detailed = TRUE

if(detailed == TRUE){
qc_FR <- suppressWarnings(flowAI::flow_auto_qc(
  remove_from = "FR",
  fcsfiles = fc_frame,
  output = 3,
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

qc_FS <- suppressWarnings(flowAI::flow_auto_qc(
  remove_from = "FS",
  fcsfiles = fc_frame,
  output = 3,
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

qc_FM <- suppressWarnings(flowAI::flow_auto_qc(
  remove_from = "FM",
  fcsfiles = fc_frame,
  output = 3,
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

if(detailed == TRUE){
qc_FR <- suppressWarnings(flowAI::flow_auto_qc(
  remove_from = "FR",
  fcsfiles = fc_frame,
  output = 3,
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

qc_FS <- suppressWarnings(flowAI::flow_auto_qc(
  remove_from = "FS",
  fcsfiles = fc_frame,
  output = 3,
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

qc_FM <- suppressWarnings(flowAI::flow_auto_qc(
  remove_from = "FM",
  fcsfiles = fc_frame,
  output = 3,
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

qc_df <- data.frame(matrix(ncol=0, nrow=nrow(data)))

qc_FR_list <- as.list(qc_FR[[1]])
qc_FS_list <- as.list(qc_FS[[1]])
qc_FM_list <- as.list(qc_FM[[1]])

qc_df$FlowRate_flag <- ifelse(rownames(qc_df) %in% qc_FR_list, "fail", "pass")
qc_df$SigAcqtn_flag <- ifelse(rownames(qc_df) %in% qc_FS_list, "fail", "pass")
qc_df$DynRange_flag <- ifelse(rownames(qc_df) %in% qc_FM_list, "fail", "pass")

qc_df$QC_flag <- paste(qc_df$FlowRate_flag, qc_df$SigAcqtn_flag, qc_df$DynRange_flag, sep = "")
qc_df$QC_flag <- ifelse(qc_df$QC_flag == "passpasspass", "pass", "fail")

flowAI_QC <- cbind(qc_df, .ci = (0:(nrow(qc_df)-1)))

result <- ctx$addNamespace(flowAI_QC)
ctx$save(result)
}
}

