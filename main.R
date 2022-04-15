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

matrix2flowFrame <- function(a_matrix){ 
  
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
  flowFrame <- flowCore::flowFrame(a_matrix, params)
  
  return(flowFrame)
}


flowAI_QC<-function(fc_frame,input.pars){
  if(input.pars$Detailed == FALSE){
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
    qc_df <- cbind(flag, .ci = (0:(nrow(qc_df)-1)))
    
    return(qc_df)
  }  
  else if (input.pars$Detailed == TRUE) {
    
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
    
    return(qc_df)
    
  }
}

#################
ctx <- tercenCtx()

if(ctx$cnames[1] == "filename") {filename <- TRUE
if(ctx$cnames[2] != "Time") stop("Time not detected in the second column.")
}else{filename <- FALSE
if(ctx$cnames[1] != "Time") stop("filename or Time not detected in the top column.")
}

celldf <- ctx %>% dplyr::select(.ri, .ci) 
if(nrow(celldf) != length(table(celldf)))stop("There are multiple values in one of the cells.")

input.pars <- list(
  Detailed = ifelse(is.null(ctx$op.value('Detailed')), TRUE, as.logical(ctx$op.value('Detailed'))),
  second_fractionFR = ifelse(is.null(ctx$op.value('second_fractionFR')), 0.1, as.double(ctx$op.value('second_fractionFR'))),
  alphaFR = ifelse(is.null(ctx$op.value('alphaFR')), 0.01, as.double(ctx$op.value('alphaFR'))),
  decompFR = ifelse(is.null(ctx$op.value('decompFR')), TRUE, as.logical(ctx$op.value('decompFR'))),
  outlier_binsFS = ifelse(is.null(ctx$op.value('outlier_binsFS')), FALSE, as.logical(ctx$op.value('outlier_binsFS'))),
  pen_valueFS = ifelse(is.null(ctx$op.value('pen_valueFS')), 500, as.double(ctx$op.value('pen_valueFS'))),
  max_cptFS = ifelse(is.null(ctx$op.value('max_cptFS')), 3, as.double(ctx$op.value('max_cptFS'))),
  sideFM = ifelse(is.null(ctx$op.value('sideFM')), "both", ctx$op.value('sideFM')),
  neg_valuesFM = ifelse(is.null(ctx$op.value('neg_valuesFM')), 1, as.double(ctx$op.value('neg_valuesFM')))
)

marker_names <- ctx$rselect()
marker_list <- lapply(seq_len(nrow(marker_names)), function(i) marker_names[[i,1]])

if(filename == TRUE){
  data <- ctx$as.matrix() %>% t() 
  colnames(data) <- marker_list
  data <- data %>% cbind((ctx$cselect(ctx$cnames[[2]]))) %>% cbind((ctx$cselect(ctx$cnames[[1]])))
}
if(filename == FALSE){
  data <- ctx$as.matrix() %>% t() 
  colnames(data) <- c(marker_list)
  data <- data %>% cbind((ctx$cselect(ctx$cnames[[1]])))
  data$filename <- "singlefile"
}

filenames <- unique(data$filename)
qc_df <- data.frame(matrix(ncol=0, nrow=nrow(data)))


QC_allfiles <- do.call("rbind",lapply(filenames, function(x) {
  singlefiledata <- data[data$filename == x,]
  singlefileflowframe <- singlefiledata[1:(ncol(singlefiledata)-1)] %>% as.matrix() %>% matrix2flowFrame()
  
  QC_indice <- flowAI_QC(singlefileflowframe, input.pars)
  QC_indice["filename"]<-x
  return(QC_indice)
}))


flowai_res_QC <- cbind(QC_allfiles, .ci = (0:(nrow(qc_df)-1)))
ctx$addNamespace(flowai_res_QC) %>% ctx$save()


#####





