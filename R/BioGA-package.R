#' @keywords internal
"_PACKAGE"


## usethis namespace: start
#' @importFrom biocViews getBiocViews
#' @importFrom ggplot2 ggtitle ggplot geom_line labs
#' @importFrom sessioninfo session_info
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom graphics boxplot hist par
#' @importFrom animation saveGIF
#' @importFrom rlang local_options
#' @importFrom BiocStyle html_document
#' @import RcppParallel
#' @importFrom GEOquery getGEO
#' @importFrom caret train trainControl createDataPartition
#' @importFrom lattice xyplot bwplot histogram
#' @importFrom BiocParallel bplapply bpmapply register
#' @importFrom caretEnsemble caretList caretEnsemble
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom dplyr filter select mutate arrange summarise rename distinct
#' @importFrom iml Predictor FeatureImp
#' @importFrom lime explain
#' @importFrom pROC roc ci.auc plot.roc
#' @importFrom pheatmap pheatmap
#' @importFrom randomForest combine randomForest
#' @importFrom survival cluster Surv coxph survfit
#' @importFrom survminer ggsurvplot
#' @importFrom timeROC timeROC
#' @importFrom xgboost slice xgb.DMatrix xgb.train
#' @useDynLib BioGA, .registration = TRUE
## usethis namespace: end
NULL
