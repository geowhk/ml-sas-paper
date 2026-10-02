suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

source("04_current-statistics-comparison/R/run_current_statistics_comparison.R")

결과폴더 <- "04_current-statistics-comparison/runs/revision_v2"
객체경로 <- file.path(결과폴더, "objects", "current_statistics_comparison.rds")

if (!file.exists(객체경로)) {
  stop("저장된 분석 객체를 찾을 수 없다: ", 객체경로, call. = FALSE)
}

객체 <- readRDS(객체경로)
그림저장(객체$results, 객체$correlations, 결과폴더)
message("기존 통계량 비교 그림을 다시 생성했다: ", file.path(결과폴더, "figures"))
