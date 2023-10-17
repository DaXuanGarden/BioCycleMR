#' Retrieve Data for 731 Immune Cells
#'
#' This function retrieves data for 731 immune cells. It can either return
#' preprocessed data directly or reprocess the data using provided parameters.
#'
#' @param use_preprocessed Logical. If TRUE, the function will return the preprocessed data.
#'                         If FALSE, the function will process the data using the provided parameters.
#'                         Default is TRUE.
#' @param p1 Threshold parameter. Default is 1e-5.
#' @param p2 Second threshold parameter. Default is 5e-08.
#' @param count_try_max Maximum number of tries. Default is 50.
#' @param r2 R squared threshold parameter. Default is 0.001.
#' @param kb Distance threshold parameter in kilobases. Default is 10000.
#' @param mc_cores Number of cores to use for parallel processing. Default is 30.
#'
#' @return A data frame with extracted data for the 731 immune cells.
#' @export
#'
#' @examples
#' \dontrun{
#' # Using preprocessed data
#' immc_data_preprocessed <- get_immc()
#'
#' # Processing data with custom parameters
#' immc_data_custom <- get_immc(use_preprocessed = FALSE, p1 = 5e-5, mc_cores = 15)
#' }

get_immc <- function(use_preprocessed = TRUE, p1 = 1e-05, p2 = 5e-08, count_try_max = 50, r2 = 0.001, kb = 10000, mc_cores = 30) {

  if (use_preprocessed) {
    # Load preprocessed data
    preprocessed_data <- readRDS(system.file("data", "dx_immu_cell_raw_df.rda", package = "BioCycleMR"))
    return(preprocessed_data)
  }

  # 加载所需的库
  library(tidyverse)
  library(data.table)
  library(TwoSampleMR)
  library(rlist)

  # 从包内的数据加载imc731id
  imc731id <- readRDS(system.file("data", "imc731id_data.rda", package = "BioCycleMR"))

  online_id <- function(id, p1 = 1e-05, p2 = 5e-08, count_try_max = 50, r2 = 0.001, kb = 10000) {
    immu.cell <- data.frame()
    count_try <- 0

    repeat {
      count_try <- count_try + 1

      try({
        cat(count_try, "try", id, "\n")
        immu.cell <- extract_instruments(id, p1 = p1, p2 = p2, r2 = r2, kb = kb)
      })

      if ((exists("immu.cell") && !is.null(immu.cell) && nrow(immu.cell) > 0) || count_try == count_try_max) {
        break
      }

      Sys.sleep(0.5)
    }

    return(immu.cell)
  }

  # 并行处理
  immu_cell_raw <- parallel::mclapply(imc731id$id, function(id) {
    online_id(id, p1 = p1, p2 = p2, count_try_max = count_try_max, r2 = r2, kb = kb)
  }, mc.cores = mc_cores)

  # 遍历 immu_cell_raw 列表中的每个元素
  for (i in seq_along(immu_cell_raw)) {
    # 获取与当前列表元素对应的 Trait 和 Trait||id.exposure
    trait_id_exposure <- paste(imc731id$Trait[i], imc731id$id[i], sep=" || ")

    # 更新 exposure 列
    immu_cell_raw[[i]]$exposure <- trait_id_exposure
  }

  immu_cell_raw <- do.call(rbind, immu_cell_raw)
  return(immu_cell_raw)
}
