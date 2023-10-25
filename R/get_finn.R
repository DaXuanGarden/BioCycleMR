#' Harmonise and process FinnGen R9 data.
#'
#' This function processes and harmonises the FinnGen R9 data. The function reads each file,
#' formats and harmonises it, and then saves the processed data in a specified directory.
#' It uses parallel computation to process the list of files.
#'
#' @param finn_dir The directory where the FinnGen data files are located.
#' @param save_dir Directory where processed data should be saved. Default is "finn_r".
#' @param cores The number of cores to use for parallel computation. Default is 2.
#'
#' @examples
#' \dontrun{
#' get_finn(finn_dir = "finn", save_dir = "finn_r", cores = 30)
#' }
#' @export

get_finn <- function(finn_dir = "finn", save_dir = "finn_r", cores = 30) {
  # 加载内置数据
  data("Finn_R9_data")

  # 定义逐文件处理的子函数
  process_finn <- function(trait) {
    # 1. 读取数据文件
    finn <- data.table::fread(file = trait, data.table = F)

    # 删除重复的SNPs，保留p值最小的
    finn <- finn[order(finn$pval),]
    finn <- finn[!duplicated(finn$rsids), ]

    # 2. 格式化数据
    trait_basename <- sub(paste0(finn_dir, "/"), "", trait)  # 移除文件夹名称
    trait_basename <- sub("\\.gz", "", trait_basename)  # 移除.gz后缀

    finndata <- TwoSampleMR::format_data(finn, type = 'outcome', snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
                                         eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "ref",
                                         pval_col = "pval", gene_col = "nearest_genes", chr_col = "#chrom",
                                         pos_col = "pos")

    # 3. 提取trait相关的行信息，精准匹配
    trait_row <- finn_info[grepl(paste0("\\b", trait_basename, "\\b"), finn_info$path_https),]
    if (nrow(trait_row) == 0) {
      stop(paste("Couldn't find phenotype details for:", trait_basename))
    }

    # 4. 更新finndata的列数据
    finndata$ncase.outcome <- trait_row$num_cases
    finndata$ncontrol.outcome <- trait_row$num_controls
    finndata$samplesize.outcome <- trait_row$num_cases + trait_row$num_controls

    # 5. 更新finndata的outcome列
    outcome <- trait_row$name
    finndata$outcome <- outcome
    ##赋值id.exposure😊
    finndata$id.outcome <- finndata$outcome
    # 6. 保存数据
    save(finndata, file = paste0(save_dir, "/", "finndataR9_", trait_basename, ".rda"))

    # 7. 输出完成信息
    cat("Finished:", trait_basename, "\n")
  }

  # 获取文件列表
  finn_list <- list.files(finn_dir, pattern = "\\.gz$", full.names = TRUE, recursive = FALSE)

  # 使用并行处理文件列表
  parallel::mclapply(finn_list, process_finn, mc.cores = cores)
}
