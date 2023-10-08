#' Title: Immunce cells as exposure
#'
#' @param outcome_id
#' @param max_retries
#'
#' @return
#' @export
#'
#' @examples
#' for (outcome_id in outcome_ids){
#' write.table(analyze_outcome_immune(outcome_id=outcome_id,max_retries = 50),file = paste0(outcome_id,".txt"),sep = "\t",quote = F,row.names = F)
#' }
analyze_outcome_immune <- function(outcome_id,max_retries=50) {
  cat("!!Immune cells: Analyzing outcome ", outcome_id, "\n")
  result = data.frame()
  max_retries = max_retries  # 设置最大重试次数
  for (i in imc) {
    retry_count = 0  # 初始化重试计数器
    success1 = FALSE  # 标志，表示是否成功连接并获得结果
    success2 = FALSE
    while (retry_count < max_retries && !success1) {
      retry_count = retry_count + 1

      tryCatch({
        # 在这里执行您的连接和分析操作
        # 如果成功，将结果存储在 result_or 中，并将 success 标志设置为 TRUE
        expo_rt <- extract_instruments(outcomes = i, p1 = 1e-5,
                                       clump = T, p2 = 1e-5,
                                       r2 = 0.001, kb = 10000,access_token=NULL)
        success1 = TRUE  # 成功获得结果
      },
      error = function(e) {
        # 捕获到错误时的操作
        # 输出错误信息
        cat("Error1:", conditionMessage(e), "\n")

        if (retry_count < max_retries) {
          cat("Your network is unstable: Retry getting exposure tool variables...\n")
          Sys.sleep(1)  # 等待一段时间后重试
        } else {
          cat("Failed to complete the connection because the maximum number of retries was reached.\n")
        }
      })
    }
    while (retry_count < max_retries && !success2) {
      retry_count = retry_count + 1
      tryCatch({
        outc_rt <- extract_outcome_data(snps = expo_rt$SNP,
                                        outcomes = outcome_id,access_token = NULL)

        harm_rt <- harmonise_data(exposure_dat = expo_rt,
                                  outcome_dat = outc_rt, action = 2)

        mr_result <- mr(harm_rt,method_list = c("mr_wald_ratio","mr_ivw"))
        print(mr_result)
        # result_or <- generate_odds_ratios(mr_result)

        # 如果 p 值小于 0.05，则将结果添加到 result 数据框中
        table1 <- mr_result %>%
          filter(pval < 0.05,
                 method %in% c("Wald ratio","Inverse variance weighted"))

        table1 <- table1 %>%
          generate_odds_ratios()%>%
          mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95),
                 `P value` = scales::scientific(pval))

        result <- rbind(result, table1)

        success2 = TRUE  # 成功获得结果
      }, error = function(e) {
        # 捕获到错误时的操作
        # 输出错误信息
        cat("Error2:", conditionMessage(e), "\n")

        if (retry_count < max_retries) {
          cat("Your network is unstable: Retry getting outcome tool variables...\n")
          Sys.sleep(1)  # 等待一段时间后重试
        } else {
          cat("Failed to complete the connection because the maximum number of retries was reached.。\n")
        }
      })
    }
  }
  cat("Finished:", outcome_id, "\n")
  return(result)

}
