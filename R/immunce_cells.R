#' Title: Immunce cells as exposure
#'
#' @param outcome_id
#' @param max_retries
#'
#' @return result
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

#' Title immune_total_check
#'
#' @param exposure_id
#' @param outcome_name
#' @param outcome_source
#'
#' @export
#'
#' @exampleslibrary(dplyr)
#' library(TwoSampleMR)
#' #https://gwas.mrcieu.ac.uk/datasets/
#' outcome_name = "finn-b-N14_ENDOMETRIOSIS"
#' outcome_source = "finn"
#' iddf=read.table(paste0(outcome_name,".txt"),header = T,sep = "\t")
#' # 暴露数据的SNP，工具变量
#' exposure_id=iddf$id.exposure
#' immune_total_check(exposure_id=exposure_id,
#' outcome_name=outcome_name,
#' outcome_source=outcome_source)
immune_total_check<-function(exposure_id,
                             outcome_name,
                             outcome_source){
  max_retries <- 50  # 设置最大重试次数
  retry_count <- 0  # 初始化重试计数器
  success <- FALSE  # 初始化成功标志
  exposure_dat <- NULL  # 初始化 exposure_dat 变量
  print("Extracting instruments varibales")
  while (retry_count < max_retries && !success) {
    retry_count <- retry_count + 1

    tryCatch({
      # 在这里执行网络请求
      exposure_dat <- extract_instruments(exposure_id, p1 = 5e-6,
                                          clump = T,
                                          # p2 = 1e-5,
                                          r2 = 0.001, kb = 10000, access_token = NULL)
      success <- TRUE  # 成功获得结果
    }, error = function(e) {
      # 捕获到错误时的操作
      # 输出错误信息
      cat("Error1:", conditionMessage(e), "\n")

      if (retry_count < max_retries) {
        cat("On the ",retry_count, "attempt,network connection failed: Retrying...\n")
        Sys.sleep(5)  # 休眠5秒后重试
      } else {
        cat("Unable to complete connection, maximum number of retries reached.\n")
      }
    })
  }

  if (success) {
    # 如果成功，可以继续使用 exposure_dat
    print("calculate R2")
    R2a=2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
    R2b=2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
    R2=R2a/(R2a+R2b)
    print("calculate F")
    exposure_dat$F_statistics=R2*(exposure_dat$samplesize.exposure-2)/(1-R2)

    exposure_dat=exposure_dat[exposure_dat$F_statistics>10,]

    # library(readxl)
    # df <- read_excel("imc731id.xlsx")
    # colnames(df)
    # names(df) = c("Trait type","Panel","Statistical trait name","Trait",
    #               "GWAS Catalog Accession Number", "id.exposure" )
    # save(df,file = "data/immune_ref.rda")
    data("immune_ref")
    library(tidyverse)
    exposure_dat= left_join(exposure_dat,df, by = "id.exposure")

    exposure_dat$exposure=NULL
    exposure_dat$exposure=exposure_dat$Panel
    exposure_dat$Panel=NULL
    print("create dir")
    if (!dir.exists(paste0(outcome_name,"/"))) {
      dir.create(paste0(outcome_name,"/"))
    }
    if (!dir.exists(paste0(outcome_name,"/",outcome_source,"/"))) {
      dir.create(paste0(outcome_name,"/",outcome_source,"/"))
    }
    save(exposure_dat,file =paste0(outcome_name,"/",'exposure_dat.Rdata'))
  } else {
    cat("无法完成连接，达到最大重试次数。\n")
  }



  load(paste0(outcome_name,"/",'exposure_dat.Rdata'))
  print("Extracting outcome")
  # #提取结局数据
  # outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcome_name)
  max_retries <- 50  # 设置最大重试次数
  retry_count <- 0  # 初始化重试计数器
  success <- FALSE  # 初始化成功标志
  outcome_dat <- NULL  # 初始化 outcome_dat 变量

  while (retry_count < max_retries && !success) {
    retry_count <- retry_count + 1

    tryCatch({
      # 在这里执行网络请求
      outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcome_name,access_token = NULL)
      harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)

      success <- TRUE  # 成功获得结果
    }, error = function(e) {
      # 捕获到错误时的操作
      # 输出错误信息
      cat("Error2:", conditionMessage(e), "\n")

      if (retry_count < max_retries) {
        cat("On the ",retry_count, "attempt,network connection failed: Retrying...\n")
        Sys.sleep(5)  # 休眠5秒后重试
      } else {
        cat("Unable to complete connection, maximum number of retries reached.\n")
      }
    })
  }

  # if (success) {
  #   # 如果成功，可以继续使用 outcome_dat
  # } else {
  #   cat("无法完成连接，达到最大重试次数。\n")
  # }

  # 取交集
  exposure_dat=exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]
  ## MR分析
  print("MR analysis")
  mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
  {
    mr_res <- mr(dat)

    pve <- dat %>%
      dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>%
      dplyr::group_by(id.exposure) %>%
      dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))

    if(prop_var_explained)
    {
      mr_res <- mr_res %>%
        dplyr::left_join(pve, by = "id.exposure")
    }

    return(mr_res)
  }




  mr_res <- mr_modified(harmonised_dat, prop_var_explained = T)

  save(exposure_dat,outcome_dat,mr_res,harmonised_dat,file =paste0(outcome_name,"/",outcome_source,'/mr_input_res.Rdata'))

  load(paste0(outcome_name,"/",outcome_source,'/mr_input_res.Rdata'))
  print("mr_heterogeneity and mr_pleiotropy_test")
  #  异质性 ，用处不大
  # 当存在显著的异质性（p<0.05）时，结果不太可靠
  heter_dat=mr_heterogeneity(harmonised_dat)
  write.csv(heter_dat, file=paste0(outcome_name,"/",outcome_source,'/table.heterogeneity.csv'), row.names=F)

  # 水平多效性检验，用处不大
  #如果该截距项与0差异很大（pval < 0.05），说明存在水平多效性
  pleio_dat=mr_pleiotropy_test(harmonised_dat)
  write.csv(pleio_dat, file=paste0(outcome_name,"/",outcome_source,'/table.pleiotropy.csv'), row.names=F)


  ### table
  table1 <- mr_res %>%
    filter(pval < 0.05,
           method %in% c("Wald ratio","Inverse variance weighted")) %>%
    left_join(exposure_dat, by = "exposure")

  table1 <- table1 %>%
    generate_odds_ratios()%>%
    mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95),
           `P value` = scales::scientific(pval),
           `PVE` = paste0(sprintf("%.2f",100*pve),"%"),
           `F statistics` = sprintf("%.2f",F_statistics)) %>%
    dplyr::select(Immune = exposure, Trait = Trait, Trait_type = `Trait type`,id=id.exposure.x,
                  SNP, `Effect allele` = effect_allele.exposure,
                  `OR (95% CI)`, `P value`,
                  PVE, `F statistics`)

  ## 保存，后面常用
  save(table1,file =paste0(outcome_name,"/",outcome_source,'/table1.Rdata'))

  library(dplyr)
  library(ggplot2)
  library(ggrepel)

  volcano_plot <- function(.data,
                           number_comparasion = 1,
                           title = "Immune cells",
                           legend.position = "none") {

    p_thershold <- 0.05/number_comparasion

    p <- .data %>%
      mutate(y = -log10(pval),
             label = ifelse(pval < p_thershold, exposure, NA)) %>%
      ggplot(aes(x = b, y = y)) +
      geom_point(aes(size = pve), alpha = 0.5, color = "#0072b5") +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
      theme_classic() +
      theme(panel.grid = element_blank(),
            legend.title = element_text(size = 6.5),
            legend.text = element_text(size = 6.5),
            legend.position = legend.position) +
      labs(x = "ln(OR)",
           y = parse(text = "-log[10]*(italic(P)-value)"),
           title = title) +
      scale_size(name = "PVE",
                 breaks = c(0.2*1:3)) +
      ggrepel::geom_label_repel(aes(label = label), size = 3)
    plot(p)
  }

  # volcano plot
  p = mr_res %>%
    filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>%
    volcano_plot(number_comparasion = 1)
  ggsave(p,file = paste0(outcome_name,"/",outcome_source,'/volcano.pdf'))

  # modify mr_res
  data("immune_ref")
  library(tidyverse)
  mr_res= left_join(mr_res,df, by = "id.exposure")
  mr_res2=mr_res[mr_res$id.exposure %in% table1$id,]
  result_or=generate_odds_ratios(mr_res2)
  write.table(result_or,file = paste0(outcome_name,"/",outcome_source,'/OR.txt'),row.names = F,sep = "\t",quote = F)

  #将统计结果绘制森林图
  library(grid)
  library(forestploter)


  mydata=read.table(file = paste0(outcome_name,"/",outcome_source,'/OR.txt'),header = T,sep = "\t")
  mydata$` ` <- paste(rep(" ", 20), collapse = " ")
  mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                              mydata$or, mydata$or_lci95,
                                                              mydata$or_uci95))
  data_plot = cbind(exposure = mydata$exposure,Trait=mydata$Trait,Trait_type=mydata$Trait.type,method = mydata$method,` `=mydata$` `,nsnp=mydata$nsnp,pval=mydata$pval,`OR (95% CI)`=mydata$`OR (95% CI)`)
  p2 = forest(data_plot,
              est = mydata$or,
              lower =mydata$or_lci95,
              upper = mydata$or_uci95,
              sizes =0.3,
              ci_column =5 ,
              ref_line = 1,
              xlim = c(0.05, 3),
  )
  ggsave(p2,file = paste0(outcome_name,"/",outcome_source,'/forest.pdf'),width = 30,height = 30)
  print(paste0("finished",exposure_id))
}
