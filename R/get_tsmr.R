#' get_tsmr function
#'
#' This function is intended to be used as a part of BioCycleMR package.
#' It performs multiple complex operations on the input data.
#'
#' @param immu_cell_f_select Input exposure data, expected to be a list of data frames.
#' @param finn_r_dir Directory containing .rda files for each outcome.
#' @param cores Number of cores to use for parallel operations. Default is 64.
#'
#' @return This function does not return anything but creates multiple directories, .rda and .csv files as a part of its operation.
#' @export
#'
#' @examples
#' get_tsmr(immu_cell_f_select, finn_r_dir, cores = 64)
#'
ivw_Speed <- function(dat,prop_var_explained = T,sample=NULL) {
  mr_res <- mr(dat,
               method_list =c("mr_wald_ratio","mr_ivw"))
  if (is.null(sample) || is.na(sample)) {
    pve <- dat %>%
      dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>%
      dplyr::group_by(id.exposure) %>%
      dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  }else{
    pve <- dat %>%
      dplyr::select(id.exposure, beta.exposure, se.exposure) %>%
      dplyr::group_by(id.exposure) %>%
      dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + sample*se.exposure^2)))
  }


  if(prop_var_explained)
  {
    mr_res <- mr_res %>%
      dplyr::left_join(pve, by = "id.exposure")
  }

  return(mr_res)
}
mr_Speed <- function(dat = harmonised_dat, prop_var_explained = T,sample=NULL)
{
  mr_res <- mr(dat)
  if (is.null(sample) || is.na(sample)) {
    pve <- dat %>%
      dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>%
      dplyr::group_by(id.exposure) %>%
      dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  }else{
    pve <- dat %>%
      dplyr::select(id.exposure, beta.exposure, se.exposure) %>%
      dplyr::group_by(id.exposure) %>%
      dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + sample*se.exposure^2)))
  }
  if(prop_var_explained)
  {
    mr_res <- mr_res %>%
      dplyr::left_join(pve, by = "id.exposure")
  }

  return(mr_res)
}
Presso_Speed <- function(dat) {
  res=TwoSampleMR::run_mr_presso(
    dat)
  return(res)
}

FDR_correct_immune <- function(dat=iddf,sample=731,outcome_id){
  if(is.null(dat$pval[1])==T || is.na(dat$pval[1])==T){print("数据不包含pvalue，无法计算FDR")
    return(dat)
  }else{
    data_sorted = dat
    # 使用 order() 函数按照 pvalue 列升序排列数据框
    # data_sorted <- dat[order(dat$pval), ]
    # data_sorted$rank = c(1:dim(data_sorted)[1])
    data_sorted$FDR =p.adjust(data_sorted$pval,method = "BH",n = sample)
    # data_sorted[data_sorted[,"FDR"]>1,"FDR"] <-1
    # table1 <- data_sorted %>%
    #   filter(FDR < 0.05)
    # if (!dir.exists(outcome_id)){
    #   dir.create(outcome_id)
    # }
    # if (dim(table1)[1] == 0){
    #   cat("no significant!")
    # }else{
    #   write.csv(table1,file =paste0(outcome_id,"/hub.csv"))
    #   write.csv(data_sorted,file =paste0(outcome_id,"/FDR.csv"))
    # }
    return(data_sorted)
  }
}
p.adjust <- function(p, method = p.adjust.methods, n = length(p))
{
  ## Methods 'Hommel', 'BH', 'BY' and speed improvements
  ## contributed by Gordon Smyth
  method <- match.arg(method)
  if(method == "DB") method <- "bonferroni" # back compatibility
  if(method == "fdr") method <- "BH"	# back compatibility
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if(all(nna <- !is.na(p))) nna <- TRUE
  p <- p[nna]
  lp <- length(p)
  stopifnot(n >= lp)
  if (n <= 1) return(p0)
  if (n == 2 && method == "hommel") method <- "hochberg"

  p0[nna] <-
    switch(method,
           bonferroni = pmin(1, n * p),
           holm = {
             i <- seq_len(lp)
             o <- order(p)
             ro <- order(o)
             pmin(1, cummax( (n+1L - i) * p[o] ))[ro]
           },
           hommel = { ## needs n-1 >= 2 in for() below
             if(n > lp) p <- c(p, rep.int(1, n-lp))
             i <- seq_len(n)
             o <- order(p)
             p <- p[o]
             ro <- order(o)
             q <- pa <- rep.int( min(n*p/i), n)
             for (j in (n-1L):2L) {
               ij <- seq_len(n-j+1L)
               i2 <- (n-j+2L):n
               q1 <- min(j*p[i2]/(2L:j))
               q[ij] <- pmin(j*p[ij], q1)
               q[i2] <- q[n-j+1L]
               pa <- pmax(pa, q)
             }
             pmax(pa, p)[if(lp < n) ro[1L:lp] else ro]
           },
           hochberg = {
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             pmin(1, cummin( (n+1L - i) * p[o] ))[ro]
           },
           BH = {
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             pmin(1, cummin( n / i * p[o] ))[ro]
           },
           BY = {
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             q <- sum(1/(1L:n))
             pmin(1, cummin(q * n / i * p[o]))[ro]
           },
           none = p)
  p0
}

get_tsmr <- function(immu_cell_f_select, finn_r_dir, cores = 64) {
  # Combine all exposure data frames into one
  exposure_dat = do.call(rbind, immu_cell_f_select)

  # Get list of all .rda files in finn_r_dir
  finn_list <- list.files(finn_r_dir, full.names = TRUE, recursive = FALSE)
  finn_name_list <- sapply(list.files(finn_r_dir, full.names = FALSE, recursive = FALSE), FUN = function(x){sub("\\.rda", "", x)})

  # Create directories for each outcome
  for (i in finn_name_list){
    if (!dir.exists(i)){
      dir.create(i)
    }
  }

  for (i in c(1:length(finn_name_list))){
    filename = finn_list[i]
    outcome_name = finn_name_list[i]

    # Load .rda file
    load(filename)
    cat("load",filename,"\n")

    outcome_dat = finndata
    cat("harmonise data for ",outcome_name,"\n")

    # Harmonise data
    dat<-harmonise_data(exposure_dat=exposure_dat, outcome_dat=outcome_dat)
    har_true=dat[dat$mr_keep=="TRUE",]
    write.csv(har_true, file=paste0(outcome_name,"/SNP.csv"), row.names=F)
    save(har_true,file = paste0(outcome_name,"/SNP.rda"))

    # Split harmonised data by exposure
    dataList = split(har_true, har_true$id.exposure)
    cat("run heterogeneity for ",outcome_name,"\n")

    # Run heterogeneity test
    res_heterogeneity <- parallel::mclapply(dataList, mr_heterogeneity, mc.cores = cores/2)
    save(res_heterogeneity,file = paste0(outcome_name,"/heterogeneity.rda"))
    heterogeneity = do.call(rbind,res_heterogeneity)
    write.csv(heterogeneity, file=paste0(outcome_name,"/heterogeneity.csv"), row.names=F)

    #多效性检验
    cat("run pleiotropy for ",outcome_name,"\n")
    res_pleiotropy<- parallel::mclapply(dataList, mr_pleiotropy_test, mc.cores = 32)
    save(res_pleiotropy,file = paste0(outcome_name,"/pleiotropy.rda"))
    pleiotropy = do.call(rbind,res_pleiotropy)
    write.csv(pleiotropy, file=paste0(outcome_name,"/pleiotropy.csv"), row.names=F)

    # 进行水平多效性
    filtered_dataList <- lapply(dataList, function(df) {
      if (nrow(df) >= 4) {
        return(df)
      } else {
        return(NULL)
      }
    })
    # 清除列表中的 NULL 元素
    filtered_dataList <- filtered_dataList[!sapply(filtered_dataList, is.null)]
    # filtered_data<-do.call(rbind, filtered_dataList[1:2])
    # MRPRESSO <- run_mr_presso(dat = filtered_data)
    cat("run presso for ",outcome_name,"\n")
    MRPRESSO <-parallel::mclapply(filtered_dataList, run_mr_presso, mc.cores = 32)
    presso_result_list = lapply(MRPRESSO,FUN = function(mr){
      result =mr[[1]][["Main MR results"]]
      ebi_id = names(mr)[1]
      result$id.exposure = ebi_id
      return(result)
    })
    res_MRPRESSO = do.call(rbind, presso_result_list)
    save(MRPRESSO,file = paste0(outcome_name,"/presso.rda"))
    write.csv(res_MRPRESSO, file=paste0(outcome_name,"/presso.csv"), row.names=F)

    # mr result
    cat("run mr for ",outcome_name,"\n")
    list_mr = parallel::mclapply(dataList, mr_Speed, mc.cores = 32)
    res_df <- do.call(rbind, list_mr)
    save(list_mr,file = paste0(outcome_name,"/mr.rda"))
    write.csv(res_df, file=paste0(outcome_name,"/mr.csv"), row.names=F)


    # pajust result
    cat("extract effect exposure for ",outcome_name,"\n")
    list_ivw = parallel::mclapply(dataList, ivw_Speed, mc.cores = 32)
    effect_df <- do.call(rbind, list_ivw)
    data_sorted <- FDR_correct_immune(effect_df,sample = 731)
    # data_sorted$p.ajust = p.adjust(p = data_sorted$pval,n=731,method = "fdr")
    table1 <- data_sorted %>%
      filter(FDR <= 0.05,
             method %in% c("Wald ratio","Inverse variance weighted"))
    save(data_sorted,table1,file = paste0(outcome_name,"/effect.rda"))
    write.csv(data_sorted, file=paste0(outcome_name,"/effect.csv"), row.names=F)
    # write.csv(table1, file=paste0(outcome_name,"/fdr.csv"), row.names=F)
    cat("finished: ",outcome_name,"\n")
    sig_n = sig_n+dim(table1)[1]
  }
}
