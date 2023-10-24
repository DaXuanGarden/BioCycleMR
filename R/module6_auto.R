#' Title auto_find_risks
#' @description
#' find the risk factor of outcome uesd the open GWAS id.
#'
#' @param outcome_id the id of outcome
#' @param outcome_file
#' @param outcome_full
#' @param outcome_name
#' @param num1
#' @param num2
#' @param pval
#' @param clump_r2
#' @param clump_kb
#' @param proxies
#' @param pop
#' @param renew_id
#' @param action
#' @param save_path
#'
#' @return risk_data
#' @export
#'
#' @examples
Find_risks_auto<-function (outcome_id = NULL,
                           outcome_file = NULL,
                           outcome_full = NULL,
                           outcome_name,
                           num1 = 1,
                           num2 = 18115,
                           pval = 5e-08,
                           clump_r2 = 0.001,
                           clump_kb = 10000,
                           proxies = TRUE,
                           pop = "European",
                           renew_id = FALSE,
                           action = 3,
                           save_path = "./test/")
{
  if (renew_id) {
    ao <- TwoSampleMR::available_outcomes(access_token = NULL)
  }
  else {
    data("ao")
  }
  mrbase_id <- ao %>% dplyr::arrange(., id) %>% subset(.,population ==pop)
  eqtl_id <- mrbase_id[grep("eqtl-a-*", mrbase_id$id), ]
  mrbase_id <- mrbase_id[!mrbase_id$id %in% eqtl_id$id, ]
  if (!is.null(outcome_id)) {
    mrbase_id <- mrbase_id[!mrbase_id$id %in% outcome_id,
    ]
  }
  trait <- mrbase_id %>% dplyr::arrange(., id) %>% dplyr::select(id,trait) %>%
    dplyr::filter(id != c("ieu-a-819"))%>%
    dplyr::select(trait)
  id <- mrbase_id %>% dplyr::arrange(., id) %>% dplyr::select(id) %>%
    dplyr::filter(id != c("ieu-a-819"))

  risk_data = data.frame()
  for (i in num1:num2) {
    # exp <- TwoSampleMR::extract_instruments(outcomes = id$id[i],
    #                                         p1 = pval, r2 = clump_r2,
    #                                         kb = clump_kb,
    #                                         clump = TRUE,
    #                                         access_token = NULL)
    exp = online_id_exposure(id$id[i], p1 = pval,  count_try_max = 15, r2 = clump_r2, kb = clump_kb)
    print(paste0(i, ": ", id$id[i]))
    if (!is.null(exp)) {
      # exp <- TwoSampleMR::clump_data(dat = exp, clump_kb = clump_kb,
      #                                clump_r2 = clump_r2)
      if (!is.null(outcome_id)) {
        out <- online_id_outcome(snps = exp$SNP,
                                                 outcomes = outcome_id, proxies = proxies, access_token = NULL)
      }
      else if (!is.null(outcome_file)) {
        outcome_full <- readRDS(outcome_file)
        out <- merge(outcome_full, exp, by = "SNP") %>%
          dplyr::select(SNP, grep("outcome", colnames(.)))
      }
      else {
        out <- merge(outcome_full, exp, by = "SNP") %>%
          dplyr::select(SNP, grep("outcome", colnames(.)))
      }
      if (!is.null(out)) {
        dat <- TwoSampleMR::harmonise_data(exposure_dat = exp,
                                           outcome_dat = out, action = action)
        dat <- subset(dat, dat$mr_keep == TRUE)
        if (nrow(dat) != 0) {
          library(TwoSampleMR)
          res <- GagnonMR::primary_MR_analysis(dat = dat)
          if (nrow(res) != 0) {
            print(res)
            if (res$pval < 0.05) {
              print(res)
              print(i)
              utils::write.table(res, file = paste0(save_path,
                                                    "/", id$id[i], "_MRbase_id.txt"),
                                 quote = FALSE, col.names = FALSE, append = TRUE,
                                 row.names = FALSE, sep = "\t")
              risk_data = rbind(risk_data,cbind(id$id[i],trait$trait[i]))

            }
          }
        }
      }
    }
  }
  return(risk_data)
}


#' Title
#'
#' @param snps
#' @param outcomes
#' @param proxies
#' @param access_token
#'
#' @return
#' @export
#'
#' @examples
online_id_outcome <- function(snps = NULL,
                              outcomes = NULL,
                              proxies = NULL,
                              access_token = NULL,
                              count_try_max = 20) {
  immu.cell <- data.frame()
  count_try <- 0

  repeat {
    count_try <- count_try + 1

    try({
      cat(count_try, "try", outcomes, "\n")
      immu.cell <- TwoSampleMR::extract_outcome_data(snps = snps,
                                                            outcomes = outcomes, proxies = proxies, access_token = access_token)
    })

    if ((exists("immu.cell") && !is.null(immu.cell) && !nrow(immu.cell)==0  )|| count_try == count_try_max) {
      break
    }

    Sys.sleep(0.5)
  }

  return(immu.cell)

}


#' Title online_id
#'
#' @paraidm
#' @param p1
#' @param p2
#' @param count_try_max
#' @param r2
#' @param kb
#'
#' @return
#' @export
#'
#' @examples
online_id_exposure <- function(id, p1 = 5e-08, p2 = 5e-08, count_try_max = 20, r2 = 0.001, kb = 10000) {
  immu.cell <- data.frame()
  count_try <- 0

  repeat {
    count_try <- count_try + 1

    try({
      cat(count_try, "try", id, "\n")
      immu.cell <- TwoSampleMR::extract_instruments(id, p1 = p1, p2 = p2, r2 = r2, kb = kb,access_token = NULL)
    })

    if ((exists("immu.cell") && !is.null(immu.cell) )|| count_try == count_try_max) {
      break
    }

    Sys.sleep(0.5)
  }

  return(immu.cell)

}
#' Title
#'
#' @param exposure_id
#' @param exposure_file
#' @param exposure_data
#' @param exposure_name
#' @param num1
#' @param num2
#' @param pval
#' @param clump_r2
#' @param clump_kb
#' @param proxies
#' @param pop
#' @param renew_id
#' @param action
#' @param save_path
#'
#' @return
#' @export
#'
#' @examples
auto_find_outcome<-function (exposure_id = NULL, exposure_file = NULL, exposure_data = NULL,
                             exposure_name, num1 = 1, num2 = 18115, pval = 5e-08, clump_r2 = 0.001,
                             clump_kb = 10000, proxies = TRUE, pop = "European", renew_id = FALSE,
                             action = 3, save_path = "D:")
{
  if (renew_id) {
    ao <- TwoSampleMR::available_outcomes()
  }
  else {
    data("ao")
  }
  mrbase_id <- ao %>% dplyr::arrange(., id) %>% subset(., population ==
                                                         pop)
  eqtl_id <- mrbase_id[grep("eqtl-a-*", mrbase_id$id), ]
  mrbase_id <- mrbase_id[!mrbase_id$id %in% eqtl_id$id, ]
  id <- mrbase_id %>% dplyr::arrange(., id) %>% dplyr::select(id) %>%
    dplyr::filter(!id %in% exposure_id)
  if (!is.null(exposure_id)) {
    exp <- TwoSampleMR::extract_instruments(outcomes = exposure_id,
                                            p1 = pval, r2 = clump_r2, kb = clump_kb, clump = TRUE)
    exp <- TwoSampleMR::clump_data(exp, clump_kb = clump_kb,
                                   clump_r2 = clump_r2)
  }
  else if (!is.null(exposure_file)) {
    exp <- vroom::vroom(exposure_file)
  }
  else if (!is.null(exposure_data)) {
    exp <- exposure_data
  }
  if (!is.null(exp) || nrow(exp) != 0) {
    for (i in num1:num2) {
      print(paste0(i, ": ", id$id[i]))
      out <- TwoSampleMR::extract_outcome_data(snps = exp$SNP,
                                               outcomes = id$id[i], proxies = proxies)
      if (!is.null(out)) {
        dat <- TwoSampleMR::harmonise_data(exposure_dat = exp,
                                           outcome_dat = out, action = action)
        dat <- subset(dat, dat$mr_keep == TRUE)
        if (nrow(dat) != 0) {
          library(TwoSampleMR)
          res <- GagnonMR::primary_MR_analysis(dat = dat)
          if (nrow(res) != 0) {
            print(res)
            if (any(res$pval < 0.05)) {
              print(res)
              print(i)
              utils::write.table(res, file = paste0(save_path,
                                                    "/", exposure_name, "_MRbase_id.txt"),
                                 quote = FALSE, col.names = FALSE, append = TRUE,
                                 row.names = FALSE, sep = "\t")
            }
          }
        }
      }
    }
  }
}


#' Title
#'
#' @param exposure_id
#' @param exposure_file
#' @param exposure_name
#' @param outcome_id
#' @param outcome_file
#' @param outcome_name
#' @param pop
#' @param num1
#' @param num2
#' @param clump_kb
#' @param clump_r2
#' @param pval
#' @param action
#' @param palindromic
#' @param proxies_exp
#' @param proxies_med
#' @param proxies_out
#' @param remove_HLA
#' @param two_step_MR
#' @param remove_exp_med_intersect
#' @param mvmr_mediation
#' @param mvmr_med_sig
#' @param mvmr_exp_med_sig_all
#' @param save_path
#'
#' @return
#' @export
#'
#' @examples
auto_find_med2<-function (exposure_id = NULL, exposure_file = NULL, exposure_name,
                          outcome_id = NULL, outcome_file = NULL, outcome_name, pop = "European",
                          num1 = 1, num2 = 18113, clump_kb = 10000, clump_r2 = 0.001,
                          pval = 5e-08, action = 3, palindromic = TRUE, proxies_exp = TRUE,
                          proxies_med = TRUE, proxies_out = TRUE, remove_HLA = TRUE,
                          two_step_MR = FALSE, remove_exp_med_intersect = FALSE, mvmr_mediation = FALSE,
                          mvmr_med_sig = FALSE, mvmr_exp_med_sig_all = FALSE, save_path)
{
  data("ao")
  mrbase_id <- ao %>% dplyr::arrange(., id) %>% subset(., population ==
                                                         pop)
  eqtl_id <- mrbase_id[grep("eqtl-a-*", mrbase_id$id), ]
  mrbase_id <- mrbase_id[!mrbase_id$id %in% eqtl_id$id, ]
  if (!file.exists(paste0(save_path, "/", "mrbse_mediator_id.csv"))) {
    utils::write.csv(x = mrbase_id, file = paste0(save_path,
                                                  "/", "mrbse_mediator_id.csv"))
  }
  if (!is.null(exposure_id)) {
    mrbase_id <- subset(mrbase_id, mrbase_id$id != exposure_id)
    if (!is.null(outcome_id)) {
      mrbase_id <- subset(mrbase_id, mrbase_id$id != outcome_id)
      for (i in num1:num2) {
        mediator_id <- mrbase_id$id[i]
        print(paste0(i, ": ", mediator_id))
        mediation_MR(exp_name = exposure_name, out_name = outcome_name,
                     mediation_name = mrbase_id$trait[i], id_exposure = exposure_id,
                     id_mediation = mediator_id, id_outcome = outcome_id,
                     clump_kb = clump_kb, clump_r2 = clump_r2, pval = pval,
                     action = action, proxies_exp = proxies_exp,
                     proxies_med = proxies_med, proxies_out = proxies_out,
                     palindromic = palindromic, remove_HLA = remove_HLA,
                     two_step_MR = two_step_MR, remove_exp_med_intersect = remove_exp_med_intersect,
                     mvmr_mediation = mvmr_mediation, mvmr_med_sig = mvmr_med_sig,
                     mvmr_exp_med_sig_all = mvmr_exp_med_sig_all,
                     save_path = save_path)
      }
    }
    else if (!is.null(outcome_file)) {
      for (i in num1:num2) {
        mediator_id <- mrbase_id$id[i]
        print(paste0(i, ": ", mediator_id))
        mediation_MR(exp_name = exposure_name, out_name = outcome_name,
                     mediation_name = mrbase_id$trait[i], id_exposure = exposure_id,
                     id_mediation = mediator_id, out_full = outcome_file,
                     clump_kb = clump_kb, clump_r2 = clump_r2, pval = pval,
                     action = action, proxies_exp = proxies_exp,
                     proxies_med = proxies_med, proxies_out = proxies_out,
                     palindromic = palindromic, remove_HLA = remove_HLA,
                     two_step_MR = two_step_MR, remove_exp_med_intersect = remove_exp_med_intersect,
                     mvmr_mediation = mvmr_mediation, mvmr_med_sig = mvmr_med_sig,
                     mvmr_exp_med_sig_all = mvmr_exp_med_sig_all,
                     save_path = save_path)
      }
    }
  }
  else if (!is.null(exposure_file)) {
    if (!is.null(outcome_id)) {
      mrbase_id <- subset(mrbase_id, mrbase_id$id != outcome_id)
      for (i in num1:num2) {
        mediator_id <- mrbase_id$id[i]
        print(paste0(i, ": ", mediator_id))
        mediation_MR(exp_name = exposure_name, out_name = outcome_name,
                     mediation_name = mrbase_id$trait[i], exp_full = exposure_file,
                     id_mediation = mediator_id, id_outcome = outcome_id,
                     clump_kb = clump_kb, clump_r2 = clump_r2, pval = pval,
                     action = action, proxies_exp = proxies_exp,
                     proxies_med = proxies_med, proxies_out = proxies_out,
                     palindromic = palindromic, remove_HLA = remove_HLA,
                     two_step_MR = two_step_MR, remove_exp_med_intersect = remove_exp_med_intersect,
                     mvmr_mediation = mvmr_mediation, mvmr_med_sig = mvmr_med_sig,
                     mvmr_exp_med_sig_all = mvmr_exp_med_sig_all,
                     save_path = save_path)
      }
    }
    else if (!is.null(outcome_file)) {
      for (i in num1:num2) {
        mediator_id <- mrbase_id$id[i]
        print(paste0(i, ": ", mediator_id))
        mediation_MR(exp_name = exposure_name, out_name = outcome_name,
                     mediation_name = mrbase_id$trait[i], exp_full = exposure_file,
                     id_mediation = mediator_id, out_full = outcome_file,
                     clump_kb = clump_kb, clump_r2 = clump_r2, pval = pval,
                     action = action, proxies_exp = proxies_exp,
                     proxies_med = proxies_med, proxies_out = proxies_out,
                     palindromic = palindromic, remove_HLA = remove_HLA,
                     two_step_MR = two_step_MR, remove_exp_med_intersect = remove_exp_med_intersect,
                     mvmr_mediation = mvmr_mediation, mvmr_med_sig = mvmr_med_sig,
                     mvmr_exp_med_sig_all = mvmr_exp_med_sig_all,
                     save_path = save_path)
      }
    }
  }
}


