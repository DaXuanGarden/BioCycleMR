rm(list=ls())
gc()
getwd()
#setwd("D:\文档\GitHub\BioCycleMR")
# Load the devtools package
library(usethis)
library(utils)
library(devtools)
# install.packages(c("devtools", "roxygen2", "usethis", "testthat"))
#detach(package:DXMarkers)
# rm(list = c("annotate_markers"))
# 在R控制台中运行以下命令
library(roxygen2)
roxygen2::roxygenize()
library(remotes)
install_github("MRCIEU/TwoSampleMR")
devtools::install_github("renkun-ken/rlist")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
devtools::install_github("MRCIEU/gwasglue")


devtools::document()
# Check the package
check()


# 假设你的R包项目目录是"DXMarkers"，请将路径替换为你的实际路径
getwd()
#setwd("/home/data/t050446/DX_Package/DXMarkers")
# 使用devtools包进行打包
devtools::build("D:/文档/GitHub/BioCycleMR")
library(devtools)
devtools::install_local("/home/data/t050446/DX_Package/DXMarkers_1.0.tar.gz", dependencies = TRUE, upgrade = FALSE)

library(DXMarkers)
setwd("/home/data/t050446/01 Single Cell Project/MR&ScRNA")

DXMarkers_heatmap(scRNA, DXMarkers_result, "00DXMarkers.png")

DXMarkers_dotplots(scRNA, "SCINA_markers.csv", "DXMarkers")

DXMarkers_dotplots(scRNA, "~/SCINA_markers.csv", "~/DXMarkers")



#######生成并查看说明书#####
# 使用usethis包添加一个名为"my-vignette"的vignette
usethis::use_vignette("DXMarkers")

# Load devtools package
library(devtools)

# Build your vignettes
build_vignettes()

# 查看使用说明书
browseVignettes("DXMarkers")

