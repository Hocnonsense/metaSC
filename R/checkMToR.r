###
#* @Date: 2020-10-05 16:27:53
#* @LastEditors: Hwrn
#* @LastEditTime: 2020-10-05 16:28:17
#* @FilePath: /HScripts/R/checkMToR.r
#* @Description: help to binning with checkm
###
# ---------------------------------------------------------- If you are initing, go from HERE @Hwrn
library(mmgenome2)
setwd("D:/Files/TODU/2020-09-MgAffact")
cov_name = "Analyze/origin-depth._cov"
cov = read.table(cov_name, header = T)

# ------------------------------------------------------------- change here for several rolls @Hwrn
bin_dir = "F-06-MAG/03_modify/5_m5/modify/"
mks_dir = "F-06-MAG/03_modify/5_m5/modify/"
out_dir = "F-06-MAG/03_modify/5_m5/rmmed/"

bin = ""
# --------------------------------------------------------------- dealing with given bin HERE @Hwrn
bin = "maxbin2_107.029_sub112"
bin = "maxbin2_107.115_sub3"
bin = "maxbin2_107.047_sub11"
i = 0
bin_name <- paste0(bin_dir, bin, ".fa")
bin_name = "F-06-MAG/03_modify/5_m5/rmmed/maxbin2_107.029_sub1121.fa"
marksets <- paste0(mks_dir, bin, ".marker_gene.csv")  # ! output by checkMarkToR.py
rm("assembly")
mm <- mmload(
    assembly = bin_name,
    coverage = cov,
    additional = read.csv(marksets)
)
# kmer_pca =TRUE,
# verbose = TRUE,  default
# kmer_size = 4L,  default
# kmer_BH_tSNE = TRUE,  加上这一个就跑不通

# ------------------------ To choose several different subbins HERE | choose another out_name @Hwrn
i = i + 1
# select the bins, colored_by can in [cluster_id, x0-9]
se = mmplot(mm,
            y = "cov_totalAvgDepth",
            x = "gc",
            color_by = 'cluster_id',
            min_length = 1000,
            locator=TRUE,
            color_vector = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0))
)$selection

# ---------------------------------------------------------------- continue and save .fa HERE @Hwrn
bin1=mmextract(mm,selection = se)
out_name <- paste0(out_dir, bin, i, ".fa")
print(out_name)
mmexport(bin1,assembly = assembly,file = out_name)
# remind what you choose
