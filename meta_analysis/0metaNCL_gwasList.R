## to generate the list for looping

require("data.table")
root.dir <- "/home/nwl15/WORKING_DATA/HRCimpvData/"
sum.dir <- "ResultSummary"

res.dir <- c(paste0(root.dir, "NCL1_2/", sum.dir),
             paste0(root.dir, "NCL3/", sum.dir),
             paste0(root.dir, "NCL4_GT/", sum.dir)
)

dataprefix <- c("NCL1_2_finalQCed",
                "NCL3_finalQCed",
                "NCL4_gtQCed")


write.table(data.table(res.dir, dataprefix), file=paste0("NCL1_4_folder.lst"), quote=FALSE,
            sep="\t", row.names=FALSE, col.names=FALSE
)
