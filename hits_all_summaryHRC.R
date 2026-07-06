args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
} else {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

dataprefix

require("data.table")
rds.list <- list.files(path="./", pattern="_hitsALL.rds$")
kkpattern <- paste0(paste0(dataprefix, "_"), "|_hitsALL.rds", collapse="")
#phenos <- gsub("NCL1_2_finalQCed_|_hitsALL.rds", "", rds.list)
phenos <- gsub(kkpattern, "", rds.list)
phenos.res <- vector("list", length= length(phenos))
names(phenos.res) <- phenos
for (n in seq_along(phenos)){
    #phenos.res[[n]] <- readRDS(paste0("NCL1_2_finalQCed_", phenos[n], "_hitsALL.rds"))
    phenos.res[[n]] <- readRDS(paste0(dataprefix, "_", phenos[n], "_hitsALL.rds"))
}

phenos.res <- rbindlist(phenos.res, use.name=TRUE, fill=TRUE, idcol=TRUE)
#phenos.res1 <- rbindlist(phenos.res, use.name=TRUE, fill=TRUE, idcol=TRUE)
write.table(phenos.res, file=paste0(dataprefix, "_hits_10^5_HRC.txt"), sep="\t",
            quote=FALSE, col.names=TRUE, row.names=FALSE)

## get rid of mulitple allelic variants
if ( phenos.res[, .N, by=rsid][N!=10, .N] !=0 ){
    cat(phenos.res[, .N, by=rsid][N!=10, rsid], file=paste0(dataprefix, "_", "MultipleAllelic.lst"), sep="\n")
}

phenos.res1 <- phenos.res[!(rsid %chin% phenos.res[, .N, by=rsid][N!=10, rsid]), ]

chr.uni <- phenos.res[, unique(chromosome)]
chr.dt <- data.table(rsid=rep(paste0("CHROM ", chr.uni), times=1, each=10),
                     chromosome=as.numeric(rep(paste0(chr.uni), times=1, each=10)),
                     position=0,
                     .id=rep(c("status", "Normal", "t15_17", "Complex", "CBF", "del5_7", "Trans", "Trisomies", "Any.Monosomy", "monosomal.karyotype"), times=length(chr.uni))
                     )
phenos.res1 <- rbindlist(list(chr.dt, phenos.res1), use.names=TRUE, fill=TRUE)

all.markers <- unique(phenos.res1[, .(rsid, chromosome, position)])
setorder(all.markers, chromosome, position)[, numRec:=1:.N][, batch:=(numRec+49) %/% 50]


data.plot <- phenos.res1[all.markers, batch:=i.batch, on=c(rsid="rsid") ]
setnames(data.plot, c(".id"), c("pheno"))
data.plot[, pheno:=factor(pheno, levels=c("status", "Normal", "t15_17", "Complex", "CBF", "del5_7", "Trans", "Trisomies", "Any.Monosomy", "monosomal.karyotype"))][, rsid:=factor(rsid, levels=all.markers[,rsid])]
setorder(data.plot, rsid, pheno)
#data.plot <- data.inReg[,.SD[which.max(logP)],by=c("region", "variable")]

require("ggplot2")
require("grid")
require("gridExtra")

batch.order <- data.plot[,unique(batch)]

switch(dataprefix,
       NCL1_2_finalQCed={main.title <- paste0("Newcastle1_2", " HRC")},
       NCL3_finalQCed={main.title <- paste0("Newcastle3", " HRC")},
       NCL4_gtQCed={main.title <- paste0("Newcastle4", " HRC")},
       NCL5_finalQCed_1821={main.title <- paste0("Newcastle5_1821", " HRC")},
       NCL5_finalQced4imp={main.title <- paste0("Newcastle5_ukb500k", " HRC")}
       )
main.title
#main.title <- "Newcastle1_2"

base_size <- 9

pdf(paste0(dataprefix, "_Res_summary_hits_HRCdata",".pdf"), width=12, height=8)
i = 1
plot = list()
for (n in seq_along(batch.order)){
    ## subset batch
    plotdsn <- data.plot[batch==batch.order[n], .(rsid, pheno, logP)]
    x.labels <- copy(sort(unique(plotdsn[, rsid])))
    if (plotdsn[, .N] !=500){
        no.empty <- ((500-plotdsn[,.N])/10) -1 
        last.rsid <- plotdsn[,max(as.numeric(rsid))] + 1
        faked.kk <- data.table(rsid=rep(seq(from=last.rsid, to=last.rsid+no.empty ), times=1, each=10),
                               pheno=rep(c("status", "Normal", "t15_17", "Complex", "CBF", "del5_7", "Trans", "Trisomies", "Any.Monosomy", "monosomal.karyotype"),
                                         times=(no.empty+1))
                               )
        faked.kk[, logP:=NA]
        plotdsn <- rbindlist(list(plotdsn, faked.kk))
        empty.str <- rep("", (no.empty+1))
        x.breaks <-  as.character(copy(sort(unique(plotdsn[, rsid]))))
        x.labels <- c(as.character(x.labels), empty.str)
    }
    x.breaks <- as.character(x.labels)
    x.labels <- as.character(x.labels)

    
    ## check if calls are not made for batch
### process data for plotting here ####
    plot[[i]] = ggplot(plotdsn, aes(rsid, pheno)) + geom_tile(aes(fill = logP), colour="white") +
                     scale_fill_gradient(low = "yellow", high = "red", na.value="white") +
        geom_text(aes(fill=logP, label=round(logP, 1)), size=2) +
                                theme_grey(base_size = base_size) +
                            labs(x = "", y = "") +
                                scale_x_discrete(breaks=x.breaks, labels=x.labels, expand = c(0, 0)) +
                                    scale_y_discrete(expand = c(0, 0)) + theme(#legend.position = "none",
                    # axis.ticks = element_blank(),
                                                         axis.text.x = element_text(size = base_size *0.7,
                                                             angle = 90, hjust = 0.5, colour = "black"))
                      
     #axis.line = element_line(colour = "black") 
    if (i %% 3 == 0) { ## print 6 plots on a page
        argstest.list <- c(plot,list(nrow=3,ncol=1,top=textGrob(paste0(main.title), gp=gpar(cex=1.2), just="top")))
        print (do.call(grid.arrange,  argstest.list))
        plot = list() # reset plot 
        i = 0 # reset index
    }
    i = i + 1
}
if (length(plot) != 0) { 
    argstest.list <- c(plot, list(nrow=3, ncol=1, top=textGrob(paste0(main.title), gp=gpar(cex=1.2), just="top")))
    print (do.call(grid.arrange,  argstest.list))
}
dev.off()
