# make SL table
# make expression table

library(data.table)
library(VennDiagram)
library(pbapply)

load("kw.RData")
# patients with both exp and str: exp, repeats

minpts=5
nexp=ncol(exp);nrepeats=ncol(repeats);pvalues = c()
headers <- data.frame(coor=character(),gene=character(), region=character(),repeats=character(),isoform=character(),stringsAsFactors=FALSE)
### get p and statistics: case 1,2,3
#r square

# ss vs ll
group <- c('SS','LL')
association <- pbapply(repeats,1,function(x) {
    repeattype=as.factor(x[-c(1:4)])
    rpkm=unlist(subset(exp,TargetID==x[2]))[-c(1:2)]
    comm=intersect(names(repeattype)[!is.na(repeattype)],names(rpkm)[!is.na(rpkm)])
    repeattype=repeattype[comm];rpkm=as.numeric(rpkm[comm]);names(rpkm) <- comm
    repeattype_group <- repeattype[repeattype %in% group]
    rpkm_group <- rpkm[repeattype %in% group]
    if(sum(repeattype_group == group[1])>=minpts & sum(repeattype_group == group[2])>=minpts){
        wtest <- wilcox.test(rpkm_group ~ repeattype_group,conf.int=T)
        p <- signif(wtest$p.value,3)#,subset = repeattype %in% c(1,3)
        wc <- wtest$estimates
    }else{p <- NA; wc <- NA}
    return(list(p,wc))
})
 
# ss vs sl vs ll
group <- c('SS','SL')
association_sssl <- pbapply(repeats,1,function(x) {
    repeattype=as.factor(x[-c(1:4)])
    rpkm=unlist(subset(exp,TargetID==x[2]))[-c(1:2)]
    comm=intersect(names(repeattype)[!is.na(repeattype)],names(rpkm)[!is.na(rpkm)])
    repeattype=repeattype[comm];rpkm=as.numeric(rpkm[comm]);names(rpkm) <- comm
    repeattype_group <- repeattype[repeattype %in% group]
    rpkm_group <- rpkm[repeattype %in% group]
    if(sum(repeattype_group == group[1])>=minpts & sum(repeattype_group == group[2])>=minpts){
        wtest <- wilcox.test(rpkm_group ~ repeattype_group,conf.int=T)
        p <- signif(wtest$p.value,3)#,subset = repeattype %in% c(1,3)
        wc <- wtest$estimates
    }else{p <- NA; wc <- NA}
    return(list(p,wc))
})
association_sssl <- association_sl

group <- c('SL','LL')
association_slll <- pbapply(repeats,1,function(x) {
    repeattype=as.factor(x[-c(1:4)])
    rpkm=unlist(subset(exp,TargetID==x[2]))[-c(1:2)]
    comm=intersect(names(repeattype)[!is.na(repeattype)],names(rpkm)[!is.na(rpkm)])
    repeattype=repeattype[comm];rpkm=as.numeric(rpkm[comm]);names(rpkm) <- comm
    repeattype_group <- repeattype[repeattype %in% group]
    rpkm_group <- rpkm[repeattype %in% group]
    if(sum(repeattype_group == group[1])>=minpts & sum(repeattype_group == group[2])>=minpts){
        wtest <- wilcox.test(rpkm_group ~ repeattype_group,conf.int=T)
        p <- signif(wtest$p.value,3)#,subset = repeattype %in% c(1,3)
        wc <- wtest$estimates
    }else{p <- NA; wc <- NA}
    return(list(p,wc))
})
### case done: association

# adjust p

# add trend

# count num for sig 4 level

# intersect with Gymrek and count

# venn plot

# draw picture for top 3

# wetlab check
