
get.gsea.gp<-function(gp.p,sig.gp,sig.gene){
        gm=unlist(strsplit(names(gp.p),split="_"),use.names=FALSE)
        all.gene<-unique(gm)

        #significant gene groups and significant genes
#        if(!is.null(nTopPG)){
#                sig.gp=names(gp.p)[order(gp.p,decreasing = FALSE)][1:nTopPG]}else{
#                sig.gp<-names(gp.p)[p.adjust(gp.p,method="fdr")<FDRthre]}
#        sig.gene<-unique(unlist(strsplit(sig.gp,split="_"),use.names=FALSE))

        n.sig.gp<-length(sig.gp)
        n.gp<-length(gp.p)

        # pwf for gene group
        gene.freq<-c(table(gm)[all.gene])
        gene.prob<-dhyper(0,gene.freq,n.gp-gene.freq,n.sig.gp)
	gene.prob[gene.prob==0]=min(gene.prob[!gene.prob==0])
        data.frame(entrezid=all.gene,pwf=unname(gene.prob),score=ifelse(all.gene %in% sig.gene,1,0))
}



# a Monte-Carlo GSEA function
p.mc<-function(gs,gm.mc,gm.mc.cs,sig.idx,all.gene,mc.sig){
        gs.idx=which(all.gene %in% gs)
        stat<-sum(gs.idx %in% sig.idx)/length(sig.idx)
        gs2=rep(0,nrow(gm.mc));gs2[gs.idx]=1
        temp=as.vector(gs2 %*% gm.mc)
        stat.null<-"if"(mc.sig,temp/gm.mc.cs,(length(gs)-temp)/(length(all.gene)-gm.mc.cs))
        return(sum(stat.null>=stat)/ncol(gm.mc))
}



#main function
gsPG<-function(probe.p,Data4Cor=NULL,FDRthre=0.05,nTopPG=NULL,MonteCarlo=FALSE,
		 GeneProbeTable=NULL,arrayType=NULL,gSetName=NULL,geneSet=NULL,
		 species="Human",combpMethod="fisher",combpAdjust = "nyholt",
		 outfile="gsPG",ncore=5){

#gene-probe flat table
if(is.null(Data4Cor)){combpAdjust = "none"}
if(is.null(GeneProbeTable)){
        if(!(arrayType %in% c("450K","EPIC"))){stop("Specify arrayType as 450K or EPIC, or providing GeneProbeTable")}
        GeneProbeTable<-getIlluminaAnnotation(arrayType=arrayType)
}

names(probe.p)[which(names(probe.p)=="probe")]="Name"
probe.p$Name=as.character(as.vector(probe.p$Name))

if(combpAdjust=="none"){Data4Cor=NULL}else{
		   tmp=sum(!probe.p$Name %in% rownames(Data4Cor))
		   if(tmp>0){
			   probe.p=probe.p[probe.p$Name %in% rownames(Data4Cor),]
			   message(paste0("Warning: ",tmp," probes were missing in Data4Cor, and thus were excluded from analysis"))
		   }
	}

GeneProbeTable$p=probe.p$p[match(GeneProbeTable$Name, probe.p$Name)]
GeneProbeTable<-GeneProbeTable[complete.cases(GeneProbeTable),]
GeneProbeTable<-GeneProbeTable[order(GeneProbeTable$entrezid,GeneProbeTable$p),]
GeneProbeTable$Name=as.character(as.vector(GeneProbeTable$Name))
GeneProbeTable$entrezid=as.character(as.vector(GeneProbeTable$entrezid))
GeneProbeTable$p=as.numeric(as.vector(GeneProbeTable$p))

#get probe group p value
cpg.gp<-aggregate(entrezid~Name,data=GeneProbeTable,FUN=function(x){paste(sort(x),collapse="_")})
cpg.gp$p<-probe.p$p[match(cpg.gp$Name,probe.p$Name)]

genelist<-unique(cpg.gp$entrezid)
genep=mclapply(genelist,combinep,Data4Cor=Data4Cor,GeneProbeTablep=cpg.gp,
	       combpMethod=combpMethod,combpAdjust=combpAdjust,mc.cores=ncore)
genep=unlist(genep)
names(genep)=genelist
genep=genep[!is.na(genep)]

if(is.null(nTopPG)){
	sig.gp<-names(genep)[p.adjust(genep,method="fdr")<FDRthre]
}else{sig.gp=names(genep)[order(genep,decreasing = FALSE)][1:nTopPG]}
sig.gene<-unique(unlist(strsplit(sig.gp,split="_"),use.names=FALSE))
if(length(sig.gene)<10){stop("Number of significant gene is less than 10, 
			     please adjust threshold or use ranking method in gsGene()")}
sgName=suppressMessages(geneID2geneName(sig.gene,species))

gSetName12=ParseGsetName(geneSet,gSetName)
gSetName1=gSetName12$gSetName1
gSetName2=gSetName12$gSetName2

for(gset1 in gSetName1){
	for(gset2 in gSetName2){
	if(is.null(geneSet)){
	if(gset1=="KEGG"){geneSet=getKEGG(species)
	}else if(gset1=="GO"){geneSet=getGO(gset2,species)
	}else if(gset1=="MSigDB"){geneSet=getMSigDB(gset2,species)
	}else if(gset1=="Reactome"){geneSet=getReactome(species)
	}else{stop("Please specify gSetName as KEGG, GO, MSigDB, or Reactome, or provide a geneSet list")}
	}
# Check geneSet, remove genes that are not in test universe
universe<-unique(GeneProbeTable$entrezid)
if(!is.list(geneSet))stop("geneSet should be a list object")
geneSet <- lapply(geneSet,function(x){x=as.character(as.vector(x)); x=x[!is.na(x)];x=x[x %in% universe];unique(x)})
inUniv <- sapply(geneSet, function(x) length(x) > 0)
geneSet <- geneSet[inUniv]

if(MonteCarlo){
#	if(!requireNamespace("dqrng", quietly = TRUE)){stop("R package dqrng required, but not installed")}
	ss<-1e5
	n.sig.gp<-length(sig.gp)
	n.gp<-length(genep)
	gm<-strsplit(names(genep),split="_")
	all.gene=unique(unlist(gm,use.names=FALSE))
	sig.idx.gp<-which(all.gene %in% sig.gene)
	if(n.sig.gp<n.gp/2){
		mc.sig=TRUE
		gm.mc<-mclapply(1:ss,function(i){idx<-dqrng::dqsample.int(n.gp,n.sig.gp)
                                      which(all.gene %in% unique(unlist(gm[idx],use.names=FALSE)))},mc.cores=ncore)
	}else{
		mc.sig=FALSE
		gm.mc<-mclapply(1:ss,function(i){idx<-dqrng::dqsample.int(n.gp,n.sig.gp)
                                      which(!all.gene %in% unique(unlist(gm[idx],use.names=FALSE)))},mc.cores=ncore)
	}
	gm.mc.gp<-Matrix::sparseMatrix(i=unlist(gm.mc,use.names=FALSE),j=rep(1:ss,times=sapply(gm.mc,length)),
				       x=1,dims=c(length(all.gene),ss))
	rm(gm.mc)
	gm.mc.cs.gp<-colSums(gm.mc.gp)
	p=unlist(mclapply(geneSet,p.mc,gm.mc=gm.mc.gp,gm.mc.cs=gm.mc.cs.gp,sig.idx=sig.idx.gp,all.gene=all.gene,mc.sig=mc.sig,mc.cores=ncore))
	p[p==0]=min(p[p>0])/2
	nGene=sapply(geneSet,length)
	nSigGene=sapply(geneSet,function(x)sum(sig.gene %in% x))
	Sig.Gene=sapply(geneSet,function(x)paste(sgName[sig.gene[sig.gene %in% x]],collapse=";"))
	resu=data.frame(ID=names(geneSet),nGene=nGene,nSigGene=nSigGene,p=p,FDR=p.adjust(p,method="fdr"),Sig.Gene=Sig.Gene)
}else{
gsea.gp<-get.gsea.gp(genep,sig.gp,sig.gene)
resu1=mclapply(geneSet,function(x){p.gs(x,gsea.gp)},mc.cores=ncore)
resu <- do.call(rbind, lapply(resu1, unlist))
resu=data.frame(ID=names(geneSet),resu)
names(resu)=c("ID","nGene","nSigGene","p")
resu$Sig.Gene=sapply(geneSet,function(x)paste(sgName[sig.gene[sig.gene %in% x]],collapse=";"))
resu$FDR=p.adjust(resu$p,method="fdr")
}
#resu=resu[resu$Sig.Gene>0,]
outfn=paste0(outfile,"_",gset1,"_",gset2,".csv")
write.csv(resu,outfn,row.names=FALSE)
}}}



