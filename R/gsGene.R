#es: Enrichment Score;gs: a gene set; rgl:ranked gene test result
es<-function(gsi,rgl){
        ng=length(rgl)
        idx=names(rgl) %in% gsi
        x=rep(0,ng);x[idx]=1;y=1-x
        cumsum(x*rgl/sum(rgl[idx])-y/(ng-length(gsi)))
}

fitBeta<-function(i,dat){
        m <- mean(dat[,i])
        v <- var(dat[,i])
        shape1_mom <- m * ((m * (1 - m) / v) - 1)
        shape2_mom <- (1 - m) * ((m * (1 - m) / v) - 1)
        return(c(shape1_mom,shape2_mom))
}

permuTest2<-function(gSet,rgl,nSim=10000,es0,ncore){
        gsSize0=sapply(gSet,length)
        ng=length(rgl)
        gsSize=unique(gsSize0)
        maxGsSize=max(gsSize)
        resu=mclapply(1:nSim,function(k){
                rgs=sample.int(ng,maxGsSize)
                sapply(gsSize,function(gsSize_i){
                 rgs.i<-sort(rgs[1:gsSize_i])
                 rglsumi=sum(rgl[rgs.i])
                 ngi=ng-gsSize_i
                 x<-rgl[rgs.i]
                 y<-rgs.i-c(1:gsSize_i)
                 max(cumsum(x)/rglsumi-y/ngi)
                 })
        },mc.cores=ncore)
        mat.null=do.call(rbind,lapply(resu,unlist))

#get permutation P
        colnames(mat.null)=gsSize
        gsSize0=as.character(as.vector(gsSize0))
        permuP=sapply(1:length(es0),function(i){
               mean(mat.null[,gsSize0[i]]>= es0[i])})
        permuP[permuP==0]=1/(nSim*2)

#get betafit P
        b_para=lapply(1:ncol(mat.null),fitBeta,dat=mat.null)
        b_para=do.call(rbind,lapply(b_para,unlist))
        rownames(b_para)=gsSize
        b_para=b_para[as.character(as.vector(gsSize0)),]
        p=sapply(1:length(es0),function(i){
                exp(pbeta(es0[i], shape1 = b_para[i,1], shape2 =b_para[i,2],log.p =TRUE, lower.tail = FALSE))})
        return(list(p=p,permuP=permuP))
}


p.gs.gene<-function(gs,sig.gene,nGene,nSig){
	gs.ng=length(gs)
	gs.nSig=sum(gs %in% sig.gene)-1
	phyper(gs.nSig,gs.ng,nGene-gs.ng,nSig,lower.tail=FALSE)
}


#main function
gsGene<-function(probe.p,Data4Cor=NULL,method="Threshold",FDRthre=0.05,nTopGene=NULL,
		   GeneProbeTable=NULL,arrayType=NULL,gSetName="KEGG",geneSet=NULL,
		   species="Human",combpMethod="fisher",combpAdjust = "nyholt",
		   outfile="gsGene",outGenep=FALSE,gseaParam=1,nperm=1e4,ncore=5){

if(is.null(Data4Cor)){combpAdjust = "none"}
if(is.null(GeneProbeTable)){
        if(!(arrayType %in% c("450K","EPIC"))){stop("Specify arrayType as 450K or EPIC, or providing GeneProbeTable")}
        GeneProbeTable<-getIlluminaAnnotation(arrayType=arrayType)
}
names(probe.p)[which(names(probe.p)=="probe")]="Name"
names(probe.p)[which(names(probe.p)=="P")]="p"
probe.p$Name=as.character(as.vector(probe.p$Name))

if(method=="Ranking"){nTopGene=NULL}
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

#get genep
genelist<-unique(GeneProbeTable$entrezid)
genep=mclapply(genelist,combinep,Data4Cor=Data4Cor,GeneProbeTablep=GeneProbeTable,
	       combpMethod=combpMethod,combpAdjust=combpAdjust,mc.cores=ncore)
genep=unlist(genep)
names(genep)=genelist
genep=genep[!is.na(genep)]
if(outGenep){
        tmpgn=suppressMessages(geneID2geneName(names(genep),species))
        tmp=data.frame(geneid=names(genep),gene=tmpgn,p=genep)
        write.csv(tmp,"gsGene_genep.csv",row.names=FALSE)
}

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

if(method=="Ranking"){
	genep[genep==0]=min(genep[!genep==0])
	genep[genep==1]=max(genep[!genep==1])
	genep=sort(-log(genep),decreasing=TRUE)
	genep=genep^gseaParam
        es0=sapply(geneSet,function(gsi)max(es(gsi,genep)))
	p=permuTest2(gSet=geneSet,rgl=genep,nSim=nperm,es0=es0,ncore=ncore)
	nGene=sapply(geneSet,length)
	resu=data.frame(ID=names(geneSet),nGene=nGene,EnrichScore=es0,p=p$p,permuP=p$permuP)
	resu$FDR=p.adjust(resu$p,method="fdr")
}else if(method=="Threshold"){
	all.gene<-unique(names(genep))
        if(!is.null(nTopGene)){
                sig.gene=names(genep)[order(genep,decreasing = FALSE)][1:nTopGene]
	}else{
                sig.gene<-names(genep)[p.adjust(genep,method="fdr")<FDRthre]
	}
	if(length(sig.gene)<10){stop("Number of significant gene is less than 10, 
				     please adjust threshold or use ranking method in gsGene()")}
	sgName=suppressMessages(geneID2geneName(sig.gene,species))
	nSigGene=length(sig.gene)
	nGene=length(all.gene)
	p=sapply(geneSet,p.gs.gene,sig.gene=sig.gene,nGene=nGene,nSig=nSigGene)
        nGene=sapply(geneSet,length)
        nSigGene=sapply(geneSet,function(x)sum(sig.gene %in% x))
        Sig.Gene=sapply(geneSet,function(x)paste(sgName[sig.gene[sig.gene %in% x]],collapse=";"))
        resu=data.frame(ID=names(geneSet),nGene=nGene,nSigGene=nSigGene,p=p,FDR=p.adjust(p,method="fdr"),Sig.Gene=Sig.Gene)
#	resu=resu[resu$Sig.Gene>0,]
}else{stop("method need to be Threshold or Ranking")}

outfn=paste0(outfile,"_",gset1,"_",gset2,".csv")
write.csv(resu,outfn,row.names=FALSE)
}}}



