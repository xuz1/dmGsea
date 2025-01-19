#es: Enrichment Score;gs: a gene set; rgl:ranked gene test result
es<-function(gsi,rgl){
        ng=length(rgl)
        idx=names(rgl) %in% gsi
        x=rep(0,ng);x[idx]=1;y=1-x
        cumsum(x*rgl/sum(rgl[idx])-y/(ng-length(gsi)))
}

fitBeta<-function(i,dat){
	dd=dat[,i]
        m <- mean(dd)
        v <- var(dd)
        shape1 <- m * ((m * (1 - m) / v) - 1)
        shape2 <- (1 - m) * ((m * (1 - m) / v) - 1)
        return(c(shape1,shape2))
}

fitBeta2<-function(i,dat){
        dd=dat[,i]
#dat1=dat2=abs(dd)
	dat1=dd[dd>0]
        dat2=abs(dd[dd<=0])
	#approximate dist if size is too small
	if(length(dat1)<10)dat1=abs(dat)
	if(length(dat2)<10)dat2=abs(dat)
        m <- mean(dat1)
        v <- var(dat1)
        shape1_pos <- m * ((m * (1 - m) / v) - 1)
        shape2_pos <- (1 - m) * ((m * (1 - m) / v) - 1)
        m <- mean(dat2)
        v <- var(dat2)
        shape1_neg <- m * ((m * (1 - m) / v) - 1)
        shape2_neg <- (1 - m) * ((m * (1 - m) / v) - 1)
        return(c(shape1_pos,shape2_pos,shape1_neg,shape2_neg))
}

betaFitPermu<-function(gSet,rgl,nSim=10000,es0,ncore,scoreType){
        gsSize0=sapply(gSet,length)
        ng=length(rgl)
if(1 %in% gsSize0){
	idx1=which(gsSize0==1)
	idx2=which(gsSize0>1)
	es00=es0
	es1=es0[idx1]
	es0=es0[idx2]
	gsSize0=gsSize0[idx2]
	null.p<-1-(seq_along(rgl)-1)/(length(rgl)-1)
	null.n<-(-(seq_along(rgl)-1)/(length(rgl)-1))
	null.s<-ifelse(null.p>abs(null.n),null.p,null.n)
	null.s.p=matrix(null.s[null.s>0],ncol=1)
	null.s.n=matrix(null.s[null.s<=0],ncol=1)
	if(scoreType=="std"){
		p1=sapply(es1,function(x){
				  if(x>0){mean(null.s.p>=x)}else{mean(null.s.n<=x)}})
	}else if(scoreType=="pos"){
		p1=sapply(es1,function(x)mean(null.p>=x))
	}else if(scoreType=="neg"){
		 p1=sapply(es1,function(x)mean(null.n<=x))
	}else{stop("scoreType should be std, pos or neg")}
}
	
	gsSize=unique(gsSize0)
        maxGsSize=max(gsSize)
        resu=mclapply(1:nSim,function(k){
                rgs=sample.int(ng,maxGsSize)
                sapply(gsSize,function(gsSize_i){
                 rgs.i<-sort(rgs[1:gsSize_i])
                 ngi=ng-gsSize_i
                 x<-rgl[rgs.i]
                 sumx=sum(x)
                 y<-rgs.i-c(1:gsSize_i)
		 tmp=cumsum(x)/sumx-y/ngi
		 es.p=max(tmp)
		 es.n=min(c(0,tmp-x/sum(x)))
		 es.s=ifelse(es.p>abs(es.n),es.p,ifelse(es.p!=abs(es.n),es.n,0))
		 switch(scoreType,"std"=es.s,"pos"=es.p,"neg"=es.n)
                 })
        },mc.cores=ncore)
        mat.null=do.call(rbind,lapply(resu,unlist))
	rm(resu)

#get permutation P	
        colnames(mat.null)=gsSize
        gsSize0=as.character(as.vector(gsSize0))
        if(scoreType == "pos"){
                permuP=sapply(1:length(es0),function(i){
                mean(mat.null[,gsSize0[i]]>= es0[i])})
        }else if(scoreType == "neg"){
                permuP=sapply(1:length(es0),function(i){
                mean(mat.null[,gsSize0[i]]<= es0[i])})
	}else if(scoreType == "std"){
                permuP=sapply(1:length(es0),function(i){
                        null=mat.null[,gsSize0[i]]
                        if(es0[i]>0){
                                null=null[null>0]
				mean(null>=es0[i])
#				(sum(null>=es0[i])+1)/(length(null)+1)
                        }else{
                                null=null[null<=0]
				mean(null<=es0[i])
#				(sum(null<=es0[i])+1)/(length(null)+1)
                        }
                })
        }else{stop("scoreType should be pos, neg or std")}
        permuP[permuP==0]=1/(nSim*2)

#get betafit P
        if(scoreType=="std"){
		b_para=lapply(1:ncol(mat.null),fitBeta2,dat=mat.null)
	        b_para=do.call(rbind,lapply(b_para,unlist))
	        rownames(b_para)=gsSize

	        b_para=b_para[as.character(as.vector(gsSize0)),]
	        p=sapply(1:length(es0),function(i){
			 if(es0[i]>0){
                exp(stats::pbeta(es0[i], shape1 = b_para[i,1], shape2 =b_para[i,2],log.p =TRUE, lower.tail = FALSE))
			 }else{
                exp(stats::pbeta(abs(es0[i]), shape1 = b_para[i,3], shape2 =b_para[i,4],log.p =TRUE, lower.tail = FALSE))
			 }
		})
	}else if(scoreType %in% c("pos","neg")){
		mat.null=abs(mat.null)
		b_para=lapply(1:ncol(mat.null),fitBeta,dat=mat.null)
	        b_para=do.call(rbind,lapply(b_para,unlist))
	        rownames(b_para)=gsSize
	        b_para=b_para[as.character(as.vector(gsSize0)),]
		es0=abs(es0)
	        p=sapply(1:length(es0),function(i){
                exp(stats::pbeta(es0[i], shape1 = b_para[i,1], shape2 =b_para[i,2],log.p =TRUE, lower.tail = FALSE))
			 })
	}else{stop("scoreType should be pos, neg or std")}
	p[p==0]=min(p[p!=0])

	if(scoreType=="std"){
        p[p==0]=min(p[!p==0])
        p[p==1]=max(p[!p==1])
        permuP[permuP==0]=min(permuP[!permuP==0])
        permuP[permuP==1]=max(permuP[!permuP==1])
     
	logP=-log10(p)
	pfit=data.frame(p=p,permuP=permuP,es0=es0)
	pfit=na.omit(pfit)
	pfit=pfit[pfit$permuP>=1/nSim,]
#	pfit=pfit[pfit$p>=1/nSim,]
	pfit$p=-log10(pfit$p)
	pfit$permuP=-log10(pfit$permuP)

	permuPx=pfit$permuP[pfit$es0>0]
	px=pfit$p[pfit$es0>0]
	model=lm(permuPx~px-1)
#	wt <- 1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
#	model=lm(permuPx~px-1,weights=wt)
	p[es0>0]=predict(model,data.frame(px=logP[es0>0]))

	permuPx=pfit$permuP[pfit$es0<=0]
	px=pfit$p[pfit$es0<=0]
	model=lm(permuPx~px-1)
#	wt <- 1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
#	model=lm(permuPx~px-1,weights=wt)
	p[es0<=0]=predict(model,data.frame(px=logP[es0<=0]))

	p=1/10^p
	}

	if(exists("p1")){
		p2=p
		p=rep(NA,length(es00))
		p[idx1]=p1
		p[idx2]=p2
		p2=permuP
		permuP=rep(NA,length(es00))
		permuP[idx1]=p1
		permuP[idx2]=p2
	}
        return(data.frame(p=p,permuP=permuP))
}


#main function
gsRank<-function(stats,outfile="gsRank",scoreType="std",gSetName=NULL,
		   geneSet=NULL,gseaParam=1,species="Human",nperm=1e4,ncore=5){

	genep=stats
	genep=genep[!is.na(genep)]
	genep=sort(genep,decreasing=TRUE)
	genep=abs(genep)
	genep=genep^gseaParam

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
	}else{stop("Please specify gSetName as KEGG, GO, MSigDB, 
		   or Reactome, or provide a geneSet list")}
	}

	# Check geneSet, remove genes that are not in test universe
	universe<-names(genep)
	if(!is.list(geneSet))stop("geneSet should be a list object")
	geneSet <- lapply(geneSet,function(x){x=as.character(as.vector(x)) 
			  x=x[!is.na(x)];x=x[x %in% universe];unique(x)})
	inUniv <- sapply(geneSet, function(x) length(x) > 0)
	geneSet <- geneSet[inUniv]

        es0=sapply(geneSet,function(gsi){
			   es1=es(gsi,genep)
			   switch(scoreType,"std"={
   				 m1 <- max(es1); m2 <- min(es1)
				 ifelse(m1 > abs(m2), m1, ifelse(m1 == abs(m2), 0, m2))
				  },"pos"=max(es1),"neg"=min(es1))})
	p=betaFitPermu(gSet=geneSet,rgl=genep,nSim=nperm,es0=es0,ncore=ncore,scoreType=scoreType)
	nGene=sapply(geneSet,length)
	resu=data.frame(ID=names(geneSet),nGene=nGene,EnrichScore=es0,p=p$p,permuP=p$permuP)
	resu$FDR=p.adjust(resu$p,method="fdr")

	outfn=paste0(outfile,"_",gset1,"_",gset2,".csv")
	write.csv(resu,outfn,row.names=FALSE)
}}}


