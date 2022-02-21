utils::globalVariables(c('create.MK.AL_gene_buffer',
                         'G_gene_buffer_surround',
                         'LD.filter','surround.region'))

#' Data example for knockoff generation of gene buffer region
#'
#' This example dataset contains a genotype matrix of gene buffer surrounding region and the positions for gene buffer region. The example dataset is extracted from the 'plinkforGRM_1000samples_10kMarkers' plink data in the extdata folder of SAIGE package, which contains 1000 samples.
#'
#' @name Example.Knock.generation.gene.buffer
#' @docType data
#' @keywords data
#' @usage data("Knock.generation.gene.buffer.example")
#' @examples
#' data('Knock.generation.gene.buffer.example')
#' G_gene_buffer_surround=Knock.generation.gene.buffer.example$G_gene_buffer_surround
#' pos_gene_buffer=Knock.generation.gene.buffer.example$pos_gene_buffer
"Knock.generation.gene.buffer.example"

#' Data example for knockoff generation of enhancer
#'
#' This example dataset contains a genotype matrix of enhancer surrounding region and the positions for enhancer. The example dataset is extracted from the 'plinkforGRM_1000samples_10kMarkers' plink data in the extdata folder of SAIGE package, which contains 1000 samples.
#'
#' @name Example.Knock.generation.enhancer
#' @docType data
#' @keywords data
#' @usage data("Knock.generation.enhancer.example")
#' @examples
#' data('Knock.generation.enhancer.example')
#' G_enhancer_surround=Knock.generation.enhancer.example$G_enhancer_surround
#' pos_enhancer=Knock.generation.enhancer.example$pos_enhancer
"Knock.generation.enhancer.example"


#' Knockoff Generation of Gene buffer region
#'
#' This function generates multiple knockoff genotypes for gene buffer region. The knockoff generations are optimized using shrinkage leveraging algorithm.
#'
#' @param G_gene_buffer_surround The genotype matrix of the surrounding region for gene buffer region.
#' @param gene_buffer_start The start position of gene buffer region.
#' @param gene_buffer_end The end position of gene buffer region.
#' @param surround.region Surrounding region for gene buffer region, default is +-100kb.
#' @param LD.filter The correlation threshold for hierarchical clustering. Default LD filter at correlation 0.75
#' @param M Numer of multiple knockoffs.
#' @return \item{G_gene_buffer_knockoff}{A list file contains M knockoff genotypes for gene buffer region.}
#' @examples
#'library(Matrix)
#'data('Knock.generation.gene.buffer.example')
#'G_gene_buffer_surround=Matrix(Knock.generation.gene.buffer.example$G_gene_buffer_surround)
#'pos_gene_buffer=Knock.generation.gene.buffer.example$pos_gene_buffer
#'gene_buffer_start=min(pos_gene_buffer)
#'gene_buffer_end=max(pos_gene_buffer)
#'
#'G_gene_buffer_knockoff=Knockoffgeneration.gene.buffer(G_gene_buffer_surround=G_gene_buffer_surround,
#'                                                gene_buffer_start=gene_buffer_start,
#'                                                gene_buffer_end=gene_buffer_end,
#'                                                M=5,surround.region=100000,LD.filter=0.75)
#'#M=5 knockoffs
#'G_gene_buffer_knockoff1=G_gene_buffer_knockoff[1,,]
#'G_gene_buffer_knockoff2=G_gene_buffer_knockoff[2,,]
#'G_gene_buffer_knockoff3=G_gene_buffer_knockoff[3,,]
#'G_gene_buffer_knockoff4=G_gene_buffer_knockoff[4,,]
#'G_gene_buffer_knockoff5=G_gene_buffer_knockoff[5,,]
#'
#' @import SKAT
#' @import Matrix
#' @import SPAtest
#' @import CompQuadForm
#' @import irlba
#' @export
Knockoffgeneration.gene.buffer=function(G_gene_buffer_surround=G_gene_buffer_surround,
                                        gene_buffer_start=gene_buffer_start,
                                        gene_buffer_end=gene_buffer_end,
                                        M=5,surround.region=100000,LD.filter=0.75){

  #missing genotype imputation
  G_gene_buffer_surround[G_gene_buffer_surround<0 | G_gene_buffer_surround>2]<-NA
  N_MISS<-sum(is.na(G_gene_buffer_surround))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G_gene_buffer_surround)/ncol(G_gene_buffer_surround))
    #print(msg,call.=F)
    colmean<-colMeans(x = G_gene_buffer_surround, na.rm = T)
    index <- which(is.na(G_gene_buffer_surround), arr.ind=TRUE)
    G_gene_buffer_surround[index] <- colmean[index[,2]]
  }

  #sparse matrix operation
  MAF<-colMeans(G_gene_buffer_surround)/2;MAC<-colSums(G_gene_buffer_surround)
  MAF[MAF>0.5]<-1-MAF[MAF>0.5]
  MAC[MAF>0.5]<-nrow(G_gene_buffer_surround)*2-MAC[MAF>0.5]
  s<-colMeans(G_gene_buffer_surround^2)-colMeans(G_gene_buffer_surround)^2
  SNP.index<-which(MAF>0 & MAC>=25 & s!=0 & !is.na(MAF))

  if(length(SNP.index)<=1 ){
    msg<-'Number of variants with missing rate <=10% in the specified range is <=1'
    #print(msg,call.=F)
    stop
  }
  G_gene_buffer_surround<-G_gene_buffer_surround[,SNP.index,drop=F]

  #get positions and reorder G_gene_buffer_surround
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_gene_buffer_surround)))
  G_gene_buffer_surround<-G_gene_buffer_surround[,order(pos),drop=F]

  MAF<-colMeans(G_gene_buffer_surround)/2
  G_gene_buffer_surround<-as.matrix(G_gene_buffer_surround)
  G_gene_buffer_surround[,MAF>0.5 & !is.na(MAF)]<-2-G_gene_buffer_surround[,MAF>0.5 & !is.na(MAF)]
  MAF<-colMeans(G_gene_buffer_surround)/2;MAC<-colSums(G_gene_buffer_surround)

  G_gene_buffer_surround<-Matrix(G_gene_buffer_surround,sparse=T)
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_gene_buffer_surround)))
  n=dim(G_gene_buffer_surround)[1]

  max.corr=1
  while(max.corr>=LD.filter){ #max corr < 0.75
    #clustering and filtering
    G_gene_buffer_surround=G_gene_buffer_surround
    sparse.fit<-sparse.cor(G_gene_buffer_surround)
    cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov
    range(c(cor.X)[round(c(cor.X),digits = 2)!=1.00])
    max.corr=max(abs(c(cor.X)[round(c(cor.X),digits = 2)!=1.00]))

    Sigma.distance = as.dist(1 - abs(cor.X))
    if(ncol(G_gene_buffer_surround)>1){
      fit = hclust(Sigma.distance, method="complete")
      corr_max = 0.75
      clusters = cutree(fit, h=1-corr_max)
    }else{clusters<-1}

    ##apply the LD filter before knockoff generation
    #One variant is randomly selected as the representative per cluster.
    #If a cluster is inside the gene-buffer region, we prioritize to keep one variant inside the gene buffer region instead of outsides
    gene_buffer_ind=(pos>=gene_buffer_start&pos<=gene_buffer_end)

    set.seed(12345)
    temp.index.gene_buffer<-sample(sum(gene_buffer_ind))
    temp.index.gene_buffer<-temp.index.gene_buffer[match(unique(clusters[gene_buffer_ind]),clusters[gene_buffer_ind][temp.index.gene_buffer])]
    if(length(temp.index.gene_buffer)<=1 ){
      msg<-'Number of variants after LD filtering in the gene buffer is <=1'
      warning(msg,call.=F)
      break
    }
    gene_buffer.index=which(gene_buffer_ind)[temp.index.gene_buffer]

    ##Then filter other variants in +-100kb surrounding region
    temp.index.surround<-sample(length(pos)-sum(gene_buffer_ind))
    temp.index.surround<-temp.index.surround[match(unique(clusters[!gene_buffer_ind]),clusters[!gene_buffer_ind][temp.index.surround])]
    surround.index=which(!gene_buffer_ind)[temp.index.surround]
    surround.index=surround.index[!clusters[which(!gene_buffer_ind)[temp.index.surround]]%in%unique(clusters[gene_buffer_ind])]

    temp.index=unique(c(gene_buffer.index,surround.index))

    G_gene_buffer_surround<-G_gene_buffer_surround[,temp.index,drop=F]
    pos=pos[temp.index]
  }

  #print('generating knockoffs of gene buffer region') #knockoff-AL for gene buffer, adapt the code of KnockoffScreen-AL
  set.seed(12345)
  G_gene_buffer_knockoff=create.MK.AL_gene_buffer(X=G_gene_buffer_surround,pos=pos,
                                                  gene_buffer_start=gene_buffer_start,gene_buffer_end=gene_buffer_end,M=M,
                                                  corr_max=LD.filter,maxN.neighbor=Inf,
                                                  maxBP.neighbor=surround.region,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                                  thres.ultrarare=25,R2.thres=LD.filter)

  return(G_gene_buffer_knockoff)
}

#' Knockoff Generation of enhancer
#'
#' This function generates multiple knockoff genotypes for enhancer. The knockoff generations are optimized using shrinkage leveraging algorithm.
#'
#' @param G_enhancer_surround The genotype matrix of the surrounding region for enhancer.
#' @param enhancer_start The start position of enhancer.
#' @param enhancer_end The end position of enhancer.
#' @param surround.region Surrounding region for enhancer, default is +-50kb.
#' @param LD.filter The correlation threshold for hierarchical clustering. Default LD filter at correlation 0.75
#' @param M Numer of multiple knockoffs.
#' @return {A list file contains original genotype and M knockoff genotypes for enhancer.}
#' @examples
#'library(Matrix)
#'data('Knock.generation.enhancer.example')
#'G_enhancer_surround=Matrix(Knock.generation.enhancer.example$G_enhancer_surround)
#'pos_enhancer=Knock.generation.enhancer.example$pos_enhancer
#'enhancer_start=min(pos_enhancer)
#'enhancer_end=max(pos_enhancer)
#'
#'G_enhancer_knockoff=Knockoffgeneration.enhancer(G_enhancer_surround=G_enhancer_surround,
#'                                              enhancer_start=enhancer_start,
#'                                              enhancer_end=enhancer_end,
#'                                              M=5,surround.region=50000,LD.filter=0.75)
#'#M=5 knockoffs
#'G_enhancer_knockoff1=G_enhancer_knockoff[1,,]
#'G_enhancer_knockoff2=G_enhancer_knockoff[2,,]
#'G_enhancer_knockoff3=G_enhancer_knockoff[3,,]
#'G_enhancer_knockoff4=G_enhancer_knockoff[4,,]
#'G_enhancer_knockoff5=G_enhancer_knockoff[5,,]
#'
#' @import SKAT
#' @import Matrix
#' @import SPAtest
#' @import CompQuadForm
#' @import irlba
#' @export
Knockoffgeneration.enhancer=function(G_enhancer_surround=G_enhancer_surround,
                                     enhancer_start=enhancer_start,
                                     enhancer_end=enhancer_start,
                                     M=5,surround.region=50000,LD.filter=0.75){

  #missing genotype imputation
  G_enhancer_surround[G_enhancer_surround<0 | G_enhancer_surround>2]<-NA
  N_MISS<-sum(is.na(G_enhancer_surround))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G_enhancer_surround)/ncol(G_enhancer_surround))
    print(msg,call.=F)
    colmean<-colMeans(x = G_enhancer_surround, na.rm = T)
    index <- which(is.na(G_enhancer_surround), arr.ind=TRUE)
    G_enhancer_surround[index] <- colmean[index[,2]]
  }

  #sparse matrix operation
  MAF<-colMeans(G_enhancer_surround)/2;MAC<-colSums(G_enhancer_surround)
  MAF[MAF>0.5]<-1-MAF[MAF>0.5]
  MAC[MAF>0.5]<-nrow(G_enhancer_surround)*2-MAC[MAF>0.5]
  s<-colMeans(G_enhancer_surround^2)-colMeans(G_enhancer_surround)^2
  SNP.index<-which(MAF>0 & MAC>=25 & s!=0 & !is.na(MAF))

  if(length(SNP.index)<=1 ){
    msg<-'Number of variants with missing rate <=10% in the specified range is <=1'
    print(msg,call.=F)
    stop
  }
  G_enhancer_surround<-G_enhancer_surround[,SNP.index,drop=F]

  #get positions and reorder G_enhancer_surround
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_enhancer_surround)))
  G_enhancer_surround<-G_enhancer_surround[,order(pos),drop=F]

  MAF<-colMeans(G_enhancer_surround)/2
  G_enhancer_surround<-as.matrix(G_enhancer_surround)
  G_enhancer_surround[,MAF>0.5 & !is.na(MAF)]<-2-G_enhancer_surround[,MAF>0.5 & !is.na(MAF)]
  MAF<-colMeans(G_enhancer_surround)/2;MAC<-colSums(G_enhancer_surround)

  G_enhancer_surround<-Matrix(G_enhancer_surround,sparse=T)
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_enhancer_surround)))
  n=dim(G_enhancer_surround)[1]

  max.corr=1
  while(max.corr>=LD.filter){ #max corr < 0.75
    #clustering and filtering
    G_enhancer_surround=G_enhancer_surround
    sparse.fit<-sparse.cor(G_enhancer_surround)
    cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov
    range(c(cor.X)[round(c(cor.X),digits = 2)!=1.00])
    max.corr=max(abs(c(cor.X)[round(c(cor.X),digits = 2)!=1.00]))

    Sigma.distance = as.dist(1 - abs(cor.X))
    if(ncol(G_enhancer_surround)>1){
      fit = hclust(Sigma.distance, method="complete")
      corr_max = 0.75
      clusters = cutree(fit, h=1-corr_max)
    }else{clusters<-1}

    ##apply the LD filter before knockoff generation
    #One variant is randomly selected as the representative per cluster.
    #If a cluster is inside the gene-buffer region, we prioritize to keep one variant inside the gene buffer region instead of outsides
    enhancer_ind=(pos>=enhancer_start&pos<=enhancer_end)

    set.seed(12345)
    temp.index.enhancer<-sample(sum(enhancer_ind))
    temp.index.enhancer<-temp.index.enhancer[match(unique(clusters[enhancer_ind]),clusters[enhancer_ind][temp.index.enhancer])]
    if(length(temp.index.enhancer)<=1 ){
      msg<-'Number of variants after LD filtering in the gene buffer is <=1'
      warning(msg,call.=F)
      break
    }
    enhancer.index=which(enhancer_ind)[temp.index.enhancer]

    ##Then filter other variants in +-100kb surrounding region
    temp.index.surround<-sample(length(pos)-sum(enhancer_ind))
    temp.index.surround<-temp.index.surround[match(unique(clusters[!enhancer_ind]),clusters[!enhancer_ind][temp.index.surround])]
    surround.index=which(!enhancer_ind)[temp.index.surround]
    surround.index=surround.index[!clusters[which(!enhancer_ind)[temp.index.surround]]%in%unique(clusters[enhancer_ind])]

    temp.index=unique(c(enhancer.index,surround.index))

    G_enhancer_surround<-G_enhancer_surround[,temp.index,drop=F]
    pos=pos[temp.index]
  }

  print('generating knockoffs of enhancer')
  set.seed(12345)
  G_enhancer_knockoff=create.MK.AL_enhancer(X=G_enhancer_surround,pos=pos,
                                            enhancer_start=enhancer_start,enhancer_end=enhancer_end,M=M,
                                            corr_max=LD.filter,maxN.neighbor=Inf,
                                            maxBP.neighbor=surround.region,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                            thres.ultrarare=25,R2.thres=LD.filter)

  return(G_enhancer_knockoff)
}


######### Other functions #########
#Optimize create.MK.AL function provided by Zihuai
#Knockoff generation for gene buffer regions
create.MK.AL_gene_buffer <- function(X=G_gene_buffer_surround,pos,gene_buffer_start,gene_buffer_end,M,corr_max=LD.filter,maxN.neighbor=Inf,
                                     maxBP.neighbor=surround.region,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                     thres.ultrarare=25,R2.thres=LD.filter) {

  method='shrinkage'
  sparse.fit<-sparse.cor(X)
  cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov  #correlation

  #svd to get leverage score, can be optimized;update: tried fast leveraging, but the R matrix is singular possibly because X is sparse.
  #Fast Truncated Singular Value Decomposition
  if(method=='shrinkage'){
    svd.X.u<-irlba(X,nv=floor(sqrt(ncol(X)*log(ncol(X)))))$u #U is the orthogonal singular vectors
    h1<-rowSums(svd.X.u^2)
    h2<-rep(1,nrow(X))
    prob1<-h1/sum(h1)
    prob2<-h2/sum(h2)
    prob<-0.5*prob1+0.5*prob2 #shrinkage leveraging estimator, probability weights for sampling
  }

  index.AL<-sample(1:nrow(X),min(n.AL,nrow(X)),replace = FALSE,prob=prob) #sampling r samples from n samples, using shrinkage leveraging estimator
  w<-1/sqrt(n.AL*prob[index.AL])
  rm(svd.X.u) #remove temp file

  X.AL<-w*X[index.AL,] #n.AL samples
  sum(is.na(X.AL)) #0

  sparse.fit<-sparse.cor(X.AL)
  cor.X.AL<-sparse.fit$cor;cov.X.AL<-sparse.fit$cov
  skip.index<-colSums(X.AL!=0)<=thres.ultrarare #skip features that are ultra sparse, permutation will be directly applied to generate knockoffs

  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single") #hierarchical clustering
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)  #variants from two different clusters do not have a correlation greater than 0.75.
  }else{clusters<-1}

  gc()
  X_k<-list()
  for(k in 1:M){
    X_k[[k]]<-matrix(0,nrow=nrow(X),ncol=ncol(X))
    #X_k[[k]]<-big.matrix(nrow=nrow(X),ncol=ncol(X),init=0,shared=FALSE)
  }

  ##only run snps within gene buffer
  snps_ind=which(pos<=gene_buffer_end&pos>=gene_buffer_start)

  index.exist<-c()
  for (k in unique(clusters[snps_ind])){
    #print(paste0('cluster',k))
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)[which(clusters==k)%in%snps_ind]){
      #print(i)
      rate<-1;R2<-1;temp.maxN.neighbor<-maxN.neighbor
      while(R2>=R2.thres){ #avoid over-fitting
        temp.maxN.neighbor<-floor(temp.maxN.neighbor/rate)
        snp.pos=as.numeric(gsub("^.*\\:","",names(clusters[i])))
        #+-100kb surrounding region
        index.pos<-which(pos>=max(snp.pos-maxBP.neighbor,pos[1]) & pos<=min(snp.pos+maxBP.neighbor,pos[length(pos)]))
        #correlation between this snp with other snps in +-100kb surrounding region
        temp<-abs(cor.X[i,])
        temp[which(clusters==k)]<-0 #exclude variants if they are in the same cluster as the target variant
        temp[-index.pos]<-0 #only focus on +-100kb surrounding region
        temp[which(temp<=corr_base)]<-0
        index<-order(temp,decreasing=T)
        if(sum(temp!=0,na.rm=T)==0 | temp.maxN.neighbor==0){index<-NULL}else{
          index<-setdiff(index[1:min(length(index),floor((nrow(X))^(1/3)),temp.maxN.neighbor,sum(temp!=0,na.rm=T))],i)
        } #top K snps up to K=n^1/3=75

        y<-X[,i] #n samples
        if(length(index)==0){fitted.values<-0}
        if(i %in% skip.index){fitted.values<-0}
        if(!(i %in% skip.index |length(index)==0)){
          x.AL<-X.AL[,index,drop=F]; #n.AL by K
          n.exist<-length(intersect(index,index.exist))
          x.exist.AL<-matrix(0,nrow=nrow(X.AL),ncol=n.exist*M)
          if(length(intersect(index,index.exist))!=0){
            for(j in 1:M){ # this is the most time-consuming part
              x.exist.AL[,((j-1)*n.exist+1):(j*n.exist)]<-w*X_k[[j]][index.AL,intersect(index,index.exist),drop=F]
            }
          }
          y.AL<-w*X[index.AL,i]; #n.AL

          temp.xy<-rbind(mean(y.AL),crossprod(x.AL,y.AL)/length(y.AL)-colMeans(x.AL)*mean(y.AL))
          temp.xy<-rbind(temp.xy,crossprod(x.exist.AL,y.AL)/length(y.AL)-colMeans(x.exist.AL)*mean(y.AL))
          temp.cov.cross<-sparse.cov.cross(x.AL,x.exist.AL)$cov
          temp.cov<-sparse.cor(x.exist.AL)$cov
          temp.xx<-cov.X.AL[index,index]
          temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))
          temp.xx<-cbind(0,temp.xx)
          temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)

          svd.fit<-svd(temp.xx)
          v<-svd.fit$v
          cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
          n.svd<-which(cump>=0.999)[1]
          svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
          temp.inv<-v[,svd.index,drop=F]%*%(svd.fit$d[svd.index]^(-1)*t(v[,svd.index,drop=F]))
          temp.beta<-temp.inv%*%temp.xy #least square estimate for regression coefficient, alpha and beta_k

          x<-X[,index,drop=F]
          temp.j<-1
          fitted.values<-temp.beta[1]+x%*%temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F]-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])
          length(fitted.values) #n samples

          if(length(intersect(index,index.exist))!=0){
            temp.j<-temp.j+ncol(x)
            for(j in 1:M){
              temp.x<-X_k[[j]][,intersect(index,index.exist),drop=F]
              if(ncol(temp.x)>=1){
                fitted.values<-fitted.values+temp.x%*%temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F]-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
              }
              temp.j<-temp.j+ncol(temp.x)
            }
          }
        }
        residuals<-as.numeric(y-fitted.values)
        #overfitted model
        R2<-1-var(residuals,na.rm=T)/var(y,na.rm=T)
        rate<-rate*2;temp.maxN.neighbor<-length(index)
      }
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)
      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[[j]][,which(clusters==k)]<-round(cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F],digits=1)
    }
  }

  #save knockoffs of gene buffer region
  print('saving knockoffs of gene buffer region')
  #G_gene_buffer=X[,snps_ind]
  G_gene_buffer_knockoff <- array(0, dim = c(M, nrow(X), length(snps_ind)))
  for (j in 1:M) {
    G_gene_buffer_knockoff[j, ,] <-X_k[[j]][,snps_ind]
  }
  rm(X_k)

  #G_gene_buffer_knockoff=list(G_gene_buffer=G_gene_buffer,G_gene_buffer_knockoff=G_gene_buffer_knockoff)
  return(G_gene_buffer_knockoff)
}

create.MK.AL_enhancer <- function(X=G_enhancer_surround,pos,enhancer_start,enhancer_end,M,corr_max=0.75,maxN.neighbor=Inf,
                                  maxBP.neighbor=50000,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                  thres.ultrarare=25,R2.thres=0.75) {

  method='shrinkage'
  sparse.fit<-sparse.cor(X)
  cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov

  #svd to get leverage score, can be optimized;update: tried fast leveraging, but the R matrix is singular possibly because X is sparse.
  if(method=='shrinkage'){
    svd.X.u<-irlba(X,nv=floor(sqrt(ncol(X)*log(ncol(X)))))$u
    h1<-rowSums(svd.X.u^2)
    h2<-rep(1,nrow(X))
    prob1<-h1/sum(h1)
    prob2<-h2/sum(h2)
    prob<-0.5*prob1+0.5*prob2
  }

  index.AL<-sample(1:nrow(X),min(n.AL,nrow(X)),replace = FALSE,prob=prob)
  w<-1/sqrt(n.AL*prob[index.AL])
  rm(svd.X.u) #remove temp file

  X.AL<-w*X[index.AL,]
  sparse.fit<-sparse.cor(X.AL)
  cor.X.AL<-sparse.fit$cor;cov.X.AL<-sparse.fit$cov
  skip.index<-colSums(X.AL!=0)<=thres.ultrarare #skip features that are ultra sparse, permutation will be directly applied to generate knockoffs

  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single")
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)
  }else{clusters<-1}

  X_k<-list()
  ##only focus on snps within gene buffer
  for(k in 1:M){
    #X_k[[k]]<-big.matrix(nrow=nrow(X),ncol=ncol(X),init=0,shared=FALSE)
    X_k[[k]]<-matrix(0,nrow=nrow(X),ncol=ncol(X))
  }

  snps_ind=which(pos<=enhancer_end&pos>=enhancer_start)

  index.exist<-c()
  for (k in unique(clusters[snps_ind])){
    #print(paste0('cluster',k))
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)[which(clusters==k)%in%snps_ind]){
      #print(i)
      rate<-1;R2<-1;temp.maxN.neighbor<-maxN.neighbor

      while(R2>=R2.thres){

        temp.maxN.neighbor<-floor(temp.maxN.neighbor/rate)
        snp.pos=as.numeric(gsub("^.*\\:","",names(clusters[i])))
        index.pos<-which(pos>=max(snp.pos-maxBP.neighbor,pos[1]) & pos<=min(snp.pos+maxBP.neighbor,pos[length(pos)]))

        temp<-abs(cor.X[i,]);temp[which(clusters==k)]<-0;temp[-index.pos]<-0
        temp[which(temp<=corr_base)]<-0

        index<-order(temp,decreasing=T)
        if(sum(temp!=0,na.rm=T)==0 | temp.maxN.neighbor==0){index<-NULL}else{
          index<-setdiff(index[1:min(length(index),floor((nrow(X))^(1/3)),temp.maxN.neighbor,sum(temp!=0,na.rm=T))],i)
        }

        y<-X[,i]
        if(length(index)==0){fitted.values<-0}
        if(i %in% skip.index){fitted.values<-0}
        if(!(i %in% skip.index |length(index)==0)){

          x.AL<-X.AL[,index,drop=F];
          n.exist<-length(intersect(index,index.exist))
          x.exist.AL<-matrix(0,nrow=nrow(X.AL),ncol=n.exist*M)
          if(length(intersect(index,index.exist))!=0){
            for(j in 1:M){ # this is the most time-consuming part
              x.exist.AL[,((j-1)*n.exist+1):(j*n.exist)]<-w*X_k[[j]][index.AL,intersect(index,index.exist),drop=F]
            }
          }
          y.AL<-w*X[index.AL,i];

          temp.xy<-rbind(mean(y.AL),crossprod(x.AL,y.AL)/length(y.AL)-colMeans(x.AL)*mean(y.AL))
          temp.xy<-rbind(temp.xy,crossprod(x.exist.AL,y.AL)/length(y.AL)-colMeans(x.exist.AL)*mean(y.AL))
          temp.cov.cross<-sparse.cov.cross(x.AL,x.exist.AL)$cov
          temp.cov<-sparse.cor(x.exist.AL)$cov
          temp.xx<-cov.X.AL[index,index]
          temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))
          temp.xx<-cbind(0,temp.xx)
          temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)

          svd.fit<-svd(temp.xx)
          v<-svd.fit$v
          cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
          n.svd<-which(cump>=0.999)[1]
          svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
          temp.inv<-v[,svd.index,drop=F]%*%(svd.fit$d[svd.index]^(-1)*t(v[,svd.index,drop=F]))
          temp.beta<-temp.inv%*%temp.xy

          x<-X[,index,drop=F]
          temp.j<-1
          fitted.values<-temp.beta[1]+x%*%temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F]-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])

          if(length(intersect(index,index.exist))!=0){
            temp.j<-temp.j+ncol(x)
            for(j in 1:M){
              temp.x<-X_k[[j]][,intersect(index,index.exist),drop=F]
              if(ncol(temp.x)>=1){
                fitted.values<-fitted.values+temp.x%*%temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F]-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
              }
              temp.j<-temp.j+ncol(temp.x)
            }
          }
        }
        residuals<-as.numeric(y-fitted.values)
        #overfitted model
        R2<-1-var(residuals,na.rm=T)/var(y,na.rm=T)
        rate<-rate*2;temp.maxN.neighbor<-length(index)
      }
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)
      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[[j]][,which(clusters==k)]<-round(cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F],digits=1)
    }
  }

  #save knockoffs of enhancer
  print('saving knockoffs of enhancer')
  G_enhancer_knockoff <- array(0, dim = c(M, nrow(X), length(snps_ind)))
  for (j in 1:M) {
    G_enhancer_knockoff[j, ,] <-X_k[[j]][,snps_ind]
  }
  rm(X_k)

  return(G_enhancer_knockoff)
}
