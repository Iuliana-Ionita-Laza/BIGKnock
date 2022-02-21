#' @importFrom stats as.dist cutree dbeta hclust pcauchy pchisq qchisq rbinom sd var
utils::globalVariables(c('G_gene_buffer','G_EnhancerAll','p_EnhancerAll','pos_gene_buffer','G_Enhancer','n',
                         'G_enhancer_surround','pos_enhancer','G_gene_buffer_surround'))

#' Data example for GeneScan3D UKBB gene-based testing based on GLMM.
#'
#' This example dataset contains a fitted null GLMM results, variance ratio and sparse Sigma matrix, which are obtained by SAIGE/SAIGE-Gene package. Besides, it also contains the original and knockoff genotype for gene and buffer region and enhancer, as well as the positions of gene buffer region.
#'
#' @name Example.GeneScan3D.UKB.GLMM
#' @docType data
#' @keywords data
#' @usage data("GeneScan3D.UKB.GLMM.example")
#' @examples
#' data("GeneScan3D.UKB.GLMM.example")
#'
#'result.null.model.GLMM=GeneScan3D.UKB.GLMM.example$result.null.model.GLMM
#'ratio=GeneScan3D.UKB.GLMM.example$ratio
#'sparseSigma=GeneScan3D.UKB.GLMM.example$sparseSigma
#'G_gene_buffer=GeneScan3D.UKB.GLMM.example$G_gene_buffer
#'pos_gene_buffer=GeneScan3D.UKB.GLMM.example$pos_gene_buffer
#'G_gene_buffer_knockoff1=GeneScan3D.UKB.GLMM.example$G_gene_buffer_knockoff1
#'G_enhancer=GeneScan3D.UKB.GLMM.example$G_enhancer
#'G_enhancer_knockoff1=GeneScan3D.UKB.GLMM.example$G_enhancer_knockoff1
#'Gsub.id=result.null.model.GLMM$sampleID
"GeneScan3D.UKB.GLMM.example"

#' Conduct GeneScan3D analysis on UKBB data using fitted null GLMM results.
#'
#' This function perform the gene-based test for each gene using the fitted null GLMM and estimated variance ratio obtained from SAIGE/SAIGE-Gene package. For binary traits, we conduct SPA gene-based tests to deal with imbalance case-control issues.
#'
#' @param G The genotype matrix in the gene buffer region, which is a n*p matrix where n is the number of individuals and p is the number of genetic variants in the gene buffer region.
#' @param G.EnhancerAll The genotype matrix for R enhancers, by combining the genotype matrix of each enhancer by columns.
#' @param R Number of enhancers.
#' @param p_Enhancer Number of variants in R enhancers, which is a 1*R vector.
#' @param window.size The 1-D window sizes in base pairs to scan the gene buffer region. The recommended window sizes are c(1000,5000,10000).
#' @param pos  The positions of genetic variants in the gene buffer region, an p dimensional vector. Each position corresponds to a column in the genotype matrix G and a row in the functional annotation matrix Z.
#' @param MAC.threshold Threshold for minor allele count. Variants below MAC.threshold are ultra-rare variants. The recommended level is 10.
#' @param MAF.threshold Threshold for minor allele frequency. Variants below MAF.threshold are rare variants. The recommended level is 0.01.
#' @param Gsub.id The subject id corresponding to the genotype matrix, an n dimensional vector. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param result.null.model.GLMM The fitted null GLMM results obtained from SAIGE/SAIGE-Gene package.
#' @param outcome 'C' for quantitative trait, 'D' for binary trait.
#' @param sparseSigma n by n sparse Sigma matrix obtained from SAIGE/SAIGE-Gene package.
#' @param ratio Variance ratio to calibrate test statistics,obtained from SAIGE/SAIGE-Gene package.
#' @return \item{GeneScan3D.Cauchy.pvalue}{Cauchy combination p-values of all, common and rare variants for GeneScan3D analysis.}
#' @return \item{M}{Number of 1D scanning windows.}
#' @return \item{minp}{Minimum p-values of all, common and rare variants for 3D windows.}
#' @return \item{RE_minp}{The regulartory elements in the 3D windows corresponding to the minimum p-values, for all, common and rare variants. 0 represents promoter and a number from 1 to R represents promoter and r-th enhancer.}
#' @examples
#'data("GeneScan3D.UKB.GLMM.example")
#'result.null.model.GLMM=GeneScan3D.UKB.GLMM.example$result.null.model.GLMM
#'ratio=GeneScan3D.UKB.GLMM.example$ratio
#'sparseSigma=GeneScan3D.UKB.GLMM.example$sparseSigma
#'G_gene_buffer=GeneScan3D.UKB.GLMM.example$G_gene_buffer
#'pos_gene_buffer=GeneScan3D.UKB.GLMM.example$pos_gene_buffer
#'G_gene_buffer_knockoff1=GeneScan3D.UKB.GLMM.example$G_gene_buffer_knockoff1
#'G_enhancer=GeneScan3D.UKB.GLMM.example$G_enhancer
#'G_enhancer_knockoff1=GeneScan3D.UKB.GLMM.example$G_enhancer_knockoff1
#'Gsub.id=result.null.model.GLMM$sampleID
#'
#'G_EnhancerAll=G_enhancer
#'p_EnhancerAll=dim(G_enhancer)[2]
#'G_EnhancerAll_knockoff1=G_enhancer_knockoff1
#'p_EnhancerAll_knockoff1=dim(G_enhancer_knockoff1)[2]
#'R=1
#'
#'result.GeneScan3D_orginal=GeneScan3D.UKB.GLMM(G=G_gene_buffer,
#'                                              G.EnhancerAll=G_EnhancerAll,
#'                                              R=R,
#'                                              p_Enhancer=p_EnhancerAll,
#'                                              window.size=c(1000,5000,10000),
#'                                              pos=pos_gene_buffer,
#'                                              MAC.threshold=10,
#'                                              MAF.threshold=0.01,
#'                                              Gsub.id=Gsub.id,
#'                                              result.null.model.GLMM,
#'                                              outcome='C',
#'                                             sparseSigma=sparseSigma,
#'                                             ratio=ratio)
#'result.GeneScan3D_orginal$GeneScan3D.Cauchy.pvalue[1]
#'
#' @import SKAT
#' @import Matrix
#' @import SPAtest
#' @import CompQuadForm
#' @export
GeneScan3D.UKB.GLMM<-function(G=G_gene_buffer,G.EnhancerAll=G_EnhancerAll,R=length(p_EnhancerAll),
                              p_Enhancer=p_EnhancerAll,window.size=c(1000,5000,10000),pos=pos_gene_buffer,
                              MAC.threshold=10,MAF.threshold=0.01,Gsub.id=Gsub.id,
                              result.null.model.GLMM=result.null.model.GLMM,outcome='C',
                              sparseSigma=sparseSigma,ratio=ratio){

  #load preliminary features
  mu<-as.vector(result.null.model.GLMM$fitted.values)
  Y.res<-as.vector(result.null.model.GLMM$residuals)
  X<-result.null.model.GLMM$X #covariates include intercept

  invSigma_X<-solve(sparseSigma, X, sparse=T)
  C<-solve(t(X)%*%invSigma_X)

  #genotype filtering/checking/missing values imputation
  G_filter=Genotype_filter(G,pos,impute.method='fixed')
  G=G_filter$G
  pos=G_filter$pos

  #match phenotype id (phecode) and genotype id
  if(length(Gsub.id)==0){match.index<-match(as.numeric(result.null.model.GLMM$sampleID),1:nrow(G))}else{
    match.index<-match(result.null.model.GLMM$sampleID,Gsub.id)
  }
  if(mean(is.na(match.index))>0){
    msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
    warning(msg,call.=F)
  }

  #individuals ids are matched with genotype
  G=Matrix(G[match.index,])

  #generate window matrix to specify the variants in each window
  window.matrix0_gene_buffer<-c()
  for(size in window.size){
    if (size==1){next}
    pos.tag<-seq(min(pos),max(pos),by=size*1/2)
    pos.tag<-sapply(pos.tag,function(x)pos[which.min(abs(x-pos))])
    window.matrix0_gene_buffer<-cbind(window.matrix0_gene_buffer,sapply(pos.tag,function(x)as.numeric(pos>=x & pos<x+size)))
  }

  window.string_gene_buffer<-apply(window.matrix0_gene_buffer,2,function(x)paste(as.character(x),collapse = ""))
  window.matrix_gene_buffer<-Matrix(window.matrix0_gene_buffer[,match(unique(window.string_gene_buffer),window.string_gene_buffer)])
  #Number of 1-D windows to scan the gene buffer region
  M_gene_buffer=dim(window.matrix_gene_buffer)[2]

  ##single variant score tests, related samples using SAIGE null GLMM
  if(outcome=='D'){v=as.numeric((mu*(1-mu)))}
  if(outcome=='C'){v=1/result.null.model.GLMM$theta[1]} #phi is residual variance

  #covariate adjusted genotypes
  G_tilde=G-X%*%solve(t(X)%*%(v*X))%*%(t(X)%*%(v*G))

  #variance-adjusted score statistics
  #as.vector(t(G_tilde)%*%Y.res)==as.vector(t(G)%*%Y.res)
  S=as.vector(t(G_tilde)%*%Y.res)/result.null.model.GLMM$theta[1]

  ##GLMM
  #adjusted score statistics, without SPA
  if(outcome=='C'){
    invSigma_G_tilde<-solve(sparseSigma, G_tilde, sparse=T)
    V=t(G_tilde)%*%invSigma_G_tilde
    p.single=pchisq(S^2/(ratio*diag(V)),df=1,lower.tail=F)
  }

  #with SPA
  if(outcome=='D'){
    qtilde =S/sqrt(ratio) +as.vector(t(G_tilde)%*%mu)
    #The term as.vector(t(G_tilde)%*%mu) would be removed in SPAtest:::Saddle_Prob
    #keep the ratio to estimate variance of scores in Saddle_Prob
    p.single=rep(NA,ncol(G))
    for (p in 1:ncol(G)){
      p.single[p]=Saddle_Prob(q=as.vector(qtilde)[p], mu = mu, g = G_tilde[,p])$p.value
    }
  }

  GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
  #Burden test: for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically; for binary traits, use SPA gene- or region-based score test
  #SKAT test: for continuous traits, compute p-value use Davies; for binary traits, use SPA gene- or region-based score test
  for (m in 1:M_gene_buffer){
    print(paste0('1D-window',m))

    #Create index for each window
    index.window<-(window.matrix_gene_buffer[,m]==1)
    G.window=G[,index.window]
    G.window=Matrix(G.window)

    #if there is no variant in this window, then do not conduct combined test in this window, move to the next one
    if(dim(G.window)[2]==0){
      next
    }

    MAF.window<-apply(G.window,2,mean)/2
    MAC.window<-apply(G.window,2,sum)
    weight.beta_125<-dbeta(MAF.window,1,25)
    weight.beta_1<-dbeta(MAF.window,1,1)

    weight.matrix<-cbind(MAC.window<MAC.threshold,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*weight.beta_125,(MAF.window>=MAF.threshold)*weight.beta_1)
    #ultra-rare variants, rare and common variants
    colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
    weight.matrix<-Matrix(weight.matrix)

    #Single variant score test for all variants in the window, SPA p-values for binary traits
    p.single.window<-p.single[index.window]

    #approximation the covariance matrix for GLMM: t(G) P_S G
    #G.window=G.window-X%*%solve(t(X)%*%(v*X))%*%(t(X)%*%(v*G.window))
    invSigma_G.window<-solve(sparseSigma, G.window, sparse=T)

    A<-t(G.window)%*%invSigma_G.window
    B<-t(X)%*%invSigma_G.window
    K_S=A-t(B)%*%C%*%B
    #adjusted covariance matrix
    K=K_S*ratio

    #SPA gene-based tests
    if(outcome=='D'){
      V=diag(K)
      #adjusted variance
      v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
      #adjusted covariance matrix
      K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
    }

    #Burden test: for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically
    #for binary traits, calculate the SPA gene-based p-value of Burden
    p.burden<-matrix(NA,1,ncol(weight.matrix))
    for (j in 1:ncol(weight.matrix)){
      if (sum(weight.matrix[,j]!=0)>1){
        #only conduct Burden test for at least 1 variants
        temp.window.matrix<-weight.matrix[,j]
        G.window2<-as.matrix(G.window%*%temp.window.matrix)
        weights=as.vector(weight.matrix[,j])
        if(outcome=='D'){ #SPA-adjusted
          p.burden[,j]<-pchisq(as.numeric((t(G.window2)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) ;
        }else{
          #continuous
          p.burden[,j]<-pchisq(as.numeric((t(G.window2)%*%Y.res/result.null.model.GLMM$theta[1])^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) ;
        }
      }
    }

    score<-as.vector(S)[index.window]
    p.dispersion<-matrix(NA,1,ncol(weight.matrix))
    #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling based moment matching
    weight.matrix0=(MAC.window>=MAC.threshold)*weight.matrix
    for (j in 2:ncol(weight.matrix)){
      if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
        if(outcome=='D'){
          #binary
          p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j])
        }else{
          #continuous
          p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j])
        }
      }
    }

    p.individual1<-Get.cauchy.scan(p.single.window,as.matrix((MAC.window>=MAC.threshold & MAF.window<MAF.threshold))) #rare variants
    p.individual2<-Get.cauchy.scan(p.single.window,as.matrix((MAF.window>=MAF.threshold))) #common and low frequency variants
    p.individual<-cbind(p.burden,p.dispersion,p.individual1,p.individual2);
    colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')

    p.Cauchy<-as.matrix(apply(p.individual,1,Get.cauchy))
    #aggregated Cauchy association test
    test.common<-grep('MAF>=MAF.threshold',colnames(p.individual))
    p.Cauchy.common<-as.matrix(apply(p.individual[,test.common,drop=FALSE],1,Get.cauchy))
    p.Cauchy.rare<-as.matrix(apply(p.individual[,-test.common,drop=FALSE],1,Get.cauchy))
    GeneScan1D.Cauchy.window[m,]=c(p.Cauchy,p.Cauchy.common,p.Cauchy.rare)
  }
  GeneScan1D.Cauchy=c(Get.cauchy(GeneScan1D.Cauchy.window[,1]),Get.cauchy(GeneScan1D.Cauchy.window[,2]),Get.cauchy(GeneScan1D.Cauchy.window[,3]))

  ###Obtain p-values for R enhancers
  GeneScan3D.Cauchy.EnhancerAll=c()
  if(R!=0){
    for (r in 1:R){ #Loop for each enhancer
      print(paste0('Enhancer',r))
      if (r==1){
        G.Enhancer=as.matrix(G.EnhancerAll[,1:cumsum(p_Enhancer)[r]])
      }else{
        G.Enhancer=as.matrix(G.EnhancerAll[,(cumsum(p_Enhancer)[r-1]+1):cumsum(p_Enhancer)[r]])
      }

      G.Enhancer=Genotype_filter_Enhancer(G.Enhancer=G.Enhancer,impute.method='fixed')

      #individuals ids are matched with genotype
      G.window.Enhancer=Matrix(G.Enhancer[match.index,])
      MAF.window.Enhancer<-apply(G.window.Enhancer,2,mean)/2
      MAC.window.Enhancer<-apply(G.window.Enhancer,2,sum)

      weight.beta_125<-dbeta(MAF.window.Enhancer,1,25)
      weight.beta_1<-dbeta(MAF.window.Enhancer,1,1)
      weight.matrix<-cbind(MAC.window.Enhancer<MAC.threshold,(MAF.window.Enhancer<MAF.threshold&MAC.window.Enhancer>=MAC.threshold)*weight.beta_125,(MAF.window.Enhancer>=MAF.threshold)*weight.beta_1)
      colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
      weight.matrix<-Matrix(weight.matrix)

      #Single variant score test for all variants in the enhancer
      G_tilde.Enhancer=G.window.Enhancer-X%*%solve(t(X)%*%(v*X))%*%(t(X)%*%(v*G.window.Enhancer))

      S.Enhancer=as.vector(t(G_tilde.Enhancer)%*%Y.res)/result.null.model.GLMM$theta[1]

      ##GLMM
      #adjusted score statistics, without SPA
      if(outcome=='C'){
        invSigma_G_tilde.Enhancer<-solve(sparseSigma, G_tilde.Enhancer, sparse=T)
        V.Enhancer=t(G_tilde.Enhancer)%*%invSigma_G_tilde.Enhancer
        p.single.Enhancer=pchisq(S.Enhancer^2/(ratio*diag(V.Enhancer)),df=1,lower.tail=F)
      }
      #with SPA
      if(outcome=='D'){
        #Observed test statistic
        qtilde.Enhancer =as.vector(S.Enhancer)/sqrt(ratio) +as.vector(t(G_tilde.Enhancer)%*%mu)
        #The term as.vector(t(G_tilde.Enhancer)%*%mu) would be removed in SPAtest:::Saddle_Prob
        #keep the ratio to estimate variance of scores in Saddle_Prob
        p.single.Enhancer=rep(NA,ncol(G_tilde.Enhancer))
        for (p in 1:ncol(G_tilde.Enhancer)){
          p.single.Enhancer[p]=Saddle_Prob(q=as.vector(qtilde.Enhancer)[p], mu = mu, g = G_tilde.Enhancer[,p])$p.value
        }
      }

      p.burden.Enhancer<-matrix(NA,1,ncol(weight.matrix))
      p.dispersion.Enhancer<-matrix(NA,1,ncol(weight.matrix))

      if(length(p.single.Enhancer)>1){
        #enhancer have more than 1 variant, then conduct SKAT and burden; otherwise only conduct single variant score test
        #approximation the covariance matrix for GLMM
        #t(G) P_S G
        invSigma_G.window.Enhancer<-solve(sparseSigma, G.window.Enhancer, sparse=T)
        A<-t(G.window.Enhancer)%*%invSigma_G.window.Enhancer
        B<-t(X)%*%invSigma_G.window.Enhancer
        K_S=A-t(B)%*%C%*%B
        #adjusted covariance matrix
        K=K_S*ratio

        #SPA gene-based tests
        if(outcome=='D'){
          V=diag(K)
          #adjusted variance
          v_tilde=as.vector(S.Enhancer^2)/qchisq(p.single.Enhancer,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
          #adjusted covariance matrix
          K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
        }
      }

      #Burden
      for (j in 1:ncol(weight.matrix)){
        if (sum(weight.matrix[,j]!=0)>1){
          #only conduct Burden test for at least 1 variants
          temp.window.matrix<-weight.matrix[,j]
          G.window.Enhancer2<-as.matrix(G.window.Enhancer%*%temp.window.matrix)
          weights=as.vector(weight.matrix[,j])
          if(outcome=='D'){ #SPA-adjusted
            p.burden.Enhancer[,j]<-pchisq(as.numeric((t(G.window.Enhancer2)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F)
          }else{
            #continuous
            p.burden.Enhancer[,j]<-pchisq(as.numeric((t(G.window.Enhancer2)%*%Y.res/result.null.model.GLMM$theta[1])^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F)
          }
        }
      }

      #SKAT
      #For extremely rare variants, do not conduct SKAT
      for (j in 2:ncol(weight.matrix)){
        if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
          if(outcome=='D'){
            #binary
            p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(S.Enhancer,K=K_tilde,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j])
          }else{
            #continuous
            p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(S.Enhancer,K=K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j])
          }
        }
      }

      p.individual1.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAC.window.Enhancer>=MAC.threshold & MAF.window.Enhancer<MAF.threshold))) #rare variants
      p.individual2.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAF.window.Enhancer>=MAF.threshold))) #common and low frequency variants
      p.individual.Enhancer<-cbind(p.burden.Enhancer ,p.dispersion.Enhancer,p.individual1.Enhancer,p.individual2.Enhancer);
      colnames(p.individual.Enhancer)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')

      #aggregated Cauchy association test
      p.Cauchy.Enhancer<-as.matrix(apply(p.individual.Enhancer,1,Get.cauchy))
      test.common<-grep('MAF>=MAF.threshold',colnames(p.individual.Enhancer))
      p.Cauchy.common.Enhancer<-as.matrix(apply(p.individual.Enhancer[,test.common,drop=FALSE],1,Get.cauchy))
      p.Cauchy.rare.Enhancer<-as.matrix(apply(p.individual.Enhancer[,-test.common,drop=FALSE],1,Get.cauchy))
      GeneScan3D.Cauchy.Enhancer=c(p.Cauchy.Enhancer,p.Cauchy.common.Enhancer,p.Cauchy.rare.Enhancer)

      GeneScan3D.Cauchy.EnhancerAll=rbind(GeneScan3D.Cauchy.EnhancerAll,GeneScan3D.Cauchy.Enhancer)
    }  #end of the loop of R enhancers
  }
  ##Obtain 3D windows and p-values
  #do not add promoter
  #M 1D windows + Enhancer r, r=1, ..., R
  GeneScan3D.window.EnhancerAll=c()
  if(R!=0){
    for (r in 1:dim(GeneScan3D.Cauchy.EnhancerAll)[1]){
      GeneScan3D.window.enhancer=data.frame(apply(cbind(GeneScan1D.Cauchy.window[,1],GeneScan3D.Cauchy.EnhancerAll[r,1]),1,Get.cauchy),
                                            apply(cbind(GeneScan1D.Cauchy.window[,2],GeneScan3D.Cauchy.EnhancerAll[r,2]),1,Get.cauchy),
                                            apply(cbind(GeneScan1D.Cauchy.window[,3],GeneScan3D.Cauchy.EnhancerAll[r,3]),1,Get.cauchy))
      colnames(GeneScan3D.window.enhancer)=c('all','common','rare')
      GeneScan3D.window.EnhancerAll=rbind(GeneScan3D.window.EnhancerAll,GeneScan3D.window.enhancer)
    }
  }else{
    GeneScan3D.window.enhancer=data.frame(Get.cauchy(GeneScan1D.Cauchy.window[,1]),
                                          Get.cauchy(GeneScan1D.Cauchy.window[,2]),
                                          Get.cauchy(GeneScan1D.Cauchy.window[,3]))
    colnames(GeneScan3D.window.enhancer)=c('all','common','rare')
    GeneScan3D.window.EnhancerAll=rbind(GeneScan3D.window.EnhancerAll,GeneScan3D.window.enhancer)
  }

  GeneScan3D.Cauchy.RE=GeneScan3D.window.EnhancerAll
  GeneScan3D.Cauchy=c(Get.cauchy(GeneScan3D.Cauchy.RE[,1]), Get.cauchy(GeneScan3D.Cauchy.RE[,2]), Get.cauchy(GeneScan3D.Cauchy.RE[,3]))

  ###min-p and RE with min-p
  RE_minp.all=NA;RE_minp.common=NA;RE_minp.rare=NA
  if(R!=0){
    RE.indicator=c(rep(1:R,each=M_gene_buffer))

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))){
      RE_minp.all=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,1]==min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))])
    }

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))){
      RE_minp.common=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,2]==min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))])
    }

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))){
      RE_minp.rare=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,3]==min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))])
    }
  }

  #best enhancer
  return(list(GeneScan3D.Cauchy.pvalue=GeneScan3D.Cauchy,M=M_gene_buffer,R=R,
              minp=c(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE)),
              RE_minp=cbind(RE_minp.all,RE_minp.common,RE_minp.rare)))  #GeneScan1D.Cauchy.pvalue=GeneScan1D.Cauchy,
}


######### Other functions #########
Get.p.SKAT_noMA<-function(score,K,window.matrix,weight){

  Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2) #SKAT statistics
  K.temp<-weight*t(weight*K)

  temp<-K.temp[window.matrix[,1]!=0,window.matrix[,1]!=0]
  if(sum(temp^2)==0){p<-NA}else{
    lambda=eigen(temp,symmetric=T,only.values=T)$values #eigenvalues, mixture of chi-square
    temp.p<-SKAT_davies(Q,lambda,acc=10^(-6))$Qq

    if(length(temp.p)==0 || temp.p > 1 || temp.p <= 0){
      temp.p<-Get_Liu_PVal.MOD.Lambda(Q,lambda)
    }
    p<-temp.p
  }
  return(p)
}
SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  r <- length(lambda)
  if (length(h) != r) warning("lambda and h should have the same length!")
  if (length(delta) != r) warning("lambda and delta should have the same length!")
  #out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
  out=davies(q, lambda, h = rep(1, length(lambda)), delta = rep(0,length(lambda)), sigma = 0, lim = 10000, acc = 0.0001)
  out$res <- 1 - out$res
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
}
Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){
  param<-Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  return(p.value)
}
Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation

  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }

  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get.cauchy.scan<-function(p,window.matrix){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16)
  temp<-rep(0,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[!is.small]<-as.numeric(tan((0.5-p[!is.small])*pi))

  cct.stat<-as.numeric(t(temp)%*%window.matrix/apply(window.matrix,2,sum))
  is.large<-cct.stat>1e+15 & !is.na(cct.stat)
  is.regular<-cct.stat<=1e+15 & !is.na(cct.stat)
  pval<-rep(NA,length(cct.stat))
  pval[is.large]<-(1/cct.stat[is.large])/pi
  pval[is.regular]<-1-pcauchy(cct.stat[is.regular])
  return(pval)
}
Get.cauchy<-function(p){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))

  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}
Impute<-function(Z, impute.method){
  p<-dim(Z)[2]
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  return(as.matrix(Z))
}

Genotype_filter=function(G,pos,impute.method='fixed'){

  if(ncol(G)==0|ncol(G)==1){
    stop('Number of variants in the gene buffer region is 0 or 1')
  }

  #missing genotype imputation
  G[G==-9 | G==9]=NA
  N_MISS=sum(is.na(G))
  MISS.freq=apply(is.na(G),2,mean)

  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G=Impute(G,impute.method)
  }

  #MAF filtering
  MAF<-apply(G,2,mean)/2 #MAF of nonfiltered variants
  G[,MAF>0.5 & !is.na(MAF)]<-2-G[,MAF>0.5 & !is.na(MAF)]
  MAF<-apply(G,2,mean)/2
  s<-apply(G,2,sd)
  SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF))

  check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
  if(length(check.index)<=1 ){
    stop('Number of variants with missing rate <=10% in the gene plus buffer region is <=1')
  }

  G<-Matrix(G[,SNP.index])
  pos=pos[SNP.index]
  genotype_filter=list(G=G,pos=pos)
  return(genotype_filter)
}
Genotype_filter_Enhancer=function(G.Enhancer,impute.method='fixed'){

  #missing genotype imputation
  G.Enhancer[G.Enhancer==-9 | G.Enhancer==9]=NA
  N_MISS.Enhancer=sum(is.na(G.Enhancer))
  MISS.freq.Enhancer=apply(is.na(G.Enhancer),2,mean)
  if(N_MISS.Enhancer>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS.Enhancer/nrow(G.Enhancer)/ncol(G.Enhancer))
    warning(msg,call.=F)
    G.Enhancer=Impute(G.Enhancer,impute.method)
  }

  #MAF filtering
  MAF.Enhancer<-apply(G.Enhancer,2,mean)/2 #MAF of nonfiltered variants
  G.Enhancer[,MAF.Enhancer>0.5 & !is.na(MAF.Enhancer)]<-2-G.Enhancer[,MAF.Enhancer>0.5 & !is.na(MAF.Enhancer)]
  MAF.Enhancer<-apply(G.Enhancer,2,mean)/2
  s.Enhancer<-apply(G.Enhancer,2,sd)
  SNP.index.Enhancer<-which(MAF.Enhancer>0 & s.Enhancer!=0 & !is.na(MAF.Enhancer))

  G.Enhancer<-Matrix(G.Enhancer[,SNP.index.Enhancer])
  return(G.Enhancer)
}

#####knockoff AL functions
#percentage notation
percent <- function(x, digits = 3, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}
sparse.cov.cross <- function(x,y){
  n <- nrow(x)
  cMeans.x <- colMeans(x);cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x,y)) - n*tcrossprod(cMeans.x,cMeans.y))/(n-1)
  list(cov=covmat)
}
max.nth<-function(x,n){return(sort(x,partial=length(x)-(n-1))[length(x)-(n-1)])}


