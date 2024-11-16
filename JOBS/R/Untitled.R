library(MASS)
jobs.eqtls.new<-function(beta,se,weight,COR){
  if(is.null(COR)){
    COR=F
  }

  df<-beta
  df_s<-se

  cat(paste0("Check NA in sc-eqtl \n"))
  na.list<-list()
  for(r in 3:ncol(df)){
    snps<-df[which(is.na(df[,r])),1]
    names<-colnames(df)[r]
    na.list[[names]]<-snps
    cat(paste0(length(snps)," NA in ",names," \n"))
  }

  df[is.na(df)]<-0
  df_s[is.na(df_s)]<-1

  x_old<-df[,3:ncol(df)]

  cat(paste0("Check if all sc-eqtl effects are 0 \n"))
  cat(paste0(nrow(x_old), " gene-snp pairs in total \n"))

  sum1<-rowSums(x_old)

  length(sum1)
  length(which(sum1!=0))

  df<-df[which(sum1!=0),]

  df_s<-df_s[which(sum1!=0),]

  cat(paste0(nrow(df), " gene-snp pairs passed the filtering \n"))

  weight<-t(as.matrix(weight))
  weight_total<-t(matrix(as.numeric(weight),nrow=ncol(x_old),ncol=nrow(x_old)))
  pred<-as.matrix(x_old)%*%t(weight)
  sum1<-rowSums(x_old)
  se<-df_s[,3:ncol(df)]


  if(COR){
    cat(paste0("Estimating correlation across cell type \n"))
    z2<-(x_old/se)^2
    z2<-apply( z2,2,as.numeric)
    z2[is.na(z2)]<-0
    pval<-pchisq(z2,lower.tail = F,df=1)
    pval_min<-apply(pval,1,min)
    insig<-which(pval_min>0.05)
    cor_m<-cor(x_old[insig,])
  }
  cor_m<-as.data.frame(cor_m)


  se_new<-se
  x_new<-x_old


  cat(paste0("Step2: Start analyzing ",nrow(x_old)," gene-snp pairs \n"))

  if(!COR){
    ###jointly model eqtls
    for(j in 1:ncol(x_old)){
      #print(j)
      ###estimate se by fisher information
      se_new[,j]<-sqrt(df_s[,2]^2*se[,j]^2/(df_s[,2]^2+weight_total[,j]^2*se[,j]^2))
      ###estimate effect size by maximize likelihood jointly
      x_new[,j]<-(x_old[,j]*df_s[,2]^2+df[,2]*weight_total[,j]*se[,j]^2-(pred-x_old[,j]*weight_total[,j])*weight_total[,j]*se[,j]^2)/(df_s[,2]^2+se[,j]^2*weight_total[,j]^2)
    }
  }else{

    jobs.with.cor<-function(x){
      A<-rbind( as.numeric(weight_total[x,]),diag(ncol(x_new)))
      cov_m<-as.numeric(df_s[x,3:ncol(df_s)])*cor_m*as.numeric(df_s[x,3:ncol(df_s)])
      omega<-as.matrix(rbind(c(df_s[x,2]^2,rep(0,nrow(cor_m))),cbind(0,cov_m)))
      omega_inv<-ginv(omega)
      v<-ginv(t(A)%*%omega_inv%*%A)
      s<-v%*%t(A)%*%omega_inv%*%t(as.matrix(df[x,2:ncol(df)]))
      return(c(diag(v),s))
    }

    s1<-sapply(1:nrow(x_old),function(x)jobs.with.cor(x))
    se_new<-t(s1[1:ncol(x_old),]^0.5)
    x_new<-t(s1[(ncol(x_old)+1):(2*ncol(x_old)),])
  }

  ##check results
  df_new<-as.matrix(cbind(df[,1:2],x_new))
  df_s_new<-as.matrix(cbind(df_s[,1:2],se_new))


  colnames(df_new)<-colnames(df)
  colnames(df_s_new)<-colnames(df_s)

  ###replace NA list with NA
  for(r in 3:ncol(df_s)){
    names<-colnames(df)[r]
    snps<-na.list[[names]]
    df_new[na.omit(match(snps,df_new)),r]<-NA
    df_s_new[na.omit(match(snps,df_new)),r]<-NA
  }

  cat("Finished! \n")

  res<-list("eqtls_new"=as.data.frame(df_new),"eqtls_se_new"=as.data.frame(df_s_new))

  return(res)
}
