# Distribution calculation for sum of chi-square random variables

------
 - The problem reduces to compute the tail probability of weighted sum of 1-DF central chi-square rvs
   - ![f1]
   - ![f2] and X_i^2 are iid 1-DF central chi-square rvs
   - need to calculate Pr(Q>t)
 - A non-central chi-square dist approximation is implemented.
   - match (n-1)-th and n-th moments 
   - Ref: Wu, B., Pankow, J.S., 2017. On computing the tail probability of non-negative definite quadratic forms in central normal variables. *tech report*
   - The current implementation requires the R *minqa* package.
     - http://cran.r-project.org/web/packages/minqa/index.html
     - can use other numerical optimization prog. 
 - R codes
   - Liu method ref
     - Liu, H., Tang, Y., Zhang, H.H., 2009. A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables. Computational Statistics & Data Analysis 53, 853–856.
   - Implementation of various methods including Davies etc
     - Duchesne, P., Lafaye De Micheaux, P., 2010. Computing the distribution of quadratic forms: Further comparisons between the Liu–Tang–Zhang approximation and exact methods. Computational Statistics & Data Analysis 54, 858–862.
   - Lee method ref
     - Lee, S., Wu, M.C., Lin, X., 2012. Optimal tests for rare variant effects in sequencing association studies. Biostat 13, 762–775.
   - Satterwaite method ref
     - Kwee, L.C., Liu, D., Lin, X., Ghosh, D., Epstein, M.P., 2008. A Powerful and Flexible Multilocus Association Test for Quantitative Traits. The American Journal of Human Genetics 82, 386–397. 
```r
  library(minqa)
  #####
  cum2mnc = function(kappa){
  ### convert cumulants to non-central moments
  ###    recursive formula produces as many cumulants as moments
  ###    References: Kenneth Lange: Numerical Analysis for Statisticians, 2nd ed. Page 15
    N = length(kappa)+1
    mc = rep(0, N); mc[1] = 1
    for(k in 1:(N-1)){
      mc[k+1] = sum(choose(k-1, 0:(k-1))*kappa[k:1]*mc[1:k])
    }
    return(mc[-1])
  }
  mnc2mc = function(mnc){
  ### convert non-central to central moments, uses recursive formula
    N = length(mnc)
    mc = rep(0,N); mc[1] = 0
    s1 = rep(c(1,-1), N)
    mnc = c(1,mnc)
    for(k in 1:(N-1)){
      mc[k+1] = sum( choose(k+1, 0:(k+1))*s1[(k+2):1]*mnc[1:(k+2)]*mnc[2]^((k+1):0) )
    }
    return(mc)
  }
  #### non-central chi-square $\chi_k^2(\lambda)$ cumulants 
  chisq.cum = function(k, lam, N){
  ### k: DF; lam: ncp 
    ik = 1:N
    2^(ik-1)*gamma(ik)*(k+ik*lam)
  }
  ## 1-DF central chisq mix cumulants
  chi1sqm.cum = function(lam, N){
  ### lam: weight coef
    ik = 1:N
    a1 = 2^(ik-1)*gamma(ik)
    cl = rep(0, N)
    for(i in 1:N) cl[i] = a1[i]*sum(lam^i)
    cl
  }

  ## Satterwaite approximation
  satter.lambda = function(lam){
    E = sum(lam)
    V = 2*sum(lam^2)
    dta = V/E/2
    d = 2*E^2/V
    return(list(df=d,scale=dta))
  }  
  satter.pval = function(Qq,lam){
    E = sum(lam)
    V = 2*sum(lam^2)
    dta = V/E/2
    d = 2*E^2/V
    pchisq(Qq/dta,d,lower=FALSE)
  }
  ### Liu approach
  liu.lambda = function(lambda){
    c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lambda^i) }
    muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
    s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
    if(s1^2 > s2){
      a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
    } else {
      l = 1/s1^2;  a = sqrt(l);  d = 0
    }
    muX = l+d;  sigmaX = sqrt(2)*a
    list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  }
  liu.pval = function(Q.all, lambda){
    param = liu.lambda(lambda)
    Q.Norm = (Q.all - param$muQ)/param$sigmaQ
    Q.Norm1 = Q.Norm*param$sigmaX + param$muX
    pchisq(Q.Norm1, df = param$l,ncp=param$d, lower.tail=FALSE)
  }
  ### Lee SKATO approach
  lee.lambda = function(lambda){
    c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lambda^i) }
    muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
    s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
    if(s1^2 > s2){
      a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
    } else {
      l = 1/s2;  a = sqrt(l);  d = 0
    }
    muX = l+d;  sigmaX = sqrt(2)*a
    list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  }
  lee.pval = function(Q.all, lambda){
    param = lee.lambda(lambda)
    Q.Norm = (Q.all - param$muQ)/param$sigmaQ
    Q.Norm1 = Q.Norm*param$sigmaX + param$muX
    pchisq(Q.Norm1, df = param$l,ncp=param$d, lower.tail=FALSE)
  }
  ## match higher moments
  wu.lambda = function(lam, N=12){
    cl = chi1sqm.cum(lam, N)
    muQ = cl[1]; sigmaQ = sqrt(cl[2])
    a1 = mnc2mc(cum2mnc(cl))
    a1 = a1/sqrt(a1[2])^(1:N)  
    f1 = function(xpar){
      k = exp(xpar[1])
      v = xpar[2]
      a2 = mnc2mc(cum2mnc(chisq.cum(k,v,N)))
      a2 = a2/sqrt(a2[2])^(1:N)  
      (a1[N-1]-a2[N-1])^2 + (a1[N]-a2[N])^2
    }
    tmp = bobyqa(c(0,1), f1, lower=c(-Inf,0),upper=c(Inf,Inf))
    xpar = tmp$par
    l = exp(xpar[1])
    d = xpar[2]
    if(f1(c(xpar[1],0))<=tmp$fval){
      d=0
      f.1 = function(xpar) f1(c(xpar,0))
      l = exp(bobyqa(xpar[1], f.1)$par)
    }
    muX = l+d; sigmaX = sqrt(chisq.cum(l,d,N=2)[2])
    list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  }
  wu.pval = function(Q.all, lambda, N=12){
    param = wu.lambda(lambda,N)
    Q.Norm = (Q.all - param$muQ)/param$sigmaQ
    Q.Norm1 = Q.Norm*param$sigmaX + param$muX
    pchisq(Q.Norm1, df = param$l,ncp=param$d, lower.tail=FALSE)
  }
```
  - Compute inflation factor (IF): reproduce Table 2 in the paper.
```r
  library(CompQuadForm)
  ### IF comp
  satter.IF = function(pval, lampar1,lam){
    davies(qchisq(pval,lampar1$df, lower=FALSE)*lampar1$scale, lam, acc=1e-12, lim=1e8)$Qq/pval
  }
  wu.IF = function(pval, lampar,lam){
    Q1 = qchisq(pval,df=lampar$l,ncp=lampar$d,lower=FALSE)
    davies( (Q1-lampar$muX)/lampar$sigmaX*lampar$sigmaQ+lampar$muQ, lam, acc=1e-12, lim=1e8)$Qq/pval
  }
  liu.IF = function(pval, lampar,lam){
    Q1 = qchisq(pval,df=lampar$l,ncp=lampar$d,lower=FALSE)
    davies( (Q1-lampar$muX)/lampar$sigmaX*lampar$sigmaQ+lampar$muQ, lam, acc=1e-12, lim=1e8)$Qq/pval
  }
  lee.IF = function(pval, lampar,lam){
    Q1 = qchisq(pval,df=lampar$l,ncp=lampar$d,lower=FALSE)
    davies( (Q1-lampar$muX)/lampar$sigmaX*lampar$sigmaQ+lampar$muQ, lam, acc=1e-12, lim=1e8)$Qq/pval
  }

  pval = c(1e-4,1e-5,1e-6,1e-7)
  lam1 = c(0.1,0.5,0.8)
  lam2 = c(0.3,7)
  lam3 = c(1.68, 1.45, 1.01, 1.006, 1.003, 1.001, 1.0005, 1.00023, 1.00022, 0.98, 0.94, 0.54, 0.39)
  lam4 = c(7.05, 4.28, 2.35, 1.94, 1.05, 0.93, 0.76, 0.47, 0.38, 0.30,0.18, 0.10, 0.07, 0.06, 0.04, 0.03, 0.02, 0.003, 0.001, 0.0001)

  lam = lam1 ## assign lam1, lam2, lam3, lam4 to lam
  for(a in pval) cat(satter.IF(a, satter.lambda(lam), lam), '\n')
  for(a in pval) cat(liu.IF(a, liu.lambda(lam), lam), '\n')
  for(a in pval) cat(lee.IF(a, lee.lambda(lam), lam), '\n')
  for(a in pval) cat(wu.IF(a, wu.lambda(lam), lam), '\n')
```



[f1]: http://chart.apis.google.com/chart?cht=tx&chl=Q=\sum_{i=1}^m\lambda_iX_i^2
[f2]: http://chart.apis.google.com/chart?cht=tx&chl=\lambda_1\geq\ldots\geq\lambda_m\geq 0
