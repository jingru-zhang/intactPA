## The function implements the fast MFPCA method

fastmfpca <- function(Y, id, visit = NULL, twoway = TRUE, weight = "obs", argvals = NULL,
                       pve = 0.99, npc = NULL, pve2 = NULL, npc2 = NULL, p = 3, m = 2, 
                       knots = 35, lambda.Gt=NULL,lambda.Gw=NULL, silent = TRUE){

  pspline.setting.mfpca <- function(x,knots=35,p=3,m=2,weight=NULL,type="full",
                                    knots.option="equally-spaced"){
    # design matrix
    K = length(knots)-2*p-1
    B = spline.des(knots=knots, x=x, ord=p+1, outer.ok=TRUE,sparse=TRUE)$design
    bs = "ps"
    if(knots.option == "quantile"){
      bs = "bs"
    }
    s.object = s(x=x, bs=bs, k=K+p, m=c(p-1,2), sp=NULL)
    object  = smooth.construct(s.object,data = data.frame(x=x),knots=list(x=knots))
    P = object$S[[1]]

    if(knots.option == "quantile") P = P / max(abs(P))*10 # rescaling

    if(is.null(weight)) weight <- rep(1,length(x))

    if(type=="full"){
      Sig = crossprod(matrix.multiply.mfpca(B,weight,option=2),B)
      eSig = eigen(Sig)
      V = eSig$vectors
      E = eSig$values
      if(min(E)<=0.0000001) {
        E <- E + 0.000001
      }
      Sigi_sqrt = matrix.multiply.mfpca(V,1/sqrt(E))%*%t(V)
      tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
      Esig = eigen(tUPU,symmetric=TRUE)
      U = Esig$vectors
      s = Esig$values
      s[(K+p-m+1):(K+p)]=0
      A = B%*%(Sigi_sqrt%*%U)
    }

    if(type=="simple"){
      A = NULL
      s = NULL
      Sigi_sqrt = NULL
      U = NULL
    }

    List = list("A" = A, "B" = B, "s" = s, "Sigi.sqrt" = Sigi_sqrt, "U" = U, "P" = P)
    return(List)
  }

  quadWeights.mfpca <- function(argvals, method = "trapezoidal"){
    ret <- switch(method,
                  trapezoidal = {D <- length(argvals)
                  1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])},
                  midpoint = c(0,diff(argvals)),
                  stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
    return(ret)
  }

  ## The function implements the face algorithm 

  facecov <- function(Y, argvals, An, s, Cov=TRUE, pve=1, npc=NULL, lambda=NULL, alpha=0.7, 
                             search.grid=TRUE, search.length=100, lower=-20, upper=20){
    
    ######## precalculation for missing data ########
    imputation <- FALSE
    Niter.miss <- 1
    L <- ncol(Y)
    n <- nrow(Y)
    
    Index.miss <- is.na(Y)
    if(sum(Index.miss)>0){
      num.miss <- rowSums(is.na(Y))
      for(i in 1:n){
        if(num.miss[i]>0){
          y <- Y[i,]
          seq <- (1:L)[!is.na(y)]
          seq2 <-(1:L)[is.na(y)]
          t1 <- argvals[seq]
          t2 <- argvals[seq2]
          fit <- smooth.spline(t1,y[seq])
          temp <- predict(fit,t2,all.knots=TRUE)$y
          if(max(t2)>max(t1)) temp[t2>max(t1)] <- mean(y[seq])
          if(min(t2)<min(t1)) temp[t2<min(t1)] <- mean(y[seq])
          Y[i,seq2] <- temp
        }
      }
      imputation <- TRUE
      Niter.miss <- 100
    }
    convergence.vector <- rep(0,Niter.miss)
    iter.miss <- 1
    lambda.input <- lambda
    totalmiss <- mean(Index.miss)
    
    
    while(iter.miss <= Niter.miss&&convergence.vector[iter.miss]==0) {
      ###################################################
      ######## Transform the Data           #############
      ###################################################
      Ytilde <- t(as.matrix(Y%*%An))
      C_diag <- rowSums(Ytilde^2)
      
      ###################################################
      ########  Select Smoothing Parameters #############
      ###################################################
      Y_square <- sum(Y^2)
      Ytilde_square <- sum(Ytilde^2)
      face_gcv <- function(x) {
        lambda <- exp(x)
        lambda_s <- (lambda*s)^2/(1 + lambda*s)^2
        gcv <- sum(C_diag*lambda_s) - Ytilde_square + Y_square
        trace <- sum(1/(1+lambda*s))
        gcv <- gcv/(1-alpha*trace/L/(1-totalmiss))^2
        return(gcv)
      }
      
      
      if(is.null(lambda.input) && iter.miss<=2) {
        if(!search.grid){
          fit <- optim(0,face_gcv,lower=lower,upper=upper)
          if(fit$convergence>0) {
            expression <- paste("Smoothing failed! The code is:",fit$convergence)
            print(expression)
          }
          lambda <- exp(fit$par)
        } else {
          Lambda <- seq(lower,upper,length=search.length)
          Length <- length(Lambda)
          Gcv <- rep(0,Length)
          for(i in 1:Length)
            Gcv[i] <- face_gcv(Lambda[i])
          i0 <- which.min(Gcv)
          lambda <- exp(Lambda[i0])
        }
      }
      YS <- matrix.multiply.mfpca(Ytilde,1/(1+lambda*s),2)
      
      ###################################################
      ####  Eigendecomposition of Smoothed Data #########
      ###################################################
      temp0 <- YS%*%t(YS)/n
      Eigen <- eigen(temp0,symmetric=TRUE)
      Phi = An %*% Eigen$vectors
      Sigma = Eigen$values
      
      
      if(iter.miss>1&&iter.miss< Niter.miss) {
        diff <- norm(YS-YS.temp,"F")/norm(YS,"F")
        if(diff <= 0.02)
          convergence.vector[iter.miss+1] <- 1
      }
      
      YS.temp <- YS
      iter.miss <- iter.miss + 1
      N <- min(n, ncol(An))
      d <- Sigma[1:N]
      d <- d[d>0]
      per <- cumsum(d)/sum(d)
      N <- ifelse (is.null(npc), min(which(per>=pve)), min(npc, length(d)))
      
      #########################################
      #######     Principal  Scores   #########
      ########   data imputation      #########
      #########################################
      if(imputation) {
        Phi.N <- as.matrix(Phi[,1:N, drop = FALSE])
        d <- Sigma[1:N]
        sigmahat2  <-  max(mean(Y[!Index.miss]^2) -sum(Sigma)/nrow(An),0)
        if(N>1){
          Xi <- solve(diag(N) + diag(sigmahat2/d)) %*% t(as.matrix(Y%*%Phi.N))
        } else{
          Xi <- solve(diag(N) + sigmahat2/d) %*% t(as.matrix(Y%*%Phi.N))
        }
        Yhat <- t(Phi.N %*% Xi)
        Y <- Y*(1-Index.miss) + Yhat*Index.miss
        if(sum(is.na(Y))>0) print("error")
      }
      
    } ## end of while loop
    
    Phi.N <- as.matrix(Phi[,1:N, drop = FALSE])
    evalues <- Sigma[1:N]
    Ktilde <- NULL
    if(Cov) {
      Ktilde <- Phi.N %*%  matrix.multiply.mfpca(t(Phi.N),evalues,2)
    }
    
    return(list(Yhat=Y, decom=temp0, Ktilde=Ktilde, evalues=evalues, efunctions=Phi.N,lambda=lambda))
  }

  matrix.multiply.mfpca <- function(A,s,option=1){
    if(option==2)
      return(A*(s%*%t(rep(1,dim(A)[2]))))
    if(option==1)
      return(A*(rep(1,dim(A)[1])%*%t(s)))
  }



  ##################################################################################
  ## Organize the input
  ##################################################################################
  if(silent == FALSE) print("Organize the input")

  stopifnot((!is.null(Y) & !is.null(id)))
  stopifnot(is.matrix(Y))

  ## specify visit variable if not provided
  if (!is.null(visit)){
    visit <- as.factor(visit)
  }else{ ## if visit is not provided, assume the visit id are 1,2,... for each subject
    visit <- as.factor(ave(id, id, FUN=seq_along))
  }
  id <- as.factor(id) ## convert id into a factor

  ## organize data into one data frame
  df <- data.frame(id = id, visit = visit, Y = I(Y))
  rm(id, visit, Y)

  ## derive several variables that will be used later
  J <- length(levels(df$visit)) ## number of visits
  L <- ncol(df$Y) ## number of observations along the domain
  nVisits <- data.frame(table(df$id)[unique(df$id)])  ## calculate number of visits for each subject
  colnames(nVisits) = c("id", "numVisits")
  ID = sort(unique(df$id)) ## id of each subject
  I <- length(ID) ## number of subjects
  ## assume observations are equally-spaced on [0,1] if not specified
  if (is.null(argvals))  argvals <- seq(0, 1, length.out=L)


  ##################################################################################
  ## Estimate population mean function (mu) and visit-specific mean function (eta)
  ##################################################################################
  if(silent == FALSE) print("Estimate population and visit-specific mean functions")

  meanY <- colMeans(df$Y, na.rm = TRUE)
  fit_mu <- gam(meanY ~ s(argvals))
  mu <- as.vector(predict(fit_mu, newdata = data.frame(argvals = argvals)))
  rm(meanY, fit_mu)

  mueta = matrix(0, L, J)
  eta = matrix(0, L, J) ## matrix to store visit-specific means
  colnames(mueta) <- colnames(eta) <- levels(df$visit)
  Ytilde <- matrix(NA, nrow = nrow(df$Y), ncol = ncol(df$Y))
  if(twoway==TRUE) {
    for(j in 1:J) {
      ind_j <- which(df$visit == levels(df$visit)[j])
      if(length(ind_j) > 1){
        meanYj <- colMeans(df$Y[ind_j,], na.rm=TRUE)
      }else{
        meanYj <- df$Y[ind_j,]
      }
      fit_mueta <- gam(meanYj ~ s(argvals))
      mueta[,j] <- predict(fit_mueta, newdata = data.frame(argvals = argvals))
      eta[,j] <- mueta[,j] - mu
      Ytilde[ind_j,] <- df$Y[ind_j,] - matrix(mueta[,j], nrow = length(ind_j), ncol = L, byrow = TRUE)
    }
    rm(meanYj, fit_mueta, ind_j, j)
  } else{
    Ytilde <- df$Y - matrix(mu, nrow = nrow(df$Y), ncol = L, byrow = TRUE)
  }
  df$Ytilde <- I(Ytilde) ## Ytilde is the centered multilevel functional data
  rm(Ytilde)


  ##################################################################################
  ## FACE preparation: see Xiao et al. (2016) for details
  ##################################################################################
  if(silent == FALSE) print("Prepare ingredients for FACE")

  ## Specify the knots of B-spline basis
  if(length(knots)==1){
    if(knots+p>=L) cat("Too many knots!\n")
    stopifnot(knots+p<L)

    K.p <- knots
    knots <- seq(-p, K.p+p, length=K.p+1+2*p)/K.p
    knots <- knots*(max(argvals)-min(argvals)) + min(argvals)
  }
  if(length(knots)>1) K.p <- length(knots)-2*p-1
  if(K.p>=L) cat("Too many knots!\n")
  stopifnot(K.p < L)
  c.p <- K.p + p # the number of B-spline basis functions

  ## Precalculation for smoothing
  List <- pspline.setting.mfpca(argvals, knots, p, m)
  B <- List$B # B is the J by c design matrix
  Sigi.sqrt <- List$Sigi.sqrt # (t(B)B)^(-1/2)
  s <- List$s # eigenvalues of Sigi_sqrt%*%(P%*%Sigi_sqrt)
  U <- List$U # eigenvectors of Sigi_sqrt%*%(P%*%Sigi_sqrt)
  A0 <- Sigi.sqrt %*% U
  An = B%*%A0 # the transformed basis functions
  
  ##################################################################################
  ## Impute missing data of Y using FACE and estimate the total covariance (Kt)
  ##################################################################################
  if(silent == FALSE) print("Estimate the total covariance (Kt)")

  Ji <- as.numeric(table(df$id)[unique(df$id)])
  diagD <- rep(Ji, Ji)
  smooth.Gt = facecov(Y=unclass(df$Ytilde), argvals, An, s,lambda=lambda.Gt)
  ## impute missing data of Y using FACE approach
  if(sum(is.na(df$Ytilde))>0){
    df$Ytilde[which(is.na(df$Ytilde))] <- smooth.Gt$Yhat[which(is.na(df$Ytilde))]
    df$Y <- df$Ytilde + matrix(mu, nrow = nrow(df$Y), ncol = L, byrow = TRUE)
  }
  if(weight=="subj"){
    YH <- unclass(df$Ytilde)*sqrt(nrow(df$Ytilde)/(I*diagD))
    smooth.Gt <- facecov(Y=YH, argvals, An, s,lambda=lambda.Gt)
  }
  diag_Gt <- colMeans(df$Ytilde^2)


  ##################################################################################
  ## Estimate principal components of the within covariance (Kw)
  ##################################################################################
  if(silent == FALSE) print("Estimate principal components of the within covariance (Kw)")

  inx_row_ls <- split(1:nrow(df$Ytilde), f=factor(df$id, levels=unique(df$id)))
  Ysubm <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x,,drop=FALSE],na.rm=TRUE), numeric(L)))
  if(weight=="obs"){
    weights <- sqrt(nrow(df$Ytilde)/(sum(diagD) - nrow(df$Ytilde)))
    YR <-  do.call("rbind",lapply(1:I, function(x) {
      weights*sqrt(Ji[x]) * t(t(df$Ytilde[inx_row_ls[[x]],,drop=FALSE]) - Ysubm[x,]/Ji[x])
    }))
  }
  if(weight=="subj"){
    weights <- sqrt(nrow(df$Ytilde)/sum(Ji>1))
    YR <-  do.call("rbind",lapply(1:I, function(x) {
      if(Ji[x]>1) return((weights/sqrt(Ji[x]-1)) * t(t(df$Ytilde[inx_row_ls[[x]],,drop=FALSE]) - Ysubm[x,]/Ji[x]))
    }))
  }
  smooth.Gw <- facecov(Y=YR, argvals, An, s,lambda=lambda.Gw)
  sigma.Gw <- smooth.Gw$evalues # raw eigenvalues of Gw
  per <- cumsum(sigma.Gw)/sum(sigma.Gw)
  if (is.null(pve2)) pve2 <- pve
  if (!is.null(npc)){
    if (is.null(npc2)) npc2 <- npc
  }
  N.Gw <- ifelse (is.null(npc2), min(which(per>=pve2)), min(npc2, length(sigma.Gw)))
  pve2 <- per[N.Gw]
  smooth.Gw$efunctions <- smooth.Gw$efunctions[,1:N.Gw]
  smooth.Gw$evalues <- smooth.Gw$evalues[1:N.Gw]
 
  ##################################################################################
  ## Estimate principal components of the between covariance (Kb)
  ##################################################################################
  if(silent == FALSE) print("Estimate principal components of the between covariance (Kb)")

  temp = smooth.Gt$decom - smooth.Gw$decom
  Eigen <- eigen(temp,symmetric=TRUE)
  Sigma <- Eigen$values
  d <- Sigma[1:c.p]
  d <- d[d>0]
  per <- cumsum(d)/sum(d)
  N.Gb <- ifelse (is.null(npc), min(which(per>=pve)), min(npc, length(d)))
  smooth.Gb <- list(evalues=Sigma[1:N.Gb], efunctions=An%*%Eigen$vectors[,1:N.Gb], decom=temp)
  pve1  <- per[N.Gb]
 
  ###########################################################################################
  ## Estimate eigenvalues and eigenfunctions at two levels by calling the
  ## eigenfunction (in R "base" package) on discretized covariance matrices.
  ###########################################################################################
  if(silent == FALSE) print("Estimate eigenvalues and eigenfunctions at two levels")

  ddpos = function(dd){
    e1 = eigen(dd)
    if(sum(e1$values<=0)>0){
      s1 = e1$values
      u1 = e1$vectors
      u1 = u1[,s1>0]
      s1 = s1[s1>0]
      dpos = u1%*%diag(s1)%*%t(u1)
    }else{
      dpos = dd
    }
    return(dpos)
  }

  efunctions <- list(level1=as.matrix(smooth.Gb$efunctions), level2=as.matrix(smooth.Gw$efunctions))
  evalues <- list(level1=smooth.Gb$evalues,level2=smooth.Gw$evalues)
  decom <- list(level1=ddpos(smooth.Gb$decom),level2=ddpos(smooth.Gw$decom),
    levelall=ddpos(smooth.Gt$decom))
  npc <- list(level1=length(evalues[[1]]), level2=length(evalues[[2]]))
  pve <- list(level1 = pve1, level2 = pve2)
  lambda_face <- list(lambda.Gt=smooth.Gt$lambda,lambda.Gw=smooth.Gw$lambda)


  ###################################################################
  # Estimate the principal component scores
  ###################################################################
  if(silent == FALSE) print("Estimate principal component scores")
    cumday = c(0,cumsum(Ji))
    score.x = matrix(0,length(Ji),ncol(efunctions[[1]]))
    score.u = matrix(0,sum(Ji),ncol(efunctions[[2]]))
    Xhat = Xhat.subject = matrix(0,sum(Ji),ncol(df$Ytilde))
    for (i in 1:length(Ji)) {
        Yi = df$Ytilde[(cumday[i]+1):cumday[i+1],]
        if(Ji[i]>1){
            cYi = colMeans(Yi)
        }else{
            cYi = Yi
        }
        score.x[i,] = as.vector(diag(1/sqrt(evalues[[1]]))%*%t(efunctions[[1]])%*%cYi)
        Xhat.subject[(cumday[i]+1):cumday[i+1],] = matrix(efunctions[[1]]%*%diag(sqrt(evalues[[1]]))%*%score.x[i,],Ji[i],ncol(df$Ytilde),byrow=T)
        temp = matrix(Yi - Xhat.subject[(cumday[i]+1):cumday[i+1],],Ji[i],ncol(df$Ytilde)) # nday*p
        score.u[(cumday[i]+1):cumday[i+1],] = t(diag(1/sqrt(evalues[[2]]))%*%t(efunctions[[2]])%*%t(temp))
        Xhat[(cumday[i]+1):cumday[i+1],] = Xhat.subject[(cumday[i]+1):cumday[i+1],]+score.u[(cumday[i]+1):cumday[i+1],]%*%diag(sqrt(evalues[[2]]))%*%t(efunctions[[2]])
    }
    scores <- list(level1 = score.x, level2 = score.u)
  
  
  ###################################################################
  # Organize the results
  ###################################################################
  if(silent == FALSE) print("Organize the results")

  res <- list(Yhat = Xhat,Y.df = df$Y,Thetamat=decom,funbasis=An,lambda_face=lambda_face)

  return(res)
}
