INTACT <- function(mydf0,Yraw,formula,needlog=0,comparisons=NULL,act=1,wpve=0.5,dd=NULL,
    loga0=seq(-1.5,1.5,0.5),cw=.01,ro=c(10,.9),itermax=100,nfold=1,endpara=1e-3,th.value=0.01){
    matnorm2 <- function(X){
      X <- scale(X, center = TRUE, scale = FALSE)
        if(sum(is.na(X))>0){       
            fit <- softImpute(X, rank.max = 1, lambda = 0) 
            sigma1 <- fit$d[1]      
        }else{
            sigma1 <- irlba(X, nv = 1, nu = 0)$d[1]
        }
        norm2 <- sigma1^2 / (nrow(X) - 1)
        return(norm2)
    }
    XWtoY <- function(X,W,daynum){
        subnum <- length(daynum)
        obsnum <- nrow(W)
        p <- ncol(W)
        cumday <- c(0,cumsum(daynum))
        Y <- matrix(0,obsnum,p)
        for(sub in 1:subnum){
          subrow <- (cumday[sub]+1):cumday[sub+1]
          Y[subrow,] <- matrix(X[sub,],daynum[sub],p,byrow = TRUE) + 
            W[subrow,]
        }
        return(Y)
    }

    start_time <- proc.time()
    mydf0$batch <- as.factor(mydf0$batch)
    daynumall <- table(mydf0$subid)[unique(mydf0$subid)]
    if(is.null(colnames(Yraw))){
      colnames(Yraw) <- paste("t",1:ncol(Yraw),sep="")
    }
    if(needlog==1){
        Yraw <- log(Yraw+1)
    }
    M <- length(unique(mydf0$batch))
    Yrawm <- featuresm <- daynum <- vector("list",M)
    yn2_sq <- rep(0,M)
    for(i in 1:M){
        tmp = which(mydf0$batch==i-1)
        Yrawm[[i]] <- Yraw[tmp,]
        yn2_sq[i] <- sqrt(matnorm2(Yrawm[[i]]))
        featuresm[[i]] <- Yrawm[[i]]/yn2_sq[i]
        daynum[[i]] <- table(mydf0[tmp,]$subid)[unique(mydf0[tmp,]$subid)]
    }
    features <- do.call(rbind, featuresm)
    multipar <- sqrt(matnorm2(Yraw)/matnorm2(features))
    mydf <- data.frame(mydf0,features)
    featurenames <- colnames(features)
    intactXW <- INTACTforXW(mydf,featurenames,daynum,formula,loga0,cw,ro,itermax,nfold,dd,endpara,th.value,wpve)
    Xint <- intactXW$Xint
    Wint <- intactXW$Wint
    Sx <- intactXW$Sx
    Sw <- intactXW$Sw
    eps <- intactXW$eps
    predicted <- intactXW$predicted
    funbasis <- intactXW$funbasis
    
    Eta_intact0 <- Eta_intact <- errint <- vector("list",M)
    Y_intact0 <- Y_intact <- vector("list",M)
    norm2Yxm <- norm2Ywm <- norm2errm <- rep(0,M)
    mSx <- mSw <- 0
    epsall <- NULL
    for(m in 1:M){
        mSx <- mSx + Sx[[m]]
        mSw <- mSw + Sw[[m]]
        epsall <- rbind(epsall,scale(eps[[m]],center = TRUE, scale = FALSE)/sqrt(nrow(eps[[m]])-1))
        norm2Yxm[m] <- matnorm2(Xint[[m]])
        norm2Ywm[m] <- matnorm2(Wint[[m]])
        norm2errm[m] <- matnorm2(eps[[m]])
    }
    mSx <- mSx/M
    mSw <- mSw/M
    norm2Yx <- irlba(mSx, nv = 1, nu = 0)$d[1]
    norm2Yw <- irlba(mSw, nv = 1, nu = 0)$d[1]
    norm2err <- (irlba(as.matrix(epsall), nv = 1, nu = 0)$d[1])^2/M
    c1 <- sqrt(norm2Yx/norm2Yxm)
    c2 <- sqrt(norm2Yw/norm2Ywm)
    c3 <- sqrt(norm2err/norm2errm)
    for(m in 1:M){
        Xint[[m]] <- c1[m]*Xint[[m]]
        Wint[[m]] <- c2[m]*Wint[[m]]
        errint[[m]] <- c3[m]*eps[[m]]
    }
    for(m in 1:M){
        Eta_intact0[[m]] <- XWtoY(Xint[[m]],Wint[[m]],daynum[[m]])%*%t(funbasis)
        Eta_intact[[m]] <- as.matrix(Eta_intact0[[m]]) + errint[[m]]
       
        Y_intact[[m]] <- Eta_intact[[m]] + predicted[[m]]
        Y_intact[[m]] <- Y_intact[[m]]*multipar
        Y_intact0[[m]] <- Eta_intact0[[m]] + predicted[[m]]
        Y_intact0[[m]] <- Y_intact0[[m]]*multipar
        if(act==1){
            Y_intact[[m]][Y_intact[[m]]<0] <- 0
            Y_intact0[[m]][Y_intact0[[m]]<0] <- 0
        }
    }
    for(i in 1:M){
        Y_intact[[i]] <- as.matrix(Y_intact[[i]])
        Y_intact0[[i]] <- as.matrix(Y_intact0[[i]])
        Eta_intact[[i]] <- as.matrix(Eta_intact[[i]])
        Eta_intact0[[i]] <- as.matrix(Eta_intact0[[i]])
        rownames(Y_intact[[i]]) <- rownames(Y_intact0[[i]]) <- rownames(Eta_intact[[i]]) <- rownames(Eta_intact0[[i]]) <- rownames(featuresm[[i]])
    }
    end_time <- proc.time()
    elapsed <- end_time - start_time

    time.intact <- intactXW$elapsed.int
    time.intact[[4]] <- elapsed['elapsed']
    names(time.intact)[4] <- 'elapsed.total'

    reintact <- list(Y=Y_intact,Eta=Eta_intact,multipar=multipar,elapsed=time.intact)
    reintact0 <- list(Y=Y_intact0,Eta=Eta_intact0)
    reall <- list(reintact=reintact,reintact0=reintact0,
        d=intactXW$d,Xi.x=intactXW$Xi.x,Xi.w=intactXW$Xi.w,fixeffs=intactXW$fixeffs,batch_effects_adjusted=intactXW$batch_effects_adjusted)

    if(comparisons==1){
        ### raw
        start_time.raw <- proc.time()
        Y_raw <- Eta_raw <- vector("list",M)
        for (m in 1:M) {
            Eta_raw[[m]] <- as.matrix(intactXW$Eta[[m]]) 
            if(sum(is.na(featuresm[[m]]))>0){
                id <- mydf$subid[which(mydf$batch==m-1)]
                featuresm[[m]] <- fastmfpca(Y = as.matrix(featuresm[[m]]), id = id,twoway=FALSE,pve=1)$Y.df
            }
            Y_raw[[m]] <- as.matrix(featuresm[[m]])*multipar 
        }
        features <- do.call(rbind, featuresm)
        mydf <- data.frame(mydf0,features)
        for(i in 1:M){
            Y_raw[[i]] <- as.matrix(Y_raw[[i]])
            Eta_raw[[i]] <- as.matrix(Eta_raw[[i]])
            rownames(Y_raw[[i]]) <- rownames(Eta_raw[[i]]) <- rownames(featuresm[[i]])
        }
        end_time.raw <- proc.time()
        elapsed.raw <- end_time.raw - start_time.raw
        reraw <- list(Y=Y_raw,Eta=Eta_raw,elapsed=elapsed.raw["elapsed"])

        ### simple standardization
        start_time.sim <- proc.time()
        p <- ncol(Eta_raw[[1]])
        sig <- del <- matrix(0,M,p)
        Y_sim <- Eta_sim <- vector("list",M)
        sigerr <- apply(do.call(rbind, Eta_raw),2,var)
        for (m in 1:M) {
            sig[m,] <- apply(Eta_raw[[m]],2,var)
            del[m,] <- sqrt(sig[m,]/sigerr)
            Eta_sim[[m]] <- Eta_raw[[m]]/matrix(del[m,],nrow(Eta_raw[[m]]),p,byrow=T)
            Y_sim[[m]] <- Eta_sim[[m]] + predicted[[m]]
            Y_sim[[m]] <- Y_sim[[m]]*multipar
            if(act==1){
                Y_sim[[m]][Y_sim[[m]]<0] <- 0
            }
        }
        for(i in 1:M){
            Y_sim[[i]] <- as.matrix(Y_sim[[i]])
            Eta_sim[[i]] <- as.matrix(Eta_sim[[i]])
            rownames(Y_sim[[i]]) <- rownames(Eta_sim[[i]]) <- rownames(featuresm[[i]])
        }
        end_time.sim <- proc.time()
        elapsed.sim <- end_time.sim - start_time.sim
        elapsed.sim <- elapsed.sim["elapsed"]+time.intact[[1]]
        resim <- list(Y=Y_sim,Eta=Eta_sim,elapsed=elapsed.sim)

        ### quantile normalization
        start_time.quan <- proc.time()
        Y_quan <- vector("list",M)
        Eta_quan <- Eta_raw
        ni <- sapply(Eta_raw, nrow)
        for(pid in 1:p){
            aligned_mat <- lapply(Eta_raw, function(mat) mat[, pid])
            maxni <- max(ni)
            common_grid <- seq(0, 1, length.out = maxni)
            aligned_mat <- sapply(aligned_mat, function(x) {
                    approx(seq(0, 1, length.out = length(x)), x, xout = common_grid)$y
                })
            temp <- normalize(aligned_mat)
            Eta_quan_pid <- lapply(1:ncol(temp), function(j) {
                grid_out <- seq(0, 1, length.out = ni[j])
                approx(x = seq(0, 1, length.out = nrow(temp)), y = temp[,j], xout = grid_out)$y
            })
            for (i in 1:M) {
                Eta_quan[[i]][,pid] = Eta_quan_pid[[i]]
            }
        }
        for (i in 1:M) {
            Y_quan[[i]] <- Eta_quan[[i]] + predicted[[i]]
            Y_quan[[i]] <- Y_quan[[i]]*multipar
            if(act==1){
                Y_quan[[i]][Y_quan[[i]]<0] <- 0
            }
        }
        for(i in 1:M){
            Y_quan[[i]] <- as.matrix(Y_quan[[i]])
            Eta_quan[[i]] <- as.matrix(Eta_quan[[i]])
            rownames(Y_quan[[i]]) <- rownames(Eta_quan[[i]]) <- rownames(featuresm[[i]])
        }
        end_time.quan <- proc.time()
        elapsed.quan <- end_time.quan - start_time.quan
        elapsed.quan <- elapsed.quan["elapsed"]+time.intact[[1]]
        requan <- list(Y=Y_quan,Eta=Eta_quan,elapsed=elapsed.quan)

        ### longitudinal ComBat
        start_time.lcombat <- proc.time()
        mydf$time <- unlist(sapply(daynumall,function(x) 1:x))-1
        Ycombat <- longCombat(idvar='subid',timevar='time',batchvar='batch',features=colnames(features), 
            formula=formula,ranef='(1|subid)',data=mydf,verbose=FALSE)
        Y_lcombat <- Eta_lcombat <- vector("list",M)
        for (i in 1:M) {
            Y_lcombat[[i]] <- Ycombat$data_combat[which(mydf0$batch==i-1),-(1:3)] 
            Y_lcombat[[i]] <- as.matrix(Y_lcombat[[i]])
            Eta_lcombat[[i]] <- Y_lcombat[[i]] - predicted[[i]]
            if(act==1){
                Y_lcombat[[i]][Y_lcombat[[i]]<0] <- 0
            }
            Y_lcombat[[i]] <- Y_lcombat[[i]]*multipar
        }
        for(i in 1:M){
            Y_lcombat[[i]] <- as.matrix(Y_lcombat[[i]])
            Eta_lcombat[[i]] <- as.matrix(Eta_lcombat[[i]])
            rownames(Y_lcombat[[i]]) <- rownames(Eta_lcombat[[i]]) <- rownames(featuresm[[i]])
        }
        end_time.lcombat <- proc.time()
        elapsed.lcombat <- end_time.lcombat - start_time.lcombat
        relcombat <- list(Y=Y_lcombat,Eta=Eta_lcombat,elapsed=elapsed.lcombat["elapsed"])
        
        reall <- list(reintact=reintact,reintact0=reintact0,reraw=reraw,resim=resim,requan=requan,relcombat=relcombat,
            d=intactXW$d,Xi.x=intactXW$Xi.x,Xi.w=intactXW$Xi.w,fixeffs=intactXW$fixeffs,batch_effects_adjusted=intactXW$batch_effects_adjusted)
    }

    return(reall)
}

INTACTforXW <- function(mydf,featurenames,daynum,formula,loga0=seq(-1.5,1.5,0.5),cw=.01,ro=c(10,.9),itermax=100,nfold=1,dd=NULL,endpara=1e-3,th.value=0.01,wpve=0.5){
    
    preadj <- function(data,formula,featurenames,batchvar,ranef,ifsm=1){
    # make batch a factor if not already
    data[[batchvar]] <- as.factor(data[[batchvar]])
    batch <- droplevels(data[[batchvar]])
    M <- nlevels(batch)
    batches <- lapply(levels(batch), function(x) which(batch==x))
    ni <- sapply(batches, length)
    V <- length(featurenames)
    L <- nrow(data)

    varname <- unlist(strsplit(formula,"[+]"))
    fixnum <- length(varname) + M
    fixeffs <- matrix(nrow=fixnum, ncol=V)
    batchnames <- paste(batchvar,1:(M-1),sep="")
    rownames(fixeffs) <- c("intercept",varname,batchnames)
    
    for(v in 1:V){
        lme_formula <- as.formula(paste0(featurenames[v], '~', formula, '+' , batchvar, '+', ranef))
        lme_fit <- lmer(lme_formula, data=data, REML=TRUE, control=lmerControl(optimizer='bobyqa'))
         # save fixed effects
        fixeffs[,v] <- fixef(lme_fit)
    }

    if(ifsm == 1){
        # smooth fixed effects
        argvals <- seq(0, 1, length.out=V) 
        for(ff in 1:fixnum){
            fit_fixeffs <- gam(fixeffs[ff,] ~ s(argvals))
            fixeffs[ff,] <- as.vector(predict(fit_fixeffs, newdata = data.frame(argvals = argvals)))
        }
    }

    # smooth batch effects
    batch_effects <- matrix(fixeffs[-(1:(1+length(varname))),],M-1,V)

    gamma1hat <- -(ni[2:M] %*% batch_effects)/L
  
    batch_effects_adjusted <- sweep(batch_effects, 2, gamma1hat, FUN='+')
    batch_effects_adjusted <- rbind(gamma1hat, batch_effects_adjusted)
    batch_effects_expanded <- matrix(nrow=L, ncol=V)
    for(i in 1:M){ # begin loop over batches
        batch_effects_expanded[batches[[i]],] <- matrix(
      rep(batch_effects_adjusted[i,],length(batches[[i]])),
      ncol=V, byrow=TRUE) 
    } # end loop over batches

    intercept_adjusted <- fixeffs[1,] - gamma1hat
    
    return(list(fixeffs=fixeffs,batch_effects_expanded=batch_effects_expanded,
        batch_effects_adjusted=batch_effects_adjusted,intercept_adjusted=intercept_adjusted))
    }

    newscore.N2 <- function(YY,dm,evalues.l1,evalues.l2,efunctions.l1,efunctions.l2){
        cumday <- c(0,cumsum(dm))
        score.x <- matrix(0,length(dm),ncol(efunctions.l1))
        score.u <- matrix(0,sum(dm),ncol(efunctions.l2))

        for (i in 1:length(dm)) {
            Yi <- YY[(cumday[i]+1):cumday[i+1],]
            if(dm[i]>1){
                cYi <- colMeans(Yi)
            }else{
                cYi <- Yi
            }
            score.x[i,] <- as.vector(diag(1/sqrt(evalues.l1),length(evalues.l1),length(evalues.l1))%*%t(efunctions.l1)%*%cYi)
            temp <- Yi - matrix(efunctions.l1%*%diag(sqrt(evalues.l1),length(evalues.l1),length(evalues.l1))%*%score.x[i,],dm[i],ncol(YY),byrow=T) # nday*p
            score.u[(cumday[i]+1):cumday[i+1],] <- as.matrix(t(diag(1/sqrt(evalues.l2),length(evalues.l2),length(evalues.l2))%*%t(efunctions.l2)%*%t(temp)))
        }

        return(list(score.x=score.x,score.u=score.u))
    }

    estintact <- function(S,d,funbasis,alpha1,alpha2,cw,ro,itermax,endpara=1e-3){
    initial <- function(S,d){
        M <- length(S)
        Smean <- 0
        for (j in 1:M) {
            Smean <- Smean + S[[j]]
        }
        Smean <- Smean/M
        eigS <- eigen(Smean)
        Xi <- as.matrix(eigS$vectors[,1:d])
        Lam <- diag(eigS$values[1:d],d,d)
        D <- U <- vector("list",M)
        sumdiagD <- 0
        matE <- diag(1,d,d)

        Xim.ini <- vector("list",M)
        for (j in 1:M) {
            eigSj <- eigen(S[[j]])
            Xim.ini[[j]] <- as.matrix(eigSj$vectors[,1:d])
            for(k in 1:d){
                if(sum(Xi[,k]*Xim.ini[[j]][,k])<0){
                    Xim.ini[[j]][,k] <- -Xim.ini[[j]][,k]
                }
            }
            temp <- t(Xi)%*%Xim.ini[[j]]
            svdtemp <- svd(temp)
            U[[j]] <- svdtemp$u%*%t(svdtemp$v)
            for(k in 1:d){
                if(sum((U[[j]][,k]-matE[,k])^2)>sum((-U[[j]][,k]-matE[,k])^2)){
                    U[[j]][,k] <- -U[[j]][,k]
                }
            }

            D[[j]] <- diag(eigSj$values[1:d]/diag(Lam),d,d)
            sumdiagD <- sumdiagD + D[[j]][1,1]/M
        }

        b <- rep(0,M)
        for (j in 1:M) {
            D[[j]] <- D[[j]]/sumdiagD
            b[j] <- D[[j]][1,1]
        }
        Lam <- Lam*sumdiagD
        lam <- diag(Lam)
        lam <- lam/lam[d]
        nv <- lam

        return(list(Xi=Xi,Lam=Lam,D=D,U=U,b=b,lam=lam,nv=nv,Xim.ini=Xim.ini))
    }

    M <- length(S)
    aa <- rep(1/M,M)

    # initialization ...
    ini <- initial(S,d)
    Xi <- ini$Xi
    Lam <- ini$Lam
    D <- ini$D
    U <- ini$U
    b <- ini$b
    lam <- ini$lam*alpha1
    nv <- ini$nv*alpha2
    p <- nrow(funbasis)
    D2 <- diff(diag(p), differences = 2)
    D2B <- D2%*%funbasis # ((p-2)*p)*(p*q)=(p-2)*q
    B0 <- t(D2B)%*%D2B
    omega <- cw/sum(diag(t(Xi)%*%B0%*%Xi))

    crho <- ro[1]
    brho <- ro[2]
    for (i in 1:itermax) {
        rho <- crho*brho^i
        bold <- b
        # Xi-update
        A <- 0
        for(j in 1:M){
            A <- A + S[[j]]%*%Xi%*%U[[j]]%*%Lam%*%D[[j]]%*%t(U[[j]])
        }
        A <- A - omega*B0%*%Xi/2 + Xi*rho/2 
        svdA <- svd(A)
        Xi <- svdA$u%*%t(svdA$v)

        # Lambda-update
        R <- E <- 0
        for (j in 1:M) {
            R <- R + D[[j]]%*%t(U[[j]])%*%t(Xi)%*%S[[j]]%*%Xi%*%U[[j]]
            E <- E + D[[j]]^2
        }
        R <- R + rho*Lam
        E <- E + rho*diag(1,d,d)
        tLam <- diag(R)/diag(E)
        Lam <- diag(pmax(tLam,0),d,d)

        # D-update
        Fb <- diag(M)*sum(nv)
        psi <- rep(0,M)
        for (j in 1:M) {
            F <- t(U[[j]])%*%t(Xi)%*%S[[j]]%*%Xi%*%U[[j]]%*%Lam+b[j]*diag(nv,d,d)+rho*D[[j]]
            G <- Lam^2 + diag(nv,d,d) + rho*diag(1,d,d)
            tD <- diag(F)/diag(G)
            D[[j]] <- diag(pmax(tD,0),d,d)
            psi[j] <- sum(diag(nv,d,d)*D[[j]])
        }

        # b-update
        quad_D <- Fb+rho*diag(M)
        quad_d <- psi+rho*b
        quad_A <- t(rbind(aa,diag(M)))
        quad_b0 <- c(1,rep(0,M))
        sc <- norm(quad_D,"2")
        b <- solve.QP(quad_D/sc, quad_d/sc, quad_A, quad_b0, meq=1, factorized=FALSE)$solution

        # U-update
        for(j in 1:M){
            H <- t(Xi)%*%S[[j]]%*%Xi%*%U[[j]]%*%Lam%*%D[[j]]+U[[j]]*rho/2+diag(lam,d,d)/2
            svdH <- svd(H)
            U[[j]] <- svdH$u%*%t(svdH$v)
        }

        if(mean((b-bold)^2)<endpara){
            break
        }
    }
    
    return(list(Xi=Xi,Lam=Lam,D=D,U=U,b=b,ini=ini,iter=i))
    }

    cvest <- function(mydf,featurenames,daynum,formula,funbasis,d,loga0,cw,ro,itermax,nfold,endpara){
    estintact.UD <- function(S,d,ini,Xi,Lam,funbasis,alpha1,alpha2,ro,itermax,endpara=1e-3){
        M <- length(S)
        aa <- rep(1/M,M)
        # initialization ...
        D <- ini$D
        U <- ini$U
        b <- ini$b
        lam <- ini$lam*alpha1
        nv <- ini$nv*alpha2

        p <- nrow(funbasis)
        D2 <- diff(diag(p), differences <- 2)
        D2B <- D2%*%funbasis # ((p-2)*p)*(p*q)=(p-2)*q
        B0 <- t(D2B)%*%D2B

        crho <- ro[1]
        brho <- ro[2]
        for (i in 1:itermax) {
            rho <- crho*brho^i
            bold <- b

            # D-update
            Fb <- diag(M)*sum(nv)
            psi <- rep(0,M)
            for (j in 1:M) {
                F <- t(U[[j]])%*%t(Xi)%*%S[[j]]%*%Xi%*%U[[j]]%*%Lam+b[j]*diag(nv,d,d)+rho*D[[j]]
                G <- Lam^2 + diag(nv,d,d) + rho*diag(1,d,d)
                tD <- diag(F)/diag(G)
                D[[j]] <- diag(pmax(tD,0),d,d)
                psi[j] <- sum(diag(nv,d,d)*D[[j]])
            }

            # b-update
            quad_D <- Fb+rho*diag(M)
            quad_d <- psi+rho*b
            quad_A <- t(rbind(aa,diag(M)))
            quad_b0 <- c(1,rep(0,M))
            sc <- norm(quad_D,"2")
            b <- solve.QP(quad_D/sc, quad_d/sc, quad_A, quad_b0, meq=1, factorized=FALSE)$solution

            # U-update
            for(j in 1:M){
                H <- t(Xi)%*%S[[j]]%*%Xi%*%U[[j]]%*%Lam%*%D[[j]]+U[[j]]*rho/2+diag(lam,d,d)/2
                svdH <- svd(H)
                U[[j]] <- svdH$u%*%t(svdH$v)
            }

            if(mean((b-bold)^2)<endpara){
                break
            }
        }
        return(list(Xi=Xi,Lam=Lam,D=D,U=U,b=b))
    }

    EstTheta <- function(fixeffs,batch_effects_adjusted,dfuse,dfleft,formula,featurenames){
        M <- length(dfuse)
        if(!is.null(formula)){
            varname <- unlist(strsplit(formula,"[+]"))
        }

        Etause <- Etaleft <- vector("list",M)
        Thetause.l1 <- Thetause.l2 <-  Thetaleft.l1 <- Thetaleft.l2 <-  vector("list",M)
        for (i in 1:M) {
            if(is.null(formula)){
                myvab <- matrix(1,nrow(dfuse[[i]]),1)
            }else{
                myvab <- cbind(1,matrix(as.numeric(as.matrix(dfuse[[i]][,varname])),nrow(dfuse[[i]]),length(varname)))
            }
            predicted <- myvab%*%fixeffs
            Etause[[i]] <- dfuse[[i]][,featurenames] - predicted - matrix(batch_effects_adjusted[i,],nrow(predicted),ncol(predicted),byrow=T)
            if(is.null(formula)){
                myvab <- matrix(1,nrow(dfleft[[i]]),1)
            }else{
                myvab <- cbind(1,matrix(as.numeric(as.matrix(dfleft[[i]][,varname])),nrow(dfleft[[i]]),length(varname)))
            }
            predicted <- myvab%*%fixeffs
            Etaleft[[i]] <- dfleft[[i]][,featurenames] - predicted - matrix(batch_effects_adjusted[i,],nrow(predicted),ncol(predicted),byrow=T)
            
            fpcause <- fastmfpca(Y = as.matrix(Etause[[i]]), id = dfuse[[i]]$subid,twoway=FALSE,pve=1)
            Thetause.l1[[i]] <- fpcause$Thetamat[[1]]
            Thetause.l2[[i]] <- fpcause$Thetamat[[2]]
            lambdause <- fpcause$lambda_face
            fpcaleft <- fastmfpca(Y = as.matrix(Etaleft[[i]]), id = dfleft[[i]]$subid,twoway=FALSE,pve=1,
                lambda.Gt=lambdause$lambda.Gt,lambda.Gw=lambdause$lambda.Gw)
            Thetaleft.l1[[i]] <- fpcaleft$Thetamat[[1]]
            Thetaleft.l2[[i]] <- fpcaleft$Thetamat[[2]]
            lambdaleft <- fpcaleft$lambda_face
        }

        return(list(Thetause.l1=Thetause.l1,Thetause.l2=Thetause.l2,
            Thetaleft.l1=Thetaleft.l1,Thetaleft.l2=Thetaleft.l2,
            lambdaleft=lambdaleft,lambdause=lambdause))
    }

    M <- length(unique(mydf$batch))

    Thetause.l1 <- Thetause.l2 <- Thetaleft.l1 <- Thetaleft.l2 <- vector("list",M)

    alpha <- exp(loga0)
    res.x <- res.w <- matrix(0,length(alpha),nfold)

    subnum <- nleft <- rep(0,M)
    cumday <- dfM <- vector("list",M)
    for(j in 1:M){
        subnum[j] <- length(daynum[[j]]) 
        cumday[[j]] <- c(0,cumsum(daynum[[j]]))
        dfM[[j]] <- mydf[which(mydf$batch==j-1),]
        if(nfold>1){
            nleft[j] <- round(subnum[j]/nfold) 
        }else{ # nfold=1
            nleft[j] <- round(subnum[j]/10) 
        }
    }

    for(k in 1:nfold){   
        dfuse <- dfleft <- vector("list",M) 
        for (j in 1:M) {
            if(nfold>1){
                if(k!=nfold){
                    leftsub <- (1+(k-1)*nleft[j]):(k*nleft[j])
                }else{
                    leftsub <- (1+(k-1)*nleft[j]):subnum[j]
                }
            }else{# k=nfold=1
                leftsub <- 1:nleft[j]
            }
            leftrow <- (cumday[[j]][leftsub[1]]+1):cumday[[j]][leftsub[length(leftsub)]+1]

            dfuse[[j]] <- dfM[[j]][-leftrow,]
            dfleft[[j]] <- dfM[[j]][leftrow,]
        }

        combined_dfuse <- do.call(rbind, dfuse)
        if(is.null(formula)){
            batch_effects_adjusted <- matrix(0,M,length(featurenames))
            Ni <- rep(0,M)
            for(i in 1:M){
                Ni[i] <- nrow(dfuse[[i]])
                batch_effects_adjusted[i,] <- colMeans(dfuse[[i]][,featurenames],na.rm=T)
            }
            fixeffs <- t(Ni)%*%batch_effects_adjusted/sum(Ni)
            for (i in 1:M) {
                batch_effects_adjusted[i,] <- batch_effects_adjusted[i,]-fixeffs
            }
            fixeffs <- matrix(fixeffs,1,length(fixeffs),byrow=T)           
        }else{
            mypre <- preadj(data=combined_dfuse,formula=formula,featurenames=featurenames,
                batchvar='batch',ranef='(1|subid)',ifsm=1)
            fixeffs <- mypre$fixeffs # dim*T
            fixeffs <- fixeffs[2:(nrow(fixeffs)-M+1),]
            intercept_adjusted <- mypre$intercept_adjusted
            fixeffs <- rbind(intercept_adjusted,fixeffs)
            batch_effects_adjusted <- mypre$batch_effects_adjusted
        }
        Theta <- EstTheta(fixeffs,batch_effects_adjusted,dfuse,dfleft,formula,featurenames=featurenames)
        Thetause.l1 <- Theta$Thetause.l1
        Thetause.l2 <- Theta$Thetause.l2
        Thetaleft.l1 <- Theta$Thetaleft.l1
        Thetaleft.l2 <- Theta$Thetaleft.l2

        for(ai in 1:length(alpha)){
            alpha1 <- alpha[ai]
            alpha2 <- alpha[ai]
            re.x <- estintact(Thetause.l1,d[1],funbasis,alpha1,alpha2,cw,ro,itermax,endpara)
            re.w <- estintact(Thetause.l2,d[2],funbasis,alpha1,alpha2,cw,ro,itermax,endpara)

            Xi.x <- re.x$Xi
            Lam.x <- re.x$Lam
            ini.x <- re.x$ini
            releft.x <- estintact.UD(Thetaleft.l1,d[1],ini.x,Xi.x,Lam.x,funbasis,alpha1,alpha2,ro,itermax,endpara)
            Dleft.x <- releft.x$D
            Uleft.x <- releft.x$U

            Xi.w <- re.w$Xi
            Lam.w <- re.w$Lam
            ini.w <- re.w$ini
            releft.w <- estintact.UD(Thetaleft.l2,d[2],ini.w,Xi.w,Lam.w,funbasis,alpha1,alpha2,ro,itermax,endpara)
            Dleft.w <- releft.w$D
            Uleft.w <- releft.w$U

            for(j in 1:M){
                res.x[ai,k] <- res.x[ai,k] + sum((Thetaleft.l1[[j]]-Xi.x%*%Uleft.x[[j]]%*%Lam.x%*%Dleft.x[[j]]%*%t(Uleft.x[[j]])%*%t(Xi.x))^2)
                res.w[ai,k] <- res.w[ai,k] + sum((Thetaleft.l2[[j]]-Xi.w%*%Uleft.w[[j]]%*%Lam.w%*%Dleft.w[[j]]%*%t(Uleft.w[[j]])%*%t(Xi.w))^2)
            }
        }
    }

    return(list(res.x = res.x, res.w = res.w,Theta=Theta))
    }

    tag <- 0
    start_time1 <- proc.time()
    M <- length(unique(mydf$batch))
    mydf$batch <- as.factor(mydf$batch)

    if(is.null(formula)){
        batch_effects_adjusted = matrix(0,M,length(featurenames))
        Ni = rep(0,M)
        argvals <- seq(0, 1, length.out=length(featurenames)) 
        for(i in 1:M){
            Ni[i] = length(which(mydf$batch==(i-1)))
            batch_effects_adjusted[i,] = colMeans(mydf[which(mydf$batch==(i-1)),featurenames],na.rm=T)

            fit_fixeffs <- gam(batch_effects_adjusted[i,] ~ s(argvals))
            batch_effects_adjusted[i,] <- as.vector(predict(fit_fixeffs, newdata = data.frame(argvals = argvals)))
        }
        fixeffs <- t(Ni)%*%batch_effects_adjusted/sum(Ni)
        for (i in 1:M) {
            batch_effects_adjusted[i,] <- batch_effects_adjusted[i,]-fixeffs
        }
        fixeffs <- matrix(fixeffs,1,length(fixeffs),byrow=T)           
    }else{
        varname <- unlist(strsplit(formula,"[+]"))
        mypre <- preadj(data=mydf,formula=formula,featurenames,
                batchvar='batch',ranef='(1|subid)',ifsm=1)
        fixeffs <- mypre$fixeffs 
        fixeffs <- fixeffs[2:(nrow(fixeffs)-M+1),]
        intercept_adjusted <- mypre$intercept_adjusted
        fixeffs <- rbind(intercept_adjusted,fixeffs)
        batch_effects_adjusted <- mypre$batch_effects_adjusted
    }

    dfM <- Eta <- predicted <- vector("list",M)
    for(i in 1:M){
        dfM[[i]] <- mydf[which(mydf$batch==i-1),]
        if(is.null(formula)){
            myvab <- matrix(1,nrow(dfM[[i]]),1)
        }else{
            myvab <- cbind(1,matrix(as.numeric(as.matrix(dfM[[i]][,varname])),nrow(dfM[[i]]),length(varname)))
        }
        predicted[[i]] <- myvab%*%fixeffs
        Eta[[i]] <- dfM[[i]][,featurenames] - predicted[[i]] - matrix(batch_effects_adjusted[i,],nrow(predicted[[i]]),ncol(predicted[[i]]),byrow=T)
    }

    eps <- Theta.l1 <- Theta.l2 <- vector("list",M)
    for (i in 1:M) {
        id <- mydf$subid[which(mydf$batch==i-1)]
        fpca1 <- fastmfpca(Y = as.matrix(Eta[[i]]), id = id,twoway=FALSE,pve=1)
        if(sum(is.na(Eta[[i]]))>0){
          Eta[[i]] <- fpca1$Y.df
        }
        eps[[i]] <- Eta[[i]] - fpca1$Yhat
        Theta.l1[[i]] <- fpca1$Thetamat[[1]]
        Theta.l2[[i]] <- fpca1$Thetamat[[2]]
    }
    
    funbasis <- as.matrix(fpca1$funbasis)
    efun.l1 <- evalue.l1 <- efun.l2 <- evalue.l2 <- vector("list",M)
    pve.l1 <- pve.l2 <- matrix(0,M,nrow(Theta.l1[[1]]))
    ev.l1 <- ev.l2 <- rep(0,M)
    for (i in 1:M) {
      temp <- eigen(Theta.l1[[i]])
      efun.l1[[i]] <- temp$vectors
      evalue.l1[[i]] <- temp$values
      evalue.l1[[i]][which(evalue.l1[[i]]<0)] <- 0
      pve.l1[i,] <- cumsum(evalue.l1[[i]])/sum(evalue.l1[[i]])
      if(length(which(evalue.l1[[i]]/sum(evalue.l1[[i]])>th.value))>0){
        ev.l1[i] <- max(which(evalue.l1[[i]]/sum(evalue.l1[[i]])>th.value))
      }else{
        ev.l1[i] <- 2
        tag <- 1
      }

      temp <- eigen(Theta.l2[[i]])
      efun.l2[[i]] <- temp$vectors
      evalue.l2[[i]] <- temp$values
      evalue.l2[[i]][which(evalue.l2[[i]]<0)] <- 0
      pve.l2[i,] <- cumsum(evalue.l2[[i]])/sum(evalue.l2[[i]])
      if(length(which(evalue.l2[[i]]/sum(evalue.l2[[i]])>th.value))>0){
        ev.l2[i] <- max(which(evalue.l2[[i]]/sum(evalue.l2[[i]])>th.value))
      }else{
        ev.l2[i] <- 2
        tag <- 1
      }
    }
    ev.l1.min <- min(ev.l1)
    ev.l2.min <- min(ev.l2)
    pve.l1 <- colMeans(pve.l1)
    pve.l2 <- colMeans(pve.l2)
    tempsim.l1 <- tempsim.l2 <- matrix(0,M*(M-1)/2,nrow(Theta.l1[[1]]))
    k <- 1
    for (i in 1:(M-1)) {
        for (j in (i+1):M) {
            tempsim.l1[k,] <- abs(colSums(efun.l1[[i]]*efun.l1[[j]]))
            tempsim.l2[k,] <- abs(colSums(efun.l2[[i]]*efun.l2[[j]]))
            k <- k+1
        }
    }
    sim.l1 <- apply(tempsim.l1,2,min)
    sim.l1 <- cummin(sim.l1)
    if(length(which(sim.l1>th.value))>0){
        kmax.l1 <- min(ev.l1.min,max(which(sim.l1>th.value)))
    }else{
        kmax.l1 <- min(ev.l1.min,2)
        tag <- 1
    }
    sim.l2 <- apply(tempsim.l2,2,min)
    sim.l2 <- cummin(sim.l2)
    if(length(which(sim.l2>th.value))>0){
        kmax.l2 <- min(ev.l2.min,max(which(sim.l2>th.value)))
    }else{
        kmax.l2 <- min(ev.l2.min,2)
        tag <- 1
    }
    
    npc <- which.max(((1-wpve)*sim.l1+wpve*pve.l1)[2:kmax.l1])+(2-1)
    npc2 <- which.max(((1-wpve)*sim.l2+wpve*pve.l2)[2:kmax.l2])+(2-1)
    
    d <- c(npc,npc2)
    if(!is.null(dd)){
        d <- dd
    }
    end_time1 <- proc.time()
    elapsed1 <- end_time1 - start_time1

    start_time2 <- proc.time()
    if(length(loga0)>1){
      re.cv <- cvest(mydf,featurenames,daynum,formula,funbasis,d,loga0,cw,ro,itermax,nfold,endpara)
      res.x <- re.cv$res.x
      res.w <- re.cv$res.w
      ai.x <- which.min(rowSums(res.x))
      ai.w <- which.min(rowSums(res.w))
    }else{
      re.cv <- NULL
      ai.x <- ai.w <- 1
    }
    alpha <- exp(loga0)
    end_time2 <- proc.time()
    elapsed2 <- end_time2 - start_time2

    start_time3 <- proc.time()
    alpha1 <- alpha[ai.x]
    alpha2 <- alpha[ai.x]
    result.x <- estintact(Theta.l1,d[1],funbasis,alpha1,alpha2,cw,ro,itermax,endpara)
    alpha1 <- alpha[ai.w]
    alpha2 <- alpha[ai.w]
    result.w <- estintact(Theta.l2,d[2],funbasis,alpha1,alpha2,cw,ro,itermax,endpara)

    Xim.ini.x <- result.x$ini$Xim.ini
    Xim.ini.w <- result.w$ini$Xim.ini

    score.l1 <- score.l2 <- vector("list",M)
    for(i in 1:M){
        rescore <- newscore.N2(YY=as.matrix(Eta[[i]]),dm=daynum[[i]],evalues.l1=evalue.l1[[i]][1:d[1]],
            evalues.l2=evalue.l2[[i]][1:d[2]],efunctions.l1=funbasis%*%Xim.ini.x[[i]],efunctions.l2=funbasis%*%Xim.ini.w[[i]])
        score.l1[[i]] <- rescore$score.x
        score.l2[[i]] <- rescore$score.u

        if(ncol(score.l1[[i]])<d[1]){
            score.l1[[i]] <- cbind(score.l1[[i]],matrix(0,nrow(score.l1[[i]]),d[1]-ncol(score.l1[[i]])))
        }
        if(ncol(score.l2[[i]])<d[2]){
            score.l2[[i]] <- cbind(score.l2[[i]],matrix(0,nrow(score.l2[[i]]),d[2]-ncol(score.l2[[i]])))
        }
    }
   
    Xi.x <- result.x$Xi
    Lam.x <- result.x$Lam
    Xi.w <- result.w$Xi
    Lam.w <- result.w$Lam

    Xint <- Wint <- vector("list",M)
    for (j in 1:M) {
        Xint[[j]] <- score.l1[[j]]%*%sqrt(Lam.x)%*%t(Xi.x) # subnum*q
        Wint[[j]] <- score.l2[[j]]%*%sqrt(Lam.w)%*%t(Xi.w) # obsnum*q
    }
    end_time3 <- proc.time()
    elapsed3 <- end_time3 - start_time3
    elapsed.int <- list(elapsed.fixedeff=elapsed1["elapsed"],elapsed.cv=elapsed2["elapsed"],elapsed.intact=elapsed3["elapsed"])
    Xi.x <- as.matrix(funbasis%*%Xi.x)
    Xi.w <- as.matrix(funbasis%*%Xi.w)
    intactXW = list(Xint=Xint,Wint=Wint,Sx=Theta.l1,Sw=Theta.l2,eps=eps,predicted=predicted,
        funbasis=funbasis,Eta=Eta,elapsed.int=elapsed.int,d=d,Xi.x=Xi.x,Xi.w=Xi.w,
        fixeffs=fixeffs,batch_effects_adjusted=batch_effects_adjusted)
    return(intactXW)
}

getMVDV <- function(Phixfun,Lamx,Phiwfun,Lamw,Y,multipar){
    d0 <- c(ncol(Phixfun),ncol(Phiwfun))
    var.l1.true <- Phixfun%*%Lamx%*%t(Phixfun)
    var.l2.true <- Phiwfun%*%Lamw%*%t(Phiwfun)
    M <- length(Y)
    evar.l1 <- evar.l2 <- rep(0,M)
    var.l1 <- var.l2 <- vector("list",M)
    for(m in 1:M){
        id <- rownames(Y[[m]])
        fpca <- fastmfpca(Y = as.matrix(Y[[m]]*multipar), id = id,twoway=FALSE,pve=1,npc=d0[1],npc2=d0[2])
        funbasis <- as.matrix(fpca$funbasis)
        var.l1[[m]] <- funbasis%*%fpca$Thetamat[[1]]%*%t(funbasis)
        var.l2[[m]] <- funbasis%*%fpca$Thetamat[[2]]%*%t(funbasis)
        evar.l1[m] <- sum((var.l1[[m]]-var.l1.true)^2)
        evar.l2[m] <- sum((var.l2[[m]]-var.l2.true)^2)
    }
    
    MV.X <- mean(evar.l1)
    MV.W <- mean(evar.l2)
    DV.X <- DV.W <- 0
    for(m1 in 1:(M-1)){
        for(m2 in (m1+1):M){
            DV.X <- DV.X + sum((var.l1[[m1]]-var.l1[[m2]])^2)
            DV.W <- DV.W + sum((var.l2[[m1]]-var.l2[[m2]])^2)
        }
    }
    DV.X <- DV.X/(M*(M-1)/2)
    DV.W <- DV.W/(M*(M-1)/2)

    return(list(MV.X=MV.X,MV.W=MV.W,DV.X=DV.X,DV.W=DV.W,var.X=var.l1,var.W=var.l2))
}


