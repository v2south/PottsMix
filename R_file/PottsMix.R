    ########################################################################### 
    # R code implementation for manuscript "The Analysis of Face Perception MEG 
    # and EEG Data Using a Potts-Mixture Spatiotemporal Joint Model" by  XYZ (anonymized). 
    # This R code will run the Iterative Conditional Mode algorithm. 
    ########################################################################### 
    # Last modified 20/09/2017 by XYZ (anonymized).
    
    # Get the current directory and set it to be the woring directory.
    rm(list = ls())
    path <- getwd()
    setwd(dir = path)
    
    library(rgl)
    library(R.matlab)
    library(scatterplot3d)
    library(MASS)
    library(PottsUtils)
    library(MCMCpack)
    library(bigRR)
    library(nnet)
    library(glmnet)
    library(fields)
    
    
    # load the Data Set
    # There are 4 different dataset
    
    # The first dataset Data1.RData contains simulated data with K=3 and two Gaussian signals.
    #load("./data/Data1/Data1.RData")
    
    # The second dataset Data3.RData contains simulated data with K=4 and three Gaussian signals. 
    # load("./data/Data2/Data2.RData")
    
    # The third dataset Data3.RData contains simulated data with K=4 and one sinusoid signal and two Gaussian signals. 
    #load("./data/Data3/Data3.RData")
    
    # The fourth dataset Data4.RData contains real data for application. 
    load("./data/Data4/Data4.RData")
  
    div <- TRUE 
    
    while (div)
    {
      # Load the vertices
      data <- readMat("./Data/vert.mat")
      vert.data <- data.frame(x = data$vert[,1],y = data$vert[,2],z = data$vert[,3])
     
      
      sub.vert <- vert.data
      op_sub.vert <- sub.vert
      
      # Built-in data 
      n_M  <- dim(Y_M)[1] #NUMBER OF MEG SENSORS
      n_E  <- dim(Y_E)[1] #NUMBER OF EEG SENSORS
      T <-  dim(Y_M)[2] #NUMBER OF TIME POINTS
      P <- dim(sub.vert)[1] #NUMBER OF BRAIN LOCATIONS
      K <- 10 #NUMBER OF MIXTURE COMPONENTS (MESO-STATES)
     
      if(sim){
      true.S <- S ; 
      }
      
      S <- NULL #STORE THE TRUE VALUES
      sigma2_E <- NULL
      sigma2_M <- NULL
     
      H_M <- diag(1, nrow = n_M, ncol = n_M)
      H_E <- diag(1, nrow = n_E, ncol = n_E)
      
      
      
      row.names(sub.vert) <- seq(1,P,1)
  
      
      
      # Number of brain locations before clustering.
      O_P <- dim(sub.vert)[1] 
      
      
      #Scale forward operator X for MEG and EEG.
      X_E <-  X_E / sqrt((1 / n_E)* sum(diag(X_E %*% t(X_E))))
      X_M <-  X_M / sqrt((1 / n_M)* sum(diag(X_M %*% t(X_M))))
      
      
      
      
      #Use Kmeans to cluster the vertices into 1000 clusters if number of vertices is greater than 1000.
      
      sub_cluster <- kmeans(sub.vert, centers = 250, iter.max = 100, nstart = 100)
      sub.vert <- sub_cluster$centers
      colnames(sub.vert) <- c("x", "y","z")
      sub.vert <- data.frame(sub.vert)
      P <- dim(sub.vert)[1] #Number of brain locations after clustering.
  
      CX_E <- matrix(NA, nrow = n_E, ncol = P)
      CX_M <- matrix(NA, nrow = n_M, ncol = P)
      
      for ( c in 1:P)
      {
        c_idx <- which( as.vector(sub_cluster$cluster) %in% c)
        if (length(c_idx) > 1){
          CX_M[,c] <- rowSums(X_M[,c_idx])
          CX_E[,c] <- rowSums(X_E[,c_idx])
        } else{
          CX_M[,c] <- X_M[,c_idx]
          CX_E[,c] <- X_E[,c_idx]
        }
        
      }
      
      X_M <- CX_M
      X_E <- CX_E
      
      sd_X <- apply(rbind(X_M, X_E), 2, sd)
      
      X_M  <-  X_M %*% diag(1/sd_X)  
      X_E  <-  X_E %*% diag(1/sd_X) 
      
      
      
      # Initial Number of voxels - this input will be modified
      n.v <- 400
      
      d.length <- ((range(sub.vert$x)[2] - range(sub.vert$x)[1])*(range(sub.vert$y)[2] - range(sub.vert$y)[1])*(range(sub.vert$z)[2] - range(sub.vert$z)[1]) / n.v)^(1/3)
      len.cubes <- d.length 
      
      x.cut <- seq(floor(range(sub.vert$x)[1]),range(sub.vert$x)[2]+len.cubes,len.cubes)
      # Numeber of intervals on x axis
      n.x <- length(x.cut) - 1
      
      y.cut <- seq(floor(range(sub.vert$y)[1]),range(sub.vert$y)[2]+len.cubes,len.cubes)
      # Numeber of intervals on y axis
      n.y <- length(y.cut) -1
      
      z.cut <- seq(floor(range(sub.vert$z)[1]),range(sub.vert$z)[2]+len.cubes,len.cubes)
      # Numeber of intervals on z axis
      n.z <- length(z.cut) - 1
      
      # Actual total number of voxels:
      n.v <- n.x*n.y*n.z
      
      # For each vertex, finding which intervals its x, y, z in. 
      vl <- cbind(findInterval(sub.vert$x,x.cut),findInterval(sub.vert$y,y.cut),findInterval(sub.vert$z,z.cut))
      
      # Mapping the indices into the labelling of each cube. 
      vert.Z <- rep(NA, P)
      for(i in 1:P)
      {
        vert.Z[i] <- vl[i,1] + (vl[i,2] -1)*n.x + (vl[i,3] -1)*(n.x*n.y)
      }
      
      
      
      ####################################################################
      ####################################################################
      # Scale the measurement of MEG and EEG.
      # The following two lines are used for REAL but not simulated data
      # For simulated data the model generates Y_M and Y_E that are already scaled
      
      if(sim){
        Y_M 
        Y_E
      }else{
      Y_M <- Y_M / sqrt((1 / n_M)* sum(diag(Y_M %*% t(Y_M))))
      Y_E <- Y_E / sqrt((1 / n_E)* sum(diag(Y_E %*% t(Y_E))))
      }
      
      
    
      
      ################   LASSSO #####################################
    
      Y.start<-rbind(Y_E,Y_M)
      X.start<-rbind(X_E,X_M)
      mfit <- glmnet(X.start, Y.start, family = "mgaussian", alpha= 1)
      s_median <- median(mfit$lambda)
      non_zero_coef <- coef(mfit, s = s_median)
      
      S <- matrix(0, ncol = T, nrow = P)
      for ( t in 1:T){
        S[,t] <- as.vector(non_zero_coef[[t]])[-1]
      }
      
     
      
      index.active <- which(rowSums(abs(S))!=0)
      index.inactive <-  which(!(rowSums(abs(S))!=0))
      
      S_smooth <- S
      
      for ( i in index.active){
        time <- 1:T
        act_source <- S[i,] 
        S_smooth[i,] <- (loess(act_source ~ time,span = 0.2))$fitted
      }
      
      S <- S_smooth
      S_initial <- S
      
      S_original_scale <- diag(1/sd_X) %*% S
  
      # Map the S_original_scale here to origanl vertice size.
      S_trans  <- matrix(NA,nrow = O_P,ncol = T)
      for ( c in 1:O_P)
      {
        p_idx <-  as.vector(as.vector(sub_cluster$cluster))[c]
        S_trans[c, ] <- S_original_scale[p_idx,]
      } 
      
  
      
      Z.state <- rep(0,P)
      Z.state[index.inactive] <- 1
      
      if (K ==2){
        Z.state[index.active] <- 2
      }else{
        S_kms <- S
        max_S_active <- apply(abs(S[index.active,]), 1, max)
        S_kms[index.active, ] <-  S[index.active,] * (1/max_S_active)
        
        vert.cluster <- kmeans(S_kms[index.active,],K-1,nstart = 100)
        
        for (k in 1:(K-1)){
          c.index <- index.active[which(vert.cluster$cluster==k)]
          Z.state[c.index] <- k+1
        }
      }
    
      # State variable for each voxel 
      cube.state <- rep(0,n.v)
      
      # Empty voxels are assigned to be in state 1, which is the inactive state.
      v.index.nz <- as.numeric(names(table(vert.Z)))
      cube.state[-v.index.nz] <- 1
      
      for (v in v.index.nz)
      {
        idx <- which(vert.Z==v) 
        p.Z.state.v <- (as.vector( table(c(Z.state[idx],1:K))) - 1) / as.vector(table(Z.state))
        
        cube.state[v] <- which.is.max(p.Z.state.v)
      }
      
      # State variable for each cluster of vertices.
      Z.state <- cube.state[vert.Z]
  
      ################   LASSSO #####################################
      
      
      ####################################################################
      ####################################################################
      ## Hyper-parameters set-up
      a_E <- 0.1; b_E <- 0.1
      a_M <- 0.1; b_M <- 0.1
      a_a <- 0.1; b_a <- 0.1
      a_alpha <- 0.1; b_alpha <- 0.1
      sigma2_A <- 0.25
      sigma2_mu1 <- 1 # check sensitivity
      beta_u <- 0.5*2/3*log(0.5*(sqrt(2) + sqrt(4*K - 2))) 
  
      ###################################################################
      ## Initializing values
      beta <- runif(1,0,beta_u)
      
      if (K ==2){
        A <- 0.95
      }else
      {A <- diag(sample(seq(0.9,0.95,0.01),K-1,replace = TRUE))}
  
      # 3D array specifying cubes.
      mask <- array(1, dim = c(n.x, n.y, n.z))
      
      # Define the neighborhood structure(First order) and get the neighbor matrix.
      neiStruc <- matrix(c(2,2,0,0,
                           0,2,0,0,
                           0,0,0,0), nrow=3, byrow=TRUE)
      # Neighborhood matrix
      neighbors <- getNeighbors(mask, neiStruc)
      
      # Blocks for Chequeboard Updating
      blocks <- getBlocks(mask, nblock=2)
  
      # Get the index for all the white blocks
      white_blocks <- blocks[[1]]
      
      # Get the index for all the black blocks
      black_blocks <- blocks[[2]]
      
      
      # For the variance components 
      sigma2_a <- b_a/(a_a+1) #initialize at prior mode
      
      alpha <- rep(b_alpha/(a_alpha+1),K)
  
      # Initial values for mu
      mu <- matrix(0, nrow = K, ncol = T)
      
      for (k in 2:K)
      {
        if( sum(Z.state==k) ==1)
        {mu[k,] <- S[Z.state==k,] }
        else
        {mu[k,] <- colMeans(S[Z.state==k,]) }
        
      }
  
      
      sigma2_E <- mean(apply(Y_E - X_E %*% S, 1, var))
      sigma2_M <- mean(apply(Y_M - X_M %*% S, 1, var))
      
      
      
      
      ## ICM updating algorithm and intial set-up.
      R  <- 100 #NUMBER OF ICM ITERATIONS TO RUN
      
      ##Stored Values
      sigma2_M_star <- rep(0,R)
      sigma2_M_star[1]  <- sigma2_M
      
      sigma2_E_star <- rep(0,R)
      sigma2_E_star[1]  <- sigma2_E
      
      sigma2_a_star <- rep(0,R)
      sigma2_a_star <- sigma2_a
      
      vec_A_star <- matrix(0,nrow = (K-1)^2, ncol = R)
      vec_A_star[,1] <- as.vector(A)
      
      alpha_star <- matrix(0, nrow = K, ncol = R)
      alpha_star[,1] <- alpha
      
      mu_star <- array(0,dim = c(K,T,R))
      mu_star[, , 1] <- mu
      
      Z.nv <- matrix(0,nrow = n.v, ncol = R)
      Z.nv[,1] <- cube.state
      
      beta_star <- rep(0,R)
      beta_star[1] <- beta
      
      inv_H_M <- solve(H_M)
      inv_H_E <- solve(H_E)
      
      W_1j_C1 <- rep(0,P)
      W_1j_C2 <- rep(0,P)
      
      HX.M <- matrix(0,nrow = n_M,ncol = P)
      HX.E <- matrix(0,nrow = n_E,ncol = P)
      
      W2j_C1 <- matrix(0, nrow = T, ncol = P )
      W2j_C2 <- matrix(0, nrow = T, ncol = P )
      
      
      for ( j in 1:P)
      {
        HX.M[,j] <- inv_H_M %*% X_M[,j]
        HX.E[,j] <- inv_H_E %*% X_E[,j]
        
        W_1j_C1[j] <- crossprod(X_M[,j],HX.M[,j])
        W_1j_C2[j] <- crossprod(X_E[,j],HX.E[,j])
        
        for(t in 1:T){
          W2j_C1[t,j] <- -2 * crossprod(Y_M[,t], HX.M[,j])
          W2j_C2[t,j] <- -2 * crossprod(Y_E[,t], HX.E[,j])
        }
      }
      
      # Define function of -H(beta) and -H'(beta)
      
      H_beta <- function(beta,cube.state,n.v,neighbors)
      {
        test <- matrix(cube.state[neighbors],ncol = 6,byrow = FALSE)
        term1 <- 2*beta*sum(rowSums(test == cube.state,na.rm = TRUE))
        nn.count<-apply(test,1,tabulate,nbins=K)
        term2 <-sum(log(colSums(exp(2*beta*nn.count))))
        # Return the negative of H_beta 
        result <- -1*(term1 - term2)
        result
      }
      
      HP_beta <- function(beta,cube.state,n.v,neighbors)
      {
        test <- matrix(cube.state[neighbors],ncol = 6,byrow = FALSE)
        term1 <- 2*sum(rowSums(test == cube.state,na.rm = TRUE))
        nn.count<-apply(test,1,tabulate,nbins=K)
        term2.1 <- 1 / colSums(exp(2*beta*nn.count))
        term2.2 <- 2*colSums(exp(2*beta*nn.count) * nn.count)
        term2 <- sum(term2.1 * term2.2)
        result <- -1*(term1 - term2)
        result
      }
      
      f_norm <- rep(NA, R)
      
      f_norm[1] <- sqrt(sum(S*S))
      
      count <- 0
      failure <- 0
      div <- 0
      
      r <- 1
      
      while (r < R) {
        time.iter<-proc.time()
        # Update the sigma2_M
        a_M_star <- a_M + T*n_M / 2
        temp <- Y_M - X_M %*% S
        b_M_star <- 1/2 * sum( diag(crossprod(temp,inv_H_M %*%temp ))) + b_M
        sigma2_M_star[r+1] <- b_M_star / (a_M_star + 1)
        sigma2_M <- b_M_star / (a_M_star + 1)
        #
        # Update the sigma2_E
        a_E_star <- a_E + T*n_E /2
        temp <- Y_E - X_E %*% S
        b_E_star <- 1/2 * sum( diag(crossprod(temp,inv_H_E%*%temp))) + b_E
        sigma2_E_star[r+1] <- b_E_star / (a_E_star + 1)
        sigma2_E <- b_E_star / (a_E_star + 1)
        #
        # Update the  sigma2_a
        a_a_star <- a_a + (T -1) * (K - 1) / 2
        temp <- mu[2:K, 2:T]- A%*%mu[2:K,1:T-1]
        b_a_star <- 1/2 * sum( diag( crossprod(temp))) + b_a
        sigma2_a_star[r+1] <- b_a_star / (a_a_star + 1)
        sigma2_a <- b_a_star / (a_a_star + 1)
        
        # 
        # Update the vec(A)
        sKr_t <- 0
        vc <- 0
        for (t in 2:T)
        {
          Kr_t <- kronecker( t(mu[2:K,t - 1]), diag(1, K-1))
          sKr_t <-  sKr_t + crossprod(Kr_t)#t(Kr_t) %*% Kr_t
          vc <- vc  + crossprod(mu[2:K,t],Kr_t)#t(mu[2:K,t]) %*% Kr_t
        }
        
        C_1 <- (1 / sigma2_A) * diag(1,(K-1)^2) + (1/sigma2_a) * sKr_t
        V_1 <- t( (1/sigma2_a) * vc %*% solve(C_1))
        A <- matrix(V_1, nrow = K-1)
        vec_A_star[,r+1] <- as.vector(A)
        
        # Update for each alpha_l, for each alpha_l, it's a inverse-gamma distribution.
        #  t.temp<-proc.time()
        a_alpha_l_star <- rep(0,K)
        b_alpha_l_star <- rep(0,K)
        for ( l in 1:K)
        {
          a_alpha_l_star[l] <- a_alpha + 1/2 * T * sum(Z.state == l)
          if (length(which(Z.state == l))==1)
          {
            b_alpha_l_star[l] <- b_alpha + 1/2 * sum ( (sweep(t(S[which(Z.state == l),]), 2, mu[l,]))^2 )
          }
          else
          {
            b_alpha_l_star[l] <- b_alpha + 1/2 * sum ( (sweep(S[which(Z.state == l),], 2, mu[l,]))^2 )
          }
          alpha[l] <- b_alpha_l_star[l] / ( a_alpha_l_star[l] + 1)
        }
        alpha_star[,r+1] <- alpha
        #  t.temp<-proc.time()-t.temp
        # 
        # Update mu_l(t=1) for all l =1, ... K.
        D_j <- array(0, dim = c(K-1, K-1,P))
        STD_j <- matrix(0,P,K-1)
        for (j in 1:P)
        {
          if (Z.state[j] != 1)
          {
            if(K ==2){
              D_j[, , j] <- 1 / alpha[Z.state[j]]
            }else{
              D_j[, , j][Z.state[j]-1,Z.state[j]-1] <- 1 / alpha[Z.state[j]]}
            
          }
          STD_j[j,] <- t(rep(S[j,1],K-1)) %*% D_j[, , j]
        }
        
        SD_j <- apply(D_j, 1:2,sum)
        B_1 <- SD_j + 1/sigma2_a * t(A)%*%A + 1 / sigma2_mu1 * diag(1,nrow = K-1, ncol = K-1)
        M_1 <- t( ( t(apply(STD_j,2,sum)) + 1/sigma2_a*t(mu[2:K,2]) %*% A) %*% solve(B_1))
        mu[2:K,1] <- M_1
        
        # Update mu_l(t) for all l =2 , ... K when 1 < t < T,
        B_2 <- SD_j + 1 / sigma2_a * (t(A)%*%A + diag(1,K-1,K-1))
        inv_B_2 <- solve(B_2)
        STD_jt <- matrix(0,P,K-1)
        time_interval <- seq(2,T-1)
        for ( t in time_interval)
        {
          for (j in 1:P)
          {
            STD_jt[j,] <- t(rep(S[j,t],K-1)) %*% D_j[, , j]
          }
          SSTD <- t(apply(STD_jt,2,sum))
          M_2_1 <- 1/sigma2_a*t(mu[2:K,t+1])%*%A
          M_2_2 <- 1/sigma2_a*t(mu[2:K,t-1]) %*% t(A)
          M_2 <- t((SSTD + M_2_1 + M_2_2) %*% inv_B_2)
          mu[,t] <- rbind(0,M_2)
        }
        
        #Update mu_l(T) for all l=2, ..., K, when t = T.
        B_3 <- SD_j +1 / sigma2_a*diag(1,K-1,K-1)
        inv_B_3 <- solve(B_3)
        STD_jT <- matrix(0,P,K-1)
        for (j in 1:P)
        {
          STD_jT[j,] <- t(rep(S[j,T],K-1)) %*% D_j[, , j]
        }
        # SD_j <- apply(D_j, 1:2,sum)
        M_3 <- t( ( t(apply(STD_jT,2,sum)) + 1/sigma2_a*t(mu[2:K,T-1]) %*% t(A)) %*% inv_B_3)
        mu[2:K,T] <- M_3
        
        mu_star[,,r+1] <- mu
        
        #UPDATE FOR BETA NEEDS TO BE ADDED
        beta.opt <-nlminb(start = beta, objective = H_beta,gradient = HP_beta, cube.state = cube.state,n.v= n.v, neighbors = neighbors,lower = 0,upper = beta_u)
        if(beta.opt$convergence != 0){
          cat("Warning, nlminb didn't converge for beta, iteration = ", r,"\n")
        }
        beta <- beta.opt$par
        beta_star[r+1] <- beta.opt$par
        
        
        # Replace the full conditional updating with cheuqboard updating. 
        # Updating all the white blocks first 
        
        P.Z <- matrix(0, nrow = n.v,ncol = K)
        n_r <- rep(0,n.v)
        n_r.index <- as.numeric(names(table(vert.Z)))
        n_r.values <- as.vector(table(vert.Z))
        n_r[n_r.index] <- n_r.values
        #for (v in v.index.nz)
        for (v in white_blocks)
        {
          for (h in 1:K)
          {
            log.term1 <- (-T*n_r[v]/2)*log(alpha[h])
            
            if (n_r[v] == 0)
            {
              log.term2 <- 0
            }
            else
            {
              if(n_r[v] ==1)
              {
                S.term2 <- t(S[which(vert.Z == v),])
              }
              else
              {
                S.term2 <- S[which(vert.Z == v),]
              }
              log.term2 <- -sum ((sweep(S.term2, 2, mu[h,]))^2)/ (2*alpha[h])
            }
            term3.neighbors <-neighbors[v,neighbors[v,] != (n.v+1)]
            log.term3 <- 2*beta*sum(cube.state[term3.neighbors] ==h)
            P.Z[v,h] <-  log.term1 + log.term2+ log.term3
            
          }
          P.Z[v,] <-  P.Z[v,] -  max(P.Z[v,]) #ADDED FOR NUMRECIAL STABILITY
          P.Z[v,] <- exp(P.Z[v,])
          P.Z[v,] <-  P.Z[v,] / sum( P.Z[v,])
          cube.state[v] <- which.is.max(P.Z[v,]) #UPDATE THE STATE FOR VOXEL
          
        }
        
        # Now, we will update all the black blocks
        for (v in black_blocks)
        {
          for (h in 1:K)
          {
            log.term1 <- (-T*n_r[v]/2)*log(alpha[h])
            
            if (n_r[v] == 0)
            {
              log.term2 <- 0
            }
            else
            {
              if(n_r[v] ==1)
              {
                S.term2 <- t(S[which(vert.Z == v),])
              }
              else
              {
                S.term2 <- S[which(vert.Z == v),]
              }
              log.term2 <- -sum ((sweep(S.term2, 2, mu[h,]))^2)/ (2*alpha[h])
            }
            term3.neighbors <-neighbors[v,neighbors[v,] != (n.v+1)]
            log.term3 <- 2*beta*sum(cube.state[term3.neighbors] == h)
            P.Z[v,h] <-  log.term1 + log.term2+ log.term3
            
          }
          P.Z[v,] <-  P.Z[v,] -  max(P.Z[v,]) #ADDED FOR NUMRECIAL STABILITY
          P.Z[v,] <- exp(P.Z[v,])
          P.Z[v,] <-  P.Z[v,] / sum( P.Z[v,])
          cube.state[v] <- which.is.max(P.Z[v,]) #UPDATE THE STATE FOR VOXEL
        }
        
        Z.state <- cube.state[vert.Z]
        Z.nv[,r+1] <- cube.state #STORE UPDATED STATES FROM THIS SWEEP
        
        
        # Updating for S
        W_1j <- 1/sigma2_M * W_1j_C1 + 1/sigma2_E*W_1j_C2+ 1/alpha[Z.state]
        for (j in 1:P)
        {
          # j = 1
          # W_2j <- 1/sigma2_M*( W2j_C1[,j] + 2*crossprod(X_M[,-j] %*% S[-j,], HX.M[,j]) ) + 1/sigma2_E*( W2j_C2[,j] + 2*crossprod(X_E[,-j] %*% S[-j,], HX.E[,j])) -2*rowSums( crossprod(mu, diag(1/alpha,nrow = K,ncol = K)))
          W_2j <- 1/sigma2_M*( W2j_C1[,j] + 2*crossprod(S[-j,], t(X_M[,-j])%*% HX.M[,j])) + 1/sigma2_E*( W2j_C2[,j] + 2*crossprod(S[-j,],t(X_E[,-j])%*% HX.E[,j])) -2*rowSums( crossprod(mu, diag(1/alpha,nrow = K,ncol = K)))
          Sigma_S_j <- (1/as.numeric(W_1j[j]))*diag(1,nrow = T,ncol = T)
          mu_sj <- -1/2* Sigma_S_j %*% W_2j
          S[j,] <- mu_sj
        }
        
        
        S_original_scale <- diag(1/sd_X) %*% S
  
        # Map the S_original_scale here to origanl vertice size.
        S_trans <- matrix(NA,nrow = O_P,ncol = T)
        new_Z.state  <- rep(NA, O_P)
        for ( c in 1:O_P)
        {
          p_idx <-  as.vector(as.vector(sub_cluster$cluster))[c]
          S_trans[c, ] <- S_original_scale[p_idx,]
          new_Z.state[c] <- Z.state[p_idx]
        } 
        
        
        f_norm[r+1] <- sqrt(sum(S_trans*S_trans))
        
        f_ratio <- (abs(f_norm[r+1] - f_norm[r]) / f_norm[r]) *100
        
        if (f_ratio < 0.001){
          count <- count +1
        }else if(f_ratio > 20) {
          failure <- failure +1
        } else{
          count <- 0
          failure <- 0
        }
        
        
  #       time.iter<-proc.time()-time.iter
  #       cor.iter <- cor(c(true.S),c(S_trans))
  #       rmse <- sqrt(sum( ( true.S - S_trans)^2) /(O_P*T))
  
        r <- r+1
        
        if (count == 10)
        {break}
        
        if (failure == 10)
        {
          cat("Diverge! Re-run the algorithm!", "\n")
          div <- TRUE
          break
        }
        
      } 
      
      
      # Kernel Smoothing for active state.
      index.active <- which(Z.state!=1)
      
      S_smooth <- S
      
      for ( i in index.active){
        time <- 1:T
        act_source <- S[i,] 
        S_smooth[i,] <- (loess(act_source ~ time,span = 0.2))$fitted
      }
      
      S <- S_smooth
      
      S_original_scale <- diag(1/sd_X) %*% S
  
      # Map the S_original_scale here to origanl vertice size.
      S_trans  <- matrix(NA,nrow = O_P,ncol = T)
      for ( c in 1:O_P)
      {
        #c_idx <- which( as.vector(sub_cluster$cluster) %in% c)
        p_idx <-  as.vector(as.vector(sub_cluster$cluster))[c]
        S_trans[c, ] <- S_original_scale[p_idx,]
      } 
      
      #close 
    }
    
    
    color_vec <- c("black", "blue", "green","red", "yellow", "pink", "purple","orange")
    K_hat_icm <- length(names(table(Z.state)))
    active_labels <- as.numeric(names(table(Z.state))[-1])
    signal_max  <- matrix(NA,nrow = length(active_labels), ncol = T)
    #cat("Algorithm Complete. The estimated number of mixture component is ",  K_hat_icm, sep = "")
    
    est_Z_state <- new_Z.state
    est_Z <- as.numeric(names(table(new_Z.state)))
    for (k in 1: length(est_Z)){
          new_Z.state[which(new_Z.state == est_Z[k])] <- color_vec[k]
    }
    
    # Estimated mixture components allocation 
    par3d(windowRect = c(20, 30, 600, 600))
    plot3d(op_sub.vert, col = new_Z.state, pch=19)
    bgplot3d({ 
      plot.new() 
      title(main = 'ICM Estimated Partition: Z_hat', line = 1, cex.main=1.5)})
    
    x11()
    matplot(t(S_trans), col = new_Z.state, ylab = "Estimated Source",  xlab = "Time Points",type = 'l')
    
   
    if(sim)
    {
    open3d()
    tr_Z <- as.numeric(names(table(true_Z)))
    for (k in 1: length(tr_Z)){
      true_Z[which(true_Z == tr_Z[k])] <- color_vec[k]
    }
     # True mixture component allocation
    par3d(windowRect = c(20, 30, 600, 600))
    plot3d(op_sub.vert, col = true_Z, pch=19)
    bgplot3d({ 
    plot.new() 
    title(main = 'True Partition:Z_true', line = 1, cex.main=1.5)})
    }
    
    # Compute the total power and get the clusters have the highest power.   
    vt_max <- list()
    Power_sub <-  rowSums(S_original_scale^2)
    N <- 15
    ndx <- order(Power_sub, decreasing = T)[1:N]
    top_power <- Power_sub[ndx]
    
    for(k in 1:N){
      vt_max[[k]] <- which(as.vector(sub_cluster$cluster) ==ndx[k])
    }
    
    
    open3d()
    par3d(windowRect = c(20, 30, 600, 600))
    act_v  <-  NULL
    for ( i in 1:length(vt_max)){
      act_v <- c(act_v,vt_max[[i]])
    }
    
    plot3d(op_sub.vert[-c(act_v),],pch=19)
    
    
    for (k in 1:N)
    {
      points3d(op_sub.vert[vt_max[[k]],], pch=19,col="red")
    }
    
    bgplot3d({ 
    plot.new()
    title(main = 'Highest Source Power', line = 1, cex.main=1.5)})
    
     est_S <- S_trans
     
     # Save the estimated resutls to results directory.
     save(est_S,est_Z_state,beta,mu,A, sigma2_a, sigma2_E, sigma2_M, alpha, file = "./results/estimated_result.RData")

     cat("\n\n\n")
     cat("Algorithm Complete. The estimated number of mixture components is K_ICM = ",  K_hat_icm, sep = "")
