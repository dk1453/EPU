### function ###
GFEVD = function(Phi, Sigma, n.ahead=10,normalize=TRUE,standardize=TRUE) {
  tvp.Phi = function (x, nstep = 10, ...) {
    nstep = abs(as.integer(nstep))
    K=nrow(x)
    p=floor(ncol(x)/K)
    A = array(0, c(K,K,nstep))
    for (i in 1:p){
      A[,,i]=x[,((i-1)*K+1):(i*K)]#因为这行horizon不能为1
    }
    
    Phi = array(0, dim = c(K, K, nstep + 1))
    Phi[, , 1] = diag(K)
    Phi[, , 2] = Phi[, , 1] %*% A[, , 1]
    if (nstep > 1) {
      for (i in 3:(nstep + 1)) {
        tmp1 = Phi[, , 1] %*% A[, , i - 1]
        tmp2 = matrix(0, nrow = K, ncol = K)
        idx = (i - 2):1
        for (j in 1:(i - 2)) {
          tmp2 = tmp2 + Phi[, , j + 1] %*% A[, , idx[j]]
        }
        Phi[, , i] = tmp1 + tmp2
      }
    }
    return(Phi)
  }
  A = tvp.Phi(Phi, (n.ahead-1))
  Sigma = Sigma
  gi = array(0, dim(A))
  sigmas = sqrt(diag(Sigma))
  for (j in 1:dim(A)[3]) {
    gi[,,j] = t(A[,,j]%*%Sigma%*%MASS::ginv(diag(sqrt(diag(Sigma)))))
  }
  if (standardize==TRUE){
    girf=array(NA, c(dim(gi)[1],dim(gi)[2], (dim(gi)[3])))
    for (i in 1:dim(gi)[3]){
      girf[,,i]=((gi[,,i])%*%MASS::ginv(diag(diag(gi[,,1])))) #generate a diag matrix with only elements on diag
    }
    gi=girf
  }
  
  num = apply(gi^2,1:2,sum)
  den = c(apply(num,1,sum))
  fevd = t(num)/den
  nfevd = fevd
  if (normalize==TRUE) {
    fevd=(fevd/apply(fevd, 1, sum))
  } else {
    fevd=(fevd)
  }
  return = list(GFEVD=fevd, GIRF=gi)
}

cnct_analysis <- function(x){
  Self <- diag(x)
  From <- (rowSums(x) - Self) 
  To <- (colSums(x) - Self) 
  Net <- To - From
  Sum <- (sum(x) - sum(diag(x))) 
  return <- list(Self = Self, From = From, To = To, Net = Net, Sum = Sum)
}

cnct_pairwise_analysis <- function(x, i, j){
  pairwise <- (x[j, i] - x[i, j])
}


give_weight_rw <- function(x, y){
  cnct_a_rw_weighted <- array(NA, dim = c(dim(x)[1], dim(x)[2], dim(x)[3]))
  for (t in 1:dim(x)[3])
  {
    weighted_matrix_rw <- matrix(nrow = nrow(x), ncol = 0)
    for (j in 1:ncol(x))
    {
      weighted_vector_rw <- x[,j,t] * as.numeric(y[t,j])
      weighted_matrix_rw <- cbind(weighted_matrix_rw,weighted_vector_rw)
    }
    weighted_matrix_rw <- (weighted_matrix_rw / apply(weighted_matrix_rw, 1, sum)) *100
    cnct_a_rw_weighted[,,t] <- weighted_matrix_rw
  } 
  return(cnct_a_rw_weighted)
}

### NEW FUNCTION
cnct_analysis_updated <- function(x){
  Self <- diag(x)
  Sum <- (sum(x) - sum(diag(x)))
  From <- ((rowSums(x) - Self)/Sum) 
  To <- ((colSums(x) - Self) /Sum)
  Net <- To - From
  return <- list(Self = Self, From = From, To = To, Net = Net, Sum = Sum)
}

#目标是做出一个24*24的weight matrix，然后就可以直接和GVD matrix相乘相加
Calculate_To  <- function(x, y){
  matrixofweight <- matrix(NA, ncol = ncol(x), nrow = nrow(y))
  for (i in 1: nrow(y))
  {
    for (j in 1:nrow(y))
    {
      if(j==i){matrixofweight[j,i] <- 0}
      else{matrixofweight[j,i] <- y[j,]/(1-y[i,])}
    } 
  }
  matrix_for_To <- x * matrixofweight
  return <- matrix(colSums(matrix_for_To))*(ncol(x)-1)
}

give_weight_updated <- function(x,y,i){
  weighted_matrix <- matrix(nrow = 0, ncol = ncol(x))
  for (j in 1: nrow(x))
  {    
    weighted_vector <- x[j,] * as.numeric(y[j,i])
    weighted_matrix <- rbind(weighted_matrix,weighted_vector)
  }
  return(weighted_matrix)
}

give_weight <- function(x, y, i){
  #这里需要把矩阵里的元素as.numeric，不然r会使用矩阵运算，而不是数字乘矩阵。
  weighted_matrix <- matrix(nrow = 0, ncol = ncol(x))
  for (j in 1: nrow(x))
  {
    weighted_vector <- x[j,] * as.numeric(y[j,i])
    weighted_matrix <- rbind(weighted_matrix,weighted_vector)
  }
  weighted_matrix_st <- cnct_analysis_updated(weighted_matrix) # use cnct_analysis_updated from cnct_analysis from the ols.R
  weighted_matrix <- cbind(weighted_matrix, weighted_matrix_st$From)
  weighted_matrix <- rbind(weighted_matrix, append(weighted_matrix_st$To, 0), append(weighted_matrix_st$Net, 0))
  weighted_matrix[nrow(weighted_matrix),ncol(weighted_matrix)] <- weighted_matrix_st$Sum
  return(weighted_matrix)
}

weight_moving_average <- function(x, t){
  weight_MA <- matrix(NA, nrow(x) - (t-1), ncol(x))
  for (i in 1:(nrow(x)-t+1)){
    weight_MA[i,] <- colMeans(x[i:i+t-1,])
  }
  return(weight_MA)
}

give_weight_rw <- function(x, y){
  cnct_a_rw_weighted <- array(NA, dim = c(dim(x)[1], dim(x)[2], dim(x)[3]))
  for (t in 1:dim(x)[3])
  {
    weighted_matrix_rw <- matrix(nrow = 0, ncol = ncol(x))
    for (j in 1:nrow(x))
    {
      weighted_vector_rw <- x[j,,t] * as.numeric(y[t,j]) * dim(x)[1]
      weighted_matrix_rw <- rbind(weighted_matrix_rw,weighted_vector_rw)
    }
    cnct_a_rw_weighted[,,t] <- weighted_matrix_rw
  } 
  return(cnct_a_rw_weighted)
}

Rolling_Calculate_To <- function(x, y)
{
  array_for_To <- array(NA, dim = c(dim(y)[1], dim(y)[2], dim(y)[3]))
  for (t in 1:dim(x)[3])
  {
    array_for_To[,,t] <- Calculate_To(x[,,t],matrix(y[,,t]))
  }
  return(array_for_To)
}

