library(MASS)
library(ggplot2)
library(sf)
library(matrixStats)
library(boot)

## plotAreaCol
# This functions plots the values and saves the figure to the specified file.
# It is advisable to not plot directly in R since it will take a long time.
# The arguments are:
#   fNamme: file name for saving figure
#   width: width of figure in inches
#   height: height of figure in inches
#   estVal: the k values to be plotted on geoMap
#   geoMap: the map containing k regions
#   leg: name to use on top of legend
#   colLim: control lower limit and upper limit of color scale (colLim = c(lowVal, highVal))
plotAreaCol = function(fName, width, height, estVal, geoMap, leg, colLim = NULL){
  if(is.null(colLim)){
    colLim = range(estVal)
  }
  
  # Set up data object for plotting
  nigeriaMapTmp = geoMap
  nigeriaMapTmp$MCV1 = estVal
  
  # Plot
  map = ggplot() +
    geom_sf(data = nigeriaMapTmp,
            aes(fill = MCV1),
            color = 'gray', size = .2)+
    scale_fill_viridis_c(direction = 1,
                         begin = 1,
                         end = 0,
                         limit = colLim,
                         name = leg) + 
    theme(text = element_text(size=30),
          legend.key.height = unit(1, 'cm'),
          legend.key.width  = unit(0.4, 'cm'))
  ggsave(filename = fName,
         plot = map,
         width = width, 
         height = height)
}


setwd("C:/Users/marcu/OneDrive/Dokumenter/8. semester/Romlig statistikk/Project 3/")
admin1 <- read.table("admin1.txt", header = TRUE)
admin1 <- as.matrix(admin1)
admin2 <- read.table("admin2.txt", header = TRUE)
admin2 <- as.matrix(admin2)
#dimension 37x37 admin1, 775x775 admin2
#precision matrix -> diagonal elements equal to the number of elements, 
createPrecisionMatrix <- function(neighbour){
  n <- dim(neighbour)[1]
  precMatrix <- matrix(rep(0,n*n ), ncol = n)
  for (i in 1:n){
    precMatrix[i,i] <- sum(neighbour[i,])
    for (j in 1:n){
      if (neighbour[i,j] == 1 && i != j){
        precMatrix[i,j] <-  -1     
      }
    }
  }
  return(precMatrix)
}


admin1_precMatrix <- createPrecisionMatrix(admin1)
admin1_precMatrix
admin2_precMatrix <- createPrecisionMatrix(admin2)
admin2_precMatrix
#number of non-zero elements in the precision matrices
admin1_precMatrix_nonzero <-sum(colSums(admin1_precMatrix != 0))
admin2_precMatrix_nonzero <-sum(colSums(admin2_precMatrix != 0))
#proportion of non-zero elements in each of the matrices
admin1_prop_nonzero <- admin1_precMatrix_nonzero/(37*37)
admin2_prop_nonzero <- admin2_precMatrix_nonzero/(775*775)
admin1_prop_nonzero
admin2_prop_nonzero
admin1_sparsity <- admin1_precMatrix
admin2_sparsity <- admin2_precMatrix
for (i in 1:37){
  for (j in 1:37){
    if(admin1_sparsity[i,j] != 0){
      admin1_sparsity[i,j] <- 1
    }
  }
}
for (i in 1:775){
  for (j in 1:775){
    if(admin2_sparsity[i,j] != 0){
      admin2_sparsity[i,j] <- 1
    }
  }
}
admin1_sparsity <- admin1_sparsity[,c(37:1), drop = FALSE]
admin2_sparsity <- admin2_sparsity[,c(775:1), drop = FALSE]
par(pty = "s")
image(admin1_sparsity)
image(admin2_sparsity)

#sampling from besag model
epsilon <- 1e-8
Q_admin1 <- admin1_precMatrix
Q_admin2 <- admin2_precMatrix
besag_samples <- function(Q, epsilon){
  n <- dim(Q)[1]
  Q_tilde <- Q + epsilon*diag(n) #full rank matrix
  L_tilde <- t(chol(Q_tilde))
  z <- mvrnorm(1, mu = rep(0,n), diag(n))
  v <- solve(t(L_tilde),z)
  return (v-mean(v)*rep(1,n))
}
besag_1_1 <- besag_samples(Q_admin1, epsilon)
besag_1_2 <- besag_samples(Q_admin1, epsilon) 
mvr1_1 <- mvrnorm(n = 1, mu = rep(0, 37), Sigma = diag(37)) 
mvr1_2 <- mvrnorm(n = 1, mu = rep(0, 37), Sigma = diag(37)) 
besag_2_1 <- besag_samples(Q_admin2, epsilon)
besag_2_2 <- besag_samples(Q_admin2, epsilon)
mvr2_1 <- mvrnorm(n = 1, mu = rep(0, 775), Sigma = diag(775)) 
mvr2_2 <- mvrnorm(n = 1, mu = rep(0, 775), Sigma = diag(775)) 
load("Admin1Geography.RData")
load("Admin2Geography.RData")
#plots
plotAreaCol("besag_1_1_file.pdf", 12,7,besag_1_1,nigeriaAdm1, leg = "Besag simulation\n admin1", colLim = c(-4,4)) 
plotAreaCol("mvr1_1_file.pdf", 12,7,mvr1_1,nigeriaAdm1, leg = "MVN simulation\n admin1", colLim = c(-4,4)) 
plotAreaCol("besag_1_2_file.pdf", 12,7,besag_1_2, nigeriaAdm1,leg = "Besag simulation\n admin1", colLim = c(-4,4))
plotAreaCol("mvr1_2_file.pdf", 12,7,mvr1_2,nigeriaAdm1, leg = "MVN simulation\n admin1", colLim = c(-4,4)) 

plotAreaCol("besag_2_1_file.pdf", 12,7,besag_2_1,nigeriaAdm2, leg = "Besag simulation\n admin2", colLim = c(-4,4)) 
plotAreaCol("mvr2_1_file.pdf", 12,7,mvr2_1,nigeriaAdm2, leg = "MVN simulation\n admin2", colLim = c(-4,4)) 
plotAreaCol("besag_2_2_file.pdf", 12,7,besag_2_2, nigeriaAdm2,leg = "Besag simulation\n admin2", colLim = c(-4,4))
plotAreaCol("mvr2_2_file.pdf", 12,7,mvr2_2,nigeriaAdm2, leg = "MVN simulation\n admin2", colLim = c(-4,4)) 

#computing 100 realizations with the Besag model on the admin2 graph 
set.seed(42)
admin2_realizations <- matrix(rep(0, 775*100), ncol = 775, nrow = 100)
for (i in 1:100){
  realization <- besag_samples(Q_admin2, epsilon)
  admin2_realizations[i,] <- realization
}
variances <- colVars(admin2_realizations)
plotAreaCol("emp_variances_besag.pdf", 12,7, variances, nigeriaAdm2, leg = "Besag empirical \n variances")
#looking at Gubio index 150
gubio <- admin2_realizations[,150]
correlations <- cor(gubio, admin2_realizations)
plotAreaCol("emp_correlation_gubio_besag.pdf", 12,7, t(correlations), nigeriaAdm2, leg = "Besag empirical \n correlations")


#Problem 2
directEstimates <- read.table("DirectEstimates.txt", header = TRUE)
plotAreaCol("observed_proportions_admin1.pdf", 12,7, inv.logit(directEstimates$Observation), nigeriaAdm1, leg = "Observed \nproportions",
            c(0,1))

y <- directEstimates$Observation
V <- directEstimates$StdDev^2
n <- length(y)
variance <- 100^2
V_inv <- diag(n)*(1/V)
mu_b <- solve(variance^(-1)*V + diag(n))%*%y
Q_b <- variance^(-1)*diag(n)+V_inv
P_b <- matrix(rep(0, n*100), nrow = 100)
L_b <- t(chol(Q_b))
for (i in 1:100){
  z <- rnorm(n)
  x <- solve(t(L_b), z) + mu_b
  P_b[i,] <- inv.logit(x)
}
medians <- colMedians(P_b)
coefficientOfVariations <- sqrt(colVars(P_b))/colMeans(P_b)
plotAreaCol("medians_2b_admin1.pdf", 12,7, medians, nigeriaAdm1, leg = "Medians", colLim = c(0,1))
plotAreaCol("CV_2b_admin1.pdf", 12,7, coefficientOfVariations, nigeriaAdm1, leg = "CV", colLim = c(0,0.5))
#c
mu_c <- solve(Q_admin1+V_inv)%*%V_inv%*%y
Q_c <- Q_admin1+V_inv
L_c <- t(chol(Q_c))

P_c <- matrix(rep(0, n*100), nrow = 100)
for (i in 1:100){
  z <- rnorm(n)
  x <- solve(t(L_c), z)+mu_c
  P_c[i,] <- inv.logit(x)
}
medians_c <- colMedians(P_c)
coefficientOfVariations_c <- sqrt(colVars(P_c))/colMeans(P_c)
plotAreaCol("medians_2c_admin1.pdf", 12,7, medians_c, nigeriaAdm1, leg = "Medians", c(0,1))
plotAreaCol("CV_2c_admin1.pdf", 12,7, coefficientOfVariations_c, nigeriaAdm1, leg = "CV", c(0,0.5))

#Kaduna at index 19
extra_row <- rep(0,n)
extra_row[19] <- 1
M <- rbind(diag(n),extra_row)
V_tilde <- cbind(rbind(diag(n)*V, rep(0,n)),rep(0,n+1))
V_tilde[38,38] <- 0.1^2
V_tilde_inv <- solve(V_tilde)
y_tilde <- c(y,0.5)
mu_d <- solve(Q_admin1+t(M)%*%V_tilde_inv%*%M)%*%t(M)%*%V_tilde_inv%*%y_tilde
Q_d <- Q_admin1+t(M)%*%V_tilde_inv%*%M
L_d <- t(chol(Q_d))
P_d <- matrix(rep(0,n*100), nrow = 100)
for (i in 1:100){
  z <- rnorm(n)
  x <- solve(t(L_d), z)+mu_d
  P_d[i,] <- inv.logit(x)
}
medians_d <- colMedians(P_d)
coefficientOfVariations_d <- sqrt(colVars(P_d))/colMeans(P_d)
plotAreaCol("medians_2d_admin1.pdf", 12,7, medians_d, nigeriaAdm1, leg = "Medians", c(0,1))
plotAreaCol("CV_2d_admin1.pdf", 12,7, coefficientOfVariations_d, nigeriaAdm1, leg = "CV", c(0,0.5))

#e
tau_first <- 0.1
tau_second <- 10
Q_admin1_first <- tau_first*Q_admin1
Q_admin1_second <- tau_second*Q_admin1
#using tau = 0.1
mu_e_first <- solve(Q_admin1_first+V_inv)%*%V_inv%*%y
Q_e_first <- Q_admin1_first+V_inv
L_e_first <- t(chol(Q_e_first))

P_e_first <- matrix(rep(0, n*100), nrow = 100)
for (i in 1:100){
  z <- rnorm(n)
  x <- solve(t(L_e_first), z)+mu_e_first
  P_e_first[i,] <- inv.logit(x)
}
medians_e_first <- colMedians(P_e_first)
coefficientOfVariations_e_first <- sqrt(colVars(P_e_first))/colMeans(P_e_first)
plotAreaCol("medians_2e_first_admin1.pdf", 12,7, medians_e_first, nigeriaAdm1, leg = "Medians", c(0,1))
plotAreaCol("CV_2e_first_admin1.pdf", 12,7, coefficientOfVariations_e_first, nigeriaAdm1, leg = "CV", c(0,0.5))

#using tau = 10
mu_e_second <- solve(Q_admin1_second+V_inv)%*%V_inv%*%y
Q_e_second <- Q_admin1_second+V_inv
L_e_second <- t(chol(Q_e_second))

P_e_second <- matrix(rep(0, n*100), nrow = 100)
for (i in 1:100){
  z <- rnorm(n)
  x <- solve(t(L_e_second), z)+mu_e_second
  P_e_second[i,] <- inv.logit(x)
}
medians_e_second <- colMedians(P_e_second)
coefficientOfVariations_e_second <- sqrt(colVars(P_e_second))/colMeans(P_e_second)
plotAreaCol("medians_2e_second_admin1.pdf", 12,7, medians_e_second, nigeriaAdm1, leg = "Medians", c(0,1))
plotAreaCol("CV_2e_second_admin1.pdf", 12,7, coefficientOfVariations_e_second, nigeriaAdm1, leg = "CV", c(0,0.5))

#f
logit <- function(p){
  return(log(p/(1-p)))
}
logit(0.000000001)
logit(0.999999)
x <- rep(-20, 37)
mu_tau <- function(tau){
  return(solve(tau*Q_admin1 + V_inv)%*%V_inv%*%y)
}
Q_tau <- function(tau){
  return(tau*Q_admin1 + V_inv)
}

loglikmodel_c <- function(tau){
  return((-1)*(18*log(tau)-tau/2*t(x)%*%Q_admin1%*%x-
           1/2*log(det(Q_tau(tau)))+1/2*t(x-mu_tau(tau))%*%Q_tau(tau)%*%(x-mu_tau(tau))))
}
tau_hat <- optimize(loglikmodel_c, c(0.001,10))
tau_hat
tau_hat$minimum
#making 100 realizations with tau = tau_hat
L_f <- t(chol(Q_tau(tau_hat$minimum)))
P_f <- matrix(rep(0, n*100), nrow = 100)
mu_f <- mu_tau(tau_hat$minimum)
for (i in 1:100){
  z <- rnorm(n)
  x <- solve(t(L_f), z)+ mu_f
  P_f[i,] <- inv.logit(x)
}
P_f
medians_f <- colMedians(P_f)
coefficientOfVariations_f <- sqrt(colVars(P_f))/colMeans(P_f)
plotAreaCol("medians_2f_admin1.pdf", 12,7, medians_f, nigeriaAdm1, leg = "Medians" , c(0,1))
plotAreaCol("CV_2f_admin1.pdf", 12,7, coefficientOfVariations_f, nigeriaAdm1, leg = "CV", c(0,0.5))

