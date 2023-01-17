# loading the packages
library(extRemes)
library(tidyverse)
library(geodist)
library(foreach)
library(mvnfast)
library(mvtnorm)

#### coordinates and geodesic distance ####
geocoord <- read.csv("coordinates.csv", header = TRUE)

ind_sites <- which(between(geocoord$geolat, 35, 67) & 
                     between(geocoord$geolon, -25, 8))
geocoordinates <- data.frame(geocoord$geolon[ind_sites],                                                                                                                  
                             geocoord$geolat[ind_sites])
names(geocoordinates) <- c("lon", "lat")

dist  <- geodist(geocoordinates, measure = "geodesic")
dist <- dist/100000
ind_dist <- which(upper.tri(dist, diag = FALSE), arr.ind = TRUE)
geodist <-  dist[ind_dist[order(ind_dist[,1]),]]

#### data ####
data <- read.csv("data.csv", header = TRUE)

# obtaining only the nine sites of interest
df_sites <- data[, ind_sites]

# removing rows where there is at least one missing value
df_final <- df_sites[complete.cases(df_sites),]

# preparing the pairwise datasets
pairlist <- as.list(df_final)
pairwise <- as.list(combn(pairlist, 2))
ind_col_1 <- seq(1, length(pairwise) - 1, by = 2)
ind_col_2 <- ind_col_1 + 1
col_1 <- purrr::map(ind_col_1, ~pairwise[[.x]])
col_2 <- purrr::map(ind_col_2, ~pairwise[[.x]])
length(col_1) == length(col_2)

# function to obtain the dataframes in a pairwise fashion
make_df <- function(col1, col2){
  new_data <- data.frame(col1, col2)
  return(new_data)
}

# complete dataframes in a list
complete_dfs <- purrr::map2(.x = col_1, .y = col_2, make_df) 


#### transforming to uniform margins ####
# function to transform the margins of each dataframe to uniform
unif_margins <- function(data){
  data1 <- data[,1]
  data2 <- data[,2]
  
  vec_unif1 <- numeric(length = length(data1))
  vec_unif2 <- numeric(length = length(data2))
  
  q <- 0.95
  th1 <- quantile(data1, q)
  th2 <- quantile(data2, q)
  
  overth1 <- data1 >= th1
  underth1 <- data1 < th1
  overth2 <- data2 >= th2
  underth2 <- data2 < th2
  
  fitgpd1 <- fevd(data1,  method = "MLE", type = "GP", threshold = th1)
  par1 <- fitgpd1$results$par
  fitgpd2 <- fevd(data2,  method = "MLE", type = "GP", threshold = th2)
  par2 <- fitgpd2$results$par
  
  ranks1 <- rank(data1)/(length(data1)+1)
  ranks2 <- rank(data2)/(length(data2)+1)
  
  pgpd_new <- function(data, scale, shape, lambda, threshold){
    p <- pmax(1 + (shape*(data - threshold))/scale, 0)
    p <- 1 - lambda*p^(-1/shape)
    return(p)
  }
  
  unifgpd1 <-  pgpd_new(data1, scale = par1[1], shape = par1[2], 
                        lambda = 1 - q, threshold = th1)
  unifgpd2 <-  pgpd_new(data2, scale = par2[1], shape = par2[2], 
                        lambda = 1 - q, threshold = th2)
  vec_unif1[underth1] <- ranks1[underth1]
  vec_unif1[overth1] <- unifgpd1[overth1]
  vec_unif2[underth2] <- ranks2[underth2]
  vec_unif2[overth2] <- unifgpd2[overth2]
  unif_mar <- data.frame(vec_unif1, vec_unif2)
  return(unif_mar)
}

# applying the function above to the list of dataframes
unif_df <- purrr::map(complete_dfs, unif_margins)

# censored Gaussian negative log-likelihood
nll_Gaussian <- function(rho, datU, thresh){
  
  z <- matrix(sapply(datU, qnorm), ncol = ncol(datU), nrow = nrow(datU))   
  uz <- qnorm(thresh)
  if (length(uz) == 1) {
    uz <- rep(uz, dim(z)[2])
  }
  else if (length(uz) < dim(z)[2]) {
    stop("Invalid censoring threshold")
  }
  
  if (rho < -0.9999 || rho > 0.9999) {
    return(1e+11)
  }
  
  Sig <- matrix(c(1, rho, rho, 1), ncol = 2)
  
  
  if (!exists(".Random.seed", mode = "numeric", envir = globalenv())) 
    sample(NA)
  oldSeed <- get(".Random.seed", mode = "numeric", envir = globalenv())
  
  ind <-  which(apply(z, 1, function(x) {sum(!is.na(x)) == dim(z)[2]}))
  
  z <- z[ind,]
  
  tmp <- apply(z, 1, function(t) {(sum(t > uz))})
  
  ind_part_cens <- c(1:dim(z)[1])[tmp > 0 & tmp < dim(z)[2]]
  ind_full_cens <- c(1:dim(z)[1])[tmp == 0]
  ind_no_cens <- c(1:dim(z)[1])[tmp == dim(z)[2]]
  
  if(length(ind_no_cens) > 0){
    ll1 <- -sum(mvnfast::dmvn(z[ind_no_cens, ], 
                              mu = rep(0, ncol(z)), 
                              sigma = Sig, 
                              log = TRUE)) +
      sum(dnorm(z[ind_no_cens,], log = TRUE)) 
  } else{ll1 <- 0}
  
  
  ll2 <- foreach::foreach(j = ind_part_cens, .combine = 'c') %dopar% {
    cens <- which(z[j, ] <= uz)
    nocens <- which(z[j, ] > uz)
    Sig11 <- Sig[cens, cens] - Sig[cens, nocens] %*% 
      (solve(Sig[nocens, nocens]) %*% Sig[nocens, cens])
    Sig11 <- as.matrix(Sig11)
    if(!isSymmetric.matrix(Sig11)){
      Sig11 <- (Sig11 + t(Sig11))/2
    }
    mu11 <- c(Sig[cens, nocens] %*% 
                (solve(Sig[nocens, nocens]) %*% 
                   z[j, nocens]))
    set.seed(123)
    mvnfast::dmvn(z[j, nocens], 
                  mu = rep(0, length(nocens)),
                  sigma = as.matrix(Sig[nocens, nocens]), 
                  log = TRUE) + 
      log(mvtnorm::pmvnorm(upper = uz[cens], 
                           mean = mu11, sigma = Sig11)[1]) - 
      sum(dnorm(z[j,nocens], log = TRUE))
  }
  ll2 <- -sum(ll2)
  
  set.seed(123)
  ll3 <- -length(ind_full_cens)*log(mvtnorm::pmvnorm(upper = uz, 
                                                     sigma = Sig)[1])
  
  assign(".Random.seed", oldSeed, envir = globalenv())
  return(ll1 + ll2 + ll3)
}

#### optmisation ####
mle <- matrix(rep(NA, 2), nrow = 1, ncol = 2)

comb <- function(...) { 
  mapply('rbind', ..., SIMPLIFY = FALSE)
}

# loop to obtain the MLEs and negative log-likelihood estimates 
# for the Gaussian copula 
loop <- foreach::foreach(i = 1:length(unif_df),
                         .combine = 'comb',
                         .multicombine = TRUE) %dopar% {
                           data <- unif_df[[i]]
                           data <- as.matrix(data)
                           fit <- optim(0.6,
                                        nll_Gaussian,
                                        method = "BFGS",
                                        thresh = 0.95,
                                        datU = data)
                           
                           mle[, 1] <- fit$par
                           mle[, 2] <- fit$value
                           list(mle)
                         }

mles <- data.frame(loop)
names(mles) <- c("rho", "nll_Gaussian")

# estimates for each model
#> print(mles)
#rho nll_Gaussian
#1  0.9863241  -473.198683
#2  0.9652405  -283.995622
#3  0.9857896  -435.782406
#4  0.8613851   120.186416
#5  0.7097369   312.610705
#6  0.9852277  -425.114979
#7  0.9866764  -430.010106
#8  0.9714839  -303.335174
#9  0.9648448  -276.455859
#10 0.9799662  -366.290839
#11 0.8702448   107.717808
#12 0.7319095   290.910204
#13 0.9795275  -368.535953
#14 0.9728336  -291.429854
#15 0.9687396  -273.059565
#16 0.9794567  -392.146757
#17 0.9131043     9.867028
#18 0.7461897   273.535119
#19 0.9646458  -275.517008
#20 0.9483830  -154.114057
#21 0.9821864  -395.823026
#22 0.8788641    86.939471
#23 0.7298333   292.799895
#24 0.9866249  -448.133902
#25 0.9772714  -319.506604
#26 0.9836648  -421.980756
#27 0.8159237   185.917185
#28 0.8591979   124.030247
#29 0.8366267   169.222046
#30 0.9059431     9.707709
#31 0.7095187   312.061879
#32 0.6863212   334.474167
#33 0.7474842   272.522087
#34 0.9811539  -345.494980
#35 0.9701906  -281.419806
#36 0.9578660  -191.293670