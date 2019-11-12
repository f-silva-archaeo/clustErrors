#' Error-Based Clustering for known number of clusters
#'
#' Generalization of k-means clustering algorithm for error-based clustering, developed by Kumar and Patel (2007).
#' @param data Array of measurements.
#' @param errors Array of errors. Must have same length as \emph{data}.
#' @param k Number of clusters.
#' @param verbose (Optional) Boolean to trigger printing to screen. Defaults to \emph{TRUE}.
#' @export
#' @seealso \code{\link[clustErrors]{hError}}
#' @examples
#' # Comparison with k-means when measurements all have same error
#' data <- c(rnorm(15, 3, 2), rnorm(15, 0, 5))
#' error <- rep(2, 30)
#' test0 <- kmeans(data, 2)
#' test1 <- kError(data, error, 2)
#'
#' # results are the same, since errors are the same
#' test0$centers
#' test1$metadata$cluster
#'
#' test0$cluster
#' test1$metadata$classification
#'
#' # Comparison when measurements have varying errors
#' error <- runif(30, 0, 5)
#' test2 <- kError(data, error, 2)
#'
#' # results are not the same, since errors are not the same
#' test0$centers
#' test2$metadata$cluster
#'
#' test0$cluster
#' test2$metadata$classification
kError = function(data, errors, k, verbose=T) {
  if (NROW(data) != NROW(errors)) { stop('data and errors must have the same length.') }
  n <- NROW(data)

  # 1. random partition of data into k clusters
  c <- sample(1:k, n, replace=T)

  iter <- 0
  repeat{
    # 2. compute cluster centers from equation 5
    theta <- c(); omega <- c()
    for (i in 1:k) {
      ind <- which(c==i);
      theta[i] <- Hmisc::wtd.mean(data[ind], 1/(errors[ind]^2))
      # omega[i] <- 1/sum(1/(errors[ind]^2), na.rm=T)
      omega[i] <- sum(1/errors[ind]^2, na.rm=T) / (sum(1/errors[ind], na.rm=T))^2 * stats::sd(data[ind])^2

    }

    # 3. reassign points to nearest cluster based on equation 33
    d <- matrix(NA, nrow = k, ncol = NROW(data))
    for (i in 1:k) {
      d[i,] <- (data - theta[i])*1/(errors^2)*(data - theta[i])
    }

    c.new <- apply(d,2,which.min)

    diff <- sum((c.new-c)^2)
    c <- c.new
    iter <- iter + 1

    # 4. repeat 2 and 3 above until convergence
    if(diff==0) { break}
  }


  # 5. output cluster means, and classification
  if (verbose) { cat(paste0("Number of iterations: ", iter, sep="")) }
  n <- c()
  for (i in 1:k) {
    n[i] <- length(which(c == i))
  }

  cl <- data.frame(mean=theta, sd=sqrt(omega))

  cluster <- c()
  cluster$metadata$method <- 'kError'
  cluster$metadata$k <- k
  cluster$metadata$n <- table(c)
  cluster$metadata$cluster <- cl
  cluster$metadata$classification <- c
  cluster$data <- data.frame(data=data, errors=errors)
  class(cluster) <- 'clustErrors'
  return(cluster)
}


#' Hierarchical Error-Based Clustering
#'
#' Generalization of Ward's hierarchical clustering algorithm for error-based clustering, developed by Kumar and Patel (2007).
#' @param data Array of measurements.
#' @param errors Array of errors. Must have same length as \emph{data}.
#' @param siglevel (Optional) Significance level for Chi-Squared test. Defaults to 0.05.
#' @param verbose (Optional) Boolean to trigger printing to screen. Defaults to \emph{TRUE}.
#' @export
#' @seealso \code{\link[clustErrors]{kError}}
#' @examples
#' data <- c(rnorm(50, 3, 2), rnorm(50, 0, 5))
#' error <- c(rep(2, 50), rep(5, 50))
#' test <- hError(data, error)
#'
#' test$metadata$cluster
#' test$metadata$classification
hError = function(data, errors, siglevel=0.05, verbose=T) {
  if (NROW(data) != NROW(errors)) { stop('data and errors must have the same length.') }

  # 1. initialization
  n <- NROW(data)
  k <- n
  clusters <- seq(1,k,1)
  classif <- clusters; dim(classif) <- c(1,n)
  theta <- data; psi <- errors^2

  dist <- matrix(NA, nrow=k, ncol=k)
  for (i in 1:k) {
    for (j in 1:k) {
      dist[i,j] <- (theta[i]-theta[j]) * 1/(psi[i] + psi[j]) * (theta[i] - theta[j])  # equation 5
    }
  }
  dist[dist==0] <- NA


  closest <- matrix(NA, nrow=2, ncol=k)
  for (i in 1:k) {
    aux <- which.min(dist[i,])
    if (length(aux)>0) {
      closest[1,i] <- aux
      closest[2,i] <- dist[i, aux]
    }
  }

  # 2. cycle
  iter <- 0
  if (verbose) { pb <- utils::txtProgressBar(max = n-1, style = 3) }

  repeat{
    # calculate Z2G statistic
    Z2G <- 0
    for (i in 1:k) {
      ind <- which(classif[iter+1,]==clusters[i])
      for (j in ind) {
        Z2G <-  Z2G + (data[j] - theta[which(clusters==clusters[i])]) * (1/errors[j]^2) * (data[j] - theta[which(clusters==clusters[i])])
      }
    }

    Chi <- stats::qchisq(siglevel, df = n-k, lower.tail = F)
    if (Z2G > Chi) { break }

    closestCl <- which.min(closest[2,]); closestCl[2] <- closest[1,closestCl]

    union.theta <- Hmisc::wtd.mean(theta[closestCl], 1/psi[closestCl])
    union.psi <- 1/sum(1/psi[closestCl], na.rm=T)
    theta[closestCl[1]] <- union.theta; theta <- theta[-closestCl[2]]
    psi[closestCl[1]] <- union.psi; psi <- psi[-closestCl[2]]
    classif <- rbind(classif, classif[iter+1,])
    classif[iter+2, which(classif[iter+2,]==clusters[closestCl[2]])] <- clusters[closestCl[1]]
    clusters <- clusters[-closestCl[2]]
    k <- k-1
    if (k==1) { break }

    # recalculate distances between classif
    dist <- matrix(NA, nrow=k, ncol=k)
    for (i in 1:k) {
      for (j in 1:k) {
        dist[i,j] <- (theta[i]-theta[j]) * 1/(psi[i] + psi[j]) * (theta[i] - theta[j])
      }
    }
    dist[dist==0] <- NA

    # update closest cluster
    closest <- matrix(NA, nrow=2, ncol=k)
    for (i in 1:k) {
      aux <- which.min(dist[i,])
      if (length(aux)>0) {
        closest[1,i] <- aux
        closest[2,i] <- dist[i, aux]
      }
    }

    iter <- iter + 1
    if (verbose) { utils::setTxtProgressBar(pb, iter) }
  }
  if (verbose) { utils::setTxtProgressBar(pb, n-1); cat('\nDone.') }


  c.relabel <- classif[dim(classif)[1],]
  un <- unique(c.relabel)
  for (i in 1:length(un)){
    c.relabel[which(c.relabel==un[i])] <- i
  }

  cl <- data.frame(mean=theta, error=sqrt(psi))

  cluster <- c()
  cluster$metadata$method <- 'hError'
  cluster$metadata$k <- k
  cluster$metadata$n <- table(c.relabel)
  cluster$metadata$cluster <- cl
  cluster$metadata$classification <- c.relabel
  cluster$metadata$history <- rbind(classif, as.numeric(c.relabel))
  cluster$data <- data.frame(data=data, errors=errors)
  class(cluster) <- 'clustErrors'
  return(cluster)
}

