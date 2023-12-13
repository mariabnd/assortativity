#' Calculate Q statistic
#'
#' @description Calculates the Q from by \insertCite{Aral1999-10.2105/AJPH.89.6.825;textual}{assortativity}
#' @param mat An n x n transmission matrix
#' @returns The value of Q (takes a value between - 1 / (n - 1) and 1)
#' @examples
#' mat <- matrix(data = c(0.25, 0.1, 0.15,
#' 0.05, 0.1, 0.07,
#' 0.08, 0.09, 0.11),
#' ncol = 3, nrow = 3, byrow = TRUE)
#' calculate_Q_stat(mat)
#' # -0.27
#' @export
#' @references
#' \insertAllCited{}
calculate_Q_stat <- function(mat){
  if(sum(mat) > 1){
    warning("mat may be given by counts rather than proportions")
  }
  nom <- sum(diag(mat)) - 1
  denom <- nrow(mat) - 1
  return(nom / denom)
}

#' Calculate r statistic
#'
#' @description Calculates the r from \insertCite{Newman2003-10.1103/PhysRevE.67.026126;textual}{assortativity}
#' @param mat An n x n transmission matrix
#' @returns The value of r (takes a value between 0 and 1)
#' @examples
#' mat <- matrix(data = c(0.25, 0.1, 0.15,
#' 0.05, 0.1, 0.07,
#' 0.08, 0.09, 0.11),
#' ncol = 3, nrow = 3, byrow = TRUE)
#' calculate_r_stat(mat)
#' # 0.1740593
#' @seealso [assortnet::assortment.discrete()]
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}
calculate_r_stat <- function(mat){
  if(sum(mat) > 1){
    warning("mat may be given by counts rather than proportions")
  }
  nom <- sum(diag(mat)) - sum(colSums(mat) * rowSums(mat))
  denom <- 1 - sum(colSums(mat) * rowSums(mat))
  return(nom / denom)
}

#' Calculate q statistic
#'
#' @description Calculates the q from \insertCite{KeelingRohani2008-10.2307/j.ctvcm4gk0.6;textual}{assortativity}
#' @param mat An n x n transmission matrix
#' @returns The value of q
#' @examples
#' mat <- matrix(data = c(0.25, 0.1, 0.15,
#' 0.05, 0.1, 0.07,
#' 0.08, 0.09, 0.11),
#' ncol = 3, nrow = 3, byrow = TRUE)
#' calculate_ev_ratio(mat)
#' # 3.994067
#' @export
#' @references
#' \insertAllCited{}
calculate_ev_ratio <- function(mat){
  if(sum(mat) > 1){
    warning("mat may be given by counts rather than proportions")
    }
  vals <- eigen(mat, only.values = TRUE)$values
  return(vals[1] / vals[2])
}

#' Calculate index of disassortativity
#'
#' @description Calculates the I-squared from \insertCite{Farrington2009-10.1002/bimj.200800160;textual}{assortativity}
#' @param pop_frac A vector of length n with the population distribution
#' @param mat An n x n transmission matrix
#' @param classes A vector of length n with the population sizes
#' @returns A list with three variants of the index
#' \describe{
#'   \item{I_stat_raw}{The index of disassortativity}
#'   \item{I_stat_standardised}{The standardised index of disassortativity (`I_stat_raw` standardised to the variance)}
#'   \item{I_stat_margin}{The index of assortativity standardised with respect to the marginal distribution}
#' }
#' @examples
#' mat <- matrix(data = c(0.25, 0.1, 0.15,
#' 0.05, 0.1, 0.07,
#' 0.08, 0.09, 0.11),
#' ncol = 3, nrow = 3, byrow = TRUE)
#' pop_frac <- c(0.35, 0.56, 0.21)
#' classes <- c(30, 27, 45)
#' calculate_I_stat(pop_frac, mat, classes)
#' # $I_stat_raw
#' # 40.76486
#' @export
#' @references
#' \insertAllCited{}
calculate_I_stat <- function(pop_frac, mat, classes = NULL){
  if(is.null(classes)){
    names(pop_frac) <- colnames(mat) <- rownames(mat) <- rep(1 : length(pop_frac))
    warning("It is assumed each population group consists of one (1) member\n
            To change this behaviour please include a classes argument")
  } else {
    names(pop_frac) <- colnames(mat) <- rownames(mat) <- classes
  }
  tmp_mat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
  for(i in 1 : nrow(mat)){
    for(j in 1 : ncol(mat)){
      tmp_mat[i, j] <- pop_frac[i] * mat[i, j] * pop_frac[j]
    }
  }
  f_st <- function(x, y){
    return(tmp_mat[x, y] / sum(tmp_mat))
  }
  f_s <- function(x){
    nom <- pop_frac[x] * sum(mat[x, ] * pop_frac)
    return(nom / sum(tmp_mat))
  }
  f_t <- function(x){
    nom <- sum(pop_frac * mat[, x]) * pop_frac[x]
    return(nom / sum(tmp_mat))
  }
  mn <- sum(as.numeric(names(pop_frac)) * pop_frac)
  vr <- sum((as.numeric(names(pop_frac)) - mn) ** 2 * pop_frac)
  marg_pop_frac <- sapply(1 : nrow(mat), f_s)
  marg_mn <- sum(as.numeric(names(pop_frac)) * marg_pop_frac)
  marg_vr <- sum((as.numeric(names(pop_frac)) - marg_mn) ** 2 * marg_pop_frac)
  # Calculate the index
  idx <- expand.grid(as.numeric(names(pop_frac)), as.numeric(names(pop_frac)))
  dist <- idx$Var1 - idx$Var2
  dist <- dist ** 2
  # Below required else indexing of matrix is (potentially) out of bounds
  idx <- expand.grid(1 : nrow(mat), 1 : ncol(mat))
  I_stat_parts <- dist * sapply(idx$Var1, f_s) * sapply(idx$Var2, f_t)
  I_stat_raw <- sum(I_stat_parts) / 2
  return(list(mean = mn,
              var = vr,
              f_s = sapply(1 : nrow(mat), f_s),
              marg_mn = marg_mn,
              marg_vr = marg_vr,
              f_t = sapply(1 : ncol(mat), f_t),
              I_stat_parts = I_stat_parts,
              I_stat_raw = I_stat_raw,
              I_stat_standardised = I_stat_raw / vr,
              I_stat_margin = I_stat_raw / marg_vr))
}

#' Plot the statistic on a number line
#'
#' @description Create a number line plot in the style of \insertCite{Garnett1996-10.1097/00007435-199605000-00015;textual}{assortativity}
#' @seealso [calculate_Q_stat()]
#' #' mat <- matrix(data = c(0.25, 0.1, 0.15,
#' 0.05, 0.1, 0.07,
#' 0.08, 0.09, 0.11),
#' ncol = 3, nrow = 3, byrow = TRUE)
#' plot_stat(mat)
#' @export
#' @references
#' \insertAllCited{}
plot_stat <- function(mat){
  n <- nrow(mat)
  Q <- calculate_Q_stat(mat)
  plot(NA,
       axes = FALSE,
       ann = FALSE,
       yaxs = "i",
       xlim = c(- 1 / (n - 1), 1),
       ylim = c(0, 2))
  axis(1, at = c(- 1 / (n - 1), 0, 1))
  points(x = Q, y = 0, pch = 16, col = 4)
  axis(1, line = 2, at = c(Q), col = 4)
}

#' Plot the matrix as an annotated heatmap
#'
#' @seealso [hhh4contacts::plotC()]
#' @examples
#' mat <- matrix(data = c(0.25, 0.1, 0.15,
#' 0.05, 0.1, 0.07,
#' 0.08, 0.09, 0.11),
#' ncol = 3, nrow = 3, byrow = TRUE)
#' plot_mat(mat)
#' @export
#' @references
#' \insertAllCited{}
plot_mat <- function(mat){
  n <- nrow(mat)
  grd <- expand.grid(seq(from = 0, to = 1, length.out = n),
                     seq(from = 0, to = 1, length.out = n))
  image(mat, axes = FALSE)
  text(x = grd$Var1, y = grd$Var2, labels = mat)
  }
