############################
## PreRequisite Functions ##
############################

## TODO these could really be cleaned up a lot

## * Tools
## ** Diagonal Scaling Matrix
   #' Diagonal Factors of Sparse Matrix
   #'
   #' Takes a sparse matrix
   #' and returns a diagonal matrix such that each column contains
   #' either 1 / colsum() or 0. (D^-1 in the papers)
   #' multiplication with this matrix would have all columns
   #' sum to 1 (or zero if that column already summed to zero)
   #'
   #' @param mat Input Matrix, may either be a matrix, sparse matrix
   #' dgCMatrix or probably even a data frame, as long as `nrow` and `colSum`
   #' work on the data as expected it should work, output is always a
   #' sparse matrix of class dgCMatrix
   #'
   #' This should take the an adjacency matrix as input and the output
   #' can be multiplied by the original matrix to scale it to 1.
   #'
   #' A directed adjacency matrix should be such that A[i,j]
   #' describes the the weight
   #' of an edge from vertex j to vertex i (this is, unfourtunately
   #' the transpose of what igraph gives back.)
   #'
   #' @return A diagonal matrix such that multiplication by the original
   #' matrix will lead to all columns adding
   #' to 1 (or 0 in the case of a column of zeros)
   #' The output class is sparse matrix "dgCMatrix" from the Matrix package
   #' use `as.matrix()` in order to get a matrix type, I didn't put that
   #' Inside the function because it would have been confusing.
   #'
   #'
   #' @examples
   #' (mat <- matrix(1:9, nrow = 3))
   #' mat %*% create_sparse_diag_scaling_mat(mat) %>% colSum()
   #' mat %*% create_sparse_diag_scaling_mat(mat)
   #'
   #' > [1] 1 1 1 1 1
   #' @export
create_sparse_diag_scaling_mat<- function(mat) {
  if ("matrix" %in% class(mat)) {
    write("Warning: input is matrix, output will be dgcMatrix")
  }

   ## Get the Dimensions
   n <- nrow(mat)

   ## Make a Diagonal Matrix of Column Su m
  D <- Matrix::sparseMatrix(i = 1:n,
                            j = 1:n,
                            x = Matrix::colSums(mat),
                            dims = c(n,n))

   ## Throw away explicit Zeroes
   D <- Matrix::drop0(D)

   ## Inverse the non-zero Values
   D@x <- 1/D@x

   ## Return the Diagonal Matrix
   return(D)

}

## ** Diagonal Scaling Matrix Inv
   #' Diagonal Factors of Sparse Matrix Inv
   #'
   #' Takes a sparse matrix
   #' and returns a diagonal matrix such that each column contains
   #' either colsum() or 0. (D in the Papers)
   #' multiplication with this matrix would have all columns
   #' sum to 1 (or zero if that column already summed to zero)
   #'
   #' @param mat Input Matrix, may either be a matrix, sparse matrix
   #' dgCMatrix or probably even a data frame, as long as `nrow` and `colSum`
   #' work on the data as expected it should work, output is always a
   #' sparse matrix of class dgCMatrix
   #'
   #' This should take the an adjacency matrix as input and the output
   #' can be multiplied by the original matrix to scale it to 1.
   #'
   #' A directed adjacency matrix should be such that A[i,j]
   #' describes the the weight
   #' of an edge from vertex j to vertex i (this is, unfourtunately
   #' the transpose of what igraph gives back.)
   #'
   #' @return A diagonal matrix such that multiplication by the original
   #' matrix will lead to all columns adding
   #' to 1 (or 0 in the case of a column of zeros)
   #' The output class is sparse matrix "dgCMatrix" from the Matrix package
   #' use `as.matrix()` in order to get a matrix type, I didn't put that
   #' Inside the function because it would have been confusing.
   #'
   #'
   #' @examples
   #' (mat <- matrix(1:9, nrow = 3))
   #' mat %*% create_sparse_diag_scaling_mat(mat) %>% colSum()
   #' mat %*% create_sparse_diag_scaling_mat(mat)
   #'
   #' > [1] 1 1 1 1 1
   #' @export
create_sparse_diag_sc_inv_mat<- function(mat) {
  if ("matrix" %in% class(mat)) {
    write("Warning: input is matrix, output will be dgcMatrix")
  }

   ## Get the Dimensions
   n <- nrow(mat)

   ## Make a Diagonal Matrix of Column Su m
  D <- Matrix::sparseMatrix(i = 1:n,
                            j = 1:n,
                            x = Matrix::colSums(mat),
                            dims = c(n,n))

   ## Throw away explicit Zeroes
   D <- Matrix::drop0(D)

   ## Return the Diagonal Matrix
   return(D)

}

## * Ergodic Graph
## ** Probability Transition Matrix

#' Adjacency to Probability Transition Matrix
#'
#' Takes an Adjacency matrix and scales each column to 1 or 0.
#'
#' The returned matrix will be such that each entry A[i,j] describes the
#' probability of travelling from vertex j to vertex i during a random
#' walk. (Note that column -> row is the transpose of what igraph returns)
#' which is row to column)
#'
#' @return the function returns a matrix of the form dgCMatrix from
#' from the Matrix package, wrap in as.matrix() if necessary
#'
#'
#' @param mat A matrix like object (either a matrix, sparse matrix or dataframe)
#'
#' @export
#'
#' @examples adj_to_probTrans(matrix(1:3, 3))
#'
#'
adj_to_probTrans <- function(mat) {
  D_inv  <-  create_sparse_diag_scaling_mat(mat)
      T  <- mat %*% D_inv
  return(T)
}

## ** TODO Stationary Distribution of Ergodic Graph

## * Power Walk
## ** Power Walk Probability Transition Matrix
#' Power Walk Probability Transition Matrix
#'
#' Takes an adjacency matrix and returns the transition
#' probability matrix using the Power Walk Method
#'
#' @section Power Walk Method
#'
#' The power walk method may take a non ergodic adjacency matrix with
#' any weights in the real numbers.
#'
#' @param beta is the ratio of probability between following an edge or
#' making a random jump in a graph. If beta is 0 the probability of
#' following an edge is no more likely than a randum jump.
#' As beta is made larger, the probability of making random jumps becomes smaller.
#' TODO if beta is zero the trans mat should be uniform.
#' Value defaults to 10
#'
#' @mat is a weighted adjacency matrix with real values.
#' may be of the class matrix or dgCmatrix
#'
#' @return The output matrix will be of the class dgCmatrix
#' @export
power_walk_prob_trans <- function(A, beta = 10) {
    n       <- nrow(A)
    B       <- beta^A

    DB_inv  <- create_sparse_diag_scaling_mat(B)
    T     <- B %*% DB_inv
    return(T)
}

## ** Power Walk Stationary Point

#' Power Walk Stationary Point
#'
#' Takes an adjacency matrix and returns the stationary
#' and returns the stationary distribution using the power walk method
#'
#' The returned vector will describe the probability of landing on any
#' given vertex of a graph during a random walk with a smoothing parameter.
#'
#' @section Power Walk Method
#' See \code{\link{power_walk_prob_trans}}
#'
#' @param beta
#' See \code{\link{power_walk_prob_trans}}
#'
#' probability of travelling from vertex j to vertex i during a random
#' walk. (Note that column -> row is the transpose of what igraph returns)
#' which is row to column)
#'
#' @return the function returns a matrix of the form dgCMatrix from
#' from the Matrix package, wrap in as.matrix() if necessary
#'
#'
#' @param mat A matrix like object (either a matrix, sparse matrix or dataframe)
#'
#' @export
#'
#' @examples adj_to_probTrans(matrix(1:3, 3))
#'
#'
power_walk_stationary_point  <- function(A, beta=10, eta) {
  ## I can't just rese the power_walk_prob_trans function because
  ## I need δ in both places, performing ColSums multiple times and having
  ## nested functions could also get very slow if this was inside a loop so
  ## it's simpler to just redefine it.
    A      <- Matrix::Matrix(A, sparse = TRUE)
    n      <- nrow(A)
    B      <- A
    B      <- beta^A   # Element Wise exponentiation

    Bo     <- A
    # These two approaches are equivalent
    Bo@x   <- beta^(A@x) -1   # This in theory would be faster
    # Bo     <- β^(A) -1
    # Bo     <- drop0(Bo)

## Create the Scaling Matrix to make row sums 1
##
  δB   <- 1/(Matrix::colSums(Bo)+n) # = 1/(Matrix::colSums(B))
  δBt  <- t(δB)
  DB   <- diag(δB)

## Create the Trans Prob Mat using Power Walk
  T <- Bo %*% DB

  n      <- nrow(A)
  p_new  <- rep(1/n, n)  # Uniform
  p      <- rep(0, n)    # Zero
# Implement the Loop

 while (sum(abs(p_new - p)) > eta) {
    (p <- as.vector(p_new)) # P should remain a vector
    sum(p <- as.vector(p_new)) # P should remain a vector
     p_new  <- T %*% p + rep(t(δB) %*% p, n)
  }

  return(as.vector(p))

}
## * Random Surfer
## ** Probability Transition Matrix
## This is the ordinary Transition Matrix because the RS is:
## RS = ɑT + (1-ɑ)B
## ** Stationary Distribution
## To keep everything sparse this has to be done differently
## than just iterating T <- TD
#' Random Surfer Stationary Distribution
#'
#' Takes an Adjacency matrix and returns the Stationary Distribution using Random Surfer
#'
#' The returned vector will describe the probability of landing at each
#' vertex during a random walk
#' Using the Random Surfer model.
#'
#' @param alpha the value of the smoothing parameter,
#'
#'   * 0 would indicate that every vertex is teleported to during a random walk
#'   * 1 would indicate that teleporting can't happen
#' @param A An adjacancy matrix such that A[i,j] indicates travel j -> i
#' @param dp how many decimal points of accuracy to use TODO make eta
#'
#' @export
#'
#' @examples adj_to_probTrans(matrix(1:3, 3))
#'
#'
random_surfer_stationary <- function(A, alpha=0.85, dp = 2) {
  T     <- adj_to_probTrans(A)
  A     <- Matrix::Matrix(A, sparse = TRUE)
  N     <- nrow(A)
  F     <- rep((1-alpha)/N, nrow(A))  ## A nx1 vector of (1-alpha)/N

  ## Solve using the power method
  p     <- rep(0, length.out = ncol(T)); p[1] <- 1
  p_new <- alpha*T %*% p + F

  ## use a Counter to debug
## **** TODO use eta not DP
  i <- 0
  while (sum(round(p, dp) != round(p_new, dp))) {
      p     <- p_new
      p_new <- alpha*T %*% p + F
  }

  return(as.vector(p))
}
