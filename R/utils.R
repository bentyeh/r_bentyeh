#' Evaluate an object or expression while suppressing output.
#'
#' @param x An object or expression
#' @param quiet_messages logical. Suppress messages.
#' @param quiet_warnings logical. Suppress warnings.
#'
#' @references Based on Hadley Wickham's suggestion in
#'   \url{https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html}.
#'
#' @export
quiet <- function(x, quiet_messages = TRUE, quiet_warnings = TRUE) {
  wrapper_messages <- identity
  if (quiet_messages) {
    wrapper_messages <- suppressMessages
  }
  wrapper_warnings <- identity
  if (quiet_warnings) {
    wrapper_warnings <- suppressWarnings
  }
  sink(tempfile())
  on.exit(sink())
  invisible(wrapper_warnings(wrapper_messages(force(x))))
}

#' Wait on child processes.
#'
#' Currently only supports POSIX-compliant systems.
#'
#' @export
wait <- if (.Platform$OS.type == "unix") {
  Rcpp::cppFunction(
    "void wait() {int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};}",
    includes = "#include <sys/wait.h>"
  )
} else {
  function() {
    stop("wait() not implemented for Windows platforms.")
  }
}

#' Get the number of CPUs the current process can use.
#'
#' @param silent \code{logical}. Warn if unable to accurately detect number of available CPUs.
#'
#' @return \code{integer} or \code{NA} if the number of CPUs in the system is undetermined.
#'
#' @seealso \link[parallel]{mcaffinity}, \link[parallel]{detectCores}.
#'
#' @export
get_available_cpus <- function(silent = FALSE) {
  try(
    {
      getFromNamespace("mcaffinity", "parallel")
      return(length(parallel::mcaffinity()))
    },
    silent = TRUE
  )
  if (identical(silent, FALSE)) {
    warning("Unable to detect number of available CPUs. Returning number of total CPUs.")
  }
  return(parallel::detectCores())
}

#' Convert an array to a data frame.
#'
#' @param ar An array.
#'
#' @return \code{data.frame. dim = c(prod(dim(ar)), length(dim(ar)) + 1)}.
#' \itemize{
#'   \item Rows: One row per element of the array. Rows are "entered" into the data frame in
#'     \href{https://en.wikipedia.org/wiki/Row-_and_column-major_order}{column major order}.
#'   \item Columns: One column per dimension of the array, plus one column for the values in the
#'     array. \itemize{
#'     \item Column names are retained from names of array dimensions (\code{names(dimnames(ar)}).
#'       Where non-existant, column names are given as "d#", where # is the index of the dimension.
#'     \item For each column, if the corresponding array dimension was named, then the column is of
#'       class \code{character}. Otherwise, the column is of class \code{numeric}.
#'     \item The last column (the values from the array) is named "value".
#'   }
#' }
#'
#' @seealso \code{\link[tidyr]{pivot_longer}}, \url{https://stackoverflow.com/a/42810479}.
#'
#' @export
arrayToDf <- function(ar) {

  dims <- dim(ar)
  ndims <- length(dims)

  if (is.null(dimnames(ar))) {
    nullDims <- 1:ndims
    nullDimNames <- 1:ndims
  } else {
    nullDims <- which(sapply(dimnames(ar), is.null))
    if (is.null(names(dimnames(ar)))) {
      nullDimNames <- 1:ndims
    } else {
      nullDimNames <- which(sapply(names(dimnames(ar)), function(x) identical(x, "")))
    }
  }
  namedDims <- setdiff(1:ndims, nullDims)

  df <- as.data.frame.table(ar)
  df[namedDims] <- lapply(df[namedDims], as.character)
  df[nullDims] <- lapply(df[nullDims], as.numeric)

  colnames(df)[nullDimNames] <- paste0("d", nullDimNames)
  colnames(df)[ncol(df)] <- "value"

  return(df)
}

#' Convert a data frame to an array.
#'
#' @param df \code{data.frame. dim = c(nr, nc)}. The first (nc - 1) columns represent dimensions.
#'   The last column gives values in the array.
#' @param dimOrders \code{list} of \code{vector}. Ordering (i.e., factor levels) of each dimension
#'   in the array. Must be a fully named list (arbitrary length) or fully unnamed list (\code{length
#'   = nc}). If a column ordering is given for a column, any values present in that column but not
#'   in the ordering assumes a value of NA. Any columns for which an ordering is not given is sorted
#'   lexicographically. The class of each vector should match the class of the corresponding column.
#'
#' @return \code{array. mode = mode(df[[length(df)]])}. Mode of array is the same as the mode of the
#'   values column (i.e., last column) in the data frame.
#'
#' @references \url{https://stackoverflow.com/a/9617424}, \url{https://stackoverflow.com/a/46129338}.
#'
#' @examples
#' dfToArray(CO2) # CO2 data frame exported from datasets package
#'
#' @export
dfToArray <- function(df, dimOrders = NULL) {

  nDim <- ncol(df) - 1
  stopifnot(is.null(dimOrders) || !is.null(names(dimOrders)) || length(dimOrders) == nDim)

  df <- unique(df)

  # convert columns in df to factors
  if (is.null(dimOrders) || !is.null(names(dimOrders))) {
    namedCols <- intersect(names(df[1:nDim]), names(dimOrders))
    unnamedCols <- setdiff(names(df[1:nDim]), namedCols)
    for (c in namedCols) {
      df[[c]] <- factor(df[[c]], levels = dimOrders[[c]])
    }
    for (c in unnamedCols) {
      df[[c]] <- factor(df[[c]])
    }
  } else {
    for (c in 1:nDim) {
      df[[c]] <- factor(df[[c]], levels = dimOrders[[c]])
    }
  }

  # initialize array
  ar <- array(
    dim = sapply(df[1:nDim], function(c) length(levels(c))),
    dimnames = sapply(df[1:nDim], levels)
  )

  # input values into array
  ar[do.call(cbind, df[1:nDim])] <- df[[nDim + 1]]

  return(ar)
}

#' Order each dimension of array by hierarchical clustering.
#'
#' @param ar \code{array}.
#' @param dims \code{vector, integer}. Dimensions to order. If NULL, all dims are ordered.
#' @param metric \code{character}. Distance metric by which to compare hyperplanes along dimension.
#'   Hyperplanes are flattened into vectors for comparison. Support \code{"cor"} ((1 - r) / 2) or
#'   any \link[stats]{dist} method.
#' @param method \code{character}. Any \link[stats]{hclust} agglomeration method.
#' @param cor_use \code{character}. Only applicable if \code{metric = "cor"}. \code{use} argument of
#'   \link[stats]{cor}.
#' @param cor_method \code{character}. Only applicable if \code{metric = "cor"}. \code{method}
#'   argument of \link[stats]{cor}.
#' @param return_hclust \code{logical}. Return list of hclust objects for ordered dimensions.
#'
#' @return \code{list. length = length(dim(ar))}. \itemize{
#'   \item If \code{return_hclust = FALSE}: Each element is a vector giving the permutation of the
#'     corresponding array dimension.
#'   \item If \code{return_hclust = TRUE}: Each element is an \code{hclust} object for the
#'     corresponding dimension, or NULL if that dimension was not ordered.
#' }
#'
#' @seealso \link[stats]{cor}, \link[stats]{dist}, \link[stats]{hclust}.
#'
#' @examples
#' ar <- matrix(c(1, 1, 1, 2, 2, 2, 3, 4, 5), nrow = 3)
#' orderArray(ar)
#'
#' @export
orderArray <- function(
  ar,
  dims = NULL,
  metric = "cor",
  method = "complete",
  cor_use = "pairwise.complete.obs",
  cor_method = "pearson",
  return_hclust = FALSE
) {
  nDims <- length(dim(ar))
  if (is.null(dims)) {
    dims <- 1:nDims
  }

  # validate that dims to order are present in array
  stopifnot(all(dims %in% 1:nDims))

  # validate return_hclust
  stopifnot(is.logical(return_hclust) && length(return_hclust) == 1)

  # create named dimOrder list
  if (!return_hclust) {
    dimOrder <- vector("list", nDims)
    if (!is.null(names(dimnames(ar)))) {
      names(dimOrder) <- sapply(names(dimnames(ar)), function(x) ifelse(identical(x, ""), NULL, x))
    }
  }
  hclusts <- vector("list", nDims)

  # order specified dims
  for (d in dims) {
    # construct matrix to pass into dist()
    # - each row represents a hyperplane along dimension d
    nValues <- dim(ar)[d]
    mat <- matrix(nrow = nValues, ncol = prod(dim(ar)[-d]))

    # flatten each hyperplane in dimension d into a row vector
    for (value in 1:nValues) {
      mat[value, ] <- as.vector(indexArray(ar, d, value))
    }

    if (metric == "cor") {
      distMat <- (1 - cor(t(mat), use = cor_use, method = cor_method))/2
      distMat[is.na(distMat)] <- 1
      distMat <- as.dist(distMat)
    } else {
      distMat <- dist(mat, method = metric)
    }
    hclusts[[d]] <- hclust(distMat, method = method)
    if (!return_hclust) {
      dimOrder[[d]] <- hclusts[[d]]$order
      # alternatively, order.dendrogram(as.dendrogram(hclusts[[d]]))
    }
  }

  if (return_hclust) {
    return(hclusts)
  } else {
    # maintain original orders for other dims
    unOrderedDims <- setdiff(1:nDims, dims)
    dimOrder[unOrderedDims] <- lapply(unOrderedDims, function(x) 1:(dim(ar)[x]))
    return(dimOrder)
  }
}

#' Select along one dimension in multidimensional array
#'
#' @param x \code{array}.
#' @param dim \code{integer}. Dimension along which to select one hyperplane.
#' @param value \code{integer}. Hyperplane to select from given dimension.
#' @param drop \code{logical}. For matrices and arrays. If TRUE the result is coerced to the lowest
#'   possible dimension.
#' @param verbose \code{logical}. Print out equivalent command using `[` operator.
#'
#' @return \code{array} or \code{vector}.
#'
#' @references \url{https://stackoverflow.com/a/14502298}.
#'
#' @examples
#' A <- array(data = 1:24, dim = c(2, 3, 4), dimnames = list(c("x1", "x2"), c("y1", "y2", "y3"), NULL))
#' identical(indexArray(A, 1, 1), A[1, , ])
#' identical(indexArray(A, 1, 1, drop = FALSE), A[1, , , drop = FALSE])
#' identical(indexArray(A, 1, 2), A[2, , ])
#' identical(indexArray(A, 2, 3), A[, 3, ])
#'
#' @export
indexArray <- function(x, dim, value, drop = TRUE, verbose = FALSE) {

  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- value

  # Generate the call to [
  call <- as.call(c(
    list(as.name("["), quote(x)),
    indices,
    list(drop = drop)))

  if (verbose) {
    print(sub("^x", deparse(substitute(x)), capture.output(print(call))))
  }

  eval(call)
}
