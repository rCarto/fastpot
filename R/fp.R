#' @title Fast Potential Calulation
#' @description stewart with cuttoff (limit) and parrallel.
#' @param x set of observations to compute the potentials from, sf points.
#' @param y set of points for which the potentials are computed.
#' @param var names of the variables in \code{x} from which potentials are computed.
#' Quantitative variable with no negative values.
#' @param fun spatial interaction function. Options are "p"
#' (pareto, power law) or "e" (exponential).
#' For pareto the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' For "exponential" the interaction is defined as:
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user
#' (\code{beta} and \code{span}).
#' @param span distance where the density of probability of the spatial
#' interaction function equals 0.5.
#' @param beta impedance factor for the spatial interaction function.
#' @param limit maximum distance used to retrieved \code{x} points, in map units.
#' @param ncl number of clusters. \code{ncl} is set to
#' \code{parallel::detectCores() - 1} by default.
#' @param size \code{fp_fastpot} splits \code{y} in smaller chunks and
#' dispatches the computation in \code{ncl} cores, \code{size} indicates the
#' size of each chunks.
#' @return If only one variable is computed a vector is returned, if more than
#' one variable is computed a matrix is returned.
#' @export
#' @importFrom sf st_buffer st_centroid st_geometry st_intersects
#' @examples
#' \dontrun{
#' library(sf)
#' library(cartography)
#' library(SpatialPosition)
#'
#' mtq <- st_read(system.file("gpkg/mtq.gpkg", package="cartography"))
#' grid <- st_make_grid(mtq, cellsize = 500)
#' grid <- st_sf(id = 1:length(grid), geom = grid)
#' set.seed(1)
#' size = 50000
#' pt <- st_sample(grid, size = size)
#' pt <- st_sf(id = 1:length(pt), geom = pt)
#' cc <- st_coordinates(pt)
#' pt$v <- runif(n = size, 10000, 20000)
#' pt$v <- pt$v + cc[,1] + cc[,2] / 1000
#'
#' system.time(
#'   pot2 <- mcStewart(pt, grid, varname = "v", span = 2000, beta = 2,
#'                     returnclass = "sf", cl = 3, size  = 100)
#' )
#' system.time(
#'   pot <- fp_fastpot(pt, grid, var = c('v'), fun="e",beta = 2, limit = 10000,
#'                     ncl = 3, size = 100)
#' )
#'
#' plot(pot,pot2$OUTPUT, asp = 1)
#' lm(pot ~pot2$OUTPUT)$coeff
#'
#' grid$pot <- pot
#' par(mfrow = c(1,2))
#' choroLayer(pot2, var = "OUTPUT", border = NA)
#' choroLayer(grid, var = "pot", border = NA)
#' }
fp_fastpot <- function(x, y, var = "v", fun = "e",
                       span = 2000, beta = 2,
                       limit = 10000, ncl = 3, size = 500){

  # launch multiple cores
  if (missing(ncl)){
    ncl <- parallel::detectCores(all.tests = FALSE, logical = FALSE) - 1
  }
  cl <- parallel::makeCluster(ncl)
  doParallel::registerDoParallel(cl)

  # data simplification
  xsfc <- st_geometry(x)
  kgeom <- matrix(unlist(xsfc), ncol = 2, byrow = TRUE)



  v <- as.matrix(x= x[,var,drop = TRUE])
  print(dim(v))
  ysfc <- st_centroid(st_geometry(y))

  # sequence to split unknowpts
  ny <- nrow(y)
  sequence <- unique(c(seq(1, ny, size), ny + 1))
  lseq <- length(sequence) - 1

  # split unknownpts and put it on a list
  ml <- list()
  for  (i in 1:lseq){
    ml[[i]] <- ysfc[(sequence[i]):(sequence[i+1]-1)]
  }

  # dispatch
  pot <- foreach::`%dopar%`(
    foreach::foreach(ysfc = ml,
                     .packages = 'sf',
                     .combine = c,
                     .inorder = FALSE),
    {
      # FUNS
      eucledian_simple <- function(from, to){
        sqrt( (from[1] - to[1])^2 + (from[ 2] - to[ 2])^2 )
      }
      if(fun =="e"){
        alpha  <- log(2) / span ^ beta
        fric <- function(alpha,matdist,beta){
          exp(- alpha * matdist ^ beta)
        }
      }
      if(fun == "p"){
        alpha  <- (2 ^ (1 / beta) - 1) / span
        fric <- function(alpha,matdist,beta){
          (1 + alpha * matdist) ^ (-beta)
        }
      }

      # Buffer limit
      gbuf <- st_buffer(ysfc, limit)
      inter <- st_intersects(gbuf, xsfc, prepared = TRUE)

      # data transformation
      ugeom <- matrix(unlist(ysfc), ncol = 2, byrow = TRUE)

      # go through each y
      l <- vector("list", nrow(ugeom))
      for (i in seq_along(l)){
        kindex <- unlist(inter[i])
        kn <- kgeom[kindex,,drop=FALSE]
        un <- ugeom[i,]
        matdist <- apply(kn, 1, eucledian_simple, un)
        un <- apply(X = v[kindex, , drop = FALSE], MARGIN = 2,
                    FUN = function(x){sum(x * fric(alpha, matdist, beta), na.rm = TRUE)})
        l[[i]] <- un
      }
      unlist(l)
    }
  )
  # stop parralel
  parallel::stopCluster(cl)
  print("kkk")
  if(length(var)==1){
    pot <- as.numeric(pot)
  }else{
    pot <- matrix(pot, ncol = length(var), byrow = T, dimnames = list(NULL, var))
  }
  return(pot)
}

