#' Fast Potential Calulation
#'
#' @description stewart with cuttoff (limit) and parrallel
#' @param knownpts sf point, known pts
#' @param unknownpts sf point, known points
#' @param var variable name
#' @param fun friction function, "e" (expon) or "p" (pareto)
#' @param span span
#' @param beta beta
#' @param limit limit
#' @param ncl number of core
#' @param size size chunk
#' @return a vector
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
fp_fastpot <- function(knownpts, unknownpts, var = "v", fun = "e",
                       span = 2000, beta = 2,
                       limit = 10000, ncl = 3, size = 500){


  #   # knownpts = pt
  #   # unknownpts = grid
  #   # var = "v"
  #   # fun = "e"
  #   # span = 2000
  #   #
  #   # beta = 2
  #   # limit = 5000
  #   # size = 250
  # launch multiple cores
  if (missing(ncl)){
    ncl <- parallel::detectCores(all.tests = FALSE, logical = FALSE) - 1
  }
  cl <- parallel::makeCluster(ncl)
  doParallel::registerDoParallel(cl)


  # sequence to split unknowpts
  sequence <- unique(c(seq(1,nrow(unknownpts), size),nrow(unknownpts)+1))
  lseq <- length(sequence)-1

  # split unknownpts and put it on a list
  ml <- list()
  for  (i in 1:lseq){
    ml[[i]] <- unknownpts[(sequence[i]):(sequence[i+1]-1),]
  }

  pot <- foreach::`%dopar%`(foreach::foreach(unknownpts = ml,
                                             .packages = c('sf'),
                                             .combine = c, .inorder = FALSE),
                            {
                              gbuf <- st_buffer(st_centroid(st_geometry(unknownpts)), limit)
                              x <- st_intersects(gbuf, knownpts, prepared = TRUE)
                              kgeom <- st_geometry(knownpts)
                              ugeom <- st_centroid(st_geometry(unknownpts))
                              kgeom <- matrix(unlist(kgeom), ncol = 2, byrow = TRUE)
                              ugeom <- matrix(unlist(ugeom), ncol = 2, byrow = TRUE)
                              eucledian_simple <- function(from, to){
                                sqrt( (from[1] - to[1])^2 + (from[ 2] - to[ 2])^2 )
                              }
                              if(fun =="e"){
                                alpha  <- log(2) / span ^ beta
                                fric <- function(alpha,matdist,beta){exp(- alpha * matdist ^ beta)}
                              }
                              if(fun == "p"){
                                alpha  <- (2 ^ (1 / beta) - 1) / span
                                fric <- function(alpha,matdist,beta){(1 + alpha * matdist) ^ (-beta)}
                              }

                              v <- knownpts[[var]]

                              l <- vector("list", nrow(ugeom))

                              for (i in seq_along(l)){
                                kn <- kgeom[unlist(x[i]),,drop=FALSE]
                                un <- ugeom[i,]
                                matdist <- apply(kn, 1, eucledian_simple, un)
                                un <- sum(v[unlist(x[i])] * fric(alpha, matdist, beta), na.rm = TRUE)
                                l[[i]] <- un
                              }
                              unlist(l)
                            })
  parallel::stopCluster(cl)
  return(pot)
}
