#' ST-DBSCAN
#' @article{birant2007st,
#'  title={ST-DBSCAN: An algorithm for clustering spatial--temporal data},
#'  author={Birant, Derya and Kut, Alp},
#'  journal={Data \& knowledge engineering},
#'  volume={60},
#'  number={1},
#'  pages={208--221},
#'  year={2007},
#'  publisher={Elsevier}
#'}
#' https://www.sciencedirect.com/science/article/pii/S0169023X06000218
#' 
#'
#' @importFrom sf st_read
#' @importFrom sf st_transform
#' @importFrom sf st_centroid
#' @importFrom sf st_geometry
#' @importFrom spdep poly2nb
#' @importFrom geosphere distm
#' @importFrom compiler::cmpfun
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end

library(geosphere)


#' calc ST-DBSCAN clustering.
#' @article{birant2007st,
#'  title={ST-DBSCAN: An algorithm for clustering spatial--temporal data},
#'  author={Birant, Derya and Kut, Alp},
#'  journal={Data \& knowledge engineering},
#'  volume={60},
#'  number={1},
#'  pages={208--221},
#'  year={2007},
#'  publisher={Elsevier}
#'}
#' https://www.sciencedirect.com/science/article/pii/S0169023X06000218
#' nb: neighbor list
#' minPts: minimum number of points
#' vals: as delta_eps list
#' type:
#'   "GridCell": For Grid data. Adjust nb size before calc.
#'   "Random": For normal(random point) data. not Adjustment.
st_dbscan <- function(nb = NULL,
                      vals = list(list(D = NULL, ## val list
                                       delta_eps = 5)),
                      type = c("GridCell", "Random"),
                      minPts = 10) {
  type <- match.arg(type)
  stopifnot(is.list(vals))
  stopifnot(!is.null(nb))
  
  ## set cluster number
  cluster <- 0

  ## set cluster column
  label <- rep(NA, times = length(nb))
  
  ## st-dbscan
  message("\nStart Clustering:  ", date())
  for (i in seq(1, length(nb))) {
    if (length(nb[[i]]) == 0) {
      label[i] <- "Noise"
    } else if (is.na(label[i]) || label[i] == "Noise") {
      ## 格子データを用いる場合
      if (type == "GridCell") {
        strClust <- do.call(c, lapply(vals, function(x) {
          nb[[i]][which(abs(x$D[i] - x$D[nb[[i]]]) <= x$delta_eps)]
        }))
      } else if (type == "Random") {
        strClust <- nb[[i]]
      }

      ## minPtsより隣接数が少ない場合は "Noise" 判定
      if (length(strClust) <= length(nb[[i]]) && length(nb[[i]]) < minPts) {
        label[i] <- "Noise"
      } else if (minPts <= length(strClust)  && length(strClust) <= length(nb[[i]])) {
        ## up to cluster number
        cluster <- cluster + 1
        label[i] <- as.character(cluster)
        message("\n --- Create Cluster : ", as.character(cluster))

        ## check cluster
        if (type == "GridCell") {
          label[strClust] <- as.character(cluster)
          label[setdiff(nb[[i]], strClust)] <- "Noise"
        } else if (type == "Random") {
          label[nb[[i]]] <- as.character(cluster)
        }
        
        ## check node and Add Cluster.
        num <- nb[[i]]
        while (length(num) != 0) {
          nodeList <- lapply(vals, function(x) {
            ## search cluster value
            m <- x$D[which(label == as.character(cluster))]

            ## check nb =====
            linkNb <- unique(do.call(c, nb[num]))
            if (type == "GridCell") {
              fitNb <- do.call(c, lapply(linkNb, function(l) {
                xd <- x$D[nb[[l]]]
                mu <- mean(c(m, xd), na.rm = TRUE)
                nbs <- nb[[l]][which(abs(mu - xd) <= x$delta_eps)]
                if (minPts <= length(nbs)) {
                  return(nbs)
                } 
              })) |> unique()
            } else if (type == "Random") {
              fitNb <- linkNb
            }
            checkNb <- setdiff(fitNb, num)
            rm(linkNb)
            rm(fitNb)
            
            ## Fit Cluster:: search node =====
            empNode <- checkNb[is.na(label[checkNb])]
            add <- c()
            for (n in empNode) {
              mu <- mean(c(m, x$D[add]), na.rm = TRUE)
              if (abs(mu - x$D[n]) <= x$delta_eps) {
                add <- c(add, n)
              }
              rm(mu)
            }
            return(add)
          })
          if (length(nodeList) != 0) {
            node <- duplicated2(nodeList)
            rm(nodeList)
          } else {
            node <- c()
          }
          if (length(node) != 0) {
            label[node] <- as.character(cluster)
            num <- node
          } else {
            num <- c()
          }
        }
      } else {
        label[i] <- "Noise"
      }
    }
  }
  message("\n", date(), "  Completed.")
  return(label)
}


#' create neighbor list
stnb <- function(x = double(),
                 y = double(),
                 time = NULL,
                 eps0 = 0,
                 eps1 = NULL,
                 eps2 = NULL,
                 method = c("geo", "euclidian"),
                 type = c("GridCell", "Random")) {
  stopifnot(!is.null(x))
  stopifnot(!is.null(y))
  stopifnot(!is.null(time))
  stopifnot(!is.null(eps1))
  stopifnot(!is.null(eps2))
  method <- match.arg(method)
  type <- match.arg(type)
  
  if (type == "GridCell") {
    len <- length(x)
    timediff <- lapply(time, function(t) {
      which(abs(time - t) <= eps2)
    })

    g <- lapply(seq(1, length(x), by = 1), function(i) {
      d <- dist2(x, y, x[i], y[i], method = method)
      which(eps0 < d & d <= eps1)
    })
    do.call(c, lapply(timediff, function(num) {
      lapply(g, function(gd) {
        unique(do.call(c, Map(`+`, list(gd), (num - 1) * len)))
      })
    }))
  } else if (type == "Random") {
    lapply(seq(1, length(x), by = 1), function(i) {
      d <- dist2(x, y, x[i], y[i], method = method)
      td <- abs(time - time[i])
      intersect(which(eps0 < d & d <= eps1), which(0 < td & td <= eps2))
    })
  }
}
dist2 <- function(x1, y1, x2, y2, method) {
  if (method == "euclidian") {
    sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
  } else if (method == "geo") {
    as.vector(geosphere::distm(cbind(x1, y1), cbind(x2, y2)) / 1000)
  }
}

#' 重複している値のみを反す。
duplicated2 <- function(data = list()) {
  if (length(data) == 1) {
    data[[1]]
  } else if (length(data) == 2) {
    data[[1]][data[[1]] %in% data[[2]]]
  } else {
    data[[1]][data[[1]] %in% duplicated2(data[-1])]
  }
}

geoSailing <- function(Lon1, Lat1, Lon2, Lat2,
                       A = 6378137.0, F = 1/298.257222101,
                       Unit = c("km", "m"), digit = 1) {
  Unit <- match.arg(Unit)
  if (Lon1 == Lon2 && Lat1 == Lat2) {
    return(0.0)
  }

  ## 扁平率
  B <- A * (1 - F)
  
  ## ラジアン角度に変換
  deg2rad <- function(x) {
    pi * x / 180
  }
  lonA <- deg2rad(Lon1)
  latA <- deg2rad(Lat1)
  lonB <- deg2rad(Lon2)
  latB <- deg2rad(Lat2)

  ## 化成緯度に変換
  phiA <- atan(B / A * tan(latA))
  phiB <- atan(B / A * tan(latB))

  ## 球面上の距離
  chi <- acos(sin(phiA) * sin(phiB) + cos(phiA) * cos(phiB) * cos(lonA - lonB))

  ## Lanbert-Andoyer補正
  delta_rho <- F / 8 * ((sin(chi) - chi) * (sin(phiA) + sin(phiB)) ^ 2 / cos(chi / 2) ^ 2 - (sin(chi) + chi) * (sin(phiA) - sin(phiB)) ^ 2 / sin(chi / 2) ^ 2)
  if (is.na(delta_rho)) {
    delta_rho <- 0.0
  }

  ## 測地線長
  rho <- A * (chi + delta_rho)

  if (Unit == "km") {
    return(round(rho / 1000, digits = digit))
  } else if (Unit == "m") {
    return(round(rho, digits = digit))
  }
}
