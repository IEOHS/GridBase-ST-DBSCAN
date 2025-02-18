#' https://rion778.hatenablog.com/entry/20121203/1354461231
#'
#' @importFrom units set_units

#' create airmass data
#' wpdata <- with(wp, {
#'   do.call(rbind, lapply(names(data), function(x) {
#'     cbind(date = as.Date(x), data[[x]])
#'   }))
#' })
#' wpdata <- with(wpdata, {
#'   ## potential temprature [K]
#'   pt <- local({
#'     pface <- which_press(pface)
#'     eval(parse(text = paste0(
#'       "apply(cbind(",
#'       paste0("potentialTemp(t = P_", pface,
#'              "_temperature, p = ", pface, ")", collapse = ", "),
#'       ") , 1, mean, na.rm = TRUE)"
#'     )))
#'   })
#'   ## specific-humidity [kg-water/kg-water]
#'   q <- local({
#'     pface <- which_press(pface)
#'     eval(parse(text = paste0(
#'       "apply(cbind(",
#'       paste0("TVR(t = P_", pface, "_temperature, ",
#'              "RH = P_", pface, "_relative_humidity, ",
#'              "P = ", pface, ", output = \"specific-humidity\")",
#'              collapse = ", "),
#'       "), 1, mean, na.rm = TRUE)"
#'     )))
#'   })
#'   airmass <- airmass(sfc = wp$geo,
#'                      ## clustering by st-dbscan
#'                      cl = st_dbscan(
#'                        ## calc neighbor list
#'                        nb = stnb(type = "GridCell",
#'                                  poly = wp$geo,
#'                                  times = unique(date),
#'                                  stNeigh = 1),
#'                        values = list(list(v = pt,
#'                                           delta_eps = 2),
#'                                      list(v = q,
#'                                           delta_eps = 0.0002)),
#'                        minPts = 18, type = "GridCell"),
#'                      lower = sqrt(100000 ^ 2 + 100000 ^ 2),
#'                      upper = 40000000)
#'   data.frame(date, pt, q, airmass)
#' })
#' ## return
#' wp$data <- local({
#'   d <- lapply(names(wp$data), function(x) {
#'     subset(wpdata, date == x)
#'   })
#'   names(d) <- as.character(names(wp$data))
#'   return(d)
#' })
#' wp$units <-  list(date = set_one_units("day"),
#'                   potential_temperature = set_one_units("K"),
#'                   specific_humidity = set_one_units("kg / kg"),
#'                   airmass = NA)
#' rm(wpdata)

#' get airmass data set
#' d <- readRDS("../data/WeatherData/MSM-P_700,850hPa_daily_2020-05-31.rds")
#' ds <- readRDS("../data/WeatherData/MSM_daily_2020-05-31.rds")
#' df <- do.call(rbind, lapply(names(d$data), function(x) {
#'   t <- ISOdate(year = substring(x, 1, 4),
#'                month = substring(x, 6, 7),
#'                day = substring(x, 9, 10),
#'                hour = substring(x, 12, 13),
#'                tz = "Asia/Tokyo")
#'   within(sf::st_as_sf(d$data[[x]]["P_700_temperature"],
#'                       SP = ds$data[[x]]["sea_level_pressure"],
#'                       geometry = d$geo), {
#'     Date <- t
#'     PT <- potentialTemp(t = P_700_temperature, p = 700)
#'   })
#' }))



#' Goff–Gratch equation 飽和水蒸気圧
#' https://www.metsoc.jp/tenki/pdf/1988/1988_02_0115.pdf
## ggfunc <- function(t = 20, T1 = 273.16) {
##     T <- t + 273.15
##     e <- 10.79574 * (1 - T1 / T) -
##         5.02800 * log(T / T1, base = 10) +
##         1.50475 * 10 ^ (-4) * (1 - 10 ^ (-8.2969 * (T / T1 - 1))) +
##         0.42873 * 10 ^ (-3) * (10 ^ (4.76955 * (1 - T1 / T)) - 1) +
##         0.78614
##     units::set_units(10 ^ e, hPa)
## }
## #' tetensの式による飽和水蒸気圧の計算
## tetens <- function(t = 20,
##                    a = 7.5,
##                    b = 237.3) {
##     units::set_units(6.1078 * 10 ^ (a * t / (t + 273.15)),
##                      hPa)
## }


#' check airmass.
#' sfc = sf::sfc object.
#' cl = cluster label.
#' lower = smallest check geospatial size. (default size 100000 [m])
#' upper = biggest geospatial size. (default size 40000000 [m])
#' if cluster size(bbox) under your set lower then unset airmass.
airmass <- function(sfc = NULL,
                    cl = c(),
                    lower = sqrt(100000^2 + 100000^2),
                    upper = 40000000) {
  label <- unique(cl)
  label <- label[label != "Noise"]
  nsfc <- length(sfc)
  force(lower)

  ## check airmass.
  message("----------------------------------------")
  message("Check airmass.", "\n")
  for (x in label) {
    n <- which(cl == x)
    bbox <- sf::st_bbox(sfc[ifelse(n %% nsfc == 0, nsfc, n %% nsfc)])
    clsize <- geosphere::distGeo(bbox[1:2], bbox[3:4])
    if (clsize < lower) {
      message("Cluster: ", x, " < ", lower, "[m].  Delete cluster.")
      cl[n] <- "Noise"
    } else if (clsize > upper) {
      message("Cluster: ", x, " > ", upper, "[m].  Delete cluster.")
      cl[n] <- "Noise"
    }
  }
  ## rename cluster number
  message("----------------------------------------")
  message("Create airmass number.", "\n")
  label <- unique(cl)
  label <- label[label != "Noise"]
  for (i in seq(1, length(label))) {
    cl[which(cl == label[i])] <- paste0("airmass_", i)
    message("Create airmass. number: ", i)
  }
  return(cl)
}

goffgratch <- function(t = 20,      ## 気温
                       T1 = 373.16, ## 沸点
                       est = 1013.25   ## 基準気圧[hPa]
                       ) {
    T <- t + 273.15
    e <- -7.90298 * (T1 / T - 1) +
        5.02808 * log(T1 / T, base = 10) -
        1.3816 * 10 ^ (-7) * (10 ^ (11.344 * (1 - T / T1)) - 1) +
        8.1328 * 10 ^ (-3) * (10 ^ (-3.49149 * (T1 / T - 1)) - 1) +
        log(est, base = 10)
    units::set_units(10 ^ e, hPa)
}

#' 水蒸気量の計算
wv <- function(t = 20, ## 気温
               RH = 50,## 相対湿度
               et      ## 飽和水蒸気圧
               ) {
    ## 飽和水蒸気量
    w <- 217 * as.numeric(et) / (t + 273.15)
    units::set_units(w * RH / 100,
                     g * m^(-3))
}

#' 混合比の計算
mixRatio <- function(t = 20,       ## 気温（℃）
                     RH = 50,      ## 相対湿度(%)
                     P = 1013.25,  ## 気圧
                     W = NULL      ## 飽和水蒸気圧
                     ) {
    T <- t + 273.15
    W <- as.numeric(W)
    RH <- RH / 100
    qw <- 622 * (RH * W) / (P - RH * W)
    attr(qw, "units") <- "g-water / kg-dryair"
    qw
}

#' 比湿
spHum <- function(qw = NULL) {
    sh <- qw / (1 + qw)
    attr(sh, "units") <- "g-water / kg-wetair"
    sh
}

#' 仮温度
vt <- function(t = 20, qw = NULL) {
    T <- t + 273.15
    units::set_units((1 + 0.61 * qw) * T,
                     K)
}

#' 仮温位
vpt <- function(vt = NULL,
                p = NULL,
                p0 = 1000,
                Rd = 287,
                Cp = 1004
                ) {
    ## tv * (p0 / p) ^ (Rd / Cp)
    units::set_units(vt * (p0 / p) ^ (Rd / Cp),
                     K)
}

## 温位
potentialTemp <- function(t = double(), # temp℃
                          p = double(), # pressure hPa
                          p0 = 1000,    # standard pressure
                          Rd = 287,     # 乾燥空気の気体定数 J/K/kg
                          Cp = 1004    # J/K/kg
                          ) {
  ## units::set_units((t + 273.15) * (p / p0)^(-Rd / Cp), "K")
  ret <- (t + 273.15) * (p / p0)^(-Rd / Cp)
  attr(ret, "unit") <- "K"
  return(ret)
}

## temp-vapor ratio
TVR <- function(t = 20,
                RH = 50,
                P = 1013,
                output = c("saturation-vapor-pressure",
                           "mixing-ratio",
                           "virtual-temp",
                           "specific-humidity",
                           "temp-water_vapor-ratio",
                           "all")
                ) {
    output <- match.arg(output)
    
    ## 飽和水蒸気圧の計算
    et <- goffgratch(t = t, est = P)

    ## 水蒸気量の計算
    ##w <- wv(t = t, RH = RH, et = et)

    ## 混合比の計算
    mr <- mixRatio(t = t, RH = RH, P = P, W = et)

    ## 比湿の計算
    sh <- spHum(qw = mr)
    
    ## 仮温度の計算
    vt <- vt(t = t, qw = as.numeric(mr))

    ## output
    switch(output,
           "saturation-vapor-pressure" = et,
           "mixing-ratio" = mr,            
           "virtual-temp" = vt,            
           "specific-humidity" = sh,             
           "temp-water_vapor-ratio" = vt / w,
           "all" = list(saturation_vapor_pressure = et,
                        mixing_ratio = mr,            
                        virtual_temp = vt,            
                        specific_humidity = sh)
           )
}

