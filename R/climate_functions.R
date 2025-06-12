#' Altcorrect_point
#' apply altitudinal correction to point climate data
#' code adapted from microclimf:::.runmodel2Cpp
#'
#' @param climdata climate dataframe
#' @param dtmc elevation of input macroclimate data
#' @param altcor numeric; 1 = fixed lapse rate, 2 = humidity-dependent lapse rate
#'
#' @return altitduinally corrected climate data (adjusts temp and pressure) for input into
#' micropoint::runmicropoint


altcorrect_point = function(climdata, dtmc_p, elev_p, altcor = 2) {
  # Altitudinal correction applied
  dtmc_p[is.na(dtmc_p)]<-0
  h = nrow(climdata)
  pk = climdata$pres
  es = microclimf:::.satvap(climdata$temp)
  ea = es*climdata$relhum/100
  tdew = microclimf:::.dewpoint(ea, climdata$temp)
  psl<-pk/(((293-0.0065*rep(dtmc_p,h))/293)^5.26)
  pk = psl*(((293-0.0065*rep(elev_p,h))/293)^5.26)
  elevd = rep(dtmc_p - elev_p, h)

  if (altcor==1) {  # Fixed lapse rate
    tcdif<-elevd*(5/1000)
  } else { # Humidity-dependent lapse rate
    lr<-microclimf:::.lapserate(climdata$temp,ea,pk)
    tcdif<-lr*elevd
  }
  tc<-tcdif+climdata$temp
  # tc = new temp, pk = new pressure
  return(list(tc = tc, pk = pk))

}


#' climdata distance weight
#' Apply distance weighting to point climate data extracted from era5 grid
#'
#' @param lon decimal longitude
#' @param lat decimal latitude
#' @param climr list of spatrasters of climate variables (output of mcera5::extract_clima)
#' @param tme POSIX vector of hourly times represented in climr
#'
#' @return dataframe of climate variables for input to micropoint::runmicropoint
clim_distweight = function(lon, lat, climr, tme) {
  # distance to nearest era5 points
  focal = mcera5:::focal_dist(lon, lat)
  focal_collect = list()

  nms = names(climr)

  # get climate data at four cells nearest focal point
  clim_point = foreach(c = 1:length(climr), .combine = "cbind") %do% {
    df = terra::extract(climr[[c]], focal[,c('x', 'y')])
    df = df %>%  pivot_longer(2:(length(tme)+1), names_to = "obs_time", values_to = nms[c])
    if(c > 1) {
      df = df[,3]
    }
    df
  }

  iw = data.frame(ID = 1:4, inverse_weight = focal$inverse_weight)
  clim_point = clim_point %>% left_join(iw, by = "ID") %>%
    rename(neighbour = ID) %>%
    mutate(obs_time = rep(tme, times = 4),
           obs_time = as.POSIXlt(obs_time))

  clim_point <- clim_point %>%
    dplyr::group_by(obs_time)%>%
    dplyr::summarise_at(dplyr::vars(temp, relhum, pres, swdown, difrad,
                                    lwdown, windspeed, winddir, precip),
                        weighted.mean, w = dplyr::quo(inverse_weight)) %>%
    dplyr::mutate(timezone = lubridate::tz(obs_time))

  # convert negative precip to 0
  if(sum(clim_point$precip < 0) > 0) {
    clim_point$precip[clim_point$precip<0] = 0
    warning("Negative precipitation values converted to zero")
  }


  return(clim_point)

}
