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
  clim_point = foreach::foreach(c = 1:length(climr), .combine = "cbind") %do% {
    df = terra::extract(climr[[c]], focal[,c('x', 'y')])
    df = df %>%  tidyr::pivot_longer(2:(length(tme)+1), names_to = "obs_time", values_to = nms[c])
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


#' Process era5 data using extract_clima
#' @description
#' Forked from mcera5::extract_clima. Adjustments include, faster cropping to extent
#' and addition of land-sea mask separately from the era5.nc file
#'
#' @param nc character vector containing the path to the nc file. Use the
#' `build_era5_request` and `request_era5` functions to acquire an nc file with
#' the correct set of variables. Data within nc file must span the period
#' defined by start_time and end_time.
#' @param long_min minimum longitude of the grid for which data are required (decimal
#' @param long_max maximum longitude of the grid for which data are required (decimal
#' degrees, -ve west of Greenwich Meridian).
#' @param lat_min minimum latitude of the grid for which data are required (decimal
#' @param lat_max maximum latitude of the grid for which data are required (decimal
#' degrees, -ve south of the equator).
#' @param start_time a POSIXlt or POSIXct object indicating the first day or hour
#' for which data are required. Encouraged to specify desired timezone as UTC (ERA5
#' data are in UTC by default), but any timezone is accepted.
#' @param end_time a POSIXlt or POSIXct object indicating the last day or hour for
#' which data are required. Encouraged to specify desired timezone as UTC (ERA5
#' data are in UTC by default), but any timezone is accepted.
#' @param land_sea_mask spatraster; era5 lsm layer
#' @param dtr_cor logical value indicating whether to apply a diurnal temperature
#' range correction to air temperature values. Default = `TRUE`.
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction. Default = 1.285, based on calibration against UK Met Office
#' observations.
#' @param format specifies what microclimate package extracted climate data will
#' be used for. Data will be formatted accordingly. Default is "microclimf".
#' Options: "microclima", "NicheMapR", "microclimc", "microclimf", "micropoint".
#' Note: of these options, only "microclimf" accepts as input an array of climate
#' variables. For all other models you will need to iterate through each spatial
#' point to run the model
#'
extract_clima_2 = function(nc, long_min, long_max, lat_min, lat_max, start_time, end_time,
                        land_sea_mask, dtr_cor = TRUE, dtr_cor_fac = 1.285,
                        format = "microclimf") {

  # Open nc file for error trapping
  nc_dat = ncdf4::nc_open(nc)

  ## Error trapping ---------------

  # Specify the base date-time, which differs between the CDS versions, and the
  # first and last timesteps from timeseries, which has different names across
  # versions
  # Extract time dimension from data queried from either old or new CDS
  timedim <- mcera5:::extract_timedim(nc_dat)
  # Find basetime from units
  base_datetime <- as.POSIXct(gsub(".*since ", "", timedim$units), tz = "UTC")
  # Extract time values
  nc_datetimes <- c(timedim$vals)
  # If units in hours, multiply by 3600 to convert to seconds
  nc_datetimes <- nc_datetimes * ifelse(
    grepl("hours", timedim$units), 3600, 1
  )
  # Find first timestep
  first_timestep <- nc_datetimes[1]
  # Find last timestep
  last_timestep <- utils::tail(nc_datetimes, n = 1)

  # Confirm that start_time and end_time are date-time objects
  if (any(!class(start_time) %in% c("Date", "POSIXct", "POSIXt", "POSIXlt")) |
      any(!class(end_time) %in% c("Date", "POSIXct", "POSIXt", "POSIXlt"))) {
    stop("`start_time` and `end_time` must be provided as date-time objects.")
  }
  # Confirm that start_time and end_time are same class of date-time objects
  if (any(class(start_time) != class(end_time))) {
    stop("`start_time` and `end_time` must be of the same date-time class.")
  }

  # Check if start_time is after first time observation
  start <- base_datetime + first_timestep
  if (start_time < start) {
    stop("Requested start time is before the beginning of time series of the ERA5 netCDF.")
  }

  # Check if end_time is before last time observation
  end <- base_datetime + last_timestep
  if (end_time > end) {
    stop("Requested end time is after the end of time series of the ERA5 netCDF.")
  }

  if (long_max <= long_min) {
    stop("Maximum longitude must be greater than minimum longitude.")
  }

  if (lat_max <= lat_min) {
    stop("Maximum longitude must be greater than minimum longitude.")
  }

  if (abs(long_min) > 180 | abs(long_max) > 180 |
      abs(lat_min) > 90 | abs(lat_max) > 90) {
    stop("Coordinates must be provided in decimal degrees (longitude between -180 and 180, latitude between -90 and 90).")
  }

  # Check if requested coordinates are in spatial grid
  if (long_min < (min(nc_dat$dim$longitude$vals) - 0.125) |
      long_max > (max(nc_dat$dim$longitude$vals) + 0.125)
  ) {
    long_out <- TRUE
  } else {
    long_out <- FALSE
  }

  if (lat_min < (min(nc_dat$dim$latitude$vals) - 0.125) |
      lat_min > (max(nc_dat$dim$latitude$vals) + 0.125)
  ) {
    lat_out <- TRUE
  } else {
    lat_out <- FALSE
  }
  # close nc file
  ncdf4::nc_close(nc_dat)
  if(long_out & lat_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (both longitude and latitude out of range).")
  }
  if(long_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (longitude out of range).")
  }
  if(lat_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (latitude out of range).")
  }

  if (lubridate::tz(start_time) != lubridate::tz(end_time)) {
    stop("start_time and end_time are not in the same timezone.")
  }

  if (lubridate::tz(start_time) != "UTC" | lubridate::tz(end_time) != "UTC") {
    warning("provided times (start_time and end_time) are not in timezone UTC (default timezone of ERA5 data). Output will be provided in timezone UTC however.")
  }

  # Specify hour of end_time as last hour of day, if not specified
  if (lubridate::hour(end_time) == 0) {
    end_time <- as.POSIXlt(paste0(lubridate::year(end_time), "-",
                                  lubridate::month(end_time), "-",
                                  lubridate::day(end_time),
                                  " 23:00"), tz = lubridate::tz(end_time))
  }

  # Check that `format` is an accepted value
  if (!format %in% c("NicheMapR", "microclima", "microclimc", "micropoint", "microclimf")) {
    stop("Argument `format` must be one of the following values: `NicheMapR`, `microclima`, `microclimc`, `micropoint`, `microclimf`")
  }

  if (dtr_cor == TRUE & !is.numeric(dtr_cor_fac)) {
    stop("Invalid diurnal temperature range correction value provided.")
  }

  tme <- as.POSIXct(seq(start_time,
                        end_time, by = 3600), tz = lubridate::tz(end_time))

  ## Load in and subset netCDF variables --------------

  # Rename mean surface rad variables, which changed names on the CDS in late Jan 2025
  nc_dat <- ncdf4::nc_open(nc)
  varnames <- names(nc_dat$var)

  if ("msnlwrf" %in% varnames) {
    varname_list <- c("t2m", "d2m", "sp", "u10" , "v10",  "tp", "tcc", "msnlwrf",
                      "msdwlwrf", "fdir", "ssrd", "lsm")
  } else {
    varname_list <- c("t2m", "d2m", "sp", "u10" , "v10",  "tp", "tcc", "avg_snlwrf",
                      "avg_sdlwrf", "fdir", "ssrd", "lsm")
  }


  crop_fast = function(r, e) {

    colmin = terra::colFromX(r, e[1])
    colmax = terra::colFromX(r, e[2])
    rowmin = terra::rowFromY(r, e[4]) # 4 because spatrasters indexed from topleft corner
    rowmax = terra::rowFromY(r, e[3])

    rcrop = r[rowmin:rowmax, colmin:colmax, drop = F]
    return(rcrop)

  }

  # Remove land-sea mask from varname list if using pre-downloaded layers
  # We will pull in the land-sea mask separately
  fun = function(nc, v, e) {
    r <- terra::rast(nc, subds = v)
    r <- r[[as.POSIXct(nc_datetimes, tz = "UTC", origin = "1900-01-01") %in% tme]]
    names(r) <- tme
    r <- crop_fast(r, e)
    return(r)
  }
  if(!("lsm" %in% varnames)) {
    varname_list = grep("lsm", varname_list, value = T, invert = T)

    var_list = list()

    e = terra::ext(long_min, long_max, lat_min, lat_max)
    var_list[["t2m"]] = fun(nc, "t2m", e)
    var_list[["d2m"]] = fun(nc, "d2m", e)
    var_list[["sp"]] = fun(nc, "sp", e)
    var_list[["u10"]] = fun(nc, "u10", e)
    var_list[["v10"]] = fun(nc, "v10", e)
    var_list[["tp"]] = fun(nc, "tp", e)
    var_list[["tcc"]] = fun(nc, "tcc", e)
    var_list[["msnlwrf"]] = fun(nc, "msnlwrf", e)
    var_list[["msdwlwrf"]] = fun(nc, "msdwlwrf", e)
    var_list[["fdir"]] = fun(nc, "fdir", e)
    var_list[["ssrd"]] = fun(nc, "ssrd", e)
    var_list[["lsm"]] = fun(nc, "lsm", e)



    # for(v in varname_list) {
    #   # subset down to desired time period
    #   # terra::time() not identifying time data of ERA5 data from new CDS, so
    #   # use nc_datetimes
    #   r <- terra::rast(nc, subds = v)
    #
    #   r <- r[[as.POSIXct(nc_datetimes, tz = "UTC", origin = "1900-01-01") %in% tme]]
    #   # Name layers as timesteps
    #   names(r) <- tme
    #
    #   # Subset down to desired spatial extent
    #   # r <- terra::crop(r, terra::ext(long_min, long_max, lat_min, lat_max))
    #   r <- crop_fast(r, terra::ext(long_min, long_max, lat_min, lat_max))
    #   var_list[[v]] = r
    # }

    # var_list = lapply(varname_list, function(v) {
    #   # subset down to desired time period
    #   # terra::time() not identifying time data of ERA5 data from new CDS, so
    #   # use nc_datetimes
    #   r <- terra::rast(nc, subds = v)
    #
    #   # Subset down to desired spatial extent
    #   r <- terra::crop(r, terra::ext(long_min, long_max, lat_min, lat_max))
    #   # r <- crop_fast(r, terra::ext(long_min, long_max, lat_min, lat_max))
    #
    #   r <- r[[as.POSIXct(nc_datetimes, tz = "UTC", origin = "1900-01-01") %in% tme]]
    #   # Name layers as timesteps
    #   names(r) <- tme
    #
    #   return(r)
    # })


    # Add land-sea mask into the raster list
    lsm = crop_fast(land_sea_mask, var_list[[1]])
    lsm = terra::resample(lsm, var_list[[1]][[1]], method = "average")
    var_list[[12]] = lsm

    names(var_list) <- c(varname_list, "lsm")

  } else {
    var_list = list()
    for(v in varname_list) {
      if (v == "lsm") {
        # only need one timestep for land-sea mask
        r <- terra::rast(nc, subds = v)[[1]]
      } else {
        # For all others, subset down to desired time period
        # terra::time() not identifying time data of ERA5 data from new CDS, so
        # use nc_datetimes
        r <- terra::rast(nc, subds = v)
        r <- r[[as.POSIXct(nc_datetimes, tz = "UTC") %in% tme]]
        # Name layers as timesteps
        names(r) <- tme
      }

      # Subset down to desired spatial extent
      #r <- terra::crop(r, terra::ext(long_min, long_max, lat_min, lat_max))
      r <- crop_fast(r, terra::ext(long_min, long_max, lat_min, lat_max))
      var_list[[v]] = r
    }
    # var_list <- lapply(varname_list, function(v) {
    #   if (v == "lsm") {
    #     # only need one timestep for land-sea mask
    #     r <- terra::rast(nc, subds = v)[[1]]
    #   } else {
    #     # For all others, subset down to desired time period
    #     # terra::time() not identifying time data of ERA5 data from new CDS, so
    #     # use nc_datetimes
    #     r <- terra::rast(nc, subds = v)
    #     r <- r[[as.POSIXct(nc_datetimes, tz = "UTC") %in% tme]]
    #     # Name layers as timesteps
    #     names(r) <- tme
    #   }
    #
    #   # Subset down to desired spatial extent
    #   #r <- terra::crop(r, terra::ext(long_min, long_max, lat_min, lat_max))
    #   r <- crop_fast(r, terra::ext(long_min, long_max, lat_min, lat_max))
    #   return(r)
    # })

    names(var_list) <- c(varname_list)
  }

  if ("msnlwrf" %in% names(var_list)) {
    var_list$avg_snlwrf <- var_list$msnlwrf
    var_list$avg_sdlwrf <- var_list$msdwlwrf
  }

  t2m <- var_list$t2m
  d2m <- var_list$d2m
  sp <- var_list$sp
  u10 <- var_list$u10
  v10 <- var_list$v10
  tcc <- var_list$tcc
  avg_snlwrf <- var_list$avg_snlwrf
  avg_sdlwrf <- var_list$avg_sdlwrf
  fdir <- var_list$fdir
  ssrd <- var_list$ssrd
  prec <- var_list$tp * 1000 # convert from mm to metres
  lsm <- var_list$lsm
  temperature <- t2m - 273.15 # kelvin to celcius
  ## Coastal correction ----------

  # see land-sea mask here:
  # https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

  # Only conduct if there are cells with proximity to water
  if (any(terra::values(lsm) < 1)) {
    # Calculate daily average
    # Indices to associate each layer with its yday
    ind <- rep(1:(dim(temperature)[3]/24), each = 24)
    # Average across days
    tmean <- terra::tapp(temperature, ind, fun = mean, na.rm = T)
    # Repeat the stack 24 times to expand back out to original timeseries
    tmean <- rep(tmean, 24)
    # Sort according to names so that the stack is now in correct order: each
    # daily mean, repeated 24 times
    # your pasted command properly sorts the names, X1 to X365 (or X366)
    tmean <- tmean[[paste0("X", sort(rep(seq(1:(dim(temperature)[3]/24)), 24)))]]
    m <- (1 - lsm) * dtr_cor_fac + 1
    tdif <- (temperature - tmean) * m
    temperature <- tmean + tdif
  }

  humidity <- mcera5:::humfromdew(d2m - 273.15,
                                  temperature,
                                  sp)
  windspeed = sqrt(u10^2 + v10^2)
  windspeed = mcera5::windheight(windspeed, 10, 2)
  winddir = (terra::atan2(u10, v10) * 180/pi + 180)%%360
  cloudcover = tcc * 100
  netlong = abs(avg_snlwrf) * 0.0036 # Convert to MJ/m^2/hr
  downlong = avg_sdlwrf * 0.0036 # Convert to MJ/m^2/hr
  uplong = netlong + downlong
  emissivity = downlong/uplong
  jd = mcera5::julday(lubridate::year(tme),
              lubridate::month(tme),
              lubridate::day(tme))
  rad_dni = fdir * 0.000001 # Convert form J/m^2 to MJ/m^2
  rad_glbl = ssrd * 0.000001 # Convert form J/m^2 to MJ/m^2
  ## si processing -----------------
  # use t2m as template of dimensions for iterating through
  si <- t2m
  foo <- t2m[[1]]
  coords <- as.data.frame(terra::crds(t2m[[1]]))
  # # Get a vect (points) of the scene coords
  # foo <- coordinates(raster(foo))
  # points <- terra::vect(foo, crs = terra::crs(t2m, proj = T), type = "points")
  # # And coerce points to dataframe
  # points <- terra::geom(points, df = TRUE)
  # Create a template with dimensions (x * y) and length(tme)
  out <- array(NA, dim = c(nrow(coords), length(tme)))
  for (i in 1:nrow(coords)) {
    out[i,] <- mcera5::siflat(lubridate::hour(tme),
                      lat = coords$y[i],
                      long = coords$x[i],
                      jd)
  }
  si <- terra::setValues(si, out)

  # Calc rad_dif
  rad_dif = rad_glbl - rad_dni * si
  ## szenith processing -----------------
  # use t2m as template of dimensions for iterating through
  szenith <- t2m
  # Create a template with dimensions (x * y) and length(tme)
  out <- array(NA, dim = c(nrow(coords), length(tme)))
  for (i in 1:nrow(coords)) {
    out[i,] <- mcera5::solalt(lubridate::hour(tme),
                      lat = coords$y[i],
                      long = coords$x[i],
                      julian = jd)
  }
  szenith <- terra::setValues(szenith, out)

  ## Add timesteps back to names ---------------
  # Only necessary for temperature at the moment, all other variables retain info
  names(temperature) <- names(t2m)
  terra::time(temperature) <- tme

  ## Format of output use ----------
  ## Equivalent of hourlyncep_convert
  if (format %in% c("microclimc", "micropoint", "microclimf")) {
    pres <- sp / 1000
    ## Convert humidity from specific to relative
    relhum <- humidity
    terra::values(relhum) <- mcera5::converthumidity(h = terra::as.array(humidity),
                                             intype = "specific",
                                             tc  = terra::as.array(temperature),
                                             pk = terra::as.array(pres))$relative
    relhum[relhum > 100] <- 100

    ## NOTE: extract_clima()  DOES NOT CURRENTLY CONVERT ACCUMULATED MEASUREMENTS
    # ACROSS THE HOUR TO A MEAN ON THE HOUR, AS IS PERFORMED IN rad_calc(). THIS
    # YIELDS SMALL BUT NON-MARGINABLE DIFFERENCES IN VALUES
    message("note: extract_clima() does not currently convert accumulated measurements across the hour to a mean on the hour, as is performed in extract_clim(). This yields small but non-marginable differences in values.")

    raddr <- (rad_dni * si)/0.0036 # convert back to W/m^2
    difrad <- rad_dif/0.0036 # convert from MJ/hr to W/m^2
    swrad <- raddr + difrad
    downlong <- downlong / 0.0036 # convert from MJ/hr to W/m^2
  }
  # Return list - SpatRasters now wrapped as won't store as list if saved otherwise'
  if (format %in% c("micropoint", "microclimf")) {
    return(list(
      obs_time = tme,
      temp = terra::wrap(temperature),
      relhum = terra::wrap(relhum),
      pres = terra::wrap(pres),
      swdown = terra::wrap(swrad),
      difrad = terra::wrap(difrad),
      lwdown = terra::wrap(downlong),
      windspeed = terra::wrap(windspeed),
      winddir = terra::wrap(winddir),
      precip = terra::wrap(prec)
    ))
  }
  if (format == "microclimc") {
    return(list(
      obs_time = tme,
      temp = terra::wrap(temperature),
      relhum = terra::wrap(relhum),
      pres = terra::wrap(pres),
      swrad = terra::wrap(swrad),
      difrad = terra::wrap(difrad),
      skyem = terra::wrap(emissivity),
      windspeed = terra::wrap(windspeed),
      winddir = terra::wrap(winddir),
      precip = terra::wrap(prec)
    ))
  }
  if (format %in% c("microclima", "NicheMapR")) {
    return(list(
      obs_time = tme,
      temperature = terra::wrap(temperature),
      humidity = terra::wrap(humidity),
      pressure = terra::wrap(sp),
      windspeed = terra::wrap(windspeed),
      winddir = terra::wrap(winddir),
      emissivity = terra::wrap(emissivity),
      netlong = terra::wrap(netlong),
      uplong = terra::wrap(uplong),
      downlong = terra::wrap(downlong),
      rad_dni = terra::wrap(rad_dni),
      rad_dif = terra::wrap(rad_dif),
      szenith = terra::wrap(szenith),
      cloudcover = terra::wrap(tcc)
    ))
  }

}
