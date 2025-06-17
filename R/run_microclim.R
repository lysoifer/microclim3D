#' Prep_micropoint
#' prepare input to micropoint model
#'
#' @param e spatextent of area to model
#' @param start start date character "YYYY-MM-DD"
#' @param end end date; character "YYYY-MM-DD"
#' @param era_path file path to downloaded era5 data (should be a .nc file from get_era5)
#' @param landcover_path path to landcover data
#' @param soilpath path to soil data
#' @param elev_path path to elevation data
#' @param dtmc_path path to era5 elevation raster
#'
prep_micropoint = function(e, start, end, era_path, landcover_path, soilpath, elev_path, dtmc_path) {
  # extend requested extent (so distance weighting works)
  e_extend = extend(e, 1)
  # extract climate data
  climr = extract_clima(era_path,
                        start_time = as.POSIXlt(start, tz = "UTC"),
                        end_time = as.POSIXlt(end, tz = "UTC"),
                        e_extend[1], e_extend[2], e_extend[3], e_extend[4],
                        format = "micropoint")
  climr = lapply(climr, terra::unwrap)

  # VEGETATION PARAMETERS
  # landcover
  lc = rast(landcover_path)
  # reclassify landcover to categories used in microclimf
  lc = reclass_landcover(lc[[1]])
  # reproject landcover to match projection of gedi data ("epsg:4326)
  lc = .reproj_bigr(lc, "epsg:4326", e, method = "near")

  # need to deal with PAI
  # will have to estimate monthly change in lai based on modis or landsat and apply offset to gedi

  # EXTENT SHOULD BE PER 1 X 1 DEGREE TILE
  # get habitat-based vegetation data
  load("scripts/package_devel/data/globclim.rda") # global climate variables (see documentation in microclimf)
  tme = climr$obs_time
  vegp = microclimf::vegpfromhab(
    habitats = lc,
    clump0 = F,
    tme = as.POSIXlt(tme),
    lat = mean(e[3], e[4]),
    long = mean(e[1], e[2])
  )
  vegp$pai = rast(vegp$pai, crs = crs(vegp$hgt), extent = ext(vegp$hgt))
  #vegp = lapply(vegp, unwrap)
  vegp$clump = rast(vegp$clump, crs = crs(vegp$hgt), extent = ext(vegp$hgt))

  # get leaf reflectance - can work on this for improving the model??
  # albedo = rast(list.files(albedopath, pattern = ".tif$", full.names = T), lyrs = seq(1, 23, 2))
  # eproj_alb = project(e, "epsg:4326", crs(albedo))
  # albedo = crop(albedo, eproj_alb)
  # albedo = project(albedo, "epsg:4326")
  # albedo = resample(albedo, vegp$hgt)
  # albedo = mean(albedo, na.rm = T)

  # can the model take monthly leafr values?
  # or use same months for pai and albedo - see example in package - errors here
  #leafr = leafrfromalb(pai = pai$m_1, x = vegp$x, alb = albedo)

  #vegp$leafr = leafr

  # GROUND PARAMETERS
  soil = rast(soilpath)
  soil = .reproj_bigr(soil, "epsg:4326", e, method = "bilinear")
  soil = resample(soil, vegp$hgt)
  soil = soil/10

  # extend elevation to remove edge effects when calculating aspect and slope
  elev = rast(elev_path)
  eproj = project(extend(e, c(1,1,1,1)), "epsg:4326", crs(elev))
  elev = crop(elev, eproj)
  elev = project(elev, "epsg:4326")
  elev = project(elev, vegp$hgt, align_only = T)
  #elev = resample(elev, vegp$hgt)

  asp = terra::terrain(elev, "aspect", unit = "degrees")
  slp = terra::terrain(elev, "slope", unit = "degrees")

  # era5 elevation for altitudinal correction
  dtmc = rast(dtmc_path)
  crs(dtmc) = "epsg:4326"

  return(list(climr = climr, vegp = vegp, soil = soil, elev = elev, asp = asp, slp = slp, dtmc = dtmc))
}


#' Run micropoint
#' Iteratively run micropoint model for all points in gedi
#'
#' @param tme vector of POSIXlt representing times of climdata
#' @param gedi datatable of gedi shots to model microclimates
#' @param climr climate rasters - output of mcera5::extract_clima
#' @param vegp output of microclimf::vegpfromhab
#' @param soil 3 layer soil raster with sand, silt, and clay from soilgrids database
#' @param elev elevation raster
#' @param asp aspect raster
#' @param slp slope raster
#' @param dtmc raster of elevations of era5 climate data (in climr)
#' @param method either "temporal_year", "temporal_month, or "vertprof";
#' temporal_year models microclimate over entire year at specified reqhgt and estimates monthly variation in PAI based on MODIS LAI data.
#' The estimate of monthly variation in based on offsets from MODIS in the month when the gedi shot was taken and assumes a constant offset throughout the year.
#' temporal_month models hourly microclimates at specified reqhgts for the month when the gedi shot was taken. If reqhgt = NA, then
#' the reqhgts are taken to be 0.15, 2, and then 5 to floor(top of canopy) at 5 m intervals.
#' vertprof models vertical microclimate profile at specified hour and for specified climate variable over n canopy layers.
#' There are still some minor bugs with temporal_year and vertprof. The estimate of monthly variation in vertical profiles is unlikely to be reliable
#' because greenup happens at different times in the understory vs. in the canopy. Vertprof ends up estimating microclimate above the height given to the model.
#' I'm not sure why this is happening.
#' @param plotout one of "tair", "relhum", or "tleaf" - the climate variable to be output for the vertical profile model
#' only needed if method = "vertprof"
#' @param n the number of canopy nodes used to estimate the vertical profile
#' only needed if method = "vertprof"
#' @param fout .csv file to save model outputs
#' @param modis_path file path to modis lai data
#' @param reqhgts vector of hgts above ground at which to model microclimates.
#' default NA models microclimate at 0.15, 2, and then 5 to floor(canopy height) at 5 m intervals
#' @param vertpai_method determines the method used to calculate the vertical distribution of foliage.
#' Can be "pai" or "pavd". If "pai", then vertical foliage profiles are calculated
#' by finding the difference between cumulative PAI at different heights in the canopy. For example,
#' the pai in the 0-1m voxel would be calculated as the difference between the pai at 0m (cumulative pai)
#' and the pai at 1m (where this represents the PAI from the top of the canopy down to 1m). If "pavd", then
#' vertical pai is calculated as proportion of total pai where proportions are determined by the contribution
#' of pavd in a given layer to the total pavd. These two methods will return slightly different vertical profiles.

#'
run_micropoint = function(tme, gedi, climr, vegp, soil, elev, asp, slp, dtmc,
                          method, plotout, n, maxiter, fout, modis_path, reqhgts = NA, vertpai_method = "pai") {
  # if(method == "vertprof") {
  #   vertmat = matrix(ncol = 7, nrow = 0)
  #   colnames(vertmat) = c("z", "tair", "canopy_height", "elev", "lon", "lat", "doy")
  #   write.table(vertmat, fout)
  # }

  # Unpack spatrasters
  # if we have multiple years of climate data
  if(is(climr[[1]], "list")) {
    for(a in 1:length(climr)) {
      climr[[a]] = lapply(climr[[a]], unwrap)
    }
  } else { # if there is only one year of climate data
    climr = unwrap(climr)
  }

  # unpack vegp
  vegp = lapply(vegp, unwrap)
  soil = unwrap(microin$soil)
  elev = unwrap(microin$elev)
  asp = unwrap(microin$asp)
  slp = unwrap(microin$slp)
  dtmc = unwrap(microin$dtmc)

  # Loop through GEDI points and calculate vertical microclimate gradient for each
  v = foreach(i = 1:nrow(gedi), .combine = rbind) %do% {
    print(i)

    # If we are modelling for only the month in which the GEDI point was observed
    # Then we need to pull out the year and month of the observation
    # and select the climate data for the correct year
    if(method == "temporal_month") {
      gedi[, year:= year(as.Date(date))]
      pyear = gedi[i, year]
      pmonth = gedi[i, month]
      climri = climr[[as.character(pyear)]]
      tmei = climri$obs_time
      climri = climri[2:10]
    }


    # extract climate data from raster at requested point
    # and apply distance weighting correction

    clim_point = clim_distweight(lon = as.numeric(gedi[i, "lon_lm_a0"]), lat = as.numeric(gedi[i,"lat_lm_a0"]),
                                 climri, tmei) %>%
      dplyr::select(!timezone) %>%
      mutate(obs_time = format(obs_time, format = "%Y-%m-%d %H:%M:%S")) %>%
      as.data.frame()

    # extract dtmc and elevation at point
    dtmc_p = terra::extract(dtmc, gedi[i,c("lon_lm_a0", "lat_lm_a0")], ID = F)[1,1]
    elev_p = as.numeric(gedi[i, "elev_lm_a0"])

    # apply altitudinal correction
    altcor = altcorrect_point(clim_point, dtmc_p, elev_p, altcor = 2)
    clim_point$temp = altcor[["tc"]]
    clim_point$pres = altcor[["pk"]]

    # vegetation parameters
    h = as.numeric(gedi[i, "rh_100_a0"]) # height
    pai = as.numeric(gedi[i, "pai_a0"]) #pai

    # get vertical pai profile
    pai_z = as.numeric(gedi[i, pai_a0:pai_l21])
    pavd = as.numeric(gedi[i,pavd_0_5:pavd_95_100]) # pavd

    paii = .pai_vertprofile(pai_z, h, vertpai_method)

    # If we want to model climate over an entire year, we need monthly estimates of pai
    # MODIS provides temporally resolved PAI, but GEDI is a single time point
    # The way I have approached it here is to model temporal variation in GEDI PAI
    # assuming a constant offset with MODIS
    # I additionally assume that the temporal variation is constant across the vertical profile
    # which is probably not a good assumption as understory greenup is typically offset from
    # canopy greenup
    # I also should be estimating temporal variability in clumpiness, which I am not currently doing

    # get vertical profile by season
    if(method == "temporal_year") {
      modis = rast(modis_path)
      pai_z_season = pai_seasonality(gedi[i,], xy_crs = "epsg:4326", modis)
      pai = pai_z_season[[1]]
      paii = pai_z_season[[2]]
    }

    # test if height is greater than 5m (errors when there is only one paiz value - not sure why)
    #print(pai_z)
    print(h)

    # extract values from vegp raster to get veg parameters at point
    # pai and height are derived from GEDI
    # remaining is derived from auto-generated params based on habitat type and lat/lon
    vegp_p = vegp_point(vegp, lon = as.numeric(gedi[i, "lon_lm_a0"]),
                        lat = as.numeric(gedi[i, "lat_lm_a0"]),
                        pai, h)

    # soil parameters
    # get soil parameters
    asp_p = terra::extract(asp$aspect, gedi[i, c("lon_lm_a0", "lat_lm_a0")])[1,2]
    slp_p = terra::extract(slp$slope, gedi[i, c("lon_lm_a0", "lat_lm_a0")])[1,2]

    # get ground parameters based on soilgrids and average values from micropoint::soilparams
    ground_p = .ground_point(soil, asp_p, slp_p,
                             as.numeric(gedi[i, "lon_lm_a0"]),
                             as.numeric(gedi[i, "lat_lm_a0"]))

    # apply dtr_correct
    clim_point_correct = micropoint:::dtr_correct(clim_point, zref = 2,
                                                  lat = gedi[i, lat_lm_a0],
                                                  long = gedi[i, lon_lm_a0])

    # model hourly microclimate over one year
    # takes into account temporal variation in PAI by modeling by month and changing PAI with each month based on changes in MODIS PAI
    if(method == "temporal_year") {
      # NEEDS SOME WORK, SEE JUNE 9 2025 IN LAB NOTEBOOK
      mout = foreach(m=1:12, .combine = "rbind") %do% {
        vegp_p$pai = pai[m]
        hrs = which(month(clim_point$obs_time)==m)
        out = micropoint::runpointmodel(
          clim_point_correct[hrs,],
          reqhgt = reqhgt,
          vegp = vegp_p,
          paii = paii[,m],
          groundp = ground_p,
          lat = as.numeric(gedi[i, "lat_lm_a0"]),
          long = as.numeric(gedi[i, "lon_lm_a0"])
        )
      }

    }

    # model microclimates at specified reqhgts for the month of the gedi shot
    if(method == "temporal_month") {
      #m = month(as.Date(gedi[i,date]))
      hrs = which(month(clim_point$obs_time)==pmonth)
      req = reqhgts
      if(is.na(sum(reqhgts))) {req = c(0.15, 2, seq(5,h - h%%5, 5), h-0.1)}
      # pai is currently constant
      mout = foreach(hi = req, .combine = rbind) %do% {
        out = micropoint::runpointmodel(
          clim_point_correct[hrs,],
          reqhgt = hi,
          vegp = vegp_p,
          paii = paii,
          groundp = ground_p,
          lat = as.numeric(gedi[i, "lat_lm_a0"]),
          long = as.numeric(gedi[i, "lon_lm_a0"])
        )
        out$canopy_hgt = vegp_p$h
        out$reqhgt = hi
        out
      }
    }

    pai_test = pai_z[pai_z>0]
    if(method == "vertprof" & h >= 6 & length(pai_test)>1 & min(pai > 0.2)) {
      # NEEDS A BIT OF WORK SEE JUNE9, 2025 IN LAB NOTEBOOK
      hour_yr = gedi[i, doy]*24-12 # noon on day of year that gedi shot was taken
      tair = c()
      height = c()
      hour = c()
      for(i in seq(hour_yr,hour_yr + 1,1)) {
        vertprof <- micropoint::plotprofile(clim_point_correct, hr = hour_yr, plotout = plotout,
                                            vegp_p, paii = paii, ground_p,
                                            lat = as.numeric(gedi[i, "lat_lm_a0"]),
                                            long= as.numeric(gedi[i, "lon_lm_a0"]),
                                            maxiter = maxiter, n = length(paii))
        tair = c(tair, vertprof$var)
        height = c(height, vertprof$z)
        hour = c(hour, rep(i, length(vertprof$z)))
      }

      mat = cbind(z = vertprof$z, var = vertprof$var,
                  canopy_height = vegp_p$h,
                  elev = elev_p,
                  lon = gedi[i,lon_lm_a0],
                  lat = gedi[i, lat_lm_a0],
                  doy = gedi[i, doy])
      colnames(mat)[2] = plotout

      # mat

      #vertmat = rbind(vertmat, mat)
      #write.csv(vertmat, file = fout, row.names = F)
    }
    #toc()
    mout$shot_num = gedi[i, shot_num]
    mout
  }
  #return(vertmat)
  write.csv(v, file = fout, row.names = F)
  return(v)
}


# need to return canopy height and canopy pai
# also would be nice to find a way of returning vertical veg profiles used in the model
# test plots
# df = v %>%
#   mutate(mday = mday(obs_time)) %>%
#   group_by(shot_num, reqhgt, mday) %>%
#   summarise(tair.minday = min(tair)) %>%
#   group_by(shot_num, reqhgt) %>%
#   summarise(tmin_cold = mean(tair.minday))
#   # group_by(reqhgt) %>%
#   # summarise(tmin_cold.mean = mean(tmin_cold), tmin_cold.sd = sd(tmin_cold))
#   #
#
# ggplot(df, aes(tmin_cold, reqhgt, group = as.factor(shot_num))) +
#   geom_line() +
#   theme_bw()
#
#
# df2 = v2 %>%
#   mutate(mday = mday(obs_time)) %>%
#   group_by(shot_num, reqhgt, mday) %>%
#   summarise(tair.minday = min(tair)) %>%
#   group_by(shot_num, reqhgt) %>%
#   summarise(tmin_cold = mean(tair.minday))
#
# ggplot(df2, aes(tmin_cold, reqhgt, group = as.factor(shot_num))) +
#   geom_line() +
#   theme_bw()
#
#
