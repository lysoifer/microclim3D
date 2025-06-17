#' Reproject big raster
#' Crop and reproject large raster to specified extent and projection of extent
#'
#' @param r spatraster to crop and reproject
#' @param proj_e projection of spatextent
#' @param e spatextent to project to
#' @param method method for reprojection (see terra::project method for options)
.reproj_bigr = function(r, proj_e, e, method) {
  e_extend = extend(e, 0.5)
  eproj = terra::project(e_extend, proj_e, crs(r))
  r = crop(r, eproj)
  r = project(r, proj_e, method = method)
  r = crop(r, e_extend)
  r = mask(r, e_extend)

  return(r)
}

#' PAI vertical profile
#' Calculate pai vertical profile from gedi data
#' Values add up to total PAI (for input into micropoint paii argument)
#'
#' @param pai_z vector of PAI values at height intervals from 0 to top of canopy (as in gedi data)
#' @param h height of canopy
#' @param vertpai_method can be "pai" or "pavd". If "pai", then vertical foliage profiles are calculated
#' by finding the difference between cumulative PAI at different heights in the canopy. For example,
#' the pai in the 0-1m voxel would be calculated as the difference between the pai at 0m (cumulative pai)
#' and the pai at 1m (where this represents the PAI from the top of the canopy down to 1m). If "pavd", then
#' vertical pai is calculated as proportion of total pai where proportions are determined by the contribution
#' of pavd in a given layer to the total pavd. These two methods will return slightly different vertical profiles.
#'
#' @return PAI profile at 1 m intervals from ground to top of canopy. Values add up to total PAI in gedi data
.pai_vertprofile = function(pai_z, h, vertpai_method) {

  # sum of vertical pai profile must equal total pai for the model to run

  if(vertpai_method == "pai") {
    # need to include zero here because pai_a0 does not equal pai0-5
    h_int = c(0,1,seq(6,floor(h)+5, 5))
    h_int2 = seq(0,floor(h), 1)

    pai_z = pai_z[1:length(h_int)]

    # apply a monotonically decreasing spline
    monospline = splinefun(h_int, pai_z, method = "monoH.FC")
    pai_zspline = monospline(h_int2)
    pai_z2 = pai_zspline - c(pai_zspline[2:length(pai_zspline)],0)
  }


  # We add zero to the list of PAIz values because pai at the top of the canopy = 0 (see paiz description in gedi documentation)

  # Try using proportion of pavd
  if(vertpai_method == "pavd") {
    h_int = seq(1,floor(h) + 5,5)
    h_int2 = seq(1,floor(h),1)
    pavd = pavd[1:length(h_int)]
    splinepavd = splinefun(h_int, pavd, method = "monoH.FC")
    pavdspline = splinepavd(h_int2)
    pavdprop = pavdspline/sum(pavdspline)
    pai_z2 = pai*pavdprop
  }
  return(pai_z2)
}


#' Get point vegetation data
#' Get vegetation data at point based on habitat variables and manually specified height and pai
vegp_point = function(vegp, lon, lat, pai, h) {
  x = terra::extract(vegp$x, data.frame(lon, lat))[1,2]
  clump = terra::extract(vegp$clump, data.frame(lon, lat))[1,2] # can replace with value from clumpiness raster or use single month
  lref = terra::extract(vegp$leafr, data.frame(lon, lat))[1,2]
  ltra = terra::extract(vegp$leaft, data.frame(lon, lat))[1,2]
  leafd = terra::extract(vegp$leafd, data.frame(lon, lat))[1,2]
  em = 0.97 # recommended value by micropoint vignette
  gsmax = terra::extract(vegp$gsmax, data.frame(lon, lat))[1,2]
  q50 = 100 # recommended value by micropoint vignette

  vegp_p = list(h = h, pai = pai, x = x, clump = clump, lref = lref,
                ltra = ltra, leafd = leafd, em = em, gsmax = gsmax, q50 = q50)

  return(vegp_p)
}

#' Ground_point
#' Make groundparams object for input into micropoint::runmicropoint
#' @param soil 3 layer spatraster of clay, sand, and silt content (from soilgrids)
#' @param asp_p aspect (degrees) at lon, lat
#' @param slp_p slope (degrees) at lon, lat
#' @param lon longitude decimal degrees
#' @param lat latitude decimal degrees
.ground_point = function(soil, asp_p, slp_p, lon, lat) {
  soilclass = terra::extract(soil, data.frame(lon, lat))
  colnames(soilclass) = c("ID", "CLAY", "SAND", "SILT")

  # classify soil texture according to USDA classification
  soilclass = TT.points.in.classes(soilclass, class.sys = "USDA.TT")

  colnames(soilclass) = c(
    "Clay", "Silty clay", "Sandy clay",
    "Clay loam", "Silty clay loam", "Sandy clay loam",
    "Loam", "Silty loam", "Sandy loam",
    "Silt", "Loamy sand", "Sand"
  )

  soilclass = names(which(soilclass[1,]==1))
  #soilnum = which(microclimf::soilparameters$Soil.type==soilclass)

  soilp = micropoint::soilparams[soilparams$Soil.type == soilclass,]

  groundparams_p = list(
    gref = 0.15, # ground reflectance - microctools automatically sets gref to 0.15
    slope = slp_p, # slope of ground surface (deg)
    aspect = asp_p, # aspect of ground surface
    em = 0.97, # emissivity of ground surface
    rho = soilp$rho, # soil bulk density (mg/m3)
    Vm = soilp$Vm, # volumetric mineral fraction of soil
    Vq = soilp$Vq, # Volumetric quartz fraction of soil
    Mc = soilp$Mc, # Mass fraction of clay
    b = soilp$b, # Shape parameter fo campbell soil moisture
    Psie = -soilp$psi_e, # Matric potential
    Smax = soilp$Smax, # Volumetric water content
    Smin = soilp$Smin, # Residual water content
    alpha = soilp$alpha, # shape parameter of the van Genuchten model
    n = soilp$n, # pore size distribution
    Ksat = soilp$Ksat # saturated hydraulic conductivity
  )

}


#' reclass_landcover
#' Reclassify landcover types to those specified in microclimf
#'
#' @param lc spatraster with numeric landcover classifications
#' @param type character classification scheme (currently supported is LC_Type 1
#' (annual international geosphere-biosphere programme classification) from
#' MODIS MCD12Q1.)
#'
#' @description
#' type options: modis_igbp (IGBP classification from MODIS12Q1)
#'
#'
reclass_landcover = function(lc, type = "modis_igbp") {
  # currently shows all grassland as short grassland
  if(type == "modis_igbp") {
    rcl = matrix(c(11,12,
                   12,13,
                   13,14,
                   14,15,
                   15,NA,
                   17,NA,
                   255,NA),
                 byrow = T, ncol = 2)
  }

  reclass = classify(lc, rcl, others = NULL)
  return(reclass)
}

