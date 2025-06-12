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
#'
#' @return PAI profile at 1 m intervals from ground to top of canopy. Values add up to total PAI in gedi data
.pai_vertprofile = function(pai_z, pai, h) {
  pai_z = pai_z[pai_z != 0] # remove zeros above top of canopy
  # convert to 1m intervals
  # assume even distribution within 5m intervals
  #pai_z2 = spline(pai_z, n = round(h))[[2]]

  # first is pai 0-5m and is equivalent to total pai
  # get pai in each height bin so that sum of pai = total pai
  if(length(pai_z)!=1) {
    pai_z = pai_z - c(pai_z[2:length(pai_z)],0)
  }

  if(length(pai_z != 1)) {
    pai_z1 = pai_z[1:length(pai_z)-1]
    pai_z1 = rep(pai_z1/5, each = 5)
    pai_z2 = pai_z[length(pai_z)]
    x = 5 - (length(pai_z)*5 - (ceiling(h)))
    pai_z2 = rep(pai_z2/x, each = x)
    pai_z = c(pai_z1, pai_z2)
  } else {
    x = floor(h)
    pai_z = rep(pai/x, each = x)
  }


  return(pai_z)
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
#' @param type classification scheme (currently supported is LC_Type 1
#' (annual international geosphere-biosphere programme classification) from
#' MODIS MCD12Q1)
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
  reclass = classify(lc[[1]], rcl, others = NULL)
  return(reclass)
}

