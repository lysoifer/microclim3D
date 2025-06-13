#' soil_params
#' get soil parameters for microclimate model
#'
#' @param soilr spatraster of soildata with three layers named sand, silt, and clay
#' @param e extent to crop to
#'
#' @return spatraster of soiltypes numbered according to microclimf::soilparams
soil_params = function(soilr,e) {

  terra::crop(soilr, e)

  soildf = as.data.frame(soilr, xy = T)
  cl = which(grepl("clay", colnames(soildf)))
  sa = which(grepl("sand", colnames(soildf)))
  sl = which(grepl("silt", colnames(soildf)))

  colnames(soildf)[cl] = "CLAY"
  colnames(soildf)[sa] = "SAND"
  colnames(soildf)[sl] = "SILT"

  # classify soil texture according to USDA classification
  soilclass = TT.points.in.classes(soildf, class.sys = "USDA.TT", text.tol = 0.1)

  colnames(soilclass) = c(
    "Clay", "Silty clay", "Sandy clay",
    "Clay loam", "Silty clay loam", "Sandy clay loam",
    "Loam", "Silty loam", "Sandy loam",
    "Silt", "Loamy sand", "Sand"
  )

  soilclass = apply(soilclass, 1, function(x) {names(which(x>0))})
  soilclass = sapply(soilclass, function(x) ifelse(length(x)>1,sample(x,1),x)) # if on an edge or corner, between soil classes, randomly choose one of them
  soilnum = sapply(soilclass, function(x){which(microclimf::soilparameters$Soil.type == x)}) # what number is soilclass in microclimf soil parameters

  soildf$soilnum = soilnum

  r = rast(soildf, crs = crs(soilr), extent = ext(soilr))

  return(r)

}






# get soil class from soilgrids
# soilclass = terra::extract(soil, gedi_shot[,c("longitude_bin0", "latitude_bin0")])
# colnames(soilclass) = c("ID", "CLAY", "SAND", "SILT")
#
#
#
#
# if(model == "microclimf") {
#   soilrast = rast(xmin = gedi_shot$longitude_bin0 - 0.1,
#                   xmax = gedi_shot$longitude_bin0 + 0.1,
#                   ymin = gedi_shot$latitude_bin0 - 0.1,
#                   ymax = gedi_shot$latitude_bin0 + 0.1,
#                   resolution = 0.2, crs = "epsg:4326",
#                   vals = soilnum)
#   soilc = microclimf::soilcfromtype(soilrast, groundr = 0.15) # need to refine groundr
# }
#
# if(model == "micropoint") {
#   soilp = soilparams[soilparams$Soil.type == soilclass,]
#
#   groundparams = list(
#     gref = habvars$refg, # ground reflectance
#     slope = slp, # slope of ground surface (deg)
#     aspect = asp, # aspect of ground surface
#     em = em, # emissivity of ground surface
#     rho = soilp$rho, # soil bulk density (mg/m3)
#     Vm = soilp$Vm, # volumetric mineral fraction of soil
#     Vq = soilp$Vq, # Volumetric quartz fraction of soil
#     Mc = soilp$Mc, # Mass fraction of clay
#     b = soilp$b, # Shape parameter fo campbell soil moisture
#     Psie = -soilp$psi_e, # Matric potential
#     Smax = soilp$Smax, # Volumetric water content
#     Smin = soilp$Smin, # Residual water content
#     alpha = soilp$alpha, # shape parameter of the van Genuchten model
#     n = soilp$n, # pore size distribution
#     Ksat = soilp$Ksat # saturated hydraulic conductivity
#   )
# }
