#' Download GEDI tiles
#'
#'@description
#'Download cleaned tiles of NASA GEDI L2b data. Tiles are indexed by the lower left corner of the quadrant.
#'Each tile is 1x1 degree. Cleaning of GEDI data was done following methods in
#'Burns et al. 2024. Multi-resolution gredded maps of vegetation structure from GEDI
#'https://doi.org/10.1038/s41597-024-03668-4
#'
#' @param e spatraster, spatvector, or spatextent
#' @param fpath file path for downloaded GEDI tiles
#' @param unzip logical; should the downloaded file be unzipped warning if true, file sizes will be large
#' @param download logical (default TRUE) if TRUE will download gedi tiles, if false will only return file paths
#'
#' @return vector names of filepaths and downloads gedi tiles if download == TRUE
#' @export
#'
#'
get_gedi_tiles = function(e, fpath, unzip = F, download = T) {
  # download GEDI tiles from Patrick Burns' prefiltered tile footprint database
  # see filters in Burns et al. 2024
  # tiles are indexed by lower left corner of the quadrant
  # each tile is 1 by 1 degree

  if(methods::is(e, "SpatVector") | methods::is(e, "SpatRaster")) {
    e = terra::ext(e)
  }

  xmin = floor(e[1])
  ymin = floor(e[3])
  xmax = floor(e[2])
  ymax = floor(e[4])

  xs = seq(xmin, xmax, 1)
  ys = seq(ymin, ymax, 1)

  for(x in xs[1:length(xs)]) {
    for(y in ys[1:length(ys)]) {
      url = "https://rcdata.nau.edu/geode_data/GEDIv002_L0204A_20190417to20230316_proc202312/tables/gediv002_l2l4a_va_g"
      if(x<0) {
        lon = paste0(numform::f_pad_zero(abs(x),3),"w")
      } else {
        lon = paste0(numform::f_pad_zero(x, 3), "e")
      }
      if(y<0) {
        lat = paste0(numform::f_pad_zero(abs(y), 2), "s")
      } else {
        lat = paste0(numform::f_pad_zero(y, 2), "n")
      }
      url = paste0(url,lon,lat,".zip")

      fname = strsplit(url, split = "/")
      fname = fname[[1]][length(fname[[1]])]
      fname = paste0(fpath, fname)

      fname.csv = gsub("zip", "csv", fname)
      if(download & (!(file.exists(fname) | file.exists(fname.csv)))) {
        dl = try(utils::download.file(url = url, destfile = fname))
        if(dl != 0) {
          warning(paste0("Download of ", url, " failed. Tile might be in the ocean"))
        }
      }

      if(!file.exists(fname.csv) & download) {
        if(unzip & dl == 0) {
          unzip(fname, exdir = fpath)
          file.remove(fname)
        }
      }
    }
  }
  return(fname.csv)
}


#' Download ERA5 for area
#'
#'@description
#'Download ERA5 climate data using mcera5. Will not download if file is already present.
#'cds access token must be set prior to running this function. See mcera5 vignette for details
#'
#' @param st_time start time of download, character date in format YYYY-MM-DD (e.g., "2023-01-01")
#' @param en_time end time of download, character date in format YYYY-MM-DD (e.g., "2023-12-31)
#' @param uid uid for climate data store (should be your email for your account)
#' @param file_prefix prefix for .nc file
#' @param file_path path for .nc file
#' @param e spat extent to define spatial extent of the download
#' @param overwrite; logical (default FALSE), should files be overwritten
#'
#' @return NA downloads and saves data
#' @export

get_era5 = function(st_time, en_time, uid, file_prefix, file_path, e, overwrite = F) {
  # st_time: vector of start times (e.g., 1st of each month)
  # en_time: vector of end times (e.g. last day/hour of each month)
  # uid: uid for climate data store (should be your email for your account)
  # file_prefix: prefix for .nc file
  # file_path: path for .nc file
  # e spatextent to define spatial extent of the download

  # check if file exists and if not download era5 data

  req = mcera5::build_era5_request(xmin = e[1], xmax = e[2],
                                   ymin = e[3], ymax = e[4],
                                   start_time = as.POSIXlt(st_time, tz = "UTC"),
                                   end_time = as.POSIXlt(en_time, tz = "UTC"),
                                   by_month = T, outfile_name = file_prefix)

  #wd = getwd()
  #setwd(file_path)
  mcera5::request_era5(req, uid, out_path = file_path, overwrite = overwrite, combine = T)
  #setwd(wd)

}


#' get_soil
#' Download sand, silt, and clay content at 0-5 cm depth from soilgrids database
#' @param e spatvector giving extent over which to get soil data
#' @param outdir filename to save soil data (ending in .tif)
#' @param overwrite whether to overwrite file is it exists (TRUE) or not (default FALSE)
#'
#' @returns saves spatraster to outdir and returns spatraster
#'
#' @export
get_soil = function(e, outdir, overwrite = FALSE) {
  soil = soilDB::fetchSoilGrids(
    x = e,
    depth_intervals = c("0-5"),
    variables = c("clay", "sand", "silt"),
    grid = T,
    filename = outdir,
    summary_type = "mean",
    overwrite = F
  )
}
