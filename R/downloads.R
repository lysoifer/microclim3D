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
