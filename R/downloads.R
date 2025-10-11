#' Download GEDI tiles
#'
#'@description
#'Download cleaned tiles of NASA GEDI L2b data. Tiles are indexed by the lower left corner of the quadrant.
#'Each tile is 1x1 degree. Cleaning of GEDI data was done following methods in
#'Burns et al. 2024. Multi-resolution gredded maps of vegetation structure from GEDI
#'https://doi.org/10.1038/s41597-024-03668-4
#'These points all have leaf on conditions
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
#' This is an older version of the function (pre-microclimdata package)
#' Preference is to use get_soil2
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
    overwrite = overwrite
  )
}


#' get_soil2
#' Download sand, silt, and clay content at 0-5 cm depth from soilgrids database using microclimdata package
#' @param x spatVector, list of spatrasters, or dataframe of coordinates
#' @param coords if providing a dataframe, the column names for x and y coordinates
#' @param crs if providing a dataframe, the crs
#' @param tempdir_location directory to create tempdir folders
#' @param fout directory to story output
#' @param overwrite whether to overwrite file is it exists (TRUE) or not (default FALSE)
#'
#' @returns saves spatraster to outdir and returns spatraster
#'
#' @export

get_soil2 = function(x, coords = c("x", "y"), crs = "epsg:4326",
                     tempdir_location, fout, overwrite = T) {
  # define world grid for downloading data
  # use 5degree file size
  grid = .get_grid(crs = "epsg:4326", tilesize = 5)

  for(i in 1:length(br)) {
    # select the ith location object provided
    if(is(x, "SpatVector")) {
      xi = x[i]
    }
    if(is(x, "data.frame")) {
      xi = x[i,c(coords[1], coords[2])]
      xi = terra::vect(xi, geom = c(coords[1], coords[2]), crs = crs)
    }
    if(is(x, "list")) {
      xi = x[[i]]
      xi = ext(xi)
    }
    c = cells(grid, xi)[,2]
    for(cc in c) {
      # extract cell from the world grid to download
      r = grid[cc, drop = F]
      values(r) = 1

      # name for writing out tiles by extent
      ename = paste(c(ext(r)[1], ext(r)[2], ext(r)[3], ext(r)[4]), collapse = "_")
      fout = paste0(fout, "soilGrids_depthWeightedSoil_", ename, ".tif")

      # make a temporary directory to store downloads - necessary to parallelize so workers don't share tempdirs
      tempdir1= paste0(tempdir_location, i, cc)
      if(!dir.exists(tempdir1)) {dir.create(tempdir_location, recursive = T)}

      if(!file.exists(fout) | overwrite) {
        soil = microclimdata::soildata_download(r, pathdir = tempdir1, deletefiles = T)
        writeRaster(soil, fout, overwrite)
      }
      unlink(tempdir1)
    }
  }
}




#' Download MODIS data
#'
#' @description
#' Downloads MODIS data using modisfast R package. Defaults are to download MODIS landcover
#'
#' @param e list of spatExtents to download modis landcover data (e.g., list(ext(xmin, xmax, ymin, ymax)))
#' @param start character start date (e.g., "2019-01-01")
#' @param end character end date
#' @param usr earthdata username
#' @param pwd earthdata password
#' @param collection collection to download (default MCD12Q1.061) Can be any collection supported by modisfast package
#' @param variables vector of variables to download (must be variables in the collection)
#' @param regions character name of regions for which data will be downloaded (this will be used to name folders for downloaded files)
#' @param outpath file path to save downloads
#'
#' @return dataframe describing downloaded data. See documentation of modisfast::mf_download_data for variable descriptions
#' Data will be downloaded and combined into a single netcds file
#'
#' @export

get_modis_lc = function(e, start, end, usr, pwd, collection = "MCD12Q1.061",
                        variables = c("LC_Type1", "QC"), regions, outpath) {

  # convert roi to sf polygons
  roi = lapply(e, terra::as.polygons, crs = "epsg:4326")
  roi = terra::vect(roi)
  roi$id = regions
  roi = sf::st_as_sf(roi)

  # convert start and end to dates
  time_range = as.Date(c(start, end))

  # add login credentials
  log = modisfast::mf_login(credentials = c(usr, pwd))

  # get urls

  opt_param = modisfast::mf_get_opt_param(collection = collection, roi = roi,
                                          credentials = c(usr, pwd))

  urls = modisfast::mf_get_url(collection = collection,
                    variables = variables,
                    roi = roi,
                    time_range = time_range,
                    single_netcdf = T,
                    opt_param = opt_param)

  # download data
  res_dl = modisfast::mf_download_data(df_to_dl = urls, path = outpath, parallel = T)

  return(res_dl)
}


#'Find gedi orbits
#'@description Define Function to Query CMR
#'The function returns a list of links to download the files of GEDI transects intersecting the bounding box
#'sub-orbit V2 granules directly from the LP DAAC's Data Pool.
#'
#'@param product can be  'GEDI01_B.002', 'GEDI02_A.002', 'GEDI02_B.002'
#'@param bbox area of interest. bbox can be a character "LLlong,LLlat,URlong,URlat" coordinates in WGS84, a spatVector, or a spatRaster
#'
#' @return list of granules with shots in bbox
#'
#' @export
gedi_finder <- function(product, bbox) {

  # Define the base CMR granule search url, including LPDAAC provider name and max page size (2000 is the max allowed)

  cmr <- "https://cmr.earthdata.nasa.gov/search/granules.json?pretty=true&provider=LPCLOUD&page_size=2000&concept_id="

  # Set up list where key is GEDI shortname + version and value is CMR Concept ID
  concept_ids <- list('GEDI01_B.002'='C2142749196-LPCLOUD',
                      'GEDI02_A.002'='C2142771958-LPCLOUD',
                      'GEDI02_B.002'='C2142776747-LPCLOUD')

  if(is(bbox, "SpatVector") | is(bbox, "SpatRaster")) {
    bbox = .bbox_to_char(bbox, "xmin,ymin,xmax,ymax")
  }

  # CMR uses pagination for queries with more features returned than the page size
  page <- 1
  bbox <- sub(' ', '', bbox)  # Remove any white spaces
  granules <- list()          # Set up a list to store and append granule links to

  # Send GET request to CMR granule search endpoint w/ product concept ID, bbox & page number
  cmr_response <- httr::GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page))

  # Verify the request submission was successful
  if (cmr_response$status_code==200){

    # Send GET request to CMR granule search endpoint w/ product concept ID, bbox & page number, format return as a list
    cmr_url <- sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)
    cmr_response <- httr::content(httr::GET(cmr_url))$feed$entry

    # If 2000 features are returned, move to the next page and submit another request, and append to the response
    while(length(cmr_response) %% 2000 == 0){
      page <- page + 1
      cmr_url <- sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)
      cmr_response <- c(cmr_response, httr::content(httr::GET(cmr_url))$feed$entry)
    }

    # CMR returns more info than just the Data Pool links, below use for loop to grab each DP link, and add to list
    for (i in 1:length(cmr_response)) {
      granules[[i]] <- cmr_response[[i]]$links[[1]]$href
    }

    # Return the list of links
    return(granules)
  } else {

    # If the request did not complete successfully, print out the response from CMR
    print(httr::content(httr::GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)))$errors)
  }
}

