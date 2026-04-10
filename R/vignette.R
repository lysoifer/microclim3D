library(terra)

# Download GEDI data

# Make an AOI in Mark Twain National Forest
pts = data.frame(lon = -91.08, lat = 37.69, region = "mark_twain_nf")
pts = vect(pts, geom = c("lon", "lat"), crs = "epsg:4326")
pts = buffer(pts, 500)


