# crsLL = CRS("+init=epsg:4326")
crsLL = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

# CRS("+init=epsg:3857") without the nagrids stuff
crsMerc = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +no_defs")

#crsMerc =  CRS("+proj=merc +ellps=sphere +units=m")
#crsLlSphere = CRS("+proj=longlat +ellps=sphere")

crsModis <- CRS(
  "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
)
