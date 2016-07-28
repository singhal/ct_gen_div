library(raster)
library(rgdal)
data(wrld_simpl)

# downloaded from here: http://www.worldclim.org/paleo-climate
aus = wrld_simpl[wrld_simpl$NAME == "Australia",]

dir = list( '/Users/sonal/Desktop/paleomodel/LGM/cclgmbi_2-5m/',
			'/Users/sonal/Desktop/paleomodel/LGM/mrlgmbi_2-5m/',
			'/Users/sonal/Desktop/paleomodel/midHolocene/hemidbi_2-5m/',
			'/Users/sonal/Desktop/paleomodel/midHolocene/hgmidbi_2-5m/')

# get all the rasters
files = lapply(dir, list.files, full.names=T, pattern=".tif")
names(files) = dir

for (i in 1:length(files)) {
	rasters = lapply(files[[i]], raster)
	rasters = lapply(rasters, crop, aus)
	for (j in 1:length(rasters)) {
		outfile = gsub("paleomodel", "paleo_aus", files[[i]][j])
		writeRaster(rasters[[j]], outfile)
	}
}

dir = '/Users/sonal/Desktop/paleomodel/LIG'
files = list.files(dir, full.names=T, pattern=".bil")
for (i in 1:length(files)) {
	r = raster(files[i])
	
	# crop to extent
	r = crop(r, aus)
	
	# aggregate by factor 5
	r = aggregate(r, fact=5, fun=mean)
	
	# write raster as tif format for equivalency
	outfile = gsub("paleomodel", "paleo_aus", files[i])
	writeRaster(r, outfile, format="GTiff")
	}