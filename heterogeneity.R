library(sp)
library(maptools)
library(geosphere)
library(raster)
library(corrplot)
library(rgeos)

################################################
# get a measure of environmental heterogeneity 
# based on all the variables
################################################

# australia coast
aus = readShapePoly("/Users/Sonal/macroevolution/eco_IBD_oz/documentation/eco_IBD/pascal_scripts/AusCoast.shp", proj4string=CRS("+proj=longlat +datum=WGS84"))
# sample random points in australia
pts = spsample(aus, 50000, type="regular")

#  all rasters
r_files = list.files('/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/', pattern=".tif", full.names=T)
rasters = stack(r_files)
# get env data at sampled points
vals = extract(rasters, pts)

# combine with randomly sampled points
vals = cbind(vals, data.frame(pts))
# get rid of NA
vals2 = vals[complete.cases(vals),]

# do the pca
valspca = prcomp(vals2[,1:27], scale=T, center=T)

# turn the first pcas (0.35, 0.70, 0.79) into data frame
pca1 = cbind(vals2[, 28:29], valspca$x[,1])
names(pca1) = c("x", "y", "PC1")
pca2 = cbind(vals2[, 28:29], valspca$x[,2])
names(pca2) = c("x", "y", "PC2")
pca3 = cbind(vals2[, 28:29], valspca$x[,3])
names(pca3) = c("x", "y", "PC3")

# turn them into rasters
pca1 = rasterFromXYZ(pca1)
pca2 = rasterFromXYZ(pca2)
pca3 = rasterFromXYZ(pca3)
pcas = stack(pca1, pca2, pca3)

d = read.csv('~/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv', stringsAsFactors=F, na.string="")
rangedir = "~/macroevolution/eco_IBD_oz/data/geography/new_ranges/"

res_names = as.vector(na.omit(unique(d$LatinName)))
res = data.frame(	size=rep(NA, length(res_names)), 
					elev_range=rep(NA, length(res_names)),
					elev_sd=rep(NA, length(res_names)),
					bio1_range=rep(NA, length(res_names)),
					bio1_sd=rep(NA, length(res_names)),
					bio12_range=rep(NA, length(res_names)),
					bio12_sd=rep(NA, length(res_names)),
					silica_range=rep(NA, length(res_names)),
					silica_sd=rep(NA, length(res_names)),
					aridity_range=rep(NA, length(res_names)),
					aridity_sd=rep(NA, length(res_names)),
					vegetation_range=rep(NA, length(res_names)),
					vegetation_sd=rep(NA, length(res_names)),
					PC1_range=rep(NA, length(res_names)),
					PC1_sd=rep(NA, length(res_names)),
					PC2_range=rep(NA, length(res_names)),
					PC2_sd=rep(NA, length(res_names)),
					PC3_range=rep(NA, length(res_names)),
					PC3_sd=rep(NA, length(res_names)),
					veg_zones = rep(NA, length(res_names))
					)
rownames(res) = res_names

# the rasters that win
rasters = list(	'/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/aus5min_alt.tif',
				'/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/aus5min_bio1.tif',
				'/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/aus5min_bio12.tif',
				'/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/aus5min_silica.tif',
				'/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/aus5min_aridityIndex.tif',
				'/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/aus5min_vegetation.tif')
names(rasters) = c("elev", "bio1", "bio12", "silica", "aridity", "vegetation")
rasters = stack(rasters)
# vegetation zones
ogel = raster("/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/glccapl20/apogel20.tif")

# get the values
for (i in 1:nrow(res)) {
	# get range
	file = paste(rangedir, gsub("\\.? ", "_", rownames(res)[i]), ".shp", sep="")
	r = readShapePoly(file, proj4string=CRS("+proj=longlat +datum=WGS84"))
	# get range size
	res[i, "size" ] = sum(areaPolygon(r)) / 1e6
	
	# get environmental values
	vals = data.frame(extract(rasters, r)[[1]])
	for (j in 1:ncol(vals)) {
		val = names(vals)[j]
		
		res[i, paste(val, "_range", sep="")] = max(vals[, j], na.rm=T) - min(vals[, j], na.rm=T)
		res[i, paste(val, "_sd", sep="")] = sd(vals[, j], na.rm=T)
	}
	
	# get pca values
	vals2 = data.frame(extract(pcas, r)[[1]])
	for (j in 1:ncol(vals2)) {
		val = names(vals2)[j]
		
		res[i, paste(val, "_range", sep="")] = max(vals2[, j], na.rm=T) - min(vals2[, j], na.rm=T)
		res[i, paste(val, "_sd", sep="")] = sd(vals2[, j], na.rm=T)
	}
	
	# get vegetation zones
	veg = extract(ogel, r)[[1]]
	veg_counts = table(veg) / length(veg)
	# don't want to include a vegetation zone if it is not much of the zone
	res[i, 'veg_zones'] = length(veg_counts[which(veg_counts > 0.05)])
	
	cat(i, "\n")
}

# take the log values
res2 = res
for (i in 1:13) {
	res2[, i] = log(res2[, i])
}

###########################
# now for the individuals
###########################

ll = read.csv('/Users/sonal/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv', stringsAsFactors=F, na.string="")
ll = ll[ll$sample_id %in% d$sample, c("sample_id", "lon", "lat")]
ll = ll[complete.cases(ll),]
# turn all lat long to a SpatialPoints object
pts = SpatialPoints(ll[, c("lon", "lat")], proj=CRS("+proj=longlat +datum=WGS84"))

# function to calculate range extent
# ... means to ignore passed in value, which is na.rm because of logic of raster::aggregate
range_extent <- function(vals, ...) {
	rextent = max(vals, na.rm=T) - min(vals, na.rm=T)
	return(rextent)
}

# number of zones after aggregrating
num_zones <- function(vals, ...) {
	veg_counts = table(vals) / length(vals)
	# don't want to include a vegetation zone if it is not much of the zone
	counts = length(veg_counts[which(veg_counts > 0.05)])
	return(counts)
}

# get aggregated data
# factors change to make it about the same for all the rastesr
range_rasters = aggregate(rasters, fact=10, fun=range_extent)
sd_rasters = aggregate(rasters, fact=10, fun=sd)

range_pcas = aggregate(pcas, fact=7, fun=range_extent)
sd_pcas = aggregate(pcas, fact=7, fun=sd)

num_zones_ogel = aggregate(ogel, fact=38.48335, fun=num_zones)

# extract the points
# change the names
range1 = data.frame(extract(range_rasters, pts))
names(range1) = unlist(lapply(names(range1), paste, "_range", sep=""))
sd1 = data.frame(extract(sd_rasters, pts))
names(sd1) = unlist(lapply(names(sd1), paste, "_sd", sep=""))
range2 = data.frame(extract(range_pcas, pts))
names(range2) = unlist(lapply(names(range2), paste, "_range", sep=""))
sd2 = data.frame(extract(sd_pcas, pts))
names(sd2) = unlist(lapply(names(sd2), paste, "_sd", sep=""))
zones1 = data.frame(extract(num_zones_ogel, pts))
names(zones1) = c("veg_zones")

# bind the values
res2 = cbind(range1, sd1, range2, sd2, zones1, ll)

write.csv(res2, "/Users/sonal/Desktop/habitat_heterogeneity_inds.csv")