require(alphahull)
require(rgeos)
require(rgdal)
require(maps)
require(mapdata)
require(maptools)
require(rJava)
require(dismo)
require(raster)
require(sp)
data(wrld_simpl)

# location of raster files
cur_rast_dir = '/Users/sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/'
hist_rast_dir = list(	c('/Users/sonal/macroevolution/eco_IBD_oz/data/paleo_aus/LIG/', 'lig_30s_bio_'),
						c('/Users/sonal/macroevolution/eco_IBD_oz/data/paleo_aus/LGM/cclgmbi_2-5m/', 'cclgmbi'),						
						c('/Users/sonal/macroevolution/eco_IBD_oz/data/paleo_aus/midHolocene/hemidbi_2-5m/', 'hemidbi'))
names(hist_rast_dir) = c("LIG", "LGM", "midHol")
# min_pts
min_pts = 3
# threshold for ENM
# anything in the 99% of the distribution 
# for suitability will be treated as part of the range
threshold <- 0.01

# occurence data
occ_dir = '/Users/sonal/macroevolution/eco_IBD_oz/data/geography/niche_points/'
# range_dir
r_dir = '/Users/sonal/macroevolution/eco_IBD_oz/data/geography/new_ranges/'
# out dir
out_dir1 = '/Users/sonal/macroevolution/eco_IBD_oz/data/stability/niches/'
out_dir2 = '/Users/sonal/macroevolution/eco_IBD_oz/data/stability/'

# ind data 
ll = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv", stringsAsFactors=F)
# cluster data
d = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv", na.string="", stringsAsFactors=F)
d = merge(d, ll[, c("sample_id", "lon", "lat")], by.x="sample", by.y="sample_id")
d = d[, c("GMYC_RAxML2", "LatinName", "lon", "lat", "sample")]
d = d[complete.cases(d),]


get_cur_env_data <- function(rast_dir) {
	#load rasters
	setwd(rast_dir)
	env <- stack(list.files(pattern='.tif$'))
	names(env) <- gsub('aus5min_', '', names(env))
	
	# only select bioclim variables
	# don't have historical data for other variables
	env <- env[[paste('bio',1:19,sep='')]]

	#mask so that all rasters have same NA values
	for (i in 2:19) {
		env[[paste('bio',i,sep='')]] <- mask(env[[paste('bio',i,sep='')]], env[['bio1']])
	}
	
	return(env)
}


get_hist_env_data <- function(rast_dir, stub) {
	#load rasters
	setwd(rast_dir)
	env <- stack(list.files(pattern='.tif$'))
	names(env) <- gsub(stub, 'bio', names(env))
	
	# only select bioclim variables
	# don't have historical data for other variables
	env <- env[[paste('bio',1:19,sep='')]]

	#mask so that all rasters have same NA values
	for (i in 2:19) {
		env[[paste('bio',i,sep='')]] <- mask(env[[paste('bio',i,sep='')]], env[['bio1']])
	}
	
	return(env)
}


run_maxent <- function(cur_env, hist_env, occ) {
	occ = read.csv(occ)
	occ = SpatialPointsDataFrame(coords=occ[,c('origLong','origLat')], data=occ, proj4string=CRS('+proj=longlat +datum=WGS84'))
	
	# don't want to run maxent if all the coordinates
	# are in the same grid - no information there
	if (nrow(occ) < 10) {
		# number unique cells
		num_cells = length(unique(cellFromXY(cur_env[[1]], occ)))
		if (num_cells == 1) {
			# will halt things in the next step
			occ = occ[1,]
		}
	}
	
	# only run MaxEnt if you have a certain number of points
	if (nrow(occ) >= min_pts) {
		# figure out what occurrence points are associated with missing data
		x <- extract(cur_env[[c('bio1')]], occ)
		occ <- occ[which(complete.cases(x)),]
		
		# actually run maxent!
		xm <- maxent(x=cur_env, p=occ, args=c("randomseed=true"))

		# thins variables so no overfitting	
		perm <- xm@results
		# gets permutation significance results
		perm <- perm[grep('permutation', rownames(perm)),]
		names(perm) <- gsub('.permutation.importance', '', names(perm))
		# selects only those that are significant above 5
		topVar <- names(which(perm > 5))
		# cannot build a model with a single variable
		# so if there is only one top variable, then add the second
		if (length(topVar) == 1) { 
			topVar <- c(topVar, names(sort(perm, decreasing=TRUE))[2])
		}

		# build the model again with just the topvar
		xm <- maxent(x=cur_env[[topVar]], p=occ, args=c("randomseed=true"))
		# make current prediction
		cur_pred <- predict(xm, cur_env[[topVar]], progress='text')
		
		# make historical prediction
		hist_pred = vector("list", length(hist_env))
		names(hist_pred) = names(hist_env)
		for (i in 1:length(hist_pred)) {
			hist_xm = maxent(x=hist_env[[i]][[topVar]], p=occ, args=c("randomseed=true"))
			hist_pred[[i]] <- predict(hist_xm, hist_env[[i]][[topVar]], progress='text')
		}
		
		#least training presence threshold
		# presenceProbs <- extract(px, occ)
		# thresh <- quantile(presenceProbs,threshold)
		# converts 0 - 1 habitat suitability to binary presence / absence
		# px2 = px
		# px2[which(values(px) < thresh)] <- NA
		# px2[which(values(px) >= thresh)] <- 1
	
		# master <- binaryToPolygon(px, projection = proj4string(occ), buffer=10000, occ, dropNoOcc=TRUE)
		# results = list(master, xm)
		
		results = c(unlist(hist_pred), cur_pred)
		names(results) = c(names(hist_pred), "current")
		return(results)
	} else {
		# returns NA if too few points
		# or if all points in the same grid cell
		return(list(NA, NA, NA, NA))
	}
}

make_plots <- function(out_dir, results, occ_file, samp_pts, r, sp) {
	pdf(paste(out_dir, gsub('\\.? ', '_', sp), ".pdf", sep=""), width=12.5, height=2)
	par(mfrow=c(1, 5), bty="n", mai=c(0.1, 0.1, 0.1, 0.1))
	
	occ = read.csv(occ_file)
	occ = SpatialPointsDataFrame(coords=occ[,c('origLong','origLat')], data=occ, proj4string=CRS('+proj=longlat +datum=WGS84'))
	
	for (i in 1:length(results)) {
		plot(results[[i]], axes=F, legend=F, xlim=c(112, 155), ylim=c(-44.2, -9))
		title(names(results)[i])
		plot(r, add=T)	
	}
	plot(wrld_simpl[wrld_simpl$NAME == "Australia",], axes=F, xlim=c(112, 155), ylim=c(-44.2, -9), col=alpha("gray", 0.6), border=F)
	plot(r, border=F, col=alpha("forestgreen", 0.6), add=T)	
	points(occ, pch=3, cex=0.5, col=alpha("black", 0.3))
	points(samp_pts[, c("lon", "lat")], pch=3, cex=0.5, col=alpha("red3", 0.7))
	
	dev.off()
}


write_results <- function(out_dir, results, sp) {
	for (i in 1:length(results)) {
		out = paste(out_dir, gsub('\\.? ', '_', sp), "_", names(results)[i], ".tif", sep="")
		writeRaster(results[[i]], out, format="GTiff")
	}
}


gm_mean = function(x){
	x = x[complete.cases(x)]
	mean = prod(x) ^ (1/length(x))
	return(mean)
}


ind_stability <- function(indres, samp_pts, results) {
	pts = SpatialPoints(coords=samp_pts[,c('lon','lat')], proj4string=CRS('+proj=longlat +datum=WGS84'))
	
	stab = cbind(samp_pts, data.frame(extract(stack(results), pts)))
	
	indres[stab$sample, "cur_suitability"] = stab$current
	indres[stab$sample, "min_suitability"] = apply(stab[,6:9], 1, min, na.rm=T)
	indres[stab$sample, "mean_suitability"] = apply(stab[,6:9], 1, gm_mean)
	
	return(indres)
}


sp_stability <- function(spres, r, results, sp) {
	range_vals = data.frame(extract(stack(results), r)[[1]])
	range_vals = range_vals[complete.cases(range_vals),]
	
	spres[sp, 'max_suitability'] = max(range_vals[,4])
	
	range_vals = apply(range_vals, 2, sum)
	
	spres[sp, 'avg_range'] = gm_mean(range_vals)
	spres[sp, 'min_range'] = min(range_vals)
	spres[sp, 'cur_range'] = range_vals['current']
	
	return(spres)
}

cur_env = get_cur_env_data(cur_rast_dir)
hist_env = vector("list", length(hist_rast_dir)) 
names(hist_env) = names(hist_rast_dir)
for (i in 1:length(hist_rast_dir)) {
	hist_env[[i]] = get_hist_env_data(hist_rast_dir[[i]][1], hist_rast_dir[[i]][2])
}

# to run for
sps = unique(d$LatinName)

# ind results
indres = data.frame(	mean_suitability = rep(NA, nrow(d)),
						min_suitability = rep(NA, nrow(d)),
						cur_suitability = rep(NA, nrow(d)),
						max_suitability = rep(NA, nrow(d)))
rownames(indres) = d$sample

# sp results
spres = data.frame(	avg_suitability = rep(NA, length(sps)),
					max_suitability = rep(NA, length(sps)),
					avg_range = rep(NA, length(sps)),
					min_range = rep(NA, length(sps)),
					cur_range = rep(NA, length(sps)))
					
rownames(spres) = sps

for (i in 1:length(sps)) {
	# get Atlas of Living Oz occurrences
	occ_file = paste(occ_dir, gsub('\\.? ', '_', sps[i]), ".csv", sep="")
	results = run_maxent(cur_env, hist_env, occ_file)
	
	# get genetic sample points
	samp_pts = d[d$LatinName == sps[i],]
	
	# get range
	range_file = paste(r_dir, gsub('\\.? ', '_', sps[i]), ".shp", sep="")
	r = readShapePoly(range_file, proj4string=CRS('+proj=longlat +datum=WGS84'))
	
	# print plots with inferred range and starting points
	make_plots(out_dir2, results, occ_file, samp_pts, r, sps[i])
	
	# change all to same projection
	for (x in 1:3) {
		results[[x]] = projectRaster(results[[x]], results[[4]])
	}
	
	# write results to drive
	write_results(out_dir1, results, sps[i])
	
	# extract values for sampled points and take geometric mean for stability 
	indres = ind_stability(indres, samp_pts, results)
	
	# extract current range from each prediction and calculate something for range-wide stability
	spres = sp_stability(spres, r, results, sps[i])
	
	# average values for sampled points for range wide stability
	spres[sps[i], 'avg_suitability'] = gm_mean(indres[samp_pts$sample, 'mean_suitability'])
	
	# get max suitability assigned to inds
	indres[samp_pts$sample, 'max_suitability'] = spres[sps[i], 'max_suitability']
	
	### these stability values have different maximums for different species
	### don't think this means anything biologically
	### should I get rid of this?
	### could divide by current stability values?
	
	}