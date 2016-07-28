require(maptools)
require(sp)
require(rgeos)
require(geosphere)
require(plyr)
require(spatstat)
require(RColorBrewer)
data(wrld_simpl)

get_points <- function(genus) {
	d = read.csv(paste("/Users/sonal/macroevolution/eco_IBD_oz/data/occurrence_data/", genus, "_ALA_4Nov.csv", sep=""), na.string="")
	d$sp = gsub(paste(genus, " ", sep=""), "", d$Matched.Scientific.Name)
	d$sp = gsub(" \\S+$", "", d$sp)
	x = split(d, d$sp)
	return(x)
}

point_in_range <- function(row, r) {
	d = data.frame(x=row[1], y=row[2])
	pt = SpatialPoints(d, proj4string=CRS('+proj=longlat +datum=WGS84'))
	within = gIntersection(pt, r)
	if (class(within) == "NULL") {
		return(FALSE)
	} else {
		return(TRUE)
	}
}

mapdir = '/Users/sonal/macroevolution/eco_IBD_oz/data/geography/new_ranges/'
range = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/geography/ranges.csv", stringsAsFactors=F, na.string="")
spnames = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/metadata/species_names.csv", stringsAsFactors=F, na.string="")
rownames(spnames) = spnames$shortname

res = data.frame(all_orig=rep(NA, nrow(range)), all_ninds=rep(NA, nrow(range)), 
                 all_nlocs=rep(NA, nrow(range)), mus_orig=rep(NA, nrow(range)), 
                 mus_ninds=rep(NA, nrow(range)), mus_nlocs=rep(NA, nrow(range)))
rownames(res) = range$LatinName

pts = list(list(), list())
pts[[1]] = get_points("Ctenotus")
pts[[2]] = get_points("Lerista")
names(pts) = c("Ctenotus", "Lerista")

for (i in 1:nrow(range)) {
	latinname = range[i, "LatinName"]
	r = paste(mapdir, gsub("\\.? ", "_", latinname), ".shp", sep="")
	r = readShapePoly(r, proj4string=CRS('+proj=longlat +datum=WGS84'))
	r = gBuffer(r, width=1)
	
	sp_names = strsplit(range[i, "WithWhat"], ",")[[1]]
	
	all_sp_pts = data.frame(character(0), character(0), character(0))
	names(all_sp_pts) = c("lon", "lat", "type")

	for (j in 1:length(sp_names)) {
		short_name = sp_names[j]
		if (short_name != "NA") {
			genus = spnames[short_name, "genus"]
			sp = spnames[short_name, "sp"]
			
			allpts = pts[[genus]][[sp]][, c("Longitude...processed", "Latitude...processed", "Basis.Of.Record...processed")]
			names(allpts) = c("lon", "lat", "type")
			allpts = allpts[complete.cases(allpts),]
			all_sp_pts = rbind(all_sp_pts, allpts)
		}
	}
	if (nrow(all_sp_pts) > 0) {
		orig_pts = all_sp_pts
		sp_pts = ddply(all_sp_pts, .(all_sp_pts$lon, all_sp_pts$lat, all_sp_pts$type), nrow)
		names(sp_pts) = c("lon", "lat", "type", "count")
		sp_pts$lon = round(sp_pts$lon, 4)
		sp_pts$lat = round(sp_pts$lat, 4)
		
		sp_pts2 = SpatialPoints(sp_pts[, c("lon", "lat")], proj4string=CRS('+proj=longlat +datum=WGS84'))
		sp_pts3 = gIntersection(sp_pts2, r)
		
		res[latinname, "all_orig"] = sum(sp_pts$count)
		res[latinname, "mus_orig"] = sum(sp_pts[sp_pts$type == 'PreservedSpecimen',]$count)
		
		sp_pts3 = data.frame(sp_pts3)
		if (nrow(sp_pts3) > 0) {
			names(sp_pts3) = c("lon", "lat")
		
			sp_pts3 = merge(sp_pts3, sp_pts, by=c("lon", "lat"), all.x=TRUE)
		
			res[latinname, "all_nlocs"] = nrow(sp_pts3)
			res[latinname, "all_ninds"] = sum(sp_pts3$count)
			
			res[latinname, "mus_nlocs"] = nrow(sp_pts3[sp_pts3$type == 'PreservedSpecimen',])
			res[latinname, "mus_ninds"] = sum(sp_pts3[sp_pts3$type == 'PreservedSpecimen',]$count)
			}
		}
	cat(latinname, "\n")
}

### get cluster data set
cl = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv", stringsAsFactors=F, na.string="")
cl = cl[!is.na(cl$GMYC_RAxML2),]
cl = unique(cl[,c("GMYC_RAxML2", "LatinName")])

### range sizes
rangedir = '/Users/sonal/macroevolution/eco_IBD_oz/data/geography/new_ranges/'
res$range_size = rep(NA, dim(res)[1])
for (i in 1:nrow(res)) {
	sp_name = gsub('\\.? ', '_', rownames(res)[i])
	shp = paste(rangedir, sp_name, '.shp', sep='')
	range = readShapePoly(fn = shp, proj4string=CRS('+proj=longlat +datum=WGS84'))
	
	# need to sum because some ranges consist of multiple ranges
	# divide by 1e6 to convert from m^2 to km^2
	res[i, "range_size"] = sum(areaPolygon(range)) / 1e6
	}	

# add pop density
res$all_pop_density = res$all_ninds / res$range_size
write.csv(res, "/Users/sonal/macroevolution/eco_IBD_oz/data/occurrence_data/species_counts.csv")

########################################
# individual based analysis  ###########
########################################

pi = read.csv("~/macroevolution/eco_IBD_oz/data/pop_gen/all_clusters.individual_pi.csv", stringsAsFactors=F, na.string="NA")
# get rid of low quality dat
pi = pi[pi$denom > 10000,]

# add lat long
ll = read.csv("~/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv", stringsAsFactors=F, na.string="NA")
pi$lon = ll[match(pi$ind, ll$sample_id), "lon"]
pi$lat = ll[match(pi$ind, ll$sample_id), "lat"]
pi$LatinName = cl[match(pi$cluster, cl$GMYC_RAxML2), "LatinName"]
pi = pi[complete.cases(pi$lon),]

range = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/geography/ranges.csv", stringsAsFactors=F, na.string="")

pi2 = split(pi, pi$LatinName)
# second result dir
res2 = data.frame(	rep(NA, length(pi2)), 
					rep(NA, length(pi2)))
names(res2) = c("genus", "density_cv")
rownames(res2) = names(pi2)

for (i in 1:length(pi2)) {
	pi2[[i]]$density = rep(NA, nrow(pi2[[i]]))
	latinname = names(pi2)[i]
	r = paste(mapdir, gsub("\\.? ", "_", latinname), ".shp", sep="")
	r = readShapePoly(r, proj4string=CRS('+proj=longlat +datum=WGS84'))
	r = gBuffer(r, width=0)
	
	sp_names = strsplit(range[range$LatinName == latinname, "WithWhat"], ",")[[1]]
	sp_pts = data.frame(character(0), character(0))
	names(sp_pts) = c("lon", "lat")
	for (j in 1:length(sp_names)) {
		short_name = sp_names[j]
		if (short_name != "NA") {
			genus = spnames[short_name, "genus"]
			sp = spnames[short_name, "sp"]
			
			tmppts = pts[[genus]][[sp]][, c("Longitude...processed", "Latitude...processed")]
			names(tmppts) = c("lon", "lat")
			tmppts = tmppts[complete.cases(tmppts),]
			sp_pts = rbind(sp_pts, tmppts)
		}
	}
	if (nrow(sp_pts) > 0) {
		sp_pts2 = ddply(sp_pts, .(sp_pts$lon, sp_pts$lat), nrow)
		lon = c()
		lat = c()
		for (j in 1:nrow(sp_pts2)) {
			lon = c(lon, jitter(rep(sp_pts2[j, 1], sp_pts2[j, 3]), amount=0.0001))
			lat = c(lat, jitter(rep(sp_pts2[j, 2], sp_pts2[j, 3]), amount=0.0001))
		}
	
		sp_pts3 = ppp(lon, lat, range(lon), range(lat))
		dens_rast = raster(density(sp_pts3))
	
		res2[latinname, "density_cv"] = cv(sample(getValues(dens_rast), 100))
	
		samp_pts = SpatialPoints(pi2[[i]][, c("lon", "lat")], proj4string=CRS('+proj=longlat +datum=WGS84'))
		pi2[[i]]$density = extract(dens_rast, samp_pts)
	}
	cat(i, "\n")
}

pi = unsplit(pi2, f=pi$LatinName)
write.csv(pi, "/Users/sonal/macroevolution/eco_IBD_oz/data/occurrence_data/density_at_individual_pts.csv")

res2 = data.frame(ninds = rep(NA, nrow(range)), cor = rep(NA, nrow(range)), pval = rep(NA, nrow(range)))
rownames(res2) = range$LatinName
for (i in 1:length(pi2)) {
	tmp = pi2[[i]]
	tmp = tmp[complete.cases(tmp$density),]
	if (nrow(tmp) > 2) {
		test = cor.test(tmp$pi, tmp$density, method="spearman", alternative="g")
		res2[names(pi2)[i], "cor"] = test$estimate
		res2[names(pi2)[i], "pval"] = test$p.value
	}
	res2[names(pi2)[i], "ninds"] = nrow(tmp)
}

# look at the resutls with a sufficient number of inds
x = res2[res2$ninds > 9,]
# plot correlations colored by signifiance
breaks = hist(x$cor, xlab="correlation", ylab="frequency", main="")$breaks
hist(x$cor, xlab="correlation", ylab="frequency", main="", col="gray", border=NA)
hist(x[x$pval < 0.05,]$cor, breaks=breaks, add=T, col="black", border=NA)
legend("topright", c("sig", "non-sig"), fill=c("black", "gray"), bty="n")