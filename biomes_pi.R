library(raster)
library(sp)
library(maptools)
library(rgeos)
library(geosphere)
library(raster)
library(RColorBrewer)

colors = brewer.pal(6, "Set2")
# minimum part of a species range that has to be there 
min = 0.3

#############################
# individual based analysis #
#############################

# lat long data
ll = read.csv('~/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv', stringsAsFactors=F, na.string='')
cl = read.csv('~/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv', stringsAsFactors=F, na.string='')
ll$sp = cl[match(ll$sample_id, cl$sample), "LatinName"]
ll = ll[,c("sample_id", "lon", "lat", "sp")]
ll = na.omit(ll)

# the original biomes downloaded from the TNC
# http://maps.tnc.org/gis_data.html
biomes = readShapePoly("~/macroevolution/eco_IBD_oz/data/ecoregions/terr-ecoregions-aus_combined/terr_ecoregions-aus_combined.shp")
# simplify this so that it takes less time to process
x = gSimplify(biomes, 0.1, topologyPreserve=T)
x$biome = as.character(biomes$WWF_MHTNAM)
# rename the biomes
bnames1 = c("Tropical and Subtropical Moist Broadleaf Forests", "Montane Grasslands and Shrublands", "Temperate Broadleaf and Mixed Forests", "Tropical and Subtropical Grasslands, Savannas and Shrublands", "Deserts and Xeric Shrublands", "Mediterranean Forests, Woodlands and Scrub", "Temperate Grasslands, Savannas and Shrublands")
bnames2 = c("tropical forests", "montane", "temperate forests", "tropical grasslands", "desert", "Mediterranean forests", "temperature grasslands")
for (i in 1:length(bnames1)) {
	 x$biome = gsub(bnames1[i], bnames2[i],  x$biome)
}
b = extract(x, ll[,2:3])
ll = cbind(ll, b[,3])
names(ll)[5] = "biomes"
ll$biomes = as.character(ll$biomes)

# pop gen data
dfile = '/Users/sonal/macroevolution/eco_IBD_oz/data/pop_gen/all_clusters.individual_pi.csv'
d = read.csv(dfile, na.string="NA", stringsAsFactors=FALSE)
d = d[d$denom > 10000,]
ll = merge(ll, d, by.x="sample_id", by.y="ind")
ll = na.omit(ll)

# only compare biomes with enough individuals
n = 100
biomesKeep = names(which(table(ll$biomes) > n))
ll2 = ll[ll$biomes %in% biomesKeep,]
ll2$genus = rep('Ct', dim(ll2)[1])
ll2[grep("Le", ll2$cluster), "genus"] = "Le"
fit = aov(ll2$pi ~ ll2$biomes * ll2$genus)
fit
summary(fit)
boxplot(ll2$pi ~ ll2$biomes + ll2$genus)

# or should it be this? because otherwise, not really independent measures?
ll3 = aggregate(ll2$pi, by=list(ll2$biomes, ll2$sp, ll2$genus), mean)
names(ll3) = c("biome", "sp", "genus", "pi")
fit2 = aov(ll3$pi ~ ll3$biome * ll3$genus)
fit3 = aov(ll3$pi ~ ll3$biomed)

###########################
# Species level analysis  #
###########################

# the original biomes downloaded from the TNC
# http://maps.tnc.org/gis_data.html
biomes = readShapePoly("/Users/sonal/macroevolution/eco_IBD_oz/data/ecoregions/terr-ecoregions-aus_combined/terr_ecoregions-aus_combined.shp")
# simplify this so that it takes less time to process
x = gSimplify(biomes, 0.1, topologyPreserve=T)
x$biome = as.character(biomes$WWF_MHTNAM)

# range dir
rangedir = '/Users/sonal/macroevolution/eco_IBD_oz/data/geography/new_ranges/'

# rename the biomes so that it is easier to deal with
biomes = vector("list", length(x))
names(biomes) = x$biome
bnames1 = c("Tropical and Subtropical Moist Broadleaf Forests", "Montane Grasslands and Shrublands", "Temperate Broadleaf and Mixed Forests", "Tropical and Subtropical Grasslands, Savannas and Shrublands", "Deserts and Xeric Shrublands", "Mediterranean Forests, Woodlands and Scrub", "Temperate Grasslands, Savannas and Shrublands")
bnames2 = c("tropical forests", "montane", "temperate forests", "tropical grasslands", "desert", "Mediterranean forests", "temperature grasslands")
for (i in 1:length(bnames1)) {
	names(biomes) = gsub(bnames1[i], bnames2[i], names(biomes))
}
# helps deal with issues around illegal spatial geometries
for (i in 1:length(x)) {
	a = x[i, ]
	biomes[[i]] = gBuffer(SpatialPolygons(a@polygons,proj4string=a@proj4string), width=0)
}

# lat long data
ll = read.csv('~/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv', stringsAsFactors=F, na.string='')
cl = read.csv('~/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv', stringsAsFactors=F, na.string='')
ll$sp = cl[match(ll$sample_id, cl$sample), "LatinName"]
ll = ll[,c("sample_id", "lon", "lat", "sp")]
ll = na.omit(ll)

# initialize result dataframe
sps = sort(unique(cl$LatinName))
res = data.frame(matrix(NA, nrow=length(sps), ncol=2))
names(res) = c("species", "biome")
res$species = sps

# identify the biomes in which
# more than (min) of the area overlaps with the biome
classify <- function(overlaps) {
	overlaps = overlaps / sum(overlaps)
	hab = overlaps[ overlaps > min ]
	# if nothing is significant overlapping, then call it a generalist
	if (length(hab) < 1) {
		return("generalist")
	} else {
		return(paste(sort(names(hab)), collapse=","))
	}
}

for (i in 1:length(sps)) {
	# get the range
	range = readShapePoly(paste(rangedir, gsub("\\.? ", "_", sps[i]), ".shp", sep=""))
	if (class(range) == 'SpatialPolygonsDataFrame') {
		range = SpatialPolygons(range@polygons,proj4string=range@proj4string)
	}
	# calculate all the overlaps in km2
	overlaps = rep(0, length(biomes))
	names(overlaps) = names(biomes)
	for (x in 1:length(biomes)) {
		overlap = gIntersection(gBuffer(range, width=0), biomes[[x]])
		if (class(overlap) != 'NULL') {
			overlaps[x] = sum(areaPolygon(overlap)) / 1e6
			}
		}	

	res[i, "biome"] = classify(overlaps)
	
	# sampled points
	# sp = ll[ll$sp == sps[i],]
	# pts2 = extract(x, sp[,2:3])
	# res[i, "biome"] = biome1
	cat(i, "\n")
}
write.csv(res, "~/macroevolution/eco_IBD_oz/data/ecoregions/biomes.csv", row.names=F)


# get range sizes
sizes = rep(NA, length(sps))
names(sizes) = sps
for (i in 1:length(sps)) {
	# random points
	range = readShapePoly(paste(rangedir, gsub("\\.? ", "_", sps[i]), ".shp", sep=""))
	sizes[i] = sum(areaPolygon(range)) / 1e6
}

# pop gen data
dfile = '/Users/sonal/macroevolution/eco_IBD_oz/data/pop_gen/all_clusters.individual_pi.csv'
d = read.csv(dfile, na.string="NA", stringsAsFactors=FALSE)
d = merge(cl, d, by.x="sample", by.y="ind")
d$range = sizes[match(d$LatinName, names(sizes))]
d = aggregate(d[,c("pi", "range")], by=list(d$LatinName), median, na.rm=T)
names(d)[1] = "species"
d = merge(d, res)
d$genus = rep(NA, dim(d)[1])
d[grep("C. ", d$species), "genus"] = 'Ct'
d[grep("L. ", d$species), "genus"] = 'Le'

# do anova fits
fit1 = aov(d$pi ~ d$biome + log(d$range) + d$genus + d$genus * log(d$range))
summary(fit1)

# only do anova fits for biome types that are common
counts = table(res$biome)
keep = names(counts[counts > 15])
x = d[d$biome %in% keep,]
boxplot(x$pi ~ x$biome, col=colors, xlab="habitat", ylab=expression(pi))
fit2 = aov(x$pi ~ x$biome + x$genus + log(x$range))
summary(fit2)