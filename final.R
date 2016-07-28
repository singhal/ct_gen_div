library(AICcmodavg)
library(ape)
library(geosphere)
library(maptools)
library(nlme)
library(phylolm)
library(phytools)
library(pcaMethods)
library(rgeos)

##########################################
# variables
out_dir = '/Users/sonal/Desktop/activeWork/genetic_diversity_ctenotus/'
pi_val = 'median_pi'
genus = 1

###########################################
### functions #############################
###########################################


## from https://github.com/mrhelmus/phylogeny_manipulation/blob/master/AICc.r
AICc.phylolm<-function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
    LL <- logLik(mod)$logLik
    K <- logLik(mod)$df  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    return(AICc)
  }


get_tree <- function(genus, cl) {
	# function to get a dated tree
	tree = paste('/Users/sonal/macroevolution/eco_IBD_oz/data/species_tree/astrid_brlens/', genus, '_best_summary.tre', sep="")
	a = read.tree(tree)
	a = drop.tip(a, to_drop)
	
	a$tip.label = as.vector(cl[match(a$tip.label, row.names(cl)), 'LatinName'])
	
	# define the tree root; this is defined from Huateng's tree
	if (genus == 'Ct') {
		root = c("C. catenifer", "C. youngsoni", "C. impar")
		}
	else {
		root = c("L. baynesi", "L. picturata", "L. edwardsae", "L. punctatovittata", "L. emmotti", "L. apoda")
		}
	a = root(a, root, resolve.root=TRUE)
	
	# date the tree using a modified version of Sanderson's approach
	# the lambda was set by looking at a range of lambdas - this has the best likelihood
	date = chronos(a, lambda=1e-5)
	# need to write the tree and then re-read it to change the class type
	write.tree(date, file="~/Desktop/tmp.tre")
	date = read.newick("~/Desktop/tmp.tre")
	date$tip.label = gsub("_", " ", date$tip.label)
	file.remove("~/Desktop/tmp.tre")
	return(date)
}

get_correlations <- function(a, iv, name) {
	for (y in 1:(length(iv) - 1)) {
		for (z in (y + 1):length(iv)) {
			tmp = a[, c(iv[y], iv[z])]
			tmp = tmp[complete.cases(tmp),]
			
			if (nrow(tmp) > 4) {
				r = cor.test(tmp[,1], tmp[,2], method="spearman")
				if (r$p.value < 0.05) {
					pval = "***"
					}
				else {
					pval = ''
				}
				rho = round(r$estimate ^ 2, 3)
				if (rho > 0.7) {
					cat(name, ": ", iv[y], ", ", iv[z], " ---- ", rho, pval, "\n", sep="")
				}
			}	
		}
	}
}

fit_models <- function(x, pi_val, variables, models, tree) {
	# fit all additive models and store results
	fits = vector("list", length=length(models))
	for (i in 1:length(models)) {
		tmp1 = x[complete.cases(x[, models[[i]]]), ]
		# rownames(tmp1) = tmp1$LatinName
	
		drop = setdiff(tree$tip.label, rownames(tmp1))
		tmp_tree = drop.tip(tree, drop)
		tmp1 = tmp1[tmp_tree$tip.label, ]
		
		fmla <- as.formula(paste(pi_val, " ~ ", paste(models[[i]], collapse= "+")))
		fit = phylolm(fmla, data=tmp1, tmp_tree, model="lambda", lower.bound=0)
	
		fits[[i]] = fit
		if (i %% 100 == 0) {
			cat("Model ", i, " of ", length(models), " done.\n", sep="")
		}
	}

	aics = sapply(fits, AICc.phylolm)
	best = min(aics)
	raw_weights = sapply(aics, weights <- function(x) {exp((best - x)/2)})
	weights = raw_weights / sum(raw_weights)

	results = vector("list", length=length(variables))
	names(results) = variables
	
	for (i in 1:length(variables)) {
		var_models = sapply(models, y <- function(x) {variables[i] %in% x})
	
		var_fits = fits[var_models]
		var_coef = sapply(var_fits, y <- function(x) {x$coefficients[variables[i]]})
		var_pval = sapply(var_fits, y <- function(x) {summary(x)$coefficients[variables[i],'p.value']})
		var_wt = weights[var_models]
		var_raw_wt = raw_weights[var_models]
		var_rel_raw_wt = var_raw_wt / sum(var_raw_wt)
	
		coef = sum(var_coef * var_rel_raw_wt)
		pval = sum(var_pval * var_rel_raw_wt)
		rel_imp = sum(var_wt)
	
		vals = c(coef, pval, rel_imp)
		names(vals) = c("coefficient", "pvalue", "rel_importance")
		results[[i]] = vals
	}

	full_results = list(fits, aics, weights, results)
	names(full_results) = c("fits", "AIC", "weight", "results")
	return(full_results)
}

###########################################
### prep the data #########################
###########################################

### get cluster data set
cl = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised2.csv", stringsAsFactors=F, na.string="")
cl = cl[!is.na(cl$GMYC_RAxML2),]
cl = unique(cl[,c("GMYC_RAxML2", "LatinName")])
rownames(cl) = cl$GMYC_RAxML2
cl$GMYC_RAxML2 = NULL

### get and write the trees
ct_tree = get_tree("Ct", cl)
le_tree = get_tree("Le", cl)
tree = list(ct_tree, le_tree)

### median nuc pi (based on inds)
pi = read.csv("~/macroevolution/eco_IBD_oz/data/pop_gen/all_clusters.individual_pi.csv", stringsAsFactors=F, na.string="NA")
# get rid of low quality data
pi = pi[pi$denom > 10000,]
# get counts
counts = table(pi$cluster)
# get median pi for a cluster
d = aggregate(pi[,c("pi")], by=list(pi$cluster), median, na.rm=T)
names(d) = c("GMYC_RAxML2", "nuc_pi")
cl$median_pi = d[match(row.names(cl), d$GMYC_RAxML2), "nuc_pi"]

### number of inds
cl$ninds = counts[row.names(cl)]

### mt pi
pi2 = read.csv("~/macroevolution/eco_IBD_oz/data/pop_gen/species_diversity.csv", stringsAsFactors=F, na.string="nan")
cl$mt_pi = pi2[match(row.names(cl), pi2$cluster), "mt_pi"]

### species-wide nuc pi
cl$species_pi = pi2[match(row.names(cl), pi2$cluster), "pi"]

### biomes
biome = read.csv("~/macroevolution/eco_IBD_oz/data/ecoregions/biomes.csv", stringsAsFactors=F)
cl$biome = biome[match(cl$LatinName, biome$species), "biome"]

### genus
cl$genus = rep("Ct", nrow(cl))
cl[grep("Le", rownames(cl)), "genus"] = "Le"

### range size
### &
### longitudinial midpoint
rangedir = '/Users/sonal/macroevolution/eco_IBD_oz/data/geography/new_ranges/'
cl$range_size = rep(NA, nrow(cl))
cl$lat_midpoint = rep(NA, nrow(cl))
for (i in 1:nrow(cl)) {
	sp_name = gsub('\\.? ', '_', cl[i, "LatinName"])
	shp = paste(rangedir, sp_name, '.shp', sep='')
	range = readShapePoly(fn = shp, proj4string=CRS('+proj=longlat +datum=WGS84'))
	
	# need to sum because some ranges consist of multiple ranges
	# divide by 1e6 to convert from m^2 to km^2
	cl[i, "range_size"] = sum(areaPolygon(range)) / 1e6
	
	# get the latitude of the range
	cl[i, "lat_midpoint"] = gCentroid(range)@coords[2]
	}	

### habitat heterogeneity
d = read.csv("~/macroevolution/eco_IBD_oz/data/heterogeneity/habitat_heterogeneity_species.csv", stringsAsFactors=F, row.names=1)
# drop first column because it is range size
cl = cbind(cl, d[cl$LatinName, 2:ncol(d)])

### habitat stability
d = read.csv("~/macroevolution/eco_IBD_oz/data/stability/species_suitability.csv", stringsAsFactors=F, row.names=1)
cl = cbind(cl, d[cl$LatinName, 1:ncol(d)])

### occurrence raw numbers
### &
### occurrence density
d = read.csv("~/macroevolution/eco_IBD_oz/data/occurrence_data/species_counts.csv", stringsAsFactors=F, row.names=1)
names(d)[names(d) == 'all_ninds'] = "all_noccs"
names(d)[names(d) == 'mus_ninds'] = "mus_noccs"
cl = cbind(cl, d[cl$LatinName, c("all_noccs", "mus_noccs", "all_nlocs", "all_pop_density")])

### time in tree
d = read.csv("~/macroevolution/eco_IBD_oz/data/history/time_in_tree.csv", stringsAsFactors=F, row.names=2)
cl$time_in_tree = d[cl$LatinName, "length"]

### morphology
## PCs
m = read.csv('~/macroevolution/eco_IBD_oz/data/morphology/skink_indmeans_june7.csv', stringsAsFactors=F, na.string="NA")
# change the names so it plays nice with my files
m$treename = gsub("ctenotus_", "C. ", m$treename)
m$treename = gsub("lerista_", "L. ", m$treename)

# get means by species
means = aggregate(m[1:8], by=list(m$treename), mean)
names(means)[1] = "LatinName"
morpho = data.frame(character(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0))
names(morpho) = names(means)

# this is janky code to create a summary data set
sps = rownames(d)
for (i in 1:length(sps)) {
	sp = gsub(" \\d", "", sps[i])
	# we spell one species differently
	if (sp == "C. essingtoni") {
		sp = "C. essingtonae"
	}
	tmp = means[grep(gsub(" \\d", "", sp), means$LatinName, ),]
	morpho = rbind(morpho, tmp[1,])
	rownames(morpho)[i] = sps[i]
}

# add toe data
toe = read.csv("~/macroevolution/eco_IBD_oz/data/morphology/Lerista_toes.csv", stringsAsFactors=F)
morpho = cbind(morpho, toe[match(gsub('\\s\\d', "", rownames(morpho)), toe$species), 2:3])
morpho$LatinName = NULL

# add in missing data
morpho["L. haroldi", "svl"] = 40
morpho["L. kendricki", "svl"] = 67
morpho$genus = rep("Ct", nrow(morpho))
morpho[grep("L\\.", row.names(morpho)), "genus"] = "Le"
a = split(morpho, morpho$genus)

# transform the data before doing a PCA
morph_var = list(c("svl", "head_depth", "d_tymp_rost", "head_width", "antebrachium", "shank", "toe", "depth_pect"), c("svl", "head_depth", "d_tymp_rost", "head_width", "antebrachium", "shank", "toe", "depth_pect", "hand_digits", "toe_digits"))
trans = list(c("svl", "head_depth", "d_tymp_rost", "head_width", "antebrachium", "shank", "toe", "depth_pect"), c("svl", "head_depth", "d_tymp_rost", "head_width", "depth_pect"))

for (i in 1:length(a)) {
	m = a[[i]][, morph_var[[i]]]
	
	for (j in 1:length(trans[[i]])) {
		m[, trans[[i]][j]] = log(m[, trans[[i]][j]])	
	}
	
	if (i == 1) {
		x = prcomp(na.omit(m), scale=T, vector=T)
		a[[i]] = merge(a[[i]], x$x[, c("PC1", "PC2", "PC3")], by="row.names", all.x=T)
		rownames(a[[i]]) = a[[i]]$Row.names
		a[[i]]$Row.names = NULL
	} else {
		x = pca(m, scale="vector", nPcs=7, method="ppca")
		a[[i]] = cbind(a[[i]], x@scores[, c("PC1", "PC2", "PC3")])
	}
}

morpho2 = do.call("rbind", a)
morpho2$genus = NULL
rownames(morpho2) = gsub("Ct.", "", rownames(morpho2))
rownames(morpho2) = gsub("Le.", "", rownames(morpho2))
names(morpho2)[names(morpho2) == "PC1"] = "morph_PC1"
names(morpho2)[names(morpho2) == "PC2"] = "morph_PC2"
names(morpho2)[names(morpho2) == "PC3"] = "morph_PC3"

cl = cbind(cl, morpho2[cl$LatinName, 1:ncol(morpho2)])
saveRDS(cl, paste(out_dir, "Ctenotus_Lerista.pi_iv_raw_data.rds", sep=""))

### separate out the two genera
g = split(cl, cl$genus)

# identify correlations for IV within a genera
# the number 8 is hard-coded here - beware if things change in the future
iv = names(cl)[8:length(names(cl))]
to_drop = c("elev_sd", "bio1_range", "bio1_sd", "bio12_range", "bio12_sd", "silica_range", "silica_sd", "aridity_range", "aridity_sd", "vegetation_range", "vegetation_sd", "svl", "head_depth", "d_tymp_rost", "head_width", "antebrachium", "shank", "toe", "depth_pect", "hand_digits", "toe_digits", "max_suitability", "pop_density", "PC1_sd", "PC2_sd", "PC3_sd", "veg_zones", "avg_range", "min_range", "cur_range", "all_nlocs", "morph_PC3")
iv = setdiff(iv, to_drop)

get_correlations(g[[1]], iv, "Ct")

to_drop2 = c("PC2_range", "PC3_range", "all_noccs")
iv = setdiff(iv, to_drop2)

get_correlations(g[[1]], iv, "Ct")

pdf(paste(out_dir, "test_normality/testNormality_", pi_val, ".pdf", sep=""))
# test for what should be logged
test = setdiff(iv, c("lat_midpoint"))
for (i in 1:length(test)) {
	fit1 = lm(g[[genus]][, pi_val] ~ g[[genus]][, test[i]])
	fit2 = lm(g[[genus]][, pi_val] ~ log(g[[genus]][, test[i]]))
	
	par(mfrow=c(2,2))
	plot(fitted(fit1), resid(fit1), main=paste(test[i], ", untrans"))
	abline(h=0, lty=2, col="red")
	lines(lowess(fitted(fit1), resid(fit1)), col="blue")
	plot(fitted(fit2), resid(fit2), main=paste(test[i], ", trans"))
	abline(h=0, lty=2, col="red")
	lines(lowess(fitted(fit2), resid(fit2)), col="blue")
	qqnorm(resid(fit1))
	qqline(resid(fit1), col="red")
	qqnorm(resid(fit2))
	qqline(resid(fit2), col="red")
	
	aic1 = AICc(fit1)
	aic2 = AICc(fit2)
	
	if (aic2 < aic1) {
		diff = exp((aic2 - aic1) / 2)
		cat(test[i], ": ", diff, "\n")
	}
}
dev.off()
g[[genus]]$log_mus_noccs = log(g[[genus]]$mus_noccs)
g[[genus]]$log_range = log(g[[genus]]$range_size)

# not sure if I should keep PC1_sd, PC2_sd etc or PC1_range, PC2_range
iv_ct = c('log_range', 'lat_midpoint', 'elev_range', 'PC1_range', 'avg_suitability', 'log_mus_noccs', 'time_in_tree', 'morph_PC1', 'morph_PC2')
# iv_ct = c("range_size", "lat_midpoint", "elev_range", "PC1_range", "avg_suitability", "mus_noccs", "time_in_tree", "morph_PC1", "morph_PC2")

# create all possible models
# this is so many models
models = list()
for (i in 1:length(iv_ct)) {
	l = combn(iv_ct, i, simplify=FALSE)
	models = append(models, l)
}

# only use cases fit across all models
x = g[[genus]][, c(pi_val, iv_ct)]
rownames(x) = g[[genus]]$LatinName
x = x[complete.cases(x),]
all_fit = fit_models(x, pi_val, iv_ct, models, tree[[genus]])
full_fit = all_fit[["results"]]
saveRDS(all_fit, paste(out_dir, "Ctenotus.full_fit_", pi_val, ".rds", sep=""))

# nsamp = 100
# per_sample = 0.8

# cv = vector("list", length=nsamp)
# for (n in 1:nsamp) {
#	sub = x[sample(nrow(x), round(nrow(x) * per_sample)),]
#	sub_fit = fit_models(sub, iv_ct, models, ct_tree)
#	cv[[n]] = t(as.data.frame(sub_fit))
#	cat(n, "\n")
# }

# cv_sum = as.data.frame(do.call(rbind, cv))
# cv_sum$factor = rownames(cv_sum)
# rownames(cv_sum) = NULL
