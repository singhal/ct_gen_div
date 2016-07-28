library(ape)
library(phangorn)

d = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv")
d = d[,c("GMYC_RAxML2", "LatinName")]
d = unique(d)

genus = 'Le'

#### 


tree = paste('/Users/sonal/Desktop/trees/astrid_brlens/', genus, '_best_summary.tre', sep="")
a = read.tree(tree)

a$tip.label = as.vector(d[match(a$tip.label, d$GMYC_RAxML2),]$LatinName)
a$node.label = as.numeric(a$node.label)

if (genus == 'Ct') {
	root = c("C. atlas 2", "C. regius 2", "C. catenifer", "C. youngsoni", "C. impar")
	}
else {
	root = c("L. baynesi", "L. picturata", "L. edwardsae", "L. punctatovittata", "L. emmotti", "L. apoda")
}

a = root(a, root, resolve.root=TRUE)

lambdas = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10)
nrep = 5

loglk = matrix(rep(NA, length(lambdas) * nrep), nrow=nrep)
for (i in 1:length(lambdas)) {
	for (j in 1:nrep) {
		date = chronos(a, lambda=lambdas[i], quiet=TRUE)
		loglk[j, i] = attr(date, "ploglik")
		cat(lambdas[i], j, "\n")			
	}
}

date = chronos(a, lambda=1e-5)

plot(date, cex=0.5)
bs = as.numeric(date$node.label)
num_inds = length(date$tip.label)
for (i in 1:date$Nnode) {
	if (!is.na(bs[i])) {
		if (bs[i] < 0.95) {
			nodelabels(bs[i], i + num_inds, frame="rect", bg="white", cex=0.5)
		}
	}	
}