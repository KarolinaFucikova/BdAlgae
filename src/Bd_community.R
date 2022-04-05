library(vegan)

Bd <- read.csv(file = "data/Bd.csv", row.names =1)
# these data have abundance/concentrations of species per 10ul

# at some point it might need to be cleaned, leave out some taxa or samples
# clean_bd <- subset(Bd, select = -HB1967)

# the commands that compute dissimilarity require the data samples to be rows
# we will transpose 
transp_bd <- t(Bd)

# calculate Shannon index for each site
H <- diversity(transp_bd)
H

# calculate the Jaccard index using
# Jaccard is appropriate for presence/absence data
eucl <- vegdist(transp_bd,method="euclidean")

# now we'll make a hierarchical cluster object, which will be the basis of the dendrogram
bd_hc <-hclust(eucl)
# and finally the dendrogram
plot(bd_hc)

# trying principal coordinates analysis
# probably don't have enough data for it but seems most appropriate according to
# https://www.davidzeleny.net/anadat-r/doku.php/en:pcoa_nmds
# there are a lot of zeroes, so PCA and similar analyses are not ideal

PCoA.bd<-capscale(transp_bd~1,distance="bray")
PCoA.bd
summary(PCoA.bd)

plot(PCoA.bd, type = "n")
points(PCoA.bd, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(PCoA.bd, display = "sites", cex=0.7, col="blue")

# https://websites.pmc.ucsc.edu/~mclapham/Rtips/ordination.htm#:~:text=PCoA%20is%20a%20distance%2Dbased,for%20more%20information%20on%20distances).
