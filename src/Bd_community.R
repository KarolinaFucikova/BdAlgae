library(vegan)

Bd <- read.csv(file = "data/Bd_data_clean.csv", row.names =1)
# these data have abundances of species calculated from 3 passes (transects) 
# across each of two slides per sample
# first column has species name, second has group (family/order) affiliation

# taxa can be excluded if needed
# clean_bd <- subset(Bd, select = -speciesname)

# first, combine data from species in the same groups (e.g., merge all Selenastraceae)
Bd_grouped <- aggregate(Bd[2:30], by=list(Group=Bd$Group), FUN=sum)
# turn groups into row names
rownames(Bd_grouped) <- Bd_grouped[,1]
# the first column remains as group and needs to be removed prior to analyses
Bd_grouped[,1] <- NULL

# if we wish to save the cleaned, grouped data as csv:
write.csv(Bd_grouped,"data/Bd_data_grouped.csv", row.names = TRUE)

# the commands that compute dissimilarity require the data samples to be rows
# we will transpose 
transp_bd <- as.data.frame(t(Bd_grouped))

# calculate Shannon index for each sample
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
# seems most appropriate according to
# https://www.davidzeleny.net/anadat-r/doku.php/en:pcoa_nmds
# there are a lot of zeroes, so PCA and similar analyses are not ideal
# https://websites.pmc.ucsc.edu/~mclapham/Rtips/ordination.htm#:~:text=PCoA%20is%20a%20distance%2Dbased,for%20more%20information%20on%20distances).

PCoA.bd<-capscale(transp_bd~1,distance="bray")
PCoA.bd
summary(PCoA.bd)

plot(PCoA.bd, type = "n")
points(PCoA.bd, display = "sites", pch = 2, col="#56B4E9", cex = 1.5)
text(PCoA.bd, display = "sites", cex=1, col="black", adj = c(1.2,0))

# the code below is in progress, contains various options for plotting
# nicer plot with ggvegan:
library("ggvegan")
library("ggplot2")
library("viridis")
library("tidyverse")
autoplot(PCoA.bd, layers = "sites", arrows = FALSE)
#auplot <- autoplot(PCoA.bd, layers = c("sites","species"), arrows = FALSE)
#auplot
aufort <- fortify(PCoA.bd, display = c("sites","species"))
aufort
# extract site coords (points)
fort_sites <- aufort %>% 
  filter(Score == "sites")

# extract species coords
fort_species <- aufort %>% 
  filter(Score == "species")

# we can plot with coloring by amount of Bd (not available in new data), 
# one way is to arbitrarily assign categories, e.g.
# L (low) - 0-10/sample, M (medium) - 11-30, H (high) - >30
#fort_sites$bdabund <- c("L","M","H","L","L","H","M","H","H")
# or make continuous based on bdabundance
# fort_sites$Bd_abundance <- transp_bd$sperangium

# colorblind palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# or just add after geom point line to use viridis
#  scale_color_viridis(discrete = TRUE, option = "D") +
library(wesanderson)

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/


#fortify_plot <- ggplot() +
  #geom_point(data = fort_species, aes(x = MDS1, y = MDS2)) +
  #geom_point(data = fort_sites, aes(x = MDS1, y = MDS2), size = 3, shape = 2) +
#  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#  geom_point(data = fort_sites, aes(x = MDS1, y = MDS2, color = Bd_abundance),size = 7, shape = 17) +
  #scale_colour_manual(values=wes_palette("Darjeeling1", n = 3)) +
#  scale_color_viridis(discrete = TRUE, option = "D") +
#  geom_text(data = fort_sites, size = 6, aes(x = MDS1, y = MDS2, label = Label), nudge_x = 0.1, nudge_y = 0.1) +
#  theme(legend.key.size = unit(2, 'cm'), legend.title = element_text(size=18),
#        legend.text = element_text(size=20),
#        axis.text.x = element_text(size = 15),
#        axis.text.y = element_text(size = 15),
#        axis.title = element_text(size = 20)) + 
#     labs(x = "MDS1",
#       y = "MDS2")
#fortify_plot



