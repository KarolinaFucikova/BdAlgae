library(vegan)
#library(ggdendro)
library("ggplot2")

Bd <- read.csv(file = "data/Bd_data_clean.csv", row.names =1)
Treat <- read.csv(file="data/Bd_treatment.csv")
# these data have abundances of species calculated from 3 passes (transects) 
# across each of two slides per sample
# first column has species name, second has group (family/order) affiliation

# taxa can be excluded if needed
# clean_bd <- subset(Bd, select = -speciesname)

# first, combine data from species in the same groups (e.g., merge all Selenastraceae)
Bd_grouped <- aggregate(Bd[2:29], by=list(Group=Bd$Group), FUN=sum)
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

# calculate the distances based on abundances
# euclidean distances are appropriate for abundance data
eucl <- vegdist(transp_bd,method="euclidean")

# now we'll make a hierarchical cluster object, which will be the basis of the dendrogram
bd_hc <-hclust(eucl)
# and finally the dendrogram
plot(bd_hc)

# labeling - first converting fully to dendrogram
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning#plot.dendrogram-function
hcd <- as.dendrogram(bd_hc)
plot(hcd, type = "rectangle", ylab = "Height")
# dend_data <- dendro_data(hcd, type = "rectangle")


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
# or make continuous based on bdabundance
# fort_sites$Bd_abundance <- transp_bd$sperangium

# or better by treatment
fort_sites$treat <- Treat$treatment
# could add columns with other data about tanks
fort_sites$day<-Treat$day
fort_sites$tank<-Treat$tank

# colorblind palette with grey:
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# or just add after geom point line to use viridis
#  scale_color_viridis(discrete = TRUE, option = "D") +
library(wesanderson)

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/


fortify_plot <- ggplot() +
  #geom_point(data = fort_species, aes(x = MDS1, y = MDS2)) +
  geom_point(data = fort_sites, aes(x = MDS1, y = MDS2), size = 3, shape = 2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = fort_sites, aes(x = MDS1, y = MDS2, color = treat),size = 7, shape = 17) +
  scale_colour_manual(values=wes_palette("Moonrise2", n = 3)) +
  #scale_color_viridis(discrete = TRUE, option = "D") +
  geom_text(data = fort_sites, size = 6, aes(x = MDS1, y = MDS2, label = day), nudge_x = 0.05, nudge_y = 0.05) +
  theme(legend.key.size = unit(2, 'cm'), legend.title = element_text(size=18),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20)) + 
     labs(x = "MDS1",
       y = "MDS2")
fortify_plot

# other analyses that take environmental variables into account
# DO, pH, conductivity in particular
# turn groups into row names
# first removing the tank number as variable, and creating a modified Treat df
Treat1 <- Treat[,-4] 
rownames(Treat1) <- Treat1[,1]
# the first column remains as group and needs to be removed prior to analyses
Treat1[,1] <- NULL
# the ojbect Treat remains unchanged, this is important for later removing rows


#perform the CCA. Type ccamodel(species_data ~., environmental_data)
ccamodel <- cca(transp_bd ~., Treat1)
#view the results of the CCA with summary()
summary(ccamodel)
#test for the significance of the CCA model with anova(ccamodel)
anova(ccamodel)
#test for the significance of the CCA model axes with anova(ccamodel, by = "axis")
anova(ccamodel, by = "axis")

#plotting fancy
autoplot(ccamodel)
# https://github.com/gavinsimpson/ggvegan/issues/9
aufort <- fortify(ccamodel)
aufort
# extract site coords (points)
fort_sites <- aufort %>% 
  filter(Score == "sites")
#fort_sites$latitude <- environmental_variables_cca$lat

# extract species coords
fort_species <- aufort %>% 
  filter(Score == "species")

# extract arrows coords
fort_arrows <- aufort %>% 
  filter(Score == "biplot")
# remove treatment from fort_arrows
fort_arr1 <- fort_arrows[-c(2,3),]

# extract centroid coords for treatments to be separate
fort_arr2 <- aufort %>% 
  filter(Score=="centroids")

# colorblind palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# or just add after geom point line to use viridis
#  scale_color_viridis(discrete = TRUE, option = "D") +
library(wesanderson)
# pal <- wes_palette(name = "Darjeeling1", type = "continuous")

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
# https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels/#change-the-legend-font-size-color-and-face

fortify_plot <- ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = fort_arr1,
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_segment(data = fort_arr2,
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_point(data = fort_sites, aes(x = CCA1, y = CCA2, color = Treat1$treatment),size = 5, shape = 17) +
  scale_colour_manual(values=wes_palette("Moonrise2", n = 3)) +
  #geom_point(data = fort_species, aes(x = CCA1, y = CCA2),size = 3, shape = 1) +
  #  scale_color_viridis(discrete = TRUE, option = "D") +
  #geom_text(data = fort_sites, size = 4, aes(x = CCA1, y = CCA2, label = Label), nudge_x = 0.1, nudge_y = 0.1) +
  geom_text(data = fort_arr1, size = 4,
            mapping = aes(label = Label, x = CCA1, y = CCA2, hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  geom_text(data = fort_arr2, size = 5,
            mapping = aes(label = c("Bd", "No Bd", "No tadpoles"), x = CCA1, y = CCA2, hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
    theme(legend.key.size = unit(2, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20)) + 
  labs(x = "CCA1",
       y = "CCA2")
fortify_plot

# for colors in the Moonrise2 scale color =c("#798E87", "#C27D38", "#CCC591")
# using Bd prevalence as continuous predictor variable
# only available for day 18, so different subset of data
# day 18 only may be best anyway

# excluding the no tadpole data to only compare Bd infected and uninfected tanks
# first use the Treat data to find samples that are in the "no tadpole" treatment
no_tad <- subset(Treat, subset = treatment!="NoTad")
# no_tad_list <- c(no_tad$Sample)
transp_bd1 <- transp_bd[c(no_tad$Sample),]

PCoA.tadonly<-capscale(transp_bd1~1,distance="bray")
PCoA.tadonly
summary(PCoA.tadonly)

autoplot(PCoA.tadonly, layers = "sites", arrows = FALSE)
aufort.tad <- fortify(PCoA.tadonly, display = c("sites","species"))
aufort.tad
# extract site coords (points)
fort_sites.tad <- aufort.tad %>% 
  filter(Score == "sites")

# extract species coords
fort_species.tad <- aufort.tad %>% 
  filter(Score == "species")

# plot with coloring by treatment
fort_sites.tad$treat <- no_tad$treatment
# could add columns with other data about tanks
fort_sites.tad$day<-no_tad$day
fort_sites.tad$tank<-no_tad$tank

fortify_plot.tad <- ggplot() +
  #geom_point(data = fort_species, aes(x = MDS1, y = MDS2)) +
  geom_point(data = fort_sites.tad, aes(x = MDS1, y = MDS2), size = 3, shape = 2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = fort_sites.tad, aes(x = MDS1, y = MDS2, color = treat),size = 7, shape = 17) +
  scale_colour_manual(values=wes_palette("Moonrise2", n = 3)) +
  #scale_color_viridis(discrete = TRUE, option = "D") +
  geom_text(data = fort_sites.tad, size = 6, aes(x = MDS1, y = MDS2, label = day), nudge_x = 0.05, nudge_y = 0.05) +
  theme(legend.key.size = unit(2, 'cm'), legend.title = element_text(size=18),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20)) + 
  labs(x = "MDS1",
       y = "MDS2")
fortify_plot.tad

# other analyses that take environmental variables into account
# DO, pH, conductivity in particular
# turn groups into row names
# first removing the tank number as variable, and creating a modified no_tad
no_tad1 <- no_tad[,-4] 
rownames(no_tad1) <- no_tad1[,1]
# the first column remains as group and needs to be removed prior to analyses
no_tad1[,1] <- NULL


#perform the CCA. Type ccamodel(species_data ~., environmental_data)
ccamodel.tad <- cca(transp_bd1 ~., no_tad1)
#view the results of the CCA with summary()
summary(ccamodel.tad)
#test for the significance of the CCA model with anova(ccamodel)
anova(ccamodel.tad)
#test for the significance of the CCA model axes with anova(ccamodel, by = "axis")
anova(ccamodel.tad, by = "axis")

#plotting fancy
autoplot(ccamodel.tad)
# https://github.com/gavinsimpson/ggvegan/issues/9
aufort.tad1 <- fortify(ccamodel.tad)
aufort.tad1
# extract site coords (points)
fort_sites.tad1 <- aufort.tad1 %>% 
  filter(Score == "sites")

# extract species coords
fort_species.tad1 <- aufort.tad1 %>% 
  filter(Score == "species")

# extract arrows coords
fort_arrows.tad1 <- aufort.tad1 %>% 
  filter(Score == "biplot")
# remove treatment from fort_arrows
fort_arr1.tad1 <- fort_arrows.tad1[-c(2,3),]

# extract centroid coords for treatments to be separate
fort_arr2.tad1 <- aufort.tad1 %>% 
  filter(Score=="centroids")

fortify_plot.tad1 <- ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = fort_arr1.tad1,
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_segment(data = fort_arr2.tad1,
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  #geom_point(data = fort_species.tad1, aes(x = CCA1, y = CCA2),size = 3, shape = 16) +
  geom_point(data = fort_sites.tad1, aes(x = CCA1, y = CCA2, color = no_tad1$treatment),size = 5, shape = 17) +
  scale_colour_manual(values=wes_palette("Moonrise2", n = 3)) +
  #geom_point(data = fort_species, aes(x = CCA1, y = CCA2),size = 3, shape = 1) +
  #  scale_color_viridis(discrete = TRUE, option = "D") +
  #geom_text(data = fort_sites, size = 4, aes(x = CCA1, y = CCA2, label = Label), nudge_x = 0.1, nudge_y = 0.1) +
  geom_text(data = fort_arr1.tad1, size = 4,
            mapping = aes(label = Label, x = CCA1, y = CCA2, hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  geom_text(data = fort_arr2.tad1, size = 5,
            mapping = aes(label = c("Bd", "No Bd"), x = CCA1, y = CCA2, hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  theme(legend.key.size = unit(2, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20)) + 
  labs(x = "CCA1",
       y = "CCA2")
fortify_plot.tad1


