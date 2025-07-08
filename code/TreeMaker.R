# load required packages
library(ape)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(stringr)
library(gridExtra)
library(grid)
library(ggtext)

# Set complimentary color palettes for host-of-origin and clade affinities

hostcolors <- c(
  "steelblue1",
  "mediumpurple4", 
  "darkorange1", 
  "green4", 
  "mediumpurple1", 
  "lightsalmon", 
  "#E31A1C", 
  "#20B2AA", 
  "orchid", 
  "olivedrab3",
  "darkolivegreen",
  "lightgray"
  )

# assign colors to host names
# or read host list from metadata file
names(hostcolors) <- scan(text = "P.brutia P.echinata P.elliottii P.nigra P.palustris P.ponderosa P.radiata P.strobus P.sylvestris P.taeda P.thunbergii Unknown", what = "")

cladecolors <- c(
  "#20B2AA",  
  "royalblue",
  "firebrick2",
  "mediumpurple1",
  "goldenrod4", 
  "darkolivegreen", 
  "olivedrab3", 
  "palegreen",
  "darkorange1",
  "mediumpurple1",
  "mediumpurple4",
  "purple4",
  "azure4",
  "#005F5F",
  "deeppink1",
  "navyblue",
  "steelblue1",
  "orchid1", 
  "royalblue",
  "palevioletred",
  "rosybrown1",
  "maroon" 
)

# assign colors to clade names
# or read clade list from metadata file
names(cladecolors) <- scan(text = "N.USA SpainN SpainS S.USA", what = "")

# Read in tree data
  Tree1 <- read.tree("~/RAxML_bipartitions.LaJul9")

Tree1 <- root(Tree1, outgroup = "CBS133791")


# Read in metadata
tipdata <- read.table("~/LaAllTreeStrains.txt")
colnames(tipdata) <- c("strain", "cladeID", "hostID")
tipdata$strain = factor(tipdata$strain)
tipdata$hostID <- factor(tipdata$hostID)
tipdata$cladeID <- factor(tipdata$cladeID)

# create genus labels
genusLabels <- paste0("*", names(hostcolors), "*")

# generate metadata for fruit
fruitdata <- tipdata
fruitdata$clade <- fruitdata$cladeID

# Plot the tree
p1 <- ggtree(Tree1, layout = "rectangular")

p2 <- p1 %<+% tipdata + 
  geom_tiplab(aes(color = hostID), size = 3, align = T, line.size = 0.25, offset = 0.08, show.legend = F) +
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) ==100), size = 1.5, show.legend=F) +
  geom_tippoint(aes(color = hostID), stroke = 0, alpha=0.6, size = 3, position = position_nudge(x = 0.01)) +
  scale_color_manual(values=hostcolors,
                     name ="Host Genus",
                     labels = genusLabels,
                     guide=guide_legend(keywidth=0.8,
                                        keyheight=0.8,
                                        ncol=1,
                                        order=2,
                                        override.aes=list(size=4,alpha=0.6)))

p3 <- p2 +
  geom_fruit(data=fruitdata, geom=geom_tile, mapping = aes(y=strain, fill=clade), width = 0.06, offset = 0.3, alpha = 0.6) +
  scale_fill_manual(values=cladecolors,
                    name = "Clade ID",
                    na.translate = F,
                    guide=guide_legend(keywidth=0.7,
                                       keyheight=0.4,
                                       ncol=1,
                                       order=1
                    )) +
  
  geom_treescale(x = 0, y = -5, width = 0.1, color = "black", linesize = 0.7, fontsize = 0) +
        annotate("text", x = 0.05, y = -3.8, label = "0.1", size = 5) +
  
  theme(plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"), 
        legend.position = "right", 
        legend.box = "vertical", 
        legend.direction = "vertical",
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0.5),
        legend.spacing.y = unit(18, "cm"),
 #       legend.margin = margin(r = 20, l = -10),
        title = element_text(size = rel(1.5)), 
        legend.text = element_markdown(size = rel(1)),
        legend.key.size = unit(0.8, "cm"))


print(p3)

