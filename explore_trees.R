library(ape)
library(phangorn)




# Pama-Nyungan languages
all.trees <- read.nexus("PamaNyungan/prunedPN.trees")

# Sino-Tibetan languages
 all.trees <- read.nexus("SinoTibetan/SinoTibetanSubset.nex")


# Consensus tree ----------------------------------------------------------

# Consensus topology
cons <- consensus(all.trees, p = .5, rooted = T)
cons$node.label <- round(cons$node.label * 100, 0)
plot(cons, show.node.label=TRUE)

# Only display posterior probabilities for uncertain splits
cons$node.label[cons$node.label >95] <- NA
plot(cons, show.node.label=TRUE)


# Consensus topology + branch lengths. This might take a couple of minutes
ct <- consensus.edges(all.trees, consensus.tree = cons, rooted = TRUE)
ct$root.edge <- .15

plot(ct, show.node.label=FALSE, root.edge = TRUE, no.margin = TRUE, 
     font = 1, label.offset = .05)
nodelabels(c(NA, ct$node.label[-1]), frame="none", adj = c(1.15,1.25), cex = .9)


# Densitree ----------------------------------------------------------------

densiTree(all.trees, consensus = cons, alpha = .005, font = 1, 
          label.offset = .01, cex = 1, scale.bar = FALSE)


maxBT <- max(getAges(tt_scaled))
label <- rev(pretty(c(maxBT, 0)))
axis(side = 1, at = seq(0, 1.0, length.out = length(label)), 
     labels = label, line = -2.5, cex = .9)




# MCC tree --------------------------------------------------------------------------------------------------------

mcct <- maxCladeCred(all.trees)
mcct$node.label <- round(mcct$node.label * 100, 0)
mcct$node.label[mcct$node.label > 95] = NA
mcct$root.edge <- .15

plot(mcct, show.node.label=TRUE, root.edge = TRUE, no.margin = TRUE, 
     font = 1, label.offset = .05)



