library(phytools)
library(corrplot)
library(RColorBrewer)


all.trees <- read.nexus("PamaNyungan/PNY10_285.trees")

# We have Pama-Nyungan trees inferred from lexical data.
# We will map Functional Load data on the same trees.

d <- read.table("PamaNyungan/FL_VlengthC_normalised.tsv", header=T)

# Prune the trees to keep only leaves for which we have Functional Load data.
ll <- d$tip_label
idx <- which(tr$tip.label %in% ll)


# For now, work on a single tree. 
ntrees <- length(all.trees)
tr <- keep.tip(all.trees[[ntrees]], idx)

tip.sort <- match(tr$tip.label, d$tip_label)


# Compute phylogenetic signal of functional load using Blomberg's K
phylosig(tr, d$FL_Vlength[tip.sort], method = "K", test = TRUE)

# Plot continuous trait over the tree
trait <- d$FL_Vlength
names(trait) <- d$tip_label
trait <- trait[tip.sort]
contMap(tr,  trait)


# Reconstruct ancestral value
fastAnc(tr, trait, CI=T)




# Test for correlated evolution between traits
ci =<-1:6 # index of columns of interest
covm <- phyl.vcv(as.matrix(d[tip.sort, ci]), vcv(tr), 1)

corrplot(cov2cor(covm$R))

## Average over trees
covm <- array(NA, c(length(ci), length(ci), ntrees))
for(i in 1:ntrees){
  tr <- all.trees[[i]]
  idx <- which(tr$tip.label %in% ll)
  tr <- keep.tip(tr, idx)
  tip.sort <- match(tr$tip.label, d$tip_label)
  covm[,, i] <- cov2cor(phyl.vcv(as.matrix(d[tip.sort, ci]), vcv(tr), 1)$R)
}

mean.covm <- apply(covm, 1:2, mean)
rownames(mean.covm) = colnames(mean.covm) =  colnames(d)[ci]
corrplot(mean.covm)

# if we don't take into account the phylogeny, we vastly overestimate the correlation
covm.nophylo <- cor(d[ , ci])
corrplot(covm.nophylo)

# Consider a binary trait: does a language exhibit Vowel length?
# We can reconstruct ancestral values for the binary trait on a known tree
bintrait <- (trait>0)

fit <- fitMk(tr, bintrait, model="ARD", pi="fitzjohn")

anc <- ancr(fit)

plot(anc, piecol=brewer.pal(3, "Dark2"))
