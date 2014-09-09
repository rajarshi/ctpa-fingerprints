library(rcdk)
library(ggplot2)
library(parallel)

#########################################
## Reading fingerprints
#########################################
fps <- fp.read('data/cdk.fp', size=1024, lf=cdk.lf)

fps.100 <- fps[1:100]
fpm <- fp.to.matrix(fps.100)

#########################################
## Random fingerprints
#########################################
nfp <- 300
sizes <- c(64, 128, 512, 1024, 4096, 8192)
times <- sapply(sizes, function(size) {
  fps <- lapply(1:nfp, function(i) random.fingerprint(size, size * 0.35))
  system.time(junk <- fp.sim.matrix(fps))[3]
})

tmp <- data.frame(size=sizes, time=times)
ggplot(tmp, aes(x=size, y=time))+
  geom_point()+xlab("Fingerprint Length")+ylab("Time (s)")

densities <- c(0.1, 0.25, 0.5, 0.75, 0.95)
times <- sapply(densities, function(density) {
  fps <- lapply(1:nfp, function(i) random.fingerprint(1024, 1024 * density))
  system.time(junk <- fp.sim.matrix(fps))[3]
})

tmp <- data.frame(density=densities, time=times)
ggplot(tmp, aes(x=density, y=time))+
  geom_point()+xlab("Bit Density")+ylab("Time (s)")

#########################################
## Compare similarity metrics
#########################################
fps <- fp.read('data/cdk.fp', size=881, lf=cdk.lf, header=TRUE)[1:500]
s.tanimoto <- fp.sim.matrix(fps, method='tanimoto')
s.dice <- fp.sim.matrix(fps, method='dice')
d <- rbind(data.frame(method='Tanimoto', s=as.numeric(s.tanimoto)),
           data.frame(method='Dice', s=as.numeric(s.dice)))
ggplot(d, aes(x=s, fill=method))+geom_density(alpha=0.25)+
  scale_fill_discrete(name='Metric')+
  xlab("Similarity")+theme(legend.position = c(0.85, 0.8))

#########################################
## Predicting with fingerprints
#########################################
sol <- read.csv('data/solubility.csv', header=TRUE)
fps <- fp.read('data/solubility.maccs', header=FALSE, size=166, lf=function(line) {
  toks <- strsplit(line, " ")[[1]]
  title <- toks[1]
  bits <- as.numeric(toks[2:length(toks)])
  list(title, bits, list())
})
## Extract fingerprint for which we have a label
common <- which( sapply(fps, function(x) x@name) %in% sol$sid )
fps <- fps[common]
## Order the fingerprints & data
sol <- sol[order(sol$sid),]
fps <- fps[order(sapply(fps, function(x) as.integer(x@name)))]
## Make X matrix
fpm <- fp.to.matrix(fps)
## Model!
library(randomForest)
m1 <- randomForest(x=fpm, y=as.factor(sol$label))
## force equal class size
m2 <- randomForest(x=fpm, y=as.factor(sol$label), sampsize=854)

ggplot(data.frame(table(sol$label)),
       aes(x=Var1, y=Freq))+
  geom_bar(stat='identity', colour='black', fill='beige')+
  xlab("Solubility Class")+ylab("Frequency")

#########################################
## Clustering fingerprints
#########################################
fps <- fp.read('data/cdk.fp', size=881, lf=cdk.lf, header=TRUE)[1:500]
sims <- fp.sim.matrix(fps)
dmat <- as.dist(1-sims)
clus <- hclust(dmat)
par(mar=c(1,4,1,1))
plot(clus, label=FALSE, xlab='', main='')

## compare clusterins from different fingerprints
library(rcdk)
library(dendextend)
mols <- load.molecules('data/chembl.smi')
fps <- lapply(c('pubchem', 'extended', 'graph', 'maccs'), function(type) {
  lapply(mols[1:300], get.fingerprint, type)
})
dms <- lapply(fps, function(x) as.dist(1-fp.sim.matrix(x)))
cls <- lapply(dms, function(x) as.dendrogram(hclust(x)))
tanglegram(cls[[1]], cls[[2]], main_left='Pubchem 881', main_right='CDK Ext 1024',
           k_branches=5)

csim <- do.call(rbind, mclapply(cls, function(x) {
  sapply(cls, function(y) {
    cor_cophenetic(x,y)
  })
}))
rownames(csim) <- c('Pubchem', 'CDK Extended', 'CDK Graph', 'MACCS')
colnames(csim) <- rownames(csim)
image(csim)

#########################################
## Bit spectrum
#########################################
fps <- fp.read('data/cdk.fp', size=881, lf=cdk.lf, header=TRUE)

##library(rcdk)
##mols <- load.molecules('data/chembl.smi')
##fps <- lapply(mols, get.fingerprint, 'pubchem')

bs <- bit.spectrum(fps)
d <- data.frame(x=1:length(bs), y=bs)
ggplot(d, aes(x=x,y=y))+geom_line()+
  xlab('Bit Position')+ylab('Normalized Frequency')

## make two subsets and generate bit spectra for each
fps1 <- fps[ 1:(length(fps)*0.3) ]
fps2 <- fps[ (length(fps)*0.3 + 1):length(fps) ]

bs1 <- bit.spectrum(fps1)
bs2 <- bit.spectrum(fps2)

## display a difference plot along with individual spectra
bsdiff <- bs1-bs2
d <- data.frame(x=1:length(bs), y=bsdiff)
ggplot(d, aes(x=x,y=y))+geom_line()+
  xlab('Bit Position')+ylab('Normalized Frequency')+
  ylim(c(-1,1))

## Compare bit spectra for soluble & insoluble compounds
sol <- read.csv('data/solubility.csv', header=TRUE)
fps <- fp.read('data/solubility.maccs', header=FALSE, size=166, lf=function(line) {
  toks <- strsplit(line, " ")[[1]]
  title <- toks[1]
  bits <- as.numeric(toks[2:length(toks)])
  list(title, bits, list())
})
common <- which( sapply(fps, function(x) x@name) %in% sol$sid )
fps <- fps[common]
sol <- sol[order(sol$sid),]
fps <- fps[order(sapply(fps, function(x) as.integer(x@name)))]

sol.idx <- which(sol$label == 'high')
insol.idx <- which(sol$label != 'high')

sol.bs <- bit.spectrum(fps[sol.idx])
insol.bs <- bit.spectrum(fps[insol.idx])
bsdiff <- sol.bs - insol.bs
d <- data.frame(x=1:length(sol.bs), y=bsdiff)
ggplot(d, aes(x=x,y=y))+geom_line()+
  xlab('Bit Position')+ylab(expression(paste(Delta,' Normalized Frequency')))+
  ylim(c(-1,1))
