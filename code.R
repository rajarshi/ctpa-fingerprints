library(rcdk)
library(ggplot2)
library(parallel)

mols <- load.molecules('data/chembl.smi')
fps <- lapply(mols, get.fingerprint)

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
  xlab('Bit Position')+ylab('Normalized Frequency')


