## Simple demos of R embedding of HyperHMM


### simple demos of R embedding of HyperHMM

# source code for inference
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hyperhmm-r.cpp")

# source code for plots
library(ggpubr)
source("hypercube-plots.R")

# read in cross-sectional data and return a matrix
cube.read.crosssectional = function(fname) {
  data.raw = readLines(fname)
  data.mat = do.call(rbind, lapply(strsplit(data.raw, ""), as.numeric))
  return(data.mat)
}

# read in longitudinal data and return a list of two matrices
cube.read.longitudinal = function(fname) {
  data.list = list()
  data.raw = read.table(fname, header=FALSE, colClasses = "character")
  data.list$from = do.call(rbind, lapply(strsplit(data.raw[,1], ""), as.numeric))
  data.list$to = do.call(rbind, lapply(strsplit(data.raw[,2], ""), as.numeric))
  return(data.list)
}


### ovarian cancer dataset
cgh.mat = cube.read.crosssectional("Data/Cross_sectional_data.txt")
fit.cgh = HyperHMM(cgh.mat)
plot.cgh = plot.standard(fit.cgh)

### TB dataset
tb.list = cube.read.longitudinal("Data//phylogeny_data.txt")
fit.tb = HyperHMM(tb.list$to, initialstates=tb.list$from, nboot=1)
plot.tb = plot.standard(fit.tb)
#png("newplot1.png",width=1200, height=700)
# put these more specialised plots together
#ggarrange(plot.cgh, plot.tb,ncol = 2,nrow = 2, vjust = 15)
#dev.off()

plot1=plot.hypercube.flux(fit.cgh)
plot2=plot.hypercube.flux(fit.tb)

#png("cross_phylo_plot.png",width=1200, height=750)

ggarrange(plot.cgh, plot.tb,plot1, plot2,ncol = 2,nrow = 2,labels = c("A", "B", "C","D"),size=10)
#dev.off()






## Africa
cgh.afr = cube.read.crosssectional("Data/Africa.txt")
fit.afr = HyperHMM(cgh.afr)
plot.afr = plot.standard(fit.afr)
flux_afr=plot.hypercube.flux(fit.afr)
## Asia
cgh.asia = cube.read.crosssectional("Data/Asia.txt")
fit.asia = HyperHMM(cgh.asia)
plot.asia = plot.standard(fit.asia)
flux_asia=plot.hypercube.flux(fit.asia)

## Europe
cgh.europe = cube.read.crosssectional("Data/Europe.txt")
fit.europe = HyperHMM(cgh.europe)
plot.europe = plot.standard(fit.europe)
flux_europe=plot.hypercube.flux(fit.europe)

## N.america= n.a
cgh.n.a = cube.read.crosssectional("Data/North_America.txt")
fit.n.a = HyperHMM(cgh.n.a)
plot.n.a = plot.standard(fit.n.a)
flux_n.a=plot.hypercube.flux(fit.n.a)

## Oceania

cgh.oceania = cube.read.crosssectional("Data/Oceania.txt")
fit.oceania = HyperHMM(cgh.oceania)
plot.oceania = plot.standard(fit.oceania)
flux_oceania=plot.hypercube.flux(fit.oceania)

## South America =s.a
cgh.s.a = cube.read.crosssectional("Data/South_America.txt")
fit.s.a = HyperHMM(cgh.s.a)
plot.s.a = plot.standard(fit.s.a)
flux_s.a=plot.hypercube.flux(fit.s.a)


#png("Continents_plot.png",width=1200, height=750)
ggarrange(plot.afr,plot.asia,plot.europe,plot.n.a,plot.oceania,plot.s.a,ncol = 3,nrow = 3,labels = c("A", "B", "C","D","E","F"),size = 10)
#dev.off()

#png("Continents_flux.png",width=1200, height=750)

ggarrange(flux_afr, flux_asia,flux_europe, flux_n.a,flux_oceania,flux_s.a,ncol = 3,nrow =3,labels = c("A", "B", "C","D","E","F"),size = 10)
#dev.off()

