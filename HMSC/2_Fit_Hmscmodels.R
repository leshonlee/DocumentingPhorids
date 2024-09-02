setwd() 

library(Hmsc)

#Load unfitted Hmsc model created by script 1.
load(file = "unfitted_models.Rdata")

#warning, this script is extremely slow and heavy to run
nChains = 4
samples = 250
thin = 1000

for(n in 1:4)
{
  m = models[[n]]
  m = sampleMcmc(m, samples = samples, thin=thin,
                 adaptNf=rep(ceiling(0.4*samples*thin),m$nr),
                 transient = ceiling(0.5*samples*thin),
                 nChains = nChains, nParallel = nChains)
  models[[n]] = m
}

filename_out = paste("models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")


save(models, modelnames, file=filename_out)