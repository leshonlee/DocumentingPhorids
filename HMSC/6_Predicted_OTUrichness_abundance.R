setwd()

library(Hmsc)
library(ggplot2)
library(abind)

load("models_thin_1000_samples_250_chains_4.Rdata")

#Predicted mOTU occurrence and log(abundance) as a function of the continuous covariates %forest/woodland cover and mean annual temperature.
#One PDF generated per sampling season, among which that for late summer is the basis for Supplementary figure S9.
#Only those mOTUs are illustrated for which a directional trend was predicted with > 95% posterior support. 
#The numbers in the legend represent the percentages of mOTUs, of the total modelled, that were predicted to have a
#positive (red) or negative (blue) response to the focal covariate. These are marginal predictions, where values for 
#all covariates apart from the focal covariate and sampling season (late summer) were fixed at their mean value in the dataset. 

models = models[1:2]

library(khroma)
discrete_rainbow <- colour("smooth rainbow")

seasons = c("1_Late_Spring", "2_Midsummer", "3_Late_Summer", "4_Offseason")   

for(n in 1:4)
{
          pdf(paste0("mOTUresponses_Gradients_", seasons[n], "_June2023.pdf"), width = 8, height = 10)    
          
          par(mfrow = c(3, 2))
          covariates = c("ForestWood", "bio1")
          
          for(i in 1:length(covariates))
          {
            covariate = covariates[i]
            xlabel = switch(i, "Forest or Wood (%)", "Mean Annual Temperature (C)")
            for(j in 1:2)
            {
            m = switch(j, models[[1]], models[[2]])
            
            plottitle = switch(j, "OTU prevalence", "OTU abundance")
            
            ylimits = switch(j, c(0,1), c(0,8))
            ylabel = switch(j, "Predicted prevalence", "Predicted log. abundance (when present)")
            
            
            Gradient = constructGradient(m,focalVariable = covariate, non.focalVariables = 1)
            Gradient$XDataNew$Time.Period = seasons[n]
            
            predY = abind(predict(m, Gradient=Gradient, expected = TRUE), along = 3)
            
            if(j == 1)
            {
            EpredY = apply(predY, c(1,2), mean)
            }
            if(j == 2)
            {
              #predY = exp(predY)
              EpredY = apply(predY, c(1,2), mean)
            }
            
            selgradient = Gradient$XDataNew[,1]
            
            Pr = vector()
            
            for(sp in 1:ncol(EpredY))
              {
              Pr[sp] = mean(apply(predY[2:20, sp, ], 2, max)>predY[1,sp, ]) 
            }
            
            sel_taxa = which(Pr>0.975|Pr<0.025)
            Pr = Pr[which(Pr>0.975|Pr<0.025)]
            legendtext = colnames(EpredY)[sel_taxa]
            plotpredY = data.frame(EpredY[,sel_taxa])
            colours = rep("blue",length(sel_taxa))
            colours[which(Pr>0.975)] = "red"
            
            for(selsp in 1:ncol(plotpredY))
            {
                if(selsp == 1)
                  plot(selgradient, plotpredY[,selsp], xlab = xlabel, main = plottitle, ylab = ylabel, ylim = ylimits, col = colours[selsp], type = "l", lwd = 1, lty = 1)
                if(selsp > 1)
                  points(selgradient, plotpredY[,selsp], type = "l", col = colours[selsp], lwd = 1, lty = 1)
                if(selsp == 1 & length(sel_taxa>0))
                {legend("topright", legend = c(paste0(round(length(which(Pr>0.975))/ncol(EpredY)*100, 2), "%"), paste0(round(length(which(Pr<0.025))/ncol(EpredY)*100, 2), "%")), horiz = TRUE, fill = c("red", "blue"), cex=0.8)}
              }
            }
          }
          dev.off()
}


