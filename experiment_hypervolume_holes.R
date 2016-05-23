#######################################################################
### Erosion of phylogenetic and functional diversity in Amphibians  ###
### Author: "Brunno F Oliveira                                      ###
### Universidade Federal do Rio Grande do Norte, Natal, RN, Brazil. ###
### Correspondense: brunno.oliveira@me.com                          ###
#######################################################################

# Here I test what is the best bandwidth for calculating hypervolumes 
# using Amphibians traits

library(hypervolume)
library(raster)

rm(list=ls())
gc()

#setwd("/media/brunno/FAT/Chap3/R_code2")
setwd("F:/Chap3/R_code2")

load("PD+FD_erosion_getdata.RData")

traits2<-read.csv('traits2.csv')
rownames(traits2)<-traits2$X

# simulate bandwidths ranging from 0.01 to 2
bd<-seq(0,2,0.05)
bd<-bd[-1]

js <- which(Rich>10)[1]

tipnames <- names(occr)[which(occr[js,]==1)]
com <- traits2[tipnames,2:4]

Hole <- data.frame(matrix(NA,10,length(bd)))
colnames(Hole) <- bd

for(i in 1:length(bd)){
  for(j in 1:10){ cat('---------------------------------')
    write.table(c(i,j), 'run.txt')
    # MAKE HYPERVOLUME
    tmp<-hypervolume(com, bandwidth = .07, verbose = F, warnings=F, quantile = 0.05)
    # GET CONVEX EXPECTATION
    hv_convex <- expectation_convex(tmp, check_memory=F)
    # DETECT HOLES (low npoints for fast execution)
    features_annulus <- hypervolume_holes(hv_obs=tmp,hv_exp=hv_convex,set_check_memory=F)
    if(any(is.null(features_annulus),get_volume(features_annulus) < 0.001)){
      Hole[j,i] <- NA 
    }
    else{
      # CLEAN UP RESULTS
      features_segmented <- hypervolume_segment(features_annulus)
      features_segmented_pruned <- hypervolume_prune(features_segmented, minvol=0.001)
      # PLOT RETAINED HOLE(S)
      # plot(hypervolume_join(tmp, features_segmented))
      # CALCULATE THE HOLE FRACTION
      if(is.numeric(get_volume(features_segmented_pruned))){
        Hole[j,i] <- get_volume(features_segmented_pruned)/get_volume(tmp)     }
      else{
        Hole[j,i] <- NA 
      }
    }
  }
}

save.image("experiment_hypervolume.RData")


boxplot(Hole)
  
# estimated bandwidths using the method of silverman
bandw <- NA
for(i in 1:nrow(occr)){ cat(i, "from", nrow(occr), "\r")
  tipnames <- names(occr)[which(occr[i,]==1)]
  bandw[i] <- estimate_bandwidth(traits2[tipnames,2:4])
}

hist(bandw)

plot(log(Rich),bandw)



# correlate bandwidth with hypervolum and convex hulls
for(i in 1:length(bd)){ cat(i, "from", length(bd), "\r")
tmp<-hypervolume(scale(traits2[,2:4]), bandwidth = .5, warnings=F, quantile = 0.05)
plot(tmp, pairplot=F)
hyper[i] <- tmp@Volume
chull[i] <- convhulln(tmp@RandomUniformPointsThresholded, "FA")$vol
}

# Test whether hypervolume and convex are associated
plot(chull,hyper)
cor(chull,hyper) # This variables are strongly associated. 

# Test whether hypervolume and convex are associated with the bandwidth
# A larger bandwidth will potentially create a larger trait space.
plot(bd,hyper)
cor(bd,hyper) 
plot(bd,chull)
cor(bd,chull) # This variables are strongly associated and show identical results.

# Now lets test whether these indices are sensitive to the presence of holes.
# Communities with larger variation in trait values are expected to create
# larger trait spaces, and consequently high FD values. However, these FD values may be
# biased by outliers, and not representing a true high FD. 

tmp<-hypervolume(func_space[,2:4], bandwidth = .5, quantile = 0.05, warnings=F)
tmp@Volume
plot(tmp)
  
# manipulate SD in traits
# large SD should result in higher chull, but not hypervolume
hyper<-NA
chull<-NA

bw<-seq(0.01,2,0.01)

for(i in 1:length(bw)){
  func_space <- data.frame(species = paste("species", 1:30), 
                           habitat = rnorm(30, sd=bw[i]), 
                           trophic = rnorm(30, sd=bw[i]), 
                           life_history = rnorm(30, sd=bw[i])
  )
  
  tmp<-hypervolume(func_space[,2:4], bandwidth = 0.5, warnings=F, quantile = 0.05)
  hyper[i] <- tmp@Volume
  chull[i] <- convhulln(func_space[,2:4], "FA")$vol
}
  
plot(bw,hyper)
plot(bw,chull)
plot(chull,hyper)
