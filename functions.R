##### -- Calculate functional distinctness (Fdist) -- #####
### Calculate Fdist given a functional space and a set of species
calculate_Fdist <- function(Fspace = NA, 
                            sp_names = NULL){
  ### Load packages
  require(nabor)
  ### Extract desired species subset from the total functional space
  if (!is.null(sp_names)) Fspace <- Fspace[Fspace$species %in% sp_names, ] else Fspace <- Fspace
  if (nrow(Fspace) > 0){
    ### Generate empty output table
    Fdist_table <- data.frame(species = Fspace$species, 
                              Fdist_overall = NA, 
                              Fdist_Axis1 = NA, 
                              Fdist_Axis2 = NA, 
                              Fdist_Axis3 = NA)
    ### Calculate FDist for focal species[i] compared to all others
    runs <- nrow(Fdist_table)
    for (i in 1:runs){ cat(i, "of", runs, "\r")
      Fdist_table$Fdist_overall[i] <- median(knn(Fspace[, 2:4], query = Fspace[i, 2:4], k = nrow(Fspace))$nn.dists[1, ]) # overall
      Fdist_table$Fdist_Axis1[i] <- median(knn(Fspace[, 2:4]$Axis1, query = Fspace[i, 2:4], k = nrow(Fspace))$nn.dists[1, ]) # along the first dimension only
      Fdist_table$Fdist_Axis2[i] <- median(knn(Fspace[, 2:4]$Axis2, query = Fspace[i, 2:4], k = nrow(Fspace))$nn.dists[1, ]) # along the second dimension only
      Fdist_table$Fdist_Axis3[i] <- median(knn(Fspace[, 2:4]$Axis3, query = Fspace[i, 2:4], k = nrow(Fspace))$nn.dists[1, ]) # along the third dimension only
    }        
  } else Fdist_table <- NULL
  return(Fdist_table)
}

#### -- Generate 3D plot of a defined functional space -- ####
func_space_plot <- function(Fspace = NA, sp_names = NULL, Fdist = NULL, ...)
{
  ### Load packages
  require(plot3D)
  ### Extract desired species subset from the total functional space
  if (!is.null(sp_names)) Fspace <- Fspace[Fspace$species %in% sp_names, ] else 
    Fspace <- Fspace
  ### Generate 3D plot with the desired species
  scatter3D(x = Fspace$Axis1, y = Fspace$Axis2, z = Fspace$Axis3, 
            xlab = "Axis1", ylab = "Axis2", zlab = "Axis3", 
            col = grey(0.4), pch = 3, cex = 0.7)
  ### Optionally, label the most and least functionally-distinct species
  if (!is.null(Fdist)){
    points3D(x = Fspace$Axis1, 
             y = Fspace$Axis2, 
             z = Fspace$Axis3,
             colvar = Fdist$Fdist_overall,
             col = colorRampPalette(c("tomato", "deepskyblue"))(6),
             pch = 3, cex = 0.8, add = TRUE, ...)        
  }
}

#### function to random sample n species in a community and calc pd
random.drop <- function(n, com, phy){
  # species in community i
  tipnames <- names(com)[which(com[i,]==1)]
  # random sample n species from community i
  sp.tmp <- sample(1:length(tipnames), n, replace = F)
  # these are the randomly sampled species
  tipnames <- tipnames[sp.tmp]
  # calculate PD
  trx <- drop.tip(tree, tree$tip.label[-na.omit(match(tipnames, tree$tip.label))])
  return(sum(trx$edge.length))
}

#### function to random sample n tips and calc pd
random.drop.phy <- function(n, phy){
  # random sample X species from that pool
  sp.tmp <- sample(1:length(phy$tip.label), n, replace = F)
  # null.dataset
  tipnames <- phy$tip.label[-sp.tmp]
  trx <- drop.tip(phy, tipnames)
  return(sum(trx$edge.length))
}

#### function to random sample n species and calc FD
random.fd <- function(n, com, td){
  # species in community i
  tipnames <- names(com)[which(com[i,]==1)]
  # random sample n species from community i
  sp.tmp <- sample(1:length(tipnames), n, replace = F)
  # these are the randomly sampled species
  tipnames <- tipnames[sp.tmp]
  # calculate FD
  tmp <- hypervolume(td[tipnames,], bandwidth = bandw, verbose = F, warnings=F, quantile = 0.05)
  tmp <- get_volume(tmp)
}

# function to get the phylogenetic most distinctive species
distinphy <- function(phy,n,m=TRUE){
  n <- n
  names <- NA
  if(m==TRUE){
    for(i in 1:n){
      names[i] <- names(sort(colSums(cophenetic(phy)), decreasing = T)[1])
      phy <- drop.tip(phy, names[i])
    }
  }
  else {
    for(i in 1:n){
      names[i] <- names(sort(colSums(cophenetic(phy)))[1])
      phy <- drop.tip(phy, names[i])
    }
  }  
  names
}


# function to sample the n phylogenetic most distinctive species and calc PD
min.max.PD <- function(phy,n){
  n <- n
  names <- NA
  phymin <- phy
  phymax <- phy
  ### Get the most distictive species and remove
  ### This will return the min PD possible
  for(k in 1:n){
    #cat(k, "from", n, '\r')
    names[k] <- names(sort(colSums(cophenetic(phymin)), decreasing = T)[1])
    phymin <- drop.tip(phymin, names[k])
  }
  minPD=sum(phymin$edge.length)
  for(k in 1:n){
    #cat(k, "from", n, '\r')
    ### Get the less distictive species and remove
    ### This will return the min PD possible
    names[k] <- names(sort(colSums(cophenetic(phymax)))[1])
    phymax <- drop.tip(phymax, names[k])
  }
  maxPD=sum(phymax$edge.length)
  c(minPD=minPD,maxPD=maxPD)
}


## first get the node numbers of the tips
nodes<-sapply(trx$tip.label,function(x,y) which(y==x),y=trx$tip.label)
## then get the edge lengths for those nodes
edge.lengths<-setNames(trx$edge.length[sapply(nodes, function(x,y) which(y==x),y=trx$edge[,2])],names(nodes))


# function to sample the n functionaly most distinctive species and calc FD
min.max.FD <- function(td,n,m=TRUE){
  n <- n
  td <- calculate_Fdist(na.omit(td))
  ### Get the most distictive species and remove
  ### This will return the min FD possible
  minFD <- mean(td$Fdist_overall[order(test$Fdist_overall, decreasing = T)][1:n])
  ### Get the less distictive species and remove
  ### This will return the max FD possible
  maxFD <- mean(td$Fdist_overall[order(test$Fdist_overall)][1:n])
  ### Return the results
  c(minFD=minFD, maxFD=maxFD)
}