### These functions to be used for Testa and Dworkin 2016 paper at DGE. 
## the following functions are used primarily for resampling purposes 


###############
#
# Figure 2
#
##############

#returns value for SSD in wt lines
SSD_wt_func <- function()
{
  dat_temp <- temp.wings[sample(nrow(temp.wings),replace=T),]
  size_lm <- lm(Centroid.Size ~ GENOTYPE*SEX, data=dat_temp)
  as.numeric((coef(size_lm)[1]/(coef(size_lm)[1] + coef(size_lm)[3])) -1)
}

#returns value for SSD in mutant lines
SSD_mutant_func <- function()
{
  dat_temp <- temp.wings[sample(nrow(temp.wings),replace=T),]
  size_lm <- lm(Centroid.Size ~ GENOTYPE*SEX, data=dat_temp)
  as.numeric(((coef(size_lm)[1] + coef(size_lm)[2])/(sum(coef(size_lm)))) -1 )
}

#standard Procrustes Distance formula
PD <- function(x)
{
  sqrt(t(x) %*% x)
}

#returns value of SShD for wt lines
SShD_wt_func <- function()
{
  dat_temp <- temp.master[sample(nrow(temp.master),replace=T),]
  shape_PCs <- as.matrix(dat_temp[,6:19])
  PC_lm <- lm(shape_PCs ~ geno*sex, data=dat_temp)
  wt_f_shape <- coef(PC_lm)[1,]
  wt_m_shape <- coef(PC_lm)[1,] + coef(PC_lm)[3,]
  wt_SShD <- wt_f_shape - wt_m_shape
  PD(wt_SShD) # This corresponds to the distance you see in the geomorph distance for wild type SShD.
}

#returns value of SShD for mutant lines
SShD_mutant_func <- function()
{
  dat_temp <- temp.master[sample(nrow(temp.master),replace=T),]
  shape_PCs <- as.matrix(dat_temp[,6:19])
  PC_lm <- lm(shape_PCs ~ geno*sex, data=dat_temp)
  p_f_shape <- coef(PC_lm)[1,] + coef(PC_lm)[2,]
  p_m_shape <- coef(PC_lm)[1,] + coef(PC_lm)[2,] + coef(PC_lm)[3,] + coef(PC_lm)[4,]
  mt_SShD <- p_f_shape - p_m_shape
  PD(mt_SShD)
}



###############
#
# Figure 3
#
##############


##vector correlation between Sam and Ore wt backgrounds
boot_SShD_VC_wt_test <- function()
{
  repeat
  {
    sample_dat <- dat_together[sample(nrow(dat_together), replace=T),]
    
    if (length(subset(sample_dat, background=="Ore")[,1]) >= 2 & length(subset(sample_dat, background=="Sam")[,1]) >= 2 )
    {
      break
    }
  }
  shape <- as.matrix(sample_dat[,1:14])
  fitboot_PC<- lm(shape ~ sex*background, data=sample_dat)
  estimates <- coef(fitboot_PC)
  
  # wt female shape
  ore_f_shape <- estimates[1,]
  
  # wt male shape
  ore_m_shape <- estimates[1,] + estimates[2,]
  
  # mutant female shape
  sam_f_shape <- estimates[1,] + estimates[3,]
  
  # mutant male shape 
  sam_m_shape <- estimates[1,] + estimates[2,] + estimates[3,] + estimates[4,]
  
  ang.vec.abs((ore_f_shape - ore_m_shape), (sam_f_shape - sam_m_shape))
}

#returns (absolute) values for vector correlation and angle of vector
ang.vec.abs <- function(vec1, vec2)
{
  vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
  vec.angle <- acos(vec.cor)*(180/pi)
  return(c(vector.cor=vec.cor, vec.angle=vec.angle))
}
comment(ang.vec.abs) <- c(" This computes both the vector correlation, and angle, between two vectors.", " to compare to the Pearson correlation coefficient make sure to center and standardize vectors", "set it up to compute the absolute values of the vector correlation")

#function used for testing vetor correlations (deprecated)
boot_SShD_VC_test <- function()
{
  repeat
  {
    sample_dat <- dat_together[sample(nrow(dat_together), replace=T),]
    
    if (length(subset(sample_dat, sex=="M")[,1]) >= 2 & length(subset(sample_dat, sex=="F")[,1]) >= 2 )
    {
      break
    }
  }
  shape <- as.matrix(sample_dat[,1:14])
  fitboot_PC<- lm(shape ~ sex*genotype, data=sample_dat) 
  estimates <- coef(fitboot_PC)
  
  # wt female shape
  wt_f_shape <- estimates[1,]
  
  # wt male shape
  wt_m_shape <- estimates[1,] + estimates[2,]
  
  # mutant female shape
  p_f_shape <- estimates[1,] + estimates[3,]
  
  # mutant male shape 
  p_m_shape <- estimates[1,] + estimates[2,] + estimates[3,] + estimates[4,]
  
  ang.vec.abs((wt_f_shape - wt_m_shape), (p_f_shape - p_m_shape))
}

###############
#
# Figure 4
#
##############


#returns values for allometric slopes and associated confidence intervals
allometry_func <- function()
{
  dat_temp <- temp.master[sample(nrow(temp.master), replace=T),]
  shape <- as.matrix(dat_temp[,6:19])
  #   size <- log(dat_temp[,1])   #csize
  size <- log(dat_temp[,1])   #csize
  sex <- dat_temp[,2]
  genotype <- dat_temp[,5]
  
  allometry.model <- lm(shape ~ size)
  Beta <- allometry.model$coefficients[2,] 
  allometry.shape <- ShapeScore(t(Beta), shape)
  shape <- allometry.shape
  allometry <- lm(shape ~ size)
  coef(allometry)[[2]]
}

